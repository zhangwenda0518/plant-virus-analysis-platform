# -*- coding: utf-8 -*-
"""NCBI 提交准备（自 app.py 拆出）。

unified_metadata.csv 一张表驱动 GenBank + BioSample 提交产物。"""
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid

from flask import (Blueprint, abort, jsonify, render_template, request,
                   send_file, send_from_directory)

from vp.config import DIRS, PLATFORM_ROOT, db_path, engine_cmd
from vp.utils import (TaskLogger, check_path, fmt_size, run_cmd, safe_open,
                      safe_remove)
from vp.web.state import cfg, tool_runs_root as _tool_runs_root
from vp.web.tasks import tm

_log = logging.getLogger('vp.web.submit')

bp = Blueprint('submit', __name__)


@bp.route('/submit')
def page_submit():
    return render_template('submit.html')


@bp.route('/api/submit/tables')
def api_submit_tables():
    from vp.ncbi_submit import store
    return jsonify(store.list_tables())


@bp.route('/api/submit/samples')
def api_submit_samples():
    """内置示例数据集（对应原 GUI Samples 菜单）。"""
    from vp.ncbi_submit import store
    return jsonify(store.list_samples())


@bp.route('/api/submit/create', methods=['POST'])
def api_submit_create():
    from vp.ncbi_submit import store
    body = request.get_json(force=True) or {}
    name = str(body.get('name') or '').strip()
    sample = str(body.get('sample') or ('demo' if body.get('demo') else ''))
    try:
        store.create_table(name, sample=sample)
    except FileExistsError:
        abort(400, f'提交项目已存在: {name}')
    except ValueError as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'name': name})


@bp.route('/api/submit/copy', methods=['POST'])
def api_submit_copy():
    """项目另存为（原 GUI Save As）。body: {src, dst}"""
    from vp.ncbi_submit import store
    body = request.get_json(force=True) or {}
    src = str(body.get('src') or '').strip()
    dst = str(body.get('dst') or '').strip()
    try:
        store.copy_table(src, dst)
    except FileNotFoundError:
        abort(404, f'源项目不存在: {src}')
    except FileExistsError:
        abort(400, f'提交项目已存在: {dst}')
    except ValueError as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'name': dst})


@bp.route('/api/submit/upload', methods=['POST'])
def api_submit_upload():
    from vp.ncbi_submit import store
    name = str(request.form.get('name') or '').strip()
    f = request.files.get('file')
    if not f or not f.filename:
        abort(400, '缺少上传文件')
    import tempfile
    fd, tmp = tempfile.mkstemp(suffix=os.path.splitext(f.filename)[1])
    try:
        with os.fdopen(fd, 'wb') as out:
            f.save(out)
        try:
            _, rows = store.import_table(name, tmp)
        except FileExistsError:
            abort(400, f'提交项目已存在: {name}')
        except ValueError as e:
            abort(400, str(e))
        except Exception as e:
            abort(400, f'解析失败（需要 CSV/TSV/Excel）: {e}')
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass
    return jsonify({'ok': True, 'name': name, 'rows': rows})


def _submit_store(name):
    from vp.ncbi_submit import store
    try:
        store.table_dir(name)
    except ValueError as e:
        abort(400, str(e))
    return store


@bp.route('/api/submit/table/<name>')
def api_submit_table(name):
    store = _submit_store(name)
    try:
        return jsonify(store.table_payload(name))
    except FileNotFoundError:
        abort(404, f'提交项目不存在: {name}')


@bp.route('/api/submit/table/<name>/save', methods=['POST'])
def api_submit_table_save(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    cols = body.get('columns') or []
    raw = body.get('rows') or []
    rows = [dict(zip(cols, r)) for r in raw]
    try:
        n = store.save_table(name, rows)
    except ValueError as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'rows': n})


@bp.route('/api/submit/table/<name>/rows', methods=['POST'])
def api_submit_table_rows(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        if body.get('add'):
            n = store.add_rows(name, body['add'])
        elif body.get('delete') is not None:
            n = store.delete_rows(name, body['delete'])
        else:
            abort(400, '需要 add 或 delete')
    except (ValueError, IndexError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'rows': n})


@bp.route('/api/submit/table/<name>/fill', methods=['POST'])
def api_submit_table_fill(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    column = str(body.get('column') or '')
    value = str(body.get('value') or '')
    if not column or not value:
        abort(400, '需要 column 和 value')
    try:
        count = store.batch_fill(name, column, value, old_value=body.get('old'))
    except KeyError as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'count': count})


@bp.route('/api/submit/table/<name>/export', methods=['POST'])
def api_submit_table_export(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        files = store.export_files(
            name,
            assembler=str(body.get('assembler') or 'SPAdes;4.3.0;metaviral'),
            sequencer=str(body.get('sequencer') or 'Illumina NovaSeq 6000'),
            enrichment=str(body.get('enrichment') or 'rRNA depletion'))
    except Exception as e:
        abort(400, f'生成失败: {e}')
    return jsonify({'ok': True, 'files': files})


@bp.route('/api/submit/table/<name>/validate', methods=['POST'])
def api_submit_table_validate(name):
    store = _submit_store(name)
    return jsonify({'issues': store.validate_table(name)})


@bp.route('/api/submit/table/<name>/fasta', methods=['POST'])
def api_submit_table_fasta(name):
    """从平台内 FASTA（contigs 运行 / 样品结果 / 任意平台内文件）
    按表行 sequence_name 提取序列 → 项目目录 sequences.fsa + 一致性报告。"""
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    src = str(body.get('fasta') or '').strip()
    if not src:
        abort(400, '请提供序列 FASTA 路径')
    try:
        p = check_path(src if os.path.isabs(src)
                       else os.path.join(PLATFORM_ROOT, src), must_exist=True)
        rep = store.export_submission_fasta(
            name, p, min_len=int(body.get('min_len', 200) or 200))
    except (ValueError, KeyError, FileNotFoundError, OSError) as e:
        abort(400, f'FASTA 导出失败: {e}')
    return jsonify({'ok': True, **rep})


@bp.route('/api/submit/table/<name>/import_run', methods=['POST'])
def api_submit_table_import_run(name):
    """body: {run} → 从 tool_runs/<run>/virus_classification.tsv 追加表行。"""
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        res = store.import_from_run(name, str(body.get('run') or ''))
    except (ValueError, KeyError, FileNotFoundError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, **res})


@bp.route('/api/submit/contig_runs')
def api_submit_contig_runs():
    from vp.ncbi_submit import store
    return jsonify(store.list_contig_runs())


@bp.route('/api/submit/orf_runs')
def api_submit_orf_runs():
    """有 CDS 注释产物（04_orf/04b）的 orf/orfa 运行列表（关联注释用）。"""
    from vp.ncbi_submit import store
    return jsonify(store.list_orf_runs())


@bp.route('/api/submit/table/<name>/link_orf', methods=['POST'])
def api_submit_table_link_orf(name):
    """把独立跑的 orf/orfa 注释运行关联到提交项目（按 contig 名对接）。"""
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        res = store.link_orf_run(name, str(body.get('run') or ''))
    except (ValueError, FileNotFoundError, KeyError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, **res})


@bp.route('/api/submit/table/<name>/tbl', methods=['POST'])
def api_submit_table_tbl(name):
    """从已关联 orf 运行的 CDS 注释生成 featuretable.tbl（GenBank 特征表）。"""
    store = _submit_store(name)
    try:
        res = store.generate_feature_tbl(name)
    except ValueError as e:
        abort(400, str(e))
    return jsonify({'ok': True, **res})


@bp.route('/api/submit/table/<name>/source_fasta')
def api_submit_table_source_fasta(name):
    """自动推断提交序列 FASTA 路径（viral_contigs / contigs.filtered）。"""
    store = _submit_store(name)
    try:
        p = store.infer_source_fasta(name)
    except (ValueError, OSError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'fasta': p})


@bp.route('/api/submit/table/<name>/package', methods=['POST'])
def api_submit_table_package(name):
    """项目内全部提交产物 → submission_package.zip（一键下载用）。"""
    store = _submit_store(name)
    try:
        path, n = store.package_zip(name)
    except (ValueError, OSError) as e:
        abort(400, f'打包失败: {e}')
    return jsonify({'ok': True, 'files': n,
                    'download': f'/submissions/{name}/submission_package.zip'})


@bp.route('/api/submit/table/<name>/taxonomy_check', methods=['POST'])
def api_submit_table_taxonomy_check(name):
    """organism 逐个查 NCBI Taxonomy（联网）；未收录物种会提交被拒。"""
    store = _submit_store(name)
    try:
        rows = store.taxonomy_check(name)
    except KeyError as e:
        abort(400, str(e))
    except RuntimeError as e:
        abort(502, f'NCBI 查询失败: {e}')
    return jsonify({'ok': True, 'rows': rows,
                    'missing': [r['organism'] for r in rows if not r['found']]})


@bp.route('/api/submit/table/<name>/sbt', methods=['POST'])
def api_submit_table_sbt(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        tpl, path = store.generate_sbt(
            body.get('fields') or {},
            extra_authors=body.get('extra_authors') or [],
            title=body.get('title') or '', out_name=name)
    except ValueError as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'path': path,
                    'authors': 1 + len(body.get('extra_authors') or [])})


@bp.route('/api/submit/table/<name>/sqn', methods=['POST'])
def api_submit_table_sqn(name):
    """本地 .sqn 生成（suvtk features→comments→table2asn / MIUVIG）。

    body: {author:{...}, use_features:bool, miuvig:{...},
           assembler, sequencer}
    author 缺省时回退到 store.generate_sbt 需要的必填字段（last/first/affil/
    city/country/email）；use_features=True 走 suvtk BFVD 链，False 用平台
    自家 featuretable.tbl（需先 /tbl 生成）。"""
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    from vp import suvtk_submit as ss
    import pandas as pd

    # 1) 关联的 orf 运行目录（CDS 坐标来源 + dominant_genome_type 推断）
    link = store._read_link(name, 'orf')
    if not link or not link.get('run'):
        abort(400, '尚未关联 orf 注释运行——先点"关联注释"')
    run_dir = store._orf_run_dir(link['run'])

    # 2) 表行 → build_sqn 的 rows（过滤占位符 sequence_name）
    df = store.load_table(name)
    rows = []
    for _, r in df.iterrows():
        sn = str(r.get('sequence_name') or '').strip()
        if not sn or store.is_placeholder(sn):
            continue
        rows.append({c: ('' if pd.isna(v) else v) for c, v in r.items()})
    if not rows:
        abort(400, '表内没有有效 sequence_name 行')

    # 3) 作者字段
    author = body.get('author') or {}

    # 4) 其余参数
    use_features = bool(body.get('use_features', True))
    # .sqn 及中间文件生成到提交项目目录（submissions/<name>/），与④文件预览/
    # 编辑一致——用户编辑 .sbt/.src/.cmt 后可重新生成 .sqn 并复用编辑。
    out_dir = store.table_dir(name)
    try:
        res = ss.build_sqn(
            run_dir, rows=rows, author=author,
            miuvig=body.get('miuvig'),
            assembler=body.get('assembler'),
            sequencer=body.get('sequencer'),
            use_features=use_features,
            tbl_path=os.path.join(out_dir, 'featuretable.tbl')
            if not use_features else None,
            out_dir=out_dir,
            log=lambda m: _log.info('[sqn] %s', m))
    except (RuntimeError, ValueError, FileNotFoundError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, **res})


@bp.route('/api/submit/table/<name>/sqn_status')
def api_submit_table_sqn_status(name):
    """已有 .sqn 产物状态（含 .val 校验摘要），供 submit 页展示。

    .sqn 生成在提交项目目录（submissions/<name>/），与④文件预览一致。"""
    store = _submit_store(name)
    d = store.table_dir(name)
    sqn = os.path.join(d, 'sqn.sqn')
    val = os.path.join(d, 'sqn.val')
    if not os.path.isfile(sqn):
        return jsonify({'ok': True, 'exists': False, 'dir': d})
    errs = warns = infos = 0
    if os.path.isfile(val):
        with open(val, encoding='utf-8', errors='replace') as f:
            for line in f:
                s = line.strip()
                if s.startswith('Error'):
                    errs += 1
                elif s.startswith('Warning'):
                    warns += 1
                elif s:
                    infos += 1
    return jsonify({'ok': True, 'exists': True,
                    'sqn': sqn, 'val': val, 'dir': d,
                    'stats': {'error': errs, 'warning': warns, 'info': infos}})


@bp.route('/api/submit/table/<name>/biosample', methods=['POST'])
def api_submit_table_biosample(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        path, rows, skipped = store.export_biosample_tsv(
            name, skip_placeholders=not bool(body.get('include_placeholders')))
    except (ValueError, KeyError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'path': path, 'rows': rows, 'skipped': skipped})


@bp.route('/api/submit/table/<name>/file')
def api_submit_table_file(name):
    store = _submit_store(name)
    fname = request.args.get('name') or ''
    try:
        return jsonify({'content': store.read_file(name, fname)})
    except (ValueError, FileNotFoundError) as e:
        abort(400, str(e))


@bp.route('/api/submit/table/<name>/file_save', methods=['POST'])
def api_submit_table_file_save(name):
    """预览编辑保存（原 GUI Preview 标签的 Edit/Save）。body: {name, content}"""
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    fname = str(body.get('name') or '')
    try:
        store.write_file(name, fname, str(body.get('content') or ''))
    except (ValueError, FileNotFoundError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'file': fname})


@bp.route('/api/submit/table/<name>/delete', methods=['POST'])
def api_submit_table_delete(name):
    import shutil
    store = _submit_store(name)
    d = store.table_dir(name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    return jsonify({'ok': True})


@bp.route('/submissions/<name>/<path:fname>')
def submissions_download(name, fname):
    """提交产物下载（限项目目录内）。"""
    from vp.ncbi_submit import store
    try:
        p = store.read_file_path(name, fname)
    except (ValueError, FileNotFoundError) as e:
        abort(404, str(e))
    return send_file(p, as_attachment=True)

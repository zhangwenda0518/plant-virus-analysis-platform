# -*- coding: utf-8 -*-
"""参考序列与比较分析 API（自 app.py 拆出）。

NCBI Entrez 下载、ICTV 级联选参、GenBank 集合管理、共线性比较、
MSA / 进化树 / SDT 数据接口。"""
import json
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
from vp.web.common import _safe_sample
from vp.web.state import cfg, tool_runs_root as _tool_runs_root
from vp.web.tasks import tm

bp = Blueprint('refs', __name__)


@bp.route('/api/ncbi/collections')
def api_ncbi_collections():
    from vp.ncbi_download import list_collections
    return jsonify([c for c in list_collections()
                    if not c.get('name', '').startswith('_')])


@bp.route('/api/ictv/cascade')
def api_ictv_cascade():
    """级联选参下拉数据：?Realm=..&Phylum=..（已选等级作为过滤）→
    {selected, ranks:[{col, zh, options:[{name,n}]}]}。植物病毒口径。"""
    from vp import ictv_db
    if not ictv_db.available():
        abort(400, 'ICTV 库未初始化（python main.py ictv-update）')
    levels = {col: (request.args.get(col) or '').strip()
              for col, _zh in ictv_db.RANK_COLS}
    levels = {k: v for k, v in levels.items() if v}
    try:
        opts = ictv_db.rank_options(levels, plant_only=True)
    except Exception as e:
        abort(400, f'ICTV taxa 读取失败: {e}')
    return jsonify({
        'selected': levels,
        'ranks': [{'col': col, 'zh': zh, 'options': opts.get(col, [])}
                  for col, zh in ictv_db.RANK_COLS]})


def _ictv_levels(body):
    """从请求体取分类级选择 {列名: 值}（非空）+ 最深层（用于日志与定位）。"""
    from vp import ictv_db
    levels = {}
    for col, _zh in ictv_db.RANK_COLS:
        v = str((body.get(col) or '')).strip()
        if v:
            levels[col] = v
    if not levels:
        abort(400, '请至少选择一个分类等级（界/门/纲/目/科/属/种）')
    deepest_col = [c for c, _zh in ictv_db.RANK_COLS if c in levels][-1]
    return levels, deepest_col, levels[deepest_col]


@bp.route('/api/ictv/preview', methods=['POST'])
def api_ictv_preview():
    """body: {<分类级列名>: 值, ..., limit, genome} → 选参预览（不下载）。"""
    body = request.get_json(force=True) or {}
    levels, dcol, dval = _ictv_levels(body)
    genome = body.get('genome') or 'complete'
    limit = max(1, min(int(body.get('limit', 200) or 200), 2000))
    per_genus = max(0, min(int(body.get('per_genus', 0) or 0), 50))
    per_species = max(0, min(int(body.get('per_species', 0) or 0), 20))
    from vp import ictv_db
    try:
        rows, total = ictv_db.select_refs(ranks=levels, plant_only=True,
                                          genome=genome,
                                          limit=max(limit * 20, 2000))
    except Exception as e:
        abort(400, f'选参失败: {e}')
    # 抽样上限先于 limit 生效（先每属/每种封顶，再取前 limit 条）
    rows, dropped = ictv_db.cap_per_rank(rows, per_genus, per_species)
    rows = rows[:limit]
    cap_note = []
    if per_genus:
        cap_note.append(f'每属≤{per_genus}')
    if per_species:
        cap_note.append(f'每种≤{per_species}')
    return jsonify({'scope': f'{dcol}={dval}', 'total': total,
                    'n_shown': len(rows),
                    'sampling': '、'.join(cap_note) or '不限',
                    'dropped_by_cap': dropped,
                    'rows': [{'acc': r.get('Genbank', ''),
                              'source': r.get('Source', ''),
                              'species': r.get('Species', ''),
                              'coverage': r.get('Genome_Coverage', ''),
                              'genome': r.get('Genome', ''),
                              'host': r.get('Host_Source', '')}
                             for r in rows]})


@bp.route('/api/ictv/download', methods=['POST'])
def api_ictv_download():
    """body: {rank, name, coll, limit, genome} → 后台任务：
    ICTV 选参（本地优先）→ 下载整科/整属 GenBank 集合（含巡检清单）
    → 同步导出 FASTA 参考集（databases/ncbi_refs/<coll>/refs.fa，
    供样品流程建树追加参考）。"""
    body = request.get_json(force=True) or {}
    levels, dcol, dval = _ictv_levels(body)
    coll = (body.get('coll') or '').strip()
    genome = body.get('genome') or 'complete'
    if not coll:
        abort(400, '参数不完整（请填集合名）')
    limit = max(1, min(int(body.get('limit', 200) or 200), 2000))
    per_genus = max(0, min(int(body.get('per_genus', 0) or 0), 50))
    per_species = max(0, min(int(body.get('per_species', 0) or 0), 20))

    def job(log, prog, cancel):
        import shutil
        logger = TaskLogger(callback=log)
        from vp import ictv_db
        from vp.gb_collection import (download_gb_collection, gb_collection_dir,
                                      collection_records)
        from vp.ncbi_download import collection_dir as fasta_coll_dir
        from vp.utils import write_fasta_record

        prog('ictv', 0.05, f'ICTV 选参：{dcol} = {dval}（{len(levels)} 级过滤）')
        rows, total = ictv_db.select_refs(ranks=levels, plant_only=True,
                                          genome=genome,
                                          limit=max(limit * 20, 2000))
        rows, _dropped = ictv_db.cap_per_rank(rows, per_genus, per_species)
        rows = rows[:limit]
        if not rows:
            raise RuntimeError(f'ICTV {dcol} "{dval}" 无匹配参考'
                               '（换一个分类或改 genome=any）')
        accs = [r['Genbank'] for r in rows]
        logger.log(f'ICTV 命中 {total} 条，取前 {len(accs)} 条'
                   f'（{"Complete genome" if genome == "complete" else "全部"}）')

        # 重名集合替换式重建（与建库口径一致，不累积重复记录）
        shutil.rmtree(gb_collection_dir(coll), ignore_errors=True)
        prog('ictv', 0.15, f'下载 GenBank 集合 {coll}（{len(accs)} 条 accession）')
        res = download_gb_collection(coll, accessions=','.join(accs),
                                     max_records=len(accs),
                                     logger=logger, prog=prog, cancel=cancel)

        # 同步 FASTA 参考集（样品流程「NCBI 参考集合」可直接填 coll）
        prog('ictv', 0.85, '同步 FASTA 参考集（样品建树用）')
        fdir = fasta_coll_dir(coll)
        os.makedirs(fdir, exist_ok=True)
        fa = os.path.join(fdir, 'refs.fa')
        n_fa = 0
        with safe_open(fa, 'wt') as f:
            for acc, _org, seq in collection_records(coll):
                if seq:
                    write_fasta_record(f, acc, seq)
                    n_fa += 1
        logger.log(f'FASTA 参考集: {fa}（{n_fa} 条）')
        prog('ictv', 1.0, f'完成：GenBank {res.get("n_records", n_fa)} 条 + '
                    f'FASTA {n_fa} 条')
        logger.close()
        return {'collection': coll, 'scope': f'{dcol}={dval}',
                'levels': levels, 'n_records': n_fa, 'n_fasta': n_fa}

    tid = tm.start(cfg.tr(f'ICTV下载 {coll}', f'ICTV download {coll}'), job,
                   weight='light')
    return jsonify({'task': tid})


@bp.route('/api/ncbi/search', methods=['POST'])
def api_ncbi_search():
    """body: {term, db, limit} → {count, rows:[{acc,title,organism,...}]}。"""
    body = request.get_json(force=True)
    term = (body.get('term') or '').strip()
    db = body.get('db') or 'nucleotide'
    if not term or db not in ('nucleotide', 'protein'):
        abort(400, '参数不完整（term 必填，db=nucleotide/protein）')
    from vp.ncbi_download import search_preview
    try:
        return jsonify(search_preview(term, db=db,
                                      limit=int(body.get('limit', 50))))
    except RuntimeError as e:
        abort(502, str(e))


@bp.route('/api/ncbi/download', methods=['POST'])
def api_ncbi_download():
    """body: {term, name, db, max} → 后台任务下载集合。"""
    body = request.get_json(force=True)
    term = (body.get('term') or '').strip()
    name = (body.get('name') or '').strip()
    db = body.get('db') or 'nucleotide'
    if not term or not name or db not in ('nucleotide', 'protein'):
        abort(400, '参数不完整（term/name 必填，db=nucleotide/protein）')
    max_records = max(1, min(int(body.get('max', 100) or 100), 10000))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.ncbi_download import download_collection
        res = download_collection(term, name, db=db, max_records=max_records,
                                  logger=logger, prog=prog, cancel=cancel)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'NCBI下载 {name}', f'NCBI download {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


@bp.route('/api/gb/collections')
def api_gb_collections():
    """集合列表（只读概览：读已生成的清单/警告文件，不重新解析 .gb）。

    下划线前缀的集合（_内部测试/_归档）不在列表展示。
    """
    from vp.gb_collection import list_gb_collections, read_manifest
    out = []
    for c in list_gb_collections():
        if c.get('name', '').startswith('_'):
            continue
        try:
            st = read_manifest(c['name'])
            c['records'] = [{'acc': r['acc'], 'organism': r['organism'],
                             'length': r['length'], 'cds': r['cds_count']}
                            for r in st['records']]
            c['warnings'] = st['warnings']
        except (OSError, FileNotFoundError, ValueError):
            c['records'], c['warnings'] = [], []
        c['has_compare'] = os.path.isfile(
            os.path.join(c['dir'], 'compare', 'compare.html'))
        c['has_lovis4u'] = os.path.isfile(
            os.path.join(c['dir'], 'compare', 'lovis4u.pdf'))
        c['has_phylo'] = os.path.isfile(
            os.path.join(c['dir'], 'phylo', 'aln.fasta'))
        try:
            from vp.gb_collection import list_extract_files
            ex = list_extract_files(c['name'])
        except Exception:
            ex = None
        c['has_extract'] = bool(ex)
        c['n_extract_genes'] = (len([x for x in (ex or [])
                                     if x['kind'].endswith('_gene')
                                     and x['kind'].startswith('cds')]) or
                                len({x['path'].split('/')[-1]
                                     for x in (ex or [])
                                     if x['kind'] == 'cds_gene'}))
        out.append(c)
    return jsonify(out)


@bp.route('/api/gb/extract', methods=['POST'])
def api_gb_extract():
    """body: {name} → 后台任务：集合 .gb → extract/ 分类分目录提取
    （genome.fa / CDS.fa / PEP.fa / CDS/<基因>.fa / PEP/<基因>.fa / genes.tsv）。"""
    body = request.get_json(force=True)
    name = (body.get('name') or '').strip()
    if not name:
        abort(400, '参数不完整（name 必填）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gb_collection import extract_collection_features
        res = extract_collection_features(name, logger=logger, prog=prog)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'CDS/PEP 提取 {name}', f'Extract CDS/PEP {name}'),
                   job, weight='light')
    return jsonify({'task': tid})


@bp.route('/api/gb/extract_files')
def api_gb_extract_files():
    """某集合的提取产物清单（供「参考序列获取」卡展示与送下游）。"""
    name = (request.args.get('name') or '').strip()
    if not name:
        abort(400, '参数不完整（name 必填）')
    from vp.gb_collection import list_extract_files
    items = list_extract_files(name)
    if items is None:
        abort(404, '该集合尚未提取（先点「🧬 提取 CDS/PEP」）')
    return jsonify(items)


@bp.route('/api/gb/inspect', methods=['POST'])
def api_gb_inspect():
    """body: {name} → 后台任务全量巡检（重解析全部 .gb，重建清单与警告）。"""
    body = request.get_json(force=True)
    name = (body.get('name') or '').strip()
    if not name:
        abort(400, '参数不完整（name 必填）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gb_collection import inspect_collection
        st = inspect_collection(name, logger=logger)
        logger.close()
        return {'name': name, 'n_records': len(st['records']),
                'n_warnings': len(st['warnings'])}

    tid = tm.start(cfg.tr(f'GenBank巡检 {name}', f'GenBank inspect {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


@bp.route('/api/gb/phylo', methods=['POST'])
def api_gb_phylo():
    """body: {name, tree_tool, molecule, gene} → 后台任务：集合比对 + 建树。

    molecule: genome（全基因组，产物 phylo/，结果中心 gb:<name> 展示）/
    cds / pep（基因级：按 gene 关键词从 .gb 提取 CDS/蛋白，
    产物 gene_trees/<gene>_<molecule>/）。
    """
    body = request.get_json(force=True)
    name = (body.get('name') or '').strip()
    tree_tool = body.get('tree_tool') or 'fasttree'
    molecule = (body.get('molecule') or 'genome').strip().lower()
    gene = (body.get('gene') or '').strip()
    if not name:
        abort(400, '参数不完整（name 必填）')
    if tree_tool not in ('fasttree', 'iqtree', 'nj'):
        abort(400, 'tree_tool 仅支持 fasttree / iqtree / nj')
    if molecule not in ('genome', 'cds', 'pep'):
        abort(400, 'molecule 仅支持 genome / cds / pep')
    if molecule in ('cds', 'pep') and not gene:
        abort(400, '基因级建树需提供基因/产物关键词（如 coat protein）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gb_collection import build_collection_phylo
        res = build_collection_phylo(name, tree_tool=tree_tool,
                                     molecule=molecule, gene=gene or None,
                                     logger=logger, prog=prog, cancel=cancel)
        logger.close()
        return res

    tag = name if molecule == 'genome' else f'{name}·{gene}({molecule})'
    tid = tm.start(cfg.tr(f'集合建树 {tag}', f'Collection phylo {tag}'), job)
    return jsonify({'task': tid})


@bp.route('/api/gb/download', methods=['POST'])
def api_gb_download():
    """body: {name, term?/accessions?, max} → 后台任务下载 GenBank 集合。"""
    body = request.get_json(force=True)
    name = (body.get('name') or '').strip()
    term = (body.get('term') or '').strip()
    accs = (body.get('accessions') or '').strip()
    if not name or (not term and not accs):
        abort(400, '参数不完整（name 与 term/accessions 必填）')
    if term and accs:
        abort(400, '检索式与 accession 列表二选一')
    max_records = max(1, min(int(body.get('max', 50) or 50), 10000))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gb_collection import (download_gb_collection,
                                      collection_records)
        res = download_gb_collection(name, term=term or None,
                                     accessions=accs or None,
                                     max_records=max_records,
                                     logger=logger, prog=prog, cancel=cancel)
        # 同步 FASTA 参考集（样品流程「NCBI 参考集合」可直接填集合名）
        try:
            from vp.ncbi_download import collection_dir as fasta_coll_dir
            from vp.utils import write_fasta_record
            fdir = fasta_coll_dir(name)
            os.makedirs(fdir, exist_ok=True)
            fa = os.path.join(fdir, 'refs.fa')
            n_fa = 0
            with safe_open(fa, 'wt') as f:
                for acc, _org, seq in collection_records(name):
                    if seq:
                        write_fasta_record(f, acc, seq)
                        n_fa += 1
            logger.log(f'FASTA 参考集: {fa}（{n_fa} 条）')
            res['n_fasta'] = n_fa
        except Exception as e:
            logger.log(f'FASTA 参考集同步失败（不影响 GenBank 集合）: {e}',
                       'WARN')
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'GenBank下载 {name}', f'GenBank download {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


@bp.route('/api/gb/import', methods=['POST'])
def api_gb_import():
    """body: {name, files:[本机 .gb 路径]} → 后台任务导入。"""
    body = request.get_json(force=True)
    name = (body.get('name') or '').strip()
    files = body.get('files') or []
    if not name or not files:
        abort(400, '参数不完整（name 与 files 必填）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gb_collection import import_local_gb
        res = import_local_gb(name, files, logger=logger, prog=prog)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'GenBank导入 {name}', f'GenBank import {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


@bp.route('/api/compare/preview')
def api_compare_preview():
    """读取集合的 compare.html，提取 body 可视化内容（共线性 SVG + 表）供卡内内嵌。

    参数: name=集合名。返回 {ok, html}（html 为 compare.html 的 <body> 内片段）。
    """
    import re as _re
    name = (request.args.get('name') or '').strip()
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, 'invalid collection name')
    from vp.gb_collection import gb_collection_dir
    try:
        cdir = gb_collection_dir(name)
    except ValueError as e:
        abort(400, str(e))
    p = os.path.join(cdir, 'compare', 'compare.html')
    if not os.path.isfile(p):
        abort(404, '该集合还没有比较结果，先运行比较')
    with safe_open(p) as f:
        html = f.read()
    m = _re.search(r'<body[^>]*>(.*)</body>', html, _re.S | _re.I)
    body = m.group(1) if m else html
    # 去掉可能存在的 <script>（compare.html 通常无，防御）
    body = _re.sub(r'<script[^>]*>.*?</script>', '', body, flags=_re.S | _re.I)
    return jsonify({'ok': True, 'html': body})


@bp.route('/api/compare/run', methods=['POST'])
def api_compare_run():
    """body: {collection?, files?, min_ident, min_cov} → 后台同属比较任务。"""
    body = request.get_json(force=True)
    collection = (body.get('collection') or '').strip()
    files = body.get('files') or []
    if not collection and not files:
        abort(400, '参数不完整（collection 或 files 至少一项）')
    try:
        min_ident = max(0.05, min(float(body.get('min_ident', 0.30)), 0.95))
        min_cov = max(0.05, min(float(body.get('min_cov', 0.50)), 1.0))
    except (TypeError, ValueError):
        abort(400, '阈值参数不合法')
    style = body.get('style') or 'lovis'
    if style not in ('lovis', 'category'):
        abort(400, 'style 仅支持 lovis / category')
    # LoVis4u 为共线性主输出（不可用时由内置引擎回退，见 synteny 日志）
    lovis4u = bool(body.get('lovis4u', True))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.synteny import run_comparison
        res = run_comparison(name=collection or None, files=files,
                             min_ident=min_ident, min_cov=min_cov,
                             style=style,
                             lovis4u_pdf=lovis4u,
                             logger=logger, prog=prog, cancel=cancel)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'同属比较 {collection or "自备文件"}',
                         f'Synteny {collection or "adhoc files"}'), job)
    return jsonify({'task': tid})


@bp.route('/compare/<name>/')
def page_compare(name):
    """集合的共线性比较结果页（compare.html，含相对引用的 plotly.min.js）。"""
    from vp.gb_collection import gb_collection_dir
    try:
        cdir = gb_collection_dir(name)
    except ValueError as e:
        abort(400, str(e))
    p = os.path.join(cdir, 'compare', 'compare.html')
    if not os.path.isfile(p):
        abort(404, '该集合还没有比较结果，先运行比较')
    return send_file(check_path(p, must_exist=True, in_platform=True))


@bp.route('/compare/<name>/<path:filename>')
def page_compare_file(name, filename):
    from vp.gb_collection import gb_collection_dir
    if '..' in filename:
        abort(400, '非法路径')
    try:
        cdir = gb_collection_dir(name)
    except ValueError as e:
        abort(400, str(e))
    p = check_path(os.path.join(cdir, 'compare', filename),
                   must_exist=True, in_platform=True)
    return send_file(p)


def _msa_group_dir(sample, group):
    """sample+group → 组目录（含 gb: 集合伪样品解析）。"""
    if sample.startswith('gb:'):
        from vp.gb_collection import gb_collection_dir
        try:
            return check_path(os.path.join(gb_collection_dir(sample[3:]),
                                           'phylo'),
                              must_exist=True, in_platform=True)
        except ValueError as e:
            raise ValueError(f'集合名不合法: {e}')
    from vp.pipeline import _safe_sample_name
    from vp.phylo import _safe_name
    sdir = check_path(os.path.join(DIRS['results'],
                                   _safe_sample_name(sample)),
                      must_exist=True, in_platform=True)
    return check_path(os.path.join(sdir, '05_phylo', _safe_name(group)),
                      must_exist=True, in_platform=True)


def _list_phylo_groups(pdir):
    """目录 → [{group, aln, trees}]（含比对文件的分组，优先 aln.trim）。"""
    groups = []
    for g in sorted(os.listdir(pdir)):
        gdir = os.path.join(pdir, g)
        if not os.path.isdir(gdir):
            continue
        aln = os.path.join(gdir, 'aln.trim.fasta')
        kind = 'trim'
        if not os.path.isfile(aln):
            aln = os.path.join(gdir, 'aln.fasta')
            kind = 'full'
        if os.path.isfile(aln):
            groups.append({'group': g, 'aln': kind, 'trees': _tree_files(gdir)})
    return groups


@bp.route('/api/msa/samples')
def api_msa_samples():
    """列出含比对结果的样品与集合建树产物（伪样品 gb:<集合名>）。"""
    out = []
    res = check_path(DIRS['results'], must_exist=False, in_platform=True)
    if os.path.isdir(res):
        for name in sorted(os.listdir(res)):
            if name.startswith('_'):      # 归档/内部样品不进下拉
                continue
            pdir = os.path.join(res, name, '05_phylo')
            if not os.path.isdir(pdir):
                continue
            groups = _list_phylo_groups(pdir)
            if groups:
                out.append({'sample': name, 'groups': groups})
    # GenBank 集合建树产物（比较基因组：下载/导入集合 → 直接建树）
    from vp.gb_collection import list_gb_collections
    try:
        cols = list_gb_collections()
    except OSError:
        cols = []
    for c in cols:
        pdir = os.path.join(c['dir'], 'phylo')
        if not os.path.isdir(pdir):
            continue
        aln = os.path.join(pdir, 'aln.trim.fasta')
        if not os.path.isfile(aln):
            aln = os.path.join(pdir, 'aln.fasta')
        if os.path.isfile(aln):
            out.append({'sample': f'gb:{c["name"]}',
                        'label': f'🧬 {c["name"]}（集合）',
                        'groups': [{'group': c['name'], 'aln': 'full',
                                    'trees': _tree_files(pdir)}]})
    return jsonify(out)


_TREE_FILES = (('tree.nwk', 'FastTree'),
               ('nj.nwk', 'NJ 快速树（identity 距离）'),
               ('iqtree.treefile', 'IQ-TREE ML 树'),
               ('iqtree.contree', 'IQ-TREE 一致树'))


def _tree_files(gdir):
    """组目录下实际存在的树文件列表（供树查看器下拉）。"""
    return [{'file': fn, 'label': label}
            for fn, label in _TREE_FILES
            if os.path.isfile(os.path.join(gdir, fn))]


@bp.route('/api/tree/data')
def api_tree_data():
    """sample+group(+file) → 进化树 Newick 与建树摘要（前端 Archaeopteryx.js 渲染）。"""
    sample = request.args.get('sample') or ''
    group = request.args.get('group') or ''
    fname = request.args.get('file') or ''
    try:
        gdir = _msa_group_dir(sample, group)
    except ValueError as e:
        abort(400, f'参数不合法: {e}')
    choices = ([(fname, '')] if fname else list(_TREE_FILES))
    tree_file = None
    for fn, _label in choices:
        p = os.path.join(gdir, fn)
        if fn and os.path.isfile(p):
            tree_file = p
            break
    if not tree_file:
        abort(404, '该组没有树文件（tree.nwk / nj.nwk / iqtree.treefile / iqtree.contree）')
    with safe_open(tree_file) as f:
        nwk = f.read().strip()
    if not nwk:
        abort(404, '树文件为空')
    # 从 summary.json 取建树工具/模型等摘要（无则忽略）
    meta = {}
    summary = os.path.join(gdir, 'summary.json')
    if os.path.isfile(summary):
        try:
            with safe_open(summary) as f:
                sj = json.load(f)
            g = (next((x for x in (sj.get('groups') or [])
                       if x.get('group') == group), None)
                 if isinstance(sj.get('groups'), list) else None) or {}
            base = os.path.basename(tree_file)
            if base == 'nj.nwk':
                tool = 'NJ'
            elif (base.endswith('.treefile') or base.endswith('.contree')
                    or g.get('tree_tool') == 'iqtree'
                    or sj.get('tree_tool') == 'iqtree'):
                tool = 'IQ-TREE'
            else:
                tool = 'FastTree'
            meta = {'tool': tool,
                    'model': g.get('tree_model') or sj.get('model'),
                    'logl': g.get('tree_logl') or sj.get('logl')}
        except (OSError, ValueError):
            pass
    return jsonify({'newick': nwk, 'file': os.path.basename(tree_file), **meta})


@bp.route('/api/msa/data')
def api_msa_data():
    """sample+group → SNP-only 展示数据（变异列矩阵/共识/多样性）。"""
    sample = request.args.get('sample') or ''
    group = request.args.get('group') or ''
    try:
        gdir = _msa_group_dir(sample, group)
    except ValueError as e:
        abort(400, f'参数不合法: {e}')
    aln = os.path.join(gdir, 'aln.trim.fasta')
    if not os.path.isfile(aln):
        aln = os.path.join(gdir, 'aln.fasta')
    if not os.path.isfile(aln):
        abort(404, '该组无比对文件')
    from vp.msa_view import snp_view_data
    try:
        return jsonify(snp_view_data(aln))
    except ValueError as e:
        abort(400, str(e))


@bp.route('/api/sdt/data')
def api_sdt_data():
    """sample+group → SDT 口径全长成对 identity 矩阵（结果中心在线热图）。

    优先读样品流程产物 sdt_matrix.csv；无则从比对 fasta 现算
    （与 vp/phylo.pairwise_identity_matrix 同口径）。
    """
    sample = request.args.get('sample') or ''
    group = request.args.get('group') or ''
    try:
        gdir = _msa_group_dir(sample, group)
    except ValueError as e:
        abort(400, f'参数不合法: {e}')
    csv = os.path.join(gdir, 'sdt_matrix.csv')
    if os.path.isfile(csv):
        names, mat = [], []
        with safe_open(csv) as f:
            rows = [line.rstrip('\n').split(',') for line in f if line.strip()]
        if len(rows) >= 2:
            names = rows[0][1:]
            for r in rows[1:]:
                mat.append([float(x) if x else 0.0 for x in r[1:]])
            return jsonify({'names': names, 'matrix': mat,
                            'source': 'sdt_matrix.csv'})
    aln = os.path.join(gdir, 'aln.trim.fasta')
    if not os.path.isfile(aln):
        aln = os.path.join(gdir, 'aln.fasta')
    if not os.path.isfile(aln):
        abort(404, '该组没有 SDT 矩阵，也没有比对文件')
    from vp.phylo import pairwise_identity_matrix
    names, mat = pairwise_identity_matrix(aln)
    return jsonify({'names': names, 'matrix': mat, 'source': 'aln'})

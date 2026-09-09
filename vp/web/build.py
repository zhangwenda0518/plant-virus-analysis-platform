# -*- coding: utf-8 -*-
"""数据库构建与库状态（自 app.py 拆出）。

Taxonomy / 宿主库 / 病毒库构建、Kraken2 转换、kv 索引、库状态、
文件浏览对话框。"""
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
from vp.web.common import _safe_sample
from vp.web.state import cfg, tool_runs_root as _tool_runs_root
from vp.web.tasks import tm

bp = Blueprint('build', __name__)


@bp.route('/api/tools')
def api_tools():
    return jsonify({'threads': cfg.threads, 'tools': cfg.tool_status()})


@bp.route('/api/dbs')
def api_dbs():
    from vp.kunpeng import db_ready
    out = {}
    for key in ('host', 'virus'):
        d = cfg.databases[key]
        out[key] = {'path': d,
                    'ready': db_ready(d) if os.path.isdir(d) else False}
    # 通用参考库（refvirus / rvdb / 其它 kunpeng 库目录自动探测）
    for key, mapping in (('refvirus', ('virus', 'ref')),
                         ('rvdb', ('virus', 'rvdb')),
                         ('k2viral', None)):
        if mapping:
            d = db_path(*mapping)
        else:
            d = os.path.join(DIRS['databases'], 'k2viral_db')
        out[key] = {'path': d,
                    'ready': db_ready(d) if os.path.isdir(d) else False}
    from vp.taxonomy import taxonomy_ready
    out['taxonomy'] = {'ready': taxonomy_ready(), 'path': DIRS['taxonomy']}
    # 病毒分类库统一为「自备预构建」：扫描 databases/ 下的 kunpeng 库目录
    # （宿主库除外）。使用者把库目录放进来即自动识别，工具②/④的病毒库
    # 下拉与构建页列表都以此为准。准入只看 db_ready（结构完整可用），
    # 不看 .building 标记——平台意外退出会留下孤儿标记，不应隐藏可用库。
    libs = []
    base = DIRS['databases']
    for dp, dn, fn in os.walk(base):
        if not any(f.startswith('hash_') and f.endswith('.k2d') for f in fn):
            continue
        rel = os.path.relpath(dp, base).replace('\\', '/')
        if rel in ('.', ''):
            continue
        # 宿主分类库（host/classify）不列入病毒库下拉
        if rel.startswith('host/'):
            continue
        try:
            if not db_ready(dp):
                continue
        except (ValueError, OSError):
            continue
        libs.append({'name': rel, 'path': f'databases/{rel}'})
    # 外部登记的自备库（platform.json extra_virus_libs，可为平台外任意位置）
    seen = {os.path.normpath(os.path.join(DIRS['databases'], l['name']))
            for l in libs}
    for p in cfg.extra_virus_libs:
        n = os.path.normpath(os.path.abspath(p))
        if n in seen or not _db_dir_ready(n):
            continue
        libs.append({'name': os.path.basename(n), 'path': p,
                     'external': True})
    out['virus_libs'] = libs
    return jsonify(out)


def _db_dir_ready(d):
    """kunpeng 库结构完整性（不限定平台目录，用于加载/登记自备库）。"""
    if not os.path.isdir(d):
        return False
    if not all(os.path.isfile(os.path.join(d, f))
               for f in ('opts.k2d', 'taxo.k2d', 'hash_config.k2d')):
        return False
    import glob as _glob
    return bool(_glob.glob(os.path.join(d, 'hash_*.k2d')))


def _abs_arg(raw):
    d = raw if os.path.isabs(raw) else os.path.join(PLATFORM_ROOT, raw)
    return os.path.normpath(os.path.abspath(d))


def _within_allowed_roots(d):
    from vp.config import write_roots
    d = os.path.normpath(os.path.abspath(d))
    return any(d == os.path.normpath(r)
               or d.startswith(os.path.normpath(r) + os.sep)
               for r in write_roots())


@bp.route('/api/load_taxonomy', methods=['POST'])
def api_load_taxonomy():
    """加载已有 NCBI Taxonomy（含 nodes.dmp/names.dmp 的目录）→
    复制 .dmp/.pkl 进平台标准位置 databases/tax/core。"""
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少 taxonomy 目录')
    d = _abs_arg(raw)
    if not os.path.isdir(d):
        abort(400, f'目录不存在: {raw}')
    missing = [f for f in ('nodes.dmp', 'names.dmp')
               if not os.path.isfile(os.path.join(d, f))]
    if missing:
        abort(400, f'缺少 {"、".join(missing)}（不是有效的 taxonomy 目录）')
    os.makedirs(DIRS['taxonomy'], exist_ok=True)
    import shutil as _shutil
    copied = []
    for fn in sorted(os.listdir(d)):
        src = os.path.join(d, fn)
        dst = os.path.join(DIRS['taxonomy'], fn)
        if not os.path.isfile(src) or not fn.lower().endswith(('.dmp', '.pkl')):
            continue
        if os.path.normpath(os.path.abspath(src)) == os.path.normpath(dst):
            copied.append(fn)                     # 已在标准位置，跳过
            continue
        _shutil.copyfile(src, dst)
        copied.append(fn)
    from vp.taxonomy import taxonomy_ready
    return jsonify({'ok': taxonomy_ready(), 'copied': copied,
                    'dir': DIRS['taxonomy']})


@bp.route('/api/load_host_db', methods=['POST'])
def api_load_host_db():
    """加载已构建宿主库（kunpeng）登记为平台宿主库（持久化 platform.json）。

    宿主库参与建库写流程，目录须位于平台目录内（或设置页配置的数据库根）。"""
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少宿主库目录')
    d = _abs_arg(raw)
    if not _db_dir_ready(d):
        abort(400, '不是有效的 kunpeng 库目录（缺 hash_*.k2d / opts.k2d / '
                   'taxo.k2d / hash_config.k2d）: ' + raw)
    if not _within_allowed_roots(d):
        abort(400, '宿主库目录须位于平台目录内（或先在「设置 → 数据库目录」'
                   '配置自定义数据库根）')
    cfg.databases['host'] = d
    cfg.save()
    return jsonify({'ok': True, 'path': d})


@bp.route('/api/register_virus_lib', methods=['POST'])
def api_register_virus_lib():
    """登记 / 移除外部已构建病毒库目录（持久化；分类只读，允许平台外任意位置）。"""
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少病毒库目录')
    d = _abs_arg(raw)
    if body.get('remove'):
        cfg.extra_virus_libs = [p for p in cfg.extra_virus_libs
                                if os.path.normpath(os.path.abspath(p)) != d]
        cfg.save()
        return jsonify({'ok': True, 'libs': cfg.extra_virus_libs})
    if not _db_dir_ready(d):
        abort(400, '不是有效的 kunpeng 库目录（缺 hash_*.k2d / opts.k2d / '
                   'taxo.k2d / hash_config.k2d）: ' + raw)
    if d not in cfg.extra_virus_libs:
        cfg.extra_virus_libs.append(d)
        cfg.save()
    return jsonify({'ok': True, 'libs': cfg.extra_virus_libs})


@bp.route('/api/browse')
def api_browse():
    """目录浏览（只读）。平台内目录可用相对路径；也支持浏览整台电脑
    （绝对路径，含盘符列表）。选择的数据输入文件可为任意位置；
    分析输出仍严格限制在平台目录内。"""
    raw = (request.args.get('path', '') or '.').strip()
    # all=1（选目录模式）：文件不限扩展名，便于查看目录里有什么
    show_all = request.args.get('all') == '1'
    if raw in ('.', '', '/', '~', '此电脑'):
        drives = [f'{c}:\\' for c in 'CDEFGHIJKLMNOPQRSTUVWXYZ'
                  if os.path.exists(f'{c}:\\')]
        return jsonify({'cwd': '此电脑', 'parent': '', 'dirs': drives,
                        'files': []})
    p = os.path.normpath(os.path.abspath(raw))
    if not os.path.isdir(p):
        abort(400, '目录不存在')
    dirs, files = [], []
    try:
        for name in sorted(os.listdir(p)):
            child = os.path.join(p, name)
            try:
                if os.path.isdir(child):
                    dirs.append(name)
                elif os.path.isfile(child) and (
                        show_all or name.lower().endswith(
                        ('.fastq.gz', '.fq.gz', '.fastq', '.fq',
                         '.fasta', '.fa', '.fna', '.fas',
                         '.fa.gz', '.fasta.gz', '.tsv', '.sra',
                         '.dmp', '.pkl', '.map',
                         '.gb', '.gbk', '.gbff', '.genbank'))):
                    files.append({'name': name,
                                  'size': fmt_size(os.path.getsize(child))})
            except OSError:
                continue                                  # 跳过无权限项
    except PermissionError:
        abort(400, '无权限访问该目录')
    except OSError:
        abort(400, '无法访问该路径')
    parent = os.path.dirname(p)
    return jsonify({'cwd': p,
                    'parent': '' if parent == p else parent,
                    'dirs': dirs, 'files': files})


def _job_build_taxonomy(log, prog, cancel):
    from vp.taxonomy import prepare_taxonomy
    logger = TaskLogger(callback=log)
    prog('taxonomy', 0.1, '准备 NCBI taxonomy')
    prepare_taxonomy(logger=logger)
    prog('taxonomy', 1.0, 'taxonomy 就绪')
    logger.close()
    return 'ok'


def _job_build_host_db(body):
    def job(log, prog, cancel):
        from vp.kunpeng import build_host_db
        logger = TaskLogger(callback=log)
        # 建库内存峰值与线程数成正比（每线程约 1GB 缓冲），默认限 8
        threads = int(body.get('threads') or min(8, cfg.threads))
        prog('build_host', 0.05, f"注入 taxid={body['taxid']} 并建库 (线程 {threads})")
        build_host_db(body['genome'], int(body['taxid']),
                      hash_capacity=body.get('hash_capacity', '256M'),
                      threads=threads,
                      logger=logger, rebuild=bool(body.get('rebuild')),
                      clean_mid=bool(body.get('clean_mid')))
        prog('build_host', 1.0, '宿主库就绪')
        logger.close()
        return 'ok'
    return job


def _job_build_virus_db(body):
    def job(log, prog, cancel):
        from vp.kunpeng import build_virus_db
        logger = TaskLogger(callback=log)
        prog('build_virus', 0.05, '解析 info 表并建库')
        build_virus_db(body['fasta'], body['info'],
                       hash_capacity=body.get('hash_capacity', '64M'),
                       threads=int(body.get('threads', cfg.threads)),
                       logger=logger, rebuild=bool(body.get('rebuild')),
                       clean_mid=bool(body.get('clean_mid')))
        prog('build_virus', 1.0, '病毒库就绪')
        logger.close()
        return 'ok'
    return job


@bp.route('/api/build_taxonomy', methods=['POST'])
def api_build_taxonomy():
    tid = tm.start(cfg.tr('下载/准备 Taxonomy', 'Download/prepare Taxonomy'),
                   _job_build_taxonomy, weight='light')
    return jsonify({'task': tid})


@bp.route('/api/build_host_db', methods=['POST'])
def api_build_host_db():
    body = request.get_json(force=True)
    for k in ('genome', 'taxid'):
        if not body.get(k):
            abort(400, f'缺少参数 {k}')
    check_path(body['genome'], must_exist=True)
    tid = tm.start(cfg.tr(f"宿主库构建 (taxid={body['taxid']})",
                         f"Host DB build (taxid={body['taxid']})"),
                   _job_build_host_db(body))
    return jsonify({'task': tid})


def _universal_db_job(source):
    """通用病毒库建库 job 构造器。"""
    def job(log, prog, cancel):
        from vp.universal_ref import build_universal_db
        logger = TaskLogger(callback=log)
        prog('prep', 0.05, '解析 accession→taxid 映射')
        prog('build', 0.3, 'kunpeng add-library + build-db（见日志）')
        res = build_universal_db(
            source, hash_capacity='2G', logger=logger, rebuild=True)
        prog('build', 0.9, '建库完成，校验库文件')
        prog('done', 1.0, '完成')
        logger.close()
        return res
    return job


@bp.route('/api/convert_kraken2', methods=['POST'])
def api_convert_kraken2():
    """Kraken2 库包/目录 → kunpeng 分片库（kunpeng hashshard，方式 C）。

    body: {tar: "databases/k2_viral_20260626.tar.gz"（平台内路径）,
           name: "k2viral"（目标库目录名 databases/<name>）,
           hash_capacity: "1G"}
    """
    body = request.get_json(force=True) or {}
    tar = (body.get('tar') or '').strip()
    name = (body.get('name') or 'k2viral').strip()
    if not tar:
        abort(400, '缺少 Kraken2 库包路径')
    check_path(tar, must_exist=True)
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '库名仅限字母数字-_')
    db_dir = os.path.join(DIRS['databases'], name)

    def job(log, prog, cancel):
        from vp.kunpeng import convert_kraken2
        logger = TaskLogger(callback=log)
        prog('convert', 0.2, '解包 + hashshard 转换（见日志）')
        res = convert_kraken2(tar, db_dir,
                              hash_capacity=body.get('hash_capacity') or '1G',
                              logger=logger)
        prog('done', 1.0, '完成')
        logger.close()
        return {'db_dir': res}

    tid = tm.start(cfg.tr(f'Kraken2 库转换 {name}', f'Kraken2 convert {name}'),
                   job)
    return jsonify({'task': tid})


@bp.route('/api/build_kv_index', methods=['POST'])
def api_build_kv_index():
    """构建「病毒鉴定库」：参考 FASTA → minibwa / salmon 比对索引。

    body: {fasta: 参考 FASTA（平台内外均可，只读）,
           engines: ['minibwa','salmon']（默认两者）,
           name: 库名（默认 kv_index，限字母数字-_）,
           threads: 线程数,
           ref_info: 可选，参考注释 TSV，与 FASTA 同目录归档}

    产物：virus-db/<name>/ 下 minibwa.{mbw,l2b}（或 <name>/minibwa/）与
    salmon_k31/，并写一份 manifest.json 记录来源与时间，便于 t-consensus /
    kvsuite 索引复用。与「病毒分类库」（kunpeng）无关，两者互不影响。
    """
    body = request.get_json(force=True) or {}
    fasta = (body.get('fasta') or '').strip()
    if not fasta:
        abort(400, '缺少参数 fasta')
    fasta = check_path(fasta if os.path.isabs(fasta)
                       else os.path.join(PLATFORM_ROOT, fasta),
                       must_exist=True)
    if not os.path.isfile(fasta):
        abort(400, f'参考 FASTA 不是文件：{fasta}')

    engines = body.get('engines') or ['minibwa', 'salmon']
    if isinstance(engines, str):
        engines = [e.strip() for e in engines.split(',') if e.strip()]
    engines = [e for e in engines if e in ('minibwa', 'salmon')]
    if not engines:
        abort(400, 'engines 需包含 minibwa 或 salmon')

    name = (body.get('name') or 'kv_index').strip()
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '库名仅限字母数字-_')
    out_dir = check_path(os.path.join(DIRS['virus_src'], name),
                         must_exist=False, in_platform=True)

    ref_info = (body.get('ref_info') or '').strip()
    if ref_info:
        ref_info = check_path(ref_info if os.path.isabs(ref_info)
                              else os.path.join(PLATFORM_ROOT, ref_info),
                              must_exist=True)

    try:
        threads = max(1, min(int(body.get('threads') or cfg.threads),
                             (os.cpu_count() or 4) * 4))
    except (TypeError, ValueError):
        threads = cfg.threads

    def job(log, prog, cancel):
        from known_virus_suite.kv_common import ToolRegistry, setup_logger
        from known_virus_suite.kv_engines import make_engine
        os.makedirs(out_dir, exist_ok=True)
        # kv_engines 的 logger 需是标准 logging.Logger（.info/.warning/.error）。
        # tag 带时间戳避免多次建库复用同一 logger 时 handler 累积。
        logger, _logs = setup_logger(
            out_dir, tag='kv_index_build_' + time.strftime('%H%M%S'))

        class _PanelHandler(logging.Handler):
            """把 kv_engines 的日志推到任务面板"""
            def emit(self, rec):
                try:
                    log(f'{time.strftime("%H:%M:%S", time.localtime(rec.created))} '
                        f'{rec.levelname} - {rec.getMessage()}')
                except Exception:                # noqa: BLE001
                    pass

        logger.addHandler(_PanelHandler())
        reg = ToolRegistry(logger)
        reg.probe()
        logger.info(f"参考 FASTA: {fasta}")
        logger.info(f"引擎: {', '.join(engines)}  线程: {threads}")
        logger.info(f"产物目录: {out_dir}")

        made, failed = [], []
        for i, eng_name in enumerate(engines):
            if cancel is not None and cancel.is_set():
                raise RuntimeError('用户取消')
            prog(eng_name, 0.05 + 0.9 * i / len(engines), f'构建 {eng_name} 索引')
            if not reg.has(eng_name):
                failed.append(f'{eng_name}: 未找到可执行文件')
                logger.error(f'{eng_name} 不可用，跳过')
                continue
            try:
                engine = make_engine(eng_name, reg, threads=threads,
                                     logger=logger)
                # 各引擎的 build_index 落在 out_dir/<引擎子目录>
                idx = engine.build_index(fasta, out_dir, threads)
                made.append({'engine': eng_name, 'index': str(idx)})
                logger.info(f'{eng_name} 索引就绪 -> {idx}')
            except Exception as exc:                # noqa: BLE001
                failed.append(f'{eng_name}: {exc}')
                logger.error(f'{eng_name} 建索引失败：{exc}')

        if not made:
            raise RuntimeError('所有引擎建索引均失败：' + '; '.join(failed))

        manifest = {
            'name': name,
            'reference': fasta,
            'ref_info': ref_info or '',
            'engines': made,
            'failed': failed,
            'threads': threads,
            'built_at': time.strftime('%Y-%m-%d %H:%M:%S'),
        }
        try:
            with open(os.path.join(out_dir, 'manifest.json'), 'w',
                      encoding='utf-8') as fh:
                json.dump(manifest, fh, ensure_ascii=False, indent=2)
        except OSError as exc:
            logger.warning(f'manifest 写入失败：{exc}')
        # 关闭 FileHandler，否则日志文件句柄会泄漏并锁住产物目录
        for h in list(logger.handlers):
            try:
                h.close()
            except Exception:                    # noqa: BLE001
                pass
            logger.removeHandler(h)
        prog('done', 1.0, '建库完成')
        return {'out': out_dir, 'engines': [m['engine'] for m in made],
                'failed': failed}

    tid = tm.start(cfg.tr(f'病毒鉴定库构建 {name}', f'Virus index build {name}'), job)
    return jsonify({'task': tid, 'out': out_dir})


def _kv_index_probe(d, name=None):
    """探测单个目录是否是一个可用的病毒鉴定库（minibwa / salmon 索引）。

    返回 lib 字典（含 minibwa/salmon 布尔与各自路径），不可用时返回 None。
    仅探测结构，不判断目录所在位置，平台内与平台外目录共用。
    """
    lib = {'name': name or os.path.basename(d), 'path': d, 'minibwa': False,
           'salmon': False, 'manifest': None}
    for cand in (os.path.join(d, 'minibwa.mbw'),
                 os.path.join(d, 'minibwa', 'minibwa.mbw'),
                 os.path.join(d, 'harmonized_minibwa.mbw')):
        if os.path.isfile(cand):
            lib['minibwa'] = True
            lib['minibwa_path'] = cand
            break
    for cand in (os.path.join(d, 'salmon_k31', 'info.json'),
                 os.path.join(d, 'salmon', 'salmon_k31', 'info.json')):
        if os.path.isfile(cand):
            lib['salmon'] = True
            lib['salmon_path'] = os.path.dirname(cand)
            break
    mf = os.path.join(d, 'manifest.json')
    if os.path.isfile(mf):
        try:
            with open(mf, encoding='utf-8') as fh:
                lib['manifest'] = json.load(fh)
        except (OSError, ValueError):
            pass
    return lib if (lib['minibwa'] or lib['salmon']) else None


@bp.route('/api/kv_index_list')
def api_kv_index_list():
    """列出可用的鉴定库（virus-db 下的 + platform.json 登记的外部目录）与引擎可用性。"""
    base = DIRS['virus_src']
    out = {'base': base, 'libs': [], 'engines': {}}
    try:
        from known_virus_suite.kv_common import ToolRegistry
        reg = ToolRegistry()
        reg.probe()
        for e in ('minibwa', 'salmon'):
            out['engines'][e] = {'available': reg.has(e),
                                 'path': str(reg.paths.get(e) or '')}
    except Exception as exc:                    # noqa: BLE001
        out['engines_error'] = str(exc)
    seen = set()
    if os.path.isdir(base):
        for nm in sorted(os.listdir(base)):
            d = os.path.join(base, nm)
            if not os.path.isdir(d):
                continue
            lib = _kv_index_probe(d, nm)
            if lib is None and nm == 'kv_index':
                # 默认库即使索引缺失也列出（便于排查），但不计入可用
                lib = {'name': nm, 'path': d, 'minibwa': False,
                       'salmon': False, 'manifest': None}
            if lib is None:
                continue
            lib['external'] = False
            out['libs'].append(lib)
            seen.add(os.path.normpath(os.path.abspath(d)))
    # 外部登记目录（platform.json extra_kv_indexes，可为平台外任意位置，只读使用）
    for p in getattr(cfg, 'extra_kv_indexes', []) or []:
        n = os.path.normpath(os.path.abspath(p))
        if n in seen or not os.path.isdir(n):
            continue
        lib = _kv_index_probe(n)
        if lib is None:
            continue
        lib['external'] = True
        out['libs'].append(lib)
        seen.add(n)
    out['extra'] = [p for p in (getattr(cfg, 'extra_kv_indexes', []) or [])]
    return jsonify(out)


@bp.route('/api/kv_index_register', methods=['POST'])
def api_kv_index_register():
    """登记 / 移除外部病毒鉴定库目录（持久化；只读使用，允许平台外任意位置）。

    与 /api/register_virus_lib 同语义：把平台外已建好的 minibwa/salmon
    索引目录挂进平台，之后在工具卡的「鉴定库」下拉里即可选择复用。
    """
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少鉴定库目录')
    d = _abs_arg(raw)
    if body.get('remove'):
        cfg.extra_kv_indexes = [
            p for p in (getattr(cfg, 'extra_kv_indexes', []) or [])
            if os.path.normpath(os.path.abspath(p)) != d]
        cfg.save()
        return jsonify({'ok': True, 'dirs': cfg.extra_kv_indexes})
    if not os.path.isdir(d):
        abort(400, '目录不存在: ' + raw)
    lib = _kv_index_probe(d)
    if lib is None:
        abort(400, '不是有效的鉴定库目录（未找到 minibwa.mbw 或 salmon_k31/info.json）: ' + raw)
    if d not in (getattr(cfg, 'extra_kv_indexes', []) or []):
        cfg.extra_kv_indexes = list(getattr(cfg, 'extra_kv_indexes', []) or []) + [d]
        cfg.save()
    return jsonify({'ok': True, 'dirs': cfg.extra_kv_indexes,
                    'lib': {'name': lib['name'], 'path': d,
                            'minibwa': lib['minibwa'], 'salmon': lib['salmon']}})


@bp.route('/api/build_universal_db', methods=['POST'])
def api_build_universal_db():
    """构建通用病毒参考库（refvirus=NCBI RefSeq Viral / rvdb=RVDB C-RVDB）。"""
    body = request.get_json(force=True) or {}
    source = body.get('source') or 'refvirus'
    if source not in ('refvirus', 'rvdb'):
        abort(400, 'source 需为 refvirus 或 rvdb')
    name = {'refvirus': '通用病毒库构建 (RefSeq Viral)',
            'rvdb': 'RVDB 库构建 (C-RVDB)'}[source]
    tid = tm.start(cfg.tr(name, name), _universal_db_job(source))
    return jsonify({'task': tid})


@bp.route('/api/build_virus_db', methods=['POST'])
def api_build_virus_db():
    body = request.get_json(force=True)
    for k in ('fasta', 'info'):
        if not body.get(k):
            abort(400, f'缺少参数 {k}')
    check_path(body['fasta'], must_exist=True)
    check_path(body['info'], must_exist=True)
    tid = tm.start(cfg.tr('病毒库构建', 'Virus DB build'), _job_build_virus_db(body))
    return jsonify({'task': tid})

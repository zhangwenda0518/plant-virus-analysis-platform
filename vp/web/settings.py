# -*- coding: utf-8 -*-
"""设置中心 + 存储/磁盘水位（自 app.py 拆出）。

路由：
    GET  /settings              设置页
    GET  /api/settings          读取（语言/线程/邮箱/默认参数/目录/库状态）
    POST /api/settings          保存（逐项校验后持久化到 platform.json）
    GET  /api/storage           磁盘容量 + 目录占用（占用走后台扫描 + 缓存）
    POST /api/storage/rescan    强制重新统计目录占用

设计要点：磁盘剩余用 shutil.disk_usage 毫秒级实时算；目录占用要遍历
6 万+ 文件（databases 12GB+），实测数十秒，因此走后台线程 + 10 分钟缓存，
绝不放在请求路径上同步扫描。
"""
import os
import threading
import time

from flask import Blueprint, abort, jsonify, render_template, request

from vp.config import DIRS, PLATFORM_ROOT
from vp.utils import fmt_size
from vp.web.state import cfg, tool_runs_root

bp = Blueprint('settings', __name__)

# ------------------------------------------------------------------
# 存储与磁盘水位
# ------------------------------------------------------------------
STORAGE_WARN_GB = 20          # 剩余低于此值触发水位横幅
_STORAGE_TTL = 600            # 目录占用缓存有效期（秒）
_storage = {'ts': 0, 'scanning': False, 'items': [], 'error': None,
            'scanned_at': None}
_storage_lock = threading.Lock()


def _dir_size(path):
    """累加目录内文件大小与文件数（os.walk，容忍无权限/已删除的文件）。"""
    total = 0
    count = 0
    for cur, _dirs, files in os.walk(path):
        for fn in files:
            try:
                total += os.path.getsize(os.path.join(cur, fn))
            except OSError:
                continue
            count += 1
    return total, count


def _scan_targets():
    """待统计的目录清单：(显示名, 路径, 类别)。"""
    items = [
        ('databases', DIRS['databases'], 'database'),
        ('host-db', DIRS['host_src'], 'database'),
        ('virus-db', DIRS['virus_src'], 'database'),
        ('results', DIRS['results'], 'output'),
        ('downloads', DIRS['downloads'], 'output'),
        ('tool_runs', tool_runs_root(), 'output'),
        ('meta_search', DIRS['meta_search'], 'output'),
        ('logan', DIRS['logan'], 'output'),
        ('logs', DIRS['logs'], 'output'),
        ('tasks', DIRS['tasks'], 'output'),
        ('fastq', DIRS.get('fastq', ''), 'input'),
        ('dist', os.path.join(PLATFORM_ROOT, 'dist'), 'build'),
    ]
    return [(n, p, k) for n, p, k in items if p and os.path.isdir(p)]


def _do_storage_scan():
    """后台扫描：平台各目录 + databases 下逐库细分。"""
    out = []
    try:
        for name, path, kind in _scan_targets():
            size, count = _dir_size(path)
            entry = {'name': name, 'path': path, 'kind': kind,
                     'size': size, 'size_h': fmt_size(size),
                     'files': count, 'children': []}
            # 数据库根太大，按子库拆开才有参考价值
            if name == 'databases':
                try:
                    subs = sorted(d for d in os.listdir(path)
                                  if os.path.isdir(os.path.join(path, d)))
                except OSError:
                    subs = []
                for sub in subs:
                    sp = os.path.join(path, sub)
                    ssz, scnt = _dir_size(sp)
                    entry['children'].append({
                        'name': sub, 'path': sp, 'kind': 'database',
                        'size': ssz, 'size_h': fmt_size(ssz), 'files': scnt})
                entry['children'].sort(key=lambda x: -x['size'])
            out.append(entry)
        out.sort(key=lambda x: -x['size'])
        with _storage_lock:
            _storage['items'] = out
            _storage['ts'] = time.time()
            _storage['scanned_at'] = time.strftime('%H:%M:%S')
            _storage['error'] = None
    except Exception as e:                      # 扫描失败不应影响平台
        with _storage_lock:
            _storage['error'] = str(e)
    finally:
        with _storage_lock:
            _storage['scanning'] = False


def _disk_snapshot():
    """盘容量快照（实时，含数据库根所在盘——可能与平台根不同盘）。"""
    import shutil
    from vp.config import get_database_root
    roots = [('平台', PLATFORM_ROOT)]
    db_root = get_database_root()
    if db_root:
        roots.append(('数据库', db_root))
    out = []
    for label, path in roots:
        try:
            u = shutil.disk_usage(path)
            out.append({'label': label, 'path': path,
                        'total': u.total, 'free': u.free, 'used': u.used,
                        'total_h': fmt_size(u.total),
                        'free_h': fmt_size(u.free),
                        'used_h': fmt_size(u.used),
                        'used_pct': round(u.used / u.total * 100, 1) if u.total else 0,
                        'warn': u.free < STORAGE_WARN_GB * 1e9})
        except OSError as e:
            out.append({'label': label, 'path': path, 'error': str(e),
                        'warn': False})
    return out


# ------------------------------------------------------------------
# 路由
# ------------------------------------------------------------------
@bp.route('/settings')
def page_settings():
    return render_template('settings.html')


@bp.route('/api/settings', methods=['GET'])
def api_settings_get():
    """设置中心读取：语言/线程/邮箱/默认参数/目录/数据库状态。"""
    from vp.config import Config
    fields = []
    for f in Config.DEFAULT_FIELDS:
        key, typ, dflt, lz, le = f[0], f[1], f[2], f[3], f[4]
        fields.append({'key': key, 'type': typ, 'def': dflt,
                       'label_zh': lz, 'label_en': le,
                       'opts': (f[5] if len(f) > 5 else [])})
    return jsonify({
        'language': cfg.language,
        'email': cfg.email,
        'threads': cfg.threads,
        'cpu_count': os.cpu_count() or 4,
        'defaults': cfg.defaults,
        'default_fields': fields,
        'databases': cfg.databases,
        'output_root': cfg.output_root,
        'input_root': cfg.input_root,
        'database_root': cfg.database_root,
        'dirs': {'results': DIRS['results'],
                 'tool_runs': tool_runs_root(),
                 'databases': DIRS['databases'],
                 'fastq': DIRS.get('fastq')
                          or os.path.join(PLATFORM_ROOT, 'fastq'),
                 'uploads': DIRS.get('uploads')
                            or os.path.join(PLATFORM_ROOT, 'uploads'),
                 'host_src': DIRS.get('host_src')
                             or os.path.join(PLATFORM_ROOT, 'host-db'),
                 'virus_src': DIRS.get('virus_src')
                              or os.path.join(PLATFORM_ROOT, 'virus-db'),
                 'logs': DIRS['logs']},
    })


@bp.route('/api/settings', methods=['POST'])
def api_settings_set():
    """设置中心保存：language/email/threads/defaults（逐项校验后持久化）。"""
    body = request.get_json(force=True) or {}
    from vp.config import Config
    known = {k: t for k, t, *_r in Config.DEFAULT_FIELDS}
    if 'language' in body:
        if body['language'] not in ('zh', 'en'):
            abort(400, 'language 仅支持 zh / en')
        cfg.language = body['language']
    if 'email' in body:
        cfg.email = str(body['email'] or '').strip()
    if 'threads' in body:
        try:
            cfg.threads = max(1, int(body['threads']))
        except (TypeError, ValueError):
            abort(400, 'threads 必须是整数')
    if 'defaults' in body and isinstance(body['defaults'], dict):
        for k, v in body['defaults'].items():
            if k not in known:
                continue
            kind = known[k]
            try:
                if kind == 'int':
                    cfg.defaults[k] = int(v)
                elif kind == 'float':
                    cfg.defaults[k] = float(v)
                elif kind == 'bool':
                    cfg.defaults[k] = bool(v)
                else:
                    cfg.defaults[k] = str(v)
            except (TypeError, ValueError):
                abort(400, f'参数 {k} 的值不合法: {v!r}')
    if 'output_root' in body:
        try:
            cfg.set_output_root(str(body.get('output_root') or ''))
        except ValueError as e:
            abort(400, f'输出目录不合法: {e}')
        except OSError as e:
            abort(400, f'输出目录无法创建: {e}')
    if 'input_root' in body:
        try:
            cfg.set_input_root(str(body.get('input_root') or ''))
        except ValueError as e:
            abort(400, f'输入目录不合法: {e}')
        except OSError as e:
            abort(400, f'输入目录无法创建: {e}')
    if 'database_root' in body:
        try:
            cfg.set_database_root(str(body.get('database_root') or ''))
        except ValueError as e:
            abort(400, f'数据库目录不合法: {e}')
        except OSError as e:
            abort(400, f'数据库目录无法创建: {e}')
    try:
        cfg.save()
    except OSError as e:
        abort(500, f'配置保存失败: {e}')
    # 返回数据库就绪状态（前端切换 database_root 后可即时提示缺库）
    from vp.kunpeng import db_ready
    from vp.taxonomy import taxonomy_ready
    return jsonify({'ok': True, 'language': cfg.language,
                    'defaults': cfg.defaults,
                    'dbs': {'host': db_ready(cfg.databases['host']),
                            'virus': db_ready(cfg.databases['virus']),
                            'taxonomy': taxonomy_ready(),
                            'root': getattr(cfg, 'database_root', '')}})


@bp.route('/api/storage')
def api_storage():
    """磁盘与存储占用。

    ?detail=1 时返回目录占用明细；无缓存或缓存过期则后台起一次扫描，
    未就绪时 ready=false（前端可稍后再问），请求本身不阻塞。
    """
    detail = request.args.get('detail') in ('1', 'true', 'yes')
    disks = _disk_snapshot()
    resp = {'disks': disks,
            'warn': any(d.get('warn') for d in disks),
            'threshold_gb': STORAGE_WARN_GB,
            'scan': {'ready': False, 'scanning': False,
                     'scanned_at': None, 'items': [], 'error': None}}
    if not detail:
        return jsonify(resp)
    with _storage_lock:
        fresh = (time.time() - _storage['ts']) < _STORAGE_TTL
        if not fresh and not _storage['scanning']:
            _storage['scanning'] = True
            threading.Thread(target=_do_storage_scan, daemon=True,
                             name='storage-scan').start()
        resp['scan'] = {'ready': bool(_storage['items']) and fresh,
                        'scanning': _storage['scanning'],
                        'scanned_at': _storage['scanned_at'],
                        'items': _storage['items'] if fresh else [],
                        'error': _storage['error']}
    return jsonify(resp)


@bp.route('/api/storage/rescan', methods=['POST'])
def api_storage_rescan():
    """强制重新统计目录占用（忽略缓存）。"""
    with _storage_lock:
        if _storage['scanning']:
            return jsonify({'ok': False, 'msg': '统计进行中'})
        _storage['ts'] = 0
        _storage['items'] = []
        _storage['scanning'] = True
    threading.Thread(target=_do_storage_scan, daemon=True,
                     name='storage-scan').start()
    return jsonify({'ok': True})

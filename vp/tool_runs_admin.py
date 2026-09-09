# -*- coding: utf-8 -*-
"""tool_runs 目录管理：规范化、归档、清理。

目录规范（见 README 与 tool_runs/*/README.txt）：

    tool_runs/
    ├─ _archive/   固定：历史运行归档（_archive/<tool>/<run>/）
    ├─ _tmp/       固定：冒烟测试与一次性实验（随时可清）
    ├─ _scripts/   固定：运维/补丁脚本
    ├─ _reports/   固定：跨运行汇总产物
    └─ <tool>_<YYYYMMDD_HHMMSS>/   动态活动区（扁平，30 天滚动归档）

CLI：
  python main.py tool-runs status                        # 现状统计
  python main.py tool-runs organize [--dry-run]          # 杂项归位
  python main.py tool-runs archive --before 20260901     # 归档旧运行
  python main.py tool-runs clean [--active-days 30] [--archive-days 90] [--dry-run]
"""
import os
import re
import time
import shutil

FIXED_DIRS = ('_archive', '_tmp', '_scripts', '_reports')
_RUN_RE = re.compile(r'^([a-zA-Z0-9]+)_(\d{8})_(\d{6})$')

# 杂项归位规则：模式 → 目标固定目录（目录名不含标准时间戳时按此分流）
_ORGANIZE_RULES = [
    (re.compile(r'^_[a-z0-9_]*$', re.IGNORECASE), '_tmp'),        # _xxx 测试目录
    (re.compile(r'.*\.py$', re.IGNORECASE), '_scripts'),          # 误落的脚本
    (re.compile(r'.*_out$'), '_tmp'),                             # 工具日志残留目录
    (re.compile(r'.*(test|smoke).*', re.IGNORECASE), '_tmp'),     # 测试命名
]


def _root():
    from .config import DIRS
    return DIRS['tool_runs']


def ensure_fixed_dirs(root=None):
    """确保四个固定子目录存在（服务启动或 CLI 均可调用）。"""
    root = root or _root()
    for d in FIXED_DIRS:
        os.makedirs(os.path.join(root, d), exist_ok=True)
    return root


def _fmt(n):
    for unit in ('B', 'KB', 'MB', 'GB'):
        if n < 1024 or unit == 'GB':
            return f'{n:.1f} {unit}' if unit != 'B' else f'{n:.0f} B'
        n /= 1024
    return f'{n:.1f} GB'


def _dir_stat(path):
    total, n = 0, 0
    for cur, _dirs, files in os.walk(path):
        for fn in files:
            try:
                total += os.path.getsize(os.path.join(cur, fn))
            except OSError:
                continue
            n += 1
    return n, total


def list_runs(root=None):
    """活动区运行列表：[{name, tool, date, mtime, path}]（只收标准命名）。"""
    root = root or _root()
    out = []
    for name in os.listdir(root) if os.path.isdir(root) else []:
        m = _RUN_RE.match(name)
        if not m:
            continue
        p = os.path.join(root, name)
        if not os.path.isdir(p):
            continue
        try:
            mt = os.path.getmtime(p)
        except OSError:
            mt = 0
        out.append({'name': name, 'tool': m.group(1), 'date': m.group(2),
                    'mtime': mt, 'path': p})
    out.sort(key=lambda x: -x['mtime'])
    return out


def status(root=None):
    """统计：活动区/固定区/归档区的条目数与体积。"""
    root = root or _root()
    runs = list_runs(root)
    info = {'root': root, 'active': [], 'fixed': {}, 'archive': {}}
    total_sz = 0
    for r in runs:
        n, sz = _dir_stat(r['path'])
        total_sz += sz
        info['active'].append({'name': r['name'], 'tool': r['tool'],
                               'date': r['date'], 'files': n,
                               'size': sz, 'size_h': _fmt(sz)})
    for d in FIXED_DIRS:
        p = os.path.join(root, d)
        if os.path.isdir(p):
            n, sz = _dir_stat(p)
            subs = len([x for x in os.listdir(p)
                        if os.path.isdir(os.path.join(p, x))])
            info['fixed'][d] = {'files': n, 'size': sz, 'size_h': _fmt(sz),
                                'dirs': subs}
    return info


def organize(dry_run=False, logger=None):
    """把活动区里非标准命名的条目按规则归位到固定目录。"""
    root = ensure_fixed_dirs()

    def _log(m):
        if logger:
            logger.log(m)
        else:
            print(m)

    moved, skipped = [], []
    for name in sorted(os.listdir(root)):
        p = os.path.join(root, name)
        if name in FIXED_DIRS or _RUN_RE.match(name):
            continue                                # 固定目录 / 标准运行
        for pat, target in _ORGANIZE_RULES:
            if pat.match(name):
                dst_dir = os.path.join(root, target)
                os.makedirs(dst_dir, exist_ok=True)
                dst = os.path.join(dst_dir, name)
                if os.path.exists(dst):
                    skipped.append((name, f'{target}/ 已存在同名'))
                    break
                moved.append((name, target))
                if not dry_run:
                    shutil.move(p, dst)
                break
        else:
            skipped.append((name, '无法识别类型，保留原位（请人工判断）'))
    for name, target in moved:
        _log(f'  {"将移动" if dry_run else "已移动"} {name} → {target}/')
    for name, why in skipped:
        _log(f'  跳过 {name}（{why}）')
    return {'moved': len(moved), 'skipped': len(skipped), 'dry_run': dry_run}


def archive(before=None, older_days=None, dry_run=False, logger=None):
    """归档旧运行到 _archive/<tool>/<run>/。

    before: 'YYYYMMDD'（归档该日期之前的运行）；older_days: 归档早于 N 天的。
    两者至少给一个。
    """
    root = ensure_fixed_dirs()

    def _log(m):
        if logger:
            logger.log(m)
        else:
            print(m)

    if before is None and older_days is None:
        raise ValueError('需指定 --before YYYYMMDD 或 --older-days N')
    cutoff_ts = None
    if older_days is not None:
        cutoff_ts = time.time() - older_days * 86400
    todo = []
    for r in list_runs(root):
        if before and r['date'] >= before:
            continue
        if cutoff_ts is not None and r['mtime'] >= cutoff_ts:
            continue
        todo.append(r)
    n = 0
    for r in todo:
        dst_dir = os.path.join(root, '_archive', r['tool'])
        dst = os.path.join(dst_dir, r['name'])
        if os.path.exists(dst):
            _log(f'  跳过 {r["name"]}（归档区已存在）')
            continue
        _log(f'  {"将归档" if dry_run else "归档"} {r["name"]} → _archive/{r["tool"]}/')
        if not dry_run:
            os.makedirs(dst_dir, exist_ok=True)
            shutil.move(r['path'], dst)
        n += 1
    return {'archived': n, 'dry_run': dry_run}


def clean(active_days=None, archive_days=None, dry_run=False, logger=None):
    """清理：删除活动区超期运行、归档区超期归档。

    删除是不可逆的（不进回收站）——dry_run 先跑一遍确认清单再正式执行。
    """
    root = ensure_fixed_dirs()

    def _log(m):
        if logger:
            logger.log(m)
        else:
            print(m)

    removed, freed = 0, 0
    if active_days is not None:
        cutoff = time.time() - active_days * 86400
        for r in list_runs(root):
            if r['mtime'] >= cutoff:
                continue
            n, sz = _dir_stat(r['path'])
            _log(f'  {"将删除" if dry_run else "删除"} {r["name"]} '
                 f'({_fmt(sz)}, {time.strftime("%Y-%m-%d", time.localtime(r["mtime"]))})')
            if not dry_run:
                shutil.rmtree(r['path'], ignore_errors=True)
            removed += 1
            freed += sz
    if archive_days is not None:
        cutoff = time.time() - archive_days * 86400
        ar = os.path.join(root, '_archive')
        for tool in sorted(os.listdir(ar)) if os.path.isdir(ar) else []:
            tdir = os.path.join(ar, tool)
            if not os.path.isdir(tdir):
                continue
            for run in sorted(os.listdir(tdir)):
                p = os.path.join(tdir, run)
                if not os.path.isdir(p):
                    continue
                try:
                    mt = os.path.getmtime(p)
                except OSError:
                    continue
                if mt >= cutoff:
                    continue
                n, sz = _dir_stat(p)
                _log(f'  {"将删除归档" if dry_run else "删除归档"} _archive/{tool}/{run} '
                     f'({_fmt(sz)})')
                if not dry_run:
                    shutil.rmtree(p, ignore_errors=True)
                removed += 1
                freed += sz
    return {'removed': removed, 'freed': freed, 'freed_h': _fmt(freed),
            'dry_run': dry_run}

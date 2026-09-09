# -*- coding: utf-8 -*-
"""数据库迁移工具：把 databases / host-db / virus-db 复制（或搬移）到
外置位置，并把 platform.json 的 database_root 指过去，实现「软件与数据
库分离部署」。

典型场景：
- 本机磁盘吃紧：把 47GB 数据库挪到大容量盘；
- 打包分发：软件包不含数据库（体积从 60+GB 降到 ~1GB），数据库包单独拷到
  目标机器后，在设置页「数据库目录」填路径（或跑本工具 --to）即完成对接。

流程（copy 模式，推荐）：
  1. 校验目标目录：绝对路径、不能在平台根内、磁盘剩余充足（> 源大小 + 10GB）
  2. robocopy /E /MT /Z（断点续传）复制三个目录 → 校验文件数与总字节一致
  3. 通过后才写 platform.json database_root（未通过则配置不变，源目录原样）
  4. move 模式（--mode move）在复制校验通过后删除源（同 move 语义）

用法：
  python main.py db-migrate --to D:\\plant_virus_databases            # 复制+切换
  python main.py db-migrate --to D:\\plant_virus_databases --mode move # 搬移（删源）
  python main.py db-migrate --check --to D:\\plant_virus_databases     # 只校验已有
"""
import os
import sys
import shutil
import subprocess

_DB_SUBDIRS = ('databases', 'host-db', 'virus-db')


def _fmt(n):
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if n < 1024 or unit == 'TB':
            return f'{n:.1f} {unit}' if unit != 'B' else f'{n:.0f} B'
        n /= 1024
    return f'{n:.1f} TB'


def _log(logger, msg):
    if logger:
        logger.log(msg)
    else:
        print(msg)


def _dir_stat(path):
    """返回 (文件数, 总字节)；目录不存在返回 (0, 0)。"""
    total, n = 0, 0
    if not os.path.isdir(path):
        return 0, 0
    for cur, _dirs, files in os.walk(path):
        for fn in files:
            try:
                total += os.path.getsize(os.path.join(cur, fn))
            except OSError:
                continue
            n += 1
    return n, total


def _robocopy(src, dst, logger=None):
    """robocopy 复制（多线程 + 断点续传）；返回 True/False。"""
    os.makedirs(dst, exist_ok=True)
    cmd = ['robocopy', src, dst, '/E', '/MT:16', '/Z', '/R:2', '/W:3',
           '/NFL', '/NDL', '/NP', '/NJH']
    _log(logger, f'    robocopy {src} → {dst}')
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           errors='replace')
        # robocopy 返回码：0-7 均表示成功（1 = 有文件被复制）
        if r.returncode >= 8:
            tail = (r.stdout or r.stderr or '')[-500:]
            _log(logger, f'    ✘ robocopy 失败(码 {r.returncode}): {tail}')
            return False
        return True
    except OSError as e:
        _log(logger, f'    ✘ robocopy 无法启动: {e}（请确认系统自带 robocopy）')
        return False


def _copy_fallback(src, dst, logger=None):
    """shutil 兜底复制（robocopy 不可用时；无断点续传）。"""
    _log(logger, f'    (shutil 兜底复制) {src} → {dst}')
    try:
        shutil.copytree(src, dst, dirs_exist_ok=True)
        return True
    except OSError as e:
        _log(logger, f'    ✘ 复制失败: {e}')
        return False


def migrate(to_dir, mode='copy', logger=None, dry_run=False):
    """迁移数据库根。返回 {'ok': bool, 'error': str, 'detail': str}。"""
    from .config import (PLATFORM_ROOT, DIRS, get_database_root)
    import platform as _platform

    to_dir = (to_dir or '').strip().strip('"')
    if not to_dir:
        return {'ok': False, 'error': '目标目录不能为空'}
    if not os.path.isabs(to_dir):
        return {'ok': False, 'error': f'目标目录需为绝对路径: {to_dir}'}
    to_dir = os.path.normpath(os.path.abspath(to_dir))
    if (to_dir == PLATFORM_ROOT
            or to_dir.startswith(PLATFORM_ROOT + os.sep)):
        return {'ok': False, 'error': '目标目录不能位于平台目录内'
                f'（{PLATFORM_ROOT}）——分离部署要放到平台目录之外'}
    if mode not in ('copy', 'move'):
        return {'ok': False, 'error': f'模式需为 copy 或 move: {mode!r}'}

    cur_root = get_database_root()
    if cur_root and os.path.normpath(cur_root) == to_dir:
        return {'ok': False, 'error': f'数据库根已是 {to_dir}，无需迁移'}

    # 源目录存在性：至少 databases 应存在且有内容
    src_dbs = os.path.join(PLATFORM_ROOT, 'databases')
    if not os.path.isdir(src_dbs):
        return {'ok': False, 'error': f'源目录不存在: {src_dbs}'}

    # 空间校验：目标盘剩余需 > 源总量 + 10GB 裕量
    try:
        need = sum(_dir_stat(os.path.join(PLATFORM_ROOT, s))[1]
                   for s in _DB_SUBDIRS
                   if os.path.isdir(os.path.join(PLATFORM_ROOT, s)))
        os.makedirs(to_dir, exist_ok=True)
        free = shutil.disk_usage(to_dir).free
    except OSError as e:
        return {'ok': False, 'error': f'目标目录不可用: {to_dir} ({e})'}
    if free < need + 10e9:
        return {'ok': False, 'error':
                f'目标盘剩余 {_fmt(free)}，不足（需 ≥ 源 {_fmt(need)} + 10GB 裕量）'}

    _log(logger, f'数据库迁移 → {to_dir}（模式: {mode}，源约 {_fmt(need)}）')
    if dry_run:
        return {'ok': True, 'dry_run': True,
                'detail': f'预检通过：目标剩余 {_fmt(free)}，将迁移 '
                          f'{_fmt(need)}（{_DB_SUBDIRS}）'}

    # 逐个复制
    copied = 0
    for sub in _DB_SUBDIRS:
        s = os.path.join(PLATFORM_ROOT, sub)
        if not os.path.isdir(s):
            _log(logger, f'    跳过（不存在）: {sub}')
            continue
        d = os.path.join(to_dir, sub)
        # robocopy 在 Linux 不存在 → 回退 shutil（主要面向 Windows，保留此分支）
        if _platform.system() == 'Windows':
            ok = _robocopy(s, d, logger)
        else:
            ok = _copy_fallback(s, d, logger)
        if not ok:
            return {'ok': False, 'error': f'复制失败: {sub} → {d}'
                    '（源目录未动，配置未改，可重试断点续传）'}
        copied += 1

    # 校验：文件数与总字节一致（抽样全比）
    _log(logger, '  校验复制结果…')
    mism = []
    for sub in _DB_SUBDIRS:
        s = os.path.join(PLATFORM_ROOT, sub)
        d = os.path.join(to_dir, sub)
        if not os.path.isdir(s):
            continue
        sn, ss = _dir_stat(s)
        dn, ds = _dir_stat(d)
        if sn != dn or ss != ds:
            mism.append(f'{sub}: 源 {sn}文件/{_fmt(ss)} ≠ 目标 {dn}文件/{_fmt(ds)}')
    if mism:
        return {'ok': False, 'error': '校验不一致: ' + '; '.join(mism)
                + '（配置未改，可用 robocopy 断点续传重跑）'}
    _log(logger, f'  ✔ 校验一致（{copied} 个目录）')

    if mode == 'move':
        # 复制校验通过后删源
        for sub in _DB_SUBDIRS:
            s = os.path.join(PLATFORM_ROOT, sub)
            if os.path.isdir(s):
                _log(logger, f'  删除源（move 模式）: {s}')
                try:
                    shutil.rmtree(s)
                except OSError as e:
                    return {'ok': False, 'error':
                            f'复制已完成但删除源失败: {s} ({e})。'
                            '配置未切换；数据已在目标位置，可手动删源后重跑。'}

    # 写配置
    from .config import cfg
    try:
        cfg.set_database_root(to_dir)
    except (OSError, ValueError) as e:
        return {'ok': False, 'error': f'复制成功但配置切换失败: {e}'
                '（可手动到「设置 → 数据库目录」填入该路径）'}

    _log(logger, f'✔ 数据库根已切换: {to_dir}')
    _log(logger, '  现在可用，检查一下: python main.py selfcheck')
    return {'ok': True, 'detail': f'已切换数据库根到 {to_dir}'}


def check(to_dir):
    """只校验目标目录是否已有完整数据库（迁移/对接前检查）。"""
    from .config import PLATFORM_ROOT, DIRS
    to_dir = (to_dir or '').strip().strip('"')
    if not os.path.isabs(to_dir):
        return {'ok': False, 'error': '需绝对路径'}
    missing, ready = [], 0
    for sub in _DB_SUBDIRS:
        p = os.path.join(to_dir, sub)
        if os.path.isdir(p) and os.listdir(p):
            n, s = _dir_stat(p)
            print(f'  ✔ {sub}: {n} 文件 / {_fmt(s)}')
            ready += 1
        else:
            missing.append(sub)
    if missing:
        return {'ok': False, 'ready': ready,
                'error': f'缺少: {", ".join(missing)}（尚未对接，'
                '可用 db-migrate --to 复制本机库，或拷入数据库包）'}
    return {'ok': True, 'ready': ready, 'detail': '目标位置已有完整数据库，可直接对接'}

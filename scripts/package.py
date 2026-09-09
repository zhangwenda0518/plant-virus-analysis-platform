# -*- coding: utf-8 -*-
"""一键打包发布版：PyInstaller exe + 外部工具 + 前端（不含数据库，默认）。

用法:
  python scripts/package.py                          # 软件包（无数据库，~1GB）
  python scripts/package.py --with-db                # 软件包 + 预置最小病毒库/宿主库
                                                    #（开箱即用，体积 ~4GB）

软件与数据库分离：databases/（约 47GB）不进软件包，分发到目标机后用
  python main.py db-migrate --to D:\\库位置            # 从本机复制/对接
（或把另一台机器的 databases/ + host-db/ + virus-db 拷入后，在「设置 →
数据库目录」指向其父目录即可）。

产物: dist/VirusPlatform/
"""
import argparse
import os
import shutil
import subprocess
import sys

# 参数最先解析：--help 直接退出，绝不触发后面的打包/清理动作
_parser = argparse.ArgumentParser(description='打包发布版（默认不含数据库）')
_parser.add_argument('--with-db', action='store_true',
                     help='同时预置病毒库/宿主库运行文件（开箱即用）')
_parser.add_argument('--db-only', action='store_true',
                     help='跳过 exe，仅把数据库预置到已有发布目录（供重打库）')
_args, _ = _parser.parse_known_args()

# 本脚本位于 scripts/，ROOT = 平台根目录（app.py / vp 所在处）
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

# 预清理发布目录（git 打包文件只读，交给系统 rmdir 处理；PyInstaller 删不动）
import subprocess as _sp
_app_dir = os.path.join('dist', 'VirusPlatform')
if os.path.isdir(_app_dir):
    _sp.run(['cmd', '/c', 'rmdir', '/s', '/q', _app_dir], check=False)
    if os.path.isdir(_app_dir):
        sys.exit('无法清理旧发布目录（被占用？）——请关闭其中的 exe 后重试')

print('== 1/4 PyInstaller 打包 ==')
r = subprocess.run([sys.executable, '-m', 'PyInstaller',
                    os.path.join('scripts', 'VirusPlatform.spec'),
                    '--noconfirm', '--distpath', 'dist', '--workpath', 'build'])
if r.returncode != 0:
    sys.exit('打包失败')
shutil.rmtree('build', ignore_errors=True)

APP = os.path.join('dist', 'VirusPlatform')

# ---- 数据库预置（默认跳过：软件与数据库分离分发） ---------------------
if _args.db_only:
    _db_dir = _app_dir
else:
    _db_dir = _app_dir if _args.with_db else None

if _db_dir:
    print('== 2/4 预置数据库（运行文件） ==')
    import glob as _glob

    def _copy_db(src_dir, names_pat):
        dst = os.path.join(_db_dir, 'databases', os.path.basename(src_dir))
        os.makedirs(dst, exist_ok=True)
        n = 0
        pats = names_pat if isinstance(names_pat, (list, tuple)) else [names_pat]
        for pat in pats:
            for p in _glob.glob(os.path.join(src_dir, pat)):
                os.makedirs(dst, exist_ok=True)
                shutil.copy2(p, os.path.join(dst, os.path.basename(p)))
                n += 1
                print('  +', f'{os.path.basename(src_dir)}/{os.path.basename(p)}')
        return n

    # 病毒库（分类必需）
    _copy_db(os.path.join('databases', 'virus', 'plant'),
             ('hash_*.k2d', 'taxo.k2d', 'opts.k2d', 'hash_config.k2d',
              'seqid2taxid.map'))
    # 宿主库运行文件（已建好则预置；不含 library/prep 原料）
    _copy_db(os.path.join('databases', 'host', 'classify'),
             ('hash_*.k2d', 'taxo.k2d', 'opts.k2d', 'hash_config.k2d'))
    # palmdb（Open-Virome RdRP 分析参考库）
    palm_src = os.path.join('databases', 'palmdb')
    palm_dst = os.path.join(_db_dir, 'databases', 'palmdb')
    if os.path.isdir(palm_src):
        if os.path.isdir(palm_dst):
            shutil.rmtree(palm_dst)
        shutil.copytree(palm_src, palm_dst)
        print('  + palmdb/')
    # taxonomy 三件套（宿主建库/校验用；存在才拷）
    _copy_db(os.path.join('databases', 'tax', 'core'),
             ('nodes.dmp', 'names.dmp', 'merged.dmp'))
    print('  （数据库已预置 —— 注意：这是"内嵌库"分发方式，'
          '与 db-migrate 的外置库方式二选一）')
else:
    print('== 2/4 数据库：跳过（软件与数据库分离分发） ==')
    print('  软件包不含 databases/。目标机就绪方式：')
    print('    ① python main.py db-migrate --to <库目录>   （从本机复制）')
    print('    ② 或拷入 databases/host-db/virus-db 后在设置页填数据库目录')

print('== 3/4 外部工具与前端 ==')
# bin/（单文件工具整体）与 tools/（工具套件整体）——分发机无需 PATH 配置
for d in ('bin', 'tools'):
    if os.path.isdir(d):
        dst = os.path.join(APP, d)
        if os.path.isdir(dst):
            shutil.rmtree(dst)
        shutil.copytree(d, dst,
                        ignore=shutil.ignore_patterns('.git', '.github'))
        print('  +', d + '/（全部，已排除 .git）')

# orfipy.exe 放 dist 根（冻结模式下 vp/config.py 按 exe 同级目录探测）
orfipy = os.path.join(os.path.dirname(sys.executable), 'Scripts', 'orfipy.exe')
if os.path.isfile(orfipy):
    shutil.copy2(orfipy, os.path.join(APP, 'orfipy.exe'))
    print('  + orfipy.exe')

# Open-Virome 前端构建（app.py _FRONTEND_BUILD 按 <平台根>/open-virome 探测）
ov_build = os.path.join('open-virome', 'frontend', 'build')
if os.path.isdir(ov_build):
    dst = os.path.join(APP, 'open-virome', 'frontend', 'build')
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    shutil.copytree(ov_build, dst)
    print('  + open-virome/frontend/build/')

# 生成干净的 platform.json：不含任何本机绝对路径与已建库路径，
# 首次运行由 vp/config.py 自动探测（bin/ tools/ 随包；databases 由 db-migrate 对接）
import json as _json
fresh_cfg = {'threads': 0, 'tools': {},
             'defaults': {'max_heavy_tasks': 2, 'max_light_tasks': 4},
             'language': 'zh'}
if _args.db_only and os.path.isfile(os.path.join(APP, 'platform.json')):
    try:
        with open(os.path.join(APP, 'platform.json'), encoding='utf-8') as f:
            old = _json.load(f)
        old.setdefault('defaults', {}).setdefault('max_heavy_tasks', 2)
        old.setdefault('defaults', {}).setdefault('max_light_tasks', 4)
        old['tools'] = {}
        fresh_cfg = old
    except Exception:
        pass
_json.dump(fresh_cfg, open(os.path.join(APP, 'platform.json'), 'w',
                           encoding='utf-8'), ensure_ascii=False, indent=2)
print('  + platform.json（干净配置，路径零硬编码，含并发闸门缺省值）')

# 随包分发一份对接说明
_readme_txt = os.path.join(APP, '数据库对接说明.txt')
if not _args.db_only and os.path.isdir(_app_dir):
    with open(_readme_txt, 'w', encoding='utf-8') as f:
        f.write('本版本为「软件包」（不含数据库，体积小、易分发）。\n'
                '首次使用需对接数据库（二选一）：\n'
                '\n'
                '方式一（本机迁移，推荐）：在开发机上执行\n'
                '  python main.py db-migrate --to <目标库目录>\n'
                '  （把 databases/host-db/virus-db 复制到目标盘并自动配置）\n'
                '\n'
                '方式二（接收分发）：把库包解压/拷到任意目录（含三个子目录\n'
                '  databases/、host-db/、virus-db/），然后启动平台 →\n'
                '  「设置 → 数据库目录」填该目录 → 点应用（会校验库就绪状态）。\n'
                '\n'
                '验证：python main.py selfcheck\n')
    print('  + 数据库对接说明.txt')

print('== 4/4 完成 ==')
total = sum(os.path.getsize(os.path.join(dp, f))
            for dp, _, fs in os.walk(APP) for f in fs)
print(f'发布目录: {os.path.abspath(APP)}  ({total / 1e9:.2f} GB)')
if not _args.with_db and not _args.db_only:
    print('提示：本包不含数据库；需要"开箱即用"请加 --with-db 重打，'
          '或对已有发布目录跑 python scripts/package.py --db-only')

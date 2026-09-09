# -*- mode: python ; coding: utf-8 -*-
# spec 位于 scripts/，平台根目录 = SPECPATH 的上级（由 scripts/package.py 调用）
import os
ROOT = os.path.abspath(os.path.join(SPECPATH, '..'))

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

datas = [(os.path.join(ROOT, 'webapp'), 'webapp'),
         (os.path.join(ROOT, 'README.md'), '.')]
datas += collect_data_files('plotly')
datas += collect_data_files('pycirclize')
datas += collect_data_files('gbdraw')


a = Analysis(
    [os.path.join(ROOT, 'app.py')],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=['vp.public_meta.search_engine',
                   'vp.public_meta.info_engine',
                   'vp.logan_submit', 'vp.lovis4u_run',
                   'vp.ncbi_submit.unified_metadata', 'vp.engine_entry',
                   'gbdraw', 'vp.public_meta.landscape_plot', 'seaborn']
                   + collect_submodules('gbdraw'),
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='VirusPlatform',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='VirusPlatform',
)

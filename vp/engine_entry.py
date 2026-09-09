# -*- coding: utf-8 -*-
"""冻结分发（PyInstaller exe）的内部引擎入口。

exe 以 `VirusPlatform.exe --run-engine <脚本名> [args...]` 方式调用时，
app.py 顶部把控制权转交本模块，用 runpy 在进程内以 __main__ 语义执行
对应引擎（与源码模式下 `python 某引擎.py args...` 等价）。
引擎脚本必须在 VirusPlatform.spec 的 hiddenimports 中登记才会被打包。
"""
import runpy
import sys

_ENGINES = {
    'search_engine.py': 'vp.public_meta.search_engine',
    'info_engine.py': 'vp.public_meta.info_engine',
    'logan_submit.py': 'vp.logan_submit',
    'lovis4u_run.py': 'vp.lovis4u_run',
    'unified_metadata.py': 'vp.ncbi_submit.unified_metadata',
}


def run_engine(script_name):
    mod = _ENGINES.get(script_name)
    if not mod:
        print(f'未知引擎脚本: {script_name}', file=sys.stderr)
        return 2
    argv = [mod] + sys.argv[3:]
    try:
        sys.argv = argv
        runpy.run_module(mod, run_name='__main__', alter_sys=True)
        return 0
    except SystemExit as e:                      # 引擎主动退出码透传
        return e.code if isinstance(e.code, int) else 0
    except ImportError as e:                     # 可选依赖未打包（如 lovis4u）
        print(f'引擎 {script_name} 不可用（依赖未打包）: {e}', file=sys.stderr)
        return 3

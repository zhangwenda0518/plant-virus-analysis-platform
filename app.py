# -*- coding: utf-8 -*-
"""植物病毒分析平台 - 本地 Web GUI（组合根）。

启动: python app.py 或双击 启动平台.bat → 自动打开浏览器 http://127.0.0.1:8765
仅监听本机回环地址；长任务在后台线程执行，页面实时轮询进度。
所有文件写入经由 utils.safe_open（路径校验+限平台内），os.startfile 前一律 check_path。

本文件只做「组装」——创建 Flask 应用、注册 blueprint、启动服务。
业务逻辑分两层：
  - vp/        计算与流程（pipeline / kunpeng / phylo / orf_annot …）
  - vp/web/    HTTP 边界与任务工厂（各 blueprint）
**请勿在 app.py 新增路由**；新路由写进 vp/web/ 对应模块，在此登记 blueprint。
"""
import os
import sys
import threading
import time
import webbrowser

from flask import Flask

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 冻结分发：`VirusPlatform.exe --run-engine <脚本名> [args...]` → 进程内执行引擎
if getattr(sys, 'frozen', False) and len(sys.argv) > 2 and sys.argv[1] == '--run-engine':
    from vp.engine_entry import run_engine
    raise SystemExit(run_engine(sys.argv[2]))

# 冻结分发下 multiprocessing（SDT 精确矩阵的 ProcessPoolExecutor 等）子进程
# 会重新启动 exe，必须在任何业务代码前调用 freeze_support 拦截子进程引导；
# 非冻结环境是空操作
if __name__ == '__main__':
    import multiprocessing
    multiprocessing.freeze_support()

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# webapp 目录 / 配置单例 / tool_runs 根统一由 vp/web/state.py 提供，
# 各 blueprint 也从那里取，避免 app ↔ blueprint 循环导入。
from vp.web.state import _WWW  # noqa: E402

app = Flask(__name__, template_folder=os.path.join(_WWW, 'templates'),
            static_folder=os.path.join(_WWW, 'static'))
app.config['JSON_AS_ASCII'] = False
# 模板热重载：debug=False 下 Jinja2 默认缓存模板，改 .html 后需重启才生效；
# 打开后每次请求按 mtime 检测，本地平台开销可忽略，改模板即刷新可见。
app.config['TEMPLATES_AUTO_RELOAD'] = True

# ------------------------------------------------------------------
# Blueprint 注册
# 路由清单守卫：python tests/_route_inventory.py
#（锁定 179 条 rule + methods，拆分/新增后必须零差异）
# ------------------------------------------------------------------
from vp.web import (  # noqa: E402
    build as _bp_build,
    download as _bp_download,
    io_api as _bp_io,
    logan as _bp_logan,
    meta as _bp_meta,
    pages as _bp_pages,
    refs as _bp_refs,
    samples as _bp_samples,
    settings as _bp_settings,
    submit as _bp_submit,
    tasks_api as _bp_tasks,
    tool_results as _bp_tool_results,
    tools_api as _bp_tools,
    virome as _bp_virome,
)

# 顺序与拆分前的 app.py 区块顺序一致（便于对照回溯）
for _bp in (_bp_pages, _bp_tool_results, _bp_refs, _bp_logan, _bp_io,
            _bp_samples, _bp_download, _bp_submit, _bp_meta,
            _bp_settings, _bp_virome, _bp_tools, _bp_tasks, _bp_build):
    app.register_blueprint(_bp.bp)
del _bp


# ------------------------------------------------------------------
# 启动
# ------------------------------------------------------------------
def _pick_port():
    """选择可绑定的端口。

    Windows 上 Hyper-V/WSL 会动态保留端口段（netsh 可见），
    固定端口可能被系统排除导致 bind 被拒（访问权限不允许），
    因此逐个探测候选端口，全部失败再随机尝试。
    """
    import socket
    candidates = [8765, 8900, 8989, 9600, 8888, 5050, 5000]
    for p in candidates:
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('127.0.0.1', p))
                return p
        except OSError:
            continue
    import random
    for _ in range(30):
        p = random.randint(10000, 19999)
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('127.0.0.1', p))
                return p
        except OSError:
            continue
    raise RuntimeError('未找到可用的本地端口，请关闭占用端口的程序后重试')


def _open_browser(url):
    time.sleep(1.2)
    try:
        webbrowser.open(url)
    except Exception:
        pass


def main():
    try:
        port = _pick_port()
    except RuntimeError as e:
        print(f'  [错误] {e}')
        input('按回车键退出...')
        return
    url = f'http://127.0.0.1:{port}'
    print('=' * 56)
    print('  植物病毒分析平台 GUI 启动中...')
    print(f'  浏览器访问: {url}')
    if port != 8765:
        print(f'  （默认端口 8765 被系统占用/保留，已自动改用 {port}）')
    print('  关闭本窗口即退出平台')
    print('=' * 56)
    threading.Thread(target=_open_browser, args=(url,), daemon=True).start()
    try:
        app.run(host='127.0.0.1', port=port, debug=False, threaded=True)
    except OSError as e:
        print(f'\n  [错误] 端口 {port} 监听失败: {e}')
        print('  可在管理员 PowerShell 运行以下命令查看被系统保留的端口段:')
        print('     netsh interface ipv4 show excludedportrange protocol=tcp')
        input('按回车键退出...')


if __name__ == '__main__':
    main()

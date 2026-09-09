# -*- coding: utf-8 -*-
"""Web 层共享状态与路径解析——app.py 拆分的「单一真源」。

为什么需要它：app.py 里 `cfg` / `_WWW` / `_tool_runs_root` 等模块级对象
被上百个路由共用。一旦把路由搬进独立模块，若各模块反过来 `import app`
取这些对象，就会形成 app → blueprint → app 的循环导入。因此统一收到
本模块，由 app.py 与各 blueprint 双向导入。

内容：
  - cfg         配置单例（vp.config.get_config() 本身即单例）
  - _WWW        webapp 目录（模板/静态资源根）
  - webapp_dir() 目录探测（源码运行 / PyInstaller 冻结分发两种布局）
  - tool_runs_root()  tool_runs 根（跟随自定义输出根）
"""
import os
import sys

from vp.config import DIRS, PLATFORM_ROOT, get_config


def webapp_dir():
    """定位 webapp 目录。

    源码运行：平台根/webapp（app.py 与 vp/ 同级）；
    PyInstaller 冻结分发：优先内置资源 `sys._MEIPASS/webapp`。
    """
    # vp/web/state.py -> vp/web -> vp -> 平台根
    root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    src = os.path.join(root, 'webapp')
    if os.path.isdir(src):
        return src
    if getattr(sys, 'frozen', False):
        base = getattr(sys, '_MEIPASS',
                       os.path.dirname(os.path.abspath(sys.executable)))
        bundled = os.path.join(base, 'webapp')
        if os.path.isdir(bundled):
            return bundled
    return src


_WWW = webapp_dir()

# 配置单例：vp.config.get_config() 内部缓存，重复调用返回同一对象，
# 因此 app.py 与各 blueprint 拿到的是同一个 cfg（改语言/参数全局生效）。
cfg = get_config()


def tool_runs_root():
    """tool_runs 根目录（设置自定义输出根后由 DIRS 决定）。"""
    return DIRS.get('tool_runs') or os.path.join(PLATFORM_ROOT, 'tool_runs')


# ------------------------------------------------------------------
# 批处理队列单例（SampleQueue 在 app.py 中实例化，2D 批次迁出）
# 用「注册 / 取用」而非 `import app`，避免 blueprint ↔ app 循环导入。
# ------------------------------------------------------------------
_sample_queue = None


def set_sample_queue(queue):
    """由 app.py 在创建 sample_queue 后调用。"""
    global _sample_queue
    _sample_queue = queue


def get_sample_queue():
    """供 blueprint 取用（download 的「转入分析流程」用）。"""
    return _sample_queue

# -*- coding: utf-8 -*-
"""Open-Virome 公共病毒组模块（平台壳页 + 同源 iframe 内嵌 SPA）。

自 app.py 拆出。三条路由：
    GET /virome              平台壳页（顶部导航 + 内嵌应用）
    GET /virome/             SPA 本体（iframe 的 src）
    GET /virome/<path>       静态资源

目录优先级：源码运行时直接用 open-virome/frontend/build（重新构建即时
生效，无需拷贝）；打包部署回退到 webapp/static/virome 副本。
"""
import os

from flask import Blueprint, render_template, send_from_directory

from vp.config import PLATFORM_ROOT
from vp.web.state import _WWW

bp = Blueprint('virome', __name__)

_VIROME_DIST = os.path.join(_WWW, 'static', 'virome')
_FRONTEND_BUILD = os.path.join(PLATFORM_ROOT, 'open-virome',
                               'frontend', 'build')
if os.path.isdir(os.path.join(_FRONTEND_BUILD, 'static')):
    _VIROME_DIST = _FRONTEND_BUILD


@bp.route('/virome')
def page_virome():
    """平台内模块页：顶部平台导航 + 内嵌 Open-Virome 应用。

    查询参数（filters/palmprintOnly）在壳页 URL 与内嵌应用间双向同步，
    因此 /virome?filters=... 可直接收藏、分享、或由程序构造。
    """
    return render_template('virome.html')


@bp.route('/virome/')
def page_virome_app():
    """SPA 本体（iframe 的 src）。"""
    return send_from_directory(_VIROME_DIST, 'index.html')


@bp.route('/virome/<path:filename>')
def virome_asset(filename):
    return send_from_directory(_VIROME_DIST, filename)

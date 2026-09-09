# -*- coding: utf-8 -*-
"""任务状态 API（自 app.py 拆出）。

/api/tasks 列表与详情、日志、SSE 流、取消/删除/清理。"""
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid

from flask import (Blueprint, Response, abort, jsonify, render_template, request,
                   send_file, send_from_directory)

from vp.config import DIRS, PLATFORM_ROOT, db_path, engine_cmd
from vp.utils import (TaskLogger, check_path, fmt_size, run_cmd, safe_open,
                      safe_remove)
from vp.web.common import _safe_sample
from vp.web.state import cfg, tool_runs_root as _tool_runs_root
from vp.web.tasks import tm

bp = Blueprint('tasks_api', __name__)


@bp.route('/api/tasks')
def api_tasks():
    # logs=0：只要元信息不要日志（任务中心徒轮询用，日志走单独端点）
    logs = request.args.get('logs') not in ('0', 'false', 'no')
    if request.args.get('full'):
        arch, arch_total = tm._read_archived(logs=logs)
        return jsonify({'active': tm.list_all(logs=logs),
                        'archived': arch,
                        'archived_total': arch_total})
    return jsonify(tm.list_all(logs=logs))


@bp.route('/api/task/<tid>')
def api_task(tid):
    try:
        lines = min(int(request.args.get('log_lines') or 80), 1000)
    except (TypeError, ValueError):
        lines = 80
    snap = tm.snapshot(tid, log_lines=lines)
    if not snap:
        abort(404)
    return jsonify(snap)


@bp.route('/api/task/<tid>/log')
def api_task_log(tid):
    """任务全量日志（尾部 N 行）：内存任务 / 归档任务统一入口。"""
    try:
        lines = min(int(request.args.get('lines') or 400), 2000)
    except (TypeError, ValueError):
        lines = 400
    log = tm.full_log(tid, lines=lines)
    if log is None:
        abort(404, '日志不存在')
    return jsonify({'log': log})


@bp.route('/api/task/<tid>/stream')
def api_task_stream(tid):
    """SSE 实时任务进度流：快照有变化立即推送，任务结束自动关流。

    前端 EventSource 订阅（LOGAN 批量面板用），替代定时轮询。"""
    def gen():
        last = None
        idle = 0
        while True:
            snap = tm.snapshot(tid, log_lines=10)
            if snap is None:
                yield 'event: gone\ndata: {}\n\n'
                return
            data = json.dumps(snap, ensure_ascii=False)
            if data != last:
                last = data
                idle = 0
                yield f'data: {data}\n\n'
            else:
                idle += 1
                if idle >= 30:                     # ~20s 心跳注释行
                    idle = 0
                    yield ': ping\n\n'
            if snap.get('status') != 'running':
                return
            time.sleep(0.7)

    return Response(gen(), mimetype='text/event-stream',
                    headers={'Cache-Control': 'no-cache',
                             'X-Accel-Buffering': 'no'})


@bp.route('/api/task/<tid>/cancel', methods=['POST'])
def api_task_cancel(tid):
    return jsonify({'ok': tm.cancel(tid)})


@bp.route('/api/task/<tid>/delete', methods=['POST'])
def api_task_delete(tid):
    ok, msg = tm.delete(tid)
    if not ok:
        abort(400, msg)
    return jsonify({'ok': True})


@bp.route('/api/tasks/clear_finished', methods=['POST'])
def api_tasks_clear_finished():
    return jsonify({'removed': tm.clear_finished()})

# -*- coding: utf-8 -*-
"""公共数据下载中心（自 app.py 拆出）。

SRR/ERR/DRR/CRR 或 URL → aria2c 批量下载 → 一键转入分析流程。"""
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid

from flask import (Blueprint, abort, jsonify, render_template, request,
                   send_file, send_from_directory)

from vp.config import DIRS, PLATFORM_ROOT, db_path, engine_cmd
from vp.utils import (TaskLogger, check_path, fmt_size, run_cmd, safe_open,
                      safe_remove)
from vp.web.state import cfg, tool_runs_root as _tool_runs_root
from vp.web.tasks import tm
from vp.web.common import _safe_sample
from vp.web import state as _state


def _queue():
    """批处理队列单例（app.py 启动时注册，避免循环导入）。"""
    return _state.get_sample_queue()

bp = Blueprint('download', __name__)


@bp.route('/download')
def page_download():
    return render_template('download.html')


@bp.route('/api/dl/create', methods=['POST'])
def api_dl_create():
    body = request.get_json(force=True) or {}
    name = str(body.get('name') or '').strip() or 'batch'
    items = body.get('items') or []
    items = [x.strip() for x in items if x and x.strip()]
    if not items:
        abort(400, '请提供要下载的编号/URL 列表')
    if len(items) > 500:
        abort(400, '单批最多 500 条（大批次请拆分）')
    bid = get_dl_manager().create(
        name, items,
        concurrency=int(body.get('concurrency') or 2),
        convert_sra=bool(body.get('convert_sra', True)))
    return jsonify({'batch': bid})


@bp.route('/api/dl/batches')
def api_dl_batches():
    return jsonify(get_dl_manager().list_snapshots())


@bp.route('/api/dl/batch/<bid>')
def api_dl_batch(bid):
    snap = get_dl_manager().snapshot(bid)
    if not snap:
        abort(404, '批次不存在')
    return jsonify(snap)


@bp.route('/api/dl/batch/<bid>/<action>', methods=['POST'])
def api_dl_batch_action(bid, action):
    """批次操作：cancel / retry / delete（记录+文件全删）/
    delete_files（只删文件，留记录与日志）/ delete_record（只删记录，留文件）。"""
    m = get_dl_manager()
    if action == 'cancel':
        ok = m.cancel(bid)
    elif action == 'retry':
        n = m.retry_failed(bid)
        return jsonify({'ok': bool(n), 'retried': n})
    elif action == 'delete':
        ok = m.delete(bid)
    elif action == 'delete_files':
        ok = m.delete_files(bid)
    elif action == 'delete_record':
        ok = m.delete_record(bid)
    else:
        abort(400, '无效操作')
    if not ok:
        abort(400, '操作失败（批次可能不存在）')
    return jsonify({'ok': True})


def _dl_manager():
    from vp import public_data
    return public_data.get_manager()


get_dl_manager = _dl_manager


@bp.route('/downloads/<bid>/<path:filename>')
def page_dl_file(bid, filename):
    """下载批次内文件（严格限定在该批次目录内）。"""
    if not re.fullmatch(r'[A-Za-z0-9_\-.]+', bid) or '..' in filename:
        abort(400, '无效路径')
    from vp.public_data import _batch_root
    p = check_path(os.path.join(_batch_root(), bid, filename),
                   must_exist=True, in_platform=True)
    bdir = check_path(os.path.join(_batch_root(), bid),
                      must_exist=True, in_platform=True)
    if not p.startswith(bdir + os.sep):
        abort(400, '路径越界')
    return send_file(p, as_attachment=request.args.get('dl') == '1')


@bp.route('/api/dl/to_pipeline', methods=['POST'])
def api_dl_to_pipeline():
    """把下载完成的 run 转成本地样品（可选加入批处理队列）。"""
    body = request.get_json(force=True) or {}
    bid = body.get('batch') or ''
    runs = body.get('runs') or []
    project = str(body.get('project') or '') or None
    enqueue = bool(body.get('enqueue'))
    if not bid or not re.fullmatch(r'[A-Za-z0-9_\-.]+', bid):
        abort(400, '无效批次 ID')
    m = get_dl_manager()
    ready = m.ready_files(bid)
    if runs:
        ready = {k: v for k, v in ready.items() if k in runs}
    if not ready:
        abort(400, '该批次没有可用的 FASTQ 文件（可能仍在下载或需要 .sra 转换）')
    from vp.pipeline import _safe_sample_name, save_sample_input
    created, skipped = [], []
    for acc, files in ready.items():
        sample = _safe_sample_name(acc)
        sd = check_path(os.path.join(DIRS['results'], sample),
                        must_exist=False, in_platform=True)
        if os.path.isdir(sd) and any(os.scandir(sd)):
            skipped.append(sample)
            continue
        os.makedirs(sd, exist_ok=True)
        r1 = files.get('r1') or files.get('single')
        r2 = files.get('r2')
        if not r1:
            skipped.append(sample)
            continue
        save_sample_input(sd, r1, r2, sample, project=project)
        created.append(sample)
    queued = 0
    if enqueue and created:
        _queue().add(created, None, {}, project=project)
        queued = len(created)
    return jsonify({'created': created, 'skipped': skipped, 'queued': queued})

# -*- coding: utf-8 -*-
"""LOGAN 溯源（自 app.py 拆出）。

病毒序列 → Logan-Search → 物种/样本溯源；含 Selenium 批量提交。"""
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

bp = Blueprint('logan', __name__)


@bp.route('/logan')
def page_logan():
    return render_template('logan.html')


@bp.route('/logan/report/<name>/')
def page_logan_report(name):
    from vp.logan_trace import _job_dir
    p = os.path.join(_job_dir(name), 'trace_report.html')
    if not os.path.isfile(p):
        abort(404, '溯源报告尚未生成（请先导入结果文件）')
    return send_file(check_path(p, must_exist=True, in_platform=True))


@bp.route('/logan/report/<name>/<path:filename>')
def page_logan_report_file(name, filename):
    """溯源报告目录内静态文件（plotly.min.js / 原始结果表 / 查询 FASTA）。"""
    from vp.logan_trace import _job_dir
    p = check_path(os.path.join(_job_dir(name), filename),
                   must_exist=True, in_platform=True)
    return send_file(check_path(p, must_exist=True, in_platform=True))


@bp.route('/api/logan/samples')
def api_logan_samples():
    """可选来源样品：已有 ③组装 病毒 contigs 产出的样品。"""
    from vp.logan_trace import list_query_samples
    return jsonify(list_query_samples())


@bp.route('/api/logan/contigs/<sample>')
def api_logan_contigs(sample):
    from vp.logan_trace import list_virus_contigs
    return jsonify(list_virus_contigs(_safe_sample(sample)))


@bp.route('/api/logan/jobs')
def api_logan_jobs():
    from vp.logan_trace import list_jobs
    return jsonify(list_jobs())


@bp.route('/api/logan/job/<name>')
def api_logan_job(name):
    from vp.logan_trace import job_detail
    return jsonify(job_detail(name))


@bp.route('/api/logan/job/<name>/segment/<int:idx>')
def api_logan_segment(name, idx):
    """片段序列（前端复制提交用）。"""
    from vp.logan_trace import get_segment_fasta
    header, seq = get_segment_fasta(name, idx)
    return jsonify({'header': header, 'seq': seq})


@bp.route('/api/logan/job/<name>/query')
def api_logan_query_dl(name):
    """下载查询 FASTA（全部片段）。"""
    from vp.logan_trace import _job_dir
    p = os.path.join(_job_dir(name), 'query_all.fasta')
    if not os.path.isfile(p):
        abort(404, '查询 FASTA 不存在')
    return send_file(check_path(p, must_exist=True, in_platform=True),
                     as_attachment=True)


@bp.route('/api/logan/create', methods=['POST'])
def api_logan_create():
    from vp.logan_trace import create_job
    body = request.get_json(force=True)
    name = (body.get('name') or '').strip()
    if not name:
        abort(400, '请填写查询名称')
    pasted = (body.get('fasta') or '').strip()
    sample = (body.get('sample') or '').strip()
    contigs = body.get('contigs') or []
    if not pasted and not (sample and contigs):
        abort(400, '请勾选样品的病毒 contigs，或粘贴序列')
    try:
        return jsonify(create_job(name, sample=_safe_sample(sample) or None,
                                  contig_ids=contigs or None,
                                  pasted=pasted or None,
                                  n_seg=int(body.get('n_seg') or 2)))
    except ValueError as e:
        abort(400, str(e))


@bp.route('/api/logan/job/<name>/import/<int:idx>', methods=['POST'])
def api_logan_import(name, idx):
    """导入某片段的 Logan-Search 结果表（CSV/TSV），解析并自动生成报告。"""
    from vp.logan_trace import import_result, MAX_UPLOAD_BYTES
    f = request.files.get('file')
    if not f:
        abort(400, '缺少结果文件')
    raw = f.read(MAX_UPLOAD_BYTES + 1)
    if len(raw) > MAX_UPLOAD_BYTES:
        abort(400, '结果文件超过 64MB 上限')
    if not raw.strip():
        abort(400, '结果文件为空')
    try:
        return jsonify(import_result(name, idx, f.filename or '', raw))
    except ValueError as e:
        abort(400, str(e))


@bp.route('/api/logan/job/<name>/delete', methods=['POST'])
def api_logan_delete(name):
    from vp.logan_trace import delete_job
    delete_job(name)
    return jsonify({'ok': True})


@bp.route('/api/logan/batch_ready')
def api_logan_batch_ready():
    """批量模式可用性（selenium 是否已安装）。"""
    from vp.logan_trace import selenium_ready, GROUP_HINTS
    return jsonify({'selenium': selenium_ready(), 'groups': GROUP_HINTS})


@bp.route('/api/logan/job/<name>/batch', methods=['POST'])
def api_logan_batch(name):
    """后台批量提交任务的全部未导入片段，完成后自动导入并生成报告。"""
    from vp.logan_trace import batch_submit, _load_job
    _jdir, meta, segs = _load_job(name)
    if not any(not s.get('imported') for s in segs):
        abort(400, '该任务所有片段均已导入结果')
    for t in tm.list_all():
        if t.get('status') == 'running' and t.get('name') in {
            f'LOGAN 批量 {meta["name"]}',
            f'LOGAN batch {meta["name"]}',
        }:
            abort(400, f'任务 {meta["name"]} 已有批量提交在运行，请等它结束')
    body = request.get_json(force=True) or {}
    emails = [e.strip() for e in (body.get('emails') or '').split(',') if e.strip()]
    if not emails:
        abort(400, '请至少填写一个通知邮箱（用于接收 Logan-Search 完成通知）')
    group = body.get('group') or 'Fast_No_human'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        batch_submit(meta['name'], emails, group=group,
                     headless=bool(body.get('headless', True)),
                     first_wait=int(body.get('first_wait') or 300),
                     max_wait=int(body.get('max_wait') or 1800),
                     logger=logger, progress=prog, cancel=cancel)
        logger.close()
        return 'ok'

    tid = tm.start(cfg.tr(f'LOGAN 批量 {meta["name"]}',
                         f'LOGAN batch {meta["name"]}'), job,
                   weight='light')
    return jsonify({'task': tid})

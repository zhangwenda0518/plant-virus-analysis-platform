# -*- coding: utf-8 -*-
"""输入输出杂项 API（自 app.py 拆出）。

粘贴序列、拖拽上传、序列查看器、打开平台目录。"""
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
from vp.web.common import _safe_sample
from vp.web.state import cfg, tool_runs_root as _tool_runs_root
from vp.web.tasks import tm

bp = Blueprint('io_api', __name__)


@bp.route('/api/paste_input', methods=['POST'])
def api_paste_input():
    """粘贴序列文本 → 写入 uploads/paste_<ts>.<ext>，返回平台相对路径。

    body: {text: "FASTA/FASTQ 文本", ext: ".fasta" | ".fastq" | ".tsv"}
    所有工具的文件输入框均可使用返回的 path。
    """
    body = request.get_json(force=True) or {}
    text = (body.get('text') or '').strip()
    ext = (body.get('ext') or '.fasta').strip().lower()
    if not text:
        abort(400, '粘贴内容为空')
    if not re.fullmatch(r'\.(fasta|fa|fna|fas|fastq|fq|tsv|txt)', ext):
        abort(400, '不支持的文件类型')
    # 基本格式校验
    first_line = text.split('\n', 1)[0].strip()
    if ext in ('.fasta', '.fa', '.fna', '.fas') and not first_line.startswith('>'):
        abort(400, 'FASTA 格式应以 > 开头')
    if ext in ('.fastq', '.fq') and not first_line.startswith('@'):
        abort(400, 'FASTQ 格式应以 @ 开头')
    ts = time.strftime('%Y%m%d_%H%M%S')
    fname = f'paste_{ts}{ext}'
    fdir = os.path.join(PLATFORM_ROOT, 'uploads')
    os.makedirs(fdir, exist_ok=True)
    fp = os.path.join(fdir, fname)
    with safe_open(fp, 'wt') as f:
        f.write(text + '\n')
    return jsonify({'path': f'uploads/{fname}', 'size': len(text)})


@bp.route('/api/upload', methods=['POST'])
def api_upload():
    """拖拽/选择上传数据文件到平台 uploads/ 目录（流式落盘，支持大文件）。

    返回 {'path': 'uploads/<文件名>'} 供前端填入输入框。"""
    f = request.files.get('file')
    if not f:
        abort(400, '缺少上传文件')
    name = os.path.basename(f.filename or '')
    name = re.sub(r'[^\w.\-\u4e00-\u9fff]+', '_', name).strip('._')
    if not name:
        abort(400, '文件名无效')
    updir = check_path(DIRS.get('uploads')
                       or os.path.join(PLATFORM_ROOT, 'uploads'),
                       must_exist=False, in_platform=True)
    os.makedirs(updir, exist_ok=True)
    dst = check_path(os.path.join(updir, name), must_exist=False,
                     in_platform=True)
    if os.path.isfile(dst):                      # 重名不覆盖：追加序号
        base, ext = os.path.splitext(name)
        i = 1
        while os.path.isfile(check_path(
                os.path.join(updir, f'{base}_{i}{ext}'),
                must_exist=False, in_platform=True)):
            i += 1
        dst = check_path(os.path.join(updir, f'{base}_{i}{ext}'),
                         must_exist=False, in_platform=True)
    f.save(str(dst))
    up_disp = (os.path.abspath(dst) if DIRS.get('uploads')
               and not DIRS['uploads'].startswith(PLATFORM_ROOT)
               else 'uploads/' + os.path.basename(dst))
    return jsonify({'path': up_disp,
                    'size': os.path.getsize(dst)})


_DATA_EXTS = ('.fasta', '.fa', '.fna', '.fas', '.ffn')


@bp.route('/api/seqview')
def api_seqview():
    """轻量序列查看器：流式统计 + 分页预览（FASTA / FASTA.gz）。"""
    rel = request.args.get('path') or ''
    page = max(0, int(request.args.get('page', 0) or 0))
    per = 50
    # 平台内相对路径或任意绝对路径均可（只读流式统计，无写风险）；
    # 与文件浏览对话框（可选任意盘文件）行为对齐。
    if os.path.isabs(rel):
        p = check_path(rel, must_exist=True)
    else:
        p = check_path(os.path.join(PLATFORM_ROOT, rel), must_exist=True,
                       in_platform=True)
    if not p.lower().endswith(_DATA_EXTS +
                              tuple(e + '.gz' for e in _DATA_EXTS)):
        abort(400, '仅支持 FASTA / FASTA.gz')
    rows, previews = [], []
    total_bp = 0
    truncated = False
    rec_i = -1

    def _push(h, seq):
        nonlocal rec_i, total_bp, truncated
        rec_i += 1
        seq = seq.upper()
        if not seq:
            return
        n_gc = seq.count('G') + seq.count('C')
        n_deg = sum(seq.count(c) for c in 'RYKMSWBDHVN')
        total_bp += len(seq)
        if rec_i >= 20000:
            truncated = True
            return
        row = {'id': h.split()[0] if h.split() else h[:30],
               'len': len(seq),
               'gc': round(n_gc * 100.0 / len(seq), 1),
               'deg': round(n_deg * 100.0 / len(seq), 1)}
        if page * per <= rec_i < page * per + per:
            row['preview'] = seq[:300]
        rows.append(row)

    try:
        with safe_open(p) as f:
            h, buf = None, []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if h is not None:
                        _push(h, ''.join(buf))
                    h, buf = line[1:], []
                elif h is not None and line:
                    buf.append(line)
            if h is not None:
                _push(h, ''.join(buf))
    except (OSError, ValueError) as e:
        abort(400, f'读取失败: {e}')
    page_rows = [r for r in rows if 'preview' in r]
    pages = max(1, (min(rec_i + 1, 20000) + per - 1) // per)
    return jsonify({'total': rec_i + 1, 'total_bp': total_bp,
                    'page': page, 'pages': pages, 'per': per,
                    'truncated': truncated, 'rows': page_rows})


@bp.route('/api/open_platform_dir')
def api_open_platform_dir():
    os.startfile(check_path(PLATFORM_ROOT, must_exist=True, in_platform=True))
    return jsonify({'ok': True})

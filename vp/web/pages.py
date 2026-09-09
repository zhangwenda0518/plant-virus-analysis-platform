# -*- coding: utf-8 -*-
"""页面路由与全局上下文（自 app.py 拆出）。

含：NAV_GROUPS 导航单一数据源、模板上下文处理器、统一错误处理、
19 个页面端点与报告静态文件服务。"""
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid

from werkzeug.exceptions import HTTPException

from flask import (Blueprint, abort, jsonify, render_template,
                   request, send_file, send_from_directory)

from vp.config import DIRS, PLATFORM_ROOT, db_path, engine_cmd
from vp.utils import (TaskLogger, check_path, fmt_size, run_cmd, safe_open,
                      safe_remove)
from vp.web.common import _safe_sample
from vp.web.state import (  # noqa: F401
    _WWW, cfg, tool_runs_root as _tool_runs_root)
from vp.web.tasks import tm

bp = Blueprint('pages', __name__)


@bp.app_context_processor
def _inject_lang():
    """所有模板可用 {{ lang }}（服务端配置语言，前端 i18n.js 据此初始化）。"""
    return {'lang': cfg.lang}


_ASSET_V = None


def _asset_version():
    """静态资源版本号（取 app.js / i18n.js / app.css 的 mtime 最大值），模板引用带 ?v= 破缓存。"""
    global _ASSET_V
    if _ASSET_V is None:
        mt = 0.0
        for _f in ('app.js', 'i18n.js', 'app.css'):
            try:
                mt = max(mt, os.path.getmtime(os.path.join(_WWW, 'static', _f)))
            except OSError:
                pass
        _ASSET_V = str(int(mt)) if mt else '1'
    return _ASSET_V


@bp.app_context_processor
def _inject_asset_v():
    return {'asset_v': _asset_version()}


NAV_GROUPS = [
    {'id': 'resource', 'label': '数据资源', 'path': '数据资源', 'items': [
        {'href': '/meta', 'title': '公共数据检索', 'desc': '检索公共样本 / 元数据'},
        {'href': '/virome', 'title': 'Open-Virome', 'desc': '公共病毒组浏览 / 导出'},
        {'href': '/build', 'title': '数据库构建', 'desc': 'Taxonomy / 宿主库 / 病毒库'},
    ]},
    {'id': 'sample', 'label': '样本处理', 'path': '样本处理', 'items': [
        {'href': '/download', 'title': '公共数据下载', 'desc': 'Run / URL / SRA → FASTQ，衔接样品与模块'},
        {'href': '/samples', 'title': '样品创建 / 批量导入', 'desc': '样品名 + R1/R2 / TSV 批量'},
        {'id': 't-convert', 'title': '格式转换', 'desc': 'sra2fastq / sra2fasta / fastq2fasta'},
        {'id': 't-fastp', 'title': '序列质控', 'desc': 'fastp 去接头 / 过滤 / 去重'},
        {'href': '/hostremoval', 'title': '宿主去除与序列提取', 'desc': 'kunpeng 宿主库分类'},
    ]},
    {'id': 'kvsuite', 'label': '病毒定量与共识', 'path': '病毒定量与共识', 'items': [
        {'id': 't-kvchain', 'title': '⚡ 一键分析（全流程）',
         'desc': '鉴定 → 过滤 → 共识 → 深度绘图 → 变异注释，一次跑完，产物收拢到汇总目录'},
        {'id': 't-kvsuite', 'title': '已知病毒识别与定量',
         'desc': 'minibwa 比对病毒库 → 鉴定/过滤 → 共识 → 深度绘图 → 变异注释'},
        {'id': 't-consensus', 'title': '共识序列分析',
         'desc': '比对到参考 → 共识序列构建（minibwa + viral_consensus）'},
        {'id': 't-variant', 'title': '病毒变异分析',
         'desc': 'BAM / VCF 直入 → bcftools caller（可选）→ SnpEff 注释 → 8 张变异图 + SNPGenie'},
    ]},
    {'id': 'virus', 'label': '病毒识别和分类分析', 'path': '病毒识别和分类分析', 'items': [
        {'id': 't-virchain', 'title': '⚡ 一键分析（②→③→④）',
         'desc': '提取序列组装 → 组装结果再鉴定 → 候选序列验证，可选先跑①提取病毒序列'},
        {'id': 't-identify', 'title': '病毒识别分类与提取', 'desc': 'FASTQ / FASTA / 去宿主 → 分类 + 病毒序列提取'},
        {'id': 't-assemble', 'title': '提取序列组装', 'desc': 'SPAdes metaviral / contigs 过滤'},
        {'id': 't-contigs', 'title': '组装结果再鉴定', 'desc': 'contig 二次分类 + 谱系注释'},
        {'id': 't-verify', 'title': '候选序列验证', 'desc': '宿主筛选 → 长度分流 → blastx/CDD 过滤 + 类病毒 blastn'},
        {'href': '/hostpredict', 'title': '宿主预测',
         'desc': 'ICTV 宿主概率级联 + NCBI 元数据交叉（独立模块）'},
    ]},
    {'id': 'annotate', 'label': '病毒注释分析', 'path': '病毒注释分析', 'items': [
        {'href': '/orf', 'title': 'ORF 预测',
         'desc': 'pyrodigal / pyrodigal-rv 基因预测（独立模块）'},
        {'href': '/annotation', 'title': '功能注释',
         'desc': 'DIAMOND / MMseqs2 / blastp 对照 RefSeq 病毒蛋白库 + 类别投票 + ICTV 映射（独立模块）'},
        {'id': 't-cdd', 'title': '保守域与功能元件注释（CDD）',
         'desc': 'NCBI CDD 保守域搜索（6-frame 翻译取最长 ORF ≥50aa）：结构域位置 / 边界，服务于基因组图谱与功能注释'},
        {'id': 't-hom', 'title': '保守域与同源比对（CDD · BLASTN · BLASTX）',
         'desc': 'NCBI nt / nr 同源比对（限病毒）+ CDD 保守域搜索：为 GenBank 提交注释提供依据'},
        {'href': '/genome', 'title': '基因组图谱',
         'desc': 'gbdraw / DFV 圈图 + 线图（独立模块）'},
        {'href': '/primer', 'title': '引物设计',
         'desc': 'primer3 全长分窗 / 保守区设计（独立模块）'},
    ]},
    {'id': 'compare', 'label': '比较基因组分析', 'path': '比较基因组分析', 'items': [
        {'id': 't-seqprep', 'title': '参考序列获取', 'desc': 'ICTV 科/属选择 或 accession / 检索式 → 下载整科整属序列（GenBank 集合 + FASTA 参考集）'},
        {'id': 't-align', 'title': '序列比对（MAFFT + trimAl）', 'desc': 'MAFFT 比对 + trimAl 清剪；彩色比对查看器支持查看与编辑，结果直接送建树 / SDT'},
        {'id': 't-treebuild', 'title': '进化树构建（科/属级）', 'desc': 'GenBank 集合（全基因组 / CDS / PEP）或 FASTA → MAFFT 比对 + NJ / FastTree / IQ-TREE 建树；页内树查看'},
        {'id': 't-sdt', 'title': 'SDT 同一性分析（属级）', 'desc': '逐对 MAFFT 精确比对 → identity 矩阵 / 热图 / 分布图；NT+AA 模式同一性表 + 复合热图'},
        {'id': 't-synteny', 'title': '同属共线性比较', 'desc': 'LoVis4u 基因组共线性图（MMseqs2 聚类）'},
    ]},
    {'id': 'result', 'label': '结果中心', 'path': '结果中心', 'items': [
        {'href': '/results', 'title': '样品结果 / 专项结果',
         'desc': '报告 / 专项运行 / SDT / MSA / 进化树 / 序列查看'},
    ]},
    {'id': 'trace', 'label': '溯源与提交', 'path': '溯源与提交', 'items': [
        {'href': '/logan', 'title': 'LOGAN 溯源 / 批量提交', 'desc': '公共样本追踪'},
        {'href': '/submit', 'title': '数据提交', 'desc': 'NCBI 提交准备'},
    ]},
]


_PATH_TO_GROUP = {'/meta': 'resource', '/virome': 'resource',
                  '/download': 'sample', '/build': 'resource',
                  '/hostremoval': 'sample', '/samples': 'sample',
                  '/hostpredict': 'virus', '/orf': 'annotate',
                  '/annotation': 'annotate',
                  '/genome': 'annotate', '/primer': 'annotate',
                  '/logan': 'trace', '/submit': 'trace'}


@bp.app_context_processor
def _inject_group_nav():
    """组导航数据 + 当前组：子页侧栏与 tools 工作台共用同一数据源。"""
    gid = None
    if request.path == '/tools':
        gid = request.args.get('g')
    if gid not in {g['id'] for g in NAV_GROUPS}:
        gid = _PATH_TO_GROUP.get(request.path)
    group = next((g for g in NAV_GROUPS if g['id'] == gid), None)
    return {'nav_groups': NAV_GROUPS, 'current_group': group}


@bp.app_errorhandler(ValueError)
@bp.app_errorhandler(FileNotFoundError)
def _bad_request(e):
    """路径非法/不存在 → 400（含路径穿越拦截）。"""
    return jsonify({'error': str(e)}), 400


@bp.app_errorhandler(PermissionError)
def _forbidden(e):
    return jsonify({'error': str(e)}), 403


@bp.app_errorhandler(HTTPException)
def _http_err(e):
    """所有 HTTP 错误统一 JSON 返回（abort(400, msg) 等）。"""
    return jsonify({'error': e.description}), e.code


@bp.route('/')
def page_index():
    return render_template('home.html')


@bp.route('/pipeline')
def page_pipeline():
    return render_template('pipeline.html')


@bp.route('/hostremoval')
def page_host_removal():
    return render_template('host_removal.html')


@bp.route('/samples')
def page_samples():
    return render_template('samples.html')


@bp.route('/hostpredict')
def page_hostpredict():
    return render_template('host_predict.html')


@bp.route('/orf')
def page_orf():
    return render_template('orf.html')


@bp.route('/annotation')
def page_annotation():
    return render_template('annotation.html')


@bp.route('/genome')
def page_genome():
    return render_template('genome.html')


@bp.route('/primer')
def page_primer():
    return render_template('primer.html')


@bp.route('/build')
def page_build():
    return render_template('build.html')


@bp.route('/tasks')
def page_tasks():
    """全局任务中心：运行中监测 / 日志折叠 / 停止 / 重启 / 删除 / 历史归档。"""
    return render_template('tasks.html')


@bp.route('/results')
def page_results():
    samples, archived = [], []
    res_root = check_path(DIRS['results'], must_exist=True, in_platform=True)
    for name in sorted(os.listdir(res_root)):
        if name.startswith('_'):
            # _archive 归档样品单独收集；其余内部目录不展示
            if name == '_archive':
                arch = os.path.join(res_root, '_archive')
                for an in sorted(os.listdir(arch)):
                    ad = os.path.join(arch, an)
                    if not os.path.isdir(ad):
                        continue
                    archived.append({
                        'name': an,
                        'report': os.path.isfile(
                            os.path.join(ad, '07_report', 'report.html')),
                        'mtime': time.strftime(
                            '%Y-%m-%d %H:%M',
                            time.localtime(os.path.getmtime(ad)))})
            continue
        d = check_path(os.path.join(res_root, name), must_exist=False,
                       in_platform=True)
        if not os.path.isdir(d):
            continue
        rpt = check_path(os.path.join(d, '07_report', 'report.html'),
                         must_exist=False, in_platform=True)
        samples.append({'name': name, 'report': os.path.isfile(rpt),
                        'mtime': time.strftime(
                            '%Y-%m-%d %H:%M',
                            time.localtime(os.path.getmtime(d)))})
    return render_template('results.html', samples=samples, archived=archived)


@bp.route('/report/<sample>/')
def page_report(sample):
    safe = _safe_sample(sample)
    rpt = check_path(os.path.join(DIRS['results'], safe, '07_report', 'report.html'),
                     must_exist=True, in_platform=True)
    return send_file(check_path(rpt, must_exist=True, in_platform=True))


@bp.route('/report/<sample>/<path:filename>')
def page_report_file(sample, filename):
    """报告目录内静态文件（plotly.min.js 等；尾斜杠路由使相对引用可解析）。"""
    safe = _safe_sample(sample)
    p = check_path(os.path.join(DIRS['results'], safe, '07_report', filename),
                   must_exist=True, in_platform=True)
    return send_file(check_path(p, must_exist=True, in_platform=True))


@bp.route('/archive_report/<sample>/')
def page_archive_report(sample):
    """归档样品报告（results/_archive/<sample>/）。"""
    safe = _safe_sample(sample)
    rpt = check_path(os.path.join(DIRS['results'], '_archive', safe,
                                  '07_report', 'report.html'),
                     must_exist=True, in_platform=True)
    return send_file(check_path(rpt, must_exist=True, in_platform=True))


@bp.route('/help')
def page_help():
    """浏览器内阅读使用手册（不依赖系统 .md 文件关联）。"""
    cands = [os.path.join(PLATFORM_ROOT, 'README.md')]
    if getattr(sys, '_MEIPASS', None):
        cands.append(os.path.join(sys._MEIPASS, 'README.md'))
    text = ''
    for p in cands:
        if os.path.isfile(p):
            try:
                with safe_open(p) as f:
                    text = f.read()
            except OSError:
                pass
            break
    return render_template('help.html', text=text)

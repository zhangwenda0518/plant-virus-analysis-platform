# -*- coding: utf-8 -*-
"""
植物病毒分析平台 - 本地 Web GUI
启动: python app.py 或双击 启动平台.bat  →  自动打开浏览器 http://127.0.0.1:8765
仅监听本机回环地址；长任务在后台线程执行，页面实时轮询进度。
所有文件写入经由 utils.safe_open（路径校验+限平台内），os.startfile 前一律 check_path。
"""
import os
import re
import sys
import json
import time
import logging
import uuid
import threading
import subprocess
import webbrowser
from collections import deque

from flask import (Flask, jsonify, request, render_template, send_file,
                   send_from_directory, abort, Response)
from werkzeug.exceptions import HTTPException

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

from vp.config import DIRS, PLATFORM_ROOT, engine_cmd, db_path
from vp.utils import check_path, safe_open, safe_remove, TaskLogger, fmt_size

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# webapp 目录 / 配置单例 / tool_runs 根统一由 vp/web/state.py 提供：
# app.py 拆分后各 blueprint 也从此处取，避免 app ↔ blueprint 循环导入。
from vp.web.state import (  # noqa: E402
    _WWW, cfg, tool_runs_root as _tool_runs_root)

app = Flask(__name__, template_folder=os.path.join(_WWW, 'templates'),
            static_folder=os.path.join(_WWW, 'static'))
app.config['JSON_AS_ASCII'] = False
# 模板热重载：debug=False 下 Jinja2 默认缓存模板，改 .html 后需重启才生效；
# 打开后每次请求按 mtime 检测，本地平台开销可忽略，改模板即刷新可见。
app.config['TEMPLATES_AUTO_RELOAD'] = True

# ------------------------------------------------------------------
# Blueprint 注册（app.py 单体拆分的模块化路由）
# 每拆一批在此登记；路由清单守卫：python tests/_route_inventory.py
# ------------------------------------------------------------------
from vp.web import download as _bp_download    # noqa: E402
from vp.web import logan as _bp_logan          # noqa: E402
from vp.web import meta as _bp_meta            # noqa: E402
from vp.web import settings as _bp_settings    # noqa: E402
from vp.web import submit as _bp_submit        # noqa: E402
from vp.web import tools_api as _bp_tools      # noqa: E402
from vp.web import virome as _bp_virome        # noqa: E402
from vp.web.common import _safe_sample         # noqa: E402
from vp.web.tasks import tm                    # noqa: E402
from vp.web.tool_jobs import _find_contig_seq  # noqa: E402

app.register_blueprint(_bp_download.bp)
app.register_blueprint(_bp_logan.bp)
app.register_blueprint(_bp_meta.bp)
app.register_blueprint(_bp_settings.bp)
app.register_blueprint(_bp_submit.bp)
app.register_blueprint(_bp_tools.bp)
app.register_blueprint(_bp_virome.bp)


@app.context_processor
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


@app.context_processor
def _inject_asset_v():
    return {'asset_v': _asset_version()}


# 一级模块组 → 二级模块（单一数据源：tools 工作台与各二级页侧栏共用，
# 与 docs/DEVELOPMENT_NOTES.md 的模块 YAML 总表一致）
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


@app.context_processor
def _inject_group_nav():
    """组导航数据 + 当前组：子页侧栏与 tools 工作台共用同一数据源。"""
    gid = None
    if request.path == '/tools':
        gid = request.args.get('g')
    if gid not in {g['id'] for g in NAV_GROUPS}:
        gid = _PATH_TO_GROUP.get(request.path)
    group = next((g for g in NAV_GROUPS if g['id'] == gid), None)
    return {'nav_groups': NAV_GROUPS, 'current_group': group}


# ------------------------------------------------------------------
# 任务引擎 / 结果预览 → vp/web/tasks.py（TaskManager、tm、build_result_preview）
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# 页面
# ------------------------------------------------------------------
@app.errorhandler(ValueError)
@app.errorhandler(FileNotFoundError)
def _bad_request(e):
    """路径非法/不存在 → 400（含路径穿越拦截）。"""
    return jsonify({'error': str(e)}), 400


@app.errorhandler(PermissionError)
def _forbidden(e):
    return jsonify({'error': str(e)}), 403


@app.errorhandler(HTTPException)
def _http_err(e):
    """所有 HTTP 错误统一 JSON 返回（abort(400, msg) 等）。"""
    return jsonify({'error': e.description}), e.code


@app.route('/')
def page_index():
    return render_template('home.html')


@app.route('/pipeline')
def page_pipeline():
    return render_template('pipeline.html')


@app.route('/hostremoval')
def page_host_removal():
    return render_template('host_removal.html')


@app.route('/samples')
def page_samples():
    return render_template('samples.html')


@app.route('/hostpredict')
def page_hostpredict():
    return render_template('host_predict.html')


@app.route('/orf')
def page_orf():
    return render_template('orf.html')


@app.route('/annotation')
def page_annotation():
    return render_template('annotation.html')


@app.route('/genome')
def page_genome():
    return render_template('genome.html')


@app.route('/primer')
def page_primer():
    return render_template('primer.html')


@app.route('/build')
def page_build():
    return render_template('build.html')


@app.route('/tasks')
def page_tasks():
    """全局任务中心：运行中监测 / 日志折叠 / 停止 / 重启 / 删除 / 历史归档。"""
    return render_template('tasks.html')


@app.route('/results')
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


@app.route('/report/<sample>/')
def page_report(sample):
    safe = _safe_sample(sample)
    rpt = check_path(os.path.join(DIRS['results'], safe, '07_report', 'report.html'),
                     must_exist=True, in_platform=True)
    return send_file(check_path(rpt, must_exist=True, in_platform=True))


@app.route('/report/<sample>/<path:filename>')
def page_report_file(sample, filename):
    """报告目录内静态文件（plotly.min.js 等；尾斜杠路由使相对引用可解析）。"""
    safe = _safe_sample(sample)
    p = check_path(os.path.join(DIRS['results'], safe, '07_report', filename),
                   must_exist=True, in_platform=True)
    return send_file(check_path(p, must_exist=True, in_platform=True))


# ------------------------------------------------------------------
# Open-Virome / 工具箱 / 设置与存储 / 任务引擎 已拆到 vp/web/ 包：
#   virome.py    Open-Virome 壳页与静态资源
#   tool_jobs.py 22 个工具任务工厂
#   chains.py    链式一键分析（virchain / kvchain）
#   tools_api.py 工具箱 HTTP 边界（/tools、/api/tool/*）
#   settings.py  设置中心 + 存储水位
#   tasks.py     TaskManager 与结果预览
# blueprint 在文件顶部「Blueprint 注册」区统一登记。
# ------------------------------------------------------------------


def _task_matches_sample(task_name, sample):
    name = str(task_name or '')
    sample = str(sample or '')
    if not name or not sample:
        return False
    if name.startswith(sample + ' · '):
        return True
    # 批处理队列的任务名为「队列·<样品>」/「Queue·<样品>」（分隔符无空格），
    # 与上面的「<样品> · 」不是同一套写法，必须单独匹配，
    # 否则队列运行中删除/清空样品时 _sample_busy 会漏判。
    if re.match(rf'^(?:队列|Queue)·{re.escape(sample)}$', name):
        return True
    return bool(re.match(rf'^(?:样品分析|Sample analysis)\s+{re.escape(sample)}(?:\b|$)',
                         name))


# ------------------------------------------------------------------
# 工具箱注册表：每个工具声明 标题 / 是否需病毒库 / job 构造器。
# 新增工具三步：① 实现 _tool_job_<name>(ctx) 返回 job(log, prog, cancel)；
#              ② 在 TOOL_REGISTRY 登记；③ 前端 tools.html 加卡片 + MODULE_NAV。
# ctx 字段：p(参数) run_dir threads db_virus req(key,what) opt(key)
# ------------------------------------------------------------------


# 共识模块可选 reads 来源：层级目录名。
# 默认 host_removed（去宿主后全量）——未经过病毒相似度筛选，变异检测无偏。


# ------------------------------------------------------------------
# 链式一键分析（chain jobs）
# 设计：每步仍调用现有 _tool_job_*，各建自己的 run_dir，一步不改；
# 汇总目录只放「符号链接 + 清单」，把一次多步运行的结果收拢成一份。
# 汇总目录：(tool_runs_root)/<chain>_<ts>/
#   _chain.json            每步的 run 名 / 产物 / 状态
#   stepN_<name>           指向该步 run_dir 的符号链接（目录）
#   关键产物文件链接        指向该步产物文件，便于直接在汇总目录取用
# ------------------------------------------------------------------


# 轻量工具（走 light 闸门，可与 heavy 任务并行）：纯 IO 或秒级计算，
# 不占满多核/大内存。其余一律 heavy（kunpeng 分类 / SPAdes / DIAMOND /
# 建树 / SDT 等吃满线程的任务，同时最多跑 max_heavy_tasks 个）。

# 工具②/④可选的内置 kunpeng 病毒库（键名 → 目录）。
# 主库 virus 由「数据库构建」页构建；refvirus/rvdb/k2viral 由
# vp/universal_ref 或外部预构建（/api/dbs 探测就绪状态）。


def _analysis_out(run_dir, contig, action):
    """分析结果 JSON 的固定路径：run/analysis/<安全化contig>_<action>.json。

    run/contig/action 先做严格字符校验 + 规范化，杜绝 ../ 穿越。
    """
    if action not in ('blastn', 'blastx', 'cdd', 'primer'):
        abort(400, '无效的分析类型')
    if not re.fullmatch(r'[A-Za-z0-9_\-]{1,80}', os.path.basename(run_dir)):
        abort(400, '无效的运行名')
    adir = check_path(os.path.join(run_dir, 'analysis'), must_exist=False,
                      in_platform=True)
    os.makedirs(adir, exist_ok=True)
    safe_c = re.sub(r'[^A-Za-z0-9_\-.]', '_', contig)[:60]
    p = os.path.abspath(os.path.join(adir, f'{safe_c}_{action}.json'))
    if not p.startswith(os.path.abspath(adir) + os.sep):
        abort(400, '非法结果路径')
    return check_path(p, must_exist=False, in_platform=True)


@app.route('/api/tool/viral_contigs')
def api_tool_viral_contigs():
    """某次 contigs 运行的病毒 contig 分类表（metabuli 风格列）。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    tsv = check_path(os.path.join(_tool_runs_root(), run,
                                  'virus_classification.tsv'),
                     must_exist=False, in_platform=True)
    rows = []
    if os.path.isfile(tsv):
        import csv
        with safe_open(tsv) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                r['near_complete'] = r.get('near_complete') == 'True'
                rows.append(r)
    # 并入宿主预测列（若该运行已做宿主预测）
    hp = check_path(os.path.join(_tool_runs_root(), run,
                                 '08_host_analysis', 'host_prediction.tsv'),
                    must_exist=False, in_platform=True)
    host_map = {}                      # 先初始化：未做宿主预测的运行无该文件
    if os.path.isfile(hp):
        import csv
        host_map = {}
        with safe_open(hp) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                c = (r.get('contig') or '').strip()
                h = (r.get('final_host') or '').strip()
                conf = (r.get('confidence_level') or '').strip()
                host_map[c] = h + (f' ({conf})' if conf else '')
        for r in rows:
            r['host'] = host_map.get((r.get('contig') or '').strip(), '')
    return jsonify(rows)


@app.route('/api/tool/verify_result')
def api_tool_verify_result():
    """候选序列验证结果：calls.tsv 行 + summary.json 摘要。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    vdir = check_path(os.path.join(_tool_runs_root(), run, 'verify'),
                      must_exist=False, in_platform=True)
    rows, summary = [], {}
    calls_tsv = os.path.join(vdir, 'calls.tsv')
    if os.path.isfile(calls_tsv):
        import csv
        with safe_open(calls_tsv) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                r['len'] = int(r.get('len') or 0)
                rows.append(r)
    summ_json = os.path.join(vdir, 'summary.json')
    if os.path.isfile(summ_json):
        with safe_open(summ_json) as f:
            summary = json.load(f)
    return jsonify({'rows': rows, 'summary': summary})


@app.route('/api/tool/kvsuite_refs')
def api_tool_kvsuite_refs():
    """列出某次已知病毒运行中可用作回贴参考的病毒基因组。

    来源（按优先级）：
      <run>/kvsuite/virus-fasta/ref_<acc>/ref_<acc>.ref.fasta   过滤后的病毒基因组
      <run>/kvsuite/virus-annotations/<acc>.gb / .gtf         注释
    供 t-consensus 选参考、t-variant 查注释用。
    """
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    import csv as _csv
    base = check_path(os.path.join(_tool_runs_root(), run, 'kvsuite'),
                      must_exist=False, in_platform=True)

    def _read_ids(path):
        ids = []
        if not os.path.isfile(path):
            return ids
        with safe_open(path) as f:
            for line in f:
                if line.startswith('>'):
                    ids.append(line[1:].strip().split()[0])
        return ids

    fa_dir = os.path.join(base, 'virus-fasta')
    ann_dir = os.path.join(base, 'virus-annotations')
    refs = []
    if os.path.isdir(fa_dir):
        for dn in sorted(os.listdir(fa_dir)):
            if not dn.startswith('ref_'):
                continue
            sub = os.path.join(fa_dir, dn)
            if not os.path.isdir(sub):
                continue
            fasta = None
            for fn in sorted(os.listdir(sub)):
                if fn.endswith(('.fasta', '.fa', '.fna')):
                    fasta = os.path.join(sub, fn)
                    break
            if not fasta:
                continue
            acc = dn[len('ref_'):]
            try:
                size = os.path.getsize(fasta)
            except OSError:
                size = 0
            gb = os.path.join(ann_dir, f'{acc}.gb')
            gtf = os.path.join(ann_dir, f'{acc}.gtf')
            gtf_n = 0
            if os.path.isfile(gtf):
                with safe_open(gtf) as f:
                    gtf_n = sum(1 for line in f if line.strip())
            refs.append({
                'accession': acc,
                'dir': os.path.relpath(sub, base).replace('\\', '/'),
                'fasta': os.path.relpath(fasta, base).replace('\\', '/'),
                'size': fmt_size(size),
                'seq_ids': _read_ids(fasta),
                'gb': os.path.isfile(gb),
                'gtf': os.path.isfile(gtf),
                'n_cds': gtf_n,
            })

    # BAM：供 t-consensus 复用映射（映射模式=复用）
    bams = []
    bam_dir = os.path.join(base, 'bam')
    if os.path.isdir(bam_dir):
        for fn in sorted(os.listdir(bam_dir)):
            if not fn.endswith('.bam'):
                continue
            p = os.path.join(bam_dir, fn)
            try:
                nbytes = os.path.getsize(p)
            except OSError:
                continue
            bams.append({
                'name': fn,
                'path': p,
                'size': fmt_size(nbytes),
                'size_mb': round(nbytes / 1048576.0, 1),
                'indexed': os.path.isfile(p + '.bai'),
            })

    # 过滤表：给前端展示物种名（末次鉴定结果）
    filt = {}
    ft = os.path.join(base, 'filter', 'filtered.tsv')
    if os.path.isfile(ft):
        with safe_open(ft) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                acc = (r.get('Accession') or '').strip()
                if acc and acc not in filt:
                    filt[acc] = {
                        'species': (r.get('Species') or '').strip(),
                        'coverage': (r.get('Coverage(%)') or '').strip(),
                        'depth': (r.get('MeanDepth') or '').strip(),
                    }
    for r in refs:
        r['species'] = filt.get(r['accession'], {}).get('species', '')
        r['coverage'] = filt.get(r['accession'], {}).get('coverage', '')
        r['depth'] = filt.get(r['accession'], {}).get('depth', '')

    return jsonify({
        'run': run,
        'n_refs': len(refs),
        'refs': refs,
        'bams': bams,
        'has_annotations': os.path.isdir(ann_dir),
    })


@app.route('/api/tool/kvsuite_result')
def api_tool_kvsuite_result():
    """known_virus_suite 结果：鉴定 / 过滤 / 共识 QC / 变异注释 四张表 + 图文件。

    产物目录 <run>/kvsuite/{identify,filter,consensus,variant,plots}
    """
    import csv as _csv
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    base = check_path(os.path.join(_tool_runs_root(), run, 'kvsuite'),
                      must_exist=False, in_platform=True)

    def _read_tsv(path, limit=None):
        rows = []
        if not os.path.isfile(path):
            return rows
        with safe_open(path) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                rows.append(r)
                if limit and len(rows) >= limit:
                    break
        return rows

    ns = {'identify': [], 'filtered': [], 'discarded': [],
          'consensus_qc': [], 'variants': []}
    ns['identify'] = _read_tsv(os.path.join(base, 'identify',
                                            'all_viruses.summary.tsv'), 500)
    ns['filtered'] = _read_tsv(os.path.join(base, 'filter', 'filtered.tsv'), 500)
    ns['discarded'] = _read_tsv(os.path.join(base, 'filter', 'discarded.tsv'), 200)
    ns['consensus_qc'] = _read_tsv(
        os.path.join(base, 'consensus', 'consensus_qc.tsv'), 500)

    # 变异注释：把所有 per-genome ann.tsv 合起来
    # 注意：kv_variant 把产物直接写在 <out> 下（<out>/annotated/、
    # <out>/variant_summary.json），没有 variant/ 中间层。
    ann_dir = os.path.join(base, 'annotated')
    # QUAL/DP/AF 不在 ann.tsv 里，需从同名的 vcf 补上（键 = CHROM|POS|ALT）
    vcf_map = _kv_load_vcf_metrics(os.path.join(base, 'vcf'))
    if os.path.isdir(ann_dir):
        for fn in sorted(os.listdir(ann_dir)):
            if not fn.endswith('.ann.tsv'):
                continue
            gid = fn[: -len('.ann.tsv')]
            for r in _read_tsv(os.path.join(ann_dir, fn), 300):
                r['genome'] = gid
                k = '{}|{}|{}'.format(r.get('CHROM', ''), r.get('POS', ''),
                                     r.get('ALT', ''))
                m = vcf_map.get(k)
                if m:
                    r['QUAL'] = m['QUAL']
                    r['DP'] = m['DP']
                    r['AF'] = m['AF']
                ns['variants'].append(r)

    # 深度图（<out>/plots/<sample>/*.png 有子目录，必须递归）
    plots = []
    pdir = os.path.join(base, 'plots')
    if os.path.isdir(pdir):
        for root, _dirs, files in os.walk(pdir):
            for fn in sorted(files):
                if fn.lower().endswith(('.png', '.pdf')):
                    rel = os.path.relpath(os.path.join(root, fn), pdir)
                    plots.append(rel.replace(os.sep, '/'))
        plots.sort()

    # 变异图（<out>/variant_plots/，平铺目录）
    variant_plots = []
    vpdir = os.path.join(base, 'variant_plots')
    if os.path.isdir(vpdir):
        for fn in sorted(os.listdir(vpdir)):
            if fn.lower().endswith(('.png', '.pdf')):
                variant_plots.append(fn)
    # 只展示 png（页面内联），pdf 供下载
    variant_plots_png = [p for p in variant_plots if p.lower().endswith('.png')]

    # 扩展变异分析（分子谱 / 变异密度 / 群体遗传），与 variant_plots 同目录平铺
    evo_plots_png = [p for p in variant_plots
                     if p.lower().endswith('.png') and '_evo_' in p]
    evo_manifest = {}
    _eman = os.path.join(base, 'variant_evo', 'evo_manifest.json')
    if os.path.isfile(_eman):
        try:
            with safe_open(_eman) as f:
                evo_manifest = json.load(f)
        except (ValueError, OSError):
            evo_manifest = {}

    # 群体遗传滑窗表 + SNPGenie 汇总（只取头几行供页面预览）
    popgen_tables = {}
    _edir = os.path.join(base, 'variant_evo')
    if os.path.isdir(_edir):
        for fn in sorted(os.listdir(_edir)):
            if fn.endswith('_popgen_window.tsv'):
                acc = fn[: -len('_popgen_window.tsv')]
                popgen_tables[acc] = _read_tsv(os.path.join(_edir, fn), 200)

    snpgenie_products = {}
    if os.path.isdir(_edir):
        for fn in sorted(os.listdir(_edir)):
            if fn.endswith('.snpgenie'):
                acc = fn[: -len('.snpgenie')]
                # 优先读带 ORF 变异计数的注释版（区分无变异/无多态）
                for pname in ('product_results.annotated.txt',
                              'product_results.txt'):
                    prod = os.path.join(_edir, fn, pname)
                    if os.path.isfile(prod):
                        snpgenie_products[acc] = _read_tsv(prod, 200)
                        break

    summary = {}
    vsum = os.path.join(base, 'variant_summary.json')
    if os.path.isfile(vsum):
        try:
            with safe_open(vsum) as f:
                summary['variant'] = json.load(f)
            # 脱敏：results 里的 vcf/ann_vcf/ann_tsv 是绝对路径，
            # 前端渲染只用 accession/genome/n_variants/n_coding_ann，
            # 不把本机路径暴露到页面。产物文件本身保持完整（含路径）。
            var = summary.get('variant') or {}
            for r in (var.get('results') or []):
                for k in ('vcf', 'ann_vcf', 'ann_tsv'):
                    r.pop(k, None)
        except (ValueError, OSError):
            pass
    summary['n_identify'] = len(ns['identify'])
    summary['n_filtered'] = len(ns['filtered'])
    summary['n_discarded'] = len(ns['discarded'])
    summary['n_consensus'] = len(ns['consensus_qc'])
    summary['n_variants'] = len(ns['variants'])
    summary['n_variant_plots'] = len(variant_plots_png)
    summary['n_evo_plots'] = len(evo_plots_png)
    summary['n_snpgenie'] = len(evo_manifest.get('snpgenie') or [])

    return jsonify({'run': run, 'summary': summary, 'plots': plots,
                    'variant_plots': variant_plots_png,
                    'evo_plots': evo_plots_png,
                    'evo_manifest': evo_manifest,
                    'popgen_tables': popgen_tables,
                    'snpgenie_products': snpgenie_products,
                    'identify': ns['identify'], 'filtered': ns['filtered'],
                    'discarded': ns['discarded'],
                    'consensus_qc': ns['consensus_qc'],
                    'variants': ns['variants']})


def _kv_load_vcf_metrics(vcf_dir):
    """扫 <out>/vcf/*.vcf 取 QUAL/DP/AF，键为 CHROM|POS|ALT。

    AF 口径与 kv_variant 过滤一致：AD[1]/INFO/DP（无 AD 时退 AC/AN）。
    单行格式异常只跳过该行并 warning，不丢整个 VCF 剩余记录。
    """
    out = {}
    if not os.path.isdir(vcf_dir):
        return out
    for fn in sorted(os.listdir(vcf_dir)):
        if not fn.lower().endswith('.vcf') or fn.startswith('all.'):
            continue
        n_bad = 0
        try:
            with safe_open(os.path.join(vcf_dir, fn)) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    try:
                        f2 = line.rstrip('\n').split('\t')
                        if len(f2) < 8:
                            n_bad += 1
                            continue
                        chrom, pos, _id, _ref, alt, qual, _filt, info = f2[:8]
                        dp = af = None
                        ad = ac = an = None
                        for item in info.split(';'):
                            if item.startswith('DP='):
                                dp = item[3:]
                            elif item.startswith('AD='):
                                ad = item[3:]
                            elif item.startswith('AC='):
                                ac = item[3:]
                            elif item.startswith('AN='):
                                an = item[3:]
                        dp_i = None
                        try:
                            dp_i = int(dp) if dp is not None else None
                        except (ValueError, TypeError):
                            dp_i = None
                        alt_ad = None
                        if ad:
                            parts = ad.split(',')
                            if len(parts) >= 2:
                                try:
                                    alt_ad = int(parts[1])
                                except (ValueError, TypeError):
                                    alt_ad = None
                        if alt_ad is not None and dp_i:
                            af = round(alt_ad / dp_i, 4)
                        elif ac and an:
                            try:
                                af = round(float(ac) / float(an), 4)
                            except (ValueError, TypeError, ZeroDivisionError):
                                af = None
                        out['{}|{}|{}'.format(chrom, pos, alt)] = {
                            'QUAL': qual, 'DP': dp_i, 'AF': af}
                    except Exception as e:
                        n_bad += 1
                        if n_bad <= 3:
                            app.logger.warning(
                                f'VCF 行解析失败 {fn} L{line[:60]!r}: {e}')
        except OSError as e:
            app.logger.warning(f'VCF 读取失败 {fn}: {e}')
            continue
        if n_bad:
            app.logger.warning(f'VCF {fn}: {n_bad} 行解析失败已跳过')
    return out


@app.route('/api/tool/consensus_result')
def api_tool_consensus_result():
    """共识与变异结果：coverage.tsv + variants.tsv + summary.json + 分箱覆盖曲线。

    cov 为按 bin 平均的深度序列（最多 ~800 点），避免把上万个逐位深度
    全量塞给前端。
    """
    import csv as _csv
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    cdir = check_path(os.path.join(_tool_runs_root(), run, 'consensus'),
                      must_exist=False, in_platform=True)
    rows, variants, summary = [], [], {}

    cov_tsv = os.path.join(cdir, 'coverage.tsv')
    if os.path.isfile(cov_tsv):
        with safe_open(cov_tsv) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                for k in ('length', 'mapped_reads', 'covered_bp',
                          'n_variants', 'n_snv', 'n_isnv'):
                    try:
                        r[k] = int(r.get(k) or 0)
                    except ValueError:
                        r[k] = 0
                for k in ('median_depth', 'min_depth', 'max_depth'):
                    try:
                        r[k] = int(float(r.get(k) or 0))
                    except ValueError:
                        r[k] = 0
                for k in ('mean_depth', 'coverage_pct',
                          'depth_ok_pct', 'called_pct'):
                    try:
                        r[k] = float(r.get(k) or 0)
                    except ValueError:
                        r[k] = 0.0
                rows.append(r)

    var_tsv = os.path.join(cdir, 'variants.tsv')
    if os.path.isfile(var_tsv):
        with safe_open(var_tsv) as f:
            for i, r in enumerate(_csv.DictReader(f, delimiter='\t')):
                if i >= 500:
                    break
                variants.append(r)

    summ_json = os.path.join(cdir, 'summary.json')
    if os.path.isfile(summ_json):
        with safe_open(summ_json) as f:
            summary = json.load(f)

    cov = {}
    for r in rows:
        name = r.get('contig') or ''
        tag = re.sub(r'[^\w\-.]+', '_', name)[:120]
        cov[name] = _bin_coverage(os.path.join(cdir, tag + '.poscounts.tsv'),
                                  int(r.get('length') or 0))
    return jsonify({'rows': rows, 'variants': variants,
                    'summary': summary, 'cov': cov})


def _bin_coverage(pc_path, length, max_bins=800):
    """poscounts.tsv → 分箱平均深度序列。"""
    if not os.path.isfile(pc_path) or length <= 0:
        return []
    bin_sz = max(1, (length + max_bins - 1) // max_bins)
    nbin = (length + bin_sz - 1) // bin_sz
    sums = [0] * nbin
    cnts = [0] * nbin
    with safe_open(pc_path) as f:
        next(f, None)
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 8:
                continue
            try:
                pos = int(p[0])
                depth = int(p[1]) + int(p[2]) + int(p[3]) + int(p[4]) + int(p[5])
            except ValueError:
                continue
            b = pos // bin_sz
            if 0 <= b < nbin:
                sums[b] += depth
                cnts[b] += 1
    return [round(sums[i] / cnts[i], 1) if cnts[i] else 0 for i in range(nbin)]


@app.route('/api/tool/structcmp_data')
def api_tool_structcmp_data():
    """结构比较结果（JSON）：identity 矩阵 + 序列名，供页内 Plotly 热图渲染。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    js = check_path(os.path.join(_tool_runs_root(), run,
                                 'identity_matrix.json'),
                    must_exist=True, in_platform=True)
    return send_file(js, mimetype='application/json')


@app.route('/api/tool/msa_data')
def api_tool_msa_data():
    """structcmp 运行 → SNP-only MSA 展示数据（比较基因组 MSA 查看卡用）。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    aln = check_path(os.path.join(_tool_runs_root(), run, 'aln.fasta'),
                     must_exist=True, in_platform=True)
    from vp.msa_view import snp_view_data
    try:
        return jsonify(snp_view_data(aln))
    except ValueError as e:
        abort(400, str(e))


@app.route('/api/tool/quicktree_data')
def api_tool_quicktree_data():
    """quicktree 运行 → 树 Newick（进化树查看器卡片直接渲染）。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    run_dir = check_path(os.path.join(_tool_runs_root(), run),
                         must_exist=True, in_platform=True)
    tree_file = None
    for name in ('nj.nwk', 'tree.nwk'):
        p = os.path.join(run_dir, name)
        if os.path.isfile(p):
            tree_file = p
            break
    if not tree_file:
        abort(404, '该运行没有树文件（尚未完成或建树失败，看运行日志）')
    with safe_open(tree_file) as f:
        nwk = f.read().strip()
    if not nwk:
        abort(404, '树文件为空')
    base = os.path.basename(tree_file)
    return jsonify({'newick': nwk, 'file': base,
                    'tool': 'NJ' if base == 'nj.nwk' else 'FastTree'})


_ALIGN_EXTS = ('.fasta', '.fa', '.fna', '.fas', '.aln', '.txt')


def _parse_alignment_fasta(path):
    """读比对/序列 FASTA → dict（names/seqs/长度一致性/类型）。"""
    from vp.utils import iter_fasta
    from vp.sdt_exact import detect_seqtype
    names, seqs = [], []
    for h, s in iter_fasta(path):
        names.append(re.split(r'[\s|]', (h or '').strip())[0][:60] or
                     f'seq{len(names) + 1}')
        seqs.append(s.upper())
    if not seqs:
        raise ValueError('文件中没有序列')
    lens = {len(s) for s in seqs}
    return {'names': names, 'seqs': seqs, 'n': len(seqs),
            'cols': max(len(s) for s in seqs),
            'aligned': len(lens) == 1,
            'type': detect_seqtype(seqs)}


@app.route('/api/align/file')
def api_align_file():
    """比对查看器数据：path= 平台内或本机（只读）FASTA 比对文件。"""
    path = (request.args.get('path') or '').strip()
    if not path.lower().endswith(_ALIGN_EXTS):
        abort(400, '请提供 FASTA/比对文件（' + '/'.join(_ALIGN_EXTS[:4]) + '…）')
    try:
        p = check_path(path if os.path.isabs(path)
                       else os.path.join(PLATFORM_ROOT, path),
                       must_exist=True)
    except (ValueError, FileNotFoundError) as e:
        abort(400, f'路径不合法: {e}')
    try:
        d = _parse_alignment_fasta(p)
    except ValueError as e:
        abort(400, str(e))
    d['path'] = os.path.relpath(p, PLATFORM_ROOT) \
        if str(p).startswith(str(PLATFORM_ROOT)) else str(p)
    return jsonify(d)


@app.route('/api/align/save', methods=['POST'])
def api_align_save():
    """编辑保存：content 写到源文件旁 <stem>.edit.fasta（平台内）；
    源在平台外时落 tool_runs/align_edit_<ts>/edited.fa。返回新路径。"""
    body = request.get_json(force=True) or {}
    path = (body.get('path') or '').strip()
    content = str(body.get('content') or '')
    if not path or not content.strip():
        abort(400, '参数不完整（path 与 content 必填）')
    try:
        p = check_path(path if os.path.isabs(path)
                       else os.path.join(PLATFORM_ROOT, path),
                       must_exist=True)
    except (ValueError, FileNotFoundError) as e:
        abort(400, f'路径不合法: {e}')
    # 校验能按 FASTA 解析（至少 2 条、只含合法字符集宽容量）
    try:
        probe = _parse_alignment_fasta_string(content)
    except ValueError as e:
        abort(400, f'内容不是合法 FASTA: {e}')
    in_plat = str(p).startswith(str(PLATFORM_ROOT) + os.sep)
    if in_plat:
        stem, _ext = os.path.splitext(p)
        dst = stem + '.edit.fasta'
    else:
        dst = check_path(os.path.join(
            _tool_runs_root(), f"align_edit_{time.strftime('%Y%m%d_%H%M%S')}",
            'edited.fasta'), must_exist=False, in_platform=True)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
    with safe_open(dst, 'wt') as f:
        f.write(content if content.endswith('\n') else content + '\n')
    return jsonify({'saved': dst, 'aligned': probe['aligned'],
                    'n': probe['n'], 'cols': probe['cols']})


def _parse_alignment_fasta_string(content):
    """与 _parse_alignment_fasta 同口径，但直接吃文本。"""
    from vp.sdt_exact import detect_seqtype
    names, seqs, cur = [], [], None
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if cur is not None:
                seqs.append(''.join(cur))
            names.append(re.split(r'[\s|]', line[1:].strip())[0][:60] or
                         f'seq{len(names) + 1}')
            cur = []
        elif cur is not None:
            cur.append(line)
        else:
            raise ValueError('第一条记录前出现序列行')
    if cur is not None:
        seqs.append(''.join(cur))
    if len(seqs) < 1 or not names:
        raise ValueError('未解析到序列')
    joined = [s.upper() for s in seqs]
    if any(not s for s in joined):
        raise ValueError('存在空序列')
    return {'names': names, 'seqs': joined, 'n': len(joined),
            'cols': max(len(s) for s in joined),
            'aligned': len({len(s) for s in joined}) == 1,
            'type': detect_seqtype(joined)}


_TREE_EXTS = ('.nwk', '.newick', '.treefile', '.contree', '.tre', '.tree')


@app.route('/api/tree/file')
def api_tree_file():
    """树文件路径 → Newick 文本（进化树查看卡「本机树文件」输入用）。

    路径可以是平台内相对路径或本机绝对路径（只读），扩展名白名单限制。
    """
    path = (request.args.get('path') or '').strip()
    if not path.lower().endswith(_TREE_EXTS):
        abort(400, '请提供树文件（'
                    + ' / '.join(_TREE_EXTS) + '）')
    try:
        p = check_path(path if os.path.isabs(path)
                       else os.path.join(PLATFORM_ROOT, path),
                       must_exist=True)
    except (ValueError, FileNotFoundError) as e:
        abort(400, f'路径不合法: {e}')
    with safe_open(p) as f:
        nwk = f.read().strip()
    if not nwk:
        abort(404, '树文件为空')
    return jsonify({'newick': nwk, 'file': os.path.basename(p)})


@app.route('/api/tool/report')
def api_tool_report():
    """生成（带缓存）并打开工具运行的结果报告（桑基/旭日/分类表）。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    run_dir = check_path(os.path.join(_tool_runs_root(), run),
                         must_exist=True, in_platform=True)
    from vp.viz import build_contig_report
    try:
        build_contig_report(run_dir)
    except RuntimeError as e:
        abort(400, str(e))
    from flask import redirect
    return redirect(f'/tool_runs/{run}/report.html')


@app.route('/api/tool/report_data')
def api_tool_report_data():
    """分类报告数据（JSON）：供工具页内联渲染桑基/旭日/分类表/明细表。

    contigs_* 运行 → mode=contig（含 contig 明细 + 宿主列）；
    identify_* 运行 → mode=reads（含 classified sequences 表）。
    """
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    run_dir = check_path(os.path.join(_tool_runs_root(), run),
                         must_exist=True, in_platform=True)
    import glob as _glob
    krs = (_glob.glob(os.path.join(run_dir, '*.kreport2')) +
           _glob.glob(os.path.join(run_dir, '**', '*.kreport2'),
                      recursive=True))
    if not krs:
        abort(400, '运行目录中没有 kreport 分类报告')
    from vp.kunpeng import parse_kreport
    krows = parse_kreport(sorted(krs)[0])
    report_rows = [{'proportion': r['percent'], 'count': int(r['frags']),
                    'rank': r['rank'], 'taxid': r['taxid'],
                    'name': r['name'], 'depth': r['depth']}
                   for r in krows]

    import csv as _csv
    tsv = os.path.join(run_dir, 'virus_classification.tsv')
    ids_tsv = os.path.join(run_dir, 'viral_ids.tsv')
    mode, vc_rows, ids_rows = ('contig', [], [])
    if os.path.isfile(tsv):
        mode = 'contig'
        with safe_open(tsv) as f:
            vc_rows = list(_csv.DictReader(f, delimiter='\t'))
    elif os.path.isfile(ids_tsv):
        mode = 'reads'
        with safe_open(ids_tsv) as f:
            ids_rows = list(_csv.DictReader(f, delimiter='\t'))

    host_map, host_cats = {}, {}
    hp = os.path.join(run_dir, '08_host_analysis', 'host_prediction.tsv')
    if os.path.isfile(hp):
        with safe_open(hp) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                c = (r.get('contig') or '').strip()
                h = (r.get('final_host') or '').strip()
                conf = (r.get('confidence_level') or '').strip()
                if c:
                    host_map[c] = h + (f' ({conf})' if conf else '')
                if h:
                    host_cats[h] = host_cats.get(h, 0) + 1
    for r in vc_rows:
        r['host'] = host_map.get((r.get('contig') or '').strip(), '')

    n_total, n_virus = 0, 0
    if mode == 'contig':
        cfa = os.path.join(run_dir, 'contigs.filtered.fasta')
        if os.path.isfile(cfa):
            from vp.utils import count_fasta_seqs
            n_total = count_fasta_seqs(cfa)
        n_virus = len(vc_rows)
    else:
        n_total = n_virus = len(ids_rows)

    host_list = [{'cat': k, 'n': v} for k, v in
                 sorted(host_cats.items(), key=lambda x: -x[1])]

    # taxburst 旭日图（Krona 等价物，离线交互 HTML）
    from vp.viz import build_taxburst
    tb_ok, tb_unc = build_taxburst(run_dir, krows)

    dl = []
    if mode == 'contig':
        dl = [('viral_contigs.fasta?dl=1', 'Virus Sequences (FASTA)'),
              ('virus_classification.tsv?dl=1', 'Virus Classification (TSV)'),
              ('contigs.filtered.fasta?dl=1', 'All Contigs (FASTA)'),
              (os.path.basename(sorted(krs)[0]) + '?dl=1', 'Report (kreport TSV)')]
        if os.path.isfile(hp):
            dl.append(('08_host_analysis/host_prediction.tsv?dl=1',
                       'Host Prediction (TSV)'))
    else:
        import glob as _g2
        for fa in sorted(_g2.glob(os.path.join(run_dir, 'viral_sequences*.fasta'))):
            dl.append((os.path.basename(fa) + '?dl=1', 'Viral Sequences (FASTA)'))
        if os.path.isfile(ids_tsv):
            dl.append(('viral_ids.tsv?dl=1', 'Classified IDs (TSV)'))
        dl.append((os.path.basename(sorted(krs)[0]) + '?dl=1',
                   'Report (kreport TSV)'))
    downloads = [{'href': href, 'label': label} for href, label in dl
                 if os.path.isfile(os.path.join(run_dir, href.split('?')[0]))]

    return jsonify({'run': run, 'mode': mode,
                    'report': report_rows, 'vc': vc_rows, 'ids': ids_rows,
                    'host': host_list, 'host_count': len(host_map),
                    'n_total': n_total, 'n_virus': n_virus,
                    'downloads': downloads,
                    'taxburst': tb_ok})


@app.route('/api/tool/hostpredict_run', methods=['POST'])
def api_tool_hostpredict_run():
    """对既有 contigs 运行目录原地做 ICTV 宿主预测（结果并入该运行）。"""
    body = request.get_json(force=True) or {}
    run = body.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    run_dir = check_path(os.path.join(_tool_runs_root(), run),
                         must_exist=True, in_platform=True)
    tsv = os.path.join(run_dir, 'virus_classification.tsv')
    if not os.path.isfile(tsv):
        abort(400, '该运行没有 virus_classification.tsv（仅 contig 分类运行支持）')
    from vp.kunpeng import db_ready
    if not db_ready(cfg.databases['host']):
        abort(400, '宿主库未就绪，请先到「数据库构建」页构建宿主库')

    def job(log, prog, cancel):
        import shutil as _shutil
        logger = TaskLogger(callback=log)
        a_dir = check_path(os.path.join(run_dir, '03_assembly'),
                           must_exist=False, in_platform=True)
        os.makedirs(a_dir, exist_ok=True)
        prog('prep', 0.05, '整理输入')
        norm = os.path.join(a_dir, 'virus_contigs.tsv')
        with safe_open(tsv) as f:
            header = f.readline().rstrip('\r\n').split('\t')
        if 'kunpeng_taxid' in header:
            _shutil.copyfile(tsv, norm)
        else:
            cols = ['contig', 'length', 'kunpeng_flag', 'kunpeng_taxid',
                    'kunpeng_species', 'blast_top_hit', 'blast_identity(%)',
                    'blast_coverage_hsp(%)', 'blast_aln_len', 'blast_species',
                    'blast_family']
            import csv as _csv
            with safe_open(tsv) as f, safe_open(norm, 'wt') as w:
                w.write('\t'.join(cols) + '\n')
                for row in _csv.DictReader(f, delimiter='\t'):
                    kt = str(row.get('taxid') or '').strip()
                    w.write('\t'.join([
                        str(row.get('contig') or '').strip(),
                        str(row.get('length') or '').strip(),
                        'C' if kt.isdigit() and int(kt) > 0 else 'U', kt,
                        str(row.get('taxon') or '').strip(),
                        '', '', '', '', '', '']) + '\n')
        fa = os.path.join(run_dir, 'viral_contigs.fasta')
        if os.path.isfile(fa):
            _shutil.copyfile(fa, os.path.join(a_dir, 'viral_contigs.fasta'))
        from vp.host_analysis import predict_hosts
        prog('predict', 0.15, 'ICTV 宿主概率级联预测')
        res = predict_hosts(run_dir, logger=logger, force=True)
        prog('done', 1.0, f"完成：宿主判定 {res.get('n_contigs', 0)} 条")
        logger.close()
        return res

    ts = time.strftime('%Y%m%d_%H%M%S')
    tid = tm.start(cfg.tr(f'宿主预测·{run} {ts}', f'Host prediction·{run} {ts}'),
                   job, log_file=os.path.join(run_dir, 'hostpredict.log'))
    return jsonify({'task': tid})


@app.route('/api/tool/analysis')
def api_tool_analysis():
    """读取已完成的 contig 分析结果（JSON 缓存）。"""
    run = request.args.get('run') or ''
    contig = request.args.get('contig') or ''
    action = request.args.get('action') or ''
    seq = (request.args.get('seq') or '').strip()
    if not contig or action not in ('blastn', 'blastx', 'cdd', 'primer'):
        abort(400, '参数不完整')
    # 直接用字符串序列查询（seq 输入的伪 run 缓存）
    if seq:
        import hashlib as _hl
        _key = _hl.md5((contig + '|' + seq).encode('utf-8')).hexdigest()[:12]
        _safe = re.sub(r'[^A-Za-z0-9_\-.]', '_', contig)[:30]
        run_dir = os.path.join(_tool_runs_root(), '_seq_input', _safe + '_' + _key)
        run_dir = check_path(run_dir, must_exist=True, in_platform=True)
        p = _analysis_out(run_dir, _safe + '_' + _key, action)
    else:
        if not run or not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
            abort(400, '无效的运行名')
        run_dir = check_path(os.path.join(_tool_runs_root(), run),
                             must_exist=True, in_platform=True)
        p = _analysis_out(run_dir, contig, action)
    if not os.path.isfile(p):
        abort(404, '该分析尚未运行')
    with safe_open(p) as f:
        return app.response_class(f.read(), mimetype='application/json')


@app.route('/api/tool/analyze', methods=['POST'])
def api_tool_analyze():
    """单 contig 深度分析四件套：blastn / blastx / cdd / primer。

    body: {run, contig, action}
    NCBI 在线分析走 vp.contig_annot（域名白名单 + IP 校验），
    结果缓存为 run/analysis/<contig>_<action>.json。
    """
    body = request.get_json(force=True) or {}
    run = body.get('run') or ''
    contig = body.get('contig') or ''
    action = body.get('action') or ''
    seq = (body.get('seq') or '').strip()
    in_seq = seq  # 供闭包安全读取（避免 UnboundLocalError）
    if action not in ('blastn', 'blastx', 'cdd', 'primer'):
        abort(400, '无效的分析类型')
    if not run or not contig:
        abort(400, '参数不完整')
    # 直接输入序列模式（无需 run），或从运行目录取 contig 序列
    if seq:
        # 用共享伪运行缓存，避免污染具体工具运行目录
        run_dir = check_path(os.path.join(_tool_runs_root(), '_seq_input'),
                             must_exist=False, in_platform=True)
        import hashlib as _hl
        _key = _hl.md5((contig + '|' + seq).encode('utf-8')).hexdigest()[:12]
        contig = re.sub(r'[^A-Za-z0-9_\-.]', '_', contig)[:30]
        _used_for_cache = contig + '_' + _key
        run_dir = os.path.join(run_dir, _used_for_cache)
        os.makedirs(run_dir, exist_ok=True)
        out_path = _analysis_out(run_dir, _used_for_cache, action)
    else:
        if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
            abort(400, '无效的运行名')
        run_dir = check_path(os.path.join(_tool_runs_root(), run),
                             must_exist=True, in_platform=True)
        out_path = _analysis_out(run_dir, contig, action)
    if os.path.isfile(out_path):
        return jsonify({'task': None, 'cached': True})

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp import contig_annot as ca
        prog('prep', 0.1, '查找 contig 序列')
        if in_seq:
            # 清洗直接输入的序列：去掉 FASTA header 与换行，仅留核酸字母
            seq = ''.join(line.strip().upper() for line in in_seq.splitlines()
                          if line.strip() and not line.strip().startswith('>'))
            if not seq:
                raise RuntimeError('输入序列为空')
        else:
            seq = _find_contig_seq(run_dir, contig)
            if not seq:
                raise RuntimeError(f'运行目录中未找到 contig: {contig}')

        if action == 'primer':
            from vp.primer import design_primers_for_seq
            prog('primer3', 0.5, '本地 primer3 设计引物')
            pairs = design_primers_for_seq(contig, seq, num_return=5)
            result = {'action': 'primer', 'contig': contig,
                      'length': len(seq), 'primers': pairs}

        elif action == 'cdd':
            prog('orf', 0.3, '6-frame 翻译取最长 ORF')
            prot = ca.longest_orf_protein(seq)
            if len(prot) < 50:
                raise RuntimeError('未找到 ≥50aa 的 ORF，无法做 CDD 域搜索')
            prog('cdd', 0.5, '提交 NCBI CDD（数分钟，耐心等待）')
            cdsid = ca.submit_cdd(f'>{contig}\n{prot}\n')
            hits = ca.poll_cdd(cdsid, cancel=cancel)
            result = {'action': 'cdd', 'contig': contig,
                      'orf_len_aa': len(prot), 'hits': hits}

        else:  # blastn / blastx（NCBI URL API，virus-restricted）
            prog('submit', 0.3, f'提交 NCBI {action}（数分钟，耐心等待）')
            rid, _rtoe = ca.submit_blast(
                action, 'nt' if action == 'blastn' else 'nr',
                f'>{contig}\n{seq}\n', entrez_query='viruses[Organism]')
            prog('poll', 0.5, f'NCBI 任务 {rid} 运行中')
            hits = ca.poll_blast(rid, cancel=cancel)
            result = {'action': action, 'contig': contig, 'rid': rid,
                      'hits': hits}

        with safe_open(out_path, 'wt') as f:
            json.dump(result, f, ensure_ascii=False, indent=1)
        prog('done', 1.0, '完成')
        logger.close()
        return {'out': out_path}

    tid = tm.start(cfg.tr(f'Contig分析·{action} {contig[:30]}',
                         f'Contig {action} {contig[:30]}'), job,
                   log_file=os.path.join(run_dir, 'run.log'))
    return jsonify({'task': tid, 'cached': False})


# ------------------------------------------------------------------
# API: NCBI Entrez 参考序列下载（借鉴 PhyloSuite 会话式翻页设计；
# 域名白名单 eutils.ncbi.nlm.nih.gov + IP 边界校验 + 限速重试）
# ------------------------------------------------------------------
@app.route('/api/ncbi/collections')
def api_ncbi_collections():
    from vp.ncbi_download import list_collections
    return jsonify([c for c in list_collections()
                    if not c.get('name', '').startswith('_')])


# ------------------------------------------------------------------
# API: ICTV 科/属选参（databases/tax/ictv taxa.txt 为数据源；
# 本地优先选 accession → 下载整科/整属 GenBank 集合 + FASTA 参考集）
# ------------------------------------------------------------------
@app.route('/api/ictv/cascade')
def api_ictv_cascade():
    """级联选参下拉数据：?Realm=..&Phylum=..（已选等级作为过滤）→
    {selected, ranks:[{col, zh, options:[{name,n}]}]}。植物病毒口径。"""
    from vp import ictv_db
    if not ictv_db.available():
        abort(400, 'ICTV 库未初始化（python main.py ictv-update）')
    levels = {col: (request.args.get(col) or '').strip()
              for col, _zh in ictv_db.RANK_COLS}
    levels = {k: v for k, v in levels.items() if v}
    try:
        opts = ictv_db.rank_options(levels, plant_only=True)
    except Exception as e:
        abort(400, f'ICTV taxa 读取失败: {e}')
    return jsonify({
        'selected': levels,
        'ranks': [{'col': col, 'zh': zh, 'options': opts.get(col, [])}
                  for col, zh in ictv_db.RANK_COLS]})


def _ictv_levels(body):
    """从请求体取分类级选择 {列名: 值}（非空）+ 最深层（用于日志与定位）。"""
    from vp import ictv_db
    levels = {}
    for col, _zh in ictv_db.RANK_COLS:
        v = str((body.get(col) or '')).strip()
        if v:
            levels[col] = v
    if not levels:
        abort(400, '请至少选择一个分类等级（界/门/纲/目/科/属/种）')
    deepest_col = [c for c, _zh in ictv_db.RANK_COLS if c in levels][-1]
    return levels, deepest_col, levels[deepest_col]


@app.route('/api/ictv/preview', methods=['POST'])
def api_ictv_preview():
    """body: {<分类级列名>: 值, ..., limit, genome} → 选参预览（不下载）。"""
    body = request.get_json(force=True) or {}
    levels, dcol, dval = _ictv_levels(body)
    genome = body.get('genome') or 'complete'
    limit = max(1, min(int(body.get('limit', 200) or 200), 2000))
    per_genus = max(0, min(int(body.get('per_genus', 0) or 0), 50))
    per_species = max(0, min(int(body.get('per_species', 0) or 0), 20))
    from vp import ictv_db
    try:
        rows, total = ictv_db.select_refs(ranks=levels, plant_only=True,
                                          genome=genome,
                                          limit=max(limit * 20, 2000))
    except Exception as e:
        abort(400, f'选参失败: {e}')
    # 抽样上限先于 limit 生效（先每属/每种封顶，再取前 limit 条）
    rows, dropped = ictv_db.cap_per_rank(rows, per_genus, per_species)
    rows = rows[:limit]
    cap_note = []
    if per_genus:
        cap_note.append(f'每属≤{per_genus}')
    if per_species:
        cap_note.append(f'每种≤{per_species}')
    return jsonify({'scope': f'{dcol}={dval}', 'total': total,
                    'n_shown': len(rows),
                    'sampling': '、'.join(cap_note) or '不限',
                    'dropped_by_cap': dropped,
                    'rows': [{'acc': r.get('Genbank', ''),
                              'source': r.get('Source', ''),
                              'species': r.get('Species', ''),
                              'coverage': r.get('Genome_Coverage', ''),
                              'genome': r.get('Genome', ''),
                              'host': r.get('Host_Source', '')}
                             for r in rows]})


@app.route('/api/ictv/download', methods=['POST'])
def api_ictv_download():
    """body: {rank, name, coll, limit, genome} → 后台任务：
    ICTV 选参（本地优先）→ 下载整科/整属 GenBank 集合（含巡检清单）
    → 同步导出 FASTA 参考集（databases/ncbi_refs/<coll>/refs.fa，
    供样品流程建树追加参考）。"""
    body = request.get_json(force=True) or {}
    levels, dcol, dval = _ictv_levels(body)
    coll = (body.get('coll') or '').strip()
    genome = body.get('genome') or 'complete'
    if not coll:
        abort(400, '参数不完整（请填集合名）')
    limit = max(1, min(int(body.get('limit', 200) or 200), 2000))
    per_genus = max(0, min(int(body.get('per_genus', 0) or 0), 50))
    per_species = max(0, min(int(body.get('per_species', 0) or 0), 20))

    def job(log, prog, cancel):
        import shutil
        logger = TaskLogger(callback=log)
        from vp import ictv_db
        from vp.gb_collection import (download_gb_collection, gb_collection_dir,
                                      collection_records)
        from vp.ncbi_download import collection_dir as fasta_coll_dir
        from vp.utils import write_fasta_record

        prog('ictv', 0.05, f'ICTV 选参：{dcol} = {dval}（{len(levels)} 级过滤）')
        rows, total = ictv_db.select_refs(ranks=levels, plant_only=True,
                                          genome=genome,
                                          limit=max(limit * 20, 2000))
        rows, _dropped = ictv_db.cap_per_rank(rows, per_genus, per_species)
        rows = rows[:limit]
        if not rows:
            raise RuntimeError(f'ICTV {dcol} "{dval}" 无匹配参考'
                               '（换一个分类或改 genome=any）')
        accs = [r['Genbank'] for r in rows]
        logger.log(f'ICTV 命中 {total} 条，取前 {len(accs)} 条'
                   f'（{"Complete genome" if genome == "complete" else "全部"}）')

        # 重名集合替换式重建（与建库口径一致，不累积重复记录）
        shutil.rmtree(gb_collection_dir(coll), ignore_errors=True)
        prog('ictv', 0.15, f'下载 GenBank 集合 {coll}（{len(accs)} 条 accession）')
        res = download_gb_collection(coll, accessions=','.join(accs),
                                     max_records=len(accs),
                                     logger=logger, prog=prog, cancel=cancel)

        # 同步 FASTA 参考集（样品流程「NCBI 参考集合」可直接填 coll）
        prog('ictv', 0.85, '同步 FASTA 参考集（样品建树用）')
        fdir = fasta_coll_dir(coll)
        os.makedirs(fdir, exist_ok=True)
        fa = os.path.join(fdir, 'refs.fa')
        n_fa = 0
        with safe_open(fa, 'wt') as f:
            for acc, _org, seq in collection_records(coll):
                if seq:
                    write_fasta_record(f, acc, seq)
                    n_fa += 1
        logger.log(f'FASTA 参考集: {fa}（{n_fa} 条）')
        prog('ictv', 1.0, f'完成：GenBank {res.get("n_records", n_fa)} 条 + '
                    f'FASTA {n_fa} 条')
        logger.close()
        return {'collection': coll, 'scope': f'{dcol}={dval}',
                'levels': levels, 'n_records': n_fa, 'n_fasta': n_fa}

    tid = tm.start(cfg.tr(f'ICTV下载 {coll}', f'ICTV download {coll}'), job,
                   weight='light')
    return jsonify({'task': tid})


@app.route('/api/ncbi/search', methods=['POST'])
def api_ncbi_search():
    """body: {term, db, limit} → {count, rows:[{acc,title,organism,...}]}。"""
    body = request.get_json(force=True)
    term = (body.get('term') or '').strip()
    db = body.get('db') or 'nucleotide'
    if not term or db not in ('nucleotide', 'protein'):
        abort(400, '参数不完整（term 必填，db=nucleotide/protein）')
    from vp.ncbi_download import search_preview
    try:
        return jsonify(search_preview(term, db=db,
                                      limit=int(body.get('limit', 50))))
    except RuntimeError as e:
        abort(502, str(e))


@app.route('/api/ncbi/download', methods=['POST'])
def api_ncbi_download():
    """body: {term, name, db, max} → 后台任务下载集合。"""
    body = request.get_json(force=True)
    term = (body.get('term') or '').strip()
    name = (body.get('name') or '').strip()
    db = body.get('db') or 'nucleotide'
    if not term or not name or db not in ('nucleotide', 'protein'):
        abort(400, '参数不完整（term/name 必填，db=nucleotide/protein）')
    max_records = max(1, min(int(body.get('max', 100) or 100), 10000))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.ncbi_download import download_collection
        res = download_collection(term, name, db=db, max_records=max_records,
                                  logger=logger, prog=prog, cancel=cancel)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'NCBI下载 {name}', f'NCBI download {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


# ------------------------------------------------------------------
# API: GenBank 集合（同属共线性比较输入）与比较任务
# （复用 ncbi_download 的白名单/限速/重试机制，GenBank 平文原样落盘）
# ------------------------------------------------------------------
@app.route('/api/gb/collections')
def api_gb_collections():
    """集合列表（只读概览：读已生成的清单/警告文件，不重新解析 .gb）。

    下划线前缀的集合（_内部测试/_归档）不在列表展示。
    """
    from vp.gb_collection import list_gb_collections, read_manifest
    out = []
    for c in list_gb_collections():
        if c.get('name', '').startswith('_'):
            continue
        try:
            st = read_manifest(c['name'])
            c['records'] = [{'acc': r['acc'], 'organism': r['organism'],
                             'length': r['length'], 'cds': r['cds_count']}
                            for r in st['records']]
            c['warnings'] = st['warnings']
        except (OSError, FileNotFoundError, ValueError):
            c['records'], c['warnings'] = [], []
        c['has_compare'] = os.path.isfile(
            os.path.join(c['dir'], 'compare', 'compare.html'))
        c['has_lovis4u'] = os.path.isfile(
            os.path.join(c['dir'], 'compare', 'lovis4u.pdf'))
        c['has_phylo'] = os.path.isfile(
            os.path.join(c['dir'], 'phylo', 'aln.fasta'))
        try:
            from vp.gb_collection import list_extract_files
            ex = list_extract_files(c['name'])
        except Exception:
            ex = None
        c['has_extract'] = bool(ex)
        c['n_extract_genes'] = (len([x for x in (ex or [])
                                     if x['kind'].endswith('_gene')
                                     and x['kind'].startswith('cds')]) or
                                len({x['path'].split('/')[-1]
                                     for x in (ex or [])
                                     if x['kind'] == 'cds_gene'}))
        out.append(c)
    return jsonify(out)


@app.route('/api/gb/extract', methods=['POST'])
def api_gb_extract():
    """body: {name} → 后台任务：集合 .gb → extract/ 分类分目录提取
    （genome.fa / CDS.fa / PEP.fa / CDS/<基因>.fa / PEP/<基因>.fa / genes.tsv）。"""
    body = request.get_json(force=True)
    name = (body.get('name') or '').strip()
    if not name:
        abort(400, '参数不完整（name 必填）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gb_collection import extract_collection_features
        res = extract_collection_features(name, logger=logger, prog=prog)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'CDS/PEP 提取 {name}', f'Extract CDS/PEP {name}'),
                   job, weight='light')
    return jsonify({'task': tid})


@app.route('/api/gb/extract_files')
def api_gb_extract_files():
    """某集合的提取产物清单（供「参考序列获取」卡展示与送下游）。"""
    name = (request.args.get('name') or '').strip()
    if not name:
        abort(400, '参数不完整（name 必填）')
    from vp.gb_collection import list_extract_files
    items = list_extract_files(name)
    if items is None:
        abort(404, '该集合尚未提取（先点「🧬 提取 CDS/PEP」）')
    return jsonify(items)


@app.route('/api/gb/inspect', methods=['POST'])
def api_gb_inspect():
    """body: {name} → 后台任务全量巡检（重解析全部 .gb，重建清单与警告）。"""
    body = request.get_json(force=True)
    name = (body.get('name') or '').strip()
    if not name:
        abort(400, '参数不完整（name 必填）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gb_collection import inspect_collection
        st = inspect_collection(name, logger=logger)
        logger.close()
        return {'name': name, 'n_records': len(st['records']),
                'n_warnings': len(st['warnings'])}

    tid = tm.start(cfg.tr(f'GenBank巡检 {name}', f'GenBank inspect {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


@app.route('/api/gb/phylo', methods=['POST'])
def api_gb_phylo():
    """body: {name, tree_tool, molecule, gene} → 后台任务：集合比对 + 建树。

    molecule: genome（全基因组，产物 phylo/，结果中心 gb:<name> 展示）/
    cds / pep（基因级：按 gene 关键词从 .gb 提取 CDS/蛋白，
    产物 gene_trees/<gene>_<molecule>/）。
    """
    body = request.get_json(force=True)
    name = (body.get('name') or '').strip()
    tree_tool = body.get('tree_tool') or 'fasttree'
    molecule = (body.get('molecule') or 'genome').strip().lower()
    gene = (body.get('gene') or '').strip()
    if not name:
        abort(400, '参数不完整（name 必填）')
    if tree_tool not in ('fasttree', 'iqtree', 'nj'):
        abort(400, 'tree_tool 仅支持 fasttree / iqtree / nj')
    if molecule not in ('genome', 'cds', 'pep'):
        abort(400, 'molecule 仅支持 genome / cds / pep')
    if molecule in ('cds', 'pep') and not gene:
        abort(400, '基因级建树需提供基因/产物关键词（如 coat protein）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gb_collection import build_collection_phylo
        res = build_collection_phylo(name, tree_tool=tree_tool,
                                     molecule=molecule, gene=gene or None,
                                     logger=logger, prog=prog, cancel=cancel)
        logger.close()
        return res

    tag = name if molecule == 'genome' else f'{name}·{gene}({molecule})'
    tid = tm.start(cfg.tr(f'集合建树 {tag}', f'Collection phylo {tag}'), job)
    return jsonify({'task': tid})


@app.route('/api/gb/download', methods=['POST'])
def api_gb_download():
    """body: {name, term?/accessions?, max} → 后台任务下载 GenBank 集合。"""
    body = request.get_json(force=True)
    name = (body.get('name') or '').strip()
    term = (body.get('term') or '').strip()
    accs = (body.get('accessions') or '').strip()
    if not name or (not term and not accs):
        abort(400, '参数不完整（name 与 term/accessions 必填）')
    if term and accs:
        abort(400, '检索式与 accession 列表二选一')
    max_records = max(1, min(int(body.get('max', 50) or 50), 10000))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gb_collection import (download_gb_collection,
                                      collection_records)
        res = download_gb_collection(name, term=term or None,
                                     accessions=accs or None,
                                     max_records=max_records,
                                     logger=logger, prog=prog, cancel=cancel)
        # 同步 FASTA 参考集（样品流程「NCBI 参考集合」可直接填集合名）
        try:
            from vp.ncbi_download import collection_dir as fasta_coll_dir
            from vp.utils import write_fasta_record
            fdir = fasta_coll_dir(name)
            os.makedirs(fdir, exist_ok=True)
            fa = os.path.join(fdir, 'refs.fa')
            n_fa = 0
            with safe_open(fa, 'wt') as f:
                for acc, _org, seq in collection_records(name):
                    if seq:
                        write_fasta_record(f, acc, seq)
                        n_fa += 1
            logger.log(f'FASTA 参考集: {fa}（{n_fa} 条）')
            res['n_fasta'] = n_fa
        except Exception as e:
            logger.log(f'FASTA 参考集同步失败（不影响 GenBank 集合）: {e}',
                       'WARN')
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'GenBank下载 {name}', f'GenBank download {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


@app.route('/api/gb/import', methods=['POST'])
def api_gb_import():
    """body: {name, files:[本机 .gb 路径]} → 后台任务导入。"""
    body = request.get_json(force=True)
    name = (body.get('name') or '').strip()
    files = body.get('files') or []
    if not name or not files:
        abort(400, '参数不完整（name 与 files 必填）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gb_collection import import_local_gb
        res = import_local_gb(name, files, logger=logger, prog=prog)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'GenBank导入 {name}', f'GenBank import {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


@app.route('/api/compare/preview')
def api_compare_preview():
    """读取集合的 compare.html，提取 body 可视化内容（共线性 SVG + 表）供卡内内嵌。

    参数: name=集合名。返回 {ok, html}（html 为 compare.html 的 <body> 内片段）。
    """
    import re as _re
    name = (request.args.get('name') or '').strip()
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, 'invalid collection name')
    from vp.gb_collection import gb_collection_dir
    try:
        cdir = gb_collection_dir(name)
    except ValueError as e:
        abort(400, str(e))
    p = os.path.join(cdir, 'compare', 'compare.html')
    if not os.path.isfile(p):
        abort(404, '该集合还没有比较结果，先运行比较')
    with safe_open(p) as f:
        html = f.read()
    m = _re.search(r'<body[^>]*>(.*)</body>', html, _re.S | _re.I)
    body = m.group(1) if m else html
    # 去掉可能存在的 <script>（compare.html 通常无，防御）
    body = _re.sub(r'<script[^>]*>.*?</script>', '', body, flags=_re.S | _re.I)
    return jsonify({'ok': True, 'html': body})


@app.route('/api/compare/run', methods=['POST'])
def api_compare_run():
    """body: {collection?, files?, min_ident, min_cov} → 后台同属比较任务。"""
    body = request.get_json(force=True)
    collection = (body.get('collection') or '').strip()
    files = body.get('files') or []
    if not collection and not files:
        abort(400, '参数不完整（collection 或 files 至少一项）')
    try:
        min_ident = max(0.05, min(float(body.get('min_ident', 0.30)), 0.95))
        min_cov = max(0.05, min(float(body.get('min_cov', 0.50)), 1.0))
    except (TypeError, ValueError):
        abort(400, '阈值参数不合法')
    style = body.get('style') or 'lovis'
    if style not in ('lovis', 'category'):
        abort(400, 'style 仅支持 lovis / category')
    # LoVis4u 为共线性主输出（不可用时由内置引擎回退，见 synteny 日志）
    lovis4u = bool(body.get('lovis4u', True))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.synteny import run_comparison
        res = run_comparison(name=collection or None, files=files,
                             min_ident=min_ident, min_cov=min_cov,
                             style=style,
                             lovis4u_pdf=lovis4u,
                             logger=logger, prog=prog, cancel=cancel)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'同属比较 {collection or "自备文件"}',
                         f'Synteny {collection or "adhoc files"}'), job)
    return jsonify({'task': tid})


@app.route('/compare/<name>/')
def page_compare(name):
    """集合的共线性比较结果页（compare.html，含相对引用的 plotly.min.js）。"""
    from vp.gb_collection import gb_collection_dir
    try:
        cdir = gb_collection_dir(name)
    except ValueError as e:
        abort(400, str(e))
    p = os.path.join(cdir, 'compare', 'compare.html')
    if not os.path.isfile(p):
        abort(404, '该集合还没有比较结果，先运行比较')
    return send_file(check_path(p, must_exist=True, in_platform=True))


@app.route('/compare/<name>/<path:filename>')
def page_compare_file(name, filename):
    from vp.gb_collection import gb_collection_dir
    if '..' in filename:
        abort(400, '非法路径')
    try:
        cdir = gb_collection_dir(name)
    except ValueError as e:
        abort(400, str(e))
    p = check_path(os.path.join(cdir, 'compare', filename),
                   must_exist=True, in_platform=True)
    return send_file(p)


# ------------------------------------------------------------------
# API: MSA 网页查看（SNP-only 彩色热图，借鉴 PhyloSuite MSAViewer 设计）
# 样品来源两种：样品流程 results/<sample>/05_phylo/<group>，
# 以及 GenBank 集合建树产物 databases/misc/gb/<name>/phylo
# （伪样品 gb:<name>，与集合名同名的单一分组）。
# ------------------------------------------------------------------
def _msa_group_dir(sample, group):
    """sample+group → 组目录（含 gb: 集合伪样品解析）。"""
    if sample.startswith('gb:'):
        from vp.gb_collection import gb_collection_dir
        try:
            return check_path(os.path.join(gb_collection_dir(sample[3:]),
                                           'phylo'),
                              must_exist=True, in_platform=True)
        except ValueError as e:
            raise ValueError(f'集合名不合法: {e}')
    from vp.pipeline import _safe_sample_name
    from vp.phylo import _safe_name
    sdir = check_path(os.path.join(DIRS['results'],
                                   _safe_sample_name(sample)),
                      must_exist=True, in_platform=True)
    return check_path(os.path.join(sdir, '05_phylo', _safe_name(group)),
                      must_exist=True, in_platform=True)


def _list_phylo_groups(pdir):
    """目录 → [{group, aln, trees}]（含比对文件的分组，优先 aln.trim）。"""
    groups = []
    for g in sorted(os.listdir(pdir)):
        gdir = os.path.join(pdir, g)
        if not os.path.isdir(gdir):
            continue
        aln = os.path.join(gdir, 'aln.trim.fasta')
        kind = 'trim'
        if not os.path.isfile(aln):
            aln = os.path.join(gdir, 'aln.fasta')
            kind = 'full'
        if os.path.isfile(aln):
            groups.append({'group': g, 'aln': kind, 'trees': _tree_files(gdir)})
    return groups


@app.route('/api/msa/samples')
def api_msa_samples():
    """列出含比对结果的样品与集合建树产物（伪样品 gb:<集合名>）。"""
    out = []
    res = check_path(DIRS['results'], must_exist=False, in_platform=True)
    if os.path.isdir(res):
        for name in sorted(os.listdir(res)):
            if name.startswith('_'):      # 归档/内部样品不进下拉
                continue
            pdir = os.path.join(res, name, '05_phylo')
            if not os.path.isdir(pdir):
                continue
            groups = _list_phylo_groups(pdir)
            if groups:
                out.append({'sample': name, 'groups': groups})
    # GenBank 集合建树产物（比较基因组：下载/导入集合 → 直接建树）
    from vp.gb_collection import list_gb_collections
    try:
        cols = list_gb_collections()
    except OSError:
        cols = []
    for c in cols:
        pdir = os.path.join(c['dir'], 'phylo')
        if not os.path.isdir(pdir):
            continue
        aln = os.path.join(pdir, 'aln.trim.fasta')
        if not os.path.isfile(aln):
            aln = os.path.join(pdir, 'aln.fasta')
        if os.path.isfile(aln):
            out.append({'sample': f'gb:{c["name"]}',
                        'label': f'🧬 {c["name"]}（集合）',
                        'groups': [{'group': c['name'], 'aln': 'full',
                                    'trees': _tree_files(pdir)}]})
    return jsonify(out)


_TREE_FILES = (('tree.nwk', 'FastTree'),
               ('nj.nwk', 'NJ 快速树（identity 距离）'),
               ('iqtree.treefile', 'IQ-TREE ML 树'),
               ('iqtree.contree', 'IQ-TREE 一致树'))


def _tree_files(gdir):
    """组目录下实际存在的树文件列表（供树查看器下拉）。"""
    return [{'file': fn, 'label': label}
            for fn, label in _TREE_FILES
            if os.path.isfile(os.path.join(gdir, fn))]


@app.route('/api/tree/data')
def api_tree_data():
    """sample+group(+file) → 进化树 Newick 与建树摘要（前端 Archaeopteryx.js 渲染）。"""
    sample = request.args.get('sample') or ''
    group = request.args.get('group') or ''
    fname = request.args.get('file') or ''
    try:
        gdir = _msa_group_dir(sample, group)
    except ValueError as e:
        abort(400, f'参数不合法: {e}')
    choices = ([(fname, '')] if fname else list(_TREE_FILES))
    tree_file = None
    for fn, _label in choices:
        p = os.path.join(gdir, fn)
        if fn and os.path.isfile(p):
            tree_file = p
            break
    if not tree_file:
        abort(404, '该组没有树文件（tree.nwk / nj.nwk / iqtree.treefile / iqtree.contree）')
    with safe_open(tree_file) as f:
        nwk = f.read().strip()
    if not nwk:
        abort(404, '树文件为空')
    # 从 summary.json 取建树工具/模型等摘要（无则忽略）
    meta = {}
    summary = os.path.join(gdir, 'summary.json')
    if os.path.isfile(summary):
        try:
            with safe_open(summary) as f:
                sj = json.load(f)
            g = (next((x for x in (sj.get('groups') or [])
                       if x.get('group') == group), None)
                 if isinstance(sj.get('groups'), list) else None) or {}
            base = os.path.basename(tree_file)
            if base == 'nj.nwk':
                tool = 'NJ'
            elif (base.endswith('.treefile') or base.endswith('.contree')
                    or g.get('tree_tool') == 'iqtree'
                    or sj.get('tree_tool') == 'iqtree'):
                tool = 'IQ-TREE'
            else:
                tool = 'FastTree'
            meta = {'tool': tool,
                    'model': g.get('tree_model') or sj.get('model'),
                    'logl': g.get('tree_logl') or sj.get('logl')}
        except (OSError, ValueError):
            pass
    return jsonify({'newick': nwk, 'file': os.path.basename(tree_file), **meta})


@app.route('/api/msa/data')
def api_msa_data():
    """sample+group → SNP-only 展示数据（变异列矩阵/共识/多样性）。"""
    sample = request.args.get('sample') or ''
    group = request.args.get('group') or ''
    try:
        gdir = _msa_group_dir(sample, group)
    except ValueError as e:
        abort(400, f'参数不合法: {e}')
    aln = os.path.join(gdir, 'aln.trim.fasta')
    if not os.path.isfile(aln):
        aln = os.path.join(gdir, 'aln.fasta')
    if not os.path.isfile(aln):
        abort(404, '该组无比对文件')
    from vp.msa_view import snp_view_data
    try:
        return jsonify(snp_view_data(aln))
    except ValueError as e:
        abort(400, str(e))


@app.route('/api/sdt/data')
def api_sdt_data():
    """sample+group → SDT 口径全长成对 identity 矩阵（结果中心在线热图）。

    优先读样品流程产物 sdt_matrix.csv；无则从比对 fasta 现算
    （与 vp/phylo.pairwise_identity_matrix 同口径）。
    """
    sample = request.args.get('sample') or ''
    group = request.args.get('group') or ''
    try:
        gdir = _msa_group_dir(sample, group)
    except ValueError as e:
        abort(400, f'参数不合法: {e}')
    csv = os.path.join(gdir, 'sdt_matrix.csv')
    if os.path.isfile(csv):
        names, mat = [], []
        with safe_open(csv) as f:
            rows = [line.rstrip('\n').split(',') for line in f if line.strip()]
        if len(rows) >= 2:
            names = rows[0][1:]
            for r in rows[1:]:
                mat.append([float(x) if x else 0.0 for x in r[1:]])
            return jsonify({'names': names, 'matrix': mat,
                            'source': 'sdt_matrix.csv'})
    aln = os.path.join(gdir, 'aln.trim.fasta')
    if not os.path.isfile(aln):
        aln = os.path.join(gdir, 'aln.fasta')
    if not os.path.isfile(aln):
        abort(404, '该组没有 SDT 矩阵，也没有比对文件')
    from vp.phylo import pairwise_identity_matrix
    names, mat = pairwise_identity_matrix(aln)
    return jsonify({'names': names, 'matrix': mat, 'source': 'aln'})


# ------------------------------------------------------------------
# LOGAN 溯源模块（独立模块页：病毒序列 → Logan-Search → 物种/样本溯源）
# 半自动流程：平台只负责切片生成查询/解析导入的结果表并出报告；
# 提交查询由用户在自己的浏览器打开 logan-search.org 完成，
# 本模块不发任何外部网络请求。
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# 拖拽上传 + 轻量序列查看器
# ------------------------------------------------------------------
@app.route('/api/paste_input', methods=['POST'])
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


@app.route('/api/upload', methods=['POST'])
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


@app.route('/api/seqview')
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


# ------------------------------------------------------------------
# 项目批处理队列（顺序执行：同一时间只跑一个样品，防内存/磁盘打满）
# ------------------------------------------------------------------
class SampleQueue:
    """样品批处理队列：FIFO 顺序跑，每项复用与单样品运行完全相同的
    run_analysis 代码路径；状态持久化到 tasks/queue.json。"""

    QUEUE_FILE = os.path.join(DIRS['tasks'], 'queue.json')

    def __init__(self):
        self.lock = threading.RLock()
        self.items = []
        self.thread = None
        self._load()

    def _load(self):
        try:
            with safe_open(self.QUEUE_FILE) as f:
                self.items = json.load(f)
            for it in self.items:
                if it.get('status') in ('queued', 'running'):
                    it['status'] = 'queued'
                    it['task'] = None
        except (OSError, ValueError):
            self.items = []

    def _persist(self):
        with self.lock:
            try:
                os.makedirs(os.path.dirname(self.QUEUE_FILE), exist_ok=True)
                # 原子写：断电/崩溃不损坏现有队列
                tmp = self.QUEUE_FILE + '.tmp'
                with safe_open(tmp, 'wt') as f:
                    json.dump(self.items, f, ensure_ascii=False, indent=1)
                os.replace(tmp, self.QUEUE_FILE)
            except OSError:
                pass

    def add(self, samples, stages, params, project=None):
        """samples: [样品名]（须已创建并含 input.json）。"""
        added = []
        with self.lock:
            existing = {it['sample'] for it in self.items
                        if it['status'] in ('queued', 'running')}
            for s in samples:
                if s in existing:
                    continue
                it = {'id': uuid.uuid4().hex[:10], 'sample': s,
                      'stages': stages, 'params': params,
                      'project': project or '', 'status': 'queued',
                      'task': None, 'error': '', 'added': time.time()}
                self.items.append(it)
                added.append(it)
        self._persist()
        self._ensure_worker()
        return added

    def _ensure_worker(self):
        with self.lock:
            if self.thread and self.thread.is_alive():
                return
            self.thread = threading.Thread(target=self._worker, daemon=True,
                                           name='sample-queue')
            self.thread.start()

    def _worker(self):
        while True:
            with self.lock:
                nxt = next((it for it in self.items
                            if it['status'] == 'queued'), None)
            if not nxt:
                return
            nxt['status'] = 'running'
            self._persist()
            try:
                self._run_entry(nxt)
            except Exception as e:
                nxt['status'] = 'failed'
                nxt['error'] = str(e)
            self._persist()

    def _run_entry(self, it):
        from vp.pipeline import (load_sample_input, run_analysis,
                                 _safe_sample_name, pipeline_overview,
                                 STAGE_ORDER)
        s = _safe_sample_name(it['sample'])
        sd = check_path(os.path.join(DIRS['results'], s),
                        must_exist=True, in_platform=True)
        r1, r2, _proj = load_sample_input(sd)
        if not r1:
            raise RuntimeError(f'样品 {s} 缺少输入记录（input.json）')
        check_path(r1, must_exist=True)
        stages = it.get('stages') or None
        if stages:
            stages = [x for x in stages if x in STAGE_ORDER]
        if not stages:
            # 与「依次运行剩余步骤」一致：所有未完成且可用的阶段
            try:
                ov = pipeline_overview(sd)['stages']
                done = {x['stage'] for x in ov if x['status'] == 'done'}
                stages = [x['stage'] for x in ov
                          if x['stage'] in STAGE_ORDER
                          and x['stage'] not in done
                          and x['status'] != 'unavailable']
            except Exception:
                stages = []
        if not stages:
            raise RuntimeError(f'样品 {s} 没有需要运行的阶段（全部已完成）')

        def job(log, prog, cancel):
            logger = TaskLogger(callback=log)
            run_analysis(sample=s, r1=r1, r2=r2, stages=stages,
                         **_analysis_kwargs(it.get('params') or {}),
                         logger=logger, progress=prog)
            logger.close()
            return sd

        q_log_dir = check_path(os.path.join(sd, 'logs'),
                               must_exist=False, in_platform=True)
        os.makedirs(q_log_dir, exist_ok=True)
        tid = tm.start(cfg.tr(f'队列·{s}', f'Queue·{s}'), job,
                       log_file=os.path.join(
                           q_log_dir,
                           f'queue_{time.strftime("%Y%m%d_%H%M%S")}.log'))
        it['task'] = tid
        self._persist()
        while True:
            time.sleep(2.0)
            snap = tm.snapshot(tid, log_lines=0, include_result=False)
            if snap is None or snap['status'] != 'running':
                it['status'] = ('done' if snap and snap['status'] == 'done'
                                else 'failed' if snap else 'failed')
                it['error'] = (snap or {}).get('error') or ''
                return

    def snapshot(self):
        with self.lock:
            items = [dict(it) for it in self.items]
        running = any(it['status'] == 'running' for it in items)
        return {'running': running, 'items': items}

    def remove(self, entry_id):
        with self.lock:
            it = next((x for x in self.items if x['id'] == entry_id), None)
            if not it or it['status'] == 'running':
                return False
            self.items.remove(it)
        self._persist()
        return True

    def clear_finished(self):
        with self.lock:
            self.items = [it for it in self.items
                          if it['status'] in ('queued', 'running')]
        self._persist()
        return True


sample_queue = SampleQueue()
# 注册到 vp/web/state.py，供 download blueprint 的「转入分析流程」取用
# （注册模式而非 import app，避免循环导入）
from vp.web import state as _web_state  # noqa: E402

_web_state.set_sample_queue(sample_queue)


@app.route('/api/queue')
def api_queue():
    return jsonify(sample_queue.snapshot())


@app.route('/api/queue/add', methods=['POST'])
def api_queue_add():
    body = request.get_json(force=True) or {}
    samples = body.get('samples') or []
    if not samples:
        abort(400, '请提供样品列表')
    from vp.pipeline import _safe_sample_name
    valid = []
    for s in samples:
        sd = check_path(os.path.join(DIRS['results'], _safe_sample_name(s)),
                        must_exist=False, in_platform=True)
        if not os.path.isfile(os.path.join(sd, '00_prep', 'input.json')):
            abort(400, f'样品 {s} 不存在或缺少输入记录（请先创建样品）')
        valid.append(_safe_sample_name(s))
    stages = body.get('stages') or None       # None = 依次运行剩余步骤
    if stages:
        from vp.pipeline import STAGE_ORDER
        stages = [s for s in stages if s in STAGE_ORDER]
        if not stages:
            abort(400, 'stages 参数不合法')
    added = sample_queue.add(valid, stages,
                             body.get('params') or {},
                             project=str(body.get('project') or '') or None)
    return jsonify({'added': len(added)})


@app.route('/api/queue/<entry_id>/remove', methods=['POST'])
def api_queue_remove(entry_id):
    if not sample_queue.remove(entry_id):
        abort(400, '条目不存在或正在运行（请先取消任务）')
    return jsonify({'ok': True})


@app.route('/api/queue/clear_finished', methods=['POST'])
def api_queue_clear():
    sample_queue.clear_finished()
    return jsonify({'ok': True})


# ------------------------------------------------------------------
# 公共数据下载中心（SRA/ENA + GSA/NGDC → aria2c 批量下载 → 进分析流程）
# ------------------------------------------------------------------


# 兼容别名（路由函数体内使用）


# ------------------------------------------------------------------
# NCBI 提交准备（vp.ncbi_submit，源自 MMPV-RNA submission_gui 功能面）
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# 公共数据检索（vp.public_meta，源自 MMPV-RNA public_metadata_pipeline）
# 单工作台模式：/meta 只操作一个当前项目（Search / Info 双存储，对应
# MMPV-RNA metadata_gui 的 8 个标签），每次在线检索自动新建
# <物种>_<时间戳> 项目，历史项目自动归档到结果中心查看。
# ------------------------------------------------------------------
# 三张标准表：search=检索结果（S 存储）core14=统一元数据（I 存储）
# full=全字段明细（只读，下载用）

# AI 元数据清洗提示词（与 MMPV-RNA metadata_gui/controllers/ai_completer.py
# 同源：只做去污染/归位/格式化，禁止凭空捏造）


# ------------------------------------------------------------------
# LOGAN 批量自动提交（logan_submit.py，Selenium 驱动本机浏览器）
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# API: 状态
# ------------------------------------------------------------------
@app.route('/api/tools')
def api_tools():
    return jsonify({'threads': cfg.threads, 'tools': cfg.tool_status()})


@app.route('/api/dbs')
def api_dbs():
    from vp.kunpeng import db_ready
    out = {}
    for key in ('host', 'virus'):
        d = cfg.databases[key]
        out[key] = {'path': d,
                    'ready': db_ready(d) if os.path.isdir(d) else False}
    # 通用参考库（refvirus / rvdb / 其它 kunpeng 库目录自动探测）
    for key, mapping in (('refvirus', ('virus', 'ref')),
                         ('rvdb', ('virus', 'rvdb')),
                         ('k2viral', None)):
        if mapping:
            d = db_path(*mapping)
        else:
            d = os.path.join(DIRS['databases'], 'k2viral_db')
        out[key] = {'path': d,
                    'ready': db_ready(d) if os.path.isdir(d) else False}
    from vp.taxonomy import taxonomy_ready
    out['taxonomy'] = {'ready': taxonomy_ready(), 'path': DIRS['taxonomy']}
    # 病毒分类库统一为「自备预构建」：扫描 databases/ 下的 kunpeng 库目录
    # （宿主库除外）。使用者把库目录放进来即自动识别，工具②/④的病毒库
    # 下拉与构建页列表都以此为准。准入只看 db_ready（结构完整可用），
    # 不看 .building 标记——平台意外退出会留下孤儿标记，不应隐藏可用库。
    libs = []
    base = DIRS['databases']
    for dp, dn, fn in os.walk(base):
        if not any(f.startswith('hash_') and f.endswith('.k2d') for f in fn):
            continue
        rel = os.path.relpath(dp, base).replace('\\', '/')
        if rel in ('.', ''):
            continue
        # 宿主分类库（host/classify）不列入病毒库下拉
        if rel.startswith('host/'):
            continue
        try:
            if not db_ready(dp):
                continue
        except (ValueError, OSError):
            continue
        libs.append({'name': rel, 'path': f'databases/{rel}'})
    # 外部登记的自备库（platform.json extra_virus_libs，可为平台外任意位置）
    seen = {os.path.normpath(os.path.join(DIRS['databases'], l['name']))
            for l in libs}
    for p in cfg.extra_virus_libs:
        n = os.path.normpath(os.path.abspath(p))
        if n in seen or not _db_dir_ready(n):
            continue
        libs.append({'name': os.path.basename(n), 'path': p,
                     'external': True})
    out['virus_libs'] = libs
    return jsonify(out)


def _db_dir_ready(d):
    """kunpeng 库结构完整性（不限定平台目录，用于加载/登记自备库）。"""
    if not os.path.isdir(d):
        return False
    if not all(os.path.isfile(os.path.join(d, f))
               for f in ('opts.k2d', 'taxo.k2d', 'hash_config.k2d')):
        return False
    import glob as _glob
    return bool(_glob.glob(os.path.join(d, 'hash_*.k2d')))


def _abs_arg(raw):
    d = raw if os.path.isabs(raw) else os.path.join(PLATFORM_ROOT, raw)
    return os.path.normpath(os.path.abspath(d))


def _within_allowed_roots(d):
    from vp.config import write_roots
    d = os.path.normpath(os.path.abspath(d))
    return any(d == os.path.normpath(r)
               or d.startswith(os.path.normpath(r) + os.sep)
               for r in write_roots())


@app.route('/api/load_taxonomy', methods=['POST'])
def api_load_taxonomy():
    """加载已有 NCBI Taxonomy（含 nodes.dmp/names.dmp 的目录）→
    复制 .dmp/.pkl 进平台标准位置 databases/tax/core。"""
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少 taxonomy 目录')
    d = _abs_arg(raw)
    if not os.path.isdir(d):
        abort(400, f'目录不存在: {raw}')
    missing = [f for f in ('nodes.dmp', 'names.dmp')
               if not os.path.isfile(os.path.join(d, f))]
    if missing:
        abort(400, f'缺少 {"、".join(missing)}（不是有效的 taxonomy 目录）')
    os.makedirs(DIRS['taxonomy'], exist_ok=True)
    import shutil as _shutil
    copied = []
    for fn in sorted(os.listdir(d)):
        src = os.path.join(d, fn)
        dst = os.path.join(DIRS['taxonomy'], fn)
        if not os.path.isfile(src) or not fn.lower().endswith(('.dmp', '.pkl')):
            continue
        if os.path.normpath(os.path.abspath(src)) == os.path.normpath(dst):
            copied.append(fn)                     # 已在标准位置，跳过
            continue
        _shutil.copyfile(src, dst)
        copied.append(fn)
    from vp.taxonomy import taxonomy_ready
    return jsonify({'ok': taxonomy_ready(), 'copied': copied,
                    'dir': DIRS['taxonomy']})


@app.route('/api/load_host_db', methods=['POST'])
def api_load_host_db():
    """加载已构建宿主库（kunpeng）登记为平台宿主库（持久化 platform.json）。

    宿主库参与建库写流程，目录须位于平台目录内（或设置页配置的数据库根）。"""
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少宿主库目录')
    d = _abs_arg(raw)
    if not _db_dir_ready(d):
        abort(400, '不是有效的 kunpeng 库目录（缺 hash_*.k2d / opts.k2d / '
                   'taxo.k2d / hash_config.k2d）: ' + raw)
    if not _within_allowed_roots(d):
        abort(400, '宿主库目录须位于平台目录内（或先在「设置 → 数据库目录」'
                   '配置自定义数据库根）')
    cfg.databases['host'] = d
    cfg.save()
    return jsonify({'ok': True, 'path': d})


@app.route('/api/register_virus_lib', methods=['POST'])
def api_register_virus_lib():
    """登记 / 移除外部已构建病毒库目录（持久化；分类只读，允许平台外任意位置）。"""
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少病毒库目录')
    d = _abs_arg(raw)
    if body.get('remove'):
        cfg.extra_virus_libs = [p for p in cfg.extra_virus_libs
                                if os.path.normpath(os.path.abspath(p)) != d]
        cfg.save()
        return jsonify({'ok': True, 'libs': cfg.extra_virus_libs})
    if not _db_dir_ready(d):
        abort(400, '不是有效的 kunpeng 库目录（缺 hash_*.k2d / opts.k2d / '
                   'taxo.k2d / hash_config.k2d）: ' + raw)
    if d not in cfg.extra_virus_libs:
        cfg.extra_virus_libs.append(d)
        cfg.save()
    return jsonify({'ok': True, 'libs': cfg.extra_virus_libs})


@app.route('/api/browse')
def api_browse():
    """目录浏览（只读）。平台内目录可用相对路径；也支持浏览整台电脑
    （绝对路径，含盘符列表）。选择的数据输入文件可为任意位置；
    分析输出仍严格限制在平台目录内。"""
    raw = (request.args.get('path', '') or '.').strip()
    # all=1（选目录模式）：文件不限扩展名，便于查看目录里有什么
    show_all = request.args.get('all') == '1'
    if raw in ('.', '', '/', '~', '此电脑'):
        drives = [f'{c}:\\' for c in 'CDEFGHIJKLMNOPQRSTUVWXYZ'
                  if os.path.exists(f'{c}:\\')]
        return jsonify({'cwd': '此电脑', 'parent': '', 'dirs': drives,
                        'files': []})
    p = os.path.normpath(os.path.abspath(raw))
    if not os.path.isdir(p):
        abort(400, '目录不存在')
    dirs, files = [], []
    try:
        for name in sorted(os.listdir(p)):
            child = os.path.join(p, name)
            try:
                if os.path.isdir(child):
                    dirs.append(name)
                elif os.path.isfile(child) and (
                        show_all or name.lower().endswith(
                        ('.fastq.gz', '.fq.gz', '.fastq', '.fq',
                         '.fasta', '.fa', '.fna', '.fas',
                         '.fa.gz', '.fasta.gz', '.tsv', '.sra',
                         '.dmp', '.pkl', '.map',
                         '.gb', '.gbk', '.gbff', '.genbank'))):
                    files.append({'name': name,
                                  'size': fmt_size(os.path.getsize(child))})
            except OSError:
                continue                                  # 跳过无权限项
    except PermissionError:
        abort(400, '无权限访问该目录')
    except OSError:
        abort(400, '无法访问该路径')
    parent = os.path.dirname(p)
    return jsonify({'cwd': p,
                    'parent': '' if parent == p else parent,
                    'dirs': dirs, 'files': files})


# ------------------------------------------------------------------
# API: 任务
# ------------------------------------------------------------------
@app.route('/api/tasks')
def api_tasks():
    # logs=0：只要元信息不要日志（任务中心徒轮询用，日志走单独端点）
    logs = request.args.get('logs') not in ('0', 'false', 'no')
    if request.args.get('full'):
        arch, arch_total = tm._read_archived(logs=logs)
        return jsonify({'active': tm.list_all(logs=logs),
                        'archived': arch,
                        'archived_total': arch_total})
    return jsonify(tm.list_all(logs=logs))


@app.route('/api/task/<tid>')
def api_task(tid):
    try:
        lines = min(int(request.args.get('log_lines') or 80), 1000)
    except (TypeError, ValueError):
        lines = 80
    snap = tm.snapshot(tid, log_lines=lines)
    if not snap:
        abort(404)
    return jsonify(snap)


@app.route('/api/task/<tid>/log')
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


@app.route('/api/task/<tid>/stream')
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


@app.route('/api/task/<tid>/cancel', methods=['POST'])
def api_task_cancel(tid):
    return jsonify({'ok': tm.cancel(tid)})


@app.route('/api/task/<tid>/delete', methods=['POST'])
def api_task_delete(tid):
    ok, msg = tm.delete(tid)
    if not ok:
        abort(400, msg)
    return jsonify({'ok': True})


@app.route('/api/tasks/clear_finished', methods=['POST'])
def api_tasks_clear_finished():
    return jsonify({'removed': tm.clear_finished()})


# ------------------------------------------------------------------
# API: 建库
# ------------------------------------------------------------------
def _job_build_taxonomy(log, prog, cancel):
    from vp.taxonomy import prepare_taxonomy
    logger = TaskLogger(callback=log)
    prog('taxonomy', 0.1, '准备 NCBI taxonomy')
    prepare_taxonomy(logger=logger)
    prog('taxonomy', 1.0, 'taxonomy 就绪')
    logger.close()
    return 'ok'


def _job_build_host_db(body):
    def job(log, prog, cancel):
        from vp.kunpeng import build_host_db
        logger = TaskLogger(callback=log)
        # 建库内存峰值与线程数成正比（每线程约 1GB 缓冲），默认限 8
        threads = int(body.get('threads') or min(8, cfg.threads))
        prog('build_host', 0.05, f"注入 taxid={body['taxid']} 并建库 (线程 {threads})")
        build_host_db(body['genome'], int(body['taxid']),
                      hash_capacity=body.get('hash_capacity', '256M'),
                      threads=threads,
                      logger=logger, rebuild=bool(body.get('rebuild')),
                      clean_mid=bool(body.get('clean_mid')))
        prog('build_host', 1.0, '宿主库就绪')
        logger.close()
        return 'ok'
    return job


def _job_build_virus_db(body):
    def job(log, prog, cancel):
        from vp.kunpeng import build_virus_db
        logger = TaskLogger(callback=log)
        prog('build_virus', 0.05, '解析 info 表并建库')
        build_virus_db(body['fasta'], body['info'],
                       hash_capacity=body.get('hash_capacity', '64M'),
                       threads=int(body.get('threads', cfg.threads)),
                       logger=logger, rebuild=bool(body.get('rebuild')),
                       clean_mid=bool(body.get('clean_mid')))
        prog('build_virus', 1.0, '病毒库就绪')
        logger.close()
        return 'ok'
    return job


@app.route('/api/build_taxonomy', methods=['POST'])
def api_build_taxonomy():
    tid = tm.start(cfg.tr('下载/准备 Taxonomy', 'Download/prepare Taxonomy'),
                   _job_build_taxonomy, weight='light')
    return jsonify({'task': tid})


@app.route('/api/build_host_db', methods=['POST'])
def api_build_host_db():
    body = request.get_json(force=True)
    for k in ('genome', 'taxid'):
        if not body.get(k):
            abort(400, f'缺少参数 {k}')
    check_path(body['genome'], must_exist=True)
    tid = tm.start(cfg.tr(f"宿主库构建 (taxid={body['taxid']})",
                         f"Host DB build (taxid={body['taxid']})"),
                   _job_build_host_db(body))
    return jsonify({'task': tid})


def _universal_db_job(source):
    """通用病毒库建库 job 构造器。"""
    def job(log, prog, cancel):
        from vp.universal_ref import build_universal_db
        logger = TaskLogger(callback=log)
        prog('prep', 0.05, '解析 accession→taxid 映射')
        prog('build', 0.3, 'kunpeng add-library + build-db（见日志）')
        res = build_universal_db(
            source, hash_capacity='2G', logger=logger, rebuild=True)
        prog('build', 0.9, '建库完成，校验库文件')
        prog('done', 1.0, '完成')
        logger.close()
        return res
    return job


@app.route('/api/convert_kraken2', methods=['POST'])
def api_convert_kraken2():
    """Kraken2 库包/目录 → kunpeng 分片库（kunpeng hashshard，方式 C）。

    body: {tar: "databases/k2_viral_20260626.tar.gz"（平台内路径）,
           name: "k2viral"（目标库目录名 databases/<name>）,
           hash_capacity: "1G"}
    """
    body = request.get_json(force=True) or {}
    tar = (body.get('tar') or '').strip()
    name = (body.get('name') or 'k2viral').strip()
    if not tar:
        abort(400, '缺少 Kraken2 库包路径')
    check_path(tar, must_exist=True)
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '库名仅限字母数字-_')
    db_dir = os.path.join(DIRS['databases'], name)

    def job(log, prog, cancel):
        from vp.kunpeng import convert_kraken2
        logger = TaskLogger(callback=log)
        prog('convert', 0.2, '解包 + hashshard 转换（见日志）')
        res = convert_kraken2(tar, db_dir,
                              hash_capacity=body.get('hash_capacity') or '1G',
                              logger=logger)
        prog('done', 1.0, '完成')
        logger.close()
        return {'db_dir': res}

    tid = tm.start(cfg.tr(f'Kraken2 库转换 {name}', f'Kraken2 convert {name}'),
                   job)
    return jsonify({'task': tid})


@app.route('/api/build_kv_index', methods=['POST'])
def api_build_kv_index():
    """构建「病毒鉴定库」：参考 FASTA → minibwa / salmon 比对索引。

    body: {fasta: 参考 FASTA（平台内外均可，只读）,
           engines: ['minibwa','salmon']（默认两者）,
           name: 库名（默认 kv_index，限字母数字-_）,
           threads: 线程数,
           ref_info: 可选，参考注释 TSV，与 FASTA 同目录归档}

    产物：virus-db/<name>/ 下 minibwa.{mbw,l2b}（或 <name>/minibwa/）与
    salmon_k31/，并写一份 manifest.json 记录来源与时间，便于 t-consensus /
    kvsuite 索引复用。与「病毒分类库」（kunpeng）无关，两者互不影响。
    """
    body = request.get_json(force=True) or {}
    fasta = (body.get('fasta') or '').strip()
    if not fasta:
        abort(400, '缺少参数 fasta')
    fasta = check_path(fasta if os.path.isabs(fasta)
                       else os.path.join(PLATFORM_ROOT, fasta),
                       must_exist=True)
    if not os.path.isfile(fasta):
        abort(400, f'参考 FASTA 不是文件：{fasta}')

    engines = body.get('engines') or ['minibwa', 'salmon']
    if isinstance(engines, str):
        engines = [e.strip() for e in engines.split(',') if e.strip()]
    engines = [e for e in engines if e in ('minibwa', 'salmon')]
    if not engines:
        abort(400, 'engines 需包含 minibwa 或 salmon')

    name = (body.get('name') or 'kv_index').strip()
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '库名仅限字母数字-_')
    out_dir = check_path(os.path.join(DIRS['virus_src'], name),
                         must_exist=False, in_platform=True)

    ref_info = (body.get('ref_info') or '').strip()
    if ref_info:
        ref_info = check_path(ref_info if os.path.isabs(ref_info)
                              else os.path.join(PLATFORM_ROOT, ref_info),
                              must_exist=True)

    try:
        threads = max(1, min(int(body.get('threads') or cfg.threads),
                             (os.cpu_count() or 4) * 4))
    except (TypeError, ValueError):
        threads = cfg.threads

    def job(log, prog, cancel):
        from known_virus_suite.kv_common import ToolRegistry, setup_logger
        from known_virus_suite.kv_engines import make_engine
        os.makedirs(out_dir, exist_ok=True)
        # kv_engines 的 logger 需是标准 logging.Logger（.info/.warning/.error）。
        # tag 带时间戳避免多次建库复用同一 logger 时 handler 累积。
        logger, _logs = setup_logger(
            out_dir, tag='kv_index_build_' + time.strftime('%H%M%S'))

        class _PanelHandler(logging.Handler):
            """把 kv_engines 的日志推到任务面板"""
            def emit(self, rec):
                try:
                    log(f'{time.strftime("%H:%M:%S", time.localtime(rec.created))} '
                        f'{rec.levelname} - {rec.getMessage()}')
                except Exception:                # noqa: BLE001
                    pass

        logger.addHandler(_PanelHandler())
        reg = ToolRegistry(logger)
        reg.probe()
        logger.info(f"参考 FASTA: {fasta}")
        logger.info(f"引擎: {', '.join(engines)}  线程: {threads}")
        logger.info(f"产物目录: {out_dir}")

        made, failed = [], []
        for i, eng_name in enumerate(engines):
            if cancel is not None and cancel.is_set():
                raise RuntimeError('用户取消')
            prog(eng_name, 0.05 + 0.9 * i / len(engines), f'构建 {eng_name} 索引')
            if not reg.has(eng_name):
                failed.append(f'{eng_name}: 未找到可执行文件')
                logger.error(f'{eng_name} 不可用，跳过')
                continue
            try:
                engine = make_engine(eng_name, reg, threads=threads,
                                     logger=logger)
                # 各引擎的 build_index 落在 out_dir/<引擎子目录>
                idx = engine.build_index(fasta, out_dir, threads)
                made.append({'engine': eng_name, 'index': str(idx)})
                logger.info(f'{eng_name} 索引就绪 -> {idx}')
            except Exception as exc:                # noqa: BLE001
                failed.append(f'{eng_name}: {exc}')
                logger.error(f'{eng_name} 建索引失败：{exc}')

        if not made:
            raise RuntimeError('所有引擎建索引均失败：' + '; '.join(failed))

        manifest = {
            'name': name,
            'reference': fasta,
            'ref_info': ref_info or '',
            'engines': made,
            'failed': failed,
            'threads': threads,
            'built_at': time.strftime('%Y-%m-%d %H:%M:%S'),
        }
        try:
            with open(os.path.join(out_dir, 'manifest.json'), 'w',
                      encoding='utf-8') as fh:
                json.dump(manifest, fh, ensure_ascii=False, indent=2)
        except OSError as exc:
            logger.warning(f'manifest 写入失败：{exc}')
        # 关闭 FileHandler，否则日志文件句柄会泄漏并锁住产物目录
        for h in list(logger.handlers):
            try:
                h.close()
            except Exception:                    # noqa: BLE001
                pass
            logger.removeHandler(h)
        prog('done', 1.0, '建库完成')
        return {'out': out_dir, 'engines': [m['engine'] for m in made],
                'failed': failed}

    tid = tm.start(cfg.tr(f'病毒鉴定库构建 {name}', f'Virus index build {name}'), job)
    return jsonify({'task': tid, 'out': out_dir})


def _kv_index_probe(d, name=None):
    """探测单个目录是否是一个可用的病毒鉴定库（minibwa / salmon 索引）。

    返回 lib 字典（含 minibwa/salmon 布尔与各自路径），不可用时返回 None。
    仅探测结构，不判断目录所在位置，平台内与平台外目录共用。
    """
    lib = {'name': name or os.path.basename(d), 'path': d, 'minibwa': False,
           'salmon': False, 'manifest': None}
    for cand in (os.path.join(d, 'minibwa.mbw'),
                 os.path.join(d, 'minibwa', 'minibwa.mbw'),
                 os.path.join(d, 'harmonized_minibwa.mbw')):
        if os.path.isfile(cand):
            lib['minibwa'] = True
            lib['minibwa_path'] = cand
            break
    for cand in (os.path.join(d, 'salmon_k31', 'info.json'),
                 os.path.join(d, 'salmon', 'salmon_k31', 'info.json')):
        if os.path.isfile(cand):
            lib['salmon'] = True
            lib['salmon_path'] = os.path.dirname(cand)
            break
    mf = os.path.join(d, 'manifest.json')
    if os.path.isfile(mf):
        try:
            with open(mf, encoding='utf-8') as fh:
                lib['manifest'] = json.load(fh)
        except (OSError, ValueError):
            pass
    return lib if (lib['minibwa'] or lib['salmon']) else None


@app.route('/api/kv_index_list')
def api_kv_index_list():
    """列出可用的鉴定库（virus-db 下的 + platform.json 登记的外部目录）与引擎可用性。"""
    base = DIRS['virus_src']
    out = {'base': base, 'libs': [], 'engines': {}}
    try:
        from known_virus_suite.kv_common import ToolRegistry
        reg = ToolRegistry()
        reg.probe()
        for e in ('minibwa', 'salmon'):
            out['engines'][e] = {'available': reg.has(e),
                                 'path': str(reg.paths.get(e) or '')}
    except Exception as exc:                    # noqa: BLE001
        out['engines_error'] = str(exc)
    seen = set()
    if os.path.isdir(base):
        for nm in sorted(os.listdir(base)):
            d = os.path.join(base, nm)
            if not os.path.isdir(d):
                continue
            lib = _kv_index_probe(d, nm)
            if lib is None and nm == 'kv_index':
                # 默认库即使索引缺失也列出（便于排查），但不计入可用
                lib = {'name': nm, 'path': d, 'minibwa': False,
                       'salmon': False, 'manifest': None}
            if lib is None:
                continue
            lib['external'] = False
            out['libs'].append(lib)
            seen.add(os.path.normpath(os.path.abspath(d)))
    # 外部登记目录（platform.json extra_kv_indexes，可为平台外任意位置，只读使用）
    for p in getattr(cfg, 'extra_kv_indexes', []) or []:
        n = os.path.normpath(os.path.abspath(p))
        if n in seen or not os.path.isdir(n):
            continue
        lib = _kv_index_probe(n)
        if lib is None:
            continue
        lib['external'] = True
        out['libs'].append(lib)
        seen.add(n)
    out['extra'] = [p for p in (getattr(cfg, 'extra_kv_indexes', []) or [])]
    return jsonify(out)


@app.route('/api/kv_index_register', methods=['POST'])
def api_kv_index_register():
    """登记 / 移除外部病毒鉴定库目录（持久化；只读使用，允许平台外任意位置）。

    与 /api/register_virus_lib 同语义：把平台外已建好的 minibwa/salmon
    索引目录挂进平台，之后在工具卡的「鉴定库」下拉里即可选择复用。
    """
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少鉴定库目录')
    d = _abs_arg(raw)
    if body.get('remove'):
        cfg.extra_kv_indexes = [
            p for p in (getattr(cfg, 'extra_kv_indexes', []) or [])
            if os.path.normpath(os.path.abspath(p)) != d]
        cfg.save()
        return jsonify({'ok': True, 'dirs': cfg.extra_kv_indexes})
    if not os.path.isdir(d):
        abort(400, '目录不存在: ' + raw)
    lib = _kv_index_probe(d)
    if lib is None:
        abort(400, '不是有效的鉴定库目录（未找到 minibwa.mbw 或 salmon_k31/info.json）: ' + raw)
    if d not in (getattr(cfg, 'extra_kv_indexes', []) or []):
        cfg.extra_kv_indexes = list(getattr(cfg, 'extra_kv_indexes', []) or []) + [d]
        cfg.save()
    return jsonify({'ok': True, 'dirs': cfg.extra_kv_indexes,
                    'lib': {'name': lib['name'], 'path': d,
                            'minibwa': lib['minibwa'], 'salmon': lib['salmon']}})


@app.route('/api/build_universal_db', methods=['POST'])
def api_build_universal_db():
    """构建通用病毒参考库（refvirus=NCBI RefSeq Viral / rvdb=RVDB C-RVDB）。"""
    body = request.get_json(force=True) or {}
    source = body.get('source') or 'refvirus'
    if source not in ('refvirus', 'rvdb'):
        abort(400, 'source 需为 refvirus 或 rvdb')
    name = {'refvirus': '通用病毒库构建 (RefSeq Viral)',
            'rvdb': 'RVDB 库构建 (C-RVDB)'}[source]
    tid = tm.start(cfg.tr(name, name), _universal_db_job(source))
    return jsonify({'task': tid})


@app.route('/api/build_virus_db', methods=['POST'])
def api_build_virus_db():
    body = request.get_json(force=True)
    for k in ('fasta', 'info'):
        if not body.get(k):
            abort(400, f'缺少参数 {k}')
    check_path(body['fasta'], must_exist=True)
    check_path(body['info'], must_exist=True)
    tid = tm.start(cfg.tr('病毒库构建', 'Virus DB build'), _job_build_virus_db(body))
    return jsonify({'task': tid})


# ------------------------------------------------------------------
# API: 样品分析
# ------------------------------------------------------------------
def _analysis_kwargs(body):
    """analyze / pipeline run 共用的 run_analysis 参数装配。

    参数优先级：请求 body > 设置页默认值(cfg.defaults) > 代码硬编码缺省。
    """
    d = cfg.defaults

    def _val(key, fallback):
        v = body.get(key)
        if v is None or v == '':
            return d.get(key, fallback)
        return v

    def _chk(key, fallback=False):
        v = body.get(key)
        if v is None:
            return bool(d.get(key, fallback))
        return bool(v)

    def _norm_db(v):
        if not v:
            return None
        return v if os.path.isabs(str(v)) else os.path.join(PLATFORM_ROOT, v)
    def _norm_methods(v):
        """验证证据方法：接受列表/逗号串，只保留 blastx/cdd，空则默认双路。"""
        if v is None or v == '':
            v = d.get('verify_methods') or ['blastx', 'cdd']
        if isinstance(v, str):
            v = [x for x in v.replace(';', ',').split(',') if x.strip()]
        keep = tuple(x for x in (str(i).strip() for i in v)
                     if x in ('blastx', 'cdd'))
        return keep or ('blastx', 'cdd')

    return dict(
        db_host=_norm_db(body.get('db_host')) or cfg.databases['host'],
        db_virus=_norm_db(body.get('db_virus')) or cfg.databases['virus'],
        chunk_dir=_norm_db(body.get('chunk_dir')),
        fastp_dedup=_chk('fastp_dedup'),
        do_fq2fa=_chk('do_fq2fa', True),
        gbdraw_max=int(_val('gbdraw_max', 12) or 12),
        gbdraw_fasta=body.get('gbdraw_fasta') or None,
        gbdraw_ann=body.get('gbdraw_ann') or None,
        plot_engine=_val('plot_engine', 'auto'),
        threads=int(body.get('threads') or cfg.threads),
        confidence=float(_val('confidence', 0) or 0),
        assembly_mode=_val('assembly_mode', 'metaviral'),
        assembly_input=_val('assembly_input', 'virus'),
        memory_gb=int(_val('memory', 64) or 64),
        subsample=int(_val('subsample', 0) or 0),
        min_contig_len=int(_val('min_contig_len', 200) or 200),
        min_orf_aa=int(_val('min_orf_aa', 100) or 100),
        top_n_refs=int(_val('top_n_refs', 10) or 10),
        tree_tool=_val('tree_tool', 'fasttree'),
        tree_sampling=_val('tree_sampling', 'blast'),
        ncbi_refs=body.get('ncbi_refs') or None,
        primer_mode=_val('primer_mode', 'conserved'),
        do_trim=_chk('do_trim', True),
        do_specificity=_chk('specificity'),
        force=_chk('force'),
        # ③b 候选序列验证（宿主类群 / 证据方法 / 并集-交集）
        verify_host=_val('verify_host', 'all'),
        verify_methods=_norm_methods(body.get('verify_methods')),
        verify_combine=_val('verify_combine', 'union'),
    )


@app.route('/api/analyze', methods=['POST'])
def api_analyze():
    body = request.get_json(force=True)
    if not body.get('r1'):
        abort(400, '缺少 R1')
    check_path(body['r1'], must_exist=True)
    if body.get('r2'):
        check_path(body['r2'], must_exist=True)
    from vp.pipeline import DEFAULT_ANALYZE_STAGES
    stages = body.get('stages') or list(DEFAULT_ANALYZE_STAGES)
    from vp.pipeline import _safe_sample_name
    sample = _safe_sample_name(body.get('sample') or _default_sample_name(body['r1']))
    for t in tm.list_all():
        if t.get('status') == 'running' and _task_matches_sample(t.get('name'), sample):
            abort(400, f'样品 {sample} 已有任务在运行（{t.get("name")}），'
                       f'请等它结束或先取消后再试')

    def job(log, prog, cancel):
        from vp.pipeline import run_analysis
        logger = TaskLogger(callback=log)
        result_dir = run_analysis(
            sample=sample, r1=body['r1'], r2=body.get('r2'),
            stages=stages, **_analysis_kwargs(body),
            logger=logger, progress=prog)
        logger.close()
        return result_dir

    log_dir = check_path(os.path.join(DIRS['results'], sample, 'logs'),
                         must_exist=False, in_platform=True)
    os.makedirs(log_dir, exist_ok=True)
    tid = tm.start(cfg.tr(f'样品分析 {sample}', f'Sample analysis {sample}'),
                   job,
                   log_file=os.path.join(
                       log_dir,
                       f'analyze_{time.strftime("%Y%m%d_%H%M%S")}.log'))
    return jsonify({'task': tid, 'sample': sample})


# ------------------------------------------------------------------
# API: 管道（模块化逐级运行）
# ------------------------------------------------------------------
def _sample_dir(sample):
    from vp.pipeline import _safe_sample_name
    from vp.config import DIRS
    s = _safe_sample_name(sample)
    return s, check_path(os.path.join(DIRS['results'], s),
                         must_exist=True, in_platform=True)


@app.route('/api/samples')
def api_samples():
    from vp.pipeline import (pipeline_overview, load_sample_input,
                             load_project_manifest)
    from vp.config import DIRS
    out = []
    res = check_path(DIRS['results'], must_exist=False, in_platform=True)
    if os.path.isdir(res):
        for name in sorted(os.listdir(res)):
            if name.startswith('_'):      # _archive / 内部测试样品不展示
                continue
            sd = check_path(os.path.join(res, name), must_exist=False,
                            in_platform=True)
            if not os.path.isdir(sd):
                continue
            try:
                stages = pipeline_overview(sd)['stages']
            except Exception:
                stages = []
            try:
                _r1, _r2, project = load_sample_input(sd)
            except Exception:
                project = None
            # 清单总是要读：即便 input.json 已带项目名，
            # last_status / last_run 也只能从 project.json 取。
            try:
                man = load_project_manifest(sd)
            except Exception:
                man = {}
            out.append({'name': name,
                        'project': project or man.get('project') or '',
                        'last_status': man.get('last_status') or '',
                        'last_run': man.get('last_run') or '',
                        'done': sum(1 for s in stages if s['status'] == 'done'),
                        'total': len(stages) or 7})
    return jsonify(out)


@app.route('/api/samples/archived')
def api_samples_archived():
    """归档样品列表（results/_archive/ 下的历史项目；结果中心展示）。"""
    arch = check_path(os.path.join(DIRS['results'], '_archive'),
                      must_exist=False, in_platform=True)
    out = []
    if os.path.isdir(arch):
        for name in sorted(os.listdir(arch)):
            sd = os.path.join(arch, name)
            if not os.path.isdir(sd):
                continue
            has_report = os.path.isfile(os.path.join(sd, '07_report',
                                                     'report.html'))
            try:
                mt = max((os.path.getmtime(os.path.join(dp, fn))
                          for dp, _, fns in os.walk(sd) for fn in fns),
                         default=0)
            except OSError:
                mt = 0
            out.append({'name': name, 'has_report': has_report,
                        'mtime': mt})
    return jsonify(out)


@app.route('/archive_report/<sample>/')
def page_archive_report(sample):
    """归档样品报告（results/_archive/<sample>/）。"""
    safe = _safe_sample(sample)
    rpt = check_path(os.path.join(DIRS['results'], '_archive', safe,
                                  '07_report', 'report.html'),
                     must_exist=True, in_platform=True)
    return send_file(check_path(rpt, must_exist=True, in_platform=True))


def _default_sample_name(path):
    """从文件名推导默认样品名：去扩展名与常见测序后缀。

    Lycium_barbarum_R1.fastq.gz -> Lycium_barbarum
    NX-5_1.fq.gz                -> NX-5
    sample.fq                   -> sample
    """
    base = os.path.basename(str(path))
    base = re.sub(r'\.(fastq|fq|fasta|fa|fna)(\.gz)+$', '', base, flags=re.I)
    base = re.sub(r'\.(fastq|fq|fasta|fa|fna)$', '', base, flags=re.I)
    # 循环剥离技术后缀（Illumina 常见 _L001_R1_001 之类为多层）
    pat = re.compile(
        r'(?:[_\-.]?R[12]|_\d{2,3}|[_\-.][12]|[_\-.]L00\d)$', re.I)
    while True:
        stripped = pat.sub('', base)
        if stripped == base or not stripped:
            break
        base = stripped
    return base or 'sample'


@app.route('/api/pipeline/create', methods=['POST'])
def api_pipeline_create():
    body = request.get_json(force=True)
    if not body.get('r1'):
        abort(400, '缺少 R1')
    check_path(body['r1'], must_exist=True)
    if body.get('r2'):
        check_path(body['r2'], must_exist=True)
    from vp.pipeline import _safe_sample_name, save_sample_input
    from vp.config import DIRS
    sample = _safe_sample_name(body.get('sample') or _default_sample_name(body['r1']))
    sd = check_path(os.path.join(DIRS['results'], sample),
                    must_exist=False, in_platform=True)
    if os.path.isdir(sd) and any(os.scandir(sd)):
        abort(400, f'样品 {sample} 已存在，请换一个样品名')
    os.makedirs(sd, exist_ok=True)
    # 可选：创建时立即截取子样本作为管道输入（大样品先小规模验证）
    sub = int(body.get('subsample', 0) or 0)
    in1, in2 = body['r1'], body.get('r2')
    if sub > 0:
        from vp.pipeline import subsample_fastq
        in1, in2 = subsample_fastq(in1, in2, os.path.join(sd, '00_prep'), sub)
    save_sample_input(sd, in1, in2, sample,
                      project=str(body.get('project') or '') or None)
    return jsonify({'sample': sample})


@app.route('/api/pipeline/<sample>')
def api_pipeline(sample):
    from vp.pipeline import pipeline_overview, load_sample_input
    s, sd = _sample_dir(sample)
    r1, r2, project = load_sample_input(sd)
    lang = request.args.get('lang') or None
    return jsonify({'sample': s, 'r1': r1, 'r2': r2, 'project': project or '',
                    **pipeline_overview(sd, lang)})


@app.route('/api/pipeline/<sample>/run', methods=['POST'])
def api_pipeline_run(sample):
    body = request.get_json(force=True)
    s, sd = _sample_dir(sample)
    # 防重复：该样品已有运行中任务时拒绝（避免并发分类写爆磁盘）
    for t in tm.list_all():
        if t.get('status') == 'running' and _task_matches_sample(t.get('name'), s):
            abort(400, f'样品 {s} 已有任务在运行（{t.get("name")}），'
                       f'请等它结束或先取消后再试')
    from vp.pipeline import load_sample_input, run_analysis, STAGE_ORDER, STAGE_NAMES
    r1, r2, _proj = load_sample_input(sd)
    if not r1:
        abort(400, '样品缺少输入记录（input.json），请重新创建样品')
    check_path(r1, must_exist=True)
    stages = body.get('stages') or ([body['stage']] if body.get('stage') else [])
    stages = [st for st in stages if st in STAGE_ORDER]
    # fastp 未安装时剔除（可选步骤）
    if 'fastp' in stages:
        from vp.preprocess import fastp_available
        if not fastp_available():
            stages = [st for st in stages if st != 'fastp']
    if not stages:
        # 未指定阶段 = 「依次运行剩余步骤」：所有未完成且可用的阶段
        try:
            from vp.pipeline import pipeline_overview
            ov = pipeline_overview(sd)['stages']
            done = {s['stage'] for s in ov if s['status'] == 'done'}
            stages = [s['stage'] for s in ov
                      if s['stage'] in STAGE_ORDER
                      and s['stage'] not in done
                      and s['status'] != 'unavailable']
        except Exception:
            stages = []
    if not stages:
        abort(400, '没有需要运行的阶段（全部已完成；重跑请勾选「强制重跑」'
                   '并指定具体步骤）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        run_analysis(sample=s, r1=r1, r2=r2, stages=stages,
                     **_analysis_kwargs(body), logger=logger, progress=prog)
        logger.close()
        return sd

    from vp.pipeline import stage_name
    label = ' → '.join(stage_name(st, cfg.lang) for st in stages)
    log_dir = check_path(os.path.join(sd, 'logs'), must_exist=False,
                         in_platform=True)
    os.makedirs(log_dir, exist_ok=True)
    tid = tm.start(f'{s} · {label}', job,
                   log_file=os.path.join(
                       log_dir,
                       f'run_{time.strftime("%Y%m%d_%H%M%S")}.log'))
    return jsonify({'task': tid, 'sample': s})


# ------------------------------------------------------------------
# API: 其他动作
# ------------------------------------------------------------------
@app.route('/api/open_platform_dir')
def api_open_platform_dir():
    os.startfile(check_path(PLATFORM_ROOT, must_exist=True, in_platform=True))
    return jsonify({'ok': True})


@app.route('/help')
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


@app.route('/api/open_report_dir/<sample>')
def api_open_report_dir(sample):
    safe = _safe_sample(sample)
    d = check_path(os.path.join(DIRS['results'], safe), must_exist=True,
                   in_platform=True)
    os.startfile(check_path(d, must_exist=True, in_platform=True))
    return jsonify({'ok': True})


@app.route('/api/samples/<sample>/<path:rel>')
def api_sample_file(sample, rel):
    """结果文件下载（严格限制在样品目录内）。"""
    safe = _safe_sample(sample)
    p = check_path(os.path.join(DIRS['results'], safe, rel),
                   must_exist=True, in_platform=True)
    return send_file(check_path(p, must_exist=True, in_platform=True),
                     as_attachment=request.args.get('dl') == '1')


def _sample_busy(safe):
    """样品是否有正在运行的任务（删除/清除结果前拦截）。"""
    for t in tm.list_all():
        if t.get('status') == 'running' and _task_matches_sample(
                t.get('name'), safe):
            return True
    return False


@app.route('/api/samples/<sample>/clear', methods=['POST'])
def api_sample_clear(sample):
    """清除样品的全部分析结果文件，保留样品登记（00_prep/input.json）。

    与 delete 的区别：样品仍留在样品列表中，输入档案不丢，
    之后再跑管道即从第一步重新分析。样品有运行中任务时拒绝。
    """
    import shutil
    safe = _safe_sample(sample)
    if _sample_busy(safe):
        abort(400, f'样品 {safe} 有任务在运行，请先取消再清除')
    d = check_path(os.path.join(DIRS['results'], safe), must_exist=True,
                   in_platform=True)
    # 双保险：必须是 results/ 的直接子目录
    if os.path.dirname(os.path.abspath(d)) != os.path.abspath(DIRS['results']):
        abort(400, '仅允许操作样品目录')
    prep = os.path.join(d, '00_prep')
    # 顶层两个档案文件与结果无关，清除结果时保留：
    #   input.json   —— 样品输入登记（在 00_prep 内）
    #   project.json —— 项目清单（项目身份/输入/参数）
    keep_top = {'project.json'}
    # 二次检查：上面的 _sample_busy 与这里之间有读档案/建目录等磁盘操作，
    # 队列或手动任务可能刚好在这段窗口内启动。Windows 上 rmtree 遇到
    # 正在写入的文件会抛 PermissionError 并留下半删目录，所以必须重新确认。
    if _sample_busy(safe):
        abort(400, f'样品 {safe} 刚刚启动了任务，请先取消再清除')
    for name in os.listdir(d):
        p = os.path.join(d, name)
        if name in keep_top:
            continue
        if os.path.abspath(p) == os.path.abspath(prep):
            # 输入档案目录：只保留 input.json（fastp/fq2fa/子采样产物一并清）
            for f in os.listdir(prep):
                if f == 'input.json':
                    continue
                fp = os.path.join(prep, f)
                try:
                    if os.path.isdir(fp):
                        shutil.rmtree(fp)
                    else:
                        os.remove(fp)
                except PermissionError:
                    abort(400, f'文件 {f} 正在被任务占用，无法删除，'
                               f'请等待任务结束或先取消')
        elif os.path.isdir(p):
            try:
                shutil.rmtree(p)
            except PermissionError:
                abort(400, f'目录 {name} 正在被任务占用，无法删除，'
                           f'请等待任务结束或先取消')
        else:
            try:
                os.remove(p)
            except PermissionError:
                abort(400, f'文件 {name} 正在被任务占用，无法删除，'
                           f'请等待任务结束或先取消')
    # 结果已清空，清单里的结果类字段必须一并归零，
    # 否则 project.json 会声称还有已完成的阶段（与卡片状态矛盾）。
    # 仅保留项目身份：project / input / params / schema / created。
    try:
        from vp.pipeline import load_project_manifest, save_project_manifest
        _m = load_project_manifest(d)
        if _m:
            for _k in ('stages_done', 'runs', 'last_status', 'last_stage',
                       'last_duration_s', 'last_run'):
                _m.pop(_k, None)
            save_project_manifest(d, _m)
    except Exception as e:
        # 清单归零失败会让 project.json 与已清空的结果目录矛盾（卡片仍显示阶段已完成）
        app.logger.warning('样品清单归零失败 %s: %s', safe, e)
    return jsonify({'ok': True})


@app.route('/api/samples/<sample>/delete', methods=['POST'])
def api_sample_delete(sample):
    """删除样品：结果目录（含 00_prep 输入档案与全部产物）一并移除。

    仅允许平台 results/ 内的一级样品目录；样品有运行中任务时拒绝。
    只想清空结果、保留样品时用 /clear。
    """
    import shutil
    safe = _safe_sample(sample)
    if _sample_busy(safe):
        abort(400, f'样品 {safe} 有任务在运行，请先取消再删除')
    d = check_path(os.path.join(DIRS['results'], safe), must_exist=True,
                   in_platform=True)
    # 双保险：必须是 results/ 的直接子目录
    if os.path.dirname(os.path.abspath(d)) != os.path.abspath(DIRS['results']):
        abort(400, '仅允许删除样品目录')
    # 二次检查：上面的 _sample_busy 与这里之间有路径校验等操作，任务可能
    # 刚好在这段窗口内启动。Windows 上 rmtree 遇到正在写入的文件会抛
    # PermissionError 并留下半删目录。
    if _sample_busy(safe):
        abort(400, f'样品 {safe} 刚刚启动了任务，请先取消再删除')
    try:
        shutil.rmtree(d)
    except PermissionError:
        abort(400, f'样品目录 {safe} 正在被任务占用，无法删除，'
                   f'请等待任务结束或先取消')
    return jsonify({'ok': True})


@app.route('/api/sample_files/<sample>')
def api_sample_files(sample):
    """列出样品结果文件树（相对路径）。"""
    safe = _safe_sample(sample)
    root = check_path(os.path.join(DIRS['results'], safe), must_exist=True,
                      in_platform=True)
    out = []
    for dirpath, _dirs, fnames in os.walk(root):
        for fn in fnames:
            if fn.startswith('.') or fn == 'plotly.min.js':
                continue
            full = check_path(os.path.join(dirpath, fn), must_exist=False,
                              in_platform=True)
            if not os.path.isfile(full):
                continue
            relp = os.path.relpath(full, root)
            out.append({'path': relp.replace(os.sep, '/'),
                        'size': fmt_size(os.path.getsize(full))})
    out.sort(key=lambda x: x['path'])
    return jsonify({'files': out})


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

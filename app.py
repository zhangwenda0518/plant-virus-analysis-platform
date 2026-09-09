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
from vp.web import settings as _bp_settings  # noqa: E402
from vp.web import virome as _bp_virome      # noqa: E402

app.register_blueprint(_bp_settings.bp)
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
# 任务引擎
# ------------------------------------------------------------------


def _slot_limits():
    """并发闸门上限（platform.json defaults 可覆盖：max_heavy_tasks /
    max_light_tasks）；未配置时用缺省值。

    heavy = 吃满多核/GB 级内存的任务（kunpeng 分类、SPAdes 组装、
    DIAMOND 注释、建库、建树、SDT）；light = 网络 IO 或轻量计算
    （下载、格式转换、绘图、检索、Selenium 提交）。
    """
    d = getattr(cfg, 'defaults', None) or {}
    try:
        heavy = max(1, min(int(d.get('max_heavy_tasks', 2)), 8))
    except (TypeError, ValueError):
        heavy = 2
    try:
        light = max(1, min(int(d.get('max_light_tasks', 4)), 16))
    except (TypeError, ValueError):
        light = 4
    return heavy, light


# 内存保留策略：最近 _KEEP_FULL 个任务保留完整信息，更早的清空日志/结果引用
# 只留状态摘要；超过 _KEEP_TOTAL 的移出内存（磁盘 tasks/*.json 仍可回溯）。
_KEEP_FULL = 60
_KEEP_TOTAL = 300


class TaskManager:
    """后台线程任务 + 状态快照（内存为主，tasks/ 目录留一份 JSON）。

    并发控制：heavy / light 两档 Semaphore。任务先排队（status='queued'，
    对外仍呈现为 running，仅 msg 标注排队位次），拿到名额后才真正开跑。
    """

    def __init__(self):
        self.tasks = {}
        self.order = []
        self.lock = threading.Lock()
        heavy, light = _slot_limits()
        self.heavy_slots = threading.Semaphore(heavy)
        self.light_slots = threading.Semaphore(light)
        self.limits = {'heavy': heavy, 'light': light}
        self._recover_interrupted()

    def _recover_interrupted(self):
        """服务启动时把上次进程遗留的 running 任务标记为中断（防幽灵），
        并把任务状态文件修剪到最近 100 个（防 tasks/ 无限膨胀）。"""
        try:
            files = [n for n in os.listdir(DIRS['tasks']) if n.endswith('.json')]
            files.sort(reverse=True)                     # 新任务 id 排前（uuid 时间无关，
                                                         # 改按 mtime 更稳，见下）
            files.sort(key=lambda n: os.path.getmtime(
                os.path.join(DIRS['tasks'], n)), reverse=True)
            for fn in files[100:]:                       # 只保留最近 100 个
                try:
                    safe_remove(os.path.join(DIRS['tasks'], fn))
                except (OSError, ValueError):
                    pass
            files = files[:100]
            for name in files:
                p = check_path(os.path.join(DIRS['tasks'], name),
                               must_exist=False, in_platform=True)
                try:
                    with safe_open(p) as f:
                        d = json.load(f)
                except (OSError, ValueError):
                    continue
                if d.get('status') in ('running', 'queued'):
                    d['status'] = 'failed'
                    d['error'] = '服务重启导致任务中断，请重新运行'
                    d['finished'] = time.time()
                    with safe_open(p, 'wt') as f:
                        json.dump(d, f, ensure_ascii=False)
        except Exception:
            pass

    def start(self, name, fn, log_file=None, weight='heavy', link=None):
        """创建任务。log_file: 任务日志持久化文件（缺省落 logs/tasks/）。

        weight: 'heavy'（吃多核/大内存，受 heavy 闸门限流）或 'light'
        （网络 IO / 轻量计算，受 light 闸门限流）。缺省 heavy——宁可多排队，
        也不要把不认识的任务放过去打满机器。

        link: 结果/操作页跳转地址（任务卡「前往」按钮）。

        重启上下文自动捕获：任务在 POST API handler 里发起时，记录
        request.path + JSON body，任务卡可一键重新提交同一请求。
        仅存内存（不落盘），避免 API key 等敏感字段写入 tasks/*.json。

        log() 同时写内存 deque（前端实时读）与磁盘文件（重启后可查）。
        """
        tid = uuid.uuid4().hex[:12]
        restart = None
        try:
            from flask import has_request_context, request as _req
            if has_request_context() and _req.method == 'POST':
                body = _req.get_json(silent=True)
                restart = {'url': _req.path, 'body': body} if body else None
        except Exception:
            restart = None
        if not log_file:
            safe = re.sub(r'[^\w\-.]+', '_', name, flags=re.UNICODE)
            safe = safe[:40].strip('_ ') or 'task'
            log_dir = check_path(os.path.join(DIRS['logs'], 'tasks'),
                                 must_exist=False, in_platform=True)
            os.makedirs(log_dir, exist_ok=True)
            log_file = os.path.join(
                log_dir, f'{safe}_{time.strftime("%Y%m%d_%H%M%S")}_{tid}.log')
        rec = {'id': tid, 'name': name, 'status': 'queued', 'stage': '',
               'pct': 0.0, 'msg': '', 'error': '', 'log': deque(maxlen=500),
               'started': time.time(), 'finished': None, 'result': None,
               'result_preview': None, 'eta': None,
               'cancel': threading.Event(),
               'procs': [],          # 本任务启动的子进程（硬停止用）
               'weight': weight if weight in ('heavy', 'light') else 'heavy',
               'log_file': log_file, 'link': link, 'restart': restart}
        with self.lock:
            self.tasks[tid] = rec
            self.order.insert(0, tid)

        def _log(line):
            rec['log'].append(line)
            try:
                with safe_open(rec['log_file'], 'at') as f:
                    f.write(time.strftime('[%H:%M:%S] ') + line + '\n')
            except Exception:
                pass

        def _progress(stage, pct, msg, *eta):
            rec['stage'], rec['pct'], rec['msg'] = stage, pct, msg
            if eta:
                try:
                    rec['eta'] = max(float(eta[0]), 0)
                except (TypeError, ValueError):
                    rec['eta'] = None
            self._persist(rec)

        def _run():
            from vp.utils import task_bind
            slot = (self.light_slots if rec['weight'] == 'light'
                    else self.heavy_slots)
            # 排队等名额：每 0.5s 探一次，期间可被取消（不会被闸门永久卡住）
            got = False
            while not rec['cancel'].is_set():
                if slot.acquire(timeout=0.5):
                    got = True
                    break
                rec['msg'] = self._queue_msg(rec)
            if not got:                       # 排队途中被取消
                rec['status'] = 'cancelled'
                rec['error'] = '任务已停止（排队中取消）'
                rec['finished'] = time.time()
                _log('[CANCEL] ' + cfg.tr('排队中已取消',
                                          'Cancelled while queued'))
                self._persist(rec)
                return
            rec['status'] = 'running'
            rec['started'] = time.time()      # 已运行时长从真正开跑算起
            self._persist(rec)
            task_bind(rec['cancel'], rec['procs'])
            try:
                result = fn(_log, _progress, rec['cancel'])
                rec['status'] = 'cancelled' if rec['cancel'].is_set() else 'done'
                rec['result'] = result
                rec['result_preview'] = build_result_preview(result, rec['name'])
                if rec['status'] == 'done':
                    _log('[OK] ' + cfg.tr('任务完成', 'Task finished'))
            except Exception as e:
                rec['status'] = 'cancelled' if rec['cancel'].is_set() else 'failed'
                rec['error'] = ('任务已停止（用户取消）' if rec['cancel'].is_set()
                                else str(e))
                _log(f"[ERROR] {rec['error']}" if rec['cancel'].is_set()
                     else f"[ERROR] {e}")
            finally:
                task_bind(None, None)
                rec['procs'].clear()
                rec['finished'] = time.time()
                slot.release()          # 归还名额
                rec['msg'] = ''
                self._trim_memory()     # 回收旧任务内存
                self._persist(rec)

        th = threading.Thread(target=_run, daemon=True, name=f'task-{tid}')
        rec['thread'] = th
        th.start()
        self._persist(rec)
        return tid

    def _persist(self, rec):
        try:
            p = check_path(os.path.join(DIRS['tasks'], rec['id'] + '.json'),
                           must_exist=False, in_platform=True)
            data = {'id': rec['id'], 'name': rec['name'],
                    'status': rec['status'], 'stage': rec['stage'],
                    'pct': rec['pct'], 'msg': rec['msg'],
                    'log_file': rec.get('log_file'),
                    'error': rec['error'],
                    'log_tail': list(rec['log'])[-40:],
                    'started': rec['started'], 'finished': rec['finished'],
                    'eta': rec.get('eta'), 'link': rec.get('link')}
            if rec.get('result') is not None:
                try:
                    json.dumps(rec['result'], ensure_ascii=False)
                    data['result'] = rec['result']
                except (TypeError, ValueError):
                    data['result'] = str(rec['result'])
            with safe_open(p, 'wt') as f:
                json.dump(data, f, ensure_ascii=False)
        except Exception as e:
            # 状态落盘失败不影响内存态，但会让刷新后丢进度：留痕便于排查
            app.logger.warning('任务状态保存失败 %s: %s', rec.get('id'), e)

    def _queue_msg(self, rec):
        """排队提示：统计同档位里排在自己前面的任务数。"""
        ahead = 0
        with self.lock:
            for tid in self.order:
                if tid == rec['id']:
                    break
                r = self.tasks.get(tid)
                if (r and r['status'] == 'queued'
                        and r.get('weight') == rec.get('weight')):
                    ahead += 1
        if ahead:
            return cfg.tr(f'排队中 · 前面还有 {ahead} 个任务',
                          f'Queued · {ahead} task(s) ahead')
        return cfg.tr('排队中 · 等待空闲名额', 'Queued · waiting for a slot')

    def _trim_memory(self):
        """回收内存：超出 _KEEP_FULL 的终态任务清掉日志与结果引用（磁盘已
        持久化，仍可经 tasks/*.json 回溯）；超出 _KEEP_TOTAL 的移出内存。

        服务长跑时任务只增不减会持续吃内存，这里做稳态回收。
        """
        with self.lock:
            for tid in self.order[_KEEP_FULL:]:
                rec = self.tasks.get(tid)
                if not rec or rec['status'] in ('running', 'queued'):
                    continue
                if rec['log']:
                    rec['log'].clear()
                rec['result'] = None
                rec['result_preview'] = None
                rec['thread'] = None
                rec['procs'] = []
            if len(self.order) > _KEEP_TOTAL:
                for tid in self.order[_KEEP_TOTAL:]:
                    self.tasks.pop(tid, None)
                self.order = self.order[:_KEEP_TOTAL]

    def snapshot(self, tid, log_lines=80, include_result=True):
        with self.lock:
            rec = self.tasks.get(tid)
            if not rec:
                return None
            # 'queued' 对外一律呈现为 running：前端据此显示取消按钮、进度条
            # 与自动滚动；排队状态本身通过 msg 文案告知用户。
            st = 'running' if rec['status'] == 'queued' else rec['status']
            out = {'id': rec['id'], 'name': rec['name'],
                   'status': st, 'stage': rec['stage'],
                   'pct': rec['pct'], 'msg': rec['msg'],
                   'error': rec['error'], 'started': rec['started'],
                   'finished': rec['finished'], 'eta': rec.get('eta'),
                   'link': rec.get('link'), 'restart': rec.get('restart'),
                   'log': list(rec['log'])[-log_lines:]}
            if include_result:
                out['result'] = rec.get('result_preview')
        return out

    def list_all(self, limit=_KEEP_FULL * 2, logs=True):
        """任务列表。

        前端每 2.5s 轮询一次本接口，历史上对所有任务都返回 30 行日志 +
        结果预览，任务累积后 payload 持续变大。这里做瘦身：
        - 运行中/排队中的任务：完整 30 行日志（实时看进度）
        - 已结束的任务：只带最后 5 行（日志已落盘，展开任务卡可单独拉取）
        - logs=False：完全不带日志（任务中心/徒轮询只需元信息，
          日志走 /api/task/<tid>/log 单独拉取）
        - 结果预览一律保留（任务完成后的结果面板依赖它）
        """
        with self.lock:
            ids = list(self.order[:limit])
            active = {i for i in ids
                      if (self.tasks.get(i) or {}).get('status')
                      in ('running', 'queued')}
        if not logs:
            return [self.snapshot(i, log_lines=0, include_result=True)
                    for i in ids if i in self.tasks]
        return [self.snapshot(i, log_lines=(30 if i in active else 5),
                              include_result=True)
                for i in ids if i in self.tasks]

    @staticmethod
    def _kill_tree(proc):
        """杀进程树（Windows 用 taskkill /T 连子进程；其它平台 kill）。"""
        if proc.poll() is not None:
            return
        try:
            if sys.platform == 'win32':
                subprocess.run(['taskkill', '/F', '/T', '/PID', str(proc.pid)],
                               capture_output=True, timeout=15)
            else:
                proc.kill()
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass

    def cancel(self, tid):
        with self.lock:
            rec = self.tasks.get(tid)
        if rec and rec['status'] in ('running', 'queued'):
            rec['cancel'].set()
            if rec['status'] == 'queued':
                # 排队中：闸门等待循环会自行退出并归还名额
                rec['msg'] = '已取消排队'
                return True
            # 硬停止：杀掉本任务已注册的全部子进程（kunpeng/spades/mafft 等）
            for proc in list(rec.get('procs') or ()):
                self._kill_tree(proc)
            rec['msg'] = '已请求取消（当前步骤结束后停止）'
            return True
        return False


    @staticmethod
    def _tid_ok(tid):
        return bool(tid) and bool(re.fullmatch(r'[0-9a-f]{6,16}', str(tid)))

    def delete(self, tid):
        """删除任务记录（仅终态）：移出内存 + 删 tasks/<tid>.json。"""
        if not self._tid_ok(tid):
            return False, '非法任务 ID'
        with self.lock:
            rec = self.tasks.get(tid)
            if rec and rec['status'] in ('running', 'queued'):
                return False, '运行中的任务请先停止再删除'
            self.tasks.pop(tid, None)
            if tid in self.order:
                self.order.remove(tid)
        try:
            p = check_path(os.path.join(DIRS['tasks'], tid + '.json'),
                           must_exist=False, in_platform=True)
            if os.path.isfile(p):
                os.remove(p)
        except (OSError, ValueError):
            pass
        return True, 'ok'

    def clear_finished(self):
        """清空内存中的终态任务（含各自的 tasks/<tid>.json）。

        注意：不碰磁盘上的历史归档（平台重启后加载的 tasks/*.json）。
        清空归档是破坏性操作，只能由用户逐个删除。
        """
        with self.lock:
            tids = [t for t, r in self.tasks.items()
                    if r['status'] not in ('running', 'queued')]
        n = 0
        for t in tids:
            ok, _ = self.delete(t)
            if ok:
                n += 1
        return n

    def list_archived(self, limit=200):
        """归档任务（磁盘 tasks/*.json，不在内存）：平台重启后的历史回溯。"""
        items, total = self._read_archived(limit=limit, logs=True)
        return items

    def count_archived(self):
        """磁盘归档任务总数（不受 list_archived 的 limit 截断影响）。"""
        try:
            return len([f for f in os.listdir(DIRS['tasks'])
                        if f.endswith('.json') and self._tid_ok(f[:-5])])
        except OSError:
            return 0

    def _read_archived(self, limit=200, logs=True):
        items = []
        try:
            files = [f for f in os.listdir(DIRS['tasks']) if f.endswith('.json')]
        except OSError:
            return items, 0
        total = 0
        for f in files:
            tid = f[:-5]
            if tid in self.tasks or not self._tid_ok(tid):
                continue
            total += 1
            if len(items) >= limit:
                continue
            try:
                with safe_open(os.path.join(DIRS['tasks'], f)) as fh:
                    d = json.load(fh)
                if not logs:
                    d.pop('log_tail', None)
                items.append(d)
            except (OSError, ValueError):
                continue
        items.sort(key=lambda d: d.get('started') or 0, reverse=True)
        return items[:limit], total

    def full_log(self, tid, lines=400):
        """任务全量日志（尾部 lines 行）：优先内存，内存被回收/归档时读
        磁盘 log_file（json 里记录了路径），最后回退 json 内的 log_tail。"""
        if not self._tid_ok(tid):
            return None
        rec = self.tasks.get(tid)
        if rec and rec['log']:
            return list(rec['log'])[-lines:]
        log_file = (rec or {}).get('log_file')
        js = None
        if not log_file:
            try:
                with safe_open(os.path.join(DIRS['tasks'], tid + '.json')) as f:
                    js = json.load(f)
                log_file = js.get('log_file')
            except (OSError, ValueError):
                js = None
        if log_file and os.path.isfile(log_file):
            # 防御深度：归档 json 的 log_file 必须在平台目录内
            try:
                check_path(log_file, must_exist=True, in_platform=True)
            except (OSError, ValueError):
                log_file = None
        if log_file and os.path.isfile(log_file):
            tail = self._tail_log_file(log_file, lines)
            if tail is not None:
                return tail
        if js is not None:
            return (js.get('log_tail') or [])[-lines:]
        try:
            with safe_open(os.path.join(DIRS['tasks'], tid + '.json')) as f:
                return (json.load(f).get('log_tail') or [])[-lines:]
        except (OSError, ValueError):
            return None

    _logcache = {}   # {(path, size, mtime_ns, lines): [行...]} 日志尾部缓存

    @classmethod
    def _tail_log_file(cls, path, lines):
        """从日志文件尾部读 lines 行。

        长任务日志可达 MB 级，而任务中心展开时会 2s 轮询一次，
        逐行读完整个文件会持续重复磁盘 IO。这里 seek 到尾部只读
        必要字节，并以 (size, mtime_ns, lines) 作键缓存：文件未变
        时直接复用，轮询场景下完全命中。
        """
        try:
            st = os.stat(path)
        except OSError:
            return None
        key = (path, st.st_size, st.st_mtime_ns, lines)
        hit = cls._logcache.get(key)
        if hit is not None:
            return list(hit)
        try:
            size = st.st_size
            chunk = min(size, max(65536, lines * 256))
            with safe_open(path, 'rb') as f:
                if size > chunk:
                    f.seek(size - chunk)
                data = f.read()
        except OSError:
            return None
        rows = data.decode('utf-8', errors='replace').splitlines()
        if size > chunk and rows:
            rows = rows[1:]           # 首行可能被字节截断，丢弃
        tail = rows[-lines:]
        if len(cls._logcache) > 32:
            cls._logcache.pop(next(iter(cls._logcache)))
        cls._logcache[key] = tuple(tail)
        return tail


tm = TaskManager()


# ------------------------------------------------------------------
# 任务结果预览（任务完成后面板/通知展示的统一结构）
# ------------------------------------------------------------------
def _sample_result_preview(sample):
    """样品任务的完成预览：报告链接 + 关键数字 + 产物文件。"""
    from vp.pipeline import pipeline_overview, _safe_sample_name
    s = _safe_sample_name(sample)
    try:
        sd = check_path(os.path.join(DIRS['results'], s), must_exist=True,
                        in_platform=True)
    except (ValueError, FileNotFoundError):
        return {'kind': 'sample', 'sample': s}
    ov = pipeline_overview(sd)
    stages = [{'stage': x['stage'], 'name': x['name'], 'status': x['status'],
               'summary': x.get('summary') or ''} for x in ov.get('stages', [])]
    rpt = os.path.join(sd, '07_report', 'report.html')
    return {'kind': 'sample', 'sample': s,
            'report': f'/report/{s}/' if os.path.isfile(rpt) else None,
            'stages': stages}


def _tool_result_files(run_name, limit=12):
    """工具运行目录的产物清单（按修改时间倒序，取前 limit 个）。"""
    root = check_path(os.path.join(_tool_runs_root(), run_name),
                      must_exist=False, in_platform=True)
    if not os.path.isdir(root):
        return []
    items = []
    for cur, _sub, fns in os.walk(root):
        for fn in fns:
            p = os.path.join(cur, fn)
            try:
                items.append((os.path.getmtime(p), os.path.relpath(p, root),
                              os.path.getsize(p)))
            except OSError:
                continue
    items.sort(reverse=True)
    return [{'path': rel.replace(os.sep, '/'), 'size': fmt_size(sz)}
            for _mt, rel, sz in items[:limit]]


def build_result_preview(result, task_name=''):
    """把任务返回值规范化成前端可渲染的预览 dict（失败安全）。"""
    try:
        if isinstance(result, str) and result:
            # 样品任务：返回 results/<样品> 目录
            rel = os.path.relpath(result, DIRS['results'])
            if rel.startswith('..') or os.path.isabs(rel):
                return {'kind': 'raw', 'text': str(result)}
            return _sample_result_preview(rel)
        if isinstance(result, dict):
            out = {'kind': 'tool',
                   'stats': {str(k): v for k, v in result.items()
                             if isinstance(v, (str, int, float, bool))
                             or v is None}}
            run = result.get('run')
            if run:
                out['run'] = str(run)          # 前端 toolrun 需顶层 run 名
                out['files'] = _tool_result_files(str(run))
            return out
    except Exception:
        pass
    return None


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


def _safe_sample(sample):
    import re
    return re.sub(r'[^A-Za-z0-9_\-.]', '_', str(sample))


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
# Open-Virome 公共病毒组模块已拆到 vp/web/virome.py（blueprint 在文件末尾注册）
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# 工具箱：独立分析工具（不依赖样品管道，任意文件即跑即得）
# ------------------------------------------------------------------


def _load_fasta_seqs(path):
    """FASTA → {首 token id: 序列}（单条记录可达 MB 级，仅限小文件用）。"""
    import gzip
    op = gzip.open if str(path).lower().endswith('.gz') else open
    out, name, buf = {}, None, []
    with op(path, 'rt', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    out[name] = ''.join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        out[name] = ''.join(buf)
    return out


def _find_contig_seq(run_dir, contig):
    """在工具运行目录的各 FASTA 中查找 contig 序列。"""
    for cur, _sub, fns in os.walk(run_dir):
        if os.path.join('analysis') in cur:
            continue
        for fn in fns:
            if not fn.lower().endswith(('.fasta', '.fa', '.fna', '.fas')):
                continue
            try:
                seqs = _load_fasta_seqs(os.path.join(cur, fn))
            except (OSError, ValueError):
                continue
            if contig in seqs:
                return seqs[contig]
    return None


def _extract_records(src, ids, dst):
    """从 FASTA/FASTQ（支持 .gz）提取 header 首 token 命中 ids 的记录。

    kunpeng 报告里的 read ID 经过 seqkit fq2fa 转换，已去掉 Illumina 配对
    后缀 /1、/2，匹配时对原始 header 兼容带/不带后缀两种写法；
    FASTQ 输入统一转成 FASTA 写出（与 .fasta 扩展名一致）。
    返回提取条数；ids 为空集时写出空文件。
    """
    import gzip
    name = str(src).lower()
    is_fq = name.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq'))
    op = gzip.open if name.endswith('.gz') else open

    def match(rid):
        """命中返回规范 ID（配对后缀 /1、/2 已去掉，与 kunpeng 报告一致），未命中返回 None。"""
        if rid in ids:
            return rid
        base = rid.rsplit('/', 1)
        if len(base) == 2 and base[1] in ('1', '2') and base[0] in ids:
            return base[0]
        return None

    n = 0
    with op(src, 'rt', errors='replace') as f, safe_open(dst, 'wt') as w:
        if is_fq:
            while True:
                h = f.readline()
                if not h:
                    break
                seq = f.readline().rstrip('\r\n')
                f.readline()
                f.readline()
                tok = h[1:].split()
                mid = match(tok[0]) if tok else None
                if mid:
                    w.write(f'>{mid}\n')
                    for i in range(0, len(seq), 70):
                        w.write(seq[i:i + 70] + '\n')
                    n += 1
        else:
            keep = False
            for line in f:
                if line.startswith('>'):
                    tok = line[1:].split()
                    keep = bool(tok) and match(tok[0]) is not None
                    if keep:
                        n += 1
                if keep:
                    w.write(line)
    return n


@app.route('/tools')
def page_tools():
    return render_template('tools.html')


# settings / storage 已拆到 vp/web/settings.py
# （blueprint 在文件末尾统一注册）


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


@app.route('/tool_runs/<run>/<path:filename>')
def page_tool_run_file(run, filename):
    """工具运行目录内文件的下载/预览（严格限制在该运行目录内）。"""
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run) or '..' in filename:
        abort(400, '无效的运行名')
    p = check_path(os.path.join(_tool_runs_root(), run, filename),
                   must_exist=True, in_platform=True)
    resp = send_file(p)
    # .svg 必须在浏览器内联渲染：send_file 在此环境把 .svg 推断成了
    # 非标准的 image/svg（浏览器会拒绝渲染 → broken image），改回标准 MIME。
    if str(filename).lower().endswith('.svg'):
        resp.mimetype = 'image/svg+xml'
    if request.args.get('dl'):          # 显式 ?dl=1 才强制下载；
        import urllib.parse as _up      # 不带时内联渲染（report iframe 需要）
        resp.headers['Content-Disposition'] = (
            "attachment; filename*=UTF-8''"
            + _up.quote(os.path.basename(filename)))
    return resp


@app.route('/api/tool/runs')
def api_tool_runs():
    """最近的工具运行列表（按运行目录修改时间倒序，30 条）。

    按目录名排序会让字母序靠前的前缀（contigs/genoplot…）挤掉真正
    最新的运行，这里以 mtime 为准。
    """
    root = _tool_runs_root()
    out = []
    if os.path.isdir(root):
        # 固定目录（_archive/_tmp/_scripts/_reports 等 _ 前缀）不算运行
        names = [n for n in os.listdir(root)
                 if not n.startswith('_')
                 and os.path.isdir(os.path.join(root, n))]
        try:
            names.sort(key=lambda n: os.path.getmtime(os.path.join(root, n)),
                       reverse=True)
        except OSError:
            names.sort(reverse=True)
        for name in names[:30]:
            d = check_path(os.path.join(root, name), must_exist=False,
                           in_platform=True)
            if not os.path.isdir(d):
                continue
            files = []
            for cur, _sub, fns in os.walk(d):
                rel = os.path.relpath(cur, d)
                for fn in fns:
                    p = os.path.join(cur, fn)
                    try:
                        size = os.path.getsize(p)
                    except OSError:
                        # MMseqs2 临时库等系统级特殊文件 stat 不到（WinError
                        # 1920），跳过大小即可，不让整个列表接口 500
                        continue
                    files.append({
                        'path': fn if rel == '.' else os.path.join(rel, fn),
                        'size': fmt_size(size)})
            out.append({'name': name, 'files': files[:40]})
    return jsonify(out)


@app.route('/api/tool/runs/<run>/delete', methods=['POST'])
def api_tool_run_delete(run):
    """删除一个工具运行目录（模块历史区的 🗑；运行名校验 + 平台内限定）。"""
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    p = check_path(os.path.join(_tool_runs_root(), run),
                   must_exist=True, in_platform=True)
    import shutil as _sh
    _sh.rmtree(p, ignore_errors=True)
    return jsonify({'ok': True, 'run': run})


@app.route('/api/tool/open', methods=['POST'])
def api_tool_open():
    """在资源管理器中打开某次工具运行的目录。"""
    body = request.get_json(force=True) or {}
    name = body.get('name') or ''
    if not name or '/' in name or '\\' in name or '..' in name:
        abort(400, '无效的运行名')
    d = check_path(os.path.join(_tool_runs_root(), name), must_exist=True,
                   in_platform=True)
    os.startfile(d)
    return jsonify({'ok': True})


# ------------------------------------------------------------------
# 工具箱注册表：每个工具声明 标题 / 是否需病毒库 / job 构造器。
# 新增工具三步：① 实现 _tool_job_<name>(ctx) 返回 job(log, prog, cancel)；
#              ② 在 TOOL_REGISTRY 登记；③ 前端 tools.html 加卡片 + MODULE_NAV。
# ctx 字段：p(参数) run_dir threads db_virus req(key,what) opt(key)
# ------------------------------------------------------------------

def _tool_job_convert(ctx):
    """格式转换：.sra→FASTQ/FASTA（sracha）、FASTQ→FASTA（seqkit）。"""
    inp = ctx.req('input', '输入文件')
    target = ctx.p.get('target') or 'fastq'
    if target not in ('fastq', 'fasta'):
        abort(400, '无效的目标格式')
    ilower = str(inp).lower()
    if ilower.endswith('.sra'):
        mode = 'sra'
    elif ilower.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        mode = 'fastq2fasta' if target == 'fasta' else None
    elif ilower.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
        abort(400, '输入已是 FASTA，无需转换')
    else:
        abort(400, '无法识别的输入格式（支持 .sra / .fastq[.gz]）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        out_files = []
        if mode == 'sra':
            from vp.public_data import sra_convert_engine
            engine, exe = sra_convert_engine()
            if engine != 'sracha':
                raise RuntimeError('sracha.exe 不可用，无法转换 .sra')
            threads = ctx.threads or max(2, min((os.cpu_count() or 4) // 2, 8))
            prog('sra', 0.2, f'sracha 解码 .sra → {target.upper()}')
            cmd = [exe, 'fastq', inp, '-O', ctx.run_dir,
                   '-t', str(threads), '-f', '-q']
            if target == 'fasta':
                cmd.append('--fasta')
            import subprocess as _sp
            _sp.run(cmd, check=True, timeout=8 * 3600,
                    stdout=_sp.DEVNULL, stderr=_sp.PIPE)
            base = os.path.basename(inp)[:-4]
            out_files = sorted(
                os.path.join(ctx.run_dir, f_) for f_ in os.listdir(ctx.run_dir)
                if f_.startswith(base + '_')
                and f_.endswith(('.fastq.gz', '.fq.gz', '.fa.gz', '.fasta.gz')))
            if not out_files:
                raise RuntimeError('sracha 无输出')
        else:
            seqkit = cfg2.tool('seqkit')
            dst = os.path.join(ctx.run_dir,
                               os.path.basename(inp).rsplit('.', 2)[0] + '.fa.gz')
            prog('fq2fa', 0.3, 'seqkit fq2fa 转换中')
            from vp.utils import run_cmd
            run_cmd([seqkit, 'fq2fa', '-w', '0',
                     '-j', str(ctx.threads or cfg2.threads), inp,
                     '-o', dst], logger=logger)
            out_files = [dst]
        prog('done', 1.0, f'完成：{len(out_files)} 个文件')
        logger.close()
        return {'n_files': len(out_files),
                'files': [os.path.relpath(f, ctx.run_dir).replace(os.sep, '/')
                          for f in out_files]}
    return job


def _tool_job_fastp(ctx):
    """① 质控预处理（fastp，单/双端）。"""
    r1 = ctx.req('r1', 'R1 FASTQ')
    r2 = ctx.opt('r2')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.preprocess import run_fastp
        prog('fastp', 0.3, 'fastp 质控中')
        res = run_fastp(ctx.run_dir, r1, r2, threads=ctx.threads, logger=logger,
                        force=True, dedup=bool(ctx.p.get('dedup')))
        if not res:
            raise RuntimeError('未检测到 fastp.exe，无法运行质控')
        prog('fastp', 1.0, '完成')
        logger.close()
        return res
    return job


def _tool_job_hostremoval(ctx):
    """宿主去除与序列提取（kunpeng 宿主库分类，独立模块，不依赖样品管道）。

    C 行 = 宿主 read 对，剔除后保留非宿主 reads（kept_R1/R2.fastq.gz）。
    """
    r1 = ctx.req('r1', 'R1 FASTQ')
    r2 = ctx.opt('r2')
    conf = float(ctx.p.get('confidence') or 0)
    db_host = ctx.opt('db') or cfg.databases['host']

    from vp.kunpeng import db_ready
    if not db_ready(db_host):
        abort(400, '宿主库未就绪，请先到「数据库构建」页构建宿主库')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.host_removal import remove_host
        prog('classify', 0.05, 'kunpeng 宿主库分类中')
        res = remove_host(ctx.run_dir, r1, r2, db_host, threads=ctx.threads,
                          confidence=conf, logger=logger, force=True,
                          progress=lambda pct, msg: prog(
                              'classify' if pct < 0.9 else 'filter', pct, msg))
        prog('done', 1.0, '完成')
        logger.close()
        return res
    return job


def _tool_job_hostpredict(ctx):
    """宿主预测（ICTV 级联 + NCBI 元数据交叉），独立模块。

    输入 = 病毒 contig 分类表 TSV：管道③ virus_contigs.tsv 或
    工具④ virus_classification.tsv（列名自动识别，后者归一为 ③ 口径），
    可选配套 viral_contigs.fasta。
    """
    import csv
    import shutil
    tsv = ctx.req('tsv', '病毒 contig 分类表 TSV')
    fa = ctx.opt('fasta')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        a_dir = check_path(os.path.join(ctx.run_dir, '03_assembly'),
                           must_exist=False, in_platform=True)
        os.makedirs(a_dir, exist_ok=True)
        prog('prep', 0.05, '整理输入')
        norm = os.path.join(a_dir, 'virus_contigs.tsv')
        with safe_open(tsv) as f:
            header = f.readline().rstrip('\r\n').split('\t')
        if 'kunpeng_taxid' in header:
            shutil.copyfile(tsv, norm)
        else:
            cols = ['contig', 'length', 'kunpeng_flag', 'kunpeng_taxid',
                    'kunpeng_species', 'blast_top_hit', 'blast_identity(%)',
                    'blast_coverage_hsp(%)', 'blast_aln_len', 'blast_species',
                    'blast_family']
            with safe_open(tsv) as f, safe_open(norm, 'wt') as w:
                w.write('\t'.join(cols) + '\n')
                for r in csv.DictReader(f, delimiter='\t'):
                    kt = str(r.get('taxid') or r.get('kunpeng_taxid')
                             or '').strip()
                    w.write('\t'.join([
                        str(r.get('contig') or '').strip(),
                        str(r.get('length') or '').strip(),
                        'C' if kt.isdigit() and int(kt) > 0 else 'U',
                        kt,
                        str(r.get('taxon') or r.get('species') or '').strip(),
                        '', '', '', '', '', '']) + '\n')
        if fa:
            shutil.copyfile(fa, os.path.join(a_dir, 'viral_contigs.fasta'))
        from vp.host_analysis import predict_hosts
        prog('predict', 0.15, 'ICTV 宿主概率级联预测')
        res = predict_hosts(ctx.run_dir, threads=ctx.threads, logger=logger,
                            force=True)
        prog('done', 1.0, f"完成：宿主判定 {res.get('n_contigs', 0)} 条")
        logger.close()
        return res
    return job


def _tool_job_orf(ctx):
    """ORF 预测（pyrodigal / pyrodigal_rv），可选 ⑥b 功能注释。"""
    fasta = ctx.req('fasta', '输入 FASTA（核酸 contigs / 基因组）')
    min_aa = max(1, min(int(ctx.p.get('min_aa') or 100), 5000))
    annotate = bool(ctx.p.get('annotate'))
    orf_tool = (ctx.p.get('orf_tool') or '').strip()
    orf_engine = (ctx.p.get('engine') or '').strip()
    orf_db = (ctx.p.get('db') or '').strip()
    # 功能注释 CDS 模型：pyrodigal_rv（默认）/ pyrodigal
    orf_model = (ctx.p.get('model') or '').strip()
    if orf_model not in ('pyrodigal_rv', 'pyrodigal'):
        orf_model = ''

    def job(log, prog, cancel):
        import json as _json
        import shutil
        logger = TaskLogger(callback=log)
        a_dir = check_path(os.path.join(ctx.run_dir, '03_assembly'),
                           must_exist=False, in_platform=True)
        os.makedirs(a_dir, exist_ok=True)
        prog('prep', 0.03, '整理输入')
        from vp.utils import iter_fasta, write_fasta_record
        ids = []
        vfa = os.path.join(a_dir, 'viral_contigs.fasta')
        cfa = os.path.join(a_dir, 'contigs.filtered.fasta')
        with safe_open(vfa, 'wt') as w, safe_open(cfa, 'wt') as wc:
            for h, s in iter_fasta(fasta):
                cid = h.split()[0]
                ids.append(cid)
                write_fasta_record(w, cid, s)
                write_fasta_record(wc, cid, s)
        if not ids:
            raise RuntimeError('输入 FASTA 中没有序列')
        with safe_open(os.path.join(a_dir, 'summary.json'), 'wt') as f:
            _json.dump({'viral_contigs': ids, 'standalone': True}, f)
        from vp.orf import predict_orfs
        prog('orf', 0.08, 'ORF 基因预测')
        res = predict_orfs(ctx.run_dir, min_aa=min_aa, threads=ctx.threads,
                           logger=logger, force=True, tools=orf_tool or None,
                           progress=lambda p, m: prog(
                               'orf', 0.08 + p * 0.6, m))
        res = dict(res)
        if annotate:
            from vp.orf_annot import run_orf_annotation
            prog('orfa', 0.72, 'ORF 功能注释')
            res['orfa'] = run_orf_annotation(
                ctx.run_dir, threads=ctx.threads, logger=logger, force=True,
                engine=orf_engine or None, db=orf_db or None,
                model=orf_model or None,
                progress=lambda p, m: prog('orfa', 0.72 + p * 0.26, m))
        prog('done', 1.0, '完成')
        logger.close()
        return res
    return job


def _tool_job_orfa(ctx):
    """功能注释（独立模块）：对已有 orf_ 运行注释，或 FASTA 预测+注释一步完成。"""
    run = (ctx.p.get('run') or '').strip()
    orf_engine = (ctx.p.get('engine') or '').strip()
    orf_db = (ctx.p.get('db') or '').strip()
    # 功能注释 CDS 模型：pyrodigal_rv（默认）/ pyrodigal
    orf_model = (ctx.p.get('model') or '').strip()
    if orf_model not in ('pyrodigal_rv', 'pyrodigal'):
        orf_model = ''
    if run:
        if not re.fullmatch(r'[A-Za-z0-9_]+', run) or not run.startswith('orf_'):
            abort(400, f'无效的 ORF 运行名: {run}')
        target = check_path(os.path.join(_tool_runs_root(), run),
                            must_exist=True, in_platform=True)

        def job(log, prog, cancel):
            logger = TaskLogger(callback=log)
            from vp.orf_annot import run_orf_annotation
            prog('orfa', 0.15, f'对运行 {run} 做 ORF 功能注释')
            res = run_orf_annotation(target, threads=ctx.threads, logger=logger,
                                     force=True, engine=orf_engine or None,
                                     db=orf_db or None,
                                     model=orf_model or None,
                                     progress=lambda p, m: prog(
                                         'orfa', 0.15 + p * 0.8, m))
            res = dict(res)
            res['run'] = run
            prog('done', 1.0, '完成')
            logger.close()
            return res
        return job
    # 无 run → FASTA 输入：预测 + 注释一步完成
    ctx.p = dict(ctx.p)
    ctx.p['annotate'] = True
    if not (ctx.p.get('fasta') or '').strip():
        abort(400, '请选择已有 ORF 运行或输入 FASTA')
    return _tool_job_orf(ctx)


def _tool_job_genoplot(ctx):
    """基因组图谱（gbdraw 首选，缺则 DFV 顶上）。

    输入 FASTA（可选配 GFF3 注释）或 GenBank（.gb/.gbk，自带注释）。
    """
    fasta = ctx.opt('fasta')
    ann = ctx.opt('ann')
    if not fasta and not ann:
        abort(400, '请选择 FASTA 或 GenBank 输入')
    if ann and str(ann).lower().endswith(('.gb', '.gbk', '.gbff', '.genbank')):
        fasta = None          # GenBank 自带序列与注释，FASTA 忽略
    engine = ctx.p.get('engine') or 'auto'
    max_plots = max(1, min(int(ctx.p.get('max_plots') or 12), 200))
    # 绘图定制参数（透传 gbdraw CLI）：仅收集有值/True 的项
    gb_opts = {}
    for _k, _v in (ctx.p.get('gb_opts') or {}).items():
        if _v is None or _v == '' or _v is False:
            continue
        gb_opts[str(_k)] = _v
    plot_mode = (ctx.p.get('mode') or 'both')   # circular / linear / both

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gbdraw_plot import run_genome_plots, gbdraw_available
        if not gbdraw_available():
            raise RuntimeError('未检测到 gbdraw（管道不支持 DFV）')
        tag = 'gbdraw'
        prog(tag, 0.05, f'{tag} 出图中')
        res = run_genome_plots(ctx.run_dir, logger=logger, force=True,
                               max_plots=max_plots, fasta_in=fasta, ann_in=ann,
                               opts=gb_opts, mode=plot_mode,
                               progress=lambda p, m: prog(tag, 0.05 + p * 0.9, m))
        if not res.get('plots'):
            raise RuntimeError('未产出任何基因组图（检查输入文件与绘图引擎）')
        res['run'] = os.path.basename(ctx.run_dir)   # 供前端 toolrun 内联展示 SVG
        prog('done', 1.0, f"完成：{len(res.get('plots', []))} 张图")
        logger.close()
        return res
    return job


@app.route('/api/tool/genoplot_preview', methods=['POST'])
def api_tool_genoplot_preview():
    """基因组图谱交互预览：接收输入 + 绘图参数，出图到临时目录，返回 SVG 内容。

    不写 tool_runs/、不建运行目录（交互模式下仅预览不落盘）。
    参数：fasta / ann（.gb/.gff）+ gb_opts（dict，透传 gbdraw CLI）。
    返回：{ok, svg, mode}（svg 为 circular SVG 文本；失败抛 400/500）。
    """
    import tempfile, glob as _glob
    body = request.get_json(force=True) or {}
    fasta = (body.get('fasta') or '').strip()
    ann = (body.get('ann') or '').strip()
    gb_opts = body.get('gb_opts') or {}
    if not fasta and not ann:
        abort(400, '请提供 FASTA 或 GenBank 输入')
    def _pabs(v):
        v = v.strip() if v else ''
        return check_path(v if os.path.isabs(v)
                          else os.path.join(PLATFORM_ROOT, v),
                          must_exist=True)
    fasta_abs = _pabs(fasta) if fasta else None
    ann_abs = _pabs(ann) if ann else None
    if ann_abs and str(ann_abs).lower().endswith(
            ('.gb', '.gbk', '.gbff', '.genbank')):
        fasta_abs = None
    # 过滤 opts（去空/False/None）。预览只出 SVG，排除仅保存时生效的格式/多记录。
    _PREVIEW_SKIP = {'format', 'multi_record_canvas'}
    opts = {str(k): v for k, v in gb_opts.items()
            if v not in (None, '', False) and k not in _PREVIEW_SKIP}
    tmp = tempfile.mkdtemp(prefix='vp_genoplot_preview_')
    try:
        # 用 ascii 临时目录（gbdraw 需 ascii 路径）
        import shutil as _sh
        tmp_ascii = os.path.join(tempfile.gettempdir(), 'vp_geno_preview')
        os.makedirs(tmp_ascii, exist_ok=True)
        from vp.gbdraw_plot import _run_gbdraw
        fasta_p, gff_p, gbk_p = None, None, None
        if ann_abs and str(ann_abs).lower().endswith(
                ('.gb', '.gbk', '.gbff', '.genbank')):
            gbk_p = ann_abs
        else:
            fasta_p = fasta_abs
            # 可选 GFF：若 ann 是 .gff 且 fasta 也给了 → 配对
            if ann_abs and str(ann_abs).lower().endswith(('.gff', '.gff3')):
                gff_p = ann_abs
        made = _run_gbdraw(fasta=fasta_p, gff=gff_p, gbk=gbk_p,
                           out_prefix=os.path.join(tmp_ascii, 'preview'),
                           mode=body.get('mode') or 'circular', opts=opts)
        if not made:
            abort(400, 'gbdraw 未产出预览图')
        svg_path = made[0]
        with safe_open(svg_path) as f:
            svg = f.read()
        return jsonify({'ok': True, 'mode': 'circular', 'svg': svg})
    except Exception as e:
        raise
    finally:
        _sh.rmtree(tmp, ignore_errors=True) if os.path.isdir(tmp) else None


def _tool_job_primer(ctx):
    """引物设计（primer3）。plain=基因组/contigs 全长分窗；
    conserved=多序列比对 FASTA 保守区（输入需已比对，如 MAFFT aln.fasta）。"""
    import shutil
    fasta = ctx.req('fasta', '输入 FASTA')
    mode = ctx.p.get('mode') or 'plain'
    if mode not in ('conserved', 'plain'):
        abort(400, '无效的引物设计模式')
    num_return = max(1, min(int(ctx.p.get('num_return') or 3), 20))
    specificity = bool(ctx.p.get('specificity'))

    def job(log, prog, cancel):
        import json as _json
        logger = TaskLogger(callback=log)
        if mode == 'conserved':
            p_dir = check_path(os.path.join(ctx.run_dir, '05_phylo'),
                               must_exist=False, in_platform=True)
            os.makedirs(os.path.join(p_dir, 'G1'), exist_ok=True)
            shutil.copyfile(fasta, os.path.join(p_dir, 'G1', 'aln.fasta'))
            with safe_open(os.path.join(p_dir, 'summary.json'), 'wt') as f:
                _json.dump({'groups': [{'group': 'G1', 'dir': 'G1'}]}, f)
        else:
            a_dir = check_path(os.path.join(ctx.run_dir, '03_assembly'),
                               must_exist=False, in_platform=True)
            os.makedirs(a_dir, exist_ok=True)
            shutil.copyfile(fasta, os.path.join(a_dir, 'viral_contigs.fasta'))
        from vp.primer import design_primers
        prog('primer', 0.1, f'primer3 引物设计（{mode}）')
        res = design_primers(ctx.run_dir, mode=mode, num_return=num_return,
                             logger=logger, force=True,
                             do_specificity=specificity)
        prog('done', 1.0, f"引物 {res.get('n_primers', 0)} 对")
        logger.close()
        return res
    return job


def _tool_job_identify(ctx):
    """② 病毒鉴定与提取（fastq 双端/单端 或 fasta contigs）。"""
    inp = ctx.req('input', '输入文件')
    itype = ctx.p.get('input_type') or 'pe'
    if itype not in ('pe', 'single', 'fasta'):
        abort(400, '无效的输入类型')
    conf = float(ctx.p.get('confidence') or 0)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.kunpeng import classify, parse_classify_output
        inputs = [inp]
        if itype == 'pe':
            inputs.append(ctx.req('input2', 'R2 FASTQ'))
        prog('classify', 0.1, 'kunpeng 病毒库分类中')
        res = classify(ctx.db_virus, inputs, os.path.join(ctx.run_dir, 'classify'),
                       paired=(itype == 'pe'), threads=ctx.threads,
                       confidence=conf, logger=logger,
                       progress=lambda pp, mm: prog(
                           'classify', 0.1 + pp * 0.75, mm))
        ids = {}
        if res['kraken']:
            for flag, rid, taxid, _l, mapping in parse_classify_output(res['kraken']):
                if flag != 'C':
                    continue
                # kraken2 口径分值：支持该 taxid 的片段数 / 映射列总片段数
                sup = tot = 0
                for seg in (mapping or '').split():
                    t, _, c = seg.rpartition(':')
                    try:
                        cnt = int(c)
                    except ValueError:
                        continue
                    tot += cnt
                    if t == str(taxid):
                        sup += cnt
                ids[rid] = (taxid, round(sup / tot, 3) if tot else 0)
        prog('extract', 0.9, '提取病毒候选序列')
        n_ext = {}
        if ids:
            for i, s_ in enumerate(inputs, 1):
                tag = '' if len(inputs) == 1 else f'_{i}'
                dst = os.path.join(ctx.run_dir, f'viral_sequences{tag}.fasta')
                n_ext[os.path.basename(dst)] = _extract_records(s_, set(ids), dst)
                logger.log(f'提取病毒序列 {os.path.basename(dst)}: '
                           f'{n_ext[os.path.basename(dst)]} 条')
            with safe_open(os.path.join(ctx.run_dir, 'viral_ids.tsv'), 'wt') as f:
                f.write('seq_id\ttaxid\tscore\n')
                for rid, (tx, sc) in ids.items():
                    f.write(f'{rid}\t{tx}\t{sc}\n')
        prog('extract', 1.0, '完成')
        logger.close()
        # n_extracted 必须是标量：结果预览的统计条只展示标量字段
        return {'n_classified': len(ids), 'n_extracted': sum(n_ext.values()),
                'kreport': res['kreport']}
    return job


def _tool_job_assemble(ctx):
    """③ 病毒组装（SPAdes，可选模式）。"""
    r1 = ctx.req('r1', 'R1 FASTQ')
    r2 = ctx.req('r2', 'R2 FASTQ')
    mode = ctx.p.get('mode') or 'metaviral'
    mem = int(ctx.p.get('memory') or 64)
    min_len = int(ctx.p.get('min_len') or 200)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.assembly import run_spades, filter_contigs
        out = os.path.join(ctx.run_dir, 'assembly')
        prog('spades', 0.05, 'SPAdes 组装中（耗时主要步骤）')
        spades_out = run_spades(r1, r2, out, mode=mode, threads=ctx.threads,
                                memory_gb=mem, logger=logger,
                                progress=lambda pp, mm: prog(
                                    'spades', 0.05 + pp * 0.85, mm))
        # 组装产物可能是 contigs.fasta（常规/metaviral-meta）或
        # transcripts.fasta（低覆盖自动降级 rna 模式），以 run_spades 实际
        # 返回的路径为准，不能硬编码 contigs.fasta（否则降级时找不到文件）。
        prog('filter', 0.92, 'contig 长度过滤')
        filtered, n_c, total_bp = filter_contigs(
            spades_out,
            os.path.join(ctx.run_dir, 'contigs.filtered.fasta'),
            min_len=min_len, logger=logger)
        prog('filter', 1.0, '完成')
        logger.close()
        return {'n_contigs': n_c, 'total_bp': total_bp,
                'contigs': filtered}
    return job


def _tool_job_contigs(ctx):
    """④ contig 病毒分类与提取（输入 contigs fasta）。"""
    contigs = ctx.req('contigs', 'contigs FASTA')
    min_len = int(ctx.p.get('min_len') or 200)
    conf = float(ctx.p.get('confidence') or 0)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.assembly import filter_contigs
        from vp.kunpeng import classify, parse_classify_output
        from vp.utils import run_cmd
        from vp.contig_annot import classify_rows, genus_avg_map, RANKS

        prog('genus_lens', 0.03, '统计属平均基因组长度（首跑需建缓存）')
        genus_map = genus_avg_map(logger=logger)

        prog('filter', 0.05, 'contig 长度过滤')
        filtered, n_c, _bp = filter_contigs(
            contigs, os.path.join(ctx.run_dir, 'contigs.filtered.fasta'),
            min_len=min_len, logger=logger)
        if n_c == 0:
            raise RuntimeError('过滤后无 contigs（检查最小长度设置）')

        prog('classify', 0.15, 'kunpeng 病毒库分类')
        res = classify(ctx.db_virus, [filtered],
                       os.path.join(ctx.run_dir, 'classify'),
                       paired=False, threads=ctx.threads, confidence=conf,
                       logger=logger,
                       progress=lambda pp, mm: prog(
                           'classify', 0.15 + pp * 0.45, mm))
        ids = {}
        if res['kraken']:
            for flag, _rid, _tx, _l, _pa in parse_classify_output(res['kraken']):
                if flag == 'C':
                    ids[_rid] = _tx

        # metabuli 风格分类表：8 级谱系 + 属平均长度 + 近完整判定
        prog('annot', 0.8, '谱系注释与属长比整理')
        rows = []
        if res['kraken']:
            rows = classify_rows(res['kraken'], genus_map)
        tsv_path = os.path.join(ctx.run_dir, 'virus_classification.tsv')
        header = (['contig', 'taxid', 'taxon'] + RANKS
                  + ['length', 'genus_avg_len', 'ratio', 'near_complete',
                     'score', 'kmer_support', 'kmer_total'])
        with safe_open(tsv_path, 'wt') as f:
            f.write('\t'.join(header) + '\n')
            for r in rows:
                f.write('\t'.join(str(r.get(k, '')) for k in header) + '\n')

        prog('extract', 0.88, '提取病毒 contigs（带谱系 header）')
        n_viral = 0
        viral_fa = None
        if ids:
            viral_fa = os.path.join(ctx.run_dir, 'viral_contigs.fasta')
            seqs = _load_fasta_seqs(filtered)
            with safe_open(viral_fa, 'wt') as f:
                for rid, taxid in ids.items():
                    seq = seqs.get(rid)
                    if seq is None:
                        continue
                    n_viral += 1
                    row = next((r for r in rows if r['contig'] == rid), None)
                    lineage = ';'.join(row[r] for r in RANKS
                                       if row and row.get(r))
                    taxon = (row or {}).get('taxon', '')
                    f.write(f'>{rid} taxid={taxid} taxon='
                            f'{taxon.replace(" ", "_")} '
                            f'lineage={lineage}\n')
                    for i in range(0, len(seq), 70):
                        f.write(seq[i:i + 70] + '\n')

        prog('extract', 1.0, '完成')
        logger.close()
        return {'n_contigs': n_c, 'n_viral': n_viral,
                'kreport': res['kreport'],
                'classification': tsv_path, 'viral_fasta': viral_fa}
    return job


def _tool_job_structcmp(ctx):
    """结构比较：多条序列 MAFFT 全长比对 → 两两 identity 矩阵（SDT 口径）。"""
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    max_n = max(3, min(int(ctx.p.get('max_n') or 30), 200))

    def _short(h, idx):
        name = re.split(r'[\s|]', (h or '').strip())[0][:40]
        name = re.sub(r'[^A-Za-z0-9_\-.]', '_', name)
        return name or f'seq{idx}'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.phylo import _run_mafft
        from vp.utils import iter_fasta

        prog('read', 0.05, '读取与筛选序列')
        recs, seen = [], {}
        for h, s in iter_fasta(seqs_fa):
            name = _short(h, len(recs) + 1)
            if name in seen:
                seen[name] += 1
                name = f'{name}_{seen[name]}'
            else:
                seen[name] = 0
            recs.append((name, s.upper()))
        if len(recs) < 2:
            raise RuntimeError('FASTA 中少于 2 条序列，无法做结构比较')
        recs = recs[:max_n]

        capped_fa = os.path.join(ctx.run_dir, 'input.fasta')
        with safe_open(capped_fa, 'wt') as f:
            for name, s in recs:
                f.write(f'>{name}\n')
                for i in range(0, len(s), 70):
                    f.write(s[i:i + 70] + '\n')

        prog('aln', 0.15, 'MAFFT 全长比对')
        aln = _run_mafft(capped_fa, os.path.join(ctx.run_dir, 'aln.fasta'),
                         threads=ctx.threads, logger=logger)

        prog('matrix', 0.75, '计算两两 identity 矩阵')
        names, cols = [], []
        for h, s in iter_fasta(aln):
            names.append(_short(h, len(names) + 1))
            cols.append(s.upper())
        n = len(names)
        matrix = [[100.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                same = comp = 0
                for a, b in zip(cols[i], cols[j]):
                    if a == '-' and b == '-':
                        continue
                    comp += 1
                    if a == b:
                        same += 1
                pid = round(same / comp * 100, 2) if comp else 0.0
                matrix[i][j] = matrix[j][i] = pid

        tsv = os.path.join(ctx.run_dir, 'identity_matrix.tsv')
        with safe_open(tsv, 'wt') as f:
            f.write('\t'.join([''] + names) + '\n')
            for i in range(n):
                f.write('\t'.join([names[i]] +
                                  [f'{matrix[i][j]:.2f}' for j in range(n)]) + '\n')
        data = {'names': names, 'matrix': matrix, 'n': n,
                'aln_cols': len(cols[0]) if cols else 0}
        js = os.path.join(ctx.run_dir, 'identity_matrix.json')
        with safe_open(js, 'wt') as f:
            json.dump(data, f, ensure_ascii=False)

        prog('done', 1.0, '完成')
        logger.close()
        return {'n_seqs': n, 'aln': aln, 'matrix_tsv': tsv,
                'matrix_json': js, 'aln_cols': data['aln_cols']}
    return job


def _tool_job_verify(ctx):
    """候选序列验证（对齐 02b/09b）：宿主筛选 → 长度分流 → 双路过滤。

    输入可选已有 contigs 运行（run）或独立 FASTA（fasta）。
    参数：host（默认 all）、methods（blastx/cdd 组合）、combine（union/intersection）。
    """
    run_ref = (ctx.p.get('run') or '').strip()
    fasta = ctx.opt('fasta') if ctx.p.get('fasta') else None
    host = (ctx.p.get('host') or 'all').strip() or 'all'
    combine = (ctx.p.get('combine') or 'union').strip()
    methods = ctx.p.get('methods') or ['blastx', 'cdd']
    if isinstance(methods, str):
        methods = [m.strip() for m in methods.split(',') if m.strip()]
    methods = [m for m in methods if m in ('blastx', 'cdd')]

    if run_ref and not re.fullmatch(r'[A-Za-z0-9_\-]+', run_ref):
        abort(400, '无效的运行名')
    if combine not in ('union', 'intersection'):
        combine = 'union'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.verify import verify

        # 输入解析：优先 run（其 viral_contigs.fasta），否则独立 fasta
        run_dir = ctx.run_dir
        if run_ref:
            src_run = check_path(os.path.join(_tool_runs_root(), run_ref),
                                 must_exist=True, in_platform=True)
            vfa = os.path.join(src_run, 'viral_contigs.fasta')
            if not os.path.isfile(vfa):
                raise RuntimeError('该运行无 viral_contigs.fasta（先跑 contig 分类）')
            # 宿主归属依赖源运行的分类表
            for _fn in ('virus_classification.tsv',):
                s = os.path.join(src_run, _fn)
                if os.path.isfile(s):
                    import shutil
                    shutil.copyfile(s, os.path.join(run_dir, _fn))
        else:
            if not fasta:
                abort(400, '请选择 contigs 运行或输入 FASTA')
            vfa = fasta

        summary = verify(run_dir, vfa, host=host, methods=methods,
                         combine=combine, threads=ctx.threads, logger=logger,
                         progress=lambda st, fr, msg: prog(st, fr, msg))
        logger.close()
        return summary
    return job


def _tool_job_consensus(ctx):
    """共识序列与变异分析（验证模块之后）：reads 回贴 → 共识序列 + 变异谱。

    参考三条来源：
      1. kvsuite 运行选参考（virus-fasta/ref_<acc>/）—— 已知病毒基因组
      2. contig 分类运行的 viral_contigs.fasta —— 未知/组装候选
      3. 独立 FASTA
    reads 默认取**去宿主后全量**（01_host_removal）——变异检测无偏；
    02_virus_screen 已按相似度丢过一轮 reads，会系统性低估变异与准种多样性。
    默认重新比对（而非复用 kvsuite 的 BAM）：旧 BAM 有比对偏好性，
    新 BAM 覆盖更全。
    """
    run_ref = (ctx.p.get('run') or '').strip()
    fasta = ctx.opt('fasta') if ctx.p.get('fasta') else None
    kv_run = (ctx.p.get('kv_run') or '').strip()
    kv_refs = [x.strip() for x in (ctx.p.get('kv_refs') or '').split(',') if x.strip()]
    reuse_bam = (ctx.p.get('reuse_bam') or '').strip()
    reads_src = (ctx.p.get('reads_src') or 'host_removed').strip()
    ambig = (ctx.p.get('ambig') or 'N').strip()[:1] or 'N'

    def _num(v, default, cast):
        try:
            return cast(v)
        except (TypeError, ValueError):
            return default

    min_qual = _num(ctx.p.get('min_qual'), 20, int)
    min_depth = _num(ctx.p.get('min_depth'), 10, int)
    min_freq = _num(ctx.p.get('min_freq'), 0.5, float)
    min_mapq = _num(ctx.p.get('min_mapq'), 10, int)
    min_minor_freq = _num(ctx.p.get('min_minor_freq'), 0.02, float)
    min_cov_pct = _num(ctx.p.get('min_cov_pct'), 10.0, float)
    min_qual = max(0, min(min_qual, 60))
    min_depth = max(1, min(min_depth, 100000))
    min_freq = min(1.0, max(0.0, min_freq))

    if run_ref and not re.fullmatch(r'[A-Za-z0-9_\-]+', run_ref):
        abort(400, '无效的运行名')
    if kv_run and not re.fullmatch(r'[A-Za-z0-9_\-]+', kv_run):
        abort(400, '无效的 kvsuite 运行名')
    if kv_refs and not kv_run:
        abort(400, '选了参考序列但未指定 kvsuite 运行')
    for _r in kv_refs:
        if not re.fullmatch(r'[A-Za-z0-9_.\-]+', _r):
            abort(400, '无效的 accession: %s' % _r)
    if reads_src not in CONSENSUS_READS_LAYERS:
        reads_src = 'host_removed'
    if reuse_bam:
        # 只允许平台运行目录内的 BAM（防任意路径写入）
        reuse_bam = check_path(reuse_bam, must_exist=False, in_platform=True)
        if not os.path.isfile(reuse_bam):
            abort(400, '复用的 BAM 不存在：%s' % reuse_bam)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.consensus import consensus_and_variants, find_read_pairs

        # ---- 参考 ----
        # 三条来源，优先级：kvsuite 选参考 > contig 运行 > 独立 FASTA
        #   kvsuite：<run>/kvsuite/virus-fasta/ref_<acc>/ref_<acc>.ref.fasta
        #            多选按 accession 合并成临时库（序列名已保证唯一）
        src_run = None
        vfa = None
        if kv_refs:
            kv_base = check_path(os.path.join(_tool_runs_root(), kv_run,
                                              'kvsuite'),
                                 must_exist=True, in_platform=True)
            fa_dir = os.path.join(kv_base, 'virus-fasta')
            parts, missing = [], []
            for acc in kv_refs:
                sub = os.path.join(fa_dir, 'ref_%s' % acc)
                cand = os.path.join(sub, 'ref_%s.ref.fasta' % acc)
                if os.path.isfile(cand) and os.path.getsize(cand) > 0:
                    parts.append((acc, cand))
                    continue
                # 容错：目录在但文件名不同
                found = None
                if os.path.isdir(sub):
                    for fn in sorted(os.listdir(sub)):
                        if fn.endswith(('.fasta', '.fa', '.fna')):
                            found = os.path.join(sub, fn)
                            break
                if found:
                    parts.append((acc, found))
                else:
                    missing.append(acc)
            if not parts:
                raise RuntimeError(
                    '所选参考在 %s 下均无 FASTA（先跑共识段生成 virus-fasta）'
                    % fa_dir)
            if missing:
                log('跳过无 FASTA 的参考: %s' % ', '.join(missing))
            # 写到运行目录，作为本次分析的临时参考
            tmp_ref = os.path.join(ctx.run_dir, 'ref_%s.fasta' % kv_run)
            n_seq = 0
            with open(tmp_ref, 'w', encoding='utf-8', newline='\n') as out:
                for acc, path in parts:
                    with open(path, encoding='utf-8', errors='replace') as fh:
                        for line in fh:
                            if line.startswith('>'):
                                n_seq += 1
                                head = line[1:].strip().split()[0]
                                # accession 唯一化：序列名已一致时保留原样
                                if head == acc:
                                    out.write(line if line.endswith('\n')
                                              else line + '\n')
                                else:
                                    out.write('>%s\n' % acc)
                            else:
                                out.write(line if line.endswith('\n')
                                          else line + '\n')
            log('kvsuite 参考库：%d 条（%d 个 accession）-> %s'
                % (n_seq, len(parts), os.path.basename(tmp_ref)))
            vfa = tmp_ref
        elif run_ref:
            src_run = check_path(os.path.join(_tool_runs_root(), run_ref),
                                 must_exist=True, in_platform=True)
            vfa = os.path.join(src_run, 'viral_contigs.fasta')
            if not os.path.isfile(vfa):
                raise RuntimeError('该运行无 viral_contigs.fasta（先跑 contig 分类）')
        else:
            if not fasta:
                abort(400, '请选择 contigs 运行、kvsuite 参考或输入参考 FASTA')
            vfa = fasta

        # ---- reads ----
        # 默认重新比对：用去宿主后全量 reads，让新 BAM 覆盖更全。
        # kvsuite 自己的 bam/virus_reads 已按病毒筛选过，复用会引入比对
        # 偏好性并系统性低估变异，因此不作为默认来源。
        reads = []
        if not reuse_bam:
            for p in (ctx.p.get('reads') or '').split(','):
                p = p.strip()
                if not p:
                    continue
                if os.path.isdir(p):
                    reads.extend(find_read_pairs(p))
                else:
                    reads.append(p)
            if not reads:
                if not src_run:
                    # kvsuite 参考或独立 FASTA：也需要 reads
                    if kv_run:
                        abort(400, '选择了 kvsuite 参考，请显式指定回贴 reads 文件/目录')
                    abort(400, '未指定运行时必须直接给出 reads 文件路径')
                layer = os.path.join(src_run, CONSENSUS_READS_LAYERS[reads_src])
                reads = find_read_pairs(layer)
                if not reads:
                    raise RuntimeError(
                        '在 %s 下未找到 reads（期望 *_R1/*.fastq.gz 配对文件）'
                        % CONSENSUS_READS_LAYERS[reads_src])
        else:
            log('映射模式：复用 BAM %s' % os.path.basename(reuse_bam))

        summary = consensus_and_variants(
            ctx.run_dir, vfa, reads, out_subdir='consensus',
            min_qual=min_qual, min_depth=min_depth, min_freq=min_freq,
            ambig=ambig, min_mapq=min_mapq, min_minor_freq=min_minor_freq,
            min_cov_pct=min_cov_pct, preset='sr', threads=ctx.threads,
            logger=logger, progress=lambda st, fr, msg: prog(st, fr, msg),
            reuse_bam=reuse_bam or None)
        logger.close()
        return summary
    return job


# 共识模块可选 reads 来源：层级目录名。
# 默认 host_removed（去宿主后全量）——未经过病毒相似度筛选，变异检测无偏。
CONSENSUS_READS_LAYERS = {
    'host_removed': '01_host_removal',
    'viral':        '02_virus_screen',
    'clean':        '00_prep',
}


def _tool_job_kvsuite(ctx):
    """已知病毒识别与定量（known_virus_suite 五段整合）。

    完全照搬 D:/桌面/延伸基因组/MMPV-RNA/virome_analysis_pipeline 的做法：
      鉴定 → 过滤 → 共识 → 深度绘图 → 变异注释
    引擎 minibwa 替代 bowtie2，**其余一点不改**。

    caller 用 bcftools mpileup+call（freebayes/lofreq/ivar 在 Windows 上均不可得），
    参数与阈值为实测定稿，详见 known_virus_suite/POSCOUNTS_REMOVAL_PLAN.md §4.2。
    """
    import json
    import subprocess

    # ── 输入：勾选样品 → 临时 sample-sheet TSV(name,r1,r2) ──
    # 变异段（variant）不吃 reads，只需 BAM 或 VCF，因此允许 samples 为空。
    stage = (ctx.p.get('stage') or 'all').strip().lower()
    if stage not in ('all', 'index', 'identify', 'filter', 'consensus',
                     'plot', 'variant'):
        stage = 'all'
    samples_raw = (ctx.p.get('samples') or '').strip()
    sheet_in = (ctx.p.get('sample_sheet') or '').strip()
    if sheet_in:
        sheet = ctx.req('sample_sheet', '样本表 TSV')
    elif stage == 'variant':
        sheet = None          # 变异段不用样本表
    else:
        if not samples_raw:
            abort(400, '请选择样品，或提供样本表 TSV')
        names = [s.strip() for s in samples_raw.split(',') if s.strip()]
        if not names:
            abort(400, '未选择有效样品')
        from vp.pipeline import load_sample_input
        rows = []
        for nm in names:
            sd = check_path(os.path.join(DIRS['results'], nm),
                            must_exist=True, in_platform=True)
            r1, r2, _proj = load_sample_input(sd)
            if not r1 or not os.path.isfile(r1):
                abort(400, f'样品 {nm} 的 R1 不可用: {r1}')
            rows.append((nm, r1, r2 or ''))
        sheet = os.path.join(ctx.run_dir, 'sample_sheet.tsv')
        with open(sheet, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write('name\tr1\tr2\n')
            for nm, r1, r2 in rows:
                fh.write(f'{nm}\t{r1}\t{r2}\n')

    # ── 参考序列与注释 ──
    ref = (ctx.p.get('reference') or '').strip()
    if ref:
        ref = ctx.req('reference', '参考序列')
    else:
        ref = os.path.join(DIRS['virus_src'], 'final.cluster.ref.fasta')
        ref = check_path(ref, must_exist=True, in_platform=True)
    ref_info = (ctx.p.get('ref_info') or '').strip()
    if ref_info:
        ref_info = ctx.req('ref_info', '参考注释')
    else:
        cand = os.path.join(DIRS['virus_src'], 'final.cluster.ref_info.tsv')
        ref_info = cand if os.path.isfile(cand) else None

    engine = (ctx.p.get('engine') or 'minibwa').strip().lower()
    if engine not in ('minibwa', 'salmon'):
        engine = 'minibwa'
    # 索引复用：① 用户在前端选的「鉴定库」优先（build 页建好的）；
    # ② 默认参考 + 预建索引在位时自动复用（避免每个 run 重建 60 MB 索引）；
    # ③ 都不满足则回退 <out>/index 临时建。
    index_dir = None
    user_idx = (ctx.p.get('index_dir') or '').strip()
    if user_idx:
        index_dir = ctx.req('index_dir', '鉴定库目录')
    elif engine == 'minibwa':
        kv_idx = os.path.join(DIRS['virus_src'], 'kv_index')
        default_ref = os.path.join(DIRS['virus_src'], 'final.cluster.ref.fasta')
        if (os.path.isfile(os.path.join(kv_idx, 'minibwa.mbw'))
                and os.path.abspath(ref) == os.path.abspath(default_ref)):
            index_dir = kv_idx
    elif engine == 'salmon':
        kv_idx = os.path.join(DIRS['virus_src'], 'kv_index')
        default_ref = os.path.join(DIRS['virus_src'], 'final.cluster.ref.fasta')
        if (os.path.isfile(os.path.join(kv_idx, 'salmon_k31', 'info.json'))
                and os.path.abspath(ref) == os.path.abspath(default_ref)):
            index_dir = kv_idx

    def _num(v, default, cast, lo=None, hi=None):
        try:
            x = cast(v)
        except (TypeError, ValueError):
            return default
        if lo is not None:
            x = max(lo, x)
        if hi is not None:
            x = min(hi, x)
        return x

    min_cov = _num(ctx.p.get('min_cov'), 10.0, float, 0.0, 100.0)
    min_depth = _num(ctx.p.get('min_depth'), 0.5, float, 0.0)
    min_reads = _num(ctx.p.get('min_reads'), 10, int, 0)
    min_poisson = _num(ctx.p.get('min_poisson'), 0.3, float, 0.0, 1.0)
    variant_qual = _num(ctx.p.get('variant_qual'), 3.5, float, 0.0)
    min_freq = _num(ctx.p.get('min_freq'), 0.05, float, 0.0, 1.0)
    aa_label_cutoff = _num(ctx.p.get('aa_label_cutoff'), 0.50, float, 0.0, 1.01)
    max_aa_labels = _num(ctx.p.get('max_aa_labels'), 40, int, 0)

    def job(log, prog, cancel):
        suite = os.path.join(PLATFORM_ROOT, 'known_virus_suite',
                             'known_virus_suite.py')
        if not os.path.isfile(suite):
            raise RuntimeError(f'未找到 known_virus_suite.py: {suite}')

        out_dir = os.path.join(ctx.run_dir, 'kvsuite')
        os.makedirs(out_dir, exist_ok=True)
        cmd = [sys.executable, suite, stage,
               '--out', out_dir,
               '--reference', ref,
               '--engine', engine,
               '--min-cov', str(min_cov),
               '--min-depth', str(min_depth),
               '--min-reads', str(min_reads),
               '--min-poisson', str(min_poisson),
               '--variant-qual', str(variant_qual),
               '--min-freq', str(min_freq),
               '--aa-label-cutoff', str(aa_label_cutoff),
               '--max-aa-labels', str(max_aa_labels)]
        if sheet:
            cmd += ['--sample-sheet', sheet]
        # 变异段输入：外部 VCF 优先于 BAM 目录（CLI 内也是这个优先级）
        input_vcf = (ctx.p.get('input_vcf') or '').strip()
        if input_vcf:
            cmd += ['--input-vcf', ctx.req('input_vcf', '变异 VCF 文件')]
        bam_dir = (ctx.p.get('bam_dir') or '').strip()
        if bam_dir:
            cmd += ['--bam-dir', ctx.req('bam_dir', 'BAM 目录')]
        if ref_info:
            cmd += ['--ref-info', ref_info]
        if index_dir:
            cmd += ['--index-dir', index_dir]
        if ctx.threads:
            cmd += ['--threads', str(ctx.threads),
                    '--align-threads', str(ctx.threads)]
        if ctx.p.get('ncbi_email'):
            cmd += ['--ncbi-email', str(ctx.p['ncbi_email']).strip()]
        if ctx.p.get('ncbi_api_key'):
            cmd += ['--ncbi-api-key', str(ctx.p['ncbi_api_key']).strip()]
        if ctx.p.get('no_genes'):
            cmd.append('--no-genes')
        if ctx.p.get('no_variant_evo'):
            cmd.append('--no-variant-evo')
        if ctx.p.get('all_variants'):
            cmd.append('--all-variants')
        if ctx.p.get('no_snpgenie'):
            cmd.append('--no-snpgenie')

        log('$ ' + ' '.join(cmd))
        env = dict(os.environ)
        env.setdefault('PYTHONIOENCODING', 'utf-8')
        proc = subprocess.Popen(cmd, cwd=PLATFORM_ROOT,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                text=True, encoding='utf-8',
                                errors='replace', bufsize=1, env=env)
        tail = []
        for line in proc.stdout:
            line = line.rstrip('\n')
            if not line:
                continue
            log(line)
            tail.append(line)
            if len(tail) > 40:
                tail.pop(0)
            # cancel 是 threading.Event（见 TaskServer._run），必须用 is_set()，
            # 直接 cancel() 会报 'Event' object is not callable。
            if cancel is not None and cancel.is_set():
                proc.kill()
                raise RuntimeError('用户取消')
        rc = proc.wait()
        if rc != 0:
            raise RuntimeError(f'known_virus_suite 退出码 {rc}；末几行：\n'
                               + '\n'.join(tail[-8:]))

        # ── 汇总产物（供结果面板展示）──
        result = {'stage': stage, 'engine': engine, 'out': out_dir,
                  'reference': ref}
        fsum = os.path.join(out_dir, 'filter', 'filtered.tsv')
        if os.path.isfile(fsum):
            result['filtered_tsv'] = fsum
            try:
                with open(fsum, encoding='utf-8', errors='replace') as fh:
                    result['n_confirmed'] = max(0, sum(1 for _ in fh) - 1)
            except OSError:
                pass
        vsum = os.path.join(out_dir, 'variant_summary.json')
        if os.path.isfile(vsum):
            result['variant_summary'] = vsum
            try:
                with open(vsum, encoding='utf-8') as fh:
                    data = json.load(fh)
                rows = data.get('results') if isinstance(data, dict) else data
                rows = rows or []
                result['n_variant_genomes'] = len(rows)
                result['n_variants'] = sum(int(d.get('n_variants') or 0)
                                           for d in rows)
                if isinstance(data, dict):
                    result['n_variant_skipped'] = int(data.get('n_skipped') or 0)
            except (OSError, ValueError, TypeError):
                pass
        return result
    return job


def _tool_job_quicktree(ctx):
    """快速建树：多条序列 MAFFT 全长比对 → NJ（纯 Python）/ FastTree。"""
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    method = (ctx.p.get('method') or 'nj').strip().lower()
    if method not in ('nj', 'fasttree'):
        abort(400, '建树方法仅支持 nj / fasttree')
    max_n = max(3, min(int(ctx.p.get('max_n') or 100), 500))

    def _short(h, idx):
        name = re.split(r'[\s|]', (h or '').strip())[0][:40]
        name = re.sub(r'[^A-Za-z0-9_\-.]', '_', name)
        return name or f'seq{idx}'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.phylo import _run_mafft, _run_nj, _run_fasttree
        from vp.utils import iter_fasta

        prog('read', 0.05, '读取与筛选序列')
        recs, seen = [], {}
        for h, s in iter_fasta(seqs_fa):
            name = _short(h, len(recs) + 1)
            if name in seen:
                seen[name] += 1
                name = f'{name}_{seen[name]}'
            else:
                seen[name] = 0
            recs.append((name, s.upper()))
        if len(recs) < 2:
            raise RuntimeError('FASTA 中少于 2 条序列，无法建树')
        recs = recs[:max_n]

        capped_fa = os.path.join(ctx.run_dir, 'input.fasta')
        with safe_open(capped_fa, 'wt') as f:
            for name, s in recs:
                f.write(f'>{name}\n')
                for i in range(0, len(s), 70):
                    f.write(s[i:i + 70] + '\n')

        prog('aln', 0.2, 'MAFFT 全长比对')
        aln = _run_mafft(capped_fa, os.path.join(ctx.run_dir, 'aln.fasta'),
                         threads=ctx.threads, logger=logger)

        prog('tree', 0.75, 'NJ 建树' if method == 'nj' else 'FastTree 建树')
        tree_name = 'nj.nwk' if method == 'nj' else 'tree.nwk'
        if method == 'nj':
            _run_nj(aln, os.path.join(ctx.run_dir, tree_name), logger=logger)
        else:
            _run_fasttree(aln, os.path.join(ctx.run_dir, tree_name),
                          logger=logger)

        prog('done', 1.0, '完成')
        logger.close()
        return {'n_seqs': len(recs), 'method': method, 'tree': tree_name}
    return job


def _tool_job_align(ctx):
    """序列比对（独立模块）：MAFFT（auto/L-INS-i/fast）→ trimAl 清剪。

    输入 FASTA（核酸或蛋白，自动判别；≥2 条）。产物：input.fasta /
    aln.fasta / aln.trim.fasta（trimAl 关闭时无）/ summary.json，
    可在「比对查看器」彩色浏览与编辑。
    """
    from vp.utils import iter_fasta, write_fasta_record
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    strategy = (ctx.p.get('strategy') or 'auto').strip().lower()
    if strategy not in ('auto', 'linsi', 'fast'):
        abort(400, '比对策略仅支持 auto / linsi / fast')
    trimal = (ctx.p.get('trimal') or 'automated1').strip().lower()
    if trimal not in ('automated1', 'gappyout', 'strict', 'none'):
        abort(400, 'trimAl 方法仅支持 automated1 / gappyout / strict / none')
    max_n = max(3, min(int(ctx.p.get('max_n') or 200), 500))

    def job(log, prog, cancel):
        import shutil
        logger = TaskLogger(callback=log)
        from vp.phylo import _run_mafft, _run_trimal
        from vp.sdt_exact import detect_seqtype
        prog('read', 0.05, '读取与筛选序列')
        recs, seen = [], {}
        for h, s in iter_fasta(seqs_fa):
            name = re.split(r'[\s|]', (h or '').strip())[0][:60] or \
                f'seq{len(recs) + 1}'
            name = re.sub(r'[^A-Za-z0-9_\-.]', '_', name)
            if name in seen:
                seen[name] += 1
                name = f'{name}_{seen[name]}'
            else:
                seen[name] = 0
            recs.append((name, s.upper()))
        if len(recs) < 2:
            raise RuntimeError('FASTA 中少于 2 条序列，无法比对')
        recs = recs[:max_n]
        seqtype = detect_seqtype([s for _, s in recs])
        capped = os.path.join(ctx.run_dir, 'input.fasta')
        with safe_open(capped, 'wt') as f:
            for name, s in recs:
                write_fasta_record(f, name, s)
        logger.log(f'序列类型判定: {"蛋白(aa)" if seqtype == "aa" else "核酸(nt)"}'
                   f'，{len(recs)} 条参与比对')
        prog('align', 0.15, f'MAFFT 比对（{strategy}）')
        aln = _run_mafft(capped, os.path.join(ctx.run_dir, 'aln.fasta'),
                         threads=ctx.threads, logger=logger, strategy=strategy)
        trim_info = {'applied': False}
        aln_used = aln
        if trimal != 'none':
            prog('trim', 0.7, f'trimAl 清剪（{trimal}）')
            aln_used, trim_info = _run_trimal(
                aln, os.path.join(ctx.run_dir, 'aln.trim.fasta'),
                logger=logger)
        prog('done', 1.0, '完成')
        logger.close()
        return {'n_seqs': len(recs), 'seqtype': seqtype, 'strategy': strategy,
                'trimal': trimal, 'trim': trim_info,
                'aln': 'aln.fasta',
                'aln_used': os.path.basename(aln_used)}
    return job


def _tool_job_sdt(ctx):
    """SDT 精确分析（纯 Python 复刻 SDTv1.3，替代外部 SDT exe）：
    逐对 MAFFT 独立比对 → Get_Similarity 公式 → 聚类排序热图 + 分布图。
    aligned=True 时输入为已比对 MSA，跳过比对器直接按公式计算；
    seqtype='auto' 自动判别核酸/蛋白（MMPV 双轨口径），蛋白按 AA 同一性。"""
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    max_n = max(3, min(int(ctx.p.get('max_n') or 30), 200))
    orient = bool(ctx.p.get('orient', True))
    seqtype = (ctx.p.get('seqtype') or 'auto').strip().lower()
    if seqtype not in ('auto', 'nt', 'aa'):
        seqtype = 'auto'
    palette = (ctx.p.get('palette') or 'sdt').strip().lower()
    if palette not in ('sdt', 'cividis', 'viridis', 'RdYlBu', 'Spectral',
                       'YlGnBu', 'coolwarm', 'magma'):
        palette = 'sdt'
    aligned = bool(ctx.p.get('aligned'))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        mafft = cfg.tool('mafft')
        from vp.sdt_exact import run_sdt_exact
        prog('sdt', 0.02, 'SDT 精确分析（MAFFT 逐对独立比对，SDT v1.3 口径）')
        res = run_sdt_exact(seqs_fa, ctx.run_dir, mafft, max_n=max_n,
                            orient=orient, threads=ctx.threads,
                            seqtype=seqtype,
                            palette=palette, logger=logger,
                            progress=lambda p, m: prog('sdt', p, m),
                            cancel=cancel, aligned=aligned)
        logger.close()
        return res
    return job


def _tool_job_identity(ctx):
    """核苷酸+氨基酸同一性表（BioAider Sequence Identity Matrix 口径）：
    NT 矩阵 + AA 矩阵（可选输入或最长 ORF 翻译）+ 复合热图（NT 上 /
    AA 下）+ 逐对同一性长表。"""
    nt_fa = ctx.req('nt_seqs', '核苷酸 FASTA')
    aa_fa = ctx.opt('aa_seqs')
    max_n = max(3, min(int(ctx.p.get('max_n') or 30), 200))
    palette = (ctx.p.get('palette') or 'sdt').strip().lower()
    if palette not in ('sdt', 'cividis', 'viridis', 'RdYlBu', 'Spectral',
                       'YlGnBu', 'coolwarm', 'magma'):
        palette = 'sdt'
    aligned = bool(ctx.p.get('aligned'))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        mafft = cfg.tool('mafft')
        from vp.sdt_exact import run_identity_table
        prog('idty', 0.02, '核苷酸+氨基酸同一性表（BioAider 口径）')
        res = run_identity_table(nt_fa, ctx.run_dir, mafft, aa_fasta=aa_fa,
                                 max_n=max_n, aligned=aligned,
                                 threads=ctx.threads, palette=palette,
                                 logger=logger,
                                 progress=lambda p, m: prog('idty', p, m),
                                 cancel=cancel)
        logger.close()
        return res
    return job


# ------------------------------------------------------------------
# 链式一键分析（chain jobs）
# 设计：每步仍调用现有 _tool_job_*，各建自己的 run_dir，一步不改；
# 汇总目录只放「符号链接 + 清单」，把一次多步运行的结果收拢成一份。
# 汇总目录：(tool_runs_root)/<chain>_<ts>/
#   _chain.json            每步的 run 名 / 产物 / 状态
#   stepN_<name>           指向该步 run_dir 的符号链接（目录）
#   关键产物文件链接        指向该步产物文件，便于直接在汇总目录取用
# ------------------------------------------------------------------

_CHAIN_LINKABLE = ('.fasta', '.fa', '.fna', '.tsv', '.json', '.txt',
                   '.vcf', '.gz', '.html', '.pdf', '.png')


def _chain_link(dst, src, logger=None, target_is_dir=None):
    """在汇总目录建链接指向 src。

    Windows 下 os.symlink 需开发者模式/管理员；不可用时按序降级：
    目录 → 硬链接不可用则跳过并记录；文件 → 硬链接；再不行则跳过。
    返回实际使用的方式（'symlink' / 'hardlink' / None），失败不抛异常。
    """
    if not os.path.exists(src):
        return None
    if target_is_dir is None:
        target_is_dir = os.path.isdir(src)
    try:
        os.symlink(src, dst, target_is_directory=bool(target_is_dir))
        return 'symlink'
    except (OSError, NotImplementedError, AttributeError):
        pass
    if not target_is_dir:
        try:
            os.link(src, dst)
            return 'hardlink'
        except OSError:
            pass
    if logger:
        logger.log(f'链接创建跳过（源或权限不可用）: {os.path.basename(dst)}')
    return None


def _chain_write_manifest(chain_dir, meta):
    p = os.path.join(chain_dir, '_chain.json')
    with safe_open(p, 'wt') as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    return p


def _chain_run_step(chain_dir, root, step_key, step_name, job_factory,
                    p, threads, db_virus, logger, prog, prog_lo, prog_hi,
                    cancel=None, collect=None, optional=()):
    """跑链中一步：新建 run_dir、构造 ctx、调用 job，返回 (run_name, result)。

    prog 为阶段名；进度在 [prog_lo, prog_hi] 区间内线性映射。
    collect: 该步要链接进汇总目录的产物文件名列表（可选）。
    optional: 这些 key 就算下游 job 用 ctx.req 取也允许为空（返回 None），
              用于兼容现有 job 把可选输入写成必填的情况。
    """
    from types import SimpleNamespace
    ts = time.strftime('%Y%m%d_%H%M%S')
    run_name = f'{step_key}_{ts}'
    run_dir = check_path(os.path.join(root, run_name), must_exist=False,
                         in_platform=True)
    os.makedirs(run_dir, exist_ok=True)

    def _input_abs(v):
        return check_path(v if os.path.isabs(v)
                          else os.path.join(PLATFORM_ROOT, v),
                          must_exist=True)

    def _req(key, what):
        v = (p.get(key) or '').strip()
        if not v:
            if key in optional:
                return None
            abort(400, f'请选择{what}')
        return _input_abs(v)

    def _opt(key):
        v = (p.get(key) or '').strip()
        return _input_abs(v) if v else None

    ctx = SimpleNamespace(p=p, run_dir=run_dir, threads=threads,
                          db_virus=db_virus, req=_req, opt=_opt)
    job = job_factory(ctx)

    def _sub_prog(_stage, frac, msg):
        f = prog_lo + max(0.0, min(1.0, float(frac))) * (prog_hi - prog_lo)
        prog(step_key, f, msg)

    prog(step_key, prog_lo, f'{step_name}：开始')
    res = job(lambda m: logger.log(m), _sub_prog, cancel)
    logger.log(f'{step_name} 完成: run={run_name}')

    # 目录链接（整步产物）
    _chain_link(os.path.join(chain_dir, f'step_{step_key}'), run_dir,
                logger=logger, target_is_dir=True)
    # 关键文件链接
    linked = []
    for fn in (collect or ()):
        src = os.path.join(run_dir, fn)
        if os.path.isfile(src):
            how = _chain_link(os.path.join(chain_dir, f'{step_key}__{fn}'),
                              src, logger=logger, target_is_dir=False)
            if how:
                linked.append(fn)
    return run_name, res, linked


def _tool_job_virchain(ctx):
    """病毒识别和分类分析 · 一键：②提取序列组装 → ③组装结果再鉴定 → ④候选序列验证。

    输入固定为去宿主序列（FASTA，可单/双端）。①（识别分类与提取）为可选：
      跳过 ①：去宿主序列 ────────────────────→ ②组装 → ③再鉴定 → ④验证
      运行 ①：去宿主序列 → ①提取病毒序列 → ②组装 → ③再鉴定 → ④验证
    即两条路都到 ②，区别只在 ② 的输入（原始去宿主序列 vs ① 提取的病毒序列）。
    组装（SPAdes）接受 FASTA，或自动 --only-assembler 跳过纠错。
    每步各建 run_dir，汇总目录放符号链接 + _chain.json。
    """
    root = _tool_runs_root()
    inp = ctx.req('input', '去宿主序列（FASTA）')
    inp2 = ctx.opt('input2')
    # 可选：先跑 ① 提取病毒序列，作为 ② 的输入
    run_identify = bool(ctx.p.get('run_identify'))
    # ① 的输入类型：双端 / 单端 / FASTA
    if inp2:
        id_type = 'pe'
    else:
        id_type = ('fasta' if inp.lower().endswith(('.fa', '.fasta', '.fna'))
                   else 'single')

    do_verify = ctx.p.get('do_verify')
    do_verify = True if do_verify is None else bool(do_verify)
    asm_mode = (ctx.p.get('mode') or 'metaviral').strip()
    mem = int(ctx.p.get('memory') or 64)
    min_len = int(ctx.p.get('min_len') or 200)
    conf = float(ctx.p.get('confidence') or 0)
    host = (ctx.p.get('host') or 'all').strip() or 'all'
    combine = (ctx.p.get('combine') or 'union').strip()
    methods = ctx.p.get('methods') or ['blastx', 'cdd']
    if isinstance(methods, str):
        methods = [m.strip() for m in methods.split(',') if m.strip()]
    methods = [m for m in methods if m in ('blastx', 'cdd')]

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        chain_dir = ctx.run_dir
        chain_name = os.path.basename(chain_dir)
        meta = {'chain': 'virchain', 'run': chain_name,
                'label': '病毒识别和分类分析 一键流程',
                'run_identify': run_identify,
                'created': time.strftime('%Y-%m-%d %H:%M:%S'), 'steps': []}

        # 进度权重：identify(若开) 0.20 / assembly 0.35 / contigs 0.25 / verify 0.20
        w_id = 0.20 if run_identify else 0.0
        w_asm = 0.35
        w_ctg = 0.25
        w_vfy = 0.20 if do_verify else 0.0
        tot = w_id + w_asm + w_ctg + w_vfy

        def _band(lo, w):
            return (lo, lo + (w / tot if tot else 0.0))

        cur = 0.0

        # ---- ①（可选）病毒识别分类与提取：产物作为 ② 的输入 ----
        asm_in1, asm_in2 = inp, inp2
        if run_identify:
            lo, hi = _band(cur, w_id)
            p1 = {'input': inp, 'input_type': id_type, 'confidence': conf}
            if inp2:
                p1['input2'] = inp2
            run1, res1, _ = _chain_run_step(
                chain_dir, root, 'identify', '① 病毒识别分类与提取',
                _tool_job_identify, p1, ctx.threads, ctx.db_virus,
                logger, prog, lo, hi, cancel=cancel,
                collect=('viral_ids.tsv',))
            meta['steps'].append({'key': 'identify', 'run': run1,
                                  'title': '① 病毒识别分类与提取',
                                  'summary': {k: v for k, v in (res1 or {}).items()
                                              if not isinstance(v, (list, dict))}})
            # 取 ① 实际产出的病毒序列作为 ② 的输入（配对输入会带 _1/_2 后缀）
            id_dir = os.path.join(root, run1)
            vseqs = []
            if os.path.isdir(id_dir):
                for fn in sorted(os.listdir(id_dir)):
                    if fn.startswith('viral_sequences') and fn.endswith('.fasta'):
                        fp = os.path.join(id_dir, fn)
                        if os.path.getsize(fp) > 0:
                            vseqs.append(fp)
            if not vseqs:
                raise RuntimeError('① 未提取到病毒序列（viral_sequences*.fasta 为空），'
                                   '可取消「先跑 ①」开关直接用去宿主序列组装')
            # 链接进汇总目录
            for fp in vseqs:
                _chain_link(os.path.join(chain_dir,
                                         'identify__' + os.path.basename(fp)),
                            fp, logger=logger, target_is_dir=False)
            asm_in1 = vseqs[0]
            asm_in2 = vseqs[1] if len(vseqs) > 1 else None
            meta['steps'][-1]['output'] = [os.path.basename(x) for x in vseqs]
            logger.log('② 组装输入改用 ① 提取的病毒序列: '
                       + ', '.join(os.path.basename(x) for x in vseqs))
            cur = hi

        # ---- ② 提取序列组装（输入 = 去宿主序列 或 ① 的提取序列） ----
        lo, hi = _band(cur, w_asm)
        p2 = {'r1': asm_in1, 'mode': asm_mode, 'memory': mem,
              'min_len': min_len}
        if asm_in2:
            p2['r2'] = asm_in2
        run2, res2, _ = _chain_run_step(
            chain_dir, root, 'assembly', '② 提取序列组装',
            _tool_job_assemble, p2, ctx.threads, ctx.db_virus,
            logger, prog, lo, hi, cancel=cancel,
            collect=('contigs.filtered.fasta',), optional=('r2',))
        contigs = (res2 or {}).get('contigs')
        meta['steps'].append({'key': 'assembly', 'run': run2,
                              'title': '② 提取序列组装',
                              'summary': {k: v for k, v in (res2 or {}).items()
                                          if not isinstance(v, (list, dict))}})
        if not contigs or not os.path.isfile(contigs):
            raise RuntimeError('② 组装未产出 contigs，后续步骤无法继续')
        cur = hi

        # ---- ③ 组装结果再鉴定 ----
        lo, hi = _band(cur, w_ctg)
        p3 = {'contigs': contigs, 'min_len': min_len, 'confidence': conf}
        run3, res3, _ = _chain_run_step(
            chain_dir, root, 'contigs', '③ 组装结果再鉴定',
            _tool_job_contigs, p3, ctx.threads, ctx.db_virus,
            logger, prog, lo, hi, cancel=cancel,
            collect=('virus_classification.tsv', 'viral_contigs.fasta'))
        viral_fa = (res3 or {}).get('viral_fasta')
        meta['steps'].append({'key': 'contigs', 'run': run3,
                              'title': '③ 组装结果再鉴定',
                              'summary': {k: v for k, v in (res3 or {}).items()
                                          if not isinstance(v, (list, dict))}})
        cur = hi

        # ---- ④ 候选序列验证（依赖 ③ 的产物，直接读其 run_dir） ----
        if do_verify:
            if not viral_fa or not os.path.isfile(viral_fa):
                logger.log('③ 未提取到病毒 contig，跳过 ④ 候选序列验证')
                meta['steps'].append({'key': 'verify', 'run': None,
                                      'title': '④ 候选序列验证',
                                      'skipped': '③ 未产出 viral_contigs.fasta'})
            else:
                lo, hi = _band(cur, w_vfy)
                # verify 的 run 模式需要源 run 名；直接传 ③ 的 run
                p4 = {'run': run3, 'host': host, 'combine': combine,
                      'methods': methods}
                run4, res4, _ = _chain_run_step(
                    chain_dir, root, 'verify', '④ 候选序列验证',
                    _tool_job_verify, p4, ctx.threads, ctx.db_virus,
                    logger, prog, lo, hi, cancel=cancel,
                    collect=('verify_summary.json',))
                meta['steps'].append({'key': 'verify', 'run': run4,
                                      'title': '④ 候选序列验证',
                                      'summary': {k: v for k, v in (res4 or {}).items()
                                                  if not isinstance(v, (list, dict))}})
                cur = hi

        prog('chain', 1.0, '一键流程完成')
        logger.close()
        meta['done'] = True
        _chain_write_manifest(chain_dir, meta)
        return {'chain': chain_name, 'run': chain_name,
                'n_steps': len(meta['steps']), 'steps': meta['steps'],
                'chain_json': os.path.join(chain_dir, '_chain.json')}
    return job


def _tool_job_kvchain(ctx):
    """病毒定量与共识 · 一键：调用 kvsuite 全流程（鉴定→过滤→共识→绘图→变异）。

    kvsuite 卡里的 stage 下拉默认已是 all，本 job 把「一次跑完并收拢产物」
    做成显式入口：建汇总目录 + 链接 kvsuite run_dir + 关键产物。
    """
    root = _tool_runs_root()

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        chain_dir = ctx.run_dir
        chain_name = os.path.basename(chain_dir)
        # kvsuite 的 stage 固定为 all（一键语义）
        p = dict(ctx.p)
        p['stage'] = 'all'
        meta = {'chain': 'kvchain', 'run': chain_name,
                'label': '病毒定量与共识 一键流程',
                'created': time.strftime('%Y-%m-%d %H:%M:%S'), 'steps': []}

        run1, res1, _ = _chain_run_step(
            chain_dir, root, 'kvsuite', '已知病毒识别与定量（全流程）',
            _tool_job_kvsuite, p, ctx.threads, ctx.db_virus,
            logger, prog, 0.0, 1.0, cancel=cancel,
            collect=('run.log',))
        meta['steps'].append({'key': 'kvsuite', 'run': run1,
                              'title': '已知病毒识别与定量（全流程）',
                              'summary': {k: v for k, v in (res1 or {}).items()
                                          if not isinstance(v, (list, dict))}})
        prog('chain', 1.0, '一键流程完成')
        logger.close()
        meta['done'] = True
        _chain_write_manifest(chain_dir, meta)
        return {'chain': chain_name, 'run': chain_name,
                'n_steps': 1, 'steps': meta['steps'],
                'kvsuite_run': run1,
                'chain_json': os.path.join(chain_dir, '_chain.json')}
    return job


def _tool_job_chainoverview(ctx):
    """占位：链式汇总目录的读取走 /api/tool/chain_result，不注册为工具。"""
    raise NotImplementedError


TOOL_REGISTRY = {
    'convert':  {'title': '格式转换', 'title_en': 'Format convert',
                 'need_db': False, 'job': _tool_job_convert},
    'fastp':    {'title': '质控预处理', 'title_en': 'QC preprocess',
                 'need_db': False, 'job': _tool_job_fastp},
    'hostremoval': {'title': '宿主去除与序列提取', 'title_en': 'Host removal',
                    'need_db': False, 'job': _tool_job_hostremoval},
    'hostpredict': {'title': '宿主预测', 'title_en': 'Host prediction',
                    'need_db': False, 'job': _tool_job_hostpredict},
    'orf':         {'title': 'ORF 预测', 'title_en': 'ORF predict',
                    'need_db': False, 'job': _tool_job_orf},
    'orfa':        {'title': '功能注释', 'title_en': 'ORF annotate',
                    'need_db': False, 'job': _tool_job_orfa},
    'genoplot':    {'title': '基因组图谱', 'title_en': 'Genome plots',
                    'need_db': False, 'job': _tool_job_genoplot},
    'primer':      {'title': '引物设计', 'title_en': 'Primer design',
                    'need_db': False, 'job': _tool_job_primer},
    'identify': {'title': '病毒鉴定', 'title_en': 'Virus identify',
                 'need_db': True, 'job': _tool_job_identify},
    'assemble': {'title': '病毒组装', 'title_en': 'Virus assembly',
                 'need_db': False, 'job': _tool_job_assemble},
    'contigs':  {'title': 'contig分类', 'title_en': 'Contig classify',
                 'need_db': True, 'job': _tool_job_contigs},
    'sdt':      {'title': 'SDT 精确分析', 'title_en': 'SDT exact',
                 'need_db': False, 'job': _tool_job_sdt},
    'identity': {'title': 'NT+AA 同一性表', 'title_en': 'NT+AA identity',
                 'need_db': False, 'job': _tool_job_identity},
    'structcmp': {'title': '结构比较', 'title_en': 'Structure compare',
                  'need_db': False, 'job': _tool_job_structcmp},
    'verify':   {'title': '候选序列验证', 'title_en': 'Candidate verify',
                 'need_db': False, 'job': _tool_job_verify},
    'consensus': {'title': '共识序列与变异', 'title_en': 'Consensus & variants',
                  'need_db': False, 'job': _tool_job_consensus},
    'kvsuite':  {'title': '已知病毒识别与定量', 'title_en': 'Known virus suite',
                 'need_db': False, 'job': _tool_job_kvsuite},
    'align':    {'title': '序列比对', 'title_en': 'Alignment',
                 'need_db': False, 'job': _tool_job_align},
    'quicktree': {'title': '快速建树', 'title_en': 'Quick tree',
                  'need_db': False, 'job': _tool_job_quicktree},
    'virchain': {'title': '病毒识别分类·一键', 'title_en': 'Virus chain',
                 'need_db': True, 'job': _tool_job_virchain},
    'kvchain':  {'title': '病毒定量与共识·一键', 'title_en': 'KV chain',
                 'need_db': False, 'job': _tool_job_kvchain},
}

# 轻量工具（走 light 闸门，可与 heavy 任务并行）：纯 IO 或秒级计算，
# 不占满多核/大内存。其余一律 heavy（kunpeng 分类 / SPAdes / DIAMOND /
# 建树 / SDT 等吃满线程的任务，同时最多跑 max_heavy_tasks 个）。
LIGHT_TOOLS = {'convert', 'genoplot', 'primer'}

# 工具②/④可选的内置 kunpeng 病毒库（键名 → 目录）。
# 主库 virus 由「数据库构建」页构建；refvirus/rvdb/k2viral 由
# vp/universal_ref 或外部预构建（/api/dbs 探测就绪状态）。
VIRUS_DB_PRESETS = {'virus': ('virus', 'plant'),
                    'refvirus': ('virus', 'ref'),
                    'rvdb': ('virus', 'rvdb'),
                    'k2viral': None}


def _resolve_virus_db(v=None):
    """工具②/④的病毒库选择：内置键名、平台相对路径或本机绝对路径。

    缺省用平台主病毒库（cfg.databases['virus']）。返回库目录绝对路径；
    键名不识别且路径不存在时 abort(400)。
    """
    if not v:
        return cfg.databases['virus']
    v = str(v).strip()
    if v in VIRUS_DB_PRESETS:
        mapping = VIRUS_DB_PRESETS[v]
        if mapping:
            return db_path(*mapping)
        return os.path.join(DIRS['databases'], 'k2viral_db')
    d = v if os.path.isabs(v) else os.path.join(PLATFORM_ROOT, v)
    try:
        return check_path(d, must_exist=True)
    except (ValueError, FileNotFoundError):
        abort(400, f'病毒库不存在: {v}（可用内置库: '
                   f'{", ".join(VIRUS_DB_PRESETS)} 或填库目录路径）')


@app.route('/api/tool/run', methods=['POST'])
def api_tool_run():
    """运行独立工具。body: {tool, threads, params{...}}；返回 {task, run}。"""
    from types import SimpleNamespace
    body = request.get_json(force=True) or {}
    tool = body.get('tool')
    threads = body.get('threads') or None
    p = body.get('params') or {}

    spec = TOOL_REGISTRY.get(tool)
    if not spec:
        abort(400, '未知工具')

    def _input_abs(v):
        """输入数据文件：相对路径按平台根解析；绝对路径允许电脑任意位置
        （只读）。输出产物仍限制在平台目录内。"""
        return check_path(v if os.path.isabs(v)
                          else os.path.join(PLATFORM_ROOT, v),
                          must_exist=True)

    def _req(key, what):
        v = (p.get(key) or '').strip()
        if not v:
            abort(400, f'请选择{what}')
        return _input_abs(v)

    def _opt(key):
        v = (p.get(key) or '').strip()
        return _input_abs(v) if v else None

    db_virus = _resolve_virus_db(p.get('db_virus'))
    if spec['need_db']:
        from vp.kunpeng import db_ready
        if not db_ready(db_virus):
            abort(400, f'所选病毒库未就绪: {db_virus}'
                       '（先到「数据库构建」页构建，或换选其他病毒库）')

    ts = time.strftime('%Y%m%d_%H%M%S')
    run_dir = check_path(os.path.join(_tool_runs_root(), f'{tool}_{ts}'),
                         must_exist=False, in_platform=True)
    os.makedirs(run_dir, exist_ok=True)
    # 数字参数钳制：负数/异常大值不直接传给外部工具
    if threads is not None:
        threads = max(1, min(int(threads), (os.cpu_count() or 4) * 4))

    ctx = SimpleNamespace(p=p, run_dir=run_dir, threads=threads,
                          db_virus=db_virus, req=_req, opt=_opt)
    job = spec['job'](ctx)

    name = cfg.tr(f"工具·{spec['title']} {ts}", f"Tool·{spec['title_en']} {ts}")
    run_name = os.path.basename(run_dir)

    def _job_with_run(log, prog, cancel):
        """包一层：结果 dict 注入 run 名，供预览面板列出产物文件。"""
        res = job(log, prog, cancel)
        if isinstance(res, dict) and 'run' not in res:
            res['run'] = run_name
        return res

    tid = tm.start(name, _job_with_run,
                   log_file=os.path.join(run_dir, 'run.log'),
                   weight=('light' if tool in LIGHT_TOOLS else 'heavy'))
    return jsonify({'task': tid, 'run': run_name})


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


@app.route('/api/tool/chain_result')
def api_tool_chain_result():
    """一键流程汇总目录：返回 _chain.json 内容 + 链接清单。

    chain = virchain / kvchain；run 为汇总目录名（留空取最新一次）。
    """
    chain = (request.args.get('chain') or '').strip()
    if chain not in ('virchain', 'kvchain'):
        abort(400, '无效的链式运行类型')
    run = (request.args.get('run') or '').strip()
    root = _tool_runs_root()
    if run:
        if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
            abort(400, '无效的运行名')
        if not run.startswith(chain + '_'):
            abort(400, '运行名与链式类型不匹配')
    else:
        cands = []
        if os.path.isdir(root):
            for d in os.listdir(root):
                if d.startswith(chain + '_') and os.path.isfile(
                        os.path.join(root, d, '_chain.json')):
                    cands.append(d)
        if not cands:
            return jsonify({'run': None, 'meta': None, 'entries': []})
        run = sorted(cands)[-1]
    cdir = check_path(os.path.join(root, run), must_exist=False,
                      in_platform=True)
    meta = None
    cj = os.path.join(cdir, '_chain.json')
    if os.path.isfile(cj):
        try:
            with safe_open(cj) as f:
                meta = json.load(f)
        except (OSError, ValueError):
            meta = None
    entries = []
    if os.path.isdir(cdir):
        for name in sorted(os.listdir(cdir)):
            if name == '_chain.json':
                continue
            full = os.path.join(cdir, name)
            entries.append({'name': name,
                            'is_dir': os.path.isdir(full),
                            'is_link': os.path.islink(full),
                            'size': (None if os.path.isdir(full)
                                     else fmt_size(os.path.getsize(full))
                                     if os.path.exists(full) else None)})
    return jsonify({'run': run, 'meta': meta, 'entries': entries})


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
@app.route('/logan')
def page_logan():
    return render_template('logan.html')


@app.route('/logan/report/<name>/')
def page_logan_report(name):
    from vp.logan_trace import _job_dir
    p = os.path.join(_job_dir(name), 'trace_report.html')
    if not os.path.isfile(p):
        abort(404, '溯源报告尚未生成（请先导入结果文件）')
    return send_file(check_path(p, must_exist=True, in_platform=True))


@app.route('/logan/report/<name>/<path:filename>')
def page_logan_report_file(name, filename):
    """溯源报告目录内静态文件（plotly.min.js / 原始结果表 / 查询 FASTA）。"""
    from vp.logan_trace import _job_dir
    p = check_path(os.path.join(_job_dir(name), filename),
                   must_exist=True, in_platform=True)
    return send_file(check_path(p, must_exist=True, in_platform=True))


@app.route('/api/logan/samples')
def api_logan_samples():
    """可选来源样品：已有 ③组装 病毒 contigs 产出的样品。"""
    from vp.logan_trace import list_query_samples
    return jsonify(list_query_samples())


@app.route('/api/logan/contigs/<sample>')
def api_logan_contigs(sample):
    from vp.logan_trace import list_virus_contigs
    return jsonify(list_virus_contigs(_safe_sample(sample)))


@app.route('/api/logan/jobs')
def api_logan_jobs():
    from vp.logan_trace import list_jobs
    return jsonify(list_jobs())


@app.route('/api/logan/job/<name>')
def api_logan_job(name):
    from vp.logan_trace import job_detail
    return jsonify(job_detail(name))


@app.route('/api/logan/job/<name>/segment/<int:idx>')
def api_logan_segment(name, idx):
    """片段序列（前端复制提交用）。"""
    from vp.logan_trace import get_segment_fasta
    header, seq = get_segment_fasta(name, idx)
    return jsonify({'header': header, 'seq': seq})


@app.route('/api/logan/job/<name>/query')
def api_logan_query_dl(name):
    """下载查询 FASTA（全部片段）。"""
    from vp.logan_trace import _job_dir
    p = os.path.join(_job_dir(name), 'query_all.fasta')
    if not os.path.isfile(p):
        abort(404, '查询 FASTA 不存在')
    return send_file(check_path(p, must_exist=True, in_platform=True),
                     as_attachment=True)


@app.route('/api/logan/create', methods=['POST'])
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


@app.route('/api/logan/job/<name>/import/<int:idx>', methods=['POST'])
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


@app.route('/api/logan/job/<name>/delete', methods=['POST'])
def api_logan_delete(name):
    from vp.logan_trace import delete_job
    delete_job(name)
    return jsonify({'ok': True})


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
@app.route('/download')
def page_download():
    return render_template('download.html')


@app.route('/api/dl/create', methods=['POST'])
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


@app.route('/api/dl/batches')
def api_dl_batches():
    return jsonify(get_dl_manager().list_snapshots())


@app.route('/api/dl/batch/<bid>')
def api_dl_batch(bid):
    snap = get_dl_manager().snapshot(bid)
    if not snap:
        abort(404, '批次不存在')
    return jsonify(snap)


@app.route('/api/dl/batch/<bid>/<action>', methods=['POST'])
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


# 兼容别名（路由函数体内使用）
get_dl_manager = _dl_manager


@app.route('/downloads/<bid>/<path:filename>')
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


@app.route('/api/dl/to_pipeline', methods=['POST'])
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
        sample_queue.add(created, None, {}, project=project)
        queued = len(created)
    return jsonify({'created': created, 'skipped': skipped, 'queued': queued})


# ------------------------------------------------------------------
# NCBI 提交准备（vp.ncbi_submit，源自 MMPV-RNA submission_gui 功能面）
# ------------------------------------------------------------------
@app.route('/submit')
def page_submit():
    return render_template('submit.html')


@app.route('/api/submit/tables')
def api_submit_tables():
    from vp.ncbi_submit import store
    return jsonify(store.list_tables())


@app.route('/api/submit/samples')
def api_submit_samples():
    """内置示例数据集（对应原 GUI Samples 菜单）。"""
    from vp.ncbi_submit import store
    return jsonify(store.list_samples())


@app.route('/api/submit/create', methods=['POST'])
def api_submit_create():
    from vp.ncbi_submit import store
    body = request.get_json(force=True) or {}
    name = str(body.get('name') or '').strip()
    sample = str(body.get('sample') or ('demo' if body.get('demo') else ''))
    try:
        store.create_table(name, sample=sample)
    except FileExistsError:
        abort(400, f'提交项目已存在: {name}')
    except ValueError as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'name': name})


@app.route('/api/submit/copy', methods=['POST'])
def api_submit_copy():
    """项目另存为（原 GUI Save As）。body: {src, dst}"""
    from vp.ncbi_submit import store
    body = request.get_json(force=True) or {}
    src = str(body.get('src') or '').strip()
    dst = str(body.get('dst') or '').strip()
    try:
        store.copy_table(src, dst)
    except FileNotFoundError:
        abort(404, f'源项目不存在: {src}')
    except FileExistsError:
        abort(400, f'提交项目已存在: {dst}')
    except ValueError as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'name': dst})


@app.route('/api/submit/upload', methods=['POST'])
def api_submit_upload():
    from vp.ncbi_submit import store
    name = str(request.form.get('name') or '').strip()
    f = request.files.get('file')
    if not f or not f.filename:
        abort(400, '缺少上传文件')
    import tempfile
    fd, tmp = tempfile.mkstemp(suffix=os.path.splitext(f.filename)[1])
    try:
        with os.fdopen(fd, 'wb') as out:
            f.save(out)
        try:
            _, rows = store.import_table(name, tmp)
        except FileExistsError:
            abort(400, f'提交项目已存在: {name}')
        except ValueError as e:
            abort(400, str(e))
        except Exception as e:
            abort(400, f'解析失败（需要 CSV/TSV/Excel）: {e}')
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass
    return jsonify({'ok': True, 'name': name, 'rows': rows})


def _submit_store(name):
    from vp.ncbi_submit import store
    try:
        store.table_dir(name)
    except ValueError as e:
        abort(400, str(e))
    return store


@app.route('/api/submit/table/<name>')
def api_submit_table(name):
    store = _submit_store(name)
    try:
        return jsonify(store.table_payload(name))
    except FileNotFoundError:
        abort(404, f'提交项目不存在: {name}')


@app.route('/api/submit/table/<name>/save', methods=['POST'])
def api_submit_table_save(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    cols = body.get('columns') or []
    raw = body.get('rows') or []
    rows = [dict(zip(cols, r)) for r in raw]
    try:
        n = store.save_table(name, rows)
    except ValueError as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'rows': n})


@app.route('/api/submit/table/<name>/rows', methods=['POST'])
def api_submit_table_rows(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        if body.get('add'):
            n = store.add_rows(name, body['add'])
        elif body.get('delete') is not None:
            n = store.delete_rows(name, body['delete'])
        else:
            abort(400, '需要 add 或 delete')
    except (ValueError, IndexError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'rows': n})


@app.route('/api/submit/table/<name>/fill', methods=['POST'])
def api_submit_table_fill(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    column = str(body.get('column') or '')
    value = str(body.get('value') or '')
    if not column or not value:
        abort(400, '需要 column 和 value')
    try:
        count = store.batch_fill(name, column, value, old_value=body.get('old'))
    except KeyError as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'count': count})


@app.route('/api/submit/table/<name>/export', methods=['POST'])
def api_submit_table_export(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        files = store.export_files(
            name,
            assembler=str(body.get('assembler') or 'SPAdes;4.3.0;metaviral'),
            sequencer=str(body.get('sequencer') or 'Illumina NovaSeq 6000'),
            enrichment=str(body.get('enrichment') or 'rRNA depletion'))
    except Exception as e:
        abort(400, f'生成失败: {e}')
    return jsonify({'ok': True, 'files': files})


@app.route('/api/submit/table/<name>/validate', methods=['POST'])
def api_submit_table_validate(name):
    store = _submit_store(name)
    return jsonify({'issues': store.validate_table(name)})


@app.route('/api/submit/table/<name>/fasta', methods=['POST'])
def api_submit_table_fasta(name):
    """从平台内 FASTA（contigs 运行 / 样品结果 / 任意平台内文件）
    按表行 sequence_name 提取序列 → 项目目录 sequences.fsa + 一致性报告。"""
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    src = str(body.get('fasta') or '').strip()
    if not src:
        abort(400, '请提供序列 FASTA 路径')
    try:
        p = check_path(src if os.path.isabs(src)
                       else os.path.join(PLATFORM_ROOT, src), must_exist=True)
        rep = store.export_submission_fasta(
            name, p, min_len=int(body.get('min_len', 200) or 200))
    except (ValueError, KeyError, FileNotFoundError, OSError) as e:
        abort(400, f'FASTA 导出失败: {e}')
    return jsonify({'ok': True, **rep})


@app.route('/api/submit/table/<name>/import_run', methods=['POST'])
def api_submit_table_import_run(name):
    """body: {run} → 从 tool_runs/<run>/virus_classification.tsv 追加表行。"""
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        res = store.import_from_run(name, str(body.get('run') or ''))
    except (ValueError, KeyError, FileNotFoundError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, **res})


@app.route('/api/submit/contig_runs')
def api_submit_contig_runs():
    from vp.ncbi_submit import store
    return jsonify(store.list_contig_runs())


@app.route('/api/submit/orf_runs')
def api_submit_orf_runs():
    """有 CDS 注释产物（04_orf/04b）的 orf/orfa 运行列表（关联注释用）。"""
    from vp.ncbi_submit import store
    return jsonify(store.list_orf_runs())


@app.route('/api/submit/table/<name>/link_orf', methods=['POST'])
def api_submit_table_link_orf(name):
    """把独立跑的 orf/orfa 注释运行关联到提交项目（按 contig 名对接）。"""
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        res = store.link_orf_run(name, str(body.get('run') or ''))
    except (ValueError, FileNotFoundError, KeyError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, **res})


@app.route('/api/submit/table/<name>/tbl', methods=['POST'])
def api_submit_table_tbl(name):
    """从已关联 orf 运行的 CDS 注释生成 featuretable.tbl（GenBank 特征表）。"""
    store = _submit_store(name)
    try:
        res = store.generate_feature_tbl(name)
    except ValueError as e:
        abort(400, str(e))
    return jsonify({'ok': True, **res})


@app.route('/api/submit/table/<name>/source_fasta')
def api_submit_table_source_fasta(name):
    """自动推断提交序列 FASTA 路径（viral_contigs / contigs.filtered）。"""
    store = _submit_store(name)
    try:
        p = store.infer_source_fasta(name)
    except (ValueError, OSError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'fasta': p})


@app.route('/api/submit/table/<name>/package', methods=['POST'])
def api_submit_table_package(name):
    """项目内全部提交产物 → submission_package.zip（一键下载用）。"""
    store = _submit_store(name)
    try:
        path, n = store.package_zip(name)
    except (ValueError, OSError) as e:
        abort(400, f'打包失败: {e}')
    return jsonify({'ok': True, 'files': n,
                    'download': f'/submissions/{name}/submission_package.zip'})


@app.route('/api/submit/table/<name>/taxonomy_check', methods=['POST'])
def api_submit_table_taxonomy_check(name):
    """organism 逐个查 NCBI Taxonomy（联网）；未收录物种会提交被拒。"""
    store = _submit_store(name)
    try:
        rows = store.taxonomy_check(name)
    except KeyError as e:
        abort(400, str(e))
    except RuntimeError as e:
        abort(502, f'NCBI 查询失败: {e}')
    return jsonify({'ok': True, 'rows': rows,
                    'missing': [r['organism'] for r in rows if not r['found']]})


@app.route('/api/submit/table/<name>/sbt', methods=['POST'])
def api_submit_table_sbt(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        tpl, path = store.generate_sbt(
            body.get('fields') or {},
            extra_authors=body.get('extra_authors') or [],
            title=body.get('title') or '', out_name=name)
    except ValueError as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'path': path,
                    'authors': 1 + len(body.get('extra_authors') or [])})


@app.route('/api/submit/table/<name>/sqn', methods=['POST'])
def api_submit_table_sqn(name):
    """本地 .sqn 生成（suvtk features→comments→table2asn / MIUVIG）。

    body: {author:{...}, use_features:bool, miuvig:{...},
           assembler, sequencer}
    author 缺省时回退到 store.generate_sbt 需要的必填字段（last/first/affil/
    city/country/email）；use_features=True 走 suvtk BFVD 链，False 用平台
    自家 featuretable.tbl（需先 /tbl 生成）。"""
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    from vp import suvtk_submit as ss
    import pandas as pd

    # 1) 关联的 orf 运行目录（CDS 坐标来源 + dominant_genome_type 推断）
    link = store._read_link(name, 'orf')
    if not link or not link.get('run'):
        abort(400, '尚未关联 orf 注释运行——先点"关联注释"')
    run_dir = store._orf_run_dir(link['run'])

    # 2) 表行 → build_sqn 的 rows（过滤占位符 sequence_name）
    df = store.load_table(name)
    rows = []
    for _, r in df.iterrows():
        sn = str(r.get('sequence_name') or '').strip()
        if not sn or store.is_placeholder(sn):
            continue
        rows.append({c: ('' if pd.isna(v) else v) for c, v in r.items()})
    if not rows:
        abort(400, '表内没有有效 sequence_name 行')

    # 3) 作者字段
    author = body.get('author') or {}

    # 4) 其余参数
    use_features = bool(body.get('use_features', True))
    # .sqn 及中间文件生成到提交项目目录（submissions/<name>/），与④文件预览/
    # 编辑一致——用户编辑 .sbt/.src/.cmt 后可重新生成 .sqn 并复用编辑。
    out_dir = store.table_dir(name)
    try:
        res = ss.build_sqn(
            run_dir, rows=rows, author=author,
            miuvig=body.get('miuvig'),
            assembler=body.get('assembler'),
            sequencer=body.get('sequencer'),
            use_features=use_features,
            tbl_path=os.path.join(out_dir, 'featuretable.tbl')
            if not use_features else None,
            out_dir=out_dir,
            log=lambda m: app.logger.info('[sqn] %s', m))
    except (RuntimeError, ValueError, FileNotFoundError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, **res})


@app.route('/api/submit/table/<name>/sqn_status')
def api_submit_table_sqn_status(name):
    """已有 .sqn 产物状态（含 .val 校验摘要），供 submit 页展示。

    .sqn 生成在提交项目目录（submissions/<name>/），与④文件预览一致。"""
    store = _submit_store(name)
    d = store.table_dir(name)
    sqn = os.path.join(d, 'sqn.sqn')
    val = os.path.join(d, 'sqn.val')
    if not os.path.isfile(sqn):
        return jsonify({'ok': True, 'exists': False, 'dir': d})
    errs = warns = infos = 0
    if os.path.isfile(val):
        with open(val, encoding='utf-8', errors='replace') as f:
            for line in f:
                s = line.strip()
                if s.startswith('Error'):
                    errs += 1
                elif s.startswith('Warning'):
                    warns += 1
                elif s:
                    infos += 1
    return jsonify({'ok': True, 'exists': True,
                    'sqn': sqn, 'val': val, 'dir': d,
                    'stats': {'error': errs, 'warning': warns, 'info': infos}})


@app.route('/api/submit/table/<name>/biosample', methods=['POST'])
def api_submit_table_biosample(name):
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    try:
        path, rows, skipped = store.export_biosample_tsv(
            name, skip_placeholders=not bool(body.get('include_placeholders')))
    except (ValueError, KeyError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'path': path, 'rows': rows, 'skipped': skipped})


@app.route('/api/submit/table/<name>/file')
def api_submit_table_file(name):
    store = _submit_store(name)
    fname = request.args.get('name') or ''
    try:
        return jsonify({'content': store.read_file(name, fname)})
    except (ValueError, FileNotFoundError) as e:
        abort(400, str(e))


@app.route('/api/submit/table/<name>/file_save', methods=['POST'])
def api_submit_table_file_save(name):
    """预览编辑保存（原 GUI Preview 标签的 Edit/Save）。body: {name, content}"""
    store = _submit_store(name)
    body = request.get_json(force=True) or {}
    fname = str(body.get('name') or '')
    try:
        store.write_file(name, fname, str(body.get('content') or ''))
    except (ValueError, FileNotFoundError) as e:
        abort(400, str(e))
    return jsonify({'ok': True, 'file': fname})


@app.route('/api/submit/table/<name>/delete', methods=['POST'])
def api_submit_table_delete(name):
    import shutil
    store = _submit_store(name)
    d = store.table_dir(name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    return jsonify({'ok': True})


@app.route('/submissions/<name>/<path:fname>')
def submissions_download(name, fname):
    """提交产物下载（限项目目录内）。"""
    from vp.ncbi_submit import store
    try:
        p = store.read_file_path(name, fname)
    except (ValueError, FileNotFoundError) as e:
        abort(404, str(e))
    return send_file(p, as_attachment=True)


# ------------------------------------------------------------------
# 公共数据检索（vp.public_meta，源自 MMPV-RNA public_metadata_pipeline）
# 单工作台模式：/meta 只操作一个当前项目（Search / Info 双存储，对应
# MMPV-RNA metadata_gui 的 8 个标签），每次在线检索自动新建
# <物种>_<时间戳> 项目，历史项目自动归档到结果中心查看。
# ------------------------------------------------------------------
# 三张标准表：search=检索结果（S 存储）core14=统一元数据（I 存储）
# full=全字段明细（只读，下载用）
_META_TABLE_FILES = {
    'search': ('search', 'SRA_GSA_Merged_Final.csv'),
    'core14': ('info', 'Global_Unified_Metadata_Core14.csv'),
    'full': ('info', 'Global_Unified_Metadata_Full.csv'),
}
_META_EDITABLE = ('search', 'core14')
_META_NA = {'na', 'n/a', 'not_provided', 'not_provided.', 'not collected',
            'missing', 'none', 'unknown', 'nan', '', '-', ' '}
_meta_df_cache = {}          # (name, which) -> [mtime, df]
_meta_cache_lock = threading.Lock()

# AI 元数据清洗提示词（与 MMPV-RNA metadata_gui/controllers/ai_completer.py
# 同源：只做去污染/归位/格式化，禁止凭空捏造）
AI_SANITIZE_PROMPT = """你是一个极其严谨的生物信息元数据清理程序。你的任务是【去污染、归位、格式化】，绝对不是【凭空捏造】。
我将给你一条 JSON 数据。请严格按以下规则逐字段处理：

1. **全面净化**：遇到 'not collected', 'missing', 'N/A', 'Not_Provided', 空字符串等无意义占位符，替换为 'Not_Provided'。

2. **ScientificName**：如果混入了部位（如 leaf）或来源词（如 wild），将其剔除，只保留纯物种名。

3. **Tissue**：强制小写单数（leaves→leaf, fruits→fruit, flowers→flower, roots→root）。如果原值为空但 Source 恰好是真正的部位词（如 leaf, root），才移动过来。Source 是品种名（如 Ningqi No.1）则不移动。

4. **Source**：如果被填成了物种名，清空为 Not_Provided。如果是 "wild", "cultivated" 或品种名（如 Ningqi No.1），必须 100% 保持原样。

5. **Location**：根据 CenterName 推断机构地理位置，输出严格三级格式【国家, 省/州, 市/县_AI】。参照示例：
   "Beijing Forestry University" → "China, Beijing, Beijing_AI"
   "Ningxia University" → "China, Ningxia, Yinchuan_AI"
   "North Minzu University" → "China, Ningxia, Yinchuan_AI"
   "University of Tokyo" → "Japan, Tokyo, Tokyo_AI"
   "USDA-ARS" → "USA, Maryland, Beltsville_AI"
   "Royal Botanic Gardens Kew" → "United Kingdom, England, London_AI"
   "Northwest A&F University" → "China, Shaanxi, Yangling_AI"
   "Gansu Agricultural University" → "China, Gansu, Lanzhou_AI"
   "Xinjiang University" → "China, Xinjiang, Urumqi_AI"
   "Qinghai University" → "China, Qinghai, Xining_AI"
   "Inner Mongolia University" → "China, Inner Mongolia, Hohhot_AI"
   "Henan Agricultural University" → "China, Henan, Zhengzhou_AI"
   不认识的机构保持 Not_Provided，绝对不要输出 Unknown。

6. **Age_GrowthStage**：标准化年龄和发育阶段文本。规则：
   - "3 year" → "3 years"（补全复数）
   - "2 years" → 保持原样（已规范）
   - "Archeocyte stage" → "archeocyte stage"（统一为首字母不大写的全小写，除非是专有名词）
   - 多个阶段用 "|" 分隔，如 "3 years | mature fruit stage"
   - 将中文风格描述转为英文，如 "flowering stage"、"ripening stage"
   - 如果为空则保持 Not_Provided

7. **LibrarySource**：如果传入长文本，浓缩为一个词：TRANSCRIPTOMIC / GENOMIC / METAGENOMIC。

8. **CollectionDate, CenterName, BioProject**：必须 100% 保持输入原样，一个字都不许改！

【核心纪律】：只做减法（去污染）、归位（移动）、格式标准化，绝对禁止凭空捏造！Location 没线索就 Not_Provided！

对输入 JSON 中出现的以下字段逐个输出清洗结果，直接输出合法 JSON（勿带 ```json 标记），未出现的字段不要输出：
{"CollectionDate": "...", "Location": "...", "Source": "...", "Tissue": "...",
 "Age_GrowthStage": "...", "ScientificName": "...", "LibrarySource": "...",
 "CenterName": "...", "BioProject": "..."}"""

AI_FILL_COLS = ('CollectionDate', 'Location', 'Source', 'Tissue',
                'Age_GrowthStage', 'ScientificName', 'LibrarySource',
                'CenterName', 'BioProject')
_AI_EMPTY = {'', ' ', 'na', 'n/a', 'not_provided', 'not collected',
             'missing', 'none', 'unknown', 'nan', '-'}


def _meta_slug(species):
    slug = re.sub(r'[^\w\-]+', '_', str(species)).strip('_') or 'species'
    return slug


def _meta_proj(name):
    return check_path(os.path.join(DIRS['meta_search'], name),
                      must_exist=False, in_platform=True)


def _meta_table_path(name, which):
    if which not in _META_TABLE_FILES:
        return None
    sub, fname = _META_TABLE_FILES[which]
    return os.path.join(_meta_proj(name), sub, fname)


def _meta_load_df(name, which):
    """读取项目表（带 mtime 缓存，编辑后由 _meta_store_df 刷新）。"""
    import pandas as pd
    path = _meta_table_path(name, which)
    if not path or not os.path.isfile(path):
        return None
    key = (name, which)
    mtime = os.path.getmtime(path)
    with _meta_cache_lock:
        ent = _meta_df_cache.get(key)
        if ent and ent[0] == mtime:
            return ent[1]
    df = pd.read_csv(path, encoding='utf-8-sig', dtype=str,
                     keep_default_na=False)
    df = df.loc[:, ~df.columns.duplicated()]
    with _meta_cache_lock:
        _meta_df_cache[key] = [mtime, df]
    return df


def _meta_store_df(name, which, df):
    """保存项目表（csv+tsv 双格式，与引擎 save_dual_format 一致）。"""
    sub, fname = _META_TABLE_FILES[which]
    base = check_path(os.path.join(_meta_proj(name), sub),
                      must_exist=False, in_platform=True)
    os.makedirs(base, exist_ok=True)
    stem = os.path.splitext(fname)[0]
    df = df.loc[:, ~df.columns.duplicated()]
    csv_p = os.path.join(base, stem + '.csv')
    with safe_open(csv_p, 'wt') as f:
        f.write('\ufeff')
        df.to_csv(f, index=False, lineterminator='\n')
    with safe_open(os.path.join(base, stem + '.tsv'), 'wt') as f:
        f.write('\ufeff')
        df.to_csv(f, index=False, sep='\t', lineterminator='\n')
    with _meta_cache_lock:
        _meta_df_cache[(name, which)] = [os.path.getmtime(csv_p), df]


def _meta_missing_count(df):
    return int(sum(df[c].astype(str).str.strip().str.lower()
                   .isin(_META_NA).sum() for c in df.columns))


def _meta_web_cache_dirs(name):
    """项目内全部 0_web_cache 目录（info / search 引擎的网页缓存）。"""
    out = []
    root = _meta_proj(name)
    if not os.path.isdir(root):
        return out
    for cur, sub, _fns in os.walk(root):
        if os.path.basename(cur) == '0_web_cache':
            out.append(cur)
    return out


def _meta_deepseek_key():
    return os.environ.get('DEEPSEEK_API_KEY', '').strip()


def _meta_ai_client(api_key, base='https://api.deepseek.com'):
    from openai import OpenAI
    return OpenAI(api_key=api_key, base_url=base)


def _meta_ai_chat(client, model, system, user, timeout=180):
    kwargs = {'model': model, 'timeout': timeout, 'messages': [
        {'role': 'system', 'content': system},
        {'role': 'user', 'content': user}]}
    if 'v4' in model.lower():
        if 'pro' in model.lower():
            kwargs['reasoning_effort'] = 'high'
            kwargs['extra_body'] = {'thinking': {'type': 'enabled'}}
        else:
            kwargs['temperature'] = 0.1
    else:
        kwargs['temperature'] = 0.3
    resp = client.chat.completions.create(**kwargs)
    return resp.choices[0].message.content or ''


def _meta_stream_task(task_name, cmd, done_msg, log_file=None, link=None):
    """子进程跑 meta 脚本，stdout 实时进任务日志；完成返回 summary。"""
    def fn(log, progress, cancel):
        from vp.utils import task_register_proc
        env = dict(os.environ)
        env.setdefault('PYTHONUNBUFFERED', '1')   # 引擎脚本日志实时流出
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True,
                                encoding='utf-8', errors='replace', env=env)
        task_register_proc(proc)                  # 登记 → 取消时硬停止
        try:
            for line in proc.stdout:
                if cancel.is_set():
                    proc.terminate()
                    raise RuntimeError('已取消')
                if line.strip():
                    log(line.rstrip())
            rc = proc.wait()
        finally:
            try:
                proc.stdout.close()
            except Exception:
                pass
        if rc != 0:
            raise RuntimeError(f'命令退出码 {rc}（详见日志）')
        return {'kind': 'meta', 'msg': done_msg}
    return tm.start(task_name, fn, log_file=log_file, link=link)


@app.route('/meta')
def page_meta():
    return render_template('meta.html')


def _csv_data_rows(path):
    """CSV 数据行数（按物理行快算，字段不含内嵌换行的场景下准确；
    供卡片徽章展示用，避免每次轮询都 pandas 全量解析）。"""
    try:
        n = sum(1 for _ in safe_open(path))
        return max(0, n - 1)
    except OSError:
        return None


@app.route('/api/meta/collections')
def api_meta_collections():
    root = check_path(DIRS['meta_search'], must_exist=False, in_platform=True)
    out = []
    if os.path.isdir(root):
        for name in sorted(os.listdir(root)):
            d = os.path.join(root, name)
            if not os.path.isdir(d):
                continue
            search_f = os.path.join(d, 'search', 'SRA_GSA_Merged_Final.csv')
            core_f = os.path.join(d, 'info', 'Global_Unified_Metadata_Core14.csv')
            full_f = os.path.join(d, 'info', 'Global_Unified_Metadata_Full.csv')
            n_runs = _csv_data_rows(search_f) if os.path.isfile(search_f) else None
            n_core = _csv_data_rows(core_f) if os.path.isfile(core_f) else None
            n_full = _csv_data_rows(full_f) if os.path.isfile(full_f) else None
            plot_dir = os.path.join(d, 'plot')
            n_plots = 0
            if os.path.isdir(plot_dir):
                n_plots = sum(1 for f in os.listdir(plot_dir)
                              if f.lower().endswith('.png'))
            # 过期检测：元数据条数 ≠ 最新检索条数 → 提示重新提取
            stale = (n_runs is not None and n_core is not None
                     and n_core != n_runs)
            out.append({'name': name,
                        'search': os.path.isfile(search_f),
                        'core14': os.path.isfile(core_f),
                        'full': os.path.isfile(full_f),
                        'n_runs': n_runs, 'n_core14': n_core,
                        'n_full': n_full, 'n_plots': n_plots,
                        'stale': stale,
                        'mtime': int(os.path.getmtime(d))})
    out.sort(key=lambda x: x['mtime'], reverse=True)
    return jsonify(out)


@app.route('/api/meta/search', methods=['POST'])
def api_meta_search():
    body = request.get_json(force=True) or {}
    species = str(body.get('species') or '').strip()
    if not species:
        abort(400, '需要物种拉丁名')
    source = str(body.get('source') or 'TRANSCRIPTOMIC')
    db = str(body.get('db') or 'both')
    detailed = bool(body.get('detailed', True))
    workers = body.get('workers') or 5
    try:
        workers = min(max(1, int(workers)), 20)
    except (TypeError, ValueError):
        workers = 5
    # 每次检索自动新建 <物种>_<时间戳> 项目；也允许 body.collection 指定项目重检
    name = str(body.get('collection') or '').strip()
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name or ''):
        ts = time.strftime('%Y%m%d_%H%M')
        name = f'{_meta_slug(species)}_{ts}'
        n = 2
        while os.path.isdir(os.path.join(DIRS['meta_search'], name)):
            name = f'{_meta_slug(species)}_{ts}_{n}'
            n += 1
    out_dir = check_path(os.path.join(DIRS['meta_search'], name, 'search'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    cmd = engine_cmd(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'vp', 'public_meta', 'search_engine.py'),
                     '-q', species, '-s', source, '--db', db,
                     '--workers', str(workers), '-o', out_dir)
    if not detailed:
        cmd.append('--no-detailed')
    if os.environ.get('NCBI_API_KEY'):
        cmd += ['--ncbi-api', os.environ['NCBI_API_KEY']]
    tid = _meta_stream_task(
        f'公共检索 {species}', cmd, f'检索完成: {out_dir}',
        log_file=os.path.join(out_dir, 'run.log'), link='/meta')
    return jsonify({'task': tid, 'name': name})


@app.route('/api/meta/info', methods=['POST'])
def api_meta_info():
    body = request.get_json(force=True) or {}
    name = str(body.get('collection') or '')
    try:
        coll = check_path(os.path.join(DIRS['meta_search'], name),
                          must_exist=True, in_platform=True)
    except (ValueError, FileNotFoundError):
        abort(400, f'检索项目不存在: {name}')
    runs = body.get('runs') or None
    list_f = os.path.join(coll, 'info_runs.txt')
    if runs:
        if not isinstance(runs, list) or not all(re.fullmatch(r'[A-Za-z0-9_]+', str(r)) for r in runs):
            abort(400, 'runs 需为编号列表')
        with safe_open(list_f, 'wt') as f:
            f.write('\n'.join(str(r) for r in runs) + '\n')
        input_arg = list_f
    else:
        search_f = os.path.join(coll, 'search', 'SRA_GSA_Merged_Final.csv')
        if not os.path.isfile(search_f):
            abort(400, '该项目还没有检索结果（先运行检索）')
        import pandas as pd
        try:
            df = pd.read_csv(search_f, encoding='utf-8-sig', dtype=str)
        except Exception as e:
            abort(400, f'检索结果解析失败: {e}')
        if 'Run' not in df.columns:
            abort(400, '检索结果缺少 Run 列')
        with safe_open(list_f, 'wt') as f:
            f.write('\n'.join(df['Run'].dropna().astype(str).str.strip()) + '\n')
        input_arg = list_f
    info_dir = check_path(os.path.join(coll, 'info'),
                          must_exist=False, in_platform=True)
    os.makedirs(info_dir, exist_ok=True)
    # 未提供 AI Key 时自动降级 local 模式（纯规则解析，不调用 LLM）
    # key 优先级：请求体 > 环境变量 DEEPSEEK_API_KEY > local 模式。
    # 环境变量避免密钥落在本机任意进程可见的命令行参数里。
    ds_key = (str(body.get('deepseek_api') or '').strip()
              or os.environ.get('DEEPSEEK_API_KEY', '').strip())
    mode = 'both' if ds_key else 'local'
    cmd = engine_cmd(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'vp', 'public_meta', 'info_engine.py'),
                     '-i', input_arg, '-o', info_dir, '-t', '4', '-m', mode)
    if ds_key:
        cmd += ['--deepseek-api', ds_key]
    if body.get('fill_date'):
        cmd.append('--fill-date')
    return jsonify({'task': _meta_stream_task(
        f'元数据提取 {name}', cmd, f'元数据完成: {info_dir}',
        log_file=os.path.join(info_dir, 'run.log'), link='/meta')})


@app.route('/api/meta/plot', methods=['POST'])
def api_meta_plot():
    """body: {name} → 后台任务：按项目元数据绘制 SCI 风格图
    （发布时间 / 机构 / 组织 / 地理 / 测序平台等，PNG 600dpi + PDF 矢量）。
    输入优先 Core14 元数据，缺失时回退 Full。"""
    body = request.get_json(force=True) or {}
    name = (body.get('name') or '').strip()
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name or ''):
        abort(400, '无效的检索项目名')
    coll = check_path(os.path.join(DIRS['meta_search'], name),
                      must_exist=True, in_platform=True)
    csv = None
    for cand in ('Global_Unified_Metadata_Core14.csv',
                 'Global_Unified_Metadata_Full.csv'):
        f = os.path.join(coll, 'info', cand)
        if os.path.isfile(f):
            csv = f
            break
    if not csv:
        abort(400, '该项目还没有元数据，请先「提取元数据」')
    out_dir = os.path.join(coll, 'plot')
    ts = time.strftime('%Y%m%d_%H%M%S')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.public_meta.landscape_plot import plot_sci_landscape
        prog('plot', 0.15, f'读取 {os.path.basename(csv)} 并绘制')
        plot_sci_landscape(csv, out_dir)
        n = len([f for f in os.listdir(out_dir)
                 if f.lower().endswith(('.png', '.pdf'))])
        prog('plot', 1.0, f'完成 {n} 个图表文件')
        logger.close()
        return {'kind': 'tool', 'stats': {'项目': name, '图表文件': n}}

    tid = tm.start(cfg.tr(f'检索绘图·{name} {ts}', f'Meta plot·{name} {ts}'), job,
                   weight='light', link='/meta')
    return jsonify({'task': tid})


@app.route('/api/meta/plot/files')
def api_meta_plot_files():
    """某检索项目的已生成图表列表（含可预览 URL）。"""
    name = request.args.get('name') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    plot_dir = check_path(os.path.join(DIRS['meta_search'], name, 'plot'),
                          must_exist=False, in_platform=True)
    files = []
    if os.path.isdir(plot_dir):
        for f in sorted(os.listdir(plot_dir)):
            if f.lower().endswith(('.png', '.pdf')):
                files.append({'name': f,
                              'url': f'/api/meta/plot/file?name={name}&file={f}',
                              'size': fmt_size(os.path.getsize(
                                  os.path.join(plot_dir, f)))})
    return jsonify({'files': files})


@app.route('/api/meta/plot/file')
def api_meta_plot_file():
    """图表文件服务（严格限制在该项目 plot 目录内，防穿越）。"""
    name = request.args.get('name') or ''
    fname = request.args.get('file') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name) or '..' in fname:
        abort(400, '无效参数')
    plot_dir = check_path(os.path.join(DIRS['meta_search'], name, 'plot'),
                          must_exist=False, in_platform=True)
    p = check_path(os.path.join(plot_dir, fname), must_exist=True,
                   in_platform=True)
    if not os.path.isfile(p):
        abort(404)
    return send_file(p)


@app.route('/api/meta/table/<name>')
def api_meta_table(name):
    """项目表数据（分页 + 排序 + 关键字/字段筛选 + 绝对行号供编辑）。
    q=全列关键字  qcol=限列关键字  f.<列>=值（可多个，AND 模糊匹配）"""
    which = request.args.get('which') or 'search'
    if which not in _META_TABLE_FILES:
        abort(400, 'which 需为 search/core14/full')
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name or ''):
        abort(400, '无效的检索项目名')
    df = _meta_load_df(name, which)
    if df is None:
        sub, _ = _META_TABLE_FILES[which]
        if which == 'search':
            abort(404, '还没有检索结果')
        abort(404, '还没有元数据结果（先提取元数据）')
    cols = list(df.columns)
    total = len(df)

    # 字段筛选（本地元数据筛选器的下拉组合，AND + 不区分大小写包含）
    for key, val in request.args.items():
        if key.startswith('f.'):
            col = key[2:]
            val = val.strip()
            if col in cols and val and val.lower() != '(any)':
                df = df[df[col].astype(str).str.lower()
                        .str.contains(val.lower(), regex=False)]
    # 快速关键字（全列或限单列）
    q = (request.args.get('q') or '').strip()
    if q:
        ql = q.lower()
        qcol = request.args.get('qcol') or ''
        if qcol and qcol in cols:
            df = df[df[qcol].astype(str).str.lower()
                    .str.contains(ql, regex=False)]
        else:
            mask = None
            for c in cols:
                m = df[c].astype(str).str.lower().str.contains(ql, regex=False)
                mask = m if mask is None else (mask | m)
            df = df[mask.fillna(False)] if mask is not None else df
    n_filtered = len(df)

    # 排序：数值列按数值、文本列按字典序；稳定排序保持原相对顺序
    sort_col = request.args.get('sort') or ''
    sort_dir = (request.args.get('dir') or 'asc').lower()
    if sort_col in cols:
        import pandas as pd
        num = pd.to_numeric(df[sort_col], errors='coerce')
        if num.notna().all():
            df = df.assign(_k=num).sort_values(
                '_k', ascending=sort_dir != 'desc').drop(columns='_k')
        else:
            df = df.sort_values(sort_col, ascending=sort_dir != 'desc',
                                kind='stable')
    try:
        page = max(1, int(request.args.get('page') or 1))
    except ValueError:
        page = 1
    try:
        per_page = min(max(1, int(request.args.get('per_page') or 20)), 500)
    except ValueError:
        per_page = 20
    pages = max(1, (n_filtered + per_page - 1) // per_page)
    page = min(page, pages)
    view = df.iloc[(page - 1) * per_page: page * per_page]
    rows = view.values.tolist()
    row_ids = [int(i) for i in view.index]      # 绝对行号 → 编辑定位
    run_col = cols.index('Run') if 'Run' in cols else None
    return jsonify({
        'columns': cols, 'rows': rows, 'row_ids': row_ids,
        'total': total, 'filtered': n_filtered,
        'page': page, 'per_page': per_page, 'pages': pages,
        'run_col': run_col,
        'n_missing': _meta_missing_count(_meta_load_df(name, which)),
        'download': f'/meta_files/{name}/{_META_TABLE_FILES[which][0]}/'
                    f'{_META_TABLE_FILES[which][1]}'})


@app.route('/api/meta/unique/<name>')
def api_meta_unique(name):
    """列去重取值（本地筛选下拉，过滤占位符，最多 200 项）。"""
    which = request.args.get('which') or 'search'
    if which not in _META_TABLE_FILES:
        abort(400, 'which 需为 search/core14/full')
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name or ''):
        abort(400, '无效的检索项目名')
    df = _meta_load_df(name, which)
    if df is None:
        return jsonify({'values': {}})
    cols_req = [c for c in (request.args.get('cols') or '').split(',') if c]
    out = {}
    for col in cols_req:
        if col not in df.columns:
            continue
        s = df[col].astype(str).str.strip()
        s = s[~s.str.lower().isin(_META_NA)]
        out[col] = sorted(s.unique().tolist())[:200]
    return jsonify({'values': out})


def _meta_editable_target():
    body = request.get_json(force=True) or {}
    name = str(body.get('name') or '')
    which = str(body.get('which') or '')
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    if which not in _META_EDITABLE:
        abort(400, '该表为只读明细，仅 search / core14 可编辑')
    df = _meta_load_df(name, which)
    if df is None:
        abort(404, '表不存在（先检索/提取元数据）')
    return body, name, which, df


@app.route('/api/meta/cell', methods=['POST'])
def api_meta_cell():
    """编辑单个单元格（双击编辑）。body: {name, which, row, col, value}"""
    body, name, which, df = _meta_editable_target()
    cols = list(df.columns)
    col = str(body.get('col') or '')
    try:
        row = int(body.get('row'))
    except (TypeError, ValueError):
        abort(400, 'row 需为整数')
    if col not in cols or not (0 <= row < len(df)):
        abort(400, '无效的单元格位置')
    value = '' if body.get('value') is None else str(body.get('value'))
    df.iat[row, cols.index(col)] = value
    _meta_store_df(name, which, df)
    return jsonify({'ok': True, 'value': value,
                    'n_missing': _meta_missing_count(df)})


@app.route('/api/meta/row', methods=['POST'])
def api_meta_row():
    """整行多字段保存（编辑面板）。body: {name, which, row, data:{col:value}}"""
    body, name, which, df = _meta_editable_target()
    data = body.get('data') or {}
    try:
        row = int(body.get('row'))
    except (TypeError, ValueError):
        abort(400, 'row 需为整数')
    if not (0 <= row < len(df)):
        abort(400, '无效的行号')
    applied = 0
    for col, val in data.items():
        if col in df.columns:
            df.iat[row, df.columns.get_loc(col)] = \
                '' if val is None else str(val)
            applied += 1
    _meta_store_df(name, which, df)
    return jsonify({'ok': True, 'applied': applied,
                    'n_missing': _meta_missing_count(df)})


@app.route('/api/meta/rows/delete', methods=['POST'])
def api_meta_rows_delete():
    """删除选中行。body: {name, which, rows:[绝对行号...]}"""
    body, name, which, df = _meta_editable_target()
    rows = body.get('rows') or []
    try:
        rows = sorted({int(r) for r in rows})
    except (TypeError, ValueError):
        abort(400, 'rows 需为行号数组')
    if not rows:
        abort(400, '未选择行')
    if rows[-1] >= len(df) or rows[0] < 0:
        abort(400, '行号越界')
    df = df.drop(index=rows).reset_index(drop=True)
    _meta_store_df(name, which, df)
    return jsonify({'ok': True, 'total': len(df),
                    'n_missing': _meta_missing_count(df)})


@app.route('/api/meta/record')
def api_meta_record():
    """单行完整记录（编辑面板）+ 该 Run 的网页缓存 JSON（若有）。"""
    which = request.args.get('which') or 'core14'
    if which not in _META_TABLE_FILES:
        abort(400, 'which 需为 search/core14/full')
    name = request.args.get('name') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    df = _meta_load_df(name, which)
    if df is None:
        abort(404, '表不存在')
    try:
        row = int(request.args.get('row') or 0)
    except ValueError:
        row = 0
    if not (0 <= row < len(df)):
        abort(400, '行号越界')
    rec = {c: str(df.iat[row, df.columns.get_loc(c)]) for c in df.columns}
    cache = None
    run = rec.get('Run', '').strip()
    if run:
        for d in _meta_web_cache_dirs(name):
            jp = os.path.join(d, f'{run}.json')
            if os.path.isfile(jp):
                try:
                    with safe_open(jp) as f:
                        cache = json.load(f)
                except Exception:
                    cache = None
                break
    return jsonify({'row': rec, 'total': len(df), 'row_index': row,
                    'cache': cache})


@app.route('/api/meta/fill_cache', methods=['POST'])
def api_meta_fill_cache():
    """从本地网页缓存 JSON 回填缺失字段（对应 GUI「Fill from Cache」）。"""
    body, name, which, df = _meta_editable_target()
    mapping = {'SubmissionDate': 'ReleaseDate', 'Organization': 'CenterName',
               'PRJ': 'BioProject', 'TaxID': 'TaxID',
               'ScientificName': 'ScientificName'}
    if 'Run' not in df.columns:
        abort(400, '表缺少 Run 列')
    cache_dirs = _meta_web_cache_dirs(name)
    filled = 0
    for idx in range(len(df)):
        run = str(df.iat[idx, df.columns.get_loc('Run')]).strip()
        if not run:
            continue
        data = None
        for d in cache_dirs:
            jp = os.path.join(d, f'{run}.json')
            if os.path.isfile(jp):
                try:
                    with safe_open(jp) as f:
                        data = json.load(f)
                except Exception:
                    data = None
                break
        if not data:
            continue
        for jk, col in mapping.items():
            if col not in df.columns:
                continue
            cur = str(df.iat[idx, df.columns.get_loc(col)]).strip()
            val = str(data.get(jk, '') or '').strip()
            if val and cur.lower() in _META_NA:
                df.iat[idx, df.columns.get_loc(col)] = val
                filled += 1
    if filled:
        _meta_store_df(name, which, df)
    return jsonify({'ok': True, 'filled': filled})


@app.route('/api/meta/ai_fill', methods=['POST'])
def api_meta_ai_fill():
    """AI 清洗/补全缺失元数据（DeepSeek，后台任务；绿色标记由前端追踪）。"""
    body, name, which, df = _meta_editable_target()
    key = _meta_deepseek_key()
    if not key:
        abort(400, '需要 DEEPSEEK_API_KEY 环境变量（或用 CLI --deepseek-api）')
    model = str(body.get('model') or 'deepseek-v4-flash')
    cols = [c for c in AI_FILL_COLS if c in df.columns]
    if not cols:
        abort(400, '表中没有可清洗的元数据字段')
    records = []
    for i in range(len(df)):
        row = {c: str(df.iat[i, df.columns.get_loc(c)]).strip() for c in cols}
        if any(v and v.lower() not in _AI_EMPTY for v in row.values()):
            records.append((i, row))

    def job(log, prog, cancel):
        client = _meta_ai_client(key)
        total, changed, done = len(records), 0, 0
        cells = []
        log(f'AI 清洗 {total} 条记录（模型 {model}）…')
        for n, (idx, rec) in enumerate(records):
            if cancel.is_set():
                raise RuntimeError('已取消')
            try:
                text = _meta_ai_chat(client, model, AI_SANITIZE_PROMPT,
                                     json.dumps(rec, ensure_ascii=False))
                text = re.sub(r'<think>.*?</think>', '', text, flags=re.S).strip()
                text = re.sub(r'^```(?:json)?\s*|\s*```$', '', text,
                              flags=re.I).strip()
                result = json.loads(text)
                for col in cols:
                    if col not in result:
                        continue
                    old = rec.get(col, '').strip()
                    new = str(result[col]).strip()
                    if not new or old == new or new.lower() in _AI_EMPTY:
                        continue
                    if col == 'Location' and all(
                            'unknown' in p.lower() for p in new.split(',')):
                        continue
                    df.iat[idx, df.columns.get_loc(col)] = new
                    changed += 1
                    cells.append([idx, col])
            except json.JSONDecodeError as e:
                log(f'第 {idx + 1} 行：AI 返回非 JSON：{str(e)[:120]}')
            except Exception as e:
                log(f'第 {idx + 1} 行：{str(e)[:200]}')
            done += 1
            prog('ai', done / max(total, 1),
                 f'AI 清洗 {done}/{total}（已改 {changed} 字段）')
            if done % 10 == 0:
                log(f'进度 {done}/{total}，已修改 {changed} 个字段')
        if changed:
            _meta_store_df(name, which, df)
        return {'kind': 'meta', 'stats': {'记录': total, '修改字段': changed,
                                          'cells': cells[:5000]}}

    ts = time.strftime('%H:%M:%S')
    tid = tm.start(cfg.tr(f'AI 元数据补全·{name} {ts}',
                          f'AI metadata fill·{name} {ts}'), job,
                   weight='light', link='/meta')
    return jsonify({'task': tid, 'records': len(records)})


@app.route('/api/meta/ai_summary', methods=['POST'])
def api_meta_ai_summary():
    """基于表统计生成 SCI 数据描述段落（DeepSeek，同步返回）。"""
    body = request.get_json(force=True) or {}
    name = str(body.get('name') or '')
    which = str(body.get('which') or 'core14')
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    if which not in _META_TABLE_FILES:
        abort(400, 'which 需为 search/core14/full')
    df = _meta_load_df(name, which)
    if df is None or df.empty:
        abort(404, '没有数据')
    key = (str(body.get('deepseek_api') or '').strip()
           or _meta_deepseek_key())
    if not key:
        abort(400, '需要 DEEPSEEK_API_KEY 环境变量')

    def counter(col):
        if col not in df.columns:
            from collections import Counter
            return Counter()
        s = df[col].astype(str).str.strip()
        s = s[~s.str.lower().isin(_META_NA)]
        from collections import Counter
        return Counter(s)

    def top(c, n=10):
        return ', '.join(f'{k}({v})' for k, v in c.most_common(n))

    years = set()
    date_col = ('CollectionDate' if 'CollectionDate' in df.columns
                else 'ReleaseDate' if 'ReleaseDate' in df.columns else None)
    if date_col:
        for v in df[date_col].astype(str):
            p = v.replace('/', '-').split('-')
            if p and p[0].isdigit() and len(p[0]) == 4:
                years.add(p[0])
    total_gb, gb_n = 0.0, 0
    size_col = ('FileSize_GB' if 'FileSize_GB' in df.columns
                else 'FileSize_MB' if 'FileSize_MB' in df.columns else None)
    if size_col:
        for v in df[size_col]:
            try:
                x = float(str(v).strip())
                total_gb += x if size_col.endswith('_GB') else x / 1024
                gb_n += 1
            except (TypeError, ValueError):
                pass
    vol = (f'Data volume: {total_gb:.1f} GB total '
           f'({gb_n} runs with size data, avg {total_gb / gb_n:.1f} GB/run)'
           if gb_n else '')
    db_counts = counter('Database')
    stats_text = f"""Total records: {len(df)}
{vol}
Database sources: {top(db_counts) or 'N/A'}
Species: {top(counter('ScientificName'), 8) or 'N/A'}
Tissues: {top(counter('Tissue'), 10) or 'N/A'}
Source types: {top(counter('Source'), 8) or 'N/A'}
Locations: {top(counter('Location'), 10) or 'N/A'}
Institutions: {top(counter('CenterName'), 10) or 'N/A'}
Growth stages: {top(counter('Age_GrowthStage'), 5) or 'N/A'}
Collection years: {', '.join(sorted(years)) or 'N/A'}"""

    prompt = f"""You are a scientific writer preparing a manuscript for a virome/metagenomics study.
Based on the following metadata statistics of public sequencing runs collected for analysis, write a concise paragraph (150-250 words) suitable for the "Data Collection" or "Sample Information" section of a scientific paper.

{stats_text}

Requirements:
- Write in formal scientific English, past tense
- Include total number of runs, database split (SRA/GSA), and total data volume (in GB) if available
- Include species covered, tissue types, geographic locations, and collection time span
- Mention key institutions that contributed data
- Note the sequencing type (transcriptomic/genomic) and average data volume per run when available
- End with a note on data availability

Output ONLY the paragraph, no markdown, no headings."""
    try:
        model = str(body.get('model') or 'deepseek-v4-flash')
        text = _meta_ai_chat(_meta_ai_client(key), model,
                             'You are a scientific writer. '
                             'Output only the requested paragraph.', prompt)
    except Exception as e:
        abort(500, f'AI 总结失败: {e}')
    return jsonify({'summary': text.strip()})


@app.route('/api/meta/import', methods=['POST'])
def api_meta_import():
    """导入 TSV/CSV/Excel 替换当前项目表（对应 GUI File > Import）。"""
    name = request.form.get('name') or ''
    which = request.form.get('which') or 'search'
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    if which not in _META_EDITABLE:
        abort(400, '仅 search / core14 表可导入')
    f = request.files.get('file')
    if not f or not f.filename:
        abort(400, '需要上传文件')
    import tempfile
    import pandas as pd
    suffix = os.path.splitext(f.filename)[1].lower()
    if suffix not in ('.csv', '.tsv', '.txt', '.xlsx', '.xls'):
        abort(400, '支持 .csv / .tsv / .txt / .xlsx / .xls')
    fd, tmp = tempfile.mkstemp(suffix=suffix)
    try:
        f.save(tmp)
        if suffix in ('.xlsx', '.xls'):
            df = pd.read_excel(tmp, dtype=str, keep_default_na=False)
        else:
            sep = '\t' if suffix in ('.tsv', '.txt') else ','
            try:
                df = pd.read_csv(tmp, sep=sep, dtype=str,
                                 keep_default_na=False, encoding='utf-8-sig')
            except UnicodeDecodeError:
                df = pd.read_csv(tmp, sep=sep, dtype=str,
                                 keep_default_na=False, encoding='gbk')
    except Exception as e:
        abort(400, f'解析失败: {e}')
    finally:
        try:
            os.close(fd)
            os.remove(tmp)
        except OSError:
            pass
    df = df.loc[:, ~df.columns.duplicated()].fillna('')
    _meta_store_df(name, which, df.astype(object))
    return jsonify({'ok': True, 'rows': len(df),
                    'columns': list(df.columns)})


@app.route('/api/meta/export')
def api_meta_export():
    """导出项目表（csv 原样 / tsv / xlsx）。"""
    which = request.args.get('which') or 'core14'
    fmt = (request.args.get('fmt') or 'csv').lower()
    name = request.args.get('name') or ''
    if which not in _META_TABLE_FILES or fmt not in ('csv', 'tsv', 'xlsx'):
        abort(400, '无效参数')
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    sub, fname = _META_TABLE_FILES[which]
    if fmt in ('csv', 'tsv'):
        ext = '.csv' if fmt == 'csv' else '.tsv'
        p = check_path(os.path.join(_meta_proj(name), sub,
                                    os.path.splitext(fname)[0] + ext),
                       must_exist=True, in_platform=True)
        return send_file(p, as_attachment=True,
                         download_name=f'{name}_{which}{ext}')
    import io
    import pandas as pd
    df = _meta_load_df(name, which)
    if df is None:
        abort(404, '表不存在')
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine='openpyxl') as w:
        df.to_excel(w, index=False)
    buf.seek(0)
    return send_file(buf, as_attachment=True,
                     download_name=f'{name}_{which}.xlsx',
                     mimetype='application/vnd.openxmlformats-officedocument'
                              '.spreadsheetml.sheet')


@app.route('/api/meta/overview')
def api_meta_overview():
    """图表数据包：各列 Top 取值 / 缺失率 / 年份分布 / 汇总统计。"""
    which = request.args.get('which') or 'core14'
    if which not in _META_TABLE_FILES:
        abort(400, 'which 需为 search/core14/full')
    name = request.args.get('name') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    df = _meta_load_df(name, which)
    if df is None or df.empty:
        abort(404, '表不存在或为空')
    fields = {}
    for c in df.columns:
        s = df[c].astype(str).str.strip()
        miss = int(s.str.lower().isin(_META_NA).sum())
        s2 = s[~s.str.lower().isin(_META_NA)]
        vc = s2.value_counts().head(50)
        fields[c] = {'counts': [[k, int(v)] for k, v in vc.items()],
                     'missing': miss,
                     'unique': int(s2.nunique())}
    date_col = ('CollectionDate' if 'CollectionDate' in df.columns
                and df['CollectionDate'].astype(str).str.strip()
                .str.lower().isin(_META_NA).sum() < len(df)
                else 'ReleaseDate' if 'ReleaseDate' in df.columns else None)
    years = {}
    if date_col:
        for v in df[date_col].astype(str):
            p = v.replace('/', '-').split('-')
            if p and p[0].isdigit() and len(p[0]) == 4:
                years[p[0]] = years.get(p[0], 0) + 1
    total, ncol = len(df), len(df.columns)
    filled = total * ncol - _meta_missing_count(df)
    complete = 0
    for _, row in df.iterrows():
        if not any(str(v).strip().lower() in _META_NA for v in row):
            complete += 1
    total_gb, gb_n = 0.0, 0
    size_col = ('FileSize_GB' if 'FileSize_GB' in df.columns
                else 'FileSize_MB' if 'FileSize_MB' in df.columns else None)
    if size_col:
        for v in df[size_col]:
            try:
                x = float(str(v).strip())
                total_gb += x if size_col.endswith('_GB') else x / 1024
                gb_n += 1
            except (TypeError, ValueError):
                pass
    return jsonify({
        'columns': list(df.columns), 'fields': fields,
        'years': sorted(years.items()),
        'date_col': date_col,
        'summary': {'total': total, 'cols': ncol, 'filled': filled,
                    'missing': total * ncol - filled, 'complete': complete,
                    'total_gb': round(total_gb, 1) if gb_n else None,
                    'gb_runs': gb_n}})


@app.route('/api/meta/collection/<name>/delete', methods=['POST'])
def api_meta_collection_delete(name):
    import shutil
    try:
        coll = check_path(os.path.join(DIRS['meta_search'], name),
                          must_exist=True, in_platform=True)
    except (ValueError, FileNotFoundError):
        abort(404, f'检索项目不存在: {name}')
    shutil.rmtree(coll)
    with _meta_cache_lock:
        for key in [k for k in _meta_df_cache if k[0] == name]:
            _meta_df_cache.pop(key, None)
    return jsonify({'ok': True})


@app.route('/meta_files/<name>/<path:subpath>')
def meta_files_download(name, subpath):
    """meta 产物下载（限 meta_search/<name>/ 目录内）。"""
    try:
        coll = check_path(os.path.join(DIRS['meta_search'], name),
                          must_exist=True, in_platform=True)
        p = check_path(os.path.join(coll, subpath), must_exist=True,
                       in_platform=True)
        if not p.startswith(coll + os.sep):
            abort(404)
    except (ValueError, FileNotFoundError):
        abort(404)
    return send_file(p, as_attachment=True)


# ------------------------------------------------------------------
# LOGAN 批量自动提交（logan_submit.py，Selenium 驱动本机浏览器）
# ------------------------------------------------------------------
@app.route('/api/logan/batch_ready')
def api_logan_batch_ready():
    """批量模式可用性（selenium 是否已安装）。"""
    from vp.logan_trace import selenium_ready, GROUP_HINTS
    return jsonify({'selenium': selenium_ready(), 'groups': GROUP_HINTS})


@app.route('/api/logan/job/<name>/batch', methods=['POST'])
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

# -*- coding: utf-8 -*-
"""样品与批处理队列（自 app.py 拆出）。

SampleQueue 顺序队列、样品 CRUD、管道逐级运行、分析任务提交。"""
import json
import logging
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

_log = logging.getLogger('vp.web.samples')

bp = Blueprint('samples', __name__)


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
# 注册到 state，供 download blueprint 的「转入分析流程」取用
from vp.web import state as _web_state

_web_state.set_sample_queue(sample_queue)


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


@bp.route('/api/queue')
def api_queue():
    return jsonify(sample_queue.snapshot())


@bp.route('/api/queue/add', methods=['POST'])
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


@bp.route('/api/queue/<entry_id>/remove', methods=['POST'])
def api_queue_remove(entry_id):
    if not sample_queue.remove(entry_id):
        abort(400, '条目不存在或正在运行（请先取消任务）')
    return jsonify({'ok': True})


@bp.route('/api/queue/clear_finished', methods=['POST'])
def api_queue_clear():
    sample_queue.clear_finished()
    return jsonify({'ok': True})


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


@bp.route('/api/analyze', methods=['POST'])
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


def _sample_dir(sample):
    from vp.pipeline import _safe_sample_name
    from vp.config import DIRS
    s = _safe_sample_name(sample)
    return s, check_path(os.path.join(DIRS['results'], s),
                         must_exist=True, in_platform=True)


@bp.route('/api/samples')
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


@bp.route('/api/samples/archived')
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


@bp.route('/api/pipeline/create', methods=['POST'])
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


@bp.route('/api/pipeline/<sample>')
def api_pipeline(sample):
    from vp.pipeline import pipeline_overview, load_sample_input
    s, sd = _sample_dir(sample)
    r1, r2, project = load_sample_input(sd)
    lang = request.args.get('lang') or None
    return jsonify({'sample': s, 'r1': r1, 'r2': r2, 'project': project or '',
                    **pipeline_overview(sd, lang)})


@bp.route('/api/pipeline/<sample>/run', methods=['POST'])
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


@bp.route('/api/open_report_dir/<sample>')
def api_open_report_dir(sample):
    safe = _safe_sample(sample)
    d = check_path(os.path.join(DIRS['results'], safe), must_exist=True,
                   in_platform=True)
    os.startfile(check_path(d, must_exist=True, in_platform=True))
    return jsonify({'ok': True})


@bp.route('/api/samples/<sample>/<path:rel>')
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


@bp.route('/api/samples/<sample>/clear', methods=['POST'])
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
        _log.warning('样品清单归零失败 %s: %s', safe, e)
    return jsonify({'ok': True})


@bp.route('/api/samples/<sample>/delete', methods=['POST'])
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


@bp.route('/api/sample_files/<sample>')
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

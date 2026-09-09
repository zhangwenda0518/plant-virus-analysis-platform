# -*- coding: utf-8 -*-
"""
分析流程编排器：按阶段依赖链调度，支持断点续跑、加权总进度、
预计剩余时间（历史耗时模型）与逐阶段资源预估日志（GUI 用）。
阶段: subsample → ⓪fastp → ⓪b fq2fa → ①host → ②virus → ③assembly
      → ④hostana → ⑥orf → ⑦phylo → ⑧primer
      → ⑨基因组图(gbdraw/DFV) → ⑩report
"""
import os
import re
import json
import time

from .config import DIRS
from .utils import (check_path, safe_open, iter_fastq_records,
                    log_res_plan, fmt_eta, is_step_done, mark_step_done)

STAGE_ORDER = ['subsample', 'fastp', 'fq2fa', 'host', 'virus', 'assembly',
               'verify', 'consensus', 'hostana', 'orf', 'orfa', 'phylo',
               'primer', 'gbdraw', 'report']
DEFAULT_ANALYZE_STAGES = [
    'fq2fa', 'host', 'virus', 'assembly', 'verify', 'consensus', 'hostana',
    'orf', 'orfa', 'phylo', 'primer', 'gbdraw', 'report',
]
STAGE_NAMES = {
    'subsample': '预处理(子采样)', 'fastp': '⓪ Fastp 质控',
    'fq2fa': '⓪b 序列转换(FASTQ→FASTA)',
    'host': '① 宿主去除', 'virus': '② 病毒筛查与提取',
    'assembly': '③ 组装·分类·提取', 'verify': '③b 候选序列验证',
    'consensus': '③c 共识序列与变异',
    'hostana': '④ 宿主预测(ICTV)',
    'orf': '⑥ ORF 预测',
    'orfa': '⑥b ORF 功能注释',
    'phylo': '⑦ 进化树与 SDT', 'primer': '⑧ 引物设计',
    'gbdraw': '⑨ 基因组图(gbdraw/DFV)',
    'report': '⑩ 可视化报告',
}
STAGE_NAMES_EN = {
    'subsample': 'Prep (subsample)', 'fastp': '⓪ Fastp QC',
    'fq2fa': '⓪b FASTQ→FASTA',
    'host': '① Host removal', 'virus': '② Virus screening',
    'assembly': '③ Assembly & extraction', 'verify': '③b Candidate verify',
    'consensus': '③c Consensus & variants',
    'hostana': '④ Host prediction (ICTV)',
    'orf': '⑥ ORF prediction',
    'orfa': '⑥b ORF annotation',
    'phylo': '⑦ Phylogeny & SDT', 'primer': '⑧ Primer design',
    'gbdraw': '⑨ Genome plots (gbdraw/DFV)',
    'report': '⑩ Visual report',
}
# 管道按功能模块分组显示（PhyloSuite 风格）
STAGE_GROUPS = [
    ('🧹 测序数据预处理', ['subsample', 'fastp', 'fq2fa', 'host']),
    ('🦠 病毒鉴定', ['virus']),
    ('🧬 病毒组装', ['assembly', 'verify', 'consensus']),
    ('🧲 宿主预测', ['hostana']),
    ('🔬 下游分析', ['orf', 'orfa', 'phylo', 'primer', 'gbdraw']),
    ('📄 报告', ['report']),
]
STAGE_GROUPS_EN = [
    ('🧹 Read preprocessing', ['subsample', 'fastp', 'fq2fa', 'host']),
    ('🦠 Virus identification', ['virus']),
    ('🧬 Virus assembly', ['assembly', 'verify', 'consensus']),
    ('🧲 Host prediction', ['hostana']),
    ('🔬 Downstream analysis', ['orf', 'orfa', 'phylo', 'primer', 'gbdraw']),
    ('📄 Report', ['report']),
]


def stage_name(stage, lang=None):
    """阶段显示名（lang 缺省取平台配置语言）。"""
    if lang is None:
        from .config import get_config
        lang = get_config().lang
    table = STAGE_NAMES if lang == 'zh' else STAGE_NAMES_EN
    return table.get(stage) or STAGE_NAMES.get(stage) or stage

# 阶段耗时权重（全局进度占比初值；有历史耗时后按实测覆盖）
STAGE_WEIGHTS = {
    'fastp': 3, 'fq2fa': 2, 'host': 22, 'virus': 18, 'assembly': 20,
    'verify': 8, 'consensus': 6, 'hostana': 3, 'orf': 4, 'orfa': 5, 'phylo': 6,
    'primer': 3, 'gbdraw': 2, 'report': 3,
}
# 无历史记录时的每阶段粗估：('gb', 秒/GB输入) 按输入量缩放，('plain', 秒) 固定
DEFAULT_STAGE_EST = {
    'fastp': ('gb', 30), 'fq2fa': ('gb', 15), 'host': ('gb', 210),
    'virus': ('gb', 170), 'assembly': ('plain', 600), 'hostana': ('plain', 60),
    'verify': ('plain', 240),
    'consensus': ('plain', 180),
    'orf': ('plain', 120), 'orfa': ('plain', 180),
    'phylo': ('plain', 180), 'primer': ('plain', 60), 'gbdraw': ('plain', 45),
    'report': ('plain', 30),
}
_GB_STAGES = {s for s, (k, _v) in DEFAULT_STAGE_EST.items() if k == 'gb'}
PERF_FILE = os.path.join(DIRS['logs'], 'stage_perf.json')


def _load_stage_perf():
    try:
        with safe_open(PERF_FILE) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def _update_stage_perf(perf, stage, seconds, input_bytes):
    """EMA 更新每阶段耗时模型（gb 阶段存每 GB 秒数，其余存固定秒数）。"""
    if stage in _GB_STAGES:
        gb = max(input_bytes / 1e9, 0.05)
        rate = max(seconds / gb, 1.0)
        old = perf.get(stage, {}).get('spg')
        val = rate if old is None else old * 0.6 + rate * 0.4
        perf.setdefault(stage, {})['spg'] = round(val, 2)
    else:
        old = perf.get(stage, {}).get('plain')
        val = seconds if old is None else old * 0.6 + seconds * 0.4
        perf.setdefault(stage, {})['plain'] = round(val, 1)
    try:
        with safe_open(PERF_FILE, 'wt') as f:
            json.dump(perf, f, ensure_ascii=False, indent=1)
    except OSError:
        pass


def _est_stage_sec(perf, stage, input_bytes=0):
    kind, v = DEFAULT_STAGE_EST.get(stage, ('plain', 60))
    h = perf.get(stage) or {}
    if kind == 'gb':
        return (h.get('spg') or v) * max(input_bytes / 1e9, 0.05)
    return h.get('plain') or v


def _safe_sample_name(name):
    s = re.sub(r'[^A-Za-z0-9_\-.]+', '_', str(name)).strip('._-')
    return s[:50] or 'sample'


def _load_done_summary(sample_dir, stage_dir, fname='summary.json'):
    """阶段未参与本次运行时，读取已完成的上一轮产物摘要（管道衔接）。"""
    p = os.path.join(sample_dir, stage_dir, fname)
    if not os.path.isfile(p):
        return None
    try:
        with safe_open(p) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def subsample_fastq(r1, r2, out_dir, n_pairs, logger=None, progress=None):
    """取前 n_pairs 对 reads 写子样本。返回 (sub_r1, sub_r2|None)。"""
    sub1 = os.path.join(out_dir, 'sub_R1.fastq.gz')
    n = 0
    if r2:
        sub2 = os.path.join(out_dir, 'sub_R2.fastq.gz')
        with safe_open(sub1, 'wt') as w1, safe_open(sub2, 'wt') as w2:
            for rec1, rec2 in zip(iter_fastq_records(r1), iter_fastq_records(r2)):
                w1.write(''.join(rec1))
                w2.write(''.join(rec2))
                n += 1
                if n >= n_pairs:
                    break
                if progress and n % 100_000 == 0:
                    progress(min(n / n_pairs, 0.98),
                             f'已截取 {n:,} / {n_pairs:,} 对')
        if logger:
            logger.log(f"子采样完成: {n:,} 对 reads")
        return sub1, sub2
    with safe_open(sub1, 'wt') as w1:
        for rec1 in iter_fastq_records(r1):
            w1.write(''.join(rec1))
            n += 1
            if n >= n_pairs:
                break
            if progress and n % 100_000 == 0:
                progress(min(n / n_pairs, 0.98), f'已截取 {n:,} / {n_pairs:,} 条')
    if logger:
        logger.log(f"子采样完成: {n:,} 条 reads")
    return sub1, None


class _RunState:
    """全局加权进度 + ETA（剩余估时 = 历史耗时模型 × 剩余权重）。"""

    def __init__(self, stages, input_bytes, progress, logger):
        self.perf = _load_stage_perf()
        self.stages = stages
        self.input_bytes = input_bytes
        self._emit = progress
        self._logger = logger
        self.cur = None
        self.cur_est = 0.0
        self.done_est = 0.0
        self.t0 = time.time()
        self.total_est = sum(self.est(s) for s in stages) or 1.0

    def est(self, stage):
        return _est_stage_sec(self.perf, stage,
                              self.input_bytes if stage in _GB_STAGES else 0)

    def begin(self, stage):
        self.cur = stage
        self.cur_est = max(self.est(stage), 1.0)
        self.t0 = time.time()
        if self._logger:
            self._logger.log(
                f"──── {STAGE_NAMES.get(stage, stage)} 开始"
                f"（预估 ~{fmt_eta(self.cur_est)}）────", "STAGE")
        return self.t0

    def end(self, stage, input_bytes=0):
        dt = time.time() - self.t0
        _update_stage_perf(self.perf, stage, dt,
                           input_bytes if stage in _GB_STAGES else 0)
        self.done_est += self.cur_est
        self.cur_est = 0.0
        self.cur = None
        if self._logger:
            self._logger.log(
                f"──── {STAGE_NAMES.get(stage, stage)} 完成"
                f"（实际 {fmt_eta(dt)}）────", "STAGE")
        return dt

    def prog(self, stage, pct, msg):
        pct = min(max(float(pct), 0.0), 1.0)
        done = self.done_est + (self.cur_est * pct if stage == self.cur else 0.0)
        remain = max(self.total_est - done, 0.0)
        if self._logger:
            self._logger.log(f"[{STAGE_NAMES.get(stage, stage)} {pct:.0%}] {msg}")
        if self._emit:
            g = min(max(done / self.total_est, 0.0), 1.0)   # 估时偏差防溢出
            try:
                self._emit(stage, g, msg, remain)
            except TypeError:
                try:
                    self._emit(stage, g, msg)
                except Exception:
                    pass
            except Exception:
                pass


def _fsize(p):
    try:
        return os.path.getsize(str(p)) if p else 0
    except OSError:
        return 0


class _Ctx:
    """阶段执行上下文：跨阶段流动的产物 + 全部运行参数。

    阶段函数只通过 C 读写：输入取上一阶段写入的字段（C.cur_r1 /
    C.host_stats / C.vs …），产物写回同名字段，从而与执行顺序解耦。
    新增/重排阶段只需改 STAGE_REGISTRY，不改执行器。
    """
    pass


def _adopt_qc(C):
    """质控产物优先作为后续输入（幂等：fastp 产物存在才切换）。"""
    qc = C.qc
    if qc and qc.get('fastp_r1'):
        prep = os.path.join(C.sample_dir, '00_prep')
        q1 = os.path.join(prep, qc['fastp_r1'])
        if os.path.isfile(q1):
            C.cur_r1 = q1
            q2 = os.path.join(prep, qc.get('fastp_r2') or '')
            C.cur_r2 = q2 if q2 and os.path.isfile(q2) else None


def _stage_subsample(C):
    """0. 子采样（确定性：前 N 对 reads → 断点复用，避免使 fq2fa
    转换产物的新鲜度失效）。"""
    if not (C.subsample and C.subsample > 0):
        return
    prep = os.path.join(C.sample_dir, '00_prep')
    os.makedirs(prep, exist_ok=True)
    sub1 = os.path.join(prep, 'sub_R1.fastq.gz')
    sub2 = os.path.join(prep, 'sub_R2.fastq.gz') if C.r2 else None
    reuse = (not C.force and is_step_done(prep, 'subsample')
             and os.path.isfile(sub1) and (not C.r2 or os.path.isfile(sub2)))
    if reuse:
        try:
            with safe_open(os.path.join(prep, 'subsample.json')) as f:
                prev_n = json.load(f).get('n_pairs')
            reuse = (prev_n == C.subsample)
        except (OSError, ValueError):
            reuse = False
    C.st.begin('subsample')
    if reuse:
        if C.logger:
            C.logger.log(f"子采样已完成（{C.subsample:,} 对），复用 00_prep/sub_*")
    else:
        C.st.prog('subsample', 0.05, f'截取前 {C.subsample:,} 对 reads')
        C.cur_r1, C.cur_r2 = subsample_fastq(
            C.cur_r1, C.cur_r2, prep, C.subsample, C.logger,
            progress=lambda p, m: C.st.prog('subsample', p, m))
        C.st.prog('subsample', 1.0, '子采样完成')
        with safe_open(os.path.join(prep, 'subsample.json'), 'wt') as f:
            json.dump({'n_pairs': C.subsample}, f)
        mark_step_done(prep, 'subsample')
    C.st.end('subsample')


def _stage_fastp(C):
    """⓪ Fastp 质控（可选；工具缺失时上游已剔除本阶段）。"""
    st, logger = C.st, C.logger
    st.begin('fastp')
    log_res_plan(logger, '⓪ Fastp 质控', threads=C.threads,
                 mem_gb=max(C.threads or 4, 4) * 0.3,
                 disk_gb=(_fsize(C.cur_r1) + _fsize(C.cur_r2)) * 5.0 / 1e9,
                 note="读写各一份 gzip 输出（crabz 加速时更快）")
    st.prog('fastp', 0.05, 'fastp 去接头/低质量裁剪')
    from .preprocess import run_fastp
    C.qc = run_fastp(C.sample_dir, C.cur_r1, C.cur_r2, threads=C.threads,
                     logger=logger, force=C.force, dedup=C.fastp_dedup)
    st.prog('fastp', 1.0, '质控完成' if C.qc else 'fastp 未安装，已跳过')
    st.end('fastp')


def _stage_fq2fa(C):
    """⓪b FASTQ→FASTA 预转换（分类加速；产物复用给 ①宿主去除）。"""
    st = C.st
    _adopt_qc(C)
    st.begin('fq2fa')
    from .preprocess import convert_fq2fa
    fq2 = convert_fq2fa(C.sample_dir, C.cur_r1, C.cur_r2, threads=C.threads,
                        logger=C.logger, force=C.force,
                        progress=lambda p, m: st.prog('fq2fa', p, m))
    st.prog('fq2fa', 1.0,
            '输入已是 FASTA，无需转换' if fq2 and fq2.get('skipped')
            else '转换完成')
    st.end('fq2fa')


def _stage_host(C):
    """① 宿主去除（分类输入优先复用 ⓪b 转换产物）。"""
    st, logger = C.st, C.logger
    _adopt_qc(C)
    conv1, conv2 = _fresh_converted(C.sample_dir, C.cur_r1, C.cur_r2)
    st.begin('host')
    C.host_stats = remove_host_logged(
        C.sample_dir, C.cur_r1, C.cur_r2, C.db_host, threads=C.threads,
        confidence=C.confidence, logger=logger, force=C.force,
        chunk_dir=C.chunk_dir, classify_r1=conv1, classify_r2=conv2,
        allow_convert=C.do_fq2fa or bool(conv1),
        progress=lambda p, m: st.prog('host', p, m))
    st.prog('host', 1.0, f"宿主占比 {C.host_stats['host_ratio'] * 100:.2f}%")
    st.end('host')


def _stage_virus(C):
    """② 病毒筛查（输入：去宿主后的 kept 或原始），并提取病毒 reads。"""
    st, logger = C.st, C.logger
    _adopt_qc(C)
    v_in1, v_in2 = C.cur_r1, C.cur_r2
    if C.host_stats:
        hr = os.path.join(C.sample_dir, '01_host_removal')
        v_in1 = os.path.join(hr, C.host_stats.get('kept_r1') or 'kept_R1.fastq.gz')
        kr2 = C.host_stats.get('kept_r2')
        v_in2 = os.path.join(hr, kr2) if kr2 else None
    st.begin('virus')
    from .virus_screen import screen_virus
    C.vs = screen_virus(C.sample_dir, v_in1, v_in2, C.db_virus,
                        threads=C.threads, confidence=C.confidence,
                        logger=logger, force=C.force, chunk_dir=C.chunk_dir,
                        allow_convert=C.do_fq2fa,
                        progress=lambda p, m: st.prog('virus', p, m))
    st.prog('virus', 1.0, f"检出物种 {len(C.vs['species_detected'])} 个")
    st.end('virus')


def _assembly_inputs(C):
    """③ 组装输入解析：病毒 reads 优先 → 回退去宿主 / 原始 reads。"""
    v_in1, v_in2 = C.cur_r1, C.cur_r2
    if C.host_stats:
        hr = os.path.join(C.sample_dir, '01_host_removal')
        v_in1 = os.path.join(hr, C.host_stats.get('kept_r1') or 'kept_R1.fastq.gz')
        kr2 = C.host_stats.get('kept_r2')
        v_in2 = os.path.join(hr, kr2) if kr2 else None

    a_in1, a_in2, a_src = v_in1, v_in2, 'kept' if C.host_stats else 'raw'
    if C.assembly_input == 'raw':
        a_in1, a_in2, a_src = C.cur_r1, C.cur_r2, 'raw'
    elif C.assembly_input == 'virus':
        vp1 = os.path.join(C.sample_dir, '02_virus_screen', 'viral_R1.fastq.gz')
        vp2 = os.path.join(C.sample_dir, '02_virus_screen', 'viral_R2.fastq.gz')
        n_viral = (C.vs or {}).get('viral_pairs', 0)
        if os.path.isfile(vp1) and n_viral >= C.min_viral_pairs:
            a_in1, a_in2, a_src = vp1, (vp2 if os.path.isfile(vp2) else None), 'virus'
        else:
            C.st.prog('assembly', 0.01,
                      f"病毒 reads 仅 {n_viral:,} 对(<{C.min_viral_pairs})，"
                      f"回退{'去宿主' if C.host_stats else '原始'} reads 组装")
    return a_in1, a_in2, a_src


def _stage_assembly(C):
    """③ 组装 + contig 分类提取。"""
    st, logger = C.st, C.logger
    a_in1, a_in2, a_src = _assembly_inputs(C)
    st.begin('assembly')
    log_res_plan(logger, '③ SPAdes 组装', threads=C.threads,
                 mem_gb=C.memory_gb,
                 disk_gb=max((_fsize(a_in1) + _fsize(a_in2)) * 5.5 * 4,
                             10e8) / 1e9,
                 note=f"内存上限={C.memory_gb}GB（-m 参数）；中间图目录≈输入"
                      f"解压量数倍；中文路径自动经 ASCII 中转")
    st.prog('assembly', 0.05, f'SPAdes 组装 ({C.assembly_mode}, 输入={a_src})')
    from .assembly import assemble_and_classify
    C.asm = assemble_and_classify(C.sample_dir, a_in1, a_in2, C.db_virus,
                                  mode=C.assembly_mode, threads=C.threads,
                                  memory_gb=C.memory_gb,
                                  min_contig_len=C.min_contig_len,
                                  logger=logger, force=C.force,
                                  chunk_dir=C.chunk_dir,
                                  progress=lambda p, m:
                                  st.prog('assembly', p, m))
    C.asm['assembly_input'] = a_src
    st.prog('assembly', 1.0, f"病毒 contigs {len(C.asm['viral_contigs'])} 条")
    st.end('assembly')


def _stage_verify(C):
    """③b 候选序列验证：kunpeng 初筛结果的证据判定。

    输入 03_assembly/viral_contigs.fasta + virus_classification.tsv，
    双路证据（blastx 病毒蛋白库 + CDD 结构域）→ 四档判定。
    unclassified 保留为新病毒候选池（不删除），但不进下游注释。
    """
    st, logger = C.st, C.logger
    st.begin('verify')
    a_dir = os.path.join(C.sample_dir, '03_assembly')
    vfa = os.path.join(a_dir, 'viral_contigs.fasta')
    if not os.path.isfile(vfa):
        raise RuntimeError('03_assembly 无 viral_contigs.fasta，请先完成 ③组装阶段')
    log_res_plan(logger, '③b 候选序列验证', threads=C.threads, mem_gb=2.0,
                 disk_gb=0.5,
                 note="DIAMOND blastx vs 病毒蛋白库 + mmseqs2 vs CDD；"
                      "耗时随候选数线性增长")
    st.prog('verify', 0.05, '宿主筛选 → 长度分流 → 双路证据')
    from .verify import verify as run_verify
    vs = run_verify(C.sample_dir, vfa, host=C.verify_host,
                    methods=C.verify_methods, combine=C.verify_combine,
                    threads=C.threads, logger=logger,
                    progress=lambda s_, f_, m: st.prog('verify', f_, m),
                    out_subdir='03b_verify', assembly_dir=a_dir)
    C.verify = vs
    calls = vs.get('calls') or {}
    st.prog('verify', 1.0,
            f"验证 {vs.get('n_calls', 0)} 条："
            + '，'.join(f'{k} {v}' for k, v in calls.items()))
    st.end('verify')


def _stage_consensus(C):
    """③c 共识序列与变异：reads 回贴候选 contig → 共识序列 + 变异谱。

    参考 = 03_assembly/viral_contigs.fasta；
    reads = 01_host_removal（去宿主后全量）——未经过病毒相似度筛选，
    变异与准种多样性检测无偏（02_virus_screen 会低估变异）。
    """
    st, logger = C.st, C.logger
    st.begin('consensus')
    a_dir = os.path.join(C.sample_dir, '03_assembly')
    vfa = os.path.join(a_dir, 'viral_contigs.fasta')
    if not os.path.isfile(vfa):
        raise RuntimeError('03_assembly 无 viral_contigs.fasta，请先完成 ③组装阶段')
    from .consensus import consensus_and_variants, find_read_pairs
    reads = find_read_pairs(os.path.join(C.sample_dir, '01_host_removal'))
    if not reads:
        raise RuntimeError('01_host_removal 无配对 reads，请先完成 ①宿主去除阶段')
    log_res_plan(logger, '③c 共识序列与变异', threads=C.threads, mem_gb=2.0,
                 disk_gb=1.0,
                 note="minimap2 -ax sr 回贴 + 逐位统计；"
                      "内存只与参考长度相关，与测序深度无关")
    st.prog('consensus', 0.05, 'minimap2 回贴 reads 到候选 contig')
    cs = consensus_and_variants(
        C.sample_dir, vfa, reads, out_subdir='03c_consensus',
        threads=C.threads, logger=logger,
        progress=lambda s_, f_, m: st.prog('consensus', f_, m))
    C.consensus = cs
    st.prog('consensus', 1.0,
            f"共识 {cs.get('n_refs', 0)} 条，比对 reads "
            f"{cs.get('n_mapped_reads', 0)}，变异位点 {cs.get('n_variants', 0)}"
            f"（SNV {cs.get('n_snv', 0)} / iSNV {cs.get('n_isnv', 0)}），"
            f"判定存在 {len(cs.get('present') or [])} 条")
    st.end('consensus')


def _stage_hostana(C):
    """④ 宿主预测（ICTV 级联 + NCBI 宿主元数据交叉）。"""
    st, logger = C.st, C.logger
    st.begin('hostana')
    log_res_plan(logger, '④ ICTV 宿主预测', threads=C.threads, mem_gb=1.0,
                 disk_gb=0.1, note="查表+BLAST 回退，资源占用小")
    st.prog('hostana', 0.1, 'ICTV 宿主概率级联预测')
    from .host_analysis import predict_hosts
    ha = predict_hosts(C.sample_dir, threads=C.threads, logger=logger,
                       force=C.force)
    st.prog('hostana', 1.0,
            f"宿主判定 {ha['n_contigs']} 条，"
            f"{len(ha['categories'])} 类")
    st.end('hostana')


def _stage_orf(C):
    """⑥ ORF 预测。"""
    st, logger = C.st, C.logger
    st.begin('orf')
    log_res_plan(logger, '⑥ ORF 预测', threads=C.threads, mem_gb=2.0,
                 disk_gb=0.5, note="orfipy 多进程 + pyrodigal")
    st.prog('orf', 0.1, 'orfipy + pyrodigal 基因预测')
    from .orf import predict_orfs
    predict_orfs(C.sample_dir, min_aa=C.min_orf_aa, threads=C.threads,
                 logger=logger, force=C.force,
                 progress=lambda p, m: st.prog('orf', p, m))
    st.prog('orf', 1.0, 'ORF 预测完成')
    st.end('orf')


def _stage_orfa(C):
    """⑥b ORF 功能注释（DIAMOND/MMseqs2/blastp 对照 RefSeq 病毒蛋白库）。"""
    st, logger = C.st, C.logger
    st.begin('orfa')
    log_res_plan(logger, '⑥b ORF 功能注释', threads=C.threads, mem_gb=2.0,
                 disk_gb=2.5,
                 note="对照 RefSeq 病毒蛋白库（首次自动下载 ~107MB，"
                      "之后全缓存）")
    st.prog('orfa', 0.05, '参考库与搜索引擎准备')
    from .orf_annot import run_orf_annotation
    oa = run_orf_annotation(C.sample_dir, threads=C.threads, logger=logger,
                            force=C.force,
                            progress=lambda p, m: st.prog('orfa', p, m))
    st.prog('orfa', 1.0,
            f"功能注释 {oa['n_annotated']}/{oa['n_orfs']} 个 ORF，"
            f"{oa['n_families']} 科")
    st.end('orfa')


def _stage_phylo(C):
    """⑦ 进化树 + SDT。"""
    st, logger = C.st, C.logger
    st.begin('phylo')
    log_res_plan(logger, '⑦ 进化树与 SDT', threads=C.threads, mem_gb=2.0,
                 disk_gb=0.5,
                 note=f"MAFFT+{C.tree_tool}；仅比对近缘参考小片段")
    st.prog('phylo', 0.1, '近缘参考比对建树 + SDT 矩阵')
    from .phylo import build_phylo, resolve_ncbi_refs
    extra_refs = None
    if C.ncbi_refs:
        extra_refs = resolve_ncbi_refs(
            [n.strip() for n in str(C.ncbi_refs).split(',') if n.strip()])
        if logger:
            logger.log(f"额外参考: {len(extra_refs)} 个 NCBI 集合")
    ps = build_phylo(C.sample_dir, top_n_refs=C.top_n_refs,
                     tree_tool=C.tree_tool,
                     do_trim=C.do_trim,
                     threads=C.threads, logger=logger, force=C.force,
                     extra_refs=extra_refs, sampling=C.tree_sampling,
                     db_virus=C.db_virus)
    st.prog('phylo', 1.0, f"分析组 {len(ps.get('groups', []))} 个")
    st.end('phylo')


def _stage_primer(C):
    """⑧ 引物设计。"""
    st = C.st
    st.begin('primer')
    st.prog('primer', 0.1, f'primer3 引物设计 ({C.primer_mode})')
    from .primer import design_primers
    pr = design_primers(C.sample_dir, mode=C.primer_mode, logger=C.logger,
                        force=C.force, do_specificity=C.do_specificity)
    st.prog('primer', 1.0, f"引物 {pr.get('n_primers', 0)} 对")
    st.end('primer')


def _stage_gbdraw(C):
    """⑨ 基因组图（gbdraw 首选；不可用时 dna_features_viewer 顶上）。"""
    st, logger = C.st, C.logger
    st.begin('gbdraw')
    engine = C.plot_engine if C.plot_engine in ('gbdraw', 'dfv') else \
        ('gbdraw' if _gbdraw_ok() else 'dfv')
    log_res_plan(logger, f'⑨ 基因组图 ({engine})', threads=1, mem_gb=0.6,
                 disk_gb=0.05, note="每条病毒 contig 出圈图+线图 SVG")
    if engine == 'gbdraw':
        from .gbdraw_plot import run_genome_plots
        plots = run_genome_plots(C.sample_dir, logger=logger, force=C.force,
                                 max_plots=C.gbdraw_max,
                                 fasta_in=C.gbdraw_fasta, ann_in=C.gbdraw_ann,
                                 progress=lambda p, m: st.prog('gbdraw', p, m))
    else:
        from .dfv_plot import run_dfv_plots
        plots = run_dfv_plots(C.sample_dir, logger=logger, force=C.force,
                              max_plots=C.gbdraw_max,
                              fasta_in=C.gbdraw_fasta, ann_in=C.gbdraw_ann,
                              progress=lambda p, m: st.prog('gbdraw', p, m))
    st.prog('gbdraw', 1.0, f"基因组图 {len(plots.get('plots', []))} 张（{engine}）")
    st.end('gbdraw')


def _stage_report(C):
    """⑩ 可视化报告。"""
    st = C.st
    st.begin('report')
    st.prog('report', 0.2, '生成可视化报告')
    from .viz import build_report
    rpt = build_report(C.sample_dir, logger=C.logger)
    st.prog('report', 1.0, f"报告: {rpt}")
    st.end('report')


# ------------------------------------------------------------------
# 阶段注册表：key → {deps, fn}。新增阶段登记一行即可接入；
# 调整执行顺序只需改本表（执行器按 deps 做稳定拓扑排序）。
# deps 只约束"已选阶段之间"的先后；未选中的依赖不强制补跑，
# 由阶段函数读取上一轮产物衔接（断点续跑语义，与旧实现一致）。
# ------------------------------------------------------------------
STAGE_REGISTRY = {
    'subsample': {'deps': [],                    'fn': _stage_subsample},
    'fastp':     {'deps': [],                    'fn': _stage_fastp},
    'fq2fa':     {'deps': [],                    'fn': _stage_fq2fa},
    'host':      {'deps': [],                    'fn': _stage_host},
    'virus':     {'deps': [],                    'fn': _stage_virus},
    'assembly':  {'deps': [],                    'fn': _stage_assembly},
    'verify':    {'deps': ['assembly'],          'fn': _stage_verify},
    'consensus': {'deps': ['assembly'],          'fn': _stage_consensus},
    'hostana':   {'deps': ['assembly'],          'fn': _stage_hostana},
    'orf':       {'deps': ['assembly'],          'fn': _stage_orf},
    'orfa':      {'deps': ['orf'],               'fn': _stage_orfa},
    'phylo':     {'deps': ['assembly'],          'fn': _stage_phylo},
    'primer':    {'deps': ['assembly'],          'fn': _stage_primer},
    'gbdraw':    {'deps': ['assembly'],          'fn': _stage_gbdraw},
    'report':    {'deps': ['virus'],              'fn': _stage_report},
}

# 兼容导出：依赖解锁表（管道页卡片）沿用旧名字
STAGE_DEPS = {k: list(v['deps']) for k, v in STAGE_REGISTRY.items()}


def _stage_execution_order(selected, logger=None):
    """已选阶段的稳定拓扑排序：保证依赖先跑，其余保持调用方顺序。

    未知阶段名会被忽略，但一定要留日志——否则用户拼错阶段名
    （如 README 曾经写过的 virome）会得到「什么都没跑但也没报错」。
    """
    unknown = [s for s in selected if s not in STAGE_REGISTRY]
    if unknown and logger:
        logger(f"忽略未知阶段名：{', '.join(unknown)}（"
               f"可用阶段：{', '.join(STAGE_ORDER)}）")
    sel = [s for s in selected if s in STAGE_REGISTRY]
    sel_set = set(sel)
    order, seen = [], set()

    def visit(k):
        if k in seen:
            return
        seen.add(k)
        for d in STAGE_REGISTRY[k]['deps']:
            if d in sel_set:
                visit(d)
        order.append(k)

    for k in sel:
        visit(k)
    return order


def run_analysis(sample, r1, r2, stages, db_host=None, db_virus=None,
                 threads=None, confidence=0.0, assembly_mode='metaviral',
                 memory_gb=64, subsample=0, min_contig_len=200, top_n_refs=10,
                 tree_tool='fasttree', tree_sampling='blast', ncbi_refs=None,
                 primer_mode='conserved', do_trim=True,
                 do_specificity=False, min_orf_aa=100, force=False,
                 assembly_input='virus', min_viral_pairs=500,
                 chunk_dir=None, fastp_dedup=False, do_fq2fa=True,
                 gbdraw_max=12, gbdraw_fasta=None, gbdraw_ann=None,
                 plot_engine='auto',
                 verify_host='all', verify_methods=('blastx', 'cdd'),
                 verify_combine='union',
                 logger=None, progress=None):
    """全流程编排。progress(stage, pct 0~1, msg[, eta_seconds]) 供 GUI 更新。

    执行体来自 STAGE_REGISTRY：按已选阶段的稳定拓扑排序逐个运行；
    未选中阶段的上一轮产物摘要照常衔接（断点续跑）。
    """
    stages = [s for s in STAGE_ORDER if s in stages]
    # 可选步骤在工具缺失/用户关闭时静默剔除
    if 'fastp' in stages:
        from .preprocess import fastp_available
        if not fastp_available():
            stages = [s for s in stages if s != 'fastp']
    if 'fq2fa' in stages and (not do_fq2fa or not _seqkit_ok()):
        stages = [s for s in stages if s != 'fq2fa']
    if 'gbdraw' in stages and not _plots_ok(plot_engine):
        stages = [s for s in stages if s != 'gbdraw']
    if 'orfa' in stages:
        from .orf_annot import annotation_engine
        if not annotation_engine():
            stages = [s for s in stages if s != 'orfa']
    if 'verify' in stages:
        from .verify import verify_engine
        _ve = verify_engine()
        if not (_ve.get('diamond') or _ve.get('mmseqs')):
            if logger:
                logger.log('未检测到 DIAMOND / MMseqs2，跳过 ③b 候选序列验证',
                           'WARN')
            stages = [s for s in stages if s != 'verify']

    sample = _safe_sample_name(sample)
    sample_dir = check_path(os.path.join(DIRS['results'], sample),
                            must_exist=False, in_platform=True)
    os.makedirs(sample_dir, exist_ok=True)
    st = _RunState(stages, _fsize(r1) + _fsize(r2), progress, logger)

    C = _Ctx()
    C.sample_dir = sample_dir
    C.st = st
    C.logger = logger
    C.force = force
    C.threads = threads
    C.chunk_dir = chunk_dir
    C.r1, C.r2 = r1, r2
    C.db_host, C.db_virus = db_host, db_virus
    C.confidence, C.assembly_mode = confidence, assembly_mode
    C.memory_gb, C.subsample = memory_gb, subsample
    C.min_contig_len, C.top_n_refs = min_contig_len, top_n_refs
    C.tree_tool, C.tree_sampling = tree_tool, tree_sampling
    C.ncbi_refs = ncbi_refs
    C.primer_mode, C.do_trim = primer_mode, do_trim
    C.do_specificity, C.min_orf_aa = do_specificity, min_orf_aa
    C.assembly_input, C.min_viral_pairs = assembly_input, min_viral_pairs
    C.fastp_dedup, C.do_fq2fa = fastp_dedup, do_fq2fa
    C.gbdraw_max, C.gbdraw_fasta = gbdraw_max, gbdraw_fasta
    C.gbdraw_ann, C.plot_engine = gbdraw_ann, plot_engine
    # ③b 候选序列验证参数（宿主类群 / 证据方法 / 并集-交集）
    C.verify_host = verify_host
    C.verify_methods = tuple(verify_methods or ('blastx', 'cdd'))
    C.verify_combine = verify_combine
    # 跨阶段产物
    C.cur_r1 = check_path(r1, must_exist=True)
    C.cur_r2 = check_path(r2, must_exist=True) if r2 else None
    C.qc = None
    C.host_stats = None
    C.vs = None
    C.asm = None

    selected = set(stages)
    # 未选中阶段的断点衔接：读取上一轮已完成产物摘要（与旧实现一致）
    if 'fastp' not in selected:
        prep = os.path.join(sample_dir, '00_prep')
        if os.path.isfile(os.path.join(prep, 'fastp_R1.fastq.gz')):
            from .preprocess import _summary as _fastp_summary
            C.qc = _fastp_summary(prep)
    if 'host' not in selected:
        C.host_stats = _load_done_summary(sample_dir, '01_host_removal',
                                          'stats.json')
    if 'virus' not in selected:
        C.vs = _load_done_summary(sample_dir, '02_virus_screen')

    # 项目清单：开始运行即登记本次输入/参数（复现依据），逐阶段追加进度
    _run_t0 = time.time()
    update_project_manifest(
        sample_dir, sample=sample, r1=r1, r2=r2,
        params=_manifest_params({
            'assembly_mode': assembly_mode, 'assembly_input': assembly_input,
            'confidence': confidence, 'threads': threads,
            'memory_gb': memory_gb, 'subsample': subsample,
            'min_contig_len': min_contig_len, 'top_n_refs': top_n_refs,
            'tree_tool': tree_tool, 'tree_sampling': tree_sampling,
            'primer_mode': primer_mode, 'do_trim': do_trim,
            'do_specificity': do_specificity, 'min_orf_aa': min_orf_aa,
            'min_viral_pairs': min_viral_pairs, 'fastp_dedup': fastp_dedup,
            'do_fq2fa': do_fq2fa, 'gbdraw_max': gbdraw_max,
            'plot_engine': plot_engine, 'db_host': db_host,
            'db_virus': db_virus, 'ncbi_refs': ncbi_refs,
            'chunk_dir': chunk_dir,
        }),
        status='running')

    def _manifest_sync(stages_done, status, duration_s=None):
        """把已完成阶段同步到清单，只记真正出结果的。

        「跑过但无结果」（skipped：无病毒 contigs / 无引物 / 无基因组图）
        不进 stages_done，与 pipeline_overview 卡片状态保持同一口径。
        real_done 已覆盖历史与本次的实际产物，因此不需再并入
        「本次跑过的阶段名单」——那会把 skipped 又放回来。
        阶段函数写完 summary 才返回，磁盘顺序保证本次结果已可读。
        """
        try:
            real_done = [s['stage'] for s in pipeline_overview(sample_dir)['stages']
                         if s['status'] == 'done']
        except Exception:
            # overview 不可用时保守回退：保留旧值 + 本次跑过的
            real_done = list(load_project_manifest(sample_dir).get('stages_done') or [])
            real_done += [s for s in stages_done if s not in real_done]
        update_project_manifest(sample_dir, stages_done=real_done,
                                status=status, duration_s=duration_s,
                                replace_stages=True)

    _finished = []
    try:
        for key in _stage_execution_order(stages, logger=logger):
            STAGE_REGISTRY[key]['fn'](C)
            _finished.append(key)
            _manifest_sync(_finished, 'running')
    except Exception:
        # 失败也要留下痕迹：记下已跑完的阶段与实际耗时
        _manifest_sync(_finished, 'failed',
                       duration_s=time.time() - _run_t0)
        raise
    _manifest_sync(_finished, 'ok', duration_s=time.time() - _run_t0)

    return sample_dir


def remove_host_logged(sample_dir, r1, r2, db_host, **kw):
    """①宿主去除 + 资源预估日志（磁盘按经验系数，内存细节由 classify 记录）。"""
    logger = kw.get('logger')
    from .host_removal import remove_host
    log_res_plan(logger, '① kunpeng 宿主分类',
                 threads=kw.get('threads'),
                 disk_gb=max(_fsize(r1) + _fsize(r2), 1) * 16 / 1e9,
                 note="chunk 中间盘≈输入gz×16（运行结束记录实测值）；"
                      "内存=库hash表+resolve(reads×190B)，详见分类开始行；"
                      "内存不足自动分片")
    return remove_host(sample_dir, r1, r2, db_host, **kw)


def _seqkit_ok():
    try:
        from .preprocess import seqkit_available
        return seqkit_available()
    except Exception:
        return False


def _gbdraw_ok():
    try:
        from .config import get_config
        get_config().tool('gbdraw')
        return True
    except Exception:
        return False


def _dfv_ok():
    try:
        from .dfv_plot import dfv_available
        return dfv_available()
    except Exception:
        return False


def _plots_ok(engine):
    """⑨基因组图引擎可用性：auto=任一引擎即可，指定引擎则须该引擎可用。"""
    if engine == 'gbdraw':
        return _gbdraw_ok()
    if engine == 'dfv':
        return _dfv_ok()
    return _gbdraw_ok() or _dfv_ok()


def _fresh_converted(sample_dir, cur_r1, cur_r2):
    """可安全复用的预处理 FASTA（conv_R1/R2.fa.gz）。

    转换产物 mtime 不旧于当前输入才复用（fastp/输入更新后旧转换文件失效）。
    """
    prep = check_path(os.path.join(sample_dir, '00_prep'),
                      must_exist=False, in_platform=True)
    c1 = os.path.join(prep, 'conv_R1.fa.gz')
    c2 = os.path.join(prep, 'conv_R2.fa.gz')
    out1 = out2 = None

    def _fresh(conv, src):
        if conv and os.path.isfile(conv) and src and os.path.isfile(str(src)):
            try:
                return os.path.getmtime(conv) >= os.path.getmtime(str(src)) - 5
            except OSError:
                return False
        return False

    if _fresh(c1, cur_r1):
        out1 = c1
        if _fresh(c2, cur_r2):
            out2 = c2
    return out1, out2


def run_report_only(sample, logger=None, force=False):
    sample = _safe_sample_name(sample)
    sample_dir = check_path(os.path.join(DIRS['results'], sample),
                            must_exist=True, in_platform=True)
    from .viz import build_report
    return build_report(sample_dir, logger=logger)


# ------------------------------------------------------------------
# 管道视图（GUI 模块化逐级运行）
# ------------------------------------------------------------------
PIPELINE_STAGES = [
    # (阶段键, 展示名, 目录, 摘要文件)
    ('subsample', '预处理·子采样', '00_prep', 'subsample.json'),
    ('fastp', '⓪ Fastp 质控', '00_prep', 'fastp_report.json'),
    ('fq2fa', '⓪b 序列转换 (FASTQ→FASTA)', '00_prep', 'fq2fa.json'),
    ('host', '① 宿主去除', '01_host_removal', 'stats.json'),
    ('virus', '② 病毒筛查与提取', '02_virus_screen', 'summary.json'),
    ('assembly', '③ 组装·分类·提取', '03_assembly', 'summary.json'),
    ('verify', '③b 候选序列验证', '03b_verify', 'summary.json'),
    ('hostana', '④ 宿主预测 (ICTV)', '08_host_analysis', 'summary.json'),
    ('orf', '⑥ 编码区预测 ORF', '04_orf', 'summary.json'),
    ('orfa', '⑥b ORF 功能注释', '04b_orf_annot', 'summary.json'),
    ('phylo', '⑦ 进化树与 SDT', '05_phylo', 'summary.json'),
    ('primer', '⑧ 引物设计', '06_primer', 'summary.json'),
    ('gbdraw', '⑨ 基因组图 (gbdraw/DFV)', '09_genome_plots', 'summary.json'),
    ('report', '⑩ 可视化报告', '07_report', None),
]
PIPELINE_STAGE_NAMES_EN = {
    'subsample': 'Preprocess · subsample',
    'fastp': '⓪ Fastp QC', 'fq2fa': '⓪b FASTQ→FASTA',
    'host': '① Host removal', 'virus': '② Virus screening & extraction',
    'assembly': '③ Assembly & extraction',
    'verify': '③b Candidate verify',
    'hostana': '④ Host prediction (ICTV)',
    'orf': '⑥ ORF prediction',
    'orfa': '⑥b ORF annotation', 'phylo': '⑦ Phylogeny & SDT',
    'primer': '⑧ Primer design', 'gbdraw': '⑨ Genome plots (gbdraw/DFV)',
    'report': '⑩ Visual report',
}
# 工具缺失时的提示（中/英）
_MISSING_TOOL_MSG = {
    'fastp': ('未检测到 fastp.exe（可选步骤，已自动跳过）',
              'fastp.exe not found (optional; auto-skipped)'),
    'fq2fa': ('未检测到 seqkit.exe（可选步骤；分类将直接用 FASTQ）',
              'seqkit.exe not found (optional; classify on FASTQ)'),
    'gbdraw': ('未检测到基因组图引擎（任一即可：pip install git+'
               'https://github.com/satoshikawato/gbdraw.git 或 '
               'pip install dna_features_viewer）',
               'No genome-plot engine found (install gbdraw or '
               'dna_features_viewer)'),
    'orfa': ('未检测到 DIAMOND / MMseqs2 / blastp（任一即可）',
             'No DIAMOND / MMseqs2 / blastp found (any one suffices)'),
    'verify': ('未检测到 DIAMOND / MMseqs2（验证需 blastx 与 CDD）',
               'No DIAMOND / MMseqs2 found (verify needs blastx + CDD)'),
}
# 严格上游依赖（未列出的阶段可直接用原始/默认输入运行）
# 单一数据源：直接取自 STAGE_REGISTRY 的 deps，避免两处硬编码走神。
# 注：subsample 与其余预处理阶段无硬依赖（都用原始/上轮输入）。
STAGE_DEPS = {k: list(v['deps']) for k, v in STAGE_REGISTRY.items()}
# 卡片摘要取值函数
def _sum_subsample(s):
    n = s.get('n_pairs')
    return f"子采样 {n:,} 对 reads" if n else '子采样完成'


def _sum_subsample_en(s):
    n = s.get('n_pairs')
    return f"Subsampled {n:,} read pairs" if n else 'Subsampled'


def _sum_fastp(s):
    parts = []
    if s.get('reads_before') and s.get('reads_after'):
        parts.append(f"{s['reads_before']:,}→{s['reads_after']:,} reads")
    if s.get('q30_rate') is not None:
        parts.append(f"Q30 {s['q30_rate']}%")
    return '，'.join(parts) or '质控完成'


def _sum_fq2fa(s):
    if s.get('skipped'):
        return '输入已是 FASTA，无需转换'
    files = s.get('files') or []
    return f"已转换 {len(files)} 个文件：" + '，'.join(files)


def _sum_host(s): return f"宿主占比 {s.get('host_ratio', 0) * 100:.2f}%，保留 {s.get('kept_pairs', 0):,} 对"
def _sum_virus(s): return (f"检出物种 {len(s.get('species_detected', []))} 个，"
                           f"提取病毒 reads {s.get('viral_pairs', 0):,} 对")
def _sum_asm(s): return (f"contigs {s.get('contigs', '?')} 条"
                         f"（病毒 {len(s.get('viral_contigs', []))} 条）")


def _sum_verify(s):
    calls = s.get('calls') or {}
    if not calls:
        return '验证完成，无判定记录'
    order = ['known', 'novel', 'domain_only', 'unclassified', 'viroid']
    parts = [f'{k} {calls[k]}' for k in order if k in calls]
    parts += [f'{k} {v}' for k, v in calls.items() if k not in order]
    return f"共 {s.get('n_calls', 0)} 条：" + '，'.join(parts)


def _sum_verify_en(s):
    calls = s.get('calls') or {}
    if not calls:
        return 'Verification finished, no calls'
    order = ['known', 'novel', 'domain_only', 'unclassified', 'viroid']
    parts = [f'{k} {calls[k]}' for k in order if k in calls]
    parts += [f'{k} {v}' for k, v in calls.items() if k not in order]
    return f"{s.get('n_calls', 0)} total: " + ', '.join(parts)


def _sum_orf(s):
    n = s.get('n_orfs') or (s.get('orfipy') or {}).get('n_orfs')
    return f"ORF {n if n is not None else '?'} 个"


def _sum_orfa(s):
    extra = ''
    if s.get('n_hmm_only') or s.get('n_cdd_only'):
        extra = (f"，HMM 兜底 {s.get('n_hmm_only', 0)} / "
                 f"CDD 兜底 {s.get('n_cdd_only', 0)}")
    return (f"有效注释 {s.get('n_annotated', 0)}/{s.get('n_orfs', 0)} 个，"
            f"{s.get('n_families', 0)} 病毒科（引擎 {s.get('engine', '?')}）"
            f"{extra}")

def _sum_phylo(s): return f"分析组 {len(s.get('groups', []))} 个"
def _sum_primer(s): return f"引物 {s.get('n_primers', 0)} 对"


def _sum_hostana(s):
    cats = s.get('categories', {})
    top = ', '.join(f'{k} {v}' for k, v in list(cats.items())[:3])
    return f"{s.get('n_contigs', 0)} 条 contigs：{top}"


def _sum_gbdraw(s):
    eng = s.get('engine', 'gbdraw')
    return (f"基因组图 {len(s.get('plots', []))} 张"
            f"（SVG，{eng.replace('dna_features_viewer', 'DFV')}，已嵌入报告）")

STAGE_SUMMARIES = {'subsample': _sum_subsample,
                   'fastp': _sum_fastp, 'fq2fa': _sum_fq2fa, 'host': _sum_host,
                   'virus': _sum_virus,
                   'assembly': _sum_asm, 'verify': _sum_verify,
                   'hostana': _sum_hostana,
                   'orf': _sum_orf, 'orfa': _sum_orfa,
                   'phylo': _sum_phylo, 'primer': _sum_primer,
                   'gbdraw': _sum_gbdraw}

# ---- 英文摘要（设置页切到 English 时管道卡片用） ----
def _sum_fastp_en(s):
    parts = []
    if s.get('reads_before') and s.get('reads_after'):
        parts.append(f"{s['reads_before']:,} -> {s['reads_after']:,} reads")
    if s.get('q30_rate') is not None:
        parts.append(f"Q30 {s['q30_rate']}%")
    return ', '.join(parts) or 'QC finished'


def _sum_fq2fa_en(s):
    if s.get('skipped'):
        return 'Input already FASTA; conversion skipped'
    files = s.get('files') or []
    return f"Converted {len(files)} file(s): " + ', '.join(files)


def _sum_host_en(s):
    return (f"Host fraction {s.get('host_ratio', 0) * 100:.2f}%, "
            f"{s.get('kept_pairs', 0):,} pairs kept")


def _sum_virus_en(s):
    return (f"{len(s.get('species_detected', []))} species detected, "
            f"{s.get('viral_pairs', 0):,} viral read pairs extracted")


def _sum_asm_en(s):
    return (f"{s.get('contigs', '?')} contigs "
            f"({len(s.get('viral_contigs', []))} viral)")


def _sum_orf_en(s):
    n = s.get('n_orfs') or (s.get('orfipy') or {}).get('n_orfs')
    return f"{n if n is not None else '?'} ORFs"


def _sum_orfa_en(s):
    extra = ''
    if s.get('n_hmm_only') or s.get('n_cdd_only'):
        extra = (f", HMM fallback {s.get('n_hmm_only', 0)} / "
                 f"CDD fallback {s.get('n_cdd_only', 0)}")
    return (f"{s.get('n_annotated', 0)}/{s.get('n_orfs', 0)} ORFs annotated, "
            f"{s.get('n_families', 0)} viral families "
            f"(engine {s.get('engine', '?')}){extra}")


def _sum_phylo_en(s): return f"{len(s.get('groups', []))} analysis group(s)"
def _sum_primer_en(s): return f"{s.get('n_primers', 0)} primer pair(s)"


def _sum_hostana_en(s):
    cats = s.get('categories', {})
    top = ', '.join(f'{k} {v}' for k, v in list(cats.items())[:3])
    return f"{s.get('n_contigs', 0)} contigs: {top}"


def _sum_gbdraw_en(s):
    eng = s.get('engine', 'gbdraw')
    return (f"{len(s.get('plots', []))} genome plot(s) "
            f"(SVG, {eng.replace('dna_features_viewer', 'DFV')}, in report)")

STAGE_SUMMARIES_EN = {
    'subsample': _sum_subsample_en,
    'fastp': _sum_fastp_en, 'fq2fa': _sum_fq2fa_en, 'host': _sum_host_en,
    'virus': _sum_virus_en, 'assembly': _sum_asm_en,
    'verify': _sum_verify_en,
    'hostana': _sum_hostana_en,
    'orf': _sum_orf_en, 'orfa': _sum_orfa_en, 'phylo': _sum_phylo_en,
    'primer': _sum_primer_en, 'gbdraw': _sum_gbdraw_en,
}
# 各阶段关键产物（文件名 glob）
STAGE_OUTPUTS = {
    'subsample': ['sub_R1.fastq.gz', 'sub_R2.fastq.gz', 'subsample.json'],
    'fastp': ['fastp_R1.fastq.gz', 'fastp_R2.fastq.gz', 'fastp_report.html'],
    'fq2fa': ['conv_R1.fa.gz', 'conv_R2.fa.gz'],
    'hostana': ['host_prediction.tsv', 'host_summary.tsv',
                'sankey_host.html', 'sunburst_host.html'],
    'host': ['kept_R1.fastq.gz', 'kept_R2.fastq.gz', 'stats.json'],
    'virus': ['virus_summary.tsv', 'viral_R1.fastq.gz', 'viral_R2.fastq.gz',
              'summary.json'],
    'assembly': ['contigs.filtered.fasta', 'virus_contigs.tsv',
                 'virus_classification.tsv', 'summary.json'],
    'verify': ['calls.tsv', 'summary.json'],
    'orf': ['*'],
    'orfa': ['orf_annotation.tsv', 'orf_function_summary.tsv',
             'orf_family_summary.tsv', 'contig_function_profile.tsv',
             'orf_annotation.gff3', 'genome_diagrams/*.svg', 'summary.json'],
    'phylo': ['*/*'],
    'primer': ['primers.tsv'],
    'gbdraw': ['*.svg', 'summary.json'],
    'report': ['report.html'],
}
# 卡片「查看」按钮指向的结果页文件（可在线预览的入口文件；
# 报告/图表类 HTML 直接内嵌打开，表格类走预览弹窗渲染）
STAGE_VIEW = {
    'report': 'report.html',
    'fastp': 'fastp_report.html',
    'hostana': 'sankey_host.html',
    'orfa': 'orf_annotation.tsv',
    'primer': 'primers.tsv',
    'virus': 'virus_summary.tsv',
    'assembly': 'virus_contigs.tsv',
    'verify': 'calls.tsv',
    'host': 'stats.json',
}


def _list_outputs(sample_dir, stage_key, stage_dirname, limit=10):
    import glob as _glob
    pats = STAGE_OUTPUTS.get(stage_key, ['*'])
    out, seen = [], set()
    for pat in pats:
        for p in sorted(_glob.glob(os.path.join(sample_dir, stage_dirname, pat))):
            if os.path.isfile(p) and p not in seen:
                seen.add(p)
                out.append({'name': os.path.relpath(p, sample_dir),
                            'size': os.path.getsize(p)})
    return out[:limit]


def pipeline_overview(sample_dir, lang=None):
    """样品管道总览：各阶段状态/摘要/产物/可运行性（GUI 模块化视图数据）。

    返回 {'stages': [...], 'groups': [(组名, [阶段键]), ...]}。
    fastp/seqkit/gbdraw 未安装时对应阶段状态为 unavailable（灰色禁用）。
    lang: 'zh' | 'en'（缺省取平台配置语言），影响阶段名/摘要/提示文案。
    """
    if lang is None:
        from .config import get_config
        lang = get_config().lang
    try:
        from .preprocess import fastp_available
        has_fastp = fastp_available()
    except Exception:
        has_fastp = False
    has_seqkit = _seqkit_ok()
    has_gbdraw = _gbdraw_ok()
    has_dfv = _dfv_ok()
    try:
        from .orf_annot import annotation_engine
        has_orfa = bool(annotation_engine())
    except Exception:
        has_orfa = False
    try:
        from .verify import verify_engine
        _ve = verify_engine()
        has_verify = bool(_ve.get('diamond') or _ve.get('mmseqs'))
    except Exception:
        has_verify = False
    sum_tab = STAGE_SUMMARIES if lang == 'zh' else STAGE_SUMMARIES_EN
    ran_map, rows = {}, []
    for key, name, dirname, sfile in PIPELINE_STAGES:
        disp_name = name if lang == 'zh' else \
            PIPELINE_STAGE_NAMES_EN.get(key, name)
        summary = _load_done_summary(sample_dir, dirname, sfile or 'summary.json')
        if key == 'report':
            import glob as _glob
            summary = ({'report': 'report.html'}
                       if _glob.glob(os.path.join(sample_dir, dirname,
                                                  'report.html')) else None)
        if key == 'fastp':
            # 以质控产物文件判定完成；摘要走 preprocess._summary 解析原生报告
            from .preprocess import _summary as _fastp_summary
            prep = os.path.join(sample_dir, '00_prep')
            summary = (_fastp_summary(prep)
                       if os.path.isfile(os.path.join(prep, 'fastp_R1.fastq.gz'))
                       else None)
        has_summary = summary is not None
        # 「跳过」不等于「完成」：summary 里 skipped=true（或把 skipped
        # 写成 reason 且无实际产物）时，该阶段虽然跑过但没结果，既要让
        # 用户看见真实情况，也要允许分步运行时重新发起。
        skipped_reason = ''
        if isinstance(summary, dict) and summary.get('skipped'):
            skipped_reason = str(summary.get('reason') or '')
        if key == 'primer' and isinstance(summary, dict) and \
                not summary.get('n_primers'):
            skipped_reason = skipped_reason or str(summary.get('reason') or '')
        if key == 'gbdraw' and isinstance(summary, dict) and \
                not summary.get('plots'):
            # 一条图都没画出来（无 contig / 引擎失败）也算无效跑，
            # 否则卡片会显示 done + “基因组图 0 张”误导读者。
            skipped_reason = skipped_reason or str(
                summary.get('reason') or '未生成基因组图（无可用 contig）')
        is_skipped = bool(has_summary and skipped_reason)
        # 依赖判定用的是「上游跑过」，不是「上游成功」：
        # skipped 阶段确实跑过且产物目录就绪，下游取到空数据会自行跳过；
        # 若把它当未完成，用户反而无法手动补跑下游，更糟。
        ran_map[key] = has_summary
        deps_ok = all(ran_map.get(d) for d in STAGE_DEPS.get(key, []))
        missing_tool = None
        if key == 'fastp' and not has_fastp:
            missing_tool = _MISSING_TOOL_MSG['fastp']
        elif key == 'fq2fa' and not has_seqkit:
            missing_tool = _MISSING_TOOL_MSG['fq2fa']
        elif key == 'gbdraw' and not (has_gbdraw or has_dfv):
            missing_tool = _MISSING_TOOL_MSG['gbdraw']
        elif key == 'orfa' and not has_orfa:
            missing_tool = _MISSING_TOOL_MSG['orfa']
        elif key == 'verify' and not has_verify:
            missing_tool = _MISSING_TOOL_MSG['verify']
        if missing_tool:
            status = 'unavailable'
            summary_txt = missing_tool[0] if lang == 'zh' else missing_tool[1]
        elif is_skipped:
            status = 'skipped'
            summary_txt = (f'已跳过：{skipped_reason}'
                           if skipped_reason else '已跳过') if lang == 'zh' \
                else (f'Skipped: {skipped_reason}' if skipped_reason
                      else 'Skipped')
        elif has_summary:
            status = 'done'
            fn = sum_tab.get(key)
            summary_txt = fn(summary) if fn else \
                ('报告已生成' if lang == 'zh' else 'Report generated')
        else:
            status = 'ready' if deps_ok else 'blocked'
            summary_txt = ''
        view = STAGE_VIEW.get(key)
        view_ok = bool(view) and os.path.isfile(
            os.path.join(sample_dir, dirname, view))
        card = {'stage': key, 'name': disp_name, 'dir': dirname,
                'status': status, 'summary': summary_txt,
                'view': (f'/api/samples/{os.path.basename(sample_dir)}'
                         f'/{dirname}/{view}') if view_ok else None,
                'outputs': (_list_outputs(sample_dir, key, dirname)
                            if has_summary else [])}
        rows.append(card)
    groups = STAGE_GROUPS if lang == 'zh' else STAGE_GROUPS_EN
    return {'stages': rows, 'groups': groups}


def save_sample_input(sample_dir, r1, r2, sample_name, project=None):
    """存档样品输入（管道逐级运行时自动带出）；project 为可选项目分组标签。"""
    prep = os.path.join(sample_dir, '00_prep')
    os.makedirs(prep, exist_ok=True)
    with safe_open(os.path.join(prep, 'input.json'), 'wt') as f:
        json.dump({'sample': sample_name, 'r1': str(r1), 'r2': str(r2 or ''),
                   'project': str(project or '')},
                  f, ensure_ascii=False, indent=2)


def load_sample_input(sample_dir):
    d = _load_done_summary(sample_dir, '00_prep', 'input.json') or {}
    return d.get('r1'), d.get('r2') or None, d.get('project') or None


# ------------------------------------------------------------------
# 项目清单（project.json）：样品目录即项目根，清单记录「这次的输入、
# 参数、已完成阶段、运行历史」，使一次运行可复现、可追溯。
# 位置固定在样品目录下（不随阶段目录变动），老样品无此文件时读取端
# 一律返回空 dict，由 ensure_project_manifest() 从既有产物回填。
# ------------------------------------------------------------------
MANIFEST_NAME = 'project.json'
# 阶段键 → 实际断点标记文件名（部分阶段标记名与键不同）
_DONE_ALIAS = {
    'host': '.host_removal.done',
    'virus': '.virus_screen.done',
    'hostana': '.host_analysis.done',
}
# 需要落盘记录的分析参数（键与 run_analysis 形参同名，值是它的取值）
_MANIFEST_PARAM_KEYS = (
    'assembly_mode', 'assembly_input', 'confidence', 'threads',
    'memory_gb', 'subsample', 'min_contig_len', 'top_n_refs',
    'tree_tool', 'tree_sampling', 'primer_mode', 'do_trim',
    'do_specificity', 'min_orf_aa', 'min_viral_pairs', 'fastp_dedup',
    'do_fq2fa', 'gbdraw_max', 'plot_engine',
)
# 这些参数含绝对路径，只记文件名以避免清单泄露/绑定本机路径
_MANIFEST_PATH_KEYS = ('db_host', 'db_virus', 'ncbi_refs', 'chunk_dir')


def _manifest_path(sample_dir):
    return os.path.join(sample_dir, MANIFEST_NAME)


def load_project_manifest(sample_dir):
    """读项目清单；不存在或损坏时返回空 dict（调用方自行兜底）。"""
    p = _manifest_path(sample_dir)
    if not os.path.isfile(p):
        return {}
    try:
        with safe_open(p) as f:
            d = json.load(f)
        return d if isinstance(d, dict) else {}
    except (OSError, ValueError):
        return {}


def save_project_manifest(sample_dir, manifest):
    """原子写项目清单（先写临时文件再替换，避免中途崩溃留下半截）。"""
    p = _manifest_path(sample_dir)
    tmp = p + '.tmp'
    os.makedirs(sample_dir, exist_ok=True)
    with safe_open(tmp, 'wt') as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)
    os.replace(tmp, p)
    return p


def _manifest_params(kw):
    """从 run_analysis 关键字参数提取待记录的分析参数。"""
    out = {}
    for k in _MANIFEST_PARAM_KEYS:
        if k in kw and kw[k] is not None:
            out[k] = kw[k]
    for k in _MANIFEST_PATH_KEYS:
        v = kw.get(k)
        if v:
            out[k] = os.path.basename(str(v))
    return out


def update_project_manifest(sample_dir, *, sample=None, r1=None, r2=None,
                            params=None, stage=None, stages_done=None,
                            status=None, duration_s=None,
                            replace_stages=False):
    """增量更新项目清单；已存在的字段只覆盖本次给出的部分。

    - 首次调用写入 created / project / input
    - 每次运行写 last_run 与 params（本次实际使用的参数）
    - 阶段完成后追加 stages_done（去重、按 STAGE_ORDER 排序）；
      replace_stages=True 时直接替换而不做并集（用于「跳过」回退）
    - runs 保存最近 20 次运行记录（阶段集合、耗时、状态）
    """
    m = load_project_manifest(sample_dir)
    now = time.strftime('%Y-%m-%dT%H:%M:%S')
    m.setdefault('schema', 1)
    # 项目名以 00_prep/input.json 里用户填的为准（唯一权威源）；
    # 没有时才用样品目录名兑底。sample 形参只是样品名，不用于推断项目。
    _up = None
    try:
        _up = load_sample_input(sample_dir)[2]
    except Exception:
        _up = None
    if _up:
        if m.get('project') != _up:
            m['project'] = _up
    elif not m.get('project'):
        m['project'] = (_safe_sample_name(sample) if sample
                        else os.path.basename(str(sample_dir).rstrip('\\/')))
    m.setdefault('created', now)
    if r1:
        m['input'] = {'r1': str(r1), 'r2': str(r2 or '')}
    if params:
        m['params'] = params
    if stages_done is not None:
        if replace_stages:
            prev = set(stages_done)
        else:
            prev = set(m.get('stages_done') or [])
            prev |= set(stages_done)
        m['stages_done'] = [s for s in STAGE_ORDER if s in prev]
    m['last_run'] = now
    if status:
        m['last_status'] = status
    if stage:
        m['last_stage'] = stage
    if duration_s is not None:
        m['last_duration_s'] = round(float(duration_s), 1)
    # runs 只记一次运行的终态（ok/failed）；运行中的进度快照不写历史，
    # 否则同一次运行会先记 running 再记 ok，历史里凭空多一条。
    if stages_done and (status or 'ok') != 'running':
        runs = m.get('runs') or []
        runs.append({'at': now, 'stages': [s for s in STAGE_ORDER
                                           if s in set(stages_done)],
                     'status': status or 'ok',
                     'duration_s': (round(float(duration_s), 1)
                                    if duration_s is not None else None)})
        m['runs'] = runs[-20:]
    save_project_manifest(sample_dir, m)
    return m


def backfill_manifest_from_products(sample_dir, r1=None, r2=None):
    """从既有产物回填清单（老样品无 project.json 时用）。

    只做事实推断，且状态口径与 pipeline_overview 保持一致：
      - 只有真正跑出结果的阶段才计入 stages_done；
      - 「跑过但无结果」(skipped：无病毒 contigs / 无引物 / 无基因组图)
        不计入，它们仍属待运行，与 GUI 卡片显示一致。
    输入来源优先级：00_prep/input.json > 调用方传入。
    """
    m = load_project_manifest(sample_dir)
    changed = False
    ir1 = ir2 = _proj = None
    try:
        ir1, ir2, _proj = load_sample_input(sample_dir)
    except Exception:
        pass
    # 项目名优先级：已有清单 > 00_prep/input.json 用户填的 > 目录名
    if not m.get('project'):
        m['project'] = (_proj or
                        os.path.basename(str(sample_dir).rstrip('\\/')))
        changed = True
    elif _proj and m.get('project') != _proj and m.get('backfilled'):
        # 回填生成的目录名不如用户填的真实项目名，以用户填的为准
        m['project'] = _proj
        changed = True
    if not m.get('input'):
        ir1, ir2 = ir1 or r1, ir2 or r2
        if ir1:
            m['input'] = {'r1': str(ir1), 'r2': str(ir2 or '')}
            changed = True
    # 从 pipeline_overview 取真实状态，避免自己再猜一遍
    try:
        done = [s['stage'] for s in pipeline_overview(sample_dir)['stages']
                if s['status'] == 'done']
    except Exception:
        # 抾底：至少标出 .done 标记存在的阶段
        done = [key for key, _n, dirname, _sf in PIPELINE_STAGES
                if os.path.isfile(os.path.join(
                    sample_dir, dirname, f'.{key}.done'))]
    if done != (m.get('stages_done') or []):
        m['stages_done'] = [s for s in STAGE_ORDER if s in set(done)]
        changed = True
    if changed:
        m.setdefault('schema', 1)
        m.setdefault('created', time.strftime('%Y-%m-%dT%H:%M:%S'))
        m.setdefault('backfilled', True)
        save_project_manifest(sample_dir, m)
    return m
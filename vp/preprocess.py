# -*- coding: utf-8 -*-
"""
阶段⓪ 预处理（可选）：fastp 自动去接头 + 低质量裁剪（双端 / 单端）。
fastp 未安装时该阶段自动跳过（SPAdes 自带 BayesHammer 纠错，质控非必需）。
输出: 00_prep/fastp_R1.fastq.gz(与 R2)、fastp_report.json / fastp_report.html

阶段⓪·序列转换（可选，推荐）：FASTQ→FASTA 预转换（seqkit fq2fa）。
kunpeng 分类不读质量值，转换后分类阶段少解压一半碱基行，速度更快；
转换产物同时供 ①宿主去除 的分类步骤直接使用（FASTA 输入，中间数据更小）。
输出: 00_prep/conv_R1.fa.gz(与 R2)、fq2fa.json
"""
import os

from .config import get_config
from .utils import (check_path, safe_open, run_cmd, is_step_done,
                    mark_step_done, log_res_plan, fmt_size)
import json


def fastp_available():
    try:
        return bool(get_config().tool('fastp'))
    except Exception:
        return False


# ------------------------------------------------------------------
# FASTQ → FASTA 预转换（分类加速；可选，推荐运行）
# ------------------------------------------------------------------
def seqkit_available():
    try:
        return bool(get_config().tool('seqkit'))
    except Exception:
        return False


def convert_fq2fa(sample_dir, r1, r2, threads=None, logger=None, force=False,
                  progress=None):
    """FASTQ→FASTA 预转换（kunpeng 分类用；分类不读质量值）。

    输入已是 FASTA 时无需转换，直接原样返回（分类器原生支持 FASTA）。
    返回 {'fa_r1','fa_r2','files':[...]}；seqkit 未安装返回 None（阶段跳过，
    宿主去除会直接用 FASTQ 分类）。
    """
    step = 'fq2fa'
    out_dir = check_path(os.path.join(sample_dir, '00_prep'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)

    def _is_fastq(p):
        return p and str(p).lower().endswith(('.fastq', '.fq',
                                              '.fastq.gz', '.fq.gz'))
    if not _is_fastq(r1) and not (r2 and _is_fastq(r2)):
        if logger:
            logger.log("输入已是 FASTA，无需 FASTQ→FASTA 转换")
        with safe_open(os.path.join(out_dir, 'fq2fa.json'), 'wt') as f:
            json.dump({'skipped': True, 'files': []}, f)
        mark_step_done(out_dir, step)
        return {'fa_r1': r1, 'fa_r2': r2, 'files': [], 'skipped': True}

    try:
        seqkit = get_config().tool('seqkit')
    except (FileNotFoundError, RuntimeError):
        if logger:
            logger.log("未检测到 seqkit.exe，跳过 FASTQ→FASTA 预转换"
                       "（可选步骤；宿主去除将直接用 FASTQ 分类）")
        return None

    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("FASTQ→FASTA 转换已完成，跳过")
        return _fq2fa_summary(out_dir)

    tasks = [(r1, os.path.join(out_dir, 'conv_R1.fa.gz'))]
    if r2:
        tasks.append((r2, os.path.join(out_dir, 'conv_R2.fa.gz')))

    in_bytes = sum(os.path.getsize(check_path(p, must_exist=True))
                   for p, _ in tasks)
    if logger:
        logger.log(f"阶段⓪·序列转换: FASTQ→FASTA（seqkit fq2fa，"
                   f"{'双端' if r2 else '单端'}）")
        log_res_plan(logger, 'FASTQ→FASTA 转换', threads=threads,
                     mem_gb=max(threads or 4, 4) * 0.15,
                     disk_gb=in_bytes * 4.0 / 1e9,
                     note="内存极小；磁盘为转换产物（FASTA.gz ≈输入 gz×4~5）")
    for i, (src, dst) in enumerate(tasks):
        if progress:
            progress(i / len(tasks) * 0.95,
                     f"转换 {os.path.basename(str(src))} ({i + 1}/{len(tasks)})")
        run_cmd([seqkit, 'fq2fa', '-w', '0', '-j', str(threads or 4),
                 check_path(src, must_exist=True),
                 '-o', check_path(dst, must_exist=False, in_platform=True)],
                logger=logger)
    if progress:
        progress(1.0, '转换完成')
    mark_step_done(out_dir, step)
    s = _fq2fa_summary(out_dir)
    with safe_open(os.path.join(out_dir, 'fq2fa.json'), 'wt') as f:
        json.dump({'skipped': False, 'files': s['files'],
                   'sizes': {f: os.path.getsize(os.path.join(out_dir, f))
                             for f in s['files']}}, f, ensure_ascii=False, indent=1)
    if logger:
        outs = ', '.join(f"{os.path.basename(t[1])}({fmt_size(os.path.getsize(t[1]))})"
                         for t in tasks if os.path.isfile(t[1]))
        logger.log(f"FASTQ→FASTA 完成: {outs}（分类阶段直接用 FASTA，"
                   f"不再解压质量行）")
    return s


def _fq2fa_summary(out_dir):
    out = {'fa_r1': os.path.join(out_dir, 'conv_R1.fa.gz'),
           'fa_r2': os.path.join(out_dir, 'conv_R2.fa.gz'),
           'files': [], 'skipped': False}
    for k in ('fa_r1', 'fa_r2'):
        if not os.path.isfile(out[k]):
            out[k] = None
    out['files'] = [os.path.basename(p) for p in (out['fa_r1'], out['fa_r2'])
                    if p and os.path.isfile(p)]
    return out


def run_fastp(sample_dir, r1, r2, threads=None, logger=None, force=False,
              dedup=False):
    """fastp 质控。返回产物 dict；未安装 fastp 时返回 None（跳过）。

    dedup=True 时启用 fastp --dedup（按序列去 PCR 重复，宏基因组慎用：
    高深度真实重复会被误删，病毒筛查/组装建议默认关闭）。
    """
    step = 'fastp'
    cfg = get_config()
    try:
        fastp = cfg.tool('fastp')
    except (FileNotFoundError, RuntimeError):
        fastp = None
    if not fastp:
        if logger:
            logger.log("未检测到 fastp.exe，跳过质控（可选步骤；"
                       "可将 fastp.exe 放到平台目录后重跑）")
        return None

    out_dir = check_path(os.path.join(sample_dir, '00_prep'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("Fastp 质控已完成，跳过")
        return _summary(out_dir)

    o1 = os.path.join(out_dir, 'fastp_R1.fastq.gz')
    o2 = os.path.join(out_dir, 'fastp_R2.fastq.gz') if r2 else None
    json_rep = os.path.join(out_dir, 'fastp_report.json')
    html_rep = os.path.join(out_dir, 'fastp_report.html')

    if logger:
        logger.log(f"阶段⓪ Fastp 质控（{'双端' if r2 else '单端'}"
                   f"{'，--dedup 去重' if dedup else ''}）")
    cmd = [fastp, '-i', check_path(r1, must_exist=True),
           '-o', o1, '-w', str(threads or cfg.threads),
           '-j', json_rep, '-h', html_rep,
           '--qualified_quality_phred', '20',
           '--length_required', '50']
    if dedup:
        cmd.append('--dedup')
    if r2:
        cmd += ['-I', check_path(r2, must_exist=True), '-O', o2,
                '--detect_adapter_for_pe']
    run_cmd(cmd, logger=logger)
    # 质控后 0 条 reads：后续所有阶段必然失败，这里给出直接原因
    s = _summary(out_dir)
    if s.get('reads_after') == 0:
        raise RuntimeError(
            "Fastp 质控后 reads 为 0（全部被过滤）。常见原因：输入质量值过低"
            "（当前阈值 Q20 / 最短 50bp）或输入文件为空。可跳过 ⓪ 质控阶段，"
            "或换用更高质量的输入数据")
    mark_step_done(out_dir, step)
    if logger:
        logger.log("Fastp 质控完成: " + os.path.basename(str(o1)))
    return _summary(out_dir)


def _summary(out_dir):
    """读 fastp json 摘要（读取失败不报错，只给文件路径）。

    fastp(Windows) 的 JSON 含未转义反斜杠路径，属无效 JSON——
    先标准解析，失败则把反斜杠替换为 '/' 再解析（我们只读统计字段，
    路径字段不使用）。
    """
    p = os.path.join(out_dir, 'fastp_report.json')
    s = {'fastp_r1': 'fastp_R1.fastq.gz',
         'fastp_r2': 'fastp_R2.fastq.gz',
         'report_html': 'fastp_report.html'}
    try:
        with safe_open(p) as f:
            text = f.read()
        try:
            d = json.loads(text)
        except ValueError:
            d = json.loads(text.replace('\\', '/'))
        summ = d.get('summary', {}) or {}
        bf = summ.get('before_filtering', {}) or {}
        af = summ.get('after_filtering', {}) or {}
        filt = d.get('filtering_result', {}) or {}
        s['reads_before'] = bf.get('total_reads')
        s['reads_after'] = af.get('total_reads')
        s['q30_rate'] = (round(af.get('q30_rate', 0) * 100, 2)
                         if af.get('q30_rate') is not None else None)
        s['lowq_removed'] = filt.get('low_quality_reads')
    except (OSError, ValueError):
        pass
    return s

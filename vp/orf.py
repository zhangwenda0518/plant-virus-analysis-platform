# -*- coding: utf-8 -*-
"""
阶段④ ORF 预测：
- pyrodigal（meta 模式）：基因预测，输出 faa/ffn/gff
- pyrodigal-rv（RNA/反转录病毒专用，若安装可用）
- orfipy（6 框 start-stop）已默认去除，仅显式指定时运行（历史 04b 无模型兜底）
默认同时跑 pyrodigal + pyrodigal_rv；输入默认为阶段③的病毒 contigs。
"""
import os
import csv
import json

from Bio import Seq
from Bio.Seq import Seq as BioSeq

from .config import get_config
from .utils import (check_path, safe_open, run_cmd, iter_fasta, write_fasta_record,
                    count_fasta_seqs, is_step_done, mark_step_done)


def verified_contigs(sample_dir, logger=None):
    """从 03b_verify/calls.tsv 取「验证通过」的 contig id。

    通过档位 = known（已知病毒，强证据）+ novel（新病毒候选，有证据但库无记录）。
    domain_only（仅结构域信号，无序列同源）与 unclassified（无任何证据）
    不进注释，但仍保留在 calls.tsv 作新病毒候选池。

    返回 (ids, 状态)。状态：'ok' 正常读到；'absent' 无验证产物；
    'empty' 有产物但无通过项；'error' 读取失败（已记日志）。
    """
    PASS_CALLS = ('known', 'novel')
    # 阶段产物优先，其次工具运行目录（手动跑过也可以被采纳）
    cands = [os.path.join(sample_dir, '03b_verify', 'calls.tsv'),
             os.path.join(sample_dir, 'verify', 'calls.tsv')]
    path = next((p for p in cands if os.path.isfile(p)), None)
    if not path:
        return [], 'absent'
    ids = []
    try:
        with safe_open(path) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                call = (row.get('call') or '').strip().lower()
                cid = (row.get('contig') or '').strip()
                if cid and call in PASS_CALLS:
                    ids.append(cid)
    except Exception as e:
        if logger:
            logger.log(f"读取验证结果失败（{e}），回退全量 contigs", 'WARN')
        return [], 'error'
    if logger:
        logger.log(f"验证通过 contigs: {len(ids)} 条"
                   f"（{'/'.join(PASS_CALLS)}）← {os.path.relpath(path, sample_dir)}")
    return ids, ('ok' if ids else 'empty')


def extract_viral_contigs(assembly_dir, max_n=0, only=None):
    """从 03_assembly 提取病毒 contigs 子集。返回 (fasta路径, [contig ids])。

    only：限定 id 集合（如 03b_verify 验证通过的 contig）。传空列表时
    不提取任何序列（避免静默回退成全量）。
    """
    import csv
    summary_p = os.path.join(assembly_dir, 'summary.json')
    with safe_open(summary_p) as f:
        summary = json.load(f)
    ids = list(summary.get('viral_contigs') or [])
    if not ids:
        # 无病毒 contig：退化为最长的前 10 条
        lens = sorted(((len(s), h.split()[0]) for h, s in
                       iter_fasta(os.path.join(assembly_dir, 'contigs.filtered.fasta'))),
                      reverse=True)
        ids = [c for _, c in lens[:10]]
    if only is not None:
        # 保留 verify 判定顺序，且只留确实存在于 contigs 里的
        keep_set = set(only)
        ids = [c for c in ids if c in keep_set] or \
              [c for c in only if c in keep_set]
    if max_n:
        ids = ids[:max_n]
    keep = set(ids)
    out = check_path(os.path.join(assembly_dir, 'viral_contigs.fasta'),
                     must_exist=False, in_platform=True)
    with safe_open(out, 'wt') as f:
        for h, s in iter_fasta(os.path.join(assembly_dir, 'contigs.filtered.fasta')):
            if h.split()[0] in keep:
                write_fasta_record(f, h.split()[0], s)
    return out, ids


# ------------------------------------------------------------------
# orfipy
# ------------------------------------------------------------------
def run_orfipy(in_fasta, out_dir, min_aa=100, threads=None, logger=None):
    """orfipy ORF 预测（等价 EMBOSS getORF）。返回输出文件 dict。"""
    from .config import get_config
    try:
        orfipy = get_config().tool('orfipy')
    except FileNotFoundError:
        raise FileNotFoundError("orfipy 未安装 (python -m pip install orfipy)")
    # orfipy 的 --outdir/--pep 等对相对路径敏感（相对输出 + 绝对输入会退出码 1），
    # 统一规范为绝对路径；--outdir 把运行日志收进任务目录而非输入文件旁边
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.join(out_dir, 'orfipy')
    out = {'pep': prefix + '_pep.fa', 'nt': prefix + '_nt.fa', 'bed': prefix + '.bed'}
    log_dir = os.path.join(out_dir, 'orfipy_logs')
    min_nt = int(min_aa) * 3
    if logger:
        logger.log(f"orfipy ORF 预测 (min {min_aa} aa / {min_nt} nt)")
    run_cmd([orfipy, in_fasta, '--pep', out['pep'], '--dna', out['nt'],
             '--bed', out['bed'], '--min', str(min_nt),
             '--outdir', log_dir,
             '--procs', str(threads or 4)], logger=logger)
    return out


# ------------------------------------------------------------------
# pyrodigal / pyrodigal-rv
# ------------------------------------------------------------------
def _translate(seq_str, begin, end, strand):
    sub = BioSeq(seq_str[begin:end])
    if strand == -1:
        sub = sub.reverse_complement()
    prot = str(sub.translate(table=11))
    if prot.endswith('*'):
        prot = prot[:-1]
    return prot


def run_pyrodigal(in_fasta, out_prefix, tag='pyrodigal', logger=None,
                  progress=None):
    """pyrodigal(meta 模式) 基因预测。tag 用于区分 pyrodigal / pyrodigal_rv。

    返回 {'faa', 'ffn', 'gff', 'genes': n}；依赖不可用时返回 None。
    """
    try:
        if tag == 'pyrodigal':
            import pyrodigal as pdlib
        else:
            import pyrodigal_rv as pdlib
    except ImportError:
        if logger:
            logger.log(f"{tag} 未安装，跳过", "WARN")
        return None

    # 兼容 pyrodigal v2(OrfFinder) / v3(GeneFinder) 与 pyrodigal-rv(ViralGeneFinder)
    if tag == 'pyrodigal_rv':
        finder_cls = (getattr(pdlib, 'ViralGeneFinder', None)
                      or getattr(pdlib, 'OrfFinder', None)
                      or getattr(pdlib, 'GeneFinder', None))
    else:
        finder_cls = (getattr(pdlib, 'GeneFinder', None)
                      or getattr(pdlib, 'OrfFinder', None))
    orf_finder = finder_cls(meta=True)
    faa = out_prefix + '.faa'
    ffn = out_prefix + '.ffn'
    gff = out_prefix + '.gff'
    n_genes = 0
    n_seqs = 0
    # 先数序列总数供进度分母（FASTA 很小，代价可忽略）
    try:
        n_seqs = count_fasta_seqs(in_fasta) or 1
    except (OSError, ValueError):
        n_seqs = 1
    with safe_open(faa, 'wt') as f_prot, safe_open(ffn, 'wt') as f_nt, \
            safe_open(gff, 'wt') as f_gff:
        f_gff.write("##gff-version 3\n")
        for i, (header, s) in enumerate(iter_fasta(in_fasta), 1):
            cid = header.split()[0]
            seq_str = s.upper()
            preds = list(orf_finder.find_genes(seq_str))
            for g in preds:
                begin, end, strand = g.begin, g.end, g.strand
                # pyrodigal v2/v3 坐标均为 1-based closed（Prodigal 口径）：
                # Python 切片须换算 [begin-1:end]，否则正链基因整体移码
                begin0 = begin - 1
                gid = f"{cid}_{n_genes + 1}"
                prot = _translate(seq_str, begin0, end, strand)
                nt = seq_str[begin0:end]
                if strand == -1:
                    nt = str(BioSeq(nt).reverse_complement())
                write_fasta_record(f_prot, f"{gid} #{begin} #{end} {strand}", prot)
                write_fasta_record(f_nt, gid, nt)
                f_gff.write(f"{cid}\t{tag}\tCDS\t{begin}\t{end}\t.\t"
                            f"{'+' if strand == 1 else '-'}\t0\tID={gid}\n")
                n_genes += 1
            if progress:
                progress(i / n_seqs * 0.95, f'{tag}: {i}/{n_seqs} 条 contigs')
    if logger:
        logger.log(f"{tag} 预测: {n_genes} 个基因")
    return {'faa': faa, 'ffn': ffn, 'gff': gff, 'genes': n_genes}


# ------------------------------------------------------------------
# 阶段入口
# ------------------------------------------------------------------
def predict_orfs(sample_dir, min_aa=100, threads=None, logger=None, force=False,
                 assembly_dir=None, progress=None, tools=None):
    """tools: 运行哪些 ORF 预测工具。None=默认（pyrodigal + pyrodigal_rv）；
    可选 'orfipy' / 'pyrodigal' / 'pyrodigal_rv'（逗号或列表）。orfipy 六框已
    默认去除，仅显式指定时运行（无基因模型的 04b 兜底）。"""
    step = 'orf'
    out_dir = check_path(os.path.join(sample_dir, '04_orf'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    summary_file = os.path.join(out_dir, 'summary.json')
    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("阶段④ORF 预测已完成，跳过")
        with safe_open(summary_file) as f:
            return json.load(f)

    def _sp(frm, to):
        return (lambda p, m: progress(frm + p * (to - frm), m)
                if progress else None)

    a_dir = check_path(assembly_dir or os.path.join(sample_dir, '03_assembly'),
                       must_exist=True)
    # 上游收窄：只注释 03b_verify 判定通过（known/novel）的 contig。
    # 无验证产物时优雅降级为全量（老样品 / 未选该阶段），但一定留日志。
    vid, vstate = verified_contigs(sample_dir, logger=logger)
    if vstate == 'ok':
        viral_fa, ids = extract_viral_contigs(a_dir, only=vid)
    elif vstate == 'empty':
        if logger:
            logger.log("验证结果中无通过项（known/novel），ORF 预测跳过",
                       'WARN')
        raise RuntimeError(
            '03b_verify 判定无通过项（known/novel），无可注释 contig。'
            '若确需全量注释，请改用未收窄的输入或检查验证参数')
    else:
        if logger:
            logger.log(f"无 03b_verify 验证结果（{vstate}），"
                       "ORF 预测回退到全量病毒 contigs", 'WARN')
        viral_fa, ids = extract_viral_contigs(a_dir)
    if logger:
        logger.log(f"阶段④ ORF 预测: {len(ids)} 条 contigs"
                   f"（来源 {'验证通过' if vstate == 'ok' else '全量回退'}\）")

    if not ids:
        raise RuntimeError(
            '03_assembly 中没有病毒 contigs（viral_contigs 列表为空且无 '
            'contigs.filtered.fasta 兜底）——请先完成 ③组装/分类 阶段')
    summary = {'stage': step, 'input': os.path.basename(viral_fa),
               'n_contigs': len(ids)}
    # 工具选择：默认 pyrodigal + pyrodigal_rv（orfipy 已默认去除，仅显式指定时跑）
    if tools is None:
        tools = ['pyrodigal', 'pyrodigal_rv']
    elif isinstance(tools, str):
        tools = [t.strip() for t in tools.split(',') if t.strip()]
    tools = [t for t in tools if t in ('orfipy', 'pyrodigal', 'pyrodigal_rv')]

    if 'orfipy' in tools:
        if progress:
            progress(0.05, f'orfipy 六框 ORF（{len(ids)} 条 contigs）')
        orfipy_out = run_orfipy(viral_fa, out_dir, min_aa=min_aa,
                                threads=threads, logger=logger)
        summary['orfipy'] = {k: os.path.basename(v) for k, v in orfipy_out.items()}
        try:
            summary['n_orfs'] = count_fasta_seqs(
                os.path.join(out_dir, orfipy_out['pep']))
        except (OSError, ValueError):
            summary['n_orfs'] = None
        if progress:
            progress(0.45, f"orfipy 完成：{summary['n_orfs']} 个 ORF")

    if 'pyrodigal' in tools:
        pyro = run_pyrodigal(viral_fa, os.path.join(out_dir, 'pyrodigal'),
                             tag='pyrodigal', logger=logger,
                             progress=_sp(0.45, 0.75))
        summary['pyrodigal'] = ({'faa': os.path.basename(pyro['faa']),
                                 'ffn': os.path.basename(pyro['ffn']),
                                 'gff': os.path.basename(pyro['gff']),
                                 'genes': pyro['genes']} if pyro else None)

    if 'pyrodigal_rv' in tools:
        pyro_rv = run_pyrodigal(viral_fa, os.path.join(out_dir, 'pyrodigal_rv'),
                                tag='pyrodigal_rv', logger=logger,
                                progress=_sp(0.75, 0.95))
        summary['pyrodigal_rv'] = ({'faa': os.path.basename(pyro_rv['faa']),
                                    'ffn': os.path.basename(pyro_rv['ffn']),
                                    'gff': os.path.basename(pyro_rv['gff']),
                                    'genes': pyro_rv['genes']} if pyro_rv else None)
    if progress:
        progress(1.0, 'ORF 预测完成')

    with safe_open(summary_file, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    mark_step_done(out_dir, step)
    if logger:
        logger.log("阶段④ 完成")
    return summary

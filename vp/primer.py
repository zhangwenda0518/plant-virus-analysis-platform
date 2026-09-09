# -*- coding: utf-8 -*-
"""
阶段⑥ 引物设计（primer3）：
- conserved 模式：基于 05_phylo 多比对找保守区，在一致序列上设计（跨近缘参考通用）
- plain 模式：对病毒 contigs 直接设计覆盖全长的引物对
- 可选：blastn 扩增子 vs 宿主基因组特异性检查
"""
import os
import json
from collections import Counter

import primer3

from .config import get_config, DIRS, current_host_genome
from .utils import (check_path, safe_open, run_cmd, iter_fasta, fmt_size,
                    is_step_done, mark_step_done)

P3_GLOBAL = {
    'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 24,
    'PRIMER_OPT_TM': 60.0, 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0,
    'PRIMER_MIN_GC': 40.0, 'PRIMER_MAX_GC': 60.0,
    'PRIMER_MAX_POLY_X': 4,
    'PRIMER_PRODUCT_SIZE_RANGE': [[300, 1500]],
    'PRIMER_NUM_RETURN': 3,
    'PRIMER_EXPLAIN_FLAG': 1,
    'PRIMER_PICK_LEFT_PRIMER': 1, 'PRIMER_PICK_RIGHT_PRIMER': 1,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
}


def design_primers_for_seq(name, seq, num_return=3, product_range=None):
    """对单条序列设计引物对。返回引物 dict 列表。

    序列短于默认产物区间（300bp）时自动收缩产物范围，避免 primer3 直接报错。
    """
    opts = dict(P3_GLOBAL)
    opts['PRIMER_NUM_RETURN'] = num_return
    if product_range:
        opts['PRIMER_PRODUCT_SIZE_RANGE'] = [list(product_range)]
    elif len(seq) < 300:
        opts['PRIMER_PRODUCT_SIZE_RANGE'] = [[60, max(len(seq), 80)]]
    res = primer3.bindings.designPrimers(
        {'SEQUENCE_ID': name[:40], 'SEQUENCE_TEMPLATE': seq}, opts)
    pairs = []
    n = res.get('PRIMER_PAIR_NUM_RETURNED', 0)
    for i in range(n):
        pairs.append({
            'name': f"{name}_P{i + 1}",
            'F_seq': res.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''),
            'F_pos': res.get(f'PRIMER_LEFT_{i}', (0, 0))[0] + 1,
            'F_len': res.get(f'PRIMER_LEFT_{i}', (0, 0))[1],
            'F_tm': round(res.get(f'PRIMER_LEFT_{i}_TM', 0), 1),
            'F_gc': round(res.get(f'PRIMER_LEFT_{i}_GC_PERCENT', 0), 1),
            'R_seq': res.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
            'R_pos': res.get(f'PRIMER_RIGHT_{i}', (0, 0))[0] + 1,
            'R_len': res.get(f'PRIMER_RIGHT_{i}', (0, 0))[1],
            'R_tm': round(res.get(f'PRIMER_RIGHT_{i}_TM', 0), 1),
            'R_gc': round(res.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 0), 1),
            'product': res.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0),
            'penalty': round(res.get(f'PRIMER_PAIR_{i}_PENALTY', 0), 3),
            'self_any_max': round(max(res.get(f'PRIMER_LEFT_{i}_SELF_ANY_TH', 0),
                                      res.get(f'PRIMER_RIGHT_{i}_SELF_ANY_TH', 0)), 1),
            'hairpin_max': round(max(res.get(f'PRIMER_LEFT_{i}_HAIRPIN_TH', 0),
                                     res.get(f'PRIMER_RIGHT_{i}_HAIRPIN_TH', 0)), 1),
        })
    return pairs


def conserved_regions_and_consensus(aln_fasta, min_len=400, min_ident=0.9,
                                    min_cov=0.7):
    """从比对找保守区段。返回 [(consensus_seq, aln_start, aln_end)]。"""
    names, seqs = [], []
    for h, s in iter_fasta(aln_fasta):
        names.append(h.split()[0])
        seqs.append(s.upper())
    if not seqs:
        return []
    n = len(seqs)
    L = len(seqs[0])
    cons_chars, ok_cols = [], []
    for col in zip(*seqs):
        non_gap = [c for c in col if c in 'ACGTU']
        if len(non_gap) >= n * min_cov:
            most, cnt = Counter(non_gap).most_common(1)[0]
            frac = cnt / len(non_gap)
            ok_cols.append(frac >= min_ident)
            cons_chars.append(most if frac >= min_ident else 'N')
        else:
            ok_cols.append(False)
            cons_chars.append('-')
    # 连续 ok 区段
    regions = []
    start = None
    for i in range(L):
        if ok_cols[i] and start is None:
            start = i
        elif not ok_cols[i] and start is not None:
            regions.append((start, i))
            start = None
    if start is not None:
        regions.append((start, L))
    out = []
    for s, e in regions:
        if e - s >= min_len:
            sub = ''.join(cons_chars[s:e]).replace('-', '').replace('U', 'T')
            if len(sub) >= min_len:
                out.append((sub, s, e))
    return out


def check_specificity_blast(amplicons_fasta, host_genome=None, logger=None):
    """可选：扩增子 blastn vs 宿主基因组，返回命中 contig/引物 集合。"""
    from .assembly import _ascii_work_base
    cfg = get_config()
    blastn = cfg.tool('blastn')
    makeblastdb = cfg.tool('makeblastdb')
    host = host_genome or current_host_genome()
    if not os.path.isfile(host):
        return {}
    db_dir = _ascii_work_base('vp_blast')
    prefix = os.path.join(db_dir, 'host')
    import glob as _g
    if not _g.glob(prefix + '.n??'):
        if logger:
            logger.log("构建宿主 BLAST 库（仅首次，1.8GB 基因组约需 10-30 分钟）...", "WARN")
        run_cmd([makeblastdb, '-in', host, '-dbtype', 'nucl', '-out', prefix],
                logger=logger)
    out_tsv = os.path.join(db_dir, 'amp_hits.tsv')
    run_cmd([blastn, '-query', amplicons_fasta, '-db', prefix,
             '-outfmt', '6 qseqid sseqid pident length', '-evalue', '1e-5',
             '-max_target_seqs', '1', '-out', out_tsv], logger=logger)
    hits = {}
    with safe_open(out_tsv) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4 and float(parts[3]) >= 100:
                hits.setdefault(parts[0], []).append(parts[1])
    return hits


def design_primers(sample_dir, mode='conserved', num_return=3,
                   logger=None, force=False, do_specificity=False,
                   assembly_dir=None, phylo_dir=None):
    """阶段⑥ 主入口。mode: conserved | plain"""
    step = 'primer'
    out_dir = check_path(os.path.join(sample_dir, '06_primer'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    summary_file = os.path.join(out_dir, 'summary.json')
    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("阶段⑥引物设计已完成，跳过")
        with safe_open(summary_file) as f:
            return json.load(f)

    all_pairs = []
    groups_used = []

    if mode == 'conserved':
        p_dir = check_path(phylo_dir or os.path.join(sample_dir, '05_phylo'),
                           must_exist=True)
        p_summary_p = os.path.join(p_dir, 'summary.json')
        with safe_open(p_summary_p) as f:
            p_summary = json.load(f)
        for g in p_summary.get('groups', []):
            if g.get('skipped'):
                continue
            aln = check_path(os.path.join(p_dir, g['dir'], 'aln.fasta'),
                             must_exist=True)
            regions = conserved_regions_and_consensus(aln)
            if logger:
                logger.log(f"组 {g['group']}: 保守区 {len(regions)} 个")
            for k, (cons, s, e) in enumerate(regions[:5]):
                pairs = design_primers_for_seq(f"{g['group']}_C{k + 1}", cons,
                                               num_return=num_return)
                for p in pairs:
                    p['target'] = f"{g['group']} 保守区{k + 1}(比对列 {s + 1}-{e})"
                all_pairs.extend(pairs)
            if regions:
                groups_used.append(g['group'])
    else:
        a_dir = check_path(assembly_dir or os.path.join(sample_dir, '03_assembly'),
                           must_exist=True)
        vfa = os.path.join(a_dir, 'viral_contigs.fasta')
        if not os.path.isfile(vfa):
            from .orf import extract_viral_contigs
            vfa, _ = extract_viral_contigs(a_dir)
        for h, s in iter_fasta(vfa):
            cid = h.split()[0]
            # 分窗设计覆盖全长
            win, step_len = 3000, 2500
            idx = 1
            for start in range(0, max(len(s) - 400, 1), step_len):
                sub = s[start:start + win]
                if len(sub) < 600:
                    break
                pairs = design_primers_for_seq(f"{cid}_W{idx}", sub,
                                               num_return=max(1, num_return - 1))
                for p in pairs:
                    p['target'] = f"{cid} {start + 1}-{min(start + win, len(s))}bp"
                all_pairs.extend(pairs)
                idx += 1

    if not all_pairs:
        msg = "未设计出引物（保守区不足或序列问题）"
        if logger:
            logger.log(msg, "WARN")
        summary = {'stage': step, 'n_primers': 0, 'reason': msg}
        with safe_open(summary_file, 'wt') as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)
        mark_step_done(out_dir, step)
        return summary

    # 特异性检查（可选）
    spec_hits = {}
    if do_specificity:
        amp_fa = os.path.join(out_dir, 'amplicons.fa')
        with safe_open(amp_fa, 'wt') as f:
            for p in all_pairs:
                # 近似扩增子：模板未知时用引物对连接串占位仅检查引物区不可行；
                # 此处保存引物左臂+右臂反向互补的串联做 BLAST 粗筛
                from Bio.Seq import Seq as BSeq
                seq = p['F_seq'] + str(BSeq(p['R_seq']).reverse_complement())
                from .utils import write_fasta_record
                write_fasta_record(f, p['name'], seq)
        try:
            spec_hits = check_specificity_blast(amp_fa, logger=logger)
        except Exception as e:
            if logger:
                logger.log(f"特异性检查失败: {e}", "WARN")

    tsv = os.path.join(out_dir, 'primers.tsv')
    with safe_open(tsv, 'wt') as f:
        f.write("pair\ttarget\tF_primer(5'-3')\tF_pos\tF_Tm\tF_GC(%)\t"
                "R_primer(5'-3')\tR_pos\tR_Tm\tR_GC(%)\tproduct(bp)\t"
                "penalty\tself_dimer_Th\tself_hairpin_Th\thost_hit\n")
        for p in all_pairs:
            host_hit = 'NA'
            if do_specificity:
                host_hit = 'YES' if p['name'] in spec_hits else 'NO'
            f.write(f"{p['name']}\t{p.get('target', '')}\t{p['F_seq']}\t{p['F_pos']}\t"
                    f"{p['F_tm']}\t{p['F_gc']}\t{p['R_seq']}\t{p['R_pos']}\t"
                    f"{p['R_tm']}\t{p['R_gc']}\t{p['product']}\t{p['penalty']}\t"
                    f"{p['self_any_max']}\t{p['hairpin_max']}\t{host_hit}\n")

    summary = {'stage': step, 'mode': mode, 'n_primers': len(all_pairs),
               'tsv': 'primers.tsv',
               'specificity_checked': bool(do_specificity),
               'groups': groups_used}
    with safe_open(summary_file, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    mark_step_done(out_dir, step)
    if logger:
        logger.log(f"阶段⑥ 完成: {len(all_pairs)} 对引物 -> primers.tsv")
    return summary

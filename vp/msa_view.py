# -*- coding: utf-8 -*-
"""
MSA 网页查看数据提取（借鉴 PhyloSuite MSAViewer 的 SNP-only 设计）。

比对 fasta → 只保留变异位点的展示矩阵：
- positions: 变异列号（1-based，对应原比对）
- rows:      每条序列在变异列上的字符（SNP-only 行）
- consensus: 每个变异列的多数碱基
- diversity: 每个变异列的多样性 (1 - 最高频率)，前端画变异度柱
前端按列分块渲染彩色字符热图（ACGT 经典配色）。
"""
from .utils import iter_fasta

MAX_SEQS = 80          # 展示序列数上限（防巨型比对卡死浏览器）
_GAP = '-'


def _is_variable(col_chars):
    """列内忽略大小写后是否多于一种字符。"""
    upper = {c.upper() for c in col_chars}
    return len(upper) > 1


def snp_view_data(aln_fasta, max_seqs=MAX_SEQS):
    """读取比对 fasta → SNP-only 展示数据 dict。

    返回 {names, rows, positions, consensus, diversity, aln_len,
          n_seq, n_snp, seqs_truncated}
    """
    names, seqs = [], []
    for h, s in iter_fasta(aln_fasta):
        names.append(h.split()[0])
        seqs.append(s.rstrip().upper())
        if len(seqs) >= max_seqs:
            break
    if len(seqs) < 2:
        raise ValueError('比对至少需要 2 条序列')
    aln_len = min(len(s) for s in seqs)

    positions, consensus, diversity = [], [], []
    rows = ['' for _ in seqs]
    n_seq = len(seqs)
    for i in range(aln_len):
        col = [s[i] for s in seqs]
        if not _is_variable(col):
            continue
        counts = {}
        for c in col:
            counts[c] = counts.get(c, 0) + 1
        top_char, top_n = max(counts.items(), key=lambda kv: kv[1])
        positions.append(i + 1)                    # 1-based
        consensus.append(top_char)
        diversity.append(round(1.0 - top_n / n_seq, 3))
        for k, c in enumerate(col):
            rows[k] += c

    return {'names': names, 'rows': rows, 'positions': positions,
            'consensus': consensus, 'diversity': diversity,
            'aln_len': aln_len, 'n_seq': n_seq, 'n_snp': len(positions),
            'seqs_truncated': len(names) >= max_seqs}

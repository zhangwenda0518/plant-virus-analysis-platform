# -*- coding: utf-8 -*-
"""构造合成阳性测试数据：植物病毒基因组模拟 reads + 真实样本背景。

用法: python tests/make_synthetic.py
输出: tests/syn_R1.fastq.gz / syn_R2.fastq.gz
"""
import os
import sys
import random

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.utils import iter_fasta, safe_open, write_fasta_record

random.seed(42)

# 选取的病毒（覆盖不同科）：NC accession -> 期望名
WANTED = {
    'NC_004103.2': 'Tobacco curly shoot virus',       # Geminiviridae (ssDNA 环状)
    'NC_001367.1': 'Tobacco mosaic virus',            # Virgaviridae (ssRNA 线性)
    'NC_003530.1': 'Pepper mild mottle virus',        # Virgaviridae
}

READ_LEN = 150
INSERT_MEAN, INSERT_SD = 300, 60
N_PAIRS_PER_GENOME = 8000
ERROR_RATE = 0.001

COMP = str.maketrans('ACGTN', 'TGCAN')


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate(s):
    if random.random() >= ERROR_RATE * len(s):
        return s
    l = list(s)
    for i in range(len(l)):
        if random.random() < ERROR_RATE:
            l[i] = random.choice('ACGT')
    return ''.join(l)


def fake_qual():
    # Phred+33 区间 Q22-Q39，贴近真实 Illumina 输出；绝不含空格（SPAdes
    # 解析器会 strip 行首空格）。下限须 < ASCII 59：SPAdes 判定 offset 要求
    # 质量字符中出现 < 59 的字符（Q27-39 全 ≥ 60 会报 "Failed to determine
    # offset"）；同时保持 Q≥20（fastp Q20 阈值）不被质控全滤。
    # 注意：Q0-20 的旧合成数据会被 fastp(Q20 阈值) 全量过滤、被 SPAdes
    # 判为超低质量全丢弃——真实数据无此问题，仅测试数据需高 Q。
    return ''.join(chr(random.randint(55, 72)) for _ in range(READ_LEN))


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ref_fa = os.path.join(root, 'virus-db', 'final.cluster.ref.fasta')
    out_dir = os.path.join(root, 'tests')
    os.makedirs(out_dir, exist_ok=True)
    out1 = os.path.join(out_dir, 'syn_R1.fastq.gz')
    out2 = os.path.join(out_dir, 'syn_R2.fastq.gz')

    genomes = {}
    for h, s in iter_fasta(ref_fa):
        acc = h.split()[0]
        if acc in WANTED and acc not in genomes:
            genomes[acc] = s.upper()
    print('找到病毒基因组:', {a: len(s) for a, s in genomes.items()})

    n = 0
    with safe_open(out1, 'wt') as w1, safe_open(out2, 'wt') as w2:
        for acc, g in genomes.items():
            for i in range(N_PAIRS_PER_GENOME):
                ins = max(READ_LEN * 2, int(random.gauss(INSERT_MEAN, INSERT_SD)))
                if len(g) <= ins:            # 环状/短基因组：随机起点环绕
                    start = random.randint(0, len(g) - 1)
                    frag = (g[start:] + g[:start])[:ins] if ins <= len(g) else g
                else:
                    start = random.randint(0, len(g) - ins)
                    frag = g[start:start + ins]
                r1 = mutate(frag[:READ_LEN])
                r2 = mutate(revcomp(frag[-READ_LEN:]))
                name = f'SYN_{acc}_{i}'
                w1.write(f'@{name} 1:N:0:synthetic\n{r1}\n+\n{fake_qual()}\n')
                w2.write(f'@{name} 2:N:0:synthetic\n{r2}\n+\n{fake_qual()}\n')
                n += 1
    print(f'生成 {n:,} 对模拟病毒 reads -> tests/syn_R*.fastq.gz')


if __name__ == '__main__':
    main()

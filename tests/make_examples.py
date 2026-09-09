# -*- coding: utf-8 -*-
"""生成内置示例数据（databases/examples/）。

产物（全部确定性生成，重跑幂等）：
- example_viral_contigs.fasta  已随平台内置（2 条完整植物病毒基因组，勿覆盖）
- example_virus_set.fasta      6 条烟草花叶病毒属近缘基因组（3 条真实参考 +
                               3 条受控突变衍生株），供 MSA / 结构比较 /
                               SDT / NT+AA 同一性表做示例
- example_tree.nwk             上集 MAFFT+FastTree 建的示例树（进化树查看器）
- example_genome.gb            示例 contig 1 的 GenBank（pyrodigal 基因坐标 +
                               占位 product），供基因组图谱 GenBank 模式
- example_synteny_A/B/C.gb     3 条同属小基因组（6 个同源 CDS，受控分化），
                               供「同属共线性比较」离线示例

用法: python tests/make_examples.py
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from vp.config import PLATFORM_ROOT, get_config
from vp.utils import iter_fasta, safe_open, write_fasta_record

EX_DIR = os.path.join(PLATFORM_ROOT, 'databases', 'examples')
REF_FA = os.path.join(PLATFORM_ROOT, 'databases', 'ncbi_refs',
                      '_smoke_tmv', 'refs.fa')
EX_CONTIGS = os.path.join(EX_DIR, 'example_viral_contigs.fasta')

AA = 'ACDEFGHIKLMNPQRSTVWY'
BASES = 'ACGT'
COMP = str.maketrans('ACGT', 'TGCA')


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate_seq(s, rate):
    """按替换率做确定性点突变（不引入 indel，保持长度）。"""
    out = list(s)
    for i in range(len(out)):
        if random.random() < rate:
            out[i] = random.choice(BASES)
    return ''.join(out).upper()


def make_virus_set():
    """6 条 tobamovirus 近缘基因组：3 条真实参考 + 3 条受控突变衍生。"""
    refs = list(iter_fasta(REF_FA))
    assert len(refs) == 3, f'基准参考缺失: {REF_FA}'
    random.seed(7)
    recs = [(h.split()[0], s.upper()) for h, s in refs]
    for i, rate in enumerate((0.02, 0.08, 0.15), 1):
        base_h, base_s = refs[i % 3]
        recs.append((f'{base_h.split()[0]}_variant{i}',
                     mutate_seq(base_s, rate)))
    out = os.path.join(EX_DIR, 'example_virus_set.fasta')
    with safe_open(out, 'wt') as f:
        for h, s in recs:
            write_fasta_record(f, h, s)
    print('ok example_virus_set.fasta:', [(h, len(s)) for h, s in recs])
    return out


def make_conserved_set():
    """5 条同种近缘序列（同一参考 0.1%~1.2% 受控突变）——保守区引物示例。

    保守区引物要求全序列 ≥90% 一致的连续列 ≥400（vp/primer.min_len/min_ident），
    突变率须压到同种分离株水平（~1%）才能形成足够长的连续保守段；
    跨物种的 example_virus_set.fasta 达不到，故另备近缘集。
    """
    refs = list(iter_fasta(REF_FA))
    base_h, base_s = refs[0]
    base_h = base_h.split()[0]
    random.seed(13)
    recs = [(base_h, base_s.upper())]
    for i, rate in enumerate((0.001, 0.003, 0.006, 0.012), 1):
        recs.append((f'{base_h}_isolate{i}', mutate_seq(base_s, rate)))
    out = os.path.join(EX_DIR, 'example_conserved_set.fasta')
    with safe_open(out, 'wt') as f:
        for h, s in recs:
            write_fasta_record(f, h, s)
    print('ok example_conserved_set.fasta:', [(h, len(s)) for h, s in recs])
    return out


def make_tree(set_fa):
    """示例集 → MAFFT 比对 → FastTree 示例树（复用平台链路）。"""
    from vp.phylo import _run_mafft, _run_fasttree

    class _Log:
        def log(self, m, level='INFO'):
            pass

        def close(self):
            pass

    work = os.path.join(EX_DIR, '_work')
    os.makedirs(work, exist_ok=True)
    aln = _run_mafft(set_fa, os.path.join(work, 'aln.fasta'), logger=_Log())
    nwk_path = _run_fasttree(aln, os.path.join(work, 'tree.nwk'), logger=_Log())
    with safe_open(nwk_path) as f:
        nwk = f.read().strip()
    out = os.path.join(EX_DIR, 'example_tree.nwk')
    with safe_open(out, 'wt') as f:
        f.write(nwk + '\n')
    print('ok example_tree.nwk:', nwk[:70], '...')
    return out


def make_genbank():
    """示例 contig 1 → GenBank（pyrodigal 基因坐标 + 分级 product 占位）。"""
    import numpy as np
    from pyrodigal import GeneFinder

    recs = list(iter_fasta(EX_CONTIGS))
    head, seq = recs[0]
    acc = head.split()[0]
    gf = GeneFinder(meta='meta', closed=True)
    genes = gf.find_genes(seq.encode())

    rec = SeqRecord(Seq(seq), id=acc, name=acc,
                    description='plant virus complete genome (built-in example)')
    rec.annotations['organism'] = 'Plant virus (built-in example)'
    rec.annotations['molecule_type'] = 'DNA'
    rec.annotations['topology'] = 'linear'
    rec.annotations['data_file_division'] = 'VRL'
    rec.annotations['date'] = '01-JAN-2026'
    rec.annotations['accessions'] = [acc]
    rec.annotations['comment'] = ('Built-in platform example: gene coordinates '
                                  'predicted by pyrodigal (meta mode); products '
                                  'are placeholders.')
    n_cds = 0
    for g in genes:
        strand = 1 if g.strand == 1 else -1
        loc = FeatureLocation(int(g.begin) - 1, int(g.end), strand=strand)
        f = SeqFeature(loc, type='CDS')
        f.qualifiers = {
            'locus_tag': [f'example_{n_cds + 1:02d}'],
            'product': [f'hypothetical protein {n_cds + 1}'],
            'transl_table': ['11'],
        }
        rec.features.append(f)
        n_cds += 1
    out = os.path.join(EX_DIR, 'example_genome.gb')
    SeqIO.write(rec, out, 'genbank')
    print(f'ok example_genome.gb: {len(seq)}bp, {n_cds} CDS')


def make_synteny_trio():
    """3 条同属小基因组（6 个同源 CDS，受控分化）供共线性比较离线示例。"""
    random.seed(21)
    # 6 个同源基因（~280aa），产物名贴近 Potexvirus 口径
    products = ['RNA-dependent RNA polymerase', 'Triple gene block 1 protein',
                'Triple gene block 2 protein', 'Triple gene block 3 protein',
                'Coat protein', 'Nucleic acid-binding protein']
    base = [''.join(random.choice(AA) for _ in range(280)) for _ in products]

    def derive(rate):
        out = []
        for p in base:
            p = list(p)
            for _ in range(int(len(p) * rate)):
                j = random.randrange(len(p))
                p[j] = random.choice(AA)
            out.append(''.join(p))
        return out

    spec = [('example_synteny_A', 'Potexvirus example A', derive(0.00)),
            ('example_synteny_B', 'Potexvirus example B', derive(0.08)),
            ('example_synteny_C', 'Potexvirus example C', derive(0.30))]
    for acc, org, prots in spec:
        rec = SeqRecord(Seq('N' * (300 + len(prots) * 500)), id=acc, name=acc,
                        description=f'{org}, complete genome (built-in example)')
        rec.annotations['organism'] = org
        rec.annotations['molecule_type'] = 'DNA'
        rec.annotations['topology'] = 'linear'
        rec.annotations['data_file_division'] = 'VRL'
        rec.annotations['date'] = '01-JAN-2026'
        pos = 150
        for i, p in enumerate(prots):
            f = SeqFeature(FeatureLocation(pos, pos + 280 * 3, strand=1),
                           type='CDS')
            f.qualifiers = {'product': [products[i]], 'gene': [f'ORF{i + 1}'],
                            'translation': [p]}
            rec.features.append(f)
            pos += 500
        SeqIO.write(rec, os.path.join(EX_DIR, f'{acc}.gb'), 'genbank')
    print('ok example_synteny_A/B/C.gb: 3 records x 6 CDS')


def main():
    os.makedirs(EX_DIR, exist_ok=True)
    set_fa = make_virus_set()
    make_conserved_set()
    make_tree(set_fa)
    make_genbank()
    make_synteny_trio()
    print('ALL EXAMPLES GENERATED ->', EX_DIR)


if __name__ == '__main__':
    main()

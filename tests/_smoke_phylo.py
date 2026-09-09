# -*- coding: utf-8 -*-
"""冒烟测试：MAFFT -> trimAl -> FastTree / IQ-TREE3 链路（合成小比对）。

合成序列用固定算术规则生成（非随机），保证可复现且无随机性告警。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.config import PLATFORM_ROOT, get_config  # noqa: E402
from vp.phylo import (_run_mafft, _run_trimal, _run_fasttree, _run_iqtree,  # noqa: E402
                      _parse_iqtree_log)
from vp.utils import check_path, write_fasta_record, safe_open  # noqa: E402

BASES = 'ACGT'


def mutate(base_seq, k):
    """第 k 个变体：按固定步长替换 ~80 个位点，确定性生成。"""
    s = list(base_seq)
    n = len(s)
    for j in range(80):
        p = (j * 97 + k * 311) % n
        s[p] = BASES[(BASES.index(s[p]) + 1 + (j + k) % 3) % 4]
    return ''.join(s)


work = check_path(os.path.join(PLATFORM_ROOT, 'results', '_smoke_phylo'),
                  in_platform=True)
os.makedirs(work, exist_ok=True)

base = ''.join(BASES[(i * i + i // 3) % 4] for i in range(1500))
with safe_open(os.path.join(work, 'seqs.fa'), 'wt') as f:
    write_fasta_record(f, 'WT_ref', base)
    for i in range(1, 6):
        write_fasta_record(f, f'variant_{i}', mutate(base, i))

cfg = get_config()
print('tools:', {k: cfg.tools.get(k) for k in ('mafft', 'trimal', 'gblocks',
                                              'fasttree', 'iqtree2', 'iqtree3')})

aln = check_path(os.path.join(work, 'aln.fasta'))
_run_mafft(check_path(os.path.join(work, 'seqs.fa')), aln, threads=4)
with safe_open(aln) as f:
    print('MAFFT OK, first seq line length:',
          len(f.readlines()[1].strip()))

aln_trim = check_path(os.path.join(work, 'aln.trim.fasta'))
aln_used, trim_info = _run_trimal(aln, aln_trim)
print('trimAl:', trim_info, '->', os.path.basename(aln_used))

t1 = check_path(os.path.join(work, 'tree.nwk'))
_run_fasttree(aln_used, t1)
with safe_open(t1) as f:
    print('FastTree OK:', f.read()[:80], '...')

pf = check_path(os.path.join(work, 'iqtree'))
nwk = _run_iqtree(aln_used, pf, threads=4)
with safe_open(nwk) as f:
    print('IQ-TREE OK:', f.read()[:80], '...')
print('iqtree log parse:', _parse_iqtree_log(pf + '.iqtree'))
print('ALL SMOKE TESTS PASSED')

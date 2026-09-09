# -*- coding: utf-8 -*-
"""MSA 查看模块自测：snp_view_data 单元验证 + Flask API 冒烟。"""
import os
import sys
import json
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.config import DIRS  # noqa: E402
from vp.msa_view import snp_view_data  # noqa: E402
from vp.utils import check_path, safe_open, write_fasta_record  # noqa: E402

BASES = 'ACGT'


def check(cond, msg):
    print(('  ok ' if cond else '  FAIL ') + msg)
    assert cond, msg


work = check_path(os.path.join(DIRS['results'], 'itmsa'), in_platform=True)
shutil.rmtree(work, ignore_errors=True)
gdir = os.path.join(work, '05_phylo', 'grpA')
os.makedirs(gdir, exist_ok=True)

# 长度 20 比对：位点 5/11/17（A↔G）与 6/12/18（C/T/gap）为变异列，其余全同
seqs = {
    's1': 'AAAAACAAAAACAAAAACAA',
    's2': 'AAAAGCAAAAGCAAAAACAA',
    's3': 'AAAAGTAAAAGTAAAA--AA',
    's4': 'AAAAACAAAAGCAAAAACAA',
}
aln = os.path.join(gdir, 'aln.fasta')
with safe_open(aln, 'wt') as f:
    for k, s in seqs.items():
        write_fasta_record(f, k, s)

d = snp_view_data(aln)
check(d['positions'] == [5, 6, 11, 12, 17, 18], f'变异列 = {d["positions"]}')
check(d['n_snp'] == 6 and d['aln_len'] == 20, '计数正确')
check(d['rows'][0] == 'ACACAC', f's1 SNP 行 = {d["rows"][0]!r}')
check(d['rows'][2] == 'GTGT--', f's3 SNP 行 = {d["rows"][2]!r}')
check(d['consensus'] == ['A', 'C', 'G', 'C', 'A', 'C'], f'共识 = {d["consensus"]}')  # 位点5 A/G 平票取先出现
check(d['diversity'][0] == 0.5, f'位点5 多样性 = {d["diversity"][0]}')
check(d['diversity'][1] == 0.25, f'位点6 多样性 = {d["diversity"][1]}')

# 单序列应报错
aln1 = os.path.join(work, 'one.fa')
with safe_open(aln1, 'wt') as f:
    write_fasta_record(f, 'only', 'AAAA')
try:
    snp_view_data(aln1)
    check(False, '单序列应抛 ValueError')
except ValueError:
    check(True, '单序列正确抛错')

# ---- API 冒烟（Flask test client）----
import app as appmod  # noqa: E402
c = appmod.app.test_client()
r = c.get('/api/msa/samples')
check(r.status_code == 200, f'/api/msa/samples 200')
items = r.get_json()
check(any(x['sample'] == 'itmsa' for x in items), '样品列表含 itmsa')

r = c.get('/api/msa/data?sample=itmsa&group=grpA')
check(r.status_code == 200, '/api/msa/data 200')
check(r.get_json()['n_snp'] == 6, 'API 数据正确')

r = c.get('/api/msa/data?sample=itmsa&group=nope')
check(r.status_code in (400, 404), '不存在的组返回 4xx')

r = c.get('/api/msa/data?sample=../etc&group=x')
check(r.status_code == 400, '路径注入返回 400')

shutil.rmtree(work, ignore_errors=True)
print('MSA VIEW TESTS PASSED')

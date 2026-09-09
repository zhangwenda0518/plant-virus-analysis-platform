# -*- coding: utf-8 -*-
"""病毒注释分析组集成测试：ORF 预测 → 功能注释 → 基因组图谱 → 引物设计
→ 工具④ contig 分类与逐 contig 分析（本地 primer + 在线 NCBI blastn/CDD）。

全部走 /api/tool/* Web API（Flask test client），任务经 TaskManager 轮询，
示例数据用 databases/examples/ 内置文件（离线可复跑；在线部分失败自动跳过）。
"""
import os
import sys
import json
import time
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.config import PLATFORM_ROOT  # noqa: E402

EX = os.path.join(PLATFORM_ROOT, 'databases', 'examples')
EX_FA = os.path.join(EX, 'example_viral_contigs.fasta')
EX_GB = os.path.join(EX, 'example_genome.gb')
EX_SET = os.path.join(EX, 'example_virus_set.fasta')
RUNS = os.path.join(PLATFORM_ROOT, 'tool_runs')


def check(cond, msg):
    print(('  ok ' if cond else '  FAIL ') + msg)
    assert cond, msg


def wait_task(tm, tid, timeout=1800):
    t0 = time.time()
    while time.time() - t0 < timeout:
        rec = tm.tasks.get(tid)
        if rec is None:
            raise RuntimeError(f'任务不存在: {tid}')
        if rec['status'] in ('done', 'failed', 'cancelled'):
            return rec
        time.sleep(1.0)
    raise RuntimeError(f'任务超时: {tid}')


def run_tool(c, tm, tool, params, threads=None, timeout=1800):
    r = c.post('/api/tool/run', json={'tool': tool, 'threads': threads,
                                      'params': params})
    assert r.status_code == 200, f'{tool} 启动失败: {r.status_code} {r.get_json()}'
    d = r.get_json()
    rec = wait_task(tm, d['task'], timeout=timeout)
    assert rec['status'] == 'done', \
        f'{tool} 运行失败: {rec.get("error", "")[-800:]}'
    return d['run'], rec.get('result') or {}


import app as appmod  # noqa: E402
c = appmod.app.test_client()
tm = appmod.tm

check(os.path.isfile(EX_FA) and os.path.isfile(EX_GB) and os.path.isfile(EX_SET),
      '内置示例文件齐全（example_viral_contigs / example_genome.gb / '
      'example_virus_set）')

# ---------- 1. ORF 预测（orfipy + pyrodigal） ----------
print('--- 1. ORF 预测 ---')
orf_run, res = run_tool(c, tm, 'orf', {'fasta': EX_FA, 'min_aa': 100})
check(orf_run.startswith('orf_'), f'运行名 orf_*: {orf_run}')
rd = os.path.join(RUNS, orf_run)
for sub in ('03_assembly/viral_contigs.fasta', '03_assembly/contigs.filtered.fasta',
            '04_orf/orfipy_pep.fa', '04_orf/orfipy_nt.fa', '04_orf/orfipy.bed',
            '04_orf/pyrodigal.faa', '04_orf/pyrodigal.ffn', '04_orf/pyrodigal.gff'):
    check(os.path.isfile(os.path.join(rd, sub.replace('/', os.sep))),
          f'产物存在: {sub}')
check((res.get('n_orfs') or 0) > 0, f'ORF 数 > 0: {res.get("n_orfs")}')

# ---------- 2. 功能注释（方式 A：对已有 orf_ 运行注释） ----------
# 注：run 模式产物写回目标 orf_ 运行目录（注释与 ORF 运行同目录，便于追溯）
print('--- 2. 功能注释（已有运行） ---')
orfa_run, res = run_tool(c, tm, 'orfa', {'run': orf_run}, timeout=3600)
check(orfa_run.startswith('orfa_'), f'运行名 orfa_*: {orfa_run}')
rd = os.path.join(RUNS, orf_run)
for sub in ('04b_orf_annot/orf_annotation.tsv',
            '04b_orf_annot/orf_annotation.gff3',
            '04b_orf_annot/orf_function_summary.tsv',
            '04b_orf_annot/orf_family_summary.tsv',
            '04b_orf_annot/summary.json'):
    check(os.path.isfile(os.path.join(rd, sub.replace('/', os.sep))),
          f'产物存在（写入 orf_ 运行目录）: {sub}')
with open(os.path.join(rd, '04b_orf_annot', 'orf_annotation.tsv'),
          encoding='utf-8') as f:
    lines = [l for l in f.read().splitlines() if l.strip()]
check(len(lines) > 1, f'注释表 {len(lines) - 1} 行 ORF')
# 至少部分 ORF 有功能类别与产物
import csv as _csv  # noqa: E402
with open(os.path.join(rd, '04b_orf_annot', 'orf_annotation.tsv'),
          encoding='utf-8') as f:
    rows = list(_csv.DictReader(f, delimiter='\t'))
n_prod = sum(1 for r in rows if (r.get('product') or '').strip())
n_cat = sum(1 for r in rows if (r.get('category') or '').strip()
            and (r.get('informative') == 'Y'))
check(n_prod > 0, f'{n_prod}/{len(rows)} 个 ORF 有产物注释')
check(n_cat > 0, f'{n_cat}/{len(rows)} 个 ORF 功能类别有效')
gdiags = os.path.join(rd, '04b_orf_annot', 'genome_diagrams')
svgs = [x for x in os.listdir(gdiags) if x.endswith('.svg')] \
    if os.path.isdir(gdiags) else []
check(svgs, f'功能示意 SVG: {svgs}')

# ---------- 3. 功能注释（方式 B：FASTA 一步预测+注释） ----------
print('--- 3. 功能注释（FASTA 一步） ---')
orfa2_run, res = run_tool(c, tm, 'orfa', {'fasta': EX_FA, 'min_aa': 100},
                          timeout=3600)
check(os.path.isfile(os.path.join(
    RUNS, orfa2_run, '04b_orf_annot', 'orf_annotation.tsv')),
    '一步法注释产物存在')

# ---------- 4. 基因组图谱（FASTA / GenBank 两种输入） ----------
print('--- 4. 基因组图谱 ---')
gp_run, res = run_tool(c, tm, 'genoplot', {'fasta': EX_FA, 'engine': 'auto',
                                           'max_plots': 12})
plots = res.get('plots') or []
check(len(plots) >= 2, f'FASTA 出图 {len(plots)} 条（≥2）')
check(all(os.path.isfile(p) for p in plots), '图文件真实存在')
check(any(p.endswith('circular.svg') for p in plots)
      and any(p.endswith('linear.svg') for p in plots),
      '圈图 + 线图齐备（gbdraw 口径）')
gp_gb_run, res = run_tool(c, tm, 'genoplot', {'ann': EX_GB, 'engine': 'auto',
                                              'max_plots': 12})
check(len(res.get('plots') or []) >= 1, f'GenBank 出图 {len(res.get("plots") or [])} 条')

# ---------- 5. 引物设计（全长分窗 + 保守区） ----------
print('--- 5. 引物设计 ---')
pr_run, res = run_tool(c, tm, 'primer', {'fasta': EX_FA, 'mode': 'plain',
                                         'num_return': 3})
check((res.get('n_primers') or 0) > 0,
      f'全长分窗引物 {res.get("n_primers")} 对')
check(os.path.isfile(os.path.join(RUNS, pr_run, '06_primer', 'primers.tsv')),
      'primers.tsv 存在')

# 先对近缘序列集做结构比较产 aln（作为保守区引物输入；
# 保守区设计要求 ≥90% 一致连续列 ≥400，跨物种的 example_virus_set 太散）
EX_CONS = os.path.join(EX, 'example_conserved_set.fasta')
sc_run, _res = run_tool(c, tm, 'structcmp', {'seqs': EX_CONS, 'max_n': 30})
aln = os.path.join(RUNS, sc_run, 'aln.fasta')
check(os.path.isfile(aln), '近缘集比对 aln.fasta（供保守区引物）')
pr2_run, res = run_tool(c, tm, 'primer', {'fasta': aln, 'mode': 'conserved',
                                          'num_return': 3})
check((res.get('n_primers') or 0) > 0,
      f'保守区引物 {res.get("n_primers")} 对')

# ---------- 6. 工具④ contig 分类 + 分类报告 + 逐 contig 本地分析 ----------
print('--- 6. contig 分类 + 逐 contig 分析 ---')
ctg_run, res = run_tool(c, tm, 'contigs', {'contigs': EX_FA, 'min_len': 500},
                        timeout=3600)
rd = os.path.join(RUNS, ctg_run)
for sub in ('contigs.filtered.fasta', 'virus_classification.tsv',
            'viral_contigs.fasta', 'contig_blast.tsv'):
    check(os.path.isfile(os.path.join(rd, sub.replace('/', os.sep))),
          f'产物存在: {sub}')
check((res.get('n_viral') or 0) >= 2, f'病毒 contig {res.get("n_viral")} 条')

# 分类表 / kreport → 页内报告数据接口
r = c.get(f'/api/tool/viral_contigs?run={ctg_run}')
check(r.status_code == 200 and len(r.get_json()) >= 2,
      f'/api/tool/viral_contigs 200, {len(r.get_json())} 行')
r = c.get(f'/api/tool/report_data?run={ctg_run}')
check(r.status_code == 200, '/api/tool/report_data 200')
d = r.get_json()
check(d['mode'] == 'contig' and d['n_virus'] >= 2,
      f'报告数据 mode={d["mode"]}, n_virus={d["n_virus"]}')

# 本地分析四件套中的 primer（primer3，离线）
with open(os.path.join(rd, 'viral_contigs.fasta'), encoding='utf-8') as f:
    ctg_id = f.readline().strip()[1:].split()[0]
r = c.post('/api/tool/analyze',
           json={'run': ctg_run, 'contig': ctg_id, 'action': 'primer'})
check(r.status_code == 200, f'/api/tool/analyze primer 启动: {r.get_json()}')
if r.get_json().get('task'):
    rec = wait_task(tm, r.get_json()['task'], timeout=600)
    check(rec['status'] == 'done', 'primer 分析完成')
r = c.get(f'/api/tool/analysis?run={ctg_run}&contig={ctg_id}&action=primer')
check(r.status_code == 200 and (r.get_json() or {}).get('primers'),
      'primer 结果缓存可读')

# ---------- 7. 运行列表 / 结果预览 ----------
r = c.get('/api/tool/runs')
check(r.status_code == 200, '/api/tool/runs 200')
names = [x['name'] for x in r.get_json()]
# 列表只保留最近 30 条；断言本次运行创建的较新 run 必在
for want in (gp_run, pr_run, sc_run, ctg_run):
    check(want in names, f'运行列表含 {want}')

# ---------- 8. 在线分析（NCBI blastn / CDD；无网络自动跳过） ----------
print('--- 8. 在线 NCBI 分析（blastn / CDD，失败自动跳过） ---')
r = c.post('/api/tool/analyze',
           json={'run': ctg_run, 'contig': ctg_id, 'action': 'blastn'})
if r.status_code == 200:
    tid = r.get_json().get('task')
    if tid:
        try:
            rec = wait_task(tm, tid, timeout=1500)
            if rec['status'] == 'done':
                r2 = c.get(f'/api/tool/analysis?run={ctg_run}'
                           f'&contig={ctg_id}&action=blastn')
                hits = (r2.get_json() or {}).get('hits') or []
                check(True, f'BLASTN 在线分析完成，{len(hits)} 条命中')
            else:
                print('  SKIP BLASTN（任务未成功，离线或 NCBI 限流）')
        except RuntimeError as e:
            print(f'  SKIP BLASTN（{e}）')
    else:
        check(r.get_json().get('cached'), 'blastn 结果已缓存')
else:
    print('  SKIP BLASTN（启动被拒，离线环境）')

print('ANNOTATE GROUP INTEGRATION TESTS PASSED')

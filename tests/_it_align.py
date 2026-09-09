# -*- coding: utf-8 -*-
"""比对模块 + ICTV 级联 + 基因级建树 集成测试。

覆盖：① ICTV 界门纲目科属种级联 API；② 序列比对工具（MAFFT+trimAl）
与查看/编辑保存 API；③ GenBank 集合 CDS/PEP 提取建树（coat protein）。
"""
import os
import sys
import time
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.config import PLATFORM_ROOT  # noqa: E402

EX = os.path.join(PLATFORM_ROOT, 'databases', 'examples')
EX_SET = os.path.join(EX, 'example_virus_set.fasta')
EX_SYNTENY = [os.path.join(EX, f'example_synteny_{x}.gb') for x in 'ABC']
RUNS = os.path.join(PLATFORM_ROOT, 'tool_runs')
COLL = 'it_align_coll'


def check(cond, msg):
    print(('  ok ' if cond else '  FAIL ') + msg, flush=True)
    assert cond, msg


def wait_task(tm, tid, timeout=1800):
    t0 = time.time()
    while time.time() - t0 < timeout:
        rec = tm.tasks.get(tid)
        if rec is None:
            raise RuntimeError(f'任务不存在: {tid}')
        if rec['status'] in ('done', 'failed', 'cancelled'):
            return rec
        time.sleep(1.5)
    raise RuntimeError(f'任务超时: {tid}')


def main():
    import app as appmod  # noqa: E402
    c = appmod.app.test_client()
    tm = appmod.tm

    def run_tool(tool, params, threads=None, timeout=1800):
        r = c.post('/api/tool/run', json={'tool': tool, 'threads': threads,
                                          'params': params})
        assert r.status_code == 200, \
            f'{tool} 启动失败: {r.status_code} {r.get_json()}'
        d = r.get_json()
        rec = wait_task(tm, d['task'], timeout=timeout)
        assert rec['status'] == 'done', \
            f'{tool} 运行失败: {rec.get("error", "")[-600:]}'
        return d['run'], rec.get('result') or {}

    # ---------- 1. ICTV 级联 ----------
    print('--- 1. ICTV 级联 ---', flush=True)
    r = c.get('/api/ictv/cascade')
    check(r.status_code == 200, '/api/ictv/cascade 200')
    d = r.get_json()
    ranks = {x['col']: x for x in d['ranks']}
    check(len(ranks['Family']['options']) >= 55, '科选项 ≥55')
    check(len(ranks['Genus']['options']) >= 230, '属选项 ≥230')
    check(any(o['name'] == 'Riboviria' and o['n'] > 2000
              for o in ranks['Realm']['options']), '界 Riboviria 计数 >2000')
    # 下钻：Riboviria → Tombusviridae
    r = c.get('/api/ictv/cascade?Realm=Riboviria&Family=Tombusviridae')
    d = r.get_json()
    g = next(x for x in d['ranks'] if x['col'] == 'Genus')
    sp = next(x for x in d['ranks'] if x['col'] == 'Species')
    check(len(g['options']) == 19, f'Tombusviridae 下 19 属（实得 {len(g["options"])}）')
    check(len(sp['options']) == 94, f'Tombusviridae 下 94 种（实得 {len(sp["options"])}）')
    # 预览（多级 levels + 植物口径）
    r = c.post('/api/ictv/preview',
               json={'Realm': 'Riboviria', 'Family': 'Tombusviridae',
                     'Genus': 'Dianthovirus', 'genome': 'complete',
                     'limit': 20})
    d = r.get_json()
    check(d['scope'] == 'Genus=Dianthovirus' and d['total'] == 6,
          f"级联预览 Dianthovirus total={d['total']}")
    check(all('plant' in (x['host'] or '').lower() for x in d['rows']),
          '预览行全部为植物病毒（plant 口径）')

    # ---------- 2. 序列比对工具 + 查看/编辑 ----------
    print('--- 2. 序列比对 + 查看器 + 编辑 ---', flush=True)
    al_run, res = run_tool('align', {'seqs': EX_SET, 'strategy': 'auto',
                                     'trimal': 'automated1', 'max_n': 30})
    check(res.get('n_seqs') == 6 and res.get('seqtype') == 'nt',
          f"比对 6 条 nt（{res.get('n_seqs')}/{res.get('seqtype')}）")
    rd = os.path.join(RUNS, al_run)
    check(os.path.isfile(os.path.join(rd, 'aln.fasta')), 'aln.fasta 落盘')
    check(os.path.isfile(os.path.join(rd, 'aln.trim.fasta')),
          'aln.trim.fasta 落盘（automated1）')
    rel = f'tool_runs/{al_run}/aln.trim.fasta'
    r = c.get(f'/api/align/file?path={rel}')
    check(r.status_code == 200, '/api/align/file 200')
    d = r.get_json()
    check(d['aligned'] and d['n'] == 6 and d['cols'] > 500,
          f"查看器数据 aligned={d['aligned']} n={d['n']} cols={d['cols']}")
    # 编辑：删掉一条 → 保存副本
    keep = '\n'.join(f'>{n}\n{s}' for n, s in
                     list(zip(d['names'], d['seqs']))[:5])
    r = c.post('/api/align/save', json={'path': rel, 'content': keep})
    check(r.status_code == 200 and r.get_json()['n'] == 5,
          '编辑保存副本（5 条）')
    edited = r.get_json()['saved']
    check(edited.endswith('.edit.fasta') and os.path.isfile(edited),
          f'副本落盘: {os.path.basename(edited)}')
    r = c.get('/api/align/file?path=' + edited.replace('\\', '/'))
    check(r.status_code == 200 and r.get_json()['n'] == 5, '副本可再打开')
    # 坏内容被拒
    r = c.post('/api/align/save', json={'path': rel, 'content': 'not fasta'})
    check(r.status_code == 400, '非法内容保存被拒')
    r = c.get('/api/align/file?path=../platform.json')
    check(r.status_code == 400, '查看器路径穿越被拒')

    # ---------- 3. CDS/PEP 基因级建树 ----------
    print('--- 3. 基因级建树（CDS/PEP）---', flush=True)
    from vp.gb_collection import gb_collection_dir  # noqa: E402
    shutil.rmtree(gb_collection_dir(COLL), ignore_errors=True)
    r = c.post('/api/gb/import', json={'name': COLL, 'files': EX_SYNTENY})
    tid = r.get_json().get('task')
    if tid:
        rec = wait_task(tm, tid, timeout=600)
        check(rec['status'] == 'done', '示例集合导入完成')

    r = c.post('/api/gb/phylo', json={'name': COLL, 'tree_tool': 'fasttree',
                                      'molecule': 'pep', 'gene': 'coat protein'})
    check(r.status_code == 200, 'PEP 建树启动')
    rec = wait_task(tm, r.get_json()['task'], timeout=1200)
    check(rec['status'] == 'done',
          f"PEP 建树完成: {rec.get('error', '')[-300:]}")
    gdir = os.path.join(gb_collection_dir(COLL), 'gene_trees',
                        'coat_protein_pep')
    check(os.path.isfile(os.path.join(gdir, 'tree.nwk')), 'PEP tree.nwk')
    with open(os.path.join(gdir, 'tree.nwk')) as f:
        nwk = f.read().strip()
    check(nwk.endswith(';') and nwk.count(':') >= 2, 'PEP Newick 合法')
    check(os.path.isfile(os.path.join(gdir, 'combined.fa')), 'PEP combined.fa')
    with open(os.path.join(gdir, 'combined.fa')) as f:
        heads = [l.strip() for l in f if l.startswith('>')]
    check(len(heads) == 3, f'PEP 提取 3 条（{len(heads)}）')

    r = c.post('/api/gb/phylo', json={'name': COLL, 'tree_tool': 'fasttree',
                                      'molecule': 'cds', 'gene': 'coat protein'})
    rec = wait_task(tm, r.get_json()['task'], timeout=1200)
    check(rec['status'] == 'done', f"CDS 建树完成: {rec.get('error', '')[-200:]}")
    gdir2 = os.path.join(gb_collection_dir(COLL), 'gene_trees',
                         'coat_protein_cds')
    check(os.path.isfile(os.path.join(gdir2, 'tree.nwk')), 'CDS tree.nwk')

    # 参数守卫
    r = c.post('/api/gb/phylo', json={'name': COLL, 'molecule': 'pep'})
    check(r.status_code == 400, 'PEP 无关键词被拒')
    r = c.post('/api/gb/phylo', json={'name': COLL, 'molecule': 'pep',
                                      'gene': 'nonexistent_gene_xyz'})
    rec = wait_task(tm, r.get_json()['task'], timeout=600)
    check(rec['status'] == 'failed' and '命中' in (rec.get('error') or ''),
          '无命中基因报错清晰')

    # 清理
    shutil.rmtree(gb_collection_dir(COLL), ignore_errors=True)
    print('ALIGN MODULE INTEGRATION TESTS PASSED', flush=True)


if __name__ == '__main__':
    main()

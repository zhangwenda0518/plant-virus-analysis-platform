# -*- coding: utf-8 -*-
"""比较基因组分析组集成测试：结构比较 → MSA 查看 → SDT 精确矩阵 →
NT+AA 同一性表 → 快速建树 → 树文件 API → NCBI 在线（可跳过）→
GenBank 集合示例导入 → 同属共线性比较 → 集合建树。

全部走 Web API（Flask test client），示例数据用 databases/examples/ 内置文件。
注意：sdt/identity 工具内部用 ProcessPoolExecutor，Windows spawn 会重导入
__main__——测试体必须收进 main() + __name__ 保护，否则子进程递归崩溃。
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.config import PLATFORM_ROOT  # noqa: E402

EX = os.path.join(PLATFORM_ROOT, 'databases', 'examples')
EX_SET = os.path.join(EX, 'example_virus_set.fasta')
EX_TREE = os.path.join(EX, 'example_tree.nwk')
EX_SYNTENY = [os.path.join(EX, f'example_synteny_{x}.gb')
              for x in 'ABC']
RUNS = os.path.join(PLATFORM_ROOT, 'tool_runs')
COLL = 'it_compare_coll'


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
        time.sleep(1.0)
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
            f'{tool} 运行失败: {rec.get("error", "")[-800:]}'
        return d['run'], rec.get('result') or {}

    check(all(os.path.isfile(p) for p in [EX_SET, EX_TREE] + EX_SYNTENY),
          '内置示例文件齐全（example_virus_set / example_tree.nwk / '
          'example_synteny_A-C.gb）')

    # ---------- 1. 结构比较（identity 矩阵 / 热图数据） ----------
    print('--- 1. 结构比较 ---', flush=True)
    sc_run, res = run_tool('structcmp', {'seqs': EX_SET, 'max_n': 30})
    check(res.get('n_seqs') == 6, f"参与比较 {res.get('n_seqs')} 条（6）")
    check(os.path.isfile(os.path.join(RUNS, sc_run, 'identity_matrix.json'))
          and os.path.isfile(os.path.join(RUNS, sc_run, 'identity_matrix.tsv')),
          'identity_matrix.json / .tsv 落盘')
    r = c.get(f'/api/tool/structcmp_data?run={sc_run}')
    check(r.status_code == 200, '/api/tool/structcmp_data 200')
    d = r.get_json()
    check(d['n'] == 6 and len(d['matrix']) == 6, '矩阵 6×6')
    diag_ok = all(abs(d['matrix'][i][i] - 100.0) < 1e-9 for i in range(6))
    check(diag_ok, '对角线 = 100')
    # 变种与母本应高度一致（>92%），跨参考应低一些
    names, m = d['names'], d['matrix']
    i_v1 = names.index('NC_002692_variant1')
    i_b = names.index('NC_002692')
    check(m[i_v1][i_b] > 92, f'variant1 vs 母本 identity {m[i_v1][i_b]:.1f}% > 92')
    i_t1, i_t2 = names.index('NC_001367'), names.index('NC_002692')
    check(m[i_t1][i_t2] < m[i_v1][i_b] - 3,
          f'跨参考 identity {m[i_t1][i_t2]:.1f}% 明显低于近缘对')

    # ---------- 2. MSA 查看（structcmp 运行 → SNP 数据） ----------
    print('--- 2. MSA 查看 ---', flush=True)
    r = c.get(f'/api/tool/msa_data?run={sc_run}')
    check(r.status_code == 200, '/api/tool/msa_data 200')
    d = r.get_json()
    check(d.get('n_seq') == 6 and d.get('n_snp', 0) > 20,
          f"SNP 展示数据 n_seq={d.get('n_seq')}, n_snp={d.get('n_snp')}")

    # ---------- 3. SDT 精确矩阵（NT） ----------
    print('--- 3. SDT 精确矩阵（NT） ---', flush=True)
    sdt_run, res = run_tool('sdt', {'seqs': EX_SET, 'max_n': 30,
                                    'seqtype': 'nt', 'orient': True},
                            timeout=3600)
    rd = os.path.join(RUNS, sdt_run)
    for fn in ('sdt_matrix.csv', 'sdt_matrix.json', 'sdt_heatmap.png',
               'sdt_distribution.png'):
        check(os.path.isfile(os.path.join(rd, fn)), f'产物存在: {fn}')
    check(res.get('pairs') == 15 and res.get('n') == 6,
          f"6 选 2 = {res.get('pairs')} 对")
    r = c.get(f'/tool_runs/{sdt_run}/sdt_matrix.json')
    check(r.status_code == 200, 'sdt_matrix.json 可经 /tool_runs 预览')
    d = r.get_json()
    check(all(abs(d['matrix'][i][i] - 100.0) < 1e-9
              for i in range(len(d['matrix']))), 'SDT 矩阵对角线 = 100')
    # 近缘对（变种 vs 母本）的 SDT 相似度应高于跨参考对
    nm = d.get('names') or []
    j = nm.index('NC_002692_variant1')
    k = nm.index('NC_002692')
    j2 = nm.index('NC_001367')
    check(d['matrix'][j][k] > d['matrix'][j2][k],
          f"SDT 近缘对 {d['matrix'][j][k]:.1f}% > 跨参考对 "
          f"{d['matrix'][j2][k]:.1f}%")

    # ---------- 4. NT+AA 同一性表（identity 工具） ----------
    print('--- 4. NT+AA 同一性表 ---', flush=True)
    idt_run, res = run_tool('identity', {'nt_seqs': EX_SET, 'max_n': 6},
                            timeout=3600)
    rd = os.path.join(RUNS, idt_run)
    for fn in ('identity.json', 'identity_table.csv', 'nt_matrix.csv',
               'aa_matrix.csv'):
        check(os.path.isfile(os.path.join(rd, fn)), f'产物存在: {fn}')
    r = c.get(f'/tool_runs/{idt_run}/identity.json')
    check(r.status_code == 200, 'identity.json 可预览')
    d = r.get_json()
    check(d.get('n') == 6, f"identity.json n = {d.get('n')}")

    # ---------- 5. 快速建树（NJ + FastTree）+ 树文件 API ----------
    print('--- 5. 快速建树 / 树文件 ---', flush=True)
    qt_run, res = run_tool('quicktree', {'seqs': EX_SET, 'method': 'nj',
                                         'max_n': 100})
    check(res.get('tree') == 'nj.nwk', f"NJ 树文件 = {res.get('tree')}")
    r = c.get(f'/api/tool/quicktree_data?run={qt_run}')
    check(r.status_code == 200 and r.get_json()['newick'].endswith(';'),
          '/api/tool/quicktree_data 返回 Newick')
    qt2_run, res = run_tool('quicktree', {'seqs': EX_SET,
                                          'method': 'fasttree',
                                          'max_n': 100})
    check(res.get('tree') == 'tree.nwk', f"FastTree 树文件 = {res.get('tree')}")
    r = c.get(f'/api/tool/quicktree_data?run={qt2_run}')
    d = r.get_json()
    check(d.get('tool') == 'FastTree', f"树摘要 tool = {d.get('tool')}")

    # 本机树文件 API（内置示例树 + 非法路径）
    r = c.get(f'/api/tree/file?path=databases/examples/example_tree.nwk')
    check(r.status_code == 200 and r.get_json()['newick'].endswith(';'),
          '/api/tree/file 读内置示例树')
    r = c.get('/api/tree/file?path=../platform.json')
    check(r.status_code == 400, '路径穿越被拒')

    # ---------- 6. NCBI 在线（无网络自动跳过） ----------
    print('--- 6. NCBI 在线（可跳过） ---', flush=True)
    r = c.post('/api/ncbi/search',
               json={'term': 'Tobamovirus[ORGN] AND complete genome[TITL]',
                     'db': 'nucleotide', 'limit': 10})
    if r.status_code == 200:
        d = r.get_json()
        check(d.get('count', 0) > 50 and len(d.get('rows', [])) > 0,
              f"Entrez 检索命中 {d.get('count')} 条"
              f"（预览 {len(d.get('rows', []))}）")
    else:
        print(f'  SKIP Entrez 检索（{r.status_code}，离线环境）')

    # ---------- 7. GenBank 集合示例（导入 → 共线性 → 集合建树） ----------
    print('--- 7. GenBank 集合示例（离线全流程） ---', flush=True)
    import shutil  # noqa: E402
    from vp.gb_collection import gb_collection_dir  # noqa: E402
    shutil.rmtree(gb_collection_dir(COLL), ignore_errors=True)

    r = c.post('/api/gb/import', json={'name': COLL,
                                       'files': EX_SYNTENY})
    check(r.status_code == 200, f'/api/gb/import 启动: {r.get_json()}')
    tid = r.get_json().get('task')
    if tid:
        rec = wait_task(tm, tid, timeout=600)
        check(rec['status'] == 'done', '导入完成')

    r = c.get('/api/gb/collections')
    cols = r.get_json()
    mine = next((x for x in cols if x['name'] == COLL), None)
    check(mine is not None and mine['n_records'] == 3,
          f"集合列表含 {COLL}（{mine and mine['n_records']} 条记录）")

    r = c.post('/api/compare/run', json={'collection': COLL, 'min_ident': 0.3,
                                         'min_cov': 0.5, 'style': 'lovis'})
    check(r.status_code == 200, '/api/compare/run 启动')
    tid = r.get_json().get('task')
    rec = wait_task(tm, tid, timeout=1800)
    check(rec['status'] == 'done',
          f'共线性比较完成: {rec.get("error", "")[-300:]}')
    files = (rec.get('result') or {}).get('files') or {}
    for k in ('html', 'svg', 'png', 'clusters', 'similarity', 'faa', 'm8'):
        p = files.get(k)
        check(p and os.path.isfile(p),
              f'共线性产物存在: {k}'
              f'（{os.path.basename(p) if p else "缺失"}）')

    # 集合建树 → 结果中心 gb: 伪样品可看
    r = c.post('/api/gb/phylo', json={'name': COLL, 'tree_tool': 'fasttree'})
    check(r.status_code == 200, '/api/gb/phylo 启动')
    tid = r.get_json().get('task')
    rec = wait_task(tm, tid, timeout=1800)
    check(rec['status'] == 'done',
          f'集合建树完成: {rec.get("error", "")[-300:]}')

    r = c.get('/api/msa/samples')
    gb_sample = next((x for x in r.get_json()
                      if x['sample'] == f'gb:{COLL}'), None)
    check(gb_sample is not None, '结果中心出现 gb: 伪样品')
    check(bool(gb_sample['groups'][0]['trees']), '伪样品分组含树文件')
    r = c.get(f'/api/tree/data?sample=gb:{COLL}&group={COLL}')
    check(r.status_code == 200 and r.get_json()['newick'].endswith(';'),
          '/api/tree/data (gb:) 返回 Newick')
    r = c.get(f'/api/sdt/data?sample=gb:{COLL}&group={COLL}')
    check(r.status_code == 200, '/api/sdt/data (gb:) 200')

    print('COMPARE GROUP INTEGRATION TESTS PASSED', flush=True)


if __name__ == '__main__':
    main()

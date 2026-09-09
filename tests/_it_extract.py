# -*- coding: utf-8 -*-
"""CDS/PEP 提取（PhyloSuite 布局）+ ICTV 每属/每种抽样上限 集成测试。"""
import os
import sys
import time
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.config import PLATFORM_ROOT  # noqa: E402

EX = os.path.join(PLATFORM_ROOT, 'databases', 'examples')
EX_SYNTENY = [os.path.join(EX, f'example_synteny_{x}.gb') for x in 'ABC']
COLL = 'it_extract_coll'


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

    # ---------- 1. ICTV 抽样上限 ----------
    print('--- 1. ICTV 每属/每种抽样 ---', flush=True)
    r = c.post('/api/ictv/preview',
               json={'Family': 'Tombusviridae', 'genome': 'any', 'limit': 300})
    d0 = r.get_json()
    r = c.post('/api/ictv/preview',
               json={'Family': 'Tombusviridae', 'genome': 'any', 'limit': 300,
                     'per_genus': 5})
    d5 = r.get_json()
    check(d5['n_shown'] < d0['n_shown'],
          f"每属≤5 收窄：{d0['n_shown']} → {d5['n_shown']}")
    check(d5['sampling'] == '每属≤5', '预览标注抽样口径')
    r = c.post('/api/ictv/preview',
               json={'Family': 'Tombusviridae', 'genome': 'any', 'limit': 300,
                     'per_genus': 5, 'per_species': 1})
    d51 = r.get_json()
    check(d51['n_shown'] <= d5['n_shown'],
          f"加每种≤1 进一步收窄：{d5['n_shown']} → {d51['n_shown']}")
    # 大属场景：Potyvirus 科级下载配 per_genus
    from vp import ictv_db
    rows, total = ictv_db.select_refs(family='Potyviridae', genome='any',
                                      limit=5000)
    kept, dropped = ictv_db.cap_per_rank(rows, per_genus=5)
    from collections import Counter
    gen_cnt = Counter((r_.get('Genus') or '').strip().lower() for r_ in kept)
    check(max(gen_cnt.values()) <= 5 and dropped == total - len(kept),
          f'Potyviridae 每属≤5：{total} → {len(kept)}（max/属 '
          f'{max(gen_cnt.values())}）')

    # ---------- 2. CDS/PEP 提取 ----------
    print('--- 2. CDS/PEP 提取（PhyloSuite 布局）---', flush=True)
    from vp.gb_collection import gb_collection_dir, extract_dir  # noqa: E402
    shutil.rmtree(gb_collection_dir(COLL), ignore_errors=True)
    r = c.post('/api/gb/import', json={'name': COLL, 'files': EX_SYNTENY})
    tid = r.get_json().get('task')
    if tid:
        rec = wait_task(tm, tid, timeout=600)
        check(rec['status'] == 'done', '示例集合导入')

    r = c.post('/api/gb/extract', json={'name': COLL})
    check(r.status_code == 200, '提取任务启动')
    rec = wait_task(tm, r.get_json()['task'], timeout=1200)
    check(rec['status'] == 'done',
          f"提取完成: {rec.get('error', '')[-300:]}")
    res = rec.get('result') or {}
    check(res.get('n_genome') == 3 and res.get('n_cds') == 18
          and res.get('n_pep') == 18,
          f"3 基因组 × 6 CDS（{res.get('n_genome')}/{res.get('n_cds')}/"
          f"{res.get('n_pep')}）")
    check(res.get('n_genes') == 6, f"6 个基因分目录（{res.get('n_genes')}）")

    exd = extract_dir(COLL)
    # PhyloSuite 布局断言
    for rel in ('genome.fa', 'CDS.fa', 'PEP.fa', 'genes.tsv'):
        check(os.path.isfile(os.path.join(exd, rel)), f'extract/{rel}')
    for sub in ('CDS', 'PEP'):
        d_ = os.path.join(exd, sub)
        fas = sorted(f for f in os.listdir(d_) if f.endswith('.fa'))
        check(len(fas) == 6, f'extract/{sub}/ 按基因 6 个文件（{len(fas)}）')
    # 按基因文件内容：coat protein 每基因组一份
    cds_cp = os.path.join(exd, 'CDS', 'ORF5.fa')
    if not os.path.isfile(cds_cp):
        # 基因名以集合注释为准（ORF1..ORF6）；退化为任一基因文件抽查
        d_ = os.path.join(exd, 'CDS')
        cds_cp = os.path.join(d_, sorted(os.listdir(d_))[0])
    with open(cds_cp) as f:
        n_seq = sum(1 for l in f if l.startswith('>'))
    check(n_seq == 3, f'按基因文件含 3 个基因组的同源序列（{n_seq}）')
    with open(os.path.join(exd, 'genes.tsv')) as f:
        rows_tsv = f.read().splitlines()
    check(len(rows_tsv) == 7, f'genes.tsv 6 基因 + 表头（{len(rows_tsv)} 行）')

    # 产物清单 API
    r = c.get(f'/api/gb/extract_files?name={COLL}')
    check(r.status_code == 200 and len(r.get_json()) >= 15,
          f'extract_files 清单（{len(r.get_json())} 项）')
    r = c.get('/api/gb/extract_files?name=no_such_coll')
    check(r.status_code == 404, '未提取集合 404')
    # 集合列表标记
    r = c.get('/api/gb/collections')
    mine = next(x for x in r.get_json() if x['name'] == COLL)
    check(mine.get('has_extract') is True, '集合列表 has_extract 标记')

    shutil.rmtree(gb_collection_dir(COLL), ignore_errors=True)
    print('EXTRACT + SAMPLING TESTS PASSED', flush=True)


if __name__ == '__main__':
    main()

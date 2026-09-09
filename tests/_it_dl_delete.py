# -*- coding: utf-8 -*-
"""下载批次三模式删除测试：删文件（留记录）/ 删记录（留文件）/ 全删。"""
import io
import os
import sys
import json
import time
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.config import PLATFORM_ROOT  # noqa: E402

DL = os.path.join(PLATFORM_ROOT, 'downloads')


def check(cond, msg):
    print(('  ok ' if cond else '  FAIL ') + msg, flush=True)
    assert cond, msg


def make_batch(m, bid, files=('a_R1.fastq.gz', 'a_R2.fastq.gz')):
    bdir = os.path.join(DL, bid)
    shutil.rmtree(bdir, ignore_errors=True)
    os.makedirs(bdir, exist_ok=True)
    for f in files:
        with open(os.path.join(bdir, f), 'wb') as fh:
            fh.write(b'FAKEFASTQ' * 100)
    with open(os.path.join(bdir, 'extra.txt'), 'w') as fh:
        fh.write('x')
    b = {'id': bid, 'name': 'it-del', 'status': 'completed',
         'created': time.strftime('%Y-%m-%d %H:%M'), 'dir': bdir,
         'concurrency': 1, 'converting': False, 'convert_msg': '',
         'error': '', 'runs': {},
         'files': [{'acc': 'SRRTEST', 'url': 'https://example.org/' + f,
                    'out': os.path.join(bdir, f), 'size': 900,
                    'done_bytes': 900, 'progress': 1.0, 'status': 'done',
                    'md5': '', 'verified': True, 'error': ''}
                   for f in files]}
    with m.lock:
        m.batches[bid] = b
    m._save(b)
    return bdir


def main():
    from vp import public_data
    m = public_data.get_manager()

    # ---------- ① 删文件（留记录） ----------
    print('--- ① delete_files ---', flush=True)
    bdir = make_batch(m, 'it_delfiles_test')
    check(m.delete_files('it_delfiles_test'), 'delete_files 返回 True')
    left = set(os.listdir(bdir))
    check(left == {'batch.json', 'batch.log'},
          f'目录只剩记录+日志: {sorted(left)}')
    b = m.batches.get('it_delfiles_test')
    check(b is not None and b['files_deleted']
          and all(fe['status'] == 'deleted' for fe in b['files']),
          '记录保留，文件条目标记 deleted')
    with io.open(os.path.join(bdir, 'batch.log'), encoding='utf-8') as f:
        log = f.read()
    check('已删除' in log and '保留批次记录' in log, '日志写入删除记录')
    # 列表仍显示
    snaps = m.list_snapshots()
    check(any(s['id'] == 'it_delfiles_test' for s in snaps), '列表仍显示该批次')

    # ---------- ② 删记录（留文件） ----------
    print('--- ② delete_record ---', flush=True)
    bdir2 = make_batch(m, 'it_delrecord_test')
    check(m.delete_record('it_delrecord_test'), 'delete_record 返回 True')
    left2 = set(os.listdir(bdir2))
    check('batch.json' not in left2 and 'batch.log' not in left2,
          'batch.json/batch.log 已删')
    check({'a_R1.fastq.gz', 'a_R2.fastq.gz', 'extra.txt'} <= left2,
          f'数据文件保留: {sorted(left2)}')
    check(m.batches.get('it_delrecord_test') is None, '内存记录移除')
    check(not any(s['id'] == 'it_delrecord_test'
                  for s in m.list_snapshots()), '列表不再显示')

    # ---------- ③ 全删 ----------
    print('--- ③ delete（全删）---', flush=True)
    bdir3 = make_batch(m, 'it_delall_test')
    check(m.delete('it_delall_test'), 'delete 返回 True')
    check(not os.path.isdir(bdir3), '目录整体移除')
    check(m.batches.get('it_delall_test') is None, '记录移除')

    # ---------- ④ Web API 链路 ----------
    print('--- ④ API ---', flush=True)
    import app as appmod  # noqa: E402
    c = appmod.app.test_client()
    bdir4 = make_batch(public_data.get_manager(), 'it_delapi_test')
    r = c.post('/api/dl/batch/it_delapi_test/delete_files')
    check(r.status_code == 200 and
          set(os.listdir(bdir4)) == {'batch.json', 'batch.log'},
          'API delete_files 生效')
    r = c.post('/api/dl/batch/it_delapi_test/delete_record')
    check(r.status_code == 200 and not os.path.isfile(
        os.path.join(bdir4, 'batch.json')), 'API delete_record 生效')
    r = c.post('/api/dl/batch/it_delapi_test/delete')
    check(r.status_code == 400, '记录已不存在 → 400')
    r = c.post('/api/dl/batch/xxx/bad_action')
    check(r.status_code == 400, '无效动作 400')
    shutil.rmtree(bdir2, ignore_errors=True)   # ② 留下的孤儿文件清理
    print('DL DELETE MODES TESTS PASSED', flush=True)


if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
"""LOGAN 批量自动提交测试：用假提交脚本替代真实 selenium 全链路验证。

覆盖：任务准备/映射、进度事件（PROGRESS 行 → 平台进度回调）、
SSE 实时进度流、结果回收导入、补漏场景、参数校验。"""
import os
import sys
import json as _json
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp import logan_trace as lt

SAMPLE = 'REGRESS'
JOB = '_test_batch'

# 假 logan_submit：解析输入 FASTA，逐条输出 PROGRESS 进度事件并写入
# 固定结果表（模拟提交+下载全流程）。文件名安全化 + 守卫校验与生产脚本同口径。
FAKE_SCRIPT = r'''
import os, sys, time, json
from pathlib import Path
args = sys.argv[1:]
inp, out = None, None
i = 0
while i < len(args):
    a = args[i]
    if a == "-o":
        out = Path(args[i + 1]); i += 2
    elif a == "--email":
        assert "@" in args[i + 1]; i += 2
    elif a.startswith("--"):
        i += 2
    else:
        inp = Path(a); i += 1
time.sleep(0.2)
out.mkdir(parents=True, exist_ok=True)
base = os.path.normpath(os.path.abspath(str(out)))


def safe_stem(acc):
    s = "".join(c if c.isalnum() or c in "._-" else "_" for c in str(acc))
    s = ".".join(seg for seg in s.split(".") if seg)
    return s or "query"


def out_path(stem):
    p = os.path.normpath(os.path.join(base, f"{stem}.tsv"))
    if os.path.isabs(stem) or not p.startswith(base + os.sep):
        raise SystemExit("非法结果文件名")
    return p


def prog(**kw):
    print("PROGRESS " + json.dumps(kw, ensure_ascii=False), flush=True)


tsv = ("acc\tsample_acc\torganism\tkmer_coverage\tANI_estimation\tassay_type\n"
       "DRR111\tSAM111\tNicotiana tabacum\t0.91\t98.5\tRNA-Seq\n"
       "DRR222\tSAM222\tChenopodium quinoa\t0.88\t97.1\tWGS\n")
accs = [ln[1:].strip() for ln in inp.read_text(encoding="utf-8").splitlines()
        if ln.startswith(">")]
for n, acc in enumerate(accs):
    prog(stage="submit", acc=acc, done=n, total=len(accs), msg=f"正在提交 {acc}")
    prog(stage="submitted", acc=acc, done=n, total=len(accs),
         sid=f"kmviz-{n:08d}-0000", msg="session ok")
    prog(stage="waiting", acc=acc, done=n, total=len(accs), waited=30,
         http=400, span=2100)
    prog(stage="downloading", acc=acc, done=n, total=len(accs))
    Path(out_path(safe_stem(acc))).write_text(tsv, encoding="utf-8")
    prog(stage="downloaded", acc=acc, done=n, total=len(accs))
    prog(stage="segment_done", acc=acc, done=n + 1, total=len(accs))
print("[fake] 提交完成")
'''


def check(cond, msg):
    print(('  ✔ ' if cond else '  ✘ ') + msg)
    assert cond, msg


def main():
    fake = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        '_fake_logan_submit.py')
    Path(fake).write_text(FAKE_SCRIPT, encoding='utf-8')
    lt.BATCH_SCRIPT = fake
    try:
        try:
            lt.delete_job(JOB)
        except Exception:
            pass

        print('== 批量准备 ==')
        d = lt.create_job(JOB, sample=SAMPLE, n_seg=2)
        check(d['n_segments'] == 2, f'查询任务 2 个片段（{d["name"]}）')
        in_path, out_dir, acc2seg, n = lt.batch_prepare(d['name'])
        check(n == 2 and len(acc2seg) == 2,
              f'batch_prepare: {n} 个待提交, 映射 {sorted(acc2seg)}')
        with open(in_path, encoding='utf-8') as f:
            txt = f.read()
        heads = [ln[1:] for ln in txt.splitlines() if ln.startswith('>')]
        check(sorted(heads) == sorted(acc2seg), f'输入 FASTA 头正确: {heads}')

        print('== 进度事件抓取 ==')
        events = []
        r = lt.batch_submit(d['name'], emails=['a@b.c'], group='Fast',
                            first_wait=0, max_wait=60,
                            progress=lambda st, pct, msg: events.append((st, pct, msg)))
        check(r['n_imported'] == 2, f"回收导入 {r['n_imported']}/2")
        prog_events = [e for e in events if e[2] and ('等待服务器结果' in e[2]
                                                     or '正在提交' in e[2]
                                                     or '下载' in e[2])]
        check(len(events) >= 6, f'收到 {len(events)} 次进度回调（含最终完成）')
        check(len(prog_events) >= 6, f'其中结构化进度 {len(prog_events)} 次')
        pcts = [p for _, p, _ in events]
        check(pcts == sorted(pcts), f'进度单调不回退: {[round(p, 2) for p in pcts]}')
        check(events[-1][1] == 1.0 and '回收' in events[-1][2], '最终进度 100%')
        check(any('等待服务器结果' in m and 'HTTP 400' in m for _, _, m in events),
              '等待阶段文案含已等时长与 HTTP 状态')

        print('== 任务状态与报告 ==')
        detail = lt.job_detail(d['name'])
        check(detail['status'] == 'done' and detail['has_report'],
              '全部片段已导入且报告已生成')
        rpt = os.path.join(lt.DIRS['logan'], lt.safe_job_name(d['name']),
                           'trace_report.html')
        with open(rpt, encoding='utf-8') as f:
            check('Nicotiana tabacum' in f.read(), '报告聚合了导入的物种')

        print('== 全部导入后再批量 → 应拒绝 ==')
        try:
            lt.batch_submit(d['name'], emails=['a@b.c'])
            check(False, '应拒绝重复批量')
        except ValueError as e:
            check(True, f'已拒绝: {str(e)[:40]}')

        print('== 补漏场景（仅未导入片段进入批量） ==')
        try:
            lt.delete_job(JOB + '2')
        except Exception:
            pass
        d2 = lt.create_job(JOB + '2', sample=SAMPLE, n_seg=2)
        in_path2, out_dir2, acc2seg2, n2 = lt.batch_prepare(d2['name'])
        accs = sorted(acc2seg2)
        with lt.safe_open(os.path.join(out_dir2, accs[0] + '.tsv'), 'wt') as f:
            f.write('acc\torganism\tkmer_coverage\nDRR1\tNicotiana tabacum\t0.9\n')
        check(lt.batch_collect(d2['name'], {accs[0]: acc2seg2[accs[0]]}) == 1,
              '部分回收 1 片段')
        check(lt.job_detail(d2['name'])['status'] == 'partial', '任务状态 partial')
        _p, _o, m2, nn = lt.batch_prepare(d2['name'])
        check(nn == 1 and list(m2.values()) == [acc2seg2[accs[1]]],
              '补漏批量仅含未导入片段')
        r2 = lt.batch_submit(d2['name'], emails=['a@b.c'])
        check(r2['n_pending'] == 1 and r2['n_imported'] == 1, '补漏批量回收 s2')
        check(lt.job_detail(d2['name'])['status'] == 'done', '补漏后结果齐全')

        print('== SSE 实时进度流 ==')
        import app as flaskapp
        c = flaskapp.app.test_client()
        try:
            lt.delete_job(JOB + '3')
        except Exception:
            pass
        lt.create_job(JOB + '3', sample=SAMPLE, n_seg=1)
        rr = c.post(f'/api/logan/job/{JOB + "3"}/batch', json={'emails': 'a@b.c'})
        check(rr.status_code == 200, '批量任务经 API 启动')
        tid = rr.get_json()['task']
        # 等任务跑完（假脚本 <1s），再连 SSE：应推一条终态快照后自动关流
        time.sleep(1.5)
        resp = c.get(f'/api/task/{tid}/stream')
        buf, got = b'', None
        for chunk in resp.response:
            buf += chunk if isinstance(chunk, bytes) else chunk.encode()
            if b'data:' in buf:
                line = [l for l in buf.split(b'\n') if l.startswith(b'data:')][0]
                got = _json.loads(line[5:].strip())
                break
        resp.close()
        check(got is not None and got.get('id') == tid, 'SSE 推送任务快照')
        check(got.get('status') != 'running' or got.get('pct') is not None,
              f"SSE 快照含状态/进度（{got.get('status')}）")

        print('\n批量模式 + 实时进度全部通过 ✔')
    finally:
        lt.BATCH_SCRIPT = os.path.join(os.path.dirname(
            os.path.abspath(lt.__file__)), 'logan_submit.py')   # 还原
        for nm in (JOB, JOB + '2', JOB + '3'):
            try:
                lt.delete_job(nm)
            except Exception:
                pass
        try:
            os.remove(fake)
        except OSError:
            pass
        print('（测试任务与假脚本已清理）')


if __name__ == '__main__':
    import time
    main()

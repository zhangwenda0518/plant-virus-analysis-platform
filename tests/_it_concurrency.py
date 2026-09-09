# -*- coding: utf-8 -*-
"""并发闸门 / 内存回收 / 列表瘦身 集成验证（TaskManager 层，不启 HTTP）。

跑法：  C:\\Python312\\python.exe tests\\_it_concurrency.py
判定：  全部断言通过打印 "ALL PASS"，否则抛出 AssertionError。

覆盖：
  1. heavy 闸门：同时运行的 heavy 任务不超过 max_heavy_tasks
  2. light 闸门：independent 于 heavy，不超过 max_light_tasks
  3. 排队任务可取消（不会卡在 Semaphore 上）
  4. 取消后名额归还（后续排队任务能继续拿到）
  5. 排队状态对外呈现为 running（前端取消按钮可用）
  6. list_all 瘦身：终态任务日志 <= 5 行，运行中为 30 行
"""
import os
import sys
import time
import threading

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import app  # noqa: E402

tm = app.tm
HEAVY = tm.limits['heavy']
LIGHT = tm.limits['light']
fails = []


def check(cond, msg):
    print(('  OK   ' if cond else '  FAIL ') + msg)
    if not cond:
        fails.append(msg)


def sleep_job(sec):
    def _j(log, prog, cancel):
        for i in range(int(sec * 10)):
            if cancel.is_set():
                return None
            time.sleep(0.1)
        return 'ok'
    return _j


def running_now(status=('running',)):
    return sum(1 for r in tm.tasks.values() if r['status'] in status)


print('闸门上限:', tm.limits)

# ---- 1 & 2 & 5: heavy 闸门 + 排队对外呈现 running ---------------------
print('\n[1] heavy 并发不超过 %d' % HEAVY)
n = HEAVY + 3
tids = [tm.start('T-heavy-%d' % i, sleep_job(1.5), weight='heavy')
        for i in range(n)]
time.sleep(0.8)
peak = 0
for _ in range(15):
    peak = max(peak, running_now())
    time.sleep(0.1)
check(peak <= HEAVY, '并发峰值 %d <= %d' % (peak, HEAVY))
check(peak >= 1, '并发峰值 %d >= 1（确认任务真的在跑）' % peak)

snaps = {t['id']: t for t in tm.list_all()}
queued = [s for s in snaps.values()
          if tm.tasks[s['id']]['status'] == 'queued']
check(len(queued) > 0, '存在排队任务（%d 个）' % len(queued))
check(all(s['status'] == 'running' for s in queued),
      '排队任务对外 status=running（前端取消按钮可用）')
if queued:
    check('排队' in (queued[0].get('msg') or '') or 'Queued' in (queued[0].get('msg') or ''),
          '排队任务 msg 含排队提示: %r' % (queued[0].get('msg') or '')[:40])

# ---- 3 & 4: 排队中取消 + 名额归还 -------------------------------------
print('\n[2] 排队中取消')
if queued:
    victim = queued[0]['id']
    ok = tm.cancel(victim)
    check(ok, 'cancel() 对排队任务返回 True')
    time.sleep(1.2)
    st = tm.tasks[victim]['status']
    check(st == 'cancelled', '排队任务状态 -> cancelled（实际 %s）' % st)

print('\n[3] 所有任务最终都能完成（名额正确归还，无死锁）')
deadline = time.time() + 60
while time.time() < deadline:
    if all(tm.tasks[t]['status'] in ('done', 'cancelled', 'failed')
           for t in tids):
        break
    time.sleep(0.3)
done_n = sum(1 for t in tids if tm.tasks[t]['status'] == 'done')
check(done_n == n - (1 if queued else 0),
      '完成任务数 %d == 预期 %d' % (done_n, n - (1 if queued else 0)))

# ---- light 闸门独立于 heavy -------------------------------------------
print('\n[4] light 闸门（上限 %d）与 heavy 独立' % LIGHT)
lt = [tm.start('T-light-%d' % i, sleep_job(1.2), weight='light')
      for i in range(LIGHT + 2)]
time.sleep(0.8)
lpeak = 0
for _ in range(10):
    lpeak = max(lpeak, sum(1 for r in tm.tasks.values()
                           if r['status'] == 'running'
                           and r.get('weight') == 'light'))
    time.sleep(0.1)
check(lpeak <= LIGHT, 'light 并发峰值 %d <= %d' % (lpeak, LIGHT))
for t in lt:
    tm.cancel(t)
time.sleep(0.5)

# ---- 6: list_all 瘦身 --------------------------------------------------
print('\n[5] list_all 日志瘦身')
t2 = tm.start('T-log', sleep_job(0.3), weight='light')
for i in range(40):
    tm.tasks[t2]['log'].append('line-%d' % i)
snap_running = next(s for s in tm.list_all() if s['id'] == t2)
check(len(snap_running['log']) <= 30,
      '运行中任务日志 %d 行 <= 30' % len(snap_running['log']))
time.sleep(1.5)
snap_done = next(s for s in tm.list_all() if s['id'] == t2)
check(tm.tasks[t2]['status'] == 'done', '任务已完成')
check(len(snap_done['log']) <= 5,
      '终态任务日志 %d 行 <= 5' % len(snap_done['log']))

# ---- 内存回收 ----------------------------------------------------------
print('\n[6] _trim_memory 回收')
tm._trim_memory()
check(len(tm.order) <= app._KEEP_TOTAL,
      'order 长度 %d <= %d' % (len(tm.order), app._KEEP_TOTAL))

print('\n' + '=' * 46)
if fails:
    print('FAILED: %d 项' % len(fails))
    for f in fails:
        print('  - ' + f)
    sys.exit(1)
print('ALL PASS')

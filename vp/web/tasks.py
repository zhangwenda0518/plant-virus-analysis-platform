# -*- coding: utf-8 -*-
"""后台任务引擎（自 app.py 拆出）。

TaskManager：后台线程执行 + 状态快照（内存为主，tasks/ 留一份 JSON）。
并发控制：heavy / light 两档 Semaphore，任务先排队（status='queued'）
再等名额；日志经 SSE（/api/task/<id>/stream）实时推给前端。

内存保留策略：最近 _KEEP_FULL 个任务保留完整信息，更早的清空日志/结果
引用只留状态摘要；超过 _KEEP_TOTAL 的移出内存（磁盘 tasks/*.json 仍可回溯）。

注：原实现用 `app.logger.warning`，拆出后改用本模块 logger，避免对 Flask
应用对象的模块级依赖（后台线程内 app 对象不可靠）。
"""
import json
import logging
import os
import re
import subprocess
import sys
import threading
import time
import uuid
from collections import deque

from vp.config import DIRS, PLATFORM_ROOT
from vp.utils import check_path, fmt_size, safe_open, safe_remove
from vp.web.state import cfg, tool_runs_root as _tool_runs_root

_log = logging.getLogger('vp.web.tasks')


def _slot_limits():
    """并发闸门上限（platform.json defaults 可覆盖：max_heavy_tasks /
    max_light_tasks）；未配置时用缺省值。

    heavy = 吃满多核/GB 级内存的任务（kunpeng 分类、SPAdes 组装、
    DIAMOND 注释、建库、建树、SDT）；light = 网络 IO 或轻量计算
    （下载、格式转换、绘图、检索、Selenium 提交）。
    """
    d = getattr(cfg, 'defaults', None) or {}
    try:
        heavy = max(1, min(int(d.get('max_heavy_tasks', 2)), 8))
    except (TypeError, ValueError):
        heavy = 2
    try:
        light = max(1, min(int(d.get('max_light_tasks', 4)), 16))
    except (TypeError, ValueError):
        light = 4
    return heavy, light


_KEEP_FULL = 60


_KEEP_TOTAL = 300


class TaskManager:
    """后台线程任务 + 状态快照（内存为主，tasks/ 目录留一份 JSON）。

    并发控制：heavy / light 两档 Semaphore。任务先排队（status='queued'，
    对外仍呈现为 running，仅 msg 标注排队位次），拿到名额后才真正开跑。
    """

    def __init__(self):
        self.tasks = {}
        self.order = []
        self.lock = threading.Lock()
        heavy, light = _slot_limits()
        self.heavy_slots = threading.Semaphore(heavy)
        self.light_slots = threading.Semaphore(light)
        self.limits = {'heavy': heavy, 'light': light}
        self._recover_interrupted()

    def _recover_interrupted(self):
        """服务启动时把上次进程遗留的 running 任务标记为中断（防幽灵），
        并把任务状态文件修剪到最近 100 个（防 tasks/ 无限膨胀）。"""
        try:
            files = [n for n in os.listdir(DIRS['tasks']) if n.endswith('.json')]
            files.sort(reverse=True)                     # 新任务 id 排前（uuid 时间无关，
                                                         # 改按 mtime 更稳，见下）
            files.sort(key=lambda n: os.path.getmtime(
                os.path.join(DIRS['tasks'], n)), reverse=True)
            for fn in files[100:]:                       # 只保留最近 100 个
                try:
                    safe_remove(os.path.join(DIRS['tasks'], fn))
                except (OSError, ValueError):
                    pass
            files = files[:100]
            for name in files:
                p = check_path(os.path.join(DIRS['tasks'], name),
                               must_exist=False, in_platform=True)
                try:
                    with safe_open(p) as f:
                        d = json.load(f)
                except (OSError, ValueError):
                    continue
                if d.get('status') in ('running', 'queued'):
                    d['status'] = 'failed'
                    d['error'] = '服务重启导致任务中断，请重新运行'
                    d['finished'] = time.time()
                    with safe_open(p, 'wt') as f:
                        json.dump(d, f, ensure_ascii=False)
        except Exception:
            pass

    def start(self, name, fn, log_file=None, weight='heavy', link=None):
        """创建任务。log_file: 任务日志持久化文件（缺省落 logs/tasks/）。

        weight: 'heavy'（吃多核/大内存，受 heavy 闸门限流）或 'light'
        （网络 IO / 轻量计算，受 light 闸门限流）。缺省 heavy——宁可多排队，
        也不要把不认识的任务放过去打满机器。

        link: 结果/操作页跳转地址（任务卡「前往」按钮）。

        重启上下文自动捕获：任务在 POST API handler 里发起时，记录
        request.path + JSON body，任务卡可一键重新提交同一请求。
        仅存内存（不落盘），避免 API key 等敏感字段写入 tasks/*.json。

        log() 同时写内存 deque（前端实时读）与磁盘文件（重启后可查）。
        """
        tid = uuid.uuid4().hex[:12]
        restart = None
        try:
            from flask import has_request_context, request as _req
            if has_request_context() and _req.method == 'POST':
                body = _req.get_json(silent=True)
                restart = {'url': _req.path, 'body': body} if body else None
        except Exception:
            restart = None
        if not log_file:
            safe = re.sub(r'[^\w\-.]+', '_', name, flags=re.UNICODE)
            safe = safe[:40].strip('_ ') or 'task'
            log_dir = check_path(os.path.join(DIRS['logs'], 'tasks'),
                                 must_exist=False, in_platform=True)
            os.makedirs(log_dir, exist_ok=True)
            log_file = os.path.join(
                log_dir, f'{safe}_{time.strftime("%Y%m%d_%H%M%S")}_{tid}.log')
        rec = {'id': tid, 'name': name, 'status': 'queued', 'stage': '',
               'pct': 0.0, 'msg': '', 'error': '', 'log': deque(maxlen=500),
               'started': time.time(), 'finished': None, 'result': None,
               'result_preview': None, 'eta': None,
               'cancel': threading.Event(),
               'procs': [],          # 本任务启动的子进程（硬停止用）
               'weight': weight if weight in ('heavy', 'light') else 'heavy',
               'log_file': log_file, 'link': link, 'restart': restart}
        with self.lock:
            self.tasks[tid] = rec
            self.order.insert(0, tid)

        def _log(line):
            rec['log'].append(line)
            try:
                with safe_open(rec['log_file'], 'at') as f:
                    f.write(time.strftime('[%H:%M:%S] ') + line + '\n')
            except Exception:
                pass

        def _progress(stage, pct, msg, *eta):
            rec['stage'], rec['pct'], rec['msg'] = stage, pct, msg
            if eta:
                try:
                    rec['eta'] = max(float(eta[0]), 0)
                except (TypeError, ValueError):
                    rec['eta'] = None
            self._persist(rec)

        def _run():
            from vp.utils import task_bind
            slot = (self.light_slots if rec['weight'] == 'light'
                    else self.heavy_slots)
            # 排队等名额：每 0.5s 探一次，期间可被取消（不会被闸门永久卡住）
            got = False
            while not rec['cancel'].is_set():
                if slot.acquire(timeout=0.5):
                    got = True
                    break
                rec['msg'] = self._queue_msg(rec)
            if not got:                       # 排队途中被取消
                rec['status'] = 'cancelled'
                rec['error'] = '任务已停止（排队中取消）'
                rec['finished'] = time.time()
                _log('[CANCEL] ' + cfg.tr('排队中已取消',
                                          'Cancelled while queued'))
                self._persist(rec)
                return
            rec['status'] = 'running'
            rec['started'] = time.time()      # 已运行时长从真正开跑算起
            self._persist(rec)
            task_bind(rec['cancel'], rec['procs'])
            try:
                result = fn(_log, _progress, rec['cancel'])
                rec['status'] = 'cancelled' if rec['cancel'].is_set() else 'done'
                rec['result'] = result
                rec['result_preview'] = build_result_preview(result, rec['name'])
                if rec['status'] == 'done':
                    _log('[OK] ' + cfg.tr('任务完成', 'Task finished'))
            except Exception as e:
                rec['status'] = 'cancelled' if rec['cancel'].is_set() else 'failed'
                rec['error'] = ('任务已停止（用户取消）' if rec['cancel'].is_set()
                                else str(e))
                _log(f"[ERROR] {rec['error']}" if rec['cancel'].is_set()
                     else f"[ERROR] {e}")
            finally:
                task_bind(None, None)
                rec['procs'].clear()
                rec['finished'] = time.time()
                slot.release()          # 归还名额
                rec['msg'] = ''
                self._trim_memory()     # 回收旧任务内存
                self._persist(rec)

        th = threading.Thread(target=_run, daemon=True, name=f'task-{tid}')
        rec['thread'] = th
        th.start()
        self._persist(rec)
        return tid

    def _persist(self, rec):
        try:
            p = check_path(os.path.join(DIRS['tasks'], rec['id'] + '.json'),
                           must_exist=False, in_platform=True)
            data = {'id': rec['id'], 'name': rec['name'],
                    'status': rec['status'], 'stage': rec['stage'],
                    'pct': rec['pct'], 'msg': rec['msg'],
                    'log_file': rec.get('log_file'),
                    'error': rec['error'],
                    'log_tail': list(rec['log'])[-40:],
                    'started': rec['started'], 'finished': rec['finished'],
                    'eta': rec.get('eta'), 'link': rec.get('link')}
            if rec.get('result') is not None:
                try:
                    json.dumps(rec['result'], ensure_ascii=False)
                    data['result'] = rec['result']
                except (TypeError, ValueError):
                    data['result'] = str(rec['result'])
            with safe_open(p, 'wt') as f:
                json.dump(data, f, ensure_ascii=False)
        except Exception as e:
            # 状态落盘失败不影响内存态，但会让刷新后丢进度：留痕便于排查
            _log.warning('任务状态保存失败 %s: %s', rec.get('id'), e)

    def _queue_msg(self, rec):
        """排队提示：统计同档位里排在自己前面的任务数。"""
        ahead = 0
        with self.lock:
            for tid in self.order:
                if tid == rec['id']:
                    break
                r = self.tasks.get(tid)
                if (r and r['status'] == 'queued'
                        and r.get('weight') == rec.get('weight')):
                    ahead += 1
        if ahead:
            return cfg.tr(f'排队中 · 前面还有 {ahead} 个任务',
                          f'Queued · {ahead} task(s) ahead')
        return cfg.tr('排队中 · 等待空闲名额', 'Queued · waiting for a slot')

    def _trim_memory(self):
        """回收内存：超出 _KEEP_FULL 的终态任务清掉日志与结果引用（磁盘已
        持久化，仍可经 tasks/*.json 回溯）；超出 _KEEP_TOTAL 的移出内存。

        服务长跑时任务只增不减会持续吃内存，这里做稳态回收。
        """
        with self.lock:
            for tid in self.order[_KEEP_FULL:]:
                rec = self.tasks.get(tid)
                if not rec or rec['status'] in ('running', 'queued'):
                    continue
                if rec['log']:
                    rec['log'].clear()
                rec['result'] = None
                rec['result_preview'] = None
                rec['thread'] = None
                rec['procs'] = []
            if len(self.order) > _KEEP_TOTAL:
                for tid in self.order[_KEEP_TOTAL:]:
                    self.tasks.pop(tid, None)
                self.order = self.order[:_KEEP_TOTAL]

    def snapshot(self, tid, log_lines=80, include_result=True):
        with self.lock:
            rec = self.tasks.get(tid)
            if not rec:
                return None
            # 'queued' 对外一律呈现为 running：前端据此显示取消按钮、进度条
            # 与自动滚动；排队状态本身通过 msg 文案告知用户。
            st = 'running' if rec['status'] == 'queued' else rec['status']
            out = {'id': rec['id'], 'name': rec['name'],
                   'status': st, 'stage': rec['stage'],
                   'pct': rec['pct'], 'msg': rec['msg'],
                   'error': rec['error'], 'started': rec['started'],
                   'finished': rec['finished'], 'eta': rec.get('eta'),
                   'link': rec.get('link'), 'restart': rec.get('restart'),
                   'log': list(rec['log'])[-log_lines:]}
            if include_result:
                out['result'] = rec.get('result_preview')
        return out

    def list_all(self, limit=_KEEP_FULL * 2, logs=True):
        """任务列表。

        前端每 2.5s 轮询一次本接口，历史上对所有任务都返回 30 行日志 +
        结果预览，任务累积后 payload 持续变大。这里做瘦身：
        - 运行中/排队中的任务：完整 30 行日志（实时看进度）
        - 已结束的任务：只带最后 5 行（日志已落盘，展开任务卡可单独拉取）
        - logs=False：完全不带日志（任务中心/徒轮询只需元信息，
          日志走 /api/task/<tid>/log 单独拉取）
        - 结果预览一律保留（任务完成后的结果面板依赖它）
        """
        with self.lock:
            ids = list(self.order[:limit])
            active = {i for i in ids
                      if (self.tasks.get(i) or {}).get('status')
                      in ('running', 'queued')}
        if not logs:
            return [self.snapshot(i, log_lines=0, include_result=True)
                    for i in ids if i in self.tasks]
        return [self.snapshot(i, log_lines=(30 if i in active else 5),
                              include_result=True)
                for i in ids if i in self.tasks]

    @staticmethod
    def _kill_tree(proc):
        """杀进程树（Windows 用 taskkill /T 连子进程；其它平台 kill）。"""
        if proc.poll() is not None:
            return
        try:
            if sys.platform == 'win32':
                subprocess.run(['taskkill', '/F', '/T', '/PID', str(proc.pid)],
                               capture_output=True, timeout=15)
            else:
                proc.kill()
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass

    def cancel(self, tid):
        with self.lock:
            rec = self.tasks.get(tid)
        if rec and rec['status'] in ('running', 'queued'):
            rec['cancel'].set()
            if rec['status'] == 'queued':
                # 排队中：闸门等待循环会自行退出并归还名额
                rec['msg'] = '已取消排队'
                return True
            # 硬停止：杀掉本任务已注册的全部子进程（kunpeng/spades/mafft 等）
            for proc in list(rec.get('procs') or ()):
                self._kill_tree(proc)
            rec['msg'] = '已请求取消（当前步骤结束后停止）'
            return True
        return False


    @staticmethod
    def _tid_ok(tid):
        return bool(tid) and bool(re.fullmatch(r'[0-9a-f]{6,16}', str(tid)))

    def delete(self, tid):
        """删除任务记录（仅终态）：移出内存 + 删 tasks/<tid>.json。"""
        if not self._tid_ok(tid):
            return False, '非法任务 ID'
        with self.lock:
            rec = self.tasks.get(tid)
            if rec and rec['status'] in ('running', 'queued'):
                return False, '运行中的任务请先停止再删除'
            self.tasks.pop(tid, None)
            if tid in self.order:
                self.order.remove(tid)
        try:
            p = check_path(os.path.join(DIRS['tasks'], tid + '.json'),
                           must_exist=False, in_platform=True)
            if os.path.isfile(p):
                os.remove(p)
        except (OSError, ValueError):
            pass
        return True, 'ok'

    def clear_finished(self):
        """清空内存中的终态任务（含各自的 tasks/<tid>.json）。

        注意：不碰磁盘上的历史归档（平台重启后加载的 tasks/*.json）。
        清空归档是破坏性操作，只能由用户逐个删除。
        """
        with self.lock:
            tids = [t for t, r in self.tasks.items()
                    if r['status'] not in ('running', 'queued')]
        n = 0
        for t in tids:
            ok, _ = self.delete(t)
            if ok:
                n += 1
        return n

    def list_archived(self, limit=200):
        """归档任务（磁盘 tasks/*.json，不在内存）：平台重启后的历史回溯。"""
        items, total = self._read_archived(limit=limit, logs=True)
        return items

    def count_archived(self):
        """磁盘归档任务总数（不受 list_archived 的 limit 截断影响）。"""
        try:
            return len([f for f in os.listdir(DIRS['tasks'])
                        if f.endswith('.json') and self._tid_ok(f[:-5])])
        except OSError:
            return 0

    def _read_archived(self, limit=200, logs=True):
        items = []
        try:
            files = [f for f in os.listdir(DIRS['tasks']) if f.endswith('.json')]
        except OSError:
            return items, 0
        total = 0
        for f in files:
            tid = f[:-5]
            if tid in self.tasks or not self._tid_ok(tid):
                continue
            total += 1
            if len(items) >= limit:
                continue
            try:
                with safe_open(os.path.join(DIRS['tasks'], f)) as fh:
                    d = json.load(fh)
                if not logs:
                    d.pop('log_tail', None)
                items.append(d)
            except (OSError, ValueError):
                continue
        items.sort(key=lambda d: d.get('started') or 0, reverse=True)
        return items[:limit], total

    def full_log(self, tid, lines=400):
        """任务全量日志（尾部 lines 行）：优先内存，内存被回收/归档时读
        磁盘 log_file（json 里记录了路径），最后回退 json 内的 log_tail。"""
        if not self._tid_ok(tid):
            return None
        rec = self.tasks.get(tid)
        if rec and rec['log']:
            return list(rec['log'])[-lines:]
        log_file = (rec or {}).get('log_file')
        js = None
        if not log_file:
            try:
                with safe_open(os.path.join(DIRS['tasks'], tid + '.json')) as f:
                    js = json.load(f)
                log_file = js.get('log_file')
            except (OSError, ValueError):
                js = None
        if log_file and os.path.isfile(log_file):
            # 防御深度：归档 json 的 log_file 必须在平台目录内
            try:
                check_path(log_file, must_exist=True, in_platform=True)
            except (OSError, ValueError):
                log_file = None
        if log_file and os.path.isfile(log_file):
            tail = self._tail_log_file(log_file, lines)
            if tail is not None:
                return tail
        if js is not None:
            return (js.get('log_tail') or [])[-lines:]
        try:
            with safe_open(os.path.join(DIRS['tasks'], tid + '.json')) as f:
                return (json.load(f).get('log_tail') or [])[-lines:]
        except (OSError, ValueError):
            return None

    _logcache = {}   # {(path, size, mtime_ns, lines): [行...]} 日志尾部缓存

    @classmethod
    def _tail_log_file(cls, path, lines):
        """从日志文件尾部读 lines 行。

        长任务日志可达 MB 级，而任务中心展开时会 2s 轮询一次，
        逐行读完整个文件会持续重复磁盘 IO。这里 seek 到尾部只读
        必要字节，并以 (size, mtime_ns, lines) 作键缓存：文件未变
        时直接复用，轮询场景下完全命中。
        """
        try:
            st = os.stat(path)
        except OSError:
            return None
        key = (path, st.st_size, st.st_mtime_ns, lines)
        hit = cls._logcache.get(key)
        if hit is not None:
            return list(hit)
        try:
            size = st.st_size
            chunk = min(size, max(65536, lines * 256))
            with safe_open(path, 'rb') as f:
                if size > chunk:
                    f.seek(size - chunk)
                data = f.read()
        except OSError:
            return None
        rows = data.decode('utf-8', errors='replace').splitlines()
        if size > chunk and rows:
            rows = rows[1:]           # 首行可能被字节截断，丢弃
        tail = rows[-lines:]
        if len(cls._logcache) > 32:
            cls._logcache.pop(next(iter(cls._logcache)))
        cls._logcache[key] = tuple(tail)
        return tail


tm = TaskManager()


def _sample_result_preview(sample):
    """样品任务的完成预览：报告链接 + 关键数字 + 产物文件。"""
    from vp.pipeline import pipeline_overview, _safe_sample_name
    s = _safe_sample_name(sample)
    try:
        sd = check_path(os.path.join(DIRS['results'], s), must_exist=True,
                        in_platform=True)
    except (ValueError, FileNotFoundError):
        return {'kind': 'sample', 'sample': s}
    ov = pipeline_overview(sd)
    stages = [{'stage': x['stage'], 'name': x['name'], 'status': x['status'],
               'summary': x.get('summary') or ''} for x in ov.get('stages', [])]
    rpt = os.path.join(sd, '07_report', 'report.html')
    return {'kind': 'sample', 'sample': s,
            'report': f'/report/{s}/' if os.path.isfile(rpt) else None,
            'stages': stages}


def _tool_result_files(run_name, limit=12):
    """工具运行目录的产物清单（按修改时间倒序，取前 limit 个）。"""
    root = check_path(os.path.join(_tool_runs_root(), run_name),
                      must_exist=False, in_platform=True)
    if not os.path.isdir(root):
        return []
    items = []
    for cur, _sub, fns in os.walk(root):
        for fn in fns:
            p = os.path.join(cur, fn)
            try:
                items.append((os.path.getmtime(p), os.path.relpath(p, root),
                              os.path.getsize(p)))
            except OSError:
                continue
    items.sort(reverse=True)
    return [{'path': rel.replace(os.sep, '/'), 'size': fmt_size(sz)}
            for _mt, rel, sz in items[:limit]]


def build_result_preview(result, task_name=''):
    """把任务返回值规范化成前端可渲染的预览 dict（失败安全）。"""
    try:
        if isinstance(result, str) and result:
            # 样品任务：返回 results/<样品> 目录
            rel = os.path.relpath(result, DIRS['results'])
            if rel.startswith('..') or os.path.isabs(rel):
                return {'kind': 'raw', 'text': str(result)}
            return _sample_result_preview(rel)
        if isinstance(result, dict):
            out = {'kind': 'tool',
                   'stats': {str(k): v for k, v in result.items()
                             if isinstance(v, (str, int, float, bool))
                             or v is None}}
            run = result.get('run')
            if run:
                out['run'] = str(run)          # 前端 toolrun 需顶层 run 名
                out['files'] = _tool_result_files(str(run))
            return out
    except Exception:
        pass
    return None

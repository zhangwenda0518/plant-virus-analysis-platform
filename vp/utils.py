# -*- coding: utf-8 -*-
"""
通用工具：路径安全校验、命令执行、日志、FASTA/FASTQ 流式处理、步骤断点标记。
所有 subprocess 调用一律使用参数列表（shell=False），不拼接命令字符串。
所有文件访问统一经由 safe_open()/safe_remove()：先 check_path 规范化校验
（拒绝 .. 段、写模式限定平台根内、读模式要求存在），再执行 I/O。
"""
import os
import re
import sys
import gzip
import time
import threading
import subprocess

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

from .config import PLATFORM_ROOT


def decode_output(raw):
    """子进程输出解码：优先 UTF-8，失败退 GBK，最后才 replace。

    Windows 上不少工具（bcftools/minimap2/samtools 等）按系统 ANSI 代码页
    (936/GBK) 打印命令行，其中含中文路径。若直接用 encoding='utf-8',
    errors='replace' 解码，每字节变 U+FFFD 且不可逆。先按字节读入再逐
    编码试解，可无损保留中文。
    """
    if isinstance(raw, str):
        return raw
    for enc in ('utf-8', 'gbk'):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode('utf-8', errors='replace')


# ------------------------------------------------------------------
# 路径安全
# ------------------------------------------------------------------
def check_path(path, must_exist=False, in_platform=False):
    """规范化并校验路径。

    - 拒绝路径中含 '..' 段（防穿越）
    - in_platform=True 时强制路径位于平台根目录内（用于输出/中间文件）
    - must_exist=True 时要求文件/目录已存在
    返回规范化后的绝对路径。
    """
    raw = str(path)
    for seg in raw.replace('/', os.sep).split(os.sep):
        if seg == '..':
            raise ValueError(f"非法路径（含 ..）: {raw}")
    p = os.path.normpath(os.path.abspath(raw))
    if in_platform:
        # 平台根 + 自定义输出根（若有）均视为合法写入区
        from .config import write_roots
        for root in write_roots():
            root = os.path.normpath(root)
            if p == root or p.startswith(root + os.sep):
                break
        else:
            raise ValueError(f"路径越界（仅允许平台或输出目录内）: {raw}")
    if must_exist and not os.path.exists(p):
        raise FileNotFoundError(f"路径不存在: {p}")
    return p


_STEP_RE = re.compile(r'^[A-Za-z0-9_\-]+$')


def check_step(step):
    """步骤名白名单校验：仅字母/数字/下划线/连字符，防止拼入文件名做路径注入。"""
    s = str(step)
    if not _STEP_RE.match(s):
        raise ValueError(f"非法步骤名: {s!r}")
    return s


# ---- 统一文件 I/O 入口（open/os.remove 仅出现在以下两个函数内） ----
def safe_open(path, mode='rt', compress_level=6):
    """安全打开文件。

    读模式（mode 含 r）：路径必须已存在；
    写模式（mode 含 w/a/x）：路径必须位于平台根目录内，父目录自动创建。
    自动识别 .gz（读与写均支持；写 gz 需显式以 .gz 结尾）。
    .gz 文件在 crabz.exe 可用时改走 crabz 多线程管道（比 Python gzip
    快数倍，写默认级别 6 而非 gzip 模块的 9），失败自动回退。
    compress_level 仅写 .gz 时生效（crabz -l；回退 gzip.compresslevel）。
    """
    is_read = 'r' in mode
    is_write = any(c in mode for c in ('w', 'a', 'x'))
    is_bin = 'b' in mode
    if is_read and not is_write:
        p = check_path(path, must_exist=True)
        if p.endswith('.gz'):
            cz = crabz_path()
            if cz:
                try:
                    return _GzReadPipe(p, mode)
                except OSError:
                    pass     # crabz 启动失败 → 回退 gzip
            if is_bin:
                return gzip.open(p, mode)
            return gzip.open(p, mode, encoding='utf-8', errors='replace',
                             newline=None)
        if is_bin:
            return open(p, mode)
        return open(p, mode, encoding='utf-8', errors='replace', newline=None)
    # 写模式：限定平台内（追加模式 gz 不支持 crabz 管道，回退 gzip）
    p = check_path(path, in_platform=True)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    if p.endswith('.gz'):
        cz = crabz_path()
        if cz and 'a' not in mode:
            try:
                return _GzWritePipe(p, mode, level=compress_level)
            except OSError:
                pass
        if is_bin:
            return gzip.open(p, mode, compresslevel=compress_level)
        # fastq/fasta 一律 LF 换行（CRLF 会导致 SPAdes 等解析失败）
        return gzip.open(p, mode, encoding='utf-8', errors='replace', newline='\n',
                         compresslevel=compress_level)
    if is_bin:
        return open(p, mode)
    return open(p, mode, encoding='utf-8', errors='replace', newline='\n')


def safe_remove(path, in_platform=True):
    """安全删除文件（默认限定平台内）。"""
    p = check_path(path, must_exist=True, in_platform=in_platform)
    os.remove(p)


# 兼容别名（历史调用点）
def open_maybe_gzip(path, mode='rt'):
    return safe_open(path, mode)


def open_write(path):
    return safe_open(path, 'wt')


def open_write_bin(path):
    return safe_open(path, 'wb')


# ------------------------------------------------------------------
# 日志
# ------------------------------------------------------------------
class TaskLogger:
    """线程安全日志器：写 UTF-8 文件 + 推送到回调（GUI 实时显示）。"""

    def __init__(self, log_file=None, callback=None, echo=False):
        self.log_file = log_file
        self.callback = callback          # callable(str) or None
        self.echo = echo
        self._lock = threading.Lock()
        self._fh = None
        if log_file:
            self._fh = safe_open(log_file, 'at')

    def log(self, msg, level="INFO"):
        line = f"[{time.strftime('%H:%M:%S')}] [{level}] {msg}"
        with self._lock:
            if self._fh:
                self._fh.write(line + "\n")
                self._fh.flush()
            if self.echo:
                try:
                    print(line)
                except Exception:
                    pass
        if self.callback:
            try:
                self.callback(line)
            except Exception:
                pass

    def close(self):
        with self._lock:
            if self._fh:
                self._fh.close()
                self._fh = None


# ------------------------------------------------------------------
# crabz 加速的 gzip 流（可选；未装 crabz 时回退 Python gzip）
# ------------------------------------------------------------------
_CRABZ_PATH = None
_CRABZ_CHECKED = False


def crabz_path():
    """crabz.exe 路径（未安装返回 None）。Rust 多线程 gzip，替代 Python gzip。"""
    global _CRABZ_PATH, _CRABZ_CHECKED
    if not _CRABZ_CHECKED:
        _CRABZ_CHECKED = True
        try:
            from .config import get_config
            try:
                _CRABZ_PATH = get_config().tool('crabz')
            except (FileNotFoundError, RuntimeError):
                _CRABZ_PATH = None
        except Exception:
            _CRABZ_PATH = None
    return _CRABZ_PATH


def _crabz_threads():
    try:
        from .config import get_config
        return max(2, min(4, get_config().threads))
    except Exception:
        return 4


def _popen_hidden(cmd):
    flags = 0
    if sys.platform == 'win32':
        flags = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            stdin=subprocess.PIPE, creationflags=flags)


class _GzReadPipe:
    """crabz -d 解压流。用法与 gzip.open(path,'rb'/'rt') 等价的上下文对象。

    safe_open 的调用方直接使用返回对象（不经 with），故代理内部文件
    对象的读写方法（readline/read/迭代/...），close() 附加进程回收。
    """

    def __init__(self, path, mode):
        proc = _popen_hidden([crabz_path(), '-d', '-Q', '-p', '2', str(path)])
        self._proc = proc
        if 't' in mode:
            import io
            self._fh = io.TextIOWrapper(
                proc.stdout, encoding='utf-8', errors='replace', newline=None)
        else:
            self._fh = proc.stdout
        self._closed = False

    def __enter__(self):
        return self._fh

    def __exit__(self, *exc):
        self.close()

    def __iter__(self):
        return iter(self._fh)

    def __getattr__(self, name):
        return getattr(self._fh, name)

    def close(self):
        if self._closed:
            return
        self._closed = True
        try:
            self._fh.close()
        except OSError:
            pass
        rc = self._proc.wait(timeout=60)
        if rc not in (0, None):
            raise RuntimeError(f"crabz 解压失败 (退出码 {rc})")


class _GzWritePipe:
    """crabz 压缩流（stdin 进，文件出）。等价 gzip.open(path,'wb'/'wt')。"""

    def __init__(self, path, mode, level=6):
        p = check_path(path, in_platform=True)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        proc = _popen_hidden([crabz_path(), '-Q', '-l', str(level),
                              '-p', str(_crabz_threads()), '-o', p])
        self._proc = proc
        self._path = p
        if 't' in mode:
            import io
            self._fh = io.TextIOWrapper(
                proc.stdin, encoding='utf-8', errors='replace', newline='\n')
        else:
            self._fh = proc.stdin
        self._closed = False

    def __enter__(self):
        return self._fh

    def __exit__(self, *exc):
        self.close()

    def __iter__(self):
        return iter(self._fh)

    def __getattr__(self, name):
        return getattr(self._fh, name)

    def close(self):
        if self._closed:
            return
        self._closed = True
        try:
            self._fh.close()
        except (OSError, ValueError):
            pass
        rc = self._proc.wait(timeout=600)
        if rc not in (0, None):
            try:
                os.remove(self._path)
            except OSError:
                pass
            raise RuntimeError(f"crabz 压缩失败 (退出码 {rc})")


def dir_size(path):
    """目录总字节数（元数据遍历，用于实测 chunk 磁盘占用）。"""
    total = 0
    for dp, _dn, fns in os.walk(path):
        for fn in fns:
            try:
                total += os.path.getsize(os.path.join(dp, fn))
            except OSError:
                pass
    return total


def est_decompressed(paths):
    """估算一组（可能 .gz）输入解压后的总字节数。
    FASTQ 压缩比经验值 ~5.5:1，FASTA ~4:1（无质量行），普通文件 1:1。"""
    total = 0
    for p in paths:
        try:
            sz = os.path.getsize(str(p))
        except OSError:
            continue
        low = str(p).lower()
        if low.endswith(('.fastq.gz', '.fq.gz')):
            total += sz * 5.5
        elif low.endswith('.gz'):
            total += sz * 4.0
        else:
            total += sz
    return total


def log_res_plan(logger, title, threads=None, mem_gb=None, disk_gb=None,
                 note=''):
    """阶段资源预估日志（统一格式，便于核对各阶段内存/线程/磁盘）。"""
    parts = [f"线程 {threads if threads is not None else '自动'}"]
    if mem_gb is not None:
        parts.append(f"内存 ~{mem_gb:.1f} GB")
    if disk_gb is not None:
        parts.append(f"磁盘 ~{disk_gb:.1f} GB")
    if note:
        parts.append(note)
    if logger:
        logger.log(f"📊 资源预估 [{title}] " + "；".join(parts), "PLAN")


def fmt_eta(seconds):
    """秒 → 人读时长（2小时3分 / 5分10秒 / 40秒）。"""
    if seconds is None or seconds < 0:
        return ''
    s = int(seconds)
    if s >= 3600:
        return f"{s // 3600}小时{(s % 3600) // 60:02d}分"
    if s >= 60:
        return f"{s // 60}分{s % 60:02d}秒"
    return f"{s}秒"


# ------------------------------------------------------------------
# 命令执行
# ------------------------------------------------------------------
def total_ram():
    """物理内存总量（字节）；失败返回 0。"""
    import ctypes
    import sys as _sys
    if not _sys.platform.startswith('win'):
        try:
            with open('/proc/meminfo') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        return int(line.split()[1]) * 1024
        except OSError:
            return 0
        return 0

    class MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [('dwLength', ctypes.c_ulong),
                    ('dwMemoryLoad', ctypes.c_ulong),
                    ('ullTotalPhys', ctypes.c_ulonglong),
                    ('ullAvailPhys', ctypes.c_ulonglong),
                    ('ullTotalPageFile', ctypes.c_ulonglong),
                    ('ullAvailPageFile', ctypes.c_ulonglong),
                    ('ullTotalVirtual', ctypes.c_ulonglong),
                    ('ullAvailVirtual', ctypes.c_ulonglong),
                    ('ullAvailExtendedVirtual', ctypes.c_ulonglong)]

    stat = MEMORYSTATUSEX()
    stat.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
    if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat)):
        return int(stat.ullTotalPhys)
    return 0


# ------------------------------------------------------------------
# 任务级取消上下文：TaskManager 在每个任务线程绑定 cancel 事件与子进程
# 注册表；run_cmd / run_cmd_redirect 自动注册所启动的进程并在取消时
# 由 TaskManager 杀进程树（硬停止）。非任务线程（CLI 等）绑定为 None，
# 行为与从前完全一致。
# ------------------------------------------------------------------
_TASK_STATE = threading.local()


def task_bind(cancel=None, procs=None):
    """绑定当前线程（任务）的取消事件与子进程注册表；不传即解绑。"""
    _TASK_STATE.cancel = cancel
    _TASK_STATE.procs = procs


def task_check_cancel():
    """任务取消事件置位时抛 RuntimeError（长循环内可随时调用）。"""
    ev = getattr(_TASK_STATE, 'cancel', None)
    if ev is not None and ev.is_set():
        raise RuntimeError('任务已停止（用户取消）')


def _task_register(proc):
    procs = getattr(_TASK_STATE, 'procs', None)
    if procs is not None:
        procs.append(proc)


def task_register_proc(proc):
    """公开入口：run_cmd 之外自行 Popen 的代码也把子进程登记进任务，
    供 TaskManager.cancel 硬停止。非任务线程内调用是空操作。"""
    _task_register(proc)


def run_cmd(cmd, logger=None, timeout=None, cwd=None, check=True,
            env=None, silent=False, monitor_fn=None, monitor_interval=5.0):
    """执行外部命令（参数列表形式）。stdout/stderr 实时写入日志。

    monitor_fn: 可选回调，命令运行期间每 monitor_interval 秒调用一次
    （用于长命令的进度看门狗，如 chunk 目录增长、spades 日志解析）。
    返回 returncode；check=True 时非零退出抛 RuntimeError。
    """
    cmd = [str(c) for c in cmd]
    cwd = check_path(cwd, must_exist=True) if cwd else None
    if logger:
        logger.log("$ " + subprocess.list2cmdline(cmd), "CMD")
    task_check_cancel()
    creationflags = 0
    if sys.platform == 'win32':
        creationflags = getattr(subprocess, 'CREATE_NO_WINDOW', 0)

    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=False,
        cwd=cwd, env=env, creationflags=creationflags)
    _task_register(proc)

    stop_mon = threading.Event()
    mon_thread = None
    if monitor_fn:
        def _watch():
            while not stop_mon.wait(monitor_interval):
                try:
                    monitor_fn()
                except Exception:
                    pass
        mon_thread = threading.Thread(target=_watch, daemon=True)
        mon_thread.start()

    # timeout 看门狗：阻塞读 stdout 时 wait(timeout) 永远轮不到，
    # 由看门狗线程在超时后强杀进程，读循环随即拿到 EOF。
    if timeout:
        def _kill_watchdog():
            if not stop_mon.wait(timeout):
                try:
                    proc.kill()
                except Exception:
                    pass
        _wd = threading.Thread(target=_kill_watchdog, daemon=True)
        _wd.start()

    t0 = time.time()
    try:
        for raw_line in proc.stdout:
            # 按字节读入后 lossless 解码：工具若按 GBK 打印中文路径可完整保留
            line = decode_output(raw_line).rstrip()
            if line and not silent:
                if logger:
                    logger.log(line, "OUT")
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
        raise RuntimeError(f"命令超时被终止: {cmd[0]}")
    except KeyboardInterrupt:
        proc.kill()
        proc.wait()
        raise
    finally:
        if mon_thread:
            stop_mon.set()
            mon_thread.join(timeout=2)

    rc = proc.returncode
    dt = time.time() - t0
    if logger:
        logger.log(f"[退出码 {rc}] 耗时 {dt:.1f}s", "OUT")
    task_check_cancel()
    if check and rc != 0:
        raise RuntimeError(f"命令失败 (退出码 {rc}): {subprocess.list2cmdline(cmd)}")
    return rc


def run_cmd_redirect(cmd, out_path, logger=None, timeout=None, cwd=None):
    """执行外部命令并将 stdout 重定向到平台内文件（mafft/fasttree 等输出流）。

    返回 returncode；非零抛 RuntimeError。
    """
    cmd = [str(c) for c in cmd]
    out = check_path(out_path, must_exist=False, in_platform=True)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cwd = check_path(cwd, must_exist=True) if cwd else None
    if logger:
        logger.log(f"$ {subprocess.list2cmdline(cmd)} > {os.path.basename(out)}", "CMD")
    task_check_cancel()
    creationflags = 0
    if sys.platform == 'win32':
        creationflags = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
    t0 = time.time()
    with safe_open(out, 'wb') as fo:
        proc = subprocess.Popen(
            cmd, stdout=fo, stderr=subprocess.PIPE,
            cwd=cwd, creationflags=creationflags)
        _task_register(proc)
        err_lines = []
        try:
            for line in iter(proc.stderr.readline, b''):
                txt = line.decode('utf-8', errors='replace').rstrip()
                if txt:
                    err_lines.append(txt)
                    if logger:
                        logger.log(txt, "OUT")
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
            raise RuntimeError(f"命令超时被终止: {cmd[0]}")
    task_check_cancel()
    rc = proc.returncode
    if logger:
        logger.log(f"[退出码 {rc}] 耗时 {time.time() - t0:.1f}s", "OUT")
    if rc != 0:
        raise RuntimeError(
            f"命令失败 (退出码 {rc}): {subprocess.list2cmdline(cmd)}\n"
            + "\n".join(err_lines[-10:]))
    return rc


# ------------------------------------------------------------------
# 步骤断点标记（work_dir 均位于平台 results/ 内）
# ------------------------------------------------------------------
def step_done_file(work_dir, step):
    step = check_step(step)
    wd = check_path(work_dir, in_platform=True, must_exist=True)
    return os.path.join(wd, f".{step}.done")


def is_step_done(work_dir, step):
    return os.path.isfile(step_done_file(work_dir, step))


def mark_step_done(work_dir, step):
    with safe_open(step_done_file(work_dir, step), 'wt') as f:
        f.write(time.strftime('%Y-%m-%d %H:%M:%S'))


def clear_step_done(work_dir, step=None):
    wd = check_path(work_dir, in_platform=True, must_exist=True)
    targets = []
    if step is None:
        for fn in os.listdir(wd):
            if fn.startswith('.') and fn.endswith('.done'):
                targets.append(check_path(os.path.join(wd, fn), in_platform=True))
    else:
        step = check_step(step)
        p = os.path.join(wd, f".{step}.done")
        if os.path.isfile(p):
            targets.append(check_path(p, in_platform=True))
    for p in targets:
        try:
            os.remove(p)
        except OSError:
            pass


# ------------------------------------------------------------------
# FASTA / FASTQ 流式处理
# ------------------------------------------------------------------
def iter_fasta(path):
    """迭代 FASTA：yield (header_no_arrow, seq)"""
    header, chunks = None, []
    with safe_open(path) as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, ''.join(chunks)


def write_fasta_record(fh, name, seq, width=60):
    fh.write(f">{name}\n")
    if width and len(seq) > width:
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")
    else:
        fh.write(seq + "\n")


def iter_fastq(path):
    """迭代 FASTQ：yield (name_no_at, seq, plus, qual)。"""
    with safe_open(path) as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            if not l4:
                raise ValueError(f"FASTQ 记录数不是 4 的倍数（文件截断?）: {path}")
            yield l1.rstrip('\n\r'), l2.rstrip('\n\r'), l3.rstrip('\n\r'), l4.rstrip('\n\r')


def iter_fastq_records(path):
    """迭代 FASTQ：yield 4 行原始文本（原样写出用）。

    path 可为文件路径或已打开的文本文件对象（调用方负责关闭）。
    """
    if hasattr(path, 'readline'):
        f = path
        outer = None
    else:
        f = outer = safe_open(path)
    try:
        while True:
            lines = [f.readline() for _ in range(4)]
            if not lines[0]:
                break
            if not lines[3]:
                raise ValueError(f"FASTQ 记录数不是 4 的倍数（文件截断?）: {path}")
            yield lines
    finally:
        if outer is not None:
            outer.close()


def count_fasta_seqs(path):
    n = 0
    with safe_open(path) as f:
        for line in f:
            if line.startswith('>'):
                n += 1
    return n


# ------------------------------------------------------------------
# FASTA 头 taxid 注入
# ------------------------------------------------------------------
def _clean_seq_id(header, fallback_no):
    """从 FASTA 头提取规范 seq_id：取首空白前 token、幂等去除历史
    kraken:taxid 标签；空 ID 自动编号 seq_1/seq_2/...（kunpeng 要求非空）。
    返回 (seq_id, 其余描述部分)。"""
    tokens = header.split(None, 1)
    seq_id = tokens[0] if tokens else ''
    desc = tokens[1] if len(tokens) > 1 else ''
    if '|kraken:taxid|' in seq_id:
        seq_id = seq_id.split('|kraken:taxid|')[0]
    if not seq_id:
        seq_id = f'seq_{fallback_no}'
    return seq_id, (' ' + desc if desc else '')


def inject_taxid_to_fasta(in_fasta, out_fasta, taxid, logger=None):
    """为 FASTA 每条序列头注入 taxid 标签（kunpeng add-library 要求）。

    kunpeng 解析要求标签位于第一个空白分隔 token 内: >ID|kraken:taxid|N
    """
    tag = f"|kraken:taxid|{taxid}"
    n = 0
    with safe_open(in_fasta) as fin, safe_open(out_fasta, 'wt') as fout:
        for line in fin:
            if line.startswith('>'):
                header = line[1:].rstrip('\n\r')
                n += 1
                if 'kraken:taxid|' in header:
                    fout.write('>' + header + '\n')   # 已有标签，保留
                    continue
                seq_id, rest = _clean_seq_id(header, n)
                fout.write(f'>{seq_id}{tag}{rest}\n'.rstrip() + '\n')
            else:
                fout.write(line if line.endswith('\n') else line + '\n')
    if logger:
        logger.log(f"taxid 注入完成: {n} 条序列 -> {os.path.basename(out_fasta)} (taxid={taxid})")
    return n


def inject_taxid_chunked(in_fasta, out_fasta, taxid, chunk_bp=1000000,
                         overlap_bp=34, logger=None):
    """注入 taxid 并把每条序列切成 ≤chunk_bp bp 的片段（kunpeng add-library 用）。

    kunpeng convert 阶段按 60 条/批整批载入序列且无字节上限：多条大染色体
    （如玉米 chr01 176MB）同批时内存随批内序列总长暴涨，宿主库 1.77GB 基因组
    实测触发 >20GB 分配失败。切成小片段后批内存有界，且 kmer 分类只依赖
    k-mer→taxid 映射，片段边界仅丢失每切口 k-1 个 k-mer，对去宿主无影响。

    overlap_bp: 相邻片段的重叠碱基数，默认 34 = k-1（kunpeng 默认 k=35），
    使跨切口 k-mer 全部保留、kmer 集合与不切割时完全一致（零损失）；0 关闭。

    片段标题: >{orig_id}.p{序号}|kraken:taxid|N {原描述}
    返回 (原始序列数, 输出片段数)。
    """
    tag = f"|kraken:taxid|{taxid}"
    n_seq = n_frag = 0
    seq_id = desc = ''
    part = 0
    acc = []          # 当前未写出序列的纯碱基块（无换行）
    acc_len = 0
    seen_ids = set()  # 记录已用 ID，重复 ID 自动加序号去重

    def flush():
        nonlocal n_frag, acc, acc_len
        if not acc_len:
            return
        # 尾部残余 ≤ overlap_bp 时内容已完全含于上一片段，无需重复写出
        if part > 0 and acc_len <= overlap_bp:
            acc, acc_len = [], 0
            return
        _write_piece(fout, seq_id, desc, tag, part, ''.join(acc))
        n_frag += 1
        acc = []
        acc_len = 0

    with safe_open(in_fasta) as fin, safe_open(out_fasta, 'wt') as fout:
        for line in fin:
            if line.startswith('>'):
                flush()
                header = line[1:].rstrip('\n\r')
                n_seq += 1
                seq_id, desc = _clean_seq_id(header, n_seq)
                if seq_id in seen_ids:          # 重复 ID 去重（如 draft 同名 scaffold）
                    seq_id = f'{seq_id}_r{n_seq}'
                seen_ids.add(seq_id)
                part = 0
                continue
            s = line.strip()
            if not s:
                continue
            acc.append(s)
            acc_len += len(s)
            if acc_len >= chunk_bp:
                # 依次切出整片段；下一片段从切点前 overlap_bp 处开始，
                # 跨切口的 k-mer 因此在相邻两片中各存在一份
                seq = ''.join(acc)
                while len(seq) >= chunk_bp:
                    _write_piece(fout, seq_id, desc, tag, part,
                                 seq[:chunk_bp])
                    n_frag += 1
                    part += 1
                    keep = chunk_bp - overlap_bp if overlap_bp else chunk_bp
                    seq = seq[keep:]
                acc = [seq]
                acc_len = len(seq)
        flush()
    if logger:
        logger.log(f"taxid 注入+分块完成: {n_seq} 条序列 -> {n_frag} 个片段"
                   f"（≤{chunk_bp // 1000}kb/片段, 重叠{overlap_bp}bp, "
                   f"taxid={taxid}）")
    return n_seq, n_frag


def _write_piece(fout, seq_id, desc, tag, part, seq):
    head = f"{seq_id}.p{part:05d}{tag}"
    fout.write(f">{head}{desc}\n".rstrip() + "\n")
    for i in range(0, len(seq), 60):
        fout.write(seq[i:i + 60] + "\n")


def inject_taxid_map(in_fasta, out_fasta, acc2taxid, logger=None, missing_fatal=True):
    """按 accession→taxid 映射逐条注入 kraken:taxid|N。

    acc2taxid: {accession: taxid(int)}；FASTA 头第一个空格前的 ID 用于查询。
    返回 (注入数, 未匹配列表)。
    """
    n, missing = 0, []
    with safe_open(in_fasta) as fin, safe_open(out_fasta, 'wt') as fout:
        for line in fin:
            if line.startswith('>'):
                header = line[1:].rstrip('\n\r')
                seq_id = header.split()[0] if header.split() else ''
                rest = header[len(seq_id):]
                taxid = acc2taxid.get(seq_id)
                if taxid is None:
                    missing.append(seq_id)
                    if missing_fatal:
                        raise KeyError(
                            f"序列 {seq_id} 在 info 表中未找到 Taxid（全部缺失项已记录）")
                    fout.write('>' + header + '\n')
                else:
                    fout.write(f'>{seq_id}|kraken:taxid|{taxid}{rest}\n'.rstrip() + '\n')
                n += 1
            else:
                fout.write(line if line.endswith('\n') else line + '\n')
    if logger:
        logger.log(f"taxid 映射注入: {n} 条, 未匹配 {len(missing)} 条")
        if missing:
            logger.log("未匹配示例: " + ", ".join(missing[:10]), "WARN")
    if missing and missing_fatal:
        raise KeyError(f"共 {len(missing)} 条序列缺少 Taxid 映射")
    return n, missing


# ------------------------------------------------------------------
# 杂项
# ------------------------------------------------------------------
def fmt_size(nbytes):
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if nbytes < 1024:
            return f"{nbytes:.1f}{unit}"
        nbytes /= 1024
    return f"{nbytes:.1f}PB"


def sample_name_from_fastq(path):
    """NX-5_S2_L001_R1_001.fastq.gz -> NX-5"""
    base = os.path.basename(str(path))
    for suffix in ('.fastq.gz', '.fq.gz', '.fastq', '.fq', '.fasta.gz', '.fa.gz', '.fasta', '.fa', '.fna'):
        if base.lower().endswith(suffix):
            base = base[:-len(suffix)]
            break
    # 去 Illumina 命名尾部
    for pat in ('_R1_001', '_R2_001', '_R1', '_R2', '_1', '_2'):
        if base.endswith(pat):
            base = base[:-len(pat)]
            break
    base = base.split('_L00')[0]
    return base or 'sample'

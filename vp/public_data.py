# -*- coding: utf-8 -*-
"""
公共数据下载：SRA/ENA (SRR/ERR/DRR) 与 GSA/NGDC (CRR) 的链接解析与
批量下载管理（aria2c 优先、内置 HTTP 回退）。

安全约束（强制）：
- 仅允许 http/https；主机白名单（ENA 镜像 / NGDC 域）；所有出站请求前
  逐次校验：拒绝环回/私有/保留/链路本地地址，重定向每一跳重新校验
  （防 SSRF / DNS rebinding）。
- 批次目录与输出文件在创建时统一经 vp.utils.check_path 校验
  （拒绝 .. 段、限定平台 downloads/ 内），文件名严格白名单清洗。
- 完整性：ENA 提供的字节数比对 + gzip 魔数抽样。
- aria2c 低并发分段（-x4 -s4），NGDC 场景防 IP 封控（经验吸收自
  MMPV-RNA public_metadata_pipeline/gsa_sra.down.py）。
"""
import os
import re
import json
import time
import gzip
import shutil
import socket
import subprocess
import threading
import urllib.request
import urllib.error

from .config import DIRS, PLATFORM_ROOT
from .utils import check_path, safe_open

# ------------------------------------------------------------------
# 常量与安全
# ------------------------------------------------------------------
ALLOWED_HOSTS = {
    'ftp.sra.ebi.ac.uk',            # ENA 读文件 HTTPS 镜像
    'www.ebi.ac.uk',                # ENA API（解析用）
    'ngdc.cncb.ac.cn',              # NGDC 搜索页（解析用）
    'download.cncb.ac.cn',          # GSA 下载
    'download.big.ac.cn',           # GSA 下载（big 镜像）
    'ddbj.nig.ac.jp',               # DDBJ（DRR 编号的原始 archive）
    'eutils.ncbi.nlm.nih.gov',      # NCBI eutils（DRR→DRA/DRX 元数据解析）
}
_IP_RE = re.compile(r'^\d{1,3}(\.\d{1,3}){3}$')
_UA = ('Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 '
       '(KHTML, like Gecko) Chrome/124.0 Safari/537.36')

DATA_FILE_RE = re.compile(
    r'\.(fastq\.gz|fq\.gz|fastq|fq|sra|fasta\.gz|fa\.gz|fasta|fa|fastq\.bz2)$', re.I)


def _is_private_ip(ip):
    parts = ip.split('.')
    if len(parts) != 4:
        return True
    try:
        a, b = int(parts[0]), int(parts[1])
    except ValueError:
        return True
    if a in (0, 10, 127) or a >= 224:            # 0.x / 10.x / loopback / 组播保留
        return True
    if a == 169 and b == 254:                     # link-local
        return True
    if a == 172 and 16 <= b <= 31:                # 172.16/12
        return True
    if a == 192 and b == 168:                     # 192.168/16
        return True
    return False


def check_url_host(url):
    """出站 URL 安全校验：仅 http(s)、主机在白名单、且不得解析到
    私有/保留地址。返回 (ok, msg)。所有网络请求前必须调用。"""
    m = re.match(r'^(https?)://([^/:?#]+)(?::(\d+))?([/?#].*)?$',
                 str(url or '').strip(), re.I)
    if not m:
        return False, f'仅支持 http/https URL: {url!r}'
    scheme, host = m.group(1).lower(), m.group(2).lower().rstrip('.')
    if _IP_RE.match(host) and _is_private_ip(host):
        return False, f'拒绝私有/保留 IP: {host}'
    if host == 'localhost' or host.endswith(('.local', '.internal')):
        return False, f'拒绝本机/内部主机名: {host}'
    if host not in ALLOWED_HOSTS:
        return False, (f'主机不在白名单内: {host}'
                       f'（允许: {", ".join(sorted(ALLOWED_HOSTS))}）')
    try:
        infos = socket.getaddrinfo(host, 443 if scheme == 'https' else 80,
                                   proto=socket.IPPROTO_TCP)
    except OSError as e:
        return False, f'主机解析失败: {host} ({e})'
    for info in infos:
        if _is_private_ip(info[4][0]):
            return False, f'{host} 解析到私有/保留地址: {info[4][0]}'
    return True, ''


class _SafeRedirectHandler(urllib.request.HTTPRedirectHandler):
    """重定向每一跳重新过白名单与私网校验。"""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        ok, why = check_url_host(newurl)
        if not ok:
            raise urllib.error.URLError(f'重定向目标被拒绝: {why}')
        return super().redirect_request(req, fp, code, msg, headers, newurl)


_OPENER = urllib.request.build_opener(_SafeRedirectHandler)


def _http_get(url, timeout=40):
    """受限 GET：请求前逐次校验（白名单 + 解析 IP + 重定向逐跳校验）。"""
    ok, why = check_url_host(url)
    if not ok:
        raise ValueError(f'出站请求被拒绝: {why}')
    req = urllib.request.Request(url,
                                 headers={'User-Agent': _UA, 'Accept': '*/*'})
    with _OPENER.open(req, timeout=timeout) as r:
        return r.read().decode('utf-8', 'replace')


def aria2c_path():
    """内置 aria2c.exe（bin/，兼容旧版根目录布局）；找不到再试 PATH。"""
    for p in (os.path.join(PLATFORM_ROOT, 'bin', 'aria2c.exe'),
              os.path.join(PLATFORM_ROOT, 'aria2c.exe')):
        if os.path.isfile(p):
            return p
    return shutil.which('aria2c') or shutil.which('aria2c.exe')


def sracha_path():
    """内置 sracha.exe（bin/，官方 0.7.0；兼容旧版根目录布局）；
    找不到再试 PATH。"""
    for p in (os.path.join(PLATFORM_ROOT, 'bin', 'sracha.exe'),
              os.path.join(PLATFORM_ROOT, 'sracha.exe')):
        if os.path.isfile(p):
            return p
    return shutil.which('sracha') or shutil.which('sracha.exe')


def fasterq_path():
    for p in (os.path.join(PLATFORM_ROOT, 'bin', 'fasterq-dump.exe'),
              os.path.join(PLATFORM_ROOT, 'fasterq-dump.exe')):
        if os.path.isfile(p):
            return p
    return shutil.which('fasterq-dump') or shutil.which('fasterq-dump.exe')


def sra_convert_engine():
    """选定 .sra→FASTQ 转换引擎：优先 sracha（纯 Rust 单文件，比
    fasterq-dump 快 5-13 倍，直接输出 gzip），回退 fasterq-dump。
    返回 (引擎名, 可执行路径) 或 (None, None)。"""
    p = sracha_path()
    if p and os.path.isfile(p):
        return 'sracha', p
    w = shutil.which('sracha') or shutil.which('sracha.exe')
    if w:
        return 'sracha', w
    p = fasterq_path()
    if p:
        return 'fasterq-dump', p
    return None, None


def _safe_filename(url):
    """从 URL 取文件名并强制白名单清洗（杜绝路径段与保留名）。"""
    fn = os.path.basename(str(url).split('?')[0].split('#')[0])
    fn = re.sub(r'[^A-Za-z0-9._\-]+', '_', fn).strip('._')
    if not fn or fn in {'.', '..'} or not DATA_FILE_RE.search(fn):
        raise ValueError(f'URL 文件名不合规（需数据文件后缀）: {url!r}')
    return fn[:120]


# ------------------------------------------------------------------
# Accession → 文件 URL 解析
# ------------------------------------------------------------------
_ENA_API = ('https://www.ebi.ac.uk/ena/portal/api/filereport'
            '?accession={acc}&result=read_run'
            '&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes,'
            'scientific_name,sample_accession&format=tsv')


def _dedupe_urls(urls):
    out, seen = [], set()
    for u in urls:
        u = u.strip()
        if u.startswith('//'):
            u = 'https:' + u
        u = re.sub(r'(?<=//SRR)', '/SRR', u)       # NGDC 页面常见的 //SRR 双斜杠
        u = re.sub(r'(?<=//ERR)', '/ERR', u)
        u = re.sub(r'(?<=//CRR)', '/CRR', u)
        if u not in seen:
            seen.add(u)
            out.append(u)
    return out


def _resolve_ena(acc):
    """SRR/ERR/DRR → ENA filereport → [{url, bytes}]。acc 已过严格正则。"""
    txt = _http_get(_ENA_API.format(acc=acc))
    lines = [l for l in txt.splitlines() if l.strip()]
    if len(lines) < 2:
        raise ValueError(f'ENA 未返回 {acc} 的数据行（run 可能不存在或未公开）')
    header = lines[0].split('\t')
    row = lines[1].split('\t')
    rec = dict(zip(header, row))
    ftps = [x for x in (rec.get('fastq_ftp') or '').split(';') if x]
    byts = [x for x in (rec.get('fastq_bytes') or '').split(';') if x]
    if not ftps:
        # ENA 未生成 FASTQ.GZ（常见于刚公布 / 发布者仅提供 SRA 原始数据）。
        # 只要 NCBI SRA 有数据（INSDC 三库元数据互通），就用 sracha
        # 从 NCBI .sra 拉取（--prefer-ena 会回落 NCBI）。
        # 不再直接判“受控数据”判死，避免“原来能下”的 run 因 ENA 无 fastq
        # 而整批失败。
        if re.fullmatch(r'[SED]RR\d+', acc) and sracha_path():
            return {'db': 'ENA→SRA', 'organism': rec.get('scientific_name', ''),
                    'bio_sample': rec.get('sample_accession', ''),
                    'files': [{'url': '', 'sracha': True}]}
        raise ValueError(f'ENA 无 {acc} 的 fastq 文件（可能为受控数据）')
    files = []
    for i, f in enumerate(ftps):
        url = 'https://' + f
        ok, msg = check_url_host(url)
        if not ok:
            raise ValueError(f'ENA 返回了不允许的地址: {msg}')
        try:
            size = int(byts[i]) if i < len(byts) else 0
        except ValueError:
            size = 0
        files.append({'url': url, 'bytes': size})
    # 单端库有时同时列出 .fastq.gz 与 _1/_2 —— _1/_2 齐全时只保留配对
    paired = [x for x in files if re.search(r'_1\.', x['url'])
              or re.search(r'_2\.', x['url'])]
    if len(paired) >= 2:
        files = paired
    return {'db': 'ENA', 'organism': rec.get('scientific_name', ''),
            'bio_sample': rec.get('sample_accession', ''), 'files': files}


def _resolve_ddbj(acc):
    """DRR → DDBJ 公共 FASTQ（.fastq.bz2）。

    DDBJ 与 ENA/NCBI 只共享元数据，原始文件留在提交方 archive，
    因此 DRR 编号必须走 DDBJ。布局：
      dra/fastq/<DRA前6>/<DRA>/<DRX>/<acc>.fastq.bz2        （单端）
      dra/fastq/<DRA前6>/<DRA>/<DRX>/<acc>_1/_2.fastq.bz2   （双端）
    DRA/DRX 号经 NCBI eutils 元数据解析（INSDC 三库元数据互通）。
    """
    import urllib.parse as _up
    term = _up.quote(acc)
    txt = _http_get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
                    f'esearch.fcgi?db=sra&term={term}')
    m = re.search(r'<Id>(\d+)</Id>', txt)
    if not m:
        raise ValueError(f'NCBI 元数据中未找到 {acc}（可能未公开）')
    xml = _http_get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
                    f'efetch.fcgi?db=sra&id={m.group(1)}&rettype=xml')
    dra = re.search(r'(DRA\d{6})', xml)
    drx = re.search(r'(DRX\d{6})', xml)
    if not (dra and drx):
        raise ValueError(f'NCBI 元数据缺少 {acc} 的 DRA/DRX 号')
    base = ('https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/'
            f'{dra.group(1)[:6]}/{dra.group(1)}/{drx.group(1)}/')
    try:
        idx = _http_get(base)
    except Exception as e:
        raise ValueError(f'DDBJ 目录不可达: {e}') from e
    names = sorted(set(re.findall(
        r'href="(' + re.escape(acc) + r'[^"/]*\.fastq\.bz2)"', idx)))
    if not names:
        raise ValueError(f'DDBJ 暂无 {acc} 的 FASTQ（新提交的 bz2 生成有延迟，'
                         f'可稍后重试）')
    return {'db': 'DDBJ',
            'organism': '',
            'files': [{'url': base + n, 'bytes': 0, 'md5': ''}
                      for n in names]}


def _resolve_ngdc(acc):
    """CRR → NGDC 搜索页 → browse 详情页 → 下载 URL（吸收 gsa_sra.down.py 流程）。"""
    html = _http_get(f'https://ngdc.cncb.ac.cn/gsa/search?searchTerm={acc}')
    urls = [u for u in _dedupe_urls(re.findall(
        r'(?:(?:ftp|https?):)?//download[0-9]*\.(?:cncb|big)\.ac\.cn/[^"\'<>\s]+',
        html)) if acc in u]
    m = re.search(rf'href="((?:/gsa/)?browse/[^"]*{acc}[^"]*)"', html)
    if not m:
        html2 = _http_get(f'https://ngdc.cncb.ac.cn/gsa/search?db=GSA&term={acc}')
        m2 = re.search(rf'href="((?:/gsa/)?browse/[^"]*{acc}[^"]*)"', html2)
        path = m2.group(1) if m2 else None
    else:
        path = m.group(1)
    if path and not urls:
        detail = ('https://ngdc.cncb.ac.cn' + path
                  if path.startswith('/gsa/')
                  else 'https://ngdc.cncb.ac.cn/gsa/' + path.lstrip('/'))
        html3 = _http_get(detail)
        urls = [u for u in _dedupe_urls(re.findall(
            r'(?:(?:ftp|https?):)?//download[0-9]*\.(?:cncb|big)\.ac\.cn/[^"\'<>\s]+',
            html3)) if acc in u]
    files = []
    for u in urls:
        u = re.sub(r'^ftp://', 'https://', u)      # 统一 HTTPS（协议白名单）
        if not DATA_FILE_RE.search(u):
            continue
        ok, msg = check_url_host(u)
        if not ok:
            raise ValueError(f'NGDC 返回了不允许的地址: {msg}')
        if all(f['url'] != u for f in files):
            files.append({'url': u, 'bytes': 0})
    # cncb/big 双镜像去重：同文件名优先 cncb 主站
    by_name = {}
    for f in files:
        bn = os.path.basename(f['url'])
        cur = by_name.get(bn)
        if cur is None or ('cncb.ac.cn' in f['url'] and 'cncb.ac.cn' not in cur['url']):
            by_name[bn] = f
    files = list(by_name.values())
    if not files:
        raise ValueError(f'NGDC 未找到 {acc} 的公开下载文件（可能受控/未公开）')
    return {'db': 'GSA', 'organism': '', 'bio_sample': '', 'files': files}


_ACC_RE_ENA = re.compile(r'[SED]RR\d+')
_ACC_RE_CRR = re.compile(r'CRR\d+')


def resolve_accession(acc):
    """单条 accession → {db, organism, files:[{url,bytes}]}；
    不支持的编号抛 ValueError（CRA 项目请先取得 CRR 运行号列表）。"""
    acc = str(acc).strip().upper()
    if _ACC_RE_ENA.fullmatch(acc):
        if acc.startswith('DRR'):
            # DRR 归 DDBJ 管理，ENA 通常无数据；先走 DDBJ，失败回退 ENA
            try:
                return _resolve_ddbj(acc)
            except Exception:
                return _resolve_ena(acc)
        return _resolve_ena(acc)
    if _ACC_RE_CRR.fullmatch(acc):
        return _resolve_ngdc(acc)
    raise ValueError(f'不支持的编号格式: {acc}（支持 SRR/ERR/DRR/CRR）')


# ------------------------------------------------------------------
# 下载管理器（写入统一走 safe_open：内部再做一次路径与平台边界校验）
# ------------------------------------------------------------------
def _batch_root():
    return DIRS.get('downloads') or os.path.join(PLATFORM_ROOT, 'downloads')


def _sanitize(name):
    s = re.sub(r'[^A-Za-z0-9_\-.]+', '_', str(name)).strip('._-')
    return s[:40] or 'batch'


def _gzip_ok(path):
    """gzip 魔数 + 抽样解压：快速完整性检查。"""
    try:
        with open(path, 'rb') as f:
            if f.read(2) != b'\x1f\x8b':
                return False
        with gzip.open(path, 'rb') as f:
            f.read(1 << 16)
        return True
    except (OSError, EOFError):
        return False


class DownloadManager:
    """批量下载：文件队列 + N 并发 aria2c（内置 HTTP 回退）+ 进度状态。

    batch = {id, name, status, created, dir(已校验), n_files, done, failed,
             overall, files: [{acc, url, out, size, done_bytes, progress,
                               status, md5, verified, error, db, organism}]}
    状态: resolving → downloading → (converting) → completed / failed /
          cancelled / interrupted
    """

    def __init__(self):
        self.lock = threading.RLock()
        self.batches = {}
        self.procs = {}
        self._spd = {}          # out -> [last_bytes, last_ts, ema_bytes_per_s]
        self._load_all()

    # ---- 下载速率（字节/秒，EMA 平滑）----
    # 设计要点：不依赖 aria2c 的输出格式。任何一路下载只要周期性更新
    # done_bytes，就能算出速率 —— 内置 HTTP 回退路径同样生效。
    _SPD_ALPHA = 0.3

    @staticmethod
    def _aria_bytes(s):
        """aria2c 容量串（'12MiB' / '1.2GiB' / '0B'）→ 字节数。"""
        m = re.match(r'^\s*(\d+(?:\.\d+)?)\s*([KMGT]?)i?B\s*$', s or '', re.I)
        if not m:
            return 0
        mult = {'': 1, 'K': 1024, 'M': 1024 ** 2,
                'G': 1024 ** 3, 'T': 1024 ** 4}
        return int(float(m.group(1)) * mult[m.group(2).upper()])

    def _tick_speed(self, out, done):
        """记录一次采样，返回平滑后的速率（B/s）。

        采样间隔 <0.5s 时不更新（避免抖动），但也不清零 —— 保留上一次
        平滑值，前端显示不会闪成 0。
        """
        now = time.time()
        prev = self._spd.get(out)
        if prev is None:
            self._spd[out] = [done, now, 0.0]
            return 0.0
        last_bytes, last_ts, ema = prev
        dt = now - last_ts
        if dt >= 0.5:
            if done >= last_bytes:
                inst = (done - last_bytes) / dt
                ema = inst if ema <= 0 else \
                    (self._SPD_ALPHA * inst + (1 - self._SPD_ALPHA) * ema)
            prev[0], prev[1], prev[2] = done, now, ema
        return prev[2]

    def _prime_speed(self, out, done):
        """同步 EMA 基准（采用外部速率时调用，避免回退时瞬间跳变）。"""
        prev = self._spd.get(out)
        now = time.time()
        if prev is None:
            self._spd[out] = [done, now, 0.0]
        else:
            prev[0], prev[1] = done, now

    def _drop_speed(self, out):
        self._spd.pop(out, None)

    # ---- 持久化 ----
    def _save(self, b):
        with self.lock:
            try:
                meta = check_path(b['dir'] + os.sep + 'batch.json',
                                  must_exist=False, in_platform=True)
                os.makedirs(b['dir'], exist_ok=True)
                with safe_open(meta, 'wt') as f:
                    json.dump(b, f, ensure_ascii=False, indent=1)
            except (OSError, ValueError):
                pass

    def _blog(self, b, msg):
        """追加一条批次过程日志到 <批次目录>/batch.log。

        任何写入失败都静默吞掉，绝不影响下载流程本身。"""
        try:
            d = b.get('dir') or ''
            if not d:
                return
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, 'batch.log'), 'a',
                      encoding='utf-8') as f:
                f.write(f'[{time.strftime("%H:%M:%S")}] {msg}\n')
        except OSError:
            pass

    def _load_all(self):
        root = _batch_root()
        try:
            os.makedirs(root, exist_ok=True)
        except OSError:
            return
        for name in sorted(os.listdir(root)):
            try:
                meta = check_path(os.path.join(root, name, 'batch.json'),
                                  must_exist=False, in_platform=True)
            except ValueError:
                continue
            if not os.path.isfile(meta):
                continue
            try:
                with open(meta, encoding='utf-8') as f:
                    b = json.load(f)
                if b.get('status') in ('downloading', 'resolving', 'converting'):
                    b['status'] = 'interrupted'
                    for fe in b.get('files', []):
                        if fe.get('status') in ('downloading', 'pending',
                                                'verifying'):
                            fe['status'] = 'pending'
                self.batches[b['id']] = b
            except (OSError, ValueError):
                continue

    # ---- 创建 ----
    def create(self, name, accessions, concurrency=2, convert_sra=True):
        bid = f'{_sanitize(name)}_{time.strftime("%Y%m%d_%H%M%S")}'
        bdir = check_path(os.path.join(_batch_root(), bid),
                          must_exist=False, in_platform=True)
        b = {'id': bid, 'name': _sanitize(name), 'status': 'resolving',
             'created': time.strftime('%Y-%m-%d %H:%M'),
             'concurrency': max(1, min(int(concurrency or 2), 4)),
             'convert_sra': bool(convert_sra), 'dir': bdir,
             'files': [], 'runs': {}, 'error': '',
             'converting': False, 'convert_msg': ''}
        with self.lock:
            self.batches[bid] = b
        self._save(b)
        self._blog(b, f'创建批次 {b["name"]}：共 {len(accessions)} 个条目，'
                      f'输出目录 {os.path.basename(bdir)}')

        def resolver():
            ok_n = 0
            for raw in accessions:
                raw = str(raw).strip()
                if not raw or raw.startswith('#'):
                    continue
                try:
                    if re.match(r'^https?://', raw, re.I):
                        # 完整数据文件 URL：主机白名单 + 文件名白名单
                        ok, why = check_url_host(raw)
                        if not ok:
                            raise ValueError(why)
                        fn = _safe_filename(raw)
                        out = check_path(bdir + os.sep + fn,
                                         must_exist=False, in_platform=True)
                        if not out.startswith(bdir + os.sep):
                            raise ValueError('输出文件名越界（已拦截）')
                        with self.lock:
                            acc = fn.split('.')[0] or raw[:24]
                            b['runs'][acc] = {'db': 'URL', 'organism': ''}
                            b['files'].append({
                                'acc': acc, 'url': raw, 'out': out,
                                'size': 0, 'done_bytes': 0, 'progress': 0.0,
                                'status': 'pending', 'md5': '',
                                'verified': False, 'error': '', 'db': 'URL',
                                'organism': ''})
                        self._blog(b, f'URL 条目: {raw[:80]}')
                        ok_n += 1
                        continue
                    acc = raw.upper()
                    info = resolve_accession(acc)
                    with self.lock:
                        b['runs'][acc] = {'db': info['db'],
                                          'organism': info.get('organism', '')}
                        for fe in info['files']:
                            # sracha 条目：无 ENA url，输出名用 acc.sra（sracha 落地）
                            if fe.get('sracha'):
                                fn = acc + '.sra'
                            else:
                                fn = _safe_filename(fe['url'])
                            out = check_path(bdir + os.sep + fn,
                                             must_exist=False,
                                             in_platform=True)
                            if not out.startswith(bdir + os.sep):
                                raise ValueError('输出文件名越界（已拦截）')
                            # size        = 展示/进度用总量，下载中会被
                            #               aria2c 实时值覆盖（可能按整数
                            #               MiB 缩写，仅供进度条，不可校验）
                            # expect_bytes = 校验用精确总量，只来自 ENA/
                            #               DDBJ 元数据，全程不被覆盖
                            b['files'].append({
                                'acc': acc, 'url': fe['url'], 'out': out,
                                'size': int(fe.get('bytes') or 0),
                                'expect_bytes': int(fe.get('bytes') or 0),
                                'done_bytes': 0, 'progress': 0.0,
                                'status': 'pending',
                                'md5': str(fe.get('md5', ''))[:64],
                                'verified': False, 'error': '',
                                'db': info['db'],
                                'organism': info.get('organism', '')})
                    self._blog(b, f'解析 {acc} -> {info["db"]} '
                                  f'{len(info["files"])} 个文件')
                    ok_n += 1
                except Exception as e:
                    with self.lock:
                        b['files'].append({
                            'acc': raw.upper()[:40], 'url': '', 'out': '',
                            'size': 0, 'done_bytes': 0, 'progress': 0.0,
                            'status': 'failed', 'md5': '', 'verified': False,
                            'error': f'解析失败: {e}', 'db': '?',
                            'organism': ''})
                    self._blog(b, f'解析失败 {raw}: {e}')
                self._save(b)
            with self.lock:
                # cancel 可能落在解析期间：仅在仍是 resolving 时才转入下载
                if b['status'] == 'resolving':
                    b['status'] = 'downloading' if ok_n else 'failed'
                    if not ok_n:
                        b['error'] = '没有可下载的文件（全部解析失败）'
            if ok_n:
                self._blog(b, f'解析完成 {ok_n} 个条目，开始下载')
            else:
                self._blog(b, '解析全部失败')
            self._save(b)
            if ok_n and b['status'] == 'downloading':
                self._start_workers(bid)

        threading.Thread(target=resolver, daemon=True,
                         name=f'dl-resolve-{bid}').start()
        return bid

    # ---- 下载执行 ----
    def _start_workers(self, bid):
        threading.Thread(target=self._worker_loop, args=(bid,),
                         daemon=True, name=f'dl-{bid}').start()

    def _worker_loop(self, bid):
        with self.lock:
            b = self.batches.get(bid)
            conc = max(1, min(int((b or {}).get('concurrency', 2)), 4))
        threads = [threading.Thread(target=self._worker, args=(bid,),
                                    daemon=True) for _ in range(conc)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self._finalize(bid)

    def _next_pending(self, b):
        with self.lock:
            if b['status'] != 'downloading':
                return None
            for fe in b['files']:
                if fe['status'] == 'pending':
                    fe['status'] = 'downloading'
                    return fe
        return None

    def _worker(self, bid):
        while True:
            with self.lock:
                b = self.batches.get(bid)
                if not b or b['status'] != 'downloading':
                    return
            fe = self._next_pending(b)
            if not fe:
                return
            self._blog(b, f'开始下载 {fe["acc"]} ({fe["url"]})')
            try:
                self._download_file(b, fe)
            except Exception as e:                  # 单文件失败不拖垮批次
                with self.lock:
                    fe['status'] = 'failed'
                    fe['error'] = str(e)
                self._blog(b, f'失败 {fe["acc"]}: {e}')

    def _download_file(self, b, fe):
        try:
            out = check_path(fe['out'], must_exist=False, in_platform=True)
        except ValueError as e:
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = f'路径校验失败: {e}'
            self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
            return
        # sracha 条目：无 ENA url，直接从 NCBI .sra 拉取（--prefer-ena 自动回落）
        if fe.get('sracha') or (not fe.get('url') and re.fullmatch(r'[SED]RR\d+', fe.get('acc') or '') and sracha_path()):
            os.makedirs(b['dir'], exist_ok=True)
            ok_dl = self._dl_sracha(b, fe)
            if ok_dl:
                return
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = fe.get('error') or 'sracha 拉取失败'
            self._blog(b, f'失败 {fe["acc"]}: {fe.get("error")}')
            return
        ok, why = check_url_host(fe['url'])
        if not ok:
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = f'URL 校验失败: {why}'
            self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
            return
        if not out.startswith(b['dir'] + os.sep):   # 写盘前包含性检查
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = '输出路径越界（已拦截）'
            self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
            return
        os.makedirs(b['dir'], exist_ok=True)
        a2 = aria2c_path()
        ok_dl = self._dl_aria2(a2, b, fe) if a2 else self._dl_python(fe)
        if not ok_dl:
            # 三级回退：INSDC run 编号（SRR/ERR/DRR）可用 sracha fetch
            # （--prefer-ena 优先 ENA FASTQ.GZ，失败自动回落 NCBI .sra，
            # 下载即断点续传 + MD5 校验；.sra 由收尾转换链自动转 FASTQ）。
            # GSA CRR / 任意 URL 不适用（NCBI 无此数据），维持内置 HTTP。
            if re.fullmatch(r'[SED]RR\d+', fe.get('acc') or '') \
                    and sracha_path():
                with self.lock:
                    fe['status'] = 'downloading'
                    fe['error'] = ''
                ok_dl = self._dl_sracha(b, fe)
        if not ok_dl:
            with self.lock:
                if fe['status'] != 'cancelled':
                    fe['status'] = 'failed'
                    fe['error'] = fe.get('error') or '下载失败'
            if fe['status'] == 'cancelled':
                self._blog(b, f'取消 {fe["acc"]}')
            else:
                self._blog(b, f'失败 {fe["acc"]}: {fe.get("error")}')
            return
        # 完整性：字节数比对（仅当元数据给出精确值时才做强校验）
        # ⚠️ 不能用下载中 aria2c 实时报的 size 做精确比对 —— aria2c 对
        #    ≥10MiB 的数按整数 MiB 缩写（36.84MiB 显示成 '36MiB'），
        #    误差可达 0.5MiB，用它卡字节数会把正常文件判成损坏。
        #    实测：DDBJ DRR003724.fastq.bz2 真实 38638145，aria2c 报
        #    37748736（=36.00MiB），差值 889409 → 误判 failed。
        #    expect_bytes 为空时退化为 ±1MiB 容差的粗校验。
        exp = int(fe.get('expect_bytes') or 0)
        try:
            actual = os.path.getsize(out)
        except OSError as e:
            with self.lock:
                fe['error'] = f'完整性检查出错: {e}'
            actual = 0
        if exp and actual != exp:
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = (f'字节数不一致（期望 {exp}, '
                               f'实际 {actual}）')
            self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
            return
        if not exp and fe.get('size') and actual:
            if abs(actual - int(fe['size'])) > 1048576:
                with self.lock:
                    fe['status'] = 'failed'
                    fe['error'] = (f'字节数偏差过大（约 {fe["size"]}, '
                                   f'实际 {actual}）')
                self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
                return
        if out.lower().endswith('.gz') and not _gzip_ok(out):
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = 'gzip 文件损坏（魔数/抽样解压失败）'
            self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
            return
        with self.lock:
            fe['status'] = 'done'
            fe['verified'] = True
            fe['progress'] = 1.0
        self._blog(b, f'完成 {fe["acc"]}')
        self._save(b)

    def _dl_aria2(self, a2, b, fe):
        out = fe['out']
        cmd = [a2, '-x', '4', '-s', '4', '-k', '1M', '-c',
               '--file-allocation=none', '--console-log-level=warn',
               '--summary-interval=2', '--retry-wait=3', '--max-tries=5',
               '-d', b['dir'], '-o', os.path.basename(out), fe['url']]
        try:
            proc = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding='utf-8', errors='replace')
        except OSError:
            return self._dl_python(fe)
        with self.lock:
            self.procs[fe['out']] = proc
        # aria2c 实测进度行（--summary-interval=2 每 2 秒一帧）：
        #   总量已知  [#289963 4.7MiB/0.9GiB(0%) CN:4 DL:1.5MiB ETA:9m53s]
        #   总量未知  [#289963 0B/0B CN:1 DL:0B]        ← 无百分比，size 保持 0
        # 解析要点：
        #   · DONE/TOTAL 直接取 aria2c 认定的总量。原先 `size=max(size,
        #     done_bytes)` 会让总量跟着已下载量一起涨，进度条永远贴着 100% 假动。
        #   · DL: 是 aria2c 自报速率，首帧即有值，比落盘差分快一拍；
        #     拿不到时（格式变更 / 内置 HTTP 回退）才用 _tick_speed 差分。
        pct_re = re.compile(r'\((\d+)%\)')
        sz_re = re.compile(r'\[#\w+\s+([\d.]+[KMGT]?i?B)/([\d.]+[KMGT]?i?B)')
        dl_re = re.compile(r'\bDL:\s*([\d.]+[KMGT]?i?B)')
        try:
            for line in proc.stdout:
                m = pct_re.search(line)
                ms = sz_re.search(line)
                md = dl_re.search(line)
                if not (m or ms or md):
                    continue
                with self.lock:
                    if ms:
                        total = self._aria_bytes(ms.group(2))
                        if total > 0:
                            fe['size'] = total
                    try:
                        fe['done_bytes'] = os.path.getsize(out)
                    except OSError:
                        pass
                    # 有精确字节数就自己算（比整数百分比更细），否则退回 (NN%)
                    if fe.get('size') and fe.get('done_bytes'):
                        fe['progress'] = min(
                            fe['done_bytes'] / fe['size'], 1.0)
                    elif m:
                        fe['progress'] = min(int(m.group(1)), 100) / 100.0
                    if md:
                        fe['speed'] = self._aria_bytes(md.group(1))
                        self._prime_speed(out, fe.get('done_bytes') or 0)
                    elif fe.get('done_bytes'):
                        sp = self._tick_speed(out, fe['done_bytes'])
                        if sp:
                            fe['speed'] = sp
        finally:
            rc = proc.wait()
            with self.lock:
                self.procs.pop(fe['out'], None)
                self._drop_speed(out)
                fe.pop('speed', None)
        if rc != 0:
            with self.lock:
                fe['error'] = f'aria2c 退出码 {rc}'
            return False
        if not os.path.isfile(out) or os.path.getsize(out) == 0:
            with self.lock:
                fe['error'] = 'aria2c 完成但文件缺失/为空'
            return False
        with self.lock:
            try:
                fe['done_bytes'] = os.path.getsize(out)
            except OSError:
                pass
        return True

    def _dl_sracha(self, b, fe):
        """第三级回退：sracha fetch（仅 SRR/ERR/DRR）。

        --prefer-ena 优先下载 ENA 预生成的 FASTQ.GZ（与主路径同源同内容），
        ENA 无 FASTQ 时自动回落 NCBI 的 .sra（后续收尾链自动转 FASTQ）。
        自带断点续传与 MD5 校验。"""
        exe = sracha_path()
        acc = fe['acc']
        cmd = [exe, 'fetch', acc, '-O', b['dir'], '-f', '-q', '--no-progress',
               '--prefer-ena']
        try:
            rc = subprocess.run(
                cmd, timeout=24 * 3600,
                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE).returncode
        except (OSError, subprocess.TimeoutExpired) as e:
            with self.lock:
                fe['error'] = f'sracha 回退失败: {e}'
            return False
        if rc != 0:
            with self.lock:
                fe['error'] = f'sracha 退出码 {rc}'
            return False
        # 识别产物：acc.sra（NCBI 路径）或 acc[_N].fastq.gz / acc.fastq.gz（ENA 路径）
        names = set(os.listdir(b['dir']))
        sra_p = b['dir'] + os.sep + acc + '.sra'
        fastqs = sorted(
            n for n in names
            if n.startswith(acc)
            and (n == acc + '.fastq.gz'
                 or re.search(r'_1\.f(ast)?q\.gz$', n, re.I)
                 or re.search(r'_2\.f(ast)?q\.gz$', n, re.I)))
        with self.lock:
            if os.path.isfile(sra_p):
                # 交给收尾转换链（条目 out 指向 .sra，_finalize 会触发转换）
                fe['out'] = sra_p
                fe['size'] = os.path.getsize(sra_p)
                fe['done_bytes'] = fe['size']
                fe['status'] = 'done'
                fe['verified'] = True                # sracha 已做 MD5 校验
                fe['progress'] = 1.0
                return True
            if fastqs:
                produced_paths = []
                for n in fastqs:
                    produced_paths.append(
                        check_path(b['dir'] + os.sep + n, must_exist=True,
                                   in_platform=True))
                new_entries = []
                for e in b['files']:
                    if e['out'] == fe['out']:
                        for pp in produced_paths:
                            new_entries.append({
                                'acc': acc, 'url': '', 'out': pp,
                                'size': os.path.getsize(pp),
                                'done_bytes': os.path.getsize(pp),
                                'progress': 1.0, 'status': 'done',
                                'md5': '', 'verified': True, 'error': '',
                                'db': fe.get('db', 'ENA'),
                                'organism': fe.get('organism', '')})
                    else:
                        new_entries.append(e)
                b['files'] = new_entries
                return True
        with self.lock:
            fe['error'] = 'sracha 完成但未找到产物文件'
        return False

    def _dl_python(self, fe):
        """内置 HTTP 下载回退（Range 断点续传；URL 已在外层校验）。

        先写 .aria_part 临时文件（无 .gz 后缀 → safe_open 原样二进制写，
        不会对已压缩数据二次 gzip），完成后改名到位。"""
        out, url = fe['out'], fe['url']
        part = check_path(out + '.aria_part', must_exist=False,
                          in_platform=True)
        pos = os.path.getsize(part) if os.path.isfile(part) else 0
        headers = {'User-Agent': _UA}
        if pos:
            headers['Range'] = f'bytes={pos}-'
        try:
            req = urllib.request.Request(url, headers=headers)
            with _OPENER.open(req, timeout=60) as r:
                if r.status == 416:                 # 本地已完整
                    shutil.move(part, out)
                    return True
                total = int(r.headers.get('content-length') or 0) + \
                    (pos if r.status == 206 else 0)
                if r.status != 206:
                    pos = 0
                with self.lock:
                    fe['size'] = total or fe.get('size') or 0
                with safe_open(part, 'wb') as f:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        f.write(chunk)
                        pos += len(chunk)
                        with self.lock:
                            fe['done_bytes'] = pos
                            fe['progress'] = (pos / total) if total else 0.0
                            sp = self._tick_speed(out, pos)
                            if sp:
                                fe['speed'] = sp
            shutil.move(part, out)
            return True
        except (urllib.error.URLError, OSError) as e:
            with self.lock:
                fe['error'] = f'HTTP 下载失败: {e}'
            return False
        finally:
            with self.lock:
                self._drop_speed(out)
                fe.pop('speed', None)

    # ---- 收尾：.sra 转换 ----
    def _finalize(self, bid):
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return
            done = sum(1 for f in b['files'] if f['status'] == 'done')
            failed = sum(1 for f in b['files'] if f['status'] == 'failed')
            b['done'], b['failed'] = done, failed
            sra_left = [f for f in b['files']
                        if f['status'] == 'done'
                        and f['out'].lower().endswith('.sra')]
            bz2_left = [f for f in b['files']
                        if f['status'] == 'done'
                        and f['out'].lower().endswith('.fastq.bz2')]
            engine, _path = sra_convert_engine()
            b['status'] = ('converting'
                           if ((sra_left and b.get('convert_sra') and engine)
                               or bz2_left)
                           else ('completed' if done else 'failed'))
        if b['status'] == 'converting':
            self._blog(b, f'批次转换中（{len(sra_left)} 个 .sra 待转换）')
        elif b['status'] == 'completed':
            self._blog(b, f'批次完成 {b.get("done", 0)}/'
                          f'{len(b.get("files", []))} 成功')
        self._save(b)
        if b['status'] == 'converting':
            self._convert_sras(b)

    def _convert_sras(self, b):
        engine, exe = sra_convert_engine()
        with self.lock:
            if not b or not engine:
                if b:
                    b['status'] = 'completed'
                    b['convert_msg'] = ('无 .sra 需转换' if engine is None and b
                                        else b.get('convert_msg', ''))
                    self._save(b)
                    self._blog(b, f'批次完成 {b.get("done", 0)}/'
                                  f'{len(b.get("files", []))} 成功')
                return
            b['converting'] = True
            sras = [f['out'] for f in b['files']
                    if f['status'] == 'done'
                    and f['out'].lower().endswith('.sra')]
            threads = max(2, min((os.cpu_count() or 4) // 2, 8))
        msgs = []
        for sra in sras:
            try:
                if engine == 'sracha':
                    # 单文件静态引擎：直接输出 gzip 的 _1/_2（split-3 命名
                    # 与 fasterq-dump/ENA 一致），无需再压缩
                    subprocess.run(
                        [exe, 'fastq', sra, '-O', b['dir'],
                         '-t', str(threads), '-f', '-q'],
                        check=True, timeout=8 * 3600,
                        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                    base = os.path.basename(sra)[:-4]
                    produced = sorted(
                        p for p in os.listdir(b['dir'])
                        if p.startswith(base + '_')
                        and p.endswith(('.fastq.gz', '.fq.gz')))
                    if not produced:
                        raise RuntimeError('sracha 无输出')
                    os.remove(sra)
                    msgs.append(f'✔ {os.path.basename(sra)} → '
                                f'{len(produced)} FASTQ.GZ（sracha）')
                else:
                    subprocess.run(
                        [exe, '--split-files', '-e', str(threads), '-O',
                         b['dir'], sra],
                        check=True, timeout=4 * 3600,
                        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                    produced = sorted(
                        p for p in os.listdir(b['dir'])
                        if p.startswith(os.path.basename(sra) + '_')
                        and p.endswith('.fastq'))
                    if not produced:
                        raise RuntimeError('fasterq-dump 无输出')
                    for p in produced:
                        src = check_path(b['dir'] + os.sep + p,
                                         must_exist=True, in_platform=True)
                        with open(src, 'rb') as fin, \
                                safe_open(src + '.gz', 'wb', compress_level=6) \
                                as fout:
                            shutil.copyfileobj(fin, fout, 1 << 20)
                        os.remove(src)
                    os.remove(sra)
                    msgs.append(f'✔ {os.path.basename(sra)} → '
                                f'{len(produced)} FASTQ.GZ（fasterq-dump）')
                # 把 .sra 条目替换为产物 FASTQ 条目（供 ready_files 配对进流程）
                produced_paths = [check_path(b['dir'] + os.sep + p,
                                             must_exist=True,
                                             in_platform=True)
                                  for p in produced]
                with self.lock:
                    new_entries = []
                    for fe in b['files']:
                        if fe['out'] == sra:
                            acc = fe['acc']
                            for pp in produced_paths:
                                new_entries.append({
                                    'acc': acc, 'url': '', 'out': pp,
                                    'size': os.path.getsize(pp),
                                    'done_bytes': os.path.getsize(pp),
                                    'progress': 1.0, 'status': 'done',
                                    'md5': '', 'verified': True,
                                    'error': '', 'db': fe.get('db', 'URL'),
                                    'organism': fe.get('organism', '')})
                        else:
                            new_entries.append(fe)
                    b['files'] = new_entries
            except Exception as e:
                msgs.append(f'✘ {os.path.basename(sra)}: {e}')
        # DDBJ .fastq.bz2 条目（与 .sra 转换互不影响）
        bz2_entries = [f for f in b['files']
                       if f['status'] == 'done'
                       and f['out'].lower().endswith('.fastq.bz2')]
        if bz2_entries:
            self._blog(b, f'批次转换中（{len(bz2_entries)} 个 .fastq.bz2 待转换）')
            self._convert_bz2(b, bz2_entries)
        with self.lock:
            b['converting'] = False
            done = sum(1 for f in b['files'] if f['status'] == 'done')
            failed = sum(1 for f in b['files'] if f['status'] == 'failed')
            b['done'], b['failed'] = done, failed
            b['convert_msg'] = '；'.join(msgs) or '无 .sra'
            b['status'] = 'completed' if done or not b['files'] else 'failed'
        self._blog(b, f'批次完成 {b.get("done", 0)}/'
                      f'{len(b.get("files", []))} 成功')
        self._save(b)

    def _convert_bz2(self, b, entries):
        """DDBJ .fastq.bz2 → .fastq.gz（流式转压缩，gzip level 1）。
        单端单文件 → <acc>.fastq.gz；双端 _1/_2 → 同名 _1/_2.fastq.gz。
        完成后条目替换为 FASTQ（供 ready_files 配对进流程）。"""
        import bz2 as _bz2
        import gzip as _gzip
        for fe in entries:
            src = fe['out']
            try:
                base = os.path.basename(src)[:-len('.fastq.bz2')]
                outs = []
                with _bz2.open(src, 'rb') as fin:
                    if base.endswith('_1') or base.endswith('_2'):
                        dst = src[:-len('.fastq.bz2')] + '.fastq.gz'
                        with safe_open(dst, 'wb', compress_level=1) as fout:
                            shutil.copyfileobj(fin, fout, 1 << 20)
                        outs.append(dst)
                    else:
                        # 单端：直接转 gz，命名与 ENA/sracha 产物一致
                        dst = os.path.join(os.path.dirname(src),
                                           base + '.fastq.gz')
                        with safe_open(dst, 'wb', compress_level=1) as fout:
                            shutil.copyfileobj(fin, fout, 1 << 20)
                        outs.append(dst)
                os.remove(src)
                with self.lock:
                    new_entries = []
                    for e in b['files']:
                        if e is fe:
                            for pp in outs:
                                new_entries.append({
                                    'acc': e['acc'], 'url': '', 'out': pp,
                                    'size': os.path.getsize(pp),
                                    'done_bytes': os.path.getsize(pp),
                                    'progress': 1.0, 'status': 'done',
                                    'md5': '', 'verified': True,
                                    'error': '', 'db': e.get('db', 'DDBJ'),
                                    'organism': e.get('organism', '')})
                        else:
                            new_entries.append(e)
                    b['files'] = new_entries
                self._blog(b, f'✔ {os.path.basename(src)} → FASTQ.GZ（DDBJ）')
            except Exception as e:
                with self.lock:
                    fe['status'] = 'failed'
                    fe['error'] = f'bz2 转换失败: {e}'
                self._blog(b, f'✘ {os.path.basename(src)}: {e}')

    # ---- 状态 / 控制 ----
    def snapshot(self, bid):
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return None
            files = [{k: f.get(k) for k in
                      ('acc', 'url', 'out', 'size', 'done_bytes', 'progress',
                       'status', 'md5', 'verified', 'error', 'db',
                       'organism', 'speed')}
                     for f in b['files']]
            total = len(files)
            done = sum(1 for f in files if f['status'] == 'done')
            failed = sum(1 for f in files if f['status'] == 'failed')
            overall = (sum((f.get('progress') or 0) for f in files) / total
                       if total else 0.0)
            # 批次速率 = 正在下载条目的速率之和；ETA 由剩余字节 / 速率推算
            active = ('pending', 'downloading', 'verifying')
            speed = sum(float(f.get('speed') or 0) for f in files
                        if f.get('status') == 'downloading')
            done_bytes = sum(int(f.get('done_bytes') or 0) for f in files)
            # 总量取 max(size, done)：元数据缺失时 size 来自 aria2c 的
            # 整数 MiB 缩写值，可能略小于实际，避免出现「已下载 36.8MB
            # / 总量 36.0MB」这种倒挂显示。校验仍走 expect_bytes。
            total_bytes = sum(max(int(f.get('size') or 0),
                                  int(f.get('done_bytes') or 0))
                              for f in files)
            remain = sum(max(int(f.get('size') or 0)
                             - int(f.get('done_bytes') or 0), 0)
                         for f in files if f.get('status') in active)
            eta = int(remain / speed) if speed > 1024 else 0
            runs = b.get('runs', {})
            return {'id': b['id'], 'name': b['name'], 'status': b['status'],
                    'created': b['created'], 'error': b.get('error', ''),
                    'concurrency': b.get('concurrency', 2),
                    'converting': b.get('converting', False),
                    'convert_msg': b.get('convert_msg', ''),
                    'n_files': total, 'done': done, 'failed': failed,
                    'n_runs': len(runs),
                    'speed': round(speed, 1), 'eta': eta,
                    'done_bytes': done_bytes, 'total_bytes': total_bytes,
                    'overall': round(min(overall, 1.0), 4), 'files': files}

    def _read_log_tail(self, b, n=40):
        """读取批次过程日志 batch.log 的最后 n 行（读失败返回空列表）。"""
        try:
            d = b.get('dir') or ''
            if not d:
                return []
            with open(os.path.join(d, 'batch.log'), encoding='utf-8',
                      errors='replace') as f:
                return f.read().splitlines()[-n:]
        except OSError:
            return []

    def list_snapshots(self):
        with self.lock:
            ids = list(self.batches.keys())
        outs = []
        for bid in sorted(ids, reverse=True)[:40]:
            s = self.snapshot(bid)
            if s:
                s['files'] = None
                s['log_tail'] = self._read_log_tail(self.batches.get(bid) or {})
                outs.append(s)
        return outs

    def cancel(self, bid):
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return False
            bdir = b['dir']
            for out, proc in list(self.procs.items()):
                if os.path.normpath(out).startswith(
                        os.path.normpath(bdir) + os.sep):
                    try:
                        proc.terminate()
                    except OSError:
                        pass
            for fe in b['files']:
                if fe['status'] in ('pending', 'downloading', 'verifying'):
                    fe['status'] = 'cancelled'
            b['status'] = 'cancelled'
        self._blog(b, '批次已取消')
        self._save(b)
        return True

    def retry_failed(self, bid):
        n = 0
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return 0
            for fe in b['files']:
                if fe['status'] == 'failed' and fe.get('url'):
                    fe['status'] = 'pending'
                    fe['error'] = ''
                    fe['progress'] = 0.0
                    n += 1
            if n:
                b['status'] = 'downloading'
        self._save(b)
        if n:
            self._blog(b, f'重试 {n} 个失败文件')
            self._start_workers(bid)
        return n

    _RECORD_FILES = ('batch.json', 'batch.log')
    _ACTIVE = ('resolving', 'downloading', 'converting')

    def delete(self, bid):
        """全删：批次记录 + 日志 + 已下载文件（目录在创建时已校验）。"""
        self.cancel(bid)
        with self.lock:
            b = self.batches.pop(bid, None)
        if not b:
            return False
        if b.get('dir'):
            shutil.rmtree(b['dir'], ignore_errors=True)
        return True

    def delete_files(self, bid):
        """只删文件：清掉批次目录内的数据文件，保留记录(batch.json)与
        日志(batch.log)；文件条目标记 deleted（列表仍可见、可追溯）。"""
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return False
            if b.get('status') in self._ACTIVE:
                active = True
            else:
                active = False
        if active:
            self.cancel(bid)          # 运行中先停（cancel 会改状态并保存）
        bdir = b.get('dir') or ''
        n = 0
        if bdir and os.path.isdir(bdir):
            for name in os.listdir(bdir):
                if name in self._RECORD_FILES:
                    continue
                p_ = os.path.join(bdir, name)
                try:
                    if os.path.isdir(p_) and not os.path.islink(p_):
                        shutil.rmtree(p_, ignore_errors=True)
                    else:
                        os.remove(p_)
                    n += 1
                except OSError:
                    pass
        with self.lock:
            b = self.batches.get(bid)
            if b:
                for fe in b['files']:
                    if fe.get('status') in ('done', 'failed', 'cancelled',
                                            'pending'):
                        fe['status'] = 'deleted'
                        fe['progress'] = 0.0
                b['files_deleted'] = True
                self._save(b)
        self._blog(b, f'已删除 {n} 个数据文件（保留批次记录与日志）')
        return True

    def delete_record(self, bid):
        """只删记录：移除 batch.json/batch.log（批次从列表消失），已下载
        数据文件保留在 downloads/<批次>/ 磁盘原处。"""
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return False
            active = b.get('status') in self._ACTIVE
        if active:
            self.cancel(bid)
        with self.lock:
            b = self.batches.pop(bid, None)
        if not b:
            return False
        if b.get('dir'):
            for name in self._RECORD_FILES:
                try:
                    os.remove(os.path.join(b['dir'], name))
                except OSError:
                    pass
        return True

    # ---- 已完成文件的清单（进分析流程用） ----
    def ready_files(self, bid):
        snap = self.snapshot(bid)
        if not snap:
            return {}
        by_acc = {}
        for f in snap['files']:
            if f['status'] != 'done' or not f['out']:
                continue
            fn = os.path.basename(f['out'])
            if re.search(r'_1\.f(ast)?q\.gz$', fn, re.I):
                by_acc.setdefault(f['acc'], {})['r1'] = f['out']
            elif re.search(r'_2\.f(ast)?q\.gz$', fn, re.I):
                by_acc.setdefault(f['acc'], {})['r2'] = f['out']
            elif re.search(r'\.f(ast)?q\.gz$', fn, re.I):
                by_acc.setdefault(f['acc'], {}).setdefault('single', f['out'])
        return by_acc


_manager = None


def get_manager():
    global _manager
    if _manager is None:
        _manager = DownloadManager()
    return _manager

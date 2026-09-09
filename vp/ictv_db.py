# -*- coding: utf-8 -*-
"""
ICTV VMR 参考库（databases/tax/ictv/）——官方元数据驱动的按需参考池。

数据流（分层）:
  ① 解析层: VMR 官方 xlsx（ictv.global 下载或手动放入本目录）→ taxa.txt。
     沿用 acvirus_db/taxa.txt 的谱系列格式（acvirus.load_taxa 可直接读），
     另附 Subrealm/Subkingdom/Subphylum/Subclass 完整谱系 + 病毒名/缩写/
     分离物/基因组完整性/Baltimore/宿主组扩展列；多 accession（分段病毒）
     拆成一行一条。
  ② 选参层: 按属/科/种过滤 taxa.txt 挑选参考——本地已有（acvirus_db fasta
     命中或 gb 缓存已下）优先，缺的才走 ③；同层内 Complete genome 优先。
  ③ 下载层: 本地没有的 accession 经 NCBI efetch 拉 GenBank 平文（复用
     gb_collection 的批量/解析机制）落盘 gb_cache/，汇总 gb_refs.fa 供
     ⑦ 建树参考提取。
  ④ 兜底: 网络不可用/下载失败只告警——调用方（phylo ⑦）回退现有
     acvirus_db/植物库参考池，流程不中断。

目录布局:
  VMR_MSL*.xlsx       官方 VMR（ictv-update 下载或手动放入；多份取最新）
  taxa.txt            解析产物（xlsx 更新后自动重建）
  DATA_VERSION        MSL 版本锚点（MSL= / SOURCE= / XLSX= / PARSED=）
  acvirus_accs.txt    acvirus fasta accession 索引缓存（判"本地已有"）
  gb_cache/<acc>.gb   按需下载的 GenBank 平文（幂等缓存）
  gb_refs.fa          缓存序列汇总（⑦ 参考提取数据源之一）
  refs_meta.tsv       缓存记录元数据（acc/物种/宿主/长度等）
"""
import os
import re
import time
import socket
import ipaddress
import urllib.parse
import urllib.request
import urllib.error

from .config import DIRS, db_path
from .utils import check_path, safe_open, write_fasta_record

_DIR = db_path('tax', 'ictv')

# 谱系列（前 16 列；含 acvirus_db/taxa.txt 12 列的全部字段名），其后为扩展列
LINEAGE_FIELDS = ['Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum',
                  'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder',
                  'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species',
                  'Genbank']
EXT_FIELDS = ['Virus_Name', 'Abbreviation', 'Isolate', 'Genome_Coverage',
              'Genome', 'Host_Source', 'MSL']
TAXA_FIELDS = LINEAGE_FIELDS + EXT_FIELDS

# VMR xlsx 表头 → taxa 列名（按表头映射，与列序无关）
_VM_MAP = {
    'Realm': 'Realm', 'Subrealm': 'Subrealm', 'Kingdom': 'Kingdom',
    'Subkingdom': 'Subkingdom', 'Phylum': 'Phylum', 'Subphylum': 'Subphylum',
    'Class': 'Class', 'Subclass': 'Subclass', 'Order': 'Order',
    'Suborder': 'Suborder', 'Family': 'Family', 'Subfamily': 'Subfamily',
    'Genus': 'Genus', 'Subgenus': 'Subgenus', 'Species': 'Species',
    'Virus name(s)': 'Virus_Name',
    'Virus name abbreviation(s)': 'Abbreviation',
    'Virus isolate designation': 'Isolate',
    'Virus GENBANK accession': 'Genbank',
    'Genome coverage': 'Genome_Coverage', 'Genome': 'Genome',
    'Host source': 'Host_Source',
}
# GenBank nucleotide 号（含 .v）：
#   经典  V01234 / MH447526 / NC_010393（1-2 字母 + 可选下划线 + 5-6 数字）
#   WGS   AKVG01000002 / CAJDJZ010000002（4-6 字母 + 2 位版本码 + 6-9 数字）
_ACC_TOKEN_RE = re.compile(
    r'^[A-Z]{1,2}_?\d{5,6}(\.\d+)?$|^[A-Z]{4,6}\d{2}\d{6,9}(\.\d+)?$')
_VMR_XLSX_RE = re.compile(r'^VMR_.+\.xlsx$', re.IGNORECASE)
_MSL_RE = re.compile(r'VMR\s*(MSL\d+)', re.IGNORECASE)

ICTV_HOSTS = {'ictv.global', 'www.ictv.global'}
VMR_URL = 'https://ictv.global/vmr/current?fid=15873'
_MAX_VMR_BYTES = 128 * 1024 * 1024     # VMR xlsx 大小上限（当前 ~4 MB）
_GB_CACHE_CAP = 2000                   # gb 缓存条数上限（防御性）


def db_dir():
    return _DIR if os.path.isdir(_DIR) else None


def _ensure_dir():
    d = check_path(_DIR, must_exist=False, in_platform=True)
    os.makedirs(d, exist_ok=True)
    return d


def vmr_xlsx():
    """目录内最新一份 VMR xlsx（按修改时间取最新）。"""
    if not os.path.isdir(_DIR):
        return None
    cands = [os.path.join(_DIR, f) for f in os.listdir(_DIR)
             if _VMR_XLSX_RE.match(f)]
    if not cands:
        return None
    return max(cands, key=os.path.getmtime)


def taxa_path():
    p = os.path.join(_DIR, 'taxa.txt') if os.path.isdir(_DIR) else None
    return p if p and os.path.isfile(p) else None


def available():
    return taxa_path() is not None


# ------------------------------------------------------------------
# ① 解析层：VMR xlsx → taxa.txt + DATA_VERSION
# ------------------------------------------------------------------
def _msl_version(xlsx):
    """从 xlsx 工作表名/文件名识别 MSL 版本（如 MSL41）。"""
    try:
        import openpyxl
        wb = openpyxl.load_workbook(xlsx, read_only=True)
        for name in wb.sheetnames:
            m = _MSL_RE.search(name or '')
            if m:
                wb.close()
                return m.group(1).upper()
        wb.close()
    except Exception:
        pass
    m = _MSL_RE.search(os.path.basename(xlsx) or '')
    return m.group(1).upper() if m else ''


def _split_accessions(raw):
    """VMR accession 字段 → 合法 GenBank nucleotide 号列表（去重保序）。

    分段/多基因组行一个字段含多个号（空格/逗号/分号分隔），按 token 正则
    过滤——只保留经典 nucleotide 格式，WGS master/注释性文字自动剔除。
    """
    out = []
    for tok in re.split(r'[\s,;]+', str(raw or '').strip()):
        tok = tok.strip().upper().rstrip('.')
        if tok and _ACC_TOKEN_RE.fullmatch(tok) and tok not in out:
            out.append(tok)
    return out


def parse_vmr(xlsx=None, logger=None):
    """解析 VMR xlsx → taxa.txt + DATA_VERSION（换版/更新后调用）。"""
    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    xlsx = xlsx or vmr_xlsx()
    if not xlsx or not os.path.isfile(xlsx):
        raise FileNotFoundError(
            'ictv_db 内没有 VMR xlsx（用 ictv-update 下载或手动放入）')
    import openpyxl
    msl = _msl_version(xlsx)
    log(f'解析 VMR（{msl or "未知版本"}）: {os.path.basename(xlsx)}')

    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    ws = None
    for name in wb.sheetnames:
        if name.lower().startswith('vmr'):
            ws = wb[name]
            break
    if ws is None:
        wb.close()
        raise RuntimeError('VMR xlsx 中未找到 VMR 数据工作表')

    # 表头行定位：含 'Realm' 与 accession 列名的首行
    header = None
    rows = ws.iter_rows(values_only=True)
    for r in rows:
        vals = [str(c).strip() if c is not None else '' for c in (r or [])]
        if 'Realm' in vals and 'Virus GENBANK accession' in vals:
            header = vals
            break
    if header is None:
        wb.close()
        raise RuntimeError('VMR 表头行未识别（缺 Realm/accession 列）')
    col = {h: i for i, h in enumerate(header)}
    rev = {v: k for k, v in _VM_MAP.items()}

    def cell(r, dst):
        i = col.get(rev.get(dst, dst))
        if i is None or i >= len(r) or r[i] is None:
            return ''
        return str(r[i]).strip()

    out_rows, dups, no_acc = [], 0, 0
    seen_acc = set()
    for r in rows:
        if not r or col['Realm'] >= len(r) or not r[col['Realm']]:
            continue
        accs = _split_accessions(cell(r, 'Genbank'))
        if not accs:
            no_acc += 1
            continue
        base = {k: cell(r, k) for k in TAXA_FIELDS if k != 'Genbank'}
        for a in accs:
            if a in seen_acc:
                dups += 1
                continue
            seen_acc.add(a)
            row = dict(base)
            row['Genbank'] = a
            row['MSL'] = msl
            out_rows.append(row)
    wb.close()

    taxa_p = os.path.join(_ensure_dir(), 'taxa.txt')
    with safe_open(taxa_p, 'wt') as f:
        f.write('\t'.join(TAXA_FIELDS) + '\n')
        for row in out_rows:
            f.write('\t'.join((row.get(k) or '').replace('\t', ' ')
                              for k in TAXA_FIELDS) + '\n')
    _write_version(msl, xlsx)
    log(f'taxa.txt 已生成: {len(out_rows)} 条 accession'
        + (f'（跳过重复 {dups}，无有效 accession 行 {no_acc}）'
           if dups or no_acc else ''))
    return taxa_p


def _write_version(msl, xlsx):
    try:
        st = os.stat(xlsx)
        stamp = f'{st.st_size}|{int(st.st_mtime)}'
    except OSError:
        stamp = ''
    with safe_open(os.path.join(_ensure_dir(), 'DATA_VERSION'), 'wt') as f:
        f.write(f'MSL={msl}\n')
        f.write(f'SOURCE={VMR_URL}\n')
        f.write(f'XLSX={os.path.basename(xlsx)}\n')
        f.write(f'XLSX_STAMP={stamp}\n')
        f.write(f'PARSED={time.strftime("%Y-%m-%d %H:%M:%S")}\n')


def _read_version():
    p = os.path.join(_DIR, 'DATA_VERSION')
    out = {}
    if p and os.path.isfile(p):
        with safe_open(p) as f:
            for line in f:
                line = line.strip()
                if '=' in line and not line.startswith('#'):
                    k, v = line.split('=', 1)
                    out[k.strip()] = v.strip()
    return out


def _xlsx_stamp(xlsx):
    try:
        st = os.stat(xlsx)
        return f'{st.st_size}|{int(st.st_mtime)}'
    except OSError:
        return ''


def ensure_taxa(logger=None):
    """taxa.txt 缺失或 xlsx 已更新时自动重建，返回 taxa.txt 路径。"""
    xlsx = vmr_xlsx()
    if not xlsx:
        return taxa_path()
    v = _read_version()
    if taxa_path() and v.get('XLSX_STAMP') == _xlsx_stamp(xlsx):
        return taxa_path()
    return parse_vmr(xlsx, logger=logger)


_taxa_cache = {}


def load_rows(force=False):
    """taxa.txt → [dict]（行序即文件序）。"""
    if _taxa_cache.get('rows') and not force:
        return _taxa_cache['rows']
    p = taxa_path()
    rows = []
    if p:
        import csv
        with safe_open(p) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                if (row.get('Genbank') or '').strip():
                    rows.append(row)
    _taxa_cache['rows'] = rows
    return rows


def load_taxa(force=False):
    """{无版本 accession: 谱系+扩展 dict}——与 acvirus.load_taxa 同构。"""
    if _taxa_cache.get('map') and not force:
        return _taxa_cache['map']
    out = {}
    for row in load_rows(force=force):
        base = row['Genbank'].split('.')[0].upper()
        out[base] = row
    _taxa_cache['map'] = out
    return out


def lineage_for(accession):
    base = (accession or '').split('.')[0].upper()
    return load_taxa().get(base)


def species_for(accession):
    return (lineage_for(accession) or {}).get('Species', '')


# ------------------------------------------------------------------
# ② 选参层：属/科过滤 + 本地优先
# ------------------------------------------------------------------
def acvirus_acc_index(force=False, logger=None):
    """acvirus_db/all_virus.fasta 的 accession 索引（去版本 base 集）。

    首次扫描 586 MB fasta 头部（数十秒级），结果缓存 acvirus_accs.txt，
    由 fasta 大小+时间戳控制失效——选参层据此判"本地已有，零下载"。
    """
    from . import acvirus
    ref = acvirus.fasta_path()
    if not ref:
        return set()
    cache_p = os.path.join(_ensure_dir(), 'acvirus_accs.txt')
    stamp = acvirus._stamp()
    marker = f'#stamp={stamp}'
    if not force and os.path.isfile(cache_p):
        try:
            with safe_open(cache_p) as f:
                if f.readline().strip() == marker:
                    return {line.strip() for line in f if line.strip()}
        except (OSError, ValueError):
            pass
    accs = set()
    with safe_open(ref) as f:
        for line in f:
            if line.startswith('>'):
                accs.add(line[1:].split()[0].split('.')[0].upper())
    with safe_open(cache_p, 'wt') as f:
        f.write(marker + '\n')
        for a in sorted(accs):
            f.write(a + '\n')
    if logger:
        logger.log(f'acvirus accession 索引已缓存: {len(accs)} 条')
    return accs


def _gb_cached_accs():
    cache = gb_cache_dir()
    # 文件名形如 OQ943946.1.gb，统一按无版本 base 归一
    return {f[:-3].split('.')[0].upper() for f in os.listdir(cache)
            if f.upper().endswith('.GB')}


# 级联选参支持的主分类级（列名 → 中文等级）
RANK_COLS = [('Realm', '界'), ('Kingdom', '界下 Kingdom'), ('Phylum', '门'),
             ('Class', '纲'), ('Order', '目'), ('Family', '科'),
             ('Subfamily', '亚科'), ('Genus', '属'), ('Species', '种')]


def rank_options(levels=None, plant_only=True):
    """级联下拉数据：给定已选等级 {列名: 值}，返回每个等级的
    [{name, n}]（该等级的独立值 + accession 条数，n 降序）。

    plant_only=True 限定宿主含 plant 的记录（平台定位植物病毒）。
    """
    levels = {k: str(v).strip().lower()
              for k, v in (levels or {}).items() if str(v or '').strip()}
    out = {}
    for r in load_rows():
        if plant_only and 'plant' not in (r.get('Host_Source') or '').lower():
            continue
        ok = True
        for col, val in levels.items():
            if (r.get(col) or '').strip().lower() != val:
                ok = False
                break
        if not ok:
            continue
        for col, _zh in RANK_COLS:
            v = (r.get(col) or '').strip()
            if v:
                out.setdefault(col, {})
                out[col][v] = out[col].get(v, 0) + 1
    return {col: [{'name': k, 'n': n} for k, n in
                  sorted(d.items(), key=lambda x: (-x[1], x[0]))]
            for col, d in out.items()}


def cap_per_rank(rows, per_genus=0, per_species=0):
    """按属/物种上限抽样（保持传入顺序 = 本地优先序）。

    per_genus / per_species: 每属/每种最多保留的 accession 数，0=不限。
    高等级（科/目…）整包下载时用它防单个大属（如 Begomovirus 846 条）
    挤占 limit 配额、把集合淹成单属集。返回 (kept_rows, dropped)。
    """
    per_genus = max(0, int(per_genus or 0))
    per_species = max(0, int(per_species or 0))
    if not per_genus and not per_species:
        return list(rows), 0
    g_cnt, s_cnt, kept, dropped = {}, {}, [], 0
    for r in rows:
        g = (r.get('Genus') or '').strip().lower() or '(无属)'
        s = (r.get('Species') or '').strip().lower() or '(无种)'
        if per_species and s_cnt.get(s, 0) >= per_species:
            dropped += 1
            continue
        if per_genus and g_cnt.get(g, 0) >= per_genus:
            dropped += 1
            continue
        g_cnt[g] = g_cnt.get(g, 0) + 1
        s_cnt[s] = s_cnt.get(s, 0) + 1
        kept.append(r)
    return kept, dropped


def select_refs(genus=None, family=None, species=None, limit=20,
                genome='complete', logger=None, ranks=None, plant_only=False):
    """按分类级挑选参考（选参层主入口）。

    兼容旧参 genus/family/species；ranks={列名: 值} 支持界门纲目科属种
    任意组合过滤（与 genus/family/species 可并用，取交集）。
    plant_only=True 只取宿主含 plant 的记录（Web 端级联选参口径）。
    返回 (rows, total)：rows 按优先级排序后取前 limit 条，每行附
    Source（acvirus_db/gb_cache/ncbi——ncbi 表示需下载）；排序键：
    本地已有 > Complete genome > accession 字典序（结果确定可复现）。
    genome: 'complete' 只取 Complete genome；'any' 不过滤。
    """
    g = (genus or '').strip().lower()
    fam = (family or '').strip().lower()
    sp = (species or '').strip().lower()
    rank_f = {k: str(v).strip().lower()
              for k, v in (ranks or {}).items() if str(v or '').strip()}
    if not (g or fam or sp or rank_f):
        raise ValueError('至少指定一个分类级')
    rows = load_rows()
    local = acvirus_acc_index(logger=logger)
    cached = _gb_cached_accs()
    total = 0
    sel = []
    for r in rows:
        if plant_only and 'plant' not in (r.get('Host_Source') or '').lower():
            continue
        if g and (r.get('Genus') or '').strip().lower() != g:
            continue
        if fam and (r.get('Family') or '').strip().lower() != fam:
            continue
        if sp and (r.get('Species') or '').strip().lower() != sp:
            continue
        if rank_f:
            skip = False
            for col, val in rank_f.items():
                if (r.get(col) or '').strip().lower() != val:
                    skip = True
                    break
            if skip:
                continue
        total += 1
        cov = (r.get('Genome_Coverage') or '').lower()
        if genome == 'complete' and cov and not cov.startswith('complete'):
            continue
        base = r['Genbank'].split('.')[0].upper()
        in_acv = base in local
        in_cache = base in cached
        src = 'acvirus_db' if in_acv else ('gb_cache' if in_cache else 'ncbi')
        sel.append({**r, 'Source': src,
                    'Priority': (0 if in_acv else 1 if in_cache else 2,
                                 0 if cov.startswith('complete') else 1,
                                 r['Genbank'])})
    sel.sort(key=lambda r: r['Priority'])
    return sel[:max(1, int(limit))], total


# ------------------------------------------------------------------
# ③ 下载层：缺的 accession 经 NCBI efetch 落盘缓存
# ------------------------------------------------------------------
def gb_cache_dir():
    d = os.path.join(_ensure_dir(), 'gb_cache')
    os.makedirs(d, exist_ok=True)
    return d


def _refs_fa_path():
    return os.path.join(_DIR, 'gb_refs.fa')


def _gb_fasta_stamp_path():
    return os.path.join(_DIR, 'gb_refs.stamp')


def _rebuild_gb_fasta(logger=None):
    """gb_cache 全部 .gb → gb_refs.fa + refs_meta.tsv（条数变化才重建）。"""
    from Bio import SeqIO
    from .gb_collection import _record_summary, _META_FIELDS
    cache = gb_cache_dir()
    gbs = sorted(f for f in os.listdir(cache) if f.upper().endswith('.GB'))
    if len(gbs) > _GB_CACHE_CAP:
        raise RuntimeError(f'gb 缓存 {len(gbs)} 条超过上限 {_GB_CACHE_CAP}，'
                           f'请清理 {cache}')
    stamp = str(len(gbs))
    fa_p = _refs_fa_path()
    if not gbs:
        return None
    if os.path.isfile(fa_p) and os.path.isfile(_gb_fasta_stamp_path()):
        try:
            with safe_open(_gb_fasta_stamp_path()) as f:
                if f.read().strip() == stamp:
                    return fa_p
        except (OSError, ValueError):
            pass
    meta = []
    with safe_open(fa_p, 'wt') as out:
        for fn in gbs:
            path = os.path.join(cache, fn)
            try:
                for rec in SeqIO.parse(path, 'genbank'):
                    if not rec.seq:
                        continue
                    m = _record_summary(rec)
                    org = (m.get('organism') or 'unknown').replace(' ', '_')
                    write_fasta_record(out, f"{m['acc']}|{org}",
                                       str(rec.seq).upper())
                    meta.append({k: m.get(k, '') for k in _META_FIELDS})
            except Exception as e:
                if logger:
                    logger.log(f'{fn} 解析失败，跳过: {e}', 'WARN')
    with safe_open(_gb_fasta_stamp_path(), 'wt') as f:
        f.write(stamp)
    with safe_open(os.path.join(_DIR, 'refs_meta.tsv'), 'wt') as f:
        f.write('\t'.join(_META_FIELDS) + '\n')
        for m in meta:
            f.write('\t'.join(str(m.get(k, '')) for k in _META_FIELDS) + '\n')
    if logger:
        logger.log(f'gb 缓存 fasta 重建: {len(meta)} 条 -> '
                   f'{os.path.basename(fa_p)}')
    return fa_p


def ensure_gb(accs, logger=None):
    """把 accession 列表补齐到本地：已有（acvirus/gb 缓存）跳过，缺的下载。

    返回 {downloaded, skipped, missing, fasta}；fasta 非 None 时可直接
    加入 ⑦ 的参考来源列表。无网络/下载失败时抛异常，由调用方兜底。
    """
    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    from .gb_collection import _fetch_accessions_gb, _safe_acc_name
    cache = gb_cache_dir()
    local = acvirus_acc_index()
    have = _gb_cached_accs()
    need, skipped = [], 0
    for a in accs:
        base = (a or '').split('.')[0].upper()
        if not base:
            continue
        if base in local or base in have:
            skipped += 1
            continue
        if base not in need:
            need.append(base)
    fa = _rebuild_gb_fasta(logger=logger)
    if not need:
        return {'downloaded': 0, 'skipped': skipped, 'missing': [],
                'fasta': fa}
    log(f'ICTV 参考按需下载: {len(need)} 条（本地已有跳过 {skipped}）')
    pairs, missing = _fetch_accessions_gb(need, lambda m: log(m))
    downloaded = 0
    for idx, (chunk, m) in enumerate(pairs):
        base = m['acc'].split('.')[0].upper()
        if base in local or base in have:
            continue
        with safe_open(os.path.join(cache, _safe_acc_name(m['acc'], idx)),
                       'wt') as f:
            f.write(chunk)
        have.add(base)
        downloaded += 1
    if missing:
        log(f'NCBI 未找到 {len(missing)} 条: {", ".join(missing[:10])}'
            + ('…' if len(missing) > 10 else ''), 'WARN')
        with safe_open(os.path.join(_DIR, 'missing_accs.txt'), 'at') as f:
            f.write(time.strftime('# %Y-%m-%d %H:%M:%S\n'))
            for a in missing:
                f.write(a + '\n')
    if downloaded:
        fa = _rebuild_gb_fasta(logger=logger)
    log(f'下载落盘 {downloaded} 条，缓存合计 {len(have)} 条')
    return {'downloaded': downloaded, 'skipped': skipped,
            'missing': missing, 'fasta': fa}


# ------------------------------------------------------------------
# VMR 在线更新（ictv.global，域名白名单 + IP 校验 + 大小上限）
# ------------------------------------------------------------------
def _check_host_ip(hostname):
    """域名不得解析到私网/环回等保留地址（防 SSRF，同 ncbi_download）。"""
    for info in socket.getaddrinfo(hostname, None):
        ip = ipaddress.ip_address(info[4][0])
        if (ip.is_private or ip.is_loopback or ip.is_reserved
                or ip.is_link_local or ip.is_multicast):
            raise RuntimeError(f'域名解析到保留地址，已拒绝: {ip}')


def _open_ictv(url, timeout=120, retries=3):
    """白名单域名的 HTTPS GET（自动重定向，但最终必须仍在 ictv.global）。

    返回 (字节内容, 最终 URL)。响应必须是 xlsx，超过大小上限即中止。
    """
    u = urllib.parse.urlparse(url)
    if u.scheme != 'https' or (u.hostname or '').lower() not in ICTV_HOSTS:
        raise RuntimeError(f'仅允许 https://ictv.global: {url}')
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            _check_host_ip(u.hostname.lower())
            req = urllib.request.Request(url, headers={
                'User-Agent': 'PlantVirusPlatform/1.0 (+ictv-db-updater)'})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                final = urllib.parse.urlparse(r.geturl() or url)
                if (final.hostname or '').lower() not in ICTV_HOSTS:
                    raise RuntimeError(f'重定向跳出白名单域: {r.geturl()}')
                ctype = (r.headers.get('Content-Type') or '').lower()
                data = b''
                while True:
                    chunk = r.read(1024 * 1024)
                    if not chunk:
                        break
                    data += chunk
                    if len(data) > _MAX_VMR_BYTES:
                        raise RuntimeError(
                            f'VMR 下载超过 {_MAX_VMR_BYTES} 字节上限')
            if 'spreadsheetml' not in ctype and 'octet-stream' not in ctype:
                raise RuntimeError(f'响应不是 xlsx（Content-Type={ctype}），'
                                   f'请到 ictv.global 手动下载后放入 {_DIR}')
            return data, r.geturl()
        except (urllib.error.URLError, socket.timeout, OSError) as e:
            last_err = str(e)
            if attempt < retries:
                time.sleep(2 * attempt)
    raise RuntimeError(f'VMR 下载失败（重试 {retries} 次）: {last_err}')


def download_current_vmr(logger=None):
    """下载 ictv.global 当前版 VMR xlsx 到 ictv_db（保留最新，旧的删除）。"""
    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    data, final_url = _open_ictv(VMR_URL)
    name = os.path.basename(urllib.parse.urlparse(final_url).path)
    if not _VMR_XLSX_RE.match(name):
        name = 'VMR_current.xlsx'
    d = _ensure_dir()
    path = check_path(os.path.join(d, name), must_exist=False,
                      in_platform=True)
    with open(path, 'wb') as f:
        f.write(data)
    log(f'VMR 已下载: {name}（{len(data) / 1048576:.1f} MB）')
    # 只保留最新一份，避免 vmr_xlsx() 取最新时产生歧义
    for old in os.listdir(d):
        p = os.path.join(d, old)
        if _VMR_XLSX_RE.match(old) \
                and os.path.abspath(p) != os.path.abspath(path):
            try:
                os.remove(p)
                log(f'已移除旧版: {old}')
            except OSError:
                pass
    return path


def update(logger=None, xlsx=None):
    """更新主入口：可选在线下载 → 解析 taxa.txt → 返回 taxa.txt 路径。"""
    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    if xlsx:
        src = check_path(xlsx, must_exist=True, in_platform=False)
        d = _ensure_dir()
        dst = os.path.join(d, os.path.basename(src))
        if os.path.abspath(src) != os.path.abspath(dst):
            with open(src, 'rb') as fi, open(dst, 'wb') as fo:
                fo.write(fi.read())
        taxa_p = parse_vmr(dst, logger=logger)
    else:
        try:
            download_current_vmr(logger=logger)
        except Exception as e:
            log(f'在线下载失败（改用本地已有 xlsx 继续解析）: {e}', 'WARN')
        taxa_p = ensure_taxa(logger=logger)
        if not taxa_p:
            raise FileNotFoundError('ictv_db 内没有可解析的 VMR xlsx')
    n_gen, n_fam = _taxa_stats()
    log(f'ICTV 参考库就绪: {_read_version().get("MSL", "?")}, '
        f'{len(load_rows(force=True))} 条 accession, '
        f'{n_gen} 属 / {n_fam} 科')
    return taxa_p


def _taxa_stats():
    gens, fams = set(), set()
    for r in load_rows():
        if r.get('Genus'):
            gens.add(r['Genus'])
        if r.get('Family'):
            fams.add(r['Family'])
    return len(gens), len(fams)


def status():
    """库状态摘要（GUI/CLI 用）。"""
    v = _read_version()
    n_gb = len(_gb_cached_accs()) if os.path.isdir(gb_cache_dir()) else 0
    n_gen, n_fam = _taxa_stats()
    n_local = 0
    if taxa_path():
        try:
            local = acvirus_acc_index()
            n_local = sum(1 for r in load_rows()
                          if r['Genbank'].split('.')[0].upper() in local)
        except Exception:
            pass
    return {
        'dir': db_dir(),
        'msl': v.get('MSL', ''),
        'xlsx': v.get('XLSX', ''),
        'parsed': v.get('PARSED', ''),
        'n_rows': len(load_rows()) if taxa_path() else 0,
        'n_genus': n_gen, 'n_family': n_fam,
        'n_in_acvirus': n_local,
        'n_gb_cache': n_gb,
        'taxa_ready': available(),
        'has_xlsx': vmr_xlsx() is not None,
    }

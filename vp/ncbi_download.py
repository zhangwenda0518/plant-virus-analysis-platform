# -*- coding: utf-8 -*-
"""
NCBI Entrez 参考序列批量下载（借鉴 PhyloSuite 的会话式翻页设计）。

流程：esearch(usehistory=y) 拿 WebEnv/QueryKey 会话 → esummary 预览元数据
→ efetch 按 retstart/retmax 分批拉取（不需重传 ID 列表）→ GenBank 解析
提取序列+宿主等元数据 → 按 accession.version 去重落盘。

与 PhyloSuite 差异（按平台安全模型重写）：
- 域名白名单（仅 eutils.ncbi.nlm.nih.gov）+ 解析 IP 阻断私网/环回
- 请求间隔限速（无 api_key 3 req/s）+ 指数退避重试（PhyloSuite 无重试）
- 断点续传：目录内已有 refs.fa 时按 accession 跳过，重复运行幂等
- 元数据 TSV（物种/宿主/分离物/国家/日期/taxid）供建树命名与溯源

集合目录：databases/ncbi_refs/<name>/（refs.fa + refs_meta.tsv + query.json）
"""
import io
import os
import re
import json
import time
import socket
import ipaddress
import urllib.parse
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET

from .config import DIRS, PLATFORM_ROOT, get_config
from .utils import check_path, safe_open, write_fasta_record

EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
_EUTILS_HOST = 'eutils.ncbi.nlm.nih.gov'
_RETRYABLE_HTTP = {429, 500, 502, 503, 504}
_MAX_RETRIES = 5
_MIN_INTERVAL = 0.35          # 请求最小间隔秒（≤3 req/s，NCBI 无 key 限速）
FETCH_BATCH = 200             # efetch 每批记录数
SUMMARY_BATCH = 100           # esummary 每批 ID 数
_MAX_XML_BYTES = 64 * 1024 * 1024   # 单响应 XML 大小上限

_last_request = [0.0]         # 模块级限速时钟


def _open_eutils(url, params, timeout=120):
    """带白名单/IP 校验/限速/重试的 EUtils 请求，返回响应文本。"""
    u = urllib.parse.urlparse(url)
    if u.scheme != 'https' or (u.hostname or '').lower() != _EUTILS_HOST:
        raise RuntimeError(f'仅允许 https://{_EUTILS_HOST}: {url}')
    for info in socket.getaddrinfo(_EUTILS_HOST, None):
        ip = ipaddress.ip_address(info[4][0])
        if (ip.is_private or ip.is_loopback or ip.is_reserved
                or ip.is_link_local or ip.is_multicast):
            raise RuntimeError(f'域名解析到保留地址，已拒绝: {ip}')

    common = {'tool': 'PlantVirusPlatform', 'usehistory': 'y'}
    email = getattr(get_config(), 'email', '')
    if email:
        common['email'] = email
    data = urllib.parse.urlencode({**common, **params}).encode('ascii')

    delay = 1.0
    last_err = None
    for attempt in range(1, _MAX_RETRIES + 1):
        wait = _MIN_INTERVAL - (time.monotonic() - _last_request[0])
        if wait > 0:
            time.sleep(wait)
        _last_request[0] = time.monotonic()
        req = urllib.request.Request(url, data=data)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode('utf-8', errors='replace')
        except urllib.error.HTTPError as e:
            last_err = f'HTTP {e.code}'
            if e.code not in _RETRYABLE_HTTP:
                raise RuntimeError(f'EUtils HTTP {e.code}: {e.reason}')
        except (urllib.error.URLError, socket.timeout, OSError) as e:
            last_err = str(e)
        if attempt < _MAX_RETRIES:
            time.sleep(delay)
            delay = min(delay * 2, 30)
    raise RuntimeError(f'EUtils 请求失败（重试 {_MAX_RETRIES} 次）: {last_err}')


def _safe_xml(text, what='EUtils'):
    """安全 XML 解析：优先 defusedxml（默认禁止 DTD/实体扩展），大小上限兜底。

    defusedxml 不可用时回退 stdlib：拒绝 <!ENTITY 内部定义，并剥离响应中
    NCBI 合法的外部 DOCTYPE 声明（ElementTree 不加载外部 DTD）后再解析；
    仅对白名单域名 eutils.ncbi.nlm.nih.gov 的响应用此回退路径。
    """
    if len(text) > _MAX_XML_BYTES:
        raise RuntimeError(f'{what} 响应超过 {_MAX_XML_BYTES} 字节上限')
    if '<!ENTITY' in text:
        raise RuntimeError(f'{what} 响应含内部实体定义，已拒绝解析')
    try:
        from defusedxml.ElementTree import fromstring as safe_fromstring
        return safe_fromstring(text)
    except ImportError:
        if '<!DOCTYPE' in text[:4096]:
            stripped = _DOCTYPE_RE.sub('', text, count=1)
            if '<!DOCTYPE' in stripped or '<!ENTITY' in stripped:
                raise RuntimeError(f'{what} 响应 DOCTYPE 异常，已拒绝解析')
            text = stripped
        return ET.fromstring(text)


_DOCTYPE_RE = re.compile(r'<!DOCTYPE[^>\[]*(\[[^\]]*\])?[^>]*>', re.IGNORECASE)


# ------------------------------------------------------------------
# esearch / esummary
# ------------------------------------------------------------------
def esearch(term, db='nucleotide', retmax=0):
    """Entrez 检索（usehistory=y 会话）。

    返回 {count, webenv, query_key, ids}。count=总命中数；
    retmax=0 时不取 ID 列表（仅会话，供 efetch 翻页用）。
    """
    params = {'db': db, 'term': term, 'retmax': str(retmax)}
    xml_text = _open_eutils(f'{EUTILS_BASE}/esearch.fcgi', params)
    root = _safe_xml(xml_text, 'esearch')
    err = root.findtext('ErrorList/PhraseNotFound')
    if root.find('ErrorList') is not None and not root.findall('IdList/Id'):
        raise RuntimeError(f'检索无结果或查询词错误: {err or term}')
    return {
        'count': int(root.findtext('Count') or 0),
        'webenv': root.findtext('WebEnv') or '',
        'query_key': root.findtext('QueryKey') or '',
        'ids': [e.text for e in root.findall('IdList/Id')],
    }


def summarize_ids(ids, db='nucleotide'):
    """按 ID 批量取摘要 → [{acc, title, organism, taxid, length, updated}]。"""
    rows = []
    for i in range(0, len(ids), SUMMARY_BATCH):
        batch = ids[i:i + SUMMARY_BATCH]
        xml_text = _open_eutils(f'{EUTILS_BASE}/esummary.fcgi',
                                {'db': db, 'id': ','.join(batch),
                                 'version': '2.0'})   # 2.0 = DocumentSummarySet 格式
        root = _safe_xml(xml_text, 'esummary')
        for doc in root.findall('DocumentSummarySet/DocumentSummary'):
            rows.append({
                'acc': doc.findtext('Caption') or '',
                'title': doc.findtext('Title') or '',
                'organism': doc.findtext('Organism') or '',
                'taxid': doc.findtext('TaxId') or '',
                'length': doc.findtext('Slen') or doc.findtext('Length') or '',
                'updated': doc.findtext('UpdateDate') or '',
            })
    return rows


def search_preview(term, db='nucleotide', limit=50):
    """检索 + 前 limit 条元数据预览（GUI 表格用）。"""
    limit = max(1, min(int(limit), 500))
    res = esearch(term, db=db, retmax=limit)
    rows = summarize_ids(res['ids'], db=db) if res['ids'] else []
    return {'count': res['count'], 'rows': rows}


# ------------------------------------------------------------------
# efetch 下载 + GenBank 解析
# ------------------------------------------------------------------
def _fetch_batch_gb(session, db, retstart, retmax):
    return _open_eutils(f'{EUTILS_BASE}/efetch.fcgi', {
        'db': db, 'rettype': 'gb', 'retmode': 'text',
        'query_key': session['query_key'], 'WebEnv': session['webenv'],
        'retstart': str(retstart), 'retmax': str(retmax)})


def parse_gb_records(gb_text):
    """解析 GenBank flatfile → 元数据 dict 列表（Bio.SeqIO）。

    提取：accession.version、definition、organism、taxid、length、
    host、isolate、country、collection_date（source feature qualifiers）。
    """
    from Bio import SeqIO
    out = []
    for rec in SeqIO.parse(io.StringIO(gb_text), 'genbank'):
        src = {}
        taxid = ''
        for feat in rec.features:
            if feat.type == 'source':
                for k in ('host', 'isolate', 'country', 'collection_date'):
                    v = feat.qualifiers.get(k, [''])[0]
                    if v:
                        src.setdefault(k, v)
                for xref in feat.qualifiers.get('db_xref', []):
                    if xref.startswith('taxon:'):
                        taxid = xref.split(':', 1)[1]
                        break
                break
        out.append({
            'acc': rec.id,
            'title': rec.description,
            'organism': rec.annotations.get('organism', ''),
            'taxid': taxid,
            'length': len(rec.seq),
            'host': src.get('host', ''),
            'isolate': src.get('isolate', ''),
            'country': src.get('country', ''),
            'date': rec.annotations.get('date', ''),
            'seq': str(rec.seq).upper(),
        })
    return out


META_FIELDS = ['acc', 'title', 'organism', 'taxid', 'length', 'host',
               'isolate', 'country', 'date']


def _collections_root():
    d = check_path(os.path.join(DIRS['databases'], 'ncbi_refs'), in_platform=True)
    os.makedirs(d, exist_ok=True)
    return d


def collection_dir(name):
    """集合目录（名称白名单校验，防路径注入）。"""
    import re as _re
    if not _re.fullmatch(r'[A-Za-z0-9_\-\u4e00-\u9fff]{1,64}', str(name)):
        raise ValueError(f'集合名仅允许中英文/数字/_/-（≤64字符）: {name!r}')
    return check_path(os.path.join(_collections_root(), name), in_platform=True)


def list_collections():
    """已下载集合列表 [{name, n_seqs, query, date, dir}]。"""
    out = []
    root = _collections_root()
    for name in sorted(os.listdir(root)):
        cdir = os.path.join(root, name)
        if not os.path.isdir(cdir):
            continue
        fa = os.path.join(cdir, 'refs.fa')
        n = 0
        if os.path.isfile(fa):
            with safe_open(fa) as f:
                n = sum(1 for line in f if line.startswith('>'))
        info = {}
        qfile = os.path.join(cdir, 'query.json')
        if os.path.isfile(qfile):
            try:
                with safe_open(qfile) as f:
                    info = json.load(f)
            except (OSError, ValueError):
                pass
        out.append({'name': name, 'n_seqs': n, 'query': info.get('term', ''),
                    'db': info.get('db', ''), 'date': info.get('date', ''),
                    'dir': cdir})
    return out


def _existing_accessions(fa_path):
    accs = set()
    if os.path.isfile(fa_path):
        with safe_open(fa_path) as f:
            for line in f:
                if line.startswith('>'):
                    accs.add(line[1:].split('|')[0].split()[0])
    return accs


def download_collection(term, name, db='nucleotide', max_records=100,
                        logger=None, prog=None, cancel=None):
    """检索并下载一个参考序列集合（幂等：重复运行按 accession 跳过已有）。

    返回 {name, dir, downloaded, skipped, total_hits}。
    下载 GenBank 并解析出宿主等元数据；集合内已有 refs.fa 时续传。
    """
    if db not in ('nucleotide', 'protein'):
        raise ValueError(f'不支持的数据库: {db}')
    max_records = max(1, min(int(max_records), 10000))
    cdir = collection_dir(name)
    fa_path = os.path.join(cdir, 'refs.fa')
    meta_path = os.path.join(cdir, 'refs_meta.tsv')

    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    if prog:
        prog('search', 0.05, 'Entrez 检索中')
    session = esearch(term, db=db)
    total_hits = session['count']
    n_target = min(total_hits, max_records)
    log(f"检索 [{term}] 命中 {total_hits} 条，下载前 {n_target} 条 (db={db})")

    existing = _existing_accessions(fa_path)
    if existing:
        log(f"续传模式: 已有 {len(existing)} 条，重复项将跳过")

    with safe_open(os.path.join(cdir, 'query.json'), 'wt') as f:
        json.dump({'term': term, 'db': db, 'max': max_records,
                   'date': time.strftime('%Y-%m-%d %H:%M:%S')}, f,
                  ensure_ascii=False, indent=2)

    downloaded = skipped = 0
    new_meta = []
    retstart = 0
    while retstart < n_target:
        if cancel is not None and getattr(cancel, 'is_set', lambda: False)():
            log('用户取消下载', 'WARN')
            break
        batch_max = min(FETCH_BATCH, n_target - retstart)
        if prog:
            prog('fetch', 0.1 + 0.85 * retstart / n_target,
                 f'下载 {min(retstart + batch_max, n_target)}/{n_target}')
        gb_text = _fetch_batch_gb(session, db, retstart, batch_max)
        records = parse_gb_records(gb_text)
        if not records:
            log(f'批次 retstart={retstart} 解析到 0 条记录，停止', 'WARN')
            break
        # 追加写 fasta + 收集元数据（按 accession 去重）
        with safe_open(fa_path, 'at') as ffa:
            for r in records:
                acc = r['acc'].split('.')[0]
                if r['acc'] in existing or acc in existing:
                    skipped += 1
                    continue
                org = r['organism'].replace(' ', '_') or 'unknown'
                write_fasta_record(ffa, f"{r['acc']}|{org}", r['seq'])
                existing.add(r['acc'])
                existing.add(acc)
                downloaded += 1
                new_meta.append({k: r.get(k, '') for k in META_FIELDS})
        retstart += batch_max

    # 元数据 TSV（追加模式，与 fasta 同步增长）
    if new_meta:
        write_header = not os.path.isfile(meta_path)
        with safe_open(meta_path, 'at') as f:
            if write_header:
                f.write('\t'.join(META_FIELDS) + '\n')
            for m in new_meta:
                f.write('\t'.join(str(m.get(k, '')) for k in META_FIELDS) + '\n')

    if prog:
        prog('done', 1.0, f'完成: 新增 {downloaded}，跳过 {skipped}')
    log(f"下载完成: 新增 {downloaded} 条，跳过(已有/重复) {skipped} 条 "
        f"-> {os.path.basename(fa_path)}")
    return {'name': name, 'dir': cdir, 'downloaded': downloaded,
            'skipped': skipped, 'total_hits': total_hits}

# -*- coding: utf-8 -*-
"""
Contig 深度注释（参考 plant_virus_db_pipeline/9.metabuli 的结果整理方式）：
  1) NCBI Taxonomy 8 级谱系统列（realm→kingdom→phylum→class→order→
     family→genus→species，从 nodes.dmp/names.dmp 实时走父链）；
  2) 属平均基因组长度（genus_lens）：由病毒参考库自算并缓存，
     length/genus_avg ∈ [0.8, 1.2] 判定为「近完整基因组」；
  3) NCBI 在线接口：BLASTN/BLASTX（URL API，virus-restricted）与
     CDD 保守域（Batch CD-Search，6-frame 最长 ORF 翻译）。

对外网络仅访问 NCBI 官方域名，且逐次校验解析出的 IP 非私有/保留地址。
"""
import os
import re
import glob
import gzip
import time
import json
import socket
import ipaddress
import urllib.parse
import urllib.request

from .config import DIRS, get_config, db_path
from .utils import check_path, safe_open

RANKS = ['realm', 'kingdom', 'phylum', 'class',
         'order', 'family', 'genus', 'species']

NCBI_BLAST_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
NCBI_CDD_URL = 'https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi'
NCBI_BLAST_SLEEP = 10
NCBI_CDD_SLEEP = 10
BLAST_TIMEOUT = 2400
CDD_TIMEOUT = 1800
BLAST_HITLIST = 20
_NCBI_HOSTS = {'blast.ncbi.nlm.nih.gov', 'www.ncbi.nlm.nih.gov',
               'eutils.ncbi.nlm.nih.gov'}


# ------------------------------------------------------------------
# 对外请求安全校验：仅 http(s) + 域名白名单 + 解析 IP 非私有/保留
# ------------------------------------------------------------------
def open_checked(url, data=None, headers=None, timeout=60):
    """带安全校验的 urllib 打开（仅允许 NCBI 白名单域名）。"""
    u = urllib.parse.urlparse(url)
    if u.scheme not in ('http', 'https'):
        raise RuntimeError(f'仅允许 http/https: {u.scheme}')
    host = (u.hostname or '').lower()
    if host not in _NCBI_HOSTS:
        raise RuntimeError(f'非白名单域名: {host}')
    for info in socket.getaddrinfo(host, None):
        ip = ipaddress.ip_address(info[4][0])
        if (ip.is_private or ip.is_loopback or ip.is_reserved
                or ip.is_link_local or ip.is_multicast):
            raise RuntimeError(f'域名解析到保留地址，已拒绝: {host} -> {ip}')
    req = urllib.request.Request(url, data=data)
    for k, v in (headers or {}).items():
        req.add_header(k, v)
    return urllib.request.urlopen(req, timeout=timeout)


# ------------------------------------------------------------------
# NCBI Taxonomy：nodes.dmp/names.dmp → 8 级谱系
# ------------------------------------------------------------------
_tax_cache = {}


def _taxonomy_dir():
    return check_path(DIRS['taxonomy'], must_exist=True, in_platform=True)


def _load_taxonomy(force=False):
    """解析 nodes/names dmp → {taxid: (parent, rank)}, {taxid: name}。"""
    if _tax_cache.get('nodes') and not force:
        return _tax_cache['nodes'], _tax_cache['names']
    tdir = _taxonomy_dir()
    nodes, names = {}, {}
    nd = os.path.join(tdir, 'nodes.dmp')
    nm = os.path.join(tdir, 'names.dmp')
    with open(nd, encoding='utf-8', errors='replace') as f:
        for line in f:
            parts = [x.strip() for x in line.split('|')]
            if len(parts) >= 3 and parts[0]:
                nodes[parts[0]] = (parts[1], parts[2])
    with open(nm, encoding='utf-8', errors='replace') as f:
        for line in f:
            parts = [x.strip() for x in line.split('|')]
            if len(parts) >= 4 and parts[3] == 'scientific name':
                names[parts[0]] = parts[1]
    _tax_cache['nodes'], _tax_cache['names'] = nodes, names
    return nodes, names


def lineage_ranks(taxid):
    """taxid → {rank: name}（只保留 8 个目标级别）+ sciname。"""
    nodes, names = _load_taxonomy()
    out = {r: '' for r in RANKS}
    cur = str(taxid)
    hops = 0
    while cur and cur != '0' and cur in nodes and hops < 80:
        parent, rank = nodes[cur]
        if rank in out and not out[rank]:
            out[rank] = names.get(cur, '')
        cur = parent
        hops += 1
    out['sciname'] = names.get(str(taxid), '')
    return out


# ------------------------------------------------------------------
# 属平均基因组长度（由病毒参考库自算，缓存 databases/genus_lens.tsv）
# ------------------------------------------------------------------
def _seqid2taxid():
    """virus_db/seqid2taxid.map → {seq_id: taxid}。

    map 键形如 `NC_xxx|kraken:taxid|N`；同时登记 `|` 前的裸 accession，
    便于与原始参考 FASTA 头匹配。
    """
    m = {}
    p = os.path.join(db_path('virus', 'plant'), 'seqid2taxid.map')
    if not os.path.isfile(p):
        return m
    with safe_open(p) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                key, taxid = parts[0], parts[1]
                m[key] = taxid
                base = key.split('|')[0]
                if base and base not in m:
                    m[base] = taxid
    return m


def _iter_fasta(path):
    op = gzip.open if str(path).lower().endswith('.gz') else open
    name, buf, n = None, [], 0
    with op(path, 'rt', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, n
                name = line[1:].split()[0]
                buf, n = [], 0
            else:
                n += len(line.strip())
        if name is not None:
            yield name, n


def genus_lens_path():
    return os.path.join(DIRS['databases'], 'genus_lens.tsv')


def build_genus_lens(force=False, logger=None):
    """按病毒参考库统计 属→平均基因组长度，缓存 TSV。返回 {genus: avg}。"""
    cache = genus_lens_path()
    if os.path.isfile(cache) and not force:
        out = {}
        with safe_open(cache) as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) == 2 and parts[0] != 'genus':
                    try:
                        out[parts[0]] = float(parts[1])
                    except ValueError:
                        pass
        return out

    from .assembly import find_virus_ref_fasta
    ref = find_virus_ref_fasta()
    s2t = _seqid2taxid()
    nodes, _names = _load_taxonomy()
    acc = {}
    for sid, length in _iter_fasta(ref):
        taxid = s2t.get(sid)
        if not taxid:
            continue
        cur, hops = str(taxid), 0
        while cur and cur in nodes and hops < 80:
            parent, rank = nodes[cur]
            if rank == 'genus':
                g = _tax_cache['names'].get(cur, cur)
                acc.setdefault(g, []).append(length)
                break
            cur = parent
            hops += 1

    out = {g: sum(v) / len(v) for g, v in acc.items() if v}
    with safe_open(cache, 'wt') as f:
        f.write('genus\ttotal\n')
        for g, avg in sorted(out.items()):
            f.write(f'{g}\t{avg:.2f}\n')
    if logger:
        logger.log(f'属平均长度表已生成: {len(out)} 属 -> {cache}')
    return out


def genus_avg_map(logger=None):
    """带缓存的属平均长度查询表（不足时自动构建）。"""
    m = build_genus_lens(logger=logger)
    if not m:
        m = build_genus_lens(force=True, logger=logger)
    return m


# ------------------------------------------------------------------
# 分类表整理（metabuli 风格）
# ------------------------------------------------------------------
def classify_rows(kraken_txt, genus_map):
    """kunpeng classify 输出 → 病毒 contig 分类行（8 级谱系 + 属长比 + 分值）。

    score 为 kraken2 口径的可靠性分值：支持该 taxid 的 k-mer 片段数 /
    映射列总片段数（0~1）。-T 阈值即要求 score ≥ T 才判 C。
    返回 rows 列表；ratio = length/genus_avg，0.8~1.2 视为近完整基因组。
    """
    rows = []
    for flag, rid, taxid, length, mapping in _iter_kraken(kraken_txt):
        if flag != 'C' or not taxid:
            continue
        lin = lineage_ranks(taxid)
        genus = lin.get('genus', '')
        avg = genus_map.get(genus, 0)
        try:
            seqlen = int(length)
        except (TypeError, ValueError):
            seqlen = 0
        ratio = round(seqlen / avg, 3) if avg else 0
        support = total = 0
        for seg in (mapping or '').split():
            _t, _, c = seg.rpartition(':')
            try:
                cnt = int(c)
            except ValueError:
                continue
            total += cnt
            if _t == str(taxid):
                support += cnt
        score = round(support / total, 3) if total else 0
        rows.append({
            'contig': rid, 'taxid': taxid, 'taxon': lin.get('sciname', ''),
            **{r: lin.get(r, '') for r in RANKS},
            'length': seqlen, 'genus_avg_len': round(avg, 1) if avg else 0,
            'ratio': ratio,
            'near_complete': bool(avg and 0.8 <= ratio <= 1.2),
            'score': score,
            'kmer_support': support, 'kmer_total': total,
        })
    return rows


def _iter_kraken(txt):
    with safe_open(txt) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4 or parts[0] != 'C':
                continue
            try:
                tx = int(parts[2]) if parts[2].isdigit() else 0
            except ValueError:
                tx = 0
            yield 'C', parts[1], tx, parts[3], (parts[4] if len(parts) > 4 else '')


# ------------------------------------------------------------------
# 6-frame 翻译 + 最长 ORF（CDD 查询用）
# ------------------------------------------------------------------
_CODON = {}


def _codon_table():
    if _CODON:
        return _CODON
    bases = 'TCAG'
    aas = ('FFLLSSSSYY**CC*W' 'LLLLPPPPHHQQRRRR'
           'IIIMTTTTNNKKSSRR' 'VVVVAAAADDEEGGGG')
    i = 0
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                _CODON[b1 + b2 + b3] = aas[i]
                i += 1
    return _CODON


def _revcomp(dna):
    return dna.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def longest_orf_protein(dna, min_aa=50):
    """6-frame 翻译取最长 ORF 蛋白（≥min_aa），无则返回 ''。"""
    table = _codon_table()
    dna = re.sub(r'[^ACGTN]', 'N', dna.upper())
    best = ''
    for strand in (dna, _revcomp(dna)):
        for frame in range(3):
            sub = strand[frame:]
            prot = ''.join(table.get(sub[i:i + 3], 'X')
                           for i in range(0, len(sub) - 2, 3))
            for seg in prot.split('*'):
                seg = seg.strip('X')
                idx = seg.find('M')
                cand = seg[idx:] if idx >= 0 else seg
                if len(cand) >= min_aa and len(cand) > len(best):
                    best = cand
    return best


# ------------------------------------------------------------------
# NCBI BLAST (URL API) / CDD (Batch CD-Search)
# ------------------------------------------------------------------
def submit_blast(program, database, query_fasta, entrez_query=''):
    """提交 BLAST → (rid, rtoe 秒)。program: blastn/blastx。"""
    params = {
        'CMD': 'Put', 'PROGRAM': program, 'DATABASE': database,
        'QUERY': query_fasta, 'FORMAT_TYPE': 'JSON2_S',
        'HITLIST_SIZE': str(BLAST_HITLIST),
        'ALIGNMENTS': str(BLAST_HITLIST),
        'DESCRIPTIONS': str(BLAST_HITLIST),
    }
    if entrez_query:
        params['ENTREZ_QUERY'] = entrez_query
    if program == 'blastn':
        params['MEGABLAST'] = 'on'
    data = urllib.parse.urlencode(params).encode('ascii')
    with open_checked(NCBI_BLAST_URL, data=data) as resp:
        body = resp.read().decode('utf-8', errors='replace')
    rid = re.search(r'QBlastInfoBegin\s+RID\s*=\s*(\S+)', body)
    rtoe = re.search(r'QBlastInfoBegin\s+RTOE\s*=\s*(\d+)', body)
    if not rid:
        err = re.search(r'QBlastInfoBegin\s+Error\s*=\s*(.*)', body)
        raise RuntimeError(err.group(1) if err
                           else f'BLAST 提交失败: {body[:300]}')
    return rid.group(1), int(rtoe.group(1)) if rtoe else 30


def poll_blast(rid, max_wait=BLAST_TIMEOUT, sleep=None, cancel=None):
    """轮询 BLAST → 解析 hits 列表；cancel 为 threading.Event 可中断。"""
    sleep = sleep or NCBI_BLAST_SLEEP
    deadline = time.time() + max_wait
    while time.time() < deadline:
        if cancel is not None and cancel.is_set():
            raise RuntimeError('已取消')
        url = NCBI_BLAST_URL + '?' + urllib.parse.urlencode(
            {'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'JSON2_S'})
        with open_checked(url) as resp:
            body = resp.read().decode('utf-8', errors='replace')
        try:
            return parse_blast_json(json.loads(body))
        except (json.JSONDecodeError, ValueError):
            pass
        if 'Status=WAITING' in body or 'Status=UNKNOWN' in body:
            time.sleep(sleep)
            continue
        if 'Status=FAILED' in body:
            raise RuntimeError('NCBI BLAST 任务失败')
        if 'Status=READY' in body:
            time.sleep(sleep)
            continue
        if 'QBlastInfo' in body:
            err = re.search(r'Error\s*=\s*(.*)', body)
            raise RuntimeError(f'BLAST 错误: {err.group(1) if err else body[:200]}')
        time.sleep(sleep)
    raise TimeoutError(f'BLAST 轮询超时（{max_wait}s）')


def parse_blast_json(result):
    """BlastOutput2 JSON → hits 列表。"""
    hits = []
    for output in result.get('BlastOutput2', []) or []:
        search = (((output.get('report') or {}).get('results') or {})
                  .get('search') or {})
        for hit in search.get('hits', []) or []:
            desc = (hit.get('description') or [{}])[0]
            hsps = hit.get('hsps') or []
            best = hsps[0] if hsps else {}
            aln = best.get('align_len', 1) or 1
            hits.append({
                'query': search.get('query_title', 'query'),
                'accession': desc.get('accession', ''),
                'title': desc.get('title', ''),
                'sciname': desc.get('sciname', ''),
                'identity': round((best.get('identity', 0) or 0)
                                  / aln * 100, 1) if aln else 0,
                'align_len': aln,
                'evalue': best.get('evalue', ''),
                'bit_score': best.get('bit_score', 0),
            })
    return hits


def submit_cdd(query_protein_fasta):
    """提交 NCBI Batch CD-Search → cdsid。"""
    form = {
        'queries': query_protein_fasta, 'smode': 'auto', 'db': 'cdd',
        'evalue': '0.01', 'maxhit': '500', 'useid1': 'true',
        'filter': 'false', 'compbasedadj': '1', 'tdata': 'hits',
        'dmode': 'rep', 'qdefl': 'true', 'cddefl': 'true',
    }
    data = urllib.parse.urlencode(form).encode('ascii')
    with open_checked(NCBI_CDD_URL, data=data, headers={
            'Content-Type': 'application/x-www-form-urlencoded',
            'User-Agent': 'VirusPlatform/1.0'}) as resp:
        text = resp.read().decode('utf-8', errors='replace')
    m = re.search(r'^#cdsid\s+(\S+)', text, re.MULTILINE)
    if not m:
        raise RuntimeError(f'CDD 未返回 Search ID: {text[:300]}')
    return m.group(1)


def poll_cdd(cdsid, max_wait=CDD_TIMEOUT, sleep=None, cancel=None):
    """轮询 CDD → hits 列表。status: 0/4 完成, 1/2/5 失败, 3 运行中。"""
    sleep = sleep or NCBI_CDD_SLEEP
    deadline = time.time() + max_wait
    form = urllib.parse.urlencode({
        'cdsid': cdsid, 'tdata': 'hits', 'dmode': 'rep',
        'qdefl': 'true', 'cddefl': 'true'}).encode('ascii')
    while time.time() < deadline:
        if cancel is not None and cancel.is_set():
            raise RuntimeError('已取消')
        with open_checked(NCBI_CDD_URL, data=form, headers={
                'Content-Type': 'application/x-www-form-urlencoded',
                'User-Agent': 'VirusPlatform/1.0'}) as resp:
            text = resp.read().decode('utf-8', errors='replace')
        m = re.search(r'^#status\s+(\S+)', text, re.MULTILINE)
        status = m.group(1) if m else '3'
        if status in ('0', '4', 'success'):
            return parse_cdd_hits(text)
        if status in ('1', '2', '5'):
            raise RuntimeError(f'CDD 任务失败: status={status}')
        time.sleep(sleep)
    raise TimeoutError(f'CDD 轮询超时（{max_wait}s）')


def parse_cdd_hits(text):
    """Batch CD-Search TSV → hits 列表。"""
    hits = []
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.startswith('Query\t') and 'Hit type' in line:
            start = i
            break
    if start is None:
        return hits
    headers = [h.strip() for h in lines[start].split('\t')]
    for line in lines[start + 1:]:
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.split('\t')
        row = {h: (parts[j].strip() if j < len(parts) else '')
               for j, h in enumerate(headers)}
        hits.append({
            'query': row.get('Query', ''),
            'hit_type': row.get('Hit type', ''),
            'accession': row.get('Accession', ''),
            'short_name': row.get('Short name', ''),
            'evalue': row.get('E-Value', ''),
            'bit_score': row.get('Bit score', ''),
            'from': row.get('From', ''),
            'to': row.get('To', ''),
            'superfamily': row.get('Superfamily', ''),
        })
    return hits

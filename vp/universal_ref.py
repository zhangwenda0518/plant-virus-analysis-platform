# -*- coding: utf-8 -*-
"""通用病毒参考库（ref-virus / RVDB）建库支持。

数据源（databases/virus_ref/，支持多个候选文件名，自动匹配）：
- RefSeq Viral FASTA（viral*.fna.gz / *.genomic.fna.gz / *.fasta.gz / 任意 .fa/.fasta）
- RVDB FASTA      （C-RVDB*.fasta.gz / rvdb*.fasta.gz / 任意 .fa/.fasta）
- RVDB_Taxon_Current.tab.gz  accession → taxid（可选，缺失时回退 NCBI 官方映射）

也支持用户自备同名/近名文件，序列头格式（RefSeq 标准头 / RVDB acc|DB|ACC|desc /
ref| / gi| / 大小写混合）均由 _extract_accession 统一解析。

产出：info.tsv（accession→taxid 两列）+ 打好 kraken:taxid 标签的 FASTA，
走现有 vp.kunpeng.build_db 建独立库（databases/virus/ref / rvdb）。
"""
import os
import gzip
import re

from .config import DIRS, get_config, db_path
from .utils import check_path, safe_open

_REF_DIR = os.path.join(DIRS['databases'], 'virus_ref')
_TAXON_GZ = os.path.join(_REF_DIR, 'RVDB_Taxon_Current.tab.gz')
_UNI_DIR = os.path.join(_REF_DIR, 'universal')
_NCBI_ACC2TAX = os.path.join(_UNI_DIR, 'nucl_gb.accession2taxid.gz')
_NCBI_FTP = ('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/'
             'nucl_gb.accession2taxid.gz')

# 数据库名段（管道格式中需跳过的非 accession 段，小写比较）
_DB_TOKENS = frozenset((
    'ref', 'gb', 'emb', 'dbj', 'sp', 'tr', 'pir', 'prf', 'pdb',
    'tpg', 'tpe', 'tpd', 'gpp', 'gnl', 'pat', 'acc', 'genbank',
    'refseq', 'nr', 'nt', 'gi',
))


def _norm_acc(acc):
    """accession 归一化：取首个 token、去版本号、转大写。查表统一用。"""
    if not acc:
        return ''
    a = acc.strip().split()[0] if acc.strip() else ''
    a = a.split('|kraken:taxid|')[0]
    return a.split('.')[0].upper()


def _is_accession_like(s):
    """NCBI accession 模式：字母(1-6)+可选字母+数字(>=4位)(+版本号)。大小写无关。"""
    return bool(re.fullmatch(r'[A-Za-z]{1,6}_?[A-Za-z]{0,4}\d{4,}(?:\.\d+)?', s))


def extract_accession(header):
    """从 FASTA 头提取 accession（保留原大小写与版本号）。

    覆盖：RefSeq 标准头、RVDB `acc|GenBank|ACC|desc`、`ref|ACC|`、
    `gi|N|ref|ACC|`、已注入 `ACC|kraken:taxid|N`、大小写混合、Tab/CRLF。
    无法识别时退化为首个非数字 token；空头返回 None。
    """
    if not header:
        return None
    h = header.strip()
    if not h:
        return None
    # 去已注入标记，避免重复注入时污染
    h = h.split('|kraken:taxid|')[0]

    if '|' in h:
        parts = [p.strip() for p in h.split('|')]
        for p in parts:                      # 优先找像 accession 的段
            if not p or p.isdigit():
                continue
            if p.lower() in _DB_TOKENS:
                continue
            if _is_accession_like(p):
                return p
        for p in parts:                      # 退化：首个非空非数字段
            if p and not p.isdigit() and p.lower() not in _DB_TOKENS:
                return p
        return None

    tok = h.split()[0] if h.split() else h
    if _is_accession_like(tok):
        return tok
    m = re.match(r'([A-Za-z]{1,6}_?[A-Za-z]{0,4}\d{4,}(?:\.\d+)?)', h)
    return m.group(1) if m else (tok or None)


def _pick_taxid(cols):
    """从一行剩余列中挑出 taxid：首个纯数字且 >=2 的值（跳过 gi/坐标/0）。"""
    for c in cols:
        c = c.strip()
        if c.isdigit():
            v = int(c)
            if v >= 2:
                return v
    return None


def _find_ref_file(patterns):
    """在 _REF_DIR 下按通配候选名找文件，返回首个命中的绝对路径（无则 None）。

    patterns 按优先级排列，每个都是 glob 模式（小写扩展名已包含 .gz 与不压缩）。
    """
    import glob
    if not os.path.isdir(_REF_DIR):
        return None
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(_REF_DIR, pat)))
        # 排除索引/映射类文件
        hits = [p for p in hits
                if not os.path.basename(p).lower().startswith(('rvdb_taxon',
                                                               'nucl_'))]
        if hits:
            return hits[0]
    return None


def _ref_fasta(source):
    """定位源 FASTA：先精确名，再通配兜底（支持用户自备不同文件名）。"""
    exact = ('viral.1.1.genomic.fna.gz' if source == 'refvirus'
             else 'C-RVDBv32.1.fasta.gz')
    p = os.path.join(_REF_DIR, exact)
    if os.path.isfile(p):
        return p
    if source == 'refvirus':
        pats = ['viral*.genomic.fna*', 'viral*.fna*', 'viral*.fasta*',
                '*refseq*viral*.fa*', '*.fna.gz', '*.fna', '*.fasta.gz',
                '*.fasta', '*.fa.gz', '*.fa', '*.fas', '*.fsa*']
    else:
        pats = ['C-RVDB*.fasta*', 'c-rvdb*.fasta*', 'rvdb*.fasta*',
                'rvdb*.fa*', '*rvdb*.fa*', '*RVDB*.fa*', '*.fasta.gz',
                '*.fasta', '*.fa.gz', '*.fa', '*.fas', '*.fsa*']
    return _find_ref_file(pats)


def _ref_path(name):
    return check_path(os.path.join(_REF_DIR, name), must_exist=True,
                      in_platform=True)


def _fasta_accessions(source):
    """源 FASTA 内全部 accession（保留原形式，查表时统一归一化）。"""
    fa = _ref_fasta(source)
    if not fa:
        raise FileNotFoundError(
            f"未找到 {source} 源 FASTA（目录 {_REF_DIR}）。"
            f"请放入 RefSeq Viral 或 RVDB 的 FASTA 文件（.fa/.fasta，可 .gz）")
    want = set()
    opener = gzip.open if fa.endswith('.gz') else open
    with opener(fa, 'rt', errors='replace') as f:
        for line in f:
            if not line.startswith('>'):
                continue
            acc = extract_accession(line[1:])
            if acc:
                want.add(acc)
    return want


def _k2_seqid_map():
    """从 Kraken2 库包（databases/k2_viral_*.tar.gz）提取 seqid2taxid.map
    → {accession: taxid}。包不存在返回空 dict。"""
    import glob
    import tarfile
    hits = sorted(glob.glob(os.path.join(_REF_DIR, 'k2_viral_*.tar.gz')) +
                  glob.glob(os.path.join(DIRS['databases'],
                                         'k2_viral_*.tar.gz')))
    if not hits:
        return {}
    out = {}
    with tarfile.open(hits[0], 'r:*') as tf:
        try:
            fp = tf.extractfile('seqid2taxid.map')
        except KeyError:
            return {}
        for raw in fp:
            line = raw.decode().strip()
            parts = line.split('\t')
            if len(parts) == 2 and parts[1].isdigit():
                acc = parts[0].rsplit('|', 1)[-1]
                out[acc] = int(parts[1])
                out[acc.split('.')[0]] = int(parts[1])
    return out


def build_universal_meta(source='refvirus', logger=None):
    """构建 accession→taxid 映射 TSV（双源合并）。

    源1: RVDB_Taxon_Current.tab.gz（RVDB 收录集）
    源2: Kraken2 库包 seqid2taxid.map（RefSeq 官方构建集）
    并集覆盖（RefSeq 病毒实测 96.4%）。同源已存在则复用。
    """
    out_dir = check_path(os.path.join(_REF_DIR, 'universal'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    info_tsv = os.path.join(out_dir, f'{source}_acc2taxid.tsv')
    if os.path.isfile(info_tsv):
        if logger:
            logger.log(f"复用已有映射: {info_tsv}")
        return info_tsv

    want = _fasta_accessions(source)
    if logger:
        logger.log(f"FASTA 中 accession {len(want):,} 条")

    # 源0（权威）：NCBI 官方 nucl_gb.accession2taxid（覆盖全部 RefSeq/GenBank 核酸）
    ncbi = _load_ncbi_acc2tax(want, logger)

    # 源1：RVDB_Taxon（可选：用户只自备了 FASTA 时可能没有该表）
    # 列位置自适应：RVDB 各版本表头列序不统一，不假定 taxid 在第几列，
    # 而是在每行中找「首个落在合法 taxid 范围的纯数字列」（>=2，排除 gi/坐标）。
    rvdb = {}
    taxon_gz = os.path.join(_REF_DIR, 'RVDB_Taxon_Current.tab.gz')
    if os.path.isfile(taxon_gz):
        if logger:
            logger.log("解析 RVDB_Taxon 表（1100 万行，首次约 1-3 分钟）...")
        want_norm = {_norm_acc(a) for a in want}
        with gzip.open(taxon_gz, 'rt', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 2:
                    continue
                acc = parts[0].strip()
                if not acc or _norm_acc(acc) not in want_norm:
                    continue
                taxid = _pick_taxid(parts[1:])
                if taxid is not None:
                    rvdb[_norm_acc(acc)] = taxid
    elif logger:
        logger.log("未找到 RVDB_Taxon 表（可选），仅用 NCBI/Kraken2 映射源",
              "WARN")

    # 源2：Kraken2 seqid2taxid.map（补 RVDB_Taxon 未收录部分）
    k2 = _k2_seqid_map()
    if logger:
        logger.log(f"RVDB_Taxon 命中 {len(rvdb):,}；"
                   f"Kraken2 seqid2taxid 补充 {len(k2):,} 候选")

    n = ncbi_hits = 0
    with safe_open(info_tsv, 'wt') as w:
        w.write('accession\ttaxid\n')
        seen = set()
        for acc in want:
            key = _norm_acc(acc)
            t = (ncbi.get(acc) or ncbi.get(key) or
                 rvdb.get(key) or k2.get(acc) or k2.get(key))
            if t and acc not in seen:
                w.write(f'{acc}\t{t}\n')
                seen.add(acc)
                n += 1
                if acc in ncbi or key in ncbi:
                    ncbi_hits += 1
    if not n:
        raise RuntimeError(
            f'所有映射源均未命中（FASTA 内 {len(want):,} 条 accession）。'
            f'请确认源 FASTA 序列头含标准 accession，或提供 '
            f'RVDB_Taxon_Current.tab.gz / NCBI accession2taxid 映射')
    if logger:
        logger.log(f"映射完成: {n:,} 条（NCBI 官方源命中 {ncbi_hits:,}，"
                   f"双源兜底 {n - ncbi_hits:,}）→ {info_tsv}")
    return info_tsv


def build_universal_db(source='refvirus', db_dir=None, hash_capacity='2G',
                       threads=None, logger=None, rebuild=False,
                       progress=None):
    """建通用病毒参考库（refvirus / rvdb）。

    source: 'refvirus'（NCBI RefSeq Viral）或 'rvdb'（RVDB C-RVDB）。
    走 kunpeng build_db（与植物病毒库同管线），产出独立库目录。
    """
    from .kunpeng import build_db

    db_dir = db_dir or (db_path('virus', 'ref') if source == 'refvirus'
                        else db_path('virus', 'rvdb'))
    db_dir = check_path(db_dir, must_exist=False, in_platform=True)
    os.makedirs(db_dir, exist_ok=True)

    info_tsv = build_universal_meta(source, logger=logger)

    # 读映射（同时建归一化索引，兼容大小写/版本号差异）
    acc2taxid = {}
    with safe_open(info_tsv) as f:
        f.readline()
        for line in f:
            a, t = line.rstrip('\n').split('\t')
            acc2taxid[a] = int(t)
    norm2taxid = {}
    for a, t in acc2taxid.items():
        norm2taxid.setdefault(_norm_acc(a), t)

    def _lookup(acc):
        if not acc:
            return None
        return (acc2taxid.get(acc) or acc2taxid.get(acc.split('.')[0])
                or norm2taxid.get(_norm_acc(acc)))

    if logger:
        logger.log(f"映射加载: {len(acc2taxid):,} 条")

    work = os.path.join(db_dir, 'prep')
    os.makedirs(work, exist_ok=True)
    src = _ref_fasta(source)
    if not src:
        raise FileNotFoundError(
            f"未找到 {source} 源 FASTA（目录 {_REF_DIR}）。"
            f"RefSeq Viral / RVDB 的 FASTA 均可（.fa/.fasta，可 .gz）")
    if logger:
        logger.log(f"源 FASTA: {os.path.basename(src)}")
    tagged = os.path.join(work, f'{source}_tagged.fa')

    # 流式注入（源可 .gz 或不压缩）
    opener = gzip.open if src.endswith('.gz') else open
    n = skip = no_acc = 0
    skipped_samples = []
    with opener(src, 'rt', errors='replace') as f, \
            safe_open(tagged, 'wt') as w:
        for line in f:
            if line.startswith('>'):
                h = line[1:].strip()
                acc = extract_accession(h)
                if not acc:
                    no_acc += 1
                # 序列 ID 用干净的 accession（不含管道/gi/db 名段），
                # 原始头整段作为描述保留便于溯源
                sid = acc or 'seq_no_id'
                taxid = _lookup(acc)
                if taxid is None:
                    skip += 1
                    if len(skipped_samples) < 5 and acc:
                        skipped_samples.append(acc)
                    w.write(f'>{sid}|kraken:taxid|0 {h}\n')
                else:
                    n += 1
                    w.write(f'>{sid}|kraken:taxid|{taxid} {h}\n')
            else:
                w.write(line)
    if logger:
        logger.log(f"taxid 注入: {n:,} 条（跳过无 taxid {skip:,} 条"
                   + (f"，含无 accession 头 {no_acc:,} 条" if no_acc else "")
                   + "）")
        if skipped_samples:
            logger.log("跳过样本（前5例，请核对序列头格式）: "
                       + ', '.join(skipped_samples), "WARN")
    if not n:
        raise RuntimeError(
            f'没有任何序列成功注入 taxid（共跳过 {skip:,} 条）。'
            f'请确认源 FASTA 序列头含标准 accession，或提供映射表')

    # build_db 无 progress 参数（add-library/build-db 两阶段由日志体现）
    return build_db(db_dir, tagged, hash_capacity=hash_capacity,
                    threads=threads, logger=logger, rebuild=rebuild)


def _load_ncbi_acc2tax(want, logger=None):
    """NCBI 官方 nucl_gb.accession2taxid → {accession: taxid}（仅 want 内）。

    本地无该文件时自动从 NCBI FTP 下载（2.7GB gz，一次性，之后复用缓存）。
    下载失败返回空 dict（回退双源）。
    """
    if not os.path.isdir(_UNI_DIR):
        os.makedirs(_UNI_DIR, exist_ok=True)

    if not os.path.isfile(_NCBI_ACC2TAX):
        if logger:
            logger.log(f"下载 NCBI 官方映射（{_NCBI_FTP}，2.7GB，一次性）...")
        import urllib.request
        try:
            req = urllib.request.Request(_NCBI_FTP)
            with urllib.request.urlopen(req, timeout=60) as resp,                     open(_NCBI_ACC2TAX + '.part', 'wb') as out:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            os.replace(_NCBI_ACC2TAX + '.part', _NCBI_ACC2TAX)
        except Exception as e:
            if logger:
                logger.log(f"NCBI 映射下载失败（回退双源）: {e}", "WARN")
            return {}
    if logger:
        logger.log("解析 NCBI nucl_gb.accession2taxid（全量约 4 亿行，"
                   "流式过滤）...")
    out = {}
    want_base = {a.split('.')[0] for a in want}
    with gzip.open(_NCBI_ACC2TAX, 'rt', errors='replace') as f:
        f.readline()                      # header: accession	accession.version	taxid	gi
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            acc = parts[1].strip()        # 带版本号（如 NC_003214.2）
            if acc in want:
                out[acc] = int(parts[2])
            elif parts[0].strip() in want_base:
                out[parts[0].strip()] = int(parts[2])
    if logger:
        logger.log(f"NCBI 官方映射命中 {len(out):,} 条")
    return out


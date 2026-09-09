# -*- coding: utf-8 -*-
"""
GenBank 参考集合（同属共线性比较的输入层）。

与 ncbi_download（进化树参考集合，只留 FASTA+元数据）的关键差异：
本模块把 NCBI GenBank **平文原样落盘**（每条记录一个 .gb 文件），
保留 CDS/mat_peptide 等特征注释，供 vp/synteny 做基因级同属比较直接消费。

三种集合来源：
- Entrez 检索式：esearch 会话 + efetch 翻页（复用 ncbi_download 机制）
- Accession 列表：efetch 按 id 直取（逗号/空格/换行分隔，支持带版本号）
- 本机 .gb/.gbk/.gbff/.genbank 导入（含 .gz，多记录自动拆成单记录文件）

集合目录：databases/misc/gb/<name>/
  <acc>.gb          单条 GenBank 记录（比较输入，文件名即 accession）
  collection.tsv    清单（accession/物种/长度/CDS 数/宿主等，自动重建）
  query.json        来源参数（检索式或 accession 列表，供续传与溯源）

安全模型与 ncbi_download 一致：域名白名单 + IP 边界校验 + 限速重试；
文件读写经 check_path 限制在平台目录内。
"""
import io
import os
import re
import json
import time

from .config import DIRS, get_config, db_path
from .utils import check_path, safe_open
from .ncbi_download import (_open_eutils, esearch, _fetch_batch_gb,
                            EUTILS_BASE, FETCH_BATCH)

ACC_BATCH = 100               # efetch 按 accession 直取每批条数
_GB_EXTS = ('.gb', '.gbk', '.gbff', '.genbank')
_NAME_RE = re.compile(r'[A-Za-z0-9_\-\u4e00-\u9fff]{1,64}')
_ACC_RE = re.compile(r'^[A-Za-z0-9_.]{4,32}$')
_META_FIELDS = ['file', 'acc', 'title', 'organism', 'taxid', 'length',
                'molecule', 'cds_count', 'mat_peptides', 'has_translation',
                'host', 'isolate', 'country', 'date']


# ------------------------------------------------------------------
# 集合目录
# ------------------------------------------------------------------
def _gb_root():
    d = check_path(db_path('misc', 'gb'),
                   in_platform=True)
    os.makedirs(d, exist_ok=True)
    return d


def gb_collection_dir(name):
    """集合目录（名称白名单校验，防路径注入）。"""
    if not _NAME_RE.fullmatch(str(name or '')):
        raise ValueError(f'集合名仅允许中英文/数字/_/-（≤64字符）: {name!r}')
    return check_path(os.path.join(_gb_root(), name), in_platform=True)


def list_gb_collections():
    """已下载集合列表 [{name, n_records, source, date, dir}]。"""
    out = []
    root = _gb_root()
    for name in sorted(os.listdir(root)):
        cdir = os.path.join(root, name)
        if not os.path.isdir(cdir):
            continue
        n = sum(1 for f in os.listdir(cdir) if f.lower().endswith(_GB_EXTS))
        info = {}
        qfile = os.path.join(cdir, 'query.json')
        if os.path.isfile(qfile):
            try:
                with safe_open(qfile) as f:
                    info = json.load(f)
            except (OSError, ValueError):
                pass
        out.append({'name': name, 'n_records': n,
                    'source': info.get('source', ''),
                    'term': info.get('term', ''),
                    'accessions': len(info.get('accessions', [])),
                    'date': info.get('date', ''), 'dir': cdir})
    return out


# ------------------------------------------------------------------
# GenBank 平文拆分与解析
# ------------------------------------------------------------------
def split_flatfile(gb_text):
    """多记录 GenBank 平文 → 单记录文本列表（补回 // 终止符，供 SeqIO 逐条读）。"""
    chunks = []
    for c in re.split(r'(?m)^//\s*$', gb_text or ''):
        c = c.strip('\ufeff\r\n')
        if 'LOCUS' in c[:4096] and any(l.startswith('FEATURES') for l in c.splitlines()[:80]):
            chunks.append(c + '\n//')
    return chunks


def _record_summary(rec):
    """Bio.SeqIO 记录 → 元数据 dict（含 CDS/mat_peptide 统计）。"""
    src, taxid = {}, ''
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
    cds = [f for f in rec.features
           if f.type == 'CDS' and not f.qualifiers.get('pseudo')]
    mat = [f for f in rec.features if f.type == 'mat_peptide']
    return {
        'acc': rec.id,
        'title': rec.description,
        'organism': rec.annotations.get('organism', ''),
        'taxid': taxid,
        'length': len(rec.seq),
        'molecule': rec.annotations.get('molecule_type', ''),
        'cds_count': len(cds),
        'mat_peptides': len(mat),
        'has_translation': any(f.qualifiers.get('translation') for f in cds),
        'host': src.get('host', ''),
        'isolate': src.get('isolate', ''),
        'country': src.get('country', ''),
        'date': rec.annotations.get('date', ''),
    }


def parse_flatfile(gb_text):
    """平文 → [(单记录文本, 元数据 dict)]，两列表按序一一对应。"""
    chunks = split_flatfile(gb_text)
    from Bio import SeqIO
    recs = list(SeqIO.parse(io.StringIO('\n'.join(chunks)), 'genbank'))
    if len(recs) != len(chunks):          # 拆分与解析条数不一致时逐条解析兜底
        recs = []
        for c in chunks:
            r = next(iter(SeqIO.parse(io.StringIO(c), 'genbank')), None)
            if r is not None:
                recs.append(r)
        chunks = chunks[:len(recs)]
    return list(zip(chunks, [_record_summary(r) for r in recs]))


def _safe_acc_name(acc, idx):
    acc = re.sub(r'[^A-Za-z0-9_.\-]', '_', acc or '').rstrip('.')
    return f'{acc}.gb' if _ACC_RE.match(acc or 'x') else f'record_{idx + 1}.gb'


def _norm_acc(a):
    """accession 规范化（去空白/前缀，统一大写，保留版本号）。"""
    a = re.sub(r'(?i)^(gb|ref|ncbi)[:|]', '', str(a or '').strip())
    return re.sub(r'[^A-Za-z0-9_.]', '', a).upper()


def _versionless(a):
    return (a or '').split('.')[0].upper()


def _write_manifest(cdir, logger=None):
    """扫描集合内全部 .gb 重建 collection.tsv（导入/续传后清单始终一致）。"""
    rows, warnings = [], []
    for fn in sorted(os.listdir(cdir)):
        if not fn.lower().endswith(_GB_EXTS):
            continue
        path = os.path.join(cdir, fn)
        try:
            with safe_open(path) as f:
                pairs = parse_flatfile(f.read())
        except (OSError, ValueError) as e:
            warnings.append(f'{fn}: 解析失败（{e}）')
            if logger:
                logger.log(f'{fn} 解析失败，已跳过: {e}', 'WARN')
            continue
        for _, m in pairs:
            row = {**m, 'file': fn}
            if m['cds_count'] == 0:
                w = f"{m['acc']}: 无 CDS 注释，基因级比较会跳过该记录"
                warnings.append(w)
                if logger:
                    logger.log(w, 'WARN')
            elif not m['has_translation']:
                w = f"{m['acc']}: CDS 缺 /translation，将按位置翻译核苷酸"
                warnings.append(w)
                if logger:
                    logger.log(w, 'WARN')
            rows.append(row)
    meta_path = os.path.join(cdir, 'collection.tsv')
    with safe_open(meta_path, 'wt') as f:
        f.write('\t'.join(_META_FIELDS) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(k, '')) for k in _META_FIELDS) + '\n')
    # 警告随清单持久化：列表页（只读概览）无需重解析即可展示
    warn_path = os.path.join(cdir, 'warnings.txt')
    with safe_open(warn_path, 'wt') as f:
        for w in warnings:
            f.write(w + '\n')
    return rows, warnings


def inspect_collection(name, logger=None):
    """集合巡检：逐记录摘要 + 警告（不改动 .gb，只重建清单与警告文件）。"""
    cdir = gb_collection_dir(name)
    if not any(f.lower().endswith(_GB_EXTS) for f in os.listdir(cdir)):
        raise FileNotFoundError(f'集合 [{name}] 内没有 GenBank 文件')
    rows, warnings = _write_manifest(cdir, logger=logger)
    return {'name': name, 'dir': cdir, 'records': rows, 'warnings': warnings}


def read_manifest(name):
    """集合快速概览：读已生成的 collection.tsv / warnings.txt（不重解析 .gb）。

    供列表页 GET 使用；清单缺失（旧集合 / 手工放入的 .gb）时回退全量巡检。
    """
    cdir = gb_collection_dir(name)
    tsv = os.path.join(cdir, 'collection.tsv')
    if not os.path.isfile(tsv):
        return inspect_collection(name)
    with safe_open(tsv) as f:
        header = f.readline().rstrip('\n').split('\t')
        records = [dict(zip(header, line.rstrip('\n').split('\t')))
                   for line in f if line.strip()]
    warn_path = os.path.join(cdir, 'warnings.txt')
    warnings = []
    if os.path.isfile(warn_path):
        with safe_open(warn_path) as f:
            warnings = [line.rstrip('\n') for line in f if line.strip()]
    return {'name': name, 'dir': cdir, 'records': records,
            'warnings': warnings}


# ------------------------------------------------------------------
# 下载（检索式 / accession 列表）
# ------------------------------------------------------------------
def _fetch_accessions_gb(accs, log):
    """accession 列表分批 efetch → [(记录文本, 元数据)], 未找到列表。"""
    pairs, missing = [], []
    for i in range(0, len(accs), ACC_BATCH):
        batch = accs[i:i + ACC_BATCH]
        log(f'efetch 直取 {len(pairs) + 1}-{len(pairs) + len(batch)} / {len(accs)}')
        text = _open_eutils(f'{EUTILS_BASE}/efetch.fcgi', {
            'db': 'nucleotide', 'id': ','.join(batch),
            'rettype': 'gb', 'retmode': 'text'})
        got = parse_flatfile(text)
        seen = {_versionless(m['acc']) for _, m in got}
        pairs.extend(got)
        for a in batch:
            if _versionless(a) not in seen:
                missing.append(a)
    return pairs, missing


def download_gb_collection(name, term=None, accessions=None, max_records=50,
                           logger=None, prog=None, cancel=None):
    """下载一个 GenBank 集合（幂等：按 accession 跳过已有文件，可续传）。

    term（Entrez 检索式）与 accessions（列表）二选一。
    返回 {name, dir, downloaded, skipped, missing, total_hits}。
    """
    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    if bool(term) == bool(accessions):
        raise ValueError('term（检索式）与 accessions（列表）必须二选一')
    cdir = gb_collection_dir(name)
    os.makedirs(cdir, exist_ok=True)
    total_hits = n_target = None

    if term:
        if prog:
            prog('search', 0.05, 'Entrez 检索中')
        session = esearch(term, db='nucleotide')
        total_hits = session['count']
        n_target = min(total_hits, max(1, min(int(max_records), 10000)))
        log(f"检索 [{term}] 命中 {total_hits} 条，下载前 {n_target} 条")
        source, qinfo = 'query', {'term': term, 'max': n_target}
    else:
        accs = []
        for a in accessions if isinstance(accessions, (list, tuple)) \
                else re.split(r'[\s,;]+', str(accessions)):
            a = _norm_acc(a)
            if a and a not in accs:
                accs.append(a)
        if not accs:
            raise ValueError('accession 列表为空')
        for a in accs:
            if not _ACC_RE.match(a):
                raise ValueError(f'accession 格式异常: {a}')
        total_hits, n_target = len(accs), len(accs)
        source, qinfo = 'accessions', {'accessions': accs}
        log(f'按 accession 直取 {n_target} 条')

    existing = {f[:-3].upper() for f in os.listdir(cdir)
                if f.lower().endswith(_GB_EXTS)}
    if existing:
        log(f'续传模式: 集合内已有 {len(existing)} 条，重复项将跳过')

    with safe_open(os.path.join(cdir, 'query.json'), 'wt') as f:
        json.dump({'source': source, **qinfo,
                   'date': time.strftime('%Y-%m-%d %H:%M:%S')},
                  f, ensure_ascii=False, indent=2)

    pairs, missing = [], []
    if term:
        retstart = 0
        while retstart < n_target:
            if cancel is not None and getattr(cancel, 'is_set', lambda: False)():
                log('用户取消下载', 'WARN')
                break
            batch_max = min(FETCH_BATCH, n_target - retstart)
            if prog:
                prog('fetch', 0.1 + 0.85 * retstart / n_target,
                     f'下载 {min(retstart + batch_max, n_target)}/{n_target}')
            text = _fetch_batch_gb(session, 'nucleotide', retstart, batch_max)
            pairs.extend(parse_flatfile(text))
            retstart += batch_max
    else:
        if prog:
            prog('fetch', 0.15, f'efetch 直取 {n_target} 条')
        pairs, missing = _fetch_accessions_gb(qinfo['accessions'], log)
        if missing:
            log(f'未找到（请检查 accession）: {", ".join(missing)}', 'WARN')

    downloaded = skipped = 0
    for idx, (chunk, m) in enumerate(pairs):
        fname = _safe_acc_name(m['acc'], idx)
        if fname[:-3].upper() in existing:
            skipped += 1
            continue
        with safe_open(os.path.join(cdir, fname), 'wt') as f:
            f.write(chunk)
        existing.add(fname[:-3].upper())
        downloaded += 1
    log(f'落盘 {downloaded} 条 GenBank 平文，跳过(已有/重复) {skipped} 条')

    rows, _warns = _write_manifest(cdir, logger=logger)
    if prog:
        prog('done', 1.0, f'完成: 新增 {downloaded}，跳过 {skipped}')
    log(f'集合 [{name}] 就绪: {len(rows)} 条记录（含历史），'
        f'{sum(r["cds_count"] for r in rows)} 个 CDS')
    return {'name': name, 'dir': cdir, 'downloaded': downloaded,
            'skipped': skipped, 'missing': missing,
            'total_hits': total_hits or 0}


# ------------------------------------------------------------------
# 本机导入
# ------------------------------------------------------------------
def import_local_gb(name, paths, logger=None, prog=None):
    """导入本机 GenBank 文件（.gb/.gbk/.gbff/.genbank[.gz]）。

    多记录文件自动拆成单记录文件；同名 accession 跳过（幂等）。
    返回 {name, dir, imported, skipped, warnings}。
    """
    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    if isinstance(paths, str):
        paths = [p.strip() for p in re.split(r'[,\n;]+', paths) if p.strip()]
    if not paths:
        raise ValueError('未提供要导入的文件')
    cdir = gb_collection_dir(name)
    os.makedirs(cdir, exist_ok=True)
    with safe_open(os.path.join(cdir, 'query.json'), 'wt') as f:
        json.dump({'source': 'local', 'files': [os.path.abspath(p) for p in paths],
                   'date': time.strftime('%Y-%m-%d %H:%M:%S')},
                  f, ensure_ascii=False, indent=2)

    existing = {f[:-3].upper() for f in os.listdir(cdir)
                if f.lower().endswith(_GB_EXTS)}
    imported = skipped = 0
    for i, p in enumerate(paths):
        if prog:
            prog('import', 0.05 + 0.9 * i / len(paths), os.path.basename(p))
        base = p or ''
        if not base.lower().split('.gz')[0].endswith(_GB_EXTS):
            raise ValueError(f'仅支持 GenBank 文件（{", ".join(_GB_EXTS)}）: {base}')
        path = check_path(base, must_exist=True, in_platform=True)
        from .utils import open_maybe_gzip
        with open_maybe_gzip(path) as f:
            text = f.read()
        for idx, (chunk, m) in enumerate(parse_flatfile(text)):
            fname = _safe_acc_name(m['acc'], idx)
            if fname[:-3].upper() in existing:
                skipped += 1
                continue
            with safe_open(os.path.join(cdir, fname), 'wt') as f:
                f.write(chunk)
            existing.add(fname[:-3].upper())
            imported += 1
        log(f'已导入 {os.path.basename(path)}')
    rows, _ = _write_manifest(cdir, logger=logger)
    if prog:
        prog('done', 1.0, f'完成: 导入 {imported}，跳过 {skipped}')
    log(f'集合 [{name}]: 新增 {imported} 条，跳过 {skipped} 条，'
        f'共 {len(rows)} 条记录')
    return {'name': name, 'dir': cdir, 'imported': imported,
            'skipped': skipped}


def resolve_gb_files(name=None, files=None):
    """比较输入解析：集合名 → 其下全部 .gb；files → 本机 .gb 路径列表。

    两者都给时合并（集合在前）；都为空报错。返回排序后的路径列表。
    files 可为单个路径字符串（逗号/换行分隔多个）。
    """
    out = []
    if isinstance(files, str):
        files = [p.strip() for p in re.split(r'[,\n;]+', files) if p.strip()]
    if name:
        cdir = gb_collection_dir(name)
        got = [os.path.join(cdir, f) for f in sorted(os.listdir(cdir))
               if f.lower().endswith(_GB_EXTS)]
        if not got:
            raise FileNotFoundError(f'集合 [{name}] 内没有 GenBank 文件'
                                    f'（先用 gb-dl 下载或 gb-import 导入）')
        out.extend(got)
    for p in (files or []):
        if isinstance(p, str):
            p = p.strip()
        if not p:
            continue
        if not p.lower().split('.gz')[0].endswith(_GB_EXTS):
            raise ValueError(f'仅支持 GenBank 文件: {p}')
        out.append(check_path(p, must_exist=True, in_platform=True))
    if len(out) < 2:
        raise ValueError('同属比较至少需要 2 条 GenBank 记录')
    return out


# ------------------------------------------------------------------
# 集合直接建树（MSA / 进化树 / SDT 查看器的集合侧数据源）
# ------------------------------------------------------------------
def collection_records(name):
    """集合内全部记录 → [(acc, organism, seq)]（重复 accession 保留首条）。"""
    from Bio import SeqIO
    cdir = gb_collection_dir(name)
    seen, out = set(), []
    for fn in sorted(os.listdir(cdir)):
        if not fn.lower().endswith(_GB_EXTS):
            continue
        for rec in SeqIO.parse(os.path.join(cdir, fn), 'genbank'):
            if rec.id in seen or not rec.seq:
                continue
            seen.add(rec.id)
            out.append((rec.id, rec.annotations.get('organism', '')
                        or rec.description, str(rec.seq).upper()))
    if len(out) < 2:
        raise ValueError(f'集合 [{name}] 可用记录不足 2 条，无法比对建树')
    return out


def _gene_feature_seqs(name, gene, molecule):
    """从集合 .gb 记录提取匹配基因的 CDS(nt) / PEP(aa) 序列。

    gene 关键词对 CDS 的 product / gene / note 做不区分大小写子串匹配
    （每条记录取首个命中，避免多拷贝重复）；PEP 优先用 translation
    限定符，缺则按 transl_table 翻译。返回 [(id, organism, seq, product)]。
    """
    from Bio import SeqIO
    from Bio.Seq import Seq as BioSeq
    kw = (gene or '').strip().lower()
    hits = []
    cdir = gb_collection_dir(name)
    for fn in sorted(os.listdir(cdir)):
        if not fn.lower().endswith(_GB_EXTS):
            continue
        for rec in SeqIO.parse(os.path.join(cdir, fn), 'genbank'):
            for feat in rec.features:
                if feat.type != 'CDS':
                    continue
                q = feat.qualifiers
                text = ' '.join(str(v) for k in ('product', 'gene', 'note')
                                for v in q.get(k, [])).lower()
                if kw and kw not in text:
                    continue
                product = next((str(v) for v in q.get('product', [''])), '')
                gene_tag = next((str(v) for v in q.get('gene', [''])), '')
                if molecule == 'pep':
                    tr = next((str(v) for v in q.get('translation', [])), None)
                    if tr:
                        seq = tr.upper()
                    else:
                        table = next((str(v) for v in q.get('transl_table',
                                                           ['11'])), '11')
                        nt = str(feat.extract(rec.seq)).upper()
                        seq = str(BioSeq(nt).translate(table=table,
                                                       cds=False)).rstrip('*')
                    if len(seq) < 30:
                        continue
                else:                                  # cds
                    seq = str(feat.extract(rec.seq)).upper()
                    if len(seq) < 90:
                        continue
                tag = gene_tag or (product.split()[0].lower()
                                   if product else 'gene')
                hits.append((f"{rec.id}|{tag}", rec.annotations.get(
                    'organism', '') or rec.description, seq, product))
                break                                 # 每记录取首个命中
    return hits


def extract_dir(name):
    """集合特征提取产物目录 extract/。"""
    return check_path(os.path.join(gb_collection_dir(name), 'extract'),
                      must_exist=False, in_platform=True)


def _gene_tag_of(q, product):
    tag = next((str(v) for v in q.get('gene', [])), '').strip()
    if not tag and product:
        tag = product.split()[0]
    return re.sub(r'[^A-Za-z0-9_\-.]+', '_', tag)[:40] or 'gene'


def extract_collection_features(name, logger=None, prog=None):
    """集合 .gb → 分类分目录特征提取（PhyloSuite 布局）。

    extract/
      genome.fa            全基因组核酸
      CDS.fa / PEP.fa      全部 CDS 核酸 / 蛋白（汇总）
      CDS/<基因>.fa        按基因拆分的核酸（同源基因各基因组一份）
      PEP/<基因>.fa        按基因拆分的蛋白
      genes.tsv            基因 × 基因组矩阵摘要（条数/长度/产物示例）
    序列 id = <accession>|<基因>；PEP 优先 translation，缺则 transl_table 翻译。
    返回 summary dict（文件清单 + 基因统计）。
    """
    from Bio import SeqIO
    from Bio.Seq import Seq as BioSeq
    from .utils import write_fasta_record

    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    cdir = gb_collection_dir(name)
    out = extract_dir(name)
    cds_dir = os.path.join(out, 'CDS')
    pep_dir = os.path.join(out, 'PEP')
    for d in (out, cds_dir, pep_dir):
        os.makedirs(d, exist_ok=True)

    files = {'genome': os.path.join(out, 'genome.fa'),
             'cds_all': os.path.join(out, 'CDS.fa'),
             'pep_all': os.path.join(out, 'PEP.fa'),
             'genes_tsv': os.path.join(out, 'genes.tsv')}
    n_genome = n_cds = n_pep = 0
    gene_stats = {}          # tag -> {'n': int, 'lens_nt': [], 'prods': set()}
    cds_by_gene, pep_by_gene = {}, {}

    gbs = [fn_ for fn_ in sorted(os.listdir(cdir))
           if fn_.lower().endswith(_GB_EXTS)]
    if prog:
        prog('extract', 0.05, f'解析 {len(gbs)} 个 GenBank 记录')
    with safe_open(files['genome'], 'wt') as fg, \
            safe_open(files['cds_all'], 'wt') as fc, \
            safe_open(files['pep_all'], 'wt') as fp:
        for i, fn_ in enumerate(gbs, 1):
            for rec in SeqIO.parse(os.path.join(cdir, fn_), 'genbank'):
                write_fasta_record(fg, rec.id, str(rec.seq).upper())
                n_genome += 1
                for feat in rec.features:
                    if feat.type != 'CDS':
                        continue
                    q = feat.qualifiers
                    product = next((str(v) for v in q.get('product', [''])),
                                   '').strip()
                    tag = _gene_tag_of(q, product)
                    sid = f"{rec.id}|{tag}"
                    nt = str(feat.extract(rec.seq)).upper()
                    tr = next((str(v) for v in q.get('translation', [])), None)
                    if tr:
                        aa = tr.upper()
                    else:
                        table = next((str(v) for v in q.get('transl_table',
                                                           ['11'])), '11')
                        aa = str(BioSeq(nt).translate(table=table,
                                                      cds=False)).rstrip('*')
                    if len(nt) >= 90:
                        write_fasta_record(fc, sid, nt)
                        n_cds += 1
                        cds_by_gene.setdefault(tag, []).append((sid, nt))
                        st = gene_stats.setdefault(
                            tag, {'n': 0, 'lens_nt': [], 'prods': set()})
                        st['n'] += 1
                        st['lens_nt'].append(len(nt))
                        if product:
                            st['prods'].add(product[:60])
                    if len(aa) >= 30:
                        write_fasta_record(fp, sid, aa)
                        n_pep += 1
                        pep_by_gene.setdefault(tag, []).append((sid, aa))
            if prog:
                prog('extract', 0.05 + 0.75 * i / max(len(gbs), 1),
                     f'解析 {i}/{len(gbs)}')

    if prog:
        prog('write', 0.85, '按基因分目录写出')
    for tag, recs in sorted(cds_by_gene.items()):
        with safe_open(os.path.join(cds_dir, f'{tag}.fa'), 'wt') as f:
            for sid, s in recs:
                write_fasta_record(f, sid, s)
    for tag, recs in sorted(pep_by_gene.items()):
        with safe_open(os.path.join(pep_dir, f'{tag}.fa'), 'wt') as f:
            for sid, s in recs:
                write_fasta_record(f, sid, s)

    with safe_open(files['genes_tsv'], 'wt') as f:
        f.write('gene\tn_genomes\tavg_len_nt\tn_cds\tn_pep\tproducts\n')
        for tag, st in sorted(gene_stats.items()):
            prods = '; '.join(sorted(st['prods']))[:120]
            f.write(f"{tag}\t{st['n']}\t"
                    f"{sum(st['lens_nt']) // max(len(st['lens_nt']), 1)}\t"
                    f"{len(cds_by_gene.get(tag, []))}\t"
                    f"{len(pep_by_gene.get(tag, []))}\t{prods}\n")

    n_gene_files = len(set(cds_by_gene) | set(pep_by_gene))
    log(f'提取完成: 基因组 {n_genome} · CDS {n_cds} · PEP {n_pep} · '
        f'基因 {n_gene_files} 个（extract/ 分类分目录）')
    if prog:
        prog('done', 1.0, f'完成：{n_gene_files} 个基因')
    return {'collection': name, 'n_genome': n_genome, 'n_cds': n_cds,
            'n_pep': n_pep, 'n_genes': n_gene_files,
            'files': {k: os.path.basename(v) for k, v in files.items()},
            'gene_dirs': {'cds': 'CDS/', 'pep': 'PEP/'}}


def list_extract_files(name):
    """提取产物清单（未提取时返回 None）。"""
    out = extract_dir(name)
    if not os.path.isdir(out):
        return None
    items = []
    for rel in ('genome.fa', 'CDS.fa', 'PEP.fa', 'genes.tsv'):
        p_ = os.path.join(out, rel)
        if os.path.isfile(p_):
            items.append({'path': f'databases/misc/gb/{name}/extract/'
                          + rel, 'kind': rel})
    for sub in ('CDS', 'PEP'):
        d = os.path.join(out, sub)
        if os.path.isdir(d):
            for fn_ in sorted(os.listdir(d)):
                if fn_.endswith('.fa'):
                    items.append({
                        'path': f'databases/misc/gb/{name}/extract/'
                                f'{sub}/{fn_}',
                        'kind': f'{sub.lower()}_gene'})
    return items


def build_collection_phylo(name, tree_tool='fasttree', trim=True, threads=None,
                           logger=None, prog=None, cancel=None,
                           molecule='genome', gene=None):
    """集合序列 → MAFFT 比对（可选 trimAl）→ FastTree / NJ / IQ-TREE。

    molecule: 'genome'（全基因组，默认）/ 'cds'（基因 CDS 核酸）/
    'pep'（基因蛋白）——后两者需 gene 关键词（对 CDS product/gene/note
    子串匹配，如 coat protein / RdRp），输出进
    databases/misc/gb/<name>/gene_trees/<gene>_<molecule>/；
    genome 输出 phylo/。结果中心以伪样品 gb:<name>（genome）展示。
    """
    from .phylo import (_run_fasttree, _run_iqtree, _parse_iqtree_log,
                        _run_mafft, _run_nj, _run_trimal)
    from .utils import write_fasta_record

    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    def cancelled():
        return cancel is not None and getattr(cancel, 'is_set',
                                              lambda: False)()

    molecule = (molecule or 'genome').strip().lower()
    if molecule not in ('genome', 'cds', 'pep'):
        raise ValueError("molecule 仅支持 genome / cds / pep")
    if molecule in ('cds', 'pep') and not (gene or '').strip():
        raise ValueError('基因级建树需提供基因/产物关键词（如 coat protein）')

    if molecule == 'genome':
        recs = [(a, o, s) for a, o, s in collection_records(name)]
        out_dir = check_path(os.path.join(gb_collection_dir(name), 'phylo'),
                             must_exist=False, in_platform=True)
    else:
        hits = _gene_feature_seqs(name, gene, molecule)
        if len(hits) < 2:
            raise ValueError(
                f'基因 "{gene}" 在集合 [{name}] 中命中 {len(hits)} 条记录'
                '（需 ≥2 才能建树；换关键词或用全基因组模式）')
        recs = [(a, o, s) for a, o, s, _p in hits]
        slug = re.sub(r'[^A-Za-z0-9]+', '_', (gene or '').strip())[:40]
        out_dir = check_path(os.path.join(gb_collection_dir(name),
                                          'gene_trees',
                                          f'{slug}_{molecule}'),
                             must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    threads = threads or get_config().threads

    combined = os.path.join(out_dir, 'combined.fa')
    with safe_open(combined, 'wt') as f:
        for acc, org, seq in recs:
            write_fasta_record(f, acc, seq)
    unit = 'aa' if molecule == 'pep' else 'bp'
    log(f'集合 [{name}] {molecule}: {len(recs)} 条序列，'
        f'长度 {min(len(s) for _, _, s in recs)}–{max(len(s) for _, _, s in recs)} {unit}')

    if cancelled():
        raise RuntimeError('用户取消')
    if prog:
        prog('align', 0.15, 'MAFFT 全基因组比对')
    aln = os.path.join(out_dir, 'aln.fasta')
    _run_mafft(combined, aln, threads=threads, logger=logger)

    if trim:
        aln_used, trim_info = _run_trimal(
            aln, os.path.join(out_dir, 'aln.trim.fasta'), logger=logger)
    else:
        aln_used, trim_info = aln, {'applied': False}

    if cancelled():
        raise RuntimeError('用户取消')
    if prog:
        prog('tree', 0.6, '建树中' + ('（IQ-TREE 较慢）' if tree_tool == 'iqtree' else ''))
    tree_extra = {}
    if tree_tool == 'iqtree':
        tree_file = _run_iqtree(aln_used, os.path.join(out_dir, 'iqtree'),
                                threads=threads, logger=logger)
        tree_extra = _parse_iqtree_log(os.path.join(out_dir,
                                                    'iqtree.iqtree'))
    elif tree_tool == 'nj':
        tree_file = _run_nj(aln_used, os.path.join(out_dir, 'nj.nwk'),
                            logger=logger)
    else:
        if molecule == 'pep':              # 蛋白树：FastTree 蛋白模式（无 -nt）
            from .utils import run_cmd_redirect
            from .config import get_config as _cfg
            ft = _cfg().tool('fasttree')
            log('FastTree 蛋白建树（默认 LG 模型）')
            run_cmd_redirect([ft, '-gamma',
                              check_path(aln_used, must_exist=True)],
                             os.path.join(out_dir, 'tree.nwk'),
                             logger=logger)
            tree_file = check_path(os.path.join(out_dir, 'tree.nwk'),
                                   must_exist=True)
        else:
            tree_file = _run_fasttree(aln_used, os.path.join(out_dir, 'tree.nwk'),
                                      logger=logger)

    summary = {'stage': 'gb_phylo', 'collection': name,
               'n_seqs': len(recs), 'tree_tool': tree_tool,
               'molecule': molecule,
               **({'gene': gene} if molecule != 'genome' else {}),
               'trim': trim_info,
               'tree': os.path.basename(tree_file), **tree_extra,
               'records': [{'acc': a, 'organism': o, 'length': len(s)}
                           for a, o, s in recs]}
    with safe_open(os.path.join(out_dir, 'summary.json'), 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    if prog:
        prog('done', 1.0, f'完成: {len(recs)} 条序列建树')
    log(f'集合 [{name}] 建树完成 → {tree_file}')
    return summary

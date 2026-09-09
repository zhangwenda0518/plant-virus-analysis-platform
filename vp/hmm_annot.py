# -*- coding: utf-8 -*-
"""
⑥b 层2：HMM / 结构域注释（序列同源层之后的功能兜底）。

- pyhmmer 流式扫描三库（内存恒定 ~50MB，与库大小无关）：
    VOG   databases/hmm/vogdb/vog_all.hmm        49,116 profiles
    RVDB  databases/rvdb/rvdb_v32_prot.hmm       13,679 profiles
    vFam  databases/hmm/vfam/vFam-B_2014.hmm     5,585 profiles
- 命中过滤采用 rosekantor/viral_fams 口径双门槛：
    域级 i-Evalue ≤ 1e-3  且  HMM 模型覆盖率 ≥ 0.5
- 功能/分类归属查预计算元数据：
    VOG  → vog.annotations.tsv.gz（功能类别码+共识描述）
           vog.lca.tsv.gz（成员基因组谱系 LCA，尾段即科/亚科级）
    RVDB → annot/<FAM>.txt（KEYWORDS 词频 + LCA）
    vFam → annot/vFam_NNNN_annotations.txt（FAMILIES/GENERA 计数 + 成员标题）
- 结构域层：mmseqs2 搜索 NCBI Cdd（替代 RPS-BLAST，思路同 Cenote-Taker3 的
  `mmseqs databases CDD`），命中按 viral_cdds_and_pfams_191028.txt（1,580 条
  精选病毒域列表，源自 Cenote-Taker2/3）标记病毒相关性。

对外接口：
  hmm_available()            pyhmmer 是否可用
  annotate_orfs_hmm(...)     三库扫描 → {orf_id: {'hits': [...], 'best': hit}}
  cdd_search_orfs(...)       CDD 结构域搜索 → {orf_id: [hit, ...]}
  merge_hmm_cdd(...)         与层1结果合并 → 每_ORF 的 product/category/family 补全
"""
import os
import re
import gzip
import glob

from .config import DIRS, get_config, db_path

# i-Evalue / 覆盖率门槛（viral_fams 口径）
DOM_IEVAL_MAX = 1e-3
HMM_COV_MIN = 0.5
# CDD 命中过滤
CDD_EVALUE_MAX = 1e-3
CDD_MAX_HITS = 3

HMM_DIR = db_path('annot', 'hmm')
CDD_DIR = db_path('annot', 'cdd')
HMM_LIBS = [
    ('vog',  os.path.join(HMM_DIR, 'vogdb', 'vog_all.hmm')),
    ('rvdb', os.path.join(DIRS['databases'], 'rvdb', 'rvdb_v32_prot.hmm')),
    ('vfam', os.path.join(HMM_DIR, 'vfam', 'vFam-B_2014.hmm')),
]
VOG_ANNOT_TSV = os.path.join(HMM_DIR, 'vogdb', 'vog.annotations.tsv.gz')
VOG_LCA_TSV = os.path.join(HMM_DIR, 'vogdb', 'vog.lca.tsv.gz')
RVDB_ANNOT_DIR = os.path.join(DIRS['databases'], 'rvdb', 'annot')
RVDB_META_GZ = os.path.join(DIRS['databases'], 'rvdb', 'rvdb_annotations.tsv.gz')
VFAM_ANNOT_DIR = os.path.join(HMM_DIR, 'vfam', 'annot')
VIRAL_CDD_LIST = os.path.join(CDD_DIR, 'viral_cdds_and_pfams_191028.txt')
CDD_ID_TBL = os.path.join(CDD_DIR, 'cddid_all.tbl')

# VOG 功能类别码 → 平台类别（vogdb.functional_categories.txt 口径；
# 具体描述仍优先走 classify_product 关键词归类，此处为兜底）
VOG_CATEGORY_FALLBACK = {
    'Xr': '聚合酶/复制相关', 'Xs': '结构蛋白',
    'Xh': '宿主互作/致病', 'Xp': '其他功能蛋白', 'Xu': '假想蛋白（功能未定）',
}

_RVDB_META_CACHE = None
_VOG_META_CACHE = None


def hmm_available():
    """pyhmmer 可导入即视为 HMM 层可用（库文件缺失时逐库跳过）。"""
    try:
        import pyhmmer  # noqa: F401
        return True
    except ImportError:
        return False


def available_libs():
    """返回 (库名, 路径) 列表（仅存在的库）。"""
    return [(n, p) for n, p in HMM_LIBS if os.path.isfile(p)]


# ------------------------------------------------------------------
# VOG / RVDB / Pfam 元数据
# ------------------------------------------------------------------
def _strip(s):
    return (s or '').strip()


def load_vog_meta():
    """VOG id → (类别码, 共识描述, 谱系尾段)。全量载入并缓存（同 RVDB）。"""
    global _VOG_META_CACHE
    if _VOG_META_CACHE is not None:
        return _VOG_META_CACHE
    meta = {}
    try:
        with gzip.open(VOG_ANNOT_TSV, 'rt', encoding='utf-8',
                       errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                p = line.rstrip('\n').split('\t')
                if len(p) >= 5:
                    meta[p[0]] = {'cat': p[3].strip(), 'desc': p[4].strip()}
    except OSError:
        pass
    try:
        with gzip.open(VOG_LCA_TSV, 'rt', encoding='utf-8',
                       errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                p = line.rstrip('\n').split('\t')
                if len(p) >= 5 and p[0] in meta:
                    lineage = p[3].strip()
                    meta[p[0]]['lca'] = lineage
                    meta[p[0]]['family'] = lineage.split(';')[-1].strip()
    except OSError:
        pass
    _VOG_META_CACHE = meta
    return meta


_KW_STOP = {'viral', 'virus', 'viruses', 'protein', 'proteins', 'like',
            'putative', 'hypothetical', 'associated', 'non', 'none',
            'unknown', 'uncharacterized', 'containing', 'characterized'}


def _rvdb_consensus(keywords):
    """KEYWORDS 词频（已按权重降序）→ 拼共识描述：滤停用词/纯数字取前3。"""
    words = [w for w in keywords
             if w and not w.isdigit() and w.lower() not in _KW_STOP]
    return ' '.join(words[:3])


def _lineage_family_genus(lineage):
    """'Viruses::...::Closteroviridae::Ampelovirus' → (科, 属)。
    ICTV 命名码：科 *viridae、亚科 *virinae，属在其后一位。"""
    parts = [p.strip() for p in lineage.split('::') if p.strip()]
    fam_i = -1
    family = genus = ''
    for i, p in enumerate(parts):
        low = p.lower()
        if low.endswith('viridae') or low.endswith('virinae'):
            family = p
            fam_i = i
            break
    if fam_i >= 0 and fam_i + 1 < len(parts) \
            and parts[fam_i + 1].lower().endswith('virus'):
        genus = parts[fam_i + 1]
    return family, genus


def build_rvdb_meta_table(logger=None):
    """把 13,679 个 annot/<n>.txt 一次性解析成 VOG 同构预计算表
    rvdb_annotations.tsv.gz（GroupName|ProteinCount|LCA_Lineage|Family|
    Genus|FunctionalCategory|ConsensusDescription），之后元数据查询走单文件。"""
    from .orf_annot import classify_product   # 延迟导入避免环
    import glob as _glob
    files = _glob.glob(os.path.join(RVDB_ANNOT_DIR, '*.txt'))
    if not files:
        return None
    rows = []
    for i, path in enumerate(files):
        n = os.path.basename(path)[:-4]
        fam = f'FAM{int(n):06d}'
        length = nbseq = 0
        lca = ''
        kw = []
        section = None
        try:
            with open(path, encoding='utf-8', errors='replace') as f:
                for line in f:
                    line = line.rstrip('\n')
                    if line.startswith('LENGTH\t'):
                        length = line.split('\t', 1)[1].strip()
                    elif line.startswith('LCA\t'):
                        lca = line.split('\t', 1)[1].strip()
                    elif line.startswith('NBSEQ\t'):
                        nbseq = line.split('\t', 1)[1].strip()
                    elif line.startswith('KEYWORDS:'):
                        section = 'kw'
                    elif line.startswith('KEYWORDS FROM SEQUENCES'):
                        break          # 后段为原始标题词频，噪音大，不要
                    elif section == 'kw':
                        p = line.split('\t')
                        if len(p) == 2:
                            kw.append(p[0])
                        elif not line.strip():
                            section = None
        except OSError:
            continue
        lineage = lca.replace('::', ';')
        family, genus = _lineage_family_genus(lca)
        kws = _rvdb_consensus(kw)
        cat = classify_product(kws) if kws else ''
        rows.append([fam, nbseq, lineage, family, genus, cat, kws])
        if logger and (i + 1) % 2000 == 0:
            logger.log(f"  RVDB 元数据 {i + 1}/{len(files)}")
    rows.sort(key=lambda r: r[0])
    with gzip.open(RVDB_META_GZ, 'wt', encoding='utf-8') as f:
        f.write('#GroupName\tProteinCount\tLCA_Lineage\tFamily\tGenus\t'
                'FunctionalCategory\tConsensusDescription\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')
    if logger:
        logger.log(f"RVDB 元数据表生成: {len(rows)} 个家族 → {RVDB_META_GZ}")
    return RVDB_META_GZ


def load_rvdb_meta(fam_ids=None):
    """RVDB 家族元数据。优先读预计算 rvdb_annotations.tsv.gz（VOG 同构，
    一次载入并缓存）；表缺失时回退逐文件（注意 annot 文件按数字命名，
    FAM000001 ↔ 1.txt）。返回 {fam: {...}}。"""
    global _RVDB_META_CACHE
    if _RVDB_META_CACHE is not None:
        return (_RVDB_META_CACHE if fam_ids is None
                else {k: _RVDB_META_CACHE[k] for k in fam_ids
                      if k in _RVDB_META_CACHE})
    if os.path.isfile(RVDB_META_GZ):
        meta = {}
        with gzip.open(RVDB_META_GZ, 'rt', encoding='utf-8',
                       errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                p = line.rstrip('\n').split('\t')
                if len(p) >= 7:
                    meta[p[0]] = {'nbseq': p[1], 'lca': p[2].replace(';', '::'),
                                  'family': p[3], 'genus': p[4],
                                  'category': p[5], 'desc': p[6]}
        _RVDB_META_CACHE = meta
        return (meta if fam_ids is None
                else {k: meta[k] for k in fam_ids if k in meta})
    meta = {}
    for fam in (fam_ids or []):
        try:
            n = int(str(fam).replace('FAM', ''))
        except ValueError:
            continue
        d = {}
        try:
            from .orf_annot import classify_product
            with open(os.path.join(RVDB_ANNOT_DIR, f'{n}.txt'),
                      encoding='utf-8', errors='replace') as f:
                in_kw = False
                for line in f:
                    line = line.rstrip('\n')
                    if line.startswith('LCA\t'):
                        lca = line.split('\t', 1)[1].strip()
                        d['lca'] = lca
                        fam_, gen_ = _lineage_family_genus(lca)
                        d['family'], d['genus'] = fam_, gen_
                    elif line.startswith('KEYWORDS:'):
                        in_kw = True
                    elif line.startswith('KEYWORDS FROM SEQUENCES'):
                        break
                    elif in_kw:
                        p = line.split('\t')
                        if len(p) == 2:
                            d.setdefault('keywords', []).append(p[0])
                        elif not line.strip():
                            in_kw = False
        except OSError:
            pass
        kws = _rvdb_consensus(d.get('keywords', []))
        if kws:
            d['desc'] = kws
            d['category'] = classify_product(kws)
        meta[fam] = d
    return meta


def load_vfam_meta(fam_ids):
    """按需读取 vFam annot/vFam_NNNN_annotations.txt →
    {fam: {'family': 主科, 'genus': 主属, 'desc': 代表产物名}}。
    格式：FAMILIES/GENERA 计数字典 + 成员 FASTA 标题。"""
    import ast
    meta = {}
    for fam in fam_ids:
        path = os.path.join(VFAM_ANNOT_DIR, f'{fam}_annotations.txt')
        d = {}
        try:
            with open(path, encoding='utf-8', errors='replace') as f:
                titles = []
                in_titles = False
                for line in f:
                    line = line.rstrip('\n')
                    if line.startswith('FAMILIES\t'):
                        try:
                            fams = ast.literal_eval(
                                line.split('\t', 1)[1].strip())
                            d['family'] = (max(fams.items(),
                                               key=lambda x: x[1])[0]
                                           if fams else '')
                        except (ValueError, SyntaxError):
                            pass
                    elif line.startswith('GENERA\t'):
                        try:
                            gens = ast.literal_eval(
                                line.split('\t', 1)[1].strip())
                            d['genus'] = (max(gens.items(),
                                              key=lambda x: x[1])[0]
                                          if gens else '')
                        except (ValueError, SyntaxError):
                            pass
                    elif line.startswith('FASTA SEQUENCE TITLES'):
                        in_titles = True
                    elif in_titles and '|' in line:
                        titles.append(line)
                for t in titles:
                    # 标题形如 gi|..|ref|YP_xxx.1|vFam_1000| 产物 [物种]
                    body = t.split('|', 5)[-1] if t.count('|') >= 5 else t
                    prod = body.split(' [')[0].strip()
                    low = prod.lower()
                    if prod and 'hypothetical' not in low \
                            and 'unknown' not in low:
                        d['desc'] = prod
                        break
                if 'desc' not in d and titles:
                    first = titles[0].split('|', 5)[-1]
                    d['desc'] = first.split(' [')[0].strip()
        except OSError:
            pass
        meta[fam] = d
    return meta


# ------------------------------------------------------------------
# pyhmmer 扫描
# ------------------------------------------------------------------
_ASCII_LINK = None


def _ascii_path(path):
    """Easel(ANSI fopen) 打不开含中文的路径：经 %TEMP%\\vp_ascii_root 目录
    junction（无需管理员，建一次）把平台内路径转成纯 ASCII 别名。"""
    try:
        path.encode('ascii')
        return path
    except (UnicodeEncodeError, AttributeError):
        pass
    global _ASCII_LINK
    if _ASCII_LINK is None:
        import tempfile
        import subprocess
        from .config import PLATFORM_ROOT
        link = os.path.join(tempfile.gettempdir(), 'vp_ascii_root')
        if not os.path.isdir(link):
            subprocess.run(['cmd', '/c', 'mklink', '/J', link, PLATFORM_ROOT],
                           capture_output=True)
        _ASCII_LINK = link if os.path.isdir(link) else ''
    if not _ASCII_LINK:
        return path
    from .config import PLATFORM_ROOT
    try:
        rel = os.path.relpath(path, PLATFORM_ROOT)
    except ValueError:
        return path
    cand = os.path.join(_ASCII_LINK, rel)
    return cand if str(cand).isascii() else path


def _press_in_progress(lib_path):
    """hmmpress 预压进行中（h3m 已建但索引未写完）。预压期间同文件的其他
    读取会报格式错误，故跳过该库待下次运行自动启用。"""
    h3i = lib_path + '.h3i'
    return (os.path.isfile(lib_path + '.h3m')
            and (not os.path.isfile(h3i) or os.path.getsize(h3i) == 0))


def _scan_library(faa, lib_name, lib_path, threads, logger=None):
    """单库扫描（流式，内存恒定）。返回 {orf_id: [hit, ...]}，
    hit: {lib,target,acc,desc,evalue,score,cov,hmm_from,hmm_to}。"""
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    lib_path = _ascii_path(lib_path)
    faa = _ascii_path(faa)
    if _press_in_progress(lib_path):
        if logger:
            logger.log(f"  [{lib_name}] hmmpress 预压进行中，本次跳过"
                       f"（完成后自动启用，更快）", "WARN")
        return {}
    pressed = os.path.isfile(lib_path + '.h3m') and \
        os.path.getsize(lib_path + '.h3m') > 1024
    by_q = {}
    with pyhmmer.easel.SequenceFile(faa, digital=True, alphabet=alpha) as sf:
        seqs = list(sf)
    if not seqs:
        return by_q
    with pyhmmer.plan7.HMMFile(lib_path) as hf:
        targets = hf.optimized_profiles() if pressed else hf
        for hits in pyhmmer.hmmer.hmmscan(seqs, targets, cpus=threads,
                                          E=DOM_IEVAL_MAX):
            qname = hits.query.name
            qname = qname.decode() if isinstance(qname, bytes) else qname
            found = []
            for hit in hits:
                dom = hit.best_domain
                if dom is None:
                    continue
                ie = dom.i_evalue
                aln = dom.alignment
                hmm_len = aln.hmm_length or 0
                cov = ((aln.hmm_to - aln.hmm_from + 1) / hmm_len
                       if hmm_len else 0.0)
                # viral_fams 双门槛：域级 i-Evalue + 模型覆盖率
                if ie > DOM_IEVAL_MAX or cov < HMM_COV_MIN:
                    continue
                tname = hit.name
                tname = tname.decode() if isinstance(tname, bytes) else tname
                found.append({
                    'lib': lib_name, 'target': tname,
                    'acc': (getattr(hit, 'accession', '') or ''),
                    'desc': (getattr(hit, 'description', '') or ''),
                    'evalue': ie, 'score': dom.score, 'cov': round(cov, 2),
                    'hmm_from': aln.hmm_from, 'hmm_to': aln.hmm_to,
                })
            if found:
                found.sort(key=lambda h: h['evalue'])
                by_q[qname] = found[:5]
    if logger:
        n_orf = len(by_q)
        logger.log(f"  [{lib_name}] {n_orf} 个 ORF 获得合规 HMM 命中")
    return by_q


def annotate_orfs_hmm(faa, threads=None, logger=None, progress=None,
                      frac_from=0.0, frac_to=1.0):
    """三库顺序扫描。返回 ({orf_id: {'hits': [...], 'best': hit}}, [库名])。"""
    libs = available_libs()
    if not libs:
        return {}, []
    # RVDB 元数据表（VOG 同构）缺失时构建一次（约 1 分钟，之后全缓存）
    if any(n == 'rvdb' for n, _p in libs) and not os.path.isfile(RVDB_META_GZ):
        if glob.glob(os.path.join(RVDB_ANNOT_DIR, '*.txt')):
            build_rvdb_meta_table(logger=logger)
    merged = {}
    done_libs = []
    for i, (lib_name, path) in enumerate(libs):
        if progress:
            progress(frac_from + (frac_to - frac_from) * i / len(libs),
                     f'HMM 扫描 {lib_name}')
        try:
            by_q = _scan_library(faa, lib_name, path,
                                 threads or get_config().threads, logger)
        except Exception as e:
            if logger:
                logger.log(f"  [{lib_name}] 扫描失败: {e}", "WARN")
            continue
        done_libs.append(lib_name)
        for q, hits in by_q.items():
            ent = merged.setdefault(q, {'hits': [], 'best': None})
            ent['hits'].extend(hits)
            cur = ent['best']
            if cur is None or hits[0]['evalue'] < cur['evalue']:
                ent['best'] = hits[0]
        if progress:
            progress(frac_from + (frac_to - frac_from) * (i + 1) / len(libs),
                     f'HMM 扫描 {lib_name} 完成')
    return merged, done_libs


# ------------------------------------------------------------------
# CDD 结构域层（mmseqs2 替代 RPS-BLAST）
# ------------------------------------------------------------------
def cdd_db_prefix():
    """mmseqs 格式 CDD 库前缀。优先级：platform.json databases.cdd 覆盖 >
    病毒子集 cdd_virus_db（自服务器拷入，搜索更快）> 全量 cdd_db。"""
    cfg = get_config()
    p = cfg.databases.get('cdd')
    if p and os.path.isfile(p):
        return p
    for name in ('cdd_virus_db', 'cdd_db'):
        prefix = os.path.join(CDD_DIR, name)
        if os.path.isfile(prefix):
            return prefix
    return None


def _load_viral_cdd_list():
    ids = set()
    if os.path.isfile(VIRAL_CDD_LIST):
        with open(VIRAL_CDD_LIST, encoding='utf-8', errors='replace') as f:
            for line in f:
                v = line.strip()
                if v:
                    ids.add(v)
    return ids


_CDD_NAMES_CACHE = None


def _load_cdd_names():
    """cddid_all.tbl（可选）→ {CDD-ID: ShortName}。全量载入并缓存。"""
    global _CDD_NAMES_CACHE
    if _CDD_NAMES_CACHE is not None:
        return _CDD_NAMES_CACHE
    names = {}
    if os.path.isfile(CDD_ID_TBL):
        with open(CDD_ID_TBL, encoding='utf-8', errors='replace') as f:
            for line in f:
                p = line.rstrip('\n').split('\t')
                if len(p) >= 3:
                    # 列序: PSSM-Id  CDD-ID  ShortName  Description  Length
                    names[p[1].strip()] = p[2].strip()
    _CDD_NAMES_CACHE = names
    return names


def cdd_search_orfs(faa, out_dir, threads=None, logger=None):
    """mmseqs easy-search ORF 蛋白 vs CDD 库。返回 {orf_id: [hit, ...]} 或
    None（库/工具不可用）。"""
    cfg = get_config()
    try:
        mmseqs = cfg.tool('mmseqs')
    except FileNotFoundError:
        if logger:
            logger.log("mmseqs2 未安装，跳过 CDD 结构域层", "WARN")
        return None
    db = cdd_db_prefix()
    if not db:
        if logger:
            logger.log("未找到 mmseqs 格式 CDD 库（databases/cdd/cdd_db*），"
                       "跳过结构域层", "WARN")
        return None
    if logger:
        logger.log(f"  [cdd] 使用库: {os.path.basename(db)} "
                   f"({os.path.getsize(db) / 1e6:.0f}MB)")
    env = os.environ.copy()
    env['PATH'] = os.path.dirname(mmseqs) + os.pathsep + env.get('PATH', '')
    from .utils import run_cmd
    out_tsv = os.path.join(out_dir, 'cdd_hits.raw.tsv')
    tmp = out_tsv + '.mmseqs_tmp'
    run_cmd([mmseqs, 'easy-search', faa, db, out_tsv, tmp,
             '--format-output', 'query,target,evalue,pident,qlen,qstart,'
                                'qend,tstart,tend,bits',
             '-e', str(CDD_EVALUE_MAX), '--max-seqs', '5',
             '--threads', str(threads or cfg.threads)],
            logger=logger, env=env)
    viral_ids = _load_viral_cdd_list()
    names = _load_cdd_names()
    by_q = {}
    with open(out_tsv, encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 10:
                continue
            q, acc = p[0].split()[0], p[1].split()[0]
            try:
                ev = float(p[2])
                qlen, qs, qe = int(p[4]), int(p[5]), int(p[6])
            except ValueError:
                continue
            if ev > CDD_EVALUE_MAX:
                continue
            qcov = round((qe - qs + 1) / max(qlen, 1), 2)
            by_q.setdefault(q, []).append({
                'lib': 'cdd', 'target': acc,
                'desc': names.get(acc, ''),
                'evalue': ev, 'cov': qcov,
                'viral_list': acc in viral_ids,
            })
    for q in by_q:
        by_q[q] = sorted(by_q[q], key=lambda h: h['evalue'])[:CDD_MAX_HITS]
    if logger:
        logger.log(f"  [cdd] {len(by_q)} 个 ORF 获得结构域命中")
    return by_q


# ------------------------------------------------------------------
# 与层1合并
# ------------------------------------------------------------------
def format_hit_cell(hits):
    """命中列表 → 紧凑单元格文本：VOG00001(2e-40,cov0.93)|PF00680(...)。
    CDD 命中代号后附 ShortName 便于识读（如 *pfam00946 Mononeg_RNA_pol(...)）。"""
    parts = []
    for h in hits:
        mark = '*' if h.get('viral_list') or h.get('lib') == 'vog' else ''
        name = ''
        if h.get('lib') == 'cdd' and h.get('desc'):
            name = f" {h['desc']}"
        parts.append(f"{mark}{h['target']}{name}({h['evalue']:.0e},cov{h['cov']})")
    return '|'.join(parts[:5])


def merge_hmm_cdd(rows, hmm_by_q, cdd_by_q, classify_product):
    """把层2结果并回层1行。返回 (n_hmm_only, n_cdd_only)。

    优先级：层1 序列命中（informative）> HMM（VOG/RVDB/Pfam 元数据）> CDD 名称。
    rows 为 orf_annotation 行 dict（就地更新 product/organism/family/category/
    informative/evidence，并新增 hmm_hits/cdd_hits 列数据）。"""
    n_hmm_only = n_cdd_only = 0
    for r in rows:
        q = r['orf_id']
        hh = (hmm_by_q or {}).get(q)
        cc = (cdd_by_q or {}).get(q)
        r['hmm_hits'] = format_hit_cell(hh['hits']) if hh else ''
        r['cdd_hits'] = format_hit_cell(cc) if cc else ''
        if r.get('informative') == 'Y':
            r['evidence'] = 'seq'
            continue
        if hh:
            best = hh['best']
            lib, target = best['lib'], best['target']
            desc, family, cat = '', '', ''
            if lib == 'vog':
                meta = load_vog_meta().get(target, {})
                desc = meta.get('desc', '')
                family = meta.get('family', '')
                cat = classify_product(desc)
                if cat == '假想蛋白（功能未定）' or not desc:
                    cat = VOG_CATEGORY_FALLBACK.get(meta.get('cat', ''),
                                                    cat) or cat
            elif lib == 'rvdb':
                meta = load_rvdb_meta([target]).get(target, {})
                kws = meta.get('desc', '')
                desc = f"RVDB {target}" + (f" ({kws})" if kws else '')
                family = meta.get('family', '')
                genus = meta.get('genus', '')
                if genus and not r.get('genus'):
                    r['genus'] = genus
                cat = meta.get('category') or (
                    classify_product(kws) if kws else '其他功能蛋白')
            else:  # vfam
                meta = load_vfam_meta([target]).get(target, {})
                desc = meta.get('desc') or f"vFam {target.split('_')[-1]}"
                family = meta.get('family', '')
                cat = classify_product(desc)
            r['product'] = desc or target
            r['category'] = cat
            if not r.get('family'):
                r['family'] = family
            r['informative'] = 'Y'
            r['evidence'] = 'hmm'
            n_hmm_only += 1
            continue
        if cc:
            best = cc[0]
            name = best.get('desc') or best['target']
            r['product'] = f"CDD {best['target']}: {name}" if best.get('desc') \
                else f"CDD {best['target']}"
            r['category'] = classify_product(name)
            r['informative'] = 'Y'
            r['evidence'] = 'cdd'
            n_cdd_only += 1
            continue
        r['evidence'] = 'none'
    return n_hmm_only, n_cdd_only

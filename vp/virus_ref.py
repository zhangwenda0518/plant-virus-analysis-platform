# -*- coding: utf-8 -*-
"""
非冗余病毒参考库（databases/virus_ref/）——元数据中心。

数据件套（源自 plant_virus_db_pipeline，更新时整体替换）:
  Plant_Virus_Ref.fasta / Plant_Virus_Ref.Info.tsv      8,465 条非冗余代表集（98% ANI 聚类）
  Plant_Virus.complete_ref.fasta / _info.tsv            5,773 条完整基因组子集（BLAST 确认用）
  Plant_Virus_Full.fasta / Plant_Virus_Full.Info.tsv    ~199K 全量序列（高灵敏深扫，可选）
  DATA_VERSION                                          版本锚点（VERSION / LAST_INCREMENTAL / ...）

原始 Info 表的谱系列覆盖不全（VMR_Family/Genus 仅 ~55%，VMR_Species 仅 ~10%，
Segment 写法混乱如 "DNA A"/"DNA-A"）。本模块负责:
  1) Segment 归一化（norm_segment）；
  2) 用本地 taxonomy（nodes/names.dmp）沿 Taxid 走谱系，补全
     Order/Family/Genus/Species → 生成规范化元数据缓存 ref_meta.tsv；
  3) 对外提供 ③BLAST 库取数（complete_ref）、④宿主交叉验证
     （物种级已知宿主集）、⑦参考挑选（RefSeq/complete 优先）所需的查询接口。
"""
import os
import csv
import glob

from .config import DIRS
from .utils import check_path, safe_open

_REF_DIR = os.path.join(DIRS['databases'], 'virus_ref')

# ref_meta.tsv 输出列（缓存文件，可用 build_meta(force=True) 重建）
META_FIELDS = [
    'Accession', 'Taxid', 'Species_NCBI', 'Species_ICTV', 'Species',
    'Segment', 'Segment_Norm', 'Sequence_Type', 'Category', 'Topology',
    'Molecule_type', 'Length', 'Nuc_Completeness', 'Host',
    'Isolation_Source', 'Geo_Location',
    'VMR_Species', 'VMR_Genus', 'VMR_Family',
    'Lin_Order', 'Lin_Family', 'Lin_Genus', 'Lin_Species',
    'In_Complete_Ref',
]


def ref_dir():
    """virus_ref 目录（不存在返回 None，调用方回退旧数据源）。"""
    if os.path.isdir(_REF_DIR):
        return _REF_DIR
    return None


def _ref_file(name):
    d = ref_dir()
    if not d:
        return None
    p = os.path.join(d, name)
    return p if os.path.isfile(p) else None


def ref_fasta():
    """非冗余代表集 FASTA（8,465 条）。"""
    return _ref_file('Plant_Virus_Ref.fasta')


def complete_fasta():
    """完整基因组子集 FASTA（5,773 条；③BLAST 确认库的建库源）。"""
    return _ref_file('Plant_Virus.complete_ref.fasta')


def full_fasta():
    """全量序列 FASTA（~199K 条；仅高灵敏深扫用）。"""
    return _ref_file('Plant_Virus_Full.fasta')


def ref_info_tsv():
    return _ref_file('Plant_Virus_Ref.Info.tsv')


def available():
    return ref_fasta() is not None and ref_info_tsv() is not None


def read_version():
    """DATA_VERSION → {KEY: VALUE}（# 注释行跳过）。"""
    p = _ref_file('DATA_VERSION')
    out = {}
    if not p:
        return out
    with safe_open(p) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or '=' not in line:
                continue
            k, v = line.split('=', 1)
            out[k.strip()] = v.strip()
    return out


def version_stamp():
    """缓存失效锚点：版本关键字段 + 两个 FASTA 大小。"""
    v = read_version()
    parts = [v.get('VERSION', ''), v.get('LAST_INCREMENTAL', '')]
    for p in (ref_fasta(), complete_fasta()):
        try:
            parts.append(str(os.path.getsize(p)) if p else '-')
        except OSError:
            parts.append('-')
    return '|'.join(parts)


# ------------------------------------------------------------------
# Segment 归一化
# ------------------------------------------------------------------
def norm_segment(s):
    """'DNA A'/'dna a'/'DNA-A.' → 'DNA-A'；'rna1' → 'RNA-1'；空 → ''。

    规则：大写、空白与下划线转连字符、去首尾标点、连续连字符合并。
    数字与字母之间无连字符时补连字符（RNA1→RNA-1, DNAB→DNA-B）。
    """
    s = (s or '').strip().upper().replace('_', ' ')
    s = s.replace('-', ' ')
    s = '-'.join(p for p in s.split(' ') if p)
    s = s.strip('.').rstrip('-')
    # RNA1 / DNA2 / SATRNA 等字母数字边界补连字符
    for prefix in ('RNA', 'DNA', 'CDNA', 'SEGMENT'):
        if s.startswith(prefix) and len(s) > len(prefix) \
                and s[len(prefix)].isdigit():
            s = prefix + '-' + s[len(prefix):]
    return s


# ------------------------------------------------------------------
# 规范化元数据缓存（谱系补全）
# ------------------------------------------------------------------
def meta_cache_path():
    d = ref_dir()
    return os.path.join(d, 'ref_meta.tsv') if d else None


def _stamp_path():
    d = ref_dir()
    return os.path.join(d, 'ref_meta.stamp') if d else None


def build_meta(force=False, logger=None):
    """读 Ref.Info.tsv + taxonomy 谱系 → 写 ref_meta.tsv 缓存。

    补全规则：
      Species      = Species_ICTV > Species_NCBI > 谱系 species
      VMR_Genus/Family = 原值 > 谱系值
      VMR_Species  = 原值 > Species（即用 taxonomy 补 10% 覆盖的缺口）
    版本锚点（ref_meta.stamp）变化或 force 时重建。
    """
    out_p = meta_cache_path()
    stamp_p = _stamp_path()
    if not out_p:
        return None
    stamp = version_stamp()
    if not force and os.path.isfile(out_p) and os.path.isfile(stamp_p):
        try:
            with safe_open(stamp_p) as f:
                if f.read().strip() == stamp:
                    return out_p
        except OSError:
            pass

    from .contig_annot import lineage_ranks

    # 完整基因组子集 accession 集
    complete_accs = set()
    cp = _ref_file('Plant_Virus.complete_ref_info.tsv')
    if cp:
        with safe_open(cp) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                a = (row.get('Accession') or '').strip()
                if a:
                    complete_accs.add(a)

    info_p = ref_info_tsv()
    if not info_p:
        return None
    n = 0
    with safe_open(out_p, 'wt') as out:
        out.write('\t'.join(META_FIELDS) + '\n')
        with safe_open(info_p) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                taxid = (row.get('Taxid') or '').strip()
                lin = lineage_ranks(taxid) if taxid.isdigit() else {}
                species = ((row.get('Species_ICTV') or '').strip()
                           or (row.get('Species_NCBI') or '').strip()
                           or lin.get('species', ''))
                rec = {
                    'Accession': (row.get('Accession') or '').strip(),
                    'Taxid': taxid,
                    'Species_NCBI': (row.get('Species_NCBI') or '').strip(),
                    'Species_ICTV': (row.get('Species_ICTV') or '').strip(),
                    'Species': species,
                    'Segment': (row.get('Segment') or '').strip(),
                    'Segment_Norm': norm_segment(row.get('Segment')),
                    'Sequence_Type': (row.get('Sequence_Type') or '').strip(),
                    'Category': (row.get('Category') or '').strip(),
                    'Topology': (row.get('Topology') or '').strip(),
                    'Molecule_type': (row.get('Molecule_type') or '').strip(),
                    'Length': (row.get('Length') or '').strip(),
                    'Nuc_Completeness': (row.get('Nuc_Completeness') or '').strip(),
                    'Host': (row.get('Host') or '').strip(),
                    'Isolation_Source': (row.get('Isolation_Source') or '').strip(),
                    'Geo_Location': (row.get('Geo_Location') or '').strip(),
                    'VMR_Species': (row.get('VMR_Species') or '').strip() or species,
                    'VMR_Genus': (row.get('VMR_Genus') or '').strip()
                                 or lin.get('genus', ''),
                    'VMR_Family': (row.get('VMR_Family') or '').strip()
                                  or lin.get('family', ''),
                    'Lin_Order': lin.get('order', ''),
                    'Lin_Family': lin.get('family', ''),
                    'Lin_Genus': lin.get('genus', ''),
                    'Lin_Species': lin.get('species', ''),
                    'In_Complete_Ref': '1' if (row.get('Accession') or '').strip()
                                       in complete_accs else '',
                }
                out.write('\t'.join(rec[k] for k in META_FIELDS) + '\n')
                n += 1
    with safe_open(stamp_p, 'wt') as f:
        f.write(stamp)
    if logger:
        logger.log(f"ref_meta 规范化元数据已生成: {n} 条 -> {out_p}")
    return out_p


def load_meta(logger=None, force=False):
    """加载规范化元数据 → [dict]；缓存缺失时自动构建。"""
    p = build_meta(logger=logger, force=force)
    if not p or not os.path.isfile(p):
        return []
    rows = []
    with safe_open(p) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row.get('Accession'):
                rows.append(row)
    return rows


def meta_by_accession(logger=None, force=False):
    """{Accession: dict}——③/④/⑦ 按 accession 并表的主索引。"""
    return {r['Accession']: r for r in load_meta(logger=logger, force=force)}


def species_host_map(logger=None, force=False):
    """{Species 小写: [已知宿主列表]}——④步物种级已知宿主交叉验证。

    物种键同时收 Species_ICTV / Species_NCBI / VMR_Species 三种写法，
    便于不同来源的物种名都能命中。
    """
    out = {}
    for r in load_meta(logger=logger, force=force):
        host = r.get('Host', '').strip()
        if not host:
            continue
        keys = {r.get('Species_ICTV', '').strip().lower(),
                r.get('Species_NCBI', '').strip().lower(),
                r.get('VMR_Species', '').strip().lower()} - {''}
        for k in keys:
            lst = out.setdefault(k, [])
            if host not in lst:
                lst.append(host)
    return out


def status():
    """库状态摘要（GUI/CLI 用）。"""
    v = read_version()
    n_ref = n_complete = 0
    for h, _s in _iter_fasta_len(ref_fasta() or ''):
        n_ref += 1
    for h, _s in _iter_fasta_len(complete_fasta() or ''):
        n_complete += 1
    return {
        'dir': ref_dir(),
        'version': v.get('VERSION', ''),
        'last_incremental': v.get('LAST_INCREMENTAL', ''),
        'source_ictv': v.get('SOURCE_ICTV', ''),
        'source_ncbi': v.get('SOURCE_NCBI', ''),
        'n_ref': n_ref,
        'n_complete': n_complete,
        'meta_ready': os.path.isfile(meta_cache_path() or ''),
    }


def _iter_fasta_len(path):
    """轻量 FASTA 扫描（仅计数/长度）。"""
    if not path or not os.path.isfile(path):
        return
    name, n = None, 0
    with safe_open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        yield name, n

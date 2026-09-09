# -*- coding: utf-8 -*-
"""
全病毒参考库（databases/acvirus_db/，来自 246 服务器 acvirus_db）。

文件:
  all_virus.fasta   ~20K 条全病毒基因组（不止植物；核酸，建树参考池）
  all_virus.faa     对应蛋白（Prodigal 预测；⑥b 注释可选用，当前未接入）
  taxa.txt          accession → ICTV 全谱系（Realm…Species，GENBANK accession 无版本号）
  taxon_min_coverage.csv  各分类级最低覆盖度阈值（分类判定用，⑦暂未使用）
  species.tsv       服务器侧工作聚类表（54 行，平台未使用）

⑦建树接入：病毒 contigs 除对比植物 complete_ref 外，再对比 all_virus 库，
两组 hits 合并挑选参考——非植物病毒（或植物库缺 representativa）也能
选到近缘完整基因组建树；参考的 ICTV 谱系由 taxa.txt 提供。
"""
import os
import csv
import glob

from .config import DIRS
from .utils import check_path, safe_open, run_cmd

_DIR = os.path.join(DIRS['databases'], 'acvirus_db')

_TAXA_FIELDS = ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Suborder',
                'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species',
                'Genbank']


def db_dir():
    return _DIR if os.path.isdir(_DIR) else None


def fasta_path():
    p = os.path.join(_DIR, 'all_virus.fasta') if db_dir() else None
    return p if p and os.path.isfile(p) else None


def faa_path():
    p = os.path.join(_DIR, 'all_virus.faa') if db_dir() else None
    return p if p and os.path.isfile(p) else None


def taxa_path():
    p = os.path.join(_DIR, 'taxa.txt') if db_dir() else None
    return p if p and os.path.isfile(p) else None


def available():
    return fasta_path() is not None and taxa_path() is not None


def _stamp():
    p = fasta_path()
    try:
        st = os.stat(p)
        return f'{st.st_size}|{int(st.st_mtime)}'
    except OSError:
        return ''


# ------------------------------------------------------------------
# BLAST 库（ASCII 路径，同 ③ 的病毒库处理；版本戳存 acvirus_db 内）
# ------------------------------------------------------------------
def ensure_blast_db(logger=None):
    """构建/复用 all_virus 核酸 BLAST 库，返回 db 前缀。"""
    from .config import get_config
    cfg = get_config()
    makeblastdb = cfg.tool('makeblastdb')
    ref = fasta_path()
    if not ref:
        raise FileNotFoundError('acvirus_db/all_virus.fasta 不存在')
    from .assembly import _ascii_work_base
    base = os.path.normpath(os.path.abspath(_ascii_work_base('vp_blast')))
    prefix = os.path.join(base, 'acvirus')
    # 前缀拼装后校验仍位于工作目录内
    if not os.path.normpath(os.path.abspath(prefix)).startswith(base + os.sep):
        raise RuntimeError(f'BLAST 库路径越界: {prefix}')
    stamp = _stamp()
    stamp_file = check_path(os.path.join(_DIR, 'blast_db.stamp'),
                            must_exist=False, in_platform=True)
    if glob.glob(prefix + '.n??'):
        cur = ''
        try:
            with safe_open(stamp_file) as f:
                cur = f.read().strip()
        except (OSError, ValueError):
            cur = ''
        if cur == stamp:
            return prefix
        for p in glob.glob(prefix + '.n??'):
            if os.path.normpath(os.path.abspath(p)).startswith(base + os.sep):
                try:
                    os.remove(p)
                except OSError:
                    pass
        if logger:
            logger.log('acvirus 库已更新，重建 BLAST 库')
    if logger:
        logger.log(f'构建 acvirus BLAST 库: {ref} -> {prefix}')
    run_cmd([makeblastdb, '-in', ref, '-dbtype', 'nucl', '-out', prefix,
             '-title', 'acvirus_all_virus'], logger=logger)
    with safe_open(stamp_file, 'wt') as f:
        f.write(stamp)
    return prefix


# ------------------------------------------------------------------
# taxa.txt 谱系
# ------------------------------------------------------------------
_taxa_cache = {}


def load_taxa(force=False):
    """taxa.txt → {无版本 accession: {rank: name, ...}}（含 Genbank）。"""
    if _taxa_cache.get('map') and not force:
        return _taxa_cache['map']
    p = taxa_path()
    out = {}
    if p:
        with safe_open(p) as f:
            for row in csv.DictReader(f):
                acc = (row.get('Virus GENBANK accession')
                       or '').strip().split('.')[0].upper()
                if not acc:
                    continue
                out[acc] = {k: (row.get(k) or '').strip()
                            for k in _TAXA_FIELDS}
    _taxa_cache['map'] = out
    return out


def lineage_for(accession):
    """accession（含/不含版本号）→ ICTV 谱系 dict；未命中返回 None。"""
    base = accession.split('.')[0].upper()
    return load_taxa().get(base)


def species_for(accession):
    lin = lineage_for(accession)
    return (lin or {}).get('Species', '')

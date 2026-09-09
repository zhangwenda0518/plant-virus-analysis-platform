#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unified_metadata.py — 统一元数据模板生成器 v1.0
================================================

借鉴 SeqSender 的 "一张表覆盖 GenBank + BioSample" 设计:
  一次填写 ~25 个核心字段, 自动同时生成:
    - source.src   (GenBank source modifiers)
    - miuvig.tsv    (MIUVIG 全局参数)
    - assembly.tsv  (Assembly 注释)
    - BioSample CSV (可选, 用于 NCBI BioSample 批量注册)
    - authorset.sbt (可选, 作者信息)

字段前缀约定 (同 SeqSender):
  src-*  → source.src 列 (如 src-Isolate, src-geo_loc_name)
  cmt-*  → .cmt 注释 (如 cmt-Assembly_Method)
  bs-*   → BioSample 字段 (如 bs-isolate, bs-geo_loc_name)
  无前缀 → 通用字段 (organism, collection_date, bioproject)

用法:
  # 从 taxonomy + 公共元数据生成统一模板
  python unified_metadata.py \\
      --taxonomy taxonomy.tsv \\
      --metadata Global_Unified_Metadata_Core14.tsv \\
      --run-title my_project \\
      -o ./submission/

输出:
  submission/
  ├── unified_metadata.csv      ← 一张表 (可 Excel 编辑)
  ├── source.src                ← 自动从 CSV 生成
  ├── source_individual/        ← 每个病毒独立 source.src
  ├── biosample_template.tsv    ← NCBI BioSample 批量模板
  ├── miuvig.tsv
  └── assembly.tsv

（移植自 MMPV-RNA virome_submission_pipeline/unified_metadata.py；
 文件写出统一改走平台 vp.utils.safe_open，默认工具链/提交指引改为本平台）
"""

import argparse
import os
import sys
import csv
import re
import logging
from pathlib import Path
from datetime import datetime
from collections import OrderedDict

import pandas as pd
from tqdm import tqdm

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from vp.utils import safe_open
else:
    from ..utils import safe_open


# ══════════════════════════════════════════════════════════════
# 统一元数据表定义 (精简自 SeqSender, 保留核心 25 字段)
# ══════════════════════════════════════════════════════════════

UNIFIED_COLUMNS = OrderedDict([
    # === 通用 (Required) ===
    ("organism",               {"required": True,  "group": "general", "desc": "NCBI Taxonomy 物种名"}),
    ("sequence_name",          {"required": True,  "group": "general", "desc": "FASTA 中的序列 ID"}),
    ("authors",                {"required": True,  "group": "general", "desc": "引用作者: Last, First; ..."}),
    ("collection_date",        {"required": True,  "group": "general", "desc": "采集日期 (YYYY-MM-DD 或 YYYY-MM 或 YYYY)"}),
    ("bioproject",             {"required": True,  "group": "general", "desc": "NCBI BioProject ID (PRJNA...)"}),

    # === source.src 字段 (src-*) ===
    ("src-Isolate",            {"required": True,  "group": "src", "desc": "唯一分离株标识 (分段病毒各片段必须相同)"}),
    ("src-geo_loc_name",       {"required": True,  "group": "src", "desc": "采集地点 (格式 Country:Region, 如 China:Ningxia)"}),
    ("src-Lat_Lon",            {"required": True,  "group": "src", "desc": "经纬度 (如 38.47 N 106.27 E)"}),
    ("src-Host",               {"required": False, "group": "src", "desc": "宿主物种名"}),
    ("src-Segment",            {"required": False, "group": "src", "desc": "分段病毒的片段编号 (非分段留空)"}),
    ("src-Isolation-source",   {"required": False, "group": "src", "desc": "分离来源描述"}),
    ("src-Note",               {"required": False, "group": "src", "desc": "额外备注"}),
    ("src-Tissue_type",        {"required": False, "group": "src", "desc": "组织类型"}),
    ("src-Collected_by",       {"required": False, "group": "src", "desc": "采集人 (Core14 CenterName)"}),
    ("src-Cultivar",           {"required": False, "group": "src", "desc": "栽培品种 (Core14 Source)"}),
    ("src-Dev_stage",          {"required": False, "group": "src", "desc": "发育阶段 (Core14 Age_GrowthStage)"}),

    # === GenBank 提交字段 ===
    ("gb-sample_name",         {"required": True,  "group": "gb", "desc": "GenBank 记录名 (≤50字符)"}),
    ("gb-title",               {"required": False, "group": "gb", "desc": "提交标题 (NCBI 门户显示)"}),
    ("sra",                    {"required": False, "group": "gb", "desc": "SRA 登录号 (SRR...)"}),
    ("biosample",              {"required": True,  "group": "gb", "desc": "BioSample ID (SAMN...)"}),

    # === 结构化注释字段 (cmt-*) ===
    ("cmt-Assembly_Method",    {"required": True,  "group": "cmt", "desc": "组装方法 (如 SPAdes v4.3.0 metaviral)"}),
    ("cmt-Sequencing_Technology", {"required": True, "group": "cmt", "desc": "测序平台 (如 Illumina NovaSeq 6000)"}),
    ("cmt-Genome_Coverage",    {"required": False, "group": "cmt", "desc": "基因组覆盖度 (如 42.5x)"}),
    ("cmt-Annotation_Pipeline", {"required": False, "group": "cmt", "desc": "注释流程 (如 植物病毒分析平台 ⑥b ORF 功能注释)"}),

    # === BioSample 字段 (bs-*, 如需自动注册) ===
    ("bs-isolate",             {"required": False, "group": "bs", "desc": "BioSample: 分离株标识"}),
    ("bs-geo_loc_name",        {"required": False, "group": "bs", "desc": "BioSample: 采集地点"}),
    ("bs-host",                {"required": False, "group": "bs", "desc": "BioSample: 宿主"}),
    ("bs-isolation_source",    {"required": False, "group": "bs", "desc": "BioSample: 分离来源"}),
])

REQUIRED_COLS = [k for k, v in UNIFIED_COLUMNS.items() if v["required"]]


# ══════════════════════════════════════════════════════════════

def setup_logger(out_dir):
    logger = logging.getLogger("UnifiedMeta")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)
    os.makedirs(out_dir, exist_ok=True)
    return logger


def sanitize_filename(name, max_len=50):
    name = re.sub(r'[^\w\s-]', '', name)
    name = re.sub(r'\s+', '_', name)
    return name[:max_len]


def load_taxonomy(tax_tsv):
    """加载 taxonomy.tsv"""
    df = pd.read_csv(tax_tsv, sep='\t')
    col_map = {}
    for c in df.columns:
        if c.lower() in ('contig', 'seq_id', 'sequence_id'):
            col_map['contig'] = c
        elif c.lower() in ('taxonomy', 'tax'):
            col_map['taxonomy'] = c
    if 'contig' not in col_map:
        col_map['contig'] = df.columns[0]
    if 'taxonomy' not in col_map:
        col_map['taxonomy'] = df.columns[1] if len(df.columns) > 1 else df.columns[0]

    seqs = []
    for _, row in df.iterrows():
        seqs.append({
            'contig': str(row[col_map['contig']]),
            'taxonomy': str(row[col_map['taxonomy']]) if pd.notna(row[col_map['taxonomy']]) else 'unclassified viruses',
        })
    return seqs


def load_metadata(meta_file, log=None):
    """加载元数据表 — 自动优先 Global_Unified_Metadata_Full.tsv。

    若传入的是 Core14.tsv, 自动在同目录查找 Full.tsv 替代。
    Full 列 (34): Run, BioSample, Platform, CollectionDate, Location, Source,
                  Tissue, Age_GrowthStage, ScientificName, TaxID, CenterName,
                  BioProject, PMID, SRAStudy, Sample, LibraryLayout, ...
    Core14 列 (14): Run, ReleaseDate, CollectionDate, Location, Source, Tissue,
                    Age_GrowthStage, ScientificName, TaxID, LibrarySource,
                    CenterName, BioProject, BioSample, PMID
    """
    if not meta_file or not os.path.exists(meta_file):
        return {}

    # ── 自动优先 Full ──
    actual_file = meta_file
    if 'Core14' in str(meta_file):
        full_candidate = str(meta_file).replace('Core14', 'Full')
        if os.path.exists(full_candidate):
            actual_file = full_candidate
            if log:
                log.info("  → 自动使用 Full 元数据: %s",
                         os.path.basename(full_candidate))

    df = pd.read_csv(actual_file, sep=None, engine='python')
    # Strip BOM from column names
    df.columns = [c.lstrip('\ufeff') for c in df.columns]

    has_biosample = 'BioSample' in df.columns
    has_platform = 'Platform' in df.columns

    lookup = {}
    for _, row in df.iterrows():
        run = str(row.get('Run', '')).strip()
        if not run or run.lower() in ('nan', 'not_provided', ''):
            continue

        # Location 可能是空的 (GSA 数据), 用 Source + CenterName 推断
        loc = str(row.get('Location', '')) if pd.notna(row.get('Location')) else ''
        if not loc or loc.lower() in ('nan', 'not_provided', ''):
            center = str(row.get('CenterName', ''))
            if 'Beijing' in center or 'China' in center:
                loc = 'China'
            elif center and center.lower() not in ('nan', 'not_provided', ''):
                loc = center

        lookup[run] = {
            'collection_date': str(row.get('CollectionDate', '')) if pd.notna(row.get('CollectionDate')) else '',
            'geo_loc_name': loc,
            'bioproject': str(row.get('BioProject', '')) if pd.notna(row.get('BioProject')) else '',
            'biosample': str(row.get('BioSample', '')) if has_biosample and pd.notna(row.get('BioSample')) else '',
            'platform': str(row.get('Platform', '')) if has_platform and pd.notna(row.get('Platform')) else '',
            'host': str(row.get('ScientificName', '')) if pd.notna(row.get('ScientificName')) else '',
            'tissue': str(row.get('Tissue', '')) if pd.notna(row.get('Tissue')) else '',
            'source': str(row.get('Source', '')) if pd.notna(row.get('Source')) else '',
            'center': str(row.get('CenterName', '')) if pd.notna(row.get('CenterName')) else '',
        }

    if log:
        n_total = len(lookup)
        n_with_date = sum(1 for v in lookup.values() if v['collection_date'] and v['collection_date'].lower() not in ('nan', 'not_provided', ''))
        n_with_bs = sum(1 for v in lookup.values() if v.get('biosample') and v['biosample'].lower() not in ('nan', 'not_provided', ''))
        n_with_tissue = sum(1 for v in lookup.values() if v['tissue'] and v['tissue'].lower() not in ('nan', 'not_provided', ''))
        log.info("  元数据: %d SRA/CRR, 含日期=%d, BioSample=%d, 组织=%d",
                 n_total, n_with_date, n_with_bs, n_with_tissue)

    return lookup


def extract_sra(contig):
    m = re.match(r'([SC]RR\d+)', contig)
    return m.group(1) if m else 'UNKNOWN'


def geo_to_latlon(geo_loc):
    """从城市/地区名推断经纬度 (约值, NCBI 接受)"""
    if not geo_loc or geo_loc.lower() in ('country:region', 'nan', 'not_provided', ''):
        return ''
    loc_lower = geo_loc.lower()
    # 中国主要城市/地区 → 经纬度
    china_cities = {
        'yinchuan':     '38.47 N 106.27 E',
        'ningxia':      '37.48 N 105.68 E',
        'zhongning':    '37.48 N 105.68 E',
        'beijing':      '39.90 N 116.40 E',
        'shanghai':     '31.23 N 121.47 E',
        'guangzhou':    '23.13 N 113.26 E',
        'wuhan':        '30.59 N 114.31 E',
        'nanjing':      '32.06 N 118.79 E',
        'hangzhou':     '30.27 N 120.15 E',
        'chengdu':      '30.57 N 104.07 E',
        'xian':         '34.26 N 108.94 E',
        'kunming':      '25.04 N 102.68 E',
        'harbin':       '45.80 N 126.53 E',
        'zhengzhou':    '34.75 N 113.62 E',
        'jinan':        '36.65 N 116.98 E',
        'taiyuan':      '37.87 N 112.55 E',
        'changsha':     '28.23 N 112.94 E',
        'fuzhou':       '26.07 N 119.30 E',
        'guiyang':      '26.65 N 106.63 E',
        'lanzhou':      '36.06 N 103.79 E',
        'xining':       '36.62 N 101.77 E',
        'urumqi':       '43.79 N 87.58 E',
        'lhasa':        '29.65 N 91.10 E',
        'shenyang':     '41.80 N 123.43 E',
        'dalian':       '38.91 N 121.61 E',
        'qingdao':      '36.07 N 120.38 E',
        'suzhou':       '31.30 N 120.62 E',
        'shenzhen':     '22.54 N 114.06 E',
        'chongqing':    '29.56 N 106.55 E',
        'tianjin':      '39.13 N 117.18 E',
    }
    for city, latlon in china_cities.items():
        if city in loc_lower:
            return latlon
    # 如果 Location 包含 "China" 但没有具体城市, 用中国中心点
    if 'china' in loc_lower:
        return '35.86 N 104.19 E'
    return ''


def infer_host_from_name(taxonomy):
    """从病毒物种名推断宿主 (如 Betacytorhabdovirus lycii → Lycium)"""
    # 常见病毒-宿主对应模式
    host_patterns = [
        (r'\blycii\b', 'Lycium barbarum (goji)'),
        (r'\becsmenthae\b', 'Mentha (mint)'),
        (r'\bsolani\b', 'Solanum (nightshade)'),
        (r'\btritici\b', 'Triticum (wheat)'),
        (r'\bhordei\b', 'Hordeum (barley)'),
        (r'\bo\'ryzae\b', 'Oryza (rice)'),
        (r'\bvitis\b', 'Vitis (grape)'),
        (r'\bcitri\b', 'Citrus'),
        (r'\bnicotianae\b', 'Nicotiana (tobacco)'),
        (r'\bcapsici\b', 'Capsicum (pepper)'),
        (r'\bcucumeris\b', 'Cucumis (cucumber)'),
        (r'\bsojae\b', 'Glycine (soybean)'),
        (r'\bmaydis\b', 'Zea (corn)'),
    ]
    for pattern, host in host_patterns:
        if re.search(pattern, taxonomy, re.IGNORECASE):
            return host
    return ''


def load_ref_info(work_dir):
    """从 ref_info.tsv 加载物种级分类 (用于宿主推断)"""
    for fname in ['ref_info.tsv', 'integrated_summary.tsv']:
        path = Path(work_dir) / fname
        if path.exists():
            try:
                df = pd.read_csv(path, sep='\t')
                for col in ['Species', 'best_species', 'suvtk_species']:
                    if col in df.columns and len(df.columns) > 1:
                        contig_col = df.columns[0]
                        return dict(zip(df[contig_col].astype(str), df[col].astype(str)))
            except Exception:
                continue
    return {}


def generate_metadata_csv(seqs, meta_lookup, species_map, args, log):
    """生成统一元数据 CSV

    公共转录组数据的处理逻辑:
      - 采样信息 (日期/地点/组织/BioSample/BioProject) → 从 metadata (优先 Full.tsv) 自动填充
      - 测序平台 → 优先 metadata.Platform, 回退 --sequencer
      - BioProject/BioSample → 优先用 --submission-bioproject (你自己的)
      - 原始 SRA 的 BioProject 保存到 Note 列备查
    """
    sub_bp = args.submission_bioproject or args.bioproject or ''
    sub_bs = args.submission_biosample_prefix or ''

    rows = []
    for s in seqs:
        contig = s['contig']
        taxonomy = s['taxonomy']
        sra = extract_sra(contig)
        ref = meta_lookup.get(sra, {})

        # ── 采样信息: 从公共元数据自动填 ──
        collection_date = ref.get('collection_date', '')
        if collection_date and collection_date.lower() in ('nan', 'not_provided', ''):
            collection_date = ''

        geo_loc = ref.get('geo_loc_name', '')
        if geo_loc and geo_loc.lower() in ('nan', 'not_provided', ''):
            geo_loc = ''

        # 宿主: Core14 ScientificName > --host > auto-host (taxonomy + ref_info)
        host = ref.get('host', '')
        if not host or host.lower() in ('nan', 'not_provided', ''):
            host = args.host or ''
        if not host and hasattr(args, 'auto_host') and args.auto_host:
            species_name = species_map.get(contig, taxonomy)
            host = infer_host_from_name(species_name)
        tissue = ref.get('tissue', args.tissue or '')
        # Core14 额外字段 → NCBI source modifiers
        cultivar = ref.get('source', '')    # Source = "Ningqi No.1" 品种
        dev_stage = ref.get('age_stage', '') # Age_GrowthStage = "3 year"
        # 收集人: CenterName (如 Beijing Forestry University)
        collected_by = args.collected_by or ref.get('center', '')

        # ── BioProject: 优先用你自己的提交号, 原始号放 Note ──
        ref_bp = ref.get('bioproject', '')
        bioproject = sub_bp or ref_bp or 'PRJNAXXXXXX'

        ref_bs = ref.get('biosample', '')
        # 如果给了前缀, 自动编号
        if sub_bs and contig:
            biosample = f"{sub_bs}_{sra}"[:50]
        elif ref_bs and ref_bs.lower() not in ('nan', 'samnxxxxxxxx', ''):
            biosample = ref_bs
        else:
            biosample = 'SAMNXXXXXXXX'

        # ── Note: 记录原始数据的来源 ──
        note_parts = []
        if ref_bp and ref_bp != bioproject:
            note_parts.append(f"original_BioProject={ref_bp}")
        if ref_bs and ref_bs != biosample:
            note_parts.append(f"original_BioSample={ref_bs}")
        note = '; '.join(note_parts) if note_parts else ''

        isolate = f"{sanitize_filename(taxonomy)}_{sra}"

        rows.append({
            "organism": taxonomy,
            "sequence_name": contig,
            "authors": args.authors or "Last, First",
            "collection_date": collection_date or "YYYY-MM-DD",
            "bioproject": bioproject,

            "src-Isolate": isolate,
            "src-geo_loc_name": geo_loc or "Country:Region",
            "src-Lat_Lon": args.lat_lon or geo_to_latlon(geo_loc) or "XX.XX N XXX.XX E",
            "src-Host": host,
            "src-Segment": "",
            "src-Isolation-source": args.isolation_source or '',
            "src-Note": note,
            "src-Tissue_type": tissue,
            "src-Collected_by": collected_by,
            "src-Cultivar": cultivar,        # ← Core14 Source 字段
            "src-Dev_stage": dev_stage,       # ← Core14 Age_GrowthStage 字段

            "gb-sample_name": isolate[:50],
            "gb-title": args.title or f"{taxonomy} genome sequencing",
            "sra": sra,
            "biosample": biosample,

            "cmt-Assembly_Method": args.assembler or "SPAdes v4.3.0 (metaviral)",
            "cmt-Sequencing_Technology": ref.get('platform') or args.sequencer or "Illumina NovaSeq 6000",
            "cmt-Genome_Coverage": args.coverage or "",
            "cmt-Annotation_Pipeline": args.pipeline or "植物病毒分析平台 ⑥b ORF 功能注释 (DIAMOND/MMseqs2 × RefSeq 病毒蛋白库)",

            "bs-isolate": isolate,
            "bs-geo_loc_name": geo_loc or "Country:Region",
            "bs-host": host,
            "bs-isolation_source": args.metagenome_source or "plant virome",
        })

    df = pd.DataFrame(rows, columns=list(UNIFIED_COLUMNS.keys()))
    return df


def export_source_src(df, out_dir, log):
    """从统一 CSV 生成 source.src"""
    src_dir = Path(out_dir) / "source_individual"
    src_dir.mkdir(parents=True, exist_ok=True)

    # 完整 source.src (12 必填 + 5 推荐 NCBI modifiers)
    extra_cols = {
        'src-Host': 'Host',
        'src-Tissue_type': 'Tissue_type',
        'src-Cultivar': 'Cultivar',
        'src-Dev_stage': 'Dev_stage',
        'src-Collected_by': 'Collected_by',
    }
    # 检查哪些列有非空值
    filled_extras = []
    for csv_col, ncbi_name in extra_cols.items():
        if csv_col in df.columns and df[csv_col].notna().any() and (df[csv_col] != '').any():
            filled_extras.append((csv_col, ncbi_name))

    src_all = Path(out_dir) / "source.src"
    with safe_open(src_all, 'wt') as f:
        base_header = ("Sequence_ID\tOrganism\tIsolate\tCollection_date\tgeo_loc_name\tLat_Lon\t"
                       "Bioproject\tBiosample\tSRA\tMetagenomic\tMetagenome_source\tSegment")
        if filled_extras:
            base_header += "\t" + "\t".join(n for _, n in filled_extras)
        f.write(base_header + "\n")

        for _, row in df.iterrows():
            base_line = (f"{row['sequence_name']}\t{row['organism']}\t{row['src-Isolate']}\t"
                         f"{row['collection_date']}\t{row['src-geo_loc_name']}\t{row['src-Lat_Lon']}\t"
                         f"{row['bioproject']}\t{row['biosample']}\t{row['sra']}\tTRUE\t"
                         f"{row['src-Isolation-source'] or 'plant virome'}\t{row['src-Segment']}")
            if filled_extras:
                base_line += "\t" + "\t".join(str(row.get(c, '')) for c, _ in filled_extras)
            f.write(base_line + "\n")
    log.info("  → %s (%d 条)", src_all, len(df))

    # 按 virus 拆分
    by_org = df.groupby('organism')
    for org, group in by_org:
        fname = f"source_{sanitize_filename(org)}.src"
        src_virus = src_dir / fname
        with safe_open(src_virus, 'wt') as f:
            f.write("Sequence_ID\tOrganism\tIsolate\tCollection_date\tgeo_loc_name\tLat_Lon\t"
                    "Bioproject\tBiosample\tSRA\tMetagenomic\tMetagenome_source\tSegment\n")
            for _, row in group.iterrows():
                f.write(f"{row['sequence_name']}\t{row['organism']}\t{row['src-Isolate']}\t"
                        f"{row['collection_date']}\t{row['src-geo_loc_name']}\t{row['src-Lat_Lon']}\t"
                        f"{row['bioproject']}\t{row['biosample']}\t{row['sra']}\tTRUE\t"
                        f"{row['src-Isolation-source'] or 'plant virome'}\t{row['src-Segment']}\n")
        log.info("  → %s (%d 条)", src_virus, len(group))


def export_biosample_csv(df, out_dir, log):
    """生成 NCBI BioSample 批量提交模板 (Pathogen.cl.1.0 包)"""
    # BioSample Pathogen.cl.1.0 包的必填字段
    bs_mapping = {
        'sample_name':      'gb-sample_name',   # BioSample 记录名
        'organism':         'organism',          # 物种名
        'collection_date':  'collection_date',   # 采集日期
        'geo_loc_name':     'bs-geo_loc_name',   # 采集地点
        'isolation_source': 'bs-isolation_source', # 分离来源
        'isolate':          'bs-isolate',        # 分离株
        'host':             'bs-host',           # 宿主
        'bioproject':       'bioproject',        # BioProject
    }
    bs_df = pd.DataFrame()
    for bs_field, src_col in bs_mapping.items():
        if src_col in df.columns:
            bs_df[bs_field] = df[src_col]
    bs_path = Path(out_dir) / "biosample_template.tsv"
    with safe_open(bs_path, 'wt') as f:
        bs_df.to_csv(f, sep='\t', index=False)
    log.info("  → %s (%d 条, Pathogen.cl.1.0 格式)", bs_path, len(bs_df))


def create_authorset_sbt(config_dict, out_dir, metadata_df, log):
    """从配置生成 authorset.sbt (SeqSender 风格)"""
    sbt_path = Path(out_dir) / "authorset.sbt"

    submitter = config_dict.get('Submitter', {})
    org = config_dict.get('Organization', {})
    addr = config_dict.get('Address', {})
    authors_str = metadata_df['authors'].iloc[0] if 'authors' in metadata_df.columns else 'Author, First'

    # 解析作者列表
    author_list = [a.strip() for a in authors_str.replace(';', ',').split(',') if a.strip()]

    with safe_open(sbt_path, 'wt') as f:
        f.write("Submit-block ::= {\n")
        f.write("  contact {\n")
        f.write("    contact {\n")
        f.write("      name name {\n")
        f.write(f"        last \"{submitter.get('Last', 'LastName')}\",\n")
        f.write(f"        first \"{submitter.get('First', 'FirstName')}\"\n")
        f.write("      },\n")
        f.write("      affil std {\n")
        f.write(f"        affil \"{addr.get('Affil', 'Institution')}\",\n")
        f.write(f"        div \"{addr.get('Div', 'Department')}\",\n")
        f.write(f"        city \"{addr.get('City', 'City')}\",\n")
        f.write(f"        sub \"{addr.get('Sub', 'State')}\",\n")
        f.write(f"        country \"{addr.get('Country', 'Country')}\",\n")
        f.write(f"        street \"{addr.get('Street', '')}\",\n")
        f.write(f"        email \"{submitter.get('Email', 'email@example.com')}\",\n")
        f.write(f"        postal-code \"{addr.get('Postal_Code', '00000')}\"\n")
        f.write("      }\n")
        f.write("    }\n")
        f.write("  },\n")
        f.write("  cit {\n")
        f.write("    authors {\n")
        f.write("      names std {\n")
        for i, author in enumerate(author_list, 1):
            parts = author.strip().split()
            last = parts[-1] if parts else "Author"
            first = parts[0] if len(parts) >= 2 else ""
            f.write("        {\n")
            f.write("          name name {\n")
            f.write(f"            last \"{last}\",\n")
            f.write(f"            first \"{first}\"\n")
            f.write("          }\n")
            if i == len(author_list):
                f.write("        }\n")
            else:
                f.write("        },\n")
        f.write("      },\n")
        f.write("      affil std {\n")
        f.write(f"        affil \"{addr.get('Affil', 'Institution')}\",\n")
        f.write(f"        div \"{addr.get('Div', 'Department')}\",\n")
        f.write(f"        city \"{addr.get('City', 'City')}\",\n")
        f.write(f"        sub \"{addr.get('Sub', 'State')}\",\n")
        f.write(f"        country \"{addr.get('Country', 'Country')}\",\n")
        f.write(f"        street \"{addr.get('Street', '')}\"\n")
        f.write("      }\n")
        f.write("    },\n")
        f.write(f"    title \"{config_dict.get('Publication_Title', 'Viral genome sequencing and assembly')}\",\n")
        f.write(f"    status \"{config_dict.get('Publication_Status', 'Unpublished')}\"\n")
        f.write("  }\n")
        f.write("}\n")
    log.info("  → %s", sbt_path)


def create_submission_log(out_dir, run_title, log):
    """创建提交追踪日志 (submission_log.csv)"""
    log_path = Path(out_dir) / "submission_log.csv"
    exists = log_path.exists()
    with safe_open(log_path, 'at' if exists else 'wt') as f:
        writer = csv.writer(f)
        if not exists:
            writer.writerow(['timestamp', 'submission_name', 'status', 'n_sequences', 'biosamples', 'sqn_file', 'notes'])
        writer.writerow([
            datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            run_title,
            'FILES_GENERATED',
            '', '', '', 'Metadata files ready. Edit unified_metadata.csv, then re-run with --validate to check.'
        ])
    log.info("  → %s (更新)", log_path)


def load_config(config_path):
    """加载提交配置文件 (YAML 或 JSON)"""
    config = {
        'Submitter': {'First': 'FirstName', 'Last': 'LastName', 'Email': 'email@example.com'},
        'Address': {'Affil': 'Institution', 'Div': 'Department', 'City': 'City',
                     'Sub': 'State', 'Country': 'Country', 'Street': '', 'Postal_Code': '00000'},
        'Publication_Title': 'Viral genome sequencing and assembly',
        'Publication_Status': 'Unpublished',
        'Spuid_Namespace': '',
        'GenBank_Auto_Remove_Failed_Samples': True,
        'Specified_Release_Date': '',
    }
    if config_path and os.path.exists(config_path):
        if config_path.endswith('.yaml') or config_path.endswith('.yml'):
            import yaml
            with open(config_path) as f:
                user_config = yaml.safe_load(f)
        elif config_path.endswith('.json'):
            import json
            with open(config_path) as f:
                user_config = json.load(f)
        else:
            return config
        # 扁平化嵌套结构
        if 'Submission' in user_config and 'NCBI' in user_config['Submission']:
            ncbi = user_config['Submission']['NCBI']
        else:
            ncbi = user_config
        if 'Description' in ncbi:
            desc = ncbi['Description']
            if 'Organization' in desc:
                org = desc['Organization']
                config['Submitter'] = org.get('Submitter', config['Submitter'])
                config['Address'] = org.get('Address', config['Address'])
                config['Address']['Affil'] = org.get('Name', config['Address']['Affil'])
            config['Publication_Title'] = ncbi.get('Publication_Title', config['Publication_Title'])
            config['Publication_Status'] = ncbi.get('Publication_Status', config['Publication_Status'])
            config['Specified_Release_Date'] = ncbi.get('Specified_Release_Date', '')
            config['Spuid_Namespace'] = ncbi.get('Spuid_Namespace', '')
    return config


def inject_interactive_metadata(out_dir, df, log):
    """交互式补全缺失元数据 (VAPiD 风格)"""
    import shutil
    # 备份原始 CSV
    csv_path = Path(out_dir) / "unified_metadata.csv"
    bak_path = Path(out_dir) / "unified_metadata.csv.bak"
    if not bak_path.exists():
        shutil.copy(csv_path, bak_path)
        log.info("  备份: %s", bak_path)

    modified = False
    for col in REQUIRED_COLS:
        if col not in df.columns:
            continue
        mask = df[col].isna() | (df[col].astype(str).str.strip() == '') | \
               df[col].astype(str).str.contains('XXXX|YYYY|PRJNAXXXX|Country:Region|SAMNXXXXXXXX|Author')
        if not mask.any():
            continue

        n_missing = mask.sum()
        unique_vals = df.loc[mask, col].unique() if hasattr(df.loc[mask, col], 'unique') else []

        print(f"\n{'─'*50}")
        print(f"  [{col}] — {UNIFIED_COLUMNS.get(col, {}).get('desc', '')}")
        print(f"  缺失 {n_missing}/{len(df)} 条, 当前值: {unique_vals[:3]}")
        print(f"\n  输入新值 (回车保留原值, '{col}=ALL' 批填所有):")
        val = input(f"  > ").strip()

        if val:
            if f"{col}=ALL" in val.replace(' ', ''):
                # 批量填充
                all_val = val.split('=', 1)[1].strip()
                df.loc[mask, col] = all_val
                log.info("  批量填充 %s = %s (%d 条)", col, all_val, n_missing)
            else:
                # 只更新有具体序列的 — 这里交互式逐条填太慢, 简化为一键全填
                df.loc[mask, col] = val
                log.info("  填充 %s = %s (%d 条)", col, val, n_missing)
            modified = True

    if modified:
        with safe_open(csv_path, 'wt') as f:
            df.to_csv(f, index=False, encoding='utf-8-sig')
        log.info("  → 已保存更新: %s", csv_path)
    else:
        log.info("  无需修改")
    return df


def export_miuvig_assembly(out_dir, log, assembler, sequencer, enrichment,
                          env_broad_scale=None, env_local_scale=None, env_medium=None,
                          samp_taxon_id=None, project_name=None, source_uvig=None, mol_type='cRNA'):
    """生成 miuvig.tsv 和 assembly.tsv (覆盖 MIUVIG 全部必填字段)"""
    miuvig_path = Path(out_dir) / "miuvig.tsv"

    # MIUVIG 必填字段清单 (v6.3.0)
    miuvig_fields = [
        # === 必填 (Required) ===
        ("source_uvig",           source_uvig or "metatranscriptome (not viral targeted)"),  # 用户指定
        ("virus_enrich_appr",     enrichment or "none"),                       # 病毒富集方法
        ("assembly_software",     assembler),                                  # 组装软件;版本;参数
        ("assembly_qual",         "Genome fragment(s)"),                       # 组装质量
        ("vir_ident_software",    "kun_peng;0.7.12;BLAST+DIAMOND/MMseqs2"),    # 病毒识别工具
        ("pred_genome_type",      "uncharacterized"),                          # 由 taxonomy 覆盖
        ("pred_genome_struc",     "undetermined"),                             # 由 taxonomy 覆盖
        ("detec_type",            "independent sequence (UViG)"),              # 检测类型
        ("number_contig",         "1"),                                        # 由 comments 覆盖
        # === 必填 (环境) ===
        ("env_broad_scale",       env_broad_scale or "anthropogenic terrestrial biome [ENVO:01000219]"),
        ("env_local_scale",       env_local_scale or "agricultural field [ENVO:00000114]"),
        ("env_medium",            env_medium or "plant-associated soil [ENVO:00005789]"),
        ("geo_loc_name",          "placeholder — see source.src"),             # 每条序列不同
        ("collection_date",       "placeholder — see source.src"),
        ("lat_lon",               "placeholder — see source.src"),
        # === 必填 (测序 + 样本) ===
        ("samp_taxon_id",         samp_taxon_id or "plant metagenome [NCBITaxon:1297885]"),
        ("project_name",          project_name or "plant virus analysis platform metagenomics"),
        ("seq_meth",              sequencer),                                  # 测序仪
        ("samp_name",             "placeholder — see source.src"),             # 每条序列不同
        # === 推荐 (Recommended) ===
        ("feat_pred",             "orfipy;pyrodigal;pyrodigal-rv"),
        ("ref_db",                "RefSeq viral proteins;DIAMOND/MMseqs2 + HMM"),
        ("sim_search_meth",       "DIAMOND/MMseqs2;platform;default parameters"),
        ("tax_class",             "kunpeng LCA + BLAST;platform;ICTV MSR;default parameters"),
        ("otu_class_appr",        "95% ANI;85% AF;MMseqs2 easy-cluster"),
        ("compl_score",           ""),  # CheckV 填充
        ("compl_software",        ""),  # CheckV 填充
    ]

    # 去重
    seen = set()
    unique_fields = []
    for k, v in miuvig_fields:
        if k not in seen:
            seen.add(k)
            unique_fields.append((k, v))

    with safe_open(miuvig_path, 'wt') as f:
        f.write("MIUVIG_parameter\tvalue\n")
        for param, val in unique_fields:
            f.write(f"{param}\t{val}\n")
    log.info("  → %s (%d MIUVIG 字段)", miuvig_path, len(unique_fields))

    asm_path = Path(out_dir) / "assembly.tsv"
    with safe_open(asm_path, 'wt') as f:
        f.write("Assembly_parameter\tvalue\n")
        f.write("StructuredCommentPrefix\tAssembly-Data\n")
        f.write(f"Assembly Method\t{assembler}\n")
        f.write(f"Sequencing Technology\t{sequencer}\n")
    log.info("  → %s", asm_path)


def validate_csv(df, out_dir, log):
    """检查必填字段, 输出验证报告"""
    report_path = Path(out_dir) / "validation_report.txt"
    issues = []
    for col in REQUIRED_COLS:
        if col not in df.columns:
            issues.append(f"  [MISSING] 缺少必填列: {col}")
            continue
        missing = df[col].isna() | (df[col].astype(str).str.strip() == '') | \
                  df[col].astype(str).str.contains('XXXX|YYYY-MM-DD|Country:Region|PRJNAXXXX|SAMNXXXXXXXX')
        if missing.any():
            seqs = df.loc[missing, 'sequence_name'].tolist()
            issues.append(f"  [EMPTY] {col}: {len(seqs)} 条序列未填 ({seqs[:3]}...)")

    with safe_open(report_path, 'wt') as f:
        f.write("元数据验证报告\n")
        f.write(f"生成: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"总序列: {len(df)}\n")
        f.write(f"必填列: {len(REQUIRED_COLS)}\n")
        f.write("─" * 50 + "\n")
        if issues:
            f.write(f"⚠️  发现 {len(issues)} 个问题:\n")
            for i in issues:
                f.write(i + "\n")
            f.write("\n提交前请修复以上问题!\n")
        else:
            f.write("✓ 所有必填字段已填写, 可以提交\n")

    if issues:
        log.warning("  验证: %d 个问题 → %s", len(issues), report_path)
    else:
        log.info("  验证: 全部通过 ✓ → %s", report_path)


def main():
    parser = argparse.ArgumentParser(
        description="unified_metadata.py — 统一元数据模板生成器 v1.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--taxonomy', required=True, help='taxonomy.tsv (contig + taxonomy)')
    parser.add_argument('--metadata', help='Global_Unified_Metadata.tsv 或 Core14.tsv (自动优先同目录 Full.tsv)')
    parser.add_argument('--run-title', default='viral_submission', help='运行标题')
    parser.add_argument('-o', '--output', required=True, help='输出目录（限平台目录内）')
    parser.add_argument('--assembler', default='SPAdes;4.3.0;metaviral')
    parser.add_argument('--sequencer', default='Illumina NovaSeq 6000')
    parser.add_argument('--enrichment', default='rRNA depletion')
    parser.add_argument('--pipeline', default='植物病毒分析平台 ⑥b ORF 功能注释 (DIAMOND/MMseqs2 × RefSeq 病毒蛋白库)')
    parser.add_argument('--authors', help='引用作者')
    parser.add_argument('--bioproject', help='BioProject ID')
    parser.add_argument('--lat-lon', help='经纬度')
    parser.add_argument('--isolation-source', help='分离来源')
    parser.add_argument('--collected-by', help='采集人')
    parser.add_argument('--tissue', help='组织类型')
    parser.add_argument('--coverage', help='覆盖度')
    parser.add_argument('--title', help='提交标题')
    parser.add_argument('--metagenome-source', help='宏基因组来源')
    parser.add_argument('--config', help='提交配置文件 (YAML), 如 SeqSender 的 seqsender_config.yaml')
    parser.add_argument('--interactive', action='store_true',
                        help='交互式补全缺失元数据 (VAPiD 风格)')
    parser.add_argument('--skip-validate', action='store_true', help='跳过验证')
    parser.add_argument('--reference-metadata', dest='ref_meta',
                        help='public_metadata_pipeline 输出的 Core14.tsv (参考元数据: 日期/地点/组织)')
    parser.add_argument('--submission-bioproject',
                        help='你自己的 BioProject ID (PRJNA...) — 公共数据的 PRJNA 仅供参考')
    parser.add_argument('--submission-biosample-prefix',
                        help='你自己的 BioSample 前缀 (SAMN...) — 若未注册则留空, 用模板注册')
    parser.add_argument('--env-broad-scale',
                        help='MIUVIG env_broad_scale (如 anthropogenic terrestrial biome [ENVO:01000219])')
    parser.add_argument('--env-local-scale',
                        help='MIUVIG env_local_scale (如 agricultural field [ENVO:00000114])')
    parser.add_argument('--env-medium',
                        help='MIUVIG env_medium (如 leaf [PO:0025034] / fruit)')
    parser.add_argument('--samp-taxon-id',
                        help='MIUVIG samp_taxon_id (如 Lycium barbarum [NCBITaxon:112863])')
    parser.add_argument('--project-name',
                        help='MIUVIG project_name (如 plant virome survey)')
    parser.add_argument('--source-uvig',
                        default='metatranscriptome (not viral targeted)',
                        help='MIUVIG source_uvig: 公共数据=metatranscriptome, 自己测的富集=viral fraction metagenome (virome)')
    parser.add_argument('--host', help='所有序列的宿主 (如 Lycium barbarum [NCBITaxon:112863])')
    parser.add_argument('--auto-host', action='store_true',
                        help='从病毒名自动推断宿主 (Betacytorhabdovirus lycii → Lycium)')
    parser.add_argument('--coverage-file',
                        help='覆盖度文件 (如 all_viruses.best.summary.tsv, 自动取 Rep_MeanDepth)')
    args = parser.parse_args()

    out_dir = Path(args.output)
    log = setup_logger(str(out_dir))

    log.info("=" * 60)
    log.info("统一元数据模板生成: %s", args.run_title)

    # 加载数据
    seqs = load_taxonomy(args.taxonomy)
    meta_lookup = load_metadata(args.metadata) if args.metadata else {}
    log.info("  序列: %d, 元数据映射: %d", len(seqs), len(meta_lookup))

    # 生成统一 CSV
    species_map = {}
    if hasattr(args, 'auto_host') and args.auto_host:
        # 从 taxonomy.tsv 所在目录往上找 ref_info.tsv
        tax_dir = Path(args.taxonomy).parent
        species_map = load_ref_info(str(tax_dir))
        if not species_map:
            species_map = load_ref_info(str(tax_dir.parent))  # 再往上一级
    # 自动检测覆盖度
    if not args.coverage and args.coverage_file:
        try:
            cov_df = pd.read_csv(args.coverage_file, sep='\t')
            for col in ['Rep_MeanDepth', 'Rep_Coverage(%)', 'Covered%']:
                if col in cov_df.columns:
                    args.coverage = f"{cov_df[col].dropna().mean():.1f}x"; break
        except Exception: pass

    df = generate_metadata_csv(seqs, meta_lookup, species_map, args, log)
    csv_path = out_dir / "unified_metadata.csv"
    with safe_open(csv_path, 'wt') as f:
        df.to_csv(f, index=False, encoding='utf-8-sig')
    log.info("  → %s (%d 行 × %d 列)", csv_path, len(df), len(df.columns))

    # 导出 source.src
    export_source_src(df, out_dir, log)

    # 导出 BioSample 模板
    export_biosample_csv(df, out_dir, log)

    # 导出 miuvig/assembly (含全部 MIUVIG 必填字段)
    # 从 Core14 自动推断 samp_taxon_id
    samp_taxon = args.samp_taxon_id or ''
    if not samp_taxon and meta_lookup:
        # 取第一条有 ScientificName 的记录
        for v in meta_lookup.values():
            if v.get('host') and v['host'].lower() not in ('nan', 'not_provided', ''):
                taxid = v.get('taxid', '')
                samp_taxon = f"{v['host']}"
                if taxid:
                    samp_taxon += f" [NCBITaxon:{taxid}]"
                break

    export_miuvig_assembly(out_dir, log, args.assembler, args.sequencer, args.enrichment,
                           args.env_broad_scale, args.env_local_scale, args.env_medium,
                           samp_taxon, args.project_name or args.run_title,
                           source_uvig=args.source_uvig)

    # 加载配置 → 生成 authorset.sbt
    config = load_config(args.config)
    create_authorset_sbt(config, out_dir, df, log)

    # 提交追踪日志
    create_submission_log(out_dir, args.run_title, log)

    # 交互式补全 (VAPiD 风格)
    if args.interactive:
        df = inject_interactive_metadata(out_dir, df, log)

    # 验证
    if not args.skip_validate:
        validate_csv(df, out_dir, log)

    # 打印结果
    auto_filled = sum(1 for c in REQUIRED_COLS if c in df.columns and not df[c].astype(str).str.contains('XXXX|YYYY|PRJNA|Country:Region').any())
    print(f"\n{'='*60}")
    print(f"  统一元数据文件已生成")
    print(f"  {'='*60}")
    print(f"  输出: {out_dir}/")
    print(f"    unified_metadata.csv     ← Excel 编辑 (一张表驱动一切)")
    print(f"    source.src               ← GenBank 全部序列")
    print(f"    source_individual/       ← 每个病毒独立")
    print(f"    biosample_template.tsv   ← NCBI BioSample 批量注册")
    print(f"    authorset.sbt            ← 作者信息 (从 config 生成)")
    print(f"    submission_log.csv       ← 提交追踪日志")
    print(f"    miuvig.tsv / assembly.tsv")
    print(f"    validation_report.txt    ← 字段校验")
    print(f"  {'='*60}")
    print(f"  序列数: {len(df)}")
    print(f"  必填字段自动填充: {auto_filled}/{len(REQUIRED_COLS)}")
    print(f"  {'='*60}")
    print(f"")
    print(f"  提交流程:")
    print(f"  1. 编辑 unified_metadata.csv 占位符 (或用 --interactive 交互补全)")
    print(f"  2. 在 NCBI 注册 BioProject → 填入 bioproject 列")
    print(f"  3. 用 biosample_template.tsv 批量注册 BioSample → 填入 biosample 列")
    print(f"  4. NCBI BankIt 上传 .fsa + .tbl (或 Sequin 导入生成 .sqn)")


if __name__ == '__main__':
    main()

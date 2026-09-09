#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_plot.py — 深度绘制段
======================
per-virus 测序深度曲线 + GenBank 基因轨道

设计来源: virome_analysis_pipeline/batch_plot_virus_depth.py --mode depth
（复制后重构，原文件不动）

保留的原管线逻辑:
  - smooth() 滑动平均 + fill_between 填充
  - MeanDepth 红色虚线基准线
  - 左上角统计信息框（species/coverage/depth/reads）
  - GenBank 基因轨道：CDS/gene 坐标 + 链方向箭头多边形
  - 出版级样式：pdf.fonttype 42、白底、Arial

移植差异（大王指定：不用 sh/awk/sed/grep/gzip/rg）:
  - 原管线用 `gzip -dc | rg -f patterns` 抽取深度行 → 改为 Python 逐行过滤
  - 原管线用 `efetch`（NCBI edirect）下载 GenBank → 改为 Biopython Bio.Entrez
  - 输入表列名兼容 Accession / Rep_Accession → Virus（同原管线 _load_summary）
  - 深度表来源：pandepth 的 *.SiteDepth.gz（chr position depth 三列）
"""

import gc
import gzip
import os
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

plt.rcParams.update({
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.family': 'sans-serif',
    'axes.unicode_minus': False,
    'figure.facecolor': 'white',
    'savefig.facecolor': 'white',
})


# ── 原管线同名工具函数（逐字保留）──────────────────────────────
def smooth(y, box_pts):
    """滑动平均，边缘用 edge padding 避免收缩"""
    if box_pts <= 1 or len(y) == 0:
        return np.asarray(y, dtype=float)
    box = np.ones(box_pts) / box_pts
    y_pad = np.pad(y, (box_pts // 2, box_pts - 1 - box_pts // 2), mode='edge')
    return np.convolve(y_pad, box, mode='valid')


def safe_name(s):
    """文件名净化（与原管线一致）"""
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in str(s))


def create_info_text(v_stat):
    """左上角统计框文本（字段名照原管线）"""
    info = []
    tax = v_stat.get('taxonomy', 'Unknown')
    if tax and tax != 'Unknown':
        info.append(f"Species: {str(tax)[:40]}")
    info.append(f"Coverage: {v_stat.get('Coverage(%)', 0):.1f}%")
    info.append(f"MeanDepth: {v_stat.get('MeanDepth', 0):.2f}x")
    info.append(f"Reads: {int(v_stat.get('MappedReads', v_stat.get('Uniq_Reads', 0)))}")
    if 'TPM' in v_stat and pd.notna(v_stat.get('TPM')):
        info.append(f"TPM: {v_stat.get('TPM', 0):.2f}")
    return "\n".join(info)


# ── 数据读取（原管线用 gzip+rg 管道，这里改纯 Python）──────────
def load_depth_table(depth_file, keep_chroms, logger=None):
    """读 per-position 深度表 -> DataFrame(chr, position, depth)

    原管线: `gzip -dc x.SiteDepth.gz | rg -F -f patterns`
    这里逐行过滤，语义等价：只保留 chr 在 keep_chroms 内的行。
    """
    keep = set(str(c) for c in keep_chroms)
    rows = []
    opener = gzip.open if str(depth_file).endswith('.gz') else open
    with opener(depth_file, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                continue
            chrom = parts[0].strip()
            if chrom not in keep:
                continue
            try:
                rows.append((chrom, int(parts[1]), int(parts[2])))
            except ValueError:
                continue
    df = pd.DataFrame(rows, columns=['chr', 'position', 'depth'])
    if logger:
        logger.debug(f"    深度表 {os.path.basename(str(depth_file))}: "
                     f"{len(df)} 行 / 保留 {len(keep)} 个参考")
    return df


def find_depth_file(depth_dir, sample):
    """按原管线的候选顺序找深度文件"""
    candidates = [
        f"{sample}.SiteDepth", f"{sample}.chr.SiteDepth",
        f"{sample}.SiteDepth.gz", f"{sample}.chr.SiteDepth.gz",
    ]
    for name in candidates:
        p = os.path.join(str(depth_dir), name)
        if os.path.exists(p):
            return p
    # 兜底：pandepth 带前缀的产物名，如 S1_pandepth.pandepth.SiteDepth.gz
    d = Path(depth_dir)
    if d.is_dir():
        hits = sorted(d.glob(f"{sample}*SiteDepth*"))
        if hits:
            return str(hits[0])
    return None


# ── GenBank 注释（原管线用 efetch，这里用 Bio.Entrez）──────────
def fetch_gb(accession, save_dir, logger=None, email=None, api_key=None):
    """下载 GenBank 记录到缓存目录，已存在则复用。

    原管线: `efetch -db nuccore -id <acc> -format gb`
    这里用 Bio.Entrez.efetch，同一 NCBI 接口。
    下载失败返回 None（不抛异常，保持原管线行为）。
    """
    os.makedirs(str(save_dir), exist_ok=True)
    out_file = os.path.join(str(save_dir), f"{accession}.gb")
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return out_file
    try:
        from Bio import Entrez
    except ImportError:
        if logger:
            logger.warning("  Biopython 缺失，无法下载 GenBank 注释")
        return None
    Entrez.email = email or os.environ.get('NCBI_EMAIL',
                                           'anonymous@example.com')
    if api_key or os.environ.get('NCBI_API_KEY'):
        Entrez.api_key = api_key or os.environ.get('NCBI_API_KEY')
    try:
        with Entrez.efetch(db='nuccore', id=str(accession).split('.')[0],
                           rettype='gb', retmode='text') as handle:
            data = handle.read()
        if isinstance(data, bytes):
            data = data.decode('utf-8', 'replace')
        if not data or 'LOCUS' not in data:
            return None
        with open(out_file, 'w', encoding='utf-8') as f:
            f.write(data)
        if logger:
            logger.debug(f"    GenBank 下载成功: {accession} ({len(data)}B)")
        return out_file
    except Exception as e:
        if logger:
            logger.warning(f"    GenBank 下载失败 {accession}: {e}")
        if os.path.exists(out_file):
            os.remove(out_file)
        return None


def parse_gb_features(gb_path, logger=None):
    """解析 CDS/gene 坐标（逐字照抄原管线 parse_gb_with_biopython）"""
    try:
        from Bio import SeqIO
    except ImportError:
        if logger:
            logger.warning('  Biopython 缺失，无法解析 GenBank 注释')
        return []
    features = []
    try:
        for record in SeqIO.parse(gb_path, "genbank"):
            for feat in record.features:
                if feat.type in ["CDS", "gene"]:
                    name = feat.qualifiers.get(
                        "gene", feat.qualifiers.get("product", ["Unknown"]))[0]
                    features.append({
                        'type': feat.type,
                        'start': int(feat.location.start),
                        'end': int(feat.location.end),
                        'strand': feat.location.strand,
                        'name': name,
                    })
    except Exception as e:  # noqa: BLE001
        # 解析失败与“该参考真的没有 CDS”（viroid）是两回事，必须区分
        if logger:
            logger.warning(f'    GenBank 解析失败 {os.path.basename(str(gb_path))}: {e}')
        return []
    # 去重：同一坐标优先 CDS，其次有名字的
    seen = {}
    for f in features:
        if f['name'] == 'Unknown':
            continue
        coord = (f['start'], f['end'])
        if coord not in seen:
            seen[coord] = f
        else:
            if f['type'] == 'CDS' and seen[coord]['type'] == 'gene':
                seen[coord] = f
            elif seen[coord]['name'] == 'Unknown' and f['name'] != 'Unknown':
                seen[coord]['name'] = f['name']
    final = list(seen.values())
    final.sort(key=lambda x: x['start'])
    return final


# ── 单图绘制（绘图代码照抄原管线 process_single_sample）────────
GENE_COLORS = ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462',
               '#b3de69', '#fccde5', '#d9d9d9']


def plot_one_virus(sample, virus, chrom_df, v_stat, virus_genes, out_dir,
                   window=10, fontsize=9):
    """画单个病毒的深度图，返回 (pdf_path, png_path) 或 None"""
    if chrom_df.empty:
        return None
    x = chrom_df['position'].values
    y = chrom_df['depth'].values
    mean_depth = float(v_stat.get('MeanDepth', 0) or 0)
    tax = str(v_stat.get('taxonomy', 'Unknown'))
    info_text = create_info_text(v_stat)
    title = (f"Sample: {sample}\n"
             f"{textwrap.shorten(tax, width=65, placeholder='...')} ({virus})")
    safe_tax = safe_name(tax)
    safe_vname = safe_name(virus)
    # 未注释参考的 taxonomy 会等于 accession，避免文件名重复
    stem = safe_vname if safe_tax == safe_vname else f"{safe_tax}_{safe_vname}"
    pdf_path = os.path.join(out_dir, f"{sample}_{stem}_depth.pdf")
    png_path = os.path.join(out_dir, f"{sample}_{stem}_depth.png")

    y_smooth = smooth(y, window)
    max_x = (max(x) if len(x) > 0 else
             (virus_genes[-1]['end'] if virus_genes else 1000))

    if virus_genes:
        fig, (ax, ax_g) = plt.subplots(
            2, 1, figsize=(12, 7.5), sharex=True,
            gridspec_kw={'height_ratios': [5, 1], 'hspace': 0.08})
    else:
        fig, ax = plt.subplots(figsize=(12, 6))
        ax_g = None

    ax.fill_between(x, 0, y_smooth, alpha=0.3, color='#1f77b4')
    ax.plot(x, y_smooth, color='#1f77b4', linewidth=1.5,
            label=f'Smoothed Depth (Window={window})')
    ax.axhline(y=mean_depth, color='#d62728', linestyle='--', linewidth=2,
               label=f'Mean Depth: {mean_depth:.2f}x')
    ax.annotate(info_text, xy=(0.02, 0.96), xycoords='axes fraction',
                bbox=dict(boxstyle='round,pad=0.6', fc='#f8f9fa',
                          ec='#ced4da', alpha=0.9),
                fontsize=fontsize, family='monospace', ha='left', va='top')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.set_ylabel('Sequencing Depth (x)', fontsize=12)
    ax.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98), framealpha=0.9)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.set_ylim(bottom=0)
    ax.set_xlim(0, max_x)

    if ax_g is not None:
        y_center, height = 0.5, 0.4
        for i, gene in enumerate(virus_genes):
            start, end = gene['start'], gene['end']
            head_length = min(max_x * 0.015, (end - start) * 0.4)
            y_top, y_bottom = y_center + height / 2, y_center - height / 2
            if gene['strand'] >= 0:
                x_coords = [start, end - head_length, end, end - head_length, start]
            else:
                x_coords = [end, start + head_length, start, start + head_length, end]
            poly = patches.Polygon(
                np.column_stack((x_coords,
                                 [y_bottom, y_bottom, y_center, y_top, y_top])),
                closed=True, facecolor=GENE_COLORS[i % len(GENE_COLORS)],
                edgecolor='#444444', alpha=0.9)
            ax_g.add_patch(poly)
            short_name = (gene['name'][:10] + '..'
                          if len(gene['name']) > 12 else gene['name'])
            ax_g.text(start + (end - start) / 2, y_center, short_name,
                      ha='center', va='center', fontsize=8, fontweight='bold')
        ax_g.set_ylim(0, 1)
        ax_g.set_yticks([])
        for spine in ['top', 'right', 'left']:
            ax_g.spines[spine].set_visible(False)
        ax_g.set_xlabel('Genome Position (bp)', fontsize=12)
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel('Genome Position (bp)', fontsize=12)

    plt.tight_layout()
    try:
        plt.savefig(pdf_path, bbox_inches='tight', dpi=300)
        plt.savefig(png_path, bbox_inches='tight', dpi=300)
    finally:
        # 成百上千张图跑 bulk 时，savefig 抛异常不做 close 会累积泄漏 figure
        plt.close(fig)
    return pdf_path, png_path


class PlotStage:
    """深度绘制段：读鉴定/过滤产物 + pandepth 深度表，出 per-virus 图"""

    def __init__(self, args, tools, logger, out_dir, cfg=None):
        self.args = args
        self.tools = tools
        self.logger = logger
        self.out_dir = Path(out_dir)
        self.cfg = cfg or {}

    def run(self):
        logger = self.logger
        cfg = self.cfg or {}
        logger.info('=' * 60)
        logger.info('【绘图段】照搬 batch_plot_virus_depth.py --mode depth')
        logger.info('=' * 60)

        window = int(cfg.get('window', 10))
        fontsize = int(cfg.get('fontsize', 9))
        add_genes = bool(cfg.get('add_genes', True))
        gbk_dir = Path(cfg.get('gbk_dir') or (self.out_dir / 'gbk_files'))
        email = cfg.get('ncbi_email') or os.environ.get('NCBI_EMAIL')
        api_key = cfg.get('ncbi_api_key') or os.environ.get('NCBI_API_KEY')

        # 输入表：优先用过滤产物（已确诊），退化到全量表
        summary = None
        if getattr(self.args, 'summary', None):
            summary = Path(self.args.summary)
        if not summary or not summary.exists():
            for cand in (self.out_dir / 'filter' / 'filtered.tsv',
                         self.out_dir / 'identify' / 'all_viruses.best.summary.tsv',
                         self.out_dir / 'identify' / 'all_viruses.summary.tsv'):
                if cand.exists():
                    summary = cand
                    break
        if not summary or not summary.exists():
            logger.warning('  未找到鉴定/过滤结果表，跳绘图段')
            return {'plots': 0}

        df = pd.read_csv(summary, sep='\t')
        df.columns = [c.strip().replace('\ufeff', '') for c in df.columns]
        # 列名归一（同原管线 _load_summary 的别名规则）
        for src, tgt in (('Accession', 'Virus'), ('Rep_Accession', 'Virus'),
                         ('Species', 'taxonomy'), ('Adjusted_Species', 'taxonomy'),
                         ('Species_ICTV', 'taxonomy')):
            if src in df.columns and tgt not in df.columns:
                df = df.rename(columns={src: tgt})
        if 'Virus' not in df.columns:
            logger.error('  结果表缺少 Virus/Accession 列，跳绘图段')
            return {'plots': 0}

        depth_dir = getattr(self.args, 'depth_dir', None) or (self.out_dir / 'align')
        logger.info(f'  结果表: {summary}  ({len(df)} 行)')
        logger.info(f'  深度表目录: {depth_dir}')

        # GenBank 注释
        genes_dict = {}
        if add_genes:
            uniq = sorted(set(df['Virus'].astype(str).str.strip()))
            logger.info(f'  准备 GenBank 注释 {len(uniq)} 条 -> {gbk_dir}')
            ok = 0
            for acc in uniq:
                gb = fetch_gb(acc, gbk_dir, logger, email, api_key)
                if gb:
                    feats = parse_gb_features(gb, logger)
                    if feats:
                        genes_dict[acc] = feats
                        genes_dict[acc.split('.')[0]] = feats
                        ok += 1
            logger.info(f'  注释就绪: {ok}/{len(uniq)} 条')

        samples = sorted(set(df['Sample'].astype(str)))
        n_plots = 0
        for sample in samples:
            sub = df[df['Sample'].astype(str) == sample]
            depth_file = find_depth_file(depth_dir, sample)
            if not depth_file:
                logger.warning(f'  [{sample}] 找不到深度文件，跳过')
                continue
            chroms = sub['Virus'].astype(str).str.strip().tolist()
            # pandepth 深度表的 chr 名可能带/不带版本号，两种都收
            keep = set(chroms)
            keep |= {c.split('.')[0] for c in chroms}
            depth_df = load_depth_table(depth_file, keep, logger)
            if depth_df.empty:
                logger.warning(f'  [{sample}] 深度表无匹配坐标行，跳过')
                continue
            out_sample = self.out_dir / 'plots' / sample
            out_sample.mkdir(parents=True, exist_ok=True)
            for rec in sub.to_dict('records'):
                virus = str(rec.get('Virus', '')).strip()
                chrom_df = depth_df[
                    (depth_df['chr'] == virus) |
                    (depth_df['chr'] == virus.split('.')[0])]
                if chrom_df.empty:
                    logger.debug(f'  [{sample}] {virus}: 无坐标行')
                    continue
                genes = genes_dict.get(virus,
                                       genes_dict.get(virus.split('.')[0], []))
                res = plot_one_virus(sample, virus, chrom_df, rec, genes,
                                     str(out_sample), window, fontsize)
                if res:
                    n_plots += 1
                    logger.info(f'  [{sample}] {virus}: '
                                f'{len(chrom_df)} 位点'
                                f'{" + " + str(len(genes)) + " 基因" if genes else ""}'
                                f' -> {os.path.basename(res[0])}')
            del depth_df
            gc.collect()

        logger.info(f'  深度图完成: {n_plots} 张 -> {self.out_dir / "plots"}')
        return {'plots': n_plots}

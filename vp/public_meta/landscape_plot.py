#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
📊 SCI 级宏观组学数据可视化引擎 (精修美化终极版)
核心优化：
1. 折线图 (Panel A) 审美重构：科学蓝主色调，白底空心标记点，通透感增强。
2. 环形图 (Panel B) 标签重构：强制将 数据库名、比例、数量 三合一居中印在彩色圆环上。
3. 文本换行与防遮挡设计完美保留。
4. 除 3x2 合并总览图外，同时输出 A~F 六张独立图（PNG 600dpi + PDF 矢量）。

（移植自 MMPV-RNA public_metadata_pipeline/gsa_sra.plot.py；已用 Agg 非交互后端）
"""

import os
import sys
if sys.platform == 'win32' and sys.stdout is not None:
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

import re
import math
import argparse
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, avoids Tk memory issues
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

# ==========================================
# 0. SCI 期刊全局格式设置与安全回退
# ==========================================
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"]
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False

NULL_WORDS = ['not_provided', 'nan', 'none', 'unknown', 'missing', '', 'not applicable', 'not collected']

COLORS_DB = {'SRA': '#4C72B0', 'GSA': '#C44E52'}

def clean_series(series):
    return series.dropna()[~series.dropna().astype(str).str.strip().str.lower().isin(NULL_WORDS)]

def wrap_labels(labels, width=35):
    return [textwrap.fill(str(label), width=width) for label in labels]

def add_bar_labels(ax, values, is_horizontal=True):
    max_val = max(values) if len(values) > 0 else 1
    offset = max_val * 0.02
    for i, v in enumerate(values):
        if is_horizontal:
            ax.text(v + offset, i, str(v), va='center', ha='left', fontsize=12, color='black')
        else:
            ax.text(i, v + offset, str(v), va='bottom', ha='center', fontsize=12, color='black')

def _prep(df):
    """列名兼容与深度清洗（时间/机构/组织/地理/发育时期），返回副本。"""
    df = df.copy()
    # 兼容不同来源的列名 (info Core14 / search 旧版 / search 新版)
    if 'Database' not in df.columns and 'Run' in df.columns:
        df['Database'] = df['Run'].astype(str).str.extract(r'^([A-Za-z]+)')[0].map(
            {'SRR': 'SRA', 'ERR': 'SRA', 'DRR': 'SRA'}).fillna('GSA')
    if 'CenterName' in df.columns:
        df['Organization_CenterName'] = df['CenterName']
    elif 'Organization_CenterName' not in df.columns:
        df['Organization_CenterName'] = pd.NA

    # 1. 深度数据聚合与纠错清洗
    df['ReleaseDate'] = pd.to_datetime(df['ReleaseDate'], errors='coerce')
    df['Year'] = df['ReleaseDate'].dt.year.fillna(0).astype(int)

    if 'Organization_CenterName' in df.columns:
        orgs = df['Organization_CenterName'].astype(str).str.strip().str.title()
        orgs = orgs.str.replace('Unversity', 'University', flags=re.IGNORECASE)
        orgs = orgs.str.replace('&Amp;', '&', flags=re.IGNORECASE)
        orgs = orgs.str.replace(r'\s+', ' ', regex=True)  # collapse whitespace
        df['Organization_CenterName'] = orgs

    if 'Tissue' in df.columns:
        df['Tissue'] = df['Tissue'].astype(str).str.strip().str.title()

    if 'Location' in df.columns:
        def clean_loc(val):
            val = str(val).strip(' "')
            if val.lower() in NULL_WORDS or 'missing' in val.lower():
                return pd.NA
            val = val.replace(':', ', ')
            val = re.sub(r'_Ai$', '', val, flags=re.IGNORECASE)
            val = re.sub(r'\s+', ' ', val)
            return val.title().strip()
        df['Location_Clean'] = df['Location'].apply(clean_loc)

    if 'Age_GrowthStage' in df.columns:
        def clean_age(val):
            if pd.isna(val): return pd.NA
            parts = [p.strip().title() for p in str(val).split('|') if p.strip().lower() not in NULL_WORDS]
            return " | ".join(parts) if parts else pd.NA
        df['Age_GrowthStage_Clean'] = df['Age_GrowthStage'].apply(clean_age)
    return df

# ==========================================
# 1. 六个面板的独立绘制函数（合并图与单图共用）
# ==========================================
def _pnl_a_temporal(ax, df):
    """A. 时间分布趋势。"""
    df_valid_years = df[df['Year'] > 2000]
    year_counts = df_valid_years['Year'].value_counts().sort_index()
    if year_counts.empty:
        ax.text(0.5, 0.5, 'No valid release dates', ha='center', va='center',
                transform=ax.transAxes, fontsize=13, color='#888888')
        ax.set_title('A. Temporal Distribution of Sequencing Data', loc='left',
                     fontsize=18, fontweight='bold')
        return
    ax.plot(year_counts.index, year_counts.values, marker='o', linestyle='-', linewidth=3.5,
            markersize=10, color='#0072B2', markerfacecolor='white', markeredgewidth=2.5, zorder=3)
    ax.fill_between(year_counts.index, year_counts.values, color='#0072B2', alpha=0.15, zorder=2)
    ax.grid(axis='y', linestyle='--', alpha=0.6, color='#D3D3D3', zorder=1)
    for x, y in zip(year_counts.index, year_counts.values):
        ax.text(x, y + (max(year_counts.values)*0.03), str(y), ha='center', va='bottom',
                fontsize=13, fontweight='bold', color='#333333')
    ax.set_title('A. Temporal Distribution of Sequencing Data', loc='left', fontsize=18, fontweight='bold')
    ax.set_xlabel('Release Year', fontsize=14)
    ax.set_ylabel('Number of Runs', fontsize=14, labelpad=12,
                  rotation=0, ha='right', va='center')
    ax.yaxis.set_label_coords(-0.12, 0.5)
    ax.tick_params(axis='x', rotation=45, labelsize=12)
    ax.set_xticks(year_counts.index)

def _pnl_b_database(ax, df):
    """B. 数据来源比例环形图。"""
    db_counts = df['Database'].value_counts()
    total = db_counts.sum()
    wedges, _ = ax.pie(
        db_counts,
        labels=None,
        startangle=140,
        colors=[COLORS_DB.get(x, '#555555') for x in db_counts.index],
        wedgeprops=dict(width=0.55, edgecolor='w', linewidth=3),
        radius=1.0,
    )
    ring_center = 1.0 - 0.55 / 2
    label_r = ring_center * 0.95
    for i, wedge in enumerate(wedges):
        ang = math.radians((wedge.theta1 + wedge.theta2) / 2)
        arc_deg = abs(wedge.theta2 - wedge.theta1)
        fs = max(9, min(14, arc_deg * 0.06))
        name = db_counts.index[i]
        val = db_counts.values[i]
        pct = val / total * 100
        label = f"{name}\n{pct:.1f}%  (n={val})"
        ax.text(label_r * math.cos(ang), label_r * math.sin(ang),
                label, ha='center', va='center',
                fontsize=fs, fontweight='bold', color='white',
                linespacing=1.2)
    ax.set_title('B. Proportion of Data Origin', loc='left', fontsize=18, fontweight='bold', pad=30)

def _pnl_c_organizations(ax, df):
    """C. 核心贡献机构 Top10。"""
    org_counts = clean_series(df['Organization_CenterName']).value_counts().head(10).sort_values(ascending=True)
    n = len(org_counts)
    y_pos = range(n)
    ax.barh(y_pos, org_counts.values, height=0.65, color='#4C72B0')
    for i, v in enumerate(org_counts.values):
        ax.text(v + max(org_counts.values)*0.02, i, str(v), va='center', fontsize=11)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(org_counts.index, fontsize=10, linespacing=1.4)
    ax.set_title('C. Top 10 Contributing Organizations', loc='left', fontsize=16, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=14)
    ax.invert_yaxis()

def _pnl_d_tissues(ax, df):
    """D. 研究部位偏好 Top10。"""
    if 'Tissue' in df.columns:
        tissue_counts = clean_series(df['Tissue']).value_counts().head(10).sort_values(ascending=True)
        if not tissue_counts.empty:
            n = len(tissue_counts)
            y_pos = range(n)
            ax.barh(y_pos, tissue_counts.values, height=0.65, color='#4C72B0', edgecolor='white')
            for i, v in enumerate(tissue_counts.values):
                ax.text(v + max(tissue_counts.values)*0.02, i, str(v), va='center', fontsize=11)
            ax.set_yticks(list(y_pos))
            ax.set_yticklabels(tissue_counts.index, fontsize=10)
            ax.invert_yaxis()
    ax.set_title('D. Top Investigated Biological Tissues', loc='left', fontsize=16, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=14)

def _pnl_e_regions(ax, df):
    """E. 采样地理分布 Top10。"""
    if 'Location_Clean' in df.columns:
        loc_counts = clean_series(df['Location_Clean']).value_counts().head(10).sort_values(ascending=True)
        if not loc_counts.empty:
            n = len(loc_counts)
            y_pos = range(n)
            ax.barh(y_pos, loc_counts.values, height=0.65, color='#C44E52')
            for i, v in enumerate(loc_counts.values):
                ax.text(v + max(loc_counts.values)*0.02, i, str(v), va='center', fontsize=11)
            ax.set_yticks(list(y_pos))
            ax.set_yticklabels(loc_counts.index, fontsize=10)
            ax.invert_yaxis()
    ax.set_title('E. Top Sample Collection Regions', loc='left', fontsize=16, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=14)

def _pnl_f_stages(ax, df):
    """F. 发育时期偏好 Top10。"""
    if 'Age_GrowthStage_Clean' in df.columns:
        stage_counts = clean_series(df['Age_GrowthStage_Clean']).value_counts().head(10).sort_values(ascending=True)
        if not stage_counts.empty:
            n = len(stage_counts)
            y_pos = range(n)
            ax.barh(y_pos, stage_counts.values, height=0.65, color='#55A868', edgecolor='white')
            for i, v in enumerate(stage_counts.values):
                ax.text(v + max(stage_counts.values)*0.02, i, str(v), va='center', fontsize=11)
            ax.set_yticks(list(y_pos))
            ax.set_yticklabels(stage_counts.index, fontsize=10, linespacing=1.2)
            ax.invert_yaxis()
    ax.set_title('F. Top Developmental Stages / Ages', loc='left', fontsize=16, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=14)

# 面板注册表：单图文件名（ASCII）+ 单图尺寸
PANELS = [
    ('A_temporal',      _pnl_a_temporal,      (11, 5.5)),
    ('B_database',      _pnl_b_database,      (8, 8)),
    ('C_organizations', _pnl_c_organizations, (11, 5.5)),
    ('D_tissues',       _pnl_d_tissues,       (11, 5.5)),
    ('E_regions',       _pnl_e_regions,       (11, 5.5)),
    ('F_stages',        _pnl_f_stages,        (11, 5.5)),
]

def plot_sci_landscape(csv_path, output_dir="SCI_Figures_Output"):
    """主入口：输出 3x2 合并总览图 + A~F 六张独立图（均 PNG 600dpi + PDF 矢量）。
    独立图文件名: Panel_<key>.png/.pdf（如 Panel_A_temporal.png）。"""
    os.makedirs(output_dir, exist_ok=True)
    print(f"📥 正在读取数据: {csv_path}")
    df = _prep(pd.read_csv(csv_path))

    # ==========================================
    # 2. 合并总画布 (3x2 完美布局)
    #    显示 dpi 150; 保存时 PNG 600 / PDF 矢量 (出版要求 ≥300dpi)
    # ==========================================
    fig, axes = plt.subplots(3, 2, figsize=(16, 20), dpi=150)
    ax_list = axes.flatten()
    for i, (_key, draw, _fs) in enumerate(PANELS):
        print(f"📊 绘制 {chr(65 + i)}: ...")
        draw(ax_list[i], df)
    plt.tight_layout(pad=2.0, h_pad=3.0, w_pad=3.0)

    # 4. 输出总拼图 (PDF 矢量 + PNG 600dpi, 满足期刊 ≥300dpi 要求)
    out_pdf = os.path.join(output_dir, "Combined_Landscape_Full.pdf")
    out_png = os.path.join(output_dir, "Combined_Landscape_Full.png")
    fig.savefig(out_pdf, format='pdf', bbox_inches='tight')
    fig.savefig(out_png, format='png', dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"✅ 合并总览图: {out_png}")

    # 5. A~F 六张独立图（GUI 分开展示 / 单独引用）
    for key, draw, fsize in PANELS:
        sfig, sax = plt.subplots(figsize=fsize, dpi=150)
        draw(sax, df)
        sfig.tight_layout(pad=1.6)
        spdf = os.path.join(output_dir, f"Panel_{key}.pdf")
        spng = os.path.join(output_dir, f"Panel_{key}.png")
        sfig.savefig(spdf, format='pdf', bbox_inches='tight')
        sfig.savefig(spng, format='png', dpi=600, bbox_inches='tight')
        plt.close(sfig)
        print(f"✅ 独立图: {spng}")

    print(f"\n🎉 完美收工！所有图表已存放至文件夹: [ {output_dir}/ ]")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SCI 级宏观组学数据可视化引擎 (精修版)")
    parser.add_argument("-i", "--input", required=True, help="输入的 SRA_GSA_Merged_Final.csv 文件路径")
    parser.add_argument("-o", "--outdir", default="SCI_Figures_Output", help="图表输出的文件夹名称")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"❌ 找不到文件: {args.input}")
    else:
        plot_sci_landscape(args.input, args.outdir)

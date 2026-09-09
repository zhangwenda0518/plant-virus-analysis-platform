# -*- coding: utf-8 -*-
"""
阶段⑦ 可视化与汇总报告：
- plotly 交互图：桑基图(reads→宿主/病毒科属种流向)、旭日图(分类组成)、
  丰度柱状图、SDT identity 热图（全部离线，内嵌本地 plotly.min.js）
- matplotlib 静态图：病毒基因组圈图(pycirclize，含 ORF/GC)、进化树(Bio.PhylO)
- report.html 单文件汇总（图 + 表 + 统计卡片）
"""
import os
import io
import re
import json
import shutil
import base64

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

import plotly.graph_objects as go

from .config import DIRS
from .utils import check_path, safe_open, iter_fasta

# 中文字体 & 离线 plotly.min.js 复制
_JS_COPIED = False


def _plotly_js_path(report_dir):
    """把 plotly.min.js 复制到报告目录（离线渲染），返回文件名。"""
    global _JS_COPIED
    dst = os.path.join(report_dir, 'plotly.min.js')
    if not os.path.isfile(dst):
        import plotly
        src = os.path.join(os.path.dirname(plotly.__file__), 'package_data',
                           'plotly.min.js')
        if not os.path.isfile(src):
            src = os.path.join(os.path.dirname(plotly.__file__), 'offline',
                               'plotly.min.js')
        if os.path.isfile(src):
            shutil.copyfile(check_path(src, must_exist=True),
                            check_path(dst, must_exist=False, in_platform=True))
            _JS_COPIED = True
    return 'plotly.min.js' if os.path.isfile(dst) else None


def _save_fig_png(fig, path):
    """matplotlib 图保存 PNG（写盘走 safe_open）。"""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    with safe_open(path, 'wb') as f:
        f.write(buf.getvalue())


def png_b64(path):
    with open(check_path(path, must_exist=True), 'rb') as f:
        return base64.b64encode(f.read()).decode()


# ------------------------------------------------------------------
# 数据加载
# ------------------------------------------------------------------
def _read_json(p):
    with safe_open(p) as f:
        return json.load(f)


def _read_tsv(p):
    rows = []
    with safe_open(p) as f:
        header = None
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
                continue
            rows.append(dict(zip(header, parts)))
    return rows


# ------------------------------------------------------------------
# plotly 图
# ------------------------------------------------------------------
def _kreport_tree(sample_dir, virus_summary):
    """从 02_virus_screen 的 kreport 构建轻量分类树（避免全量 taxonomy）。"""
    from .kunpeng import KreportTree, parse_kreport
    if virus_summary and virus_summary.get('kreport'):
        p = os.path.join(sample_dir, '02_virus_screen', virus_summary['kreport'])
        if os.path.isfile(p):
            return KreportTree(parse_kreport(p))
    return KreportTree([])


def _limit_per_parent(rows, parent_of, top_n):
    """metabuli 式限流（移植自 9.metabuli renderSankey/limitSpeciesPerGenus）：
    每个 parent 分类单元下按 reads 取 top_n 个 child，其余不进入桑基图，
    避免大数据样品（数千分类单元）节点爆炸不可读。"""
    by_parent = {}
    for r in rows:
        try:
            p = parent_of(int(r['taxid']))
        except (ValueError, TypeError):
            p = None
        by_parent.setdefault(p, []).append(r)
    keep = []
    for rs in by_parent.values():
        rs.sort(key=lambda r: -int(r['reads']))
        keep.extend(rs[:top_n])
    return keep


def fig_sankey(host_stats, virus_rows, virus_total_unclassified, tree):
    """reads 流向桑基图：总reads → 宿主/非宿主 → 病毒界→科→属→种 → 未分类。

    每属仅展示 reads 前 5 的种、每科前 10 的属（低丰度不绘制，聚合层
    仍保留全量数值），保证大样品可读性。
    """
    total = host_stats['total_pairs']
    host = host_stats['dropped_pairs']
    kept = host_stats['kept_pairs']

    labels = [f'总 reads 对\n{total:,}', f'宿主\n{host:,}', f'非宿主\n{kept:,}']
    src, tgt, val = [], [], []

    def add(s, t, v):
        if v > 0:
            src.append(s)
            tgt.append(t)
            val.append(v)

    add(0, 1, host)
    add(0, 2, kept)
    # 病毒分类层级：科(family)->属(genus)->种(species)
    fam_rows = [r for r in virus_rows if r['rank'] == 'family' and int(r['reads']) > 0]
    gen_rows = [r for r in virus_rows if r['rank'] == 'genus' and int(r['reads']) > 0]
    spe_rows = [r for r in virus_rows if r['rank'] == 'species' and int(r['reads']) > 0]
    id2label = {}

    def node_for(taxid, reads):
        if taxid not in id2label:
            nm = tree.name(taxid)
            idx = len(labels)
            labels.append(f"{nm[:28]}\n{reads:,}")
            id2label[taxid] = idx
        return id2label[taxid]

    fam_reads = {int(r['taxid']): int(r['reads']) for r in fam_rows}
    gen_reads = {int(r['taxid']): int(r['reads']) for r in gen_rows}
    spe_reads = {int(r['taxid']): int(r['reads']) for r in spe_rows}

    # metabuli 式限流加强版：科 top12、每科 top8 属、每属 top5 种，
    # 另加全局上限（属 60 / 种 150）；落选科聚合为「其他科」节点。
    top_fams = sorted(fam_reads.items(), key=lambda x: -x[1])[:12]
    kept_fams = {t for t, _ in top_fams}
    other_fam_reads = sum(v for t, v in fam_reads.items() if t not in kept_fams)
    fam_reads = dict(top_fams)
    gen_rows = _limit_per_parent(
        gen_rows, lambda t: tree.get_ancestor_at_rank(t, 'family'), 8)
    gen_rows = sorted(gen_rows, key=lambda r: -int(r['reads']))[:60]
    kept_gen = {int(r['taxid']) for r in gen_rows}
    gen_reads = {t: v for t, v in gen_reads.items()
                 if t in kept_gen
                 and tree.get_ancestor_at_rank(t, 'family') in kept_fams}
    spe_rows = _limit_per_parent(
        spe_rows, lambda t: tree.get_ancestor_at_rank(t, 'genus'), 5)
    spe_rows = sorted(spe_rows, key=lambda r: -int(r['reads']))[:150]
    kept_spe = {int(r['taxid']) for r in spe_rows}
    spe_reads = {t: v for t, v in spe_reads.items()
                 if t in kept_spe
                 and tree.get_ancestor_at_rank(t, 'family') in kept_fams}

    for ftx, fr in fam_reads.items():
        add(2, node_for(ftx, fr), fr)
    if other_fam_reads > 0:
        oi = len(labels)
        labels.append(f'其他科（低丰度）\n{other_fam_reads:,}')
        add(2, oi, other_fam_reads)
    for gtx, gr in gen_reads.items():
        parent = tree.get_ancestor_at_rank(gtx, 'family')
        if parent in fam_reads:
            add(node_for(parent, fam_reads[parent]), node_for(gtx, gr), gr)
    for stx, sr in spe_reads.items():
        parent = tree.get_ancestor_at_rank(stx, 'genus')
        if parent in gen_reads:
            add(node_for(parent, gen_reads[parent]), node_for(stx, sr), sr)
    ui = len(labels)
    labels.append(f'未分类病毒\n{virus_total_unclassified:,}')
    add(2, ui, virus_total_unclassified)

    fig = go.Figure(go.Sankey(
        node=dict(pad=12, thickness=14, line=dict(color='#333', width=0.6),
                  label=labels),
        link=dict(source=src, target=tgt, value=val,
                  color='rgba(66,133,244,0.35)')))
    fig.update_layout(title='Reads 分类流向（宿主去除 + 病毒分类；'
                            '每科展示前 10 属 / 每属前 5 种）',
                      font=dict(size=11), height=560)
    return fig


def fig_sunburst(virus_rows, tree):
    """病毒分类组成旭日图：科→属→种（每属限 reads 前 5 的种）。"""
    fam_rows = [r for r in virus_rows if r['rank'] == 'family' and int(r['reads']) > 0]
    gen_rows = [r for r in virus_rows if r['rank'] == 'genus' and int(r['reads']) > 0]
    spe_rows = [r for r in virus_rows if r['rank'] == 'species' and int(r['reads']) > 0]
    spe_rows = _limit_per_parent(
        spe_rows, lambda t: tree.get_ancestor_at_rank(t, 'genus'), 5)
    ids, labels, parents, values = [], [], [], []

    def short(n):
        return (n[:24] + '..') if len(n) > 26 else n

    for r in fam_rows:
        ids.append(f"f{r['taxid']}")
        labels.append(short(r['name']))
        parents.append('')
        values.append(int(r['reads']))
    for r in gen_rows:
        pf = tree.get_ancestor_at_rank(int(r['taxid']), 'family')
        parents.append(f"f{pf}" if pf else '')
        ids.append(f"g{r['taxid']}")
        labels.append(short(r['name']))
        values.append(0)
    for r in spe_rows:
        pg = tree.get_ancestor_at_rank(int(r['taxid']), 'genus')
        parents.append(f"g{pg}" if pg else '')
        ids.append(f"s{r['taxid']}")
        labels.append(short(r['name']))
        values.append(int(r['reads']))
    if not ids:
        return None
    fig = go.Figure(go.Sunburst(ids=ids, labels=labels, parents=parents,
                                values=values, branchvalues='remainder',
                                maxdepth=3))
    fig.update_layout(title='病毒分类组成（内→外：科 → 属 → 种）', height=560)
    return fig


def fig_host_sankey(host_rows):
    """宿主预测桑基图：病毒科 → 宿主类别（contigs 数为流量）。"""
    fams, links = {}, []
    cats = {}
    for r in host_rows:
        fam = r.get('family') or '（未定科）'
        cat = r.get('final_host') or 'Unknown'
        fams.setdefault(fam, len(fams))
        cats[cat] = cats.get(cat, 0) + 1
        links.append((fam, cat))
    agg = {}
    for fam, cat in links:
        agg[(fam, cat)] = agg.get((fam, cat), 0) + 1
    if not fams:
        return None
    labels = list(fams.keys()) + list(cats.keys())
    idx = {k: len(fams) + i for i, k in enumerate(cats)}
    fig = go.Figure(go.Sankey(
        node=dict(label=labels, pad=12, thickness=16,
                  color=['#2c7a4b'] * len(fams) + ['#3b82c4'] * len(cats)),
        link=dict(
            source=[fams[f] for (f, _c) in agg],
            target=[idx[c] for (_f, c) in agg],
            value=[v for v in agg.values()])))
    fig.update_layout(title=f'病毒科 → 宿主类别（{sum(cats.values())} 条 contigs）',
                      height=520)
    return fig


def fig_host_sunburst(host_rows):
    """宿主预测旭日图：病毒科 → 属 → 种（contigs 数），外圈颜色随宿主。"""
    ids, labels, parents, values = [], [], [], []

    def short(n):
        return (n[:24] + '..') if len(n) > 26 else n

    fam_cnt, gen_cnt, sp_cnt = {}, {}, {}
    for r in host_rows:
        fam = r.get('family') or '（未定科）'
        gen = r.get('genus') or '（未定属）'
        sp = r.get('species') or r.get('contig') or '（未定种）'
        fam_cnt[fam] = fam_cnt.get(fam, 0) + 1
        gen_cnt[(fam, gen)] = gen_cnt.get((fam, gen), 0) + 1
        sp_cnt[(fam, gen, sp)] = sp_cnt.get((fam, gen, sp), 0) + 1
    if not fam_cnt:
        return None
    for fam, c in fam_cnt.items():
        ids.append(f'f|{fam}'); labels.append(short(fam))
        parents.append(''); values.append(c)
    for (fam, gen), c in gen_cnt.items():
        ids.append(f'g|{fam}|{gen}'); labels.append(short(gen))
        parents.append(f'f|{fam}'); values.append(c)
    for (fam, gen, sp), c in sp_cnt.items():
        ids.append(f's|{fam}|{gen}|{sp}'); labels.append(short(sp))
        parents.append(f'g|{fam}|{gen}'); values.append(c)
    fig = go.Figure(go.Sunburst(ids=ids, labels=labels, parents=parents,
                                values=values, branchvalues='remainder',
                                maxdepth=3))
    fig.update_layout(title='病毒 ICTV 分类旭日图（内→外：科 → 属 → 种）',
                      height=560)
    return fig


def fig_abundance_bar(virus_rows, top_n=15):
    """top 物种 reads 条形图。"""
    spe = sorted([r for r in virus_rows if r['rank'] == 'species'],
                 key=lambda r: -int(r['reads']))[:top_n]
    if not spe:
        return None
    spe = spe[::-1]
    fig = go.Figure(go.Bar(
        x=[int(r['reads']) for r in spe],
        y=[r['name'][:34] for r in spe],
        orientation='h',
        text=[f"{int(r['reads']):,} ({r.get('percent(%)', r.get('percent', ''))}%)" for r in spe],
        textposition='auto',
        marker_color='#2ca02c'))
    fig.update_layout(title=f'Top {len(spe)} 病毒物种 reads 数',
                      height=420, margin=dict(l=220))
    return fig


def fig_sdt_heatmap(names, mat):
    """SDT identity 热图（下三角）。"""
    n = len(names)
    z = [[mat[i][j] if j <= i else None for j in range(n)] for i in range(n)]
    fig = go.Figure(go.Heatmap(
        z=z, x=names, y=names, colorscale='RdYlGn',
        zmin=60, zmax=100,
        text=[[f"{mat[i][j]:.1f}%" if j <= i else '' for j in range(n)]
              for i in range(n)],
        texttemplate='%{text}', textfont=dict(size=8)))
    fig.update_layout(title='SDT 全长成对 identity 矩阵 (%)', height=620,
                      xaxis=dict(tickangle=-45, tickfont=dict(size=8)),
                      yaxis=dict(tickfont=dict(size=8)))
    return fig


# ------------------------------------------------------------------
# matplotlib 静态图
# ------------------------------------------------------------------
def plot_tree(nwk_path, out_png, title=''):
    """进化树渲染（Bio.Phylo）。"""
    from Bio import Phylo
    tree = Phylo.read(check_path(nwk_path, must_exist=True), 'newick')
    fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * tree.count_terminals())))
    Phylo.draw(tree, axes=ax, do_show=False,
               show_confidence=True, label_func=lambda c: None)
    # 叶子标签
    for tip in tree.get_terminals():
        pass
    ax.set_title(title or '系统进化树')
    ax.axis('off')
    _save_fig_png(fig, out_png)


def plot_tree_labeled(nwk_path, out_png, title=''):
    """带叶标签的进化树。"""
    from Bio import Phylo
    tree = Phylo.read(check_path(nwk_path, must_exist=True), 'newick')
    n_tips = tree.count_terminals()
    fig, ax = plt.subplots(figsize=(10, max(3.5, 0.32 * n_tips)))
    Phylo.draw(tree, axes=ax, do_show=False,
               label_func=lambda c: (c.name or '')[:40] if c.is_terminal() else '')
    ax.set_title(title or '系统进化树（FastTree/IQ-TREE）')
    _save_fig_png(fig, out_png)


def plot_genome_circos(contig_id, seq, bed_records, out_png, blast_cov=None):
    """病毒 contig 基因组圈图（pycirclize）：ORF 弧 + GC 含量曲线。"""
    try:
        from pycirclize import Circos
        sectors = {contig_id[:20]: len(seq)}
        circos = Circos(sectors, space=8)
        sector = circos.sectors[0]
        sector.text(f"{contig_id[:20]} ({len(seq):,} bp)",
                    r=105, deg=90, size=11)

        # ORF 轨道
        track = sector.add_track((72, 100))
        track.axis(fc='#EEEEEE')
        strand_colors = {1: '#E64B35', -1: '#4DBBD5'}
        for rec in bed_records:
            start, end = int(rec['start']), int(rec['end'])
            strand = -1 if rec.get('strand', '+') == '-' else 1
            track.rect(start, end, r_lim=(84, 100) if strand > 0 else (72, 88),
                       fc=strand_colors[strand], ec='none', alpha=0.95)
        # GC 含量曲线（窗口）
        track2 = sector.add_track((58, 70))
        track2.axis(fc='#FAFAFA')
        win = max(len(seq) // 360, 50)
        xs, ys = [], []
        for i in range(0, len(seq) - win, win):
            sub = seq[i:i + win].upper()
            gc = (sub.count('G') + sub.count('C')) / max(len(sub), 1)
            xs.append(i + win // 2)
            ys.append(gc)
        if xs:
            track2.line(xs, ys, color='#3C5488', lw=1.2)
            track2.yticks([0, 0.5, 1.0], labels=['0', '50%', '100%'])
        fig = circos.plotfig()
        _save_fig_png(fig, out_png)
        return True
    except Exception:
        return plot_genome_linear(contig_id, seq, bed_records, out_png)


def plot_genome_linear(contig_id, seq, bed_records, out_png):
    """圈图失败时的线性基因组图 fallback。"""
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 3.6), gridspec_kw={'height_ratios': [3, 1]})
    L = len(seq)
    ax1.set_xlim(0, L)
    for rec in bed_records:
        start, end = int(rec['start']), int(rec['end'])
        strand = -1 if rec.get('strand', '+') == '-' else 1
        color = '#E64B35' if strand > 0 else '#4DBBD5'
        ax1.arrow(start, 0.55 if strand > 0 else 0.25,
                  end - start, 0, width=0.16, head_width=0.3,
                  head_length=min((end - start) * 0.2, L * 0.02),
                  fc=color, ec='none', length_includes_head=True)
    ax1.set_yticks([0.4, 0.8])
    ax1.set_yticklabels(['-', '+'])
    ax1.set_title(f"{contig_id} ({L:,} bp) ORF 图谱")
    win = max(L // 400, 50)
    xs, ys = [], []
    for i in range(0, L - win, win):
        sub = seq[i:i + win].upper()
        xs.append(i)
        ys.append((sub.count('G') + sub.count('C')) / max(len(sub), 1) * 100)
    ax2.fill_between(xs, ys, color='#3C5488', alpha=0.6)
    ax2.set_ylabel('GC%')
    ax2.set_xlim(0, L)
    _save_fig_png(fig, out_png)
    return True


# ------------------------------------------------------------------
# 主报告
# ------------------------------------------------------------------
def build_report(sample_dir, logger=None):
    """生成 07_report/report.html。返回其路径。"""
    sample_dir = check_path(sample_dir, must_exist=True, in_platform=True)
    out_dir = check_path(os.path.join(sample_dir, '07_report'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)

    host_stats = None
    virus_summary = None
    asm = orf_sum = phylo_sum = primer_sum = None

    p = os.path.join(sample_dir, '01_host_removal', 'stats.json')
    if os.path.isfile(p):
        host_stats = _read_json(p)
    p = os.path.join(sample_dir, '02_virus_screen', 'summary.json')
    if os.path.isfile(p):
        virus_summary = _read_json(p)
    p = os.path.join(sample_dir, '03_assembly', 'summary.json')
    if os.path.isfile(p):
        asm = _read_json(p)
    p = os.path.join(sample_dir, '04_orf', 'summary.json')
    if os.path.isfile(p):
        orf_sum = _read_json(p)
    p = os.path.join(sample_dir, '05_phylo', 'summary.json')
    if os.path.isfile(p):
        phylo_sum = _read_json(p)
    p = os.path.join(sample_dir, '06_primer', 'summary.json')
    if os.path.isfile(p):
        primer_sum = _read_json(p)

    figures = []          # (title, html_div or img_tag)

    # 桑基 + 旭日 + 柱状
    if virus_summary:
        tsv_rows = _read_tsv(os.path.join(sample_dir, '02_virus_screen',
                                          'virus_summary.tsv'))
        tree = _kreport_tree(sample_dir, virus_summary)
        if host_stats:
            fig = fig_sankey(host_stats, tsv_rows,
                             virus_summary.get('unclassified_records', 0), tree)
            figures.append(('Reads 分类流向桑基图', fig.to_html(
                include_plotlyjs=False, full_html=False, div_id='sankey')))
        sb = fig_sunburst(tsv_rows, tree)
        if sb:
            figures.append(('病毒分类组成旭日图', sb.to_html(
                include_plotlyjs=False, full_html=False, div_id='sunburst')))
        bar = fig_abundance_bar(tsv_rows)
        if bar:
            figures.append(('病毒丰度柱状图', bar.to_html(
                include_plotlyjs=False, full_html=False, div_id='bar')))

    # 宿主预测（ICTV 级联）：桑基 + 旭日
    p_host = os.path.join(sample_dir, '08_host_analysis', 'summary.json')
    if os.path.isfile(p_host):
        try:
            host_sum = _read_json(p_host)
            pred_tsv = os.path.join(sample_dir, '08_host_analysis',
                                    host_sum.get('tsv', 'host_prediction.tsv'))
            host_rows = _read_tsv(pred_tsv)
            sk = fig_host_sankey(host_rows)
            if sk:
                figures.append(('宿主预测桑基图（病毒科 → 宿主类别）',
                                sk.to_html(include_plotlyjs=False,
                                           full_html=False,
                                           div_id='host_sankey')))
            sb = fig_host_sunburst(host_rows)
            if sb:
                figures.append(('宿主预测旭日图（科 → 属 → 种）',
                                sb.to_html(include_plotlyjs=False,
                                           full_html=False,
                                           div_id='host_sunburst')))
            cats = host_sum.get('categories', {})
            if cats:
                cat_rows = '\n'.join(f'<tr><td>{k}</td><td>{v}</td></tr>'
                                     for k, v in cats.items())
                figures.append(('宿主类别汇总',
                                f'<table class="ptable" style="margin:auto">'
                                f'<tr><th>宿主类别</th><th>contigs</th></tr>'
                                f'{cat_rows}</table>'))
        except Exception as e:
            if logger:
                logger.log(f"宿主分析图表失败: {e}", "WARN")

    # ⑥b ORF 功能注释：类别/科分布 + Top 注释表
    p_orfa = os.path.join(sample_dir, '04b_orf_annot', 'summary.json')
    if os.path.isfile(p_orfa):
        try:
            orfa = _read_json(p_orfa)
            from .orf_annot import fig_category_bar, fig_family_bar
            cb = fig_category_bar(orfa)
            if cb:
                figures.append(('ORF 功能类别分布',
                                cb.to_html(include_plotlyjs=False,
                                           full_html=False,
                                           div_id='orfa_cat')))
            fb = fig_family_bar(orfa)
            if fb:
                figures.append(('ORF 病毒科分布',
                                fb.to_html(include_plotlyjs=False,
                                           full_html=False,
                                           div_id='orfa_fam')))
            prev = orfa.get('rows_preview') or []
            if prev:
                trs = '\n'.join(
                    f'<tr><td>{r.get("orf_id", "")}</td>'
                    f'<td>{r.get("product", "")}</td>'
                    f'<td>{r.get("organism", "")}</td>'
                    f'<td>{r.get("family", "")}</td>'
                    f'<td>{r.get("category", "")}</td>'
                    f'<td>{r.get("pident", "")}</td></tr>'
                    for r in prev[:15])
                figures.append((
                    f"Top ORF 功能注释（{orfa.get('n_annotated', 0)}/"
                    f"{orfa.get('n_orfs', 0)} 个 ORF 有效注释，"
                    f"引擎 {orfa.get('engine', '?')}，完整表见 orf_annotation.tsv）",
                    f'<table class="ptable" style="margin:auto">'
                    f'<tr><th>ORF</th><th>产物</th><th>物种</th><th>科</th>'
                    f'<th>类别</th><th>identity%</th></tr>{trs}</table>'))
        except Exception as e:
            if logger:
                logger.log(f"ORF 功能注释图表失败: {e}", "WARN")

    # 各组：树 + SDT 热图 + 圈图
    if phylo_sum and not phylo_sum.get('skipped'):
        from .phylo import pairwise_identity_matrix
        for g in phylo_sum.get('groups', []):
            gdir = check_path(os.path.join(sample_dir, '05_phylo', g['dir']),
                              must_exist=True)
            if g.get('tree'):
                nwk = os.path.join(gdir, g['tree'])
                if os.path.isfile(nwk):
                    png = os.path.join(gdir, 'tree.png')
                    try:
                        plot_tree_labeled(nwk, png, title=f"{g['group']} 进化树")
                        figures.append((f"{g['group']} 系统进化树（与最近参考 identity "
                                        f"{g.get('best_identity_to_ref', 0)}%）",
                                        f'<img style="max-width:100%" src="data:image/png;'
                                        f'base64,{png_b64(png)}">'))
                    except Exception as e:
                        if logger:
                            logger.log(f"树渲染失败 {g['group']}: {e}", "WARN")
            aln = os.path.join(gdir, 'aln.fasta')
            if os.path.isfile(aln):
                try:
                    names, mat = pairwise_identity_matrix(aln)
                    hm = fig_sdt_heatmap(names, mat)
                    figures.append((f"{g['group']} SDT identity 矩阵", hm.to_html(
                        include_plotlyjs=False, full_html=False,
                        div_id=f'heat_{g["group"]}')))
                except Exception as e:
                    if logger:
                        logger.log(f"SDT 热图失败 {g['group']}: {e}", "WARN")

    # 病毒基因组图（gbdraw，阶段⑨产物 SVG 内嵌）
    p_gb = os.path.join(sample_dir, '09_genome_plots', 'summary.json')
    if os.path.isfile(p_gb):
        try:
            gb = _read_json(p_gb)
            n_svg = 0
            for rel in gb.get('plots', []):
                svg_p = os.path.join(sample_dir, '09_genome_plots', rel)
                if not os.path.isfile(svg_p):
                    continue
                with safe_open(svg_p) as f:
                    svg = f.read()
                # 去掉 XML 声明并让 SVG 自适应容器宽度
                svg = re.sub(r'<\?xml[^>]*\?>', '', svg, count=1)
                svg = re.sub(r'(<svg[^>]*?)\swidth="[^"]*"\sheight="[^"]*"',
                             r'\1 width="100%"', svg, count=1)
                base = os.path.basename(rel)
                figures.append((f"病毒基因组图 — {base.replace('.svg', '')}",
                                f'<div style="overflow:auto">{svg}</div>'))
                n_svg += 1
            if not n_svg and logger:
                logger.log("09_genome_plots 无可用 SVG，跳过基因组图小节", "WARN")
        except Exception as e:
            if logger:
                logger.log(f"gbdraw 嵌入失败: {e}", "WARN")

    # 病毒 contig 圈图（每组最长 contig）
    if asm:
        contigs_fa = os.path.join(sample_dir, '03_assembly', 'contigs.filtered.fasta')
        bed_file = os.path.join(sample_dir, '04_orf', 'orfipy.bed')
        bed_records = []
        if os.path.isfile(bed_file):
            # orfipy BED 无表头：标准 6+ 列（chrom start end name score strand ...）
            with safe_open(bed_file) as f:
                for line in f:
                    if line.startswith('#') or line.startswith('track') or not line.strip():
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) >= 6:
                        bed_records.append({'chrom': parts[0],
                                            'start': parts[1], 'end': parts[2],
                                            'strand': parts[5]})
        if os.path.isfile(contigs_fa) and bed_records:
            contigs = sorted(iter_fasta(contigs_fa), key=lambda x: -len(x[1]))[:5]
            for cid, s in contigs:
                recs = [r for r in bed_records if r['chrom'] == cid]
                if not recs:
                    continue
                png = os.path.join(out_dir, f'genome_{cid[:30]}.png')
                try:
                    if plot_genome_circos(cid, s, recs, png):
                        figures.append((f"病毒 contig {cid} 基因组圈图 (ORF + GC)",
                                        f'<img style="max-width:100%" src="data:image/png;'
                                        f'base64,{png_b64(png)}">'))
                except Exception as e:
                    if logger:
                        logger.log(f"圈图失败 {cid}: {e}", "WARN")

    # HTML 组装
    js = _plotly_js_path(out_dir)
    js_tag = (f'<script src="{js}"></script>' if js
              else '<script>window.Plotly||document.write("plotly.min.js 缺失，交互图不可用")</script>')

    cards = []
    if host_stats:
        cards.append(('总 reads 对', f"{host_stats['total_pairs']:,}"))
        cards.append(('宿主占比', f"{host_stats['host_ratio'] * 100:.2f}%"))
        cards.append(('非宿主 reads 对', f"{host_stats['kept_pairs']:,}"))
    if virus_summary:
        cards.append(('检出病毒物种', str(len(virus_summary['species_detected']))))
        cards.append(('病毒分类 reads', f"{virus_summary['total_classified_records']:,}"))
    if asm:
        cards.append(('contigs(≥min)', f"{asm['contigs']:,}"))
        cards.append(('病毒 contigs', str(len(asm['viral_contigs']))))
    if orf_sum:
        pyro = orf_sum.get('pyrodigal') or {}
        cards.append(('预测基因(pyrodigal)', str(pyro.get('genes', '-'))))
    if primer_sum:
        cards.append(('设计引物对', str(primer_sum.get('n_primers', 0))))

    # 表格：病毒物种 + 病毒 contigs + 引物
    tables = []
    if virus_summary and virus_summary['species_detected']:
        rows = ''.join(
            f"<tr><td>{s['name']}</td><td>{s['taxid']}</td><td>{s['reads']:,}</td>"
            f"<td>{s['percent']}%</td></tr>"
            for s in virus_summary['species_detected'][:30])
        tables.append(('检出病毒物种', 
            '<table><tr><th>物种</th><th>taxid</th><th>reads</th><th>占比%</th></tr>'
            + rows + '</table>'))
    vtsv = os.path.join(sample_dir, '03_assembly', 'virus_contigs.tsv')
    if os.path.isfile(vtsv):
        rows = _read_tsv(vtsv)
        trs = ''.join(
            f"<tr><td>{r['contig']}</td><td>{r['length']}</td>"
            f"<td>{r.get('blast_species') or r.get('kunpeng_species', '')}</td>"
            f"<td>{r.get('blast_identity(%)', '')}</td>"
            f"<td>{r.get('blast_top_hit', '')}</td></tr>" for r in rows[:50])
        tables.append(('病毒 contigs 注释',
            '<table><tr><th>contig</th><th>长度bp</th><th>物种</th>'
            f'<th>BLAST identity%</th><th>最近参考</th></tr>{trs}</table>'))
    ptsv = os.path.join(sample_dir, '06_primer', 'primers.tsv')
    if os.path.isfile(ptsv):
        rows = _read_tsv(ptsv)
        trs = ''.join(
            f"<tr><td>{r['pair']}</td><td>{r['target']}</td>"
            f"<td><code>{r['F_primer(5\'-3\')']}</code></td>"
            f"<td>{r['F_Tm']}</td><td><code>{r['R_primer(5\'-3\')']}</code></td>"
            f"<td>{r['R_Tm']}</td><td>{r['product(bp)']}</td></tr>"
            for r in rows[:60])
        tables.append(('引物设计结果',
            '<table><tr><th>引物对</th><th>靶标</th><th>正向引物</th><th>Tm</th>'
            f'<th>反向引物</th><th>Tm</th><th>产物bp</th></tr>{trs}</table>'))

    import time as _t
    html = f"""<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8">
<title>植物病毒分析报告 - {os.path.basename(sample_dir)}</title>
{js_tag}
<style>
 body{{font-family:"Microsoft YaHei",sans-serif;margin:0;background:#f5f6f8;color:#222}}
 .wrap{{max-width:1200px;margin:0 auto;padding:24px}}
 h1{{font-size:26px;border-bottom:3px solid #2c7a4b;padding-bottom:10px}}
 h2{{font-size:19px;margin-top:34px;border-left:5px solid #2c7a4b;padding-left:10px}}
 .cards{{display:flex;flex-wrap:wrap;gap:12px;margin:18px 0}}
 .card{{background:#fff;border-radius:10px;padding:14px 20px;box-shadow:0 1px 4px rgba(0,0,0,.08);min-width:150px}}
 .card .v{{font-size:22px;font-weight:700;color:#2c7a4b}}
 .card .k{{font-size:12px;color:#777}}
 .fig{{background:#fff;border-radius:10px;padding:14px;margin:14px 0;box-shadow:0 1px 4px rgba(0,0,0,.08)}}
 table{{border-collapse:collapse;width:100%;font-size:13px;background:#fff}}
 th,td{{border:1px solid #ddd;padding:6px 9px;text-align:left}}
 th{{background:#eef5ef}}
 code{{background:#f0f0f0;padding:1px 5px;border-radius:3px;font-size:12px}}
 footer{{color:#999;font-size:12px;margin:30px 0;text-align:center}}
</style></head><body><div class="wrap">
<h1>植物病毒分析报告 — {os.path.basename(sample_dir)}</h1>
<div style="color:#888">生成时间：{_t.strftime('%Y-%m-%d %H:%M:%S')}</div>
<div class="cards">{''.join(f'<div class="card"><div class="v">{v}</div><div class="k">{k}</div></div>' for k, v in cards)}</div>
{''.join(f'<h2>{t}</h2><div class="fig">{b}</div>' for t, b in figures)}
{''.join(f'<h2>{t}</h2><div style="overflow:auto">{b}</div>' for t, b in tables)}
<footer>植物病毒分析平台 · kunpeng / SPAdes / BLAST / MAFFT / FastTree / primer3</footer>
</div></body></html>"""

    report = os.path.join(out_dir, 'report.html')
    with safe_open(report, 'wt') as f:
        f.write(html)
    if logger:
        logger.log(f"报告生成: {report}")
    return report

def build_contig_report(run_dir, logger=None):
    """工具运行分类报告（1:1 复刻 <DEMO-IP>/metabuli 结果页）。

    contigs_* 运行: kreport + virus_classification.tsv（contig 谱系表）。
    identify_* 运行: kreport + viral_ids.tsv（测序 reads 分类）。
    页面结构/交互与 32-server 一致：Taxonomy Sankey（每属 top5 开关）→
    Classification Table（缩进树形表）→ Virus Sequence Classification
    （可排序 + Analyze 按钮 + 近完整绿色）→ Krona（旭日替代，显示未分类
    开关）→ 底部下载按钮。宿主预测结果（若有）并入 contig 表为 host 列
    并绘制宿主统计图。写入 <run_dir>/report.html。
    """
    import glob
    import csv as _csv
    import time as _time
    from .kunpeng import parse_kreport
    from .utils import count_fasta_seqs

    run_dir = check_path(run_dir, must_exist=True, in_platform=True)
    krs = sorted(glob.glob(os.path.join(run_dir, '*.kreport2')) +
                 glob.glob(os.path.join(run_dir, '**', '*.kreport2'),
                           recursive=True))
    if not krs:
        raise RuntimeError('运行目录中没有 kreport 分类报告')
    krows = parse_kreport(krs[0])
    _plotly_js_path(run_dir)
    tb_ok, tb_unc = build_taxburst(run_dir, krows, logger=logger)

    # 报告行（32-server /metabuli/report 同构：proportion/count/rank/taxid/name）
    report_rows = [{'proportion': r['percent'], 'count': int(r['frags']),
                    'rank': r['rank'], 'taxid': r['taxid'],
                    'name': r['name'], 'depth': r['depth']}
                   for r in krows]

    tsv = os.path.join(run_dir, 'virus_classification.tsv')
    ids_tsv = os.path.join(run_dir, 'viral_ids.tsv')
    mode = 'contig' if os.path.isfile(tsv) else 'reads'

    # ---- 宿主预测（仅 contig 运行，可选） ----
    hp_tsv = os.path.join(run_dir, '08_host_analysis', 'host_prediction.tsv')
    host_map = {}
    host_cats = {}
    if mode == 'contig' and os.path.isfile(hp_tsv):
        with safe_open(hp_tsv) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                c = (r.get('contig') or '').strip()
                h = (r.get('final_host') or '').strip()
                conf = (r.get('confidence_level') or '').strip()
                if c:
                    host_map[c] = h + (f' ({conf})' if conf else '')
                if h:
                    host_cats[h] = host_cats.get(h, 0) + 1

    vc_rows, ids_rows, n_total, n_virus = [], [], 0, 0
    if mode == 'contig':
        with safe_open(tsv) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                r['host'] = host_map.get((r.get('contig') or '').strip(), '')
                vc_rows.append(r)
        cfa = os.path.join(run_dir, 'contigs.filtered.fasta')
        if os.path.isfile(cfa):
            from .utils import count_fasta_seqs
            n_total = count_fasta_seqs(cfa)
        n_virus = len(vc_rows)
    else:
        if os.path.isfile(ids_tsv):
            with safe_open(ids_tsv) as f:
                for r in _csv.DictReader(f, delimiter='\t'):
                    ids_rows.append(r)
        n_total = len(ids_rows)
        n_virus = len(ids_rows)

    fastas = sorted(glob.glob(os.path.join(run_dir, 'viral_sequences*.fasta')) +
                    glob.glob(os.path.join(run_dir, 'viral_sequences*.fa')))
    dl = []
    if mode == 'contig':
        dl = [('viral_contigs.fasta?dl=1', 'Virus Sequences (FASTA)'),
              ('virus_classification.tsv?dl=1', 'Virus Classification (TSV)'),
              ('contigs.filtered.fasta?dl=1', 'All Contigs (FASTA)'),
              (os.path.basename(krs[0]) + '?dl=1', 'Report (TSV)')]
        if os.path.isfile(hp_tsv):
            dl.append(('08_host_analysis/host_prediction.tsv?dl=1',
                       'Host Prediction (TSV)'))
    else:
        for fa in fastas:
            dl.append((os.path.basename(fa) + '?dl=1',
                       'Viral Sequences (FASTA) · ' + os.path.basename(fa)))
        if os.path.isfile(ids_tsv):
            dl.append(('viral_ids.tsv?dl=1', 'Classified IDs (TSV)'))
        dl.append((os.path.basename(krs[0]) + '?dl=1', 'Report (kreport TSV)'))
    # 配色：首颗主输出=填充绿；其余按类型描边（分类/序列=绿，历史=蓝，报告=灰）
    def _dl_class(label, idx):
        if idx == 0:
            return 'dl primary'
        low = label.lower()
        if 'report' in low:
            return 'dl gray'
        if 'host' in low or 'classified' in low:
            return 'dl blue'
        return 'dl'
    dl_items = [(href, label, _dl_class(label, i))
                for i, (href, label) in enumerate(dl)]
    dl_html = ''.join(
        f'<a class="{cls}" href="{href}" download>⬇ {label}</a>'
        for href, label, cls in dl_items)
    # Krona（旭日图交互 HTML）：有 taxburst 时追加下载按钮
    if tb_ok and os.path.isfile(os.path.join(run_dir, 'taxburst.html')):
        dl_html += ('<a class="dl" href="taxburst.html?dl=1" download '
                    'title="交互式旭日图（Krona 等价物）">⬇ Download Krona</a>')

    n_contig = len(vc_rows)
    stamp = _time.strftime('%Y-%m-%d %H:%M')
    report_json = json.dumps(report_rows, ensure_ascii=False)
    vc_json = json.dumps(vc_rows, ensure_ascii=False)
    host_json = json.dumps(
        [{'cat': k, 'n': v} for k, v in
         sorted(host_cats.items(), key=lambda x: -x[1])],
        ensure_ascii=False)
    ids_json = json.dumps(ids_rows, ensure_ascii=False)

    if host_map:
        host_section = (f'<h4>Host Prediction (宿主预测统计 · '
                        f'{sum(host_cats.values())} 条 contigs)</h4>'
                        f'<div id="hostBar" style="height:420px"></div>')
    else:
        host_section = ''

    head = f"""<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">
<title>Classification Report</title>
<style>
body{{font-family:'Microsoft YaHei',sans-serif;margin:24px auto;max-width:1200px;
     padding:0 16px;color:#1f2d3d;background:#fff}}
h1{{font-size:20px;margin:0 0 4px}}
h4{{font-size:14px;color:#1a5276;margin:16px 0 8px}}
.ok{{color:#27ae60;font-weight:600;margin-bottom:12px}}
.hint{{color:#789;font-size:12px}}
table.tb{{font-size:12px;width:100%;border-collapse:collapse}}
table.tb th{{background:#f5f5f5;padding:7px;position:sticky;top:0;
  text-align:left;white-space:nowrap;border-bottom:2px solid #ddd}}
table.tb td{{padding:6px 7px;border-bottom:1px solid #f0f0f0;
  vertical-align:top}}
table.tb td.num{{text-align:right;font-family:'Consolas','Monaco',monospace;
  font-variant-numeric:tabular-nums}}
table.tb th.num{{text-align:right}}
table.tb tr:hover{{background:#f8fbff}}
.scrollbox{{max-height:440px;overflow:auto;border:1px solid #eee;
  border-radius:4px}}
.scrollbox.vc{{max-height:380px}}
label.ck{{font-size:12px;color:#555;display:inline-block;margin-bottom:6px;
  cursor:pointer}}
.dl{{display:inline-flex;align-items:center;gap:5px;padding:6px 14px;margin:3px 6px 3px 0;
  border:1px solid #1d7a4f;border-radius:6px;color:#1d7a4f;text-decoration:none;
  font-size:12.5px;font-weight:600;white-space:nowrap}}
.dl:hover{{background:#1d7a4f;color:#fff}}
.dl.primary{{background:#1d7a4f;color:#fff;border-color:#166844}}
.dl.primary:hover{{background:#166844}}
.dl.blue{{border-color:#2e86c1;color:#2e86c1}}
.dl.blue:hover{{background:#2e86c1;color:#fff}}
.dl.gray{{border-color:#95a5a6;color:#7f8c8d}}
.dl.gray:hover{{background:#7f8c8d;color:#fff}}
.dl-row{{display:flex;flex-wrap:wrap;align-items:center}}
th.sortable{{cursor:pointer;user-select:none}}
th.sortable:hover{{background:#e8e8e8}}
</style></head><body>
<h1>Virus Sequence Analysis · 分类报告</h1>
<div class="ok">Classification complete &middot; 报告生成于 {stamp}</div>
<p class="hint">run: {os.path.basename(run_dir)} ·
   {'contig 分类（工具④）' if mode == 'contig' else '测序 reads 分类（工具②）'}</p>
<h4>Taxonomy Sankey</h4>
<label class="ck"><input type="checkbox" id="skTop5" checked
  onchange="renderSankey()"> Only top 5 species per genus</label>
<div id="sk" style="width:100%;height:500px"></div>
<h4>Classification Table</h4>
<div class="scrollbox">
<table class="tb"><thead><tr><th>Rank</th><th>Taxon</th><th>TaxID</th>
<th>%</th><th>Reads</th></tr></thead><tbody id="ctbody"></tbody></table></div>
"""

    if tb_ok:
        krona_section = (
            '<label class="ck"><input type="checkbox" id="kronaUnc" checked'
            ' onchange="toggleKronaHtml()"> Show unclassified</label>'
            '<iframe id="kronaFrame" src="taxburst.html"'
            ' style="width:100%;height:500px;border:1px solid #eee;'
            'border-radius:4px;background:#fff"></iframe>'
            '<script>function toggleKronaHtml(){var f=document.getElementById'
            "('kronaFrame');f.src=f.src.includes('taxburst_unc.html')?"
            "'taxburst.html':'taxburst_unc.html';}</script>")
    else:
        krona_section = (
            '<label class="ck"><input type="checkbox" id="kronaUnc"'
            ' onchange="renderSun()"> Show unclassified</label>'
            '<div id="krona" style="width:100%;height:500px"></div>')

    body_mid = f"""

<h4>Virus Sequence Classification</h4>
<div class="hint" style="margin-bottom:8px">Total <b>{n_total}</b> {'contigs' if mode == 'contig' else 'sequences'}, <b style="color:#27ae60">{n_virus}</b> classified as Virus</div>
<div class="scrollbox vc">
<table class="tb" id="vtable"><thead><tr id="vthead"></tr></thead>
<tbody id="vtbody"></tbody></table></div>
{host_section}
{krona_section}
<h4>Downloads</h4>
<div class="dl-row">{dl_html or '<span class="hint">（无产物文件）</span>'}</div>
<script src="plotly.min.js"></script>
<script>
var MODE = {json.dumps(mode)};
var REPORT = {report_json};
var VC = {vc_json};
var IDS = {ids_json};
var HOST = {host_json};
var N_TOTAL = {n_total};
var N_VIRUS = {n_virus};
"""

    tail = r"""
function escapeHtml(t) {
  return String(t == null ? '' : t)
    .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}
function limitSpeciesPerGenus(rows, topN) {
  var st = [], byGenus = {};
  rows.forEach(function(p, idx) {
    while (st.length > 0 && st[st.length - 1].d >= p.depth) st.pop();
    if (p.rank === 'species') {
      var g = '(no genus)';
      for (var k = st.length - 1; k >= 0; k--) {
        if (st[k].rank === 'genus') { g = st[k].name; break; }
      }
      (byGenus[g] = byGenus[g] || []).push({idx: idx, count: p.count});
    }
    st.push({d: p.depth, name: p.name, rank: p.rank});
  });
  var keep = {};
  Object.keys(byGenus).forEach(function(g) {
    byGenus[g].slice().sort(function(a, b) { return b.count - a.count; })
      .slice(0, topN).forEach(function(x) { keep[x.idx] = true; });
  });
  var out = [], dropDepth = null;
  rows.forEach(function(p, idx) {
    if (dropDepth !== null) {
      if (p.depth > dropDepth) return;
      dropDepth = null;
    }
    if (p.rank === 'species') {
      if (keep[idx]) out.push(p); else dropDepth = p.depth;
    } else out.push(p);
  });
  return out;
}
function sankeyRows() {
  return REPORT.filter(function(r) {
    return r.rank !== 'no rank' || r.name === 'root' || r.name === 'unclassified';
  });
}
var CM = {'root': '#d4d4d4', 'Viruses': '#e41a1c', 'Riboviria': '#377eb8',
  'Orthornavirae': '#4daf4a', 'Kitrinoviricota': '#984ea3',
  'Alsuviricetes': '#ff7f00', 'Martellivirales': '#a65628',
  'Virgaviridae': '#f781bf'};
function renderSankey() {
  var rows = sankeyRows();
  if (document.getElementById('skTop5').checked) rows = limitSpeciesPerGenus(rows, 5);
  var labels = [], srcs = [], tgts = [], vals = [], n2i = {}, stack = [];
  for (var i = 0; i < rows.length; i++) {
    var n = rows[i];
    while (stack.length > 0 && stack[stack.length - 1].d >= n.depth) stack.pop();
    if (!(n.name in n2i)) { n2i[n.name] = labels.length; labels.push(n.name); }
    var ci = n2i[n.name];
    if (stack.length > 0) {
      var pi = stack[stack.length - 1].i;
      if (pi !== ci) {
        var dup = false;
        for (var k = 0; k < srcs.length; k++) {
          if (srcs[k] === pi && tgts[k] === ci) { vals[k] += n.count; dup = true; break; }
        }
        if (!dup) { srcs.push(pi); tgts.push(ci); vals.push(n.count); }
      }
    }
    stack.push({d: n.depth, i: ci});
  }
  if (labels.length > 1 && srcs.length > 0) {
    Plotly.newPlot('sk', [{type: 'sankey', orientation: 'h',
      node: {pad: 20, thickness: 25, line: {color: 'black', width: 0.5},
             label: labels,
             color: labels.map(function(l) { return CM[l] || '#2e86c1'; })},
      link: {source: srcs, target: tgts, value: vals,
             color: 'rgba(46,134,193,0.15)'}}],
      {margin: {t: 10, b: 10, l: 10, r: 10}, font: {size: 11}},
      {responsive: true});
  } else {
    document.getElementById('sk').innerHTML =
      '<div style="color:#999;padding:20px">Insufficient hierarchy levels for Sankey diagram</div>';
  }
}

var STD = {'superkingdom': 1, 'acellular root': 1, 'realm': 1, 'subrealm': 1,
  'kingdom': 1, 'subkingdom': 1, 'phylum': 1, 'subphylum': 1, 'class': 1,
  'subclass': 1, 'order': 1, 'suborder': 1, 'family': 1, 'subfamily': 1,
  'genus': 1, 'subgenus': 1, 'species': 1};
function renderClassTable() {
  var tb = document.getElementById('ctbody');
  var html = '';
  REPORT.forEach(function(r) {
    var rankColor = r.rank === 'no rank' ? '#aaa' : '#1a5276';
    var isStd = !!STD[r.rank];
    html += '<tr><td style="color:' + rankColor + ';white-space:nowrap">' +
      escapeHtml(r.rank) + '</td>' +
      '<td><span style="padding-left:' + (r.depth * 8) + 'px' +
      (isStd ? ';font-weight:600' : '') + '">' + escapeHtml(r.name) + '</span></td>' +
      '<td class="num" style="color:#888">' + escapeHtml(r.taxid) + '</td>' +
      '<td class="num">' + escapeHtml(r.proportion) + '</td><td class="num">' +
      escapeHtml(r.count.toLocaleString()) + '</td></tr>';
  });
  tb.innerHTML = html;
}

/* ---- Virus Sequence Classification（contig）或 Classified Sequences（reads） ---- */
var RK = ['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'];
var HAS_HOST = false;
function nearComplete(len, avg) {
  if (!len || !avg) return false;
  var l = Number(len), a = Number(avg);
  if (isNaN(l) || isNaN(a) || a <= 0) return false;
  var ratio = l / a;
  return ratio >= 0.8 && ratio <= 1.2;
}
function vcTbody(rows) {
  var RKX = RK.concat(HAS_HOST ? ['host'] : []);
  var tbody = '';
  rows.forEach(function(r) {
    var nc = nearComplete(r.length, r.genus_avg_len);
    var cs = r.contig && r.contig.length > 35 ? r.contig.substring(0, 32) + '...' : (r.contig || '');
    tbody += '<tr' + (nc ? ' style="background:#e8f5e9"' : '') + '><td style="font-family:monospace;font-size:11px;max-width:220px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="' + escapeHtml(r.contig) + '">' + escapeHtml(cs) + '</td>';
    RKX.forEach(function(k) {
      if (k === 'host') {
        tbody += '<td style="white-space:nowrap">' + (r.host ? escapeHtml(r.host) : '<span style="color:#ccc">-</span>') + '</td>';
        return;
      }
      var v = r[k] || '';
      var strong = (k === 'genus' || k === 'species') ? ';font-weight:600' : '';
      tbody += '<td style="white-space:nowrap' + strong + '">' + (v ? escapeHtml(v) : '<span style="color:#ccc">-</span>') + '</td>';
    });
    tbody += '<td class="num"' + (nc ? ' style="font-weight:700;color:#1b5e20"' : '') + '>' + escapeHtml(r.length) + '</td>';
    var sc = r.score;
    tbody += '<td class="num">' + (sc != null && sc !== '' ? Number(sc).toFixed(4) : '<span style="color:#ccc">-</span>') + '</td>';
    var ev = r.e_value;
    var evD = '-';
    if (ev != null && ev !== '' && ev !== '-1') {
      var en = Number(ev);
      if (en === 0) evD = '0';
      else if (en < 0.01) evD = en.toExponential(1);
      else evD = en.toFixed(3);
    }
    tbody += '<td class="num">' + evD + '</td>';
    var gal = r.genus_avg_len;
    tbody += '<td class="num">' + (gal != null && gal !== '' ? Number(gal).toFixed(0) : '<span style="color:#ccc">-</span>') + '</td>';
    tbody += '<td style="white-space:nowrap">'
      + '<button class="btn small" style="font-size:10px;margin:1px" onclick="runAnalysisFor(\'cdd\', \'' + escapeHtml(r.contig) + '\')">CDD</button>'
      + '<button class="btn small" style="font-size:10px;margin:1px" onclick="runAnalysisFor(\'blastn\', \'' + escapeHtml(r.contig) + '\')">BLASTN</button>'
      + '<button class="btn small" style="font-size:10px;margin:1px" onclick="runAnalysisFor(\'blastx\', \'' + escapeHtml(r.contig) + '\')">BLASTX</button>'
      + '<button class="btn small" style="font-size:10px;margin:1px" onclick="runAnalysisFor(\'primer\', \'' + escapeHtml(r.contig) + '\')">Primer</button>'
      + '</td></tr>';
  });
  return tbody;
}
function idsTbody(rows) {
  var tbody = '';
  rows.forEach(function(r) {
    tbody += '<tr><td style="font-family:monospace;font-size:11px">' + escapeHtml(r.seq_id || '') + '</td>' +
      '<td>' + escapeHtml(r.taxid || '') + '</td>' +
      '<td style="text-align:right;font-family:monospace;font-size:11px">' + escapeHtml(r.score || '') + '</td></tr>';
  });
  return tbody;
}
var _vSort = { col: -1, dir: -1 };
function sortVirusTable(colIdx) {
  var key = (MODE === 'contig')
    ? RK.concat(['length', 'score', 'e_value', 'genus_avg_len'].concat(HAS_HOST ? ['host'] : []))[colIdx]
    : ['seq_id', 'taxid', 'score'][colIdx];
  if (!key) return;
  if (_vSort.col === colIdx) _vSort.dir = -_vSort.dir;
  else { _vSort.col = colIdx; _vSort.dir = -1; }
  var isNum = ['length', 'score', 'genus_avg_len', 'taxid'].indexOf(key) >= 0;
  var src = (MODE === 'contig') ? VC : IDS;
  var sorted = src.slice().sort(function(a, b) {
    var va = a[key], vb = b[key];
    if (va == null || va === '') va = isNum ? -Infinity : '';
    if (vb == null || vb === '') vb = isNum ? -Infinity : '';
    if (isNum) { va = Number(va); vb = Number(vb); }
    else { va = String(va).toLowerCase(); vb = String(vb).toLowerCase(); }
    if (va < vb) return -1 * _vSort.dir;
    if (va > vb) return 1 * _vSort.dir;
    return 0;
  });
  document.getElementById('vtbody').innerHTML = (MODE === 'contig') ? vcTbody(sorted) : idsTbody(sorted);
}
function renderVirusTable() {
  var th = document.getElementById('vthead');
  var labels, keys;
  if (MODE === 'contig') {
    keys = RK.concat(['length', 'score', 'e_value', 'genus_avg_len'].concat(HAS_HOST ? ['host'] : []));
    labels = ['Contig', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species',
              'Length', 'Score', 'E-value', 'Genus Avg Len'].concat(HAS_HOST ? ['宿主(预测)'] : []);
  } else {
    keys = ['seq_id', 'taxid', 'score'];
    labels = ['Seq ID', 'TaxID', 'Score'];
  }
  th.innerHTML = '<th style="width:' + (MODE === 'contig' ? 240 : 60) + 'px">' +
    (MODE === 'contig' ? 'Contig' : 'Seq ID') + '</th>' +
    keys.slice(1).map(function(k, i) {
      var idx = i + 1;
      var numCls = ['length', 'score', 'e_value', 'genus_avg_len', 'taxid'].indexOf(k) >= 0 ? ' class="num"' : '';
      return '<th' + numCls + ' class="sortable" style="cursor:pointer;user-select:none;white-space:nowrap" onclick="sortVirusTable(' + idx + ')">' +
        escapeHtml(labels[idx]) + '</th>';
    }).join('') + '<th style="width:240px">Analyze</th>';
  var tb = document.getElementById('vtbody');
  tb.innerHTML = (MODE === 'contig') ? vcTbody(VC) : idsTbody(IDS);
}
function runAnalysisFor(tool, contig) {
  if (window.parent && typeof window.parent.runAnalyze === 'function') {
    try { window.parent.runAnalyze(tool, contig); return; } catch (e) {}
  }
  alert('请在「病毒识别和分类分析」页打开本报告以使用分析功能');
}

/* ---- 宿主统计图（若有宿主预测） ---- */
function renderHostBar() {
  if (!HOST.length) return;
  var h4s = document.querySelectorAll('h4');
  for (var i = 0; i < h4s.length; i++) {
    if (h4s[i].textContent.indexOf('宿主预测统计') === 0) {
      Plotly.newPlot('hostBar', [{type: 'bar',
        x: HOST.map(function(h) { return h.cat; }),
        y: HOST.map(function(h) { return h.n; }),
        marker: {color: '#1d7a4f'}}],
        {margin: {t: 10, b: 60, l: 50, r: 10}, yaxis: {title: 'contigs'}},
        {responsive: true});
      break;
    }
  }
}

/* ---- Krona 替代：旭日图 ---- */
function renderSun() {
  var showUnc = document.getElementById('kronaUnc').checked;
  var rows = REPORT.filter(function(r) {
    if (r.rank === 'unclassified') return showUnc;
    if (r.rank === 'no rank') return r.name === 'root';
    return true;
  });
  if (document.getElementById('skTop5').checked) rows = limitSpeciesPerGenus(rows, 5);
  var ids = [], lab = [], par = [], val = [], stack = [];
  for (var i = 0; i < rows.length; i++) {
    var n = rows[i];
    while (stack.length > 0 && stack[stack.length - 1].d >= n.depth) stack.pop();
    var sid = 'n' + i;
    ids.push(sid);
    lab.push(n.name.length > 26 ? n.name.slice(0, 24) + '..' : n.name);
    par.push(stack.length > 0 ? stack[stack.length - 1].id : '');
    val.push(n.count);
    stack.push({d: n.depth, id: sid});
  }
  if (!ids.length) { document.getElementById('krona').innerHTML = '<div style="color:#999;padding:20px">无数据</div>'; return; }
  Plotly.newPlot('krona', [{type: 'sunburst', ids: ids, labels: lab,
    parents: par, values: val, branchvalues: 'remainder', maxdepth: 6}],
    {margin: {t: 10, b: 10, l: 10, r: 10}}, {responsive: true});
}

renderSankey();
renderClassTable();
renderVirusTable();
renderSun();
renderHostBar();
</script></body></html>"""

    html = head + body_mid + tail
    out = check_path(os.path.join(run_dir, 'report.html'), must_exist=False,
                     in_platform=True)
    with safe_open(out, 'wt') as f:
        f.write(html)
    if logger:
        logger.log(f"contig 分类报告生成: {out}")
    return out


def esc_html(v):
    return (str(v).replace('&', '&amp;').replace('<', '&lt;')
            .replace('>', '&gt;').replace('"', '&quot;'))


def _taxburst_nodes(krows):
    """kreport 行 → taxburst 节点树（name/rank/count/children）。

    count 用 kreport 累计片段数（父 ≥ 子之和），符合 taxburst 预期。
    """
    top, stack = [], []
    for r in krows:
        frags = int(r['frags'])
        node = {'name': r['name'], 'rank': r['rank'], 'count': frags,
                'children': []}
        while stack and stack[-1][0] >= r['depth']:
            stack.pop()
        if stack:
            stack[-1][1]['children'].append(node)
        else:
            top.append(node)
        stack.append((r['depth'], node))
    return top


def build_taxburst(run_dir, krows, logger=None):
    """生成 taxburst 交互 HTML（Krona 等价物，离线 d3 内嵌）。

    产出 taxburst.html（仅已分类）与 taxburst_unc.html（含未分类）。
    taxburst 未安装时返回 (False, False)，调用方回退 plotly 旭日图。
    """
    try:
        from taxburst.output import generate_html as _tb_generate
    except ImportError:
        if logger:
            logger.log("taxburst 未安装，Krona 区块改用 plotly 旭日图", "WARN")
        return False, False
    top_all = _taxburst_nodes(krows)
    top_cls = [n for n in top_all if n['rank'] != 'unclassified']
    if not top_cls:
        return False, False

    def _write(html, name):
        with safe_open(check_path(os.path.join(run_dir, name),
                                  must_exist=False, in_platform=True),
                       'wt') as f:
            f.write(html)

    _write(_tb_generate(top_cls, name='Viruses'), 'taxburst.html')
    unc = any(n['rank'] == 'unclassified' for n in top_all)
    if unc:
        _write(_tb_generate(top_all, name='all'), 'taxburst_unc.html')
    if logger:
        logger.log("taxburst 旭日图生成: taxburst.html"
                   + (" (+taxburst_unc.html)" if unc else ""))
    return True, unc


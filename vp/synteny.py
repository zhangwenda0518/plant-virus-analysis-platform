# -*- coding: utf-8 -*-
"""
同属病毒共线性比较（lovis4u / clinker 的平台内置等价实现，Windows 可用）。

输入：GenBank 集合（vp/gb_collection 下载/导入）或任意本机 .gb 文件。
流程：
  1. 基因模型提取：CDS /translation 优先，缺失则按位置翻译核苷酸；
     单一 polyprotein + ≥3 mat_peptide 的多聚蛋白基因组（如马铃薯 Y 病毒属）
     自动改用成熟肽段做基因，避免"一条长 ORF"把所有同属基因组连成一团；
  2. MMseqs2 easy-search 全对全蛋白比较（内置 mmseqs/bin/mmseqs.exe）；
  3. 单 linkage 聚类成基因家族（identity / 双侧覆盖度 / E-value 三阈值，
     口径同 clinker 默认：0.3 / 0.5 / 1e-3）；
  4. 出版级共线性图（matplotlib，风格对齐 clinker/PHROG 类文献图）：
     基因画粗块箭头（按产物关键词归功能类别着色，未知基因灰色）、
     相邻基因组同源基因间浅灰连线带、轨道端点坐标、比例尺、底部类别图例；
     输出 SVG（矢量，可投稿排版）+ PNG + 网页版（SVG 内嵌，
     鼠标悬停基因显示产物/位置/家族信息）。

输出（默认 <集合目录>/compare/）：
  compare.html          网页版共线性图（离线可开，SVG 悬浮提示）
  compare.svg           矢量图（出版排版用）
  compare.png           位图（200dpi）
  gene_clusters.tsv     基因家族清单（家族号/基因数/代表产物/成员）
  genome_similarity.tsv 成对共享基因比例与平均一致性
  cds_proteins.faa      提取的全部 CDS 蛋白（可复用于其他分析）
  cds_allvsall.m8       MMseqs2 全对全原始命中（透明可复查）
  summary.json          参数、基因组清单、警告
"""
import os
import re
import sys
import json
import time
import statistics
from collections import Counter, defaultdict

from .config import DIRS, get_config
from .utils import check_path, safe_open, run_cmd
from .gb_collection import resolve_gb_files, gb_collection_dir

_MIN_PROT_LEN = 10           # 短于该长度(aa)的基因不参与比较
_PAIR_TABLE_CAP = 40         # HTML 成对表最多行数
_CLUSTER_TABLE_CAP = 200


# ------------------------------------------------------------------
# 1. 基因模型提取
# ------------------------------------------------------------------
def _short_label(product, gene):
    return (product or gene or '').strip() or 'hypothetical protein'


def _extract_genome(rec, mat_peptide=True):
    """单条 SeqIO 记录 → genome dict（genes 含蛋白序列）；无可比基因返回 None。"""
    cds = [f for f in rec.features
           if f.type == 'CDS' and not f.qualifiers.get('pseudo')
           and f.location is not None]
    mat = [f for f in rec.features
           if f.type == 'mat_peptide' and f.location is not None]
    feats, model = cds, 'CDS'
    if mat_peptide and len(cds) == 1 and len(mat) >= 3:
        feats, model = mat, 'mat_peptide'
    genes = []
    for i, f in enumerate(feats):
        q = f.qualifiers
        product = _short_label((q.get('product') or [''])[0],
                               (q.get('gene') or [''])[0])
        prot = (q.get('translation') or [''])[0]
        if not prot:
            nt = f.extract(rec.seq)
            if len(nt) < 3:
                continue
            prot = str(nt.translate(table=1))
        prot = prot.rstrip('*').replace('*', 'X').upper()
        if len(prot) < _MIN_PROT_LEN:
            continue
        genes.append({'gid': f'{rec.id}__{i}', 'acc': rec.id,
                      'start': int(f.location.start),
                      'end': int(f.location.end),
                      'strand': '-' if f.location.strand == -1 else '+',
                      'product': product, 'prot': prot})
    if not genes:
        return None
    return {'acc': rec.id,
            'organism': rec.annotations.get('organism', '') or '',
            'length': len(rec.seq), 'model': model, 'genes': genes}


def extract_genomes(gb_files, mat_peptide=True, logger=None):
    """逐个 GenBank 文件提取基因组模型（每文件可含多条记录）。

    返回 (genomes, warnings)；重复 accession 只保留首个。
    """
    from Bio import SeqIO

    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    genomes, seen, warnings = [], set(), []
    for path in gb_files:
        try:
            recs = list(SeqIO.parse(path, 'genbank'))
        except (OSError, ValueError) as e:
            warnings.append(f'{os.path.basename(path)}: 解析失败（{e}）')
            log(f'{os.path.basename(path)} 解析失败，已跳过: {e}', 'WARN')
            continue
        if not recs:
            warnings.append(f'{os.path.basename(path)}: 无 GenBank 记录')
            continue
        for rec in recs:
            if rec.id in seen:
                log(f'{rec.id} 重复，跳过 {os.path.basename(path)}', 'WARN')
                continue
            seen.add(rec.id)
            g = _extract_genome(rec, mat_peptide=mat_peptide)
            if g is None:
                w = f'{rec.id}: 无可比基因（缺 CDS 注释？），已跳过'
                warnings.append(w)
                log(w, 'WARN')
                continue
            g['file'] = path
            genomes.append(g)
            log(f"{rec.id}: {g['length']} bp, {len(g['genes'])} 基因 "
                f"[{g['model']}]")
    return genomes, warnings


# ------------------------------------------------------------------
# 2. MMseqs2 全对全
# ------------------------------------------------------------------
def _run_mmseqs_allvsall(faa, m8, threads, logger):
    exe = get_config().tools.get('mmseqs')
    if not exe:
        raise RuntimeError('未检测到 MMseqs2（期望 mmseqs/bin/mmseqs.exe）')
    env = os.environ.copy()
    env['PATH'] = os.path.dirname(exe) + os.pathsep + env.get('PATH', '')
    run_cmd([exe, 'easy-search', faa, faa, m8, m8 + '.tmp',
             '--format-output', 'query,target,fident,alnlen,evalue',
             '-e', '1e-3', '--max-seqs', '10000', '-s', '7.5',
             '--threads', str(threads)], logger=logger, env=env)


def _parse_hits(m8, gene_len):
    """mmseqs 输出 → 去重成对命中列表（q<t 归一，去自命中）。"""
    hits, seen = [], set()
    with safe_open(m8) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 5:
                continue
            q, t = p[0], p[1]
            if q == t:
                continue
            key = (q, t) if q < t else (t, q)
            if key[0] not in gene_len or key[1] not in gene_len:
                continue
            if key in seen:
                continue
            seen.add(key)
            try:
                # mmseqs fident 为 0-1 小数（BLAST 的 pident 才是百分数）
                hits.append({'q': key[0], 't': key[1],
                             'fident': float(p[2]),
                             'alnlen': int(float(p[3])),
                             'evalue': float(p[4])})
            except ValueError:
                continue
    for h in hits:
        h['cov_q'] = h['alnlen'] / gene_len[h['q']]
        h['cov_t'] = h['alnlen'] / gene_len[h['t']]
    return hits


# ------------------------------------------------------------------
# 3. 聚类（单 linkage，阈值口径同 clinker 默认）
# ------------------------------------------------------------------
class _UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def cluster_genes(genes, hits, min_ident, min_cov, max_evalue):
    """阈值过滤 + 单 linkage → (groups{root:[gid]}, kept_hits)。"""
    uf = _UF()
    kept = [h for h in hits
            if h['fident'] >= min_ident and h['cov_q'] >= min_cov
            and h['cov_t'] >= min_cov and h['evalue'] <= max_evalue]
    for h in kept:
        uf.union(h['q'], h['t'])
    groups = defaultdict(list)
    for g in genes:
        groups[uf.find(g['gid'])].append(g['gid'])
    return dict(groups), kept


def _represent(members):
    """家族代表产物：出现最多的 product（hypothetical 兜底）。"""
    c = Counter(m['product'] for m in members)
    if 'hypothetical protein' in c and len(c) > 1:
        c.pop('hypothetical protein')
    return c.most_common(1)[0][0]


def _assign(genomes, groups):
    """聚类结果落到基因上：cid（singleton 归为一类）+ 家族元数据。"""
    gid2root = {}
    for idx, (root, members) in enumerate(
            sorted(groups.items(), key=lambda kv: -len(kv[1]))):
        if len(members) < 2:
            continue                      # 单基因家族 → singleton
        cid = f'GC{idx + 1:05d}'
        for gid in members:
            gid2root[gid] = cid
    for g in genomes:
        for gene in g['genes']:
            gene['cid'] = gid2root.get(gene['gid'], 'singleton')
    row_of = {gene['gid']: gi for gi, g in enumerate(genomes)
              for gene in g['genes']}
    by_cid = defaultdict(list)
    for g in genomes:
        for gene in g['genes']:
            by_cid[gene['cid']].append(gene)
    meta = {}
    for cid, members in by_cid.items():
        rep = _represent(members)
        meta[cid] = {'rep': rep, 'n_genes': len(members),
                     'n_genomes': len({row_of[m['gid']] for m in members})}
    return row_of, by_cid, meta


# ------------------------------------------------------------------
# 4. 共线性图（matplotlib，风格对齐 clinker/PHROG 类文献图）
# ------------------------------------------------------------------
# 产物关键词 → 功能类别（植物病毒口径；顺序即优先级，未命中已知关键词
# 但有产物名的归"其他已知"，hypothetical/无名 → 灰色"未知功能"）
_CAT_RULES = [
    ('复制与聚合', ('replicase', 'replication', 'polymerase', 'rdrp',
                    'rna-dependent rna', 'helicase', 'methyltransfer',
                    'primase', 'rna replic')),
    ('外壳蛋白', ('coat protein', 'capsid', 'coat', 'envelope protein',
                  'nucleocapsid', 'structural protein')),
    ('运动蛋白', ('movement',)),
    ('蛋白酶', ('protease', 'peptidase', 'hc-pro')),
    ('VPg', ('vpg',)),
    ('传播与辅助', ('helper component', 'transmission', 'aphid',
                   'vector protein')),
    ('沉默抑制子', ('suppressor', 'silencing', 'silencer')),
]
_CAT_OTHER = '其他已知'
_CAT_UNK = '未知功能'
_CAT_COLORS = {'复制与聚合': '#4c78a8', '外壳蛋白': '#e45756',
               '运动蛋白': '#54a24b', '蛋白酶': '#f58518', 'VPg': '#b279a2',
               '传播与辅助': '#72b7b2', '沉默抑制子': '#9d755d',
               _CAT_OTHER: '#edc949', _CAT_UNK: '#c8c8c8'}


def _categorize(product):
    p = (product or '').lower()
    for name, kws in _CAT_RULES:
        for kw in kws:
            if kw in p:
                return name
    if p and p != 'hypothetical protein':
        return _CAT_OTHER
    return _CAT_UNK


def _arrow_verts(gene, y, hh, max_len):
    """基因块箭头多边形顶点（数据坐标），指向按链方向。"""
    s, e = gene['start'], gene['end']
    ah = min((e - s) * 0.45, max_len * 0.018)
    if gene['strand'] != '-':
        xs = [s, e - ah, e, e - ah, s]
    else:
        xs = [e, s + ah, s, s + ah, e]
    return [(x, y + d) for x, d in zip(xs, (-hh, -hh, 0, hh, hh))]


def _cluster_palette(n):
    """n 个家族颜色：优先 distinctipy（lovis4u 同款），回退定性色板。"""
    try:
        import distinctipy
        from matplotlib.colors import to_hex
        excl = [(1, 1, 1), (0.85, 0.85, 0.85), (0.95, 0.95, 0.95)]
        cols = distinctipy.get_colors(max(1, n), exclude=excl,
                                      n_attempts=5000)
        return [to_hex(c) for c in cols]
    except Exception:
        import plotly.colors as _pc
        pal = _pc.qualitative.Alphabet + _pc.qualitative.Safe
        return [pal[i % len(pal)] for i in range(max(1, n))]


def _draw_figure(genomes, row_of, by_cid, meta, link_ident, max_len,
                 show_labels=True, style='lovis'):
    """出版级共线性主图 → (matplotlib Figure, {gid: 悬浮文本})。

    style='lovis'：对齐 LoVis4u 画廊风格——同源家族逐个着色
    （distinctipy 调色盘，singleton 灰）、浅灰半透明同源连线带、
    左侧三行标签（坐标区间 / 斜体物种名 / 灰 accession）、无图例。
    style='category'：按产物功能类别着色 + 轨道端点坐标 + 比例尺 +
    底部类别图例（clinker/PHROG 文献图风格）。
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms
    from matplotlib.patches import Polygon, Patch
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei',
                                       'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['svg.fonttype'] = 'none'

    n = len(genomes)
    gh, rh = 0.56, 0.26                     # 基因高 / 连线带半高
    row = 1.0                               # 轨道间距
    fig_w = 13.0
    fig_h = 1.35 + n * 0.92
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.185, right=0.985, top=0.995, bottom=0.105)
    ax.set_xlim(-0.012 * max_len, max_len * 1.035)
    ax.set_ylim(-(n - 1) - 1.05, 0.62)
    ax.axis('off')
    tr_lbl = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)

    # 家族颜色：lovis=逐家族 distinctipy；category=功能类别
    fam_cids = [c for c in by_cid if c != 'singleton']
    fam_cids.sort(key=lambda c: -len(by_cid[c]))
    if style == 'lovis':
        fam_colors = {cid: col for cid, col in
                      zip(fam_cids, _cluster_palette(len(fam_cids)))}
        gene_color = lambda gene: fam_colors.get(gene['cid'], '#c8c8c8')
    else:
        gene_color = lambda gene: _CAT_COLORS[_categorize(gene['product'])]
    link_alpha = 0.30 if style == 'lovis' else 1.0

    # 相邻基因组同源连线带（浅灰，底层）
    for gi in range(n - 1):
        ya, yb = -gi * row, -(gi + 1) * row
        ga, gb = genomes[gi], genomes[gi + 1]
        sides = defaultdict(lambda: ([], []))
        for gene in ga['genes']:
            if gene['cid'] != 'singleton':
                sides[gene['cid']][0].append(gene)
        for gene in gb['genes']:
            if gene['cid'] in sides:
                sides[gene['cid']][1].append(gene)
        for cid, (ma, mb) in sides.items():
            ma = sorted(ma, key=lambda z: z['start'])
            mb = sorted(mb, key=lambda z: z['start'])
            if not ma or not mb:
                continue
            if len(ma) * len(mb) <= 12:
                pairs = [(a, b) for a in ma for b in mb]
            else:                            # 多拷贝家族按顺序就近连线
                pairs = [(a, mb[int(round(j * (len(mb) - 1) /
                                          max(1, len(ma) - 1)))])
                         for j, a in enumerate(ma)]
            for a, b in pairs:
                ax.add_patch(Polygon(
                    [(a['start'], ya), (a['end'], ya),
                     (b['end'], yb), (b['start'], yb)],
                    closed=True, facecolor='#c9c9c9', edgecolor='none',
                    alpha=link_alpha, zorder=1))

    # 骨架线（category 风格另加两端坐标）
    for gi, g in enumerate(genomes):
        y = -gi * row
        ax.plot([0, g['length']], [y, y],
                color='#666666' if style == 'lovis' else '#b8b8b8',
                lw=0.8, zorder=2)
        if style != 'lovis':
            ax.text(0, y - 0.44, '1', ha='left', va='center',
                    fontsize=7, color='#8a8a8a', zorder=2)
            ax.text(g['length'], y - 0.44, f"{g['length']:,}", ha='right',
                    va='center', fontsize=7, color='#8a8a8a', zorder=2)

    # 基因块箭头（lovis=逐家族着色；category=功能类别着色）
    title_dict = {}
    for gi, g in enumerate(genomes):
        y = -gi * row
        for gene in g['genes']:
            cat = _categorize(gene['product'])
            patch = Polygon(_arrow_verts(gene, y, gh / 2, max_len),
                            closed=True, facecolor=gene_color(gene),
                            edgecolor='#303030', linewidth=0.7,
                            joinstyle='miter', zorder=3)
            patch.set_gid(f"g_{gene['gid']}")
            ax.add_patch(patch)
            cid = gene['cid']
            fam = ('无同源' if cid == 'singleton'
                   else f"家族 {cid} {meta[cid]['rep'][:20]}"
                        f"（{meta[cid]['n_genes']} 基因/{meta[cid]['n_genomes']}"
                        "基因组）")
            title_dict[gene['gid']] = (
                f"{gene['product']} · {cat}\n{gene['acc']}:"
                f"{gene['start'] + 1}-{gene['end']} ({gene['strand']}) · "
                f"{len(gene['prot'])} aa\n{fam}")

    # 左侧标签：lovis=三行（坐标区间/斜体物种名/灰accession）；
    # category=两行（斜体物种名/灰accession）
    for gi, g in enumerate(genomes):
        y = -gi * row
        org = g['organism'] or g['acc']
        kw = {'fontfamily': 'DejaVu Sans', 'fontstyle': 'italic'} \
            if org.isascii() else {}
        if style == 'lovis':
            ax.text(-0.012, y + 0.28, f"1 - {g['length']:,}",
                    transform=tr_lbl, ha='right', va='center', fontsize=7,
                    color='#8a8a8a', clip_on=False)
            ax.text(-0.012, y + 0.04, org, transform=tr_lbl, ha='right',
                    va='center', fontsize=8.5, color='#1f1f1f',
                    clip_on=False, **kw)
            ax.text(-0.012, y - 0.20, g['acc'], transform=tr_lbl,
                    ha='right', va='center', fontsize=7.5, color='#6b6b6b',
                    clip_on=False)
        else:
            ax.text(-0.012, y + 0.17, org, transform=tr_lbl, ha='right',
                    va='center', fontsize=8.5, color='#1f1f1f',
                    clip_on=False, **kw)
            ax.text(-0.012, y - 0.11, g['acc'], transform=tr_lbl,
                    ha='right', va='center', fontsize=8, color='#6b6b6b',
                    clip_on=False)

    # 基因名标注（防碰撞：估宽后贪心放到 2 个高度层，放不下不标；
    # 左右两端避开轨道端点坐标区）
    total = sum(len(g['genes']) for g in genomes)
    if show_labels and total <= 120:
        usable_px, char_px = 1000.0, 4.8
        px2bp = max_len / usable_px
        levels, step = 2, 0.21
        for g in genomes:
            level_end = [-1e18] * levels
            for gene in sorted(g['genes'], key=lambda z: z['start']):
                if gene['cid'] == 'singleton':
                    continue
                text = gene['product']
                if len(text) > 24:
                    text = text[:23] + '…'
                half_bp = len(text) * char_px * px2bp / 2
                x0 = (gene['start'] + gene['end']) / 2
                max_right = g['length'] - 0.034 * max_len
                min_left = 0.014 * max_len
                if x0 + half_bp > max_right:
                    x0 = max_right - half_bp
                if x0 + half_bp > max_len:
                    continue          # 轨道太短放不下完整标签
                if x0 - half_bp < min_left:
                    x0 = min_left + half_bp
                placed = next((li for li in range(levels)
                               if x0 - half_bp > level_end[li]), None)
                if placed is None:
                    continue
                level_end[placed] = x0 + half_bp
                ax.text(x0, -row_of[gene['gid']] + gh / 2 + 0.08
                        + placed * step, text, ha='center', va='bottom',
                        fontsize=7.5, color='#1f1f1f', zorder=4)

    if style != 'lovis':
        # 比例尺（左下，长度取 ~1/6 基因组长的整值）
        nice = next((c for c in (1000, 2000, 5000, 10000, 20000, 50000,
                                 100000, 200000, 500000, 1000000)
                     if c >= max_len / 6), 1000000)
        y_sb = -(n - 1) - 0.82
        ax.plot([0, nice], [y_sb, y_sb], color='#1f1f1f', lw=1.4,
                solid_capstyle='butt', zorder=2)
        for x in (0, nice):
            ax.plot([x, x], [y_sb - 0.07, y_sb + 0.07], color='#1f1f1f',
                    lw=1.4, zorder=2)
        ax.text(nice / 2, y_sb - 0.24, f"{nice / 1000:g} kb", ha='center',
                va='top', fontsize=8, color='#1f1f1f')

        # 底部功能类别图例（仅画出现过的类别）
        present = {_categorize(gene['product'])
                   for g in genomes for gene in g['genes']}
        order = [name for name, _ in _CAT_RULES] + [_CAT_OTHER, _CAT_UNK]
        handles = [Patch(facecolor=_CAT_COLORS[c], edgecolor='#303030',
                         linewidth=0.7, label=c)
                   for c in order if c in present]
        fig.legend(handles=handles, loc='lower center',
                   bbox_to_anchor=(0.55, 0.004), ncol=min(len(handles), 9),
                   frameon=False, fontsize=8.5, handlelength=1.4,
                   handleheight=0.9, columnspacing=1.4)
    return fig, title_dict


def _inject_svg_titles(svg_path, title_dict):
    """SVG 基因路径注入 <title>，浏览器内嵌后悬浮即显示基因信息。"""
    import html as _html
    with safe_open(svg_path) as f:
        svg = f.read()
    n_inj = 0
    for gid, text in title_dict.items():
        pat = f'id="g_{gid}"'
        idx = svg.find(pat)
        if idx < 0:
            continue
        end = svg.find('/>', idx)
        if end < 0:
            continue
        t = _html.escape(text.replace('\n', ' · '), quote=False)
        svg = svg[:end] + f'><title>{t}</title></path' + svg[end + 2:]
        n_inj += 1
    with safe_open(svg_path, 'wt') as f:
        f.write(svg)
    return n_inj


_PAGE = """<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8">
<title>{title} · 同属共线性比较</title>
<style>
 body{{font-family:"Microsoft YaHei",system-ui,sans-serif;margin:22px;color:#1f2937}}
 h1{{font-size:20px;margin:0 0 6px}} h2{{font-size:16px;margin-top:22px}}
 .hint{{color:#6b7280;font-size:12.5px;margin:4px 0 10px}}
 table{{border-collapse:collapse;font-size:12.5px;margin-top:8px}}
 th,td{{border:1px solid #e5e7eb;padding:3px 8px;text-align:left}}
 th{{background:#f3f6fa}}
 .figbox{{overflow:auto;border:1px solid #e5e7eb;border-radius:6px;
          padding:10px;background:#fff}}
 .figbox svg{{max-width:100%;height:auto}}
</style></head>
<body>
<h1>同属病毒共线性比较 — {title}</h1>
<p class="hint">{n_genomes} 个基因组 · {n_genes} 个基因 · {n_clusters} 个同源家族 ·
阈值 identity≥{min_ident:.0%} / 双侧覆盖度≥{min_cov:.0%} / E≤{max_evalue:g} · 生成于 {date}<br>
基因按产物功能类别着色（灰=未知/无同源），浅灰连线带=相邻基因组同源家族；
<b>鼠标悬停基因</b>可查看产物/坐标/家族；矢量图 compare.svg 可直接用于论文排版。</p>
<div class="figbox">
{fig_html}
</div>
<h2>成对共享基因比例</h2>
<p class="hint">A→B 列为 A 中落在两基因组共享家族里的基因占比（≥80% 提示可能为同种/亚种分化，
40–80% 同属内近缘，明显更低则属内分化较远——经验参考，非 ICTV 划界标准）。</p>
{pair_html}
<h2>基因家族明细</h2>
{cluster_html}
</body></html>"""


# ------------------------------------------------------------------
# 主流程
# ------------------------------------------------------------------
def _run_lovis4u(gb_files, out_dir, logger):
    """可选外接引擎：调用 LoVis4u（pip install lovis4u）出原生 PDF 图。

    lovis4u 0.2.0 在 Windows 上有两处兼容 bug（临时文件句柄导致
    PermissionError、-smp 路径塞进 re.sub 模板导致 bad escape），
    由 vp/lovis4u_run.py 子进程修补；mmseqs 用平台内置 exe。
    返回 PDF 路径；不可用/失败返回 None（不影响内置图）。
    """
    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    import importlib.util
    if importlib.util.find_spec('lovis4u') is None:
        log('未安装 lovis4u（pip install lovis4u），跳过 LoVis4u PDF 输出',
            'WARN')
        return None
    exe = get_config().tools.get('mmseqs')
    if not exe:
        log('未检测到 MMseqs2，LoVis4u 无法聚类，跳过', 'WARN')
        return None
    list_file = os.path.join(out_dir, 'lovis4u_input.txt')
    with safe_open(list_file, 'wt') as f:
        f.write('\n'.join(gb_files) + '\n')
    runner = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'lovis4u_run.py')
    try:
        from vp.config import engine_cmd
        run_cmd(engine_cmd(runner, list_file, out_dir, 'lovis4u.pdf',
                 exe), logger=logger)
    except Exception as e:
        log(f'LoVis4u 运行失败（内置图不受影响）: {e}', 'WARN')
        return None
    pdf_path = os.path.join(out_dir, 'lovis4u.pdf')
    if os.path.isfile(pdf_path):
        log(f'LoVis4u 原生 PDF 已生成 → {pdf_path}')
        return pdf_path
    log('LoVis4u 未生成 PDF', 'WARN')
    return None


def run_comparison(name=None, files=None, out_dir=None, *,
                   min_ident=0.30, min_cov=0.50, max_evalue=1e-3,
                   mat_peptide=True, order=None, labels=True,
                   style='lovis', lovis4u_pdf=False,
                   threads=None, logger=None, prog=None, cancel=None):
    """同属比较主入口。name=GenBank 集合名 / files=.gb 路径列表（可并用）。

    返回 summary dict（含各产物路径）。
    """
    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    def cancelled():
        return cancel is not None and getattr(cancel, 'is_set',
                                              lambda: False)()

    gb_files = resolve_gb_files(name=name, files=files)
    if out_dir is None:
        if name:
            out_dir = os.path.join(gb_collection_dir(name), 'compare')
        else:
            out_dir = os.path.join(DIRS['results'], '_synteny',
                                   time.strftime('%Y%m%d_%H%M%S'))
    out_dir = check_path(out_dir, must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    threads = threads or get_config().threads

    if prog:
        prog('parse', 0.03, '解析 GenBank 基因模型')
    genomes, warnings = extract_genomes(gb_files, mat_peptide=mat_peptide,
                                        logger=logger)
    if len(genomes) < 2:
        raise RuntimeError('可用基因组不足 2 条（其余缺 CDS 注释被跳过），'
                           '无法比较；请检查集合内记录是否有 CDS 特征')
    if order:                              # 用户指定基因组顺序
        picked, rest = [], list(genomes)
        for acc in order:
            for g in list(rest):
                if g['acc'].upper().startswith(str(acc).upper()):
                    picked.append(g)
                    rest.remove(g)
        genomes = picked + rest
    log(f'共 {len(genomes)} 个基因组参与比较')

    if cancelled():
        raise RuntimeError('用户取消')

    faa = os.path.join(out_dir, 'cds_proteins.faa')
    with safe_open(faa, 'wt') as f:
        for g in genomes:
            for gene in g['genes']:
                f.write(f">{gene['gid']} {gene['acc']}|{gene['product']}\n"
                        f"{gene['prot']}\n")

    if prog:
        prog('search', 0.15, 'MMseqs2 全对全蛋白比较')
    m8 = os.path.join(out_dir, 'cds_allvsall.m8')
    _run_mmseqs_allvsall(faa, m8, threads, logger)
    if cancelled():
        raise RuntimeError('用户取消')

    gene_len = {gene['gid']: len(gene['prot'])
                for g in genomes for gene in g['genes']}
    hits = _parse_hits(m8, gene_len)
    log(f'全对全命中 {len(hits)} 对')
    groups, kept = cluster_genes([gene for g in genomes
                                  for gene in g['genes']],
                                 hits, min_ident, min_cov, max_evalue)
    if not any(len(m) >= 2 for m in groups.values()):
        log('阈值下无共享基因家族——这些基因组之间可能同源性过低，'
            '或阈值过严（可试 --min-ident 0.25 / --min-cov 0.4）', 'WARN')

    row_of, by_cid, meta = _assign(genomes, groups)
    link_ident = {(h['q'], h['t']): h['fident'] for h in kept}
    n_clusters = sum(1 for c in by_cid if c != 'singleton')
    n_single = len(by_cid.get('singleton', []))
    log(f'聚类: {n_clusters} 个同源家族（≥2基因），{n_single} 个无同源基因')

    if prog:
        prog('plot', 0.85, '绘制共线性图')
    fig, title_dict = _draw_figure(genomes, row_of, by_cid, meta,
                                   link_ident,
                                   max(g['length'] for g in genomes),
                                   show_labels=labels, style=style)
    svg_path = os.path.join(out_dir, 'compare.svg')
    png_path = os.path.join(out_dir, 'compare.png')
    fig.savefig(svg_path, bbox_inches='tight', pad_inches=0.15)
    fig.savefig(png_path, dpi=200, bbox_inches='tight', pad_inches=0.15)
    import matplotlib.pyplot as plt
    plt.close(fig)
    n_inj = _inject_svg_titles(svg_path, title_dict)
    log(f'SVG 悬浮提示注入 {n_inj}/{len(title_dict)} 个基因')
    with safe_open(svg_path) as f:
        svg_content = f.read()

    # 成对共享基因比例
    pairs = []
    for i in range(len(genomes)):
        for j in range(i + 1, len(genomes)):
            ci = {x['cid'] for x in genomes[i]['genes']}
            cj = {x['cid'] for x in genomes[j]['genes']}
            shared = {c for c in ci & cj if c != 'singleton'}
            ni, nj = len(genomes[i]['genes']), len(genomes[j]['genes'])
            si = sum(1 for x in genomes[i]['genes'] if x['cid'] in shared)
            sj = sum(1 for x in genomes[j]['genes'] if x['cid'] in shared)
            idents = [h['fident'] for h in kept
                      if {row_of[h['q']], row_of[h['t']]} == {i, j}]
            pairs.append({'a': genomes[i]['acc'], 'b': genomes[j]['acc'],
                          'shared': len(shared), 'frac_a': si / ni if ni else 0,
                          'frac_b': sj / nj if nj else 0,
                          'ident': (statistics.mean(idents) * 100)
                          if idents else None})
    with safe_open(os.path.join(out_dir, 'genome_similarity.tsv'), 'wt') as f:
        f.write('genome_a\tgenome_b\tshared_families\tshared_frac_a\t'
                'shared_frac_b\tmean_identity_pct\n')
        for p in pairs:
            ident = f"{p['ident']:.1f}" if p['ident'] is not None else ''
            f.write(f"{p['a']}\t{p['b']}\t{p['shared']}\t"
                    f"{p['frac_a']:.3f}\t{p['frac_b']:.3f}\t{ident}\n")

    with safe_open(os.path.join(out_dir, 'gene_clusters.tsv'), 'wt') as f:
        f.write('cluster\trepresentative\tn_genes\tn_genomes\tmembers\n')
        for cid in sorted((c for c in by_cid if c != 'singleton'),
                          key=lambda c: -len(by_cid[c])):
            members = by_cid[cid]
            mem = '; '.join(f"{m['acc']}|{m['product']}" for m in members)
            f.write(f"{cid}\t{meta[cid]['rep']}\t{len(members)}\t"
                    f"{meta[cid]['n_genomes']}\t{mem}\n")
        for m in by_cid.get('singleton', []):
            f.write(f"singleton\t{m['product']}\t1\t1\t"
                    f"{m['acc']}|{m['product']}\n")

    # HTML 明细表
    def _ident_cell(p):
        return f"{p['ident']:.1f}" if p['ident'] is not None else '—'

    pair_rows = ''.join(
        f"<tr><td>{p['a']}</td><td>{p['b']}</td><td>{p['shared']}</td>"
        f"<td>{p['frac_a']:.0%}</td><td>{p['frac_b']:.0%}</td>"
        f"<td>{_ident_cell(p)}</td></tr>"
        for p in pairs[:_PAIR_TABLE_CAP])
    pair_html = ('<table><tr><th>基因组 A</th><th>基因组 B</th>'
                 '<th>共享家族</th><th>A→B</th><th>B→A</th><th>平均一致性%</th>'
                 '</tr>' + pair_rows + '</table>')
    cids = sorted((c for c in by_cid if c != 'singleton'),
                  key=lambda c: -len(by_cid[c]))[:_CLUSTER_TABLE_CAP]
    cluster_html = ('<table><tr><th>家族</th><th>代表产物</th><th>基因数</th>'
                    '<th>基因组数</th><th>成员</th></tr>' + ''.join(
        f"<tr><td>{c}</td><td>{meta[c]['rep']}</td><td>{len(by_cid[c])}</td>"
        f"<td>{meta[c]['n_genomes']}</td><td>{'; '.join(m['acc'] for m in by_cid[c])}</td></tr>"
        for c in cids) + '</table>')

    html = _PAGE.format(
        title=name or '自备文件', n_genomes=len(genomes),
        n_genes=sum(len(g['genes']) for g in genomes), n_clusters=n_clusters,
        min_ident=min_ident, min_cov=min_cov, max_evalue=max_evalue,
        date=time.strftime('%Y-%m-%d %H:%M'), fig_html=svg_content,
        pair_html=pair_html, cluster_html=cluster_html)
    html_path = os.path.join(out_dir, 'compare.html')
    with safe_open(html_path, 'wt') as f:
        f.write(html)

    lovis4u_path = None
    if lovis4u_pdf:
        if prog:
            prog('lovis4u', 0.92, 'LoVis4u 原生 PDF 出图')
        lovis4u_path = _run_lovis4u(gb_files, out_dir, logger)

    summary = {
        'name': name, 'dir': out_dir, 'files': {
            'html': html_path, 'svg': svg_path, 'png': png_path,
            'lovis4u_pdf': lovis4u_path,
            'clusters': os.path.join(out_dir, 'gene_clusters.tsv'),
            'similarity': os.path.join(out_dir, 'genome_similarity.tsv'),
            'faa': faa, 'm8': m8},
        'n_genomes': len(genomes),
        'n_genes': sum(len(g['genes']) for g in genomes),
        'n_clusters': n_clusters, 'n_singletons': n_single,
        'genomes': [{'acc': g['acc'], 'organism': g['organism'],
                     'length': g['length'], 'n_genes': len(g['genes']),
                     'model': g['model']} for g in genomes],
        'pairs': pairs, 'warnings': warnings,
        'params': {'min_ident': min_ident, 'min_cov': min_cov,
                   'max_evalue': max_evalue, 'mat_peptide': mat_peptide,
                   'style': style, 'lovis4u': lovis4u_pdf},
    }
    with safe_open(os.path.join(out_dir, 'summary.json'), 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    if prog:
        prog('done', 1.0,
             f'完成: {len(genomes)} 基因组 · {n_clusters} 家族')
    log(f"比较完成 → {html_path}")
    return summary

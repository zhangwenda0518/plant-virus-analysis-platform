# -*- coding: utf-8 -*-
"""
⑥b 基因组示意图（littlegenomes 风格适配版）。

吸收 biolumber/littlegenomes 的画法：基线 + 500nt 刻度比例尺 + 特征矩形
（垂直偏移分层防重叠）；在此适配：
- ORF 画成**链方向箭头**（littlegenomes 为无方向矩形），按功能类别着色
- HMM/CDD 结构域画成**窄条**（s 宽度思路）叠放于 ORF 外侧，鼠标悬停可见
  域名（SVG <title>）
- 纯字符串拼 SVG，零第三方依赖（打包友好）
"""
import os
import re
import html

# 功能类别 → 颜色（与报告图表一致的口径）
CATEGORY_COLORS = {
    '聚合酶/复制相关': '#2e86c1',
    '聚合酶辅因子': '#5dade2',
    '蛋白酶/加工': '#8e44ad',
    '结构蛋白': '#e67e22',
    '运动蛋白': '#27ae60',
    '宿主互作/致病': '#c0392b',
    '基因组复制/维持': '#16a085',
    '转录/翻译调控': '#d4ac0d',
    '其他功能蛋白': '#7f8c8d',
    '假想蛋白（功能未定）': '#b2babb',
}
DOMAIN_COLOR = '#2471a3'
DOMAIN_COLOR_ALT = '#7fb3d5'

TICK = 500          # 比例尺刻度（nt）
MAX_WIDTH = 900.0   # 图最大宽（px）
ORF_H = 20.0        # ORF 矩形高
DOM_H = 7.0         # 结构域窄条高
ARROW = 12.0        # 箭头部长度
NAME_GAP = 10.0


def _esc(s):
    return html.escape(str(s or ''), quote=True)


def _fmt_int(n):
    return f"{int(n):,}"


def _safe_stem(name):
    """contig 名 → 安全文件名主干（防 ../ 与非法字符导致写出目录逃逸）。"""
    stem = re.sub(r'[^A-Za-z0-9_\-.]', '_', str(name or ''))
    stem = re.sub(r'\.{2,}', '_', stem).strip('.')
    return stem or 'contig'


def render_contig_svg(contig, length, orfs, out_path):
    """单条 contig 的线性基因组图。

    orfs: [{'start','end','strand','product','category','hmm_hits','cdd_hits'}, ...]
    （1 起始闭区间坐标，start/end 已由 orf_annotation.tsv 保证）
    """
    length = max(int(length or 0), 1)
    xscale = MAX_WIDTH / max(length, 1)
    left, top = 20.0, 58.0
    w = length * xscale
    line_y = top + 60.0

    parts = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w + 260:.0f}" '
        f'height="{line_y + 110:.0f}" viewBox="0 0 {w + 260:.0f} '
        f'{line_y + 110:.0f}" font-family="Arial,sans-serif">')
    parts.append(f'<text x="{left}" y="26" font-size="14" font-weight="bold">'
                 f'{_esc(contig)}</text>')
    parts.append(f'<text x="{left}" y="44" font-size="11" fill="#555">'
                 f'长度 {_esc(_fmt_int(length))} nt · '
                 f'{len(orfs)} 个 ORF</text>')

    # 基线
    parts.append(f'<line x1="{left}" y1="{line_y}" x2="{left + w:.1f}" '
                 f'y2="{line_y}" stroke="#333" stroke-width="1.5"/>')

    # 比例尺（线上方，逢双刻度标数字）
    sb_y = line_y - 34.0
    parts.append(f'<line x1="{left}" y1="{sb_y}" x2="{left + w:.1f}" '
                 f'y2="{sb_y}" stroke="#333"/>')
    nt = 0
    tick_i = 0
    while nt <= length + 0.001:
        x = left + nt * xscale
        emph = tick_i % 2 == 0
        h = 7.0 if emph else 4.0
        parts.append(f'<line x1="{x:.1f}" y1="{sb_y}" x2="{x:.1f}" '
                     f'y2="{sb_y - h}" stroke="#333"/>')
        if nt and (emph or nt >= length):
            parts.append(f'<text x="{x:.1f}" y="{sb_y - 10:.1f}" '
                         f'font-size="9" text-anchor="middle" fill="#333">'
                         f'{_esc(_fmt_int(nt))}</text>')
        nt += TICK
        tick_i += 1

    # ORF（按类别着色的方向箭头；+ 链在线上方，- 链在线下方）
    for i, o in enumerate(sorted(orfs, key=lambda x: (x.get('start') or 0))):
        s, e = int(o.get('start') or 0), int(o.get('end') or 0)
        if e <= s:
            continue
        strand = o.get('strand') or '+'
        x0 = left + (s - 1) * xscale
        bw = max((e - s + 1) * xscale, 4.0)
        color = CATEGORY_COLORS.get(o.get('category'), '#95a5a6')
        h = ORF_H
        if strand == '+':
            y = line_y - h
            pts = (f"{x0:.1f},{y:.1f} {x0 + max(bw - ARROW, bw * 0.85):.1f},"
                   f"{y:.1f} {x0 + bw:.1f},{y + h / 2:.1f} "
                   f"{x0 + max(bw - ARROW, bw * 0.85):.1f},{y + h:.1f} "
                   f"{x0:.1f},{y + h:.1f}")
        else:
            y = line_y
            xa = x0 + min(ARROW, bw * 0.15)
            pts = (f"{xa:.1f},{y:.1f} {x0 + bw:.1f},{y:.1f} "
                   f"{x0 + bw:.1f},{y + h:.1f} {xa:.1f},{y + h:.1f} "
                   f"{x0:.1f},{y + h / 2:.1f}")
        tip = (f"{_esc(o.get('product') or '(未注释)')} · "
               f"{_esc(o.get('category') or '')} · {s}-{e}({strand})")
        for hh in (o.get('hmm_hits'), o.get('cdd_hits')):
            if hh:
                tip += f"\n{_esc(hh)}"
        parts.append(f'<polygon points="{pts}" fill="{color}" '
                     f'fill-opacity="0.92" stroke="#2c3e50" '
                     f'stroke-width="0.6"><title>{tip}</title></polygon>')
        label = o.get('product') or ''
        if label and bw > 52:
            short = label if len(label) <= 22 else label[:20] + '…'
            parts.append(f'<text x="{x0 + bw / 2:.1f}" '
                         f'y="{y + h / 2 + 3.4:.1f}" font-size="8.5" '
                         f'text-anchor="middle" fill="#ffffff">'
                         f'{_esc(short)}</text>')

    # 结构域窄条（叠放于基线上下外侧，交替层；littlegenomes 的 s 宽度思路）
    for i, o in enumerate(sorted(orfs, key=lambda x: (x.get('start') or 0))):
        hh = o.get('hmm_hits') or o.get('cdd_hits')
        if not hh:
            continue
        s, e = int(o.get('start') or 0), int(o.get('end') or 0)
        if e <= s:
            continue
        x0 = left + (s - 1) * xscale
        bw = max((e - s + 1) * xscale, 3.0)
        above = i % 2 == 0
        y = (line_y - ORF_H - 10 - DOM_H) if above else (line_y + ORF_H + 10)
        color = DOMAIN_COLOR if above else DOMAIN_COLOR_ALT
        parts.append(f'<rect x="{x0:.1f}" y="{y:.1f}" width="{bw:.1f}" '
                     f'height="{DOM_H}" fill="{color}" stroke="#1a5276" '
                     f'stroke-width="0.4"><title>{_esc(hh)}</title></rect>')

    parts.append('</svg>')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))
    return out_path


def render_sample_diagrams(rows, out_dir, contig_lengths=None, max_plots=12,
                           logger=None):
    """按 contig 分组渲染示意图。rows 为 orf_annotation 行 dict 列表。
    返回 [相对文件名, ...]。"""
    contig_lengths = dict(contig_lengths or {})
    by_contig = {}
    for r in rows:
        cid = r.get('contig')
        if not cid:
            continue
        by_contig.setdefault(cid, []).append(r)
    # 长度缺省用该 contig 最大 end
    for cid, rs in by_contig.items():
        if cid not in contig_lengths:
            contig_lengths[cid] = max((r.get('end') or 0) for r in rs)
    # 出图顺序：ORF 数多者优先，限 max_plots
    ordered = sorted(by_contig, key=lambda c: -len(by_contig[c]))[:max_plots]
    out = []
    for cid in ordered:
        try:
            rel = os.path.join('genome_diagrams', f'{_safe_stem(cid)}.svg')
            render_contig_svg(cid, contig_lengths.get(cid),
                              by_contig[cid], os.path.join(out_dir, rel))
            out.append(rel)
        except Exception as e:
            if logger:
                logger.log(f"基因组示意图失败 {cid}: {e}", "WARN")
    if logger:
        logger.log(f"基因组示意图 {len(out)} 张（{out_dir}）")
    return out

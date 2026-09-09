# -*- coding: utf-8 -*-
"""
阶段⑨ 病毒基因组图（dna_features_viewer 引擎，gbdraw 的纯 Python 备选引擎）：
- gbdraw 圈图更精美，仍是默认首选；本引擎在其不可用时自动顶上（无需外部 CLI，
  pip install dna_features_viewer 即可），也可由 plot_engine='dfv' 强制选用
- 输入约定与 gbdraw 一致：③病毒 contigs FASTA + GFF 注释；注释优先取 ⑥b 的
  orf_annotation.gff3（带 product/category/family，可按功能类别着色），回退
  ⑥ 的 pyrodigal.gff
- standalone 模式：fasta_in + ann_in（GFF3），或 GenBank .gb/.gbk
  （BiopythonTranslator 翻译特征，无需 fasta）
- 输出与 gbdraw 同目录同 summary 格式（09_genome_plots/summary.json，plots 列表
  供 ⑩报告内嵌）：每条 contig 一张 <contig>.dfv_linear.svg + 一张 .dfv_circle.svg
"""
import os
import json

from .utils import (check_path, safe_open, iter_fasta, is_step_done,
                    mark_step_done)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from dna_features_viewer import (GraphicFeature, GraphicRecord,
                                     CircularGraphicRecord)
    _DFV_IMPORT_OK = True
except Exception:
    _DFV_IMPORT_OK = False

# 功能类别 → 颜色（与 ⑥b orf_annotation.gff3 的 category 属性对应）
CATEGORY_COLORS = {
    '结构蛋白': '#7ec8e3',
    '聚合酶/复制相关': '#ffb347',
    '运动蛋白': '#b5e7a0',
}
_COLOR_FALLBACK = '#d3d3d3'


def dfv_available():
    return _DFV_IMPORT_OK


# ---------------------------------------------------------------- GFF 解析
def parse_gff(path, seqid=None):
    """平台 GFF3 → {contig: [feat_dict, ...]}；feat 坐标为 1-based 闭区间。

    seqid 给定时只保留该 contig 的记录。
    """
    by_ctg = {}
    with safe_open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            attrs = {}
            for kv in cols[8].split(';'):
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    attrs[k.strip()] = v
            feat = {
                'type': cols[2],
                'start': int(cols[3]), 'end': int(cols[4]),
                'strand': 1 if cols[6] != '-' else -1,
                'product': attrs.get('product', ''),
                'category': attrs.get('category', ''),
                'organism': attrs.get('organism', ''),
                'id': attrs.get('ID', ''),
            }
            by_ctg.setdefault(cols[0], []).append(feat)
    if seqid is not None:
        return {seqid: by_ctg.get(seqid, [])}
    return by_ctg


def _pick_inputs(sample_dir, logger=None):
    """默认输入：③ 的病毒 contigs FASTA + ⑥b/⑥ 的 GFF 注释。"""
    fa = os.path.join(sample_dir, '03_assembly', 'viral_contigs.fasta')
    if not os.path.isfile(fa):
        fa = os.path.join(sample_dir, '03_assembly', 'contigs.filtered.fasta')
    if not os.path.isfile(fa):
        raise FileNotFoundError(
            "③组装 目录缺少 viral_contigs.fasta / contigs.filtered.fasta，"
            "请先运行组装步骤（或在卡片参数中指定自备 FASTA 路径）")
    cands = [os.path.join(sample_dir, '04b_orf_annot', 'orf_annotation.gff3'),
             os.path.join(sample_dir, '04_orf', 'pyrodigal.gff')]
    gff = next((p for p in cands if os.path.isfile(p)), None)
    return fa, gff


# ---------------------------------------------------------------- 绘图
def _label(feat):
    if feat.get('product'):
        lab = feat['product']
        if feat.get('category'):
            lab += ' [%s]' % feat['category']
        return lab
    return feat.get('id') or feat.get('type') or 'feature'


def _graphic_features(feats):
    return [GraphicFeature(
        start=max(f['start'] - 1, 0), end=f['end'], strand=f['strand'],
        label=_label(f),
        color=CATEGORY_COLORS.get(f.get('category'), _COLOR_FALLBACK))
        for f in sorted(feats, key=lambda x: x['start'])]


def _setup_fonts():
    """中文标签字体（平台 Windows 环境用雅黑/黑体，其余回退 Noto/Arial）。"""
    old = (plt.rcParams['font.sans-serif'], plt.rcParams['axes.unicode_minus'])
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei',
                                       'Noto Sans CJK SC', 'Arial']
    plt.rcParams['axes.unicode_minus'] = False
    return old


def plot_contig_svgs(seqid, feats, seq_len, out_prefix, organism=''):
    """单条 contig 出线性 + 环形两张 SVG，返回生成文件列表。"""
    if not _DFV_IMPORT_OK:
        raise RuntimeError(
            "未检测到 dna_features_viewer（安装: python -m pip install "
            "dna_features_viewer）")
    old = _setup_fonts()
    made = []
    try:
        gf = _graphic_features(feats)
        title = '%s%s' % (seqid, (' · %s · %d bp · %d 特征'
                                  % (organism, seq_len, len(gf))) if organism
                          else ' · %d bp' % seq_len)
        for circular in (False, True):
            cls = CircularGraphicRecord if circular else GraphicRecord
            record = cls(sequence_length=max(seq_len, 1), features=gf)
            fig, ax = plt.subplots(
                figsize=(8, 8) if circular else (11, 2 + min(len(gf), 6) * 0.6))
            record.plot(ax=ax, with_ruler=not circular,
                        draw_line=True, annotate_inline=True)
            ax.set_title(title, fontsize=10)
            svg = '%s.%s' % (out_prefix,
                             'dfv_circle' if circular else 'dfv_linear') + '.svg'
            fig.savefig(svg, bbox_inches='tight')
            plt.close(fig)
            made.append(svg)
    finally:
        plt.rcParams['font.sans-serif'], plt.rcParams['axes.unicode_minus'] = old
    return made


def _gb_record_svgs(record, out_prefix):
    """单条 GenBank record → 线性 + 环形两张 SVG。"""
    old = _setup_fonts()
    made = []
    try:
        for circular in (False, True):
            cls = CircularGraphicRecord if circular else GraphicRecord
            from dna_features_viewer import BiopythonTranslator
            translator = BiopythonTranslator()
            translator.ignored_features_types = ('source', 'gene')
            rec = translator.translate_record(record, record_class=cls)
            fig, ax = plt.subplots(
                figsize=(8, 8) if circular else (11, 5))
            rec.plot(ax=ax, with_ruler=not circular,
                     annotate_inline=True)
            svg = '%s.%s' % (out_prefix,
                             'dfv_circle' if circular else 'dfv_linear') + '.svg'
            fig.savefig(svg, bbox_inches='tight')
            plt.close(fig)
            made.append(svg)
    finally:
        plt.rcParams['font.sans-serif'], plt.rcParams['axes.unicode_minus'] = old
    return made


def _safe_name(name):
    from .gbdraw_plot import _safe_fn
    return _safe_fn(name)


# ---------------------------------------------------------------- 主入口
def run_dfv_plots(sample_dir, logger=None, force=False, max_plots=12,
                  fasta_in=None, ann_in=None, progress=None):
    """dna_features_viewer 出图主入口（参数与返回值同 gbdraw_plot.run_genome_plots）。

    fasta_in/ann_in 给定时走 standalone（自备文件）模式：ann_in 为 GFF3
    （与 fasta_in 配对）或 GenBank .gb/.gbk（可多条记录，逐 record 出图）。
    返回 summary dict（plots: 生成文件相对路径列表，嵌入报告用）。
    """
    step = 'gbdraw'
    out_dir = check_path(os.path.join(sample_dir, '09_genome_plots'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    summary_file = os.path.join(out_dir, 'summary.json')
    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("阶段⑨基因组图已完成，跳过")
        with safe_open(summary_file) as f:
            return json.load(f)
    if not dfv_available():
        raise RuntimeError(
            "未检测到 dna_features_viewer（安装: python -m pip install "
            "dna_features_viewer）")

    plots = []
    if ann_in and str(ann_in).lower().endswith(
            ('.gb', '.gbk', '.gbff', '.genbank')):
        # GenBank 直接出图（逐 record）
        from Bio import SeqIO
        n = 0
        for rec in SeqIO.parse(check_path(ann_in, must_exist=True), 'genbank'):
            if max_plots and n >= max_plots:
                break
            prefix = os.path.join(out_dir, _safe_name(rec.id or rec.name))
            if progress:
                progress(n / max(max_plots, 1) * 0.95, f'出图 {rec.id}')
            plots += _gb_record_svgs(rec, prefix)
            n += 1
        summary = {'stage': step, 'engine': 'dna_features_viewer',
                   'standalone': True,
                   'input': os.path.basename(str(ann_in)), 'plots': plots,
                   'n_seqs': n}
    else:
        fa = check_path(fasta_in, must_exist=True) if fasta_in else \
            _pick_inputs(sample_dir)[0]
        gff = check_path(ann_in, must_exist=True) if ann_in else \
            (None if fasta_in else _pick_inputs(sample_dir)[1])
        if logger and gff:
            logger.log(f"注释文件: {gff}")
        by_ctg = parse_gff(gff) if gff else {}
        if fasta_in and not gff and logger:
            logger.log("未提供注释文件：画裸骨架图（无特征）")
        seqs = sorted(((h.split()[0], s) for h, s in iter_fasta(fa)),
                      key=lambda x: -len(x[1]))[:max_plots] if max_plots else \
            sorted(((h.split()[0], s) for h, s in iter_fasta(fa)),
                   key=lambda x: -len(x[1]))
        for i, (cid, seq) in enumerate(seqs):
            if progress:
                progress(i / max(len(seqs), 1) * 0.95,
                         f'出图 {cid}（{i + 1}/{len(seqs)}）')
            feats = by_ctg.get(cid, [])
            organism = next((f['organism'] for f in feats if f['organism']), '')
            prefix = os.path.join(out_dir, _safe_name(cid))
            try:
                plots += plot_contig_svgs(cid, feats, len(seq), prefix,
                                          organism=organism)
            except Exception as e:
                if logger:
                    logger.log(f"{cid} 出图失败（跳过）: {e}", "WARN")
        summary = {'stage': step, 'engine': 'dna_features_viewer',
                   'standalone': bool(fasta_in or ann_in),
                   'input': os.path.basename(str(fa)), 'plots': plots,
                   'n_seqs': len(seqs)}
    # 清掉列表里实际不存在的文件（出图失败被跳过时）
    plots[:] = [p for p in plots if os.path.isfile(p)]
    with safe_open(summary_file, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)
    mark_step_done(out_dir, step)
    if logger:
        logger.log(f"阶段⑨ 完成（dna_features_viewer）: {len(plots)} 张"
                   f"基因组图 -> 09_genome_plots/")
    return summary

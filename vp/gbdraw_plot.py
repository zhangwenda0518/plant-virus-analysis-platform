# -*- coding: utf-8 -*-
"""
阶段⑨ 病毒基因组图（gbdraw）：
- 默认对 ③组装 的病毒 contigs 出图：每条 contig 一张圈图 + 一张线图（SVG）
- 注释优先用 ⑥ORF 的 pyrodigal GFF3（按 seqid 拆分）；无 ORF 时画裸骨架
- standalone 模式：fasta_in + ann_in（GFF3 或 GenBank .gb/.gbk）直接出图，
  支持用户自备文件（fasta + gb 即可）
输出: 09_genome_plots/<contig>.circle.svg / .linear.svg + summary.json
报告嵌入: SVG 内嵌 ⑪可视化报告
"""
import os
import sys
import re
import json
import shutil
import tempfile

from .config import get_config
from .utils import (check_path, safe_open, iter_fasta, write_fasta_record,
                    is_step_done, mark_step_done, log_res_plan)

_FA_EXT = ('.fasta', '.fa', '.fna', '.fas')
_ANN_EXT = ('.gff', '.gff3', '.gb', '.gbk', '.gbff', '.genbank')


def gbdraw_available():
    try:
        get_config().tool('gbdraw')
        return True
    except Exception:
        return False


def _write_skeleton_gff(fasta, out_dir=None):
    """纯 FASTA（无注释）输入：合成只含 region 特征的 GFF3。

    gbdraw CLI 要求 --fasta 必须配对 --gff，且 GFF 无任何 feature 记录时
    报 "No valid records"；一条跨全长的 region 特征即可让它画出裸骨架。
    GFF 写入 out_dir（缺省 = fasta 旁边；输入在平台外时应显式给运行目录）。
    """
    gff = os.path.join(out_dir or os.path.dirname(fasta),
                       _safe_fn(os.path.basename(fasta)) + '.skeleton.gff')
    with safe_open(gff, 'wt') as f:
        f.write('##gff-version 3\n')
        for header, seq in iter_fasta(fasta):
            f.write(f"{header.split()[0]}\tinput\tregion\t1\t{len(seq)}\t"
                    f".\t+\t.\tID=skeleton\n")
    return gff


def _run_gbdraw(fasta=None, gff=None, gbk=None, out_prefix=None, logger=None,
                mode='both', opts=None):
    """调用 gbdraw 子命令。mode: circular / linear / both。返回生成的文件列表。

    CLI 用法: gbdraw circular|linear [--gbk f.gb | (--gff f.gff --fasta f.fa)]
              [--overwrite] [-o 前缀]（子命令必须在最前）
    gbdraw 输出固定为 <前缀>.svg（不带模式后缀）——同前缀两模式会互相覆盖，
    因此内部对每个模式追加 .circular/.linear 后缀再调用。

    opts: dict {gbdraw CLI 参数名(如 'palette'/'track_type'/'separate_strands'): 值}。
    值规则：True→只输出 flag（--key）；str/int/False→--key 值；None/''→跳过。
    """
    gbdraw = get_config().tool('gbdraw')
    out_dir = os.path.dirname(out_prefix)
    os.makedirs(out_dir, exist_ok=True)
    made = []
    tail = ['--overwrite']
    if gbk:
        tail += ['--gbk', gbk]
    else:
        if not gff:
            gff = _write_skeleton_gff(fasta, out_dir=out_dir)
        tail += ['--gff', gff, '--fasta', fasta]
    # 透传绘图定制参数（调色板/布局/形态/标签/GC 等）
    if opts:
        for key, val in opts.items():
            if val is None or val == '':
                continue
            if val is True:
                tail.append('--' + key)
            else:
                tail += ['--' + key, str(val)]
    bundled = (gbdraw == 'bundled:gbdraw')
    for sub in ([mode] if mode in ('circular', 'linear')
                else ['circular', 'linear']):
        sub_prefix = f'{out_prefix}.{sub}'
        if bundled:
            _gbdraw_inproc([sub, '-o', sub_prefix] + tail)
        else:
            from .utils import run_cmd
            run_cmd([gbdraw, sub, '-o', sub_prefix] + tail, logger=logger)
        # 产物扩展名随 -f format 变化（svg/png/pdf/ps/eps），用 glob 收集任意格式
        for _ext in ('.svg', '.png', '.pdf', '.ps', '.eps'):
            _f = sub_prefix + _ext
            if os.path.isfile(_f):
                made.append(_f)
                break
    return made


def _gbdraw_inproc(argv):
    """冻结分发：进程内执行 gbdraw.cli（随包引擎，无独立 exe）。
    需 VirusPlatform.spec 中 collect_submodules('gbdraw') 与其数据文件。"""
    import contextlib
    import io as _io
    from gbdraw.cli import main
    old_argv = sys.argv
    sys.argv = ['gbdraw'] + [str(a) for a in argv]
    try:
        with contextlib.redirect_stdout(_io.StringIO()) as _out,                 contextlib.redirect_stderr(_io.StringIO()) as _err:
            main()
    except SystemExit as e:                      # CLI 正常退出路径
        if e.code not in (0, None):
            raise RuntimeError(f'gbdraw 退出码 {e.code}') from None
    except Exception as e:
        raise RuntimeError(f'gbdraw 运行失败: {e}') from e
    finally:
        sys.argv = old_argv


def _split_by_seqid(src_fasta, src_gff, out_dir, seqid):
    """从多条序列的 FASTA+GFF3 中拆出单条 contig 的小文件（gbdraw 单图输入）。"""
    fa_out = os.path.join(out_dir, f'{_safe_fn(seqid)}.fa')
    gff_out = os.path.join(out_dir, f'{_safe_fn(seqid)}.gff')
    found = False
    with safe_open(fa_out, 'wt') as f:
        for header, seq in iter_fasta(src_fasta):
            if header.split()[0] == seqid:
                write_fasta_record(f, seqid, seq)
                found = True
                break
    if not found:
        return None, None
    if src_gff and os.path.isfile(src_gff):
        with safe_open(gff_out, 'wt') as f:
            f.write('##gff-version 3\n')
            with safe_open(src_gff) as g:
                for line in g:
                    if line.startswith('#'):
                        continue
                    cols = line.rstrip('\n').split('\t')
                    if len(cols) > 3 and cols[0] == seqid:
                        f.write(line if line.endswith('\n') else line + '\n')
        return fa_out, gff_out
    return fa_out, None


def _safe_fn(name):
    return re.sub(r'[^A-Za-z0-9_\-.]+', '_', str(name))[:80] or 'seq'


def _pick_inputs(sample_dir, logger=None):
    """默认输入：③ 的病毒 contigs FASTA + ⑥ 的 pyrodigal GFF 注释。"""
    fa = os.path.join(sample_dir, '03_assembly', 'viral_contigs.fasta')
    if not os.path.isfile(fa):
        fa = os.path.join(sample_dir, '03_assembly', 'contigs.filtered.fasta')
    if not os.path.isfile(fa):
        raise FileNotFoundError(
            "③组装 目录缺少 viral_contigs.fasta / contigs.filtered.fasta，"
            "请先运行组装步骤（或在卡片参数中指定自备 FASTA 路径）")
    gff = os.path.join(sample_dir, '04_orf', 'pyrodigal.gff')
    if not os.path.isfile(gff):
        gff = None
    return fa, gff


def run_genome_plots(sample_dir, logger=None, force=False, max_plots=12,
                     fasta_in=None, ann_in=None, progress=None, opts=None,
                     mode='both'):
    """出图主入口。fasta_in/ann_in 给定时走 standalone（自备文件）模式。

    ann_in 可为 GFF3（与 fasta_in 配对）或 GenBank（.gb/.gbk，无需 fasta）。
    opts: dict {gbdraw CLI 参数名: 值}，透传给 _run_gbdraw（调色板/布局/形态/标签等）。
    mode: 'circular' / 'linear' / 'both'（默认 both 出圈图+线图）。
    返回 summary dict（plots: 生成文件列表，嵌入报告用）。
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
    if not gbdraw_available():
        raise RuntimeError(
            "未检测到 gbdraw（安装: python -m pip install git+"
            "https://github.com/satoshikawato/gbdraw.git）")

    if mode not in ('circular', 'linear'):
        mode = 'both'
    modes = (('circular', '圈图'), ('linear', '线图')) if mode == 'both' \
        else ((mode, '圈图' if mode == 'circular' else '线图'),)

    standalone = bool(fasta_in or ann_in)
    if standalone and ann_in and str(ann_in).lower().endswith(
            ('.gb', '.gbk', '.gbff', '.genbank')):
        # GenBank 直接出图（可多条记录，gbdraw 自动多面板）
        plots = []
        for mode, tag in modes:
            made = _run_gbdraw(gbk=check_path(ann_in, must_exist=True),
                               out_prefix=os.path.join(out_dir, 'genome'),
                               logger=logger, mode=mode, opts=opts)
            plots += made
        seqs = []
        summary = {'stage': step, 'standalone': True,
                   'input': os.path.basename(str(ann_in)), 'plots': plots,
                   'n_seqs': 0}
    else:
        fa = check_path(fasta_in, must_exist=True) if fasta_in else \
            _pick_inputs(sample_dir, logger=logger)[0]
        gff = check_path(ann_in, must_exist=True) if ann_in else \
            (None if fasta_in else _pick_inputs(sample_dir, logger=logger)[1])
        if fasta_in and not gff:
            if logger:
                logger.log("未提供注释文件：画裸骨架图（无 CDS 特征）")
        # 每条 contig 一组图；按长度降序最多 max_plots 条
        seqs = sorted(((h.split()[0], s) for h, s in iter_fasta(fa)),
                      key=lambda x: -len(x[1]))[:max_plots]
        tmp = check_path(os.path.join(out_dir, '_tmp'),
                         must_exist=False, in_platform=True)
        os.makedirs(tmp, exist_ok=True)
        plots = []
        try:
            for i, (cid, _seq) in enumerate(seqs):
                if progress:
                    progress(i / max(len(seqs), 1) * 0.95,
                             f'出图 {cid}（{i + 1}/{len(seqs)}）')
                fa1, gff1 = (fa, gff) if standalone and len(seqs) == 1 else \
                    _split_by_seqid(fa, gff, tmp, cid)
                if not fa1:
                    continue
                prefix = os.path.join(out_dir, _safe_fn(cid))
                try:
                    plots += _run_gbdraw(fasta=fa1, gff=gff1,
                                         out_prefix=prefix, logger=logger,
                                         opts=opts, mode=mode)
                except RuntimeError as e:
                    if logger:
                        logger.log(f"{cid} 出图失败（跳过）: {e}", "WARN")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        summary = {'stage': step, 'standalone': bool(fasta_in or ann_in),
                   'input': os.path.basename(str(fa)), 'plots': plots,
                   'n_seqs': len(seqs)}
    with safe_open(summary_file, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)
    mark_step_done(out_dir, step)
    if logger:
        logger.log(f"阶段⑨ 完成: {len(plots)} 张基因组图 -> 09_genome_plots/")
    return summary

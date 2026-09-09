# -*- coding: utf-8 -*-
"""
LoVis4u 子进程运行器（vp.synteny 的可选外接引擎）。

lovis4u 0.2.0 可直接 pip install（Windows 也有轮子），但其
DataProcessing.mmseqs_cluster() 存在 Windows 临时文件兼容 bug：
tempfile.NamedTemporaryFile() 保持句柄打开，紧接着 SeqIO.write 重新
打开同一路径 → PermissionError（POSIX 允许重开已打开文件，故官方
只在 macOS/Linux 测过）。本运行器用最小补丁（创建后立即关闭句柄，
delete=False）修正该处行为，其余完全走官方 CLI 流程；mmseqs 二进制
通过 -smp 指向平台内置 mmseqs/bin/mmseqs.exe。

用法: python vp/lovis4u_run.py <gb清单文件> <out_dir> <pdf_name> <mmseqs_exe>
      <gb清单文件> 为文本文件，每行一个 .gb/.gbk 路径。
退出码 0=成功；非 0=失败（stderr 含 traceback）。
"""
import os
import sys
import shutil
import glob
import tempfile
import traceback


def _patch_tempfile():
    """关闭 NamedTemporaryFile 的句柄（lovis4u 只用 .name，Windows 下
    打开的句柄会阻止 SeqIO/mmseqs 再打开同一路径）。
    注意：DataProcessing.py:1384 的 bedgraph 临时文件同样受益；本运行器
    不使用 bedgraph/bigwig 覆盖度轨道。"""
    _ntf = tempfile.NamedTemporaryFile

    def _patched(*args, **kwargs):
        kwargs.setdefault('delete', False)
        f = _ntf(*args, **kwargs)
        f.close()
        return f

    tempfile.NamedTemporaryFile = _patched


def main(argv):
    list_file, out_dir, pdf_name, mmseqs = argv[1:5]
    # 第 5 个位置之后的参数为可选 lovis4u CLI flags（如 --set-category-colour、
    # --run-hmmscan、-gc、-gc_skew、-align、--set-group-colour-for、-fv-off），
    # 与 gbdraw 透传方式一致，原样追加到 sys.argv。
    extra_flags = argv[5:]
    _patch_tempfile()

    work = os.path.join(out_dir, 'lovis4u_in')
    os.makedirs(work, exist_ok=True)
    n = 0
    with open(list_file, encoding='utf-8') as f:
        for line in f:
            src = line.strip()
            if not src:
                continue
            if not src.lower().split('.gz')[0].endswith(
                    ('.gb', '.gbk', '.gbff', '.genbank')):
                continue
            shutil.copy(src, work)
            n += 1
    if n < 2:
        raise RuntimeError(f'可用 GenBank 文件不足 2 个: {list_file}')

    import lovis4u
    sys.argv = ['lovis4u', '-gb', work, '-o', out_dir,
                '--pdf-name', pdf_name,
                '-hl',                     # 同源连线带
                '-fw', '260',              # 画布宽 260mm（默认太窄，标签被裁）
                '--debug'] + extra_flags
    parameters = lovis4u.Manager.Parameters()
    parameters.parse_cmd_arguments()
    parameters.load_config(parameters.cmd_arguments['config_file'])
    # 不走 -smp：其实现把路径塞进 re.sub 替换模板，Windows 反斜杠路径会
    # 触发 bad escape（lovis4u 的第二个 Windows 兼容问题）；配置加载后
    # 直接覆盖最终参数即可。
    parameters.args['mmseqs_binary'] = mmseqs

    loci = lovis4u.DataProcessing.Loci(parameters=parameters)
    loci.load_loci_from_gb(parameters.args['gb'])

    if parameters.args['mmseqs']:
        results = loci.mmseqs_cluster()
        loci.define_feature_groups(results)
        if parameters.args['clust_loci']:
            loci.cluster_sequences(results,
                                   one_cluster=parameters.args['one_cluster'])
        if parameters.args['find-variable']:
            loci.find_variable_feature_groups(results)
    if parameters.args['align_loci']:
        loci.auto_align_loci()
    if parameters.args['reorient_loci']:
        loci.reorient_loci()
    if parameters.args['run_hmmscan_search']:
        loci.pyhmmer_annotation()
    if parameters.args['set-group-colour']:
        loci.set_feature_colours_based_on_groups()
    if parameters.args['set-category-colour']:
        loci.set_category_colours()
    loci.define_labels_to_be_shown()
    loci.save_feature_annotation_table()
    loci.save_locus_annotation_table()

    canvas = lovis4u.Manager.CanvasManager(parameters)
    canvas.define_layout(loci)
    canvas.add_loci_tracks(loci)
    if parameters.args['draw_scale_line_track']:
        canvas.add_scale_line_track()
    if parameters.args['set-category-colour']:
        canvas.add_categories_colour_legend_track(loci)
    if parameters.args['homology-track']:
        canvas.add_homology_track()
    canvas.plot(filename=parameters.args['pdf-name'])
    pdf_path = os.path.join(out_dir, parameters.args['pdf-name'])
    if not os.path.isfile(pdf_path):
        raise RuntimeError('lovis4u 未生成 PDF')


if __name__ == '__main__':
    try:
        main(sys.argv)
    except Exception:
        traceback.print_exc()
        sys.exit(1)

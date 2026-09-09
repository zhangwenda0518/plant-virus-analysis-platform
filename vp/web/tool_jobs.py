# -*- coding: utf-8 -*-
"""工具箱：22 个独立分析任务工厂（自 app.py 拆出）。

每个 _tool_job_<name>(ctx) 返回 job(log, prog, cancel) 可调用；
由 vp/web/tools_api.py 的 /api/tool/run 取用。
新增工具三步：① 写 _tool_job_<name> ② 注册 TOOL_REGISTRY
（在 tools_api.py）③ webapp/templates/tools.html 加卡片。
"""
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid

from flask import abort, jsonify, render_template, request, send_file, \
    send_from_directory

from vp.config import DIRS, PLATFORM_ROOT, db_path, engine_cmd
from vp.utils import (TaskLogger, check_path, fmt_size, run_cmd, safe_open,
                      safe_remove)
from vp.web.state import cfg, tool_runs_root as _tool_runs_root


CONSENSUS_READS_LAYERS = {
    'host_removed': '01_host_removal',
    'viral':        '02_virus_screen',
    'clean':        '00_prep',
}


def _tool_job_convert(ctx):
    """格式转换：.sra→FASTQ/FASTA（sracha）、FASTQ→FASTA（seqkit）。"""
    inp = ctx.req('input', '输入文件')
    target = ctx.p.get('target') or 'fastq'
    if target not in ('fastq', 'fasta'):
        abort(400, '无效的目标格式')
    ilower = str(inp).lower()
    if ilower.endswith('.sra'):
        mode = 'sra'
    elif ilower.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        mode = 'fastq2fasta' if target == 'fasta' else None
    elif ilower.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
        abort(400, '输入已是 FASTA，无需转换')
    else:
        abort(400, '无法识别的输入格式（支持 .sra / .fastq[.gz]）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        out_files = []
        if mode == 'sra':
            from vp.public_data import sra_convert_engine
            engine, exe = sra_convert_engine()
            if engine != 'sracha':
                raise RuntimeError('sracha.exe 不可用，无法转换 .sra')
            threads = ctx.threads or max(2, min((os.cpu_count() or 4) // 2, 8))
            prog('sra', 0.2, f'sracha 解码 .sra → {target.upper()}')
            cmd = [exe, 'fastq', inp, '-O', ctx.run_dir,
                   '-t', str(threads), '-f', '-q']
            if target == 'fasta':
                cmd.append('--fasta')
            import subprocess as _sp
            _sp.run(cmd, check=True, timeout=8 * 3600,
                    stdout=_sp.DEVNULL, stderr=_sp.PIPE)
            base = os.path.basename(inp)[:-4]
            out_files = sorted(
                os.path.join(ctx.run_dir, f_) for f_ in os.listdir(ctx.run_dir)
                if f_.startswith(base + '_')
                and f_.endswith(('.fastq.gz', '.fq.gz', '.fa.gz', '.fasta.gz')))
            if not out_files:
                raise RuntimeError('sracha 无输出')
        else:
            seqkit = cfg.tool('seqkit')
            dst = os.path.join(ctx.run_dir,
                               os.path.basename(inp).rsplit('.', 2)[0] + '.fa.gz')
            prog('fq2fa', 0.3, 'seqkit fq2fa 转换中')
            from vp.utils import run_cmd
            run_cmd([seqkit, 'fq2fa', '-w', '0',
                     '-j', str(ctx.threads or cfg.threads), inp,
                     '-o', dst], logger=logger)
            out_files = [dst]
        prog('done', 1.0, f'完成：{len(out_files)} 个文件')
        logger.close()
        return {'n_files': len(out_files),
                'files': [os.path.relpath(f, ctx.run_dir).replace(os.sep, '/')
                          for f in out_files]}
    return job


def _tool_job_fastp(ctx):
    """① 质控预处理（fastp，单/双端）。"""
    r1 = ctx.req('r1', 'R1 FASTQ')
    r2 = ctx.opt('r2')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.preprocess import run_fastp
        prog('fastp', 0.3, 'fastp 质控中')
        res = run_fastp(ctx.run_dir, r1, r2, threads=ctx.threads, logger=logger,
                        force=True, dedup=bool(ctx.p.get('dedup')))
        if not res:
            raise RuntimeError('未检测到 fastp.exe，无法运行质控')
        prog('fastp', 1.0, '完成')
        logger.close()
        return res
    return job


def _tool_job_hostremoval(ctx):
    """宿主去除与序列提取（kunpeng 宿主库分类，独立模块，不依赖样品管道）。

    C 行 = 宿主 read 对，剔除后保留非宿主 reads（kept_R1/R2.fastq.gz）。
    """
    r1 = ctx.req('r1', 'R1 FASTQ')
    r2 = ctx.opt('r2')
    conf = float(ctx.p.get('confidence') or 0)
    db_host = ctx.opt('db') or cfg.databases['host']

    from vp.kunpeng import db_ready
    if not db_ready(db_host):
        abort(400, '宿主库未就绪，请先到「数据库构建」页构建宿主库')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.host_removal import remove_host
        prog('classify', 0.05, 'kunpeng 宿主库分类中')
        res = remove_host(ctx.run_dir, r1, r2, db_host, threads=ctx.threads,
                          confidence=conf, logger=logger, force=True,
                          progress=lambda pct, msg: prog(
                              'classify' if pct < 0.9 else 'filter', pct, msg))
        prog('done', 1.0, '完成')
        logger.close()
        return res
    return job


def _tool_job_hostpredict(ctx):
    """宿主预测（ICTV 级联 + NCBI 元数据交叉），独立模块。

    输入 = 病毒 contig 分类表 TSV：管道③ virus_contigs.tsv 或
    工具④ virus_classification.tsv（列名自动识别，后者归一为 ③ 口径），
    可选配套 viral_contigs.fasta。
    """
    import csv
    import shutil
    tsv = ctx.req('tsv', '病毒 contig 分类表 TSV')
    fa = ctx.opt('fasta')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        a_dir = check_path(os.path.join(ctx.run_dir, '03_assembly'),
                           must_exist=False, in_platform=True)
        os.makedirs(a_dir, exist_ok=True)
        prog('prep', 0.05, '整理输入')
        norm = os.path.join(a_dir, 'virus_contigs.tsv')
        with safe_open(tsv) as f:
            header = f.readline().rstrip('\r\n').split('\t')
        if 'kunpeng_taxid' in header:
            shutil.copyfile(tsv, norm)
        else:
            cols = ['contig', 'length', 'kunpeng_flag', 'kunpeng_taxid',
                    'kunpeng_species', 'blast_top_hit', 'blast_identity(%)',
                    'blast_coverage_hsp(%)', 'blast_aln_len', 'blast_species',
                    'blast_family']
            with safe_open(tsv) as f, safe_open(norm, 'wt') as w:
                w.write('\t'.join(cols) + '\n')
                for r in csv.DictReader(f, delimiter='\t'):
                    kt = str(r.get('taxid') or r.get('kunpeng_taxid')
                             or '').strip()
                    w.write('\t'.join([
                        str(r.get('contig') or '').strip(),
                        str(r.get('length') or '').strip(),
                        'C' if kt.isdigit() and int(kt) > 0 else 'U',
                        kt,
                        str(r.get('taxon') or r.get('species') or '').strip(),
                        '', '', '', '', '', '']) + '\n')
        if fa:
            shutil.copyfile(fa, os.path.join(a_dir, 'viral_contigs.fasta'))
        from vp.host_analysis import predict_hosts
        prog('predict', 0.15, 'ICTV 宿主概率级联预测')
        res = predict_hosts(ctx.run_dir, threads=ctx.threads, logger=logger,
                            force=True)
        prog('done', 1.0, f"完成：宿主判定 {res.get('n_contigs', 0)} 条")
        logger.close()
        return res
    return job


def _tool_job_orf(ctx):
    """ORF 预测（pyrodigal / pyrodigal_rv），可选 ⑥b 功能注释。"""
    fasta = ctx.req('fasta', '输入 FASTA（核酸 contigs / 基因组）')
    min_aa = max(1, min(int(ctx.p.get('min_aa') or 100), 5000))
    annotate = bool(ctx.p.get('annotate'))
    orf_tool = (ctx.p.get('orf_tool') or '').strip()
    orf_engine = (ctx.p.get('engine') or '').strip()
    orf_db = (ctx.p.get('db') or '').strip()
    # 功能注释 CDS 模型：pyrodigal_rv（默认）/ pyrodigal
    orf_model = (ctx.p.get('model') or '').strip()
    if orf_model not in ('pyrodigal_rv', 'pyrodigal'):
        orf_model = ''

    def job(log, prog, cancel):
        import json as _json
        import shutil
        logger = TaskLogger(callback=log)
        a_dir = check_path(os.path.join(ctx.run_dir, '03_assembly'),
                           must_exist=False, in_platform=True)
        os.makedirs(a_dir, exist_ok=True)
        prog('prep', 0.03, '整理输入')
        from vp.utils import iter_fasta, write_fasta_record
        ids = []
        vfa = os.path.join(a_dir, 'viral_contigs.fasta')
        cfa = os.path.join(a_dir, 'contigs.filtered.fasta')
        with safe_open(vfa, 'wt') as w, safe_open(cfa, 'wt') as wc:
            for h, s in iter_fasta(fasta):
                cid = h.split()[0]
                ids.append(cid)
                write_fasta_record(w, cid, s)
                write_fasta_record(wc, cid, s)
        if not ids:
            raise RuntimeError('输入 FASTA 中没有序列')
        with safe_open(os.path.join(a_dir, 'summary.json'), 'wt') as f:
            _json.dump({'viral_contigs': ids, 'standalone': True}, f)
        from vp.orf import predict_orfs
        prog('orf', 0.08, 'ORF 基因预测')
        res = predict_orfs(ctx.run_dir, min_aa=min_aa, threads=ctx.threads,
                           logger=logger, force=True, tools=orf_tool or None,
                           progress=lambda p, m: prog(
                               'orf', 0.08 + p * 0.6, m))
        res = dict(res)
        if annotate:
            from vp.orf_annot import run_orf_annotation
            prog('orfa', 0.72, 'ORF 功能注释')
            res['orfa'] = run_orf_annotation(
                ctx.run_dir, threads=ctx.threads, logger=logger, force=True,
                engine=orf_engine or None, db=orf_db or None,
                model=orf_model or None,
                progress=lambda p, m: prog('orfa', 0.72 + p * 0.26, m))
        prog('done', 1.0, '完成')
        logger.close()
        return res
    return job


def _tool_job_orfa(ctx):
    """功能注释（独立模块）：对已有 orf_ 运行注释，或 FASTA 预测+注释一步完成。"""
    run = (ctx.p.get('run') or '').strip()
    orf_engine = (ctx.p.get('engine') or '').strip()
    orf_db = (ctx.p.get('db') or '').strip()
    # 功能注释 CDS 模型：pyrodigal_rv（默认）/ pyrodigal
    orf_model = (ctx.p.get('model') or '').strip()
    if orf_model not in ('pyrodigal_rv', 'pyrodigal'):
        orf_model = ''
    if run:
        if not re.fullmatch(r'[A-Za-z0-9_]+', run) or not run.startswith('orf_'):
            abort(400, f'无效的 ORF 运行名: {run}')
        target = check_path(os.path.join(_tool_runs_root(), run),
                            must_exist=True, in_platform=True)

        def job(log, prog, cancel):
            logger = TaskLogger(callback=log)
            from vp.orf_annot import run_orf_annotation
            prog('orfa', 0.15, f'对运行 {run} 做 ORF 功能注释')
            res = run_orf_annotation(target, threads=ctx.threads, logger=logger,
                                     force=True, engine=orf_engine or None,
                                     db=orf_db or None,
                                     model=orf_model or None,
                                     progress=lambda p, m: prog(
                                         'orfa', 0.15 + p * 0.8, m))
            res = dict(res)
            res['run'] = run
            prog('done', 1.0, '完成')
            logger.close()
            return res
        return job
    # 无 run → FASTA 输入：预测 + 注释一步完成
    ctx.p = dict(ctx.p)
    ctx.p['annotate'] = True
    if not (ctx.p.get('fasta') or '').strip():
        abort(400, '请选择已有 ORF 运行或输入 FASTA')
    return _tool_job_orf(ctx)


def _tool_job_genoplot(ctx):
    """基因组图谱（gbdraw 首选，缺则 DFV 顶上）。

    输入 FASTA（可选配 GFF3 注释）或 GenBank（.gb/.gbk，自带注释）。
    """
    fasta = ctx.opt('fasta')
    ann = ctx.opt('ann')
    if not fasta and not ann:
        abort(400, '请选择 FASTA 或 GenBank 输入')
    if ann and str(ann).lower().endswith(('.gb', '.gbk', '.gbff', '.genbank')):
        fasta = None          # GenBank 自带序列与注释，FASTA 忽略
    engine = ctx.p.get('engine') or 'auto'
    max_plots = max(1, min(int(ctx.p.get('max_plots') or 12), 200))
    # 绘图定制参数（透传 gbdraw CLI）：仅收集有值/True 的项
    gb_opts = {}
    for _k, _v in (ctx.p.get('gb_opts') or {}).items():
        if _v is None or _v == '' or _v is False:
            continue
        gb_opts[str(_k)] = _v
    plot_mode = (ctx.p.get('mode') or 'both')   # circular / linear / both

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.gbdraw_plot import run_genome_plots, gbdraw_available
        if not gbdraw_available():
            raise RuntimeError('未检测到 gbdraw（管道不支持 DFV）')
        tag = 'gbdraw'
        prog(tag, 0.05, f'{tag} 出图中')
        res = run_genome_plots(ctx.run_dir, logger=logger, force=True,
                               max_plots=max_plots, fasta_in=fasta, ann_in=ann,
                               opts=gb_opts, mode=plot_mode,
                               progress=lambda p, m: prog(tag, 0.05 + p * 0.9, m))
        if not res.get('plots'):
            raise RuntimeError('未产出任何基因组图（检查输入文件与绘图引擎）')
        res['run'] = os.path.basename(ctx.run_dir)   # 供前端 toolrun 内联展示 SVG
        prog('done', 1.0, f"完成：{len(res.get('plots', []))} 张图")
        logger.close()
        return res
    return job


def _tool_job_primer(ctx):
    """引物设计（primer3）。plain=基因组/contigs 全长分窗；
    conserved=多序列比对 FASTA 保守区（输入需已比对，如 MAFFT aln.fasta）。"""
    import shutil
    fasta = ctx.req('fasta', '输入 FASTA')
    mode = ctx.p.get('mode') or 'plain'
    if mode not in ('conserved', 'plain'):
        abort(400, '无效的引物设计模式')
    num_return = max(1, min(int(ctx.p.get('num_return') or 3), 20))
    specificity = bool(ctx.p.get('specificity'))

    def job(log, prog, cancel):
        import json as _json
        logger = TaskLogger(callback=log)
        if mode == 'conserved':
            p_dir = check_path(os.path.join(ctx.run_dir, '05_phylo'),
                               must_exist=False, in_platform=True)
            os.makedirs(os.path.join(p_dir, 'G1'), exist_ok=True)
            shutil.copyfile(fasta, os.path.join(p_dir, 'G1', 'aln.fasta'))
            with safe_open(os.path.join(p_dir, 'summary.json'), 'wt') as f:
                _json.dump({'groups': [{'group': 'G1', 'dir': 'G1'}]}, f)
        else:
            a_dir = check_path(os.path.join(ctx.run_dir, '03_assembly'),
                               must_exist=False, in_platform=True)
            os.makedirs(a_dir, exist_ok=True)
            shutil.copyfile(fasta, os.path.join(a_dir, 'viral_contigs.fasta'))
        from vp.primer import design_primers
        prog('primer', 0.1, f'primer3 引物设计（{mode}）')
        res = design_primers(ctx.run_dir, mode=mode, num_return=num_return,
                             logger=logger, force=True,
                             do_specificity=specificity)
        prog('done', 1.0, f"引物 {res.get('n_primers', 0)} 对")
        logger.close()
        return res
    return job


def _tool_job_identify(ctx):
    """② 病毒鉴定与提取（fastq 双端/单端 或 fasta contigs）。"""
    inp = ctx.req('input', '输入文件')
    itype = ctx.p.get('input_type') or 'pe'
    if itype not in ('pe', 'single', 'fasta'):
        abort(400, '无效的输入类型')
    conf = float(ctx.p.get('confidence') or 0)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.kunpeng import classify, parse_classify_output
        inputs = [inp]
        if itype == 'pe':
            inputs.append(ctx.req('input2', 'R2 FASTQ'))
        prog('classify', 0.1, 'kunpeng 病毒库分类中')
        res = classify(ctx.db_virus, inputs, os.path.join(ctx.run_dir, 'classify'),
                       paired=(itype == 'pe'), threads=ctx.threads,
                       confidence=conf, logger=logger,
                       progress=lambda pp, mm: prog(
                           'classify', 0.1 + pp * 0.75, mm))
        ids = {}
        if res['kraken']:
            for flag, rid, taxid, _l, mapping in parse_classify_output(res['kraken']):
                if flag != 'C':
                    continue
                # kraken2 口径分值：支持该 taxid 的片段数 / 映射列总片段数
                sup = tot = 0
                for seg in (mapping or '').split():
                    t, _, c = seg.rpartition(':')
                    try:
                        cnt = int(c)
                    except ValueError:
                        continue
                    tot += cnt
                    if t == str(taxid):
                        sup += cnt
                ids[rid] = (taxid, round(sup / tot, 3) if tot else 0)
        prog('extract', 0.9, '提取病毒候选序列')
        n_ext = {}
        if ids:
            for i, s_ in enumerate(inputs, 1):
                tag = '' if len(inputs) == 1 else f'_{i}'
                dst = os.path.join(ctx.run_dir, f'viral_sequences{tag}.fasta')
                n_ext[os.path.basename(dst)] = _extract_records(s_, set(ids), dst)
                logger.log(f'提取病毒序列 {os.path.basename(dst)}: '
                           f'{n_ext[os.path.basename(dst)]} 条')
            with safe_open(os.path.join(ctx.run_dir, 'viral_ids.tsv'), 'wt') as f:
                f.write('seq_id\ttaxid\tscore\n')
                for rid, (tx, sc) in ids.items():
                    f.write(f'{rid}\t{tx}\t{sc}\n')
        prog('extract', 1.0, '完成')
        logger.close()
        # n_extracted 必须是标量：结果预览的统计条只展示标量字段
        return {'n_classified': len(ids), 'n_extracted': sum(n_ext.values()),
                'kreport': res['kreport']}
    return job


def _tool_job_assemble(ctx):
    """③ 病毒组装（SPAdes，可选模式）。"""
    r1 = ctx.req('r1', 'R1 FASTQ')
    r2 = ctx.req('r2', 'R2 FASTQ')
    mode = ctx.p.get('mode') or 'metaviral'
    mem = int(ctx.p.get('memory') or 64)
    min_len = int(ctx.p.get('min_len') or 200)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.assembly import run_spades, filter_contigs
        out = os.path.join(ctx.run_dir, 'assembly')
        prog('spades', 0.05, 'SPAdes 组装中（耗时主要步骤）')
        spades_out = run_spades(r1, r2, out, mode=mode, threads=ctx.threads,
                                memory_gb=mem, logger=logger,
                                progress=lambda pp, mm: prog(
                                    'spades', 0.05 + pp * 0.85, mm))
        # 组装产物可能是 contigs.fasta（常规/metaviral-meta）或
        # transcripts.fasta（低覆盖自动降级 rna 模式），以 run_spades 实际
        # 返回的路径为准，不能硬编码 contigs.fasta（否则降级时找不到文件）。
        prog('filter', 0.92, 'contig 长度过滤')
        filtered, n_c, total_bp = filter_contigs(
            spades_out,
            os.path.join(ctx.run_dir, 'contigs.filtered.fasta'),
            min_len=min_len, logger=logger)
        prog('filter', 1.0, '完成')
        logger.close()
        return {'n_contigs': n_c, 'total_bp': total_bp,
                'contigs': filtered}
    return job


def _tool_job_contigs(ctx):
    """④ contig 病毒分类与提取（输入 contigs fasta）。"""
    contigs = ctx.req('contigs', 'contigs FASTA')
    min_len = int(ctx.p.get('min_len') or 200)
    conf = float(ctx.p.get('confidence') or 0)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.assembly import filter_contigs
        from vp.kunpeng import classify, parse_classify_output
        from vp.utils import run_cmd
        from vp.contig_annot import classify_rows, genus_avg_map, RANKS

        prog('genus_lens', 0.03, '统计属平均基因组长度（首跑需建缓存）')
        genus_map = genus_avg_map(logger=logger)

        prog('filter', 0.05, 'contig 长度过滤')
        filtered, n_c, _bp = filter_contigs(
            contigs, os.path.join(ctx.run_dir, 'contigs.filtered.fasta'),
            min_len=min_len, logger=logger)
        if n_c == 0:
            raise RuntimeError('过滤后无 contigs（检查最小长度设置）')

        prog('classify', 0.15, 'kunpeng 病毒库分类')
        res = classify(ctx.db_virus, [filtered],
                       os.path.join(ctx.run_dir, 'classify'),
                       paired=False, threads=ctx.threads, confidence=conf,
                       logger=logger,
                       progress=lambda pp, mm: prog(
                           'classify', 0.15 + pp * 0.45, mm))
        ids = {}
        if res['kraken']:
            for flag, _rid, _tx, _l, _pa in parse_classify_output(res['kraken']):
                if flag == 'C':
                    ids[_rid] = _tx

        # metabuli 风格分类表：8 级谱系 + 属平均长度 + 近完整判定
        prog('annot', 0.8, '谱系注释与属长比整理')
        rows = []
        if res['kraken']:
            rows = classify_rows(res['kraken'], genus_map)
        tsv_path = os.path.join(ctx.run_dir, 'virus_classification.tsv')
        header = (['contig', 'taxid', 'taxon'] + RANKS
                  + ['length', 'genus_avg_len', 'ratio', 'near_complete',
                     'score', 'kmer_support', 'kmer_total'])
        with safe_open(tsv_path, 'wt') as f:
            f.write('\t'.join(header) + '\n')
            for r in rows:
                f.write('\t'.join(str(r.get(k, '')) for k in header) + '\n')

        prog('extract', 0.88, '提取病毒 contigs（带谱系 header）')
        n_viral = 0
        viral_fa = None
        if ids:
            viral_fa = os.path.join(ctx.run_dir, 'viral_contigs.fasta')
            seqs = _load_fasta_seqs(filtered)
            with safe_open(viral_fa, 'wt') as f:
                for rid, taxid in ids.items():
                    seq = seqs.get(rid)
                    if seq is None:
                        continue
                    n_viral += 1
                    row = next((r for r in rows if r['contig'] == rid), None)
                    lineage = ';'.join(row[r] for r in RANKS
                                       if row and row.get(r))
                    taxon = (row or {}).get('taxon', '')
                    f.write(f'>{rid} taxid={taxid} taxon='
                            f'{taxon.replace(" ", "_")} '
                            f'lineage={lineage}\n')
                    for i in range(0, len(seq), 70):
                        f.write(seq[i:i + 70] + '\n')

        prog('extract', 1.0, '完成')
        logger.close()
        return {'n_contigs': n_c, 'n_viral': n_viral,
                'kreport': res['kreport'],
                'classification': tsv_path, 'viral_fasta': viral_fa}
    return job


def _tool_job_structcmp(ctx):
    """结构比较：多条序列 MAFFT 全长比对 → 两两 identity 矩阵（SDT 口径）。"""
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    max_n = max(3, min(int(ctx.p.get('max_n') or 30), 200))

    def _short(h, idx):
        name = re.split(r'[\s|]', (h or '').strip())[0][:40]
        name = re.sub(r'[^A-Za-z0-9_\-.]', '_', name)
        return name or f'seq{idx}'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.phylo import _run_mafft
        from vp.utils import iter_fasta

        prog('read', 0.05, '读取与筛选序列')
        recs, seen = [], {}
        for h, s in iter_fasta(seqs_fa):
            name = _short(h, len(recs) + 1)
            if name in seen:
                seen[name] += 1
                name = f'{name}_{seen[name]}'
            else:
                seen[name] = 0
            recs.append((name, s.upper()))
        if len(recs) < 2:
            raise RuntimeError('FASTA 中少于 2 条序列，无法做结构比较')
        recs = recs[:max_n]

        capped_fa = os.path.join(ctx.run_dir, 'input.fasta')
        with safe_open(capped_fa, 'wt') as f:
            for name, s in recs:
                f.write(f'>{name}\n')
                for i in range(0, len(s), 70):
                    f.write(s[i:i + 70] + '\n')

        prog('aln', 0.15, 'MAFFT 全长比对')
        aln = _run_mafft(capped_fa, os.path.join(ctx.run_dir, 'aln.fasta'),
                         threads=ctx.threads, logger=logger)

        prog('matrix', 0.75, '计算两两 identity 矩阵')
        names, cols = [], []
        for h, s in iter_fasta(aln):
            names.append(_short(h, len(names) + 1))
            cols.append(s.upper())
        n = len(names)
        matrix = [[100.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                same = comp = 0
                for a, b in zip(cols[i], cols[j]):
                    if a == '-' and b == '-':
                        continue
                    comp += 1
                    if a == b:
                        same += 1
                pid = round(same / comp * 100, 2) if comp else 0.0
                matrix[i][j] = matrix[j][i] = pid

        tsv = os.path.join(ctx.run_dir, 'identity_matrix.tsv')
        with safe_open(tsv, 'wt') as f:
            f.write('\t'.join([''] + names) + '\n')
            for i in range(n):
                f.write('\t'.join([names[i]] +
                                  [f'{matrix[i][j]:.2f}' for j in range(n)]) + '\n')
        data = {'names': names, 'matrix': matrix, 'n': n,
                'aln_cols': len(cols[0]) if cols else 0}
        js = os.path.join(ctx.run_dir, 'identity_matrix.json')
        with safe_open(js, 'wt') as f:
            json.dump(data, f, ensure_ascii=False)

        prog('done', 1.0, '完成')
        logger.close()
        return {'n_seqs': n, 'aln': aln, 'matrix_tsv': tsv,
                'matrix_json': js, 'aln_cols': data['aln_cols']}
    return job


def _tool_job_verify(ctx):
    """候选序列验证（对齐 02b/09b）：宿主筛选 → 长度分流 → 双路过滤。

    输入可选已有 contigs 运行（run）或独立 FASTA（fasta）。
    参数：host（默认 all）、methods（blastx/cdd 组合）、combine（union/intersection）。
    """
    run_ref = (ctx.p.get('run') or '').strip()
    fasta = ctx.opt('fasta') if ctx.p.get('fasta') else None
    host = (ctx.p.get('host') or 'all').strip() or 'all'
    combine = (ctx.p.get('combine') or 'union').strip()
    methods = ctx.p.get('methods') or ['blastx', 'cdd']
    if isinstance(methods, str):
        methods = [m.strip() for m in methods.split(',') if m.strip()]
    methods = [m for m in methods if m in ('blastx', 'cdd')]

    if run_ref and not re.fullmatch(r'[A-Za-z0-9_\-]+', run_ref):
        abort(400, '无效的运行名')
    if combine not in ('union', 'intersection'):
        combine = 'union'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.verify import verify

        # 输入解析：优先 run（其 viral_contigs.fasta），否则独立 fasta
        run_dir = ctx.run_dir
        if run_ref:
            src_run = check_path(os.path.join(_tool_runs_root(), run_ref),
                                 must_exist=True, in_platform=True)
            vfa = os.path.join(src_run, 'viral_contigs.fasta')
            if not os.path.isfile(vfa):
                raise RuntimeError('该运行无 viral_contigs.fasta（先跑 contig 分类）')
            # 宿主归属依赖源运行的分类表
            for _fn in ('virus_classification.tsv',):
                s = os.path.join(src_run, _fn)
                if os.path.isfile(s):
                    import shutil
                    shutil.copyfile(s, os.path.join(run_dir, _fn))
        else:
            if not fasta:
                abort(400, '请选择 contigs 运行或输入 FASTA')
            vfa = fasta

        summary = verify(run_dir, vfa, host=host, methods=methods,
                         combine=combine, threads=ctx.threads, logger=logger,
                         progress=lambda st, fr, msg: prog(st, fr, msg))
        logger.close()
        return summary
    return job


def _tool_job_consensus(ctx):
    """共识序列与变异分析（验证模块之后）：reads 回贴 → 共识序列 + 变异谱。

    参考三条来源：
      1. kvsuite 运行选参考（virus-fasta/ref_<acc>/）—— 已知病毒基因组
      2. contig 分类运行的 viral_contigs.fasta —— 未知/组装候选
      3. 独立 FASTA
    reads 默认取**去宿主后全量**（01_host_removal）——变异检测无偏；
    02_virus_screen 已按相似度丢过一轮 reads，会系统性低估变异与准种多样性。
    默认重新比对（而非复用 kvsuite 的 BAM）：旧 BAM 有比对偏好性，
    新 BAM 覆盖更全。
    """
    run_ref = (ctx.p.get('run') or '').strip()
    fasta = ctx.opt('fasta') if ctx.p.get('fasta') else None
    kv_run = (ctx.p.get('kv_run') or '').strip()
    kv_refs = [x.strip() for x in (ctx.p.get('kv_refs') or '').split(',') if x.strip()]
    reuse_bam = (ctx.p.get('reuse_bam') or '').strip()
    reads_src = (ctx.p.get('reads_src') or 'host_removed').strip()
    ambig = (ctx.p.get('ambig') or 'N').strip()[:1] or 'N'

    def _num(v, default, cast):
        try:
            return cast(v)
        except (TypeError, ValueError):
            return default

    min_qual = _num(ctx.p.get('min_qual'), 20, int)
    min_depth = _num(ctx.p.get('min_depth'), 10, int)
    min_freq = _num(ctx.p.get('min_freq'), 0.5, float)
    min_mapq = _num(ctx.p.get('min_mapq'), 10, int)
    min_minor_freq = _num(ctx.p.get('min_minor_freq'), 0.02, float)
    min_cov_pct = _num(ctx.p.get('min_cov_pct'), 10.0, float)
    min_qual = max(0, min(min_qual, 60))
    min_depth = max(1, min(min_depth, 100000))
    min_freq = min(1.0, max(0.0, min_freq))

    if run_ref and not re.fullmatch(r'[A-Za-z0-9_\-]+', run_ref):
        abort(400, '无效的运行名')
    if kv_run and not re.fullmatch(r'[A-Za-z0-9_\-]+', kv_run):
        abort(400, '无效的 kvsuite 运行名')
    if kv_refs and not kv_run:
        abort(400, '选了参考序列但未指定 kvsuite 运行')
    for _r in kv_refs:
        if not re.fullmatch(r'[A-Za-z0-9_.\-]+', _r):
            abort(400, '无效的 accession: %s' % _r)
    if reads_src not in CONSENSUS_READS_LAYERS:
        reads_src = 'host_removed'
    if reuse_bam:
        # 只允许平台运行目录内的 BAM（防任意路径写入）
        reuse_bam = check_path(reuse_bam, must_exist=False, in_platform=True)
        if not os.path.isfile(reuse_bam):
            abort(400, '复用的 BAM 不存在：%s' % reuse_bam)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.consensus import consensus_and_variants, find_read_pairs

        # ---- 参考 ----
        # 三条来源，优先级：kvsuite 选参考 > contig 运行 > 独立 FASTA
        #   kvsuite：<run>/kvsuite/virus-fasta/ref_<acc>/ref_<acc>.ref.fasta
        #            多选按 accession 合并成临时库（序列名已保证唯一）
        src_run = None
        vfa = None
        if kv_refs:
            kv_base = check_path(os.path.join(_tool_runs_root(), kv_run,
                                              'kvsuite'),
                                 must_exist=True, in_platform=True)
            fa_dir = os.path.join(kv_base, 'virus-fasta')
            parts, missing = [], []
            for acc in kv_refs:
                sub = os.path.join(fa_dir, 'ref_%s' % acc)
                cand = os.path.join(sub, 'ref_%s.ref.fasta' % acc)
                if os.path.isfile(cand) and os.path.getsize(cand) > 0:
                    parts.append((acc, cand))
                    continue
                # 容错：目录在但文件名不同
                found = None
                if os.path.isdir(sub):
                    for fn in sorted(os.listdir(sub)):
                        if fn.endswith(('.fasta', '.fa', '.fna')):
                            found = os.path.join(sub, fn)
                            break
                if found:
                    parts.append((acc, found))
                else:
                    missing.append(acc)
            if not parts:
                raise RuntimeError(
                    '所选参考在 %s 下均无 FASTA（先跑共识段生成 virus-fasta）'
                    % fa_dir)
            if missing:
                log('跳过无 FASTA 的参考: %s' % ', '.join(missing))
            # 写到运行目录，作为本次分析的临时参考
            tmp_ref = os.path.join(ctx.run_dir, 'ref_%s.fasta' % kv_run)
            n_seq = 0
            with open(tmp_ref, 'w', encoding='utf-8', newline='\n') as out:
                for acc, path in parts:
                    with open(path, encoding='utf-8', errors='replace') as fh:
                        for line in fh:
                            if line.startswith('>'):
                                n_seq += 1
                                head = line[1:].strip().split()[0]
                                # accession 唯一化：序列名已一致时保留原样
                                if head == acc:
                                    out.write(line if line.endswith('\n')
                                              else line + '\n')
                                else:
                                    out.write('>%s\n' % acc)
                            else:
                                out.write(line if line.endswith('\n')
                                          else line + '\n')
            log('kvsuite 参考库：%d 条（%d 个 accession）-> %s'
                % (n_seq, len(parts), os.path.basename(tmp_ref)))
            vfa = tmp_ref
        elif run_ref:
            src_run = check_path(os.path.join(_tool_runs_root(), run_ref),
                                 must_exist=True, in_platform=True)
            vfa = os.path.join(src_run, 'viral_contigs.fasta')
            if not os.path.isfile(vfa):
                raise RuntimeError('该运行无 viral_contigs.fasta（先跑 contig 分类）')
        else:
            if not fasta:
                abort(400, '请选择 contigs 运行、kvsuite 参考或输入参考 FASTA')
            vfa = fasta

        # ---- reads ----
        # 默认重新比对：用去宿主后全量 reads，让新 BAM 覆盖更全。
        # kvsuite 自己的 bam/virus_reads 已按病毒筛选过，复用会引入比对
        # 偏好性并系统性低估变异，因此不作为默认来源。
        reads = []
        if not reuse_bam:
            for p in (ctx.p.get('reads') or '').split(','):
                p = p.strip()
                if not p:
                    continue
                if os.path.isdir(p):
                    reads.extend(find_read_pairs(p))
                else:
                    reads.append(p)
            if not reads:
                if not src_run:
                    # kvsuite 参考或独立 FASTA：也需要 reads
                    if kv_run:
                        abort(400, '选择了 kvsuite 参考，请显式指定回贴 reads 文件/目录')
                    abort(400, '未指定运行时必须直接给出 reads 文件路径')
                layer = os.path.join(src_run, CONSENSUS_READS_LAYERS[reads_src])
                reads = find_read_pairs(layer)
                if not reads:
                    raise RuntimeError(
                        '在 %s 下未找到 reads（期望 *_R1/*.fastq.gz 配对文件）'
                        % CONSENSUS_READS_LAYERS[reads_src])
        else:
            log('映射模式：复用 BAM %s' % os.path.basename(reuse_bam))

        summary = consensus_and_variants(
            ctx.run_dir, vfa, reads, out_subdir='consensus',
            min_qual=min_qual, min_depth=min_depth, min_freq=min_freq,
            ambig=ambig, min_mapq=min_mapq, min_minor_freq=min_minor_freq,
            min_cov_pct=min_cov_pct, preset='sr', threads=ctx.threads,
            logger=logger, progress=lambda st, fr, msg: prog(st, fr, msg),
            reuse_bam=reuse_bam or None)
        logger.close()
        return summary
    return job


def _tool_job_kvsuite(ctx):
    """已知病毒识别与定量（known_virus_suite 五段整合）。

    完全照搬 D:/桌面/延伸基因组/MMPV-RNA/virome_analysis_pipeline 的做法：
      鉴定 → 过滤 → 共识 → 深度绘图 → 变异注释
    引擎 minibwa 替代 bowtie2，**其余一点不改**。

    caller 用 bcftools mpileup+call（freebayes/lofreq/ivar 在 Windows 上均不可得），
    参数与阈值为实测定稿，详见 known_virus_suite/POSCOUNTS_REMOVAL_PLAN.md §4.2。
    """
    import json
    import subprocess

    # ── 输入：勾选样品 → 临时 sample-sheet TSV(name,r1,r2) ──
    # 变异段（variant）不吃 reads，只需 BAM 或 VCF，因此允许 samples 为空。
    stage = (ctx.p.get('stage') or 'all').strip().lower()
    if stage not in ('all', 'index', 'identify', 'filter', 'consensus',
                     'plot', 'variant'):
        stage = 'all'
    samples_raw = (ctx.p.get('samples') or '').strip()
    sheet_in = (ctx.p.get('sample_sheet') or '').strip()
    if sheet_in:
        sheet = ctx.req('sample_sheet', '样本表 TSV')
    elif stage == 'variant':
        sheet = None          # 变异段不用样本表
    else:
        if not samples_raw:
            abort(400, '请选择样品，或提供样本表 TSV')
        names = [s.strip() for s in samples_raw.split(',') if s.strip()]
        if not names:
            abort(400, '未选择有效样品')
        from vp.pipeline import load_sample_input
        rows = []
        for nm in names:
            sd = check_path(os.path.join(DIRS['results'], nm),
                            must_exist=True, in_platform=True)
            r1, r2, _proj = load_sample_input(sd)
            if not r1 or not os.path.isfile(r1):
                abort(400, f'样品 {nm} 的 R1 不可用: {r1}')
            rows.append((nm, r1, r2 or ''))
        sheet = os.path.join(ctx.run_dir, 'sample_sheet.tsv')
        with open(sheet, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write('name\tr1\tr2\n')
            for nm, r1, r2 in rows:
                fh.write(f'{nm}\t{r1}\t{r2}\n')

    # ── 参考序列与注释 ──
    ref = (ctx.p.get('reference') or '').strip()
    if ref:
        ref = ctx.req('reference', '参考序列')
    else:
        ref = os.path.join(DIRS['virus_src'], 'final.cluster.ref.fasta')
        ref = check_path(ref, must_exist=True, in_platform=True)
    ref_info = (ctx.p.get('ref_info') or '').strip()
    if ref_info:
        ref_info = ctx.req('ref_info', '参考注释')
    else:
        cand = os.path.join(DIRS['virus_src'], 'final.cluster.ref_info.tsv')
        ref_info = cand if os.path.isfile(cand) else None

    engine = (ctx.p.get('engine') or 'minibwa').strip().lower()
    if engine not in ('minibwa', 'salmon'):
        engine = 'minibwa'
    # 索引复用：① 用户在前端选的「鉴定库」优先（build 页建好的）；
    # ② 默认参考 + 预建索引在位时自动复用（避免每个 run 重建 60 MB 索引）；
    # ③ 都不满足则回退 <out>/index 临时建。
    index_dir = None
    user_idx = (ctx.p.get('index_dir') or '').strip()
    if user_idx:
        index_dir = ctx.req('index_dir', '鉴定库目录')
    elif engine == 'minibwa':
        kv_idx = os.path.join(DIRS['virus_src'], 'kv_index')
        default_ref = os.path.join(DIRS['virus_src'], 'final.cluster.ref.fasta')
        if (os.path.isfile(os.path.join(kv_idx, 'minibwa.mbw'))
                and os.path.abspath(ref) == os.path.abspath(default_ref)):
            index_dir = kv_idx
    elif engine == 'salmon':
        kv_idx = os.path.join(DIRS['virus_src'], 'kv_index')
        default_ref = os.path.join(DIRS['virus_src'], 'final.cluster.ref.fasta')
        if (os.path.isfile(os.path.join(kv_idx, 'salmon_k31', 'info.json'))
                and os.path.abspath(ref) == os.path.abspath(default_ref)):
            index_dir = kv_idx

    def _num(v, default, cast, lo=None, hi=None):
        try:
            x = cast(v)
        except (TypeError, ValueError):
            return default
        if lo is not None:
            x = max(lo, x)
        if hi is not None:
            x = min(hi, x)
        return x

    min_cov = _num(ctx.p.get('min_cov'), 10.0, float, 0.0, 100.0)
    min_depth = _num(ctx.p.get('min_depth'), 0.5, float, 0.0)
    min_reads = _num(ctx.p.get('min_reads'), 10, int, 0)
    min_poisson = _num(ctx.p.get('min_poisson'), 0.3, float, 0.0, 1.0)
    variant_qual = _num(ctx.p.get('variant_qual'), 3.5, float, 0.0)
    min_freq = _num(ctx.p.get('min_freq'), 0.05, float, 0.0, 1.0)
    aa_label_cutoff = _num(ctx.p.get('aa_label_cutoff'), 0.50, float, 0.0, 1.01)
    max_aa_labels = _num(ctx.p.get('max_aa_labels'), 40, int, 0)

    def job(log, prog, cancel):
        suite = os.path.join(PLATFORM_ROOT, 'known_virus_suite',
                             'known_virus_suite.py')
        if not os.path.isfile(suite):
            raise RuntimeError(f'未找到 known_virus_suite.py: {suite}')

        out_dir = os.path.join(ctx.run_dir, 'kvsuite')
        os.makedirs(out_dir, exist_ok=True)
        cmd = [sys.executable, suite, stage,
               '--out', out_dir,
               '--reference', ref,
               '--engine', engine,
               '--min-cov', str(min_cov),
               '--min-depth', str(min_depth),
               '--min-reads', str(min_reads),
               '--min-poisson', str(min_poisson),
               '--variant-qual', str(variant_qual),
               '--min-freq', str(min_freq),
               '--aa-label-cutoff', str(aa_label_cutoff),
               '--max-aa-labels', str(max_aa_labels)]
        if sheet:
            cmd += ['--sample-sheet', sheet]
        # 变异段输入：外部 VCF 优先于 BAM 目录（CLI 内也是这个优先级）
        input_vcf = (ctx.p.get('input_vcf') or '').strip()
        if input_vcf:
            cmd += ['--input-vcf', ctx.req('input_vcf', '变异 VCF 文件')]
        bam_dir = (ctx.p.get('bam_dir') or '').strip()
        if bam_dir:
            cmd += ['--bam-dir', ctx.req('bam_dir', 'BAM 目录')]
        if ref_info:
            cmd += ['--ref-info', ref_info]
        if index_dir:
            cmd += ['--index-dir', index_dir]
        if ctx.threads:
            cmd += ['--threads', str(ctx.threads),
                    '--align-threads', str(ctx.threads)]
        if ctx.p.get('ncbi_email'):
            cmd += ['--ncbi-email', str(ctx.p['ncbi_email']).strip()]
        if ctx.p.get('ncbi_api_key'):
            cmd += ['--ncbi-api-key', str(ctx.p['ncbi_api_key']).strip()]
        if ctx.p.get('no_genes'):
            cmd.append('--no-genes')
        if ctx.p.get('no_variant_evo'):
            cmd.append('--no-variant-evo')
        if ctx.p.get('all_variants'):
            cmd.append('--all-variants')
        if ctx.p.get('no_snpgenie'):
            cmd.append('--no-snpgenie')

        log('$ ' + ' '.join(cmd))
        env = dict(os.environ)
        env.setdefault('PYTHONIOENCODING', 'utf-8')
        proc = subprocess.Popen(cmd, cwd=PLATFORM_ROOT,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                text=True, encoding='utf-8',
                                errors='replace', bufsize=1, env=env)
        tail = []
        for line in proc.stdout:
            line = line.rstrip('\n')
            if not line:
                continue
            log(line)
            tail.append(line)
            if len(tail) > 40:
                tail.pop(0)
            # cancel 是 threading.Event（见 TaskServer._run），必须用 is_set()，
            # 直接 cancel() 会报 'Event' object is not callable。
            if cancel is not None and cancel.is_set():
                proc.kill()
                raise RuntimeError('用户取消')
        rc = proc.wait()
        if rc != 0:
            raise RuntimeError(f'known_virus_suite 退出码 {rc}；末几行：\n'
                               + '\n'.join(tail[-8:]))

        # ── 汇总产物（供结果面板展示）──
        result = {'stage': stage, 'engine': engine, 'out': out_dir,
                  'reference': ref}
        fsum = os.path.join(out_dir, 'filter', 'filtered.tsv')
        if os.path.isfile(fsum):
            result['filtered_tsv'] = fsum
            try:
                with open(fsum, encoding='utf-8', errors='replace') as fh:
                    result['n_confirmed'] = max(0, sum(1 for _ in fh) - 1)
            except OSError:
                pass
        vsum = os.path.join(out_dir, 'variant_summary.json')
        if os.path.isfile(vsum):
            result['variant_summary'] = vsum
            try:
                with open(vsum, encoding='utf-8') as fh:
                    data = json.load(fh)
                rows = data.get('results') if isinstance(data, dict) else data
                rows = rows or []
                result['n_variant_genomes'] = len(rows)
                result['n_variants'] = sum(int(d.get('n_variants') or 0)
                                           for d in rows)
                if isinstance(data, dict):
                    result['n_variant_skipped'] = int(data.get('n_skipped') or 0)
            except (OSError, ValueError, TypeError):
                pass
        return result
    return job


def _tool_job_quicktree(ctx):
    """快速建树：多条序列 MAFFT 全长比对 → NJ（纯 Python）/ FastTree。"""
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    method = (ctx.p.get('method') or 'nj').strip().lower()
    if method not in ('nj', 'fasttree'):
        abort(400, '建树方法仅支持 nj / fasttree')
    max_n = max(3, min(int(ctx.p.get('max_n') or 100), 500))

    def _short(h, idx):
        name = re.split(r'[\s|]', (h or '').strip())[0][:40]
        name = re.sub(r'[^A-Za-z0-9_\-.]', '_', name)
        return name or f'seq{idx}'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.phylo import _run_mafft, _run_nj, _run_fasttree
        from vp.utils import iter_fasta

        prog('read', 0.05, '读取与筛选序列')
        recs, seen = [], {}
        for h, s in iter_fasta(seqs_fa):
            name = _short(h, len(recs) + 1)
            if name in seen:
                seen[name] += 1
                name = f'{name}_{seen[name]}'
            else:
                seen[name] = 0
            recs.append((name, s.upper()))
        if len(recs) < 2:
            raise RuntimeError('FASTA 中少于 2 条序列，无法建树')
        recs = recs[:max_n]

        capped_fa = os.path.join(ctx.run_dir, 'input.fasta')
        with safe_open(capped_fa, 'wt') as f:
            for name, s in recs:
                f.write(f'>{name}\n')
                for i in range(0, len(s), 70):
                    f.write(s[i:i + 70] + '\n')

        prog('aln', 0.2, 'MAFFT 全长比对')
        aln = _run_mafft(capped_fa, os.path.join(ctx.run_dir, 'aln.fasta'),
                         threads=ctx.threads, logger=logger)

        prog('tree', 0.75, 'NJ 建树' if method == 'nj' else 'FastTree 建树')
        tree_name = 'nj.nwk' if method == 'nj' else 'tree.nwk'
        if method == 'nj':
            _run_nj(aln, os.path.join(ctx.run_dir, tree_name), logger=logger)
        else:
            _run_fasttree(aln, os.path.join(ctx.run_dir, tree_name),
                          logger=logger)

        prog('done', 1.0, '完成')
        logger.close()
        return {'n_seqs': len(recs), 'method': method, 'tree': tree_name}
    return job


def _tool_job_align(ctx):
    """序列比对（独立模块）：MAFFT（auto/L-INS-i/fast）→ trimAl 清剪。

    输入 FASTA（核酸或蛋白，自动判别；≥2 条）。产物：input.fasta /
    aln.fasta / aln.trim.fasta（trimAl 关闭时无）/ summary.json，
    可在「比对查看器」彩色浏览与编辑。
    """
    from vp.utils import iter_fasta, write_fasta_record
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    strategy = (ctx.p.get('strategy') or 'auto').strip().lower()
    if strategy not in ('auto', 'linsi', 'fast'):
        abort(400, '比对策略仅支持 auto / linsi / fast')
    trimal = (ctx.p.get('trimal') or 'automated1').strip().lower()
    if trimal not in ('automated1', 'gappyout', 'strict', 'none'):
        abort(400, 'trimAl 方法仅支持 automated1 / gappyout / strict / none')
    max_n = max(3, min(int(ctx.p.get('max_n') or 200), 500))

    def job(log, prog, cancel):
        import shutil
        logger = TaskLogger(callback=log)
        from vp.phylo import _run_mafft, _run_trimal
        from vp.sdt_exact import detect_seqtype
        prog('read', 0.05, '读取与筛选序列')
        recs, seen = [], {}
        for h, s in iter_fasta(seqs_fa):
            name = re.split(r'[\s|]', (h or '').strip())[0][:60] or \
                f'seq{len(recs) + 1}'
            name = re.sub(r'[^A-Za-z0-9_\-.]', '_', name)
            if name in seen:
                seen[name] += 1
                name = f'{name}_{seen[name]}'
            else:
                seen[name] = 0
            recs.append((name, s.upper()))
        if len(recs) < 2:
            raise RuntimeError('FASTA 中少于 2 条序列，无法比对')
        recs = recs[:max_n]
        seqtype = detect_seqtype([s for _, s in recs])
        capped = os.path.join(ctx.run_dir, 'input.fasta')
        with safe_open(capped, 'wt') as f:
            for name, s in recs:
                write_fasta_record(f, name, s)
        logger.log(f'序列类型判定: {"蛋白(aa)" if seqtype == "aa" else "核酸(nt)"}'
                   f'，{len(recs)} 条参与比对')
        prog('align', 0.15, f'MAFFT 比对（{strategy}）')
        aln = _run_mafft(capped, os.path.join(ctx.run_dir, 'aln.fasta'),
                         threads=ctx.threads, logger=logger, strategy=strategy)
        trim_info = {'applied': False}
        aln_used = aln
        if trimal != 'none':
            prog('trim', 0.7, f'trimAl 清剪（{trimal}）')
            aln_used, trim_info = _run_trimal(
                aln, os.path.join(ctx.run_dir, 'aln.trim.fasta'),
                logger=logger)
        prog('done', 1.0, '完成')
        logger.close()
        return {'n_seqs': len(recs), 'seqtype': seqtype, 'strategy': strategy,
                'trimal': trimal, 'trim': trim_info,
                'aln': 'aln.fasta',
                'aln_used': os.path.basename(aln_used)}
    return job


def _tool_job_sdt(ctx):
    """SDT 精确分析（纯 Python 复刻 SDTv1.3，替代外部 SDT exe）：
    逐对 MAFFT 独立比对 → Get_Similarity 公式 → 聚类排序热图 + 分布图。
    aligned=True 时输入为已比对 MSA，跳过比对器直接按公式计算；
    seqtype='auto' 自动判别核酸/蛋白（MMPV 双轨口径），蛋白按 AA 同一性。"""
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    max_n = max(3, min(int(ctx.p.get('max_n') or 30), 200))
    orient = bool(ctx.p.get('orient', True))
    seqtype = (ctx.p.get('seqtype') or 'auto').strip().lower()
    if seqtype not in ('auto', 'nt', 'aa'):
        seqtype = 'auto'
    palette = (ctx.p.get('palette') or 'sdt').strip().lower()
    if palette not in ('sdt', 'cividis', 'viridis', 'RdYlBu', 'Spectral',
                       'YlGnBu', 'coolwarm', 'magma'):
        palette = 'sdt'
    aligned = bool(ctx.p.get('aligned'))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        mafft = cfg.tool('mafft')
        from vp.sdt_exact import run_sdt_exact
        prog('sdt', 0.02, 'SDT 精确分析（MAFFT 逐对独立比对，SDT v1.3 口径）')
        res = run_sdt_exact(seqs_fa, ctx.run_dir, mafft, max_n=max_n,
                            orient=orient, threads=ctx.threads,
                            seqtype=seqtype,
                            palette=palette, logger=logger,
                            progress=lambda p, m: prog('sdt', p, m),
                            cancel=cancel, aligned=aligned)
        logger.close()
        return res
    return job


def _tool_job_identity(ctx):
    """核苷酸+氨基酸同一性表（BioAider Sequence Identity Matrix 口径）：
    NT 矩阵 + AA 矩阵（可选输入或最长 ORF 翻译）+ 复合热图（NT 上 /
    AA 下）+ 逐对同一性长表。"""
    nt_fa = ctx.req('nt_seqs', '核苷酸 FASTA')
    aa_fa = ctx.opt('aa_seqs')
    max_n = max(3, min(int(ctx.p.get('max_n') or 30), 200))
    palette = (ctx.p.get('palette') or 'sdt').strip().lower()
    if palette not in ('sdt', 'cividis', 'viridis', 'RdYlBu', 'Spectral',
                       'YlGnBu', 'coolwarm', 'magma'):
        palette = 'sdt'
    aligned = bool(ctx.p.get('aligned'))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        mafft = cfg.tool('mafft')
        from vp.sdt_exact import run_identity_table
        prog('idty', 0.02, '核苷酸+氨基酸同一性表（BioAider 口径）')
        res = run_identity_table(nt_fa, ctx.run_dir, mafft, aa_fasta=aa_fa,
                                 max_n=max_n, aligned=aligned,
                                 threads=ctx.threads, palette=palette,
                                 logger=logger,
                                 progress=lambda p, m: prog('idty', p, m),
                                 cancel=cancel)
        logger.close()
        return res
    return job


def _load_fasta_seqs(path):
    """FASTA → {首 token id: 序列}（单条记录可达 MB 级，仅限小文件用）。"""
    import gzip
    op = gzip.open if str(path).lower().endswith('.gz') else open
    out, name, buf = {}, None, []
    with op(path, 'rt', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    out[name] = ''.join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        out[name] = ''.join(buf)
    return out


def _find_contig_seq(run_dir, contig):
    """在工具运行目录的各 FASTA 中查找 contig 序列。"""
    for cur, _sub, fns in os.walk(run_dir):
        if os.path.join('analysis') in cur:
            continue
        for fn in fns:
            if not fn.lower().endswith(('.fasta', '.fa', '.fna', '.fas')):
                continue
            try:
                seqs = _load_fasta_seqs(os.path.join(cur, fn))
            except (OSError, ValueError):
                continue
            if contig in seqs:
                return seqs[contig]
    return None


def _extract_records(src, ids, dst):
    """从 FASTA/FASTQ（支持 .gz）提取 header 首 token 命中 ids 的记录。

    kunpeng 报告里的 read ID 经过 seqkit fq2fa 转换，已去掉 Illumina 配对
    后缀 /1、/2，匹配时对原始 header 兼容带/不带后缀两种写法；
    FASTQ 输入统一转成 FASTA 写出（与 .fasta 扩展名一致）。
    返回提取条数；ids 为空集时写出空文件。
    """
    import gzip
    name = str(src).lower()
    is_fq = name.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq'))
    op = gzip.open if name.endswith('.gz') else open

    def match(rid):
        """命中返回规范 ID（配对后缀 /1、/2 已去掉，与 kunpeng 报告一致），未命中返回 None。"""
        if rid in ids:
            return rid
        base = rid.rsplit('/', 1)
        if len(base) == 2 and base[1] in ('1', '2') and base[0] in ids:
            return base[0]
        return None

    n = 0
    with op(src, 'rt', errors='replace') as f, safe_open(dst, 'wt') as w:
        if is_fq:
            while True:
                h = f.readline()
                if not h:
                    break
                seq = f.readline().rstrip('\r\n')
                f.readline()
                f.readline()
                tok = h[1:].split()
                mid = match(tok[0]) if tok else None
                if mid:
                    w.write(f'>{mid}\n')
                    for i in range(0, len(seq), 70):
                        w.write(seq[i:i + 70] + '\n')
                    n += 1
        else:
            keep = False
            for line in f:
                if line.startswith('>'):
                    tok = line[1:].split()
                    keep = bool(tok) and match(tok[0]) is not None
                    if keep:
                        n += 1
                if keep:
                    w.write(line)
    return n

# -*- coding: utf-8 -*-
"""
阶段③ 组装与 contig 分类：
- SPAdes 组装去宿主 reads（--metaviral 默认，可选 --meta / --rna / --isolate）
- contig 长度过滤
- kunpeng 病毒库分类 contigs
- 汇总病毒 contig 表
"""
import os
import re
import csv
import glob
import shutil

from .config import get_config, DIRS
from .utils import (check_path, safe_open, run_cmd, iter_fasta, write_fasta_record,
                    count_fasta_seqs, is_step_done, mark_step_done)
from .kunpeng import classify, parse_classify_output, parse_kreport, KreportTree

def _ascii_work_base(name='vp_blast'):
    """BLAST(LMDB) 不支持中文路径：返回纯 ASCII 持久工作目录。"""
    import tempfile
    cand = os.path.abspath(tempfile.gettempdir())
    if not (os.path.isdir(cand) and all(ord(c) < 128 for c in cand)):
        cand = r'C:\Windows\Temp'
    d = os.path.join(cand, name)
    os.makedirs(d, exist_ok=True)
    return d

def find_virus_ref_fasta():
    """定位病毒参考 FASTA（virus-db 下第一个 .fasta/.fa/.fna）。"""
    src = check_path(DIRS['virus_src'], must_exist=True)
    for pat in ('*.fasta', '*.fa', '*.fna', '*.fas'):
        hits = sorted(glob.glob(os.path.join(src, pat)))
        if hits:
            return check_path(hits[0], must_exist=True)
    raise FileNotFoundError(f"{src} 下未找到病毒参考 FASTA")

def ensure_virus_blast_db(logger=None):
    """构建/复用病毒参考 BLAST 库（ASCII 路径，绕开 LMDB 中文路径缺陷）。"""
    cfg = get_config()
    makeblastdb = cfg.tool('makeblastdb')
    ref = find_virus_ref_fasta()
    db_dir = _ascii_work_base('vp_blast')
    prefix = os.path.join(db_dir, 'virus')
    # 版本戳：记录参考 FASTA 的 mtime+大小；源文件更新后自动重建 BLAST 库
    stamp = os.path.join(db_dir, 'virus.stamp')
    cur_sig = f'{os.path.getmtime(ref):.0f}:{os.path.getsize(ref)}'
    if glob.glob(prefix + '.n??'):
        old_sig = None
        try:
            with open(stamp) as sf:
                old_sig = sf.read().strip()
        except OSError:
            pass
        if old_sig == cur_sig:
            return prefix
        if logger:
            logger.log("病毒参考 FASTA 已更新，重建 BLAST 库")
        for f in glob.glob(prefix + '.*'):
            try:
                os.remove(f)
            except OSError:
                pass
    try:
        with open(stamp, 'w') as sf:
            sf.write(cur_sig)
    except OSError:
        pass
    if logger:
        logger.log(f"构建病毒 BLAST 库: {ref} -> {prefix}")
    run_cmd([makeblastdb, '-in', ref, '-dbtype', 'nucl', '-out', prefix,
             '-title', 'virus_ref'], logger=logger)
    return prefix

def _is_ascii(p):
    return all(ord(c) < 128 for c in str(p))

def _seq_file_kind(path):
    """判断序列文件的真实类型，用于给 SPAdes 中转文件起正确扩展名。

    返回 (is_gzip, kind)，kind ∈ {'fastq', 'fasta', 'other'}。
    - is_gzip：由魔数（$\x1f\x8b）或 .gz 扩展名判断。
    - kind：优先读前几字节内容嗅探（'@'开头+质量行 → fastq；'>'开头 → fasta），
      嗅探不到时回退到扩展名。
    之所以必须做：旧版把任意输入强制复制成 in_R1.fastq.gz，当输入是 FASTA
    时 SPAdes 会用 gzip 解压 FASTA 文本 → BadGzipFile（用户实报的错误）。
    """
    is_gz = str(path).lower().endswith('.gz')
    try:
        with open(path, 'rb') as fh:
            magic = fh.read(2)
        if magic == b'\x1f\x8b':
            is_gz = True
    except OSError:
        pass

    kind = 'other'
    low = str(path).lower()
    # 先按扩展名给出候选
    if is_gz:
        base = low[:-3]
    else:
        base = low
    if base.endswith(('.fastq', '.fq')):
        kind = 'fastq'
    elif base.endswith(('.fasta', '.fa', '.fna', '.fas')):
        kind = 'fasta'

    # 若非 gz，尝试从内容嗅探修正（扩展名可能被改错）
    if not is_gz:
        try:
            with open(path, 'r', encoding='utf-8', errors='replace') as fh:
                first = ''
                for ln in fh:
                    s = ln.strip()
                    if s:
                        first = s
                        break
            if first.startswith('>'):
                kind = 'fasta'
            elif first.startswith('@'):
                kind = 'fastq'
        except OSError:
            pass
    return is_gz, kind

def _shim_name(src, role):
    """给中转输入文件构成正确的目标名，保留真实格式与压缩状态。

    role: 'R1' / 'R2'。生成的扩展名与源文件实际格式一致，避免
    FASTA 被误命名为 .fastq.gz 而触发 SPAdes 的 gzip 解析崩溃。
    """
    is_gz, kind = _seq_file_kind(src)
    if kind == 'fastq':
        ext = '.fastq.gz' if is_gz else '.fastq'
    elif kind == 'fasta':
        ext = '.fasta.gz' if is_gz else '.fasta'
    else:
        # 未知类型：保留源扩展名，仅在确为 gz 时补 .gz
        low = str(src).lower()
        if is_gz and not low.endswith('.gz'):
            ext = low[low.rfind('.'):] + '.gz' if '.' in low else '.gz'
        else:
            ext = os.path.splitext(src)[1] or '.txt'
    return f'in_{role}{ext}'

def _spades_tmp_base():
    """SPAdes Windows 版不支持非 ASCII 路径：返回纯 ASCII 临时工作根目录。

    仅允许 %TEMP%（通常为 C:\\Users\\<user>\\AppData\\Local\\Temp）；
    若 TEMP 含非 ASCII 则回退固定盘根 C:\\vp_spades_tmp。
    """
    import tempfile
    cand = os.path.abspath(tempfile.gettempdir())
    if os.path.isdir(cand) and _is_ascii(cand):
        return os.path.join(cand, 'vp_spades')
    return r'C:\vp_spades_tmp'

def _spades_monitor(spades_dir, mode, progress):
    """SPAdes 看门狗：解析 spades.log 的 Assembling k=xx 行数估算进度。

    返回 monitor_fn（传给 run_cmd）；progress 为 None 时返回 None。
    """
    if progress is None:
        return None
    k_total = {'metaviral': 4, 'meta': 5, 'rna': 3, 'isolate': 4}.get(mode, 4)
    log_p = os.path.join(spades_dir, 'spades.log')
    seen = set()

    def _watch():
        try:
            with open(log_p, 'r', encoding='utf-8', errors='replace') as f:
                tail = f.read()[-8000:]
            for m in re.finditer(r'Assembling\s+(?:k\s*=\s*)?(\d+)', tail):
                seen.add(m.group(1))
            k_done = len(seen)
            if k_done:
                progress(min(0.05 + k_done / k_total * 0.85, 0.92),
                         f'SPAdes 组装中：k 阶段 {k_done}/{k_total}')
        except OSError:
            pass
    return _watch

def run_spades(r1, r2, out_dir, mode='metaviral', threads=None, memory_gb=64,
               logger=None, progress=None):
    """运行 SPAdes。返回 contigs fasta 路径（rna 模式为 transcripts.fasta）。

    metaviral/meta 模式未拼出 contigs 时（低深度样品常见——metaviralSPAdes
    的"染色体外元件"筛选门槛较严），自动清空输出并用 rna 模式重试一次。
    平台路径含中文时自动经纯 ASCII 临时目录中转（输入复制过去、结果拷回）。
    """
    import shutil
    import uuid
    try:
        return _run_spades_once(r1, r2, out_dir, mode, threads, memory_gb,
                                logger, progress)
    except RuntimeError as e:
        if mode in ('metaviral', 'meta') and '未找到 contigs' in str(e):
            if logger:
                logger.log(f"{mode} 模式未拼出 contigs（该模式对低深度病毒样品"
                           f"筛选较严），自动改用 rna 模式重试组装 ...")
            spades_dir = check_path(os.path.join(out_dir, 'spades'),
                                    must_exist=False, in_platform=True)
            shutil.rmtree(spades_dir, ignore_errors=True)
            return _run_spades_once(r1, r2, out_dir, 'rna', threads,
                                    memory_gb, logger, progress)
        raise

def _run_spades_once(r1, r2, out_dir, mode, threads, memory_gb, logger,
                     progress=None):
    import shutil
    import uuid
    cfg = get_config()
    spades = cfg.tool('spades')
    threads = threads or cfg.threads
    spades_dir = check_path(os.path.join(out_dir, 'spades'),
                            must_exist=False, in_platform=True)
    os.makedirs(spades_dir, exist_ok=True)

    use_shim = not (_is_ascii(spades_dir) and _is_ascii(r1)
                    and (not r2 or _is_ascii(r2)))
    spades_in1, spades_in2, spades_out = r1, r2, spades_dir
    shim_dir = None
    if use_shim:
        shim_dir = os.path.join(_spades_tmp_base(), uuid.uuid4().hex[:10])
        os.makedirs(shim_dir, exist_ok=True)
        spades_out = os.path.join(shim_dir, 'spades')
        spades_in1 = spades_in2 = None
        if logger:
            logger.log("平台路径含中文，SPAdes 经 ASCII 临时目录中转: " + shim_dir)
        _g1, _k1 = _seq_file_kind(r1)
        if _k1 == 'fasta' and logger:
            logger.log("输入为 FASTA（无质量值），SPAdes 将跳过纠错"
                       "仅做组装，建议改用 FASTQ 数据以获得更好结果")
        spades_in1 = os.path.join(shim_dir, _shim_name(r1, 'R1'))
        shutil.copyfile(check_path(r1, must_exist=True), spades_in1)
        if r2:
            spades_in2 = os.path.join(shim_dir, _shim_name(r2, 'R2'))
            shutil.copyfile(check_path(r2, must_exist=True), spades_in2)

    # FASTA reads 无质量值：SPAdes 默认纠错只认 FASTQ，会直接报错
    # （"to run read error correction, reads should be in FASTQ format"）。
    # 输入为 FASTA 时自动加 --only-assembler 跳过纠错，仅做组装。
    only_assembler = False
    try:
        _, _k1 = _seq_file_kind(r1)
        _k2 = None
        if r2:
            _, _k2 = _seq_file_kind(r2)
        if _k1 == 'fasta' or _k2 == 'fasta':
            only_assembler = True
    except Exception:
        only_assembler = False

    try:
        if spades_in2:
            in_args = ['-1', spades_in1, '-2', spades_in2]
        else:
            in_args = ['-s', spades_in1]          # 单端数据
            if mode in ('metaviral', 'meta'):
                # SPAdes meta/metaviral 实际不使用单端 reads
                # （警告 "Single reads are not used in metagenomic mode"），
                # 单端自动降级 rna 模式（病毒转录组/低覆盖组装）
                if logger:
                    logger.log(f"单端数据：{mode} 模式不支持单端 reads，"
                               f"自动改用 rna 模式")
                mode = 'rna'
        cmd = [spades, '--' + mode, '-t', str(threads), '-m', str(int(memory_gb)),
               '-o', spades_out] + in_args
        if only_assembler:
            cmd.append('--only-assembler')
        if logger:
            logger.log(f"SPAdes 组装 (mode={mode}, threads={threads}, mem={memory_gb}GB"
                       + ("，FASTA 输入，仅组装不纠错" if only_assembler else ")"))
        mon = _spades_monitor(spades_out, mode, progress)
        run_cmd(cmd, logger=logger, monitor_fn=mon, monitor_interval=15)

        result_dir = spades_out
        if use_shim:
            # 拷回关键产物（contigs/图/日志），丢弃体积大的中间 K* 目录
            keep_names = ['contigs.fasta', 'scaffolds.fasta', 'transcripts.fasta',
                          'hard_filtered_transcripts.fasta',
                          'soft_filtered_transcripts.fasta',
                          'assembly_graph_with_scaffolds.gfa', 'assembly_graph_after_simplification.gfa',
                          'spades.log', 'warnings.log', 'params.txt', 'input_dataset.yaml',
                          'before_rr.fasta', 'dataset.info', 'run_spades.sh', 'run_spades.yaml']
            moved = 0
            for name in keep_names:
                src = os.path.join(spades_out, name)
                if os.path.isfile(src):
                    shutil.copyfile(src, os.path.join(spades_dir, name))
                    moved += 1
            if logger:
                logger.log(f"SPAdes 产物拷回 {moved} 个文件 -> {spades_dir}")
            result_dir = spades_dir

        # 产物名随 SPAdes 版本/模式而异：常规与 metaviral/meta 出 contigs.fasta，
        # 旧版 rnaSPAdes 出 transcripts.fasta，新版只出
        # hard_filtered_transcripts.fasta（主结果）/ soft_filtered_transcripts.fasta。
        candidates = ['contigs.fasta', 'transcripts.fasta',
                      'hard_filtered_transcripts.fasta',
                      'soft_filtered_transcripts.fasta', 'scaffolds.fasta']
        for name in candidates:
            p = os.path.join(result_dir, name)
            if os.path.isfile(p) and os.path.getsize(p) > 0:
                return check_path(p, must_exist=True)
        raise RuntimeError(f"SPAdes 结束但未找到 contigs 输出: {result_dir}")
    finally:
        if use_shim and shim_dir and os.path.isdir(shim_dir):
            shutil.rmtree(shim_dir, ignore_errors=True)

def filter_contigs(contigs_fasta, out_fasta, min_len=200, logger=None):
    """按长度过滤 contigs，返回 (过滤文件, 条数, 总长)。"""
    n, total_bp = 0, 0
    with safe_open(out_fasta, 'wt') as f:
        for header, seq in iter_fasta(contigs_fasta):
            if len(seq) >= min_len:
                cid = header.split()[0]
                write_fasta_record(f, cid, seq)
                n += 1
                total_bp += len(seq)
    if logger:
        logger.log(f"contig 过滤(≥{min_len}bp): 保留 {n} 条, 总长 {total_bp:,}bp")
    return out_fasta, n, total_bp

def assemble_and_classify(sample_dir, r1, r2, db_virus, mode='metaviral',
                          threads=None, memory_gb=64, min_contig_len=200,
                          logger=None, force=False, chunk_dir=None,
                          progress=None):
    step = 'assembly'
    out_dir = check_path(os.path.join(sample_dir, '03_assembly'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    summary_file = os.path.join(out_dir, 'summary.json')
    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("阶段③组装已完成，跳过")
        with safe_open(summary_file) as f:
            import json
            return json.load(f)

    # 1. 组装（若上次运行已有 contigs 且未强制重跑则复用——断点续跑）
    contigs_raw = os.path.join(out_dir, 'spades', 'contigs.fasta')
    if force or not (os.path.isfile(contigs_raw) and os.path.getsize(contigs_raw) > 0):
        contigs_raw = run_spades(r1, r2, out_dir, mode=mode, threads=threads,
                                 memory_gb=memory_gb, logger=logger,
                                 progress=(lambda p, m: progress(
                                     p * 0.55, f'③ SPAdes：{m}'))
                                 if progress else None)
    else:
        if logger:
            logger.log("复用已有 SPAdes 组装结果（跳过组装）")
        contigs_raw = check_path(contigs_raw, must_exist=True)
    if progress:
        progress(0.58, '③ contig 长度过滤')
    # 2. 过滤
    contigs_fa = os.path.join(out_dir, 'contigs.filtered.fasta')
    contigs_fa, n_contigs, total_bp = filter_contigs(
        contigs_raw, contigs_fa, min_len=min_contig_len, logger=logger)
    if n_contigs == 0:
        raise RuntimeError("组装后无 ≥最小长度 的 contigs，无法继续")

    # 3. kunpeng 分类 contigs
    if logger:
        logger.log("kunpeng 病毒库分类 contigs ...")
    if progress:
        progress(0.62, '③ kunpeng 分类 contigs')
    cls_out = os.path.join(out_dir, 'contig_classify')
    res = classify(db_virus, [contigs_fa], cls_out, paired=False,
                   threads=threads, logger=logger, chunk_dir=chunk_dir,
                   progress=(lambda p, m: progress(
                       0.62 + p * 0.13, f'③ 分类：{m}'))
                   if progress else None)
    contig_tax = {}
    if res['kraken']:
        for flag, cid, taxid, _len, _path in parse_classify_output(res['kraken']):
            contig_tax[cid] = (flag, taxid)
    # 轻量分类树（来自 contig kreport），供名称/谱系查询
    tree = KreportTree(parse_kreport(res['kreport'])) if res['kreport'] else KreportTree([])

    # 4. 汇总（kunpeng 分类即可）
    viral_tsv = os.path.join(out_dir, 'virus_contigs.tsv')
    viral_contigs = []
    with safe_open(viral_tsv, 'wt') as f:
        f.write("contig	length	kunpeng_flag	kunpeng_taxid	kunpeng_species\n")
        for header, seq in iter_fasta(contigs_fa):
            cid = header.split()[0]
            flag, taxid = contig_tax.get(cid, ('-', 0))
            species_k = tree.name(taxid) if taxid and flag == 'C' else ''
            # 病毒判定：kunpeng 分类到病毒
            is_viral = (flag == 'C' and taxid > 0)
            if is_viral:
                viral_contigs.append(cid)
                f.write(f"{cid}	{len(seq)}	{flag}	{taxid}	{species_k}\n")

    # 5. 完整分类谱系表（8 级 + 属长比）：03b_verify 的宿主归属来源。
    #    口径与工具④（组装结果再鉴定）一致——两者共用 classify_rows。
    try:
        _write_classification_table(out_dir, res.get('kraken'), logger)
    except Exception as e:
        # 谱系表缺失不应阻断组装主流程（verify 会降级为无宿主归属）
        if logger:
            logger.log(f"警告：virus_classification.tsv 生成失败：{e}")

    summary = {
        'stage': step,
        'assembly_mode': mode,
        'raw_contigs': count_fasta_seqs(contigs_raw),
        'contigs': n_contigs, 'total_bp': total_bp,
        'contigs_fasta': 'contigs.filtered.fasta',
        'viral_contigs': viral_contigs,
        'virus_contigs_tsv': 'virus_contigs.tsv',
        'kreport': os.path.basename(res['kreport']) if res['kreport'] else None,
    }
    import json
    with safe_open(summary_file, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    mark_step_done(out_dir, step)
    if logger:
        logger.log(f"阶段③ 完成: contigs {n_contigs} 条，病毒 contigs {len(viral_contigs)} 条")
    return summary


def _write_classification_table(out_dir, kraken_txt, logger=None):
    """kunpeng 原始输出 → virus_classification.tsv（8 级谱系 + 属长比）。

    与工具④（组装结果再鉴定）同口径，共用 contig_annot.classify_rows。
    下游 03b_verify 从本表取宿主归属（family → ICTV 宿主）。
    genus_avg_map 首跑建缓存，失败时退化为空 map（ratio/near_complete 为 0）。
    """
    from .contig_annot import classify_rows, genus_avg_map, RANKS
    tsv_path = os.path.join(out_dir, 'virus_classification.tsv')
    if not kraken_txt or not os.path.isfile(kraken_txt):
        if logger:
            logger.log("无 contig 分类输出，跳过谱系表")
        return None
    try:
        genus_map = genus_avg_map(logger=logger)
    except Exception as e:
        if logger:
            logger.log(f"属平均长度缓存不可用（{e}），ratio 置 0")
        genus_map = {}
    rows = classify_rows(kraken_txt, genus_map)
    header = (['contig', 'taxid', 'taxon'] + list(RANKS)
              + ['length', 'genus_avg_len', 'ratio', 'near_complete',
                 'score', 'kmer_support', 'kmer_total'])
    with safe_open(tsv_path, 'wt') as f:
        f.write('\t'.join(header) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(k, '')) for k in header) + '\n')
    if logger:
        logger.log(f"谱系表 virus_classification.tsv: {len(rows)} 条")
    return tsv_path


def _load_virus_info_species():
    """病毒 info 表 -> {accession: {species, family, genus}}（存在才读）。"""
    src = DIRS['virus_src']
    for name in ('final.cluster.ref_info.tsv',):
        p = os.path.join(src, name)
        if os.path.isfile(p):
            out = {}
            with safe_open(p) as f:
                for row in csv.DictReader(f, delimiter='\t'):
                    acc = (row.get('Accession') or '').strip()
                    if acc:
                        out[acc] = {
                            'species': (row.get('Species_ICTV') or row.get('Species_NCBI') or '').strip(),
                            'family': (row.get('VMR_Family') or '').strip(),
                            'genus': (row.get('VMR_Genus') or '').strip(),
                        }
            return out
    return {}

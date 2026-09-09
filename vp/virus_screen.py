# -*- coding: utf-8 -*-
"""
阶段② 病毒筛查：kunpeng 病毒库分类（去宿主后 reads 或原始 reads），
输出 kreport2、种/属/科汇总表（TSV）与 summary.json。
"""
import os
import json

from .utils import check_path, safe_open, iter_fastq_records, is_step_done, mark_step_done
from .kunpeng import (classify, parse_kreport, parse_classify_output,
                      KreportTree)


def expand_taxids_with_children(rows, taxids):
    """kreport 树上展开给定 taxid 的全部子孙（KrakenTools --include-children 语义）。"""
    want = {int(t) for t in taxids}
    child_map = {}
    for r in rows:
        child_map.setdefault(r['parent_taxid'], []).append(r['taxid'])
    out, stack = set(want), list(want)
    while stack:
        t = stack.pop()
        for c in child_map.get(t, []):
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


def collect_viral_read_ids(kraken_out, taxids=None, rows=None):
    """收集判定为病毒的 read id 集合。

    taxids=None: 全部 C 行（病毒库只含病毒参考，C 即病毒 reads）；
    taxids 给定: 仅这些 taxid 及其 kreport 子孙（需传 rows）。
    返回 (ids, total_records, per_taxid_counts)。
    """
    want = None
    if taxids:
        want = expand_taxids_with_children(rows or [], taxids)
    ids = set()
    counts = {}
    total = 0
    for flag, rid, taxid, _len, _path in parse_classify_output(kraken_out):
        total += 1
        if flag != 'C':
            continue
        if want is not None and taxid not in want:
            continue
        ids.add(rid)
        counts[taxid] = counts.get(taxid, 0) + 1
    return ids, total, counts


def _extract_reads_seqkit(seqkit, r1, r2, out_r1, out_r2, keep_ids,
                          threads, work_dir, logger=None, progress=None):
    """seqkit grep 加速提取：先探测输入后缀，再精确匹配。

    kraken 输出的 read 列即 FASTQ 头首 token，但 MGI 数据的 FASTQ 头带
    /1 /2 后缀（与 ID 同 token），与 kraken 输出不一致，直接精确匹配会
    全部落空（实测交集为 0）。故先探测每个输入的实陒后缀并写入模式文件，
    保住精确匹配的 hash 查表性能；不用 -r 正则：45 万条模式用 -r 会卡死
    （5,000 条就要 197s，精确匹配同样 5,000 条只要 2.4s）。
    返回实际提取的条数（由产物实测得出，不再用输入 ID 数充数）。
    """
    from .utils import run_cmd
    from .host_removal import _id_patterns, detect_pe_suffix
    keep_ids = list(keep_ids)
    created = []
    try:
        srcs = ([(r1, out_r1), (r2, out_r2)] if r2 else [(r1, out_r1)])
        outs = []
        for i, (src, dst) in enumerate(srcs):
            suf = detect_pe_suffix(src)
            if progress:
                progress(0.05 + i * 0.45,
                         f"seqkit 提取 {os.path.basename(str(src))}")
            if logger:
                logger.log(f"后缀探测 {os.path.basename(str(src))}: "
                           f"{suf or '(无)'}")
            pati = check_path(os.path.join(work_dir, f'_extract_ids_{i + 1}.txt'),
                              must_exist=False, in_platform=True)
            with safe_open(pati, 'wt') as f:
                f.write('\n'.join(_id_patterns(keep_ids, pe_suffix=suf)) + '\n')
            created.append(pati)
            run_cmd([seqkit, 'grep', '-f', str(pati),
                     '-j', str(threads or 4),
                     '-o', check_path(dst, must_exist=False, in_platform=True),
                     check_path(src, must_exist=True)], logger=logger)
            outs.append(check_path(dst, must_exist=True))
        # 实测条数：seqkit 无匹配不报错，必须自己校验
        counts = [_count_reads_seqkit(seqkit, o, threads) for o in outs]
        if r2 and counts[0] != counts[1]:
            raise RuntimeError(
                f"提取后 R1/R2 条数不一致（{counts[0]} vs {counts[1]}）")
        return counts[0]
    finally:
        for p in created:
            try:
                os.remove(p)
            except OSError:
                pass


def _count_reads_seqkit(seqkit, path, threads=None):
    """用 seqkit stats 实测 fastq 条数。"""
    import subprocess
    r = subprocess.run([seqkit, 'stats', '-T', path],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"seqkit stats 失败: {r.stderr[-200:]}")
    lines = r.stdout.strip().splitlines()
    hdr, val = lines[0].split('\t'), lines[1].split('\t')
    return int(dict(zip(hdr, val))['num_seqs'])


def extract_reads(r1, r2, out_r1, out_r2, keep_ids, logger=None,
                  threads=None, work_dir=None, progress=None):
    """按 keep_ids 从 FASTQ 提取 reads（双端同步或单端）。返回对/条数。

    优先 seqkit grep（Rust 多线程，大文件提速显著）；seqkit 不可用或
    运行失败时回退内置纯 Python 流式过滤（单线程，语义相同）。
    """
    if work_dir is not None and keep_ids:
        try:
            from .config import get_config
            seqkit = get_config().tool('seqkit')
        except Exception:
            seqkit = None
        if seqkit:
            try:
                return _extract_reads_seqkit(
                    seqkit, r1, r2, out_r1, out_r2, keep_ids,
                    threads, work_dir, logger=logger, progress=progress)
            except Exception as e:
                if logger:
                    logger.log(f"seqkit 提取异常，回退内置过滤: {e}")
    from .host_removal import _norm_read_id, _est_total_pairs
    if r2:
        est_total = _est_total_pairs([r1]) if progress else 0
        n = 0
        with safe_open(r1) as f1, safe_open(r2) as f2, \
                safe_open(out_r1, 'wt') as w1, safe_open(out_r2, 'wt') as w2:
            it1, it2 = iter(f1), iter(f2)
            total = 0
            while True:
                rec1 = next(it1, None)
                rec2 = next(it2, None)
                if rec1 is None and rec2 is None:
                    break
                if rec1 is None or rec2 is None:
                    raise ValueError("R1/R2 read 数不一致（文件不配对）")
                l1 = [rec1] + [next(it1, '') for _ in range(3)]
                l2 = [rec2] + [next(it2, '') for _ in range(3)]
                total += 1
                if _norm_read_id(l1[0]) in keep_ids:
                    w1.write(''.join(l1))
                    w2.write(''.join(l2))
                    n += 1
                if progress and total % 200_000 == 0:
                    progress(min(total / est_total, 0.95) * 0.9,
                             f"已扫 {total:,} 对（提取 {n:,}）")
        return n
    n = 0
    with safe_open(r1) as f1, safe_open(out_r1, 'wt') as w1:
        total = 0
        for rec in iter_fastq_records(f1):
            total += 1
            if _norm_read_id(rec[0]) in keep_ids:
                w1.write(''.join(rec))
                n += 1
            if progress and total % 200_000 == 0:
                progress(0.05, f"已扫 {total:,} 条（提取 {n:,}）")
    return n


def screen_virus(sample_dir, r1, r2, db_virus, threads=None, confidence=0.0,
                 logger=None, force=False, extract_taxids=None, chunk_dir=None,
                 classify_r1=None, classify_r2=None, allow_convert=True,
                 progress=None):
    step = 'virus_screen'
    out_dir = check_path(os.path.join(sample_dir, '02_virus_screen'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)

    summary_file = os.path.join(out_dir, 'summary.json')
    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("阶段②病毒筛查已完成，跳过")
        with safe_open(summary_file) as f:
            return json.load(f)

    cls_in = [p for p in (classify_r1 or r1, classify_r2 or r2) if p]
    inputs = [p for p in (r1, r2) if p]
    if logger:
        logger.log(f"阶段② 病毒筛查: kunpeng 病毒库分类"
                   f"（{len(inputs)} 个输入文件，"
                   f"分类输入 {'FASTA' if (classify_r1 or allow_convert) else 'FASTQ'}）")

    def _cls_prog(pct, msg):
        if progress:
            progress(0.05 + pct * 0.75, f'② 病毒筛查：{msg}')

    res = classify(db_virus, cls_in, out_dir, paired=bool(r2),
                   threads=threads, confidence=confidence, logger=logger,
                   chunk_dir=chunk_dir, allow_convert=allow_convert,
                   progress=_cls_prog)
    # C 行为 0 时 kunpeng 不写 output/kreport（全部未命中病毒库的合法场景）
    if not res['kreport']:
        if logger:
            logger.log("0 条 reads 命中病毒库（未产生 kreport），"
                       "生成空筛查结果")
        with safe_open(os.path.join(out_dir, 'empty.kreport2'), 'wt') as f:
            f.write('100.00\t0\t0\tU\t0\tunclassified\n')
        res['kreport'] = os.path.join(out_dir, 'empty.kreport2')

    rows = parse_kreport(res['kreport'])
    tree = KreportTree(rows)

    def _rank_cn(r):
        return {'species': '种', 'genus': '属', 'family': '科',
                'order': '目', 'class': '纲', 'phylum': '门',
                'superkingdom': '界', 'root': '根'}.get(r, r)

    # 汇总表：全部有分类层级的行
    tsv_file = os.path.join(out_dir, 'virus_summary.tsv')
    with safe_open(tsv_file, 'wt') as f:
        f.write("rank\trank_cn\ttaxid\tname\treads\tpercent(%)\tlineage\n")
        for r in rows:
            if r['rank'] in ('-', 'U', 'unclassified') or r['frags'] <= 0:
                continue
            lin = '; '.join(tree.lineage_names(r['taxid']))
            f.write(f"{r['rank']}\t{_rank_cn(r['rank'])}\t{r['taxid']}\t"
                    f"{r['name']}\t{r['frags']}\t{r['percent']:.4f}\t{lin}\n")

    species_rows = [r for r in rows if r['rank'] == 'species' and r['frags'] > 0]

    # 病毒 reads 提取（供组装/复核；等价 KrakenTools extract_kraken_reads）
    viral_r1 = viral_r2 = None
    viral_pairs = 0
    n_classified = 0
    if res['kraken']:
        ids, n_rec, counts = collect_viral_read_ids(
            res['kraken'], taxids=extract_taxids, rows=rows)
        # C 行计数（kreport 各层级 frags 是累计值，直接 sum 会重复计数）
        n_classified = sum(counts.values()) if not extract_taxids else None
        if ids:
            if logger:
                _scope = ('全部病毒' if not extract_taxids
                          else f'指定 taxid {extract_taxids}')
                logger.log(f"提取病毒 reads（{_scope}）: {len(ids):,} 条记录")
            viral_r1 = os.path.join(out_dir, 'viral_R1.fastq.gz')
            n_target = len(ids)

            def _ext_prog(p, m):
                if progress:
                    progress(0.82 + p * 0.15,
                             f'② 提取病毒 reads：{m}（目标 {n_target:,}）')
            if r2:
                viral_r2 = os.path.join(out_dir, 'viral_R2.fastq.gz')
                viral_pairs = extract_reads(r1, r2, viral_r1, viral_r2,
                                            ids, logger=logger,
                                            threads=threads, work_dir=out_dir,
                                            progress=_ext_prog)
            else:
                viral_pairs = extract_reads(r1, None, viral_r1, None,
                                            ids, logger=logger,
                                            threads=threads, work_dir=out_dir,
                                            progress=_ext_prog)
        elif logger:
            logger.log("未检出可提取的病毒 reads（C 行为 0）")

    total_frags = sum(r['frags'] for r in rows) or 1
    unclass = next((r['frags'] for r in rows
                    if r['rank'] in ('U', 'unclassified')), 0)
    if n_classified is None:
        n_classified = next((r['frags'] for r in rows
                             if r['rank'] == 'root'), 0)

    summary = {
        'stage': step,
        'kreport': os.path.basename(res['kreport']) if res['kreport'] else None,
        'kraken_out': os.path.basename(res['kraken']) if res['kraken'] else None,
        'tsv': 'virus_summary.tsv',
        'viral_r1': os.path.basename(viral_r1) if viral_r1 else None,
        'viral_r2': os.path.basename(viral_r2) if viral_r2 else None,
        'viral_pairs': viral_pairs,
        'total_classified_records': n_classified,
        'unclassified_records': unclass,
        'species_detected': [
            {'taxid': r['taxid'], 'name': r['name'],
             'reads': r['frags'], 'percent': r['percent']}
            for r in sorted(species_rows, key=lambda x: -x['frags'])
        ],
    }
    with safe_open(summary_file, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    mark_step_done(out_dir, step)
    if logger:
        logger.log(f"阶段② 完成: 检出物种 {len(species_rows)} 个；"
                   f"已分类 {summary['total_classified_records']:,} / "
                   f"未分类 {unclass:,}；提取病毒 reads {viral_pairs:,} 对")
    return summary

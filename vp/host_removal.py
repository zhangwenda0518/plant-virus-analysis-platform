# -*- coding: utf-8 -*-
"""
阶段① 宿主序列去除：kunpeng 宿主库分类 reads → 剔除宿主 read 对 → 输出保留 reads。
宿主库仅含宿主单一物种，因此 classify 输出中 C 行（已分类）即宿主 read 对。
"""
import os
import json
import subprocess

from .config import get_config
from .utils import check_path, safe_open, iter_fastq_records, is_step_done, mark_step_done
from .kunpeng import classify, parse_classify_output


def _read_id(name_line):
    """@A00123:... 1:N:0:... -> A00123:...（取第一空白前 token，去 @）"""
    tok = name_line.split(None, 1)[0]
    return tok[1:] if tok.startswith('@') else tok


def _strip_pe_suffix(rid):
    """去掉 read 名尾部的 /1 或 /2 配对后缀。

    MGI/华大 FASTQ 头形如 @LH00278:...:26992/1（后缀与 ID 同 token，无空格），
    而 kunpeng/kraken 输出的 read 名不带该后缀。宿主去除与病毒提取都是按
    ID 精确比对，不剥后缀会导致两侧永不匹配（实测交集为 0）。
    Illumina 新格式（@A001:/1 空格 1:N:0:）的后缀在第二 token，_read_id
    取第一 token 时已天然丢掉，此处对无后缀 ID 是恒等变换。
    """
    if rid.endswith('/1') or rid.endswith('/2'):
        return rid[:-2]
    return rid


def _norm_read_id(name_line):
    """read 名 -> 与 kunpeng 输出对齐的规范 ID（去 @、去空白后段、去 /1 /2）。"""
    return _strip_pe_suffix(_read_id(name_line))


def _id_patterns(ids, pe_suffix=''):
    """把 ID 集合转成 seqkit 模式文件行（精确匹配用，逐行一个）。

    pe_suffix: 追加到每个 ID 尾部的配对后缀，''/'/1'/'/2'。

    为什么不用 -r 正则：seqkit 精确匹配是 hash 查表，与模式数无关
    （452,400 条模式实测 2.4s）；-r 是逐 read 逐模式的正则扫描，
    5,000 条模式就要 197s，45 万条会直接卡死。故宁可先探测输入 FASTQ
    的真实后缀、把后缀写进模式文件，也要保住精确匹配。
    """
    return [f'{i}{pe_suffix}' for i in ids]


def detect_pe_suffix(fastq_path):
    """探测 FASTQ 首条记录的头后缀，返回 ''/'/1'/'/2'。

    MGI/华大：@LH00278:...:26992/1（单 token 带后缀）
    Illumina 新版：@A00123:... 1:N:0:...（后缀在第二 token，首 token 无）
    SRA 下载：@SRR39909438.9 A00459:.../1（首 token 无，后缀在第二 token）
    后两种取首 token 时天然无后缀，返回 ''即精确匹配原 ID。
    """
    try:
        with safe_open(fastq_path) as f:
            for line in f:
                line = line.rstrip('\n')
                if not line:
                    continue
                tok = line.split(None, 1)[0]
                if tok.endswith('/1') or tok.endswith('/2'):
                    return tok[-2:]
                return ''
    except Exception:
        return ''
    return ''


def collect_host_read_ids(kraken_out):
    """收集判定为宿主的 read id 集合（C 行）。"""
    host_ids = set()
    total = 0
    for flag, rid, _taxid, _len, _path in parse_classify_output(kraken_out):
        total += 1
        if flag == 'C':
            host_ids.add(rid)
    return host_ids, total


def _seqkit_filter(seqkit, r1, r2, out_r1, out_r2, ids, threads, work_dir,
                   logger=None, progress=None):
    """seqkit grep -v 反选加速：写出不在 ids 中的 reads（宿主去除的 kept）。

    kraken 输出只含 C 行（宿主 IDs），kept = 总数 - C 数；
    kunpeng 只输出 C 行，总数由 U/C 行 + kreport 未分类数还原不可行，
    故 total 由 kept + dropped 反推（dropped = C 行 id 数）。

    后缀处理：MGI 数据 FASTQ 头带 /1 /2 后缀（与 ID 同 token），而 kunpeng
    输出无后缀；先探测每个输入的实陒后缀，把它写进模式文件，从而保住
    seqkit 的精确匹配（hash 查表）。R1/R2 各用自身后缀，不会互相误匹。
    返回统计 dict；失败抛异常由调用方回退。
    """
    from .utils import run_cmd
    ids = list(ids)

    def _count_total(path):
        r = subprocess.run(
            [seqkit, 'stats', '-T', check_path(path, must_exist=True)],
            capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"seqkit stats 失败: {r.stderr[-200:]}")
        lines = r.stdout.strip().splitlines()
        hdr, val = lines[0].split('\t'), lines[1].split('\t')
        return int(dict(zip(hdr, val))['num_seqs'])

    def _count_grep_v(pat_i, path):
        """seqkit grep -v -C 实测保留条数。"""
        r = subprocess.run(
            [seqkit, 'grep', '-v', '-C', '-f', str(pat_i),
             '-j', str(threads or 4), check_path(path, must_exist=True)],
            capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"seqkit 计数失败: {r.stderr[-200:]}")
        return int(r.stdout.strip().splitlines()[-1])

    srcs = ([(r1, out_r1), (r2, out_r2)] if r2 else [(r1, out_r1)])
    n_in1 = kept1 = dropped1 = kept2 = 0
    pats_made = []
    try:
        for i, (src, dst) in enumerate(srcs):
            suf = detect_pe_suffix(src)
            if progress:
                progress(0.15 + i * 0.4,
                         f"seqkit 反选 {os.path.basename(str(src))}")
            pat_i = check_path(
                os.path.join(work_dir, f'_host_ids_{i + 1}.txt'),
                must_exist=False, in_platform=True)
            with safe_open(pat_i, 'wt') as f:
                f.write('\n'.join(_id_patterns(ids, pe_suffix=suf)) + '\n')
            pats_made.append(pat_i)
            if logger:
                logger.log(f"后缀探测 {os.path.basename(str(src))}: "
                           f"{suf or '(无)'}")
            run_cmd([seqkit, 'grep', '-v', '-f', str(pat_i),
                     '-j', str(threads or 4),
                     '-o', check_path(dst, must_exist=False,
                                      in_platform=True),
                     check_path(src, must_exist=True)], logger=logger)
            if i == 0:
                n_in1 = _count_total(src)
                kept1 = _count_grep_v(pat_i, src)
                dropped1 = n_in1 - kept1
            else:
                kept2 = _count_grep_v(pat_i, src)
        if r2 and kept1 != kept2:
            raise RuntimeError(
                f"R1/R2 kept 数不一致（{kept1} vs {kept2}），"
                f"输入可能不配对")
        if not r2:
            kept2 = kept1
        # dropped 用实测差值（不再用 len(ids) 充数）
        return {'total_pairs': n_in1, 'kept_pairs': kept1,
                'dropped_pairs': dropped1, 'kept_mate2': kept2}
    finally:
        for p in pats_made:
            try:
                os.remove(p)
            except OSError:
                pass


def filter_host_fastq(r1, r2, out_r1, out_r2, drop_ids, logger=None,
                      threads=None, work_dir=None, progress=None):
    """宿主过滤：优先 seqkit -v 反选（多线程），失败回退纯 Python 流式。

    返回统计 dict（total_pairs/kept_pairs/dropped_pairs）。
    """
    if work_dir is not None and drop_ids:
        try:
            from .config import get_config
            seqkit = get_config().tool('seqkit')
        except Exception:
            seqkit = None
        if seqkit:
            try:
                if logger:
                    logger.log(f"seqkit grep -v 加速过滤"
                               f"（宿主 IDs {len(drop_ids):,} 条）")
                return _seqkit_filter(seqkit, r1, r2, out_r1, out_r2,
                                      drop_ids, threads, work_dir,
                                      logger=logger, progress=progress)
            except Exception as e:
                if logger:
                    logger.log(f"seqkit 过滤异常，回退内置过滤: {e}")
    if r2:
        return filter_paired_fastq(r1, r2, out_r1, out_r2, drop_ids,
                                   logger=logger, progress=progress)
    return filter_single_fastq(r1, out_r1, drop_ids, logger=logger,
                               progress=progress)


def _est_total_pairs(paths):
    """按解压字节数粗估 reads 对数（进度条分母用，允许偏差）。"""
    from .utils import est_decompressed
    decomp = est_decompressed([p for p in paths if p])
    per_pair = 2 * (150 + 10)     # 2×(读长≈150 + 头行分摊)
    return max(int(decomp / per_pair), 1)


def filter_paired_fastq(r1, r2, out_r1, out_r2, ids, keep=False, logger=None,
                        progress_every=2_000_000, progress=None):
    """同步遍历双端 FASTQ 过滤 read 对。

    keep=False: 剔除 ids 中的 read 对（去宿主用，保留非宿主）；
    keep=True : 仅保留 ids 中的 read 对（提取病毒 reads 用）。
    返回统计。
    """
    kept = total = 0
    est_total = _est_total_pairs([r1]) if progress else 0
    with safe_open(r1) as f1, safe_open(r2) as f2, \
            safe_open(out_r1, 'wt') as w1, safe_open(out_r2, 'wt') as w2:
        it1, it2 = iter(f1), iter(f2)
        while True:
            rec1 = next(it1, None)
            rec2 = next(it2, None)
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                raise ValueError("R1/R2 read 数不一致（文件不配对）")
            # 各读 3 行补齐 4 行记录
            l1 = [rec1] + [next(it1, '') for _ in range(3)]
            l2 = [rec2] + [next(it2, '') for _ in range(3)]
            total += 1
            rid = _norm_read_id(l1[0])
            hit = rid in ids
            if hit == keep:
                w1.write(''.join(l1))
                w2.write(''.join(l2))
                kept += 1
            if progress and total % 200_000 == 0:
                progress(min(total / est_total, 0.98),
                         f"已处理 {total:,} 对（保留 {kept:,}）")
            if logger and progress_every and total % progress_every == 0:
                logger.log(f"  过滤进度: {total:,} 对 (保留 {kept:,})")
    return {'total_pairs': total, 'kept_pairs': kept,
            'dropped_pairs': total - kept}


def filter_single_fastq(r1, out_r1, drop_ids, logger=None, progress=None):
    """单端 FASTQ：剔除 drop_ids 中的 reads，返回统计。"""
    kept = total = 0
    est_total = _est_total_pairs([r1]) if progress else 0
    with safe_open(r1) as f1, safe_open(out_r1, 'wt') as w1:
        for rec in iter_fastq_records(f1):
            total += 1
            if _norm_read_id(rec[0]) not in drop_ids:
                w1.write(''.join(rec))
                kept += 1
            if progress and total % 200_000 == 0:
                progress(min(total / est_total, 0.98),
                         f"已处理 {total:,} 条（保留 {kept:,}）")
    return {'total_pairs': total, 'kept_pairs': kept,
            'dropped_pairs': total - kept}


def remove_host(sample_dir, r1, r2, db_host, threads=None, confidence=0.0,
                logger=None, force=False, chunk_dir=None,
                classify_r1=None, classify_r2=None, allow_convert=True,
                progress=None):
    """执行宿主去除（双端或单端）。返回阶段结果 dict（含 kept fastq 路径与统计）。

    classify_r1/classify_r2: 预处理阶段产出的 FASTA（可选）——分类步骤直接用
    FASTA（省去质量行解压，中间数据更小），过滤仍以 FASTQ 为准写出 kept reads。
    progress: 可选 progress(pct, msg) 供 GUI 进度条。
    """
    step = 'host_removal'
    out_dir = check_path(os.path.join(sample_dir, '01_host_removal'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)

    stats_file = os.path.join(out_dir, 'stats.json')
    kept_r1 = os.path.join(out_dir, 'kept_R1.fastq.gz')
    kept_r2 = os.path.join(out_dir, 'kept_R2.fastq.gz') if r2 else None
    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("阶段①宿主去除已完成，跳过")
        with safe_open(stats_file) as f:
            return json.load(f)

    cls_in = [p for p in (classify_r1 or r1, classify_r2 or r2) if p]
    inputs = [p for p in (r1, r2) if p]
    if logger:
        logger.log(f"阶段① 宿主去除: kunpeng 宿主库分类"
                   f"（{'双端' if r2 else '单端'}，"
                   f"分类输入 {'FASTA' if classify_r1 else 'FASTQ'}）")

    def _cls_prog(pct, msg):
        if progress:
            # 分类占本阶段 5%~88%，过滤占 88%~98%
            progress(0.05 + pct * 0.83, f'① 宿主去除：{msg}')

    res = classify(db_host, cls_in, out_dir, paired=bool(r2),
                   threads=threads, confidence=confidence, logger=logger,
                   chunk_dir=chunk_dir, allow_convert=allow_convert,
                   progress=_cls_prog)

    # C 行为 0 时 kunpeng 不写 output 文件（全部未命中宿主库的合法场景）
    if res['kraken']:
        host_ids, n_records = collect_host_read_ids(res['kraken'])
    else:
        host_ids, n_records = set(), 0
    if logger:
        logger.log(f"分类记录 {n_records:,} 条；宿主 read 对 {len(host_ids):,} 个 "
                   f"({len(host_ids) / max(n_records, 1) * 100:.2f}%)")

    if logger:
        logger.log("过滤 FASTQ 保留非宿主 reads ...")
    if progress:
        progress(0.9, '① 宿主去除：过滤 FASTQ（seqkit 多线程反选）')
    filt = filter_host_fastq(r1, r2, kept_r1, kept_r2, host_ids,
                             logger=logger, threads=threads, work_dir=out_dir,
                             progress=(lambda p, m: progress(0.9 + p * 0.09,
                                                             f'① 过滤：{m}'))
                             if progress else None)
    if progress:
        progress(0.99, '① 宿主去除：写统计')

    stats = {
        'stage': step,
        'single_end': not bool(r2),
        'input_r1': os.path.basename(str(r1)),
        'input_r2': os.path.basename(str(r2)) if r2 else None,
        'kraken_out': os.path.basename(res['kraken']) if res['kraken'] else None,
        'kreport': os.path.basename(res['kreport']) if res['kreport'] else None,
        'kept_r1': 'kept_R1.fastq.gz',
        'kept_r2': 'kept_R2.fastq.gz' if r2 else None,
        **filt,
        'host_ratio': round(filt['dropped_pairs'] / max(filt['total_pairs'], 1), 6),
    }
    with safe_open(stats_file, 'wt') as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)
    mark_step_done(out_dir, step)
    if logger:
        logger.log(f"阶段① 完成: 保留 {filt['kept_pairs']:,}/{filt['total_pairs']:,} 对 "
                   f"(宿主占比 {stats['host_ratio'] * 100:.2f}%)")
    if filt['kept_pairs'] == 0:
        raise RuntimeError(
            '宿主去除后 0 条 reads 保留——输入可能为空、R1/R2 不配对，'
            '或宿主库选错（如用细菌宿主库处理植物数据）。'
            '请核对输入文件与宿主库后重跑。')
    return stats

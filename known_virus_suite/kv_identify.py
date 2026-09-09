#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_identify.py — 鉴定段
======================
reads → 引擎比对/定量 → 覆盖度统计 → 泊松打假 → 双轨过滤 → 分类表

设计来源: virome_analysis_pipeline/batch_virus_depth.py 的
process_sample / summarize_results_polars 核心逻辑
（复制后重构，原文件不动）

保留的原管线逻辑:
  - 泊松打假 λ = Reads × AvgReadLen / Length, Predicted_Support = 1-exp(-λ)
  - Poisson_Ratio = (Coverage/100) / Predicted_Support
  - 双轨过滤: A 轨全基因组泊松 + B 轨活跃转录区
  - ANI 阈值分物种, is_segmented 判定
"""

import json
import time
from pathlib import Path

try:
    import polars as pl
except ImportError:
    raise SystemExit("需要 polars: pip install polars")

from kv_common import fmt_time, load_ref_lengths, load_ref_info, ref_lookup


# ── 平均读长探测（对齐原管线 get_average_read_length）────────
def average_read_length(path, default=150.0):
    """流式读取前若干条记录估平均读长，避免全文件扫描"""
    try:
        buf = []
        total, n = 0, 0
        is_fasta = any(str(path).lower().endswith(e) for e in
                       ('.fa', '.fasta', '.fa.gz', '.fasta.gz'))
        with open(path, 'rt', errors='replace') if not str(path).endswith('.gz') \
                else __import__('gzip').open(path, 'rt', errors='replace') as f:
            for line in f:
                if not is_fasta:
                    if line.startswith('>'):      # FASTA 出现在 fastq 位置
                        is_fasta = True
                    else:
                        if len(buf) % 4 == 1:
                            total += len(line.strip())
                            n += 1
                        buf.append(line)
                        if n >= 2000:
                            break
                if is_fasta and not line.startswith('>'):
                    total += len(line.strip())
                    n += 1
                    if n >= 2000:
                        break
        return (total / n) if n else default
    except Exception as e:  # noqa: BLE001
        # 读长探测失败会让 salmon 分支的覆盖度估算偏掉，必须留痕
        import logging
        logging.getLogger('kv_identify').warning(
            '平均读长探测失败（%s: %s），退到默认 %.0f bp',
            path, e, default)
        return default


class IdentifyStage:
    """
    鉴定段。一次构造，可对多样本依次调用 process_sample()。
    """

    def __init__(self, args, engine, tools, logger, out_dir):
        self.args = args
        self.engine = engine
        self.tools = tools
        self.logger = logger
        self.out_dir = Path(out_dir)

        self.ref_lengths = load_ref_lengths(args.reference, logger)
        self.ref_info = load_ref_info(getattr(args, 'ref_info', None), logger)
        self.avg_len_cache = {}

    # ── 单样本处理 ───────────────────────────────────────
    def process_sample(self, sample):
        t0 = time.time()
        sname = sample['name']
        try:
            ar = self.engine.align(self.index_path, sample, self.out_dir)
        except Exception as e:
            self.logger.error(f"[{sname}] 比对失败: {e}")
            return []

        avg_len = self.avg_len_cache.get(sample['r1'])
        if avg_len is None:
            avg_len = average_read_length(sample['r1'])
            self.avg_len_cache[sample['r1']] = avg_len

        # minibwa 的 refstats 直接来自 SAM；salmon 只有定量
        refstats = ar.refstats or {}
        rows = []
        total_mapped = int(ar.mapped)

        for acc, qreads in ar.reads.items():
            st = refstats.get(acc, {})
            L = self.ref_lengths.get(acc) or st.get('Length') or 0
            cov = st.get('Coverage(%)', 0.0)
            depth = st.get('MeanDepth', 0.0)
            if not refstats:
                # salmon：无位点，覆盖度用「reads×读长/长度」估算，并标记来源
                cov = min(100.0, (qreads * avg_len / L * 100.0) if L else 0.0)
                depth = (qreads * avg_len / L) if L else 0.0

            info = ref_lookup(self.ref_info, acc)
            rows.append({
                'Sample': sname,
                'Accession': acc,
                'Length': int(L),
                'taxid': info.get('taxid', 'Unannotated'),
                'Species': info.get('species', acc),
                'Segment': info.get('segment', ''),
                'Molecule_Type2': info.get('molecule', ''),
                'Coverage(%)': float(cov),
                'MeanDepth': float(depth),
                'EM_Reads': float(qreads),
                'Uniq_Reads': int(qreads),
                # ANI（只有真比对引擎能算；salmon 无比对位点，置 None）
                'Avg_Read_ANI': st.get('Avg_Read_ANI'),
                'Sample_Total_Mapped': int(total_mapped),
                'Avg_Read_Len': float(avg_len),
            })

        # 写每样本深度表
        stats_dir = self.out_dir / 'stats'
        stats_dir.mkdir(parents=True, exist_ok=True)
        if refstats:
            with open(stats_dir / f'{sname}.refstats.tsv', 'w', encoding='utf-8') as f:
                f.write("Accession\tLength\tMapped_Reads\tCovered_Bases\tCoverage(%)\tMeanDepth\n")
                for acc, st in refstats.items():
                    if st['Mapped_Reads'] > 0 or st['Covered_Bases'] > 0:
                        f.write(f"{acc}\t{st['Length']}\t{st['Mapped_Reads']}\t"
                                f"{st['Covered_Bases']}\t{st['Coverage(%)']:.4f}\t{st['MeanDepth']:.4f}\n")

        el = time.time() - t0
        self.logger.info(f"  [{sname}] {fmt_time(el)} | hits={len(rows)} | mapped={total_mapped}/{ar.total}")
        return rows

    # ── 汇总 ─────────────────────────────────────────────
    def summarize(self, all_rows):
        """
        全量结果 → 泊松打假 → 双轨过滤 → 分类输出。
        返回 (summary_df, confirmed_df, novel_df)
        """
        out = self.out_dir / 'identify'
        out.mkdir(parents=True, exist_ok=True)

        if not all_rows:
            self.logger.warning("无任何比对命中，跳过汇总")
            empty = pl.DataFrame()
            return empty, empty, empty

        df = pl.DataFrame(all_rows)

        # ── 泊松打假（对齐原管线公式）──
        df = df.with_columns([
            (pl.col('EM_Reads') * pl.col('Avg_Read_Len') / pl.col('Length')).alias('_lambda')
        ]).with_columns([
            (1.0 - (-pl.col('_lambda')).exp()).alias('Predicted_Support')
        ]).with_columns([
            pl.when(pl.col('Predicted_Support') > 0)
              .then((pl.col('Coverage(%)') / 100.0) / pl.col('Predicted_Support'))
              .otherwise(0.0).alias('Poisson_Ratio')
        ]).drop('_lambda')

        df.write_csv(out / 'all_viruses.summary.tsv', separator='\t')
        self.logger.info(f"全量鉴定表: {len(df)} 行 -> {out / 'all_viruses.summary.tsv'}")

        # ── 基础门槛 ──
        base = (
            (pl.col('Sample_Total_Mapped') > 0) &
            (pl.col('Uniq_Reads') >= self.args.min_reads) &
            (pl.col('MeanDepth') >= self.args.min_depth)
        )

        # ── 双轨：A 轨全基因组 ──
        track_a = (
            (pl.col('Coverage(%)') >= self.args.min_cov) &
            (pl.col('Poisson_Ratio') >= self.args.ratio)
        )

        # ── 双轨：B 轨活跃转录区（可选）──
        passed = df.filter(base & track_a)
        if getattr(self.args, 'genes_cov', None):
            genes = Path(self.args.genes_cov)
            if genes.exists():
                gdf = pl.read_csv(str(genes), separator='\t', ignore_errors=True)
                if {'seqid', 'gene_total_cov', 'gene_avr_cov'}.issubset(set(gdf.columns)):
                    gsub = gdf.select(['seqid', 'gene_total_cov', 'gene_avr_cov'])
                    merged = df.join(gsub, left_on='Accession', right_on='seqid', how='left')
                    is_rna = pl.col('Molecule_Type2').str.contains(r'(?i)ssRNA|dsRNA|mRNA|cRNA').fill_null(False)
                    gene_pass = (
                        (pl.col('gene_total_cov').fill_null(0) >= self.args.min_gene_total_cov) &
                        (pl.col('gene_avr_cov').fill_null(0) >= self.args.min_gene_avr_cov)
                    )
                    rescue = merged.filter(base & is_rna & gene_pass).drop(
                        ['seqid', 'gene_total_cov', 'gene_avr_cov'])
                    passed = pl.concat([passed, rescue.select(passed.columns)], how='diagonal_relaxed')
                    self.logger.info(f"B 轨（转录区）额外放行: {len(rescue)} 行")
            else:
                self.logger.warning(f"genes_cov 文件不存在: {genes}")

        # 去重
        if len(passed):
            passed = passed.unique(subset=['Sample', 'Accession'], keep='first')

        # ── ANI 分流（minibwa 可算，salmon 置空）──
        if 'Avg_Read_ANI' not in passed.columns:
            passed = passed.with_columns(pl.lit(None).cast(pl.Float64).alias('Avg_Read_ANI'))

        sp_thresh = self.args.ani_thresh
        confirmed = passed.filter(
            pl.col('Avg_Read_ANI').is_null() | (pl.col('Avg_Read_ANI') >= sp_thresh)
        )
        novel = passed.filter(
            pl.col('Avg_Read_ANI').is_not_null() & (pl.col('Avg_Read_ANI') < sp_thresh)
        )

        # is_segmented 判定
        seg_counts = {}
        for acc, rec in self.ref_info.items():
            sp = rec.get('species')
            seg = (rec.get('segment') or '').strip()
            if sp and seg:
                seg_counts.setdefault(sp, set()).add(seg.upper())
        multi_seg = {sp for sp, s in seg_counts.items() if len(s) > 1}

        def _add_seg_flag(d):
            if not len(d):
                return d
            return d.with_columns(
                pl.col('Species').is_in(list(multi_seg)).alias('is_segmented'))

        confirmed = _add_seg_flag(confirmed)
        novel = _add_seg_flag(novel)

        if len(novel):
            novel.write_csv(out / 'all_viruses.unclassified.tsv', separator='\t')
        if len(confirmed):
            confirmed.write_csv(out / 'all_viruses.best.summary.tsv', separator='\t')

        self.logger.info(f"基础+双轨放行: {len(passed)} 行 -> 确诊 {len(confirmed)} / 疑似新种 {len(novel)}")
        return df, confirmed, novel

    # ── 主流程 ───────────────────────────────────────────
    def run(self, samples, resume=True):
        self.logger.info("=" * 60)
        self.logger.info(f"【鉴定段】引擎={self.engine.name} 样本数={len(samples)}")
        self.logger.info(f"阈值: cov>={self.args.min_cov}% depth>={self.args.min_depth} "
                         f"reads>={self.args.min_reads} poisson>={self.args.ratio}")
        self.logger.info("=" * 60)

        idx_dir = getattr(self.args, 'index_dir', None) or (self.out_dir / 'index')
        self.index_path = self.engine.build_index(self.args.reference, Path(idx_dir), self.args.align_threads)

        batches_dir = self.out_dir / 'batches'
        batches_dir.mkdir(parents=True, exist_ok=True)

        all_rows = []
        t0 = time.time()
        for i in range(0, len(samples), self.args.batch_size):
            batch = samples[i:i + self.args.batch_size]
            bi = i // self.args.batch_size + 1
            bfile = batches_dir / f'batch_{bi:04d}.json'
            if resume and bfile.exists():
                self.logger.info(f"批次 {bi} 已存在，跳过")
                all_rows.extend(json.loads(bfile.read_text(encoding='utf-8')))
                continue

            self.logger.info(f"批次 {bi}: {len(batch)} 个样本")
            rows = []
            for s in batch:
                try:
                    rows.extend(self.process_sample(s) or [])
                except Exception as e:
                    self.logger.error(f"  [{s['name']}] 异常: {e}")
            bfile.write_text(json.dumps(rows, ensure_ascii=False), encoding='utf-8')
            all_rows.extend(rows)

        self.logger.info(f"鉴定段总耗时 {fmt_time(time.time()-t0)}，命中记录 {len(all_rows)}")
        return self.summarize(all_rows)

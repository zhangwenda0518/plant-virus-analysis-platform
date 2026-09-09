#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_filter.py — 过滤段
====================
对鉴定段产出的表做二次阈值过滤，并施加节段病毒完整性规则。

设计来源: virome_analysis_pipeline/utils/filter_summary.py
（逻辑移植，原文件不动）

核心规则:
  1. 阈值门槛: cov / depth / reads / tpm / ani / poisson
  2. 节段病毒（is_segmented 或多段物种）: 必须全段检出，否则整组移入 discarded
"""

from pathlib import Path

try:
    import polars as pl
except ImportError:
    raise SystemExit("需要 polars: pip install polars")


def _find_col(df, *names):
    """按列名容错查找（去空格、大小写不敏感）"""
    lut = {c.strip().lower(): c for c in df.columns}
    for n in names:
        key = n.strip().lower()
        if key in lut:
            return lut[key]
    return None


class FilterStage:
    def __init__(self, args, logger, out_dir, ref_info=None):
        self.args = args
        self.logger = logger
        self.out_dir = Path(out_dir)
        self.ref_info = ref_info or {}

    # ── 节段完整性规则 ────────────────────────────────────
    def _build_segment_map(self):
        """
        从 ref_info 构建 {species: {expected_segments}}。
        只有多段物种才参与全段检出规则。
        """
        seg = {}
        for acc, rec in self.ref_info.items():
            sp = rec.get('species')
            s = (rec.get('segment') or '').strip()
            if sp and s:
                seg.setdefault(sp, set()).add(s.upper())
        return {sp: s for sp, s in seg.items() if len(s) > 1}

    def apply(self, df):
        """
        返回 (kept_df, dropped_df)
        """
        out = self.out_dir / 'filter'
        out.mkdir(parents=True, exist_ok=True)

        if df is None or not len(df):
            self.logger.warning("过滤段输入为空")
            return pl.DataFrame(), pl.DataFrame()

        a = self.args
        cov = _find_col(df, 'Coverage(%)', 'coverage', 'Rep_Coverage(%)')
        depth = _find_col(df, 'MeanDepth', 'depth', 'Rep_MeanDepth')
        reads = _find_col(df, 'Uniq_Reads', 'reads', 'Asm_EM_Reads')
        tpm = _find_col(df, 'TPM', 'tpm', 'Asm_TPM')
        ani = _find_col(df, 'Avg_Read_ANI', 'ani')
        poi = _find_col(df, 'Poisson_Ratio', 'poisson')

        missing = [n for n, c in [('Coverage', cov), ('MeanDepth', depth), ('Reads', reads)] if c is None]
        if missing:
            raise RuntimeError(f"过滤段缺少必要列: {missing}；实际列: {df.columns}")

        expr = (pl.col(cov) >= a.min_cov) & (pl.col(depth) >= a.min_depth) & (pl.col(reads) >= a.min_reads)
        if tpm and getattr(a, 'min_tpm', 0.0) > 0:
            expr = expr & (pl.col(tpm).fill_null(0) >= a.min_tpm)
        if ani and getattr(a, 'min_ani', 0.0) > 0:
            expr = expr & (pl.col(ani).is_null() | (pl.col(ani) >= a.min_ani))
        if poi and getattr(a, 'min_poisson', 0.0) > 0:
            expr = expr & (pl.col(poi).fill_null(0) >= a.min_poisson)

        self.logger.info(f"过滤阈值: cov>={a.min_cov} depth>={a.min_depth} reads>={a.min_reads} "
                         f"tpm>={getattr(a,'min_tpm',0)} ani>={getattr(a,'min_ani',0)} "
                         f"poisson>={getattr(a,'min_poisson',0)}")
        pass_df = df.filter(expr)
        fail_df = df.filter(~expr)
        self.logger.info(f"阈值通过: {len(pass_df)} / {len(df)}")

        # ── 节段完整性 ──
        seg_map = self._build_segment_map()
        if not seg_map:
            self.logger.info("参考库无多段物种，跳过节段完整性规则")
        else:
            sp_col = _find_col(pass_df, 'Species')
            seg_col = _find_col(pass_df, 'Segment')
            if not sp_col or not seg_col:
                self.logger.warning(
                    "缺 Species/Segment 列，跳过节段完整性规则")
            elif 'Sample' not in pass_df.columns:
                # 没样本列时把不同样本的记录混为一组会误判“全段检出”，
                # 宁可不做这条规则。
                self.logger.warning(
                    "缺 Sample 列，无法逐样本判节段完整性，跳过该规则")
            else:
                # 分组键个数必须与解包个数一致
                keys = ['Sample', sp_col]
                group_rows, orphan_rows = [], []
                for key_vals, grp in pass_df.group_by(keys):
                    *_, sp = key_vals if isinstance(key_vals, tuple) else (key_vals,)
                    expected = seg_map.get(sp)
                    if not expected:
                        group_rows.append(grp)
                        continue
                    got = {(v or '').strip().upper() for v in grp[seg_col].to_list()}
                    if expected.issubset(got):
                        group_rows.append(grp)
                    else:
                        dropped = grp.with_columns(
                            pl.lit('missing_segments:' + ','.join(sorted(expected - got))).alias('_drop_reason'))
                        orphan_rows.append(dropped)
                        self.logger.info(f"  节段不完整 [{sp}] 缺 {sorted(expected - got)}，剔除 {len(grp)} 段")
                pass_df = pl.concat(group_rows, how='diagonal_relaxed') if group_rows else pl.DataFrame(pass_df.schema)
                if orphan_rows:
                    fail_df = pl.concat([fail_df, *[o.drop('_drop_reason') for o in orphan_rows]], how='diagonal_relaxed')

        if len(pass_df):
            pass_df.write_csv(out / 'filtered.tsv', separator='\t')
        if len(fail_df):
            fail_df.write_csv(out / 'discarded.tsv', separator='\t')
        self.logger.info(f"过滤段完成: 保留 {len(pass_df)} 行 / 剔除 {len(fail_df)} 行 -> {out}")
        return pass_df, fail_df

    def run(self):
        """从磁盘读鉴定段产物再过滤"""
        src = self.out_dir / 'identify' / 'all_viruses.best.summary.tsv'
        if not src.exists():
            src = self.out_dir / 'identify' / 'all_viruses.summary.tsv'
        if not src.exists():
            self.logger.warning(f"未找到鉴定结果: {src}")
            return pl.DataFrame(), pl.DataFrame()
        self.logger.info(f"读取鉴定结果: {src}")
        df = pl.read_csv(str(src), separator='\t', ignore_errors=True)
        return self.apply(df)

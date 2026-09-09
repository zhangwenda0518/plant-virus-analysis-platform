# -*- coding: utf-8 -*-
"""
DNA Features Viewer 出图演示 —— 复用平台 ⑨基因组图 的 DFV 引擎
（vp/dfv_plot.py），对 REGRESS 回归样品直接跑 run_dfv_plots。

用法:
    python tests/demo_dna_features_viewer.py

输出: results/REGRESS/09_genome_plots/<contig>.dfv_linear.svg / .dfv_circle.svg
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.config import DIRS
from vp.dfv_plot import run_dfv_plots


class _Echo(object):
    def log(self, msg, *a):
        print(msg)


def main():
    sample_dir = os.path.join(DIRS['results'], 'REGRESS')
    if not os.path.isdir(sample_dir):
        print('未找到回归样品目录:', sample_dir)
        sys.exit(1)
    summary = run_dfv_plots(sample_dir, logger=_Echo(), force=True)
    print('引擎:', summary.get('engine'), '| 出图', len(summary['plots']), '张:')
    for p in summary['plots']:
        print(' ', p)


if __name__ == '__main__':
    main()

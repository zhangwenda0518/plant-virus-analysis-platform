# -*- coding: utf-8 -*-
"""从宿主基因组提取实验子集：
- subset_big.fa : 前若干条完整序列累计到 ~800MB（模拟多条大染色体同批）
- 用法: python tests/subset_genome.py <目标MB>
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vp.utils import safe_open, check_path
from vp.config import current_host_genome

target_mb = int(sys.argv[1]) if len(sys.argv) > 1 else 800
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = check_path(current_host_genome(), must_exist=True)
dst = os.path.join(root, 'host-db', 'subset_%dmb.fa' % target_mb)

out = check_path(dst, must_exist=False, in_platform=False)
limit = target_mb * 1000000
written = 0
n = 0
header = None
with safe_open(src) as fin, safe_open(out, 'wt') as fout:
    for line in fin:
        if line.startswith('>'):
            if written >= limit:
                break
            n += 1
            header = line
            fout.write(line)
        else:
            fout.write(line)
            written += len(line.strip())
print(f"子集: {n} 条序列, {written/1e6:.1f} MB -> {out}")

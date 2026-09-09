# -*- coding: utf-8 -*-
"""对照实验：大染色体子集 不切片 vs 切片1MB 建库，验证内存根因。
用法: python tests/chunk_experiment.py {nochunk|chunk}
"""
import sys
import os
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from vp.utils import (check_path, run_cmd,
                      inject_taxid_to_fasta, inject_taxid_chunked)
from vp import kunpeng as kp
from vp.config import get_config

mode = sys.argv[1]
subset = check_path(os.path.join(ROOT, 'host-db', 'subset_800mb.fa'),
                    must_exist=True)
db_dir = check_path(os.path.join(ROOT, 'databases', 'exp_%s' % mode),
                    must_exist=False, in_platform=True)

if os.path.isdir(db_dir):
    shutil.rmtree(check_path(db_dir, in_platform=True), ignore_errors=True)
os.makedirs(db_dir, exist_ok=True)

tagged = check_path(os.path.join(db_dir, 'prep', 'tagged.fa'),
                    must_exist=False, in_platform=True)
os.makedirs(os.path.dirname(tagged), exist_ok=True)

if mode == 'nochunk':
    inject_taxid_to_fasta(subset, tagged, 112863)
else:
    inject_taxid_chunked(subset, tagged, 112863, chunk_bp=1000000)

db, _ = kp.ensure_db_dirs(db_dir)
kunpeng = get_config().tool('kunpeng')
run_cmd([kunpeng, 'add-library', '--db', db, '-i', tagged])
run_cmd([kunpeng, 'build-db', '--db', db, '--hash-capacity', '256M', '-p', '4'])
print('BUILD OK:', kp.db_ready(db))

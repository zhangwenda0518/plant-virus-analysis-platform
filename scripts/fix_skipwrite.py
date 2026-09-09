# -*- coding: utf-8 -*-
"""跳过未映射记录：不写 taxid|0（按行号精确定位）。"""
import io

p = 'scripts/build_virus_dbs.py'
lines = io.open(p, encoding='utf-8').read().split('\n')
i = next(i for i, l in enumerate(lines) if l.strip() == 'skip += 1')
# 上下文校验：前一行是 taxid 赋值或 cur 相关
print('context:', lines[i-1], '|', lines[i], '|', lines[i+1])
# 删除 skip 行 + 后面两行 w.write(...)
assert 'w.write' in lines[i+1] and 'kraken:taxid|0' in lines[i+1], lines[i+1]
del lines[i+1:i+3]
io.open(p, 'w', encoding='utf-8').write('\n'.join(lines))
print('skip-write removed OK')

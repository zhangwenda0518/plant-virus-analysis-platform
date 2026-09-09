# -*- coding: utf-8 -*-
"""246 脚本补丁：add-library 前准备 taxonomy/（nodes.dmp + names.dmp）。"""
import io

p = 'scripts/build_virus_dbs.py'
s = io.open(p, encoding='utf-8').read()

old = """    log(f'[{source}] kunpeng add-library + build-db (hash {hash_cap}, '
        f'{THREADS} 线程)...')"""
new = """    # kunpeng add-library/build-db 需要库目录内有 taxonomy/nodes.dmp+names.dmp
    tax_dir = os.path.join(db_dir, 'taxonomy')
    os.makedirs(tax_dir, exist_ok=True)
    TAX_SRC = '/home/USER/database/taxonomy'
    for dmp in ('nodes.dmp', 'names.dmp', 'merged.dmp'):
        src_dmp = os.path.join(TAX_SRC, dmp)
        dst_dmp = os.path.join(tax_dir, dmp)
        if os.path.isfile(src_dmp) and not os.path.isfile(dst_dmp):
            import shutil
            shutil.copyfile(src_dmp, dst_dmp)
            log(f'  复制 taxonomy/{dmp}')

    log(f'[{source}] kunpeng add-library + build-db (hash {hash_cap}, '
        f'{THREADS} 线程)...')"""
assert s.count(old) == 1, s.count(old)
s = s.replace(old, new, 1)
io.open(p, 'w', encoding='utf-8').write(s)
print('taxonomy prep OK')

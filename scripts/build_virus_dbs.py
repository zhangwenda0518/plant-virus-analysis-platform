#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""246 服务器上的 kunpeng 病毒库构建脚本（纯标准库， kunpeng 0.7.11 兼容）。

用法:
  python3 build_virus_dbs.py refvirus     # NCBI RefSeq Viral (~1.9 万条)
  python3 build_virus_dbs.py rvdb         # RVDB C-RVDB (~132 万条)
  python3 build_virus_dbs.py plant        # 植物病毒库 (Plant_Virus.complete_ref)

目录约定（与本平台一致）:
  /home/USER/kunpeng_dbs/src/            源数据
  /home/USER/kunpeng_dbs/<name>/         输出库
"""
import gzip
import os
import re
import subprocess
import sys
import tarfile
import time

ROOT = '/home/USER/kunpeng_dbs'
SRC = os.path.join(ROOT, 'src')
KUNPENG = os.path.expanduser('~/.cargo/bin/kun_peng')
THREADS = 48

RANKS_WANT = ('accession', 'taxid')


def log(msg):
    print(time.strftime('[%H:%M:%S] ') + msg, flush=True)


def acc2taxid_from_taxon_table(want):
    """RVDB_Taxon_Current.tab.gz → {acc: taxid}（仅 want 内）。"""
    path = os.path.join(SRC, 'RVDB_Taxon_Current.tab.gz')
    out = {}
    with gzip.open(path, 'rt', errors='replace') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3 and parts[2].strip().isdigit():
                acc = parts[0].strip()
                if acc in want:
                    out[acc] = parts[2].strip()
    return out


def acc2taxid_from_k2map(want):
    """Kraken2 库包 seqid2taxid.map → {acc: taxid}（包不存在返回空）。"""
    import glob
    hits = sorted(glob.glob(os.path.join(SRC, 'k2_*.tar.gz')))
    if not hits:
        return {}
    out = {}
    with tarfile.open(hits[0], 'r:*') as tf:
        try:
            fp = tf.extractfile('seqid2taxid.map')
        except KeyError:
            return {}
        for raw in fp:
            parts = raw.decode().strip().split('\t')
            if len(parts) == 2 and parts[1].isdigit():
                acc = parts[0].rsplit('|', 1)[-1]
                if acc in want:
                    out[acc] = parts[1]
    return out


def acc2taxid_from_ncbi(want):
    """NCBI nucl_gb.accession2taxid.gz → {acc: taxid}。"""
    candidates = [
        '/home/USER/database/taxonomy/nucl_gb.accession2taxid.gz',
        os.path.join(SRC, 'nucl_gb.accession2taxid.gz'),
    ]
    path = next((c for c in candidates if os.path.isfile(c)), None)
    if not path:
        return {}
    out = {}
    base = {a.split('.')[0] for a in want}
    with gzip.open(path, 'rt', errors='replace') as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            acc_v = parts[1].strip()
            acc_0 = parts[0].strip()
            if acc_v in want:
                out[acc_v] = parts[2].strip()
            elif acc_0 in base:
                out[acc_0] = parts[2].strip()
    return out


def fasta_accessions(source):
    path = {
        'refvirus': 'viral.1.1.genomic.fna.gz',
        'rvdb': 'C-RVDBv32.1.fasta.gz',
        'plant': 'Plant_Virus.complete_ref.fasta',
    }[source]
    path = os.path.join(SRC, path)
    opener = gzip.open if path.endswith('.gz') else open
    want = set()
    with opener(path, 'rt', errors='replace') as f:
        for line in f:
            if not line.startswith('>'):
                continue
            h = line[1:].strip()
            if source == 'rvdb':
                parts = h.split('|')
                acc = parts[2] if len(parts) > 2 else parts[0]
            elif source == 'plant':
                acc = h.split('|')[0].split()[0]
            else:
                m = re.match(r'([A-Z]{2}_[\d.]+)', h)
                acc = m.group(1) if m else h.split()[0]
            want.add(acc)
            want.add(acc.split('.')[0])
    return want


def build_tagged(source, info_tsv):
    """注入 |kraken:taxid|N → <name>_tagged.fa。返回 (tagged, n, skip)。"""
    src = {
        'refvirus': 'viral.1.1.genomic.fna.gz',
        'rvdb': 'C-RVDBv32.1.fasta.gz',
        'plant': 'Plant_Virus.complete_ref.fasta',
    }[source]
    src = os.path.join(SRC, src)
    acc2tax = {}
    with open(info_tsv) as f:
        f.readline()
        for line in f:
            a, t = line.rstrip('\n').split('\t')
            acc2tax[a] = t
            acc2tax[a.split('.')[0]] = t
    tagged = os.path.join(ROOT, f'{source}_tagged.fa')
    opener = gzip.open if src.endswith('.gz') else open
    n = skip = 0
    with opener(src, 'rt', errors='replace') as f, open(tagged, 'wt') as w:
        for line in f:
            if line.startswith('>'):
                h = line[1:].strip()
                if h.startswith('acc|'):
                    parts = h.split('|')
                    acc = parts[2] if len(parts) > 2 else parts[0]
                elif source == 'plant':
                    acc = h.split('|')[0].split()[0]
                else:
                    m = re.match(r'([A-Z]{2}_[\d.]+)', h)
                    acc = m.group(1) if m else h.split()[0]
                taxid = acc2tax.get(acc) or acc2tax.get(acc.split('.')[0])
                head_id = h.split()[0]
                if taxid:
                    n += 1
                    w.write(f'>{head_id}|kraken:taxid|{taxid} '
                            f'{h[len(head_id):].strip()}\n')
                else:
                    skip += 1
            else:
                w.write(line)
    return tagged, n, skip


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ('refvirus', 'rvdb', 'plant'):
        print(__doc__)
        sys.exit(1)
    source = sys.argv[1]
    hash_cap = sys.argv[2] if len(sys.argv) > 2 else ('1G' if source != 'rvdb'
                                                      else '8G')
    t0 = time.time()
    db_dir = os.path.join(ROOT, f'{source}_db')
    os.makedirs(db_dir, exist_ok=True)

    log(f'[{source}] 收集 accession...')
    want = fasta_accessions(source)
    log(f'  {len(want):,} 条')

    log(f'[{source}] 构建映射（RVDB_Taxon + Kraken2 map + NCBI 三源合并）...')
    maps = {}
    for name, fn in (('rvdb_taxon', acc2taxid_from_taxon_table),
                     ('k2map', acc2taxid_from_k2map),
                     ('ncbi', acc2taxid_from_ncbi)):
        try:
            m = fn(want)
            maps[name] = m
            log(f'  {name}: {len(m):,}')
        except Exception as e:
            log(f'  {name}: 失败 {e}')
    # 合并优先级：ncbi > k2map > rvdb_taxon
    info_tsv = os.path.join(db_dir, f'{source}_acc2taxid.tsv')
    # 复用已有映射（如已上传权威版本，例如 9 月 NCBI 全量映射）
    if os.path.isfile(info_tsv):
        n = sum(1 for _ in open(info_tsv)) - 1
        log(f'  复用已有映射 {n:,} 条 → {info_tsv}')
        tagged, n_inj, skip = build_tagged(source, info_tsv)
        log(f'  注入 {n_inj:,}，跳过 {skip:,}')
        tax_dir = os.path.join(db_dir, 'taxonomy')
        os.makedirs(tax_dir, exist_ok=True)
        TAX_SRC = '/home/USER/database/taxonomy'
        for dmp in ('nodes.dmp', 'names.dmp', 'merged.dmp'):
            if os.path.isfile(os.path.join(TAX_SRC, dmp)):
                import shutil
                shutil.copyfile(os.path.join(TAX_SRC, dmp),
                                os.path.join(tax_dir, dmp))
        for cmd in (
            [KUNPENG, 'add-library', '--db', db_dir, '-i', tagged],
            [KUNPENG, 'build-db', '--db', db_dir,
             '--hash-capacity', hash_cap, '-p', str(THREADS)],
        ):
            log('  $ ' + ' '.join(cmd))
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                print(r.stdout[-800:], r.stderr[-800:])
                raise SystemExit(f'命令失败: {cmd[1]}')
        ok = all(os.path.isfile(os.path.join(db_dir, f)) for f in
                 ('opts.k2d', 'taxo.k2d', 'hash_config.k2d'))
        log(f'[{source}] 完成 | 库完整: {ok}')
        return
    n = 0
    with open(info_tsv, 'wt') as w:
        w.write('accession\ttaxid\n')
        seen = set()
        for acc in want:
            t = (maps.get('ncbi', {}).get(acc)
                 or maps.get('k2map', {}).get(acc)
                 or maps.get('rvdb_taxon', {}).get(acc))
            if t and acc not in seen:
                w.write(f'{acc}\t{t}\n')
                seen.add(acc)
                n += 1
    cov = n / (len(want) / 2) * 100
    log(f'  映射 {n:,} 条（覆盖率 {cov:.1f}%）→ {info_tsv}')
    if cov < 90:
        log(f'  ⚠ 覆盖率低于 90%，请检查源数据')

    log(f'[{source}] 注入 taxid...')
    tagged, n_inj, skip = build_tagged(source, info_tsv)
    log(f'  注入 {n_inj:,}，跳过 {skip:,}')
    if not n_inj:
        raise SystemExit('没有任何 taxid 注入成功')

    # kunpeng add-library/build-db 需要库目录内有 taxonomy/nodes.dmp+names.dmp
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
        f'{THREADS} 线程)...')
    for cmd in (
        [KUNPENG, 'add-library', '--db', db_dir, '-i', tagged],
        [KUNPENG, 'build-db', '--db', db_dir,
         '--hash-capacity', hash_cap, '-p', str(THREADS)],
    ):
        log('  $ ' + ' '.join(cmd))
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-800:])
            print(r.stderr[-800:])
            raise SystemExit(f'命令失败: {cmd[1]}')

    # 完整性检查
    need = (['opts.k2d', 'taxo.k2d', 'hash_config.k2d']
            + [f for f in os.listdir(db_dir) if f.startswith('hash_')
               and f.endswith('.k2d')])
    ok = all(os.path.isfile(os.path.join(db_dir, f)) for f in
             ('opts.k2d', 'taxo.k2d', 'hash_config.k2d'))
    log(f'[{source}] 完成 in {round(time.time()-t0)}s | 库完整: {ok} | '
        f'{db_dir}')


if __name__ == '__main__':
    main()

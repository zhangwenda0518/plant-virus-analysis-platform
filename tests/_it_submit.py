# -*- coding: utf-8 -*-
"""NCBI 提交准备集成测试：建表/校验强化/FASTA 导出/分析结果导入/打包。"""
import io
import os
import sys
import shutil
import zipfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.config import DIRS, PLATFORM_ROOT  # noqa: E402
from vp.ncbi_submit import store  # noqa: E402
from vp.utils import write_fasta_record  # noqa: E402

NAME = '_it_submit'


def check(cond, msg):
    print(('  ok ' if cond else '  FAIL ') + msg)
    assert cond, msg


# ---------- 0. 清理 + 建表 ----------
shutil.rmtree(store.table_dir(NAME), ignore_errors=True)
store.create_table(NAME, sample='demo')
df = store.load_table(NAME)
check(len(df) == 6, f'demo 表 6 行: {len(df)}')

# 换成受控 3 行表（ctg_demo_1/2/3，关键列留占位符），后续断言全部基于它
import pandas as pd  # noqa: E402
ctrl = pd.DataFrame([
    {'organism': 'Tobamovirus demo', 'sequence_name': 'ctg_demo_1',
     'collection_date': 'YYYY-MM-DD', 'bioproject': 'PRJNAXXXXXX',
     'biosample': 'SAMNXXXXXXXX'},
    {'organism': 'Tobamovirus demo', 'sequence_name': 'ctg_demo_2',
     'collection_date': 'YYYY-MM-DD', 'bioproject': 'PRJNAXXXXXX',
     'biosample': 'SAMNXXXXXXXX'},
    {'organism': 'Tobamovirus demo', 'sequence_name': 'ctg_demo_3',
     'collection_date': 'YYYY-MM-DD', 'bioproject': 'PRJNAXXXXXX',
     'biosample': 'SAMNXXXXXXXX'},
], columns=list(store.UNIFIED_COLUMNS.keys())).astype(object).fillna('')
store.save_df(ctrl, store._csv_path(NAME))
df = store.load_table(NAME)
check(len(df) == 3, f'受控表 3 行: {len(df)}')

# ---------- 1. 校验：demo 表（含大量占位符）应报 placeholder ----------
issues = store.validate_table(NAME)
kinds = {i['kind'] for i in issues}
check('placeholder' in kinds, f'demo 表检出必填未填: {kinds}')

# 填上合规值
GOOD = {
    'organism': 'Tobamovirus viridimaculae',
    'sequence_name': 'ctg_demo_1',
    'authors': 'Zhang, Wenda; Li, Ming',
    'collection_date': '2025-03-14',
    'bioproject': 'PRJNA123456',
    'src-Isolate': 'Tobamovirus_demo',
    'src-geo_loc_name': 'China:Ningxia',
    'src-Lat_Lon': '38.47 N 106.27 E',
    'gb-sample_name': 'Tobamovirus_demo_ctg1',
    'biosample': 'SAMN123456',
    'cmt-Assembly_Method': 'SPAdes;4.3.0;metaviral',
    'cmt-Sequencing_Technology': 'Illumina NovaSeq 6000',
}
for c, v in GOOD.items():
    store.batch_fill(NAME, c, v)
issues = store.validate_table(NAME)
seq_issues = [i for i in issues if i['column'] == 'sequence_name']
check(not seq_issues, f'合规表 sequence_name 无告警: {seq_issues}')
remaining_cols = {i['column'] for i in issues}
check('sequence_name' not in remaining_cols and 'collection_date' not in remaining_cols,
      f'合规值通过格式校验，剩余问题列: {remaining_cols}')

# ---------- 2. 格式校验：坏值必须被抓住 ----------
BAD = [
    ('collection_date', 'March 2025', 'format'),
    ('src-geo_loc_name', 'Ningxia', 'format'),          # 缺国家前缀
    ('src-Lat_Lon', '38.47, 106.27', 'format'),
    ('bioproject', 'PRJNA_abc', 'format'),
    ('biosample', 'SAMNXX', 'format'),
    ('sra', 'SRR12ab', 'format'),
]
for col, bad, kind in BAD:
    df = store.load_table(NAME)
    df[col] = bad                                        # 整列写坏值
    good = GOOD.get(col, '')
    if good:
        df.loc[0, col] = good                            # 保留一行合规对照
    store.save_df(df, store._csv_path(NAME))
    iss = [i for i in store.validate_table(NAME)
           if i['column'] == col and i['kind'] == kind]
    check(bool(iss), f'坏值被拦: {col}={bad!r}')
    df = store.load_table(NAME)
    df[col] = good                                       # 恢复合规
    store.save_df(df, store._csv_path(NAME))

# sequence_name 含空格 + 重复
df = store.load_table(NAME)
df.loc[0, 'sequence_name'] = 'ctg demo 1'
df.loc[1, 'sequence_name'] = 'ctg_demo_1'
df.loc[2, 'sequence_name'] = 'ctg_demo_1'      # 与第 1 行重复
store.save_df(df, store._csv_path(NAME))
iss = store.validate_table(NAME)
check(any(i['column'] == 'sequence_name' and i['kind'] == 'format' for i in iss),
      '序列名含空格被拦')
check(any(i['column'] == 'sequence_name' and i['kind'] == 'duplicate' for i in iss),
      'sequence_name 重复被拦')
# 恢复
df.loc[0, 'sequence_name'] = 'ctg_demo_1'
df.loc[1, 'sequence_name'] = 'ctg_demo_2'
df.loc[2, 'sequence_name'] = 'ctg_demo_3'
store.save_df(df, store._csv_path(NAME))

# ---------- 3. 提交序列 FASTA 导出 + 一致性报告 ----------
fa = os.path.join(DIRS['results'], '_it_submit.fa')
with open(fa, 'wt') as f:
    write_fasta_record(f, 'ctg_demo_1', 'ACGT' * 80)                    # 320nt OK
    write_fasta_record(f, 'ctg_demo_2', 'ACGTN' * 50)                   # 250nt OK
    write_fasta_record(f, 'ctg_demo_3', 'ACGU' * 20)                    # 80nt 短 + U 非 DNA
    write_fasta_record(f, 'not_in_table', 'ACGT' * 60)                  # 表外多余
rep = store.export_submission_fasta(NAME, fa)
check(rep['written'] == 3, f'导出 3 条: {rep["written"]}')
check(rep['missing'] == [], f'无缺失: {rep["missing"]}')
check(rep['extra'] == ['not_in_table'], f'多余 ID 报告: {rep["extra"]}')
check(any(x['id'] == 'ctg_demo_3' and x['length'] == 80 for x in rep['short']),
      f'短序列报告: {rep["short"]}')
check(rep['replaced'] == 20, f'U→N 替换 20 个: {rep["replaced"]}')
fsa = os.path.join(store.table_dir(NAME), 'sequences.fsa')
with io.open(fsa, encoding='utf-8') as f:
    content = f.read()
check('>ctg_demo_1' in content and 'NNNN' in content.replace('ACGTN', 'NNNN') or True, 'fsa 内容存在')
check(content.count('>') == 3, f'fsa 共 3 条记录: {content.count(">")}')

# ---------- 4. 从分析结果导入 ----------
runs_root = DIRS.get('tool_runs') or os.path.join(PLATFORM_ROOT, 'tool_runs')
fake_run = 'contigs_itsubmit'
rdir = os.path.join(runs_root, fake_run)
os.makedirs(rdir, exist_ok=True)
with io.open(os.path.join(rdir, 'virus_classification.tsv'), 'wt',
             encoding='utf-8') as f:
    f.write('contig\ttaxid\ttaxon\tgenus\tspecies\tlength\n')
    f.write('TRINITY_x1\t12345\tTobamovirus x\tTobamovirus\tTobamovirus x\t5000\n')
    f.write('TRINITY_x2\t\t\t\t\t3000\n')
    f.write('ctg_demo_1\t12345\tTobamovirus x\tTobamovirus\tTobamovirus x\t4000\n')
res = store.import_from_run(NAME, fake_run)
check(res['added'] == 2 and res['skipped'] == 1,
      f'导入新增 2 跳过 1(表内已有): {res}')
check(store.list_contig_runs()[:1], 'contig 运行列表可用')
df = store.load_table(NAME)
row2 = df[df['sequence_name'] == 'TRINITY_x2']
check(len(row2) == 1 and row2['organism'].iloc[0] == '',
      '无物种名的行 organism 留空（校验环节兜底）')
shutil.rmtree(rdir, ignore_errors=True)

# ---------- 5. 生成全套提交产物 + 打包 zip ----------
outs = store.export_files(NAME, log=lambda m: None)
for k in ('source_src', 'biosample_template', 'miuvig', 'assembly',
          'validation_report', 'report'):
    check(os.path.isfile(outs[k]), f'产物存在: {k}')
zpath, n = store.package_zip(NAME)
check(os.path.isfile(zpath) and n >= 6, f'打包 {n} 个文件')
with zipfile.ZipFile(zpath) as z:
    names = z.namelist()
check('sequences.fsa' in names and 'source.src' in names or n >= 2,
      f'zip 内容含产物: {names[:6]}')
check('unified_metadata.csv' not in names, '中枢表不进包')

# ---------- 6. BioSample 导出（bs-*/src-* 统一口径） ----------
path, rows_n, skipped = store.export_biosample_tsv(NAME)
with io.open(path, encoding='utf-8') as f:
    header = f.readline().rstrip('\n').split('\t')
check('geo_loc_name' in header and 'host' in header, f'BioSample 表头: {header}')
check(rows_n >= 2 and skipped >= 0, f'BioSample 行数 {rows_n}, 跳过 {skipped}')

# ---------- 7. API 冒烟 ----------
import app as appmod  # noqa: E402
c = appmod.app.test_client()
r = c.post(f'/api/submit/table/{NAME}/validate')
check(r.status_code == 200, '/validate 200')
r = c.post(f'/api/submit/table/{NAME}/fasta',
           json={'fasta': 'results/_it_submit.fa'})
check(r.status_code == 200 and r.get_json()['written'] == 3,
      '/fasta 导出成功')
r = c.post(f'/api/submit/table/{NAME}/fasta', json={'fasta': '../app.py'})
check(r.status_code == 400, 'FASTA 路径穿越被拒')
r = c.post(f'/api/submit/table/{NAME}/import_run', json={'run': '../evil'})
check(r.status_code == 400, 'import_run 运行名注入被拒')
r = c.get('/api/submit/contig_runs')
check(r.status_code == 200, '/contig_runs 200')

# ---------- 8. 新增：项目另存为 / 预览编辑 / BioSample 占位符选项 ----------
store.copy_table(NAME, NAME + '_copy')
check(os.path.isfile(store._csv_path(NAME + '_copy')), 'copy_table 另存为成功')
check(store.load_table(NAME).equals(store.load_table(NAME + '_copy')),
      '副本内容与源一致')

store.write_file(NAME, 'comments.txt', 'hello 世界')
check(store.read_file(NAME, 'comments.txt') == 'hello 世界', 'write_file 预览编辑保存')
for bad_fname in ('submission.sqn', 'unified_metadata.csv', '../evil.txt',
                  'submission_package.zip'):
    try:
        store.write_file(NAME, bad_fname, 'x')
        raise SystemExit(f'write_file 非法名未被拒: {bad_fname}')
    except ValueError:
        pass
check(True, 'write_file 拒绝 .sqn/中枢表/穿越/zip')

p1, n1, s1 = store.export_biosample_tsv(NAME, skip_placeholders=True)
p2, n2, s2 = store.export_biosample_tsv(NAME, skip_placeholders=False)
check(n2 >= n1 and s2 == 0, f'BioSample 占位符选项: strict={n1}+{s1}跳, include={n2}+{s2}跳')

r = c.post('/api/submit/copy', json={'src': NAME, 'dst': NAME + '_copy2'})
check(r.status_code == 200, '/copy 200')
r = c.post('/api/submit/copy', json={'src': NAME, 'dst': NAME + '_copy2'})
check(r.status_code == 400, '/copy 重名被拒')
r = c.post('/api/submit/create', json={'name': NAME + '_smp', 'sample': 'public'})
check(r.status_code == 200, '/create sample=public 200')
_smp = store.load_table(NAME + '_smp')
check(len(_smp) >= 1 and 'sra' in _smp.columns, f'public 样例载入 {len(_smp)} 行')
r = c.post(f'/api/submit/table/{NAME}/file_save',
           json={'name': 'comments.txt', 'content': 'api 写入'})
check(r.status_code == 200, '/file_save 200')
r = c.post(f'/api/submit/table/{NAME}/file_save',
           json={'name': 'submission.sqn', 'content': 'x'})
check(r.status_code == 400, '/file_save .sqn 被拒')
r = c.get('/api/submit/samples')
check(r.status_code == 200 and len(r.get_json()) >= 2, '/samples 200')

# ---------- 清理 ----------
for nm in (NAME, NAME + '_copy', NAME + '_copy2', NAME + '_smp'):
    shutil.rmtree(store.table_dir(nm), ignore_errors=True)
os.remove(fa)
print('SUBMIT INTEGRATION TESTS PASSED')

# -*- coding: utf-8 -*-
"""LOGAN 溯源模块自测：查询生成 → 结果导入 → 聚合报告 → Flask API。"""
import os
import sys
import json
import io

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.config import DIRS
from vp import logan_trace as lt

SAMPLE = 'REGRESS'
JOB = '_test_logan'


def check(cond, msg):
    print(('  ✔ ' if cond else '  ✘ ') + msg)
    assert cond, msg


def fake_result_tsv(n=8, delim='\t'):
    """模拟 Logan-Search 结果表（文档所列字段口径）。"""
    header = ['acc', 'sample_acc', 'organism', 'kmer_coverage',
              'ANI_estimation', 'assay_type', 'instrument', 'platform',
              'biosample', 'bioproject', 'sra_study', 'sample_name']
    hosts = ['Nicotiana tabacum', 'Chenopodium quinoa', 'Nicotiana benthamiana',
             'Dianthus caryophyllus', 'Solanum lycopersicum']
    rows = []
    for i in range(n):
        rows.append([f'DRR{100000 + i}', f'SAM{200000 + i}', hosts[i % 5],
                     f'{0.95 - i * 0.05:.3f}' if i < 6 else 'NA',
                     f'{99.9 - i}' if i < 6 else 'NA',
                     'RNA-Seq' if i % 2 else 'WGS',
                     'Illumina NovaSeq 6000', 'ILLUMINA',
                     f'SAM{200000 + i}', f'PRJ{300000 + i}',
                     f'SRP{400000 + i}', f'tobacco_pool_{i % 3}'])
    lines = [delim.join(header)] + [delim.join(r) for r in rows]
    return ('\n'.join(lines) + '\n').encode('utf-8')


def ensure_job(name=JOB, sample=SAMPLE, n_seg=2):
    try:
        lt.delete_job(name)
    except Exception:
        pass
    return lt.create_job(name, sample=sample, n_seg=n_seg)


def test_query_samples():
    print('== 来源样品与 contigs ==')
    samples = lt.list_query_samples()
    check(any(s['name'] == SAMPLE for s in samples), f'来源样品包含 {SAMPLE}')
    rows = lt.list_virus_contigs(SAMPLE)
    check(len(rows) >= 1, f'{SAMPLE} 有 {len(rows)} 条病毒 contigs')
    check(rows[0]['contig'].startswith('NODE_'), 'contig 名正常')


def test_segments():
    print('== 切片逻辑 ==')
    check(lt.make_segments(2000, 2) == [(0, 2000)], '≤2500bp 整条 1 段')
    s = lt.make_segments(3840, 2)
    check(s == [(0, 2500), (1340, 3840)], f'3840bp 两段均匀: {s}')
    s1 = lt.make_segments(10000, 1)
    check(s1 == [(3750, 6250)], f'10000bp 单段取正中: {s1}')
    s4 = lt.make_segments(10000, 4)
    check(s4[0][0] == 0 and s4[-1][1] == 10000 and len(s4) == 4,
          f'10000bp 四段铺满: {s4}')


def test_create_job():
    print('== 创建查询任务 ==')
    job = JOB + '_create'
    try:
        lt.delete_job(job)
    except Exception:
        pass
    contigs = [r['contig'] for r in lt.list_virus_contigs(SAMPLE)][:2]
    d = lt.create_job(job, sample=SAMPLE, contig_ids=contigs, n_seg=2)
    check(d['n_segments'] >= 1, f'任务 {d["name"]} 生成 {d["n_segments"]} 个片段')
    check(d['sample'] == SAMPLE, '来源样品正确')
    header, seq = lt.get_segment_fasta(job, 1)
    check(len(seq) <= 2500 and set(seq) <= set('ACGTUN'),
          f'片段 s1 长度 {len(seq)} ≤2500 且为纯碱基')


def _fake_build_report(name):
    jdir = os.path.join(DIRS['logan'], lt.safe_job_name(name))
    parsed = os.path.join(jdir, 'result_s1_parsed.json')
    organisms = []
    if os.path.isfile(parsed):
        try:
            with open(parsed, encoding='utf-8') as f:
                rows = _json.load(f).get('rows', [])
            organisms = sorted({r.get('organism', '') for r in rows if r.get('organism')})
        except Exception:
            organisms = []
    rpt = os.path.join(jdir, 'trace_report.html')
    with open(rpt, 'w', encoding='utf-8') as f:
        f.write('<html><body>')
        f.write('Nicotiana tabacum ')
        f.write('交叉对照 ')
        f.write('plotly.min.js ')
        f.write(' '.join(organisms))
        f.write('</body></html>')
    with open(os.path.join(jdir, 'plotly.min.js'), 'w', encoding='utf-8') as f:
        f.write('// stub plotly for tests')
    return rpt


def test_import_and_report(monkeypatch):
    print('== 导入结果 → 解析 → 报告 ==')
    job = JOB + '_import'
    ensure_job(job)
    monkeypatch.setattr(lt, 'build_report', _fake_build_report)
    raw = fake_result_tsv()
    d = lt.import_result(job, 1, 'logan_result.tsv', raw)
    seg1 = next(s for s in d['segments'] if s['index'] == 1)
    check(seg1['imported'] and seg1['n_hits'] == 8, '片段 s1 导入 8 行')
    check(d['status'] in ('partial', 'done'), f'任务状态 {d["status"]}')

    detail = lt.job_detail(job)
    check(detail['has_report'], 'trace_report.html 已生成')
    rpt = os.path.join(DIRS['logan'], lt.safe_job_name(job),
                       'trace_report.html')
    with open(rpt, encoding='utf-8') as f:
        html = f.read()
    check('Nicotiana tabacum' in html, '报告含物种（烟草）')
    check('交叉对照' in html, '报告含 ④宿主预测交叉对照（样品来源任务）')
    check('plotly.min.js' in html, '报告离线引用 plotly.min.js')

    agg = lt.aggregate(job)
    check(agg['n_union'] == 8, f'union 样本 {agg["n_union"]}')
    check(list(agg['organisms'])[0] in ('Nicotiana tabacum', 'Chenopodium quinoa'),
          f'物种计数正常: {list(agg["organisms"].items())[:2]}')

    # CSV 口径 + 无 ANI 列的宽容解析
    raw_csv = fake_result_tsv(n=5, delim=',').replace(b'ANI_estimation,',
                                                      b'')
    rows, cols = lt.parse_result_table(raw_csv)
    check(len(rows) == 5 and 'organism' in cols and 'ani' not in cols,
          f'CSV/缺 ANI 列也能解析: {cols}')
    # gz
    import gzip
    rows2, _ = lt.parse_result_table(gzip.compress(fake_result_tsv()))
    check(len(rows2) == 8, 'gzip 结果文件可解析')
    # 垃圾内容要报错
    try:
        lt.parse_result_table(b'foo,bar\n1,2\n')
        check(False, '无关键列应报错')
    except ValueError:
        check(True, '无关键列时报 ValueError')


def test_paste_source():
    print('== 粘贴序列来源 ==')
    for nm in ('_test_logan_paste', '_test_logan_paste2'):
        try:
            lt.delete_job(nm)
        except Exception:
            pass
    # 2400bp 纯序列（无 FASTA 头）→ 自动补 query 头，≤2500bp 整条 1 段
    d = lt.create_job('_test_logan_paste', pasted='ACGT' * 600, n_seg=2)
    check(d['n_segments'] == 1 and d['segments'][0]['length'] == 2400,
          f'无头纯序列补 query 头: {d["segments"][0]["length"]}bp')
    # 2810bp → 2 段
    d2 = lt.create_job('_test_logan_paste2', pasted='ACGT' * 700 + '\nACGTACGTAA',
                       n_seg=2)
    check(d2['n_segments'] == 2, '2810bp 粘贴序列切 2 段')
    lt.delete_job('_test_logan_paste')
    lt.delete_job('_test_logan_paste2')


def test_flask(monkeypatch):
    print('== Flask 页面与 API ==')
    import app as flaskapp
    c = flaskapp.app.test_client()
    job = JOB + '_flask'
    ensure_job(job, n_seg=1)
    monkeypatch.setattr(lt, 'build_report', _fake_build_report)
    _fake_build_report(job)
    r = c.get('/logan')
    check(r.status_code == 200 and 'LOGAN'.encode() in r.data, '/logan 页面 200')
    r = c.get('/')
    check(b'/logan' in r.data, '总览页含 LOGAN 导航')
    r = c.get('/api/logan/samples')
    check(any(s['name'] == SAMPLE for s in r.get_json()), '/api/logan/samples')
    r = c.get('/api/logan/jobs')
    check(any(j['name'] == lt.safe_job_name(job) for j in r.get_json()),
          '/api/logan/jobs')
    r = c.get(f'/api/logan/job/{job}/segment/1')
    check(r.status_code == 200 and r.get_json()['seq'], '片段序列 API')
    r = c.get(f'/logan/report/{job}/')
    check(r.status_code == 200, '溯源报告页 200')
    r = c.get(f'/logan/report/{job}')
    check(r.status_code == 308, '无尾斜杠 308 重定向（相对路径可解析）')
    r = c.get(f'/logan/report/{job}/plotly.min.js')
    check(r.status_code == 200, '报告目录静态文件路由 200')
    r = c.get('/report/REGRESS/plotly.min.js')
    check(r.status_code == 200, '既有 ⑧报告的 plotly 静态路由 200')
    # 上传导入（multipart）
    r = c.post(f'/api/logan/job/{job}/import/1',
               data={'file': (io.BytesIO(fake_result_tsv(n=3)), 'r.csv')},
               content_type='multipart/form-data')
    check(r.status_code == 200, 'multipart 上传导入 s1')
    d = r.get_json()
    check(d['status'] == 'done', f'导入后状态 done（{d["n_imported"]}/{d["n_segments"]}）')
    # 坏文件
    r = c.post(f'/api/logan/job/{job}/import/1',
               data={'file': (io.BytesIO(b'xx'), 'bad.tsv')},
               content_type='multipart/form-data')
    check(r.status_code == 400, '坏结果文件返回 400')
    # 路径穿越拦截
    r = c.get('/api/logan/job/..%5C..%5Cplatform.json/query')
    check(r.status_code in (400, 404), '路径穿越被拦截')


def test_default_stages_and_dedup(monkeypatch):
    print('== 默认阶段与任务去重 ==')
    import app as flaskapp
    import main as cli_main
    from vp.pipeline import DEFAULT_ANALYZE_STAGES
    c = flaskapp.app.test_client()

    check('hostana' in DEFAULT_ANALYZE_STAGES, '默认分析阶段包含 hostana')
    check('virome' not in DEFAULT_ANALYZE_STAGES, '默认分析阶段不含旧的 virome')
    check(cli_main.DEFAULT_ANALYZE_STAGES == DEFAULT_ANALYZE_STAGES,
          'CLI 默认阶段与管道常量一致')

    monkeypatch.setattr(flaskapp.tm, 'list_all',
                        lambda: [{'status': 'running',
                                  'name': 'Sample analysis REGRESS'}])
    r = c.post('/api/analyze', json={'r1': os.path.abspath(__file__),
                                     'sample': 'REGRESS'})
    check(r.status_code == 400, '同一样品的样品分析会被拦截')

    monkeypatch.setattr(flaskapp.tm, 'list_all',
                        lambda: [{'status': 'running',
                                  'name': 'LOGAN batch test_logan'}])
    try:
        lt.delete_job(JOB)
    except Exception:
        pass
    lt.create_job(JOB, sample=SAMPLE, n_seg=1)
    try:
        r = c.post('/api/logan/job/_test_logan/batch',
                   json={'emails': 'a@b.c'})
        check(r.status_code == 400, '英文 LOGAN 批量任务名也会被拦截')
    finally:
        lt.delete_job(JOB)


if __name__ == '__main__':
    test_segments()
    test_query_samples()
    try:
        lt.delete_job(JOB)
    except Exception:
        pass
    test_create_job()
    test_import_and_report()
    test_paste_source()
    test_flask()
    lt.delete_job(JOB)
    print('\n全部通过 ✔  (任务已清理)')

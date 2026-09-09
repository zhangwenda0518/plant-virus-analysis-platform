# -*- coding: utf-8 -*-
"""平台级自测：所有页面 200 + 关键 API 结构 + 页面渲染含示例按钮/示例文件
+ 既有样品报告与结果中心数据。只读，不启动分析任务。"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vp.config import PLATFORM_ROOT  # noqa: E402

EX = os.path.join(PLATFORM_ROOT, 'databases', 'examples')

PAGES = ['/', '/pipeline', '/samples', '/hostremoval', '/hostpredict',
         '/orf', '/annotation', '/genome', '/primer', '/build', '/results',
         '/tools', '/settings', '/meta', '/download', '/logan', '/submit',
         '/virome', '/tools?g=virus', '/tools?g=annotate', '/tools?g=compare']


def check(cond, msg):
    print(('  ok ' if cond else '  FAIL ') + msg, flush=True)
    assert cond, msg


def main():
    import app as appmod  # noqa: E402
    c = appmod.app.test_client()

    # ---------- 1. 页面 ----------
    print('--- 页面渲染 ---', flush=True)
    for p in PAGES:
        r = c.get(p)
        check(r.status_code == 200, f'GET {p} -> {r.status_code}')

    # 页面内含示例按钮（本次补充的示例入口都要真实渲染出来）
    html = c.get('/tools?g=compare').get_data(as_text=True)
    for marker, where in [('seqPrepExample()', '参考序列获取卡'),
                          ('ictvPreview()', 'ICTV 选参预览'),
                          ('id="spCascade"', 'ICTV 级联下拉容器'),
                          ('tbBuild(this)', '进化树构建·集合建树'),
                          ('id="tbMolecule"', '建树分子类型'),
                          ('alignRun(this)', '序列比对·运行'),
                          ('alignEditToggle()', '比对查看器·编辑模式'),
                          ('alignSave(this)', '比对查看器·保存副本'),
                          ("fillExample('al_fa', EXAMPLE_SET_FASTA)", '比对卡示例'),
                          ('fillExample(\'qt_fa\', EXAMPLE_SET_FASTA)', 'FASTA 建树'),
                          ('fillExample(\'tv_file\', EXAMPLE_TREE_NWK)', '树查看'),
                          ('fillExample(\'ms_fa\', EXAMPLE_SET_FASTA)', 'MSA 查看'),
                          ('fillExample(\'sd_fa\', EXAMPLE_SET_FASTA)', 'SDT 卡'),
                          ('fillExample(\'s_files\', EXAMPLE_SYNTENY_GBS)', '导入示例 .gb')]:
        check(marker in html, f'tools 页渲染含 {where}')
    # 比较组四模块结构：参考序列获取 / 进化树构建 / SDT / 共线性
    for gone in ('id="t-contigs-struct"', 'id="t-ncbi"', 'id="t-synteny-gb"',
                 'id="t-msa"', 'id="t-tree"', 'id="s_style"', 'id="s_lovis4u"'):
        check(gone not in html, f'旧结构已移除: {gone}')
    for present in ('id="t-seqprep"', 'id="t-treebuild"', 'id="t-sdt"',
                    'id="t-synteny"'):
        check(present in html, f'新模块存在: {present}')
    js = open(os.path.join(PLATFORM_ROOT, 'webapp', 'static', 'app.js'),
              encoding='utf-8').read()
    for fn in ('ictvCascadeRefetch', 'ictvPreview', 'ictvDownload', 'loadTbColls',
               'tbBuild', 'tbBuildFor', 'alignRun', 'alignLoad', 'alignRender',
               'alignEditToggle', 'alignSave', 'alignSend', 'tbMolChanged'):
        check(f'function {fn}' in js or f'async function {fn}' in js,
              f'app.js 定义 {fn}')
    # 旧项目归档：样品/集合列表不再出现下划线与归档项
    r = c.get('/api/samples')
    names = [x['name'] for x in r.get_json()]
    check(all(not n.startswith('_') for n in names), '样品列表无下划线项')
    r = c.get('/api/samples/archived')
    check(r.status_code == 200 and len(r.get_json()) > 0,
          f'归档样品列表可用（{len(r.get_json())} 个）')
    r = c.get('/api/gb/collections')
    check(all(not x['name'].startswith('_') for x in r.get_json()),
          '集合列表无下划线项')
    html = c.get('/annotation').get_data(as_text=True)
    check("fillExample('oa_fa')" in html, '注释页渲染含示例按钮')
    html = c.get('/genome').get_data(as_text=True)
    check('EXAMPLE_GENBANK_GB' in html, '图谱页渲染含示例 GenBank 按钮')
    html = c.get('/primer').get_data(as_text=True)
    check('EXAMPLE_CONSERVED_FASTA' in html, '引物页渲染含模式感知示例按钮')
    js = open(os.path.join(PLATFORM_ROOT, 'webapp', 'static', 'app.js'),
              encoding='utf-8').read()
    for const in ('EXAMPLE_FASTA', 'EXAMPLE_SET_FASTA', 'EXAMPLE_TREE_NWK',
                  'EXAMPLE_GENBANK_GB', 'EXAMPLE_SYNTENY_GBS',
                  'EXAMPLE_CONSERVED_FASTA'):
        check(const in js, f'app.js 定义 {const}')

    # 示例文件齐全
    import glob
    need = ['example_viral_contigs.fasta', 'example_virus_set.fasta',
            'example_conserved_set.fasta', 'example_tree.nwk',
            'example_genome.gb'] + \
        [f'example_synteny_{x}.gb' for x in 'ABC']
    for fn in need:
        check(os.path.isfile(os.path.join(EX, fn)), f'示例文件存在: {fn}')

    # ---------- 2. 关键 API ----------
    print('--- 关键 API ---', flush=True)
    for api in ('/api/tools', '/api/dbs', '/api/samples', '/api/queue',
                '/api/ictv/cascade',
                '/api/settings', '/api/tool/runs', '/api/msa/samples',
                '/api/ncbi/collections', '/api/gb/collections',
                '/api/submit/tables', '/api/submit/contig_runs',
                '/api/meta/collections', '/api/logan/jobs',
                '/api/dl/batches', '/api/tool/viral_contigs?run=x'):
        r = c.get(api)
        check(r.status_code in (200, 400), f'GET {api} -> {r.status_code}'
              + ('（空参数 400 属预期）' if r.status_code == 400 else ''))

    # ---------- 3. 既有样品报告 ----------
    print('--- 样品报告 ---', flush=True)
    r = c.get('/api/samples')
    samples = [x['name'] for x in r.get_json()]
    check(len(samples) > 0, f'样品列表 {len(samples)} 个')
    for s in samples:
        r = c.get(f'/report/{s}/')
        if r.status_code == 200:
            check(True, f'样品 {s} 报告可访问')
            break
    else:
        print('  SKIP（无样品带报告）')

    # ---------- 3a. contigs 卡：Krona 前置 / 宿主筛选与自动预测 / ID 复制 / 跳转 ----
    hv = c.get('/tools?g=virus').get_data(as_text=True)
    # Krona 旭日块已并入分类报告折叠区（summary=分类旭日图（可下钻））。
    # 报告区顺序（2026-09-09 用户确认）：分类表 → contig 明细 → 宿主预测统计 → 旭日图
    # 注：宿主预测统计在 contig 明细下方是用户明确要求（DEVELOPMENT_NOTES 五十五）；
    #     旭日图并入折叠区后排在最后，故断言为 det < host < sun。
    _i_host = hv.find('id="hostRptSec"')
    _i_det = hv.find('病毒序列分类（contig 明细）')
    _i_sun = hv.find('id="sunC"')
    check(0 < _i_det < _i_host < _i_sun,
          '报告区顺序：分类表 → contig 明细 → 宿主预测统计 → 旭日图')
    check(hv.count('<section class="card" id="t-contigs">') == 1
          and hv.count('<section class="card" id="t-assemble">') == 1,
          '④ contigs / ③ assemble 卡存在且唯一')
    check(hv.count('id="hostRptSec"') == 1, '宿主预测统计块唯一（无错插副本）')
    for marker, where in [("vEx('tmv')", '示例病毒 TMV'),
                          ("vEx('pstvd')", '示例病毒 PSTVd'),
                          ("vEx('mix')", '示例病毒 Mix All')]:
        check(marker in hv, f'病毒识别组渲染含 {where}')
    for fn in ('example_tmv.fasta', 'example_pvy.fasta', 'example_cmv.fasta',
               'example_pstvd.fasta', 'example_mix.fasta'):
        check(os.path.isfile(os.path.join(EX, fn)), f'Metabuli 示例病毒存在: {fn}')
    check('id="anaRun"' not in hv and 'id="anaRunLabel"' in hv,
          '分类表「选择运行」下拉已去除（自动跟随最新运行标签）')
    for marker, where in [('id="vcHostF"', '宿主筛选下拉'),
                          ('autoHostIfMissing', '宿主预测自动运行（无宿主列时自动触发）'),
                          ('copyContigSeq(', 'ID 点击复制序列'),
                          ('anaJump(', '四件套跳转注释分析'),
                          ('processAnaJump', '跳转自动续接')]:
        check(marker in hv, f'病毒识别组渲染含 {where}')
    ha = c.get('/tools?g=annotate').get_data(as_text=True)
    check(all(f"'{a}'" in ha for a in ('blastn', 'blastx', 'primer')),
          't-hom 卡含 BLASTN / BLASTX / Primer 分析')
    for marker in ('Best E-value', 'Conserved Domains', 'Primer design complete',
                   'ncbiCdd', 'ncbiNuc'):
        check(marker in ha, f'metabuli 风格结果渲染含 {marker}')

    # ---------- 3b. 模块历史运行组件（折叠/衔接/删除） ----------
    import os as _os
    # 工具工作台各组页渲染的是同一份 tools.html（卡片由前端按组显隐），
    # 故历史容器数直接从模板推导，避免硬编码数量随卡片增删而漂移。
    with open(_os.path.join(PLATFORM_ROOT, 'webapp', 'templates',
                            'tools.html'), encoding='utf-8') as _f:
        _n_rh = _f.read().count('class="rh" id="rh-')
    for pg, n in [('/tools?g=sample', _n_rh), ('/tools?g=compare', _n_rh),
                  ('/hostremoval', 1), ('/hostpredict', 1), ('/orf', 1),
                  ('/annotation', 1), ('/genome', 1), ('/primer', 1)]:
        h = c.get(pg).get_data(as_text=True)
        check(h.count('class="rh" id="rh-') == n, f'{pg} 历史容器 {n} 个')
    _os.makedirs(_os.path.join(PLATFORM_ROOT, 'tool_runs', 'it_rh_chk'),
                 exist_ok=True)
    r = c.post('/api/tool/runs/it_rh_chk/delete')
    check(r.status_code == 200
          and not _os.path.isdir(_os.path.join(PLATFORM_ROOT, 'tool_runs',
                                               'it_rh_chk')),
          '运行目录删除 API')
    r = c.post('/api/tool/runs/no_such_run/delete')
    check(r.status_code == 400, '删除不存在运行 400')

    # ---------- 3c. 导航项 ↔ 卡片一一对应 ----------
    # showModule() 只切换 <main> 的直接子元素（:scope > *）。若某个导航项
    # 对应的 id 不是 main 的直接子元素（例如被嵌在别的卡片里），点它就是
    # 空白页——2026-09-09 的 t-kvchain 正是如此（曾是 t-kvsuite 内的 div）。
    from html.parser import HTMLParser as _HP
    from vp.web.pages import NAV_GROUPS as _NAV

    class _MainKids(_HP):
        def __init__(self):
            super().__init__()
            self.depth = 0
            self.main_depth = None
            self.ids = set()

        def handle_starttag(self, tag, attrs):
            d = dict(attrs)
            if tag == 'main':
                self.main_depth = self.depth
            if (self.main_depth is not None
                    and self.depth == self.main_depth + 1 and d.get('id')):
                self.ids.add(d['id'])
            if tag not in ('br', 'img', 'input', 'meta', 'link', 'hr'):
                self.depth += 1

        def handle_endtag(self, tag):
            if tag not in ('br', 'img', 'input', 'meta', 'link', 'hr'):
                self.depth -= 1
            if tag == 'main':
                self.main_depth = None

    _mk = _MainKids()
    _mk.feed(hv)
    _want = {it['id'] for g in _NAV for it in g['items'] if it.get('id')}
    _missing = sorted(_want - _mk.ids)
    check(not _missing,
          f'导航项均有 <main> 直接子卡片（缺失: {_missing or "无"}）')

    # ---------- 4. 静态资源 ----------
    r = c.get('/static/app.js')
    check(r.status_code == 200, 'app.js 静态资源 200')
    r = c.get('/static/i18n.js')
    check(r.status_code == 200, 'i18n.js 静态资源 200')

    print('PLATFORM CHECKS PASSED', flush=True)


if __name__ == '__main__':
    main()

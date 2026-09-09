# -*- coding: utf-8 -*-
"""工具②报告区复选框归位到各自区域标题下方。"""
import io

def patch(path, pairs):
    s = io.open(path, encoding='utf-8').read()
    for old, new in pairs:
        assert s.count(old) == 1, f'{path}: count={s.count(old)} for {old[:70]!r}'
        s = s.replace(old, new, 1)
    io.open(path, 'w', encoding='utf-8').write(s)
    print(path, 'OK')

patch('webapp/templates/tools.html', [
    ("""    <div id="identReport" style="display:none;margin-top:12px">
      <h3 style="font-size:15px">📊 分类报告（Taxonomy Sankey / 分类表 / 旭日图）</h3>
      <label class="ck"><input type="checkbox" id="skTop5I" checked onchange="renderRptI()"> Only top 5 species per genus</label>
      <label class="ck" style="margin-left:14px"><input type="checkbox" id="skUncI" onchange="renderRptI()"> 显示未分类</label>
      <h4 style="font-size:14px;color:#1a5276;margin:16px 0 8px">Taxonomy Sankey</h4>
      <div id="sankeyI" style="width:100%;height:500px"></div>
      <h4 style="font-size:14px;color:#1a5276;margin:16px 0 8px">Classification Table（分类表 · 按层级缩进）</h4>
      <div style="max-height:440px;overflow:auto;border:1px solid #eee;border-radius:4px">
        <table class="tb" style="font-size:12px;width:100%"><thead><tr><th>Rank</th><th>Taxon</th><th>TaxID</th><th>%</th><th>Reads</th></tr></thead><tbody id="ctI"></tbody></table>
      </div>""",
     """    <div id="identReport" style="display:none;margin-top:12px">
      <h3 style="font-size:15px">📊 分类报告（Taxonomy Sankey / 分类表 / 旭日图）</h3>
      <h4 style="font-size:14px;color:#1a5276;margin:16px 0 8px">Taxonomy Sankey</h4>
      <label class="ck"><input type="checkbox" id="skTop5I" checked onchange="renderRptI()"> Only top 5 species per genus</label>
      <div id="sankeyI" style="width:100%;height:500px"></div>
      <h4 style="font-size:14px;color:#1a5276;margin:16px 0 8px">Classification Table（分类表 · 按层级缩进）</h4>
      <label class="ck"><input type="checkbox" id="skUncI" onchange="renderRptI()"> 显示未分类</label>
      <div style="max-height:440px;overflow:auto;border:1px solid #eee;border-radius:4px;margin-top:6px">
        <table class="tb" style="font-size:12px;width:100%"><thead><tr><th>Rank</th><th>Taxon</th><th>TaxID</th><th>%</th><th>Reads</th></tr></thead><tbody id="ctI"></tbody></table>
      </div>"""),
])
print('② checkboxes repositioned OK')

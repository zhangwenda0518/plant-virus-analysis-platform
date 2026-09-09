#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report_html.py — 交互式全表编辑报告 v2.0
========================================

从 unified_metadata.csv + miuvig.tsv + assembly.tsv 生成完整可编辑 HTML。
所有表格均可双击编辑、批量填充、下载导出。

用法:
  python report_html.py \
      --csv submission_metadata/unified_metadata.csv \
      --miuvig submission_metadata/miuvig.tsv \
      --assembly submission_metadata/assembly.tsv \
      --run-title my_project \
      -o submission_metadata/report.html

（移植自 MMPV-RNA virome_submission_pipeline/report_html.py；
 HTML 骨架移至 templates/submission_report.html 由 jinja2 渲染，
 HTML 写出走平台 vp.utils.safe_open，提交清单改为 Sequin/BankIt 流程）
"""

import argparse, os, json, sys
from pathlib import Path
from datetime import datetime
import pandas as pd

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from vp.utils import safe_open
    _TEMPLATE_DIR = Path(__file__).resolve().parent / 'templates'
else:
    from ..utils import safe_open
    _TEMPLATE_DIR = Path(__file__).resolve().parent / 'templates'

from jinja2 import Environment, FileSystemLoader

_ENV = None


def _jinja_env():
    global _ENV
    if _ENV is None:
        # autoescape=False：src_header/topo_body 等本身是预生成的 HTML 片段，
        # 与原始实现（字符串拼接）语义一致
        _ENV = Environment(autoescape=False,
                           loader=FileSystemLoader(str(_TEMPLATE_DIR)))
    return _ENV


CSS = """*{margin:0;padding:0;box-sizing:border-box}
body{font-family:-apple-system,Segoe UI,sans-serif;background:#f0f2f5;color:#1a1a2e}
.header{background:linear-gradient(135deg,#1a1a2e,#16213e);color:#fff;padding:22px 32px}
.header h1{font-size:19px}.header .meta{font-size:11px;opacity:.7;margin-top:3px}
.stats{display:flex;gap:20px;margin-top:12px}.stat{text-align:center}
.stat .num{font-size:26px;font-weight:700}.stat .label{font-size:10px;opacity:.7}
.stat.warn .num{color:#f9a825}.stat.ok .num{color:#4caf50}
.container{max-width:1400px;margin:16px auto;padding:0 16px}
.toolbar{background:#fff;border-radius:10px;padding:10px 16px;margin-bottom:12px;display:flex;gap:8px;align-items:center;box-shadow:0 1px 3px rgba(0,0,0,.05);flex-wrap:wrap}
.toolbar .btn{padding:7px 14px;border-radius:5px;border:none;cursor:pointer;font-size:12px;font-weight:500}
.btn:hover{opacity:.85}.btn.primary{background:#1a73e8;color:#fff}.btn.success{background:#2e7d32;color:#fff}.btn.warn{background:#e67e22;color:#fff}
.section{background:#fff;border-radius:10px;margin-bottom:12px;box-shadow:0 1px 3px rgba(0,0,0,.05);overflow:hidden}
.section-title{background:#f8f9fa;padding:10px 16px;font-weight:600;font-size:13px;border-bottom:1px solid #e0e0e0;display:flex;justify-content:space-between;align-items:center;cursor:pointer}
.section-body{padding:10px 16px;overflow-x:auto;display:block}
.section.collapsed .section-body{display:none}
.section-title::after{content:'▾';font-size:11px;color:#999}
.section.collapsed .section-title::after{content:'▸'}
table{width:100%;border-collapse:collapse;font-size:11px}
th{background:#f5f5f5;padding:6px 5px;text-align:left;border:1px solid #e0e0e0;font-size:10px;color:#666;white-space:nowrap;position:sticky;top:0}
td{padding:4px 5px;border:1px solid #e8e8e8;max-width:180px;overflow:hidden;text-overflow:ellipsis}
td .cell-content{cursor:pointer;padding:1px 0}
td.editing .cell-content{display:none}
td.editing input,td.editing select{width:100%;min-width:60px;border:1px solid #1a73e8;padding:2px 4px;font-size:11px;border-radius:3px;display:none}
td.editing input{display:block}
td.missing{background:#fff3e0}td.missing .cell-content{color:#e65100}
td.filled{background:#e8f5e9}td.filled .cell-content{color:#2e7d32}
.toast{position:fixed;top:20px;right:20px;background:#2e7d32;color:#fff;padding:10px 16px;border-radius:6px;opacity:0;transition:opacity .3s;z-index:9999}
.toast.show{opacity:1}.toast.error{background:#c62828}
.filter input{padding:5px 10px;border:1px solid #ddd;border-radius:5px;font-size:12px;width:180px}
.progress{height:5px;background:#e0e0e0;border-radius:3px;margin-top:8px}
.progress-bar{height:100%;border-radius:3px;transition:width .3s}
.legend{display:flex;gap:14px;font-size:10px;margin:6px 0}
.legend span{display:flex;align-items:center;gap:3px}
.legend .dot{width:9px;height:9px;border-radius:2px}
.dot.missing-dot{background:#fff3e0;border:1px solid #e65100}
.dot.filled-dot{background:#e8f5e9;border:1px solid #2e7d32}
.checklist li{padding:2px 0;font-size:11px}.checklist .done{color:#2e7d32}.checklist .todo{color:#e65100}
.tabs{display:flex;gap:0;margin-bottom:16px}
.tab{padding:8px 16px;background:#e0e0e0;cursor:pointer;font-size:12px;border-radius:8px 8px 0 0;margin-right:2px}
.tab.active{background:#fff;font-weight:600}
"""

JS_EDITABLE = """
const PLACEHOLDER = ['XXXX','YYYY','DD-Mmm','Country:Region','PRJNAXXXX','SAMNXXXX','XX.XX N','Author','not_provided'];

function isMissing(v){
  if(!v||v.trim()==='')return true;
  for(const p of PLACEHOLDER)if(v.includes(p))return true;
  return false;
}

function makeEditableTable(tableId){
  const tbl=document.getElementById(tableId);
  if(!tbl)return;
  const ths=tbl.querySelectorAll('thead th');
  tbl.querySelectorAll('tbody td').forEach(td=>{
    const val=td.textContent.trim();
    const row=td.parentElement.rowIndex-1;
    const col=td.cellIndex;
    td.innerHTML='';
    const span=document.createElement('span');span.className='cell-content';span.textContent=val;
    const inp=document.createElement('input');inp.value=val;
    inp.setAttribute('data-row',row);inp.setAttribute('data-col',col);
    td.appendChild(span);td.appendChild(inp);
    span.ondblclick=()=>{td.classList.add('editing');inp.focus();};
    inp.onblur=()=>{td.classList.remove('editing');span.textContent=inp.value;
      td.classList.toggle('missing',isMissing(inp.value));
      td.classList.toggle('filled',!isMissing(inp.value)&&inp.value.trim()!=='');
      updateStats();
    };
    inp.onkeydown=e=>{if(e.key==='Enter')inp.blur();if(e.key==='Escape'){inp.value=span.textContent;inp.blur();}};
    td.classList.toggle('missing',isMissing(val));
    td.classList.toggle('filled',!isMissing(val)&&val.trim()!=='');
  });
}

function collectTable(tableId){
  const tbl=document.getElementById(tableId);
  if(!tbl)return{headers:[],rows:[]};
  const headers=[];tbl.querySelectorAll('thead th').forEach(th=>headers.push(th.textContent.trim()));
  const rows=[];
  tbl.querySelectorAll('tbody tr').forEach(tr=>{
    const cols=[];tr.querySelectorAll('td').forEach(td=>{
      const inp=td.querySelector('input');
      cols.push(inp?inp.value:(td.querySelector('.cell-content')?.textContent||td.textContent||''));
    });rows.push(cols);
  });
  return{headers,rows};
}

function downloadTSV(tableId,filename){
  const{headers,rows}=collectTable(tableId);
  const lines=[headers.join('\\t')];
  rows.forEach(r=>lines.push(r.join('\\t')));
  const blob=new Blob([lines.join('\\n')],{type:'text/tab-separated-values'});
  const a=document.createElement('a');a.href=URL.createObjectURL(blob);
  a.download=filename;a.click();
  showToast(filename+' downloaded!');
}

function downloadCSV(tableId,filename){
  const{headers,rows}=collectTable(tableId);
  const lines=[headers.join(',')];
  rows.forEach(r=>lines.push(r.map(c=>'"'+c.replace(/"/g,'""')+'"').join(',')));
  const blob=new Blob([lines.join('\\n')],{type:'text/csv'});
  const a=document.createElement('a');a.href=URL.createObjectURL(blob);
  a.download=filename;a.click();
  showToast(filename+' saved!');
}

function batchFill(tableId,colName,val){
  if(!val)val=prompt('输入 '+colName+' 的值:');
  if(!val)return;
  const tbl=document.getElementById(tableId);
  const headers=[];tbl.querySelectorAll('thead th').forEach(th=>headers.push(th.textContent.trim()));
  const ci=headers.indexOf(colName);
  if(ci<0){showToast('column not found: '+colName,true);return;}
  tbl.querySelectorAll('tbody tr').forEach(tr=>{
    const td=tr.querySelectorAll('td')[ci];
    if(!td)return;
    const inp=td.querySelector('input');
    const span=td.querySelector('.cell-content');
    if(inp)inp.value=val;
    if(span)span.textContent=val;
    td.classList.remove('missing');td.classList.add('filled');
  });
  updateStats();showToast('Filled '+colName+' = '+val);
}

function updateStats(){
  let total=0,missing=0;
  document.querySelectorAll('#srcTable tbody td').forEach(td=>{
    const inp=td.querySelector('input');
    const v=inp?inp.value:(td.querySelector('.cell-content')?.textContent||'');
    if(td.classList.contains('required-col')){total++;if(isMissing(v))missing++;}
  });
  const pct=total>0?Math.round((total-missing)/total*100):100;
  document.getElementById('completeness').textContent=pct+'%';
  document.getElementById('completeBar').style.width=pct+'%';
  document.getElementById('completeBar').style.background=pct>=90?'#4caf50':pct>=70?'#f9a825':'#e65100';
  document.getElementById('missingCount').textContent=missing;
}

function showToast(msg,isErr){
  const t=document.createElement('div');t.className='toast'+(isErr?' error':'');
  t.textContent=msg;document.body.appendChild(t);
  setTimeout(()=>t.classList.add('show'),50);
  setTimeout(()=>{t.classList.remove('show');setTimeout(()=>t.remove(),300);},2200);
}

function filterTable(q){
  q=q.toLowerCase();
  document.querySelectorAll('#srcTable tbody tr').forEach(r=>{
    r.style.display=!q||r.textContent.toLowerCase().includes(q)?'':'none';
  });
}

function toggleSection(id){
  document.getElementById(id).classList.toggle('collapsed');
}

window.addEventListener('DOMContentLoaded',()=>{
  makeEditableTable('srcTable');
  makeEditableTable('cmtTable');
  makeEditableTable('bsTable');
  makeEditableTable('miuvigTable');
  makeEditableTable('asmTable');
  makeEditableTable('topoTable');
  updateStats();
});
"""


def safe_row(row, cols, max_len=80):
    """Generate a row, truncating long values"""
    cells = []
    for c in cols:
        v = str(row.get(c, '')) if pd.notna(row.get(c, '')) else ''
        if len(v) > max_len:
            v = v[:max_len-3] + '...'
        cells.append(v)
    return cells


def table_html(headers, rows):
    th = ''.join('<th>' + str(h) + '</th>' for h in headers)
    body = ''
    for row in rows:
        body += '<tr>' + ''.join('<td>' + str(c) + '</td>' for c in row) + '</tr>\n'
    return th, body


def load_table(tsv_path):
    if not tsv_path or not os.path.exists(tsv_path):
        return [], []
    df = pd.read_csv(tsv_path, sep='\t' if str(tsv_path).endswith('.tsv') else None,
                     engine='python')
    return list(df.columns), [list(df.iloc[i]) for i in range(len(df))]


def generate_html(csv_path, run_title, output_path, miuvig_path=None, asm_path=None,
                   suvtk_tax_path=None, suvtk_feat_path=None, hypo_path=None,
                   topo_path=None, gb_dir=None, log=None):
    """Generate the full interactive report"""

    df = pd.read_csv(csv_path, encoding='utf-8-sig')
    n_seqs = len(df)
    n_viruses = df['organism'].nunique() if 'organism' in df.columns else 0

    # ── 汇总统计（taxonomy/features/topology 文件可选） ──
    summary_parts = []
    if suvtk_tax_path and os.path.isfile(suvtk_tax_path):
        with open(suvtk_tax_path) as f:
            n_tax = sum(1 for _ in f) - 1
        summary_parts.append(f'taxonomy: {n_tax} sequences classified')
    if suvtk_feat_path and os.path.isfile(suvtk_feat_path):
        n_cds = sum(1 for l in open(suvtk_feat_path) if l.strip().split()[-1:] == ['CDS'])
        summary_parts.append(f'features: {n_cds} CDS annotated')
    if hypo_path and os.path.isfile(hypo_path):
        n_hypo = sum(1 for l in open(hypo_path) if 'hypothetical' in l.lower())
        n_anno = sum(1 for l in open(hypo_path) if 'inference\talignment' in l)
        summary_parts.append(f'hypothetical proteins: {n_anno} annotated, {n_hypo} remaining')
    if topo_path and os.path.isfile(topo_path):
        rows = []
        with open(topo_path) as f:
            hdr = f.readline().strip().split('\t')
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= len(hdr):
                    rows.append(dict(zip(hdr, cols)))
        n_linear = sum(1 for r in rows if 'linear' in r.get('final_topology', '').lower())
        n_circular = len(rows) - n_linear
        summary_parts.append(f'topology: {n_linear} linear / {n_circular} circular')
    if gb_dir and os.path.isdir(gb_dir):
        n_gb = len([f for f in os.listdir(gb_dir) if f.endswith('.gb')])
        if n_gb:
            summary_parts.append(f'GenBank files: {n_gb} .gb ready')
    summary_html = ' | '.join(summary_parts) if summary_parts else '<span style="color:#999">topology/features 未提供（可选）</span>'

    # --- source.src table ---
    src_cols = [c for c in df.columns if not c.startswith('cmt-') and not c.startswith('bs-')]
    src_header, src_body_str = table_html(src_cols, [safe_row(df.iloc[i], src_cols) for i in range(n_seqs)])

    # --- cmt table ---
    cmt_cols = [c for c in df.columns if c.startswith('cmt-')]
    if cmt_cols:
        cmt_header, cmt_body_str = table_html(cmt_cols, [safe_row(df.iloc[i], cmt_cols) for i in range(n_seqs)])
    else:
        cmt_header, cmt_body_str = '<th>Parameter</th><th>Value</th>', '<tr><td>Assembly_Method</td><td>SPAdes;4.3.0</td></tr>'

    # --- bs table ---
    bs_cols = [c for c in df.columns if c.startswith('bs-')]
    if bs_cols:
        bs_header, bs_body_str = table_html(bs_cols, [safe_row(df.iloc[i], bs_cols) for i in range(n_seqs)])
    else:
        bs_header, bs_body_str = '<th>Parameter</th><th>Value</th>', ''

    # --- miuvig table ---
    miuvig_headers, miuvig_rows = load_table(miuvig_path)
    miuvig_header = '<th>Parameter</th><th>Value</th>'
    miuvig_body_str = '\n'.join(f'<tr><td>{r[0]}</td><td>{r[1] if len(r)>1 else ""}</td></tr>' for r in miuvig_rows)

    # --- assembly table ---
    asm_headers, asm_rows = load_table(asm_path)
    asm_header = '<th>Parameter</th><th>Value</th>'
    asm_body_str = '\n'.join(f'<tr><td>{r[0]}</td><td>{r[1] if len(r)>1 else ""}</td></tr>' for r in asm_rows)

    # --- topology table ---
    topo_header = '<th>Contig</th><th>Taxonomy</th><th>Rule</th><th>Final</th><th>Length</th><th>DTR</th><th>Evidence</th>'
    topo_body_str = '<tr><td colspan="7" style="color:#999;text-align:center">No topology data</td></tr>'
    topo_summary = 'No data'
    if topo_path and os.path.isfile(topo_path):
        topo_rows = []
        with open(topo_path) as f:
            topo_hdr = f.readline().strip().split('\t')
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= len(topo_hdr):
                    topo_rows.append(dict(zip(topo_hdr, cols)))
        if topo_rows:
            n_linear = sum(1 for r in topo_rows if 'linear' in r.get('final_topology','').lower())
            n_circular = sum(1 for r in topo_rows if 'circular' in r.get('final_topology','').lower())
            topo_summary = f'{len(topo_rows)} seqs: {n_linear} linear / {n_circular} circular'
            topo_body_str = ''
            for row in topo_rows:
                topo = row.get('final_topology','').lower()
                cls = 'filled' if 'linear' in topo else ('missing' if 'circular' not in topo else '')
                dtr_flag = '&#10003;' if row.get('dtr_detected','False').lower() == 'true' else ''
                topo_body_str += (
                    f'<tr>'
                    f'<td>{row.get("contig","")}</td>'
                    f'<td>{row.get("taxonomy","")}</td>'
                    f'<td>{row.get("rule_topology","")}</td>'
                    f'<td class="{cls}">{row.get("final_topology","")}</td>'
                    f'<td style="text-align:right">{row.get("seq_length","")}</td>'
                    f'<td style="text-align:center">{dtr_flag}</td>'
                    f'<td style="font-size:10px;color:#888">{row.get("evidence","")}</td>'
                    f'</tr>\n')

    # --- checklist（平台流程：⑥b 注释 → 拓扑 → Sequin/BankIt 提交） ---
    chk = [
        ('done', '病毒 contigs 分类（③ 组装 + kunpeng/BLAST）'),
        ('done', 'ORF 预测与功能注释（⑥/⑥b orfipy + pyrodigal + DIAMOND/MMseqs2）'),
        ('done', '特征表 .tbl + 拓扑判定（taxonomy 规则 + DTR/ITR 检测）'),
        ('todo', '核对/补全上方表格 → 下载 source.src（或保存 CSV 备份）'),
        ('todo', '在 https://submit.ncbi.nlm.nih.gov/ 注册 BioProject'),
        ('todo', '用 biosample_template.tsv 注册 BioSample'),
        ('todo', 'NCBI BankIt 上传 .fsa + .tbl（或用 Sequin 桌面程序导入后生成 .sqn）'),
        ('todo', '提交后留意 NCBI 校对邮件（gb-sub@ncbi.nlm.nih.gov）'),
    ]
    check_mark = '\u2713'; box = '\u2610'
    chk_html = '\n'.join(f'<li class="{c}">{check_mark if c=="done" else box} {t}</li>' for c, t in chk)

    html = _jinja_env().get_template('submission_report.html').render(
        run_title=run_title,
        timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        pipeline='植物病毒分析平台 ③组装+⑥b ORF 注释 → Sequin 提交包',
        n_seqs=n_seqs, n_viruses=n_viruses,
        src_header=src_header, src_body=src_body_str,
        cmt_header=cmt_header, cmt_body=cmt_body_str,
        bs_header=bs_header, bs_body=bs_body_str,
        miuvig_header=miuvig_header, miuvig_body=miuvig_body_str,
        asm_header=asm_header, asm_body=asm_body_str,
        topo_header=topo_header, topo_body=topo_body_str,
        topo_summary=topo_summary,
        summary=summary_html,
        checklist=chk_html,
        css=CSS, js=JS_EDITABLE,
    )

    with safe_open(output_path, 'wt') as f:
        f.write(html)

    if log:
        log.info("  → %s (%d seqs)", output_path, n_seqs)
    print(f"\n  Report: {output_path} | {n_seqs} seqs, {n_viruses} viruses")
    print(f"  Double-click any cell to edit | All tables downloadable")

    return output_path


def main():
    parser = argparse.ArgumentParser(description="report_html.py v2.1 — Full interactive submission report")
    parser.add_argument('--csv', required=True, help='unified_metadata.csv')
    parser.add_argument('--miuvig', help='miuvig.tsv')
    parser.add_argument('--assembly', help='assembly.tsv')
    parser.add_argument('--suvtk-tax', help='taxonomy.tsv（分类表）')
    parser.add_argument('--suvtk-feat', help='featuretable.tbl（特征表）')
    parser.add_argument('--hypo', help='hypothetical 注释表（可选）')
    parser.add_argument('--topology', help='topology.tsv')
    parser.add_argument('--gb-dir', help='GenBank 文件目录（可选）')
    parser.add_argument('--run-title', default='viral_submission')
    parser.add_argument('-o', '--output', required=True, help='输出 HTML（限平台目录内）')
    args = parser.parse_args()
    generate_html(args.csv, args.run_title, args.output,
                  args.miuvig, args.assembly,
                  args.suvtk_tax, args.suvtk_feat, args.hypo,
                  args.topology, args.gb_dir)


if __name__ == '__main__':
    main()

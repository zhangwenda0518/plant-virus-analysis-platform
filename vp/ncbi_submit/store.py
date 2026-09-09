#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
store.py — 统一提交元数据表（unified_metadata.csv）服务端存储与操作
====================================================================

对应 MMPV-RNA submission_gui.py 中 SubmissionStore 的功能面（Web/CLI 化）：

  - 加载/保存 unified_metadata.csv（CSV/TSV/Excel，utf-8-sig）
  - 占位符/缺失检测（is_placeholder：NA 词表 + XXXX/SAMNXX/PRJNAXX… 模板残留）
  - 批量替换 / 批量填充占位符 / 增删行
  - 必填字段校验（REQUIRED_COLS）
  - NCBI BioSample 注册 TSV 导出（跳过占位符行）
  - template.sbt 生成（ASN.1 Submit-block，括号配平校验）

提交项目存放在 DIRS['submissions']/<name>/ 下，unified_metadata.csv 为中枢表；
source.src / miuvig.tsv / assembly.tsv / report.html 等产物由 export() 用
unified_metadata.py + report_html.py 的函数生成。
"""

import os
import re
import shutil
import sys
import csv as _csv
from pathlib import Path

import pandas as pd

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from vp.config import DIRS
    from vp.utils import safe_open, check_path
    from vp.ncbi_submit.unified_metadata import UNIFIED_COLUMNS, REQUIRED_COLS
else:
    from ..config import DIRS
    from ..utils import safe_open, check_path
    from .unified_metadata import UNIFIED_COLUMNS, REQUIRED_COLS

# 列描述（unified_metadata.UNIFIED_COLUMNS 的 desc 字段）
COL_DESC = {name: spec.get('desc', '') for name, spec in UNIFIED_COLUMNS.items()}

NA_PLACEHOLDERS = {"NA", "N/A", "Not_Provided", "not collected",
                   "missing", "none", "unknown", "", " "}
PLACEHOLDER_RE = re.compile(
    r'XXXX|YYYY|PRJNAXXXX|Country:Region|SAMNXXXXXXXX|XX\.\d+', re.IGNORECASE)


def is_placeholder(val) -> bool:
    """空值 / NA 占位词 / 模板残留（XXXX、SAMNXX…）都算未填。"""
    if val is None or (not isinstance(val, str) and pd.isna(val)):
        return True
    s = str(val).strip()
    if s.lower() in NA_PLACEHOLDERS or s == "":
        return True
    return bool(PLACEHOLDER_RE.search(s))


def _root():
    return DIRS['submissions']


def _table_dir(name):
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name or ''):
        raise ValueError(f"非法提交项目名: {name!r}（仅限字母/数字/_/-）")
    return check_path(os.path.join(_root(), name), must_exist=False,
                      in_platform=True)


def _csv_path(name):
    return os.path.join(_table_dir(name), 'unified_metadata.csv')


# ══════════════════════════════════════════════════════════════
# 项目管理
# ══════════════════════════════════════════════════════════════

def list_tables():
    """枚举提交项目：[{name, rows, files, mtime}]"""
    out = []
    root = _root()
    if not os.path.isdir(root):
        return out
    for name in sorted(os.listdir(root)):
        d = os.path.join(root, name)
        csv_p = os.path.join(d, 'unified_metadata.csv')
        if not os.path.isfile(csv_p):
            continue
        try:
            df = load_df(csv_p)
        except Exception:
            continue
        try:
            mt = os.path.getmtime(csv_p)
        except OSError:
            mt = 0
        out.append({'name': name, 'rows': int(len(df)),
                    'mtime': mt,
                    'files': sorted(os.listdir(d)) if os.path.isdir(d) else []})
    return out


def table_dir(name):
    return _table_dir(name)


def load_df(path):
    """读 CSV/TSV/Excel → DataFrame（dtype=str，保留空串）。"""
    path = str(path)
    if path.lower().endswith(('.xlsx', '.xls')):
        df = pd.read_excel(path, dtype=str, keep_default_na=False)
    elif path.lower().endswith('.tsv') or path.lower().endswith('.txt'):
        df = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False,
                         encoding='utf-8-sig')
    else:
        df = pd.read_csv(path, dtype=str, keep_default_na=False,
                         encoding='utf-8-sig')
    df.fillna('', inplace=True)
    return df.astype(object)


def create_table(name, df=None, sample=''):
    """新建提交项目目录 + unified_metadata.csv。

    sample: ''(空表) / 'demo'(内置示例) / 'public'(含 SRR 公共数据样例) /
    'selfseq'(自测数据样例，占位符待填)——后两者对应原 GUI Samples 菜单。
    """
    d = _table_dir(name)
    os.makedirs(d, exist_ok=True)
    csv_p = _csv_path(name)
    if os.path.isfile(csv_p):
        raise FileExistsError(f"提交项目已存在: {name}")
    if sample:
        src = os.path.join(_sample_dir(), f'sample_{sample}.csv')
        if sample == 'demo' or not os.path.isfile(src):
            cols = list(UNIFIED_COLUMNS.keys())
            df = pd.DataFrame(_DEMO_ROWS if sample == 'demo' else [],
                              columns=cols).astype(object)
        else:
            df = _align_columns(load_df(src))
    elif df is None:
        df = pd.DataFrame(columns=list(UNIFIED_COLUMNS.keys())).astype(object)
    save_df(df, csv_p)
    return csv_p


def _sample_dir():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'samples')


def list_samples():
    """可用的示例数据集（原 GUI Samples 菜单的 Web 版）。"""
    out = []
    d = _sample_dir()
    for key, label in (('public', '公共数据样例（含 SRR/BioSample 登录号）'),
                       ('selfseq', '自测数据样例（占位符待填）')):
        if os.path.isfile(os.path.join(d, f'sample_{key}.csv')):
            out.append({'key': key, 'label': label})
    return out


def copy_table(src, dst):
    """项目另存为（原 GUI Save As）：复制目录内全部文件到新项目。"""
    s = _table_dir(src)
    if not os.path.isfile(_csv_path(src)):
        raise FileNotFoundError(f"源项目不存在: {src}")
    d = _table_dir(dst)
    if os.path.isfile(_csv_path(dst)):
        raise FileExistsError(f"提交项目已存在: {dst}")
    shutil.copytree(s, d,
                    ignore=shutil.ignore_patterns('__pycache__'))
    return _csv_path(dst)


def import_table(name, src_path):
    """从外部 CSV/TSV/Excel 导入为提交项目。"""
    d = _table_dir(name)
    os.makedirs(d, exist_ok=True)
    csv_p = _csv_path(name)
    if os.path.isfile(csv_p):
        raise FileExistsError(f"提交项目已存在: {name}")
    df = load_df(src_path)
    df = _align_columns(df)
    save_df(df, csv_p)
    return csv_p, len(df)


def _align_columns(df):
    """列对齐到 UNIFIED_COLUMNS：缺的补空列，多的保留在后。"""
    cols = list(UNIFIED_COLUMNS.keys())
    for c in cols:
        if c not in df.columns:
            df[c] = ''
    extra = [c for c in df.columns if c not in cols]
    return df[cols + extra]


def save_df(df, path):
    """DataFrame → CSV（utf-8-sig，Excel 直接打开不乱码）。"""
    p = check_path(path, must_exist=False, in_platform=True)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with safe_open(p, 'wt') as f:
        f.write('\ufeff')                      # BOM（safe_open 固定 utf-8）
        df.to_csv(f, index=False, sep=',', lineterminator='\n')
    return p


def export_excel(name, df=None):
    """当前表 → unified_metadata.xlsx（openpyxl）。"""
    d = _table_dir(name)
    df = load_table(name) if df is None else df
    xlsx = os.path.join(d, 'unified_metadata.xlsx')
    df.to_excel(check_path(xlsx, must_exist=False, in_platform=True),
                index=False, engine='openpyxl')
    return xlsx


# ══════════════════════════════════════════════════════════════
# 表读写与编辑操作
# ══════════════════════════════════════════════════════════════

def load_table(name):
    csv_p = _csv_path(name)
    if not os.path.isfile(csv_p):
        raise FileNotFoundError(f"找不到 {csv_p}")
    return load_df(csv_p)


def save_table(name, rows):
    """前端编辑结果整表回写。rows: [{col: val, ...}, ...]"""
    df = load_table(name)
    cols = list(df.columns)
    out = pd.DataFrame(rows)
    for c in cols:
        if c not in out.columns:
            out[c] = ''
    out = out[cols + [c for c in out.columns if c not in cols]]
    out = out.fillna('').astype(object)
    save_df(out, _csv_path(name))
    return len(out)


def delete_rows(name, row_indices):
    df = load_table(name)
    idx = sorted({int(i) for i in row_indices}, reverse=True)
    for i in idx:
        if 0 <= i < len(df):
            df = df.drop(df.index[i])
    df = df.reset_index(drop=True)
    save_df(df, _csv_path(name))
    return len(df)


def add_rows(name, n=1):
    df = load_table(name)
    for _ in range(max(1, int(n))):
        df.loc[len(df)] = {c: '' for c in df.columns}
    save_df(df, _csv_path(name))
    return len(df)


def batch_fill(name, column, new_value, old_value=None):
    """旧值→新值替换（old_value 给定时）或填充该列所有占位符。返回修改数。"""
    df = load_table(name)
    if column not in df.columns:
        raise KeyError(f"列不存在: {column}")
    if old_value:
        mask = df[column].astype(str).str.strip() == str(old_value).strip()
    else:
        mask = df[column].apply(is_placeholder)
    count = int(mask.sum())
    if count:
        df.loc[mask, column] = str(new_value)
        save_df(df, _csv_path(name))
    return count


def table_payload(name):
    """前端表格渲染所需全部数据。"""
    df = load_table(name)
    cols = list(df.columns)
    col_meta = [{'name': c, 'required': c in REQUIRED_COLS,
                 'desc': COL_DESC.get(c, ''),
                 'group': UNIFIED_COLUMNS.get(c, {}).get('group', '')}
                for c in cols]
    rows = []
    for _, r in df.iterrows():
        rows.append([('' if pd.isna(r[c]) else str(r[c])) for c in cols])
    total_cells = len(df) * len(cols)
    filled = sum(1 for row in rows for v in row if not is_placeholder(v))
    return {'name': name, 'columns': col_meta, 'rows': rows,
            'stats': {'rows': len(df), 'cells': total_cells,
                      'filled': filled,
                      'pct': round(filled / total_cells * 100, 1) if total_cells else 0.0}}


def validate_table(name):
    """提交前校验：必填占位符 + 格式 + 唯一性（同 GUI Validate 的增强版）。

    issues 项：{column, kind, count, desc, examples}
      kind: missing_col(缺列) / placeholder(必填未填) / format(格式) /
            duplicate(重复) / blank(选填留空提示)
    """
    df = load_table(name)
    issues = []

    # ── 1. 必填列存在性 + 占位符 ──
    for col in REQUIRED_COLS:
        if col not in df.columns:
            issues.append({'column': col, 'kind': 'missing_col', 'count': -1,
                           'desc': COL_DESC.get(col, ''), 'examples': []})
            continue
        mask = df[col].apply(is_placeholder)
        n = int(mask.sum())
        if n > 0:
            examples = df.loc[mask, col].unique()[:3].tolist()
            issues.append({'column': col, 'kind': 'placeholder', 'count': n,
                           'desc': COL_DESC.get(col, ''), 'examples': examples})

    def _fmt_bad(col, pattern):
        """列内非空且不完全匹配 pattern（整串）的行索引集合。"""
        s = df[col].astype(str).str.strip()
        return s.index[s.str.len().gt(0) & ~s.str.fullmatch(pattern)]

    # ── 2. 字段格式（NCBI 拒收的常见错误）──
    _FMT = [
        # (列, 正则, 说明)
        ('collection_date', r'\d{4}(-\d{2}(-\d{2})?)?|[Nn]ot[_ ]collected|[Nn]ot[_ ]provided',
         '日期需 YYYY / YYYY-MM / YYYY-MM-DD（或 not collected）'),
        ('src-geo_loc_name', r'[A-Za-z][A-Za-z \-]+:.+',
         '需 Country:Region 格式（如 China:Ningxia，可带市级段）'),
        ('src-Lat_Lon', r'-?\d{1,3}(\.\d+)?\s*[NS]\s+-?\d{1,3}(\.\d+)?\s*[EW]',
         '需 "38.47 N 106.27 E" 格式（度 + 半球字母）'),
        ('bioproject', r'PRJ[NED][A-Z]\d+',
         'BioProject 登录号需 PRJNA/PRJNE/PRJND + 数字'),
        ('biosample', r'SAM[NED]A?\d+',
         'BioSample 登录号需 SAMN/SAME/SAMD + 数字'),
        ('sra', r'[SEDC]RR\d+',
         'SRA/CNSA 登录号需 SRR/ERR/DRR/CRR + 数字'),
    ]
    for col, pat, desc in _FMT:
        if col not in df.columns:
            continue
        bad_idx = _fmt_bad(col, pat)
        if len(bad_idx):
            issues.append({'column': col, 'kind': 'format',
                           'count': int(len(bad_idx)), 'desc': desc,
                           'examples': df.loc[bad_idx, col].unique()[:3].tolist()})

    # ── 3. sequence_name：FASTA SeqID 规范（无空白）+ 全表唯一 ──
    if 'sequence_name' in df.columns:
        s = df['sequence_name'].astype(str).str.strip()
        ws = s[s.str.contains(r'\s', na=False) & s.str.len().gt(0)]
        if len(ws):
            issues.append({'column': 'sequence_name', 'kind': 'format',
                           'count': int(len(ws)),
                           'desc': '序列名不能含空格（必须与 FASTA 头一致，NCBI SeqID 规范）',
                           'examples': ws.unique()[:3].tolist()})
        filled = s[s.str.len().gt(0) & ~s.apply(is_placeholder)]
        dup = filled[filled.duplicated(keep=False)].unique().tolist()
        if dup:
            issues.append({'column': 'sequence_name', 'kind': 'duplicate',
                           'count': len(dup),
                           'desc': 'sequence_name 重复（每行必须唯一）',
                           'examples': dup[:3]})
    return issues


# ══════════════════════════════════════════════════════════════
# BioSample 注册 TSV 导出（同 GUI 的 BioSample Export 标签）
# ══════════════════════════════════════════════════════════════

_BS_ATTR_MAP = [
    # (BioSample 列, 统一表主列 bs-*, 回退列 src-*)
    ('collection_date', None, 'collection_date'),
    ('geo_loc_name', 'bs-geo_loc_name', 'src-geo_loc_name'),
    ('host', 'bs-host', 'src-Host'),
    ('isolate', 'bs-isolate', 'src-Isolate'),
    ('isolation_source', 'bs-isolation_source', 'src-Isolation-source'),
    ('tissue_type', None, 'src-Tissue_type'),
    ('collected_by', None, 'src-Collected_by'),
    ('bioproject', None, 'bioproject'),
]
_BS_PLACEHOLDER = re.compile(r"XXXX|SAMNXX|PRJNAXX|^$", re.IGNORECASE)


def export_biosample_tsv(name, skip_placeholders=True):
    """生成 NCBI BioSample 批量注册 TSV。返回 (路径, 行数, 跳过数)。

    字段来源统一为：bs-* 列优先，缺省回退 src-*（与 export_files 生成的
    biosample_template.tsv 同口径，避免两份模板字段来源不一致）。
    """
    df = load_table(name)
    if 'sequence_name' not in df.columns:
        raise KeyError("缺少 sequence_name 列")
    attr_cols = [(a, bs, src) for a, bs, src in _BS_ATTR_MAP
                 if (bs and bs in df.columns) or (src and src in df.columns)]
    header = ['sample_name', 'organism'] + [a for a, _, _ in attr_cols]

    rows, skipped = [], 0
    for _, r in df.iterrows():
        sample = str(r.get('sequence_name', '') or '').strip()
        org = str(r.get('organism', '') or '').strip()
        if skip_placeholders and (_BS_PLACEHOLDER.search(sample)
                                  or _BS_PLACEHOLDER.search(org)):
            skipped += 1
            continue
        vals = [sample, org]
        for _a, bs, src in attr_cols:
            v = ''
            for c in (bs, src):
                if c and str(r.get(c, '') or '').strip():
                    v = str(r[c]).strip()
                    break
            vals.append('' if _BS_PLACEHOLDER.match(v) else v)
        rows.append(vals)
    if not rows:
        raise ValueError('所有行都是占位符，无可导出数据')

    out = os.path.join(_table_dir(name), 'biosample_registration.tsv')
    with safe_open(out, 'wt') as f:
        f.write('\t'.join(header) + '\n')
        for vals in rows:
            f.write('\t'.join(vals) + '\n')
    return out, len(rows), skipped


# ══════════════════════════════════════════════════════════════
# template.sbt 生成（同 GUI 的 Generate SBT 标签，ASN.1 括号配平）
# ══════════════════════════════════════════════════════════════

def generate_sbt(fields, extra_authors=None, title='', out_name=None):
    """fields: {last, first, middle, affil, div, city, sub, country,
                street, email, postal}；extra_authors: ['Li, Ming', ...]。
    返回 template.sbt 文本。"""
    for key in ('last', 'first', 'affil', 'city', 'country', 'email'):
        if not (fields.get(key) or '').strip():
            raise ValueError(f"缺少必填字段: {key}")

    def parse_author(ln):
        ln = ln.strip().rstrip(';')
        if not ln:
            return None
        if ',' in ln:
            last, _, rest = ln.partition(',')
            toks = rest.strip().split()
        else:
            toks = ln.split()
            last, toks = toks[0], toks[1:]
        return (last.strip(), toks[0] if toks else '', ' '.join(toks[1:]))

    authors = [(fields['last'].strip(), fields['first'].strip(),
                (fields.get('middle') or '').strip())]
    for ln in (extra_authors or []):
        a = parse_author(ln)
        if a:
            authors.append(a)

    def esc(s):
        return str(s).replace('"', "'")

    def nm_fields(a, d):
        return (f'{d}last "{esc(a[0])}",\n'
                f'{d}first "{esc(a[1])}",\n'
                f'{d}middle "{esc(a[2])}",\n'
                f'{d}initials "",\n'
                f'{d}suffix "",\n'
                f'{d}title ""')

    def entry(a, brace_ind):
        ni = brace_ind + '  '
        return (f'{brace_ind}' + '{\n' + ni + 'name name {\n'
                + nm_fields(a, ni + '  ') + f'\n{ni}' + '}\n' + brace_ind + '}')

    affil_c = ('      affil std {\n' + '\n'.join(
        [f'        affil "{esc(fields["affil"])}",',
         f'        div "{esc(fields.get("div") or "")}",',
         f'        city "{esc(fields["city"])}",',
         f'        sub "{esc(fields.get("sub") or "")}",',
         f'        country "{esc(fields["country"])}",',
         f'        street "{esc(fields.get("street") or "")}",',
         f'        email "{esc(fields["email"])}",',
         f'        postal-code "{esc(fields.get("postal") or "")}"']) + '\n      }')
    affil_p = ('      affil std {\n' + '\n'.join(
        [f'        affil "{esc(fields["affil"])}",',
         f'        div "{esc(fields.get("div") or "")}",',
         f'        city "{esc(fields["city"])}",',
         f'        sub "{esc(fields.get("sub") or "")}",',
         f'        country "{esc(fields["country"])}",',
         f'        street "{esc(fields.get("street") or "")}",',
         f'        postal-code "{esc(fields.get("postal") or "")}"']) + '\n      }')
    title_text = (title or '').strip() or 'Untitled submission'

    template = (
        'Submit-block ::= {\n'
        '  contact {\n'
        '    contact {\n'
        '      name name {\n'
        + nm_fields(authors[0], '        ') + '\n      },\n'
        + affil_c + '\n'
        '    }\n'
        '  },\n'
        '  cit {\n'
        '    authors {\n'
        '      names std {\n'
        + ',\n'.join(entry(a, '        ') for a in authors) + '\n'
        '      },\n'
        + affil_p + '\n'
        '    }\n'
        '  },\n'
        '  subtype new\n'
        '}\n'
        'Seqdesc ::= pub {\n'
        '  pub {\n'
        '    gen {\n'
        '      cit "unpublished",\n'
        '      authors {\n'
        '        names std {\n'
        + ',\n'.join(entry(a, '          ') for a in authors) + '\n'
        '        }\n'
        '      },\n'
        f'      title "{esc(title_text)}"\n'
        '    }\n'
        '  }\n'
        '}\n'
        'Seqdesc ::= user {\n'
        '  type str "Submission",\n'
        '  data {\n'
        '    {\n'
        '      label str "AdditionalComment",\n'
        f'      data str "ALT EMAIL:{esc(fields["email"])}"\n'
        '    }\n'
        '  }\n'
        '}\n'
        'Seqdesc ::= user {\n'
        '  type str "Submission",\n'
        '  data {\n'
        '    {\n'
        '      label str "AdditionalComment",\n'
        f'      data str "Submission Title:{esc(title_text)}"\n'
        '    }\n'
        '  }\n'
        '}\n')
    if template.count('{') != template.count('}'):
        raise RuntimeError('template.sbt 括号不配平（内部错误）')

    if out_name:
        out = os.path.join(_table_dir(out_name), 'template.sbt')
        with safe_open(out, 'wt') as f:
            f.write(template)
        return template, out
    return template, None


# ══════════════════════════════════════════════════════════════
# 产物导出（source.src / biosample / miuvig / assembly / report）
# ══════════════════════════════════════════════════════════════

def export_files(name, assembler='SPAdes;4.3.0;metaviral',
                 sequencer='Illumina NovaSeq 6000',
                 enrichment='rRNA depletion', log=None):
    """由当前 unified_metadata.csv 生成全部提交产物。返回产物路径 dict。"""
    from . import unified_metadata as um
    from . import report_html as rh

    d = _table_dir(name)
    df = load_table(name)
    log = log or (lambda msg: None)

    class _LogAdapter:
        """unified_metadata 的导出函数要求带 .info/.warning 的 logger。"""
        def __init__(self, cb):
            self._cb = cb
        def info(self, msg, *a):
            self._cb(msg % a if a else msg)
        def warning(self, msg, *a):
            self._cb('[WARN] ' + (msg % a if a else msg))

    _log = _LogAdapter(log)

    out = {}
    um.export_source_src(df, d, _log)
    out['source_src'] = os.path.join(d, 'source.src')

    um.export_biosample_csv(df, d, _log)
    out['biosample_template'] = os.path.join(d, 'biosample_template.tsv')

    um.export_miuvig_assembly(d, _log, assembler, sequencer, enrichment)
    out['miuvig'] = os.path.join(d, 'miuvig.tsv')
    out['assembly'] = os.path.join(d, 'assembly.tsv')

    um.validate_csv(df, d, _log)
    out['validation_report'] = os.path.join(d, 'validation_report.txt')

    csv_p = _csv_path(name)
    rh.generate_html(csv_p, name, os.path.join(d, 'report.html'),
                     miuvig_path=out['miuvig'], asm_path=out['assembly'],
                     log=_log)
    out['report'] = os.path.join(d, 'report.html')

    um.create_submission_log(d, name, _log)
    return out


def read_file(name, fname, limit=300000):
    """预览项目目录下的文本文件（提交目录内白名单式读取）。"""
    p = read_file_path(name, fname)
    with open(p, 'r', encoding='utf-8', errors='replace') as f:
        return f.read(limit)


# 二进制 ASN.1 与打包产物不允许在线编辑（同原 GUI：.sqn view only）
_NO_EDIT = ('.sqn', '.zip', '.xlsx', '.fsa')


def write_file(name, fname, content):
    """保存预览编辑（原 GUI Preview 标签 Edit/Save 的 Web 版）。

    unified_metadata.csv 走 save_table 整表回写，不在此入口；
    .sqn/.zip/.xlsx/.fsa 为二进制产物，拒绝编辑。
    """
    fname = str(fname or '')
    if fname.lower().endswith(_NO_EDIT):
        raise ValueError(f"{os.path.basename(fname)} 为二进制/生成产物，不支持在线编辑")
    if os.path.normpath(fname) == 'unified_metadata.csv':
        raise ValueError('unified_metadata.csv 请在表格编辑器中修改并保存')
    d = _table_dir(name)
    fname = str(fname).replace('\\', '/')
    fname = os.path.normpath(fname)
    if fname.startswith('..') or os.path.isabs(fname):
        raise ValueError(f"非法文件名: {fname}")
    p = check_path(os.path.join(d, fname), must_exist=False, in_platform=True)
    with safe_open(p, 'wt') as f:
        f.write(str(content).replace('\r\n', '\n'))
    return p


def read_file_path(name, fname):
    """校验 fname 位于项目目录内，返回绝对路径。"""
    d = _table_dir(name)
    fname = str(fname or '').replace('\\', '/')
    fname = os.path.normpath(fname)
    if fname.startswith('..') or os.path.isabs(fname):
        raise ValueError(f"非法文件名: {fname}")
    p = os.path.join(d, fname)
    return check_path(p, must_exist=True, in_platform=True)


# ══════════════════════════════════════════════════════════════
# 提交序列 FASTA（GenBank 提交三件套之一：sequences.fsa）
# ══════════════════════════════════════════════════════════════

# DNA 提交不接受 U（RNA 表示），连同其它非 IUPAC 字符统一替换为 N
_SEQ_CLEAN_RE = re.compile(r'[^ACGTYSWKMBDHVN-]', re.IGNORECASE)


def export_submission_fasta(name, fasta_path, min_len=200):
    """按表内 sequence_name 从 fasta_path 提取序列 → sequences.fsa。

    NCBI GenBank 提交三件套 = 序列 FASTA + source.src + template.sbt，
    本模块此前缺序列 FASTA 这一件。要求：
      - 表内 sequence_name 与 FASTA 头第一个 token 一一对应（顺序=表顺序）
      - 序列字符规范化为 IUPAC（其余字符替换 N）
    返回报告 dict {path, written, missing, extra, short, duplicates, replaced}。
    """
    from collections import Counter, OrderedDict
    from ..utils import iter_fasta

    df = load_table(name)
    if 'sequence_name' not in df.columns:
        raise KeyError('缺少 sequence_name 列')
    names = [str(v or '').strip() for v in df['sequence_name'].tolist()]
    names = [n for n in names if n and not is_placeholder(n)]
    if not names:
        raise ValueError('表内没有可用的 sequence_name（先填写或导入行）')
    duplicates = sorted(n for n, c in Counter(names).items() if c > 1)
    want = list(dict.fromkeys(names))          # 去重保序
    want_set = set(want)

    fa_ids, recs = [], {}
    replaced = 0
    for h, s in iter_fasta(check_path(fasta_path, must_exist=True)):
        sid = h.split()[0]
        fa_ids.append(sid)
        if sid in want_set and sid not in recs:
            clean = _SEQ_CLEAN_RE.sub('N', str(s).upper())
            replaced += sum(1 for a, b in zip(str(s).upper(), clean) if a != b)
            recs[sid] = clean

    missing = [n for n in want if n not in recs]
    extra = [sid for sid in fa_ids if sid not in want_set]
    short = [{'id': sid, 'length': len(s)}
             for sid, s in ((n, recs[n]) for n in want if n in recs)
             if len(s) < int(min_len)]

    out = os.path.join(_table_dir(name), 'sequences.fsa')
    with safe_open(out, 'wt') as f:
        for sid in want:
            if sid not in recs:
                continue
            s = recs[sid]
            f.write(f'>{sid}\n')
            for i in range(0, len(s), 70):
                f.write(s[i:i + 70] + '\n')
    return {'path': out, 'written': len(recs), 'missing': missing,
            'extra': extra, 'short': short, 'duplicates': duplicates,
            'replaced': replaced}


def fasta_id_check(name, fasta_path):
    """只做一致性检查（不生成文件）：供上传前预检。返回同 export 报告（path=None）。"""
    r = export_submission_fasta(name, fasta_path)
    r['path'] = None
    # export 已写文件；预检语义下删除生成物，保持"只检查"
    try:
        os.remove(os.path.join(_table_dir(name), 'sequences.fsa'))
    except OSError:
        pass
    return r


def infer_source_fasta(name):
    """自动推断提交序列 FASTA 路径（供前端预填，找不到返回 None）。

    优先级：
      1. 关联的 contigs 运行 viral_contigs.fasta（分类运行输出=病毒序列）
      2. 关联的 orf 运行 04_orf 上游 contigs（03_assembly/*.fasta 或输入）
      3. 关联 contigs 运行目录下其余 fasta（contigs.filtered.fasta）
    返回绝对路径或 None。
    """
    from ..config import PLATFORM_ROOT
    root = DIRS.get('tool_runs') or os.path.join(PLATFORM_ROOT, 'tool_runs')
    cands = []
    cl = _read_link(name, 'contigs')
    if cl and cl.get('run'):
        rd = os.path.join(root, cl['run'])
        for p in (os.path.join(rd, 'viral_contigs.fasta'),
                  os.path.join(rd, 'contigs.filtered.fasta')):
            if os.path.isfile(p):
                cands.append(p)
    if not cands and cl and cl.get('run'):
        for cur, _s, fns in os.walk(os.path.join(root, cl['run'])):
            for fn in sorted(fns):
                if fn.endswith('.fasta') or fn.endswith('.fna'):
                    cands.append(os.path.join(cur, fn))
            if cands:
                break
    ol = _read_link(name, 'orf')
    if ol and ol.get('run') and not cands:
        for cur, _s, fns in os.walk(os.path.join(root, ol['run'])):
            for fn in sorted(fns):
                if fn.endswith('.fasta') or fn.endswith('.fna'):
                    cands.append(os.path.join(cur, fn))
            if cands:
                break
    # 取第一个真实存在的
    for p in cands:
        if os.path.isfile(p):
            return p
    return None


# ══════════════════════════════════════════════════════════════
# 从平台分析结果导入行（打通"识别/再鉴定 → 提交"工作流）
# ══════════════════════════════════════════════════════════════

def _clean_species_name(val):
    """把分类表的物种字段规范成 GenBank organism 可用名。

    处理两类常见噪音：
      - 属名冗余重复（"Cucumovirus CMV" → 保留后者需人工确认，此处原样返回但去多余空白）
      - ICTV 样 "Cucumovirus CMV" 若 genus 段与 genus 列一致 → 取种加词部分
    返回规范化字符串；无法判断时原样返回，交 taxonomy_check 在线校验。
    """
    s = str(val or '').strip()
    return re.sub(r'\s+', ' ', s)


def _infer_organism(row):
    """分级分类 → GenBank organism。

    优先级：species(规范) → genus + ' sp.' → family + ' sp.' → taxon → ''。
    注意 NCBI 对未培养病毒 (UViG) 通常只接受属级 "Genus sp."，
    若 species 以 'virus'/'viroid' 结尾或含 sp. 则保留 species。
    """
    species = _clean_species_name(row.get('species'))
    genus = _clean_species_name(row.get('genus'))
    family = _clean_species_name(row.get('family'))
    taxon = _clean_species_name(row.get('taxon'))
    # species 有效且不像"属 种加词 缩写"混排（如 "Cucumovirus CMV"）
    if species and not is_placeholder(species):
        # 若 species == 单个词且 genus 存在 → "Genus sp."
        if ' ' not in species:
            if genus and not is_placeholder(genus):
                return f'{genus} sp.'
            return species
        return species
    if genus and not is_placeholder(genus):
        return f'{genus} sp.'
    if family and not is_placeholder(family):
        return f'{family} sp.'
    return taxon


# 导入时可自动回填的提交表列（value = 分类表同名列）
_RUN_FIELD_MAP = {
    # 提交表列          分类表列（不存在则跳过）
    'src-Host': 'host',
}


def import_from_run(name, run):
    """把 tool_runs/<run>/virus_classification.tsv 的病毒 contigs 导入为表行。

    sequence_name=contig，organism 由分级分类智能推断（species→genus sp.→family sp.）。
    表内已有同名 sequence_name 的行跳过。
    返回 {added, skipped, run, lineage_cols, rows:[{contig, organism, family, near_complete}]}
    """
    from ..config import PLATFORM_ROOT
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run or ''):
        raise ValueError(f'非法运行名: {run!r}')
    root = DIRS.get('tool_runs') or os.path.join(PLATFORM_ROOT, 'tool_runs')
    tsv = check_path(os.path.join(root, run, 'virus_classification.tsv'),
                     must_exist=True, in_platform=True)
    tdf = pd.read_csv(tsv, sep='\t', dtype=str, keep_default_na=False)
    for col in ('contig',):
        if col not in tdf.columns:
            raise KeyError(f"{run} 的分类表缺少 {col} 列")

    df = load_table(name)
    existing = {str(v or '').strip() for v in df['sequence_name'].tolist()} \
        if 'sequence_name' in df.columns else set()
    cols = list(UNIFIED_COLUMNS.keys())
    for c in cols:
        if c not in df.columns:
            df[c] = ''
    added = skipped = 0
    details = []
    for _, r in tdf.iterrows():
        contig = str(r.get('contig', '') or '').strip()
        if not contig or contig in existing:
            skipped += 1
            continue
        row = {c: '' for c in cols}
        row['sequence_name'] = contig
        row['organism'] = _infer_organism(r)
        # 可选自动回填（如宿主列存在时）
        for tbl_col, run_col in _RUN_FIELD_MAP.items():
            v = str(r.get(run_col, '') or '').strip()
            if v and not is_placeholder(v):
                row[tbl_col] = v
        df.loc[len(df)] = row
        existing.add(contig)
        added += 1
        details.append({
            'contig': contig,
            'organism': row['organism'],
            'family': str(r.get('family', '') or '').strip(),
            'near_complete': str(r.get('near_complete', '') or '').strip(),
        })
    save_df(df[_align_columns(df).columns], _csv_path(name))
    lineage_cols = [c for c in
                    ('realm', 'kingdom', 'phylum', 'class', 'order',
                     'family', 'genus', 'species')
                    if c in tdf.columns]
    # 记录来源 contigs 运行（供序列 fasta 自动推断）
    _write_link(name, 'contigs', {'run': run, 'rows': details})
    return {'added': added, 'skipped': skipped, 'run': run,
            'lineage_cols': lineage_cols, 'rows': details}


def list_contig_runs():
    """有 virus_classification.tsv 的 contigs 运行列表（导入下拉用）。"""
    from ..config import PLATFORM_ROOT
    root = DIRS.get('tool_runs') or os.path.join(PLATFORM_ROOT, 'tool_runs')
    out = []
    if os.path.isdir(root):
        for nm in sorted(os.listdir(root), reverse=True):
            if not nm.startswith('contigs_'):
                continue
            if os.path.isfile(os.path.join(root, nm,
                                           'virus_classification.tsv')):
                out.append(nm)
    return out[:30]


def list_orf_runs():
    """有 CDS 注释产物（04b_orf_annot/orf_annotation.tsv 或 04_orf/*.gff）的
    orf/orfa 运行列表（"关联注释"下拉用）。"""
    from ..config import PLATFORM_ROOT
    root = DIRS.get('tool_runs') or os.path.join(PLATFORM_ROOT, 'tool_runs')
    out = []
    if os.path.isdir(root):
        for nm in sorted(os.listdir(root), reverse=True):
            if not (nm.startswith('orf_') or nm.startswith('orfa_')):
                continue
            d = os.path.join(root, nm)
            hit = (os.path.isfile(os.path.join(d, '04b_orf_annot',
                                               'orf_annotation.tsv')) or
                   list(Path(d).glob('04_orf/*.gff')) or
                   list(Path(d).glob('04b_orf_annot/*.gff*')))
            if hit:
                out.append(nm)
    return out[:30]


def _orf_run_dir(run):
    from ..config import PLATFORM_ROOT
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run or ''):
        raise ValueError(f'非法运行名: {run!r}')
    root = DIRS.get('tool_runs') or os.path.join(PLATFORM_ROOT, 'tool_runs')
    return check_path(os.path.join(root, run), must_exist=True,
                      in_platform=True)


def _link_file(name, kind):
    """项目内记录"来自哪个运行/序列文件"的小 json（不污染主表列）。"""
    return os.path.join(_table_dir(name), f'_link_{kind}.json')


def _write_link(name, kind, payload):
    import json
    p = _link_file(name, kind)
    with safe_open(p, 'wt') as f:
        json.dump(payload, f, ensure_ascii=False, indent=1)
    return p


def _read_link(name, kind):
    import json
    p = _link_file(name, kind)
    if not os.path.isfile(p):
        return None
    try:
        with open(p, 'r', encoding='utf-8') as f:
            return json.load(f)
    except (ValueError, OSError):
        return None


def link_orf_run(name, orf_run):
    """把 orf/orfa 注释运行关联到提交项目（记录来源 + 返回可回填的 CDS/宿主信息）。

    与 contigs 分类运行是"独立跑"的——以 contig 名（sequence_name）为对接键。
    返回 {linked: True, run, contigs:[...], gffs, n_cds, host_sources, missing}
    missing = 表内 sequence_name 在 orf 注释里找不到的 contig（失配警告）。
    """
    d = _orf_run_dir(orf_run)
    ann_tsv = os.path.join(d, '04b_orf_annot', 'orf_annotation.tsv')
    cds_info = None
    if os.path.isfile(ann_tsv):
        cdf = pd.read_csv(ann_tsv, sep='\t', dtype=str, keep_default_na=False)
        if 'contig' in cdf.columns:
            cds_info = cdf
    # gff 候选（04_orf/pyrodigal*.gff 或 04b 的 gff3）
    gffs = []
    for pat in ('04_orf/*.gff', '04b_orf_annot/*.gff3',
                '04b_orf_annot/*.gff'):
        gffs += [str(p) for p in Path(d).glob(pat)]
    # 04b 的 orf_annotation.tsv 优先于 04_orf 的 pyrodigal.gff
    df = load_table(name)
    names = [str(v or '').strip() for v in df['sequence_name'].tolist()
             if str(v or '').strip() and not is_placeholder(str(v or ''))]
    want = set(names)
    n_cds = 0
    found = set()
    host_sources = {}
    profile_tsv = os.path.join(d, '04b_orf_annot',
                               'contig_function_profile.tsv')
    prof = None
    if os.path.isfile(profile_tsv):
        prof = pd.read_csv(profile_tsv, sep='\t', dtype=str,
                           keep_default_na=False)
        if 'contig' in prof.columns:
            hm = {}
            if 'Host_source' in prof.columns:
                hm = dict(zip(prof['contig'], prof['Host_source']))
            elif 'host_source' in prof.columns:
                hm = dict(zip(prof['contig'], prof['host_source']))
            host_sources = {str(k): str(v) for k, v in hm.items()
                            if str(v).strip()}
    if cds_info is not None:
        sub = cds_info[cds_info['contig'].isin(want)]
        n_cds = int(len(sub))
        found = set(sub['contig'].tolist())
        # 宿主回填（供前端展示，不自动改表——宿主需人工确认）
    missing = sorted(n for n in names if n not in found and n not in host_sources)
    _write_link(name, 'orf', {
        'run': orf_run,
        'annotation_tsv': ann_tsv if os.path.isfile(ann_tsv) else None,
        'gffs': gffs[:5],
        'n_contigs_table': len(names),
        'n_cds_matched': n_cds,
        'host_sources': host_sources,
        'missing': missing,
    })
    return {'linked': True, 'run': orf_run,
            'n_cds': n_cds, 'host_sources': host_sources,
            'missing': missing[:50], 'n_missing': len(missing),
            'gffs': gffs[:5]}


# ══════════════════════════════════════════════════════════════
# CDS 特征表生成（把 orf 运行的 CDS 注释 → GenBank feature table .tbl）
# ══════════════════════════════════════════════════════════════

_TBL_NO_CDS = re.compile(r'[^A-Za-z0-9_\-]')


def _parse_gff3_cds(gff_path):
    """gff3 → [(seqid, start, end, strand)]（仅 CDS 行；坐标 1-based 保留）。"""
    out = []
    try:
        with open(gff_path, 'r', encoding='utf-8', errors='replace') as f:
            for ln in f:
                if ln.startswith('#') or not ln.strip():
                    continue
                parts = ln.rstrip('\n').split('\t')
                if len(parts) < 8:
                    continue
                if parts[2] != 'CDS':
                    continue
                try:
                    out.append((parts[0], int(parts[3]), int(parts[4]),
                                parts[6]))
                except ValueError:
                    continue
    except OSError:
        pass
    return out


def generate_feature_tbl(name):
    """按已关联 orf 运行的 CDS 注释生成 GenBank 5 列 feature table。

    产物：项目目录 featuretable.tbl（每 contig 一个 >Feature 块，仅含表内
    sequence_name 的 contig；坐标来自 gff3，product 来自 orf_annotation.tsv，
    缺失 product 的 CDS 标 hypothetical protein）。

    返回 {path, contigs, cds_total, no_product, skipped, sources}
    """
    link = _read_link(name, 'orf')
    if not link or not link.get('run'):
        raise ValueError('尚未关联 orf 注释运行——先点"关联注释"')
    run = link['run']
    d = _orf_run_dir(run)

    # 1) CDS 坐标：优先该 orf 运行功能注释所采用的基因模型（默认 pyrodigal_rv，
    #    04b summary 记录实际 model），再退回 04b gff3 / 其余 04_orf gff（历史运行）
    import json as _json
    model = None
    summ = os.path.join(d, '04b_orf_annot', 'summary.json')
    if os.path.isfile(summ):
        try:
            with safe_open(summ) as _f:
                model = _json.load(_f).get('model')
        except (OSError, ValueError):
            model = None
    if model not in ('pyrodigal_rv', 'pyrodigal'):
        model = None
    pref = (model,) if model else ('pyrodigal_rv', 'pyrodigal')
    gffs = []
    for m in pref:
        gffs += sorted(str(p) for p in Path(d).glob(f'04_orf/{m}.gff'))
    for pat in ('04_orf/*.gff', '04_orf/*.gff3',
                '04b_orf_annot/*.gff3', '04b_orf_annot/*.gff'):
        gffs += sorted(str(p) for p in Path(d).glob(pat))
    _seen, gffs_uniq = set(), []
    for g in gffs:
        if g not in _seen:
            _seen.add(g)
            gffs_uniq.append(g)
    gffs = gffs_uniq
    cdses = []
    for g in gffs:
        cdses = _parse_gff3_cds(g)
        if cdses:
            break
    if not cdses:
        raise ValueError(f'运行 {run} 中未找到 CDS 坐标（04_orf/*.gff）')

    # 2) product 字典：(contig, start) → product（tsv 与 gff 可能差 1bp，双向容错）
    prod = {}
    ann_tsv = os.path.join(d, '04b_orf_annot', 'orf_annotation.tsv')
    if os.path.isfile(ann_tsv):
        cdf = pd.read_csv(ann_tsv, sep='\t', dtype=str, keep_default_na=False)
        if {'contig', 'start', 'end'}.issubset(cdf.columns):
            for _, r in cdf.iterrows():
                try:
                    s0 = int(float(r['start']))
                except (ValueError, TypeError):
                    continue
                p = str(r.get('product', '') or '').strip()
                if p:
                    prod.setdefault((str(r['contig']), s0), p)
                    prod.setdefault((str(r['contig']), s0 - 1), p)
                    prod.setdefault((str(r['contig']), s0 + 1), p)

    # 3) 按表行过滤 + 写 tbl
    df = load_table(name)
    names = [str(v or '').strip() for v in df['sequence_name'].tolist()
             if str(v or '').strip() and not is_placeholder(str(v or ''))]
    want = set(names)
    out = os.path.join(_table_dir(name), 'featuretable.tbl')
    written_contigs, cds_total, no_product = [], 0, 0
    with safe_open(out, 'wt') as f:
        for cid in names:                     # 保持表顺序
            rows = [(s, e, st) for (q, s, e, st) in cdses if q == cid]
            if not rows:
                continue
            written_contigs.append(cid)
            f.write(f'>Feature {cid}\n')
            for s, e, st in sorted(rows):
                product = prod.get((cid, s)) or prod.get((cid, e)) or ''
                if not product:
                    product = 'hypothetical protein'
                    no_product += 1
                f.write(f'{s}\t{e}\tCDS\n')
                # 注意：table2asn 对带斜杠+双引号的 /product="xxx" 解析失败，
                # 会丢 product 名并报 MissingProteinName×N Error。须用 suvtk
                # 同款无斜杠、无引号、tab 分隔格式（product\t名字）。
                f.write(f'\t\t\tproduct\t{product}\n')
                f.write(f'\t\t\tcodon_start\t1\n')
                f.write(f'\t\t\ttransl_table\t1\n')
                cds_total += 1
    skipped = sorted(n for n in names if n not in written_contigs)
    if not written_contigs:
        raise ValueError(
            '表内没有 contig 能在该 orf 运行中找到 CDS（全部失配）。\n'
            '原因通常是 contigs 分类运行与 orf 注释运行不是同一条序列来源——'
            '请核对：① 导入的 contigs 运行 与 ② 关联的 orf 运行 的序列名（sequence_name）'
            '必须一致；orf 运行应是用该 contigs 的 viral_contigs.fasta 跑的。\n'
            f'表内 {len(names)} 个 sequence_name 示例: {skipped[:3]} …')
    return {'path': out, 'contigs': written_contigs, 'cds_total': cds_total,
            'no_product': no_product, 'skipped': skipped,
            'sources': {'gff': gffs[0] if gffs else None,
                        'annotation': ann_tsv if os.path.isfile(ann_tsv) else None}}


# ══════════════════════════════════════════════════════════════
# 打包下载（NCBI portal / BankIt 支持整包上传）
# ══════════════════════════════════════════════════════════════

def package_zip(name):
    """项目目录内全部提交产物 → submission_package.zip（排除中枢表与旧包）。"""
    import zipfile
    d = _table_dir(name)
    out = os.path.join(d, 'submission_package.zip')
    skip = {'unified_metadata.csv', 'submission_package.zip'}
    n = 0
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as z:
        for cur, _sub, fns in os.walk(d):
            for fn in fns:
                if fn in skip:
                    continue
                p = os.path.join(cur, fn)
                z.write(p, os.path.relpath(p, d).replace('\\', '/'))
                n += 1
    return out, n


# ══════════════════════════════════════════════════════════════
# 物种名在线校验（NCBI Taxonomy；联网可选）
# ══════════════════════════════════════════════════════════════

def taxonomy_check(name):
    """逐个 unique organism 查 NCBI Taxonomy（esearch）。

    NCBI 拒绝未收录的物种名（新种需先处理），提前校验避免提交被退。
    返回 [{organism, found, taxid}]；网络失败抛 RuntimeError 由端点转 502。
    """
    from ..ncbi_download import esearch
    df = load_table(name)
    if 'organism' not in df.columns:
        raise KeyError('缺少 organism 列')
    orgs = sorted({str(v).strip() for v in df['organism'].tolist()
                   if str(v).strip() and not is_placeholder(v)})
    out = []
    for org in orgs:
        res = esearch(f'{org}[SCIN]', db='taxonomy', retmax=1)
        out.append({'organism': org,
                    'found': bool(res['count']),
                    'taxid': (res['ids'][0] if res['ids'] else '')})
    return out


# ══════════════════════════════════════════════════════════════
# 示例数据（同 GUI 的 Demo 数据，self-sequenced 场景）
# ══════════════════════════════════════════════════════════════

_DEMO_ROWS = [
    {
        "organism": "Betacytorhabdovirus lycii",
        "sequence_name": "CRR123456_Betacytorhabdovirus_lycii_contig1",
        "authors": "Zhang, Wenda; Li, Ming; Wang, Fang",
        "collection_date": "YYYY-MM-DD",
        "bioproject": "PRJNAXXXXXX",
        "src-Isolate": "Betacytorhabdovirus_lycii_CRR123456",
        "src-geo_loc_name": "China:Ningxia",
        "src-Lat_Lon": "38.47 N 106.27 E",
        "src-Host": "Lycium barbarum",
        "src-Segment": "",
        "src-Isolation-source": "plant virome",
        "src-Note": "",
        "src-Tissue_type": "root",
        "src-Collected_by": "Ningxia University",
        "src-Cultivar": "Ningqi No.5",
        "src-Dev_stage": "2 years",
        "gb-sample_name": "Betacytorhabdovirus_lycii_CRR123456",
        "gb-title": "",
        "sra": "CRR123456",
        "biosample": "SAMNXXXXXXXX",
        "cmt-Assembly_Method": "SPAdes;4.3.0;metaviral",
        "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
        "cmt-Genome_Coverage": "42.5x",
        "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
        "bs-isolate": "Betacytorhabdovirus_lycii_CRR123456",
        "bs-geo_loc_name": "China:Ningxia",
        "bs-host": "Lycium barbarum",
        "bs-isolation_source": "plant virome",
    },
    {
        "organism": "Betacytorhabdovirus lycii",
        "sequence_name": "CRR123456_Betacytorhabdovirus_lycii_contig2",
        "authors": "Zhang, Wenda; Li, Ming; Wang, Fang",
        "collection_date": "YYYY-MM-DD",
        "bioproject": "PRJNAXXXXXX",
        "src-Isolate": "Betacytorhabdovirus_lycii_CRR123456",
        "src-geo_loc_name": "China:Ningxia",
        "src-Lat_Lon": "38.47 N 106.27 E",
        "src-Host": "Lycium barbarum",
        "src-Segment": "2",
        "src-Isolation-source": "plant virome",
        "src-Note": "segment 2 of multipartite virus",
        "src-Tissue_type": "root",
        "src-Collected_by": "Ningxia University",
        "src-Cultivar": "Ningqi No.5",
        "src-Dev_stage": "2 years",
        "gb-sample_name": "Betacytorhabdovirus_lycii_CRR123456_seg2",
        "gb-title": "",
        "sra": "CRR123456",
        "biosample": "SAMNXXXXXXXX",
        "cmt-Assembly_Method": "SPAdes;4.3.0;metaviral",
        "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
        "cmt-Genome_Coverage": "38.1x",
        "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
        "bs-isolate": "Betacytorhabdovirus_lycii_CRR123456",
        "bs-geo_loc_name": "China:Ningxia",
        "bs-host": "Lycium barbarum",
        "bs-isolation_source": "plant virome",
    },
    {
        "organism": "Torradovirus Ningxiaense",
        "sequence_name": "SRR31651831_Torradovirus_contig1",
        "authors": "Zhang, Wenda; Li, Ming; Wang, Fang",
        "collection_date": "2024-07-03",
        "bioproject": "PRJNA1218117",
        "src-Isolate": "Torradovirus_Ningxiaense_SRR31651831",
        "src-geo_loc_name": "China:Ningxia:Yinchuan",
        "src-Lat_Lon": "38.47 N 106.27 E",
        "src-Host": "Lycium barbarum",
        "src-Segment": "",
        "src-Isolation-source": "plant virome",
        "src-Note": "",
        "src-Tissue_type": "leaf",
        "src-Collected_by": "Ningxia University",
        "src-Cultivar": "Ningqi No.7",
        "src-Dev_stage": "3 years",
        "gb-sample_name": "Torradovirus_Ningxiaense_SRR31651831",
        "gb-title": "Torradovirus Ningxiaense genome sequencing",
        "sra": "SRR31651831",
        "biosample": "SAMN56789012",
        "cmt-Assembly_Method": "SPAdes;4.3.0;metaviral",
        "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
        "cmt-Genome_Coverage": "56.3x",
        "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
        "bs-isolate": "Torradovirus_Ningxiaense_SRR31651831",
        "bs-geo_loc_name": "China:Ningxia:Yinchuan",
        "bs-host": "Lycium barbarum",
        "bs-isolation_source": "plant virome",
    },
    {
        "organism": "Mint virus X",
        "sequence_name": "SRR33301106_MintVirusX_contig1",
        "authors": "Zhang, Wenda; Li, Ming; Wang, Fang",
        "collection_date": "2024-09-01",
        "bioproject": "PRJNA1219886",
        "src-Isolate": "MintVirusX_SRR33301106",
        "src-geo_loc_name": "China:Inner Mongolia:Alxa",
        "src-Lat_Lon": "39.08 N 105.73 E",
        "src-Host": "Lycium ruthenicum",
        "src-Isolation-source": "plant virome",
        "src-Tissue_type": "fruit",
        "src-Collected_by": "Inner Mongolia University",
        "src-Cultivar": "wild",
        "src-Dev_stage": "3 years, mature fruit",
        "gb-sample_name": "MintVirusX_SRR33301106",
        "gb-title": "Mint virus X from Lycium ruthenicum",
        "sra": "SRR33301106",
        "biosample": "SAMN67890123",
        "cmt-Assembly_Method": "SPAdes;4.3.0;metaviral",
        "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
        "cmt-Genome_Coverage": "23.7x",
        "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
        "bs-isolate": "MintVirusX_SRR33301106",
        "bs-geo_loc_name": "China:Inner Mongolia:Alxa",
        "bs-host": "Lycium ruthenicum",
        "bs-isolation_source": "plant virome",
    },
    {
        "organism": "unclassified virus",
        "sequence_name": "SRR33389501_novel_virus_contig1",
        "authors": "Author, First",
        "collection_date": "2023-08-10",
        "bioproject": "PRJNAXXXXXX",
        "src-Isolate": "novel_virus_SRR33389501",
        "src-geo_loc_name": "Country:Region",
        "src-Lat_Lon": "XX.XX N XXX.XX E",
        "src-Isolation-source": "plant virome",
        "src-Note": "novel virus, no close reference",
        "src-Tissue_type": "leaf",
        "gb-sample_name": "novel_virus_SRR33389501",
        "gb-title": "",
        "sra": "SRR33389501",
        "biosample": "SAMNXXXXXXXX",
        "cmt-Assembly_Method": "SPAdes;4.3.0;metaviral",
        "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
        "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
        "bs-isolate": "novel_virus_SRR33389501",
        "bs-geo_loc_name": "Country:Region",
        "bs-isolation_source": "plant virome",
    },
    {
        "organism": "Potexvirus lycii",
        "sequence_name": "SRR30124789_Potexvirus_lycii_contig1",
        "authors": "Zhang, Wenda; Li, Ming; Wang, Fang",
        "collection_date": "2024-06-15",
        "bioproject": "PRJNA1189201",
        "src-Isolate": "Potexvirus_lycii_SRR30124789",
        "src-geo_loc_name": "China:Xinjiang:Urumqi",
        "src-Lat_Lon": "43.79 N 87.58 E",
        "src-Host": "Lycium barbarum",
        "src-Isolation-source": "plant virome",
        "src-Tissue_type": "fruit",
        "src-Collected_by": "Xinjiang University",
        "src-Cultivar": "cultivated",
        "src-Dev_stage": "5 years, ripening",
        "gb-sample_name": "Potexvirus_lycii_SRR30124789",
        "gb-title": "Potexvirus lycii from goji berry",
        "sra": "SRR30124789",
        "biosample": "SAMN78901234",
        "cmt-Assembly_Method": "SPAdes;4.3.0;metaviral",
        "cmt-Sequencing_Technology": "Illumina NovaSeq 6000",
        "cmt-Genome_Coverage": "61.2x",
        "cmt-Annotation_Pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1",
        "bs-isolate": "Potexvirus_lycii_SRR30124789",
        "bs-geo_loc_name": "China:Xinjiang:Urumqi",
        "bs-host": "Lycium barbarum",
        "bs-isolation_source": "plant virome",
    },
]

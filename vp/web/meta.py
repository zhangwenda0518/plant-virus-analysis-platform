# -*- coding: utf-8 -*-
"""公共数据检索（自 app.py 拆出）。

SRA+GSA 双引擎按物种检索 Run → 统一元数据 → 出版级图组；
含在线表格编辑、AI 元数据清洗（DeepSeek 可选）与项目归档。"""
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid

from flask import (Blueprint, abort, jsonify, render_template, request,
                   send_file, send_from_directory)

from vp.config import DIRS, PLATFORM_ROOT, db_path, engine_cmd
from vp.utils import (TaskLogger, check_path, fmt_size, run_cmd, safe_open,
                      safe_remove)
from vp.web.state import cfg, tool_runs_root as _tool_runs_root
from vp.web.tasks import tm

bp = Blueprint('meta', __name__)


_META_TABLE_FILES = {
    'search': ('search', 'SRA_GSA_Merged_Final.csv'),
    'core14': ('info', 'Global_Unified_Metadata_Core14.csv'),
    'full': ('info', 'Global_Unified_Metadata_Full.csv'),
}


_META_EDITABLE = ('search', 'core14')


_META_NA = {'na', 'n/a', 'not_provided', 'not_provided.', 'not collected',
            'missing', 'none', 'unknown', 'nan', '', '-', ' '}


_meta_df_cache = {}          # (name, which) -> [mtime, df]


_meta_cache_lock = threading.Lock()


AI_SANITIZE_PROMPT = """你是一个极其严谨的生物信息元数据清理程序。你的任务是【去污染、归位、格式化】，绝对不是【凭空捏造】。
我将给你一条 JSON 数据。请严格按以下规则逐字段处理：

1. **全面净化**：遇到 'not collected', 'missing', 'N/A', 'Not_Provided', 空字符串等无意义占位符，替换为 'Not_Provided'。

2. **ScientificName**：如果混入了部位（如 leaf）或来源词（如 wild），将其剔除，只保留纯物种名。

3. **Tissue**：强制小写单数（leaves→leaf, fruits→fruit, flowers→flower, roots→root）。如果原值为空但 Source 恰好是真正的部位词（如 leaf, root），才移动过来。Source 是品种名（如 Ningqi No.1）则不移动。

4. **Source**：如果被填成了物种名，清空为 Not_Provided。如果是 "wild", "cultivated" 或品种名（如 Ningqi No.1），必须 100% 保持原样。

5. **Location**：根据 CenterName 推断机构地理位置，输出严格三级格式【国家, 省/州, 市/县_AI】。参照示例：
   "Beijing Forestry University" → "China, Beijing, Beijing_AI"
   "Ningxia University" → "China, Ningxia, Yinchuan_AI"
   "North Minzu University" → "China, Ningxia, Yinchuan_AI"
   "University of Tokyo" → "Japan, Tokyo, Tokyo_AI"
   "USDA-ARS" → "USA, Maryland, Beltsville_AI"
   "Royal Botanic Gardens Kew" → "United Kingdom, England, London_AI"
   "Northwest A&F University" → "China, Shaanxi, Yangling_AI"
   "Gansu Agricultural University" → "China, Gansu, Lanzhou_AI"
   "Xinjiang University" → "China, Xinjiang, Urumqi_AI"
   "Qinghai University" → "China, Qinghai, Xining_AI"
   "Inner Mongolia University" → "China, Inner Mongolia, Hohhot_AI"
   "Henan Agricultural University" → "China, Henan, Zhengzhou_AI"
   不认识的机构保持 Not_Provided，绝对不要输出 Unknown。

6. **Age_GrowthStage**：标准化年龄和发育阶段文本。规则：
   - "3 year" → "3 years"（补全复数）
   - "2 years" → 保持原样（已规范）
   - "Archeocyte stage" → "archeocyte stage"（统一为首字母不大写的全小写，除非是专有名词）
   - 多个阶段用 "|" 分隔，如 "3 years | mature fruit stage"
   - 将中文风格描述转为英文，如 "flowering stage"、"ripening stage"
   - 如果为空则保持 Not_Provided

7. **LibrarySource**：如果传入长文本，浓缩为一个词：TRANSCRIPTOMIC / GENOMIC / METAGENOMIC。

8. **CollectionDate, CenterName, BioProject**：必须 100% 保持输入原样，一个字都不许改！

【核心纪律】：只做减法（去污染）、归位（移动）、格式标准化，绝对禁止凭空捏造！Location 没线索就 Not_Provided！

对输入 JSON 中出现的以下字段逐个输出清洗结果，直接输出合法 JSON（勿带 ```json 标记），未出现的字段不要输出：
{"CollectionDate": "...", "Location": "...", "Source": "...", "Tissue": "...",
 "Age_GrowthStage": "...", "ScientificName": "...", "LibrarySource": "...",
 "CenterName": "...", "BioProject": "..."}"""


AI_FILL_COLS = ('CollectionDate', 'Location', 'Source', 'Tissue',
                'Age_GrowthStage', 'ScientificName', 'LibrarySource',
                'CenterName', 'BioProject')


_AI_EMPTY = {'', ' ', 'na', 'n/a', 'not_provided', 'not collected',
             'missing', 'none', 'unknown', 'nan', '-'}


def _meta_slug(species):
    slug = re.sub(r'[^\w\-]+', '_', str(species)).strip('_') or 'species'
    return slug


def _meta_proj(name):
    return check_path(os.path.join(DIRS['meta_search'], name),
                      must_exist=False, in_platform=True)


def _meta_table_path(name, which):
    if which not in _META_TABLE_FILES:
        return None
    sub, fname = _META_TABLE_FILES[which]
    return os.path.join(_meta_proj(name), sub, fname)


def _meta_load_df(name, which):
    """读取项目表（带 mtime 缓存，编辑后由 _meta_store_df 刷新）。"""
    import pandas as pd
    path = _meta_table_path(name, which)
    if not path or not os.path.isfile(path):
        return None
    key = (name, which)
    mtime = os.path.getmtime(path)
    with _meta_cache_lock:
        ent = _meta_df_cache.get(key)
        if ent and ent[0] == mtime:
            return ent[1]
    df = pd.read_csv(path, encoding='utf-8-sig', dtype=str,
                     keep_default_na=False)
    df = df.loc[:, ~df.columns.duplicated()]
    with _meta_cache_lock:
        _meta_df_cache[key] = [mtime, df]
    return df


def _meta_store_df(name, which, df):
    """保存项目表（csv+tsv 双格式，与引擎 save_dual_format 一致）。"""
    sub, fname = _META_TABLE_FILES[which]
    base = check_path(os.path.join(_meta_proj(name), sub),
                      must_exist=False, in_platform=True)
    os.makedirs(base, exist_ok=True)
    stem = os.path.splitext(fname)[0]
    df = df.loc[:, ~df.columns.duplicated()]
    csv_p = os.path.join(base, stem + '.csv')
    with safe_open(csv_p, 'wt') as f:
        f.write('\ufeff')
        df.to_csv(f, index=False, lineterminator='\n')
    with safe_open(os.path.join(base, stem + '.tsv'), 'wt') as f:
        f.write('\ufeff')
        df.to_csv(f, index=False, sep='\t', lineterminator='\n')
    with _meta_cache_lock:
        _meta_df_cache[(name, which)] = [os.path.getmtime(csv_p), df]


def _meta_missing_count(df):
    return int(sum(df[c].astype(str).str.strip().str.lower()
                   .isin(_META_NA).sum() for c in df.columns))


def _meta_web_cache_dirs(name):
    """项目内全部 0_web_cache 目录（info / search 引擎的网页缓存）。"""
    out = []
    root = _meta_proj(name)
    if not os.path.isdir(root):
        return out
    for cur, sub, _fns in os.walk(root):
        if os.path.basename(cur) == '0_web_cache':
            out.append(cur)
    return out


def _meta_deepseek_key():
    return os.environ.get('DEEPSEEK_API_KEY', '').strip()


def _meta_ai_client(api_key, base='https://api.deepseek.com'):
    from openai import OpenAI
    return OpenAI(api_key=api_key, base_url=base)


def _meta_ai_chat(client, model, system, user, timeout=180):
    kwargs = {'model': model, 'timeout': timeout, 'messages': [
        {'role': 'system', 'content': system},
        {'role': 'user', 'content': user}]}
    if 'v4' in model.lower():
        if 'pro' in model.lower():
            kwargs['reasoning_effort'] = 'high'
            kwargs['extra_body'] = {'thinking': {'type': 'enabled'}}
        else:
            kwargs['temperature'] = 0.1
    else:
        kwargs['temperature'] = 0.3
    resp = client.chat.completions.create(**kwargs)
    return resp.choices[0].message.content or ''


def _meta_stream_task(task_name, cmd, done_msg, log_file=None, link=None):
    """子进程跑 meta 脚本，stdout 实时进任务日志；完成返回 summary。"""
    def fn(log, progress, cancel):
        from vp.utils import task_register_proc
        env = dict(os.environ)
        env.setdefault('PYTHONUNBUFFERED', '1')   # 引擎脚本日志实时流出
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True,
                                encoding='utf-8', errors='replace', env=env)
        task_register_proc(proc)                  # 登记 → 取消时硬停止
        try:
            for line in proc.stdout:
                if cancel.is_set():
                    proc.terminate()
                    raise RuntimeError('已取消')
                if line.strip():
                    log(line.rstrip())
            rc = proc.wait()
        finally:
            try:
                proc.stdout.close()
            except Exception:
                pass
        if rc != 0:
            raise RuntimeError(f'命令退出码 {rc}（详见日志）')
        return {'kind': 'meta', 'msg': done_msg}
    return tm.start(task_name, fn, log_file=log_file, link=link)


@bp.route('/meta')
def page_meta():
    return render_template('meta.html')


def _csv_data_rows(path):
    """CSV 数据行数（按物理行快算，字段不含内嵌换行的场景下准确；
    供卡片徽章展示用，避免每次轮询都 pandas 全量解析）。"""
    try:
        n = sum(1 for _ in safe_open(path))
        return max(0, n - 1)
    except OSError:
        return None


@bp.route('/api/meta/collections')
def api_meta_collections():
    root = check_path(DIRS['meta_search'], must_exist=False, in_platform=True)
    out = []
    if os.path.isdir(root):
        for name in sorted(os.listdir(root)):
            d = os.path.join(root, name)
            if not os.path.isdir(d):
                continue
            search_f = os.path.join(d, 'search', 'SRA_GSA_Merged_Final.csv')
            core_f = os.path.join(d, 'info', 'Global_Unified_Metadata_Core14.csv')
            full_f = os.path.join(d, 'info', 'Global_Unified_Metadata_Full.csv')
            n_runs = _csv_data_rows(search_f) if os.path.isfile(search_f) else None
            n_core = _csv_data_rows(core_f) if os.path.isfile(core_f) else None
            n_full = _csv_data_rows(full_f) if os.path.isfile(full_f) else None
            plot_dir = os.path.join(d, 'plot')
            n_plots = 0
            if os.path.isdir(plot_dir):
                n_plots = sum(1 for f in os.listdir(plot_dir)
                              if f.lower().endswith('.png'))
            # 过期检测：元数据条数 ≠ 最新检索条数 → 提示重新提取
            stale = (n_runs is not None and n_core is not None
                     and n_core != n_runs)
            out.append({'name': name,
                        'search': os.path.isfile(search_f),
                        'core14': os.path.isfile(core_f),
                        'full': os.path.isfile(full_f),
                        'n_runs': n_runs, 'n_core14': n_core,
                        'n_full': n_full, 'n_plots': n_plots,
                        'stale': stale,
                        'mtime': int(os.path.getmtime(d))})
    out.sort(key=lambda x: x['mtime'], reverse=True)
    return jsonify(out)


@bp.route('/api/meta/search', methods=['POST'])
def api_meta_search():
    body = request.get_json(force=True) or {}
    species = str(body.get('species') or '').strip()
    if not species:
        abort(400, '需要物种拉丁名')
    source = str(body.get('source') or 'TRANSCRIPTOMIC')
    db = str(body.get('db') or 'both')
    detailed = bool(body.get('detailed', True))
    workers = body.get('workers') or 5
    try:
        workers = min(max(1, int(workers)), 20)
    except (TypeError, ValueError):
        workers = 5
    # 每次检索自动新建 <物种>_<时间戳> 项目；也允许 body.collection 指定项目重检
    name = str(body.get('collection') or '').strip()
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name or ''):
        ts = time.strftime('%Y%m%d_%H%M')
        name = f'{_meta_slug(species)}_{ts}'
        n = 2
        while os.path.isdir(os.path.join(DIRS['meta_search'], name)):
            name = f'{_meta_slug(species)}_{ts}_{n}'
            n += 1
    out_dir = check_path(os.path.join(DIRS['meta_search'], name, 'search'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    cmd = engine_cmd(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'vp', 'public_meta', 'search_engine.py'),
                     '-q', species, '-s', source, '--db', db,
                     '--workers', str(workers), '-o', out_dir)
    if not detailed:
        cmd.append('--no-detailed')
    if os.environ.get('NCBI_API_KEY'):
        cmd += ['--ncbi-api', os.environ['NCBI_API_KEY']]
    tid = _meta_stream_task(
        f'公共检索 {species}', cmd, f'检索完成: {out_dir}',
        log_file=os.path.join(out_dir, 'run.log'), link='/meta')
    return jsonify({'task': tid, 'name': name})


@bp.route('/api/meta/info', methods=['POST'])
def api_meta_info():
    body = request.get_json(force=True) or {}
    name = str(body.get('collection') or '')
    try:
        coll = check_path(os.path.join(DIRS['meta_search'], name),
                          must_exist=True, in_platform=True)
    except (ValueError, FileNotFoundError):
        abort(400, f'检索项目不存在: {name}')
    runs = body.get('runs') or None
    list_f = os.path.join(coll, 'info_runs.txt')
    if runs:
        if not isinstance(runs, list) or not all(re.fullmatch(r'[A-Za-z0-9_]+', str(r)) for r in runs):
            abort(400, 'runs 需为编号列表')
        with safe_open(list_f, 'wt') as f:
            f.write('\n'.join(str(r) for r in runs) + '\n')
        input_arg = list_f
    else:
        search_f = os.path.join(coll, 'search', 'SRA_GSA_Merged_Final.csv')
        if not os.path.isfile(search_f):
            abort(400, '该项目还没有检索结果（先运行检索）')
        import pandas as pd
        try:
            df = pd.read_csv(search_f, encoding='utf-8-sig', dtype=str)
        except Exception as e:
            abort(400, f'检索结果解析失败: {e}')
        if 'Run' not in df.columns:
            abort(400, '检索结果缺少 Run 列')
        with safe_open(list_f, 'wt') as f:
            f.write('\n'.join(df['Run'].dropna().astype(str).str.strip()) + '\n')
        input_arg = list_f
    info_dir = check_path(os.path.join(coll, 'info'),
                          must_exist=False, in_platform=True)
    os.makedirs(info_dir, exist_ok=True)
    # 未提供 AI Key 时自动降级 local 模式（纯规则解析，不调用 LLM）
    # key 优先级：请求体 > 环境变量 DEEPSEEK_API_KEY > local 模式。
    # 环境变量避免密钥落在本机任意进程可见的命令行参数里。
    ds_key = (str(body.get('deepseek_api') or '').strip()
              or os.environ.get('DEEPSEEK_API_KEY', '').strip())
    mode = 'both' if ds_key else 'local'
    cmd = engine_cmd(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'vp', 'public_meta', 'info_engine.py'),
                     '-i', input_arg, '-o', info_dir, '-t', '4', '-m', mode)
    if ds_key:
        cmd += ['--deepseek-api', ds_key]
    if body.get('fill_date'):
        cmd.append('--fill-date')
    return jsonify({'task': _meta_stream_task(
        f'元数据提取 {name}', cmd, f'元数据完成: {info_dir}',
        log_file=os.path.join(info_dir, 'run.log'), link='/meta')})


@bp.route('/api/meta/plot', methods=['POST'])
def api_meta_plot():
    """body: {name} → 后台任务：按项目元数据绘制 SCI 风格图
    （发布时间 / 机构 / 组织 / 地理 / 测序平台等，PNG 600dpi + PDF 矢量）。
    输入优先 Core14 元数据，缺失时回退 Full。"""
    body = request.get_json(force=True) or {}
    name = (body.get('name') or '').strip()
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name or ''):
        abort(400, '无效的检索项目名')
    coll = check_path(os.path.join(DIRS['meta_search'], name),
                      must_exist=True, in_platform=True)
    csv = None
    for cand in ('Global_Unified_Metadata_Core14.csv',
                 'Global_Unified_Metadata_Full.csv'):
        f = os.path.join(coll, 'info', cand)
        if os.path.isfile(f):
            csv = f
            break
    if not csv:
        abort(400, '该项目还没有元数据，请先「提取元数据」')
    out_dir = os.path.join(coll, 'plot')
    ts = time.strftime('%Y%m%d_%H%M%S')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp.public_meta.landscape_plot import plot_sci_landscape
        prog('plot', 0.15, f'读取 {os.path.basename(csv)} 并绘制')
        plot_sci_landscape(csv, out_dir)
        n = len([f for f in os.listdir(out_dir)
                 if f.lower().endswith(('.png', '.pdf'))])
        prog('plot', 1.0, f'完成 {n} 个图表文件')
        logger.close()
        return {'kind': 'tool', 'stats': {'项目': name, '图表文件': n}}

    tid = tm.start(cfg.tr(f'检索绘图·{name} {ts}', f'Meta plot·{name} {ts}'), job,
                   weight='light', link='/meta')
    return jsonify({'task': tid})


@bp.route('/api/meta/plot/files')
def api_meta_plot_files():
    """某检索项目的已生成图表列表（含可预览 URL）。"""
    name = request.args.get('name') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    plot_dir = check_path(os.path.join(DIRS['meta_search'], name, 'plot'),
                          must_exist=False, in_platform=True)
    files = []
    if os.path.isdir(plot_dir):
        for f in sorted(os.listdir(plot_dir)):
            if f.lower().endswith(('.png', '.pdf')):
                files.append({'name': f,
                              'url': f'/api/meta/plot/file?name={name}&file={f}',
                              'size': fmt_size(os.path.getsize(
                                  os.path.join(plot_dir, f)))})
    return jsonify({'files': files})


@bp.route('/api/meta/plot/file')
def api_meta_plot_file():
    """图表文件服务（严格限制在该项目 plot 目录内，防穿越）。"""
    name = request.args.get('name') or ''
    fname = request.args.get('file') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name) or '..' in fname:
        abort(400, '无效参数')
    plot_dir = check_path(os.path.join(DIRS['meta_search'], name, 'plot'),
                          must_exist=False, in_platform=True)
    p = check_path(os.path.join(plot_dir, fname), must_exist=True,
                   in_platform=True)
    if not os.path.isfile(p):
        abort(404)
    return send_file(p)


@bp.route('/api/meta/table/<name>')
def api_meta_table(name):
    """项目表数据（分页 + 排序 + 关键字/字段筛选 + 绝对行号供编辑）。
    q=全列关键字  qcol=限列关键字  f.<列>=值（可多个，AND 模糊匹配）"""
    which = request.args.get('which') or 'search'
    if which not in _META_TABLE_FILES:
        abort(400, 'which 需为 search/core14/full')
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name or ''):
        abort(400, '无效的检索项目名')
    df = _meta_load_df(name, which)
    if df is None:
        sub, _ = _META_TABLE_FILES[which]
        if which == 'search':
            abort(404, '还没有检索结果')
        abort(404, '还没有元数据结果（先提取元数据）')
    cols = list(df.columns)
    total = len(df)

    # 字段筛选（本地元数据筛选器的下拉组合，AND + 不区分大小写包含）
    for key, val in request.args.items():
        if key.startswith('f.'):
            col = key[2:]
            val = val.strip()
            if col in cols and val and val.lower() != '(any)':
                df = df[df[col].astype(str).str.lower()
                        .str.contains(val.lower(), regex=False)]
    # 快速关键字（全列或限单列）
    q = (request.args.get('q') or '').strip()
    if q:
        ql = q.lower()
        qcol = request.args.get('qcol') or ''
        if qcol and qcol in cols:
            df = df[df[qcol].astype(str).str.lower()
                    .str.contains(ql, regex=False)]
        else:
            mask = None
            for c in cols:
                m = df[c].astype(str).str.lower().str.contains(ql, regex=False)
                mask = m if mask is None else (mask | m)
            df = df[mask.fillna(False)] if mask is not None else df
    n_filtered = len(df)

    # 排序：数值列按数值、文本列按字典序；稳定排序保持原相对顺序
    sort_col = request.args.get('sort') or ''
    sort_dir = (request.args.get('dir') or 'asc').lower()
    if sort_col in cols:
        import pandas as pd
        num = pd.to_numeric(df[sort_col], errors='coerce')
        if num.notna().all():
            df = df.assign(_k=num).sort_values(
                '_k', ascending=sort_dir != 'desc').drop(columns='_k')
        else:
            df = df.sort_values(sort_col, ascending=sort_dir != 'desc',
                                kind='stable')
    try:
        page = max(1, int(request.args.get('page') or 1))
    except ValueError:
        page = 1
    try:
        per_page = min(max(1, int(request.args.get('per_page') or 20)), 500)
    except ValueError:
        per_page = 20
    pages = max(1, (n_filtered + per_page - 1) // per_page)
    page = min(page, pages)
    view = df.iloc[(page - 1) * per_page: page * per_page]
    rows = view.values.tolist()
    row_ids = [int(i) for i in view.index]      # 绝对行号 → 编辑定位
    run_col = cols.index('Run') if 'Run' in cols else None
    return jsonify({
        'columns': cols, 'rows': rows, 'row_ids': row_ids,
        'total': total, 'filtered': n_filtered,
        'page': page, 'per_page': per_page, 'pages': pages,
        'run_col': run_col,
        'n_missing': _meta_missing_count(_meta_load_df(name, which)),
        'download': f'/meta_files/{name}/{_META_TABLE_FILES[which][0]}/'
                    f'{_META_TABLE_FILES[which][1]}'})


@bp.route('/api/meta/unique/<name>')
def api_meta_unique(name):
    """列去重取值（本地筛选下拉，过滤占位符，最多 200 项）。"""
    which = request.args.get('which') or 'search'
    if which not in _META_TABLE_FILES:
        abort(400, 'which 需为 search/core14/full')
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name or ''):
        abort(400, '无效的检索项目名')
    df = _meta_load_df(name, which)
    if df is None:
        return jsonify({'values': {}})
    cols_req = [c for c in (request.args.get('cols') or '').split(',') if c]
    out = {}
    for col in cols_req:
        if col not in df.columns:
            continue
        s = df[col].astype(str).str.strip()
        s = s[~s.str.lower().isin(_META_NA)]
        out[col] = sorted(s.unique().tolist())[:200]
    return jsonify({'values': out})


def _meta_editable_target():
    body = request.get_json(force=True) or {}
    name = str(body.get('name') or '')
    which = str(body.get('which') or '')
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    if which not in _META_EDITABLE:
        abort(400, '该表为只读明细，仅 search / core14 可编辑')
    df = _meta_load_df(name, which)
    if df is None:
        abort(404, '表不存在（先检索/提取元数据）')
    return body, name, which, df


@bp.route('/api/meta/cell', methods=['POST'])
def api_meta_cell():
    """编辑单个单元格（双击编辑）。body: {name, which, row, col, value}"""
    body, name, which, df = _meta_editable_target()
    cols = list(df.columns)
    col = str(body.get('col') or '')
    try:
        row = int(body.get('row'))
    except (TypeError, ValueError):
        abort(400, 'row 需为整数')
    if col not in cols or not (0 <= row < len(df)):
        abort(400, '无效的单元格位置')
    value = '' if body.get('value') is None else str(body.get('value'))
    df.iat[row, cols.index(col)] = value
    _meta_store_df(name, which, df)
    return jsonify({'ok': True, 'value': value,
                    'n_missing': _meta_missing_count(df)})


@bp.route('/api/meta/row', methods=['POST'])
def api_meta_row():
    """整行多字段保存（编辑面板）。body: {name, which, row, data:{col:value}}"""
    body, name, which, df = _meta_editable_target()
    data = body.get('data') or {}
    try:
        row = int(body.get('row'))
    except (TypeError, ValueError):
        abort(400, 'row 需为整数')
    if not (0 <= row < len(df)):
        abort(400, '无效的行号')
    applied = 0
    for col, val in data.items():
        if col in df.columns:
            df.iat[row, df.columns.get_loc(col)] = \
                '' if val is None else str(val)
            applied += 1
    _meta_store_df(name, which, df)
    return jsonify({'ok': True, 'applied': applied,
                    'n_missing': _meta_missing_count(df)})


@bp.route('/api/meta/rows/delete', methods=['POST'])
def api_meta_rows_delete():
    """删除选中行。body: {name, which, rows:[绝对行号...]}"""
    body, name, which, df = _meta_editable_target()
    rows = body.get('rows') or []
    try:
        rows = sorted({int(r) for r in rows})
    except (TypeError, ValueError):
        abort(400, 'rows 需为行号数组')
    if not rows:
        abort(400, '未选择行')
    if rows[-1] >= len(df) or rows[0] < 0:
        abort(400, '行号越界')
    df = df.drop(index=rows).reset_index(drop=True)
    _meta_store_df(name, which, df)
    return jsonify({'ok': True, 'total': len(df),
                    'n_missing': _meta_missing_count(df)})


@bp.route('/api/meta/record')
def api_meta_record():
    """单行完整记录（编辑面板）+ 该 Run 的网页缓存 JSON（若有）。"""
    which = request.args.get('which') or 'core14'
    if which not in _META_TABLE_FILES:
        abort(400, 'which 需为 search/core14/full')
    name = request.args.get('name') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    df = _meta_load_df(name, which)
    if df is None:
        abort(404, '表不存在')
    try:
        row = int(request.args.get('row') or 0)
    except ValueError:
        row = 0
    if not (0 <= row < len(df)):
        abort(400, '行号越界')
    rec = {c: str(df.iat[row, df.columns.get_loc(c)]) for c in df.columns}
    cache = None
    run = rec.get('Run', '').strip()
    if run:
        for d in _meta_web_cache_dirs(name):
            jp = os.path.join(d, f'{run}.json')
            if os.path.isfile(jp):
                try:
                    with safe_open(jp) as f:
                        cache = json.load(f)
                except Exception:
                    cache = None
                break
    return jsonify({'row': rec, 'total': len(df), 'row_index': row,
                    'cache': cache})


@bp.route('/api/meta/fill_cache', methods=['POST'])
def api_meta_fill_cache():
    """从本地网页缓存 JSON 回填缺失字段（对应 GUI「Fill from Cache」）。"""
    body, name, which, df = _meta_editable_target()
    mapping = {'SubmissionDate': 'ReleaseDate', 'Organization': 'CenterName',
               'PRJ': 'BioProject', 'TaxID': 'TaxID',
               'ScientificName': 'ScientificName'}
    if 'Run' not in df.columns:
        abort(400, '表缺少 Run 列')
    cache_dirs = _meta_web_cache_dirs(name)
    filled = 0
    for idx in range(len(df)):
        run = str(df.iat[idx, df.columns.get_loc('Run')]).strip()
        if not run:
            continue
        data = None
        for d in cache_dirs:
            jp = os.path.join(d, f'{run}.json')
            if os.path.isfile(jp):
                try:
                    with safe_open(jp) as f:
                        data = json.load(f)
                except Exception:
                    data = None
                break
        if not data:
            continue
        for jk, col in mapping.items():
            if col not in df.columns:
                continue
            cur = str(df.iat[idx, df.columns.get_loc(col)]).strip()
            val = str(data.get(jk, '') or '').strip()
            if val and cur.lower() in _META_NA:
                df.iat[idx, df.columns.get_loc(col)] = val
                filled += 1
    if filled:
        _meta_store_df(name, which, df)
    return jsonify({'ok': True, 'filled': filled})


@bp.route('/api/meta/ai_fill', methods=['POST'])
def api_meta_ai_fill():
    """AI 清洗/补全缺失元数据（DeepSeek，后台任务；绿色标记由前端追踪）。"""
    body, name, which, df = _meta_editable_target()
    key = _meta_deepseek_key()
    if not key:
        abort(400, '需要 DEEPSEEK_API_KEY 环境变量（或用 CLI --deepseek-api）')
    model = str(body.get('model') or 'deepseek-v4-flash')
    cols = [c for c in AI_FILL_COLS if c in df.columns]
    if not cols:
        abort(400, '表中没有可清洗的元数据字段')
    records = []
    for i in range(len(df)):
        row = {c: str(df.iat[i, df.columns.get_loc(c)]).strip() for c in cols}
        if any(v and v.lower() not in _AI_EMPTY for v in row.values()):
            records.append((i, row))

    def job(log, prog, cancel):
        client = _meta_ai_client(key)
        total, changed, done = len(records), 0, 0
        cells = []
        log(f'AI 清洗 {total} 条记录（模型 {model}）…')
        for n, (idx, rec) in enumerate(records):
            if cancel.is_set():
                raise RuntimeError('已取消')
            try:
                text = _meta_ai_chat(client, model, AI_SANITIZE_PROMPT,
                                     json.dumps(rec, ensure_ascii=False))
                text = re.sub(r'<think>.*?</think>', '', text, flags=re.S).strip()
                text = re.sub(r'^```(?:json)?\s*|\s*```$', '', text,
                              flags=re.I).strip()
                result = json.loads(text)
                for col in cols:
                    if col not in result:
                        continue
                    old = rec.get(col, '').strip()
                    new = str(result[col]).strip()
                    if not new or old == new or new.lower() in _AI_EMPTY:
                        continue
                    if col == 'Location' and all(
                            'unknown' in p.lower() for p in new.split(',')):
                        continue
                    df.iat[idx, df.columns.get_loc(col)] = new
                    changed += 1
                    cells.append([idx, col])
            except json.JSONDecodeError as e:
                log(f'第 {idx + 1} 行：AI 返回非 JSON：{str(e)[:120]}')
            except Exception as e:
                log(f'第 {idx + 1} 行：{str(e)[:200]}')
            done += 1
            prog('ai', done / max(total, 1),
                 f'AI 清洗 {done}/{total}（已改 {changed} 字段）')
            if done % 10 == 0:
                log(f'进度 {done}/{total}，已修改 {changed} 个字段')
        if changed:
            _meta_store_df(name, which, df)
        return {'kind': 'meta', 'stats': {'记录': total, '修改字段': changed,
                                          'cells': cells[:5000]}}

    ts = time.strftime('%H:%M:%S')
    tid = tm.start(cfg.tr(f'AI 元数据补全·{name} {ts}',
                          f'AI metadata fill·{name} {ts}'), job,
                   weight='light', link='/meta')
    return jsonify({'task': tid, 'records': len(records)})


@bp.route('/api/meta/ai_summary', methods=['POST'])
def api_meta_ai_summary():
    """基于表统计生成 SCI 数据描述段落（DeepSeek，同步返回）。"""
    body = request.get_json(force=True) or {}
    name = str(body.get('name') or '')
    which = str(body.get('which') or 'core14')
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    if which not in _META_TABLE_FILES:
        abort(400, 'which 需为 search/core14/full')
    df = _meta_load_df(name, which)
    if df is None or df.empty:
        abort(404, '没有数据')
    key = (str(body.get('deepseek_api') or '').strip()
           or _meta_deepseek_key())
    if not key:
        abort(400, '需要 DEEPSEEK_API_KEY 环境变量')

    def counter(col):
        if col not in df.columns:
            from collections import Counter
            return Counter()
        s = df[col].astype(str).str.strip()
        s = s[~s.str.lower().isin(_META_NA)]
        from collections import Counter
        return Counter(s)

    def top(c, n=10):
        return ', '.join(f'{k}({v})' for k, v in c.most_common(n))

    years = set()
    date_col = ('CollectionDate' if 'CollectionDate' in df.columns
                else 'ReleaseDate' if 'ReleaseDate' in df.columns else None)
    if date_col:
        for v in df[date_col].astype(str):
            p = v.replace('/', '-').split('-')
            if p and p[0].isdigit() and len(p[0]) == 4:
                years.add(p[0])
    total_gb, gb_n = 0.0, 0
    size_col = ('FileSize_GB' if 'FileSize_GB' in df.columns
                else 'FileSize_MB' if 'FileSize_MB' in df.columns else None)
    if size_col:
        for v in df[size_col]:
            try:
                x = float(str(v).strip())
                total_gb += x if size_col.endswith('_GB') else x / 1024
                gb_n += 1
            except (TypeError, ValueError):
                pass
    vol = (f'Data volume: {total_gb:.1f} GB total '
           f'({gb_n} runs with size data, avg {total_gb / gb_n:.1f} GB/run)'
           if gb_n else '')
    db_counts = counter('Database')
    stats_text = f"""Total records: {len(df)}
{vol}
Database sources: {top(db_counts) or 'N/A'}
Species: {top(counter('ScientificName'), 8) or 'N/A'}
Tissues: {top(counter('Tissue'), 10) or 'N/A'}
Source types: {top(counter('Source'), 8) or 'N/A'}
Locations: {top(counter('Location'), 10) or 'N/A'}
Institutions: {top(counter('CenterName'), 10) or 'N/A'}
Growth stages: {top(counter('Age_GrowthStage'), 5) or 'N/A'}
Collection years: {', '.join(sorted(years)) or 'N/A'}"""

    prompt = f"""You are a scientific writer preparing a manuscript for a virome/metagenomics study.
Based on the following metadata statistics of public sequencing runs collected for analysis, write a concise paragraph (150-250 words) suitable for the "Data Collection" or "Sample Information" section of a scientific paper.

{stats_text}

Requirements:
- Write in formal scientific English, past tense
- Include total number of runs, database split (SRA/GSA), and total data volume (in GB) if available
- Include species covered, tissue types, geographic locations, and collection time span
- Mention key institutions that contributed data
- Note the sequencing type (transcriptomic/genomic) and average data volume per run when available
- End with a note on data availability

Output ONLY the paragraph, no markdown, no headings."""
    try:
        model = str(body.get('model') or 'deepseek-v4-flash')
        text = _meta_ai_chat(_meta_ai_client(key), model,
                             'You are a scientific writer. '
                             'Output only the requested paragraph.', prompt)
    except Exception as e:
        abort(500, f'AI 总结失败: {e}')
    return jsonify({'summary': text.strip()})


@bp.route('/api/meta/import', methods=['POST'])
def api_meta_import():
    """导入 TSV/CSV/Excel 替换当前项目表（对应 GUI File > Import）。"""
    name = request.form.get('name') or ''
    which = request.form.get('which') or 'search'
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    if which not in _META_EDITABLE:
        abort(400, '仅 search / core14 表可导入')
    f = request.files.get('file')
    if not f or not f.filename:
        abort(400, '需要上传文件')
    import tempfile
    import pandas as pd
    suffix = os.path.splitext(f.filename)[1].lower()
    if suffix not in ('.csv', '.tsv', '.txt', '.xlsx', '.xls'):
        abort(400, '支持 .csv / .tsv / .txt / .xlsx / .xls')
    fd, tmp = tempfile.mkstemp(suffix=suffix)
    try:
        f.save(tmp)
        if suffix in ('.xlsx', '.xls'):
            df = pd.read_excel(tmp, dtype=str, keep_default_na=False)
        else:
            sep = '\t' if suffix in ('.tsv', '.txt') else ','
            try:
                df = pd.read_csv(tmp, sep=sep, dtype=str,
                                 keep_default_na=False, encoding='utf-8-sig')
            except UnicodeDecodeError:
                df = pd.read_csv(tmp, sep=sep, dtype=str,
                                 keep_default_na=False, encoding='gbk')
    except Exception as e:
        abort(400, f'解析失败: {e}')
    finally:
        try:
            os.close(fd)
            os.remove(tmp)
        except OSError:
            pass
    df = df.loc[:, ~df.columns.duplicated()].fillna('')
    _meta_store_df(name, which, df.astype(object))
    return jsonify({'ok': True, 'rows': len(df),
                    'columns': list(df.columns)})


@bp.route('/api/meta/export')
def api_meta_export():
    """导出项目表（csv 原样 / tsv / xlsx）。"""
    which = request.args.get('which') or 'core14'
    fmt = (request.args.get('fmt') or 'csv').lower()
    name = request.args.get('name') or ''
    if which not in _META_TABLE_FILES or fmt not in ('csv', 'tsv', 'xlsx'):
        abort(400, '无效参数')
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    sub, fname = _META_TABLE_FILES[which]
    if fmt in ('csv', 'tsv'):
        ext = '.csv' if fmt == 'csv' else '.tsv'
        p = check_path(os.path.join(_meta_proj(name), sub,
                                    os.path.splitext(fname)[0] + ext),
                       must_exist=True, in_platform=True)
        return send_file(p, as_attachment=True,
                         download_name=f'{name}_{which}{ext}')
    import io
    import pandas as pd
    df = _meta_load_df(name, which)
    if df is None:
        abort(404, '表不存在')
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine='openpyxl') as w:
        df.to_excel(w, index=False)
    buf.seek(0)
    return send_file(buf, as_attachment=True,
                     download_name=f'{name}_{which}.xlsx',
                     mimetype='application/vnd.openxmlformats-officedocument'
                              '.spreadsheetml.sheet')


@bp.route('/api/meta/overview')
def api_meta_overview():
    """图表数据包：各列 Top 取值 / 缺失率 / 年份分布 / 汇总统计。"""
    which = request.args.get('which') or 'core14'
    if which not in _META_TABLE_FILES:
        abort(400, 'which 需为 search/core14/full')
    name = request.args.get('name') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '无效的检索项目名')
    df = _meta_load_df(name, which)
    if df is None or df.empty:
        abort(404, '表不存在或为空')
    fields = {}
    for c in df.columns:
        s = df[c].astype(str).str.strip()
        miss = int(s.str.lower().isin(_META_NA).sum())
        s2 = s[~s.str.lower().isin(_META_NA)]
        vc = s2.value_counts().head(50)
        fields[c] = {'counts': [[k, int(v)] for k, v in vc.items()],
                     'missing': miss,
                     'unique': int(s2.nunique())}
    date_col = ('CollectionDate' if 'CollectionDate' in df.columns
                and df['CollectionDate'].astype(str).str.strip()
                .str.lower().isin(_META_NA).sum() < len(df)
                else 'ReleaseDate' if 'ReleaseDate' in df.columns else None)
    years = {}
    if date_col:
        for v in df[date_col].astype(str):
            p = v.replace('/', '-').split('-')
            if p and p[0].isdigit() and len(p[0]) == 4:
                years[p[0]] = years.get(p[0], 0) + 1
    total, ncol = len(df), len(df.columns)
    filled = total * ncol - _meta_missing_count(df)
    complete = 0
    for _, row in df.iterrows():
        if not any(str(v).strip().lower() in _META_NA for v in row):
            complete += 1
    total_gb, gb_n = 0.0, 0
    size_col = ('FileSize_GB' if 'FileSize_GB' in df.columns
                else 'FileSize_MB' if 'FileSize_MB' in df.columns else None)
    if size_col:
        for v in df[size_col]:
            try:
                x = float(str(v).strip())
                total_gb += x if size_col.endswith('_GB') else x / 1024
                gb_n += 1
            except (TypeError, ValueError):
                pass
    return jsonify({
        'columns': list(df.columns), 'fields': fields,
        'years': sorted(years.items()),
        'date_col': date_col,
        'summary': {'total': total, 'cols': ncol, 'filled': filled,
                    'missing': total * ncol - filled, 'complete': complete,
                    'total_gb': round(total_gb, 1) if gb_n else None,
                    'gb_runs': gb_n}})


@bp.route('/api/meta/collection/<name>/delete', methods=['POST'])
def api_meta_collection_delete(name):
    import shutil
    try:
        coll = check_path(os.path.join(DIRS['meta_search'], name),
                          must_exist=True, in_platform=True)
    except (ValueError, FileNotFoundError):
        abort(404, f'检索项目不存在: {name}')
    shutil.rmtree(coll)
    with _meta_cache_lock:
        for key in [k for k in _meta_df_cache if k[0] == name]:
            _meta_df_cache.pop(key, None)
    return jsonify({'ok': True})


@bp.route('/meta_files/<name>/<path:subpath>')
def meta_files_download(name, subpath):
    """meta 产物下载（限 meta_search/<name>/ 目录内）。"""
    try:
        coll = check_path(os.path.join(DIRS['meta_search'], name),
                          must_exist=True, in_platform=True)
        p = check_path(os.path.join(coll, subpath), must_exist=True,
                       in_platform=True)
        if not p.startswith(coll + os.sep):
            abort(404)
    except (ValueError, FileNotFoundError):
        abort(404)
    return send_file(p, as_attachment=True)

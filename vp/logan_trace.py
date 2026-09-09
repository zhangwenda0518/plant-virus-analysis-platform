# -*- coding: utf-8 -*-
"""
LOGAN 溯源模块（独立功能模块，与「公共病毒组」同级）。

把本平台鉴定的病毒 contig（或任意粘贴序列）提交到 Logan-Search
（logan-search.org，IndexThePlanet 计划：对整个 NCBI SRA 全量组装后
建立的 k-mer 索引，覆盖 ~2340 万公开样本），回答「这条病毒序列在哪些
物种和样本中真实存在」。

工作流（两种方式，平台自身不向外部发任何请求）：
  方式一 · 手动半自动：①平台切片 → ②用户复制序列到 logan-search.org 提交
           → ③导入官方结果表 → ④聚合出报告
  方式二 · 批量自动（vp/logan_submit.py）：Selenium 驱动本机 Edge/Chrome
           批量提交全部片段、自动轮询并下载结果表、回收导入（支持邮箱池
           轮换/断点续跑/补漏轮），一键完成 ②③ 两步

  ① 平台把病毒 contig 切成 ≤2.5kb 查询片段（Logan-Search 单条上限 2.5kb）
  ② 提交 Logan-Search（Groups 建议 All_No_viral_human；手动可填邮箱拿下载链接）
  ③ 结果表导入对应片段
  ④ 平台聚合出「物种分布 / 样本清单 / ANI-kmer 散点」并生成溯源报告，
     若查询来自平台样品还会与 ④ICTV 宿主预测交叉对照

查询任务存放：<平台根>/logan/<查询名>/
  meta.json            任务信息（来源样品、contigs、状态）
  segments.json        片段清单（坐标/文件/导入状态）
  query_s<i>.fasta     每片段单条 FASTA（复制提交用）
  query_all.fasta      全部片段合并（备份）
  result_s<i>.*        用户导入的原始结果文件
  result_s<i>_parsed.json  归一化后的结果行
  trace_report.html    溯源报告（离线 plotly）

引用要求：Chikhi et al. 2025, bioRxiv 10.1101/2024.07.30.605881
"""
import os
import re
import sys
import json
import time

from .config import DIRS
from .utils import (check_path, safe_open, iter_fasta, write_fasta_record)

LOGAN_SEARCH_URL = 'https://logan-search.org/dashboard'
LOGAN_CITE = ('Chikhi et al. 2025, Logan: Planetary-Scale Genome Assembly '
              'Surveys Life\'s Diversity, bioRxiv 10.1101/2024.07.30.605881')
SEG_MAX_BP = 2500          # Logan-Search 单条查询序列上限
MAX_UPLOAD_BYTES = 64 * 1024 * 1024
GROUP_HINTS = [
    ('Fast_No_human', '快速子集且排除人源（批量模式默认，最快）'),
    ('Fast', '快速子集（预检，最快）'),
    ('Fast_No_RefSeq', '快速子集，排除参考基因组'),
    ('All_No_viral_human', '全部样本但排除人源病毒（推荐，最全）'),
    ('All', '全部 2340 万样本（最全，最慢）'),
    ('Metagenomic', '仅宏基因组样本'),
    ('Metatranscriptomic', '仅宏转录组样本'),
    ('Transcriptomic', '仅转录组样本'),
    ('GenBank_RefSeq', '仅 ~4.5 万参考基因组（对照）'),
]
# logan_submit.py 支持的全部 Groups（提交校验用）
BATCH_GROUPS = [g for g, _ in GROUP_HINTS]


def logan_root():
    d = check_path(DIRS['logan'], must_exist=False, in_platform=True)
    os.makedirs(d, exist_ok=True)
    return d


def safe_job_name(name):
    n = re.sub(r'[^A-Za-z0-9_\-.]', '_', str(name or '')).strip('._')
    return n or f"query_{time.strftime('%Y%m%d_%H%M%S')}"


def _job_dir(name, must_exist=True):
    return check_path(os.path.join(logan_root(), safe_job_name(name)),
                      must_exist=must_exist, in_platform=True)


def _read_json(path):
    with safe_open(path) as f:
        return json.load(f)


# ------------------------------------------------------------------
# 查询来源
# ------------------------------------------------------------------
def list_query_samples():
    """有 ③组装 病毒 contigs 产出的样品列表（溯源查询来源）。"""
    res = check_path(DIRS['results'], must_exist=False, in_platform=True)
    out = []
    if not os.path.isdir(res):
        return out
    for name in sorted(os.listdir(res)):
        tsv = os.path.join(res, name, '03_assembly', 'virus_contigs.tsv')
        if os.path.isfile(tsv):
            try:
                out.append({'name': name, 'n_contigs': len(list_virus_contigs(name))})
            except Exception:
                continue
    return out


def list_virus_contigs(sample):
    """解析 03_assembly/virus_contigs.tsv（保持文件顺序=长度优先）。"""
    tsv = check_path(os.path.join(DIRS['results'], safe_name(sample),
                                  '03_assembly', 'virus_contigs.tsv'),
                     must_exist=True, in_platform=True)
    rows = []
    with safe_open(tsv) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {k: i for i, k in enumerate(header)}
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < len(header):
                continue
            g = lambda k: p[idx[k]] if k in idx and idx[k] < len(p) else ''
            try:
                length = int(float(g('length')))
            except ValueError:
                length = 0
            rows.append({'contig': g('contig'), 'length': length,
                         'species': g('kunpeng_species') or g('blast_species'),
                         'family': g('blast_family'),
                         'identity': g('blast_identity(%)')})
    return rows


def safe_name(sample):
    return re.sub(r'[^A-Za-z0-9_\-.]', '_', str(sample))


def _iter_pasted(pasted):
    """内存解析粘贴内容（FASTA 或纯序列文本）。yield (id, seq)。"""
    lines = [ln.strip() for ln in str(pasted).splitlines() if ln.strip()]
    if not any(ln.startswith('>') for ln in lines):
        lines = ['>query'] + lines          # 纯序列 → 补一个头
    header, chunks = None, []
    for ln in lines:
        if ln.startswith('>'):
            if header is not None and chunks:
                yield header, ''.join(chunks)
            header, chunks = ln[1:].strip() or 'query', []
        else:
            chunks.append(ln)
    if header is not None and chunks:
        yield header, ''.join(chunks)


def _load_source_fasta(sample=None, contig_ids=None, pasted=None):
    """返回 [(contig_id, seq)]。sample+contig_ids 或 pasted 二选一。"""
    recs = []
    if pasted:
        for header, seq in _iter_pasted(pasted):
            seq = re.sub(r'[^ACGTUNacgtun]', '', seq.upper())
            if seq:
                cid = re.sub(r'[^A-Za-z0-9_\-.]', '_',
                             header.split()[0]) if header.split() else 'query'
                recs.append((cid or 'query', seq))
        if not recs:
            raise ValueError('粘贴内容未解析出序列（请提供 FASTA 或纯序列）')
        return recs
    fa = check_path(os.path.join(DIRS['results'], safe_name(sample),
                                 '03_assembly', 'viral_contigs.fasta'),
                    must_exist=True, in_platform=True)
    wanted = set(contig_ids or [])
    for header, seq in iter_fasta(fa):
        cid = header.split()[0]
        if wanted and cid not in wanted:
            continue
        recs.append((cid, seq.upper()))
    if wanted and not recs:
        raise ValueError('所选 contigs 在 viral_contigs.fasta 中不存在')
    return recs


def make_segments(seqlen, n_seg):
    """把长 seqlen 的序列切成 ≤2500bp 的 n_seg 个窗口（返回 0-based 开闭区间）。

    单窗口取正中；多窗口均匀铺开，覆盖不同区域提高召回。
    """
    if seqlen <= SEG_MAX_BP:
        return [(0, seqlen)]
    n_seg = max(1, min(int(n_seg or 1), 4))
    seg = SEG_MAX_BP
    if n_seg == 1:
        s = max(0, (seqlen - seg) // 2)
        return [(s, s + seg)]
    step = (seqlen - seg) / (n_seg - 1)
    out = []
    for i in range(n_seg):
        s = int(round(i * step))
        s = min(max(0, s), seqlen - seg)
        out.append((s, s + seg))
    # 去重（极短序列多窗口可能重叠成同区间）
    seen, uniq = set(), []
    for s, e in out:
        if (s, e) not in seen:
            seen.add((s, e))
            uniq.append((s, e))
    return uniq


# ------------------------------------------------------------------
# 任务生命周期
# ------------------------------------------------------------------
def create_job(name, sample=None, contig_ids=None, pasted=None, n_seg=2):
    """创建溯源查询任务：切片 → 写 query fasta → meta/segments。"""
    jdir = _job_dir(name, must_exist=False)
    if os.path.isdir(jdir) and any(os.scandir(jdir)):
        raise ValueError(f'查询名 {safe_job_name(name)} 已存在，请换一个')
    os.makedirs(jdir, exist_ok=True)

    recs = _load_source_fasta(sample=sample, contig_ids=contig_ids,
                              pasted=pasted)
    if sample:
        contig_rows = {r['contig']: r for r in list_virus_contigs(sample)}
    else:
        contig_rows = {}

    segments, contig_infos = [], []
    for cid, seq in recs:
        wins = make_segments(len(seq), n_seg)
        info = {'contig': cid, 'length': len(seq)}
        if cid in contig_rows:
            info.update({'species': contig_rows[cid].get('species', ''),
                         'family': contig_rows[cid].get('family', '')})
        contig_infos.append(info)
        for i, (s, e) in enumerate(wins):
            idx = len(segments) + 1
            header = f"{safe_job_name(name)}|{cid}_s{idx}_{s + 1}-{e}"
            seg_fa = check_path(os.path.join(jdir, f'query_s{idx}.fasta'),
                                must_exist=False, in_platform=True)
            with safe_open(seg_fa, 'wt') as f:
                write_fasta_record(f, header, seq[s:e])
            segments.append({
                'index': idx, 'contig': cid,
                'start': s + 1, 'end': e, 'length': e - s,
                'header': header, 'file': os.path.basename(seg_fa),
                'result_file': '', 'imported': False,
                'imported_at': '', 'n_hits': 0,
            })

    meta = {
        'name': safe_job_name(name),
        'sample': safe_name(sample) if sample else '',
        'created': time.strftime('%Y-%m-%d %H:%M:%S'),
        'contigs': contig_infos,
        'n_segments': len(segments),
        'submit_url': LOGAN_SEARCH_URL,
        'group_hint': 'All_No_viral_human',
        'status': 'pending',
        'cite': LOGAN_CITE,
    }
    with safe_open(os.path.join(jdir, 'meta.json'), 'wt') as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    with safe_open(os.path.join(jdir, 'segments.json'), 'wt') as f:
        json.dump(segments, f, ensure_ascii=False, indent=2)

    all_fa = check_path(os.path.join(jdir, 'query_all.fasta'),
                        must_exist=False, in_platform=True)
    with safe_open(all_fa, 'wt') as f:
        for s in segments:
            with safe_open(os.path.join(jdir, s['file'])) as g:
                f.write(g.read())
    return job_detail(meta['name'])


def _save_json(jdir, fname, data):
    with safe_open(os.path.join(jdir, fname), 'wt') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def _load_job(name):
    jdir = _job_dir(name)
    meta = _read_json(os.path.join(jdir, 'meta.json'))
    segs = _read_json(os.path.join(jdir, 'segments.json'))
    return jdir, meta, segs


def job_detail(name):
    jdir, meta, segs = _load_job(name)
    report = os.path.join(jdir, 'trace_report.html')
    meta = dict(meta)
    meta['segments'] = segs
    meta['has_report'] = os.path.isfile(report)
    meta['n_imported'] = sum(1 for s in segs if s.get('imported'))
    return meta


def list_jobs():
    root = logan_root()
    out = []
    for d in sorted(os.listdir(root), reverse=True):
        mp = os.path.join(root, d, 'meta.json')
        if os.path.isfile(mp):
            try:
                out.append(job_detail(d))
            except Exception:
                continue
    return out


def delete_job(name):
    import shutil
    jdir = _job_dir(name)
    shutil.rmtree(jdir, ignore_errors=True)


def get_segment_fasta(name, index):
    """返回 (header, seq)。"""
    jdir, _meta, segs = _load_job(name)
    seg = next((s for s in segs if s['index'] == int(index)), None)
    if not seg:
        raise ValueError(f'片段 s{index} 不存在')
    seq = []
    with safe_open(os.path.join(jdir, seg['file'])) as f:
        for line in f:
            if not line.startswith('>'):
                seq.append(line.strip())
    return seg['header'], ''.join(seq)


# ------------------------------------------------------------------
# 结果导入与解析
# ------------------------------------------------------------------
# 结果表列名 → 平台字段（归一化后精确匹配，aliases 按优先级）
_COLMAP = {
    'acc': ('acc', 'runaccession', 'run', 'accession', 'runid'),
    'sample_acc': ('sampleacc', 'sampleaccession', 'biosampleacc'),
    'organism': ('organism', 'scientificname', 'organismname', 'speciesname'),
    'biosample': ('biosample', 'biosampleid'),
    'bioproject': ('bioproject', 'bioprojectid'),
    'experiment': ('experiment', 'experimentaccession'),
    'study': ('srastudy', 'srastudyaccession', 'study', 'studyaccession'),
    'sample_name': ('samplename', 'libraryname'),
    'kmer_cov': ('kmercoverage', 'kmercov', 'kmerscoverage', 'coverage'),
    'ani': ('aniestimation', 'ani', 'aniestimate'),
    'pvalue': ('pvalue', 'pval'),
    'evalue': ('evalue', 'eval'),
    'assay_type': ('assaytype', 'librarystrategy'),
    'instrument': ('instrument', 'instrumentmodel'),
    'platform': ('platform',),
    'location': ('location', 'country', 'geographiclocation', 'geo'),
    'lat': ('lat', 'latitude'),
    'lon': ('lon', 'long', 'longitude'),
    'disease': ('dolabel', 'disease', 'diseaselabel'),
    'tissue': ('btolabel', 'tissue', 'tissuelabel'),
}
_NUM_FIELDS = ('kmer_cov', 'ani', 'pvalue', 'evalue', 'lat', 'lon')


def _norm_header(h):
    return re.sub(r'[^a-z0-9]', '', str(h).lower())


def _sniff_delimiter(line):
    cand = {'\t': line.count('\t'), ',': line.count(','), ';': line.count(';')}
    return max(cand, key=cand.get) if max(cand.values()) > 0 else '\t'


def parse_result_table(raw):
    """解析 Logan-Search 结果表（CSV/TSV，表头驱动，宽容列名差异）。

    返回 (rows, mapped_cols)。行字段缺失补空串；数值列尽力 float。
    """
    if raw[:2] == b'\x1f\x8b':           # gzip
        import gzip
        raw = gzip.decompress(raw)
    text = raw.decode('utf-8-sig', errors='replace')
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise ValueError('结果文件为空或没有数据行')
    delim = _sniff_delimiter(lines[0])
    header = [h.strip().lstrip('\ufeff') for h in lines[0].split(delim)]
    norm = [_norm_header(h) for h in header]

    col2field, field2col = {}, {}
    for field, aliases in _COLMAP.items():
        for alias in aliases:            # 优先级顺序
            if alias in norm and norm.index(alias) not in col2field:
                col2field[norm.index(alias)] = field
                field2col[field] = norm.index(alias)
                break
    if 'organism' not in field2col and 'acc' not in field2col:
        raise ValueError('无法识别结果表列名（需至少包含 organism 或 '
                         'run accession 列）。请确认导入的是 Logan-Search '
                         '的结果表而非网页截图。')

    rows = []
    for ln in lines[1:]:
        p = ln.split(delim)
        row = {k: '' for k in _COLMAP}
        for ci, field in col2field.items():
            v = p[ci].strip() if ci < len(p) else ''
            if field in _NUM_FIELDS:
                try:
                    v = float(v) if v and v.upper() not in ('NA', 'NAN', '-') else None
                except ValueError:
                    v = None
            row[field] = v
        rows.append(row)
    return rows, sorted(field2col.keys())


def import_result(name, seg_index, filename, raw):
    """导入某片段的 Logan-Search 结果文件 → 解析入库 → 刷新任务状态。"""
    if len(raw) > MAX_UPLOAD_BYTES:
        raise ValueError('结果文件超过 64MB 上限')
    jdir, meta, segs = _load_job(name)
    seg = next((s for s in segs if s['index'] == int(seg_index)), None)
    if not seg:
        raise ValueError(f'片段 s{seg_index} 不存在')
    rows, cols = parse_result_table(raw)

    ext = os.path.splitext(filename or '')[1].lower() or '.tsv'
    raw_name = f"result_s{seg['index']}{ext if ext in ('.tsv', '.csv', '.txt') else '.tsv'}"
    with safe_open(os.path.join(jdir, raw_name), 'wb') as f:
        f.write(raw)
    seg.update({'imported': True,
                'imported_at': time.strftime('%Y-%m-%d %H:%M:%S'),
                'result_file': raw_name, 'n_hits': len(rows)})
    _save_json(jdir, 'segments.json', segs)
    _save_json(jdir, f"result_s{seg['index']}_parsed.json",
               {'cols': cols, 'n': len(rows), 'rows': rows})
    meta['status'] = ('done' if all(s.get('imported') for s in segs)
                      else 'partial')
    _save_json(jdir, 'meta.json', meta)
    build_report(meta['name'])
    return job_detail(meta['name'])


# ------------------------------------------------------------------
# 批量自动提交（方式二：vp/logan_submit.py，Selenium 驱动本机 Edge/Chrome）
# ------------------------------------------------------------------
BATCH_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'logan_submit.py')


def selenium_ready():
    import importlib.util
    return importlib.util.find_spec('selenium') is not None


def batch_acc(job_name, seg_index):
    """批量输入用 ID（= 结果文件名 stem）：{任务名}_s{片段号}。"""
    return f"{safe_job_name(job_name)}_s{int(seg_index)}"


def batch_prepare(name):
    """生成批量输入 FASTA（logan_submit.py 输入格式）与片段映射。

    返回 (in_path, out_dir, acc2seg, n_pending)。
    """
    jdir, meta, segs = _load_job(name)
    pending = [s for s in segs if not s.get('imported')]
    if not pending:
        raise ValueError('该任务所有片段均已导入结果，无需批量提交')
    if not getattr(sys, 'frozen', False) and not os.path.isfile(BATCH_SCRIPT):
        raise RuntimeError('缺少 vp/logan_submit.py（批量提交脚本）')
    if not selenium_ready():
        raise ValueError('未安装 selenium：python -m pip install selenium')
    in_path = check_path(os.path.join(jdir, 'batch_input.fa'),
                         must_exist=False, in_platform=True)
    out_dir = check_path(os.path.join(jdir, 'batch_out'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    acc2seg = {}
    with safe_open(in_path, 'wt') as f:
        for s in pending:
            acc = batch_acc(meta['name'], s['index'])
            _h, seq = get_segment_fasta(meta['name'], s['index'])
            write_fasta_record(f, acc, seq)
            acc2seg[acc] = s['index']
    return in_path, out_dir, acc2seg, len(pending)


def batch_collect(name, acc2seg=None):
    """回收 batch_out/ 下已下载的 {acc}.tsv → 导入对应片段（可重复调用）。"""
    jdir, meta, segs = _load_job(name)
    out_dir = check_path(os.path.join(jdir, 'batch_out'),
                         must_exist=True, in_platform=True)
    if acc2seg is None:
        acc2seg = {batch_acc(meta['name'], s['index']): s['index']
                   for s in segs if not s.get('imported')}
    n = 0
    for acc, idx in acc2seg.items():
        p = os.path.join(out_dir, f"{safe_stem_name(acc)}.tsv")
        if not os.path.isfile(p):
            continue
        with safe_open(p, 'rb') as f:
            raw = f.read()
        if raw.strip():
            import_result(meta['name'], idx, f"{acc}.tsv", raw)
            n += 1
    return n


def safe_stem_name(acc):
    """与 vp/logan_submit.safe_stem 同口径的文件名安全化。"""
    s = re.sub(r"[^A-Za-z0-9._\-]", "_", str(acc))
    s = re.sub(r"\.{2,}", ".", s).strip(".")
    return s or "query"


def batch_submit(name, emails, group='Fast_No_human', headless=True,
                 first_wait=300, max_wait=1800,
                 logger=None, progress=None, cancel=None):
    """一键批量：生成输入 → 子进程跑 logan_submit.py（提交+轮询+下载）
    → 回收下载好的结果表并导入片段、生成报告。

    logan_submit.py 自带断点续跑（已下载的 {acc}.tsv 跳过），
    挂起的片段稍后重跑批量即可补漏。
    """
    import subprocess
    if not emails:
        raise ValueError('至少填写一个通知邮箱')
    if group not in BATCH_GROUPS:
        raise ValueError(f'Groups 非法: {group}（可选: {", ".join(BATCH_GROUPS)}）')
    in_path, out_dir, acc2seg, n_pending = batch_prepare(name)
    if logger:
        logger.log(f'批量输入已生成: {n_pending} 个片段 '
                   f'（{", ".join(sorted(acc2seg))}）')
    from vp.config import engine_cmd
    cmd = engine_cmd(BATCH_SCRIPT, str(in_path), '-o', str(out_dir),
           '--email', ','.join(emails), '--group', group,
           '--first-wait', str(int(first_wait)), '--max-wait', str(int(max_wait)))
    if headless:
        cmd.append('--headless')
    flags = 0x08000000 if sys.platform == 'win32' else 0   # CREATE_NO_WINDOW
    if logger:
        logger.log(f'启动批量提交: Groups={group}, headless={headless}, '
                   f'邮箱 {len(emails)} 个')
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True,
                            encoding='utf-8', errors='replace',
                            creationflags=flags)
    code = None

    def _human_progress(ev):
        """PROGRESS 事件 → (pct, 中文进度文案)。

        进度条按"已完成条数"驱动（单调不回退）；当前条的细粒度状态
        （等待中/已等秒数/HTTP 状态/下载中）放在文案里实时刷新。"""
        done = int(ev.get('done') or 0)
        total = max(int(ev.get('total') or 1), 1)
        stage = ev.get('stage', '')
        acc = ev.get('acc') or ''
        pct = min(0.98, 0.02 + 0.95 * done / total)
        sid = (ev.get('sid') or '')[:13]
        if stage == 'submit':
            txt = f'正在提交 {acc}（{done}/{total}）'
        elif stage == 'submitted':
            txt = f'已提交 {acc}' + (f'（session {sid}…）' if sid else f'：{ev.get("msg", "")}')
        elif stage == 'waiting':
            txt = (f'{acc} 等待服务器结果：已等 {ev.get("waited")}s / '
                   f'约 {int((ev.get("span") or 2100) / 60)} 分钟'
                   f'（HTTP {ev.get("http")}），就绪即自动下载')
        elif stage == 'downloading':
            txt = f'{acc} 结果就绪，下载中'
        elif stage in ('downloaded', 'segment_done'):
            txt = f'已完成 {done}/{total} 条'
        else:
            txt = str(ev.get('msg') or stage)
        return pct, txt

    try:
        for line in proc.stdout:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('PROGRESS '):
                if progress:
                    try:
                        ev = json.loads(line[9:])
                    except ValueError:
                        ev = {}
                    pct, txt = _human_progress(ev)
                    progress('batch', pct, txt)
                continue
            if logger:
                logger.log(line)
            if progress:
                m = re.search(r'\[(\d+)/(\d+)\]', line)
                if m:
                    done, total = int(m.group(1)), int(m.group(2))
                    progress('batch', min(0.95, 0.05 + 0.9 * done / max(total, 1)),
                             f'批量提交 {done}/{total}')
            if cancel is not None and cancel.is_set():
                proc.terminate()
                if logger:
                    logger.log('收到取消请求，已终止批量提交进程')
                break
        code = proc.wait(timeout=60)
    finally:
        if proc.poll() is None:
            proc.kill()
    if cancel is not None and cancel.is_set():
        raise RuntimeError('批量提交已取消')
    n = batch_collect(name, acc2seg)
    if logger:
        logger.log(f'批量结束（退出码 {code}），回收导入 {n}/{n_pending} 个片段；'
                   f'未落定的片段可再次批量提交补漏')
    if progress:
        progress('batch', 1.0, f'批量完成，回收 {n}/{n_pending} 个片段')
    return {'n_pending': n_pending, 'n_imported': n, 'returncode': code}


# ------------------------------------------------------------------
# 聚合
# ------------------------------------------------------------------
def _load_rows(jdir, seg):
    p = os.path.join(jdir, f"result_s{seg['index']}_parsed.json")
    if not os.path.isfile(p):
        return []
    return _read_json(p).get('rows', [])


def _f(v, default=None):
    return v if isinstance(v, (int, float)) else default


def aggregate(name):
    """跨片段聚合：union 样本（同 accession 去重，保留最优行）、物种计数等。"""
    jdir, meta, segs = _load_job(name)
    per_seg, union = {}, {}
    for seg in segs:
        rows = _load_rows(jdir, seg)
        per_seg[seg['index']] = rows
        for r in rows:
            key = r.get('acc') or f"{seg['index']}:{len(union)}"
            cur = union.get(key)
            if cur is None or (_f(r.get('kmer_cov'), 0) > _f(cur.get('kmer_cov'), 0)):
                union[key] = dict(r, _seg=seg['index'])
    organisms, assays, platforms = {}, {}, {}
    for r in union.values():
        org = r.get('organism') or '（未知物种）'
        organisms[org] = organisms.get(org, 0) + 1
        if r.get('assay_type'):
            assays[r['assay_type']] = assays.get(r['assay_type'], 0) + 1
        if r.get('platform'):
            platforms[r['platform']] = platforms.get(r['platform'], 0) + 1
    return {
        'meta': meta, 'segments': segs, 'per_seg': per_seg,
        'union_rows': list(union.values()),
        'n_union': len(union),
        'organisms': dict(sorted(organisms.items(), key=lambda x: -x[1])),
        'assays': dict(sorted(assays.items(), key=lambda x: -x[1])),
        'platforms': dict(sorted(platforms.items(), key=lambda x: -x[1])),
    }


def _host_cross(meta, contig_ids):
    """查询来自平台样品时，取 ④宿主预测 对应行做交叉对照。"""
    if not meta.get('sample'):
        return []
    p = os.path.join(DIRS['results'], meta['sample'], '08_host_analysis',
                     'host_prediction.tsv')
    if not os.path.isfile(p):
        return []
    out = []
    with safe_open(check_path(p, must_exist=True, in_platform=True)) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {k: i for i, k in enumerate(header)}
        for line in f:
            p2 = line.rstrip('\n').split('\t')
            if len(p2) < len(header):
                continue
            g = lambda k: p2[idx[k]] if k in idx and idx[k] < len(p2) else ''
            if g('contig') in contig_ids:
                out.append({'contig': g('contig'),
                            'species': g('species'), 'family': g('family'),
                            'host_ictv': g('host_ictv'),
                            'final_host': g('final_host'),
                            'confidence': g('integrated_confidence')})
    return out


# ------------------------------------------------------------------
# 溯源报告
# ------------------------------------------------------------------
def _fmt(v, suf=''):
    if v is None or v == '':
        return '—'
    if isinstance(v, float):
        return f'{v:g}{suf}'
    return f'{v}{suf}'


def _kmer_disp(v):
    """k-mer 覆盖度展示：≤1 视为比例转百分比；>1 视为已是百分数。"""
    if not isinstance(v, (int, float)):
        return '—'
    return f'{v * 100:.1f}%' if v <= 1 else f'{v:.1f}%'


def _fig_to_div(fig, div_id):
    return fig.to_html(include_plotlyjs=False, full_html=False,
                       div_id=div_id, default_width='100%')


def build_report(name):
    """生成 trace_report.html（离线，双击或平台内打开均可）。"""
    from .viz import _plotly_js_path
    import plotly.graph_objects as go

    agg = aggregate(name)
    meta, segs = agg['meta'], agg['segments']
    jdir = _job_dir(name)
    union = sorted(agg['union_rows'],
                   key=lambda r: -(_f(r.get('kmer_cov'), 0) or 0))
    organisms = agg['organisms']

    figs = []
    # ① 物种分布（union）
    top_org = list(organisms.items())[:20][::-1]
    if top_org:
        fig = go.Figure(go.Bar(
            x=[c for _, c in top_org], y=[k for k, _ in top_org],
            orientation='h', text=[c for _, c in top_org],
            textposition='auto', marker_color='#2e7d32'))
        fig.update_layout(title=f'匹配样本物种分布（Top {len(top_org)}，'
                                f'共 {len(organisms)} 种）',
                          height=max(380, 24 * len(top_org) + 90),
                          margin=dict(l=220), xaxis_title='匹配样本数')
        figs.append(('organisms', '病毒序列出现的物种（公共样本计数）', fig))

    # ② ANI vs k-mer 覆盖度散点（有 ANI 列才画）
    pts = [r for r in union if _f(r.get('ani')) is not None
           and _f(r.get('kmer_cov')) is not None]
    if pts:
        orgs_top = [o for o, _ in list(organisms.items())[:6]]
        palette = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00',
                   '#a65628']
        cmap = {o: palette[i % len(palette)] for i, o in enumerate(orgs_top)}
        fig = go.Figure()
        shown = set()
        for r in pts:
            org = r.get('organism') or '（未知物种）'
            c = cmap.get(org, '#8a99a8')
            show = org in cmap and org not in shown
            if show:
                shown.add(org)
            fig.add_trace(go.Scatter(
                x=[r['kmer_cov']], y=[r['ani']], mode='markers',
                marker=dict(color=c, size=8, opacity=0.75),
                name=org if org in cmap else None, showlegend=bool(show),
                hovertemplate=(f"{r.get('acc') or '—'}<br>{org}<br>"
                               f"k-mer {_kmer_disp(r['kmer_cov'])}<br>"
                               f"ANI {r['ani']:.1f}%<extra></extra>")))
        fig.update_layout(title='共享 k-mer 比例 vs ANI 估计（按物种着色）',
                          xaxis_title='共享 k-mer（结果表原值）',
                          yaxis_title='ANI 估计 (%)', height=460)
        figs.append(('scatter', '近似度：k-mer 共享 × ANI', fig))

    # ③ 样本类型 / 平台分布
    if agg['assays']:
        items = list(agg['assays'].items())[:12]
        fig = go.Figure(go.Bar(
            x=[k for k, _ in items], y=[c for _, c in items],
            text=[c for _, c in items], textposition='auto',
            marker_color='#1565c0'))
        fig.update_layout(title='匹配样本的测序类型（assay_type）',
                          height=380, yaxis_title='样本数')
        figs.append(('assay', '样本类型分布', fig))

    # 交叉对照表（平台 ④宿主预测 vs LOGAN 物种）
    cross = _host_cross(meta, {s['contig'] for s in segs})
    cross_top = {s['contig']:
                 [o for o in list(organisms)[:5]] for s in segs}

    def _seg_table(seg):
        rows = sorted(agg['per_seg'].get(seg['index'], []),
                      key=lambda r: -(_f(r.get('kmer_cov'), 0) or 0))
        if not rows:
            return ('<p class="hint">该片段尚未导入结果。请到「LOGAN 溯源」'
                    '页导入 Logan-Search 返回的结果表。</p>')
        cols = [('acc', 'Run Accession'), ('organism', '物种'),
                ('kmer_cov', '共享 k-mer'), ('ani', 'ANI(%)'),
                ('assay_type', '测序类型'), ('study', 'Study'),
                ('biosample', 'BioSample'), ('location', '地点')]
        thead = ''.join(f'<th>{label}</th>' for _, label in cols)
        body = []
        for r in rows[:300]:
            tds = []
            for f, _ in cols:
                v = r.get(f)
                if f == 'kmer_cov':
                    v = _kmer_disp(v)
                tds.append(f"<td>{_fmt(v)}</td>")
            body.append('<tr>' + ''.join(tds) + '</tr>')
        accs = [r.get('acc') for r in rows if r.get('acc')]
        acc_txt = ', '.join(accs[:500])
        extra = (f"<details><summary>复制 Run 列表（{len(accs)} 条）</summary>"
                 f"<textarea readonly class='acc-list'>{acc_txt}</textarea>"
                 f"</details>") if accs else ''
        more = f"<p class='hint'>仅显示前 300 行，完整结果见导入的原始文件。</p>" \
            if len(rows) > 300 else ''
        return (f"<table class='rtable'><thead><tr>{thead}</tr></thead>"
                f"<tbody>{''.join(body)}</tbody></table>{extra}{more}")

    seg_sections = ''.join(
        f"<h3>片段 s{s['index']} · {s['contig']} "
        f"({s['start']:,}-{s['end']:,} / {s['length']:,}bp)</h3>"
        + _seg_table(s) for s in segs)

    contig_lines = ''.join(
        f"<tr><td>{c['contig']}</td><td>{c['length']:,}</td>"
        f"<td>{c.get('species') or '—'}</td><td>{c.get('family') or '—'}</td></tr>"
        for c in meta.get('contigs', []))

    cross_html = ''
    if cross:
        rows = []
        for c in cross:
            rows.append(
                f"<tr><td>{c['contig']}</td><td>{c['species'] or '—'}</td>"
                f"<td>{c['family'] or '—'}</td><td>{c['final_host'] or '—'}"
                f"（ICTV: {c['host_ictv'] or '—'}, 置信 {c['confidence'] or '—'}）</td>"
                f"<td>{'、'.join(cross_top.get(c['contig'], [])) or '—'}</td></tr>")
        cross_html = (
            "<h2 class='sec'>交叉对照：平台 ④ICTV 宿主预测 × LOGAN 实测物种</h2>"
            "<p>左侧为本平台基于分类学的宿主预测，右侧为 LOGAN 在公共样本中"
            "实测到该序列的 Top 物种（宿主物种名出现即提示其感染宿主/携带者）。</p>"
            "<table class='rtable'><thead><tr><th>Contig</th><th>病毒种</th>"
            "<th>科</th><th>平台宿主预测</th><th>LOGAN Top 物种</th></tr></thead>"
            f"<tbody>{''.join(rows)}</tbody></table>")

    js = _plotly_js_path(jdir) or ''
    # 必须放 <head>：图形内联脚本执行 Plotly.newPlot 时库须已就绪
    js_tag = (f"<script src=\"{js}\"></script>" if js else
              '<script>window.Plotly||document.write("plotly.min.js 缺失")'
              '</script>')
    fig_html = ''.join(
        f"<h2 class='sec'>{title}</h2>{_fig_to_div(fig, fid)}"
        for fid, title, fig in figs)

    status_txt = {'pending': '待提交', 'partial': '部分片段已导入',
                  'done': '结果齐全'}.get(meta.get('status'), '')
    html = f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>LOGAN 溯源报告 · {meta['name']}</title>
{js_tag}
<style>
 body {{ font-family: "Microsoft YaHei", "Segoe UI", sans-serif;
        margin: 0; background: #f5f7fa; color: #24313f; }}
 .wrap {{ max-width: 1080px; margin: 0 auto; padding: 24px 16px 60px; }}
 h1 {{ font-size: 22px; }} h2.sec {{ font-size: 17px; margin: 28px 0 10px;
        border-left: 4px solid #2e7d32; padding-left: 10px; }}
 h3 {{ font-size: 15px; margin: 20px 0 8px; color: #37505f; }}
 .cards {{ display: flex; gap: 12px; flex-wrap: wrap; margin: 14px 0; }}
 .card {{ background: #fff; border-radius: 10px; padding: 12px 18px;
         box-shadow: 0 1px 4px rgba(30,50,70,.08); min-width: 130px; }}
 .card b {{ display: block; font-size: 20px; color: #1d5c38; }}
 .card span {{ font-size: 12px; color: #6a7c8c; }}
 table.rtable {{ border-collapse: collapse; width: 100%; background: #fff;
        font-size: 13px; }}
 .rtable th, .rtable td {{ border: 1px solid #dfe6ec; padding: 5px 8px;
        text-align: left; white-space: nowrap; }}
 .rtable th {{ background: #eef3f7; position: sticky; top: 0; }}
 .rtable tbody tr:hover {{ background: #f4f9f4; }}
 .hint {{ color: #6a7c8c; font-size: 12.5px; }}
 textarea.acc-list {{ width: 100%; height: 90px; font-size: 12px;
        margin-top: 6px; }}
 .meta {{ font-size: 13px; color: #4b5d6d; line-height: 1.8; }}
 details {{ margin: 8px 0; font-size: 13px; }}
 a {{ color: #1565c0; }}
</style></head><body><div class="wrap">
<h1>🌐 LOGAN 溯源报告 · {meta['name']}</h1>
<p class="meta">
 来源：{('样品 ' + meta['sample']) if meta.get('sample') else '粘贴序列'}
 ｜ 状态：{status_txt} ｜ 创建：{meta.get('created', '')}
 ｜ 查询引擎：<a href="{meta.get('submit_url', '')}" target="_blank">Logan-Search</a>
 （NCBI SRA 全量组装 k-mer 索引）<br>
 提交片段：{meta.get('n_segments', 0)} 条（每条 ≤{SEG_MAX_BP}bp）
 ｜ 去重后匹配样本：{agg['n_union']:,}
 ｜ 物种数：{len(organisms):,}
</p>
<div class="cards">
 <div class="card"><b>{agg['n_union']:,}</b><span>匹配公共样本（去重）</span></div>
 <div class="card"><b>{len(organisms):,}</b><span>涉及物种</span></div>
 <div class="card"><b>{(_fmt(list(organisms)[0]) if organisms else '—')}</b>
   <span>样本数最多的物种</span></div>
 <div class="card"><b>{meta.get('n_segments', 0)}</b><span>提交片段</span></div>
</div>
{fig_html}
{cross_html}
<h2 class="sec">查询片段明细</h2>
<table class="rtable"><thead><tr><th>Contig</th><th>长度</th>
<th>物种（平台鉴定）</th><th>科</th></tr></thead>
<tbody>{contig_lines}</tbody></table>
<h2 class="sec">逐片段匹配样本</h2>
{seg_sections}
<h2 class="sec">方法与引用</h2>
<p class="meta">查询流程：平台切片（≤{SEG_MAX_BP}bp）→ Logan-Search
（kmindex，k=31 Bloom filter 索引，对 SRA 全量组装 unitigs 计算共享
k-mer 比例）→ 导入结果表聚合。Threshold 默认 0.25；Groups 提交时选择：
{meta.get('group_hint', '')}。<br>引用：{LOGAN_CITE}。
报告由植物病毒分析平台生成于 {time.strftime('%Y-%m-%d %H:%M:%S')}。</p>
</div></body></html>"""
    out = check_path(os.path.join(jdir, 'trace_report.html'),
                     must_exist=False, in_platform=True)
    with safe_open(out, 'wt') as f:
        f.write(html)
    return out

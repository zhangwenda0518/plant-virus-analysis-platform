# -*- coding: utf-8 -*-
"""
suvtk_submit.py — 本地 .sqn（NCBI 提交文件）生成编排
====================================================

把平台 /submit 表（unified metadata 列集）+ 关联的 orf 运行组装成
suvtk 提交链输入，调用补丁版 suvtk（features → comments → table2asn）产出
可直接投递 NCBI（gb-sub@ncbi.nlm.nih.gov / SRA 门户）的 .sqn：

    表行(sequence_name / organism / src-* / cmt-* / authors)
      + orf 运行(CDS tbl 坐标 + orfa genome_composition → pred_genome_type)
    → sequences.fsa + featuretable.tbl + source.src + template.sbt + .cmt
    → suvtk table2asn → <prefix>.sqn + <prefix>.val（NCBI 官方校验）

依赖（均已就位）：
  - suvtk 0.2.0 editable 安装（带 2 处本地补丁：features 训练恢复、table2asn 补
    import pandas；见 2026-09-08 工作日志）
  - NCBI table2asn 二进制: tools/table2asn/table2asn.exe
  - mmseqs: tools/mmseqs/bin/mmseqs.exe
  - BFVD mmseqs 库最小集: databases/misc/suvtk/（bfvd 索引族 + 两张 TSV）
"""
import os
import sys
import csv
import json
import subprocess
import tempfile
from pathlib import Path
import pandas as pd

from .config import DIRS, get_config, db_path
from .utils import safe_open, iter_fasta, write_fasta_record

PLATFORM_ROOT = str(Path(__file__).resolve().parents[1])
SUVTK_DB_DEFAULT = db_path('misc', 'suvtk')
_MMSEQS_BIN = os.path.join(PLATFORM_ROOT, 'tools', 'mmseqs', 'bin')
_T2A_BIN = os.path.join(PLATFORM_ROOT, 'tools', 'table2asn')

# NCBI MIUVIG:5.0 结构化注释字段（行级 collection_date/geo_loc_name/lat_lon 由
# taxonomy.tsv 提供，文件级参数由本表提供）。软件字段须 `软件;版本;参数` 三段格式。
DEFAULT_MIUVIG = {
    'source_uvig': 'metatranscriptome (not viral targeted)',
    'vir_ident_software': 'kun_peng;0.7.12;default parameters',
    'assembly_software': 'SPAdes;4.3.0;metaviral',
    'assembly_qual': 'Genome fragment(s)',
    'detec_type': 'independent sequence (UViG)',
    'number_contig': '1',
    'virus_enrich_appr': 'none',
    'nucl_acid_ext': 'none',
    'wga_amp_appr': 'none',
    'env_broad_scale': 'host-associated',
    'env_local_scale': 'leaf',
    'env_medium': 'plant tissue',
    'project_name': 'Lycium virome project',
    'seq_meth': 'Illumina NovaSeq 6000',
    'samp_taxon_id': 'plant metagenome [NCBITaxon:1297885]',
    'feat_pred': 'pyrodigal;3.7.1;single mode',
    'ref_db': 'BFVD;2023_02;https://bfvd.steineggerlab.workers.dev',
    'sim_search_meth': 'MMseqs2;14.7e284;default parameters',
    'size_frac': '',
}


def _logger(log):
    return (lambda msg: log(msg)) if log else (lambda msg: None)


def suvtk_db_dir():
    """BFVD 库目录：默认 databases/misc/suvtk，platform.json databases.suvtk 可覆盖。"""
    cfg = get_config()
    p = (getattr(cfg, 'databases', {}) or {}).get('suvtk')
    if p and os.path.isdir(str(p)):
        return str(p)
    return SUVTK_DB_DEFAULT if os.path.isdir(SUVTK_DB_DEFAULT) else SUVTK_DB_DEFAULT


def _suvtk_python():
    """返回能 `import suvtk` 的解释器路径。

    探测顺序：platform.json databases.suvtk_python → 环境变量 SUVTK_PYTHON →
    系统 C:\\Python312\\python.exe（若存在且可 import suvtk）→ sys.executable。
    suvtk 是 editable 安装在系统 Python 3.12 的，平台主进程（managed 3.13）
    未必有，故子进程用独立解释器。"""
    cfg = get_config()
    p = (getattr(cfg, 'databases', {}) or {}).get('suvtk_python')
    cands = []
    if p:
        cands.append(str(p))
    cands.append(os.environ.get('SUVTK_PYTHON', ''))
    cands.append(r'C:\Python312\python.exe')
    cands.append(sys.executable)
    for c in cands:
        if not c:
            continue
        if os.path.isfile(c):
            try:
                r = subprocess.run([c, '-c', 'import suvtk'], capture_output=True,
                                   timeout=30)
                if r.returncode == 0:
                    return c
            except (OSError, subprocess.TimeoutExpired):
                continue
    return sys.executable


def _check_ready(log=None):
    """环境自检：suvtk python 模块 / mmseqs / table2asn / BFVD。返回错误列表(空=OK)。"""
    errs = []
    suvtk_py = _suvtk_python()
    try:
        r = subprocess.run([suvtk_py, '-c', 'import suvtk'], capture_output=True,
                           timeout=30)
        if r.returncode != 0:
            errs.append('suvtk 未安装（需 editable 安装并带补丁：'
                        'pip install -e E:\\BandZip解压\\suvtk_db\\suvtk）')
        else:
            log = _logger(log)
            log(f'  suvtk 解释器: {suvtk_py}')
    except (OSError, subprocess.TimeoutExpired):
        errs.append(f'suvtk 解释器不可用: {suvtk_py}')
    if not os.path.isfile(os.path.join(_MMSEQS_BIN, 'mmseqs.exe')) and \
            not os.path.isfile(os.path.join(_MMSEQS_BIN, 'mmseqs')):
        errs.append(f'mmseqs 缺失: {_MMSEQS_BIN}')
    if not os.path.isfile(os.path.join(_T2A_BIN, 'table2asn.exe')):
        errs.append(f'table2asn 缺失: {_T2A_BIN}\\table2asn.exe')
    if not os.path.isfile(os.path.join(suvtk_db_dir(), 'bfvd')):
        errs.append(f'BFVD 库缺失: {suvtk_db_dir()}（需 bfvd 索引族）')
    return errs


def _env():
    env = os.environ.copy()
    env['PATH'] = _T2A_BIN + os.pathsep + _MMSEQS_BIN + os.pathsep + \
        env.get('PATH', '')
    return env


def _run_suvtk(args, cwd, log=None, timeout=None):
    """[<suvtk_python>, -m, suvtk, ...] 子进程；返回 (rc, stdout+stderr)。"""
    cmd = [_suvtk_python(), '-m', 'suvtk'] + args
    log = _logger(log)
    log('  $ ' + ' '.join(cmd))
    try:
        r = subprocess.run(cmd, cwd=cwd, env=_env(), capture_output=True,
                           text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise RuntimeError('suvtk 子命令超时: ' + ' '.join(args[:3]))
    out = (r.stdout or '') + (r.stderr or '')
    if r.returncode != 0:
        tail = '\n'.join(out.splitlines()[-12:])
        raise RuntimeError(f'suvtk {" ".join(args[:1])} 失败 (rc={r.returncode}):\n{tail}')
    return out


def run_features(fasta, out_dir, db_dir=None, threads=None, log=None,
                 work_dir=None):
    """对 fasta 跑 suvtk features（BFVD 功能 → featuretable.tbl 等）。幂等。

    work_dir: suvtk features 的工作目录（mmseqs 临时 `tmp/` 生成于此，跑完删除）。
    默认取 out_dir 父目录；当 out_dir 位于受保护的 submissions/ 下时，须把
    work_dir 指到 tool_runs 目录，避免 mmseqs 清理 tmp 触发平台批量删除保护。
    """
    log = _logger(log)
    out_dir = str(out_dir)
    tbl = os.path.join(out_dir, 'featuretable.tbl')
    if os.path.isfile(tbl):
        log(f'  suvtk features 产物已存在，跳过: {out_dir}')
    else:
        os.makedirs(out_dir, exist_ok=True)
        db_dir = db_dir or suvtk_db_dir()
        n = get_config().threads if threads is None else int(threads)
        wd = work_dir or os.path.dirname(out_dir) or os.getcwd()
        os.makedirs(wd, exist_ok=True)
        _run_suvtk(['features', '-i', os.path.abspath(fasta),
                    '-o', out_dir, '-d', db_dir, '-t', str(n)],
                   cwd=wd, log=log, timeout=1800)
        # 注：suvtk 0.2.0 已自带 transl_table/codon_start；遗传密码冲突
        # GenCodeMismatch 由 build_sqn 的 table2asn `-j "[gcode=1]"` 在
        # BioSource 层统一解决，无需给 tbl 补限定词。
    return {
        'tbl': tbl,
        'miuvig': os.path.join(out_dir, 'miuvig_features.tsv'),
        'fna': os.path.join(out_dir, 'reoriented_nucleotide_sequences.fna'),
        'faa': os.path.join(out_dir, 'proteins.faa'),
    }


def dominant_genome_type(run_dir):
    """orf 运行 04b 注释 → 主要病毒科的基因组类型（pred_genome_type）。

    读 orf_annotation.tsv 的 family 计数 + ICTV 科级 gencomp 表。无则 None。
    """
    ann = os.path.join(run_dir, '04b_orf_annot', 'orf_annotation.tsv')
    if not os.path.isfile(ann):
        return None
    fam = {}
    with open(ann, encoding='utf-8', errors='replace') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            k = (r.get('family') or '').strip()
            if k:
                fam[k] = fam.get(k, 0) + 1
    if not fam:
        return None
    top = max(fam, key=fam.get)
    gencomp = os.path.join(db_path('annot', 'prot'),
                           'ictv_family_gencomp.tsv')
    if os.path.isfile(gencomp):
        with open(gencomp, encoding='utf-8', errors='replace') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                if (r.get('family') or '').strip() == top:
                    return (r.get('genome.composition')
                            or r.get('genome_composition')
                            or r.get('gencomp') or '').strip() or None
    return None


def _clean(v):
    return str(v or '').strip()


def build_source_src(rows):
    """表行(list[dict]) → source.src 文本（列缺失容忍为空）。"""
    cols = ['sequence_name', 'organism', 'src-Isolate', 'collection_date',
            'src-geo_loc_name', 'src-Lat_Lon', 'bioproject', 'biosample',
            'sra', 'src-Isolation-source', 'src-Segment', 'src-Host',
            'src-Tissue_type', 'src-Cultivar', 'src-Dev_stage',
            'src-Collected_by']
    lines = ['Sequence_ID\tOrganism\tIsolate\tCollection_date\tgeo_loc_name\t'
             'Lat_Lon\tBioproject\tBiosample\tSRA\tMetagenomic\t'
             'Metagenome_source\tSegment\tHost\tTissue_type\tCultivar\t'
             'Dev_stage\tCollected_by']
    for r in rows:
        def g(c):
            return _clean(r.get(c))
        lines.append('\t'.join([
            g('sequence_name'), g('organism'), g('src-Isolate'),
            g('collection_date'), g('src-geo_loc_name'), g('src-Lat_Lon'),
            g('bioproject'), g('biosample'), g('sra'), 'TRUE',
            g('src-Isolation-source') or 'plant virome', g('src-Segment'),
            g('src-Host'), g('src-Tissue_type'), g('src-Cultivar'),
            g('src-Dev_stage'), g('src-Collected_by')]))
    return '\n'.join(lines) + '\n'


def build_taxonomy_tsv(rows, gtype_default=None):
    """表行 → suvtk comments -t 输入（contig 级行 + MIUVIG 行级字段）。"""
    gtype = gtype_default or 'uncharacterized'
    cols = ['contig', 'pred_genome_type', 'pred_genome_struc',
            'collection_date', 'geo_loc_name', 'lat_lon', 'host']
    lines = ['\t'.join(cols)]
    for r in rows:
        seg = _clean(r.get('src-Segment'))
        struc = 'segmented' if seg else 'non-segmented'
        lines.append('\t'.join([
            _clean(r.get('sequence_name')), gtype, struc,
            _clean(r.get('collection_date')), _clean(r.get('src-geo_loc_name')),
            _clean(r.get('src-Lat_Lon')), _clean(r.get('src-Host'))]))
    return '\n'.join(lines) + '\n'


def _platform_feat_miuvig(run_dir, out_path):
    """平台自家 orf/tbl 链的 miuvig_features 等价物（suvtk comments -f 输入）。

    从 orf 运行 04b 产物提取完整注释信息：
    - feat_pred / sim_search_meth / viral_fams_pred / cdd_pred
    - evidence_summary / family_summary / category_summary
    - ref_db / ref_db_version
    """
    ann_tsv = os.path.join(run_dir, '04b_orf_annot', 'orf_annotation.tsv')
    summ_p = os.path.join(run_dir, '04b_orf_annot', 'summary.json')

    # 基础参数（从 summary.json 取，它已有 n_orfs/n_annotated/families/categories）
    model = 'pyrodigal'
    engine = 'diamond'
    n_orfs = 0
    n_annotated = 0
    families = {}
    categories = {}
    evidence_count = {}
    hmm_total = 0
    cdd_total = 0

    if os.path.isfile(summ_p):
        try:
            with open(summ_p, encoding='utf-8') as f:
                s = json.load(f)
            model = s.get('model') or model
            engine = s.get('engine') or engine
            n_orfs = int(s.get('n_orfs') or 0)
            n_annotated = int(s.get('n_annotated') or 0)
            families = dict(s.get('families') or {})
            categories = dict(s.get('categories') or {})
        except (OSError, ValueError):
            pass

    # 从 orf_annotation.tsv 补充 evidence / hmm / cdd 详细统计
    if os.path.isfile(ann_tsv):
        try:
            cdf = pd.read_csv(ann_tsv, sep='\t', dtype=str, keep_default_na=False)
            # evidence 统计（合并 evidence + informative 列）
            for ev_col in ['evidence', 'informative']:
                if ev_col in cdf.columns:
                    for v in cdf[ev_col].tolist():
                        for part in str(v).split('|'):
                            part = part.strip()
                            if part:
                                evidence_count[part] = evidence_count.get(part, 0) + 1
            # hmm_hits 列：每行有值即计一次（非空字符串）
            if 'hmm_hits' in cdf.columns:
                hmm_total = int(cdf['hmm_hits'].ne('').sum())
            # cdd_hits 列同理
            if 'cdd_hits' in cdf.columns:
                cdd_total = int(cdf['cdd_hits'].ne('').sum())
        except (OSError, ValueError):
            pass

    # 组合参数
    # 中文 category 映射到英文（table2asn 要求 ASCII-only）
    _CAT_MAP = {
        '结构蛋白': 'structural_protein',
        '其他功能蛋白': 'other_functional_protein',
        '聚合酶/复制相关': 'polymerase_replication',
        '复制酶': 'replicase',
        '包膜蛋白': 'envelope_protein',
        '转录相关': 'transcription',
    }
    n_evidenced = len([c for c in evidence_count.values() if c > 0])
    feat_pred = f'{model};3.7.1;single mode'
    sim_search_meth = f'{engine};2.1.8;default'
    viral_fams_pred = f'HMMER;3.4;VOGdb/RVDB/vFam'
    cdd_pred = f'MMseqs2;14.7e284;CDD+Pfam'
    ref_db = 'RefSeq viral'
    ref_db_version = '2024-06'

    with safe_open(out_path, 'wt') as f:
        f.write('MIUVIG_parameter\tvalue\n')
        f.write(f'feat_pred\t{feat_pred}\n')
        f.write(f'sim_search_meth\t{sim_search_meth}\n')
        f.write(f'viral_fams_pred\t{viral_fams_pred}\n')
        f.write(f'cdd_pred\t{cdd_pred}\n')
        f.write(f'ref_db\t{ref_db}\n')
        f.write(f'ref_db_version\t{ref_db_version}\n')
        f.write(f'n_contigs\t{n_orfs}\n')
        f.write(f'n_annotated\t{n_annotated}\n')
        f.write(f'evidence_summary\tseq={n_evidenced}, '
                f'hmm={hmm_total}, cdd={cdd_total}\n')
        # 注: investigation_type 是 MIxS/MIGS 字段（prefix ##MIGS-Data-START##），
        # 不属于 MIUVIG:5.0-Data。suvtk 的 comments 命令把所有字段放在同一个
        # MIUVIG cmt 里，但 table2asn 校验器把 investigation_type 当作 MIxS
        # 必填字段检查。MIUVIG 提交时该字段可选，保留 Info 提示不影响 .sqn 产出。
        if families:
            fam_str = '; '.join(f'{k}:{v}' for k, v in sorted(families.items()))
            f.write(f'family_summary\t{fam_str}\n')
        if categories:
            # 中文 → 英文映射
            translated = {_CAT_MAP.get(k, k): v for k, v in categories.items()}
            cat_str = '; '.join(f'{k}:{v}' for k, v in sorted(translated.items()))
            f.write(f'category_summary\t{cat_str}\n')
    return out_path


def build_sqn(run_dir, *, seq_fasta=None, rows=None, author=None,
              miuvig=None, assembler=None, sequencer=None,
              use_features=True, tbl_path=None, out_dir=None, threads=None,
              log=None):
    """核心编排：给定 orf 运行（含 03_assembly / 04_orf / 04b）+ 提交表行，
    生成 .sqn/.cmt/.val。返回产物 dict。

    rows: 表行 list[dict]（sequence_name/organism/src-*/collection_date/...）。
    author: generate_sbt 字段 {last,first,middle,affil,div,city,sub,country,
            street,email,postal}。
    use_features=True: 先跑 suvtk features（BFVD）产 tbl；False: 用平台自家 tbl
    （tbl_path 须指向 store.generate_feature_tbl 产物）。
    """
    log = _logger(log)
    errs = _check_ready(log)
    if errs:
        raise RuntimeError('环境未就绪:\n' + '\n'.join('- ' + e for e in errs))

    run_dir = str(run_dir)
    if not os.path.isdir(run_dir):
        raise RuntimeError(f'orf 运行目录不存在: {run_dir}')
    rows = list(rows or [])
    if not rows:
        raise RuntimeError('没有可用的提交行（表为空或未从 contigs 导入）')

    # ---- 1) 序列 FASTA ----
    fasta = seq_fasta or os.path.join(run_dir, '03_assembly',
                                      'viral_contigs.fasta')
    if not os.path.isfile(fasta):
        raise RuntimeError(f'序列 FASTA 缺失: {fasta}')
    out_dir = out_dir or os.path.join(run_dir, 'ncbi_sqn')
    os.makedirs(out_dir, exist_ok=True)

    # ---- 1b) 提交用 FASTA（defline 带 organism + isolate 描述）----
    # NCBI table2asn 要求 defline 形如：
    #   >ACCESSION Organism name isolate XXX, complete genome
    # 裸 ID（>OR489165.1）或缺 isolate 都会让 table2asn 静默失败（exit 2）。
    org_map = {_clean(r.get('sequence_name')): _clean(r.get('organism'))
               for r in rows if _clean(r.get('sequence_name'))}
    iso_map = {_clean(r.get('sequence_name')): _clean(r.get('src-Isolate'))
               for r in rows if _clean(r.get('sequence_name'))}
    sub_fasta = os.path.join(out_dir, 'submission.fasta')
    n_seq = 0
    with safe_open(sub_fasta, 'wt') as f:
        for hdr, seq in iter_fasta(fasta):
            sid = hdr.split()[0]
            org = org_map.get(sid) or org_map.get(hdr) or ''
            iso = iso_map.get(sid) or iso_map.get(hdr) or ''
            defline = f'{sid} {org}'.rstrip()
            if iso:
                defline += f' isolate {iso}'
            defline += ', complete genome'
            write_fasta_record(f, defline, seq)
            n_seq += 1
    if n_seq == 0:
        raise RuntimeError(f'序列 FASTA 为空: {fasta}')
    log(f'  提交 FASTA 准备完成: {sub_fasta}（{n_seq} 条，defline 含 organism/isolate）')

    # ---- 2) feature table .tbl ----
    # suvtk features 会在 cwd 下建 `tmp/`（mmseqs 临时索引）跑完再删；若 out_dir
    # 位于受保护的 submissions/ 下，mmseqs 批量删除 tmp 会触发平台 SAFE_DELETE
    # 保护。故把 suvtk 的 cwd 指到 run_dir 下的安全工作目录，仅输出落到 out_dir。
    work_dir = os.path.join(run_dir, 'ncbi_sqn_work')
    if use_features:
        feats = run_features(fasta, os.path.join(out_dir, 'suvtk_features'),
                             threads=threads, log=log, work_dir=work_dir)
        tbl = feats['tbl']
        feat_miuvig = feats['miuvig']
    else:
        if not tbl_path or not os.path.isfile(str(tbl_path)):
            raise RuntimeError('自家 tbl 模式需 tbl_path 指向 featuretable.tbl'
                               '（先在 /submit 生成）')
        tbl = str(tbl_path)
        feat_miuvig = _platform_feat_miuvig(
            run_dir, os.path.join(out_dir, 'miuvig_features.tsv'))

    # ---- 3) source.src / template.sbt（编辑复用：已存在则不覆盖）----
    src_p = os.path.join(out_dir, 'source.src')
    if os.path.isfile(src_p) and os.path.getsize(src_p) > 0:
        log(f'  source.src 已存在，复用（可编辑后重新生成）: {src_p}')
    else:
        src_text = build_source_src(rows)
        with safe_open(src_p, 'wt') as f:
            f.write(src_text)
    sbt_p = os.path.join(out_dir, 'template.sbt')
    if os.path.isfile(sbt_p) and os.path.getsize(sbt_p) > 0:
        log(f'  template.sbt 已存在，复用: {sbt_p}')
    else:
        author = author or {}
        try:
            from .ncbi_submit.store import generate_sbt
            template, _ = generate_sbt(author)
        except ValueError as e:
            raise RuntimeError(f'作者模板不完整: {e}')
        with safe_open(sbt_p, 'wt') as f:
            f.write(template)

    # ---- 4) MIUVIG 文件级参数 + assembly（编辑复用）----
    mv = dict(DEFAULT_MIUVIG)
    if miuvig:
        mv.update({k: _clean(v) for k, v in miuvig.items() if _clean(v)})
    mv_p = os.path.join(out_dir, 'miuvig.tsv')
    if not (os.path.isfile(mv_p) and os.path.getsize(mv_p) > 0):
        with safe_open(mv_p, 'wt') as f:
            f.write('MIUVIG_parameter\tvalue\n')
            for k, v in mv.items():
                if v:
                    f.write(f'{k}\t{v}\n')
    asm_p = os.path.join(out_dir, 'assembly.tsv')
    if not (os.path.isfile(asm_p) and os.path.getsize(asm_p) > 0):
        asm = [('Assembly Method', assembler or 'SPAdes;4.3.0;metaviral'),
               ('Sequencing Technology', sequencer or 'Illumina NovaSeq 6000')]
        with safe_open(asm_p, 'wt') as f:
            f.write('Assembly_parameter\tvalue\n')
            for k, v in asm:
                f.write(f'{k}\t{v}\n')

    # ---- 5) taxonomy 输入（contig 行 + 行级 MIUVIG 字段）----
    gtype = dominant_genome_type(run_dir)
    tax_p = os.path.join(out_dir, 'taxonomy.tsv')
    with safe_open(tax_p, 'wt') as f:
        f.write(build_taxonomy_tsv(rows, gtype_default=gtype))

    # ---- 6) suvtk comments → .cmt（编辑复用：已存在则不覆盖）----
    cmt = os.path.join(out_dir, 'sqn')
    cmt_p = cmt + '.cmt'
    if os.path.isfile(cmt_p) and os.path.getsize(cmt_p) > 0:
        log(f'  sqn.cmt 已存在，复用: {cmt_p}')
    else:
        _run_suvtk(['comments', '-t', tax_p, '-f', feat_miuvig, '-m', mv_p,
                    '-a', asm_p, '-o', cmt],
                   cwd=out_dir, log=log, timeout=300)

    # ---- 7) table2asn → .sqn ----
    # 直接调 NCBI table2asn.exe（绕过 suvtk wrapper：wrapper 仅拼命令+读 .val，
    # 且在 .val 缺失时抛误导性 FileNotFoundError，见 table2asn.py:221 上游 bug）。
    # 关键坑：
    #  1) 输入 defline 必须含 organism + isolate（裸 ID 会被静默拒）；
    #  2) 必须用【相对路径】传参（cwd=out_dir），table2asn 对含中文的绝对路径
    #     （桌面/植物病毒分析平台）窄字符编码失败 → rc=3 静默退出；
    #  3) -o 须以 .sqn 结尾才产出 .sqn；
    #  4) -j "[gcode=1]" 给 BioSource 设标准遗传密码，消除 ssRNA(-) 病毒的
    #     GenCodeMismatch×7 Warning。
    t2a = os.path.join(_T2A_BIN, 'table2asn.exe')
    out_prefix = os.path.join(out_dir, 'sqn')
    rel = lambda p: os.path.relpath(p, out_dir)
    t2a_cmd = [t2a, '-i', rel(sub_fasta), '-o', 'sqn.sqn', '-t', rel(sbt_p),
               '-f', rel(tbl), '-src-file', rel(src_p), '-w', rel(cmt_p),
               '-j', '[gcode=1]', '-V', 'vb', '-a', 's']
    log('  $ ' + ' '.join(t2a_cmd))
    r = subprocess.run(t2a_cmd, cwd=out_dir, capture_output=True,
                       text=True, timeout=600)
    t2a_out = (r.stdout or '') + (r.stderr or '')
    if t2a_out.strip():
        log(t2a_out.strip()[-2000:])
    sqn_p = out_prefix + '.sqn'
    val_p = out_prefix + '.val'
    if not os.path.isfile(sqn_p):
        raise RuntimeError(
            f'table2asn 未产出 .sqn（rc={r.returncode}）。'
            f'table2asn 输出:\n{t2a_out.strip()[-2000:] or "(空)"}')

    # ---- 8) val 摘要 ----
    errs_v = warns_v = infos_v = 0
    val_head = []
    if os.path.isfile(val_p):
        with open(val_p, encoding='utf-8', errors='replace') as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                if s.startswith('Error'):
                    errs_v += 1
                elif s.startswith('Warning'):
                    warns_v += 1
                else:
                    infos_v += 1
                if len(val_head) < 12:
                    val_head.append(s)
    log(f'.sqn 生成完成: {sqn_p}（Error {errs_v} / Warning {warns_v} / '
        f'Info {infos_v}）')
    return {
        'sqn': sqn_p, 'cmt': cmt_p, 'val': val_p, 'tbl': tbl,
        'source_src': src_p, 'sbt': sbt_p,
        'stats': {'error': errs_v, 'warning': warns_v, 'info': infos_v},
        'val_head': val_head,
    }


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description='suvtk .sqn 编排（run_dir 模式）')
    p.add_argument('run_dir')
    p.add_argument('--no-features', action='store_true',
                   help='用平台自家 tbl（未支持, 默认 suvtk features）')
    p.add_argument('--last', default='Zhang')
    p.add_argument('--first', default='Wenda')
    p.add_argument('--affil', default='Ningxia University')
    p.add_argument('--city', default='Yinchuan')
    p.add_argument('--country', default='China')
    p.add_argument('--email', default='zhangwenda@example.com')
    a = p.parse_args()
    print(build_sqn(a.run_dir,
                    author={k: getattr(a, k) for k in
                            ('last', 'first', 'affil', 'city', 'country',
                             'email')}))

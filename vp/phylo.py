# -*- coding: utf-8 -*-
"""
阶段⑤ 进化树与 SDT 分析：
- 按病毒物种分组（BLAST top hit），每组 = 病毒 contigs + top-N 近缘参考
- MAFFT 多序列比对
- trimAl 比对清剪（automated1 自动模式；失败/过度修剪自动回退原比对）
- FastTree（快速）/ NJ（极速，纯 Python identity 距离）/ IQ-TREE
  （可选，模型选择 + UFBoot + SH-aLRT 双支持值）建树
- SDT 全长成对 identity 矩阵（pairwise gap 删除口径，与 SDT 一致）
- 输出 SDT GUI 可直接打开的比对 fasta
"""
import os
import csv
import re
import json
from collections import defaultdict

from .config import get_config
from .utils import (check_path, safe_open, run_cmd, run_cmd_redirect,
                    iter_fasta, write_fasta_record, is_step_done, mark_step_done)
from .assembly import find_virus_ref_fasta


def _safe_name(s):
    s = re.sub(r'[^A-Za-z0-9_.\-]+', '_', str(s)).strip('_')
    return s[:60] or 'group'




def _extract_refs(ref_fastas, wanted_ids, out_fasta):
    """从病毒参考库提取指定 accession 子集（多来源按序扫描，取到即止）。

    wanted_ids 允许含/不含版本号；按完整 accession 与去版本 base 双匹配
    （all_virus/taxa.txt 用无版本号，植物参考 FASTA 头带版本号）。
    """
    if isinstance(ref_fastas, str):
        ref_fastas = [ref_fastas]
    keep = {w for w in wanted_ids if w}
    keep_base = {w.split('.')[0].upper() for w in keep}
    found = []
    with safe_open(out_fasta, 'wt') as f:
        for ref in ref_fastas:
            if not keep_base or not ref or not os.path.isfile(ref):
                continue
            for h, s in iter_fasta(ref):
                acc = h.split()[0]
                base = acc.split('.')[0].upper()
                if acc in keep or base in keep_base:
                    write_fasta_record(f, acc, s)
                    found.append(acc)
                    keep.discard(acc)
                    keep_base.discard(base)
    return found


def _tax_sample_refs(top_acc, ref_tax, genus_accs, family_genera, mode,
                     cap=10, other_count=3):
    """层级抽样选参考（移植 246 服务器 acvirus_tree_pro 的抽样策略）。

    macro:   目标属取 ≤cap 条 + 同科其他属各 other_count 条 —— 科级背景树
    genus:   目标属全取(≤cap) + 同科其他属各 other_count 条 —— 属级树
    lineage: 目标种全取 + 同属其他种各 1 条 —— 种级/株系树

    top hit 在参考索引中无分类信息时返回 None，调用方回退 BLAST 模式。
    返回去版本的 accession base 列表。
    """
    info = ref_tax.get((top_acc or '').split('.')[0].upper())
    if not info:
        return None
    sp, gen, fam = info
    gen_key, fam_key, sp_key = gen.lower(), fam.lower(), sp.lower()
    picked = []
    if mode == 'lineage':
        same_sp, other_sp = [], {}
        for b in genus_accs.get(gen_key, []):
            s2 = (ref_tax[b][0] or '').strip().lower()
            if s2 == sp_key:
                same_sp.append(b)
            elif s2 and s2 not in other_sp:
                other_sp[s2] = b
        picked = same_sp + list(other_sp.values())
    else:
        picked = list(genus_accs.get(gen_key, []))[:cap]
        for g2 in sorted(family_genera.get(fam_key, set())):
            if g2 != gen_key:
                picked += genus_accs.get(g2, [])[:other_count]
    seen = set()
    out = []
    for b in picked:
        if b and b not in seen:
            seen.add(b)
            out.append(b)
    return out[:cap + other_count * 12] or None


def _ascii_stage(path, work_dir, suffix):
    """路径含非 ASCII（平台根含中文）时复制到纯 ASCII 中转目录，返回
    (暂存路径, 是否中转)。用后由调用方 _ascii_fetch 取回并清理。"""
    from .assembly import _ascii_work_base, _is_ascii
    if _is_ascii(path):
        return path, False
    stage_dir = os.path.join(_ascii_work_base('vp_mafft'), work_dir)
    os.makedirs(stage_dir, exist_ok=True)
    staged = os.path.join(stage_dir, f'in{suffix}')
    import shutil
    shutil.copyfile(path, staged)
    return staged, True


def _ascii_fetch(staged, final_path, was_staged):
    if was_staged:
        import shutil
        shutil.copyfile(staged, final_path)
        try:
            os.remove(staged)
        except OSError:
            pass


def _orient_normalize(staged_in, logger=None):
    """比对前 21-mer 定向翻正（de novo 组装的 contig 常为参考的反向互补，
    直接比对会让 identity 矩阵/树失真；平台 MAFFT v6 无 --adjustdirection）。

    以输入内最长 3 条为锚，反向共享显著高于正向的序列翻正后再比对。
    蛋白/纯短序列自动跳过。返回翻正后的 FASTA 路径（无需翻正时返回原路径）。
    staged_in 为 %TEMP% 下的 ASCII 中转文件，此处直接读写（与 _ascii_stage 同口径）。
    """
    from .sdt_exact import auto_orient, detect_seqtype
    recs = []
    try:
        with open(staged_in, 'r', encoding='utf-8', errors='replace') as f:
            cur_name = None
            cur_chunks = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if cur_name is not None:
                        recs.append([cur_name, ''.join(cur_chunks).upper()])
                    cur_name = line[1:].split()[0]
                    cur_chunks = []
                else:
                    cur_chunks.append(line)
            if cur_name is not None:
                recs.append([cur_name, ''.join(cur_chunks).upper()])
    except OSError:
        return staged_in
    if not recs or detect_seqtype([r[1] for r in recs]) != 'nt':
        return staged_in
    if any('-' in seq or '.' in seq for _, seq in recs):
        return staged_in            # 已是带 gap 的比对，翻正会打乱 gap 位置
    flips = auto_orient(recs, [r[1] for r in recs], logger=logger)
    if not any(flips.values()):
        return staged_in
    n_flipped = sum(1 for v in flips.values() if v)
    out = staged_in + '.orient.fa'
    with open(out, 'wt', encoding='utf-8', newline='') as f:
        for name, seq in recs:
            f.write(f'>{name}\n')
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + '\n')
    if logger:
        logger.log(f'比对前定向翻正: {n_flipped} 条反向互补序列已翻正', 'INFO')
    return out


def _run_mafft(in_fasta, out_fasta, threads=None, logger=None,
               strategy='auto'):
    """MAFFT 比对。strategy: auto（缺省）/ linsi（L-INS-i 精确，
    适合 <200 条）/ fast（FFT-NS-1 极速）。"""
    cfg = get_config()
    mafft = cfg.tool('mafft')
    # Windows 原生 mafft 对中文路径不可靠：输入经 ASCII 中转
    staged_in, was_staged = _ascii_stage(
        check_path(in_fasta, must_exist=True), 'mafft_in',
        os.path.splitext(in_fasta)[1] or '.fa')
    aln_input = _orient_normalize(staged_in, logger=logger)
    extra = {'auto': ['--auto'],
             'linsi': ['--localpair', '--maxiterate', '1000'],
             'fast': ['--retree', '1', '--maxiterate', '0']}.get(
                 strategy, ['--auto'])
    cmd = [mafft] + extra + ['--quiet', '--inputorder']
    if threads:
        cmd += ['--thread', str(threads)]
    cmd.append(aln_input)
    if logger:
        logger.log(f"MAFFT 比对: {os.path.basename(in_fasta)}"
                   + ("（ASCII 中转）" if was_staged else ""))
    run_cmd_redirect(cmd, out_fasta, logger=logger)
    if aln_input != staged_in:
        try:
            os.remove(aln_input)
        except OSError:
            pass
    _ascii_fetch(staged_in, in_fasta, was_staged)
    return check_path(out_fasta, must_exist=True)


def _alignment_cols(aln_fasta):
    """比对总列数（按第一条序列长度）。"""
    for _h, s in iter_fasta(aln_fasta):
        return len(s.rstrip())
    return 0


def _run_trimal(in_aln, out_aln, logger=None, min_frac=0.3, min_cols=100):
    """trimAl automated1 自动清剪。

    过度修剪保护：修剪后列数 < 原 30% 或 < min_cols 时回退原比对
    （病毒远缘 contigs 混合比对常触发激进修剪，直接建树会失真）。
    返回 (实际使用的比对路径, 修剪信息 dict)；工具缺失时静默跳过。
    """
    info = {'applied': False, 'cols_before': None, 'cols_after': None}
    cols_before = _alignment_cols(in_aln)
    info['cols_before'] = cols_before
    try:
        trimal = get_config().tool('trimal')
    except (FileNotFoundError, RuntimeError):
        if logger:
            logger.log("trimAl 未安装，跳过比对清剪", "WARN")
        return in_aln, info
    if cols_before <= min_cols:
        if logger:
            logger.log(f"比对仅 {cols_before} 列，无需清剪")
        return in_aln, info
    if logger:
        logger.log("trimAl 清剪 (automated1)")
    staged_in, staged = _ascii_stage(
        check_path(in_aln, must_exist=True), 'trimal', '.aln')
    staged_out = (os.path.join(os.path.dirname(staged_in),
                               'trimmed.aln') if staged else out_aln)
    try:
        run_cmd([trimal, '-in', staged_in,
                 '-out', staged_out, '-automated1'], logger=logger)
        if staged:
            _ascii_fetch(staged_out, out_aln, True)
        cols_after = _alignment_cols(out_aln)
        if not cols_after or not os.path.isfile(out_aln):
            raise RuntimeError("trimAl 无输出")
        if cols_after < cols_before * min_frac or cols_after < min_cols:
            if logger:
                logger.log(f"  修剪过度（{cols_before}→{cols_after} 列），"
                           f"回退未清剪比对", "WARN")
            try:
                os.remove(out_aln)
            except OSError:
                pass
            return in_aln, info
        info.update(applied=True, cols_after=cols_after)
        if logger:
            logger.log(f"  清剪完成: {cols_before} → {cols_after} 列")
        return out_aln, info
    except Exception as e:
        if logger:
            logger.log(f"  trimAl 失败({e})，使用未清剪比对", "WARN")
        info['error'] = str(e)
        return in_aln, info


def _run_fasttree(aln_fasta, tree_out, logger=None):
    cfg = get_config()
    fasttree = cfg.tool('fasttree')
    if logger:
        logger.log("FastTree 建树 (-nt -gtr -gamma)")
    # GTR+Gamma 替代默认 JC：病毒序列进化速率差异大，实测更合理（同 246 用法）
    run_cmd_redirect([fasttree, '-nt', '-gtr', '-gamma',
                      check_path(aln_fasta, must_exist=True)], tree_out, logger=logger)
    tree_out = check_path(tree_out, must_exist=True, in_platform=True)
    # FastTree 参数错误时可能退出码 0 但无输出，防御性校验
    if os.path.getsize(tree_out) == 0:
        raise RuntimeError("FastTree 未输出树（比对序列数不足或参数问题）")
    return tree_out


def _nj_distance_matrix(aln_fasta):
    """比对 FASTA → (names, numpy 距离矩阵)。

    距离 = 1 - identity，identity 用 SDT 口径（仅统计两序列均非 gap 的
    位点，精确字符相等记匹配；无可比位点的对记最大距离 1.0）。
    numpy 按行向量化，几百条序列 × 数万列秒级完成。
    """
    import numpy as np
    names, seqs = [], []
    for h, s in iter_fasta(aln_fasta):
        names.append(h.split()[0])
        seqs.append(s.upper())
    n = len(names)
    if n < 2:
        raise ValueError("比对序列不足 2 条，无法 NJ 建树")
    length = max(len(s) for s in seqs)
    # 右侧补 '.' 占位到等长（与任何碱基不等；且两侧同为 gap 时被掩码排除）
    arr = np.full((n, length), ord('.'), dtype=np.uint8)
    for i, s in enumerate(seqs):
        b = np.frombuffer(s.encode('ascii', 'ignore'), dtype=np.uint8)
        arr[i, :len(b)] = b
    gap = (arr == ord('-')) | (arr == ord('.'))
    dist = np.ones((n, n), dtype=np.float64)
    np.fill_diagonal(dist, 0.0)
    for i in range(n):
        notgap_i = ~gap[i]
        valid = ~gap[i + 1:] & notgap_i          # (n-i-1, length)
        comp = valid.sum(axis=1)
        match = ((arr[i + 1:] == arr[i]) & valid).sum(axis=1)
        ident = np.where(comp > 0, match / np.maximum(comp, 1), 0.0)
        dist[i, i + 1:] = 1.0 - ident
        dist[i + 1:, i] = 1.0 - ident
    return names, dist


def _newick_safe(name):
    """FASTA ID → Newick 叶名（无空格；转义括号/逗号等结构字符）。"""
    s = str(name).split()[0] if str(name).strip() else 'seq'
    return re.sub(r'[(),:;\[\]\'"]+', '_', s) or 'seq'


def _neighbor_joining(dist, names):
    """Saitou–Nei NJ → Newick（无内部节点名/支持值；最终两边各取半长，
    等效简易中点 rooting）。纯 Python + numpy，O(n³) 在数百条序列内秒级。"""
    import numpy as np
    d = np.array(dist, dtype=np.float64, copy=True)
    clusters = [_newick_safe(nm) for nm in names]
    while len(clusters) > 2:
        r = d.sum(axis=1)
        m = len(r)
        q = (m - 2) * d - r[:, None] - r[None, :]
        np.fill_diagonal(q, np.inf)
        flat = int(np.argmin(q))
        i, j = min(flat // m, flat % m), max(flat // m, flat % m)
        d_ij = d[i, j]
        limb_i = max(0.0, 0.5 * d_ij + (r[i] - r[j]) / (2 * (m - 2)))
        limb_j = max(0.0, d_ij - limb_i)
        to_k = 0.5 * (d[i, :] + d[j, :] - d_ij)
        merged = (f"({clusters[i]}:{limb_i:.6f},"
                  f"{clusters[j]}:{limb_j:.6f})")
        # 合并 i/j → 新簇，重建距离矩阵
        keep = [k for k in range(m) if k not in (i, j)]
        d_new = np.zeros((m - 1, m - 1), dtype=np.float64)
        for a, ka in enumerate(keep):
            d_new[a, -1] = d_new[-1, a] = max(0.0, to_k[ka])
            for b, kb in enumerate(keep):
                d_new[a, b] = d[ka, kb]
        d = d_new
        clusters = [clusters[k] for k in keep] + [merged]
    final = max(0.0, float(d[0, 1]) * 0.5)
    return f"({clusters[0]}:{final:.6f},{clusters[1]}:{final:.6f});"


def _run_nj(aln_fasta, tree_out, logger=None):
    """纯 Python NJ 快速建树（identity 距离，SDT 口径）。

    输入为 MAFFT（可选 trimAl 清剪后）比对 FASTA —— 与 FastTree /
    IQ-TREE 同一比对输入；无支持值，定位极速粗树/快速分型。
    """
    if logger:
        logger.log("NJ 建树 (identity 距离 · Saitou-Nei, 纯 Python)")
    names, dist = _nj_distance_matrix(check_path(aln_fasta, must_exist=True))
    nwk = _neighbor_joining(dist, names)
    with safe_open(tree_out, 'wt') as f:
        f.write(nwk)
    return check_path(tree_out, must_exist=False, in_platform=True)


def _parse_iqtree_log(iqtree_file):
    """从 .iqtree 报告提取最优模型与对数似然（用于摘要展示）。"""
    out = {}
    try:
        with safe_open(iqtree_file) as f:
            for line in f:
                m = re.search(r'Best-fit model according to \w+: (\S+)', line)
                if m:
                    out['model'] = m.group(1)
                m = re.search(r'Log-likelihood of the tree:\s*(-?[\d.]+)', line)
                if m:
                    out['logl'] = float(m.group(1))
    except OSError:
        pass
    return out


def _run_iqtree(aln_fasta, prefix, threads=None, logger=None, bootstrap=1000,
                alrt=1000):
    """IQ-TREE 建树：优先 v3，回退 v2。

    -m MFP 自动选模；-B UFBoot + -alrt SH-aLRT 双支持值（v2/v3 参数兼容）。
    """
    cfg = get_config()
    try:
        iqtree = cfg.tool('iqtree3')
        ver = 'v3'
    except (FileNotFoundError, RuntimeError):
        iqtree = cfg.tool('iqtree2')
        ver = 'v2'
    if logger:
        logger.log(f"IQ-TREE {ver} 建树 (MFP + UFBoot {bootstrap} + SH-aLRT {alrt})，"
                   f"耗时较长")
    pf = check_path(prefix, must_exist=False, in_platform=True)
    cmd = [iqtree, '-s', check_path(aln_fasta, must_exist=True),
           '-m', 'MFP', '-B', str(bootstrap), '-alrt', str(alrt),
           '-T', str(threads or cfg.threads), '--prefix', pf, '-redo']
    run_cmd(cmd, logger=logger)
    nwk = pf + '.treefile'
    if not os.path.isfile(nwk):
        raise RuntimeError("IQ-TREE 未生成 .treefile")
    return nwk


def pairwise_identity_matrix(aln_fasta):
    """SDT 口径成对 identity：仅统计两序列均非 gap 的位点。

    返回 (names, matrix) matrix[i][j] = 百分比。
    """
    names, seqs = [], []
    for h, s in iter_fasta(aln_fasta):
        names.append(h.split()[0])
        seqs.append(s.upper())
    n = len(names)
    mat = [[100.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a, b = seqs[i], seqs[j]
            match = comparable = 0
            for x, y in zip(a, b):
                if x != '-' and y != '-':
                    comparable += 1
                    if x == y:
                        match += 1
            ident = match * 100.0 / comparable if comparable else 0.0
            mat[i][j] = mat[j][i] = ident
    return names, mat


def write_sdt_outputs(aln_fasta, out_dir):
    """SDT 矩阵 CSV + SDT GUI 兼容 fasta。"""
    names, mat = pairwise_identity_matrix(aln_fasta)
    csv_path = os.path.join(out_dir, 'sdt_matrix.csv')
    with safe_open(csv_path, 'wt') as f:
        f.write(",".join([''] + names) + "\n")
        for i, nm in enumerate(names):
            f.write(",".join([nm] + [f"{mat[i][j]:.1f}" for j in range(len(names))]) + "\n")
    # SDT GUI 输入：未比对合并序列（SDT 内部自行比对）或比对后序列均可导入
    sdt_fas = os.path.join(out_dir, 'sdt_input.fas')
    with safe_open(sdt_fas, 'wt') as f:
        for h, s in iter_fasta(aln_fasta):
            write_fasta_record(f, h.split()[0], s.replace('-', '').upper())
    return csv_path, sdt_fas, names, mat


def resolve_ncbi_refs(names):
    """NCBI 下载集合名列表 -> fasta 路径列表（供 extra_refs 用）。

    集合位于 databases/ncbi_refs/<name>/refs.fa；名称经白名单校验。
    """
    from .ncbi_download import collection_dir
    paths = []
    for n in names or []:
        p = os.path.join(collection_dir(str(n)), 'refs.fa')
        if not os.path.isfile(p):
            raise FileNotFoundError(f"NCBI 参考集合 [{n}] 不存在（先在工具页下载）")
        paths.append(check_path(p, must_exist=True))
    return paths


def build_phylo(sample_dir, top_n_refs=10, tree_tool='fasttree', threads=None,
                db_virus=None,
                logger=None, force=False, assembly_dir=None, max_groups=5,
                extra_refs=None, sampling='blast', other_count=3,
                do_trim=True):
    """阶段⑤ 主入口。

    extra_refs: 额外参考 fasta 路径列表（如 NCBI 下载集合），追加进每组
    比对（按序列 ID 与组内已有序列去重）。
    sampling: 参考挑选方式 ——
      macro    同科建树：目标属 ≤cap + 同科其他属各 other_count 条
      genus    属级树：目标属全取 + 同科其他属各 other_count 条
      lineage  种级树：目标种全取 + 同属其他种各 1 条
      （macro/genus/lineage 移植自 246 acvirus_tree_pro；
       抽样依赖 virus_classification.tsv 的分类信息）
    """
    step = 'phylo'
    out_dir = check_path(os.path.join(sample_dir, '05_phylo'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    summary_file = os.path.join(out_dir, 'summary.json')
    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("阶段⑤进化分析已完成，跳过")
        with safe_open(summary_file) as f:
            return json.load(f)

    a_dir = check_path(assembly_dir or os.path.join(sample_dir, '03_assembly'),
                       must_exist=True)
    a_summary_p = os.path.join(a_dir, 'summary.json')
    with safe_open(a_summary_p) as f:
        a_summary = json.load(f)
    viral_ids = list(a_summary.get('viral_contigs') or [])
    if not viral_ids:
        msg = "无病毒 contigs，跳过进化分析"
        if logger:
            logger.log(msg, "WARN")
        summary = {'stage': step, 'skipped': True, 'reason': msg, 'groups': []}
        with safe_open(summary_file, 'wt') as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)
        mark_step_done(out_dir, step)
        return summary

    cfg = get_config()
    # 从 virus_classification.tsv 读 kunpeng 分类信息
    cls_map = {}
    cls_tsv = os.path.join(a_dir, 'virus_classification.tsv')
    if os.path.isfile(cls_tsv):
        with safe_open(cls_tsv) as _f:
            for _r in csv.DictReader(_f, delimiter='	'):
                _cid = (_r.get('contig') or '').strip()
                if not _cid:
                    continue
                cls_map[_cid] = {
                    'taxid': (_r.get('taxid') or '').strip(),
                    'family': (_r.get('family') or '').strip(),
                    'genus': (_r.get('genus') or '').strip(),
                    'species': (_r.get('species') or '').strip(),
                }
    if logger:
        logger.log(f"分类表读取: {len(cls_map)} 条 contig 有谱系信息")
    # 兜底：若病毒_classification.tsv 缺失（独立 FASTA 上传），用 kunpeng 分类
    if not cls_map and db_virus and os.path.isfile(contigs_fa):
        if logger:
            logger.log("无分类表，用 kunpeng 分类 contigs ...")
        from .kunpeng import classify, parse_classify_output
        cls_out = os.path.join(a_dir, 'phylo_classify')
        os.makedirs(cls_out, exist_ok=True)
        try:
            kraken_res = classify(db_virus, [contigs_fa], cls_out,
                                  paired=False, threads=threads, logger=logger)
            for flag, _rid, _tx, _len, _path in parse_classify_output(kraken_res['kraken']):
                if flag == 'C' and _tx:
                    # taxid → 8级谱系（用 contig_annot 的 taxonomy 解析）
                    from .contig_annot import lineage_ranks
                    lin = lineage_ranks(_tx)
                    cls_map[_rid] = {
                        'taxid': _tx,
                        'family': lin.get('family', ''),
                        'genus': lin.get('genus', ''),
                        'species': lin.get('species', ''),
                    }
            if logger:
                logger.log(f"kunpeng 兜底分类: {len(cls_map)} 条 contig 有谱系")
        except Exception as e:
            if logger:
                logger.log(f"kunpeng 兜底分类失败（跳过）: {e}", "WARN")
    # 参考挑选限定完整基因组：提取时优先 complete_ref、回退全参考集；
    # 元数据用于 RefSeq/complete 优先
    ref_meta = {}
    ref_sources = []
    try:
        from . import virus_ref
        if virus_ref.available():
            ref_meta = virus_ref.meta_by_accession()
            ref_sources.append(virus_ref.complete_fasta())
    except Exception:
        pass
    ref_fasta = find_virus_ref_fasta()
    ref_sources.append(ref_fasta)

    # 按 top hit 参考物种分组（用 accession 简化为 group key）
    groups = {}
    contigs_fa = os.path.join(a_dir, 'contigs.filtered.fasta')
    contig_seqs = {h.split()[0]: s for h, s in iter_fasta(contigs_fa)}

    # ICTV VMR 参考库（MSL 当前版）：谱系比 acvirus 旧表更新（种改名/
    # 科重分类），gb_cache 缓存序列并入参考池。离线安全——⑦ 内不自动
    # 联网，缺的参考用 ictv-refs --download 预先补齐。
    ictv_meta = {}
    ictv_gb_set = set()
    try:
        from . import ictv_db
        if ictv_db.ensure_taxa(logger=logger):
            ictv_meta = ictv_db.load_taxa()
            gb_fa = ictv_db._refs_fa_path()
            if gb_fa and os.path.isfile(gb_fa):
                ref_sources.append(gb_fa)
            ictv_gb_set = ictv_db._gb_cached_accs()
            if logger:
                logger.log(f"ICTV VMR（{ictv_db._read_version().get('MSL', '?')}）"
                           f"谱系接入: {len(ictv_meta)} 条")
    except Exception as e:
        if logger:
            logger.log(f"ICTV VMR 库不可用（跳过）: {e}", "WARN")

    def _ref_pref(sacc):
        """参考优选键：(RefSeq, 完整基因组)。bitscore 由调用方拼接。"""
        m = ref_meta.get(sacc) or {}
        lin = ictv_meta.get(sacc.split('.')[0].upper()) or {}
        complete = ((m.get('Nuc_Completeness') or '') == 'complete'
                    or bool(lin))          # all_virus/VMR 收录完整基因组为主
        return (1 if (m.get('Sequence_Type') or '') == 'RefSeq' else 0,
                1 if complete else 0)

    # 统一参考分类索引 {base: (species, genus, family)}：植物库元数据 +
    # ICTV VMR（优先，MSL 当前版）+ acvirus taxa.txt —— 层级抽样
    # （macro/genus/lineage）的数据基础
    ref_tax = {}
    for acc, m in ref_meta.items():
        fam = (m.get('VMR_Family') or '').strip()
        gen = (m.get('VMR_Genus') or '').strip()
        sp = (m.get('Species') or m.get('VMR_Species') or '').strip()
        if fam or gen or sp:
            ref_tax[acc.split('.')[0].upper()] = (sp, gen, fam)
    for base, row in ictv_meta.items():
        if base not in ref_tax:
            fam = (row.get('Family') or '').strip()
            gen = (row.get('Genus') or '').strip()
            sp = (row.get('Species') or '').strip()
            if fam or sp:
                ref_tax[base] = (sp, gen, fam)

    genus_accs = defaultdict(list)
    family_genera = defaultdict(set)
    for base, (_sp, gen, fam) in ref_tax.items():
        if gen:
            genus_accs[gen.lower()].append(base)
            if fam:
                family_genera[fam.lower()].add(gen.lower())

    # 按 kunpeng 分类的 (species, genus, family) 分组
    # 优先用 species，无则 genus，无则 family
    for cid in viral_ids:
        info = cls_map.get(cid, {})
        sp = (info.get('species') or '').strip()
        gen = (info.get('genus') or '').strip()
        fam = (info.get('family') or '').strip()
        # group key: 有物种名用物种，否则用 genus，否则用 family，否则 unknown
        if sp:
            top = sp
        elif gen:
            top = gen
        elif fam:
            top = fam
        else:
            top = 'unknown'
        groups.setdefault(top, []).append(cid)
    # 组按总长度排序，最多 max_groups 组
    group_items = sorted(groups.items(),
                         key=lambda kv: -sum(len(contig_seqs.get(c, '')) for c in kv[1]))
    group_items = group_items[:max_groups]

    results = []
    for top_acc, cids in group_items:
        gname = _safe_name(top_acc or 'unknown')
        gdir = check_path(os.path.join(out_dir, gname), must_exist=False, in_platform=True)
        os.makedirs(gdir, exist_ok=True)
        if logger:
            logger.log(f"组 {gname}: contigs {len(cids)} 条")

        # 参考挑选：按 (family, genus, species) 从参考库抽样
        # cls_map 提供每组 contig 的分类信息，ref_tax 提供参考序列的分类索引
        if sampling != 'blast':
            ref_bases = _tax_sample_refs(
                top_acc, ref_tax, genus_accs, family_genera, sampling,
                cap=max(top_n_refs, 10), other_count=other_count)
            if ref_bases is None and logger:
                logger.log(f"  组 {top_acc} 无分类信息，该组使用全部参考")
                ref_bases = []
        else:
            # cls_map 有 species 信息，用它来选参考；否则用全部参考
            ref_bases = []
            seen = set()
            for cid in cids:
                info = cls_map.get(cid, {})
                sp = (info.get('species') or '').strip()
                if sp and sp not in seen:
                    # 从 ref_tax 找该 species 的 accession
                    for base, (_s, _g, _f) in ref_tax.items():
                        if _s.lower() == sp.lower() and base not in seen:
                            ref_bases.append(base)
                            seen.add(base)
                            if len(ref_bases) >= max(top_n_refs, 10):
                                break
                if len(ref_bases) >= max(top_n_refs, 10):
                    break
            if not ref_bases:
                ref_bases = []
        if not ref_bases:
            # 兜底：用所有参考序列
            ref_ids = []
            seen = set()
            for base in ref_tax:
                if base not in seen:
                    seen.add(base)
                    ref_ids.append(base)
                    if len(ref_ids) >= max(top_n_refs * 3, 30):
                        break
        else:
            ref_ids = ref_bases

        refs_fa = os.path.join(gdir, 'refs.fa')
        found_refs = _extract_refs(ref_sources, ref_ids, refs_fa)

        # 兜底：若有 ref_ids 未匹配到本地序列，尝试从 NCBI 按需下载
        not_found = [rid for rid in ref_ids if rid not in found_refs]
        if not_found:
            if logger:
                logger.log(f"  参考序列缺失 {len(not_found)} 条，尝试 NCBI 下载...")
            try:
                from . import ictv_db
                result = ictv_db.ensure_gb(not_found, logger=logger)
                downloaded = result.get('downloaded', 0)
                if downloaded > 0 and logger:
                    logger.log(f"  NCBI 下载 {downloaded} 条参考序列")
                # 重新提取（gb_refs.fa 已更新）
                found_refs = _extract_refs(ref_sources, ref_ids, refs_fa)
            except Exception as e:
                if logger:
                    logger.log(f"  NCBI 下载失败（跳过）: {e}", "WARN")

        # 参考谱系标注（植物库元数据 > ICTV VMR 当前版 > acvirus taxa.txt）
        refs_info = []
        for sacc in found_refs:
            m = ref_meta.get(sacc) or {}
            base = sacc.split('.')[0].upper()
            ilin = ictv_meta.get(base) or {}
            if m:
                src = 'plant_ref'
            elif base in ictv_gb_set:
                src = 'ictv_gb'      # 序列来自 ICTV VMR 按需下载缓存
            else:
                src = 'ictv_taxa'
            refs_info.append({
                'accession': sacc,
                'source': src,
                'species': m.get('Species', '') or ilin.get('Species', ''),
                'genus': m.get('VMR_Genus', '') or ilin.get('Genus', ''),
                'family': m.get('VMR_Family', '') or ilin.get('Family', ''),
                'completeness': m.get('Nuc_Completeness', '')
                                or ilin.get('Genome_Coverage', ''),
            })
        with safe_open(os.path.join(gdir, 'refs.tsv'), 'wt') as f:
            f.write('accession\tsource\tspecies\tgenus\tfamily\tcompleteness\n')
            for r in refs_info:
                f.write('\t'.join(r[k] for k in
                                  ('accession', 'source', 'species', 'genus',
                                   'family', 'completeness')) + '\n')
        combined = os.path.join(gdir, 'combined.fa')
        with safe_open(combined, 'wt') as f:
            for cid in cids:
                write_fasta_record(f, cid, contig_seqs[cid])
            for h, s in iter_fasta(refs_fa):
                write_fasta_record(f, h.split()[0], s)
            # 额外参考（NCBI 下载集合等）：按 ID 去重后追加
            seen_ids = set(cids) | set(found_refs)
            n_extra = 0
            for ef in extra_refs or []:
                for h, s in iter_fasta(check_path(ef, must_exist=True)):
                    sid = h.split()[0]
                    base_id = sid.split('|')[0]
                    if base_id in seen_ids or sid in seen_ids:
                        continue
                    write_fasta_record(f, sid, s)
                    seen_ids.add(sid)
                    seen_ids.add(base_id)
                    found_refs.append(base_id)
                    n_extra += 1
        if logger:
            logger.log(f"  组合序列: contigs {len(cids)} + 参考 {len(found_refs)}"
                       + (f" (含额外参考 {n_extra})" if n_extra else ""))

        aln = os.path.join(gdir, 'aln.fasta')
        _run_mafft(combined, aln, threads=threads, logger=logger)

        # trimAl 清剪（可选；失败/过度修剪自动回退原比对）
        if do_trim:
            aln_used, trim_info = _run_trimal(
                aln, os.path.join(gdir, 'aln.trim.fasta'), logger=logger)
        else:
            if logger:
                logger.log("  已关闭 trimAl 清剪，使用原比对建树")
            aln_used, trim_info = aln, {'applied': False}

        tree_file = None
        tree_extra = {}
        try:
            if tree_tool == 'iqtree':
                tree_file = _run_iqtree(aln_used, os.path.join(gdir, 'iqtree'),
                                        threads=threads, logger=logger)
                tree_extra = _parse_iqtree_log(os.path.join(gdir, 'iqtree.iqtree'))
            elif tree_tool == 'nj':
                tree_file = _run_nj(aln_used, os.path.join(gdir, 'nj.nwk'),
                                    logger=logger)
            else:
                tree_file = _run_fasttree(aln_used, os.path.join(gdir, 'tree.nwk'),
                                          logger=logger)
        except Exception as e:
            if logger:
                logger.log(f"  建树失败({e})，保留比对与 SDT 结果", "WARN")

        csv_path, sdt_fas, names, mat = write_sdt_outputs(aln, gdir)
        # 代表 contig 与最近参考的 identity
        best = 0.0
        for i, nm in enumerate(names):
            if nm in cids:
                for j, nm2 in enumerate(names):
                    if nm2 in found_refs:
                        best = max(best, mat[i][j])
        results.append({'group': gname, 'contigs': cids, 'n_refs': len(found_refs),
                        'refs_species': [r['species'] for r in refs_info
                                         if r['species']][:5],
                        'dir': gname, 'aln': 'aln.fasta',
                        'trim': trim_info,
                        'tree': os.path.basename(tree_file) if tree_file else None,
                        'tree_model': tree_extra.get('model'),
                        'tree_logl': tree_extra.get('logl'),
                        'sdt_csv': 'sdt_matrix.csv', 'sdt_fas': 'sdt_input.fas',
                        'best_identity_to_ref': round(best, 2)})
        if logger:
            logger.log(f"  完成: 与最近参考 identity {best:.1f}%")

    summary = {'stage': step, 'skipped': False, 'tree_tool': tree_tool,
               'trim_tool': 'trimAl-automated1' if any(
                   g.get('trim', {}).get('applied') for g in results) else None,
               'groups': results}
    with safe_open(summary_file, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    mark_step_done(out_dir, step)
    if logger:
        logger.log(f"阶段⑤ 完成: {len(results)} 个病毒组")
    return summary

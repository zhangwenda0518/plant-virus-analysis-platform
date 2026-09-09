# -*- coding: utf-8 -*-
"""
SDT 精确引擎（SDTv1.3 的纯 Python 复刻，替代外部 SDT exe）。

移植自 MMPV-RNA virome_analysis_pipeline（sdt_genus_matrix.py +
virus_auto_pipeline.py 的 build_mat_sdt_exact / plot_matrix /
plot_distribution / get_safe_leaf_order），口径与 SDT 原版
（Brejnev & Muhire 2014）一致：

  每一对序列 → MAFFT 独立全局比对（--localpair，L-INS-i 精度）
  → SDT Get_Similarity 公式：similarity = 100 × (1 - dist / (aln_len - gaps))
    （gaps = 任一序列为 gap 的列数，即分母只统计两序列均有碱基的列）

  ProcessPoolExecutor 多进程并行；.resume_cache 原子刷盘断点续传；
  层级聚类排序（average linkage）后绘 SDT 经典红黄色阶热图
  （80/85/90/95/100 硬分界，三角模式）+ identity 分布曲线，
  PNG(300dpi)+PDF(fonttype 42) 双输出。
"""
import os
import hashlib
import atexit
import warnings
from pathlib import Path

import numpy as np

from .utils import safe_open

# matplotlib Agg（与 landscape_plot 同口径）
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

mpl = matplotlib
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# SDT 经典红黄色阶（10 档）+ 常用学术色阶（按 matplotlib cmap 采样 10 档）
SDT_COLORS = [(0.49, 0.02, 0.05), (0.66, 0.11, 0.05), (0.80, 0.24, 0.06),
              (0.90, 0.42, 0.09), (0.96, 0.60, 0.17), (0.99, 0.75, 0.32),
              (1.00, 0.86, 0.51), (1.00, 0.93, 0.70), (1.00, 0.97, 0.86),
              (1.00, 1.00, 0.94)]
PLOT_DPI = 300

PALETTE_NAMES = ['sdt', 'cividis', 'viridis', 'RdYlBu', 'Spectral',
                 'YlGnBu', 'coolwarm', 'magma']


def _cmap_colors(name):
    """任意 matplotlib colormap → 10 档 RGB 列表（供 BoundaryNorm 硬分界）。"""
    try:
        cm = plt.get_cmap(name)
        return [tuple(c[:3]) for c in cm(np.linspace(0, 1, 10))]
    except (ValueError, TypeError):
        return list(SDT_COLORS)

K_ORIENT = 21
_COMP = str.maketrans('ACGTNacgtnRYKMSWBDHVrykmswbdhv',
                      'TGCANtgcanYRMKSWVHDByrmkswvhdb')


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def detect_seqtype(seqs):
    """核酸/蛋白自动判别（MMPV 双轨口径）：统计蛋白特征残基
    E/F/I/L/P/Q 占比——这些字母不出现在核酸歧义码中，
    蛋白序列里通常 >20%，核酸序列里 ≈0%。"""
    prot_markers = 'EFILPQ'
    total = marker = 0
    for s in seqs[:5]:
        t = str(s).upper().replace('-', '')
        total += len(t)
        marker += sum(t.count(c) for c in prot_markers)
    return 'aa' if total and marker / total > 0.03 else 'nt'


def sha_seq(s):
    return hashlib.sha1(s.upper().encode()).hexdigest()


def n_percent(seq):
    seq = seq.upper()
    return 100.0 * seq.count('N') / len(seq) if len(seq) else 100.0


def kmer_set(seq, k=K_ORIENT):
    s = seq.upper()
    ok = set('ACGT')
    return {s[i:i + k] for i in range(len(s) - k + 1)
            if not (set(s[i:i + k]) - ok)}


def auto_orient(recs, ref_seqs=None, k=K_ORIENT, logger=None):
    """21-mer 共享计数自动定向（MMPV 口径）：以最长序列为全局方向锚，
    其余序列反向共享显著高于正向（>1.1×）时翻正。recs = [(name, seq)]（就地修改）。

    锚必须固定为单条：若以所有序列 k-mer 的并集为锚，反向序列自身贡献的
    k-mer 恒在锚集内（fwd≈rev），反向永远判不出。ref_seqs 提供时取其中
    最长者作锚（如外部完整基因组参考）；否则取 recs 内最长者。
    """
    if not recs:
        return {}
    if ref_seqs:
        anchor_seq = max(ref_seqs, key=len)
    else:
        anchor_seq = max((seq for _, seq in recs), key=len)
    anchor_seq = str(anchor_seq).upper()
    anchor_kmers = kmer_set(anchor_seq)
    flips = {}
    for item in recs:
        name, seq = item[0], item[1]
        if seq == anchor_seq:               # 锚自身保持原方向
            flips[name] = False
            continue
        ks = kmer_set(seq)
        if len(ks) < 5:
            flips[name] = False
            continue
        fwd = len(ks & anchor_kmers)
        rc = kmer_set(revcomp(seq))
        rev = len(rc & anchor_kmers)
        flips[name] = rev > fwd * 1.1
        if flips[name]:
            item[1] = revcomp(seq)
            if logger:
                logger.log(f"  [定向] {name} 为反向互补，已翻正 "
                           f"(fwd={fwd}, rev={rev})")
    return flips


# ------------------------------------------------------------------ worker --
_W = {}          # 子进程全局：{'mafft': 路径, 'seqs': 序列列表}


def _w_init(mafft, seqs):
    _W['mafft'] = mafft
    _W['seqs'] = seqs


def _sdt_formula(a, b):
    """SDT Get_Similarity：分母 = 非双 gap 列（含两端修剪语义一致）。"""
    dist = gaps = 0
    for x, y in zip(a, b):
        if x == '-' or y == '-':
            gaps += 1
        elif x != y:
            dist += 1
    denom = min(len(a), len(b)) - gaps
    return (1.0 - dist / denom) * 100.0 if denom > 0 else float('nan')


def _msa_pair_worker(task):
    """已比对输入：直接从 MSA 两行按 SDT 公式算（不调比对器）。"""
    i, j = task
    a, b = _W['seqs'][i], _W['seqs'][j]
    if not a or not b:
        return i, j, float('nan')
    try:
        return i, j, _sdt_formula(a.upper(), b.upper())
    except Exception:
        return i, j, float('nan')


def _sdt_pair_worker(task):
    """一对序列：写临时 FASTA → mafft --localpair 独立比对 → Get_Similarity。"""
    import subprocess as sp
    import tempfile
    i, j = task
    s1, s2 = _W['seqs'][i], _W['seqs'][j]
    if not s1 or not s2:
        return i, j, float('nan')
    tmp = None
    try:
        tmp = tempfile.mkdtemp(prefix='sdt_pair_')
        in_f = os.path.join(tmp, 'pair.fasta')
        with open(in_f, 'wt') as f:
            f.write(f'>s{i}\n')
            for k in range(0, len(s1), 70):
                f.write(s1[k:k + 70] + '\n')
            f.write(f'>s{j}\n')
            for k in range(0, len(s2), 70):
                f.write(s2[k:k + 70] + '\n')
        # mafft 输出到 stdout（--quiet；--localpair = L-INS-i 精度，
        # 与 SDT v1.3 默认 MAFFT 路径一致）
        r = sp.run([_W['mafft'], '--quiet', '--localpair', in_f],
                   capture_output=True, text=True, timeout=1800)
        if r.returncode != 0:
            return i, j, float('nan')

        # 解析两条比对序列（命名块循环逐行拼接）
        seqs, order = {}, []
        cur = None
        for line in r.stdout.splitlines():
            if line.startswith('>'):
                cur = line[1:].strip().split()[0]
                order.append(cur)
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
        if len(order) < 2:
            return i, j, float('nan')
        a = ''.join(seqs[order[0]]).upper()
        b = ''.join(seqs[order[1]]).upper()
        return i, j, _sdt_formula(a, b)
    except Exception:
        return i, j, float('nan')
    finally:
        if tmp:
            try:
                import shutil
                shutil.rmtree(tmp, ignore_errors=True)
            except Exception:
                pass


# ------------------------------------------------------------------ engine --
def build_mat_sdt_exact(names, seqs, mafft, threads=4, out_dir=None,
                        cache_prefix='sdt', resume=True,
                        logger=None, progress=None, cancel=None,
                        aligned=False):
    """逐对 SDT 相似度矩阵（n×n，上三角有效）。

    aligned=False：每对 MAFFT 独立比对（SDT v1.3 原路径）；
    aligned=True ：输入已比对（各行等长），直接按公式计算，不调比对器。
    线程 = 并行进程数；progress(frac 0~1, msg)；cancel=threading.Event。
    out_dir 下 .resume_cache/<prefix>_<n>_<hash>.npy 原子刷盘续传。
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    n = len(seqs)
    mat = np.full((n, n), np.nan)
    cache_dir = Path(out_dir) / '.resume_cache' if out_dir else None
    cache_file = None
    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)
        h = sha_seq('|'.join(seqs))[:12]
        cache_file = cache_dir / f'{cache_prefix}_{n}_{h}.npy'
        if resume and cache_file.exists():
            try:
                cached = np.load(cache_file)
                if cached.shape == (n, n):
                    if logger:
                        logger.log('  [续传] 加载 SDT 比对缓存，'
                                   '已完成对位直接跳过')
                    mat = cached
            except Exception:
                mat = np.full((n, n), np.nan)

    total_pairs = n * (n - 1) // 2
    done_pairs = sum(1 for i in range(n) for j in range(i + 1, n)
                     if not np.isnan(mat[i, j]))
    if logger:
        logger.log(f'SDT 精确计算：{n} 条序列，'
                   f'{total_pairs - done_pairs}/{total_pairs} 对待比对'
                   f'（MAFFT 逐对独立比对 + Get_Similarity）')
    if done_pairs >= total_pairs:
        return mat

    tasks = [(i, j) for i in range(n) for j in range(i + 1, n)
             if np.isnan(mat[i, j])]

    n_proc = max(1, min(int(threads or 4), len(tasks) or 1))
    worker = _msa_pair_worker if aligned else _sdt_pair_worker
    if not aligned:
        # 快速失败：比对器不可用时直接报错，避免全 NaN 结果污染续传缓存
        import shutil as _shutil
        ok = bool(mafft) and (os.path.isfile(mafft) or bool(_shutil.which(mafft)))
        if not ok:
            raise RuntimeError(f'未找到 MAFFT 比对器: {mafft}——'
                               '请到设置页检查工具路径')
    if aligned:
        _w_init(None, seqs)          # 主进程兜底；子进程仍经 initializer
        ex = ProcessPoolExecutor(max_workers=n_proc,
                                 initializer=_w_init,
                                 initargs=(None, seqs))
    else:
        ex = ProcessPoolExecutor(max_workers=n_proc,
                                 initializer=_w_init, initargs=(mafft, seqs))
    counter = 0
    n_ok = 0
    flush = max(500, len(tasks) // 20)
    try:
        futs = [ex.submit(worker, t) for t in tasks]
        for fut in as_completed(futs):
            if cancel is not None and cancel.is_set():
                raise RuntimeError('已取消')
            i, j, sim = fut.result()
            mat[i, j] = sim
            if not np.isnan(sim):
                n_ok += 1
            counter += 1
            if progress:
                progress((done_pairs + counter) / max(total_pairs, 1),
                         f'SDT 逐对比对 {done_pairs + counter}'
                         f'/{total_pairs} 对')
            if cache_file and counter % flush == 0:
                _atomic_save(cache_file, mat)
    finally:
        ex.shutdown(wait=False, cancel_futures=True)
    if not aligned and tasks and n_ok == 0:
        raise RuntimeError('所有成对比对均失败——请检查 MAFFT 可用性与序列'
                           '格式（本结果未写入缓存）')
    if cache_file:
        _atomic_save(cache_file, mat)
    return mat


def _atomic_save(cache_file, mat):
    try:
        tmp = cache_file.parent / (cache_file.stem + '.tmp.npy')
        np.save(tmp, mat)
        tmp.replace(cache_file)
    except Exception:
        pass


def cluster_order(mat, n):
    """层级聚类排序（average linkage，NaN 按 100 距离兜底）。"""
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
        if np.isnan(mat).all() or n < 3:
            return list(range(n))
        dist = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                dist[i, j] = 0.0 if i == j else (
                    100.0 if np.isnan(mat[i, j])
                    else max(0.0, 100.0 - mat[i, j]))
        z = linkage(squareform(dist), method='average')
        return leaves_list(z).tolist()
    except Exception as e:
        warnings.warn(f'聚类排序失败，退回原始顺序: {e}')
        return list(range(n))


# ------------------------------------------------------------------ plots --
def plot_heatmap(mat, labels, title, out_base, triangle_only=True,
                 colors=None):
    """色阶硬分界热图（层级聚类排序 + 三角模式；colors=10 档 RGB 列表）。"""
    n = len(labels)
    fig_size = max(8.0, min(n * 0.35, 40.0))
    font_size = min(14, max(6, 250 / max(n, 1)))

    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    ax.set_facecolor('white')

    if triangle_only:
        mask = np.triu(np.ones_like(mat, dtype=bool), k=0) | np.isnan(mat)
    else:
        mask = np.eye(n, dtype=bool) | np.isnan(mat)

    cmap = ListedColormap(colors or SDT_COLORS)
    cmap.set_bad('white', alpha=0.0)
    vals = mat[~mask]
    if len(vals) > 0:
        amin = np.nanmin(vals)
        if amin >= 95.0:
            min_v = 95.0
        elif amin >= 90.0:
            min_v = 90.0
        elif amin >= 85.0:
            min_v = 85.0
        elif amin >= 80.0:
            min_v = 80.0
        else:
            min_v = np.floor(amin / 10.0) * 10.0
    else:
        min_v = 80.0
    bounds = np.linspace(min_v, 100.0, len(SDT_COLORS) + 1)
    norm = BoundaryNorm(bounds, cmap.N)

    show = np.ma.array(mat, mask=mask)
    ax.imshow(show, cmap=cmap, norm=norm, interpolation='nearest')
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=90, fontsize=font_size)
    ax.set_yticklabels(labels, fontsize=font_size)
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)

    lw = 0.5 if n < 100 else 0.1
    for i in range(n + 1):
        if triangle_only:
            ax.plot([0, min(i, n)], [i, i], color='black', lw=lw)
            ax.plot([i, i], [i, n], color='black', lw=lw)
        else:
            ax.axhline(i, color='black', lw=lw)
            ax.axvline(i, color='black', lw=lw)

    ax.set_title(title, fontsize=max(16, font_size * 1.2),
                 fontweight='bold', pad=20)
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                        shrink=0.5, pad=0.02)
    cbar.set_label('Pairwise Identity (%)')

    fig.savefig(f'{out_base}.png', dpi=PLOT_DPI, bbox_inches='tight',
                facecolor='white')
    fig.savefig(f'{out_base}.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)


def plot_distribution(mat, title, out_base, triangle_only=True):
    """identity 分布密度曲线（SDT 原版风格，红色粗线）。"""
    plt.figure(figsize=(8, 6))
    n = mat.shape[0]
    if triangle_only:
        mask = np.triu(np.ones_like(mat, dtype=bool), k=0) | np.isnan(mat)
    else:
        mask = np.eye(n, dtype=bool) | np.isnan(mat)
    vals = mat[~mask]
    if len(vals) > 0:
        min_v = max(30.0, float(np.floor(np.nanmin(vals))))
        counts, edges = np.histogram(vals, bins=np.arange(min_v, 101, 1))
        plt.plot(0.5 * (edges[:-1] + edges[1:]), counts / len(vals),
                 color='#CC0000', lw=3.0)
        plt.ylim(0, max(counts / len(vals)) * 1.1)
    else:
        min_v = 80.0
    plt.title(f'Distribution - {title}', fontsize=16, fontweight='bold',
              pad=15)
    plt.xlabel('Identity (%)', fontsize=14)
    plt.ylabel('Proportion', fontsize=14)
    plt.xlim(min_v, 100.0)
    plt.grid(axis='y', linestyle='--', alpha=0.7, color='#CCCCCC')
    plt.savefig(f'{out_base}.png', dpi=PLOT_DPI, bbox_inches='tight',
                facecolor='white')
    plt.savefig(f'{out_base}.pdf', bbox_inches='tight', facecolor='white')
    plt.close('all')


def write_matrix_csv(names, mat, path):
    """n×n 对称矩阵 CSV（NaN 写 NaN，与 SDT 导出一致）。"""
    n = len(names)
    full = np.full((n, n), np.nan)
    for i in range(n):
        full[i, i] = 100.0
        for j in range(n):
            if i != j:
                full[i, j] = mat[i, j] if i < j else mat[j, i]
    with safe_open(path, 'wt') as f:
        f.write(',' + ','.join(names) + '\n')
        for i, nm in enumerate(names):
            row = [nm]
            for j in range(n):
                v = full[i, j]
                row.append('' if i == j else
                           ('NaN' if np.isnan(v) else f'{v:.1f}'))
            f.write(','.join(row) + '\n')
    return full


def run_sdt_exact(in_fasta, run_dir, mafft, max_n=30, orient=True,
                  threads=4, palette='sdt', logger=None, progress=None,
                  cancel=None, aligned=False, seqtype='auto'):
    """t-sdt 卡片主入口：FASTA → SDT 精确矩阵 + 聚类排序热图 + 分布图。

    aligned=True 时输入为已比对 MSA（跳过比对器与自动定向）；
    seqtype='auto' 自动判别核酸/蛋白（MMPV 双轨口径），蛋白走 AA 同一性。
    返回 dict（names / n / pairs / csv / heatmap / distribution / json）。
    """
    import json as _json
    from .utils import iter_fasta

    recs = []
    seen = {}
    for h, s in iter_fasta(in_fasta):
        name = (h or '').strip().split()[0][:40] or f'seq{len(recs) + 1}'
        name = ''.join(c if c.isalnum() or c in '_-.' else '_' for c in name)
        if name in seen:
            seen[name] += 1
            name = f'{name}_{seen[name]}'
        else:
            seen[name] = 0
        recs.append([name, s.upper()])
    if len(recs) < 2:
        raise RuntimeError('FASTA 中少于 2 条序列，无法做 SDT 分析')
    if len(recs) > max_n:
        if logger:
            logger.log(f'序列数 {len(recs)} 超过上限，截取前 {max_n} 条')
        recs = recs[:max_n]
    n = len(recs)

    if progress:
        progress(0.03, f'读取 {n} 条序列')
    st = (seqtype or 'auto').lower()
    if st not in ('nt', 'aa'):
        st = detect_seqtype([s for _, s in recs])
        if logger:
            logger.log(f'序列类型自动判别：{st.upper()}')
    if st == 'aa':
        orient = False
        if logger:
            logger.log('蛋白输入：按 AA 同一性计算（跳过 21-mer 定向）')
    if aligned:
        lens = {len(s) for _, s in recs}
        if len(lens) > 1:
            raise RuntimeError('勾选了「已比对输入」，但序列长度不一致——'
                               '请确认输入是完整的多序列比对（MSA）FASTA')
        if logger:
            logger.log('输入为已比对 MSA：跳过比对器，直接按 SDT 公式计算')
    elif orient and n >= 2:
        if progress:
            progress(0.06, '21-mer 自动定向')
        flips = auto_orient(recs, [s for _, s in recs], logger=logger)
        n_flip = sum(1 for v in flips.values() if v)
        if logger and not n_flip:
            logger.log('21-mer 定向：全部序列方向一致，无需翻正')

    names = [r[0] for r in recs]
    seqs = [r[1] for r in recs]

    if progress:
        progress(0.1, 'SDT 逐对精确比对')
    mat = build_mat_sdt_exact(names, seqs, mafft, threads=threads,
                              out_dir=run_dir, cache_prefix='sdt',
                              logger=logger, progress=lambda p, m:
                              progress(0.1 + p * 0.8, m) if progress else None,
                              cancel=cancel, aligned=aligned)

    if progress:
        progress(0.92, '写矩阵 + 聚类排序 + 出图')
    csv_path = os.path.join(run_dir, 'sdt_matrix.csv')
    full = write_matrix_csv(names, mat, csv_path)

    order = cluster_order(full, n)
    names_o = [names[i] for i in order]
    mat_o = full[np.ix_(order, order)]

    title = f'SDT Exact: {n} sequences ({st.upper()} Identity %)'
    plot_heatmap(mat_o, names_o, title, os.path.join(run_dir, 'sdt_heatmap'),
                 triangle_only=True, colors=_cmap_colors(palette))
    plot_distribution(mat_o, title, os.path.join(run_dir, 'sdt_distribution'),
                      triangle_only=True)

    data = {'names': names, 'order': order, 'n': n,
            'pairs': n * (n - 1) // 2,
            'matrix': [[None if np.isnan(full[i, j]) else round(full[i, j], 2)
                        for j in range(n)] for i in range(n)],
            'palette': palette, 'aligned': aligned, 'seqtype': st}
    js = os.path.join(run_dir, 'sdt_matrix.json')
    with safe_open(js, 'wt') as f:
        _json.dump(data, f, ensure_ascii=False)

    if progress:
        progress(1.0, '完成')
    return {'n': n, 'pairs': n * (n - 1) // 2, 'names': names,
            'seqtype': st,
            'csv': 'sdt_matrix.csv', 'heatmap': 'sdt_heatmap.png',
            'heatmap_pdf': 'sdt_heatmap.pdf',
            'distribution': 'sdt_distribution.png',
            'distribution_pdf': 'sdt_distribution.pdf',
            'json': 'sdt_matrix.json'}


# ------------------------------------------------------------------ 
# 核苷酸 + 氨基酸同一性表（BioAider「Sequence Identity Matrix」口径：
# 同一性 = 两序列比对后相同位点 / 可比较位点；NT 矩阵 + AA 矩阵 +
# 复合展示（NT 上三角 / AA 下三角）+ 逐对同一性长表）
# ------------------------------------------------------------------
def _read_fasta_recs(in_fasta, max_n, logger=None):
    from .utils import iter_fasta
    recs, seen = [], {}
    for h, s in iter_fasta(in_fasta):
        name = (h or '').strip().split()[0][:40] or f'seq{len(recs) + 1}'
        name = ''.join(c if c.isalnum() or c in '_-.' else '_' for c in name)
        if name in seen:
            seen[name] += 1
            name = f'{name}_{seen[name]}'
        else:
            seen[name] = 0
        recs.append([name, s.upper()])
    if len(recs) > max_n:
        if logger:
            logger.log(f'序列数 {len(recs)} 超过上限，截取前 {max_n} 条')
        recs = recs[:max_n]
    return recs


def plot_heatmap_composite(nt_full, aa_full, labels, title, out_base,
                           colors=None):
    """复合热图：上三角 = 核苷酸同一性，下三角 = 氨基酸同一性。"""
    n = len(labels)
    combined = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if i == j:
                combined[i, j] = 100.0
            elif i < j:
                combined[i, j] = nt_full[i, j]
            else:
                combined[i, j] = aa_full[i, j]
    mask = np.isnan(combined)
    cmap = ListedColormap(colors or SDT_COLORS)
    cmap.set_bad('white', alpha=0.0)
    vals = combined[~np.isnan(combined) &
                    ~np.eye(n, dtype=bool)]
    min_v = 80.0
    if len(vals) > 0:
        amin = float(np.nanmin(vals))
        min_v = (95.0 if amin >= 95 else 90.0 if amin >= 90 else
                 85.0 if amin >= 85 else 80.0 if amin >= 80 else
                 float(np.floor(amin / 10) * 10))
    bounds = np.linspace(min_v, 100.0, len(cmap.colors) + 1)
    norm = BoundaryNorm(bounds, cmap.N)

    fig_size = max(8.0, min(n * 0.35, 40.0))
    font_size = min(14, max(6, 250 / max(n, 1)))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    ax.set_facecolor('white')
    ax.imshow(np.ma.array(combined, mask=mask), cmap=cmap, norm=norm,
              interpolation='nearest')
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=90, fontsize=font_size)
    ax.set_yticklabels(labels, fontsize=font_size)
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    lw = 0.5 if n < 100 else 0.1
    for i in range(n + 1):
        ax.axhline(i, color='black', lw=lw)
        ax.axvline(i, color='black', lw=lw)
    ax.set_title(title, fontsize=max(16, font_size * 1.2),
                 fontweight='bold', pad=20)
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                        shrink=0.5, pad=0.02)
    cbar.set_label('Identity (%)  ·  Upper: NT / Lower: AA')
    fig.savefig(f'{out_base}.png', dpi=PLOT_DPI, bbox_inches='tight',
                facecolor='white')
    fig.savefig(f'{out_base}.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)


def run_identity_table(nt_fasta, run_dir, mafft, aa_fasta=None, max_n=30,
                       aligned=False, threads=4, palette='sdt', min_aa=50,
                       logger=None, progress=None, cancel=None):
    """核苷酸+氨基酸同一性表主入口。

    aa_fasta 提供时按名称匹配 AA 序列；否则 6-frame 最长 ORF 翻译。
    返回 dict（产物文件名 + n + n_aa + pairs）。
    """
    import json as _json
    from .contig_annot import longest_orf_protein

    if progress:
        progress(0.03, '读取核苷酸序列')
    nt_recs = _read_fasta_recs(nt_fasta, max_n, logger=logger)
    if len(nt_recs) < 2:
        raise RuntimeError('FASTA 中少于 2 条序列，无法做同一性分析')
    nt_names = [r[0] for r in nt_recs]
    nt_seqs = [r[1] for r in nt_recs]
    n = len(nt_names)
    if detect_seqtype(nt_seqs) == 'aa':
        raise RuntimeError('检测到输入是蛋白序列——请切回「AA · 蛋白同一性'
                           '矩阵」模式（NT+AA 模式的①需要核酸输入）')

    if aligned:
        lens = {len(s) for s in nt_seqs}
        if len(lens) > 1:
            raise RuntimeError('勾选了「已比对输入」，但核苷酸序列长度不一致'
                               '——请确认输入是完整 MSA FASTA')
        if logger:
            logger.log('核苷酸输入为已比对 MSA：跳过比对器直接计算')

    # ---- AA 序列来源 ----
    aa_map = {}
    if aa_fasta:
        if progress:
            progress(0.08, '读取氨基酸序列')
        for name, s in _read_fasta_recs(aa_fasta, max_n * 4, logger=logger):
            aa_map[name] = s
        missing = [x for x in nt_names if x not in aa_map]
        if logger:
            logger.log(f'AA 输入按名称匹配：命中 {n - len(missing)}/{n}'
                       + (f'，缺失 {len(missing)} 条按最长 ORF 翻译兜底'
                          if missing else ''))
    if progress:
        progress(0.1, '生成氨基酸序列' +
                 ('（6-frame 最长 ORF 翻译）' if not aa_fasta else ''))
    n_translated = 0
    for name, s in zip(nt_names, nt_seqs):
        if name not in aa_map or not aa_map[name].strip():
            p = longest_orf_protein(s, min_aa=min_aa)
            aa_map[name] = p
            if p:
                n_translated += 1
    aa_names = [x for x in nt_names if aa_map.get(x)]
    aa_seqs = [aa_map[x] for x in aa_names]
    n_aa = len(aa_names)
    if logger:
        logger.log(f'AA 序列就绪：{n_aa}/{n} 条'
                   f'（其中 {n_translated} 条来自最长 ORF 翻译）')
    if n_aa < 2:
        logger and logger.log('有效 AA 序列不足 2 条，仅输出核苷酸矩阵')

    # ---- NT 矩阵 ----
    if progress:
        progress(0.15, '核苷酸同一性矩阵')
    mat_nt = build_mat_sdt_exact(nt_names, nt_seqs, mafft, threads=threads,
                                 out_dir=run_dir, cache_prefix='idnt',
                                 logger=logger,
                                 progress=lambda p, m: progress(
                                     0.15 + p * 0.35,
                                     f'NT：{m}') if progress else None,
                                 cancel=cancel, aligned=aligned)

    # ---- AA 矩阵 ----
    mat_aa = None
    if n_aa >= 2:
        if progress:
            progress(0.55, '氨基酸同一性矩阵')
        aa_aligned = aligned   # NT 已比对时 AA 未必等长，下面校验
        if aa_aligned:
            aa_lens = {len(s) for s in aa_seqs}
            aa_aligned = len(aa_lens) == 1
        mat_aa = build_mat_sdt_exact(aa_names, aa_seqs, mafft, threads=threads,
                                     out_dir=run_dir, cache_prefix='idaa',
                                     logger=logger,
                                     progress=lambda p, m: progress(
                                         0.55 + p * 0.3,
                                         f'AA：{m}') if progress else None,
                                     cancel=cancel, aligned=aa_aligned)

    if progress:
        progress(0.9, '写同一性表 + 出图')
    full_nt = write_matrix_csv(nt_names, mat_nt,
                               os.path.join(run_dir, 'nt_matrix.csv'))
    if mat_aa is not None:
        write_matrix_csv(aa_names, mat_aa,
                         os.path.join(run_dir, 'aa_matrix.csv'))

    # BioAider 口径逐对同一性长表
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            nt_v = mat_nt[i, j]
            aa_v = None
            if mat_aa is not None:
                ii = next((k for k, x in enumerate(aa_names)
                           if x == nt_names[i]), None)
                jj = next((k for k, x in enumerate(aa_names)
                           if x == nt_names[j]), None)
                if ii is not None and jj is not None:
                    a, b = (ii, jj) if ii < jj else (jj, ii)
                    aa_v = mat_aa[a, b]
            pairs.append({'a': nt_names[i], 'b': nt_names[j],
                          'nt': None if np.isnan(nt_v) else round(nt_v, 2),
                          'aa': None if aa_v is None or np.isnan(aa_v)
                          else round(aa_v, 2)})
    pairs.sort(key=lambda p: (p['nt'] if p['nt'] is not None else -1),
               reverse=True)
    with safe_open(os.path.join(run_dir, 'identity_table.csv'), 'wt') as f:
        f.write('Comparison,NT_identity(%),AA_identity(%)\n')
        for p in pairs:
            f.write(f"{p['a']} vs {p['b']},"
                    f"{'' if p['nt'] is None else p['nt']},"
                    f"{'' if p['aa'] is None else p['aa']}\n")

    # 复合热图（NT 上三角 / AA 下三角）+ NT 分布
    colors = _cmap_colors(palette)
    if n_aa >= 2:
        full_aa = np.full((n, n), np.nan)
        idx = {x: k for k, x in enumerate(aa_names)}
        for i in range(n):
            for j in range(n):
                if i == j:
                    full_aa[i, j] = 100.0
                else:
                    ki, kj = idx.get(nt_names[i]), idx.get(nt_names[j])
                    if ki is not None and kj is not None:
                        full_aa[i, j] = (mat_aa[ki, kj] if ki < kj
                                         else mat_aa[kj, ki])
        order = cluster_order(full_nt, n)
        names_o = [nt_names[i] for i in order]
        plot_heatmap_composite(full_nt[np.ix_(order, order)],
                               full_aa[np.ix_(order, order)],
                               names_o,
                               f'Identity Matrix: {n} sequences\n'
                               'Upper: NT / Lower: AA (%)',
                               os.path.join(run_dir, 'identity_composite'),
                               colors=colors)
    plot_heatmap(full_nt, nt_names,
                 f'NT Identity: {n} sequences (%)',
                 os.path.join(run_dir, 'nt_heatmap'),
                 triangle_only=True, colors=colors)

    data = {'names': nt_names, 'n': n, 'n_aa': n_aa, 'aligned': aligned,
            'nt': [[None if np.isnan(full_nt[i, j])
                    else round(full_nt[i, j], 2) for j in range(n)]
                   for i in range(n)],
            'aa': (None if mat_aa is None else
                   [[None if np.isnan(full_aa[i, j])
                     else round(full_aa[i, j], 2) for j in range(n)]
                    for i in range(n)]),
            'pairs': pairs, 'palette': palette}
    with safe_open(os.path.join(run_dir, 'identity.json'), 'wt') as f:
        _json.dump(data, f, ensure_ascii=False)

    if progress:
        progress(1.0, '完成')
    return {'n': n, 'n_aa': n_aa, 'pairs': n * (n - 1) // 2,
            'table': 'identity_table.csv', 'nt_csv': 'nt_matrix.csv',
            'aa_csv': 'aa_matrix.csv' if mat_aa is not None else None,
            'composite': 'identity_composite.png' if n_aa >= 2 else None,
            'composite_pdf': ('identity_composite.pdf'
                              if n_aa >= 2 else None),
            'nt_heatmap': 'nt_heatmap.png', 'json': 'identity.json'}

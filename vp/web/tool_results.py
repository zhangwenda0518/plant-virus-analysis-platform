# -*- coding: utf-8 -*-
"""工具箱结果读取 API（自 app.py 拆出）。

kvsuite / consensus / verify / 结构比较 / MSA / 建树 / 比对编辑 /
报告数据 / 宿主预测 / contig 深度分析。"""
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid

from flask import (Blueprint, Response, abort, jsonify, render_template, request,
                   send_file, send_from_directory)

from vp.config import DIRS, PLATFORM_ROOT, db_path, engine_cmd
from vp.utils import (TaskLogger, check_path, fmt_size, run_cmd, safe_open,
                      safe_remove)
from vp.web.common import _safe_sample
from vp.web.state import cfg, tool_runs_root as _tool_runs_root
from vp.web.tasks import tm

_log = logging.getLogger('vp.web.tool_results')
from vp.web.tool_jobs import _find_contig_seq

bp = Blueprint('tool_results', __name__)


def _analysis_out(run_dir, contig, action):
    """分析结果 JSON 的固定路径：run/analysis/<安全化contig>_<action>.json。

    run/contig/action 先做严格字符校验 + 规范化，杜绝 ../ 穿越。
    """
    if action not in ('blastn', 'blastx', 'cdd', 'primer'):
        abort(400, '无效的分析类型')
    if not re.fullmatch(r'[A-Za-z0-9_\-]{1,80}', os.path.basename(run_dir)):
        abort(400, '无效的运行名')
    adir = check_path(os.path.join(run_dir, 'analysis'), must_exist=False,
                      in_platform=True)
    os.makedirs(adir, exist_ok=True)
    safe_c = re.sub(r'[^A-Za-z0-9_\-.]', '_', contig)[:60]
    p = os.path.abspath(os.path.join(adir, f'{safe_c}_{action}.json'))
    if not p.startswith(os.path.abspath(adir) + os.sep):
        abort(400, '非法结果路径')
    return check_path(p, must_exist=False, in_platform=True)


@bp.route('/api/tool/viral_contigs')
def api_tool_viral_contigs():
    """某次 contigs 运行的病毒 contig 分类表（metabuli 风格列）。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    tsv = check_path(os.path.join(_tool_runs_root(), run,
                                  'virus_classification.tsv'),
                     must_exist=False, in_platform=True)
    rows = []
    if os.path.isfile(tsv):
        import csv
        with safe_open(tsv) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                r['near_complete'] = r.get('near_complete') == 'True'
                rows.append(r)
    # 并入宿主预测列（若该运行已做宿主预测）
    hp = check_path(os.path.join(_tool_runs_root(), run,
                                 '08_host_analysis', 'host_prediction.tsv'),
                    must_exist=False, in_platform=True)
    host_map = {}                      # 先初始化：未做宿主预测的运行无该文件
    if os.path.isfile(hp):
        import csv
        host_map = {}
        with safe_open(hp) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                c = (r.get('contig') or '').strip()
                h = (r.get('final_host') or '').strip()
                conf = (r.get('confidence_level') or '').strip()
                host_map[c] = h + (f' ({conf})' if conf else '')
        for r in rows:
            r['host'] = host_map.get((r.get('contig') or '').strip(), '')
    return jsonify(rows)


@bp.route('/api/tool/verify_result')
def api_tool_verify_result():
    """候选序列验证结果：calls.tsv 行 + summary.json 摘要。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    vdir = check_path(os.path.join(_tool_runs_root(), run, 'verify'),
                      must_exist=False, in_platform=True)
    rows, summary = [], {}
    calls_tsv = os.path.join(vdir, 'calls.tsv')
    if os.path.isfile(calls_tsv):
        import csv
        with safe_open(calls_tsv) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                r['len'] = int(r.get('len') or 0)
                rows.append(r)
    summ_json = os.path.join(vdir, 'summary.json')
    if os.path.isfile(summ_json):
        with safe_open(summ_json) as f:
            summary = json.load(f)
    return jsonify({'rows': rows, 'summary': summary})


@bp.route('/api/tool/kvsuite_refs')
def api_tool_kvsuite_refs():
    """列出某次已知病毒运行中可用作回贴参考的病毒基因组。

    来源（按优先级）：
      <run>/kvsuite/virus-fasta/ref_<acc>/ref_<acc>.ref.fasta   过滤后的病毒基因组
      <run>/kvsuite/virus-annotations/<acc>.gb / .gtf         注释
    供 t-consensus 选参考、t-variant 查注释用。
    """
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    import csv as _csv
    base = check_path(os.path.join(_tool_runs_root(), run, 'kvsuite'),
                      must_exist=False, in_platform=True)

    def _read_ids(path):
        ids = []
        if not os.path.isfile(path):
            return ids
        with safe_open(path) as f:
            for line in f:
                if line.startswith('>'):
                    ids.append(line[1:].strip().split()[0])
        return ids

    fa_dir = os.path.join(base, 'virus-fasta')
    ann_dir = os.path.join(base, 'virus-annotations')
    refs = []
    if os.path.isdir(fa_dir):
        for dn in sorted(os.listdir(fa_dir)):
            if not dn.startswith('ref_'):
                continue
            sub = os.path.join(fa_dir, dn)
            if not os.path.isdir(sub):
                continue
            fasta = None
            for fn in sorted(os.listdir(sub)):
                if fn.endswith(('.fasta', '.fa', '.fna')):
                    fasta = os.path.join(sub, fn)
                    break
            if not fasta:
                continue
            acc = dn[len('ref_'):]
            try:
                size = os.path.getsize(fasta)
            except OSError:
                size = 0
            gb = os.path.join(ann_dir, f'{acc}.gb')
            gtf = os.path.join(ann_dir, f'{acc}.gtf')
            gtf_n = 0
            if os.path.isfile(gtf):
                with safe_open(gtf) as f:
                    gtf_n = sum(1 for line in f if line.strip())
            refs.append({
                'accession': acc,
                'dir': os.path.relpath(sub, base).replace('\\', '/'),
                'fasta': os.path.relpath(fasta, base).replace('\\', '/'),
                'size': fmt_size(size),
                'seq_ids': _read_ids(fasta),
                'gb': os.path.isfile(gb),
                'gtf': os.path.isfile(gtf),
                'n_cds': gtf_n,
            })

    # BAM：供 t-consensus 复用映射（映射模式=复用）
    bams = []
    bam_dir = os.path.join(base, 'bam')
    if os.path.isdir(bam_dir):
        for fn in sorted(os.listdir(bam_dir)):
            if not fn.endswith('.bam'):
                continue
            p = os.path.join(bam_dir, fn)
            try:
                nbytes = os.path.getsize(p)
            except OSError:
                continue
            bams.append({
                'name': fn,
                'path': p,
                'size': fmt_size(nbytes),
                'size_mb': round(nbytes / 1048576.0, 1),
                'indexed': os.path.isfile(p + '.bai'),
            })

    # 过滤表：给前端展示物种名（末次鉴定结果）
    filt = {}
    ft = os.path.join(base, 'filter', 'filtered.tsv')
    if os.path.isfile(ft):
        with safe_open(ft) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                acc = (r.get('Accession') or '').strip()
                if acc and acc not in filt:
                    filt[acc] = {
                        'species': (r.get('Species') or '').strip(),
                        'coverage': (r.get('Coverage(%)') or '').strip(),
                        'depth': (r.get('MeanDepth') or '').strip(),
                    }
    for r in refs:
        r['species'] = filt.get(r['accession'], {}).get('species', '')
        r['coverage'] = filt.get(r['accession'], {}).get('coverage', '')
        r['depth'] = filt.get(r['accession'], {}).get('depth', '')

    return jsonify({
        'run': run,
        'n_refs': len(refs),
        'refs': refs,
        'bams': bams,
        'has_annotations': os.path.isdir(ann_dir),
    })


@bp.route('/api/tool/kvsuite_result')
def api_tool_kvsuite_result():
    """known_virus_suite 结果：鉴定 / 过滤 / 共识 QC / 变异注释 四张表 + 图文件。

    产物目录 <run>/kvsuite/{identify,filter,consensus,variant,plots}
    """
    import csv as _csv
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    base = check_path(os.path.join(_tool_runs_root(), run, 'kvsuite'),
                      must_exist=False, in_platform=True)

    def _read_tsv(path, limit=None):
        rows = []
        if not os.path.isfile(path):
            return rows
        with safe_open(path) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                rows.append(r)
                if limit and len(rows) >= limit:
                    break
        return rows

    ns = {'identify': [], 'filtered': [], 'discarded': [],
          'consensus_qc': [], 'variants': []}
    ns['identify'] = _read_tsv(os.path.join(base, 'identify',
                                            'all_viruses.summary.tsv'), 500)
    ns['filtered'] = _read_tsv(os.path.join(base, 'filter', 'filtered.tsv'), 500)
    ns['discarded'] = _read_tsv(os.path.join(base, 'filter', 'discarded.tsv'), 200)
    ns['consensus_qc'] = _read_tsv(
        os.path.join(base, 'consensus', 'consensus_qc.tsv'), 500)

    # 变异注释：把所有 per-genome ann.tsv 合起来
    # 注意：kv_variant 把产物直接写在 <out> 下（<out>/annotated/、
    # <out>/variant_summary.json），没有 variant/ 中间层。
    ann_dir = os.path.join(base, 'annotated')
    # QUAL/DP/AF 不在 ann.tsv 里，需从同名的 vcf 补上（键 = CHROM|POS|ALT）
    vcf_map = _kv_load_vcf_metrics(os.path.join(base, 'vcf'))
    if os.path.isdir(ann_dir):
        for fn in sorted(os.listdir(ann_dir)):
            if not fn.endswith('.ann.tsv'):
                continue
            gid = fn[: -len('.ann.tsv')]
            for r in _read_tsv(os.path.join(ann_dir, fn), 300):
                r['genome'] = gid
                k = '{}|{}|{}'.format(r.get('CHROM', ''), r.get('POS', ''),
                                     r.get('ALT', ''))
                m = vcf_map.get(k)
                if m:
                    r['QUAL'] = m['QUAL']
                    r['DP'] = m['DP']
                    r['AF'] = m['AF']
                ns['variants'].append(r)

    # 深度图（<out>/plots/<sample>/*.png 有子目录，必须递归）
    plots = []
    pdir = os.path.join(base, 'plots')
    if os.path.isdir(pdir):
        for root, _dirs, files in os.walk(pdir):
            for fn in sorted(files):
                if fn.lower().endswith(('.png', '.pdf')):
                    rel = os.path.relpath(os.path.join(root, fn), pdir)
                    plots.append(rel.replace(os.sep, '/'))
        plots.sort()

    # 变异图（<out>/variant_plots/，平铺目录）
    variant_plots = []
    vpdir = os.path.join(base, 'variant_plots')
    if os.path.isdir(vpdir):
        for fn in sorted(os.listdir(vpdir)):
            if fn.lower().endswith(('.png', '.pdf')):
                variant_plots.append(fn)
    # 只展示 png（页面内联），pdf 供下载
    variant_plots_png = [p for p in variant_plots if p.lower().endswith('.png')]

    # 扩展变异分析（分子谱 / 变异密度 / 群体遗传），与 variant_plots 同目录平铺
    evo_plots_png = [p for p in variant_plots
                     if p.lower().endswith('.png') and '_evo_' in p]
    evo_manifest = {}
    _eman = os.path.join(base, 'variant_evo', 'evo_manifest.json')
    if os.path.isfile(_eman):
        try:
            with safe_open(_eman) as f:
                evo_manifest = json.load(f)
        except (ValueError, OSError):
            evo_manifest = {}

    # 群体遗传滑窗表 + SNPGenie 汇总（只取头几行供页面预览）
    popgen_tables = {}
    _edir = os.path.join(base, 'variant_evo')
    if os.path.isdir(_edir):
        for fn in sorted(os.listdir(_edir)):
            if fn.endswith('_popgen_window.tsv'):
                acc = fn[: -len('_popgen_window.tsv')]
                popgen_tables[acc] = _read_tsv(os.path.join(_edir, fn), 200)

    snpgenie_products = {}
    if os.path.isdir(_edir):
        for fn in sorted(os.listdir(_edir)):
            if fn.endswith('.snpgenie'):
                acc = fn[: -len('.snpgenie')]
                # 优先读带 ORF 变异计数的注释版（区分无变异/无多态）
                for pname in ('product_results.annotated.txt',
                              'product_results.txt'):
                    prod = os.path.join(_edir, fn, pname)
                    if os.path.isfile(prod):
                        snpgenie_products[acc] = _read_tsv(prod, 200)
                        break

    summary = {}
    vsum = os.path.join(base, 'variant_summary.json')
    if os.path.isfile(vsum):
        try:
            with safe_open(vsum) as f:
                summary['variant'] = json.load(f)
            # 脱敏：results 里的 vcf/ann_vcf/ann_tsv 是绝对路径，
            # 前端渲染只用 accession/genome/n_variants/n_coding_ann，
            # 不把本机路径暴露到页面。产物文件本身保持完整（含路径）。
            var = summary.get('variant') or {}
            for r in (var.get('results') or []):
                for k in ('vcf', 'ann_vcf', 'ann_tsv'):
                    r.pop(k, None)
        except (ValueError, OSError):
            pass
    summary['n_identify'] = len(ns['identify'])
    summary['n_filtered'] = len(ns['filtered'])
    summary['n_discarded'] = len(ns['discarded'])
    summary['n_consensus'] = len(ns['consensus_qc'])
    summary['n_variants'] = len(ns['variants'])
    summary['n_variant_plots'] = len(variant_plots_png)
    summary['n_evo_plots'] = len(evo_plots_png)
    summary['n_snpgenie'] = len(evo_manifest.get('snpgenie') or [])

    return jsonify({'run': run, 'summary': summary, 'plots': plots,
                    'variant_plots': variant_plots_png,
                    'evo_plots': evo_plots_png,
                    'evo_manifest': evo_manifest,
                    'popgen_tables': popgen_tables,
                    'snpgenie_products': snpgenie_products,
                    'identify': ns['identify'], 'filtered': ns['filtered'],
                    'discarded': ns['discarded'],
                    'consensus_qc': ns['consensus_qc'],
                    'variants': ns['variants']})


def _kv_load_vcf_metrics(vcf_dir):
    """扫 <out>/vcf/*.vcf 取 QUAL/DP/AF，键为 CHROM|POS|ALT。

    AF 口径与 kv_variant 过滤一致：AD[1]/INFO/DP（无 AD 时退 AC/AN）。
    单行格式异常只跳过该行并 warning，不丢整个 VCF 剩余记录。
    """
    out = {}
    if not os.path.isdir(vcf_dir):
        return out
    for fn in sorted(os.listdir(vcf_dir)):
        if not fn.lower().endswith('.vcf') or fn.startswith('all.'):
            continue
        n_bad = 0
        try:
            with safe_open(os.path.join(vcf_dir, fn)) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    try:
                        f2 = line.rstrip('\n').split('\t')
                        if len(f2) < 8:
                            n_bad += 1
                            continue
                        chrom, pos, _id, _ref, alt, qual, _filt, info = f2[:8]
                        dp = af = None
                        ad = ac = an = None
                        for item in info.split(';'):
                            if item.startswith('DP='):
                                dp = item[3:]
                            elif item.startswith('AD='):
                                ad = item[3:]
                            elif item.startswith('AC='):
                                ac = item[3:]
                            elif item.startswith('AN='):
                                an = item[3:]
                        dp_i = None
                        try:
                            dp_i = int(dp) if dp is not None else None
                        except (ValueError, TypeError):
                            dp_i = None
                        alt_ad = None
                        if ad:
                            parts = ad.split(',')
                            if len(parts) >= 2:
                                try:
                                    alt_ad = int(parts[1])
                                except (ValueError, TypeError):
                                    alt_ad = None
                        if alt_ad is not None and dp_i:
                            af = round(alt_ad / dp_i, 4)
                        elif ac and an:
                            try:
                                af = round(float(ac) / float(an), 4)
                            except (ValueError, TypeError, ZeroDivisionError):
                                af = None
                        out['{}|{}|{}'.format(chrom, pos, alt)] = {
                            'QUAL': qual, 'DP': dp_i, 'AF': af}
                    except Exception as e:
                        n_bad += 1
                        if n_bad <= 3:
                            _log.warning(
                                f'VCF 行解析失败 {fn} L{line[:60]!r}: {e}')
        except OSError as e:
            _log.warning(f'VCF 读取失败 {fn}: {e}')
            continue
        if n_bad:
            _log.warning(f'VCF {fn}: {n_bad} 行解析失败已跳过')
    return out


@bp.route('/api/tool/consensus_result')
def api_tool_consensus_result():
    """共识与变异结果：coverage.tsv + variants.tsv + summary.json + 分箱覆盖曲线。

    cov 为按 bin 平均的深度序列（最多 ~800 点），避免把上万个逐位深度
    全量塞给前端。
    """
    import csv as _csv
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    cdir = check_path(os.path.join(_tool_runs_root(), run, 'consensus'),
                      must_exist=False, in_platform=True)
    rows, variants, summary = [], [], {}

    cov_tsv = os.path.join(cdir, 'coverage.tsv')
    if os.path.isfile(cov_tsv):
        with safe_open(cov_tsv) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                for k in ('length', 'mapped_reads', 'covered_bp',
                          'n_variants', 'n_snv', 'n_isnv'):
                    try:
                        r[k] = int(r.get(k) or 0)
                    except ValueError:
                        r[k] = 0
                for k in ('median_depth', 'min_depth', 'max_depth'):
                    try:
                        r[k] = int(float(r.get(k) or 0))
                    except ValueError:
                        r[k] = 0
                for k in ('mean_depth', 'coverage_pct',
                          'depth_ok_pct', 'called_pct'):
                    try:
                        r[k] = float(r.get(k) or 0)
                    except ValueError:
                        r[k] = 0.0
                rows.append(r)

    var_tsv = os.path.join(cdir, 'variants.tsv')
    if os.path.isfile(var_tsv):
        with safe_open(var_tsv) as f:
            for i, r in enumerate(_csv.DictReader(f, delimiter='\t')):
                if i >= 500:
                    break
                variants.append(r)

    summ_json = os.path.join(cdir, 'summary.json')
    if os.path.isfile(summ_json):
        with safe_open(summ_json) as f:
            summary = json.load(f)

    cov = {}
    for r in rows:
        name = r.get('contig') or ''
        tag = re.sub(r'[^\w\-.]+', '_', name)[:120]
        cov[name] = _bin_coverage(os.path.join(cdir, tag + '.poscounts.tsv'),
                                  int(r.get('length') or 0))
    return jsonify({'rows': rows, 'variants': variants,
                    'summary': summary, 'cov': cov})


def _bin_coverage(pc_path, length, max_bins=800):
    """poscounts.tsv → 分箱平均深度序列。"""
    if not os.path.isfile(pc_path) or length <= 0:
        return []
    bin_sz = max(1, (length + max_bins - 1) // max_bins)
    nbin = (length + bin_sz - 1) // bin_sz
    sums = [0] * nbin
    cnts = [0] * nbin
    with safe_open(pc_path) as f:
        next(f, None)
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 8:
                continue
            try:
                pos = int(p[0])
                depth = int(p[1]) + int(p[2]) + int(p[3]) + int(p[4]) + int(p[5])
            except ValueError:
                continue
            b = pos // bin_sz
            if 0 <= b < nbin:
                sums[b] += depth
                cnts[b] += 1
    return [round(sums[i] / cnts[i], 1) if cnts[i] else 0 for i in range(nbin)]


@bp.route('/api/tool/structcmp_data')
def api_tool_structcmp_data():
    """结构比较结果（JSON）：identity 矩阵 + 序列名，供页内 Plotly 热图渲染。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    js = check_path(os.path.join(_tool_runs_root(), run,
                                 'identity_matrix.json'),
                    must_exist=True, in_platform=True)
    return send_file(js, mimetype='application/json')


@bp.route('/api/tool/msa_data')
def api_tool_msa_data():
    """structcmp 运行 → SNP-only MSA 展示数据（比较基因组 MSA 查看卡用）。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    aln = check_path(os.path.join(_tool_runs_root(), run, 'aln.fasta'),
                     must_exist=True, in_platform=True)
    from vp.msa_view import snp_view_data
    try:
        return jsonify(snp_view_data(aln))
    except ValueError as e:
        abort(400, str(e))


@bp.route('/api/tool/quicktree_data')
def api_tool_quicktree_data():
    """quicktree 运行 → 树 Newick（进化树查看器卡片直接渲染）。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    run_dir = check_path(os.path.join(_tool_runs_root(), run),
                         must_exist=True, in_platform=True)
    tree_file = None
    for name in ('nj.nwk', 'tree.nwk'):
        p = os.path.join(run_dir, name)
        if os.path.isfile(p):
            tree_file = p
            break
    if not tree_file:
        abort(404, '该运行没有树文件（尚未完成或建树失败，看运行日志）')
    with safe_open(tree_file) as f:
        nwk = f.read().strip()
    if not nwk:
        abort(404, '树文件为空')
    base = os.path.basename(tree_file)
    return jsonify({'newick': nwk, 'file': base,
                    'tool': 'NJ' if base == 'nj.nwk' else 'FastTree'})


_ALIGN_EXTS = ('.fasta', '.fa', '.fna', '.fas', '.aln', '.txt')


def _parse_alignment_fasta(path):
    """读比对/序列 FASTA → dict（names/seqs/长度一致性/类型）。"""
    from vp.utils import iter_fasta
    from vp.sdt_exact import detect_seqtype
    names, seqs = [], []
    for h, s in iter_fasta(path):
        names.append(re.split(r'[\s|]', (h or '').strip())[0][:60] or
                     f'seq{len(names) + 1}')
        seqs.append(s.upper())
    if not seqs:
        raise ValueError('文件中没有序列')
    lens = {len(s) for s in seqs}
    return {'names': names, 'seqs': seqs, 'n': len(seqs),
            'cols': max(len(s) for s in seqs),
            'aligned': len(lens) == 1,
            'type': detect_seqtype(seqs)}


@bp.route('/api/align/file')
def api_align_file():
    """比对查看器数据：path= 平台内或本机（只读）FASTA 比对文件。"""
    path = (request.args.get('path') or '').strip()
    if not path.lower().endswith(_ALIGN_EXTS):
        abort(400, '请提供 FASTA/比对文件（' + '/'.join(_ALIGN_EXTS[:4]) + '…）')
    try:
        p = check_path(path if os.path.isabs(path)
                       else os.path.join(PLATFORM_ROOT, path),
                       must_exist=True)
    except (ValueError, FileNotFoundError) as e:
        abort(400, f'路径不合法: {e}')
    try:
        d = _parse_alignment_fasta(p)
    except ValueError as e:
        abort(400, str(e))
    d['path'] = os.path.relpath(p, PLATFORM_ROOT) \
        if str(p).startswith(str(PLATFORM_ROOT)) else str(p)
    return jsonify(d)


@bp.route('/api/align/save', methods=['POST'])
def api_align_save():
    """编辑保存：content 写到源文件旁 <stem>.edit.fasta（平台内）；
    源在平台外时落 tool_runs/align_edit_<ts>/edited.fa。返回新路径。"""
    body = request.get_json(force=True) or {}
    path = (body.get('path') or '').strip()
    content = str(body.get('content') or '')
    if not path or not content.strip():
        abort(400, '参数不完整（path 与 content 必填）')
    try:
        p = check_path(path if os.path.isabs(path)
                       else os.path.join(PLATFORM_ROOT, path),
                       must_exist=True)
    except (ValueError, FileNotFoundError) as e:
        abort(400, f'路径不合法: {e}')
    # 校验能按 FASTA 解析（至少 2 条、只含合法字符集宽容量）
    try:
        probe = _parse_alignment_fasta_string(content)
    except ValueError as e:
        abort(400, f'内容不是合法 FASTA: {e}')
    in_plat = str(p).startswith(str(PLATFORM_ROOT) + os.sep)
    if in_plat:
        stem, _ext = os.path.splitext(p)
        dst = stem + '.edit.fasta'
    else:
        dst = check_path(os.path.join(
            _tool_runs_root(), f"align_edit_{time.strftime('%Y%m%d_%H%M%S')}",
            'edited.fasta'), must_exist=False, in_platform=True)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
    with safe_open(dst, 'wt') as f:
        f.write(content if content.endswith('\n') else content + '\n')
    return jsonify({'saved': dst, 'aligned': probe['aligned'],
                    'n': probe['n'], 'cols': probe['cols']})


def _parse_alignment_fasta_string(content):
    """与 _parse_alignment_fasta 同口径，但直接吃文本。"""
    from vp.sdt_exact import detect_seqtype
    names, seqs, cur = [], [], None
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if cur is not None:
                seqs.append(''.join(cur))
            names.append(re.split(r'[\s|]', line[1:].strip())[0][:60] or
                         f'seq{len(names) + 1}')
            cur = []
        elif cur is not None:
            cur.append(line)
        else:
            raise ValueError('第一条记录前出现序列行')
    if cur is not None:
        seqs.append(''.join(cur))
    if len(seqs) < 1 or not names:
        raise ValueError('未解析到序列')
    joined = [s.upper() for s in seqs]
    if any(not s for s in joined):
        raise ValueError('存在空序列')
    return {'names': names, 'seqs': joined, 'n': len(joined),
            'cols': max(len(s) for s in joined),
            'aligned': len({len(s) for s in joined}) == 1,
            'type': detect_seqtype(joined)}


_TREE_EXTS = ('.nwk', '.newick', '.treefile', '.contree', '.tre', '.tree')


@bp.route('/api/tree/file')
def api_tree_file():
    """树文件路径 → Newick 文本（进化树查看卡「本机树文件」输入用）。

    路径可以是平台内相对路径或本机绝对路径（只读），扩展名白名单限制。
    """
    path = (request.args.get('path') or '').strip()
    if not path.lower().endswith(_TREE_EXTS):
        abort(400, '请提供树文件（'
                    + ' / '.join(_TREE_EXTS) + '）')
    try:
        p = check_path(path if os.path.isabs(path)
                       else os.path.join(PLATFORM_ROOT, path),
                       must_exist=True)
    except (ValueError, FileNotFoundError) as e:
        abort(400, f'路径不合法: {e}')
    with safe_open(p) as f:
        nwk = f.read().strip()
    if not nwk:
        abort(404, '树文件为空')
    return jsonify({'newick': nwk, 'file': os.path.basename(p)})


@bp.route('/api/tool/report')
def api_tool_report():
    """生成（带缓存）并打开工具运行的结果报告（桑基/旭日/分类表）。"""
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    run_dir = check_path(os.path.join(_tool_runs_root(), run),
                         must_exist=True, in_platform=True)
    from vp.viz import build_contig_report
    try:
        build_contig_report(run_dir)
    except RuntimeError as e:
        abort(400, str(e))
    from flask import redirect
    return redirect(f'/tool_runs/{run}/report.html')


@bp.route('/api/tool/report_data')
def api_tool_report_data():
    """分类报告数据（JSON）：供工具页内联渲染桑基/旭日/分类表/明细表。

    contigs_* 运行 → mode=contig（含 contig 明细 + 宿主列）；
    identify_* 运行 → mode=reads（含 classified sequences 表）。
    """
    run = request.args.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    run_dir = check_path(os.path.join(_tool_runs_root(), run),
                         must_exist=True, in_platform=True)
    import glob as _glob
    krs = (_glob.glob(os.path.join(run_dir, '*.kreport2')) +
           _glob.glob(os.path.join(run_dir, '**', '*.kreport2'),
                      recursive=True))
    if not krs:
        abort(400, '运行目录中没有 kreport 分类报告')
    from vp.kunpeng import parse_kreport
    krows = parse_kreport(sorted(krs)[0])
    report_rows = [{'proportion': r['percent'], 'count': int(r['frags']),
                    'rank': r['rank'], 'taxid': r['taxid'],
                    'name': r['name'], 'depth': r['depth']}
                   for r in krows]

    import csv as _csv
    tsv = os.path.join(run_dir, 'virus_classification.tsv')
    ids_tsv = os.path.join(run_dir, 'viral_ids.tsv')
    mode, vc_rows, ids_rows = ('contig', [], [])
    if os.path.isfile(tsv):
        mode = 'contig'
        with safe_open(tsv) as f:
            vc_rows = list(_csv.DictReader(f, delimiter='\t'))
    elif os.path.isfile(ids_tsv):
        mode = 'reads'
        with safe_open(ids_tsv) as f:
            ids_rows = list(_csv.DictReader(f, delimiter='\t'))

    host_map, host_cats = {}, {}
    hp = os.path.join(run_dir, '08_host_analysis', 'host_prediction.tsv')
    if os.path.isfile(hp):
        with safe_open(hp) as f:
            for r in _csv.DictReader(f, delimiter='\t'):
                c = (r.get('contig') or '').strip()
                h = (r.get('final_host') or '').strip()
                conf = (r.get('confidence_level') or '').strip()
                if c:
                    host_map[c] = h + (f' ({conf})' if conf else '')
                if h:
                    host_cats[h] = host_cats.get(h, 0) + 1
    for r in vc_rows:
        r['host'] = host_map.get((r.get('contig') or '').strip(), '')

    n_total, n_virus = 0, 0
    if mode == 'contig':
        cfa = os.path.join(run_dir, 'contigs.filtered.fasta')
        if os.path.isfile(cfa):
            from vp.utils import count_fasta_seqs
            n_total = count_fasta_seqs(cfa)
        n_virus = len(vc_rows)
    else:
        n_total = n_virus = len(ids_rows)

    host_list = [{'cat': k, 'n': v} for k, v in
                 sorted(host_cats.items(), key=lambda x: -x[1])]

    # taxburst 旭日图（Krona 等价物，离线交互 HTML）
    from vp.viz import build_taxburst
    tb_ok, tb_unc = build_taxburst(run_dir, krows)

    dl = []
    if mode == 'contig':
        dl = [('viral_contigs.fasta?dl=1', 'Virus Sequences (FASTA)'),
              ('virus_classification.tsv?dl=1', 'Virus Classification (TSV)'),
              ('contigs.filtered.fasta?dl=1', 'All Contigs (FASTA)'),
              (os.path.basename(sorted(krs)[0]) + '?dl=1', 'Report (kreport TSV)')]
        if os.path.isfile(hp):
            dl.append(('08_host_analysis/host_prediction.tsv?dl=1',
                       'Host Prediction (TSV)'))
    else:
        import glob as _g2
        for fa in sorted(_g2.glob(os.path.join(run_dir, 'viral_sequences*.fasta'))):
            dl.append((os.path.basename(fa) + '?dl=1', 'Viral Sequences (FASTA)'))
        if os.path.isfile(ids_tsv):
            dl.append(('viral_ids.tsv?dl=1', 'Classified IDs (TSV)'))
        dl.append((os.path.basename(sorted(krs)[0]) + '?dl=1',
                   'Report (kreport TSV)'))
    downloads = [{'href': href, 'label': label} for href, label in dl
                 if os.path.isfile(os.path.join(run_dir, href.split('?')[0]))]

    return jsonify({'run': run, 'mode': mode,
                    'report': report_rows, 'vc': vc_rows, 'ids': ids_rows,
                    'host': host_list, 'host_count': len(host_map),
                    'n_total': n_total, 'n_virus': n_virus,
                    'downloads': downloads,
                    'taxburst': tb_ok})


@bp.route('/api/tool/hostpredict_run', methods=['POST'])
def api_tool_hostpredict_run():
    """对既有 contigs 运行目录原地做 ICTV 宿主预测（结果并入该运行）。"""
    body = request.get_json(force=True) or {}
    run = body.get('run') or ''
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        abort(400, '无效的运行名')
    run_dir = check_path(os.path.join(_tool_runs_root(), run),
                         must_exist=True, in_platform=True)
    tsv = os.path.join(run_dir, 'virus_classification.tsv')
    if not os.path.isfile(tsv):
        abort(400, '该运行没有 virus_classification.tsv（仅 contig 分类运行支持）')
    from vp.kunpeng import db_ready
    if not db_ready(cfg.databases['host']):
        abort(400, '宿主库未就绪，请先到「数据库构建」页构建宿主库')

    def job(log, prog, cancel):
        import shutil as _shutil
        logger = TaskLogger(callback=log)
        a_dir = check_path(os.path.join(run_dir, '03_assembly'),
                           must_exist=False, in_platform=True)
        os.makedirs(a_dir, exist_ok=True)
        prog('prep', 0.05, '整理输入')
        norm = os.path.join(a_dir, 'virus_contigs.tsv')
        with safe_open(tsv) as f:
            header = f.readline().rstrip('\r\n').split('\t')
        if 'kunpeng_taxid' in header:
            _shutil.copyfile(tsv, norm)
        else:
            cols = ['contig', 'length', 'kunpeng_flag', 'kunpeng_taxid',
                    'kunpeng_species', 'blast_top_hit', 'blast_identity(%)',
                    'blast_coverage_hsp(%)', 'blast_aln_len', 'blast_species',
                    'blast_family']
            import csv as _csv
            with safe_open(tsv) as f, safe_open(norm, 'wt') as w:
                w.write('\t'.join(cols) + '\n')
                for row in _csv.DictReader(f, delimiter='\t'):
                    kt = str(row.get('taxid') or '').strip()
                    w.write('\t'.join([
                        str(row.get('contig') or '').strip(),
                        str(row.get('length') or '').strip(),
                        'C' if kt.isdigit() and int(kt) > 0 else 'U', kt,
                        str(row.get('taxon') or '').strip(),
                        '', '', '', '', '', '']) + '\n')
        fa = os.path.join(run_dir, 'viral_contigs.fasta')
        if os.path.isfile(fa):
            _shutil.copyfile(fa, os.path.join(a_dir, 'viral_contigs.fasta'))
        from vp.host_analysis import predict_hosts
        prog('predict', 0.15, 'ICTV 宿主概率级联预测')
        res = predict_hosts(run_dir, logger=logger, force=True)
        prog('done', 1.0, f"完成：宿主判定 {res.get('n_contigs', 0)} 条")
        logger.close()
        return res

    ts = time.strftime('%Y%m%d_%H%M%S')
    tid = tm.start(cfg.tr(f'宿主预测·{run} {ts}', f'Host prediction·{run} {ts}'),
                   job, log_file=os.path.join(run_dir, 'hostpredict.log'))
    return jsonify({'task': tid})


@bp.route('/api/tool/analysis')
def api_tool_analysis():
    """读取已完成的 contig 分析结果（JSON 缓存）。"""
    run = request.args.get('run') or ''
    contig = request.args.get('contig') or ''
    action = request.args.get('action') or ''
    seq = (request.args.get('seq') or '').strip()
    if not contig or action not in ('blastn', 'blastx', 'cdd', 'primer'):
        abort(400, '参数不完整')
    # 直接用字符串序列查询（seq 输入的伪 run 缓存）
    if seq:
        import hashlib as _hl
        _key = _hl.md5((contig + '|' + seq).encode('utf-8')).hexdigest()[:12]
        _safe = re.sub(r'[^A-Za-z0-9_\-.]', '_', contig)[:30]
        run_dir = os.path.join(_tool_runs_root(), '_seq_input', _safe + '_' + _key)
        run_dir = check_path(run_dir, must_exist=True, in_platform=True)
        p = _analysis_out(run_dir, _safe + '_' + _key, action)
    else:
        if not run or not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
            abort(400, '无效的运行名')
        run_dir = check_path(os.path.join(_tool_runs_root(), run),
                             must_exist=True, in_platform=True)
        p = _analysis_out(run_dir, contig, action)
    if not os.path.isfile(p):
        abort(404, '该分析尚未运行')
    with safe_open(p) as f:
        return Response(f.read(), mimetype='application/json')


@bp.route('/api/tool/analyze', methods=['POST'])
def api_tool_analyze():
    """单 contig 深度分析四件套：blastn / blastx / cdd / primer。

    body: {run, contig, action}
    NCBI 在线分析走 vp.contig_annot（域名白名单 + IP 校验），
    结果缓存为 run/analysis/<contig>_<action>.json。
    """
    body = request.get_json(force=True) or {}
    run = body.get('run') or ''
    contig = body.get('contig') or ''
    action = body.get('action') or ''
    seq = (body.get('seq') or '').strip()
    in_seq = seq  # 供闭包安全读取（避免 UnboundLocalError）
    if action not in ('blastn', 'blastx', 'cdd', 'primer'):
        abort(400, '无效的分析类型')
    if not run or not contig:
        abort(400, '参数不完整')
    # 直接输入序列模式（无需 run），或从运行目录取 contig 序列
    if seq:
        # 用共享伪运行缓存，避免污染具体工具运行目录
        run_dir = check_path(os.path.join(_tool_runs_root(), '_seq_input'),
                             must_exist=False, in_platform=True)
        import hashlib as _hl
        _key = _hl.md5((contig + '|' + seq).encode('utf-8')).hexdigest()[:12]
        contig = re.sub(r'[^A-Za-z0-9_\-.]', '_', contig)[:30]
        _used_for_cache = contig + '_' + _key
        run_dir = os.path.join(run_dir, _used_for_cache)
        os.makedirs(run_dir, exist_ok=True)
        out_path = _analysis_out(run_dir, _used_for_cache, action)
    else:
        if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
            abort(400, '无效的运行名')
        run_dir = check_path(os.path.join(_tool_runs_root(), run),
                             must_exist=True, in_platform=True)
        out_path = _analysis_out(run_dir, contig, action)
    if os.path.isfile(out_path):
        return jsonify({'task': None, 'cached': True})

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from vp import contig_annot as ca
        prog('prep', 0.1, '查找 contig 序列')
        if in_seq:
            # 清洗直接输入的序列：去掉 FASTA header 与换行，仅留核酸字母
            seq = ''.join(line.strip().upper() for line in in_seq.splitlines()
                          if line.strip() and not line.strip().startswith('>'))
            if not seq:
                raise RuntimeError('输入序列为空')
        else:
            seq = _find_contig_seq(run_dir, contig)
            if not seq:
                raise RuntimeError(f'运行目录中未找到 contig: {contig}')

        if action == 'primer':
            from vp.primer import design_primers_for_seq
            prog('primer3', 0.5, '本地 primer3 设计引物')
            pairs = design_primers_for_seq(contig, seq, num_return=5)
            result = {'action': 'primer', 'contig': contig,
                      'length': len(seq), 'primers': pairs}

        elif action == 'cdd':
            prog('orf', 0.3, '6-frame 翻译取最长 ORF')
            prot = ca.longest_orf_protein(seq)
            if len(prot) < 50:
                raise RuntimeError('未找到 ≥50aa 的 ORF，无法做 CDD 域搜索')
            prog('cdd', 0.5, '提交 NCBI CDD（数分钟，耐心等待）')
            cdsid = ca.submit_cdd(f'>{contig}\n{prot}\n')
            hits = ca.poll_cdd(cdsid, cancel=cancel)
            result = {'action': 'cdd', 'contig': contig,
                      'orf_len_aa': len(prot), 'hits': hits}

        else:  # blastn / blastx（NCBI URL API，virus-restricted）
            prog('submit', 0.3, f'提交 NCBI {action}（数分钟，耐心等待）')
            rid, _rtoe = ca.submit_blast(
                action, 'nt' if action == 'blastn' else 'nr',
                f'>{contig}\n{seq}\n', entrez_query='viruses[Organism]')
            prog('poll', 0.5, f'NCBI 任务 {rid} 运行中')
            hits = ca.poll_blast(rid, cancel=cancel)
            result = {'action': action, 'contig': contig, 'rid': rid,
                      'hits': hits}

        with safe_open(out_path, 'wt') as f:
            json.dump(result, f, ensure_ascii=False, indent=1)
        prog('done', 1.0, '完成')
        logger.close()
        return {'out': out_path}

    tid = tm.start(cfg.tr(f'Contig分析·{action} {contig[:30]}',
                         f'Contig {action} {contig[:30]}'), job,
                   log_file=os.path.join(run_dir, 'run.log'))
    return jsonify({'task': tid, 'cached': False})

# -*- coding: utf-8 -*-
"""
宿主预测模块（吸收 MMPV-RNA virome_discovery_pipeline 的 ICTV 宿主区分方法）。

方法（C9 级联查表 + 宿主元数据交叉验证）：
  1. contig → BLAST top hit accession → 病毒库 info.tsv 的 ICTV 分类
     （VMR_Species / VMR_Genus / VMR_Family）
  2. ICTV 级联查宿主概率表（Species → Genus → Family → Order，含列错位
     容错，取分类层级最深、置信度最高者）→ Host_ICTV
  3. 交叉证据：info.tsv 的 NCBI Host 元数据（该 accession 的官方宿主记录）
     经 NCBI taxonomy 归类到宿主类别 → Host_Meta
  4. 决策：两者一致=Agree（高置信）；不一致时 Host 元数据是一手证据优先；
     仅其一则用其一；全无=Unknown

输出（08_host_analysis/）:
  host_prediction.tsv   每个病毒 contig 的宿主预测明细
  host_summary.tsv      宿主类别汇总
  {类别}.classified.fasta  按宿主拆分的病毒 contigs
  sankey_host.html      病毒科 → 宿主类别 桑基图
  sunburst_host.html    病毒 ICTV 分类旭日图（科→属→种）
  summary.json          汇总（GUI 卡片与报告用）
"""
import os
import csv
import json

from .config import DIRS, get_config, db_path
from .utils import check_path, safe_open, is_step_done, mark_step_done

# 分类层级深度（级联优先级）
_LEVEL_RANK = {'Species': 4, 'Genus': 3, 'Family': 2, 'Order': 1}
_CONF_RANK = {'High': 3, 'Medium': 2, 'Low (Singleton/Rare)': 1,
              'Low (Shared Genus)': 1, 'Unknown': 0}

# NCBI taxonomy 内的宿主类别锚点 taxid → 类别
_TAXID_CATS = [
    (9606, 'Human'),
    (50557, 'Insecta'),
    (6859, 'Arachnida'),
    (8782, 'Aves'),
    (40674, 'Animal_other'),      # Mammalia
    (33208, 'Animal_other'),      # Metazoa
    (4762, 'Oomycetes'),
    (4751, 'Fungi'),
    (33090, 'Plant'),             # Viridiplantae
    (2, 'Bacteria'),
    (2157, 'Archaea'),
    (33630, 'Protist'),           # Stramenopiles
    (2759, 'Protist'),            # Eukaryota 兜底
]


def _default_info_tsv():
    """病毒参考元数据表：优先 virus_ref 规范化元数据（谱系补全），回退旧表。"""
    try:
        from . import virus_ref
        if virus_ref.available():
            p = virus_ref.build_meta()
            if p and os.path.isfile(p):
                return check_path(p, must_exist=True, in_platform=True)
    except Exception:
        pass
    from .config import PLATFORM_ROOT
    return check_path(os.path.join(PLATFORM_ROOT, 'virus-db',
                                   'final.cluster.ref_info.tsv'),
                      must_exist=True, in_platform=True)


def _default_prob_dir():
    return check_path(db_path('misc', 'prob'),
                      must_exist=True, in_platform=True)


def load_prob_table(path, key_col):
    """加载宿主概率表 → {taxon: {...}}（与参考 C7 输出同构）。"""
    if not os.path.isfile(path):
        return {}
    lookup = {}
    with safe_open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            taxon = (row.get(key_col) or '').strip()
            if not taxon:
                continue
            try:
                conf = float(row.get('Integrated_Confidence', 0) or 0)
            except ValueError:
                conf = 0.0
            lookup[taxon] = {
                'Predicted_Host': row.get('Predicted_Host', 'Unknown'),
                'Confidence_Level': row.get('Confidence_Level', 'Unknown'),
                'Integrated_Confidence': conf,
            }
    return lookup


def cascade_lookup(order, family, genus, species, tables):
    """C9 级联查表：Species→Genus→Family→Order 优先，加列错位容错。

    返回 dict(Predicted_Host, Confidence_Level, Integrated_Confidence,
    Determination_Level)；无命中返回 Unknown。
    """
    ord_l, fam_l, gen_l, sp_l = tables
    entries = [('Order', order or '', ord_l),
               ('Family', family or '', fam_l),
               ('Genus', genus or '', gen_l),
               ('Species', species or '', sp_l)]
    candidates = []
    for level, value, table in entries:
        if not value:
            continue
        if value in table:                      # 本级精确命中
            candidates.append((level, level, table[value]))
        for lv2, _v2, t2 in entries:            # 列错位容错（值实际属于别级）
            if lv2 == level:
                continue
            if value in t2:
                candidates.append((level, lv2 + '*', t2[value]))
    if not candidates:
        return {'Predicted_Host': 'Unknown', 'Confidence_Level': 'Unknown',
                'Integrated_Confidence': 0.0, 'Determination_Level': 'None'}
    candidates.sort(key=lambda x: (_LEVEL_RANK.get(x[1].rstrip('*'), 0),
                                   x[2].get('Integrated_Confidence', 0)),
                    reverse=True)
    col, real_level, info = candidates[0]
    return {'Predicted_Host': info['Predicted_Host'],
            'Confidence_Level': info['Confidence_Level'],
            'Integrated_Confidence': info['Integrated_Confidence'],
            'Determination_Level': f'{real_level}(via {col})'}


def _load_info_map(info_tsv):
    """病毒库 info.tsv → (按 accession 索引, 按 taxid 索引)。

    每条: {species, genus, family, host_meta}（ICTV 分类 + NCBI 宿主元数据）。
    """
    by_acc, by_taxid = {}, {}
    with safe_open(info_tsv) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            acc = (row.get('Accession') or '').strip()
            taxid = (row.get('Taxid') or '').strip()
            item = {
                'species': (row.get('VMR_Species') or row.get('Species_ICTV')
                            or row.get('Species_NCBI') or '').strip(),
                'genus': (row.get('VMR_Genus') or '').strip(),
                'family': (row.get('VMR_Family') or '').strip(),
                'host_meta': (row.get('Host') or '').strip(),
            }
            if acc:
                by_acc[acc] = item
            if taxid.isdigit():
                by_taxid[int(taxid)] = item
    return by_acc, by_taxid


def _taxid_category_chain(taxid, parent):
    """沿 lineage 向上，返回 (类别, 命中 taxid)；无锚点返回 (None, None)。"""
    seen = set()
    t = taxid
    while t and t not in seen:
        seen.add(t)
        for anchor, cat in _TAXID_CATS:
            if t == anchor:
                return cat, t
        t = parent.get(t)
    return None, None


def _meta_category_map(host_names, tax_dir, cache_path):
    """NCBI 宿主名称 → 宿主类别（names.dmp/nodes.dmp 解析 + JSON 缓存）。"""
    wanted = {n.strip().lower() for n in host_names if n and n.strip()}
    cache = {}
    if os.path.isfile(cache_path):
        try:
            with safe_open(cache_path) as f:
                cache = json.load(f)
        except (OSError, ValueError):
            cache = {}
    missing = wanted - set(cache)
    if missing:
        names_dmp = check_path(os.path.join(tax_dir, 'names.dmp'),
                               must_exist=True, in_platform=True)
        # names.dmp 行格式: "taxid | 名称 | 唯一名 | 名称类别 |"
        name2tax = {}
        with safe_open(names_dmp) as f:
            for line in f:
                p = line.split('\t|\t')
                if len(p) < 4:
                    continue
                key = p[1].strip().lower()
                if key in missing and key not in name2tax:
                    try:
                        name2tax[key] = int(p[0])
                    except ValueError:
                        continue
        nodes_dmp = check_path(os.path.join(tax_dir, 'nodes.dmp'),
                               must_exist=True, in_platform=True)
        parent = {}
        with safe_open(nodes_dmp) as f:
            for line in f:
                p = line.split('\t|\t')
                if len(p) < 3:
                    continue
                try:
                    t, par = int(p[0]), int(p[1])
                except ValueError:
                    continue
                parent[t] = par
        for key, t in name2tax.items():
            cat, _ = _taxid_category_chain(t, parent)
            cache[key] = cat or 'Unknown'
        try:
            with safe_open(cache_path, 'wt') as f:
                json.dump(cache, f, ensure_ascii=False)
        except OSError:
            pass
    return cache


def _lineage_map_for(taxids_needed, tax_dir, cache_path):
    """taxid → 完整分类（order/family/genus/species 学名），nodes/names.dmp
    解析 + JSON 缓存（跨样品增量）。这是把 kunpeng C 行 taxid 转为完整
    分类的通用方法：C 行 taxid 可能是 LCA 判定的上级节点（属/科级），
    需沿 nodes.dmp 谱系上溯取各级名称，才能喂给 C9 级联查表。
    """
    cache = {}
    if os.path.isfile(cache_path):
        try:
            with safe_open(cache_path) as f:
                cache = json.load(f)
        except (OSError, ValueError):
            cache = {}
    need = {int(t) for t in taxids_needed
            if str(t).lstrip('-').isdigit() and int(t) not in cache
            and int(t) > 0}
    if need:
        # 1) nodes.dmp: taxid → (parent, rank)，沿谱系上溯定级
        parent, rank = {}, {}
        nodes_dmp = check_path(os.path.join(tax_dir, 'nodes.dmp'),
                               must_exist=True, in_platform=True)
        with safe_open(nodes_dmp) as f:
            for line in f:
                p = line.split('\t|\t')
                if len(p) < 4:
                    continue
                try:
                    t, par = int(p[0]), int(p[1])
                except ValueError:
                    continue
                parent[t] = par
                rank[t] = p[2].strip()
        # 2) 每个 taxid 上溯，记录最先遇到的 species/genus/family/order 节点
        lineage_taxids = {}
        for t in need:
            hits = {}
            cur = t
            seen = set()
            while cur and cur not in seen:
                seen.add(cur)
                rk = rank.get(cur, '')
                if rk in ('species', 'genus', 'family', 'order') \
                        and rk not in hits:
                    hits[rk] = cur
                cur = parent.get(cur, 0)
            lineage_taxids[t] = hits
        # 3) names.dmp: 取这些节点的学名
        name_taxids = {tid for hits in lineage_taxids.values()
                       for tid in hits.values()}
        names_dmp = check_path(os.path.join(tax_dir, 'names.dmp'),
                               must_exist=True, in_platform=True)
        tid2name = {}
        tid2sci = set()
        with safe_open(names_dmp) as f:
            for line in f:
                p = line.split('\t|\t')
                if len(p) < 4:
                    continue
                try:
                    t = int(p[0])
                except ValueError:
                    continue
                if t in name_taxids:
                    nm = p[1].strip()
                    if t not in tid2name:
                        tid2name[t] = nm
                    # 学名优先（同 taxid 可能有多行：学名/等名/同义词）
                    if 'scientific name' in p[3] and t not in tid2sci:
                        tid2name[t] = nm
                        tid2sci.add(t)
        # 4) 组装并写缓存
        for t in need:
            hits = lineage_taxids.get(t, {})
            entry = {rk: tid2name.get(tid, '') for rk, tid in hits.items()}
            cache[t] = {'order': entry.get('order', ''),
                        'family': entry.get('family', ''),
                        'genus': entry.get('genus', ''),
                        'species': entry.get('species', '')}
        try:
            with safe_open(cache_path, 'wt') as f:
                json.dump(cache, f, ensure_ascii=False)
        except OSError:
            pass
    return cache


def predict_hosts(sample_dir, info_tsv=None, prob_dir=None, threads=None,
                  logger=None, force=False):
    """对 ③组装 产出的病毒 contigs 做 ICTV 宿主预测。返回 summary dict。"""
    step = 'host_analysis'
    out_dir = check_path(os.path.join(sample_dir, '08_host_analysis'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    summary_file = os.path.join(out_dir, 'summary.json')
    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("宿主预测已完成，跳过")
        with safe_open(summary_file) as f:
            return json.load(f)

    info_tsv = info_tsv or _default_info_tsv()
    prob_dir = prob_dir or _default_prob_dir()
    ctg_tsv = check_path(os.path.join(sample_dir, '03_assembly',
                                      'virus_contigs.tsv'),
                         must_exist=True, in_platform=True)
    viral_fa = check_path(os.path.join(sample_dir, '03_assembly',
                                       'viral_contigs.fasta'),
                          must_exist=False, in_platform=True)

    if logger:
        logger.log("加载 ICTV 宿主概率表 ...")
    tables = (
        load_prob_table(os.path.join(prob_dir, 'order_host_probability.tsv'),
                        'Order'),
        load_prob_table(os.path.join(prob_dir, 'family_host_probability.tsv'),
                        'Family'),
        load_prob_table(os.path.join(prob_dir, 'genus_host_probability.tsv'),
                        'Genus'),
        load_prob_table(os.path.join(prob_dir, 'species_host_probability.tsv'),
                        'Species'),
    )
    if logger:
        logger.log(f"  Order {len(tables[0])} / Family {len(tables[1])} / "
                   f"Genus {len(tables[2])} / Species {len(tables[3])} 条目")

    info_acc, info_taxid = _load_info_map(info_tsv)
    if logger:
        logger.log(f"病毒库 info 映射: {len(info_acc)} 条 accession / "
                   f"{len(info_taxid)} 条 taxid")

    # 物种级已知宿主 + 逐 accession 的拓扑/分子类型/完整性（virus_ref 元数据）
    acc_meta, sp_hosts = {}, {}
    try:
        from . import virus_ref
        acc_meta = virus_ref.meta_by_accession()
        sp_hosts = virus_ref.species_host_map()
        if logger:
            logger.log(f"virus_ref 元数据: {len(acc_meta)} 条 / "
                       f"物种级已知宿主 {len(sp_hosts)} 个物种")
    except Exception as e:
        if logger:
            logger.log(f"virus_ref 元数据不可用（跳过已知宿主交叉验证）: {e}",
                       "WARN")

    # 宿主元数据类别（taxonomy 归类，带缓存）；known hosts 一并归类
    meta_names = {v['host_meta'] for v in info_acc.values() if v['host_meta']}
    meta_names |= {h.strip() for lst in sp_hosts.values()
                   for h in lst if h.strip()}
    tax_dir = DIRS['taxonomy']
    if meta_names and os.path.isfile(check_path(
            os.path.join(tax_dir, 'names.dmp'), must_exist=False,
            in_platform=True)):
        cache_path = check_path(os.path.join(prob_dir,
                                             'host_meta_category_cache.json'),
                                must_exist=False, in_platform=True)
        meta_cat = _meta_category_map(meta_names, tax_dir, cache_path)
        if logger:
            logger.log(f"NCBI 宿主元数据类别: {len(meta_names)} 个名称已归类")
    else:
        meta_cat = {}

    # 逐 contig 预测。分类信息以 kunpeng 对 contigs 的分类判定（kunpeng_taxid）
    # 为主：taxid 精确命中 info.tsv → VMR 种/属/科；未命中（LCA 判到上级节点）
    # 时经 nodes/names.dmp 谱系解析为完整分类；BLAST top hit 最后回退。
    kt_set = set()
    with safe_open(ctg_tsv) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            kt = (row.get('kunpeng_taxid') or '').strip()
            if (row.get('kunpeng_flag') or '').strip() == 'C' and kt.isdigit():
                kt_set.add(int(kt))
    info_taxid_keys = set(info_taxid.keys())
    lineage = _lineage_map_for(
        (t for t in kt_set if t not in info_taxid_keys), tax_dir,
        check_path(os.path.join(prob_dir, 'taxid_lineage_cache.json'),
                   must_exist=False, in_platform=True))
    if logger:
        n_lin = sum(1 for t in kt_set if t not in info_taxid_keys)
        logger.log(f"kunpeng taxid {len(kt_set)} 个：{n_lin} 个为上级节点"
                   f"（LCA），已谱系解析为完整分类")

    rows = []
    src_stat = {'kunpeng': 0, 'kunpeng_lineage': 0, 'blast': 0, 'none': 0}
    with safe_open(ctg_tsv) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            contig = (row.get('contig') or '').strip()
            if not contig:
                continue
            tax_src = 'none'
            tax = {}
            kt = (row.get('kunpeng_taxid') or '').strip()
            if (row.get('kunpeng_flag') or '').strip() == 'C' and kt.isdigit():
                kti = int(kt)
                tax = info_taxid.get(kti, {})
                if tax:
                    tax_src = 'kunpeng'                    # 精确命中库内 accession
                else:
                    lin = lineage.get(kti, {})
                    if any(lin.values()):
                        tax = {'species': lin.get('species', ''),
                               'genus': lin.get('genus', ''),
                               'family': lin.get('family', ''),
                               'host_meta': ''}
                        tax_src = 'kunpeng_lineage'        # LCA 上级节点 → 谱系解析
            if not tax:
                acc = (row.get('blast_top_hit') or '').strip()
                tax = info_acc.get(acc, {})
                if tax:
                    tax_src = 'blast'
                else:
                    # BLAST 也无命中：用 ③ 里的 kunpeng/blast 物种名做最后回退
                    tax = {'species': (row.get('blast_species') or '').strip(),
                           'genus': '', 'family': (row.get('blast_family') or '').strip(),
                           'host_meta': ''}
            src_stat[tax_src] += 1
            order = ''
            hit = cascade_lookup(order, tax.get('family', ''),
                                 tax.get('genus', ''), tax.get('species', ''),
                                 tables)
            meta_name = tax.get('host_meta', '')
            meta = meta_cat.get(meta_name.strip().lower(), 'Unknown') \
                if meta_name else 'Unknown'
            ictv = hit['Predicted_Host']
            ictv_conf = hit['Confidence_Level']
            # 决策
            if ictv != 'Unknown' and meta not in ('', 'Unknown'):
                if ictv == meta:
                    final, method = ictv, 'ICTV_Meta_Agree'
                else:
                    final, method = meta, 'Meta_Preferred'
            elif ictv != 'Unknown':
                final, method = ictv, 'ICTV_Prob'
            elif meta not in ('', 'Unknown'):
                final, method = meta, 'Meta_Only'
            else:
                final, method = 'Unknown', 'Unassigned'

            # 已知宿主交叉验证：virus_ref 物种级宿主集 + accession 级宿主，
            # 归类后与最终预测比对——预测宿主 ∉ 已知宿主类别 → WARN（报告标黄）
            acc = (row.get('blast_top_hit') or '').strip()
            rm = acc_meta.get(acc, {})
            sp_key = (tax.get('species') or rm.get('Species', '')
                      or '').strip().lower()
            known = list(sp_hosts.get(sp_key, []))
            acc_host = (rm.get('Host') or '').strip()
            if acc_host and acc_host not in known:
                known.insert(0, acc_host)
            known_cats = {meta_cat.get(h.strip().lower())
                          for h in known if h.strip()}
            known_cats.discard(None)
            known_cats.discard('Unknown')
            if final == 'Unknown' or not known_cats:
                check = ''
            elif final in known_cats:
                check = 'OK'
            else:
                check = 'WARN'

            rows.append({
                'contig': contig,
                'length': (row.get('length') or '').strip(),
                'kunpeng_taxid': kt if tax_src == 'kunpeng' else '',
                'class_source': tax_src,
                'accession': acc,
                'species': tax.get('species', '') or
                           (row.get('blast_species') or ''),
                'genus': tax.get('genus', ''),
                'family': tax.get('family', '') or
                          (row.get('blast_family') or ''),
                'topology': (rm.get('Topology') or '').strip(),
                'molecule_type': (rm.get('Molecule_type') or '').strip(),
                'nuc_completeness': (rm.get('Nuc_Completeness') or '').strip(),
                'host_ictv': ictv,
                'confidence_level': ictv_conf,
                'integrated_confidence': hit['Integrated_Confidence'],
                'determination': hit['Determination_Level'],
                'host_meta': meta_name,
                'host_meta_cat': meta,
                'known_hosts': '; '.join(known[:8]) +
                               ('…' if len(known) > 8 else ''),
                'final_host': final,
                'host_check': check,
                'decision_method': method,
            })

    # 汇总
    cats = {}
    for r in rows:
        cats[r['final_host']] = cats.get(r['final_host'], 0) + 1
    methods = {}
    for r in rows:
        methods[r['decision_method']] = methods.get(r['decision_method'], 0) + 1

    # 明细 TSV
    fields = list(rows[0].keys()) if rows else ['contig']
    with safe_open(os.path.join(out_dir, 'host_prediction.tsv'), 'wt') as f:
        f.write('\t'.join(fields) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(k, '')) for k in fields) + '\n')

    # 类别汇总 TSV
    with safe_open(os.path.join(out_dir, 'host_summary.tsv'), 'wt') as f:
        f.write('Host_Category\tContigs\n')
        for k in sorted(cats, key=lambda x: -cats[x]):
            f.write(f'{k}\t{cats[k]}\n')

    # 按宿主拆分病毒 contigs FASTA
    if os.path.isfile(viral_fa):
        seqs = {}
        cid = None
        with safe_open(viral_fa) as f:
            for line in f:
                if line.startswith('>'):
                    cid = line[1:].split()[0]
                    seqs[cid] = []
                elif cid:
                    seqs[cid].append(line)
        by_host = {}
        for r in rows:
            if r['contig'] in seqs and r['final_host'] != 'Unknown':
                by_host.setdefault(r['final_host'], []).append(
                    (r['contig'], ''.join(seqs[r['contig']])))
        for host, items in by_host.items():
            safe_name = ''.join(c if c.isalnum() else '_' for c in host)
            with safe_open(os.path.join(out_dir,
                                        f'{safe_name}.classified.fasta'),
                           'wt') as f:
                for c, s in items:
                    f.write(f'>{c}\n{s}\n')
        if logger:
            logger.log(f"按宿主拆分 FASTA: {len(by_host)} 个类别")
    elif logger:
        logger.log("未找到 viral_contigs.fasta，跳过按宿主拆分 FASTA")

    summary = {
        'stage': step,
        'n_contigs': len(rows),
        'n_host_warn': sum(1 for r in rows if r.get('host_check') == 'WARN'),
        'categories': dict(sorted(cats.items(), key=lambda x: -x[1])),
        'decision_methods': methods,
        'class_source': src_stat,
        'tsv': 'host_prediction.tsv',
        'summary_tsv': 'host_summary.tsv',
        'rows_preview': sorted(
            [{'contig': r['contig'], 'species': r['species'],
              'family': r['family'], 'final_host': r['final_host'],
              'confidence': r['confidence_level'],
              'host_check': r.get('host_check', ''),
              'method': r['decision_method']} for r in rows],
            key=lambda x: x['contig'])[:20],
    }
    with safe_open(summary_file, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    mark_step_done(out_dir, step)

    # 桑基图 / 旭日图（独立 HTML；报告里另嵌）
    try:
        from .viz import fig_host_sankey, fig_host_sunburst, _plotly_js_path
        _plotly_js_path(out_dir)
        sk = fig_host_sankey(rows)
        if sk:
            sk.write_html(check_path(os.path.join(out_dir, 'sankey_host.html'),
                                     must_exist=False, in_platform=True),
                          include_plotlyjs='plotly.min.js' if os.path.isfile(
                              os.path.join(out_dir, 'plotly.min.js'))
                          else 'cdn', full_html=True)
        sb = fig_host_sunburst(rows)
        if sb:
            sb.write_html(check_path(os.path.join(
                out_dir, 'sunburst_host.html'), must_exist=False,
                in_platform=True),
                include_plotlyjs=False, full_html=True)
    except Exception as e:
        if logger:
            logger.log(f"宿主桑基/旭日图生成失败: {e}", "WARN")

    if logger:
        logger.log(f"宿主预测完成: {len(rows)} 条 contigs → "
                   + ', '.join(f'{k} {v}' for k, v in
                               list(summary['categories'].items())[:6]))
    return summary

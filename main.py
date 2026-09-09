# -*- coding: utf-8 -*-
"""
植物病毒分析平台 - 命令行入口
用法:
  python main.py init-taxonomy
  python main.py build-host-db --genome host-db/112863_Lycium_barbarum/genome.fa --taxid 112863
  python main.py build-virus-db --fasta virus-db/final.cluster.ref.fasta --info virus-db/final.cluster.ref_info.tsv
  python main.py analyze --r1 a_R1.fq.gz --r2 a_R2.fq.gz --sample NX-5
  python main.py report --sample NX-5        # 重新生成可视化报告
  python main.py logan-create --name 查询名 --sample NX-5   # LOGAN 溯源查询
  python main.py gb-dl --acc "NC_001367,NC_002692" -n tobamo  # 下载 GenBank 集合
  python main.py compare -n tobamo           # 同属病毒共线性比较
图形界面: 双击 启动平台.bat
"""
import os
import re
import sys
import argparse

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from vp.config import get_config, DIRS, PLATFORM_ROOT, engine_cmd
from vp.utils import TaskLogger, check_path
from vp.pipeline import DEFAULT_ANALYZE_STAGES


def make_logger(name, echo=True):
    log_file = os.path.join(DIRS['logs'], f'{name}_{__import__("time").strftime("%Y%m%d_%H%M%S")}.log')
    return TaskLogger(log_file, echo=echo)


def cmd_init_taxonomy(args):
    from vp.taxonomy import prepare_taxonomy, taxonomy_ready
    logger = make_logger('taxonomy')
    prepare_taxonomy(logger=logger, force=args.force)
    print("\n✔ Taxonomy 就绪:", taxonomy_ready())
    logger.close()


def cmd_build_host_db(args):
    from vp.kunpeng import build_host_db, db_ready
    logger = make_logger('build_host_db')
    check_path(args.genome, must_exist=True)
    # 建库内存峰值与线程数成正比（约 1GB/线程），默认限 8
    threads = args.threads if args.threads else min(8, get_config().threads)
    build_host_db(args.genome, args.taxid,
                  hash_capacity=args.hash_capacity, threads=threads,
                  logger=logger, rebuild=args.rebuild)
    print("\n✔ 宿主库就绪" if db_ready(get_config().databases['host']) else "\n✘ 建库失败")
    logger.close()


def cmd_build_virus_db(args):
    from vp.kunpeng import build_virus_db, db_ready
    logger = make_logger('build_virus_db')
    build_virus_db(args.fasta, args.info,
                   hash_capacity=args.hash_capacity, threads=args.threads,
                   logger=logger, rebuild=args.rebuild)
    print("\n✔ 病毒库就绪" if db_ready(get_config().databases['virus']) else "\n✘ 建库失败")
    logger.close()


def cmd_analyze(args):
    from vp.pipeline import run_analysis
    check_path(args.r1, must_exist=True)
    if args.r2:
        check_path(args.r2, must_exist=True)
    logger = make_logger(f'analyze_{args.sample}')
    run_analysis(
        sample=args.sample, r1=args.r1, r2=args.r2,
        stages=args.stages.split(','),
        db_host=args.db_host, db_virus=args.db_virus,
        threads=args.threads, confidence=args.confidence,
        assembly_mode=args.assembly_mode, memory_gb=args.memory,
        subsample=args.subsample, min_contig_len=args.min_contig_len,
        top_n_refs=args.top_n_refs, tree_tool=args.tree_tool,
        tree_sampling=args.tree_sampling,
        ncbi_refs=args.ncbi_refs,
        primer_mode=args.primer_mode, do_specificity=args.specificity,
        min_orf_aa=args.min_orf_aa, force=args.force,
        do_fq2fa=not args.no_fq2fa,
        gbdraw_max=args.gbdraw_max,
        gbdraw_fasta=args.gbdraw_fasta, gbdraw_ann=args.gbdraw_ann,
        plot_engine=args.plot_engine,
        logger=logger)
    logger.close()


def cmd_ncbi_dl(args):
    from vp.ncbi_download import download_collection, list_collections
    logger = make_logger(f'ncbi_dl_{args.name}')
    res = download_collection(args.term, args.name, db=args.db,
                              max_records=args.max, logger=logger)
    print(f"\n✔ 集合 [{res['name']}]: 命中 {res['total_hits']}，"
          f"新增 {res['downloaded']}，跳过 {res['skipped']}")
    print("已下载集合：")
    for c in list_collections():
        print(f"  {c['name']:20s} {c['n_seqs']:6d} 条  {c['query'][:50]}")
    logger.close()


def cmd_ncbi_list(args):
    from vp.ncbi_download import list_collections
    cols = list_collections()
    if not cols:
        print("（无已下载集合，先用 ncbi-dl 下载）")
        return
    for c in cols:
        print(f"  {c['name']:20s} {c['n_seqs']:6d} 条  db={c['db']}  {c['date']}  "
              f"{c['query'][:60]}")


def cmd_gb_dl(args):
    from vp.gb_collection import download_gb_collection, list_gb_collections
    logger = make_logger(f'gb_dl_{args.name}')
    res = download_gb_collection(args.name, term=args.term,
                                 accessions=args.acc, max_records=args.max,
                                 logger=logger)
    tail = f"，未找到 {len(res['missing'])}" if res['missing'] else ''
    print(f"\n✔ 集合 [{res['name']}]: 命中 {res['total_hits']}，"
          f"新增 {res['downloaded']}，跳过 {res['skipped']}{tail}")
    for c in list_gb_collections():
        print(f"  {c['name']:22s} {c['n_records']:4d} 条  {c['date']}  "
              f"{c['source']}")
    logger.close()


def cmd_gb_import(args):
    from vp.gb_collection import import_local_gb
    logger = make_logger(f'gb_import_{args.name}')
    res = import_local_gb(args.name, args.files, logger=logger)
    print(f"\n✔ 集合 [{res['name']}]: 导入 {res['imported']} 条，"
          f"跳过 {res['skipped']} 条")
    logger.close()


def cmd_gb_list(args):
    from vp.gb_collection import list_gb_collections
    cols = list_gb_collections()
    if not cols:
        print('（无 GenBank 集合，先用 gb-dl 下载或 gb-import 导入）')
        return
    for c in cols:
        if c['source'] == 'query':
            src = c['term']
        elif c['source'] == 'accessions':
            src = f"{c['accessions']} 个 accession"
        else:
            src = '本机导入'
        print(f"  {c['name']:22s} {c['n_records']:4d} 条  {src[:52]:52s} "
              f"{c['date']}")


def cmd_gb_check(args):
    from vp.gb_collection import inspect_collection
    st = inspect_collection(args.name)
    print(f"集合 [{st['name']}]（{st['dir']}）: {len(st['records'])} 条记录\n")
    print(f"{'accession':18s}{'长度':>10s}{'CDS':>5s}{'成熟肽':>6s}  物种")
    for r in st['records']:
        print(f"{r['acc'][:17]:18s}{r['length']:10d}{r['cds_count']:5d}"
              f"{r['mat_peptides']:6d}  {(r['organism'] or r['title'])[:44]}")
    if st['warnings']:
        print('\n⚠ 警告:')
        for w in st['warnings']:
            print(f'  - {w}')


def cmd_compare(args):
    from vp.synteny import run_comparison
    logger = make_logger(f'compare_{args.name or "adhoc"}')
    files = None
    if args.files:
        files = [p for p in args.files.replace(';', ',').split(',')
                 if p.strip()]
    res = run_comparison(name=args.name, files=files,
                         min_ident=args.min_ident, min_cov=args.min_cov,
                         mat_peptide=not args.no_mat_peptide,
                         order=(args.order.split(',') if args.order
                                else None),
                         labels=not args.no_labels, style=args.style,
                         lovis4u_pdf=args.lovis4u, threads=args.threads,
                         logger=logger)
    print(f"\n✔ 比较: {res['n_genomes']} 基因组 · {res['n_genes']} 基因 · "
          f"{res['n_clusters']} 个同源家族"
          + (f" · {res['n_singletons']} 个无同源基因"
             if res['n_singletons'] else ''))
    for w in res['warnings']:
        print(f'  ⚠ {w}')
    print(f"共线性图 : {res['files']['html']}")
    print(f"矢量图   : {res['files']['svg']}")
    if res['files'].get('lovis4u_pdf'):
        print(f"LoVis4u  : {res['files']['lovis4u_pdf']}")
    print(f"家族清单 : {res['files']['clusters']}")
    print(f"共享比例 : {res['files']['similarity']}")
    logger.close()


def cmd_report(args):
    from vp.pipeline import run_report_only
    logger = make_logger(f'report_{args.sample}')
    run_report_only(args.sample, logger=logger, force=args.force)
    logger.close()


def cmd_host_analysis(args):
    import os
    from vp.config import DIRS
    from vp.pipeline import _safe_sample_name
    from vp.host_analysis import predict_hosts
    sample = _safe_sample_name(args.sample)
    sample_dir = os.path.join(DIRS['results'], sample)
    logger = make_logger(f'hostana_{sample}')
    s = predict_hosts(sample_dir, logger=logger, force=args.force)
    logger.close()
    print(f"宿主预测完成: {s['n_contigs']} 条 contigs")
    for cat, n in s['categories'].items():
        print(f"  {cat:15s}: {n}")


def cmd_orfa(args):
    import os
    from vp.config import DIRS
    from vp.pipeline import _safe_sample_name
    from vp.orf_annot import run_orf_annotation
    sample = _safe_sample_name(args.sample)
    sample_dir = os.path.join(DIRS['results'], sample)
    logger = make_logger(f'orfa_{sample}')
    s = run_orf_annotation(sample_dir, threads=args.threads, logger=logger,
                           force=args.force)
    logger.close()
    print(f"ORF 功能注释完成: {s['n_annotated']}/{s['n_orfs']} 个 ORF 获得"
          f"有效注释，覆盖 {s['n_families']} 个病毒科（引擎 {s['engine']}）")
    for c, n in list(s['categories'].items())[:10]:
        print(f"  {c:20s}: {n}")


def cmd_logan_create(args):
    from vp.logan_trace import create_job, get_segment_fasta
    d = create_job(args.name, sample=args.sample or None,
                   contig_ids=(args.contigs.split(',') if args.contigs else None),
                   pasted=args.fasta, n_seg=args.segments)
    print(f"查询任务已创建: {d['name']}（{d['n_segments']} 个片段）")
    print("请把每个片段的序列提交到 https://logan-search.org/dashboard"
          "（Groups 建议 All_No_viral_human，建议填邮箱拿下载链接）：")
    for s in d['segments']:
        _h, seq = get_segment_fasta(d['name'], s['index'])
        fa = os.path.join(DIRS['logan'], d['name'], s['file'])
        print(f"  s{s['index']}  {s['contig']} {s['start']}-{s['end']}  → {fa}")
    print("拿到结果表后：python main.py logan-import --name "
          f"{d['name']} --segment 1 --result 结果.tsv")


def cmd_logan_import(args):
    from vp.logan_trace import import_result
    raw = open(check_path(args.result, must_exist=True), 'rb').read()
    d = import_result(args.name, args.segment, args.result, raw)
    print(f"已导入片段 s{args.segment}（{args.result}）；"
          f"任务状态: {d['status']}（{d['n_imported']}/{d['n_segments']} 片段）")
    print(f"溯源报告: {os.path.join(DIRS['logan'], d['name'], 'trace_report.html')}")


def cmd_logan_jobs(args):
    from vp.logan_trace import list_jobs
    jobs = list_jobs()
    if not jobs:
        print("暂无 LOGAN 溯源查询任务")
        return
    for j in jobs:
        print(f"{j['name']:30s} 样品={j['sample'] or '—':12s} "
              f"片段 {j['n_imported']}/{j['n_segments']} 已导入  "
              f"{'报告✔' if j['has_report'] else ''}")


def cmd_logan_batch(args):
    from vp.logan_trace import batch_submit
    emails = [e.strip() for e in args.email.split(',') if e.strip()]
    logger = make_logger(f'logan_batch_{args.name}')
    r = batch_submit(args.name, emails, group=args.group,
                     headless=not args.show_browser,
                     first_wait=args.first_wait, max_wait=args.max_wait,
                     logger=logger)
    logger.close()
    print(f"批量完成: 待提交 {r['n_pending']} 片段, 回收导入 {r['n_imported']} 个"
          f"（未落定的片段可再次运行补漏）")


def cmd_submit_list(args):
    from vp.ncbi_submit import store
    tabs = store.list_tables()
    if not tabs:
        print("暂无提交项目（submissions/ 为空）。新建：python main.py submit-init --name demo --demo")
        return
    for t in tabs:
        extra = [f for f in t['files'] if f != 'unified_metadata.csv']
        print(f"{t['name']:30s} {t['rows']:4d} 行  产物: "
              f"{', '.join(extra) if extra else '—（submit-export 生成）'}")


def cmd_submit_init(args):
    import subprocess as _sp
    from vp.ncbi_submit import store
    from vp.ncbi_submit import unified_metadata as _um_path
    if args.demo:
        store.create_table(args.name, sample='demo')
        print(f"提交项目已创建（示例数据）: submissions/{args.name}/")
        return
    if args.import_csv:
        _, n = store.import_table(args.name, check_path(args.import_csv,
                                                        must_exist=True))
        print(f"提交项目已创建（导入 {n} 行）: submissions/{args.name}/")
        return
    if not args.taxonomy:
        raise SystemExit("需要 --taxonomy <taxonomy.tsv>（contig<TAB>taxonomy），"
                         "或 --demo / --import-csv")
    d = store.table_dir(args.name)
    os.makedirs(d, exist_ok=True)
    cmd = engine_cmd(_um_path.__file__,
           '--taxonomy', check_path(args.taxonomy, must_exist=True),
           '--run-title', args.name, '-o', d)
    for flag, val in (('--metadata', args.metadata), ('--authors', args.authors),
                      ('--title', args.title), ('--bioproject', args.bioproject),
                      ('--host', args.host), ('--lat-lon', args.lat_lon),
                      ('--sequencer', args.sequencer), ('--assembler', args.assembler),
                      ('--coverage', args.coverage)):
        if val:
            cmd += [flag, val]
    r = _sp.run(cmd)
    if r.returncode == 0:
        print(f"提交项目已创建: submissions/{args.name}/（unified_metadata.csv + source.src 等）")
        print(f"继续: python main.py submit-validate --name {args.name} && "
              f"python main.py submit-export --name {args.name}")


def cmd_submit_validate(args):
    from vp.ncbi_submit import store
    issues = store.validate_table(args.name)
    if not issues:
        print("✓ 所有必填字段已填写，可以提交")
        return
    print(f"校验：{len(issues)} 个问题")
    for it in issues:
        if it['missing_count'] == -1:
            print(f"  [缺少列] {it['column']} — {it['desc']}")
        else:
            ex = ', '.join(repr(e) for e in it['examples'][:3])
            print(f"  [{it['missing_count']} 条未填] {it['column']} — {it['desc']}\n"
                  f"      例: {ex}")
    print("修复占位符后再提交 GenBank（Web 端「提交准备」页可在线编辑）")


def cmd_submit_fill(args):
    from vp.ncbi_submit import store
    n = store.batch_fill(args.name, args.column, args.value, old_value=args.old)
    print(f"已更新 '{args.column}' 列 {n} 个单元格")


def cmd_submit_export(args):
    from vp.ncbi_submit import store
    out = store.export_files(args.name, assembler=args.assembler,
                             sequencer=args.sequencer,
                             enrichment=args.enrichment)
    print("产物已生成：")
    for k, v in out.items():
        print(f"  {k:20s} {v}")
    print("提交：NCBI BankIt 上传 source.src / .fsa+.tbl，或 Sequin 导入；"
          "BioSample 用 biosample_template.tsv 批量注册。")


def cmd_submit_sbt(args):
    from vp.ncbi_submit import store
    fields = {'last': args.last, 'first': args.first, 'middle': args.middle or '',
              'affil': args.affil, 'div': args.div or '', 'city': args.city,
              'sub': args.sub or '', 'country': args.country,
              'street': args.street or '', 'email': args.email,
              'postal': args.postal or ''}
    extra = [ln for ln in (args.extra_authors or '').splitlines() if ln.strip()] \
        if args.extra_authors else []
    tpl, path = store.generate_sbt(fields, extra_authors=extra,
                                   title=args.title or '', out_name=args.name)
    print(f"template.sbt 已生成: {path}（{1 + len(extra)} 位作者）")


def cmd_meta_search(args):
    import subprocess as _sp
    from vp.public_meta import search_engine
    cmd = engine_cmd(search_engine.__file__,
           '-q', args.species, '-s', args.source, '--db', args.db)
    if args.out:
        cmd += ['-o', args.out]
    if args.no_detailed:
        cmd.append('--no-detailed')
    if args.ncbi_api or os.environ.get('NCBI_API_KEY'):
        cmd += ['--ncbi-api', args.ncbi_api or os.environ['NCBI_API_KEY']]
    if args.deepseek_api:
        cmd += ['--deepseek-api', args.deepseek_api]
    r = _sp.run(cmd)
    if r.returncode == 0:
        print("下一步: python main.py meta-info --runs <输出目录>/sra.list")


def cmd_meta_info(args):
    import subprocess as _sp
    from vp.public_meta import info_engine
    # --runs 既可以是编号列表文件，也可以是单个 SRR/CRR 编号
    if os.path.isfile(args.runs):
        runs_arg = check_path(args.runs, must_exist=True)
    elif re.fullmatch(r'[A-Za-z]\w+', args.runs or ''):
        runs_arg = args.runs
    else:
        raise SystemExit(f"--runs 需为编号列表文件或 Run 编号: {args.runs}")
    cmd = engine_cmd(info_engine.__file__,
           '-i', runs_arg, '-m', args.mode, '-t', str(args.threads))
    if args.out:
        cmd += ['-o', args.out]
    if args.fill_date:
        cmd.append('--fill-date')
    if args.deepseek_api:
        cmd += ['--deepseek-api', args.deepseek_api]
    if args.ncbi_api:
        cmd += ['--ncbi-api', args.ncbi_api]
    _sp.run(cmd)


def cmd_meta_plot(args):
    from vp.public_meta.landscape_plot import plot_sci_landscape
    plot_sci_landscape(check_path(args.input, must_exist=True), args.out)


def cmd_host_genome(args):
    from vp.public_meta.host_genome import main as _hg_main
    sys.argv = ['host_genome.py', '--species', args.species]
    if args.out:
        sys.argv += ['--outdir', args.out]
    if args.include_organelles:
        sys.argv.append('--include-organelles')
    if args.ncbi_api:
        sys.argv += ['--ncbi-api', args.ncbi_api]
    if args.min_length:
        sys.argv += ['--min-length', str(args.min_length)]
    if args.verify_only:
        sys.argv.append('--verify-only')
    if args.skip_datasets:
        sys.argv.append('--skip-datasets')
    _hg_main()


def cmd_ref_status(args):
    from vp import virus_ref
    st = virus_ref.status()
    if not st['dir']:
        print("✘ 参考库未部署（databases/virus_ref 已移除）")
        return
    print(f"参考库目录 : {st['dir']}")
    print(f"版本       : {st['version']}（最后增量 {st['last_incremental']}）")
    print(f"数据来源   : ICTV {st['source_ictv']} / NCBI {st['source_ncbi']}")
    print(f"非冗余参考 : {st['n_ref']} 条")
    print(f"完整基因组 : {st['n_complete']} 条")
    print(f"规范化元数据: {'✔' if st['meta_ready'] else '✘（首次使用时自动构建）'}")


def cmd_ictv_update(args):
    from vp.ictv_db import update
    logger = make_logger('ictv_update')
    update(logger=logger, xlsx=args.xlsx)
    cmd_ictv_status(args)
    logger.close()


def cmd_ictv_status(args):
    from vp import ictv_db
    st = ictv_db.status()
    if not st['dir']:
        print("✘ 未找到 databases/tax/ictv/（ICTV 参考库未部署）")
        return
    print(f"ICTV 库目录 : {st['dir']}")
    if not st['taxa_ready']:
        print("  ✘ taxa.txt 未生成（用 ictv-update 更新 VMR 后自动解析）")
        return
    print(f"MSL 版本    : {st['msl'] or '?'}（VMR: {st['xlsx'] or '?'}）")
    print(f"解析时间    : {st['parsed'] or '?'}")
    print(f"accession   : {st['n_rows']} 条（{st['n_genus']} 属 / {st['n_family']} 科）")
    print(f"本地覆盖    : 本地已有 {st['n_in_acvirus']} 条，"
          f"gb 缓存已下 {st['n_gb_cache']} 条")


def cmd_ictv_refs(args):
    from vp import ictv_db
    if not any((args.genus, args.family, args.species)):
        raise SystemExit('至少指定 --genus / --family / --species 之一')
    rows, total = ictv_db.select_refs(
        genus=args.genus, family=args.family, species=args.species,
        limit=args.limit, genome='any' if args.any else 'complete')
    scope = args.genus or args.family or args.species
    print(f"[{scope}] 命中 {total} 条"
          + ('' if args.any else '（仅 Complete genome）')
          + f"，按本地优先显示前 {len(rows)} 条:\n")
    print(f"{'accession':12s}{'来源':11s}{'基因组':22s}{'物种':40s}宿主组")
    for r in rows:
        print(f"{r['Genbank']:12s}{r['Source']:11s}"
              f"{(r['Genome_Coverage'] or '')[:20]:22s}"
              f"{(r['Species'] or '')[:38]:40s}{(r['Host_Source'] or '')[:16]}")
    need = [r['Genbank'] for r in rows if r['Source'] == 'ncbi']
    if not need:
        return
    print(f"\n其中 {len(need)} 条本地没有（本地缓存未命中）")
    if not args.download:
        print("预览模式未下载；加 --download 执行按需下载（需联网）")
        return
    logger = make_logger('ictv_refs')
    res = ictv_db.ensure_gb(need, logger=logger)
    tail = f"，未找到 {len(res['missing'])}" if res['missing'] else ''
    print(f"✔ 下载 {res['downloaded']} 条，本地已有跳过 {res['skipped']} 条{tail}")
    print(f"缓存 fasta: {res['fasta']}")
    logger.close()


def cmd_tools(args):
    cfg = get_config()
    print(f"线程数: {cfg.threads}")
    for name, path in sorted(cfg.tool_status().items()):
        print(f"  {name:12s} {'✔ ' + path if path else '✘ 未找到'}")
    print(f"宿主库: {cfg.databases['host']}")
    print(f"病毒库: {cfg.databases['virus']}")


def cmd_db_migrate(args):
    """数据库迁移/对接：把 databases/host-db/virus-db 复制或搬移到外置位置，
    并切换 database_root，实现软件与数据库分离部署。"""
    from vp.db_migrate import migrate, check
    if args.check:
        r = check(args.to)
        print(('\n✔ ' if r.get('ok') else '\n✘ ') + (r.get('detail')
              or r.get('error', '')))
        return None if r.get('ok') else False
    r = migrate(args.to, mode=args.mode, dry_run=args.dry_run)
    if r.get('dry_run'):
        print('\n' + r['detail'])
        return None
    print(('\n✔ ' if r.get('ok') else '\n✘ 迁移失败: ') + (r.get('detail')
          or r.get('error', '')))
    return None if r.get('ok') else False


def cmd_tool_runs(args):
    """输出目录管理：status / organize / archive / clean 四动作。"""
    get_config()  # 确保 DIRS 已按 platform.json（含自定义输出根）就绪
    from vp import tool_runs_admin as tra
    act = args.action
    if act == 'status':
        info = tra.status()
        print(f"根目录: {info['root']}")
        act_sz = sum(a['size'] for a in info['active'])
        print(f"\n活动区（动态运行，最新在前，共 {len(info['active'])} 个，"
              f"{tra._fmt(act_sz)}）:")
        for a in info['active'][:15]:
            print(f"  {a['name']:44s} {a['size_h']:>10s}  {a['files']} 文件")
        if len(info['active']) > 15:
            print(f"  … 另有 {len(info['active']) - 15} 个，详见 GUI 运行目录页")
        print('\n固定区:')
        for d, f in info['fixed'].items():
            print(f"  {d:12s} {f['dirs']:>3d} 子目录  {f['files']:>4d} 文件  "
                  f"{f['size_h']}")
        return None
    if act == 'organize':
        r = tra.organize(dry_run=args.dry_run)
        print(f"\n{'预演' if args.dry_run else '完成'}: 归位 {r['moved']} 项，"
              f"跳过 {r['skipped']} 项")
        return None
    if act == 'archive':
        r = tra.archive(before=args.before, older_days=args.older_days,
                        dry_run=args.dry_run)
        print(f"\n{'预演' if args.dry_run else '完成'}: 归档 {r['archived']} 个运行")
        return None
    if act == 'clean':
        if args.active_days is None and args.archive_days is None:
            print('✘ clean 需至少给一个期限: --active-days N 和/或 --archive-days N')
            return False
        r = tra.clean(active_days=args.active_days,
                      archive_days=args.archive_days, dry_run=args.dry_run)
        print(f"\n{'预演' if args.dry_run else '完成'}: 删除 {r['removed']} 项，"
              f"释放 {r['freed_h']}")
        return None


def cmd_selfcheck(args):
    """环境自检：模块导入 / 外部工具 / 数据库 / 磁盘空间，一次摸底。"""
    import importlib
    import time as _t
    ok_all = True

    # 1) 核心模块导入
    print('── ① 模块导入 ──')
    mods = ['vp.config', 'vp.utils', 'vp.taxonomy', 'vp.kunpeng', 'vp.preprocess',
            'vp.host_removal', 'vp.virus_screen', 'vp.assembly', 'vp.host_analysis',
            'vp.orf', 'vp.orf_annot', 'vp.phylo', 'vp.primer',
            'vp.viz', 'vp.pipeline', 'vp.msa_view']
    t0 = _t.time()
    fails = []
    for m in mods:
        try:
            importlib.import_module(m)
        except Exception as e:
            fails.append((m, repr(e)))
    print(f'  {"✔" if not fails else "✘"} {len(mods) - len(fails)}/{len(mods)} '
          f'({_t.time() - t0:.1f}s)')
    for m, e in fails:
        print(f'    ✘ {m}: {e}')
    ok_all &= not fails

    # 2) 外部工具
    print('── ② 外部工具 ──')
    cfg = get_config()
    missing = [k for k, v in sorted(cfg.tool_status().items()) if not v]
    total = len(cfg.tool_status())
    print(f'  {"✔" if not missing else "△"} {total - len(missing)}/{total} 可用'
          + (f'（缺: {", ".join(missing)}）' if missing else ''))
    if missing:
        print('    缺失工具只影响对应步骤（可选步骤自动降级）；'

              '必装件 kunpeng/SPAdes/BLAST 缺失时对应管道无法运行。')

    # 3) 数据库
    print('── ③ 数据库 ──')
    from vp.kunpeng import db_ready
    from vp.taxonomy import taxonomy_ready
    checks = [
        ('NCBI Taxonomy', taxonomy_ready()),
        ('宿主库 host_db', db_ready(cfg.databases['host'])),
        ('病毒库 virus_db', db_ready(cfg.databases['virus'])),
    ]
    for name, ok in checks:
        print(f'  {"✔" if ok else "✘"} {name}')
        ok_all &= ok

    # 4) 磁盘空间
    print('── ④ 磁盘空间 ──')
    import shutil
    for key, d in (('平台根', PLATFORM_ROOT), ('结果', DIRS['results'])):
        try:
            free = shutil.disk_usage(d).free
            print(f'  {"✔" if free > 10e9 else "△"} {key}: 剩余 {free / 1e9:.1f} GB'
                  + ('' if free > 10e9 else '（建议保留 >10GB）'))
        except OSError as e:
            print(f'  ✘ {key}: {e}')

    # 5) 结果目录里的样品
    print('── ⑤ 样品 ──')
    try:
        samples = sorted(os.listdir(DIRS['results']))
        samples = [s for s in samples
                   if os.path.isdir(os.path.join(DIRS['results'], s))]
        print(f'  共 {len(samples)} 个样品: {", ".join(samples[:8])}'
              + ('…' if len(samples) > 8 else ''))
    except OSError:
        print('  （results/ 目录尚不存在）')

    # 6) 第三方依赖（对照 requirements.txt；缺必需依赖会影响打包分发）
    print('── ⑥ 第三方依赖 ──')
    try:
        from scripts.audit_deps import analyze as _dep_analyze
        dep = _dep_analyze()
        if not dep['hard_missing']:
            print(f"  ✔ 全部依赖就绪（扫描 {len(dep['third'])} 个，"
                  f"requirements 已覆盖）")
            if dep['missing']:
                print(f"  △ 可选依赖缺失（对应功能降级）: "
                      f"{', '.join(dep['missing'])}")
        else:
            for m in dep['hard_missing']:
                print(f'  ✘ 必需依赖未装: {m}')
            print('  安装: python -m pip install -r requirements.txt')
        ok_all &= not dep['hard_missing']
    except Exception as e:
        print(f'  ✘ 依赖审计失败: {e}')

    print('════════════════════════')
    print('✔ 自检通过' if ok_all else '△ 存在未就绪项（见上），按需处理后重跑本命令')
    return ok_all


def main():
    cfg = get_config()
    p = argparse.ArgumentParser(description="植物病毒分析平台 CLI",
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    sub = p.add_subparsers(dest='cmd', required=True)

    sp = sub.add_parser('init-taxonomy', help='下载/准备 NCBI taxonomy（首次必做）')
    sp.add_argument('--force', action='store_true')
    sp.set_defaults(func=cmd_init_taxonomy)

    sp = sub.add_parser('build-host-db', help='kunpeng 宿主数据库构建')
    sp.add_argument('--genome', required=True, help='宿主基因组 FASTA')
    sp.add_argument('--taxid', type=int, required=True, help='宿主 NCBI TaxID')
    sp.add_argument('--hash-capacity', default='256M')
    sp.add_argument('--threads', type=int, default=0,
                    help='线程数（默认自动=8；内存不足时用 4）')
    sp.add_argument('--rebuild', action='store_true')
    sp.set_defaults(func=cmd_build_host_db)

    sp = sub.add_parser('build-virus-db', help='kunpeng 病毒数据库构建')
    sp.add_argument('--fasta', required=True, help='病毒参考 FASTA')
    sp.add_argument('--info', required=True, help='Accession/Taxid 信息表 TSV')
    sp.add_argument('--hash-capacity', default='64M')
    sp.add_argument('--threads', type=int, default=cfg.threads)
    sp.add_argument('--rebuild', action='store_true')
    sp.set_defaults(func=cmd_build_virus_db)

    sp = sub.add_parser('analyze', help='样品全流程分析')
    sp.add_argument('--r1', required=True)
    sp.add_argument('--r2', default=None)
    sp.add_argument('--sample', required=True)
    sp.add_argument('--stages',
                    default=','.join(DEFAULT_ANALYZE_STAGES),
                    help='阶段: fastp,fq2fa,host,virus,assembly,hostana,'
                         'orf,orfa,phylo,primer,gbdraw,report')
    sp.add_argument('--db-host', default=cfg.databases['host'])
    sp.add_argument('--db-virus', default=cfg.databases['virus'])
    sp.add_argument('--threads', type=int, default=cfg.threads)
    sp.add_argument('--confidence', type=float, default=0.0)
    sp.add_argument('--assembly-mode', default='metaviral',
                    choices=['metaviral', 'meta', 'rna', 'isolate'])
    sp.add_argument('--subsample', type=int, default=0,
                    help='仅取前 N 对 reads（快速测试用），0=全部')
    sp.add_argument('--min-contig-len', type=int, default=500)
    sp.add_argument('--top-n-refs', type=int, default=10, help='进化分析取近缘参考数')
    sp.add_argument('--tree-tool', default='fasttree', choices=['fasttree', 'iqtree'])
    sp.add_argument('--tree-sampling', default='blast',
                    choices=['blast', 'macro', 'genus', 'lineage'],
                    help='参考挑选: blast=按比对hits(默认); macro=同科建树'
                         '(目标属+同科各属背景); genus=属级树; lineage=种级树'
                         '（后三者需分类元数据；acvirus_db/virus_ref 已移除）')
    sp.add_argument('--ncbi-refs', default=None,
                    help='NCBI 参考集合名（ncbi-dl 下载，逗号分隔可多个），'
                         '追加进进化树比对')
    sp.add_argument('--primer-mode', default='conserved', choices=['conserved', 'plain'])
    sp.add_argument('--specificity', action='store_true', help='引物宿主特异性 BLAST 检查(需建宿主BLAST库)')
    sp.add_argument('--min-orf-aa', type=int, default=100)
    sp.add_argument('--memory', type=int, default=64, help='SPAdes 内存上限 GB')
    sp.add_argument('--no-fq2fa', action='store_true',
                    help='跳过 FASTQ→FASTA 预转换（默认开启，推荐保留）')
    sp.add_argument('--gbdraw-max', type=int, default=12,
                    help='基因组图出图 contigs 上限（默认 12）')
    sp.add_argument('--gbdraw-fasta', default=None,
                    help='gbdraw 自备 FASTA（默认用③的病毒 contigs）')
    sp.add_argument('--gbdraw-ann', default=None,
                    help='gbdraw 自备注释（GFF3 或 GenBank .gb/.gbk）')
    sp.add_argument('--plot-engine', default='auto',
                    choices=['auto', 'gbdraw', 'dfv'],
                    help='⑨基因组图引擎：auto=gbdraw 优先、缺则用 '
                         'dna_features_viewer（默认 auto）')
    sp.add_argument('--force', action='store_true')
    sp.set_defaults(func=cmd_analyze)

    sp = sub.add_parser('ncbi-dl', help='NCBI Entrez 批量下载参考序列（进化树扩充）')
    sp.add_argument('term', help='Entrez 检索式，如 "Tobamovirus[ORGN] AND '
                                 'complete genome[TITL]"')
    sp.add_argument('-n', '--name', required=True, help='集合名（ analyze --ncbi-refs 引用）')
    sp.add_argument('--db', default='nucleotide', choices=['nucleotide', 'protein'])
    sp.add_argument('--max', type=int, default=100, help='最多下载条数（默认 100）')
    sp.set_defaults(func=cmd_ncbi_dl)

    sp = sub.add_parser('ncbi-list', help='列出已下载的 NCBI 参考集合')
    sp.set_defaults(func=cmd_ncbi_list)

    sp = sub.add_parser('gb-dl',
                        help='下载 GenBank 集合（同属共线性比较输入）')
    sp.add_argument('-n', '--name', required=True, help='集合名')
    sp.add_argument('--term', default=None,
                    help='Entrez 检索式，如 "Potyvirus[ORGN] AND '
                         'complete genome[TITL]"（与 --acc 二选一）')
    sp.add_argument('--acc', default=None,
                    help='accession 列表（逗号/空格分隔，支持版本号）')
    sp.add_argument('--max', type=int, default=50, help='最多下载数（检索式时）')
    sp.set_defaults(func=cmd_gb_dl)

    sp = sub.add_parser('gb-import',
                        help='导入本机 GenBank 文件为集合（多记录自动拆分）')
    sp.add_argument('-n', '--name', required=True, help='集合名')
    sp.add_argument('--files', required=True,
                    help='.gb/.gbk 文件路径（逗号分隔，支持 .gz）')
    sp.set_defaults(func=cmd_gb_import)

    sp = sub.add_parser('gb-list', help='列出 GenBank 集合')
    sp.set_defaults(func=cmd_gb_list)

    sp = sub.add_parser('gb-check',
                        help='巡检 GenBank 集合（记录数/CDS 数/警告）')
    sp.add_argument('-n', '--name', required=True, help='集合名')
    sp.set_defaults(func=cmd_gb_check)

    sp = sub.add_parser('compare',
                        help='同属病毒共线性比较（MMseqs2 全对全 + 交互式图）')
    sp.add_argument('-n', '--name', default=None,
                    help='GenBank 集合名（gb-dl/gb-import 创建的）')
    sp.add_argument('--files', default=None,
                    help='或直接给 .gb 文件（逗号分隔），可与 --name 并用')
    sp.add_argument('--min-ident', type=float, default=0.30,
                    help='家族聚类最小蛋白一致性（默认 0.30）')
    sp.add_argument('--min-cov', type=float, default=0.50,
                    help='家族聚类最小双侧覆盖度（默认 0.50）')
    sp.add_argument('--no-mat-peptide', action='store_true',
                    help='多聚蛋白基因组不改用 mat_peptide 做基因')
    sp.add_argument('--order', default=None,
                    help='基因组显示顺序（逗号分隔 accession 前缀）')
    sp.add_argument('--no-labels', action='store_true',
                    help='图上不标基因名')
    sp.add_argument('--style', default='lovis', choices=['lovis', 'category'],
                    help='着色风格：lovis=同源家族逐个着色（LoVis4u 画廊风，'
                         '默认）；category=按功能类别着色+图例')
    sp.add_argument('--lovis4u', action='store_true',
                    help='同时调用 LoVis4u（pip install lovis4u）出原生 PDF 图')
    sp.add_argument('--threads', type=int, default=None)
    sp.set_defaults(func=cmd_compare)

    sp = sub.add_parser('report', help='对已有结果重新生成可视化报告')
    sp.add_argument('--sample', required=True)
    sp.add_argument('--force', action='store_true')
    sp.set_defaults(func=cmd_report)

    sp = sub.add_parser('host-analysis',
                        help='对已有 ③组装 结果做 ICTV 宿主预测')
    sp.add_argument('--sample', required=True)
    sp.add_argument('--force', action='store_true')
    sp.set_defaults(func=cmd_host_analysis)


    sp = sub.add_parser('orfa', help='ORF 功能注释（⑥b，RefSeq 病毒蛋白搜索）')
    sp.add_argument('--sample', required=True)
    sp.add_argument('--threads', type=int, default=None)
    sp.add_argument('--force', action='store_true')
    sp.set_defaults(func=cmd_orfa)

    sp = sub.add_parser('logan-create',
                        help='LOGAN 溯源：生成待提交 Logan-Search 的查询片段')
    sp.add_argument('--name', required=True, help='查询名称')
    sp.add_argument('--sample', default=None,
                    help='来源样品（需已有 ③组装 病毒 contigs）')
    sp.add_argument('--contigs', default=None,
                    help='逗号分隔的 contig 名（默认全部）')
    sp.add_argument('--fasta', default=None, help='或直接粘贴 FASTA 文本')
    sp.add_argument('--segments', type=int, default=2,
                    help='每条 contig 切片段数 1-4（默认 2；≤2.5kb 恒为 1 段）')
    sp.set_defaults(func=cmd_logan_create)

    sp = sub.add_parser('logan-import',
                        help='LOGAN 溯源：导入 Logan-Search 结果表并生成报告')
    sp.add_argument('--name', required=True, help='查询名称')
    sp.add_argument('--segment', type=int, default=1, help='片段编号（从 1 起）')
    sp.add_argument('--result', required=True, help='结果表 CSV/TSV 文件')
    sp.set_defaults(func=cmd_logan_import)

    sp = sub.add_parser('logan-jobs', help='LOGAN 溯源：列出查询任务')
    sp.set_defaults(func=cmd_logan_jobs)

    sp = sub.add_parser('logan-batch',
                        help='LOGAN 溯源：Selenium 批量提交全部未导入片段并自动收结果')
    sp.add_argument('--name', required=True, help='查询名称')
    sp.add_argument('--email', required=True,
                    help='通知邮箱，逗号分隔多个则轮换')
    sp.add_argument('--group', default='Fast_No_human',
                    choices=['All', 'All_No_viral_human', 'Fast', 'Fast_No_human',
                             'Fast_No_RefSeq', 'Transcriptomic',
                             'Metatranscriptomic', 'Metagenomic', 'GenBank_RefSeq'])
    sp.add_argument('--show-browser', action='store_true',
                    help='显示浏览器窗口（默认后台 headless）')
    sp.add_argument('--first-wait', type=int, default=300)
    sp.add_argument('--max-wait', type=int, default=1800)
    sp.set_defaults(func=cmd_logan_batch)

    sp = sub.add_parser('submit-list', help='提交准备：列出提交项目')
    sp.set_defaults(func=cmd_submit_list)

    sp = sub.add_parser('submit-init', help='提交准备：新建提交项目（unified_metadata.csv）')
    sp.add_argument('--name', required=True, help='项目名（字母/数字/_/-）')
    sp.add_argument('--taxonomy', help='taxonomy.tsv（contig<TAB>taxonomy）→ 自动生成初表')
    sp.add_argument('--metadata', help='公共元数据表（Core14/Full，自动填日期/地点等）')
    sp.add_argument('--import-csv', help='从已有 unified CSV/TSV/Excel 导入')
    sp.add_argument('--demo', action='store_true', help='用示例数据建表')
    sp.add_argument('--authors', help='作者 "Last, First; ..."')
    sp.add_argument('--title', help='提交标题')
    sp.add_argument('--bioproject', help='BioProject ID')
    sp.add_argument('--host', help='宿主物种名')
    sp.add_argument('--lat-lon', help='经纬度 "38.47 N 106.27 E"')
    sp.add_argument('--sequencer', default='Illumina NovaSeq 6000')
    sp.add_argument('--assembler', default='SPAdes;4.3.0;metaviral')
    sp.add_argument('--coverage', help='覆盖度（如 42.5x）')
    sp.set_defaults(func=cmd_submit_init)

    sp = sub.add_parser('submit-validate', help='提交准备：必填字段校验')
    sp.add_argument('--name', required=True)
    sp.set_defaults(func=cmd_submit_validate)

    sp = sub.add_parser('submit-fill', help='提交准备：批量填充某列占位符/替换旧值')
    sp.add_argument('--name', required=True)
    sp.add_argument('--column', required=True)
    sp.add_argument('--value', required=True, help='新值')
    sp.add_argument('--old', default=None, help='旧值（缺省=填充所有占位符）')
    sp.set_defaults(func=cmd_submit_fill)

    sp = sub.add_parser('submit-export', help='提交准备：生成 source.src/miuvig/assembly/BioSample/report')
    sp.add_argument('--name', required=True)
    sp.add_argument('--assembler', default='SPAdes;4.3.0;metaviral')
    sp.add_argument('--sequencer', default='Illumina NovaSeq 6000')
    sp.add_argument('--enrichment', default='rRNA depletion')
    sp.set_defaults(func=cmd_submit_export)

    sp = sub.add_parser('submit-sbt', help='提交准备：生成 template.sbt（作者/机构）')
    sp.add_argument('--name', required=True)
    sp.add_argument('--last', required=True)
    sp.add_argument('--first', required=True)
    sp.add_argument('--middle', default='')
    sp.add_argument('--affil', required=True, help='机构')
    sp.add_argument('--div', default='', help='院系')
    sp.add_argument('--city', required=True)
    sp.add_argument('--sub', default='', help='省/州')
    sp.add_argument('--country', required=True)
    sp.add_argument('--street', default='')
    sp.add_argument('--email', required=True)
    sp.add_argument('--postal', default='')
    sp.add_argument('--title', default='', help='提交标题')
    sp.add_argument('--extra-authors', default='',
                    help='其他作者，每行一个 "Last, First"（\\n 分隔）')
    sp.set_defaults(func=cmd_submit_sbt)

    sp = sub.add_parser('meta-search', help='公共数据检索：SRA+GSA 双引擎按物种检索 Run')
    sp.add_argument('--species', required=True, help='物种拉丁名（如 "Lycium barbarum"）')
    sp.add_argument('--source', default='TRANSCRIPTOMIC',
                    help='测序类型过滤 TRANSCRIPTOMIC/GENOMIC/All（默认 TRANSCRIPTOMIC）')
    sp.add_argument('--db', default='both', choices=['sra', 'gsa', 'both'])
    sp.add_argument('--out', default=None, help='输出目录（默认 meta_search/<物种>/search）')
    sp.add_argument('--no-detailed', action='store_true', help='关闭详细模式（Tissue/Location 深提取）')
    sp.add_argument('--ncbi-api', default=None, help='NCBI API Key（也可用环境变量 NCBI_API_KEY）')
    sp.add_argument('--deepseek-api', default=None, help='DeepSeek API Key（AI 元数据清洗，可选）')
    sp.set_defaults(func=cmd_meta_search)

    sp = sub.add_parser('meta-info', help='公共数据检索：Run 列表 → Core14/Full 统一元数据')
    sp.add_argument('--runs', required=True, help='Run 编号列表文件（每行一个 SRR/CRR…）或单个编号')
    sp.add_argument('--out', default=None, help='输出目录（默认 meta_search/info/）')
    sp.add_argument('--mode', default='both', choices=['local', 'api', 'both'])
    sp.add_argument('-t', '--threads', type=int, default=4)
    sp.add_argument('--fill-date', action='store_true', help='采集日期缺失时用发布日期兜底')
    sp.add_argument('--ncbi-api', default=None)
    sp.add_argument('--deepseek-api', default=None)
    sp.set_defaults(func=cmd_meta_info)

    sp = sub.add_parser('meta-plot', help='公共数据检索：元数据 SCI 可视化（时间/机构/组织/地理）')
    sp.add_argument('--input', required=True, help='检索结果或 Core14/Full 元数据 CSV')
    sp.add_argument('--out', default='SCI_Figures_Output', help='图表输出目录')
    sp.set_defaults(func=cmd_meta_plot)

    sp = sub.add_parser('host-genome', help='宿主参考基因组下载（NCBI，供建宿主库）')
    sp.add_argument('--species', required=True, help='物种拉丁名')
    sp.add_argument('--out', default=None, help='输出目录（默认 host-db/<taxid>_<物种>/）')
    sp.add_argument('--include-organelles', action='store_true', help='同时下载叶绿体/线粒体基因组')
    sp.add_argument('--ncbi-api', default=None)
    sp.add_argument('--min-length', type=int, default=0)
    sp.add_argument('--verify-only', action='store_true', help='仅检查已有下载')
    sp.add_argument('--skip-datasets', action='store_true', help='跳过 datasets CLI，直接 E-utilities 回退')
    sp.set_defaults(func=cmd_host_genome)

    sp = sub.add_parser('ref-status', help='查看病毒参考库版本与统计')
    sp.set_defaults(func=cmd_ref_status)

    sp = sub.add_parser('ictv-update',
                        help='更新 ICTV VMR 参考库（在线下载或本地 xlsx → taxa.txt）')
    sp.add_argument('--xlsx', default=None,
                    help='本地 VMR xlsx 路径（不给则从 ictv.global 下载当前版）')
    sp.set_defaults(func=cmd_ictv_update)

    sp = sub.add_parser('ictv-status', help='查看 ICTV VMR 参考库状态')
    sp.set_defaults(func=cmd_ictv_status)

    sp = sub.add_parser('ictv-refs',
                        help='ICTV VMR 选参：按属/科/种挑参考，缺的按需下载')
    sp.add_argument('--genus', default=None, help='属名（如 Tobamovirus）')
    sp.add_argument('--family', default=None, help='科名（如 Potyviridae）')
    sp.add_argument('--species', default=None, help='种名（如 "Tobamovirus mosaic"）')
    sp.add_argument('--limit', type=int, default=20, help='最多挑多少条（默认 20）')
    sp.add_argument('--any', action='store_true',
                    help='不限定 Complete genome（含 Coding-complete 等）')
    sp.add_argument('--download', action='store_true',
                    help='对本地没有的 accession 执行 NCBI 按需下载')
    sp.set_defaults(func=cmd_ictv_refs)

    sp = sub.add_parser('tools', help='显示工具探测状态')
    sp.set_defaults(func=cmd_tools)

    sp = sub.add_parser('db-migrate',
                        help='数据库迁移/对接：软件与数据库分离部署')
    sp.add_argument('--to', required=True,
                    help='目标数据库根目录（绝对路径，平台目录之外）')
    sp.add_argument('--mode', default='copy', choices=['copy', 'move'],
                    help="copy=复制后切换（推荐，源保留）；move=复制校验后删源")
    sp.add_argument('--check', action='store_true',
                    help='只校验目标目录是否已有完整数据库（对接前检查）')
    sp.add_argument('--dry-run', action='store_true',
                    help='只做预检（空间/路径），不实际复制')
    sp.set_defaults(func=cmd_db_migrate)

    sp = sub.add_parser('tool-runs',
                        help='输出目录管理：status/organize/archive/clean')
    sp.add_argument('action', choices=['status', 'organize', 'archive', 'clean'],
                    help='status=现状统计；organize=杂项归位；archive=归档旧运行；'
                         'clean=删除超期运行/归档（不可逆，先 --dry-run）')
    sp.add_argument('--before', default=None,
                    help="archive: 归档该日期（YYYYMMDD）之前的运行")
    sp.add_argument('--older-days', type=int, default=None,
                    help='archive: 归档早于 N 天的运行')
    sp.add_argument('--active-days', type=int, default=None,
                    help='clean: 删除活动区早于 N 天的运行')
    sp.add_argument('--archive-days', type=int, default=None,
                    help='clean: 删除归档区早于 N 天的归档')
    sp.add_argument('--dry-run', action='store_true',
                    help='只预演，不实际移动/删除')
    sp.set_defaults(func=cmd_tool_runs)

    sp = sub.add_parser('selfcheck', help='环境自检（模块/工具/数据库/磁盘/依赖）')
    sp.set_defaults(func=cmd_selfcheck)

    args = p.parse_args()
    rc = args.func(args)
    # selfcheck 返回 False（存在未就绪项）时以非零码退出，供脚本判断
    if rc is False:
        sys.exit(1)


if __name__ == '__main__':
    main()

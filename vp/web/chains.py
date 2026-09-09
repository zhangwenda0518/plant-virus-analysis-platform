# -*- coding: utf-8 -*-
"""链式一键分析：把多个任务工厂串成一次运行（自 app.py 拆出）。

设计：每步仍调用 vp/web/tool_jobs.py 的现有 _tool_job_*，各建
自己的 run_dir，一步不改；汇总目录只放「符号链接 + 清单」，
把一次多步运行的结果收拢成一份。
"""
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid

from flask import abort, jsonify, render_template, request, send_file, \
    send_from_directory

from vp.config import DIRS, PLATFORM_ROOT, db_path, engine_cmd
from vp.utils import (TaskLogger, check_path, fmt_size, run_cmd, safe_open,
                      safe_remove)
from vp.web.state import cfg, tool_runs_root as _tool_runs_root
from vp.web.tool_jobs import (
    _tool_job_convert,
    _tool_job_fastp,
    _tool_job_hostremoval,
    _tool_job_hostpredict,
    _tool_job_orf,
    _tool_job_orfa,
    _tool_job_genoplot,
    _tool_job_primer,
    _tool_job_identify,
    _tool_job_assemble,
    _tool_job_contigs,
    _tool_job_structcmp,
    _tool_job_verify,
    _tool_job_consensus,
    _tool_job_kvsuite,
    _tool_job_quicktree,
    _tool_job_align,
    _tool_job_sdt,
    _tool_job_identity,
)


_CHAIN_LINKABLE = ('.fasta', '.fa', '.fna', '.tsv', '.json', '.txt',
                   '.vcf', '.gz', '.html', '.pdf', '.png')


def _chain_link(dst, src, logger=None, target_is_dir=None):
    """在汇总目录建链接指向 src。

    Windows 下 os.symlink 需开发者模式/管理员；不可用时按序降级：
    目录 → 硬链接不可用则跳过并记录；文件 → 硬链接；再不行则跳过。
    返回实际使用的方式（'symlink' / 'hardlink' / None），失败不抛异常。
    """
    if not os.path.exists(src):
        return None
    if target_is_dir is None:
        target_is_dir = os.path.isdir(src)
    try:
        os.symlink(src, dst, target_is_directory=bool(target_is_dir))
        return 'symlink'
    except (OSError, NotImplementedError, AttributeError):
        pass
    if not target_is_dir:
        try:
            os.link(src, dst)
            return 'hardlink'
        except OSError:
            pass
    if logger:
        logger.log(f'链接创建跳过（源或权限不可用）: {os.path.basename(dst)}')
    return None


def _chain_write_manifest(chain_dir, meta):
    p = os.path.join(chain_dir, '_chain.json')
    with safe_open(p, 'wt') as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    return p


def _chain_run_step(chain_dir, root, step_key, step_name, job_factory,
                    p, threads, db_virus, logger, prog, prog_lo, prog_hi,
                    cancel=None, collect=None, optional=()):
    """跑链中一步：新建 run_dir、构造 ctx、调用 job，返回 (run_name, result)。

    prog 为阶段名；进度在 [prog_lo, prog_hi] 区间内线性映射。
    collect: 该步要链接进汇总目录的产物文件名列表（可选）。
    optional: 这些 key 就算下游 job 用 ctx.req 取也允许为空（返回 None），
              用于兼容现有 job 把可选输入写成必填的情况。
    """
    from types import SimpleNamespace
    ts = time.strftime('%Y%m%d_%H%M%S')
    run_name = f'{step_key}_{ts}'
    run_dir = check_path(os.path.join(root, run_name), must_exist=False,
                         in_platform=True)
    os.makedirs(run_dir, exist_ok=True)

    def _input_abs(v):
        return check_path(v if os.path.isabs(v)
                          else os.path.join(PLATFORM_ROOT, v),
                          must_exist=True)

    def _req(key, what):
        v = (p.get(key) or '').strip()
        if not v:
            if key in optional:
                return None
            abort(400, f'请选择{what}')
        return _input_abs(v)

    def _opt(key):
        v = (p.get(key) or '').strip()
        return _input_abs(v) if v else None

    ctx = SimpleNamespace(p=p, run_dir=run_dir, threads=threads,
                          db_virus=db_virus, req=_req, opt=_opt)
    job = job_factory(ctx)

    def _sub_prog(_stage, frac, msg):
        f = prog_lo + max(0.0, min(1.0, float(frac))) * (prog_hi - prog_lo)
        prog(step_key, f, msg)

    prog(step_key, prog_lo, f'{step_name}：开始')
    res = job(lambda m: logger.log(m), _sub_prog, cancel)
    logger.log(f'{step_name} 完成: run={run_name}')

    # 目录链接（整步产物）
    _chain_link(os.path.join(chain_dir, f'step_{step_key}'), run_dir,
                logger=logger, target_is_dir=True)
    # 关键文件链接
    linked = []
    for fn in (collect or ()):
        src = os.path.join(run_dir, fn)
        if os.path.isfile(src):
            how = _chain_link(os.path.join(chain_dir, f'{step_key}__{fn}'),
                              src, logger=logger, target_is_dir=False)
            if how:
                linked.append(fn)
    return run_name, res, linked


def _tool_job_virchain(ctx):
    """病毒识别和分类分析 · 一键：②提取序列组装 → ③组装结果再鉴定 → ④候选序列验证。

    输入固定为去宿主序列（FASTA，可单/双端）。①（识别分类与提取）为可选：
      跳过 ①：去宿主序列 ────────────────────→ ②组装 → ③再鉴定 → ④验证
      运行 ①：去宿主序列 → ①提取病毒序列 → ②组装 → ③再鉴定 → ④验证
    即两条路都到 ②，区别只在 ② 的输入（原始去宿主序列 vs ① 提取的病毒序列）。
    组装（SPAdes）接受 FASTA，或自动 --only-assembler 跳过纠错。
    每步各建 run_dir，汇总目录放符号链接 + _chain.json。
    """
    root = _tool_runs_root()
    inp = ctx.req('input', '去宿主序列（FASTA）')
    inp2 = ctx.opt('input2')
    # 可选：先跑 ① 提取病毒序列，作为 ② 的输入
    run_identify = bool(ctx.p.get('run_identify'))
    # ① 的输入类型：双端 / 单端 / FASTA
    if inp2:
        id_type = 'pe'
    else:
        id_type = ('fasta' if inp.lower().endswith(('.fa', '.fasta', '.fna'))
                   else 'single')

    do_verify = ctx.p.get('do_verify')
    do_verify = True if do_verify is None else bool(do_verify)
    asm_mode = (ctx.p.get('mode') or 'metaviral').strip()
    mem = int(ctx.p.get('memory') or 64)
    min_len = int(ctx.p.get('min_len') or 200)
    conf = float(ctx.p.get('confidence') or 0)
    host = (ctx.p.get('host') or 'all').strip() or 'all'
    combine = (ctx.p.get('combine') or 'union').strip()
    methods = ctx.p.get('methods') or ['blastx', 'cdd']
    if isinstance(methods, str):
        methods = [m.strip() for m in methods.split(',') if m.strip()]
    methods = [m for m in methods if m in ('blastx', 'cdd')]

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        chain_dir = ctx.run_dir
        chain_name = os.path.basename(chain_dir)
        meta = {'chain': 'virchain', 'run': chain_name,
                'label': '病毒识别和分类分析 一键流程',
                'run_identify': run_identify,
                'created': time.strftime('%Y-%m-%d %H:%M:%S'), 'steps': []}

        # 进度权重：identify(若开) 0.20 / assembly 0.35 / contigs 0.25 / verify 0.20
        w_id = 0.20 if run_identify else 0.0
        w_asm = 0.35
        w_ctg = 0.25
        w_vfy = 0.20 if do_verify else 0.0
        tot = w_id + w_asm + w_ctg + w_vfy

        def _band(lo, w):
            return (lo, lo + (w / tot if tot else 0.0))

        cur = 0.0

        # ---- ①（可选）病毒识别分类与提取：产物作为 ② 的输入 ----
        asm_in1, asm_in2 = inp, inp2
        if run_identify:
            lo, hi = _band(cur, w_id)
            p1 = {'input': inp, 'input_type': id_type, 'confidence': conf}
            if inp2:
                p1['input2'] = inp2
            run1, res1, _ = _chain_run_step(
                chain_dir, root, 'identify', '① 病毒识别分类与提取',
                _tool_job_identify, p1, ctx.threads, ctx.db_virus,
                logger, prog, lo, hi, cancel=cancel,
                collect=('viral_ids.tsv',))
            meta['steps'].append({'key': 'identify', 'run': run1,
                                  'title': '① 病毒识别分类与提取',
                                  'summary': {k: v for k, v in (res1 or {}).items()
                                              if not isinstance(v, (list, dict))}})
            # 取 ① 实际产出的病毒序列作为 ② 的输入（配对输入会带 _1/_2 后缀）
            id_dir = os.path.join(root, run1)
            vseqs = []
            if os.path.isdir(id_dir):
                for fn in sorted(os.listdir(id_dir)):
                    if fn.startswith('viral_sequences') and fn.endswith('.fasta'):
                        fp = os.path.join(id_dir, fn)
                        if os.path.getsize(fp) > 0:
                            vseqs.append(fp)
            if not vseqs:
                raise RuntimeError('① 未提取到病毒序列（viral_sequences*.fasta 为空），'
                                   '可取消「先跑 ①」开关直接用去宿主序列组装')
            # 链接进汇总目录
            for fp in vseqs:
                _chain_link(os.path.join(chain_dir,
                                         'identify__' + os.path.basename(fp)),
                            fp, logger=logger, target_is_dir=False)
            asm_in1 = vseqs[0]
            asm_in2 = vseqs[1] if len(vseqs) > 1 else None
            meta['steps'][-1]['output'] = [os.path.basename(x) for x in vseqs]
            logger.log('② 组装输入改用 ① 提取的病毒序列: '
                       + ', '.join(os.path.basename(x) for x in vseqs))
            cur = hi

        # ---- ② 提取序列组装（输入 = 去宿主序列 或 ① 的提取序列） ----
        lo, hi = _band(cur, w_asm)
        p2 = {'r1': asm_in1, 'mode': asm_mode, 'memory': mem,
              'min_len': min_len}
        if asm_in2:
            p2['r2'] = asm_in2
        run2, res2, _ = _chain_run_step(
            chain_dir, root, 'assembly', '② 提取序列组装',
            _tool_job_assemble, p2, ctx.threads, ctx.db_virus,
            logger, prog, lo, hi, cancel=cancel,
            collect=('contigs.filtered.fasta',), optional=('r2',))
        contigs = (res2 or {}).get('contigs')
        meta['steps'].append({'key': 'assembly', 'run': run2,
                              'title': '② 提取序列组装',
                              'summary': {k: v for k, v in (res2 or {}).items()
                                          if not isinstance(v, (list, dict))}})
        if not contigs or not os.path.isfile(contigs):
            raise RuntimeError('② 组装未产出 contigs，后续步骤无法继续')
        cur = hi

        # ---- ③ 组装结果再鉴定 ----
        lo, hi = _band(cur, w_ctg)
        p3 = {'contigs': contigs, 'min_len': min_len, 'confidence': conf}
        run3, res3, _ = _chain_run_step(
            chain_dir, root, 'contigs', '③ 组装结果再鉴定',
            _tool_job_contigs, p3, ctx.threads, ctx.db_virus,
            logger, prog, lo, hi, cancel=cancel,
            collect=('virus_classification.tsv', 'viral_contigs.fasta'))
        viral_fa = (res3 or {}).get('viral_fasta')
        meta['steps'].append({'key': 'contigs', 'run': run3,
                              'title': '③ 组装结果再鉴定',
                              'summary': {k: v for k, v in (res3 or {}).items()
                                          if not isinstance(v, (list, dict))}})
        cur = hi

        # ---- ④ 候选序列验证（依赖 ③ 的产物，直接读其 run_dir） ----
        if do_verify:
            if not viral_fa or not os.path.isfile(viral_fa):
                logger.log('③ 未提取到病毒 contig，跳过 ④ 候选序列验证')
                meta['steps'].append({'key': 'verify', 'run': None,
                                      'title': '④ 候选序列验证',
                                      'skipped': '③ 未产出 viral_contigs.fasta'})
            else:
                lo, hi = _band(cur, w_vfy)
                # verify 的 run 模式需要源 run 名；直接传 ③ 的 run
                p4 = {'run': run3, 'host': host, 'combine': combine,
                      'methods': methods}
                run4, res4, _ = _chain_run_step(
                    chain_dir, root, 'verify', '④ 候选序列验证',
                    _tool_job_verify, p4, ctx.threads, ctx.db_virus,
                    logger, prog, lo, hi, cancel=cancel,
                    collect=('verify_summary.json',))
                meta['steps'].append({'key': 'verify', 'run': run4,
                                      'title': '④ 候选序列验证',
                                      'summary': {k: v for k, v in (res4 or {}).items()
                                                  if not isinstance(v, (list, dict))}})
                cur = hi

        prog('chain', 1.0, '一键流程完成')
        logger.close()
        meta['done'] = True
        _chain_write_manifest(chain_dir, meta)
        return {'chain': chain_name, 'run': chain_name,
                'n_steps': len(meta['steps']), 'steps': meta['steps'],
                'chain_json': os.path.join(chain_dir, '_chain.json')}
    return job


def _tool_job_kvchain(ctx):
    """病毒定量与共识 · 一键：调用 kvsuite 全流程（鉴定→过滤→共识→绘图→变异）。

    kvsuite 卡里的 stage 下拉默认已是 all，本 job 把「一次跑完并收拢产物」
    做成显式入口：建汇总目录 + 链接 kvsuite run_dir + 关键产物。
    """
    root = _tool_runs_root()

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        chain_dir = ctx.run_dir
        chain_name = os.path.basename(chain_dir)
        # kvsuite 的 stage 固定为 all（一键语义）
        p = dict(ctx.p)
        p['stage'] = 'all'
        meta = {'chain': 'kvchain', 'run': chain_name,
                'label': '病毒定量与共识 一键流程',
                'created': time.strftime('%Y-%m-%d %H:%M:%S'), 'steps': []}

        run1, res1, _ = _chain_run_step(
            chain_dir, root, 'kvsuite', '已知病毒识别与定量（全流程）',
            _tool_job_kvsuite, p, ctx.threads, ctx.db_virus,
            logger, prog, 0.0, 1.0, cancel=cancel,
            collect=('run.log',))
        meta['steps'].append({'key': 'kvsuite', 'run': run1,
                              'title': '已知病毒识别与定量（全流程）',
                              'summary': {k: v for k, v in (res1 or {}).items()
                                          if not isinstance(v, (list, dict))}})
        prog('chain', 1.0, '一键流程完成')
        logger.close()
        meta['done'] = True
        _chain_write_manifest(chain_dir, meta)
        return {'chain': chain_name, 'run': chain_name,
                'n_steps': 1, 'steps': meta['steps'],
                'kvsuite_run': run1,
                'chain_json': os.path.join(chain_dir, '_chain.json')}
    return job


def _tool_job_chainoverview(ctx):
    """占位：链式汇总目录的读取走 /api/tool/chain_result，不注册为工具。"""
    raise NotImplementedError

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_consensus.py — 共识段
=======================
完全照搬 virome_analysis_pipeline/batch_virus_variants.py 的共识实现，
唯一改动：把产 BAM 的比对引擎从 bowtie2 换成 minibwa。

原管线 B 版链路（batch_virus_variants.py）：
  1. extract_virus_fastas (line 700+)
       从总参考库按 '>' 头切分，每病毒写单序列 ref_{virus}.ref.fasta
       + samtools faidx
  2. worker_align (line 278)
       bowtie2 --local -p T {-f|-q} -x index {-1/-2|-U} | samtools sort -o out.bam
       索引建在全参考 harmonized fasta 上 → 一个 BAM 含多 @SQ
  3. worker_consensus (line 395)
       samtools view -h BAM virus
         | awk '/^@SQ/ && $2 != "SN:"v {next} {print}'
         | samtools view -b | samtools sort -o fixed.bam      # 瘦身成单参考
       samtools index fixed.bam
       viral_consensus -i fixed.bam -r ref.fa -o out.fa -q 20 -d 5 -f 0.5 -a N
       finally 删 fixed.bam 与 .bai
  4. 深度动态计算 (line 1171)
       d = max(1, min(10, floor(Recalc_MeanDepth/2)))  若 vc_depth == 0
       否则用 vc_depth（默认 5）

参数默认值（batch_virus_variants.py line 1240-1243）：
  -q/--vc_qual  20
  -d/--vc_depth 5
  -f/--vc_freq  0.5
  -a/--vc_ambig N
"""

import math
import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


# ══════════════════════════════════════════════════════════
# 工具函数（照抄 batch_virus_variants.py line 82）
# ══════════════════════════════════════════════════════════
def safe_name(s, max_len=100):
    s = str(s)
    s = re.sub(r'[^A-Za-z0-9\-.]', '_', s)
    s = re.sub(r'_+', '_', s)
    return s.strip('_.')[:max_len]


# # ══════════════════════════════════════════════════════════
# 外部命令执行（列表传参，不经 shell）
# ══════════════════════════════════════════════════════════
def _run(cmd, log_path=None, master_log=None, check=True, logger=None, env=None,
         stdin_data=None):
    """
    执行外部命令并记录到日志。

    ★ 与原管线的差异（Windows 移植必需）：
    原管线（batch_virus_variants.py run_cmd）用 shell=True + executable='/bin/bash'
    + 'set -o pipefail'，命令串里用单引号包路径。
    本机 Windows 无 bash，cmd.exe 不把单引号当引号字符，实测：
        单引号 -> exit=1（"The filename, directory name, or volume label syntax is incorrect."）
        双引号 -> exit=0；无引号 -> exit=0；列表 -> exit=0
    因此这里统一改为列表传参（shell=False），由 subprocess 直接处理路径含空格、
    含中文的情况，彻底避开引号解释层。原管线的 awk 过滤器改由
    _filter_sam_header() 用 Python 实现，管道改为 Python 桥接。
    """
    if isinstance(cmd, str):
        import shlex
        cmd_list = shlex.split(cmd, posix=False)
        cmd_list = [c[1:-1] if len(c) >= 2 and c[0] == c[-1] == '"' else c
                    for c in cmd_list]
    else:
        cmd_list = [str(c) for c in cmd]

    start_time = time.strftime('%Y-%m-%d %H:%M:%S')
    result = subprocess.run(cmd_list, shell=False, capture_output=True,
                            text=True, encoding='utf-8', errors='replace',
                            env=env, input=stdin_data)

    disp = ' '.join(f'"{c}"' if ' ' in str(c) else str(c) for c in cmd_list)
    log_content = f"\n[{start_time}] CMD: {disp}\nEXIT_CODE: {result.returncode}\n"
    if result.stdout:
        log_content += f"--- STDOUT ---\n{result.stdout.strip()}\n"
    if result.stderr:
        log_content += f"--- STDERR ---\n{result.stderr.strip()}\n"
    log_content += "-" * 80 + "\n"

    for lp in (log_path, master_log):
        if lp:
            lp = Path(lp)
            lp.parent.mkdir(parents=True, exist_ok=True)
            with open(lp, 'a', encoding='utf-8', errors='replace') as f:
                f.write(log_content)

    if logger and check:
        logger.debug(f"[CMD] {disp}")

    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, cmd_list, output=result.stdout,
            stderr=result.stderr)
    return result


def _write_cmd_log(log_path, master_log, disp, returncode, out=None, err=None):
    """管道式命令的日志记录（_run 的伴侣，供 Python 桥接管道使用）"""
    start_time = time.strftime('%Y-%m-%d %H:%M:%S')
    content = f"\n[{start_time}] CMD: {disp}\nEXIT_CODE: {returncode}\n"
    for tag, blob in (('STDOUT', out), ('STDERR', err)):
        if blob:
            txt = blob.decode('utf-8', 'replace') if isinstance(blob, bytes) else str(blob)
            if txt.strip():
                content += f"--- {tag} ---\n{txt.strip()}\n"
    content += "-" * 80 + "\n"
    for lp in (log_path, master_log):
        if lp:
            lp = Path(lp)
            lp.parent.mkdir(parents=True, exist_ok=True)
            with open(lp, 'a', encoding='utf-8', errors='replace') as f:
                f.write(content)
# ══════════════════════════════════════════════════════════
def build_runtime_env(tools, extra_paths=None):
    """
    构造子进程环境。

    ★ 运行时硬约束（已实测）：
    viral_consensus.exe 依赖一组 htslib DLL，而 samtools 自带同名不同版本的 DLL。
    若 PATH 中 samtools/bin 排在 MinGW bin 之前，viral_consensus 会以
    3221225477 (0xC0000005) 崩溃；若 MinGW bin 完全不在 PATH，则报
    3221225785 (0xC0000139，DLL 入口点找不到)。

    实测三种 PATH 组合：
        仅默认 PATH              -> 3221225785
        MinGW > viral_consensus > samtools -> 正常（exit 1，仅因测试文件已被清理）
        viral_consensus > samtools > MinGW -> 3221225477

    因此这里固定把 MinGW bin 插到最前，再是 viral_consensus 目录，再是 samtools。
    """
    env = dict(os.environ)
    dirs = []
    if extra_paths:
        dirs.extend(str(p) for p in extra_paths)

    # viral_consensus 与 samtools 所在目录
    for name in ('viral_consensus', 'samtools'):
        p = None
        if tools is not None:
            if hasattr(tools, 'paths'):
                p = tools.paths.get(name)
            elif isinstance(tools, dict):
                p = tools.get(name)
        if p:
            d = str(Path(str(p)).parent)
            # MinGW bin 必须最前；其余按名字次序追加
            dirs.append(d)

    # MinGW bin（viral_consensus 的同源 C++ 运行库）
    mingw_dir = None
    for exe in ('g++.exe', 'gcc.exe', 'g++', 'gcc'):
        found = shutil.which(exe)
        if found:
            mingw_dir = str(Path(found).parent)
            break
    if mingw_dir:
        dirs.insert(0, mingw_dir)

    seen, ordered = set(), []
    for d in dirs:
        if d and d not in seen and Path(d).is_dir():
            seen.add(d)
            ordered.append(d)
    if ordered:
        env['PATH'] = os.pathsep.join(ordered) + os.pathsep + env.get('PATH', '')
    return env


def _filter_sam_header(sam_text, virus):
    """
    完全等价于原管线里的 awk 过滤器：

        awk -v v='{virus}' '/^@SQ/ && $2 != "SN:"v {next} {print}'

    语义：逐行读，@SQ 行若第二列不等于 "SN:<virus>" 则丢弃，
          其余所有行（包括其他 @ 行和比对记录）原样保留。

    多参考 BAM 必须瘦身成单参考，否则 viral_consensus 会报
    "CRAB/BAM/SAM has N references, but it should have exactly 1" 并 exit(1)。
    """
    target_sn = 'SN:' + str(virus)
    out = []
    n_drop = 0
    for ln in sam_text.splitlines():
        if ln.startswith('@SQ'):
            parts = ln.split('\t')
            if len(parts) >= 2 and parts[1] != target_sn:
                n_drop += 1
                continue                      # awk 的 {next}
        out.append(ln)                        # awk 的 {print}
    return '\n'.join(out) + ('\n' if out else ''), n_drop


# ══════════════════════════════════════════════════════════
# 参考序列提取（照抄 extract_virus_fastas, line 700+）
# ══════════════════════════════════════════════════════════
def extract_virus_fastas(reference, target_set, out_dir, samtools, log_file=None,
                         master_log=None, logger=None, env=None):
    """
    从总参考库按 '>' 头切分，为每个目标病毒写单序列 ref_{name}.ref.fasta
    并 samtools faidx。返回 {virus: Path(ref_fa)}

    完全照抄 batch_virus_variants.py extract_virus_fastas。
    """
    target_set = set(target_set)
    found_map = {}
    seq_buf = []
    vid_cur = None
    d_fasta = Path(out_dir) / 'virus-fasta'
    d_fasta.mkdir(parents=True, exist_ok=True)

    def _flush():
        if vid_cur not in target_set or not seq_buf:
            return
        folder = f"ref_{safe_name(vid_cur)}"
        vdir = d_fasta / folder
        vdir.mkdir(parents=True, exist_ok=True)

        ref_fa = vdir / f"{folder}.ref.fasta"
        if not ref_fa.exists():
            with open(ref_fa, 'w', encoding='utf-8') as f:
                f.write(f">{vid_cur}\n" + "".join(seq_buf) + "\n")
        if not Path(str(ref_fa) + ".fai").exists():
            _run([samtools, 'faidx', str(ref_fa)], log_path=log_file,
                 master_log=master_log, check=False, logger=logger, env=env)

        found_map[vid_cur] = ref_fa

    with open(reference, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                _flush()
                vid_cur = line[1:].split()[0]
                seq_buf = []
            else:
                seq_buf.append(line)
    _flush()

    if logger:
        logger.info(f"  从库中提取 {len(found_map)}/{len(target_set)} 条代表株序列 -> {d_fasta}")
    return found_map


# ══════════════════════════════════════════════════════════
# 比对：minibwa 替代 bowtie2（唯一改动点）
# ══════════════════════════════════════════════════════════
def build_minibwa_index(ref_fasta, index_prefix, threads, log_file=None,
                        master_log=None, logger=None, env=None,
                        minibwa='minibwa'):
    """
    对应原管线 resolve_bam_map 里的 bowtie2-build：
        bowtie2-build --threads T '{harmonized_fasta}' '{index_prefix}'
    改为 minibwa index（索引建在全参考 harmonized fasta 上）。
    """
    if Path(str(index_prefix) + '.mbw').exists():
        if logger:
            logger.info(f"  复用已有 minibwa 索引 -> {index_prefix}")
        return

    Path(index_prefix).parent.mkdir(parents=True, exist_ok=True)
    _run([minibwa, 'index', str(ref_fasta), str(index_prefix)],
         log_path=log_file, master_log=master_log, check=True,
         logger=logger, env=env)
    if logger:
        logger.info(f"  minibwa 索引完成 -> {index_prefix}")


def align_sample(fq, index_prefix, out_bam, threads, log_file=None,
                 master_log=None, logger=None, env=None,
                 minibwa='minibwa', samtools='samtools'):
    """
    对应原管线 worker_align (line 278)：
        bowtie2 --local -p {threads} {format_flag} -x '{index_prefix}' {fq_arg}
          | samtools sort -@ {t_io} -o '{out_bam}'
    改为 minibwa map，其余（管道 samtools sort）保持不变。

    bowtie2 --local 语义 → minibwa 默认模式（-x adap，短/长读自适应）。
    原管线的 -f/-q 与 -1/-2/-U 在 minibwa 里由输入文件后缀与位置参数决定，
    minibwa map 接受 <idx> <in1> [in2]，与 bowtie2 的 -1/-2/-U 等价。
    """
    t_io = min(4, threads)
    # 注意：不加 minibwa 的 -y（copy FASTA/Q comments to output）。
    # 实测：-y 会把 FASTQ header 第二段（原始 read name）作为第 12 个字段
    # 追加到 SAM 行末尾，samtools 严格解析报 aux_parse 失败；污染只出现在
    # unmapped reads 上，而此处未加 -u（保持原管线输出 unmapped 的行为），
    # 一旦有 unmapped read 就会让整条 sort 管道失败。mapped reads 的 SAM
    # 内容与不带 -y 时完全一致，因此直接移除 -y。
    map_cmd = [minibwa, 'map', f'-t{threads}', '-o', '-',
               str(index_prefix), fq['r1']]
    if fq.get('r2'):
        map_cmd.append(fq['r2'])

    # 原管线的 `minibwa map ... | samtools sort ...`，改用 Python 桥接：
    # 逐块把 SAM 文本吸给 samtools sort 的 stdin，语义一致且不经 shell。
    p_map = subprocess.Popen(map_cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    p_sort = subprocess.Popen([samtools, 'sort', '-@', str(t_io),
                               '-o', str(out_bam), '-'], stdin=p_map.stdout,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p_map.stdout.close()                 # 交给 p_sort 后本地关闭写端
    err_sort = p_sort.communicate()[1]
    err_map = p_map.stderr.read()
    p_map.stderr.close()
    rc_map = p_map.wait()
    rc_sort = p_sort.returncode

    _write_cmd_log(log_file, master_log,
                   f"{' '.join(map_cmd)} | {samtools} sort -@ {t_io} -o {out_bam}",
                   max(rc_map, rc_sort),
                   err=(err_map or b'') + (err_sort or b''))
    if rc_map != 0 or rc_sort != 0:
        raise subprocess.CalledProcessError(
            max(rc_map, rc_sort), 'minibwa map | samtools sort',
            stderr=(err_map or b'') + (err_sort or b''))

    bai = Path(str(out_bam) + '.bai')
    if out_bam and Path(out_bam).exists() and not bai.exists():
        _run([samtools, 'index', '-@', str(t_io), str(out_bam)],
             log_path=log_file, master_log=master_log, check=False,
             logger=logger, env=env)


# ══════════════════════════════════════════════════════════
# 共识 worker（照抄 worker_consensus, line 395）
# ══════════════════════════════════════════════════════════
def worker_consensus(args):
    (sample, virus, bam_path, ref_fa, out_fa_str, fixed_bam_str,
     depth, qual, freq, ambig, threads, resume, log_file, master_log,
     minibwa, samtools, viral_consensus, env) = args

    out_fa = Path(out_fa_str)
    fixed_bam = Path(fixed_bam_str)

    if resume and out_fa.exists() and out_fa.stat().st_size > 0:
        return True

    try:
        t_io = min(4, threads)

        # ── 原管线管道（batch_virus_variants.py:395）──
        #   samtools view -h BAM virus
        #     | awk -v v=virus '/^@SQ/ && $2 != "SN:"v {next} {print}'
        #     | samtools view -b | samtools sort -o fixed.bam
        # 这里用 Python 桥接：samtools 取 SAM -> _filter_sam_header 做 awk 的活
        # -> samtools view -b 存 BAM。不用 shell，不需要 awk/sh。
        p_hdr = subprocess.Popen(
            [samtools, 'view', '-@', str(t_io), '-h', str(bam_path), str(virus)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sam_out, sam_err = p_hdr.communicate()
        if p_hdr.returncode != 0:
            raise subprocess.CalledProcessError(
                p_hdr.returncode, 'samtools view -h', stderr=sam_err)

        sam_text = sam_out.decode('utf-8', 'replace')
        filtered, n_drop = _filter_sam_header(sam_text, virus)

        p_bam = subprocess.Popen(
            [samtools, 'view', '-@', str(t_io), '-b', '-o', str(fixed_bam), '-'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, bam_err = p_bam.communicate(input=filtered.encode('utf-8'))
        _write_cmd_log(log_file, master_log,
                       f"samtools view -h {bam_path} {virus} | "
                       f"[py] 瘦身 @SQ -{n_drop} | samtools view -b -o {fixed_bam}",
                       p_bam.returncode, err=bam_err)
        if p_bam.returncode != 0:
            raise subprocess.CalledProcessError(
                p_bam.returncode, 'samtools view -b', stderr=bam_err)

        # 原管道末尾的 samtools sort（对齐原管线输出，保坐标序）
        sorted_bam = Path(str(fixed_bam) + '.sorted')
        _run([samtools, 'sort', '-@', str(t_io), '-o', str(sorted_bam),
              str(fixed_bam)], log_path=log_file, master_log=master_log, env=env)
        if sorted_bam.exists():
            Path(sorted_bam).replace(fixed_bam)

        _run([samtools, 'index', '-@', str(t_io), str(fixed_bam)],
             log_path=log_file, master_log=master_log, env=env)

        vc_cmd = [viral_consensus, '-i', str(fixed_bam), '-r', str(ref_fa),
                  '-o', str(out_fa), '-q', str(qual), '-d', str(depth),
                  '-f', str(freq), '-a', str(ambig)]
        _run(vc_cmd, log_path=log_file, master_log=master_log, env=env)
        return True
    except Exception:
        import traceback
        err = f"\n[Exception] worker_consensus: {traceback.format_exc()}\n"
        with open(log_file, 'a', encoding='utf-8') as lf:
            lf.write(err)
        return False
    finally:
        for f in (fixed_bam, Path(str(fixed_bam) + '.bai')):
            if f.exists():
                f.unlink()


# ══════════════════════════════════════════════════════════
# 蛋白 QC（保留原有实现，不属于共识算法本身）
# ══════════════════════════════════════════════════════════
def check_protein_integrity(aa_seq):
    """去掉尾部终止符后不应再有内部终止符"""
    clean = str(aa_seq).rstrip('*')
    return clean.count('*') == 0, clean.count('*')


def translate(seq, table=1):
    from Bio.Seq import Seq as _Seq
    return str(_Seq(seq).translate(table=table, cds=False))


def _load_genbank_features(gb_path):
    from Bio import SeqIO
    recs = {}
    for rec in SeqIO.parse(str(gb_path), 'genbank'):
        feats = []
        for ft in rec.features:
            if ft.type not in ('CDS', 'mat_peptide'):
                continue
            name = ft.qualifiers.get('product', [ft.type])[0]
            gene = ft.qualifiers.get('gene', ['NA'])[0]
            feats.append((int(ft.location.start), int(ft.location.end), name, gene))
        recs[rec.id] = feats
    return recs


def _read_first_sequence(fasta_file):
    """照抄原管线 read_first_sequence"""
    seq_name, seq_data = '', []
    with open(fasta_file, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if seq_name:
                    break
                seq_name = line[1:].split()[0]
            else:
                seq_data.append(line)
    return seq_name, ''.join(seq_data).upper()


# ══════════════════════════════════════════════════════════
# 共识段主流程
# ══════════════════════════════════════════════════════════
class ConsensusStage:
    def __init__(self, args, tools, logger, out_dir, engine=None, env=None):
        self.args = args
        self.tools = tools
        self.logger = logger
        self.out_dir = Path(out_dir)
        self.engine = engine
        # env 未指定时按运行时硬约束自动构造（MinGW bin 必须最前）
        self.env = env if env is not None else build_runtime_env(tools)
        self.d_consensus = self.out_dir / 'consensus'
        self.d_logs = self.out_dir / 'logs'
        self.d_bam = self.out_dir / 'bam'
        self.master_log = str(self.d_logs / 'consensus_master.log')
        for d in (self.d_consensus, self.d_logs, self.d_bam):
            d.mkdir(parents=True, exist_ok=True)

    def _tool(self, name, fallback):
        """tools 可能是 ToolRegistry 对象或普通 dict，两种都支持"""
        if self.tools is None:
            return fallback
        if hasattr(self.tools, 'paths'):
            p = self.tools.paths.get(name)
        elif isinstance(self.tools, dict):
            p = self.tools.get(name)
        else:
            p = None
        return str(p) if p else fallback

    # ── 主流程 ───────────────────────────────────────────
    def run(self, samples, passed_df, ref_info=None):
        """
        samples:   样本列表 [{'name':.., 'r1':.., 'r2':..}]
        passed_df: 过滤段保留的表，含 Sample / Accession（或 Virus）
        ref_info:  可选，用于取 Recalc_MeanDepth 动态深度
        """
        a = self.args
        logger = self.logger
        logger.info('=' * 60)
        logger.info('【共识段】照搬 batch_virus_variants.py（minibwa 替代 bowtie2）')
        logger.info('=' * 60)

        if passed_df is None or not len(passed_df):
            logger.warning('过滤结果为空，跳共识段')
            return {}

        # 列名契约：原管线用 Sample / Virus；本模块过滤产物用 Sample / Accession
        cols = passed_df.columns
        sp_col = 'Sample' if 'Sample' in cols else cols[0]
        vc_col = 'Virus' if 'Virus' in cols else ('Accession' if 'Accession' in cols else cols[1])
        tax_col = 'taxonomy' if 'taxonomy' in cols else None
        depth_col = 'Recalc_MeanDepth' if 'Recalc_MeanDepth' in cols else None

        # ── 1. 提取每病毒单序列参考 ──
        target_set = set()
        rows = []
        for r in passed_df.iter_rows(named=True):
            s = r.get(sp_col)
            v = r.get(vc_col)
            if not s or not v:
                continue
            target_set.add(v)
            rows.append(r)

        if not rows:
            logger.warning('过滤表无有效行，跳共识段')
            return {}

        ref_map = extract_virus_fastas(
            a.reference, target_set, self.out_dir,
            samtools=self._tool('samtools', 'samtools'),
            log_file=str(self.d_logs / 'extract_fasta.log'),
            master_log=self.master_log, logger=logger, env=self.env)

        # ── 2. 建索引 + 比对（minibwa 替代 bowtie2）──
        index_prefix = str(self.out_dir / 'index' / 'harmonized_minibwa')
        harmonized = self.out_dir / 'index' / 'harmonized.fasta'
        if not Path(str(index_prefix) + '.mbw').exists():
            harmonized.parent.mkdir(parents=True, exist_ok=True)
            with open(harmonized, 'w', encoding='utf-8') as out_f:
                for fp in ref_map.values():
                    with open(fp, 'r', encoding='utf-8', errors='replace') as in_f:
                        out_f.write(in_f.read() + "\n")
            build_minibwa_index(
                str(harmonized), index_prefix, a.threads,
                log_file=str(self.d_logs / 'minibwa_index.log'),
                master_log=self.master_log, logger=logger, env=self.env,
                minibwa=self._tool('minibwa', 'minibwa'))

        bam_map = {}
        for s in samples:
            sname = s['name']
            out_bam = self.d_bam / f"{sname}.sorted.bam"
            if not (a.resume and out_bam.exists() and out_bam.stat().st_size > 0):
                align_sample(
                    s, index_prefix, str(out_bam), a.threads,
                    log_file=str(self.d_logs / f"align_{sname}.log"),
                    master_log=self.master_log, logger=logger, env=self.env,
                    minibwa=self._tool('minibwa', 'minibwa'),
                    samtools=self._tool('samtools', 'samtools'))
            if out_bam.exists() and out_bam.stat().st_size > 0:
                bam_map[sname] = str(out_bam)
            else:
                logger.warning(f"  [{sname}] 比对未产出 BAM，跳过")

        # ── 3. 构造任务（照抄 line 1160-1190）──
        tasks = []
        for r in rows:
            s = r.get(sp_col)
            v = r.get(vc_col)
            if s not in bam_map or v not in ref_map:
                continue
            t = r.get(tax_col, 'Unannotated') if tax_col else 'Unannotated'

            depth_fallback = float(r.get(depth_col, 0.0) or 0.0) if depth_col else 0.0
            d = (max(1, min(10, int(math.floor(depth_fallback / 2))))
                 if a.vc_depth == 0 else a.vc_depth)

            L1 = f"{safe_name(t)}_{safe_name(v)}"
            L2 = f"{safe_name(s)}_{safe_name(v)}"
            vdir = self.d_consensus / L1 / L2
            vdir.mkdir(parents=True, exist_ok=True)
            log_file = str(self.d_logs / f"{L2}_consensus.log")

            tasks.append((
                s, v, bam_map[s], str(ref_map[v]),
                str(vdir / f"{L2}.consensus.fasta"),
                str(vdir / f"{L2}.fixed.bam"),
                d, a.vc_qual, a.vc_freq, a.vc_ambig,
                a.threads, a.resume, log_file, self.master_log,
                self._tool('minibwa', 'minibwa'),
                self._tool('samtools', 'samtools'),
                self._tool('viral_consensus', 'viral_consensus'),
                self.env,
            ))

        if not tasks:
            logger.warning('无可用共识任务')
            return {}

        logger.info(f"  共识任务 {len(tasks)} 个（Jobs:{a.jobs}, Threads/Job:{min(4, a.threads)}）")

        # ── 4. 并行执行（照抄 ProcessPoolExecutor）──
        ok = 0
        with ProcessPoolExecutor(max_workers=min(len(tasks), a.jobs)) as ex:
            futs = {ex.submit(worker_consensus, t): t for t in tasks}
            for fut in as_completed(futs):
                try:
                    if fut.result():
                        ok += 1
                except Exception as e:
                    logger.error(f"  共识任务异常: {e}")
        logger.info(f"  共识完成 {ok}/{len(tasks)}")

        # ── 5. 蛋白 QC（保留原模块行为）──
        results = {}
        for t in tasks:
            s, v = t[0], t[1]
            fa = Path(t[4])
            if not (fa.exists() and fa.stat().st_size > 0):
                logger.warning(f"  [{s}] {v} 共识未产出，不计入结果")
                continue
            results.setdefault(s, {})[v] = fa
        self._run_protein_qc(results)
        return results

    def _run_protein_qc(self, results):
        """对每个样本的共识做粗略蛋白完整性检查，输出 qc.tsv"""
        rows = []
        for sname, virs in results.items():
            for virus, fa in virs.items():
                if not fa.exists() or fa.stat().st_size == 0:
                    continue
                _, seq = _read_first_sequence(fa)
                if not seq:
                    continue
                n_ratio = seq.count('N') / len(seq) * 100 if seq else 100.0
                best_aa, best_len = '', 0
                for frame in range(3):
                    sub = seq[frame:]
                    sub = sub[:len(sub) - len(sub) % 3]
                    if not sub:
                        continue
                    try:
                        aa = translate(sub)
                    except Exception as e:  # noqa: BLE001
                        self.logger.warning(
                            f"  [{sname}] {virus} frame{frame} 翻译失败: {e}")
                        continue
                    for orf in aa.split('*'):
                        if len(orf) > best_len:
                            best_len, best_aa = len(orf), orf
                is_intact, stops = check_protein_integrity(best_aa)
                # best_aa 为空 = 三个 frame 都没找到 ORF（序列太短/全 N），
                # check_protein_integrity('') 返回 (True, 0) 会误判 PASS
                if not best_aa:
                    is_intact = False
                qc = 'PASS' if (is_intact and n_ratio < 5.0) else 'FAIL'
                rows.append({
                    'Sample': sname, 'Reference': virus,
                    'Feature': f'ORF_{best_len}aa', 'Gene': 'NA',
                    'Length': len(seq), 'N_ratio(%)': f'{n_ratio:.2f}',
                    'Stops': stops, 'QC': qc,
                })
        if not rows:
            return
        qc_tsv = self.out_dir / 'consensus' / 'consensus_qc.tsv'
        with open(qc_tsv, 'w', encoding='utf-8') as f:
            cols = ('Sample', 'Reference', 'Feature', 'Gene', 'Length',
                    'N_ratio(%)', 'Stops', 'QC')
            f.write('\t'.join(cols) + '\n')
            for r in rows:
                f.write('\t'.join(str(r[k]) for k in cols) + '\n')
        n_pass = sum(1 for r in rows if r['QC'] == 'PASS')
        self.logger.info(f"  蛋白 QC: {n_pass}/{len(rows)} 通过 -> {qc_tsv}")

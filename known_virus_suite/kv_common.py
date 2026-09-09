#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_common.py — known_virus_suite 公共层
=====================================
参考加载、日志、I/O、外部工具探测。

设计来源: virome_analysis_pipeline/batch_virus_depth.py 的
_load_reference_lengths / _load_ref_info_smart / setup_logging 等
（复制后改造，原文件不动）
"""

import gzip
import logging
import os
import shutil
import sys
import time
from pathlib import Path

# ── 平台适配 ──────────────────────────────────────────────
IS_WINDOWS = os.name == 'nt'


def smart_open(path, mode='rt'):
    """透明处理 .gz"""
    path = str(path)
    return gzip.open(path, mode, encoding=None if 'b' in mode else 'utf-8') if path.endswith('.gz') else open(path, mode, encoding=None if 'b' in mode else 'utf-8')


def which_tool(name, extra_dirs=None):
    """
    工具探测，兼顾 Windows 上的 .exe 和平台自带 tools 目录。
    返回绝对路径，找不到返回 None。
    """
    p = shutil.which(name)
    if p:
        return p
    # Windows: shutil.which 对没有 .exe 后缀的名字也能找到，但为稳妥再试一次
    if IS_WINDOWS:
        for suffix in ('.exe', '.bat', '.cmd'):
            p = shutil.which(name + suffix)
            if p:
                return p
    for d in (extra_dirs or []):
        d = Path(d)
        for cand in (d / name, d / f"{name}.exe"):
            if cand.is_file():
                return str(cand)
    return None


class ToolRegistry:
    """
    记录本模块需要的所有外部工具，支持缺件降级判断。

    用法:
        reg = ToolRegistry(logger)
        reg.probe('minimap2', extra_dirs=[...])
        if reg.has('minimap2'): ...
    """

    NEEDED = {
        'salmon': 'salmon 伪比对引擎',
        'minibwa': 'minibwa 真比对引擎',
        'minimap2': 'minimap2 比对（共识段）',
        'mafft': 'mafft 多序列比对（共识迭代）',
        'samtools': 'samtools（SAM/BAM 处理与覆盖度统计）',
        'pandepth': 'pandepth（覆盖度/深度统计，口径金标准）',
        'viral_consensus': 'viral_consensus（共识序列构建，对齐原管线）',
        'bcftools': 'bcftools（变异检出 caller，替代 freebayes/lofreq/ivar）',
    }

    # 工具可执行文件名与所在子目录（Windows 本机 tools\ 布局）
    TOOL_HINTS = {
        'salmon':   ('salmon',   r'salmon2'),
        'minibwa':  ('minibwa',  r'bin'),
        'minimap2': ('minimap2', r'minimap2'),
        'mafft':    ('mafft',    r'mafft-win\ms\bin'),
        'samtools': ('samtools', r'samtools\bin'),
        'pandepth': ('pandepth', r'pandepth'),
        'viral_consensus': ('viral_consensus', r'viral_consensus'),
        'bcftools': ('bcftools', r'bcftools\bin'),
    }

    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger('kv')
        self.paths = {}
        self.env = None

    def probe(self, extra_dirs=None):
        # tools\ 在模块上一级；部分工具（如 minibwa）直接放在平台根下
        root = Path(__file__).resolve().parent.parent
        base = root / 'tools'
        auto = []
        for name, (exe, sub) in self.TOOL_HINTS.items():
            for d in (base / sub, root / sub):
                if d.exists():
                    auto.append(str(d))
        dirs = list(dict.fromkeys(auto + [str(Path(d)) for d in (extra_dirs or [])]))
        for name in self.NEEDED:
            exe = self.TOOL_HINTS.get(name, (name, ''))[0]
            self.paths[name] = which_tool(exe, dirs)
        self.env = self._build_env(root)
        return dict(self.paths)

    def _build_env(self, root):
        """构造子进程环境：bcftools 插件目录 + 工具目录入 PATH。

        bcftools 的 +插件 与某些子命令需要 BCFTOOLS_PLUGINS 才能找到
        libexec/bcftools 下的插件；Windows 下不设会直接报插件未找到。
        PATH 里追加工具目录，让 samtools/bcftools 能互相调用。
        """
        env = dict(os.environ)
        plugins = Path(root) / 'tools' / 'bcftools' / 'libexec' / 'bcftools'
        if plugins.is_dir():
            env['BCFTOOLS_PLUGINS'] = str(plugins)
        extra = [str(Path(p).parent) for p in self.paths.values() if p]
        if extra:
            env['PATH'] = os.pathsep.join(
                list(dict.fromkeys(extra)) + [env.get('PATH', '')])
        return env

    def has(self, name):
        return bool(self.paths.get(name))

    def require(self, *names):
        """缺件直接报错，附安装提示"""
        missing = [n for n in names if not self.has(n)]
        if missing:
            lines = [f"缺少必要工具: {', '.join(missing)}"]
            for m in missing:
                lines.append(f"  - {m}  ({self.NEEDED.get(m, '')})")
            raise RuntimeError('\n'.join(lines))
        return True

    def report(self):
        rows = []
        for name, desc in self.NEEDED.items():
            p = self.paths.get(name)
            rows.append(f"  [{'OK' if p else 'MISS':>4}] {name:<12} {desc}")
            if p:
                rows[-1] += f"\n         -> {p}"
        return '\n'.join(rows)


# ── 日志 ──────────────────────────────────────────────────
def setup_logger(out_dir, verbose=False, tag='kv'):
    """文件 + 控制台双通道日志"""
    out_dir = Path(out_dir)
    logs_dir = out_dir / 'logs'
    logs_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(tag)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers = []
    logger.propagate = False

    fh = logging.FileHandler(logs_dir / f"{tag}_{time.strftime('%Y%m%d_%H%M%S')}.log",
                             mode='a', encoding='utf-8')
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s', '%H:%M:%S'))
    logger.addHandler(ch)

    return logger, logs_dir


def fmt_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    if h > 0:
        return f"{int(h)}h {int(m)}m {int(s)}s"
    if m > 0:
        return f"{int(m)}m {int(s)}s"
    return f"{s:.1f}s"


# ── 参考库加载 ────────────────────────────────────────────
def load_ref_lengths(ref_fasta, logger=None):
    """读参考 FASTA，返回 {accession: length}"""
    lengths = {}
    cur_id, cur_len = None, 0
    with smart_open(ref_fasta, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if cur_id:
                    lengths[cur_id] = cur_len
                cur_id = line.strip().split()[0][1:]
                cur_len = 0
            else:
                cur_len += len(line.strip())
        if cur_id:
            lengths[cur_id] = cur_len
    if logger:
        logger.info(f"参考库载入: {len(lengths)} 条序列, 合计 {sum(lengths.values())/1e6:.2f} Mb")
    return lengths


def load_ref_sequences(ref_fasta, only=None):
    """
    读参考 FASTA 序列本体。only 为 accession 集合时只保留这些。
    返回 {accession: sequence}
    """
    seqs = {}
    cur_id, buf = None, []
    with smart_open(ref_fasta, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if cur_id and (only is None or cur_id in only):
                    seqs[cur_id] = ''.join(buf).upper()
                cur_id = line.strip().split()[0][1:]
                buf = []
            else:
                buf.append(line.strip())
        if cur_id and (only is None or cur_id in only):
            seqs[cur_id] = ''.join(buf).upper()
    return seqs


# ref_info 列名同义词表（对齐原管线 _load_ref_info_smart）
ACC_SYNONYMS = ['Accession', 'accession', 'Virus GENBANK accession', 'ID']
TAX_SYNONYMS = ['Taxid', 'taxonomy_id', 'taxid']
SP_SYNONYMS = ['Species_NCBI', 'Species_ICTV', 'taxonomy', 'Species', 'description', 'Virus name(s)']
SEG_SYNONYMS = ['Segment', 'segment']


def load_ref_info(ref_info_tsv, logger=None):
    """
    读参考注释表，返回
      {accession: {'taxid':..,'species':..,'segment':..,'molecule':..}}
    同时建立去版本号的别名索引。
    """
    info = {}
    if not ref_info_tsv or not os.path.exists(ref_info_tsv):
        if logger:
            logger.warning("未提供 ref_info，物种名将回退为 accession")
        return info

    idx = {}
    with open(ref_info_tsv, 'r', encoding='utf-8', errors='replace') as f:
        header = None
        for line in f:
            if line.startswith('####') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if header is None:
                header = [h.strip() for h in parts]

                def find(syns):
                    return next((header.index(s) for s in syns if s in header), -1)

                idx = {
                    'acc': find(ACC_SYNONYMS), 'tax': find(TAX_SYNONYMS),
                    'sp': find(SP_SYNONYMS), 'seg': find(SEG_SYNONYMS),
                    'mol': find(['Molecule_type2', 'Molecule_Type2', 'Molecule_type']),
                }
                if idx['acc'] == -1:
                    if logger:
                        logger.error(f"ref_info 找不到 Accession 列，实际列: {header[:8]}")
                    return {}
                continue

            if len(parts) <= idx['acc']:
                continue
            acc = parts[idx['acc']].strip()
            if not acc:
                continue

            def get(key):
                i = idx.get(key, -1)
                return parts[i].strip() if i != -1 and len(parts) > i else ''

            rec = {
                'taxid': get('tax') or 'Unannotated',
                'species': get('sp') or acc,
                'segment': get('seg'),
                'molecule': get('mol'),
            }
            info[acc] = rec
            # 去版本号别名：同一 accession 多版本时后写入者胜出，
            # 会让物种名静默错配。检测到冲突就报警。
            bare = acc.split('.')[0]
            if bare in info and info[bare].get('species') != rec.get('species'):
                if logger:
                    logger.warning(f"ref_info 多版本冲突: {bare} "
                                   f"({info[bare].get('species')} vs "
                                   f"{rec.get('species')})，使用后者")
            info[bare] = rec  # 去版本号别名

    if logger:
        logger.info(f"参考注释载入: {len(info)//2} 条（含去版本号别名）")
    return info


def ref_lookup(ref_info, accession):
    """按 accession 查注释，容错版本号差异"""
    if not ref_info:
        return {'taxid': 'Unannotated', 'species': accession, 'segment': '', 'molecule': ''}
    rec = ref_info.get(accession) or ref_info.get(accession.split('.')[0])
    return rec or {'taxid': 'Unannotated', 'species': accession, 'segment': '', 'molecule': ''}

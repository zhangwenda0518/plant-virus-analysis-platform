"""变异标注段（bcftools caller + SnpEff 封装）—— known_virus_suite 第五段候选。

设计要点
--------
1. caller 用 bcftools mpileup + call。原管线的 freebayes/lofreq/ivar 在
   Windows 上均不可得（MSYS2 与 conda 都没有 win-64 包），bcftools 已实测可用。
   参数与阈值均来自实测（见 known_virus_suite/POSCOUNTS_REMOVAL_PLAN.md §4.2）：
     - 不加 --ploidy 1（会丢光低频变异）
     - 不加 -A（与 -mv 组合导致 0 变异）
     - QUAL 下限 3.5（bcftools 标度远低于 freebayes，原管线的 >20 不可移植）
     - -a 只能写 FORMAT/AD,FORMAT/DP,INFO/AD（INFO/DP4 非法）
2. 过滤表达式字段映射：原管线 freebayes 的 AO/SAF/SAR/AF 换成
   AD[1] / DP4[2] / DP4[3] / AD[1]/DP。
3. 建库输入优先复用绘图段已下载的 GenBank 缓存（<out>/gbk_files/*.gb）；
   缺失时用 Bio.Entrez.efetch 下载（与 kv_plot.fetch_gb 同一接口）。
4. SnpEff 需要 Java 21（5.4c 编译为 class 65），平台自带 tools/jdk21。
5. SnpEff 的 data.dir 走相对路径，工作目录必须切到 config 所在目录，
   否则 Windows 绝对路径里的反斜杠会被当转义吃掉。
6. 库内染色体名保留 GenBank 版本号（OR489165.1），VCF 的 CHROM 必须
   与之一致，否则报 ERROR_CHROMOSOME_NOT_FOUND。这里做自动对齐。

用法
----
python kv_variant.py --reference <ref.fasta> --out <out_dir> \
    --bam-dir <dir> [--gbk-dir <dir>] [--ncbi-email <mail>] \
    [--qual 3.5] [--min-freq 0.05]
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path

LOGGER = logging.getLogger('kv_variant')


# ── 编码兼容：外部工具产物可能是 GBK ────────────────────────────────
# Windows 中文环境下 bcftools/SnpEff 等 C/Java 工具写 header 时用的是
# 系统 ANSI 代码页（936/GBK），而不是 UTF-8。若用 encoding='utf-8'
# 加 errors='replace' 去读再写回，GBK 字节会被替换成 U+FFFD，
# 造成不可逆的信息丢失（「桌面」变成一串替换符）。
# 凡是「读取之后还要写回」的场景，都用 read_text_lossless：
# 依次尝试 utf-8 / gbk，任一成功即原样保留；都失败才退到 replace。
def read_text_lossless(path: Path) -> str:
    """读文本文件，优先 UTF-8，其次 GBK，都失败才用 replace 兜底。"""
    raw = Path(path).read_bytes()
    for enc in ('utf-8-sig', 'utf-8', 'gbk'):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode('utf-8', errors='replace')


# ── 工具定位 ────────────────────────────────────────────────────────
def _platform_root() -> Path:
    return Path(__file__).resolve().parent.parent


def find_java() -> Path | None:
    """优先平台自带 JDK 21，其次 PATH 里的 java。"""
    root = _platform_root()
    cand = sorted((root / 'tools' / 'jdk21').glob('jdk-*/bin/java.exe'))
    if cand:
        return cand[-1]
    cand = sorted((root / 'tools' / 'jdk21').glob('jdk-*/bin/java'))
    if cand:
        return cand[-1]
    w = shutil.which('java')
    return Path(w) if w else None


def find_snpeff() -> Path | None:
    root = _platform_root()
    for cand in [
        root / 'tools' / 'snpeff' / 'snpEff' / 'snpEff.jar',
        root / 'tools' / 'snpeff' / 'snpEff.jar',
    ]:
        if cand.exists():
            return cand
    return None


def java_major_version(java: Path) -> int | None:
    """探测 java 主版本号。"""
    try:
        p = subprocess.run([str(java), '-version'], capture_output=True, text=True,
                           encoding='utf-8', errors='replace', timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    text = (p.stderr or '') + (p.stdout or '')
    for token in text.replace('"', ' ').split():
        if token[:1].isdigit() and '.' in token:
            try:
                return int(token.split('.')[0])
            except ValueError:
                continue
    return None


# ── 输入收集 ────────────────────────────────────────────────────────
def read_fasta_ids(fasta: Path) -> list[str]:
    ids = []
    with open(fasta, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids


def collect_bams(bam_dir: Path) -> dict[str, Path]:
    """收集 BAM 目录下的比对文件，返回 {样本名: bam 路径}。

    多个 BAM 视为多样本，一次性喂给 bcftools mpileup（与 bam 段产物对齐）。
    """
    out: dict[str, Path] = {}
    if not bam_dir or not Path(bam_dir).is_dir():
        return out
    for p in sorted(Path(bam_dir).rglob('*.bam')):
        if p.name.endswith('.tmp.bam'):
            continue
        name = p.name[:-4]
        if name.endswith('.sorted'):      # 共识段产物 S1.sorted.bam -> S1
            name = name[: -len('.sorted')]
        out[name] = p
    return out


# ── bcftools caller ─────────────────────────────────────────────────
class BcftoolsCaller:
    """bcftools mpileup + call + 两步 filter。

    全部参数与阈值为实测定稿（POSCOUNTS_REMOVAL_PLAN.md §4.2）：
      - 不设 --ploidy 1：实测导致变异全丢（默认二倍体出 1 条，ploidy 1 出 0 条）
      - 不加 -A：与 -mv 组合实测出 0 条
      - -a 只能写 FORMAT/AD,FORMAT/DP,INFO/AD：INFO/DP4 不是合法标签
      - QUAL 下限 3.5：bcftools QUAL 标度远低于 freebayes，原管线 >20 不可移植
    """

    # 原管线 freebayes 过滤字段 -> bcftools 字段（实测定稿）
    FIELD_MAP = {
        'INFO/DP': 'INFO/DP',
        'INFO/AO': 'INFO/AD[1]',
        'INFO/AF': 'INFO/AD[1]/INFO/DP',
        'INFO/SAF': 'INFO/DP4[2]',
        'INFO/SAR': 'INFO/DP4[3]',
    }

    def __init__(self, bcftools: str, threads: int = 4,
                 logger: logging.Logger | None = None, env=None):
        self.bcftools = str(bcftools)
        self.threads = threads
        self.log = logger or LOGGER
        self.env = env

    def _run(self, args: list[str], label: str, log_file: Path | None = None) -> None:
        cmd = [self.bcftools] + args
        p = subprocess.run(cmd, capture_output=True, text=True,
                           encoding='utf-8', errors='replace',
                           env=self.env)
        if log_file is not None:
            log_file = Path(log_file)
            log_file.parent.mkdir(parents=True, exist_ok=True)
            with open(log_file, 'a', encoding='utf-8') as fh:
                fh.write('$ ' + ' '.join(cmd) + '\n')
                fh.write((p.stdout or '') + (p.stderr or '') + '\n')
        if p.returncode != 0:
            tail = ((p.stderr or p.stdout or '').strip().splitlines() or [''])[-4:]
            raise RuntimeError(f'bcftools {label} 失败(rc={p.returncode}): '
                               + ' | '.join(tail))

    @staticmethod
    def build_filter_expr(dp: int, saf: int, qual: float = 3.5,
                          min_freq: float = 0.05) -> str:
        """过滤表达式，字段名已按 BcftoolsCaller.FIELD_MAP 换过。"""
        return (f'QUAL>={qual} && INFO/DP>={dp} '
                f'&& INFO/AD[1]/INFO/DP>={min_freq} '
                f'&& (INFO/DP4[2]+INFO/DP4[3])>={saf}')

    @staticmethod
    def dynamic_thresholds(mean_depth: float) -> tuple[int, int]:
        """照原管线（batch_virus_variants.py line 461-474）按平均深度选 dp/saf。"""
        if mean_depth < 50:
            return 10, 2
        if mean_depth < 1000:
            return 20, 3
        return 100, 10

    def call(self, ref_fa: Path, bams: list[Path], out_vcf: Path,
             mean_depth: float, qual: float = 3.5,
             min_freq: float = 0.05, min_depth_mpileup: int = 100000,
             work_dir: Path | None = None,
             log_file: Path | None = None) -> tuple[Path, str]:
        """跑完整 caller 链路，返回（过滤后 VCF 路径, 过滤表达式）。

        产物：
          <out_vcf>           未经 FILTER 标记的原始变异
          <out_vcf>.filtered  仅 PASS 的变异
        """
        out_vcf = Path(out_vcf)
        out_vcf.parent.mkdir(parents=True, exist_ok=True)
        work = Path(work_dir or out_vcf.parent)
        work.mkdir(parents=True, exist_ok=True)
        stem = out_vcf.stem
        pile = work / f'{stem}.pileup.bcf'
        callb = work / f'{stem}.call.bcf'
        raw = work / f'{stem}.raw.vcf'
        soft = work / f'{stem}.soft.vcf'
        filtered = Path(str(out_vcf) + '.filtered')

        self._run(['mpileup', '-f', str(ref_fa),
                   '-d', str(min_depth_mpileup), '-q', '0', '-Q', '13', '-B',
                   '-a', 'FORMAT/AD,FORMAT/DP,INFO/AD',
                   '-Ob', '-o', str(pile)] + [str(b) for b in bams],
                  'mpileup', log_file)
        self._run(['call', '-mv', '-Ob', '-o', str(callb), str(pile)],
                  'call', log_file)
        self._run(['view', '-Ov', '-o', str(raw), str(callb)], 'view', log_file)

        dp, saf = self.dynamic_thresholds(mean_depth)
        expr = self.build_filter_expr(dp, saf, qual, min_freq)
        self._run(['filter', '-s', 'FAIL', '-i', expr, '-Ov',
                   '-o', str(soft), str(raw)], 'filter', log_file)
        self._run(['filter', '-i', 'FILTER=="PASS"', '-Ov',
                   '-o', str(filtered), str(soft)], 'filter-pass', log_file)
        # 同时保留一份未过滤 VCF 供追溯
        if read_text_lossless(raw).strip() and not out_vcf.exists():
            shutil.copy2(raw, out_vcf)
        return filtered, expr


# ── SnpEff 建库 / 注释 ──────────────────────────────────────────────
def has_coding_features(gbk: Path) -> bool:
    """GenBank 记录是否含可注释的编码 feature（CDS/gene/mRNA/tRNA/rRNA 等）。

    viroid（PSTVd、CEVd 等）在 GenBank 里只有 source 一个 feature，
    无可编码单元，SnpEff 建库会失败。这里只读 FEATURES 段的关键行，
    不解析整个记录（GBK 可能几十 KB，但只需判定有无）。

    实测对照（本平台库）：
      NC_002030.1 Potato spindle tuber viroid   CDS=0 gene=0 exon=0 -> False
      NC_000885.1 Tomato chlorotic dwarf viroid CDS=0 gene=0 exon=0 -> False
      OR489165.1  Cytorhabdovirus sp. 'lycii'   CDS=6 gene=6          -> True
      LC902918.1  Potexvirus pepini             CDS=5 gene=5          -> True
    """
    keys = ('CDS', 'gene', 'mRNA', 'tRNA', 'rRNA', 'ncRNA')
    try:
        with open(gbk, encoding='utf-8', errors='replace') as fh:
            for line in fh:
                # feature 行的形式是「5 个空格 + 关键字」，列 0 不是空格
                if not line.startswith('     ') or line[:6].strip() == '':
                    continue
                tok = line[5:21].strip().split()[0] if line[5:21].strip() else ''
                if tok in keys:
                    return True
    except OSError:
        return True  # 读不到就给 SnpEff 自己判定，不掩盖真实错误
    return False


class SnpEffRunner:
    def __init__(self, java: Path, jar: Path, work_dir: Path, data_dir: Path,
                 logger: logging.Logger | None = None):
        self.java = Path(java)
        self.jar = Path(jar)
        self.work_dir = Path(work_dir).resolve()
        self.data_dir = Path(data_dir).resolve()
        self.log = logger or LOGGER
        self.config = self.work_dir / 'snpEff.config'

    def _run(self, args: list[str], label: str) -> subprocess.CompletedProcess:
        cmd = [str(self.java), '-Dfile.encoding=UTF-8', '-jar', str(self.jar)] + args
        self.log.info('  [SnpEff] %s: %s', label, ' '.join(args[:6]))
        # 不用 text=True：SnpEff 的 header 行含路径，可能是 GBK 字节，
        # 交给调用方用 read_text_lossless 统一解码，避免在这里被 replace 掉。
        p = subprocess.run(cmd, cwd=str(self.work_dir), capture_output=True)
        out = (p.stdout or b'')
        err = (p.stderr or b'')
        for enc in ('utf-8', 'gbk'):
            try:
                out_s, err_s = out.decode(enc), err.decode(enc)
                break
            except UnicodeDecodeError:
                continue
        else:
            out_s = out.decode('utf-8', errors='replace')
            err_s = err.decode('utf-8', errors='replace')
        p.stdout, p.stderr = out_s, err_s
        if p.returncode != 0:
            tail = (p.stderr or p.stdout or '').strip().splitlines()[-5:]
            raise RuntimeError(f'SnpEff {label} 失败(rc={p.returncode}): ' + ' | '.join(tail))
        return p

    def write_config(self, genomes: dict[str, str]) -> Path:
        """写自定义 config。data.dir 用相对路径规避 Windows 反斜杠转义。"""
        lines = ['# kv_variant auto-generated', 'data.dir = ./data', '']
        for gid, desc in genomes.items():
            lines.append(f'{gid}.genome : {desc}')
        self.config.write_text('\n'.join(lines) + '\n', encoding='utf-8')
        self.data_dir.mkdir(parents=True, exist_ok=True)
        return self.config

    def build(self, genome_id: str, gbk: Path) -> Path:
        """GenBank -> snpEffectPredictor.bin。

        SnpEff 的 `build -genbank` 要求记录里存在可编码的 exon。
        viroid（PSTVd/CEVd 等）是裸 RNA，GenBank 里只有 source 一个 feature，
        没有 CDS/gene/exon，SnpEff 会报
        `FATAL ERROR: Most Exons do not have sequences!`。
        这是数据特性而非工具故障，所以先做一次预检，
        给出自解释的错误信息（调用方按“跳过该参考”处理）。
        """
        target = self.data_dir / genome_id
        target.mkdir(parents=True, exist_ok=True)
        dst = target / 'genes.gbk'
        if not dst.exists() or dst.stat().st_size != Path(gbk).stat().st_size:
            shutil.copy2(gbk, dst)

        if not has_coding_features(gbk):
            raise RuntimeError(
                f'{genome_id} 的 GenBank 记录无 CDS/gene 注释'
                f'（viroid 等非编码病毒），SnpEff 无法建库')

        self._run(['build', '-genbank', '-c', str(self.config), '-noLog', genome_id],
                  f'build {genome_id}')
        return target / 'snpEffectPredictor.bin'

    def annotate(self, genome_id: str, vcf: Path, out_vcf: Path) -> Path:
        p = self._run(['ann', '-c', str(self.config), '-noLog', genome_id, str(vcf)],
                      f'ann {genome_id}')
        out_vcf.write_text(p.stdout, encoding='utf-8')
        return out_vcf


# ── ANN 解析 ───────────────────────────────────────────────────────
ANN_KEYS = ['Allele', 'Annotation', 'Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type',
            'Feature_ID', 'BioType', 'Rank', 'HGVS_c', 'HGVS_p', 'cDNA', 'CDS', 'AA',
            'Distance', 'Note']
CODING_ANN = {
    'missense_variant', 'synonymous_variant', 'stop_gained', 'stop_lost',
    'start_lost', 'start_retained_variant', 'stop_retained_variant',
    'frameshift_variant', 'inframe_insertion', 'inframe_deletion',
    'disruptive_inframe_insertion', 'disruptive_inframe_deletion',
    'conservative_inframe_insertion', 'conservative_inframe_deletion',
    'initiator_codon_variant', 'terminator_codon_variant', 'coding_sequence_variant',
}
IMPACT_ORDER = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}


def parse_ann(vcf_path: Path, coding_only: bool = True) -> list[dict]:
    """解析 SnpEff 注释 VCF 的 ANN 字段。"""
    rows: list[dict] = []
    for line in read_text_lossless(vcf_path).splitlines():
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 8:
            continue
        chrom, pos, _id, ref, alt, _qual, filt, info = f[:8]
        ann_field = None
        for item in info.split(';'):
            if item.startswith('ANN='):
                ann_field = item[4:]
                break
        if not ann_field:
            continue
        for entry in ann_field.split(','):
            parts = entry.split('|')
            parts += [''] * (16 - len(parts))
            rec = dict(zip(ANN_KEYS, parts[:16]))
            if coding_only and rec['Annotation'] not in CODING_ANN:
                continue
            rec.update({'CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
                        'FILTER': filt})
            rows.append(rec)
    return rows


ANN_COLS = ['CHROM', 'POS', 'REF', 'ALT', 'Gene_Name', 'Annotation', 'Impact',
            'HGVS_c', 'HGVS_p', 'AA', 'CDS', 'Feature_ID', 'Note']


def write_ann_tsv(rows: list[dict], out_path: Path) -> int:
    lines = ['\t'.join(ANN_COLS)]
    for r in sorted(rows, key=lambda x: (int(x['POS']),
                                         IMPACT_ORDER.get(x['Impact'], 9))):
        lines.append('\t'.join(str(r.get(c, '')) for c in ANN_COLS))
    Path(out_path).write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return len(rows)


# ── 主流程 ──────────────────────────────────────────────────────────
def detect_contig_name(gbk: Path) -> str:
    """从 GenBank 提取序列 ID（保留版本号，与 SnpEff 库内命名一致）。"""
    acc = None
    with open(gbk, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('VERSION'):
                parts = line.split()
                if len(parts) > 1:
                    return parts[1].strip()
            if line.startswith('ACCESSION'):
                parts = line.split()
                if len(parts) > 1 and acc is None:
                    acc = parts[1].strip()
            if line.startswith('ORIGIN'):
                break
    return acc or gbk.stem


# ── 注释导出（virus-annotations/，对齐上游管线契约）──────────────────
# 照搬 virome_analysis_pipeline/batch_virus_variants.py:942-979
# build_snpeff_db 里的 GB/GTF 复原逻辑：
#   <out>/virus-annotations/<acc>.gb    逐条 GenBank
#   <out>/virus-annotations/<acc>.gtf   仅 CDS，gb2gtf 手写格式
# 上游不做长度过滤，无 CDS（如类病毒）时 GTF 为空文件但保留。
def export_annotations(gbk_map: dict, out_dir, logger) -> dict:
    """把 GenBank 缓存导出成 virus-annotations/{acc}.gb + {acc}.gtf。

    gbk_map: {accession: gb_path}
    返回 {'gb': n, 'gtf': n, 'empty_gtf': [acc...], 'dir': <path>}
    """
    ann_dir = Path(out_dir) / 'virus-annotations'
    ann_dir.mkdir(parents=True, exist_ok=True)
    n_gb = n_gtf = 0
    empty_gtf: list[str] = []
    try:
        from Bio import SeqIO
    except ImportError:
        logger.warning('  Biopython 缺失，跳过 virus-annotations 导出')
        return {'gb': 0, 'gtf': 0, 'empty_gtf': [], 'dir': str(ann_dir)}

    for acc, gb_path in gbk_map.items():
        try:
            record = next(SeqIO.parse(str(gb_path), 'genbank'))
        except Exception as e:  # noqa: BLE001
            logger.warning('  %s: GenBank 解析失败，跳过注释导出（%s）', acc, e)
            continue

        gb_file = ann_dir / f'{acc}.gb'
        gtf_file = ann_dir / f'{acc}.gtf'
        gid = record.id
        SeqIO.write(record, gb_file, 'genbank')
        n_gb += 1

        tid = 0
        with open(gtf_file, 'w', encoding='utf-8', newline='\n') as fh:
            for f in record.features:
                if f.type != 'CDS':
                    continue
                tid += 1
                keys = f.qualifiers.keys()
                gene_name = f'CDS_{tid}'
                if 'gene' in keys:
                    gene_name = f.qualifiers['gene'][0]
                elif 'product' in keys:
                    gene_name = f.qualifiers['product'][0]
                elif 'label' in keys:
                    gene_name = f.qualifiers['label'][0]
                gene_name = (str(gene_name).replace('"', '')
                             .replace(';', '_').replace(' ', '_')
                             .replace('=', '_'))
                attr = (f'gene_id "{gene_name}"; '
                        f'transcript_id "{gene_name}"; '
                        f'gene_name "{gene_name}";')
                if f.location.strand == 1:
                    strand = '+'
                elif f.location.strand == -1:
                    strand = '-'
                else:
                    strand = '.'
                start = int(f.location.start) + 1
                end = int(f.location.end)
                # 与上游同一格式：acc  gb2gtf  CDS  start  end  .  strand  0  attr
                fh.write(f'{gid}\tgb2gtf\tCDS\t{start}\t{end}\t.\t'
                         f'{strand}\t0\t{attr}\n')
        if tid:
            n_gtf += 1
        else:
            empty_gtf.append(acc)

    logger.info('virus-annotations/ 导出 %d 个 .gb，%d 个 .gtf%s',
                n_gb, n_gtf,
                f'（{len(empty_gtf)} 条无 CDS，GTF 为空）' if empty_gtf else '')
    return {'gb': n_gb, 'gtf': n_gtf, 'empty_gtf': empty_gtf,
            'dir': str(ann_dir)}


# ── GenBank 获取 ────────────────────────────────────────────────────
def fetch_gb(accession: str, gbk_dir: Path, logger, email=None,
             api_key=None) -> Path | None:
    """下载 GenBank 到缓存目录，已存在则复用（同 kv_plot.fetch_gb）。"""
    gbk_dir = Path(gbk_dir)
    gbk_dir.mkdir(parents=True, exist_ok=True)
    base = str(accession).split('.')[0]
    dst = gbk_dir / f'{accession}.gb'
    if dst.exists() and dst.stat().st_size > 0:
        return dst
    for alt in (gbk_dir / f'{base}.gb', gbk_dir / f'{base}.gbk',
                gbk_dir / f'{accession}.gbk'):
        if alt.exists() and alt.stat().st_size > 0:
            return alt
    try:
        from Bio import Entrez
    except ImportError:
        logger.warning('  Biopython 缺失，无法下载 GenBank 注释')
        return None
    Entrez.email = email or os.environ.get('NCBI_EMAIL') or 'unknown@example.org'
    key = api_key or os.environ.get('NCBI_API_KEY')
    if key:
        Entrez.api_key = key
    try:
        handle = Entrez.efetch(db='nuccore', id=base, rettype='gb', retmode='text')
        data = handle.read()
        handle.close()
        if not data or len(data) < 100:
            logger.warning('  GenBank 下载内容异常 %s', accession)
            return None
        dst.write_text(data if isinstance(data, str) else data.decode(),
                       encoding='utf-8')
        logger.debug('    GenBank 下载成功: %s (%dB)', accession, len(data))
        return dst
    except Exception as e:  # noqa: BLE001
        logger.warning('  GenBank 下载失败 %s: %s', accession, e)
        return None


# ── 段落类 ──────────────────────────────────────────────────────────
class VariantStage:
    """变异功能注释段（SnpEff 封装）。"""

    def __init__(self, args, tools, logger, out_dir, cfg=None):
        self.args = args
        self.tools = tools
        self.log = logger
        self.out_dir = Path(out_dir).resolve()
        self.cfg = cfg or {}
        # 注释导出统计（virus-annotations/），供调用方读取
        self._annotations: dict = {}

    def run(self) -> list[dict]:
        cfg = self.cfg
        out_dir = self.out_dir
        out_dir.mkdir(parents=True, exist_ok=True)

        java = find_java()
        jar = find_snpeff()
        if not java or not jar:
            self.log.error('缺少 Java 或 SnpEff：java=%s jar=%s', java, jar)
            return []
        major = java_major_version(java)
        self.log.info('Java: %s (major=%s)', java, major)
        if major is not None and major < 21:
            self.log.error('SnpEff 5.4c 需要 Java 21+，当前 major=%s', major)
            return []

        ref = Path(getattr(self.args, 'reference', None) or '').resolve()
        if not ref.exists():
            self.log.error('参考序列不存在: %s', ref)
            return []
        ref_ids = read_fasta_ids(ref)
        self.log.info('参考序列 %s：%d 条', ref.name, len(ref_ids))

        gbk_dir = Path(cfg.get('gbk_dir') or (out_dir / 'gbk_files')).resolve()

        # ── 1. 取变异 VCF ──
        # 两条入口：
        #   a) 外部 VCF 直入（cfg['input_vcf']）：跳过 caller，直接进注释
        #   b) BAM 入口（cfg['bam_dir']，默认 <out>/bam）：跑 bcftools caller
        vcf_dir = out_dir / 'vcf'
        vcf_dir.mkdir(parents=True, exist_ok=True)
        input_vcf = cfg.get('input_vcf')
        if input_vcf:
            filtered_vcf = Path(input_vcf).resolve()
            if not filtered_vcf.is_file():
                self.log.error('指定的 VCF 不存在: %s', filtered_vcf)
                return []
            self.log.info('使用外部 VCF 直接注释（跳过 caller）：%s',
                          filtered_vcf.name)
        else:
            bam_dir = cfg.get('bam_dir') or (out_dir / 'bam')
            bams = collect_bams(Path(bam_dir))
            self.log.info('BAM (%s) 命中 %d 个样本', Path(bam_dir).name, len(bams))
            if not bams:
                self.log.error('未找到 BAM，请先用 bam 段生成比对，或用 --input-vcf 直接给 VCF')
                return []

            caller = BcftoolsCaller(cfg.get('bcftools') or 'bcftools',
                                   threads=cfg.get('threads', 4), logger=self.log,
                                   env=cfg.get('env'))
            mean_depth = self._mean_depth(list(bams.values()))
            # caller.call 返回 (仅 PASS 的 VCF, 过滤表达式)；
            # all.variants.vcf 本体是未打 FILTER 的原始变异，另存供追溯。
            filtered_vcf, expr = caller.call(
                ref, list(bams.values()), vcf_dir / 'all.variants.vcf',
                mean_depth, qual=cfg.get('qual', 3.5),
                min_freq=cfg.get('min_minor_freq', 0.05),
                work_dir=vcf_dir,
                log_file=out_dir / 'logs' / 'bcftools.log')
            self.log.info('平均深度 %.1f -> 过滤表达式: %s', mean_depth, expr)

        if not Path(filtered_vcf).is_file():
            self.log.error('过滤后 VCF 不存在: %s（过滤步骤可能失败）', filtered_vcf)
            return []
        n_pass = self._count_variants(filtered_vcf)
        self.log.info('过滤后变异 %d 条 -> %s', n_pass, Path(filtered_vcf).name)
        if n_pass == 0:
            self.log.warning('过滤后无变异，跳过注释')
            return []

        # ── 2. 按 CHROM 拆 VCF（SnpEff 每库一条参考）──
        split = self._split_vcf_by_chrom(filtered_vcf, vcf_dir)
        self.log.info('拆分到 %d 条参考', len(split))
        targets = [(acc, split[acc]) for acc in ref_ids if acc in split]
        for acc in split:
            if acc not in ref_ids:
                targets.append((acc, split[acc]))

        work = out_dir / 'snpeff_work'
        data = work / 'data'
        work.mkdir(parents=True, exist_ok=True)
        runner = SnpEffRunner(java, jar, work, data, self.log)

        gbk_map: dict[str, Path] = {}
        for acc, _ in targets:
            base = acc.split('.')[0]
            found = None
            for cand in (gbk_dir / f'{acc}.gb', gbk_dir / f'{acc}.gbk',
                         gbk_dir / f'{base}.gb', gbk_dir / f'{base}.gbk'):
                if cand.exists() and cand.stat().st_size > 0:
                    found = cand
                    break
            if found is None:
                found = fetch_gb(acc, gbk_dir, self.log,
                                 cfg.get('ncbi_email'), cfg.get('ncbi_api_key'))
            if found is not None:
                gbk_map[acc] = found
        self.log.info('GenBank 就绪 %d / %d 条', len(gbk_map), len(targets))
        if not gbk_map:
            self.log.error('没有可用的 GenBank 注释，无法建库')
            return []

        # ── 2b. 导出 virus-annotations/（对齐上游管线产物契约）──
        # 上游 batch_virus_variants.py 在 build_snpeff_db 里把 GB/GTF 落盘到
        # <out>/virus-annotations/，下游（含 t-consensus 选参考后画图）直接取用，
        # 不必再从 gb 现场转换。
        try:
            ann_stat = export_annotations(gbk_map, out_dir, self.log)
            self._annotations = ann_stat
        except Exception as e:  # noqa: BLE001
            self.log.warning('virus-annotations 导出失败（不影响变异段）: %s', e)
            self._annotations = {}

        genomes = {acc.split('.')[0]: f'kv_variant {acc}' for acc in gbk_map}
        runner.write_config(genomes)

        results = []
        skipped: list[dict] = []

        for acc, gbk in gbk_map.items():
            gid = acc.split('.')[0]
            chrom = detect_contig_name(gbk)
            self.log.info('── %s (chrom=%s)', gid, chrom)

            # 拆分后的 VCF：CHROM 名可能与 GenBank 版本号不完全一致
            vcf = None
            for key in (chrom, acc, gid):
                if key in split:
                    vcf = split[key]
                    break
            if vcf is None:
                for key, val in split.items():
                    if key.startswith(gid) or key.split('.')[0] == gid:
                        vcf = val
                        break
            if vcf is None:
                self.log.info('  %s 无变异，跳过注释', gid)
                continue
            if vcf.stem != gid:
                renamed = vcf_dir / f'{gid}.vcf'
                if not renamed.exists() or renamed != vcf:
                    shutil.copy2(vcf, renamed)
                vcf = renamed
            n = self._count_variants(vcf)
            if n == 0:
                self.log.info('  %s 无变异（过滤后），跳过', gid)
                continue

            try:
                runner.build(gid, gbk)
            except RuntimeError as e:
                # viroid 等无编码基因的参考属正常情形，用 warning 不报 error
                self.log.warning('跳过注释 %s: %s', gid, e)
                skipped.append({'genome': gid, 'reason': str(e)})
                continue

            ann = out_dir / 'annotated' / f'{gid}.ann.vcf'
            ann.parent.mkdir(parents=True, exist_ok=True)
            try:
                runner.annotate(gid, vcf, ann)
            except RuntimeError as e:
                self.log.error('  注释失败 %s: %s', gid, e)
                continue

            ann_rows = parse_ann(ann, coding_only=True)
            ann_tsv = out_dir / 'annotated' / f'{gid}.ann.tsv'
            n_ann = write_ann_tsv(ann_rows, ann_tsv)
            self.log.info('  变异 %d 条，编码区注释 %d 条 -> %s',
                          n, n_ann, ann_tsv.name)
            for r in sorted(ann_rows, key=lambda x: int(x['POS']))[:20]:
                self.log.info('    %6s %s>%s %-6s %-22s %-9s %s %s',
                              r['POS'], r['REF'], r['ALT'], r['Gene_Name'],
                              r['Annotation'], r['Impact'], r['HGVS_c'],
                              r['HGVS_p'])

            results.append({'accession': acc, 'genome': gid, 'chrom': chrom,
                            'vcf': str(vcf), 'n_variants': n,
                            'n_coding_ann': n_ann,
                            'ann_vcf': str(ann), 'ann_tsv': str(ann_tsv)})

        # summary 总是写，即使全部跳过——下游（平台结果面板）靠它区分
        # 「无变异」「无编码基因被跳过」「真的失败」。
        import json
        summary = {
            'n_success': len(results),
            'n_skipped': len(skipped),
            'results': results,
            'skipped': skipped,
        }
        (out_dir / 'variant_summary.json').write_text(
            json.dumps(summary, indent=1, ensure_ascii=False),
            encoding='utf-8')
        self.log.info('变异注释完成：%d 条参考成功，%d 条跳过',
                      len(results), len(skipped))
        return results

    # ── 辅助 ─────────────────────────────────────────────────────
    def _mean_depth(self, bams: list[Path]) -> float:
        """用 samtools depth 估算平均深度（照原管线动态阈值用）。"""
        samtools = self.cfg.get('samtools') or 'samtools'
        try:
            p = subprocess.run([str(samtools), 'depth', '-a'] + [str(b) for b in bams],
                               capture_output=True, text=True, encoding='utf-8',
                               errors='replace', env=self.cfg.get('env'))
        except OSError as e:
            self.log.warning('samtools depth 失败（%s），按 dp=10/saf=2 处理', e)
            return 0.0
        vals = []
        for line in (p.stdout or '').splitlines():
            f = line.split('\t')
            if len(f) >= 3:
                try:
                    vals.append(int(f[2]))
                except ValueError:
                    continue
        if not vals:
            self.log.warning('samtools depth 无输出，按 dp=10/saf=2 处理')
            return 0.0
        return sum(vals) / len(vals)

    @staticmethod
    def _count_variants(vcf: Path) -> int:
        n = 0
        for line in read_text_lossless(vcf).splitlines():
            if line.startswith('#') or not line.strip():
                continue
            n += 1
        return n

    def _split_vcf_by_chrom(self, vcf: Path, out_dir: Path) -> dict[str, Path]:
        """按 CHROM 把 VCF 拆成每参考一份（SnpEff 每库一条参考）。

        命名用 chrom 去掉版本号（NC_000885.1 -> NC_000885.vcf），与
        annotated/{gid}.ann.tsv 的下游命名对齐。
        只有当同一次运行里两个不同 chrom 去版本号后撞名时（如 X.1 与 X.2
        同时存在）才退到用完整 chrom 名；不能拿“文件已存在”当撞名依据，
        否则断点重跑时所有产物都会漂移成带版本号的名字。
        """
        header: list[str] = []
        per: dict[str, list[str]] = {}
        for line in read_text_lossless(vcf).splitlines(keepends=True):
            if line.startswith('#'):
                header.append(line)
                continue
            if not line.strip():
                continue
            chrom = line.split('\t', 1)[0]
            per.setdefault(chrom, []).append(line)
        # 先统计去版本号后的 base 是否在同一轮内撞名
        base_count: dict[str, int] = {}
        for chrom in per:
            base_count[chrom.split('.')[0]] = \
                base_count.get(chrom.split('.')[0], 0) + 1
        out: dict[str, Path] = {}
        for chrom, rows in per.items():
            gid = chrom.split('.')[0]
            dst = (out_dir / f'{gid}.vcf' if base_count[gid] == 1
                   else out_dir / f'{chrom}.vcf')
            dst.write_text(''.join(header) + ''.join(rows), encoding='utf-8')
            out[chrom] = dst
        return out


# ── CLI ─────────────────────────────────────────────────────────────
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description='变异检出与功能注释（bcftools + SnpEff）')
    ap.add_argument('--out', required=True, help='输出目录')
    ap.add_argument('--reference', required=True, help='参考序列 FASTA')
    ap.add_argument('--bam-dir', default=None,
                    help='BAM 目录（bam 段产物；默认 <out>/bam）')
    ap.add_argument('--gbk-dir', default=None, help='GenBank 缓存目录')
    ap.add_argument('--ncbi-email', default=None, help='NCBI Entrez 邮箱')
    ap.add_argument('--ncbi-api-key', default=None, help='NCBI API key')
    ap.add_argument('--bcftools', default=None, help='bcftools 可执行文件路径')
    ap.add_argument('--samtools', default=None, help='samtools 可执行文件路径')
    ap.add_argument('--qual', type=float, default=3.5,
                    help='QUAL 下限（实测定稿 3.5，非原管线的 20）')
    ap.add_argument('--min-freq', type=float, default=0.05,
                    help='最小等位频率（AD[1]/DP）')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s %(levelname)s %(message)s', datefmt='%H:%M:%S')

    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    bcftools = args.bcftools or str(
        _platform_root() / 'tools' / 'bcftools' / 'bin' / 'bcftools.exe')
    samtools = args.samtools or str(
        _platform_root() / 'tools' / 'samtools' / 'bin' / 'samtools.exe')

    env = dict(os.environ)
    plugins = str(_platform_root() / 'tools' / 'bcftools' / 'libexec' / 'bcftools')
    if Path(plugins).is_dir():
        env['BCFTOOLS_PLUGINS'] = plugins

    cfg = {
        'bam_dir': args.bam_dir,
        'gbk_dir': args.gbk_dir,
        'ncbi_email': args.ncbi_email,
        'ncbi_api_key': args.ncbi_api_key,
        'qual': args.qual,
        'min_minor_freq': args.min_freq,
        'threads': args.threads,
        'bcftools': bcftools,
        'samtools': samtools,
        'env': env,
    }
    stage = VariantStage(args, None, LOGGER, out_dir, cfg=cfg)
    results = stage.run()
    return 0 if results else 4


if __name__ == '__main__':
    sys.exit(main())

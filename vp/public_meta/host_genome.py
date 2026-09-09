#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
📥 宿主参考基因组下载器 v1.0
============================
基于 NCBI datasets CLI，自动下载指定物种的：
  - 核基因组 (genome)
  - 注释文件 (GFF3)
  - 序列报告 (seq-report)

功能特性:
  - 自动解压与合并多文件
  - 序列名去重（处理多基因组合并时的重名问题）
  - 支持叶绿体/线粒体基因组的额外下载
  - 可选基因组大小过滤
  - 输出统一 FASTA 用于下游分析

用法:
  python host_genome.py --species "Lycium barbarum" --outdir ./host_genome
  python host_genome.py --species "Lycium barbarum" --include-organelles --ncbi-api xxx

（移植自 MMPV-RNA public_metadata_pipeline/download_host_genome.py；
 文件写入统一改走平台 vp.utils.safe_open；HTTP 统一走带域名校验的
 requests.Session；datasets CLI 不可用时自动改用 E-utilities + HTTPS 回退）
"""

import os
import sys
import re
import gzip
import argparse
import subprocess
import hashlib
from pathlib import Path
from datetime import datetime
from typing import List, Optional, Set, Tuple

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    from vp.utils import safe_open
else:
    from ..utils import safe_open

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_FTP_HOST = "ftp://ftp.ncbi.nlm.nih.gov/"
NCBI_FTP_HTTPS = "https://ftp.ncbi.nlm.nih.gov/"


def _http():
    """统一 HTTP 入口（仅允许 NCBI 官方域名，防 SSRF）。"""
    import requests

    class _Scoped:
        """只放行 eutils/ftp.ncbi 官方域名的极简 Session 包装。"""
        def __init__(self):
            self._s = requests.Session()

        @staticmethod
        def _check(url):
            host = url.split('//', 1)[-1].split('/', 1)[0]
            if host not in ('eutils.ncbi.nlm.nih.gov', 'ftp.ncbi.nlm.nih.gov'):
                raise ValueError(f"仅允许 NCBI 官方域名，拒绝: {host}")
            return url

        def get(self, url, **kw):
            kw.setdefault('allow_redirects', False)
            return self._s.get(self._check(url), **kw)

    return _Scoped()


# ==========================================
# 终端UI
# ==========================================
class UI:
    CYAN = '\033[96m'; GREEN = '\033[92m'; YELLOW = '\033[93m'
    RED = '\033[91m'; PURPLE = '\033[95m'; GRAY = '\033[90m'
    BOLD = '\033[1m'; RESET = '\033[0m'

    @staticmethod
    def ok(msg):    print(f"  {UI.GREEN}✓{UI.RESET} {msg}")
    @staticmethod
    def warn(msg):  print(f"  {UI.YELLOW}⚠{UI.RESET} {msg}")
    @staticmethod
    def err(msg):   print(f"  {UI.RED}✗{UI.RESET} {msg}")
    @staticmethod
    def info(msg):  print(f"  {UI.CYAN}→{UI.RESET} {msg}")
    @staticmethod
    def header(msg):
        print(f"\n{UI.PURPLE}{UI.BOLD}{'='*55}{UI.RESET}")
        print(f"{UI.PURPLE}{UI.BOLD} {msg}{UI.RESET}")
        print(f"{UI.PURPLE}{UI.BOLD}{'='*55}{UI.RESET}")


def check_datasets_cli() -> bool:
    """检查 NCBI datasets CLI 是否可用"""
    try:
        result = subprocess.run(
            ['datasets', '--version'],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            version = result.stdout.strip() or result.stderr.strip()
            UI.ok(f"NCBI datasets CLI 可用: {version}")
            return True
    except FileNotFoundError:
        pass
    except Exception:
        pass

    UI.warn("未找到 NCBI datasets CLI（将使用 E-utilities + HTTPS 直接下载）")
    return False


def open_fasta(path: str):
    """Context-managed opener for FASTA files (supports .gz compression)."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8')
    return open(path, 'r', encoding='utf-8')


def collect_fasta_files(directory: str) -> List[str]:
    """Recursively collect all FASTA/FA/FNA files (plain or .gz) under a directory."""
    files = []
    for root, _, filenames in os.walk(directory):
        for f in filenames:
            ext = os.path.splitext(f)[1]
            if f.endswith('.gz'):
                # 检查 .fna.gz 等双重后缀
                base = f[:-3]
                if any(base.endswith(e) for e in ['.fna', '.fasta', '.fa']):
                    files.append(os.path.join(root, f))
            elif ext in {'.fna', '.fasta', '.fa'}:
                files.append(os.path.join(root, f))
    return sorted(files)


def count_sequences(fasta_path: str) -> int:
    """Count sequences in a FASTA file by counting '>' header lines."""
    count = 0
    try:
        with open_fasta(fasta_path) as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception:
        pass
    return count


def merge_and_deduplicate(
    fasta_files: List[str],
    output_path: str,
    min_length: int = 0,
    prefix: str = ''
) -> Tuple[int, int]:
    """
    合并多个 FASTA 文件并去重序列名

    返回值: (写入序列数, 重复序列名数)
    """
    UI.header("合并与去重")

    seen_ids: Set[str] = set()
    total_written = 0
    duplicate_count = 0
    total_bp = 0

    with safe_open(output_path, 'wt') as out_f:
        out_f.write(f"# Merged host genome\n")
        out_f.write(f"# Created: {datetime.now().isoformat()}\n")
        out_f.write(f"# Source files: {len(fasta_files)}\n")
        out_f.write(f"#\n")

        for fa_path in fasta_files:
            basename = os.path.basename(fa_path)
            UI.info(f"处理: {basename}")
            local_count = 0

            try:
                with open_fasta(fa_path) as in_f:
                    current_seq = []
                    current_id = ''
                    current_full_header = ''

                    for line in in_f:
                        if line.startswith('>'):
                            # 写入上一条序列
                            if current_id and current_seq:
                                seq = ''.join(current_seq)
                                if len(seq) >= min_length:
                                    out_f.write(f">{current_full_header}\n")
                                    for i in range(0, len(seq), 60):
                                        out_f.write(seq[i:i+60] + '\n')
                                    total_written += 1
                                    total_bp += len(seq)

                            # 解析新序列头
                            full_header = line[1:].strip()
                            seq_id = full_header.split()[0]

                            # 序列名去重
                            if seq_id in seen_ids:
                                suffix = 1
                                while f"{seq_id}_dup{suffix}" in seen_ids:
                                    suffix += 1
                                new_id = f"{seq_id}_dup{suffix}"
                                seen_ids.add(new_id)
                                current_full_header = full_header.replace(seq_id, new_id, 1)
                                duplicate_count += 1
                            else:
                                seen_ids.add(seq_id)
                                current_full_header = full_header

                            current_id = seq_id
                            current_seq = []
                            local_count += 1
                        else:
                            current_seq.append(line.strip())

                    # 处理最后一条序列
                    if current_id and current_seq:
                        seq = ''.join(current_seq)
                        if len(seq) >= min_length:
                            out_f.write(f">{current_full_header}\n")
                            for i in range(0, len(seq), 60):
                                out_f.write(seq[i:i+60] + '\n')
                            total_written += 1
                            total_bp += len(seq)

                UI.ok(f"  {local_count} 条序列")

            except Exception as e:
                UI.err(f"  处理 {basename} 时出错: {e}")

    file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
    UI.ok(f"合并完成: {total_written} 条序列, {total_bp:,} bp, {file_size_mb:.1f} MB")
    if duplicate_count:
        UI.warn(f"重命名了 {duplicate_count} 个重复序列名")

    return total_written, duplicate_count


def generate_genome_report(output_dir: str, fasta_files: List[str], merged_fasta: str):
    """生成基因组下载汇总报告"""
    report_path = os.path.join(output_dir, 'genome_report.txt')

    total_seqs = count_sequences(merged_fasta)
    total_size = os.path.getsize(merged_fasta) if os.path.isfile(merged_fasta) else 0

    with safe_open(report_path, 'wt') as f:
        f.write("=" * 60 + "\n")
        f.write("  宿主参考基因组下载报告\n")
        f.write("=" * 60 + "\n")
        f.write(f"  生成时间: {datetime.now().isoformat()}\n")
        f.write(f"  合并文件: {merged_fasta}\n")
        f.write(f"  序列总数: {total_seqs}\n")
        f.write(f"  文件大小: {total_size / (1024*1024):.1f} MB\n")
        f.write(f"  源文件数: {len(fasta_files)}\n")
        f.write("\n  源文件列表:\n")
        for fa in fasta_files:
            f.write(f"    - {fa}\n")
        f.write("=" * 60 + "\n")

    UI.ok(f"报告已保存: {report_path}")


def _entrez_params(api_key=None, **extra):
    """E-utilities 参数（URL 由 requests params 序列化，不手工拼接）。"""
    p = dict(extra)
    if api_key:
        p['api_key'] = api_key
    return p


def _lookup_taxid(species: str, api_key=None):
    """用 NCBI esearch（taxonomy 库）查物种 taxid；失败返回 None。"""
    http = _http()
    term = f'"{species}"[Organism]'
    try:
        r = http.get(ESEARCH_URL,
                     params=_entrez_params(api_key, db='taxonomy', term=term,
                                           retmax=1, retmode='json'),
                     timeout=40)
        r.raise_for_status()
        ids = r.json().get('esearchresult', {}).get('idlist', [])
        return int(ids[0]) if ids else None
    except Exception:
        return None


def download_organelle_genome(species: str, organelle: str, out_dir: str,
                               ncbi_api: str = '') -> Optional[str]:
    """
    下载细胞器基因组（叶绿体/线粒体）

    参数:
        species: 物种名
        organelle: 'chloroplast' 或 'mitochondrion'
        out_dir: 输出目录
    """
    UI.info(f"正在下载{organelle}基因组...")

    org_dir = os.path.join(out_dir, organelle)
    os.makedirs(org_dir, exist_ok=True)

    if organelle == 'chloroplast':
        filter_word = 'chloroplast'
    else:
        filter_word = 'mitochondrion'
    search_term = f'"{species}"[Organism] AND {filter_word}[filter]'

    try:
        http = _http()
        resp = http.get(ESEARCH_URL,
                        params=_entrez_params(ncbi_api, db='nucleotide',
                                              term=search_term, retmode='json'),
                        timeout=30)
        data = resp.json()
        id_list = data.get('esearchresult', {}).get('idlist', [])

        if not id_list:
            UI.warn(f"  未找到 {species} 的 {organelle} 基因组")
            return None

        UI.ok(f"  找到 {len(id_list)} 条 {organelle} 记录")

        # 下载 FASTA
        ids_str = ','.join(id_list[:100])  # 限制100条
        resp = http.get(EFETCH_URL,
                        params=_entrez_params(ncbi_api, db='nucleotide',
                                              id=ids_str, rettype='fasta',
                                              retmode='text'),
                        timeout=120)
        fasta_content = resp.text

        out_name = f'{organelle}.fasta'
        out_path = os.path.join(org_dir, out_name)
        with safe_open(out_path, 'wt') as f:
            f.write(fasta_content)

        seq_count = fasta_content.count('>')
        UI.ok(f"  下载完成: {out_path} ({seq_count} 条序列)")
        return out_path

    except Exception as e:
        UI.err(f"  {organelle} 基因组下载失败: {e}")
        return None


def download_genome_https(species: str, out_dir: str, api_key=None,
                          max_assemblies: int = 5):
    """datasets CLI 不可用时的纯 Python 回退：E-utilities 搜 assembly →
    esummary 取 FTP 路径 → HTTPS 下载 RefSeq *_genomic.fna.gz 到 extracted/。"""
    http = _http()

    term = f'"{species}"[Organism] AND refseq[filter]'
    r = http.get(ESEARCH_URL,
                 params=_entrez_params(api_key, db='assembly', term=term,
                                       retmax=max_assemblies, retmode='json'),
                 timeout=40)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        UI.warn("E-utilities 未找到参考组装")
        return None
    UI.ok(f"找到 {len(ids)} 个组装: {', '.join(ids)}")

    r = http.get(ESUMMARY_URL,
                 params=_entrez_params(api_key, db='assembly',
                                       id=",".join(ids), retmode='json'),
                 timeout=40)
    r.raise_for_status()
    docs = r.json().get("result", {})

    extracted_dir = os.path.join(out_dir, "extracted")
    os.makedirs(extracted_dir, exist_ok=True)
    n_downloaded = 0
    for uid in ids:
        d = docs.get(uid) or {}
        ftp = d.get("ftppath_refseq") or d.get("ftppath_genbank") or ""
        if not ftp:
            continue
        # 仅接受 NCBI 官方 FTP 主机（换 HTTPS），_http() 内再做域名校验
        if not ftp.startswith(NCBI_FTP_HOST):
            UI.warn(f"跳过非 NCBI 路径: {ftp[:80]}")
            continue
        rel = ftp[len(NCBI_FTP_HOST):]
        dir_url = NCBI_FTP_HTTPS + rel
        try:
            lr = http.get(dir_url, timeout=60)
            lr.raise_for_status()
            fnames = re.findall(r'href="([^"]+_genomic\.fna\.gz)"', lr.text)
            if not fnames:
                continue
            fname = fnames[0]
            furl = dir_url + fname
            dest = os.path.join(extracted_dir, fname)
            if os.path.isfile(dest) and os.path.getsize(dest) > 1000:
                UI.ok(f"已存在: {fname}")
                n_downloaded += 1
                continue
            UI.info(f"下载: {furl}")
            with http.get(furl, stream=True, timeout=600) as resp:
                resp.raise_for_status()
                with safe_open(dest, 'wb') as f:
                    for chunk in resp.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
            UI.ok(f"完成: {fname}")
            n_downloaded += 1
        except Exception as e:
            UI.err(f"{dir_url} 下载失败: {e}")
    if not n_downloaded:
        UI.warn("未下载到任何基因组 FASTA")
        return None
    for fn in os.listdir(extracted_dir):
        if fn.endswith(".fna.gz"):
            src = os.path.join(extracted_dir, fn)
            dst_name = fn.replace('.fna.gz', '.fna')
            dst = os.path.join(extracted_dir, dst_name)
            if os.path.isfile(dst):
                continue
            with gzip.open(src, "rb") as fi, safe_open(dst, 'wb') as fo:
                fo.write(fi.read())
    UI.ok(f"HTTPS 回退完成: {n_downloaded} 个组装 → {extracted_dir}")
    return extracted_dir


def main():
    parser = argparse.ArgumentParser(
        description="📥 宿主参考基因组下载器 — NCBI datasets 封装",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 基础下载
  python host_genome.py --species "Lycium barbarum" --outdir ./host_genome

  # 包含细胞器基因组 + API key
  python host_genome.py --species "Lycium barbarum" \\
      --outdir ./host_genome --include-organelles --ncbi-api xxx

  # 仅验证/检查已有下载
  python host_genome.py --species "Lycium barbarum" \\
      --outdir ./host_genome --verify-only
        """
    )

    parser.add_argument('--species', required=True, help='物种拉丁学名')
    parser.add_argument('--outdir', default=None,
                        help='输出目录（默认 host-db/<taxid>_<物种>/）')
    parser.add_argument('--ncbi-api', help='NCBI API Key (提升速率)')
    parser.add_argument('--include-organelles', action='store_true',
                        help='同时下载叶绿体和线粒体基因组')
    parser.add_argument('--min-length', type=int, default=0,
                        help='最小序列长度过滤 (bp)')
    parser.add_argument('--verify-only', action='store_true',
                        help='仅验证现有下载')
    parser.add_argument('--skip-datasets', action='store_true',
                        help='跳过 NCBI datasets 步骤 (使用已有文件)')

    args = parser.parse_args()

    if not args.outdir:
        slug = re.sub(r'[^\w\-]+', '_', args.species).strip('_') or 'host'
        taxid = _lookup_taxid(args.species, args.ncbi_api)
        subdir = f'{taxid}_{slug}' if taxid else f'genome_{slug}'
        try:
            from vp.config import DIRS
            args.outdir = os.path.join(DIRS['host_src'], subdir)
        except Exception:
            args.outdir = subdir

    UI.header(f"宿主参考基因组下载: {args.species}")

    os.makedirs(args.outdir, exist_ok=True)

    # 验证模式
    if args.verify_only:
        UI.info("验证模式: 检查已有基因组文件...")
        merged_fasta = os.path.join(args.outdir, 'all.genome.uniq.fasta')
        if os.path.isfile(merged_fasta):
            n_seqs = count_sequences(merged_fasta)
            size_mb = os.path.getsize(merged_fasta) / (1024 * 1024)
            UI.ok(f"已存在合并基因组: {n_seqs} 条序列, {size_mb:.1f} MB")
        else:
            extracted_dir = os.path.join(args.outdir, 'extracted')
            fasta_files = collect_fasta_files(extracted_dir)
            if fasta_files:
                UI.ok(f"找到 {len(fasta_files)} 个源 FASTA 文件，但尚未合并")
                for f in fasta_files:
                    UI.info(f"  {f}")
            else:
                UI.warn("未找到任何基因组文件")
        return

    # 步骤 1: 使用 NCBI datasets 下载
    if not args.skip_datasets:
        if not check_datasets_cli():
            UI.info("改用 E-utilities + HTTPS 直接下载 RefSeq 基因组 ...")
            try:
                download_genome_https(args.species, args.outdir, args.ncbi_api)
            except Exception as e:
                UI.err(f"HTTPS 回退下载失败: {e}")
        else:
            genome_zip = os.path.join(args.outdir, 'genome_down.zip')

            if os.path.isfile(genome_zip) and os.path.getsize(genome_zip) > 1000:
                UI.ok(f"基因组压缩包已存在: {genome_zip}")
            else:
                UI.header("步骤 1/3: NCBI datasets 下载")
                cmd = [
                    'datasets', 'download', 'genome', 'taxon', args.species,
                    '--filename', genome_zip,
                    '--include', 'genome,gff3,seq-report'
                ]
                if args.ncbi_api:
                    cmd.extend(['--api-key', args.ncbi_api])

                UI.info(f"执行: {' '.join(cmd)}")
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

                if result.returncode != 0:
                    UI.err(f"datasets 下载失败: {result.stderr[:500]}")
                    UI.info("改用 E-utilities + HTTPS 直接下载 ...")
                    try:
                        download_genome_https(args.species, args.outdir, args.ncbi_api)
                    except Exception as e:
                        UI.err(f"HTTPS 回退下载失败: {e}")
                else:
                    UI.ok("基因组下载成功")

    # 步骤 2: 解压
    UI.header("步骤 2/3: 解压基因组")
    genome_zip = os.path.join(args.outdir, 'genome_down.zip')
    extracted_dir = os.path.join(args.outdir, 'extracted')

    if os.path.isfile(genome_zip):
        if not os.path.isdir(extracted_dir) or not os.listdir(extracted_dir):
            os.makedirs(extracted_dir, exist_ok=True)
            UI.info("正在解压...")
            result = subprocess.run(
                ['unzip', '-o', genome_zip, '-d', extracted_dir],
                capture_output=True, text=True
            )
            if result.returncode == 0:
                UI.ok("解压完成")
            else:
                UI.err(f"解压失败: {result.stderr[:200]}")
                # 尝试用 Python zipfile
                import zipfile
                try:
                    with zipfile.ZipFile(genome_zip, 'r') as zf:
                        zf.extractall(extracted_dir)
                    UI.ok("解压完成 (Python zipfile)")
                except Exception as e:
                    UI.err(f"Python 解压也失败: {e}")
        else:
            UI.ok("已解压，跳过")
    else:
        UI.warn(f"基因组压缩包不存在: {genome_zip}")

    # 步骤 3: 合并与去重
    UI.header("步骤 3/3: 合并与去重")

    fasta_files = collect_fasta_files(extracted_dir)
    if not fasta_files:
        UI.warn("未在解压目录中找到 FASTA 文件，尝试全目录搜索...")
        fasta_files = collect_fasta_files(args.outdir)

    merged_fasta = os.path.join(args.outdir, 'all.genome.uniq.fasta')
    if fasta_files:
        UI.ok(f"找到 {len(fasta_files)} 个 FASTA 文件:")
        for f in fasta_files:
            n_seqs = count_sequences(f)
            UI.info(f"  {os.path.basename(f)} — {n_seqs} 条序列")

        n_written, n_dups = merge_and_deduplicate(
            fasta_files, merged_fasta,
            min_length=args.min_length
        )
        generate_genome_report(args.outdir, fasta_files, merged_fasta)
    else:
        UI.warn("未找到任何 FASTA 文件")

    # 可选: 细胞器基因组（下载后重跑合并，把细胞器序列并入主文件）
    if args.include_organelles:
        UI.header("额外: 细胞器基因组下载")
        for org in ['chloroplast', 'mitochondrion']:
            org_fasta = download_organelle_genome(
                args.species, org, args.outdir, args.ncbi_api
            )
            if org_fasta and os.path.isfile(merged_fasta):
                UI.info(f"合并 {org} 基因组到主文件...")
                src_files = collect_fasta_files(extracted_dir)
                if org_fasta not in src_files:
                    src_files.append(org_fasta)
                merge_and_deduplicate(src_files, merged_fasta,
                                      min_length=args.min_length)
                UI.ok(f"{org} 基因组已合并")

    UI.header("下载完成")
    if os.path.isfile(merged_fasta):
        n_seqs = count_sequences(merged_fasta)
        size_mb = os.path.getsize(merged_fasta) / (1024 * 1024)
        UI.ok(f"最终输出: {merged_fasta}")
        UI.ok(f"  {n_seqs} 条序列, {size_mb:.1f} MB")
        UI.ok("下一步: 数据库页构建宿主库（python main.py build-host-db --genome ...）")
    else:
        UI.warn("最终合并文件未生成，请检查错误日志")


if __name__ == '__main__':
    main()

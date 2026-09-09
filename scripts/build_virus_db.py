#!/usr/bin/env python3
"""
Virus Database Builder
统一病毒数据库构建脚本
包含：环境依赖检查、分类学按需下载、软链接隔离、动静分离的数据库构建、表格化报告
"""

import os
import re
import sys
import time
import shutil
import logging
import argparse
import subprocess
import traceback
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass
from typing import List

try:
    from tqdm import tqdm
    import psutil
except ImportError:
    print("请安装必要的依赖: pip install tqdm psutil")
    sys.exit(1)

# ─────────────────────────── 日志配置 ───────────────────────────
def setup_logger(log_dir: Path) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"db_build_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logger = logging.getLogger("VirusDBBuilder")
    logger.setLevel(logging.DEBUG)
    
    # 文件：记录所有细节 (包含底层软件输出)
    file_fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(file_fmt)

    # 屏幕：显示 INFO 级别，配合 tqdm，保持动态但不刷屏
    console_fmt = logging.Formatter("[%(levelname)s] %(message)s")
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(console_fmt)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

# ─────────────────────────── 数据类 ───────────────────────────
@dataclass
class BuildStats:
    db_name: str
    start_time: float = 0.0
    end_time: float = 0.0
    peak_memory_mb: float = 0.0
    db_size_mb: float = 0.0
    success: bool = False
    error_msg: str = ""

    @property
    def elapsed_seconds(self) -> float:
        return self.end_time - self.start_time

# ─────────────────────────── 工具函数 ───────────────────────────
def run_cmd(cmd: str, logger: logging.Logger, check: bool = True) -> int:
    """执行命令：细节写入日志，屏幕仅显示正在执行的命令简影"""
    logger.debug(f"========== 执行命令 ==========\n{cmd}\n==============================")
    
    # 在屏幕上提示当前正在跑的核心命令（截取前80个字符，避免太长）
    short_cmd = cmd[:80] + "..." if len(cmd) > 80 else cmd
    tqdm.write(f"   [运行] {short_cmd}")
    
    result = subprocess.run(
        cmd, shell=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    
    if result.stdout:
        for line in result.stdout.splitlines():
            logger.debug(f" > {line}")
            
    if check and result.returncode != 0:
        logger.error(f"❌ 命令执行失败 (退出码: {result.returncode})")
        logger.error(f"详细报错信息请查看 log 文件。")
        raise RuntimeError(f"命令失败: {cmd}")
        
    return result.returncode

def get_dir_size_mb(path: Path) -> float:
    if not path.exists(): return 0.0
    total = sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) if path.is_dir() else path.stat().st_size
    return total / 1024 / 1024

def monitor_peak_memory(pid: int) -> float:
    try:
        proc = psutil.Process(pid)
        mem = proc.memory_info().rss
        for child in proc.children(recursive=True):
            try: mem += child.memory_info().rss
            except (psutil.NoSuchProcess, psutil.AccessDenied): pass
        return mem / 1024 / 1024
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return 0.0

def safe_symlink(src: Path, dst: Path, logger: logging.Logger):
    if dst.exists() or dst.is_symlink():
        logger.debug(f"符号链接已存在，跳过: {dst}")
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.symlink_to(src.resolve())
    logger.debug(f"创建链接: {dst} -> {src}")

# ─────────────────────────── 核心构建器 ───────────────────────────
class VirusDBBuilder:
    NUC_DBS = [
        "metabuli", "centrifuger", "kraken2", "bracken",
        "krakenuniq", "kmcp", "ganon", "sylph",
        "kun_peng", "blast", "kma", "salmon_kallisto", "lexicmap"
    ]
    PROT_DBS = ["kaiju", "kraken2x", "diamond", "mmseqs", "CAT"]
    ALL_DATABASES = NUC_DBS + PROT_DBS

    def __init__(self, args: argparse.Namespace, logger: logging.Logger):
        self.args = args
        self.logger = logger
        self.stats: list[BuildStats] = []

        if not args.nuc_fasta and not args.prot_fasta:
            self.logger.error("必须至少提供一个核苷酸序列(--nuc-fasta)或蛋白质序列(--prot-fasta)！")
            sys.exit(1)

        self.nuc_fasta  = Path(args.nuc_fasta).resolve() if args.nuc_fasta else None
        self.prot_fasta = Path(args.prot_fasta).resolve() if args.prot_fasta else None

        self.work_dir     = Path(args.work_dir).resolve()
        self.db_prefix    = args.db_prefix
        self.tax_dir      = Path(args.tax_dir).resolve() if args.tax_dir else self.work_dir / "taxonomy"
        
        self.genome_dir   = self.work_dir / f"{self.db_prefix}.genome"
        self.genome_list  = self.work_dir / "genomes.txt"

        self.nucl_gb_acc2tax  = self.tax_dir / "nucl_gb.accession2taxid"
        #self.nucl_wgs_acc2tax = self.tax_dir / "nucl_wgs.accession2taxid"
        self.prot_acc2tax     = self.tax_dir / "prot.accession2taxid"

        self.seqid2taxid     = self.work_dir / "seqid2taxid.map"
        self.protid2taxid    = self.work_dir / "protid2taxid.map"
        self.prot_is_pyrodigal = False  # 标记蛋白序列是否为 pyrodigal 基因预测格式
        self.seqid2taxid_len = self.work_dir / "seqid2taxid_len.map"
        self.protid2taxid_len= self.work_dir / "protid2taxid_len.map"
        self.combined_map    = self.work_dir / f"{self.db_prefix}.id2taxid.map"
        self.nuc_tax_txt     = self.work_dir / f"{self.db_prefix}.nuc.taxonomy.txt"
        self.prot_tax_txt    = self.work_dir / f"{self.db_prefix}.prot.taxonomy.txt"
        self.combined_tax    = self.work_dir / f"{self.db_prefix}_taxonomy.txt"

        self.threads = args.threads
        self.max_ram = args.max_ram
        self.jellyfish_bin = args.jellyfish_bin

        self.selected = self._parse_target_databases(args.databases)

        self.logger.info(f"工作目录 : {self.work_dir}")
        self.logger.info(f"分类学目录: {self.tax_dir}")
        self.logger.info(f"数据库前缀: {self.db_prefix}")
        if self.nuc_fasta: self.logger.info(f"核酸FASTA : {self.nuc_fasta}")
        if self.prot_fasta: self.logger.info(f"蛋白FASTA : {self.prot_fasta}")
        self.logger.info(f"计划构建库: {', '.join(self.selected)}")

    def _parse_target_databases(self, requested: List[str]) -> List[str]:
        selected = []
        if "all" in requested:
            if self.nuc_fasta: selected.extend(self.NUC_DBS)
            if self.prot_fasta: selected.extend(self.PROT_DBS)
        else:
            for db in requested:
                if db in self.NUC_DBS and not self.nuc_fasta:
                    self.logger.warning(f"⚠️ 构建 {db} 需要核酸序列，已忽略。")
                elif db in self.PROT_DBS and not self.prot_fasta:
                    self.logger.warning(f"⚠️ 构建 {db} 需要蛋白序列，已忽略。")
                else:
                    selected.append(db)
        return list(dict.fromkeys(selected))

    def check_dependencies(self):
        self.logger.info("-" * 60)
        self.logger.info("【预检】 验证系统环境与必要工具")
        
        global_tools = ["aria2c", "rapidgzip", "seqkit", "awk", "sed"]
        custom_scripts = ["SearchAccessionIdToTaxId.py", "db_seqid2taxid_add_legth.py", "make_ktaxonomy.py"]

        missing_tools = []
        for tool in global_tools + custom_scripts:
            if shutil.which(tool) is None: missing_tools.append(tool)

        if missing_tools:
            self.logger.error("❌ 以下工具或脚本缺失，请确认已安装并在系统 PATH 中:")
            for m in missing_tools: self.logger.error(f"  - {m}")
            sys.exit(1)
        self.logger.info("✅ 依赖检查全部通过")

    # ────────────────────────────────────────────────────────────
    # Step 1: 下载与处理 Taxonomy
    # ────────────────────────────────────────────────────────────
    def step1_taxonomy(self):
        self.logger.info("-" * 60)
        self.logger.info("【步骤 1】 下载和处理 Taxonomy 文件 (严格按需处理)")
        
        self.tax_dir.mkdir(parents=True, exist_ok=True)
        
        # 处理 taxdump
        taxdump = self.work_dir / "new_taxdump.tar.gz"
        nodes = self.tax_dir / "nodes.dmp"
        if not nodes.exists():
            if not taxdump.exists():
                tqdm.write("⏳ 正在下载 NCBI new_taxdump.tar.gz ...")
                run_cmd(f"aria2c -x 16 -s 16 -d {self.work_dir} ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz", self.logger)
            tqdm.write("⏳ 正在解压 new_taxdump.tar.gz ...")
            run_cmd(f"tar zxf {taxdump} -C {self.tax_dir}", self.logger)
        else:
            self.logger.info("✅ Taxonomy (nodes.dmp) 已存在，跳过下载与解压")

        # 动态分配需要下载的文件
        files_to_download = {}
        if self.nuc_fasta:
            files_to_download["nucl_gb.accession2taxid"] = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
            #files_to_download["nucl_wgs.accession2taxid"] = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz"
        
        if self.prot_fasta:
            files_to_download["prot.accession2taxid"] = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
        
        for out_name, url in files_to_download.items():
            out_path = self.tax_dir / out_name
            gz_name = url.split('/')[-1]
            gz_path = self.work_dir / gz_name
            tax_gz_path = self.tax_dir / gz_name  # 也检查 taxonomy 目录下的 .gz

            # 优先检查解压后的文件是否存在
            if out_path.exists():
                self.logger.info(f"✅ 解压完成的 {out_name} 已存在，跳过处理")
                continue

            # 检查压缩包: 先看 taxonomy 目录, 再看 work_dir
            actual_gz = None
            if tax_gz_path.exists():
                actual_gz = tax_gz_path
            elif gz_path.exists():
                actual_gz = gz_path

            if actual_gz:
                tqdm.write(f"⏳ 检测到已存在的 {gz_name}，直接解压...")
            else:
                tqdm.write(f"⏳ 正在下载 {gz_name} (文件极大，后台下载中)...")
                run_cmd(f"aria2c -x 16 -s 16 -d {self.work_dir} {url}", self.logger)
                actual_gz = gz_path

            tqdm.write(f"⏳ 正在解压 {gz_name} 为 {out_name} (高磁盘IO操作，请稍等)...")
            run_cmd(f"rapidgzip -d -P 0 -k {actual_gz} -o {out_path}", self.logger)
            self.logger.info(f"✅ {out_name} 处理完成！")

    # ────────────────────────────────────────────────────────────
    # Step 2: 参考基因组规范化处理
    # ────────────────────────────────────────────────────────────
    def step2_genomes(self):
        self.logger.info("-" * 60)
        self.logger.info("【步骤 2】 参考基因组规范化处理")

        if self.nuc_fasta:
            if not (self.genome_list.exists() and self.genome_dir.exists()):
                tqdm.write("⏳ 正在拆分核酸 FASTA 序列...")
                self.genome_dir.mkdir(parents=True, exist_ok=True)
                run_cmd(f"seqkit split2 -j {self.threads} -s 1 --quiet -N {self.nuc_fasta} -O {self.genome_dir}", self.logger)
                run_cmd(f"find {self.genome_dir} -name '*.fasta' > {self.genome_list}", self.logger)
            else:
                self.logger.info("✅ 核酸拆分文件已存在，跳过拆分")

            nuc_id_file = self.work_dir / f"{self.db_prefix}.nuc.id"
            if not self.seqid2taxid.exists():
                tqdm.write("⏳ 正在生成核酸 SeqID 映射...")
                run_cmd(f"grep '>' {self.nuc_fasta} | sed 's/>//' > {nuc_id_file}", self.logger)
                run_cmd(f"SearchAccessionIdToTaxId.py -j {self.threads} -b 100000 --query-file {nuc_id_file} --map-file {self.nucl_gb_acc2tax} -o {self.seqid2taxid}", self.logger)

            if not self.seqid2taxid_len.exists():
                tqdm.write("⏳ 正在生成核酸序列长度映射...")
                run_cmd(f"db_seqid2taxid_add_legth.py --input {self.seqid2taxid} --genome-file {self.nuc_fasta} --processes {self.threads} --type nucleotide --output {self.seqid2taxid_len}", self.logger)

        if self.prot_fasta:
            prot_id_file = self.work_dir / f"{self.db_prefix}.prot.id"
            prot_id_nucl = self.work_dir / f"{self.db_prefix}.prot.nucl.id"  # 提取的核酸 Accession
            # 检测蛋白 ID 格式: pyrodigal 产出的 >NC_116488.1_1 # 77 # ...
            use_nucl_map = False
            if not self.protid2taxid.exists():
                tqdm.write("⏳ 正在生成蛋白 ProtID 映射...")
                run_cmd(f"grep '>' {self.prot_fasta} | sed 's/>//;s/ .*//' > {prot_id_file}", self.logger)
                first_id = open(prot_id_file).readline().strip()
                if '#' in first_id or ('_' in first_id and '.' in first_id):
                    tqdm.write("   🔍 检测到基因预测格式蛋白 ID, 提取核酸 Accession 并用 nucl_gb 映射...")
                    # 提取 NC_116488.1_1 → NC_116488.1 (去掉最后的 _数字)
                    awk_nucl = r"""sed 's/_[0-9]*$//' """
                    run_cmd(f"{awk_nucl} {prot_id_file} | sort -u > {prot_id_nucl}", self.logger)
                    # protid2taxid: 蛋白完整ID → TaxID (下游 diamond/mmseqs 需要)
                    # 策略: 从 nucl_gb 查 TaxID, 然后每个蛋白ID继承其核酸 Accession 的 TaxID
                    tmp_map = self.work_dir / f"{self.db_prefix}.prot.nucl2taxid.map"
                    run_cmd(f"SearchAccessionIdToTaxId.py -j {self.threads} -b 100000 --query-file {prot_id_nucl} --map-file {self.nucl_gb_acc2tax} -o {tmp_map}", self.logger)
                    # 生成 prot ID → TaxID 映射: 每个蛋白ID从核酸 Accession 继承 TaxID
                    awk_join = r"""awk -F'\t' 'NR==FNR {tax[$1]=$2; next} {acc=$0; sub(/_[0-9]*$/, "", acc); if(acc in tax) print $0 "\t" tax[acc]}' """
                    run_cmd(f"{awk_join} {tmp_map} {prot_id_file} > {self.protid2taxid}", self.logger)
                    use_nucl_map = True
                    self.prot_is_pyrodigal = True
                else:
                    self.logger.info("检测到常规蛋白 ID 格式, 使用 prot.accession2taxid 映射")
                    acc2tax_map = self.prot_acc2tax
                    run_cmd(f"SearchAccessionIdToTaxId.py -j {self.threads} -b 100000 --query-file {prot_id_file} --map-file {acc2tax_map} -o {self.protid2taxid}", self.logger)

            if not self.protid2taxid_len.exists():
                tqdm.write("⏳ 正在生成蛋白序列长度映射...")
                run_cmd(f"db_seqid2taxid_add_legth.py --input {self.protid2taxid} --genome-file {self.prot_fasta} --processes {self.threads} --type protein --output {self.protid2taxid_len}", self.logger)

        nodes = self.tax_dir / "nodes.dmp"
        names = self.tax_dir / "names.dmp"
        
        tqdm.write("⏳ 正在生成供下游软件使用的 Taxonomy TXT 文件...")
        if self.nuc_fasta and not self.nuc_tax_txt.exists():
            run_cmd(f"make_ktaxonomy.py --nodes {nodes} --names {names} --seqid2taxid {self.seqid2taxid} -o {self.nuc_tax_txt}", self.logger)
        if self.prot_fasta and not self.prot_tax_txt.exists():
            run_cmd(f"make_ktaxonomy.py --nodes {nodes} --names {names} --seqid2taxid {self.protid2taxid} -o {self.prot_tax_txt}", self.logger)

        if not self.combined_map.exists():
            if self.nuc_fasta and self.prot_fasta: run_cmd(f"cat {self.seqid2taxid} {self.protid2taxid} > {self.combined_map}", self.logger)
            elif self.nuc_fasta: shutil.copy(self.seqid2taxid, self.combined_map)
            elif self.prot_fasta: shutil.copy(self.protid2taxid, self.combined_map)

        if not self.combined_tax.exists():
            run_cmd(f"make_ktaxonomy.py --nodes {nodes} --names {names} --seqid2taxid {self.combined_map} -o {self.combined_tax}", self.logger)
        
        self.logger.info("✅ 基因组规范化及映射表生成全部完成！")

    # ────────────────────────────────────────────────────────────
    # Step 3: 数据库构建调度机制 (配 tqdm 优雅显示)
    # ────────────────────────────────────────────────────────────
    def _is_pyrodigal_format(self) -> bool:
        """检测蛋白 FASTA 是否来自 pyrodigal 基因预测 (ID 含 _基因编号 或 # 分隔符)"""
        if not self.prot_fasta or not self.prot_fasta.exists():
            return False
        with open(self.prot_fasta) as f:
            first = f.readline().strip()
        # pyrodigal 格式: >NC_116488.1_1 (含 _数字后缀) 或原始 >NC_116488.1_1 # 77 # ...
        return '#' in first or bool(re.search(r'_\d+$', first.split()[0].lstrip('>')))

    def _build_with_stats(self, db_name: str, build_func) -> BuildStats:
        stat = BuildStats(db_name=db_name)
        stat.start_time = time.time()
        pid = os.getpid()

        try:
            build_func()
            stat.success = True
            tqdm.write(f"   ✅ {db_name} 数据库构建完毕。")
        except Exception as e:
            stat.success = False
            stat.error_msg = str(e)
            tqdm.write(f"   ❌ 构建 {db_name} 失败: {e} (详情见log)")
            self.logger.debug(traceback.format_exc())
        finally:
            stat.end_time = time.time()
            stat.peak_memory_mb = monitor_peak_memory(pid)
            db_dir = self.work_dir / f"{self.db_prefix}.{db_name}_db"
            stat.db_size_mb = get_dir_size_mb(db_dir)

        return stat

    def step3_build(self):
        self.logger.info("-" * 60)
        self.logger.info("【步骤 3】 数据库正式构建")

        if not self.selected:
            self.logger.warning("没有需要构建的数据库任务。")
            return

        build_map = {
            "metabuli": self.build_metabuli, "centrifuger": self.build_centrifuger,
            "kraken2": self.build_kraken2, "bracken": self.build_bracken,
            "krakenuniq": self.build_krakenuniq, "kmcp": self.build_kmcp,
            "ganon": self.build_ganon, "sylph": self.build_sylph,
            "kun_peng": self.build_kun_peng, "blast": self.build_blast,
            "kma": self.build_kma, "salmon_kallisto": self.build_salmon_kallisto,
            "lexicmap": self.build_lexicmap, "kaiju": self.build_kaiju,
            "kraken2x": self.build_kraken2x, "diamond": self.build_diamond, "mmseqs": self.build_mmseqs, "CAT": self.build_cat,
        }

        # 用 tqdm 构建时先检查哪些已完成, 支持断点续传
        skipped = []
        remaining = []
        for db_name in self.selected:
            db_dir = self.work_dir / f"{self.db_prefix}.{db_name}_db"
            if db_dir.exists() and get_dir_size_mb(db_dir) > 1:
                skipped.append(db_name)
            else:
                remaining.append(db_name)

        if skipped:
            tqdm.write(f"   ⊙ 断点续传: 以下数据库已构建, 跳过: {', '.join(skipped)}")
            for db_name in skipped:
                stat = BuildStats(db_name=db_name)
                stat.success = True
                stat.db_size_mb = get_dir_size_mb(self.work_dir / f"{self.db_prefix}.{db_name}_db")
                self.stats.append(stat)

        if not remaining:
            tqdm.write("   ✅ 全部数据库已构建完成!")
            return

        with tqdm(remaining, desc="建库总进度", unit="db", bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]") as pbar:
            for db_name in pbar:
                pbar.set_postfix_str(f"当前任务: {db_name}")
                tqdm.write(f"\n🚀 开始构建数据库: {db_name}")
                if db_name in build_map:
                    stat = self._build_with_stats(db_name, build_map[db_name])
                    self.stats.append(stat)

    # ────────────────────────────────────────────────────────────
    # 各数据库构建命令
    # ────────────────────────────────────────────────────────────
    def build_metabuli(self):
        db_dir = self.work_dir / f"{self.db_prefix}.metabuli_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        run_cmd(f"metabuli build {db_dir} {self.genome_list} {self.nucl_gb_acc2tax} --taxonomy-path {self.tax_dir} --max-ram {self.max_ram} --threads {self.threads}", self.logger)

    def build_centrifuger(self):
        db_dir = self.work_dir / f"{self.db_prefix}.centrifuger_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        prefix = db_dir / f"{self.db_prefix}.centrifuger"
        run_cmd(f"centrifuger-build -t {self.threads} --conversion-table {self.seqid2taxid} -l {self.genome_list} --taxonomy-tree {self.tax_dir}/nodes.dmp --name-table {self.tax_dir}/names.dmp -o {prefix}", self.logger)

    def build_kraken2(self):
        db_dir = self.work_dir / f"{self.db_prefix}.kraken2_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        safe_symlink(self.tax_dir, db_dir / "taxonomy", self.logger)
        
        # 【重点修复】废除软链接，直接让 kraken2-build 使用 --add-to-library 导入。
        # 这样 kraken2 内部生成的 added/ 文件夹会放在 db_dir/library/ 下，绝不会污染原 genome 目录。
        run_cmd(f"kraken2-build --add-to-library {self.nuc_fasta} --db {db_dir} --threads {self.threads}", self.logger)
        run_cmd(f"kraken2-build --build --db {db_dir} --threads {self.threads}", self.logger)

    def build_bracken(self):
        kraken2_db = self.work_dir / f"{self.db_prefix}.kraken2_db"
        if not kraken2_db.exists(): raise RuntimeError("bracken 需要先构建 kraken2 数据库")
        for read_len in [50, 75, 100, 150]:
            run_cmd(f"bracken-build -d {kraken2_db} -t {self.threads} -k 31 -l {read_len}", self.logger)

    def build_krakenuniq(self):
        db_dir = self.work_dir / f"{self.db_prefix}.krakenuniq_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        safe_symlink(self.tax_dir, db_dir / "taxonomy", self.logger)
        
        # 【重点修复】为 krakenuniq 隔离污染：新建一个干净的 library 文件夹
        # 并只把 .fasta 文件单独软链接进去，隔离 prelim_map.txt 的污染
        lib_dir = db_dir / "library"
        lib_dir.mkdir(exist_ok=True)
        for fasta_file in self.genome_dir.glob("*.fasta"):
            safe_symlink(fasta_file, lib_dir / fasta_file.name, self.logger)

        run_cmd(f"krakenuniq-build --db {db_dir} --library-dir {lib_dir} --threads {self.threads} --jellyfish-bin {self.jellyfish_bin}", self.logger)

    def build_kmcp(self):
        kmcp_genome_dir = self.work_dir / f"{self.db_prefix}.kmcp.genome"
        kmcp_genome_dir.mkdir(parents=True, exist_ok=True)
        kmcp_list = self.work_dir / "genomes.kmcp.txt"
        if not any(kmcp_genome_dir.glob("*.fasta.kmcp")):
            run_cmd(f"seqkit split2 -j {self.threads} -s 1 -N {self.nuc_fasta} -e .kmcp -O {kmcp_genome_dir}", self.logger)
        run_cmd(f"find {kmcp_genome_dir} -name '*.fasta.kmcp' > {kmcp_list}", self.logger)
        
        db_dir  = self.work_dir / f"{self.db_prefix}.kmcp_db"
        db_dir2 = self.work_dir / f"{self.db_prefix}.kmcp_db2"
        db_dir.mkdir(parents=True, exist_ok=True)
        run_cmd(f"kmcp compute --infile-list {kmcp_list} --threads {self.threads} --out-dir {db_dir}", self.logger)
        run_cmd(f"kmcp index --in-dir {db_dir} --threads {self.threads} --num-hash 1 --false-positive-rate 0.3 --out-dir {db_dir2}", self.logger)

    def build_ganon(self):
        db_dir = self.work_dir / f"{self.db_prefix}.ganon_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        ganon_map = self.work_dir / "ganon_table.map"
        if not ganon_map.exists():
            awk_script = r"""awk -F'\t' -v OFS='\t' 'NR==FNR { map[$1]=$2; next } { filepath=$0; seqid=filepath; sub(/.*\//, "", seqid); sub(/\.fasta$/, "", seqid); taxid=(seqid in map)?map[seqid]:"NA"; print filepath, seqid, taxid }' """
            run_cmd(f"{awk_script} {self.seqid2taxid} {self.genome_list} > {ganon_map}", self.logger)
        prefix = db_dir / f"{self.db_prefix}.ganon"
        run_cmd(f"ganon build-custom --input-file {ganon_map} --taxonomy-files {self.tax_dir}/nodes.dmp {self.tax_dir}/names.dmp --threads {self.threads} --skip-genome-size --db-prefix {prefix}", self.logger)

    def build_sylph(self):
        db_dir = self.work_dir / f"{self.db_prefix}.sylph_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        for c in [30, 50, 100]:
            out = db_dir / f"{self.db_prefix}.c{c}.sylph"
            if not out.exists():
                run_cmd(f"sylph sketch -c {c} -i {self.nuc_fasta} -o {out}", self.logger)

    def build_kun_peng(self):
        kraken2_db = self.work_dir / f"{self.db_prefix}.kraken2_db"
        if not kraken2_db.exists(): raise RuntimeError("kun_peng 需要先构建 kraken2 数据库")
        run_cmd(f"kun_peng hashshard --db {kraken2_db} --hash-capacity 1G", self.logger)
        db_dir = self.work_dir / f"{self.db_prefix}.kun_peng_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        for f in kraken2_db.glob("hash_*k2d"): shutil.move(str(f), db_dir / f.name)
        for f in kraken2_db.glob("*.k2d"):
            if not (db_dir / f.name).exists(): shutil.copy(str(f), db_dir / f.name)

    def build_blast(self):
        db_dir = self.work_dir / f"{self.db_prefix}.blast_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = db_dir / f"{self.db_prefix}.blast"
        run_cmd(f"makeblastdb -in {self.nuc_fasta} -dbtype nucl -out {out_prefix} -parse_seqids -taxid_map {self.nucl_gb_acc2tax}", self.logger)

    def build_kma(self):
        db_dir = self.work_dir / f"{self.db_prefix}.kma_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = db_dir / f"{self.db_prefix}.kma"
        run_cmd(f"kma index -i {self.nuc_fasta} -o {out_prefix}", self.logger)

    def build_salmon_kallisto(self):
        db_dir = self.work_dir / f"{self.db_prefix}.index_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        sal_idx = db_dir / f"{self.db_prefix}.sal_index"
        kal_idx = db_dir / f"{self.db_prefix}.kal_index"
        if not sal_idx.exists(): run_cmd(f"salmon index -k 31 -i {sal_idx} -t {self.nuc_fasta} --threads {self.threads}", self.logger)
        if not kal_idx.exists(): run_cmd(f"kallisto index -k 31 -i {kal_idx} {self.nuc_fasta} --threads {self.threads}", self.logger)

    def build_lexicmap(self):
        lmi_dir = self.work_dir / f"{self.db_prefix}.lmi"
        db_dir = self.work_dir / f"{self.db_prefix}.lexicmap_db"
        if not db_dir.exists():
            run_cmd(f"lexicmap index -I {self.genome_dir} -O {lmi_dir} --batch-size 5000 --threads {self.threads}", self.logger)
            lmi_dir.rename(db_dir)

    def build_kaiju(self):
        db_dir = self.work_dir / f"{self.db_prefix}.kaiju_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        safe_symlink(self.tax_dir, db_dir / "taxonomy", self.logger)
        kaiju_fasta = self.work_dir / f"{self.db_prefix}.prot.kaiju.fasta"
        if not kaiju_fasta.exists():
            run_cmd(f"awk 'NR==FNR {{if(NR>0) map[$1]=$2; next}} /^>/ {{id=substr($0,2); if(id in map) print \">\" id \"_\" map[id]; else print $0; next}} {{print}}' {self.protid2taxid} {self.prot_fasta} > {kaiju_fasta}", self.logger)
        kaiju_prefix = self.work_dir / f"{self.db_prefix}.kaiju"
        run_cmd(f"kaiju-mkbwt -n {self.threads} -a protein -o {kaiju_prefix} {kaiju_fasta}", self.logger)
        run_cmd(f"kaiju-mkfmi {kaiju_prefix}", self.logger)
        for ext in [".sa", ".bwt"]:
            tmp = Path(f"{kaiju_prefix}{ext}")
            if tmp.exists(): tmp.unlink()
        kaiju_fasta.unlink(missing_ok=True)
        fmi = Path(f"{kaiju_prefix}.fmi")
        if fmi.exists(): shutil.move(str(fmi), db_dir / fmi.name)

    def build_kraken2x(self):
        db_dir = self.work_dir / f"{self.db_prefix}.kraken2x_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        safe_symlink(self.tax_dir, db_dir / "taxonomy", self.logger)
        run_cmd(f"kraken2-build --add-to-library {self.prot_fasta} --db {db_dir} --threads {self.threads} --protein", self.logger)
        run_cmd(f"kraken2-build --build --db {db_dir} --threads {self.threads} --protein", self.logger)

    def build_diamond(self):
        db_dir = self.work_dir / f"{self.db_prefix}.diamond_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        dmnd = db_dir / f"{self.db_prefix}.dmnd"
        # pyrodigal 格式: protid2taxid.map 已是 蛋白ID→TaxID, 需加表头给 diamond
        if self._is_pyrodigal_format():
            diamond_taxmap = self.work_dir / f"{self.db_prefix}.diamond.taxmap"
            if not diamond_taxmap.exists():
                run_cmd(f"{{ echo -e 'accession.version\\ttaxid'; cat {self.protid2taxid}; }} > {diamond_taxmap}", self.logger)
            taxonmap = diamond_taxmap
        else:
            taxonmap = self.prot_acc2tax
        run_cmd(f"diamond makedb --threads {self.threads} --in {self.prot_fasta} --db {dmnd} --taxonmap {taxonmap} --taxonnodes {self.tax_dir}/nodes.dmp --taxonnames {self.tax_dir}/names.dmp", self.logger)

    def build_mmseqs(self):
        db_dir = self.work_dir / f"{self.db_prefix}.mmseqs_db"
        db_dir.mkdir(parents=True, exist_ok=True)
        prot_taxmap = self.work_dir / "prot.tax-map"
        # pyrodigal 格式: protid2taxid.map 已是 蛋白ID→TaxID, 直接复用
        if self._is_pyrodigal_format():
            if not prot_taxmap.exists():
                shutil.copy(self.protid2taxid, prot_taxmap)
        else:
            prot_id_file = self.work_dir / f"{self.db_prefix}.prot.id"
            if not prot_id_file.exists():
                run_cmd(f"grep '>' {self.prot_fasta} | sed 's/>//' > {prot_id_file}", self.logger)
            if not prot_taxmap.exists():
                run_cmd(f"awk 'NR==FNR {{a[$1]; next}} $1 in a {{print $1 \"\\t\" $3}}' {prot_id_file} {self.prot_acc2tax} > {prot_taxmap}", self.logger)

        mmseqs_prefix = db_dir / f"{self.db_prefix}.mmseqs"
        tmp_dir = self.work_dir / "mmseqs_tmp"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        run_cmd(f"mmseqs createdb --dbtype 1 {self.prot_fasta} {mmseqs_prefix}", self.logger)
        run_cmd(f"mmseqs createtaxdb {mmseqs_prefix} {tmp_dir} --ncbi-tax-dump {self.tax_dir} --tax-mapping-file {prot_taxmap} --threads {self.threads}", self.logger)

    # ────────────────────────────────────────────────────────────
    # 表格化 Markdown 总结报告
    # ────────────────────────────────────────────────────────────
    def build_cat(self):
        db_dir = self.work_dir / f"{self.db_prefix}.CAT-db"
        db_dir.mkdir(parents=True, exist_ok=True)
        # CAT 需要 accession.version\ttaxid 格式的头行
        cat_taxmap = self.work_dir / f"{self.db_prefix}.cat.taxmap"
        if not cat_taxmap.exists():
            run_cmd(f"sed '1i\\\accession.version\\ttaxid' {self.protid2taxid} > {cat_taxmap}", self.logger)
        run_cmd(f"CAT_pack prepare --db_fasta {self.prot_fasta} --names {self.tax_dir}/names.dmp --nodes {self.tax_dir}/nodes.dmp --acc2tax {cat_taxmap} --db_dir {db_dir}", self.logger)

    def print_report(self):
        self.logger.info("\n")
        self.logger.info("+" + "-"*76 + "+")
        self.logger.info(f"|{'构建结果汇总 (Summary Report)':^76}|")
        self.logger.info("+" + "-"*76 + "+")
        
        header = f"| {'数据库 (Database)':<20} | {'状态 (Status)':<12} | {'耗时 (Time)':<12} | {'磁盘 (Size)':<10} | {'峰值内存':<10} |"
        self.logger.info(header)
        self.logger.info("|" + "-"*22 + "+" + "-"*14 + "+" + "-"*14 + "+" + "-"*12 + "+" + "-"*12 + "|")

        total_time = sum(s.elapsed_seconds for s in self.stats)
        total_size = sum(s.db_size_mb for s in self.stats)
        success_n  = sum(1 for s in self.stats if s.success)

        for s in self.stats:
            status = "✅ 成功" if s.success else "❌ 失败"
            time_str = f"{s.elapsed_seconds:.1f}s"
            ram_str = f"{s.peak_memory_mb:.1f}MB"
            size_str = f"{s.db_size_mb:.1f}MB"
            row = f"| {s.db_name:<20} | {status:<11} | {time_str:<12} | {size_str:<10} | {ram_str:<10} |"
            self.logger.info(row)

        self.logger.info("+" + "-"*76 + "+")
        self.logger.info(f"| 统计: 成功 {success_n}/{len(self.stats)} | 总耗时: {total_time/60:.1f} min | 总磁盘占用: {total_size/1024:.2f} GB{'':<18}|")
        self.logger.info("+" + "-"*76 + "+")


# ─────────────────────────── CLI ───────────────────────────
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""
======================================================================
  🦠 病毒数据库统一构建工具 (Virus Database Builder)
======================================================================
此脚本提供一键式的自动化流程，将病毒核酸/蛋白质序列构建为多种常用
生物信息学软件的本地数据库。

【核心特性】
  1. 智能按需下载：仅在检测到对应的fasta输入时，才会下载对应的Taxonomy映射表。
  2. 软链接隔离：彻底解决 Kraken 等软件污染原 Genome 目录的问题。
  3. 日志分离：终端展示极简进度与核心操作，所有底层报错与刷屏日志写入 logs/ 文件夹。
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
【使用示例】
1. 仅构建核酸数据库 (指定特定库，并使用本地已有的 Taxonomy 目录)：
   python build_virus_db.py \\
     --nuc-fasta final.cluster.ref.fasta \\
     --db-prefix virus.ref \\
     --tax-dir ~/database/taxonomy/ \\
     --databases metabuli centrifuger kraken2

2. 同时构建所有核酸和蛋白数据库 (使用逗号或空格分隔均可)：
   python build_virus_db.py \\
     --nuc-fasta nuc.fasta --prot-fasta prot.fasta \\
     --db-prefix my_db \\
     --databases all

3. 构建 krakenuniq 并指定 jellyfish 路径：
   python build_virus_db.py \\
     --nuc-fasta nuc.fasta \\
     --databases krakenuniq \\
     --jellyfish-bin ~/.pixi/envs/krakenuniq/bin/jellyfish

4. 中断后恢复 (跳过前两个耗时阶段)：
   python build_virus_db.py \\
     --nuc-fasta ref.fasta \\
     --skip-step1 --skip-step2 \\
     --databases kraken2 blast
"""
    )

    parser.add_argument("--nuc-fasta",  default=None, help="[输入] 核酸序列 FASTA 文件 (构建核酸库时必填)")
    parser.add_argument("--prot-fasta", default=None, help="[输入] 蛋白序列 FASTA 文件 (构建蛋白库时必填)")
    parser.add_argument("--db-prefix",  default="virus_db", help="[输出] 数据库前缀，默认: virus_db")
    parser.add_argument("--work-dir",   default=".", help="[输出] 工作目录，生成文件位于此目录")
    parser.add_argument("--tax-dir",    default=None, help="[输入] 指定已有的 Taxonomy 目录 (如 ~/database/taxonomy/)")

    parser.add_argument(
        "--databases", nargs="+", default=["all"], metavar="DB",
        help=f"指定要构建的数据库(逗号或空格分隔)。输入 'all' 时自动判定。\n包含: {', '.join(VirusDBBuilder.ALL_DATABASES)}"
    )

    parser.add_argument("--threads",  type=int, default=32, help="最大线程数 (默认: 32)")
    parser.add_argument("--max-ram",  type=int, default=500, help="最大内存限制(GB) (默认: 500)")
    
    # 新增 krakenuniq 所需的 jellyfish 路径参数
    parser.add_argument("--jellyfish-bin", default="jellyfish", help="[高级] 指定 jellyfish 的绝对路径 (用于 krakenuniq 构建)")

    parser.add_argument("--skip-step1", action="store_true", help="跳过 Taxonomy 准备。前提: files(解压后) 已存在")
    parser.add_argument("--skip-step2", action="store_true", help="跳过 序列拆分和映射。前提: genome 及 map 已存在")

    args = parser.parse_args()

    raw_dbs = args.databases
    parsed_dbs = []
    for item in raw_dbs:
        parsed_dbs.extend([x.strip() for x in item.split(',') if x.strip()])
    
    valid_choices = set(VirusDBBuilder.ALL_DATABASES + ["all"])
    for db in parsed_dbs:
        if db not in valid_choices:
            parser.error(f"无效的数据库选项: '{db}'。\n请使用 -h 查看支持的数据库列表。")
    
    args.databases = parsed_dbs
    return args

def main():
    args   = parse_args()
    work   = Path(args.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(work / "logs")

    logger.info("=" * 60)
    logger.info("  🦠 病毒数据库统一构建工具 (Virus DB Builder)")
    logger.info(f"  启动时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 60)

    builder = VirusDBBuilder(args, logger)
    builder.check_dependencies()

    if not args.skip_step1: builder.step1_taxonomy()
    else: logger.info(">> 接收到 --skip-step1 指令，跳过阶段 1: Taxonomy 处理")

    if not args.skip_step2: builder.step2_genomes()
    else: logger.info(">> 接收到 --skip-step2 指令，跳过阶段 2: 基因组规范化")

    builder.step3_build()
    builder.print_report()

if __name__ == "__main__":
    main()

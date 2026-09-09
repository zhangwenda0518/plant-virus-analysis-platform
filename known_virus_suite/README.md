# known_virus_suite

已知病毒鉴定 → 二次过滤 → 共识序列 → 深度绘图 → 变异注释，五段整合成一个可独立调用、也可一键全跑的模块。

复制改造自 `D:\桌面\延伸基因组\MMPV-RNA\virome_analysis_pipeline`，**不修改原管线**。覆盖度/深度统计默认走原管线的确切路径（`samtools view -F 0x04` + `pandepth -a`），已与独立复现脚本逐位对拍通过。

---

## 一、快速开始

```powershell
cd D:\桌面\植物病毒分析平台\known_virus_suite

# 单样本全流程（鉴定 + 过滤 + 共识 + 绘图 + 变异）
python known_virus_suite.py all `
  --out ..\kv_run `
  --reference ..\virus-db\final.cluster.ref.fasta `
  --ref-info ..\virus-db\final.cluster.ref_info.tsv `
  --index-dir ..\virus-db\kv_index `
  -1 R1.fq.gz -2 R2.fq.gz
```

分步执行：

```powershell
python known_virus_suite.py index     --out ..\kv_run --reference ..\ref\ref.fasta
python known_virus_suite.py identify  --out ..\kv_run --reference ..\ref\ref.fasta -1 R1.fq.gz -2 R2.fq.gz
python known_virus_suite.py filter    --out ..\kv_run
python known_virus_suite.py consensus --out ..\kv_run --reference ..\ref\ref.fasta -1 R1.fq.gz -2 R2.fq.gz
python known_virus_suite.py plot      --out ..\kv_run
python known_virus_suite.py variant   --out ..\kv_run --reference ..\ref\ref.fasta
```

批量样本：

```powershell
# 样本清单（列：name<TAB>r1[<TAB>r2]）
python known_virus_suite.py all --out ..\kv_run --reference ..\ref\ref.fasta --sample-sheet samples.tsv
```

绘图段（只出深度图，可单独重跑）：

```powershell
# 默认读 <out>/filter/filtered.tsv + <out>/align/*.SiteDepth.gz
python known_virus_suite.py plot --out ..\kv_run

# 不联网下载 GenBank 注释
python known_virus_suite.py plot --out ..\kv_run --no-genes

# 指定缓存与邮箱
python known_virus_suite.py plot --out ..\kv_run --gbk-dir ..\gbk --ncbi-email me@example.com
```

变异段（bcftools mpileup+call → SnpEff 功能注释，可单独重跑）：

```powershell
# 输入是 <out>/bam/*.sorted.bam（共识段产物），注释复用绘图段的 gbk_files/
python known_virus_suite.py variant --out ..\kv_run --reference ..\ref\ref.fasta

# 调整 caller 参数
python known_virus_suite.py variant --out ..\kv_run --reference ..\ref\ref.fasta `
  --variant-qual 3.5 --min-freq 0.05
```

---

## 二、参数

### 引擎与后端

| 参数 | 默认 | 说明 |
|------|------|------|
| `--engine` | `minibwa` | `salmon` 伪比对 / `minibwa` 真比对 |
| `--coverage-tool` | `pandepth` | `pandepth` / `samtools` / `builtin`，见第三节 |
| `--index-dir` | `<out>/index` | 索引目录。已存在 `minibwa.mbw` 则直接复用，不重建 |
| `--threads` | 自动 | 样本级并行/批次大小 |
| `--align-threads` | 8 | 引擎内部线程数 |
| `--batch-size` | 1 | 每批样本数（断点续跑粒度） |

### 鉴定段阈值

| 参数 | 默认 | 说明 |
|------|------|------|
| `--min-reads` | 10 | 最少 uniq reads |
| `--min-tpm` | 1.0 | 最低 TPM |
| `--min-cov` | 10.0 | 覆盖率下限 % |
| `--min-depth` | 0.5 | 平均深度下限 |
| `--min-poisson` | 0.3 | 泊松比值下限（打假） |
| `--min-ani` | 0.0 | ANI 下限（salmon 下不可用） |
| `--min-mapq-identify` | 10 | **仅对 builtin/samtools 后端生效**，pandepth 后端为 0 |

### 输入

| 参数 | 说明 |
|------|------|
| `-1/--r1`、`-2/--r2` | 单样本 |
| `--sample-sheet` | TSV 清单：列 `name, r1, r2` |
| `--reference` | 参考 FASTA |
| `--ref-info` | 参考注释 TSV（Accession/Taxid/Species/Segment） |

### 绘图段

| 参数 | 默认 | 说明 |
|------|------|------|
| `--window` | 10 | 深度曲线滑动平均窗口 |
| `--plot-fontsize` | 9 | 统计信息框字号 |
| `-g/--add-genes` | 开 | 叠加 GenBank 基因轨道 |
| `--no-genes` | — | 关闭基因轨道（不联网） |
| `--gbk-dir` | `<out>/gbk_files` | GenBank 缓存目录 |
| `--depth-dir` | `<out>/align` | per-position 深度表目录 |
| `--summary` | `<out>/filter/filtered.tsv` | 结果表 |
| `--ncbi-email` | 环境变量 | NCBI Entrez 邮箱（建议填真实邮箱以免被限流） |
| `--ncbi-api-key` | 环境变量 | NCBI API key |

### 变异段

| 参数 | 默认 | 说明 |
|------|------|------|
| `--bam-dir` | `<out>/bam` | 输入 BAM 目录（共识段产物） |
| `--gbk-dir` | `<out>/gbk_files` | GenBank 缓存（复用绘图段的，无则下载） |
| `--variant-qual` | `3.5` | QUAL 下限（bcftools 标度，见下） |
| `--min-freq` | `0.05` | 最小等位频率（`AD[1]/DP`） |
| `--variant-ncbi-email` | 同 `--ncbi-email` | 变异段单独的 Entrez 邮箱 |
| `--variant-ncbi-api-key` | 同 `--ncbi-api-key` | 变异段单独的 API key |

默认参数集中在 `kv_config.yaml`。命令行优先级更高。

---

## 三、三个覆盖度后端

| 后端 | 实现 | MAPQ 过滤 | 相对速度 | 用途 |
|------|------|-----------|----------|------|
| **pandepth**（默认） | `samtools view -F 0x04` + `pandepth -a`（无 `-q`） | 无（0 起算） | 5.5s | **严格照原管线**，口径与论文一致 |
| samtools | 同一条 BAM + `samtools coverage` | 无 | 0.5s | 快速对拍/验证 |
| builtin | 纯 Python 解析 SAM/CIGAR | `--min-mapq-identify` | 0.3s | 外部工具缺失时的降级 |

三者已在 truth30 / truth3 / q06 上对拍，Coverage 与 MeanDepth **零分歧**（容差 `1e-9`）。

复现验证：

```powershell
python replicate_orig_coverage.py q06     # 照原管线独立复现
python verify_orig_parity.py q06          # 模块 vs 复现脚本，期望 PASS (bit-identical)
```

### 归一化口径

- pandepth 与 samtools/builtin 的 MAPQ 起点不同：**pandepth 为 0**，builtin 为 `--min-mapq-identify`
- 因此默认路径下 `--min-mapq-identify` 不生效，改它不会影响结果
- `hits` 有歧义：pandepth 的 `Cov>0` 条数（q06=129）≠ 鉴定表行数（q06=73，只保留有 EM reads 的条目）

---

## 四、依赖

### 必需

| 工具 | 路径 | 用途 |
|------|------|------|
| minibwa | `minibwa_win_build\minibwa.exe` | 真比对引擎 |
| samtools | `tools\samtools\bin\samtools.exe` | BAM 转换/排序/索引 |
| pandepth | `tools\pandepth\pandepth.exe` | 深度统计（默认后端） |
| viral_consensus | `tools\viral_consensus\viral_consensus.exe` | 共识序列构建（共识段必需） |
| bcftools | `tools\bcftools\bin\bcftools.exe` | 变异检出 caller（变异段必需，替代原管线 iVar/freebayes） |

### 可选

| 工具 | 路径 | 用途 |
|------|------|------|
| salmon | `tools\salmon2\salmon.exe` | 伪比对引擎（`--engine salmon`） |
| mafft | `tools\mafft-win\ms\bin\mafft` | 共识序列比对 |
| Java 21 | `tools\jdk21\jdk-21.0.12.1+1\bin\java.exe` | 跑 SnpEff（本机默认 java 是 1.8，不能用） |
| SnpEff 5.4c | `tools\snpeff\snpEff\snpEff.jar` | 变异功能注释（变异段可选，缺则跳过注释） |

### Python 包（绘图段）

| 包 | 用途 |
|----|------|
| matplotlib ≥3.11 | 绘图（300 dpi、`pdf.fonttype 42`） |
| pandas | 结果表/深度表读写 |
| numpy | 滑动平均与数组运算 |
| biopython | `Bio.Entrez.efetch` 下载 GenBank、`Bio.SeqIO` 解析特征 |

工具探测用 `python -c "from kv_common import ToolRegistry; t=ToolRegistry(); t.probe()"`。

### ⚠ pandepth 的 DLL 陷阱

pandepth 的 `hts-3.dll` 有传递依赖（libcrypto/libcurl 等）落在 `samtools\bin`。直接调 `pandepth.exe` 会报 `0xC0000135`（找不到 DLL）。模块已自动把 `pandepth` 与 `samtools\bin` 两个目录注入子进程 `PATH`，手工调用时需自己加。

### ⚠ viral_consensus 的 PATH 顺序硬约束

viral_consensus.exe 与 samtools 自带同名不同版本的 htslib DLL。实测：

| PATH 顺序 | exit |
|-----------|------|
| 仅默认 PATH | 3221225785 (0xC0000139) |
| **MinGW bin > viral_consensus > samtools\bin** | **正常** |
| viral_consensus > samtools\bin > MinGW bin | 3221225477 (0xC0000005) |

`build_runtime_env()` 已自动按正确顺序构造，用 `shutil.which('g++')` 探测 MinGW。

### Windows 上没有的东西

- **pysam**：无 Python 3.12 预编译 wheel，源码编译失败。模块用 `samtools index` 等价替代 `pysam.index()`，功能一致。
- **sh / awk / sed / grep / gzip / efetch**：全链不依赖任何 Unix 文本工具。原管线的 awk `@SQ` 瘦身由 `_filter_sam_header()` 用 Python 实现，`gzip -dc | rg` 深度表抽取由 `load_depth_table()` 逐行读实现，`efetch` 由 `Bio.Entrez.efetch` 替代。

---

## 五、产物结构

```
out_dir/
├── logs/                          运行日志
├── index/                         索引（--index-dir 指定时改用外部索引）
├── align/{sample}.sam             minibwa SAM
├── align/{sample}_pandepth.sorted.bam(.bai)  排序 BAM（覆盖度与变异段共用）
├── quant/{sample}/quant.sf        salmon 定量
├── stats/{sample}.refstats.tsv    覆盖度/深度
├── batches/batch_0001.json        鉴定段断点续跑
├── identify/
│   ├── all_viruses.summary.tsv        全量鉴定表
│   ├── all_viruses.best.summary.tsv   置信白名单
│   └── all_viruses.unclassified.tsv   疑似新种
├── filter/
│   ├── filtered.tsv               通过二次过滤
│   └── discarded.tsv              被剔除
├── align/
│   └── {sample}_pandepth.pandepth.SiteDepth.gz   逐位点深度（绘图输入）
├── plots/
│   └── {sample}/{sample}_{taxonomy}_{accession}_depth.{pdf,png}
├── gbk_files/                     GenBank 注释缓存（SnpEff 建库直接复用）
├── consensus/
│   ├── {sample}.consensus.fasta
│   ├── {sample}.qc.tsv
│   └── features/
├── bam/{sample}.sorted.bam        共识段与变异段输入 BAM
├── vcf/{accession}.vcf            每参考的过滤后变异
├── annotated/{accession}.ann.{vcf,tsv}   SnpEff 注释结果
├── snpeff_work/                   SnpEff 工作目录（config + data/）
└── variant_summary.json           变异段汇总（含 skipped 列表）
```

`variant_summary.json` 结构：

```json
{
  "n_success": 2,
  "n_skipped": 2,
  "results": [{"accession": "OR489165.1", "n_variants": 53, "n_coding_ann": 34, ...}],
  "skipped": [{"genome": "NC_002030", "reason": "无 CDS/gene 注释（viroid 等非编码病毒）"}]
}
```

`skipped` 里的参考属正常情形：viroid（PSTVd、CEVd 等）在 GenBank 里只有 `source` 一个 feature，
没有 CDS/gene/exon，SnpEff 的 `build -genbank` 无法建库。已实测对照（见 §八）。

`stats/{sample}.refstats.tsv` 列定义：

```
Accession  Length  Mapped_Reads  Covered_Bases  Coverage(%)  MeanDepth
```

这是两个引擎的统一交汇点。

---

## 六、实测基线（GQMIX.q06，8464 条参考）

| 指标 | 值 |
|------|-----|
| 比对耗时 | 5.1s（minibwa，8 线程） |
| 全流程耗时 | 7.9s |
| spliced mapped | 327 / 1051 |
| 鉴定表行数 | 73 |
| pandepth `Cov>0` 条数 | 129 |

关键参考的覆盖度（与原管线口径逐位一致）：

| Accession | Coverage(%) | MeanDepth |
|-----------|-------------|-----------|
| NC_002030.1 | 100.00 | 26.36 |
| NC_000885.1 | 71.67 | 3.11 |
| LC902918.1 | 13.02 | — |
| PP333054.1 | 9.30 | 0.41 |
| PP803563.1 | 5.32 | — |
| PP317290.1 | 3.97 | — |

---

## 七、文件说明

| 文件 | 作用 |
|------|------|
| `known_virus_suite.py` | 主入口，五段编排 |
| `kv_engines.py` | 引擎抽象层 + 三个覆盖度后端 |
| `kv_identify.py` | 鉴定段 |
| `kv_filter.py` | 过滤段 |
| `kv_consensus.py` | 共识段 |
| `kv_plot.py` | 绘图段（深度图 + GenBank 基因轨道） |
| `kv_variant.py` | 变异段（bcftools mpileup+call caller + SnpEff 注释） |
| `kv_common.py` | 工具探测、I/O、日志、参考加载 |
| `kv_config.yaml` | 默认参数 |
| `DESIGN.md` | 设计文档（含差异清单与错判更正） |
| `replicate_orig_coverage.py` | 照原管线独立复现覆盖度 |
| `verify_orig_parity.py` | 模块 vs 复现脚本对拍 |
| `cmp_coverage_backends.py` | 三后端横向对比 |
| `diag_exact.py` | 分歧位点单变量扫描 |
| `make_truth_reads.py` | 生成模拟真值 reads |
| `check_consensus.py` | 共识序列正确性校验 |

---

## 八、已知限制

- **ANI/Pi** 在 salmon 引擎下不可用（伪比对无位点）
- **未移植**：PVGA 延伸、gmcloser/abyss-sealer gap filling、环化检测、完整 ANI 实现。这些属于组装后处理，与鉴定+过滤+共识的核心目标无关
- **pandepth 启动开销** 约 5s（与数据量无关，空 BAM 也 58s 的旧测法是错的，见 DESIGN.md 第九节），样本量很大时可临时切 `--coverage-tool samtools`
- **viroid 无法做编码区注释**：viroid（PSTVd、CEVd 等）的 GenBank 记录只有 `source` 一个 feature，CDS/gene/exon 均为 0，SnpEff 建库会报 `FATAL ERROR: Most Exons do not have sequences!`。这是数据特性而非代码缺陷。实测对照：

  | 参考 | 描述 | CDS | gene | exon | SnpEff 建库 |
  |------|------|-----|------|------|------------|
  | NC_002030.1 | Potato spindle tuber viroid | 0 | 0 | 0 | ✗ |
  | NC_000885.1 | Tomato chlorotic dwarf viroid | 0 | 0 | 0 | ✗ |
  | OR489165.1 | Cytorhabdovirus sp. 'lycii' | 6 | 6 | 0 | ✓ |
  | LC902918.1 | Potexvirus pepini | 5 | 5 | 0 | ✓ |

  处理方式：`has_coding_features()` 先读 FEATURES 段预检，无编码 feature 则直接抛自解释错误，日志记 warning（不记 error），并写入 `variant_summary.json` 的 `skipped` 列表。页面汇总区以橙色徒标「跳过注释 N」展示。viroid 的变异位点仍然会被检出并写入 `vcf/`，只是不做基因功能注释。
- **bcftools 插件路径**：`ToolRegistry` 会为子进程设 `BCFTOOLS_PLUGINS=<tools>\bcftools\libexec\bcftools`，手工调用 `bcftools +counts` 等插件时需自己设。

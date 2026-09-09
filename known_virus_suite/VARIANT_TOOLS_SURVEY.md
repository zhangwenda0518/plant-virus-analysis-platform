# 变异分析工具 Windows 可行性调研（实测版）

> 调研时间：2026-09-08
> 调研方式：全部本机实测，不接受文档推断

## 一、结论速览

| 工具 | Windows 可行性 | 本机状态 | 结论 |
|---|---|---|---|
| **SnpEff 5.4c** | ✅ 已跑通 | 已部署，端到端验证通过 | **采用** |
| iVar 1.4.4 | ❌ 不可得 | 无包、无源码构建链 | 放弃 |
| SNPGenie | ⚠️ 需 Perl | 本机无 perl | 暂缓 |

---

## 二、SnpEff（已落地）

### 2.1 部署事实

| 项 | 值 |
|---|---|
| 版本 | SnpEff 5.4c (build 2026-02-23) |
| 位置 | `tools/snpeff/snpEff/snpEff.jar` (27.3 MB) |
| Java | 平台自带 Temurin **JDK 21.0.12.1 LTS** (`tools/jdk21/jdk-21.0.12.1+1/`) |
| 下载源 | `https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip` (63.6 MB) |
| JDK 源 | `https://api.adoptium.net/v3/binary/latest/21/ga/windows/x64/jdk/hotspot/normal/eclipse` (195.6 MB) |
| 磁盘占用 | jdk21 523.5 MB + snpeff 229.4 MB = 753 MB |

### 2.2 为什么必须自带 JDK 21

本机系统 `java` 是 **1.8.0_341（class 52）**，SnpEff 5.4c 编译为 **class 65（Java 21）**。
实测报错：

```
java.lang.UnsupportedClassVersionError:
org/snpeff/SnpEff has been compiled by a more recent version of the Java Runtime
(class file version 65.0), this version of the Java Runtime only recognizes
class file versions up to 52.0
```

SnpEff 官方文档也已明确「requires Java 21 or later」。旧版 4.3t 只存在于 SourceForge，
实测直链返回 HTML 重定向页（0.1 MB）、`downloads.sourceforge.net` 返回 403，
GitHub releases 资产为空数组，取用不可靠。**自带 JDK 21 是唯一稳妥路线。**

### 2.3 端到端实测结果

测试对象：`OR489165.1`（Cytorhabdovirus sp. 'lycii'，枸杞 crinkle 病毒，14812 bp，6 个 CDS：N/P/P4/M/G/L）

```
GenBank (Bio.Entrez.efetch)
   ↓ SnpEff build -genbank
snpEffectPredictor.bin + sequence.bin
   ↓ poscounts.tsv → VCF
SnpEff ann → ANN 字段 → 结构化 TSV
```

建库输出：
```
Protein check: OR489165  OK: 6  Not found: 0  Errors: 0  Error percentage: 0.0%
```

注释结果（4 条编码区变异）：

| POS | 变异 | 基因 | 效果 | 影响 | HGVS.c | HGVS.p |
|---|---|---|---|---|---|---|
| 438 | A>G | N | start_retained_variant | LOW | c.3A>G | p.Met1? |
| 8000 | T>C | L | missense_variant | MODERATE | c.439T>C | p.Ser147Pro |
| 8002 | T>A | L | synonymous_variant | LOW | c.441T>A | p.Ser147Ser |
| 8005 | T>G | L | synonymous_variant | LOW | c.444T>G | p.Leu148Leu |

单库建库耗时 **0.27 秒**。

### 2.4 三个必须记住的坑

1. **中文路径 + Java 编码**
   不加 `-Dfile.encoding=UTF-8` 时路径显示为 `D:\����\ֲ�ﲡ������ƽ̨`，
   文件读取随即失败。所有调用必须带该参数。

2. **config 里 `data.dir` 必须用相对路径**
   `data.dir = D:\...\data` 会被 SnpEff 当成转义序列吃掉，报：
   ```
   Cannot read file '...\snpeff_work/D:����ֲ�ﲡ...data/OR489165/genes.gbk'
   ```
   正确写法 `data.dir = ./data`，且**工作目录必须切到 config 所在目录**执行。

3. **染色体名必须与库内一致（保留版本号）**
   SnpEff 库内染色体名是 `OR489165.1`，VCF 写 `OR489165` 会报
   `ERROR_CHROMOSOME_NOT_FOUND`，ANN 字段降级成空白注释。
   平台端必须做自动对齐（从 GenBank `VERSION` 行取 ID）。

### 2.5 附带收获：SnpEff 会校验参考碱基

人工构造的 VCF 中 REF 与库内序列不符时，输出
`WARNING_REF_DOES_NOT_MATCH_GENOME`。这说明它不是单纯按坐标套注释，
而是真的做了参考序列一致性检查，对结果可信度是正面信号。

---

## 三、iVar（不可行，已定性）

| 项 | 实测结果 |
|---|---|
| bioconda 平台 | 仅 `linux-64` / `linux-aarch64` / `osx-64` / `osx-arm64`（61 个文件），**无 win-64** |
| 官方二进制 | 不提供 Windows 版 |
| 源码构建 | 需 autotools（autoconf/automake/libtool）+ HTSlib + GCC≥5.0 |

本机实测：gcc 16.2.0 ✅、samtools 1.24 + htslib 1.24 ✅、make ✅，
但 **autoconf / automake / libtool / perl / cpan 全部缺失**，且无 bash / Git Bash / WSL。
按官方构建路径堵死。

**功能重叠判断**：平台现有 `viral_consensus/variants.tsv` 已产出
`contig/pos/ref_base/cons_base/depth/A/C/G/T/N/gap/major_base/major_freq/minor_base/minor_freq/type`，
iVar 的主要能力（iSNV 检出）已覆盖，边际价值低。

---

## 四、SNPGenie（需 Perl）

官方描述：「command-line application written in Perl, with no additional dependencies
beyond the standard Perl package」。

- **输入**：单序列 FASTA + GTF 注释 + VCF/CLC/Geneious SNP 报告
- **输出**：πN/πS、dN/dS、gene diversity
- **三种分析**：within-pool (`snpgenie.pl`)、within-group、between-group
- **VCF 必须指定 `--vcfformat=4`**

本机 perl 缺失。若要用，需额外部署 Perl（Strawberry Perl / MSYS2）。
与 SnpEff 的分工很清楚：SnpEff 做功能注释，SNPGenie 做群体遗传学指标（πN/πS）。

---

## 五、平台集成建议

### 5.1 数据链路

```
样本 poscounts.tsv (viral_consensus 产物，已有)
   ↓ poscounts_to_vcf()          ← kv_variant.py 已实现并验证
VCF (CHROM 自动对齐库内命名)
   ↓ SnpEff build -genbank       ← GenBank 复用绘图段缓存
SnpEff 库
   ↓ SnpEff ann
ANN 注释 VCF
   ↓ parse_ann()
结构化变异注释 TSV (基因/效果/影响/HGVS.c/HGVS.p)
```

### 5.2 已验证的数字

- 平台自产 poscounts → VCF：NODE_1 得 **43 条变异**，
  与 `summary.json` 的 `n_isnv: 43` **完全一致**
- 单库建库 0.27 秒，8464 条参考库理论约 38 分钟
- 建库必须有 GenBank 注释，`virus-db/final.cluster.ref.fasta` 只有序列，
  注释需走 NCBI（复用 `kv_plot.fetch_gb`）

### 5.3 待定问题

- 批量建库策略：全量 8464 条 vs 按样本检出病毒按需建库
  （倾向按需，避免 8464 次 NCBI 请求）

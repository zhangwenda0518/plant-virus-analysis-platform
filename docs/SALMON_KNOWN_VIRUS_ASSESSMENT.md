# 已知病毒快速鉴定：salmon 接入评估

评估对象：`D:\桌面\延伸基因组\MMPV-RNA\virome_analysis_pipeline`（10 阶段已知病毒管线）
评估时间：2026-09-08
本机 salmon：`D:\桌面\植物病毒分析平台\tools\salmon2\salmon.exe`（v2.7.0，Rust 版，本机编译）

---

## 一、管线现状：salmon 已经是第一公民

这条管线在设计上**已经支持 salmon**，且是 Stage 1 的默认引擎（`auto_known_virus.py --tool` 默认 `salmon`）。

Stage 1（`batch_virus_depth.py`，820 行）的双引擎架构：

| 维度 | 伪比对（salmon / kallisto） | 传统比对（bowtie2 / bwa / minimap2 / strobealign / hisat2） |
|------|------------------------------|-------------------------------------------------------------|
| 索引 | `salmon index -k 31` | `bowtie2-build` / `bwa index` 等 |
| 定量 | `salmon quant -l A -1/-2 ...` | 比对 → `samtools sort` → `pandepth` |
| 深度来源 | `quant.sf` 的 `NumReads` + `--writeMappings` 管道出的 BAM | BAM |
| ANI/Pi | **不可用**（伪比对不产生真实比对） | 从 BAM 逐 read 算 `(aln_len - NM)/aln_len` |

`is_pseudo = tool in ['kallisto','salmon']` 这个开关决定后续走哪条分支。

**salmon 的定位是"快速筛查"，传统比对是"精细验证"。** 两者互补，不是替代关系。

---

## 二、Stage 1 的核心算法（值得保留的部分）

### 1. Poisson 打假模型

这是整条管线最有价值的创新，直接决定假阳性率：

```
λ = (Rep_Reads × Avg_Read_Len) / Rep_Length
Predicted_Support = 1 - exp(-λ)
Poisson_Ratio = Observed_Coverage / Predicted_Support
```

判读口径（代码注释原文）：

- `≈ 1.0` reads 随机分布，正常
- `< 1.0` reads 局部堆叠，可能假阳性（来自重复区域或非特异性比对）
- `> 1.0` 比随机更均匀，高置信真阳性

注意 `Rep_Reads` 用的是**参考序列上的 reads**，不是 EM 分配后的 `Asm_EM_Reads`。这个口径让 Poisson 模型在伪比对模式下依然成立。

### 2. 双轨过滤

```
filter_base = Sample_Total_Mapped > 0
            & Uniq_Reads >= min_uniq_reads
            & Asm_TPM   >= min_tpm

Track A（全基因组）
  Rep_Coverage(%) >= coverage
  & Poisson_Ratio >= ratio        # DNA 用 max(ratio, 0.5)，更严
  & Rep_MeanDepth >= meandepth

Track B（基因区，仅 RNA 病毒启用）
  gene_total_cov >= 80.0
  & gene_avr_cov >= 5.0

final = filter_base & (Track A | Track B)
```

DNA 比 RNA 的 Poisson 阈值更严，理由是 DNA 病毒突变率低、信号更可靠。这个区分有道理。

### 3. 默认阈值

| 参数 | 默认 | 含义 |
|------|------|------|
| `--coverage` | 10.0 | 代表序列全长覆盖度下限 % |
| `--ratio` | 0.3 | Poisson Ratio 下限 |
| `--meandepth` | 0.5 | 代表序列平均深度 |
| `--min_tpm` | 1.0 | 过滤痕量污染 |
| `--min_uniq_reads` | 10 | 过滤 1-2 条 reads 的假阳性 |
| `--sp_thresh` | 95.0 | 物种 ANI 识别阈值 % |
| `--min_gene_total_cov` | 80.0 | B 轨转录区总覆盖 |
| `--min_gene_avr_cov` | 5.0 | B 轨转录区平均覆盖 |

### 4. 物种归属

```
sp_candidates = ['Species_NCBI', 'Species_ICTV', 'Base_Parsed_Species']
```

按此顺序 coalesce。伪比对模式 `Avg_Read_ANI` 全为 None，所以 `df_confirmed` 的条件 `(ANI is null)` 让所有记录都进 confirmed 分支，实际等于**没有 ANI 校验**。

---

## 三、Stage 2 过滤（`utils/filter_summary.py`，236 行）

Stage 1 内部已经过滤过一轮，Stage 2 是**独立可调的二次过滤**，通过 `--filter` 启用。

列名映射（注意大小写）：

```python
KNOWN_COLS = {
    "coverage": "Rep_Coverage(%)",
    "depth":    "Rep_MeanDepth",
    "reads":    "Asm_EM_Reads",
    "tpm":      "Asm_TPM",
    "ani":      "Avg_Read_ANI",
    "poisson":  "Poisson_Ratio",
}
```

CLI：`--min_cov / --min_depth / --min_reads / --min_tpm / --min_ani / --min_poisson`，全部默认 0.0（不启用）。

列不存在时打 warning 跳过，不报错。这个行为对伪比对模式友好（ANI 列存在但为空）。

---

## 四、共识序列链路

### Stage 5 `consensus_extract.py`（231 行）

对 `virus-full.py` 产出的共识序列做 QC：

1. `calculate_assembly_stats`：contig 数、总长、N50、GC%、N 比例
2. `check_protein_integrity`：翻译后查内部终止子（尾部 `*` 不计）
3. `fill_n`（可选）：逐位用 GenBank 参考填补 N，同时记录填补数
4. `download_genbank`：`Entrez.efetch` 联网取参考，跳过 `[SCE]RR` accession
5. 输出 TSV 汇总

### Stage 4 `virus-full.py`（1120 行）12 步引擎

关键链路（Step 9 是真共识来源）：

```
Step 1&2  de novo 组装（SPAdes 等）
Step 3    Shiver 净化
Step 4    Ref Merge
Step 5    PVGA 延伸 → 可能触发"快车道"跳过 6-8
Step 6    Pre-Fusion
Step 7    rmDup 去冗余
Step 8    Divine Fusion
Step 9    迭代抛光 ← 共识序列真正生成的地方
Step 10   Gap Filling（gmcloser + abyss-sealer）
Step 11   环化检测与 Trim
Step 12   演化图谱
```

Step 9 的命令：

```bash
minimap2 -t N -a -x PRESET ref.fasta reads.fq | \
  tee >(viral_consensus -i - -r ref -o raw_cons \
        --min_qual Q --min_depth D -op position.tsv) | \
  samtools view -b -@ N > mapped.bam
```

迭代逻辑：每轮用 mafft 把共识与参考对齐，**共识的非 N 碱基保留、N 用参考回填**，作为下一轮参考。这个设计把"测序偏差纠正"和"参考回填"分开，避免单轮引入参考偏倚。

---

## 五、接入本机的硬障碍

| 依赖 | 用途 | 本机状态 |
|------|------|----------|
| `salmon` | 定量 | ✅ 已装（v2.7.0 本机编译） |
| `samtools` | BAM 排序/索引/统计 | ❌ 无 |
| `pandepth` | 深度/覆盖度统计 | ❌ 无 |
| `minimap2` | Step 9 比对 | ✅ `tools\minimap2` |
| `mafft` | Step 9 迭代对齐 | ✅ `tools\mafft-win` |
| `viral_consensus` | Step 9 共识生成 | ❌ 无（Linux 二进制） |
| `gmcloser` / `abyss-sealer` | Step 10 补洞 | ❌ 无 |
| `bowtie2` / `bwa` / `kallisto` | 传统比对/伪比对备选 | ❌ 无 |
| 病毒参考 FASTA | salmon 建索引输入 | ❌ 本机无（只有 kraken2 库和 acc2taxid 映射） |

**命令行本身还是 Linux 专用**：

```python
# L379
full_cmd = f"/usr/bin/time -v salmon quant ... --writeMappings /dev/stdout | {samtools_pipe}"
subprocess.run(full_cmd, shell=True, executable='/bin/bash', check=True, stderr=flog)
```

`/usr/bin/time -v`、`/dev/stdout`、`executable='/bin/bash'`、`pandepth`、`tee >(...)` 进程替换，这些都跑不了 Windows。

---

## 六、接入方案（三档，按投入排序）

### 方案 A：独立快速筛查脚本（推荐先做）

不碰原管线，在本机写一个 `salmon_screen.py`，复用原管线的**判定口径**但换掉执行层：

```
输入：clean reads + 病毒参考 FASTA + ref_info.tsv
  ↓
salmon index -k 31
  ↓
salmon quant -l A -1 R1 -2 R2 --writeBam mappings.bam -o quant/
  ↓
从 quant.sf 读 Name/NumReads/TPM
从 --writeBam 的 BAM 算覆盖率与深度（替代 pandepth）
  ↓
按原口径算 Poisson_Ratio、双轨过滤
  ↓
输出 summary.tsv（列名与原管线一致，下游可直接读）
```

关键收益：**salmon 2.7.0 原生支持 `--writeBam`**，比原管线用 `--writeMappings /dev/stdout | samtools` 管道干净得多，Windows 上不需要 `tee` 和进程替换。

需要新写的只有：BAM 覆盖率/深度统计（Python + pysam，或借 samtools 的 `depth`/`coverage`）。本机无 samtools，但 pysam 自带 htslib，可以纯 Python 完成。

前置缺口：**要先把病毒参考 FASTA 准备好**。本机 `databases\virus\ref` 只有 kraken2 库，没有序列。

### 方案 B：给原管线打 Windows 兼容补丁

改 `batch_virus_depth.py`，加平台分支：

```python
IS_WIN = os.name == 'nt'
TIME_PREFIX = '' if IS_WIN else '/usr/bin/time -v '
EXECUTABLE = None if IS_WIN else '/bin/bash'
```

并把 `--writeMappings /dev/stdout | samtools` 换成 `--writeBam out.bam`，把 `pandepth` 换成自研统计。

改动面大（Pipeline 类 + 两个 process 方法 + 索引分支），且原管线在 246 服务器上跑得好，为 Windows 改它有维护成本。除非确定要在本机跑全链，否则不建议。

### 方案 C：只在本机做筛查，重活仍交服务器

本机用 salmon 快速扫一遍出候选，把候选交给 246 服务器跑完整 10 阶段。这样各取所长：

- 本机 salmons 伪比对，速度是传统比对的 10-100 倍，适合大批量样本初筛
- 服务器有完整工具链（samtools/pandepth/viral_consensus/gmcloser），跑精细阶段

**这条最符合你现有的工作流**：本机开发验证 + 服务器跑量。

---

## 七、版本兼容提醒

salmon 2.7.0 与 1.x 的差异（已实测）：

| 项 | 状态 |
|----|------|
| `salmon index` / `salmon quant` 工作流 | 不变 |
| `quant.sf` 格式（Name/Length/EffectiveLength/TPM/NumReads） | 不变 |
| `--writeMappings` | 仍支持 |
| `--writeBam` | 新增，可直接出 BAM |
| 索引格式 | **变了，旧索引必须重建** |
| `-l A` 自动检测 | 不变 |
| `lib_format_counts.json` | 不变 |

上游有 `MIGRATION.md` 列完整不兼容项。原管线用的参数（`index -k 31`、`quant -l A -1/-2`、`--writeMappings`）在 2.7.0 全部可用，**这部分不用改**。

---

## 八、结论

1. 原管线对 salmon 的支持是设计好的，不是硬塞。判定逻辑（Poisson 打假 + 双轨过滤）与引擎无关，可以整体复用。
2. 真正的障碍在**执行层**（Linux 专用命令、缺失的 Linux 工具），不在算法层。
3. 伪比对模式拿不到 ANI，所以 salmon 只能做筛查。物种级确认仍需传统比对，或者接受 `Sp_Confirmed` 条件失效的现状。
4. 建议走**方案 A 或 C**：本机轻量筛查 + 服务器精细分析。直接给原管线打 Windows 补丁性价比不高。
5. 前置动作：准备病毒参考 FASTA（可从 NCBI RefSeq viral 或你服务器上的参考集拷过来）。

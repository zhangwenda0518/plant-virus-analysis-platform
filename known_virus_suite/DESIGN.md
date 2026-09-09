# known_virus_suite 设计文档

**目标**：把「已知病毒鉴定 → 二次过滤 → 共识序列 → 深度绘图」四段整合成一个模块，引擎层支持 salmon 与 minibwa。

**来源**：复制改造自 `D:\桌面\延伸基因组\MMPV-RNA\virome_analysis_pipeline`，不修改原管线。

---

## 一、原管线三段的实际数据流

| 段 | 脚本 | 输入 | 核心动作 | 输出 |
|----|------|------|----------|------|
| 鉴定 | `batch_virus_depth.py` | reads + 参考 FASTA + ref_info | 比对/定量 → 深度统计 → ANI/Pi → 双轨过滤 → 分类 | `summary/all_viruses.best.summary.tsv`<br>`summary/all_viruses.unclassified.tsv`<br>`summary/all_viruses.summary.tsv` |
| 过滤 | `utils/filter_summary.py` | 上一步的 TSV | 纯表格阈值 + 节段完整性规则 | `filtered.tsv` / `discarded.tsv` |
| 共识 | `batch_virus_variants.py` `worker_consensus` | BAM + 参考 | samtools 取 SAM → Python 瘦身 @SQ → samtools view -b/sort → **viral_consensus** | `*.consensus.fasta`<br>`consensus_qc.tsv` |

原管线里「过滤」和「共识」是分离的独立脚本，靠人工串接。

---

## 二、整合后的模块结构

```
known_virus_suite/
├── DESIGN.md                本文件
├── known_virus_suite.py     主入口（五段编排）
├── kv_engines.py            引擎抽象层：salmon / minibwa
├── kv_identify.py           鉴定段
├── kv_filter.py             过滤段
├── kv_consensus.py          共识段
├── kv_plot.py               绘图段（深度图 + GenBank 基因轨道）
├── kv_variant.py            变异段（bcftools caller + SnpEff 注释）
├── kv_common.py             公共工具（I/O、日志、参考加载）
├── kv_config.yaml           默认参数
└── README.md                使用说明
```

模块化的意义：五段可单独调用，也可一键全跑；引擎可替换。

数据流：

```
reads ──> identify ──> filter ──> consensus ──> plot
            │            │           │           ↑
            └────────────┴───────────┘           │
              stats/ + align/*.SiteDepth.gz ─────┘
                                       │
                                       └──> variant
                                            bam/ + gbk_files/
```

---

## 三、引擎抽象层契约（kv_engines.py）

所有引擎实现同一接口，上层不关心底层是伪比对还是真比对。

```python
class VirusEngine:
    name: str

    def build_index(self, ref_fasta: Path, index_dir: Path, threads: int) -> Path:
        """建索引，返回索引路径"""

    def align(self, index_path: Path, sample: dict, out_dir: Path, threads: int) -> AlignResult:
        """
        执行比对，返回：
          - sam_path   : SAM 文件路径（可为 None）
          - quant_path : 定量文件路径（可为 None）
          - reads      : {accession: 定量 reads 数}
          - mapped     : 实际比对上的记录数 / 片段数
          - total      : 输入总片段数
          - has_bam    : 是否产出可用比对记录（决定共识段能否跑）
        """

    @property
    def supports_consensus(self) -> bool:
        """是否产出真实比对位点（salmon 伪比对为 False）"""
```

### 两个引擎的实现差异

| 维度 | salmon | minibwa |
|------|--------|---------|
| 比对方式 | 伪比对（k-mer 精确匹配） | 真比对（BWA 风格，允许错配） |
| 索引 | `salmon index -k 31`，738MB | `minibwa index`，70MB |
| 定量输出 | `quant.sf`（NumReads / TPM） | 无，需从 SAM 统计 |
| 比对位点 | **无** | SAM 记录（POS + CIGAR + MD） |
| 共识支持 | **不支持**（supports_consensus=False） | 支持 |
| 远缘检测 | 弱（k-mer 必须精确） | 强（容忍错配） |
| 深度统计 | 依赖外部 pandepth | 默认 pandepth（严格照原管线），可切 samtools / 内置 |

### 关键决策：salmon 的共识兜底

salmon 不产出比对位点，而共识段必须从 BAM 出共识。**采用 fail-fast**：

- 共识段强制要求 minibwa（或其它真比对引擎）
- `tools.require('samtools', 'viral_consensus')` 在入口处检查，缺件直接报错
- 缺工具不降级、不自研代替

---

## 四、中间产物契约

```
out_dir/
├── logs/                    运行日志
├── index/                   索引文件
├── align/                   比对原始产物
│   └── {sample}.sam         minibwa 的 SAM（salmon 无）
├── quant/                   salmon 定量
│   └── {sample}/quant.sf
├── stats/                   每样本每参考的覆盖度/深度
│   └── {sample}.refstats.tsv
├── identify/
│   ├── all_viruses.summary.tsv        全量鉴定结果
│   ├── all_viruses.best.summary.tsv   置信白名单
│   └── all_viruses.unclassified.tsv   疑似新种
├── filter/
│   ├── filtered.tsv         通过二次过滤
│   └── discarded.tsv        被过滤
└── consensus/
    ├── {sample}.consensus.fasta       共识序列
    ├── {sample}.qc.tsv                逐特征 QC
    └── features/                      蛋白/核酸 FASTA
```

### stats/{sample}.refstats.tsv 列定义

```
Accession  Length  Mapped_Reads  Covered_Bases  Coverage(%)  MeanDepth
```

这个表是两个引擎的统一交汇点：salmon 用 quant.sf 的 NumReads 填 `Mapped_Reads`，minibwa 解析 SAM 得到全部字段。

---

## 五、过滤段逻辑（移植自 filter_summary.py）

### 阈值参数

| 参数 | 默认 | 作用 |
|------|------|------|
| `--min_cov` | 10.0 | 覆盖率下限 % |
| `--min_depth` | 0.5 | 平均深度下限 |
| `--min_reads` | 10 | 最少比对 reads |
| `--min_tpm` | 1.0 | 最低 TPM |
| `--min_ani` | 0.0 | 最低 ANI（salmon 下不可用） |
| `--min_poisson` | 0.3 | 最低泊松比值 |

### 双轨放行

- **A 轨**：`Coverage >= min_cov` 且 `Poisson_Ratio >= ratio` 且 `MeanDepth >= min_depth`
- **B 轨**（RNA + 提供 genes_cov 时）：`gene_total_cov >= 80` 且 `gene_avr_cov >= 5`
- 基础门槛：`Sample_Total_Mapped > 0` 且 `Uniq_Reads >= min_reads` 且 `TPM >= min_tpm`

### 节段病毒规则

某物种在参考库中若有多段，则要求全部段都检出，否则整组移到 discarded。这条规则来自 `filter_summary.py` 的 `--ref_info` 分支，保留。

---

## 六、共识段逻辑

### 主路径：完全照抄 batch_virus_variants.py

共识段严格复刻 `virome_analysis_pipeline/batch_virus_variants.py`，
唯一差异是比对引擎由 bowtie2 换成 minibwa。

原管线共识管道（`worker_consensus`, line 395）：

```bash
samtools view -h <BAM> <virus> \
  | awk -v v='<virus>' '/^@SQ/ && $2 != "SN:"v {next} {print}' \
  | samtools view -b \
  | samtools sort -o fixed.bam
samtools index fixed.bam
viral_consensus -i fixed.bam -r ref.fasta -o consensus.fasta \
                -q 20 -d 5 -f 0.5 -a N
```

移植时只做两处等价替换：

| 原管线 | 本模块 | 说明 |
|--------|--------|------|
| `bowtie2 --local` | `minibwa map` | 大王指定 |
| `awk` 过滤 `@SQ` | `_filter_sam_header()` Python 实现 | 大王指定：不用 sh/awk/sed/grep |

`_filter_sam_header` 与原 awk 语义逐行等价：@SQ 行第二列不等于 `SN:<virus>` 则丢弃，
其余行原样保留。多参考 BAM 必须瘦身成单参考，否则 viral_consensus 报
“has N references, but it should have exactly 1” 并 exit(1)。

### 管道全部走 Python 桥接

不用 `shell=True`，不依赖 sh.exe/awk.exe/sed.exe/grep.exe。
各段用 `subprocess.Popen` + stdin/stdout 对接：

```
samtools view -h  ──stdout──>  _filter_sam_header()  ──stdin──>  samtools view -b
minibwa map -o -  ──stdout──>                       ──stdin──>  samtools sort
```

### 运行时硬约束：PATH 顺序

viral_consensus.exe 与 samtools 自带同名不同版本的 htslib DLL。实测三组 PATH：

| PATH 顺序 | exit |
|-----------|------|
| 仅默认 PATH | 3221225785 (0xC0000139，DLL 入口点找不到) |
| **MinGW > viral_consensus > samtools** | **正常** |
| viral_consensus > samtools > MinGW | 3221225477 (0xC0000005) |

`build_runtime_env()` 自动构造：MinGW bin（`shutil.which('g++')` 探测）插到最前，
再是 viral_consensus 目录，再是 samtools。

### 动态深度

照抄原管线 line 1171：`--vc-depth 0` 时按 MeanDepth 动态调 `-d`。

### 蛋白完整性 QC

对齐 `consensus_extract.py`：用 GenBank 特征提取 CDS，翻译后检查内部终止符与 N 比例。

```python
# 判定
is_intact = (internal_stop_count == 0)
pass_qc   = is_intact and (n_ratio < 5.0)
```

---

## 六点五、绘图段逻辑（移植自 batch_plot_virus_depth.py --mode depth）

### 输入契约

| 输入 | 来源 | 格式 |
|------|------|------|
| 结果表 | `<out>/filter/filtered.tsv`（可 `--summary` 覆盖） | TSV，需 `Sample` + `Accession`/`Virus` 列 |
| 深度表 | `<out>/align/<sample>*.SiteDepth.gz`（可 `--depth-dir` 覆盖） | `chr	position	depth` 三列无表头 |
| GenBank | NCBI（可 `--gbk-dir` 指定缓存） | `.gb` 文本 |

深度表由鉴定段的 pandepth 后端产出（`-a` 模式会同时出 `*.chr.stat.gz` 汇总和 `*.SiteDepth.gz` 逐位点表）。

### 输出

`<out>/plots/<sample>/<sample>_<taxonomy>_<accession>_depth.{pdf,png}`，300 dpi、`pdf.fonttype 42`、白底。
未注释参考的 `taxonomy` 会等于 accession，此时文件名去重只用一次。

### 绘制内容

- 主图：`smooth(y, window)` 滑动平均深度曲线 + `fill_between` 填充，红色虚线标 MeanDepth
- 左上角统计框：Species / Coverage / MeanDepth / Reads / TPM（monospace）
- 基因轨道（`add_genes: true` 时）：CDS/gene 坐标画成带链方向箭头的多边形，9 色循环

### 移植差异

| 项 | 原管线 | 本模块 | 原因 |
|----|--------|--------|------|
| 深度表抽取 | `gzip -dc x.gz \| rg -a -F -f patterns` | Python 逐行读 + 集合过滤 | 大王指定：不用 sh/gzip/rg |
| GenBank 下载 | `efetch -db nuccore -id X -format gb` | `Bio.Entrez.efetch` | 本机无 edirect；同一 NCBI 接口 |
| 并发 | `ProcessPoolExecutor` | 串行（Windows 多进程 spawn 开销大于收益） | 单样本位点量级不大 |

绘图代码本身（`smooth`/`plot_one_virus`/基因多边形）逐字照抄。

### 无注释参考的行为

库内大量参考（如 NC_000885.1 这类 360 bp 片段记录）GenBank 里本无 CDS/gene，`parse_gb_features` 返回空列表，该图只出深度曲线，不画轨道。这是预期行为，不是失败。

---

## 六点六、变异段逻辑（kv_variant.py）

### 背景：为什么换掉原管线的 caller

原管线第五段用 iVar / freebayes / lofreq。三者在 Windows 上均不可得：iVar 无 win-64 包，
freebayes 与 lofreq 只能在 Linux。已安装并实测 **bcftools 1.24** 作为替代。
另有一条 poscounts（伪计数求变异）路线，已整体作废（见 `POSCOUNTS_REMOVAL_PLAN.md`）。

### 数据流

```
<out>/bam/*.sorted.bam  +  <out>/gbk_files/*.gb
        ↓                        ↓
  bcftools mpileup -f ref -d 100000 -q 0 -Q 13 -B
      -a FORMAT/AD,FORMAT/DP,INFO/AD -Ob
        ↓
  bcftools call -mv  →  bcftools view -Ov
        ↓
  两步过滤（QUAL + DP + AD/DP + DP4）
        ↓
  <out>/vcf/{accession}.vcf
        ↓
  SnpEff build + ann（每参考一库）
        ↓
  <out>/annotated/{accession}.ann.{vcf,tsv}
```

### 定稿参数（实测定下）

| 项 | 定稿 | 理由 |
|----|------|------|
| `--ploidy` | **不设** | 设 `1` 会把多拷贝位点压成单倍型，双端 truth set 实测丢真阳性 |
| `-A` | **不加** | 加了会保留 ADF/ADR 但 AD 语义变形，过滤表达式难对齐 |
| `-Q` | 13 | min base quality，与平台其他段一致 |
| `-B` | 开 | 关闭 BAQ，避免小参考上过度降分 |
| `-d` | 100000 | mpileup 默认深度上限太小 |
| QUAL 下限 | 3.5 | bcftools QUAL 标度远低于 freebayes，原管线的 `>20` 不可移植 |
| min-freq | 0.05 | `AD[1]/INFO/DP` |

### 过滤表达式字段映射

原管线用 freebayes 的 `AO/SAF/SAR/AF`，bcftools 下映射为：

| 原字段 | bcftools 字段 |
|--------|--------------|
| `AO` | `AD[1]` |
| `SAF` | `DP4[0]` |
| `SAR` | `DP4[1]` |
| `AF` | `AD[1]/INFO/DP` |

注意 `-a` **不能写 `INFO/DP4`**（非法），但 `call` 之后 `DP4` 会自然出现。

### SnpEff 的三个坑

1. 必须 `-Dfile.encoding=UTF-8`，否则 Windows 下中文路径与注释文本乱码
2. config 的 `data.dir` 用 `./data` 相对路径，工作目录必须切到 config 所在目录
   （绝对路径里的反斜杠会被当成转义序列吃掉）
3. 库内染色体名保留版本号（`OR489165.1`），VCF 的 CHROM 必须与之匹配，做了自动对齐

### 无编码基因参考的处理（viroid）

viroid（PSTVd、CEVd 等）GenBank 记录里只有 `source` 一个 feature，CDS/gene/exon 均为 0。
SnpEff 的 `build -genbank` 会报 `FATAL ERROR: Most Exons do not have sequences!`。
这是数据特性，不是代码缺陷。

处理方式：`has_coding_features()` 读 FEATURES 段做预检，无编码 feature 直接抛自解释错误；
调用方捕后记 **warning**（不是 error）并写入 `variant_summary.json` 的 `skipped`。
viroid 的变异位点仍会被检出并写入 `vcf/`，只是不做基因功能注释。

实测对照（本平台库）：

| 参考 | 描述 | CDS | gene | SnpEff |
|------|------|-----|------|--------|
| NC_002030.1 | Potato spindle tuber viroid | 0 | 0 | 跳过 |
| NC_000885.1 | Tomato chlorotic dwarf viroid | 0 | 0 | 跳过 |
| OR489165.1 | Cytorhabdovirus sp. 'lycii' | 6 | 6 | OK |
| LC902918.1 | Potexvirus pepini | 5 | 5 | OK |

### 性能修正：鉴定段逐参考计数

`kv_engines.count_mapped_by_ref` 原实现对 8464 条参考逐条启动 `samtools view -c`，
8364 次进程启动开销约分钟级。改为单次 `samtools idxstats` 取第 3 列，毫秒返回。
口径已对拍 12 条参考（8 非零 + 4 零值）完全一致。

### 索引复用（--index-dir）

`_index_dir(args, out_dir)` 统一取 `--index-dir`，缺省 `<out>/index`。
平台病毒库的 minibwa 索引预建在 `virus-db/kv_index/`（`minibwa.mbw` 61.9 MB + `.l2b` 7.9 MB），
仅在 engine=minibwa 且参考为平台默认时传入，避免每次跑样品重建索引。

---

## 七、与原管线的差异清单

| 项 | 原管线 | 本模块 | 原因 |
|----|--------|--------|------|
| 运行平台 | Linux（246 服务器） | Windows + Linux 双平台 | 本机是 Windows |
| 深度统计 | `pandepth -a`（无 `-q`） | **默认 pandepth，严格照原管线** | 口径一致，已逐位对拍 |
| 建索引 | `pysam.index()` | `samtools index` | Windows 无 pysam wheel；功能等价 |
| 共识器 | viral_consensus | **同一个 viral_consensus.exe**（自编译 + 补丁） | 上游 `count.cpp:192` 对 N/n 返回 -1 导致越界写 |
| 共识比对 | bowtie2 --local | minibwa map | 大王指定 |
| 管道 | bash + awk | Python 桥接 | 大王指定：不用 sh/awk/sed/grep |
| 绘图深度表 | `gzip -dc \| rg -f` | Python 逐行过滤 | 大王指定 |
| GenBank 下载 | `efetch`（edirect） | `Bio.Entrez.efetch` | 本机无 edirect |
| 共识阈值 | -q 20 -d 5 -f 0.5 -a N | 完全相同 | 照抄默认值 |
| 比对引擎 | 8 种传统 + 2 种伪比对 | salmon + minibwa（可扩展） | 本次任务范围 |
| ANI/Pi | 传统比对全支持 | minibwa 支持，salmon 置空 | 伪比对无位点 |
| 并发 | 多进程 spawn | 多进程 spawn（同） | 保持 |
| 泊松打假 | 有 | 保留 | 核心逻辑 |
| 双轨过滤 | 有 | 保留 | 核心逻辑 |

**未移植**：PVGA 延伸、gmcloser/abyss-sealer gap filling、环化检测、ANI 的完整实现（需 blast 比对全 reads）。这些属于组装后处理，与本次三段整合目标无关。

---

## 八、验证方案

用 `bench_salmon_minibwa/ref/ref.fasta`（8464 条植物病毒参考）+ GQMIX 数据端到端跑：

1. **引擎一致性**：同一份 reads，salmon 与 minibwa 分别跑鉴定段，比较检出的病毒物种集合
2. **过滤正确性**：构造边界数据（刚好卡在阈值上），验证放行/剔除符合预期
3. **共识正确性**：用模拟 reads（从参考切片段）跑共识，与原始参考比对，计算一致性
4. **端到端**：全流程跑通，产物齐全

每项都要有实测数字，不能只看「跑成功了」。

---

## 九、覆盖度后端：严格照原管线（含一次错判的更正）

### 原管线的确切做法（照搬规格）

`batch_virus_depth.py` pseudo 分支 L372-404 / traditional 分支 L493-503：

```bash
# 1. 建 BAM（只去 unmapped，不做 MAPQ 过滤）
samtools view -@ 2 -b -F 0x04 - | samtools sort -@ 2 -o sample.sorted.bam
# 2. 索引（原管线用 pysam.index，本模块用 samtools index 等价替代）
pysam.index(sample.sorted.bam)
# 3. 深度统计（无 -q，即 min MAPQ = 0）
pandepth -a -i sample.sorted.bam -o stat/sample -t <threads>
# 4. 解析 stat/sample.chr.stat.gz，取第 5 列 Cov / 第 6 列 Dep
```

要点：

- 不传 `-q`，MAPQ 0 起算，多重比对全部计入
- `-x` 默认 1796，排除 `0x4|0x100|0x200|0x800`
- `-a` 语义是「output all the site depth」，额外产逐位点文件，**不改变统计口径**

### 三条后端对拍结果（8464 条参考）

| 样本 | pandepth | samtools | 内置 | 与 pandepth 逐位一致 |
|------|----------|----------|------|----------------------|
| truth30 | 5.5s | 0.5s | 0.2s | PASS（3 条命中，0 分歧） |
| truth3 | 5.5s | 0.3s | 0.1s | PASS（2 条命中，0 分歧） |
| q06 | 5.5s | 0.5s | 0.3s | PASS（129 条命中，0 分歧） |

验证脚本：`verify_orig_parity.py` 对比模块输出与 `replicate_orig_coverage.py`（照原管线写的独立复现脚本）的 JSON 产物，容差 `1e-9`，三个样本全部 **bit-identical**。

### ❗ 一次错判的更正

我一度测得 pandepth 耗时 61.6s（q06），并据此写下「pandepth 慢 150 倍、不可用」。**这个结论是错的。**

- 错因：当时的调用路径在 `BamCoverageBackend` 里额外加了 `-q 10`，与自带封装叠加
- 实测反证：空 BAM（零记录）也要 58s，说明差异不在数据量
- 照原管线方式实测：truth30 / truth3 / q06 **均为 5.5s**

本文档保留这条记录，以标记「同一工具的不同调用方式会得出相反结论」。

### min_mapq 的作用域变化

| 后端 | MAPQ 过滤 | 说明 |
|------|-----------|------|
| pandepth（默认） | **无**（0 起算） | 严格照原管线 |
| samtools | 无 | 走同一条 `-F 0x04` BAM |
| builtin | `--min-mapq-identify`（默认 10） | 保留自研过滤，作为降级选项 |

因此 `--min-mapq-identify` 在默认路径下**不生效**，只影响 builtin 后端。

### 一个易混的陷阱：hits 的两种含义

q06 上 pandepth 报 `Cov>0` 的参考有 **129** 条，而鉴定表只有 **73** 行。两者都对：

- 129 = 有比对位点覆盖的参考数（pandepth 全量输出）
- 73 = 有 EM reads 分配进入鉴定表的条数（原管线过滤 `Reads > 0`）

对拍时不能用这两个数字互相验证。

### 内置解析器的一个真实 bug（已修）

`parse_sam_coverage` 原先把 CIGAR 的 `D`（deletion）/ `N`（intron skip）位点也累加深度，与 samtools depth 默认口径不符。

- 症状：3 条参考（PP317290.1 / PP333054.1 / PP803563.1）上有位点分歧
- 单变量扫描四个变体（去/不去 supplementary × 计/不计 deletion）：

| 变体 | PP317290.1 | PP333054.1 | PP803563.1 |
|------|-----------|-----------|-----------|
| skip-supp, 不计 D | **pos_diffs=0** | **0** | **0** |
| skip-supp, 计 D | 6 | 13 | 1 |
| keep-supp, 不计 D | 0 | 0 | 0 |
| keep-supp, 计 D | 6 | 13 | 1 |

结论：唯一真因是 `D/N` 计不计；supplementary 处理对本批数据零影响（其 MAPQ 均低于 10 被过滤）。已修为 `elif op in 'DN': rp += n`。

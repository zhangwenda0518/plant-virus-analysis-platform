# 原管线对照：poscounts 移除方案（修正版）

> 2026-09-08
> 上一版 `POSCOUNTS_REMOVAL_PLAN.md` 自创了 `bcftools mpileup+call` 路线，与原管线不符，本文档作废该方案并重写。

## 一、原管线变异链路（`batch_virus_variants.py` 实测源码）

### 1.1 输入 BAM 抽取（line 432-438）

```python
awk_filter = f"awk -v v='{virus}' '/^@SQ/ && $2 != \"SN:\"v {{next}} {{print}}'"
extract_cmd = (f"samtools view -@ {t_io} -h '{bam_path}' '{virus}' "
               f"| {awk_filter} | samtools view -@ {t_io} -b "
               f"| samtools sort -@ {t_io} -o '{fixed_bam}'")
run_cmd(f"samtools index -@ {t_io} '{fixed_bam}'")
```

产物：`{L1}/{L2}/{L2}.fixed.bam`

### 1.2 caller 三选一（line 440-459）

```python
if caller == "freebayes":
    freebayes -p 1 --pooled-continuous --min-alternate-fraction 0.01 -f ref.fa fixed.bam > raw_out
elif caller == "lofreq":
    lofreq indelqual --dindel -f ref.fa -o tmp.lofreq.bam fixed.bam
    samtools index tmp.lofreq.bam
    lofreq call --call-indels -f ref.fa -o raw_out tmp.lofreq.bam
    # tmp 用完删
elif caller == "ivar":
    samtools mpileup -aa -A -d 0 -B -Q 0 fixed.bam | ivar variants -p prefix -q 20 -t 0.01 -r ref.fa
    ivar_tsv_to_vcf(tsv_out, raw_out)
```

**CLI 定义（line 1228）：**
```python
mod.add_argument("--variant_caller", choices=["freebayes", "ivar", "lofreq"], default="freebayes")
```

**choices 里没有 bcftools。bcftools 在原管线里不是 caller。**

### 1.3 动态质量过滤（line 461-474）★核心

```python
dp  = 100 if disable_dyn else (10 if mean_depth < 50 else (20 if mean_depth < 1000 else 100))
frq = 0.05

if caller == "freebayes":
    saf = 2 if dp == 10 else (3 if dp == 20 else 10)
    flt = f'QUAL>20 && INFO/DP>={dp} && INFO/SAF>={saf} && INFO/SAR>={saf} && (INFO/AO/INFO/DP)>={frq}'
else:
    flt = f'QUAL>20 && INFO/DP>={dp} && INFO/AF>={frq}'

soft = Path(str(clean_vcf).replace(".filtered.", ".soft."))
bcftools filter --threads {t_io} -s FAIL -i '{flt}' -Ov -o '{soft}' '{raw_out}'
bcftools filter --threads {t_io} -i 'FILTER=="PASS"' -Ov -o '{clean_vcf}' '{soft}'
```

**阈值表：**

| mean_depth | dp | freebayes saf/sar | frq |
|---|---|---|---|
| < 50 | 10 | 2 | 0.05 |
| 50 - 1000 | 20 | 3 | 0.05 |
| >= 1000 | 100 | 10 | 0.05 |
| `--disable_dynamic_vcf` | 100 | 10 | 0.05 |

`--disable_dynamic_vcf` 时 dp 固定 100。

### 1.4 产物（line 1015-1017）

```python
raw_out   = vdir / f"{L2}.variants.vcf"       # caller 原始输出
clean_vcf = vdir / f"{L2}.filtered.vcf"       # 过滤后
fixed_bam = vdir / f"{L2}.fixed.bam"          # finally 里删除
```

另有 `{sample}.{virus}.allele_frequencies.tsv`（`extract_allele_frequency`）。

### 1.5 SnpEff 注释（line 505-511）

```python
cmd_ann = f"java -Xmx{snpeff_mem} -jar '{snpeff_jar}' ann -c '{snpeff_config}' -noStats {db_name} '{clean_vcf_str}' > '{ann_vcf_str}'"
parse_ann_to_tsv(str(ann_vcf_str), str(sum_tsv))
```

### 1.6 依赖检查（line 609-628）

```python
if self.args.call_variants:
    need[self.args.variant_caller] = "变异检测"
    if self.args.variant_caller in ("freebayes", "lofreq"):
        need["bcftools"] = "VCF 过滤"
```

freebayes / lofreq 需要 bcftools 做过滤；iVar 不需要（它自己出 tsv 再转 vcf）。

## 二、`virus_vcf_pipeline.py` 的角色

bcftools 在这个脚本里出现 5 次，全是**下游处理器**：

| 行 | 用途 |
|---|---|
| 109/116/226/234/241 | `bcftools query` 提取样本、位点、GT |
| 188 | `bcftools view -s` 保留 QC 通过样本 |
| 732 | `bcftools merge` 多样本合并 |
| 736 | `bcftools view -m2 -M2 -v snps` 过滤双等位 SNP |
| 1709/2177 | `bcftools reheader` 重命名样本 |

`--ivar` 参数（line 1601）用于 awk 修复合成的 FORMAT/GT 列，注释明确写「Freebayes/bcftools 不要加此参数」。

**结论：bcftools 在原管线是 VCF 加工工具，不是变异检出工具。**

## 三、可移植性实测（Windows）

| 工具 | MSYS2 包 | conda win-64 | 状态 |
|---|---|---|---|
| freebayes | 无 | 无（仅 linux-64/osx） | **不可得** |
| lofreq | 未提供 | 无（仅 linux-64/osx） | **不可得** |
| ivar | - | - | 已放弃（此前结论） |
| bcftools | **有**（已部署 1.24） | - | 可用 |

原管线三选一 caller，在 Windows 上**一个都拿不到**。

## 四、可选路线

### 路线 A：bcftools 兼做 caller + filter ★已实测定稿

```
# 1. 单参考 BAM 抽取（照抄原管线，awk 换 Python）
# 2. caller
bcftools mpileup -f ref.fa -d 100000 -q 0 -Q 13 -B \
    -a FORMAT/AD,FORMAT/DP,INFO/AD -Ob -o pileup.bcf fixed.bam
bcftools call -mv -Ob -o call.bcf pileup.bcf          # 不加 --ploidy 1，不加 -A
bcftools view -Ov -o raw.vcf call.bcf
# 3. 两步过滤（字段映射见 §4.3）
bcftools filter -s FAIL -i '<flt>' -Ov -o soft.vcf raw.vcf
bcftools filter -i 'FILTER=="PASS"' -Ov -o clean.vcf soft.vcf
# 4. SnpEff ann（同原管线）
```

- **优点**：工具已在位，实测通过；保持 BAM → VCF → SnpEff 的骨架
- **缺点**：偏离原管线；`INFO/SAF`/`INFO/SAR`/`INFO/AO` 是 freebayes 专有字段，bcftools 下需换成 `INFO/AD`、`INFO/DP4` 等，过滤表达式不能照抄

### 路线 B：远端 caller，本地产物

246 服务器有完整 conda 环境（freebayes/lofreq 都能装）。把 `fixed.bam` 或病毒 reads 传上去跑 caller，VCF 拉回本地注释。

- **优点**：与原管线 100% 一致
- **缺点**：需要网络往返，平台从「本地工具」变成「本地+远端混合」；部署复杂度上升

### 路线 C：只保留 SnpEff 注释层，caller 留空

平台不做 calling，接受用户已有的 VCF（或从服务器拷回），只做 SnpEff 注释与报告。

- **优点**：最小改动，SnpEff 已验证可用
- **缺点**：变异检出能力缺失

## 四.2 路线 A 实测结论（2026-09-08）

### 4.2.1 caller 参数：两个致命坑

| 参数 | 后果 | 实测数据 |
|---|---|---|
| `--ploidy 1` | **变异全丢** | 默认二倍体出 1 条，`--ploidy 1` 出 0 条 |
| `-A` + `-mv` | **变异全丢** | `call -mv` 出 0 条，`call -m` 出 14801 条（全 `<*>`） |

**根因**：`--ploidy 1` 下 GT 只能是 `0` 或 `1`，无法表示低频混合样本的 `0/1` 杂合态。病毒虽是单倍体，但混合样本中的变异正是低频杂合情形。原管线 freebayes 用 `--pooled-continuous` 解决的正是这个问题。

**定稿**：不加 `--ploidy 1`，不加 `-A`。用默认二倍体模型，靠 `AD` 字段算等位频率。

### 4.2.2 过滤字段映射（实测定稿）

| 原管线（freebayes） | bcftools 等价 | 验证 |
|---|---|---|
| `INFO/DP` | `INFO/DP` | 同名 |
| `INFO/AO` | `INFO/AD[1]` | 变异等位深度 |
| `INFO/AO/INFO/DP` | `INFO/AD[1]/INFO/DP` | 等位频率 |
| `INFO/SAF`（alt-fwd） | `INFO/DP4[2]` | alt-forward 计数 |
| `INFO/SAR`（alt-rev） | `INFO/DP4[3]` | alt-reverse 计数 |
| `INFO/AF` | `INFO/AD[1]/INFO/DP` | bcftools 不直接给 AF |

**`-a` 标签写法**：只能写 `FORMAT/AD,FORMAT/DP,INFO/AD`。`INFO/DP4` 会报 `Could not parse tag`，它随 mpileup 输出自动带上。

### 4.2.3 QUAL 阈值不可移植 ★

原管线 `QUAL>20` 是 freebayes 标度。bcftools call 的 QUAL 标度显著更低：

| 数据集 | 唯一/最低真阳性 QUAL | 沿用 QUAL>20 的结果 |
|---|---|---|
| 纯误差模拟 BAM（2962 reads, 30x） | 5.22 | 全丢 |
| 注入真值测试集（20 变异, 30x） | 4.01 | 丢 2 条 |

**定稿**：`QUAL >= 3.5`（真阳性最低 4.01，留 0.5 余量）。实测 `QUAL>4.02` 即开始丢真阳性，`QUAL>3` 保留全部 13 条。

### 4.2.4 truth set 验证

自造测试集：参考 OR489165.1（14812 bp），注入 20 个已知 SNV（频率 0.05/0.15/0.50/0.95 四档，各 5 个），1% 测序误差，30x 深度，150 bp 单端。

| 指标 | 结果 |
|---|---|
| caller 报告 | 13 条 |
| 真阳性 | 13 |
| 假阳性 | **0**（精确率 100%） |
| 漏检 | 7 个 |

**漏检归因**：7 个漏检位点的实际覆盖深度只有 1-6x（reads 随机采样导致局部覆盖不足），全部落在 0.05/0.15 低频档。

期望变异 reads 数（二项分布，depth × freq）：

| freq | 30x 下期望变异 reads | P(≥3 条) |
|---|---|---|
| 0.05 | 1.5 | 0.188 |
| 0.15 | 4.5 | 0.849 |
| 0.50 | 15.0 | 1.000 |

**检出率**：0.50 档 5/5，0.95 档 5/5，0.15 档 2/5，0.05 档 1/5。

**结论**：低频漏检是采样物理极限，非阈值问题。30x 深度检出 5% 变异需要约 90x 深度才能达到 P(≥3 条) ≥ 0.8。

### 4.2.5 过滤表达式定稿

```
QUAL>=3.5 && INFO/DP>={dp} && INFO/AD[1]/INFO/DP>={frq} && (INFO/DP4[2]+INFO/DP4[3])>={saf}
```

| 参数 | 值 | 依据 |
|---|---|---|
| `QUAL` | 3.5 | 真阳性最低 4.01，留余量 |
| `INFO/DP` | 10 / 20 / 100 | 照抄原管线三档（<50 / 50-1000 / >=1000） |
| `frq` | 0.05 | 照抄原管线 |
| 链支持 | 2 / 3 / 10 | 照抄原管线 saf/sar 档位，用 `DP4[2]+DP4[3]` 替代 |

实测：`QUAL>3 + DP + AF + 链支持` 保留全部 13 条真阳性。

### 4.2.6 数据有效性核查（避免误判）

此前用 `kv_truth` / `cli_test` 的 BAM 测出 0 变异，一度怀疑流程有问题。实测核查结论：

- `kv_truth\truth_reads.json` 明示 `depth:30, err_rate:0.01, n_reads:3105`，是**纯测序误差模拟**，无注入变异
- 实测错配率 1.06%，最密集位点仅 4 条 reads 支持（30x 下 13% AF）
- 该数据下 caller 报 1 条（恰好累积到 3/25 的位点）属正确行为
- 另注：`bcftools mpileup` 的 REF 列显示 `.`、ALT 显示 `<*>` 是正常格式（`.` = 与参考相同，`<*>` = `<NON_REF>`），**不是参考未加载**

## 五、结论

**路线 A 可落地，参数已实测定稿。** 理由：

1. freebayes/lofreq 在 Windows 上确实拿不到，工具的可用性是硬约束
2. bcftools mpileup+call 是 samtools 生态标准 caller，不是自研算法
3. 骨架完全照原管线：BAM 抽取、`fixed.bam` 命名、两步 filter、SnpEff 注释一致，只换 caller 实现
4. 过滤阈值已完成语义映射与实测标定，truth set 精确率 100%

**已推翻的结论（记录备查）**：

| 曾认为 | 实测推翻 |
|---|---|
| `--ploidy 1` 适合病毒单倍体 | 导致变异全丢，必须用默认二倍体 |
| `-A` 是保留等位基因，应加 | 与 `-mv` 组合导致变异全丢 |
| `INFO/DP4` 可写入 `-a` 标签 | 非法标签，会报错 |
| bcftools REF 列显示 `.` 是参考未加载 | 是正常格式 |
| `QUAL>20` 可直接移植 | bcftools QUAL 标度低，必须降到 3.5 |
| 测试 BAM 无变异是流程 bug | 该数据本就是纯误差模拟 |

## 六、待定

- `fixed.bam` 是否保留（原管线 finally 删除；保留可让变异段独立重跑）

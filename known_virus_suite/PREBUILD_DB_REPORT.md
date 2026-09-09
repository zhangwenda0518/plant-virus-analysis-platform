# 预建常用病毒库（B 策略）执行报告

日期：2026-09-08
执行：`kv_variant_test/prebuild_db.py`

## 一、清单口径

「常用」定义 = **枸杞属相关病毒 ∪ 平台历史样本高频检出病毒**

数据来源：

| 来源 | 方法 | 结果 |
|---|---|---|
| 平台样本 Kraken2 检出 | 扫 `results/*/02_virus_screen/virus_summary.tsv` 的 `species` 行（reads≥100） | 50 个 species |
| ref_info 匹配 | `Species_NCBI` / `Species_ICTV` / `VMR_Species` / `Virus name(s)` 四列索引 | **50/50 全命中** |
| 多段病毒展开 | 一个 species 对应多个 accession（segment） | 172 条 accession |
| 枸杞宿主补充 | `ref_info` 的 `Host` 含 Lycium | +2 条（NC_011803.1 枸杞病毒、MH791331.1） |
| **合计** | | **174 条** |

扫描样本目录 3 个（ERR7586041 / NX-6 / GQMIX），SRR39909438 无 `virus_summary.tsv`。

## 二、执行结果

```
目标 174 条  →  GenBank 下载 174/174（100%，约 3 分钟）
             →  SnpEff 建库成功 146 条（83.9%）
             →  失败 28 条（16.1%）
```

产物位置：

```
virus-db/snpeff_db/
├── gbk_files/                 174 个 GenBank 缓存（可复用，不必重下）
├── work/
│   ├── snpEff.config          建库配置
│   └── data/<GID>/
│       ├── genes.gbk
│       ├── sequence.bin
│       └── snpEffectPredictor.bin    ← 库本体（0.27s/条）
└── prebuild_summary.json      执行小结
```

## 三、28 条失败的三类根因（已逐条实测）

### A. 无 CDS：9 条

```
NC_001464  MZ463745  AM774357  AF298177  Y00328
NC_002030  NC_000885  PP101789  NC_003637
```

特征：GenBank 的 FEATURES 里**只有 `source`，没有任何 CDS**。

典型是类病毒（viroid，300~400 bp 裸 RNA，如 Pospiviroid exocortiscitri / Mexican papita viroid）。SnpEff 日志 `Exons created for 0 transcripts`。

**这不是缺陷**：类病毒本来就无蛋白编码基因，无功能注释可言。8 条 viroid + 1 条 PP101789（Myrica rubra tombus-like virus）。

### B. partial CDS：17 条

```
OR378271 OR378272 OR378274 AB971854 MG752888 MG752889 MH998222
MW246590 AY731821 OU862947 OU862948 OU862950 OU862951 OU862953
OU862954 OU862955 PV936563
```

特征：CDS 起点写作 `<1..3403`（5' 端不完整）。SnpEff 拿不完整序列翻译，长度非 3 倍数，日志 `Protein check: OK: 0  Errors: 1  Error percentage: 100.0%`。

这批都是**部分基因组提交**（partial genome），常见于早期提交或不完整组装。

### C. CDS 重叠：2 条

```
NC_008310  (Tobamovirus singaporense, 6485 bp)
NC_011803  (枸杞病毒, 6449 bp)
```

特征：4 个 CDS 且有重叠。以 NC_008310 为例：

```
CDS  59..4978    gp1  viral replication
CDS  59..3463    gp2  methyltransferase; helicase   ← 完全包含在 gp1 内
CDS  4968..5816  gp3  cell-to-cell movement          ← 与 gp1 重叠 10 bp
CDS  5803..6294  gp4  encapsidation
```

SnpEff 严格要求 CDS 不重叠（重叠会让变异归属产生歧义），日志 `Protein check: OK: 3  Errors: 1  Error percentage: 25.0%`。

**这是 SnpEff 的设计限制**，重叠 ORF 的病毒（Tobamovirus 的 readthrough、Potyvirus 的 PIPO 等）都会撞上。

## 四、已建成库的可用性验证

用 `OR489165`（Betacytorhabdovirus lycii，枸杞 crinkle 病毒）实测注释：

| 位点 | 变异 | 基因 | 效果 | 影响 | HGVS.p |
|---|---|---|---|---|---|
| 8000 | T>C | L | missense_variant | MODERATE | p.Ser147Pro |
| 8002 | T>A | L | synonymous_variant | LOW | p.Ser147Ser |
| 8005 | T>G | L | synonymous_variant | LOW | p.Leu148Leu |

结果与建库前的单样本测试完全一致，且额外给出 `downstream_gene_variant` 上下文。

## 五、覆盖率评估

按**检出 reads** 加权，看失败的 28 条影响了多少：

| 分类 | 条数 | 说明 |
|---|---|---|
| 无 CDS | 9 | 类病毒，本就无编码基因，不影响功能注释能力 |
| partial CDS | 17 | 部分基因组，CDS 不完整是数据本身的状态 |
| CDS 重叠 | 2 | SnpEff 设计限制 |

**对实际分析的影响**：失败条目里 reads 占比最高的是 NC_008310（220 reads）和 NC_011803（枸杞，低频）。主力高读数病毒（Coguvirus citrulli 542 万、Sobemovirus PHYRMV 406 万、Pittosporum cryptic virus-1 272 万）全部建库成功。

## 六、后续可做的改进（未执行）

1. **partial CDS**：把 `<1` 改写成 `1` 并补 `codon_start`，可强行建库，但注释会失真。**不建议**。
2. **CDS 重叠**：用 GFF3 只保留主 ORF 建库，会丢失重叠基因的注释。**按需权衡**。
3. **类病毒**：功能注释对 viroid 无意义，跳过是正确行为。

## 七、结论

- B 策略执行完毕：**174 条目标 → 146 条可用库**
- 失败 28 条全部定性为**数据本身特性 + SnpEff 设计限制**，没有实现缺陷
- 覆盖率足以支撑枸杞属病毒分析的主线（所有高读数病毒均建成）
- GenBank 缓存 174 个已留存，冷门病毒按需补建时可直接复用

## 附：失败清单完整文件

`kv_variant_test/build_fail_analysis.json`、`kv_variant_test/prebuild_plan.json`

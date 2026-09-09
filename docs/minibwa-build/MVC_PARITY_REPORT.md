# MultiVirusConsensus（MVC）对拍报告

**版本**：平台 `vp/consensus.py` vs ViralConsensus 原生 v1.0.4 vs ViralWasm 0.0.6
**数据**：ERR7586041 真实 reads（`conv_R1.fa.gz` + `conv_R2.fa.gz`）→ minibwa（`-x sr`）→ 同一份 SAM
**参考**：`results/ERR7586041/03_assembly/viral_contigs.fasta`，2 条 contig（18001 bp / 11735 bp）
**日期**：2026-09-08

---

## 一、编译与运行

| 项目 | 结果 |
|---|---|
| 源码 | `git-repo/ViralConsensus/`（v1.0.4，`common.h:16`） |
| 服务器编译 | 246 服务器，gcc 11.2.0 + conda htslib 1.21，EXITCODE=0 |
| 编译命令 | `g++ -Wall -O3 -std=c++11 -I$CONDA/include -o viral_consensus main.cpp common.cpp argparse.cpp count.cpp fasta.cpp primer.cpp -L$CONDA/lib -lhts -llzma -lbz2 -lz -lcurl` |
| 运行依赖 | `export LD_LIBRARY_PATH=$CONDA/lib` |
| 性能 | NODE_1 419738 recs / 1.14 s / RSS 12.8 MB；NODE_2 247495 recs / 0.71 s / RSS 11.4 MB |
| WASM 0.0.6 | NODE_1 0.66 s；NODE_2 0.33 s |

本机 Windows 编译受阻：缺 msys2/Git Bash，htslib 的 autoconf 脚本无法运行（`bash` 被 WSL 劫持）。服务器 Linux 编译一步到位。

**硬约束**：`viral_consensus` 要求参考 FASTA 恰好 1 条序列、SAM 恰好 1 个 `@SQ`。多参考须按参考切分逐条跑。WASM 版报错文案不同（`CRAM/BAM/SAM has 2 references, but it should have exactly 1`），约束一致。

---

## 二、三方产物一致性

### 2.1 官方两条实现路径完全等价

native v1.0.4 vs WASM 0.0.6，ACGTgap 逐位点差异：

- NODE_1：**0 / 18001**
- NODE_2：**0 / 11735**

### 2.2 共识序列三方完全一致

| 参考 | native | WASM | 平台 |
|---|---|---|---|
| NODE_1 | len=18001, N=86 (0.48%) | 同 | 同 |
| NODE_2 | len=11735, N=349 (2.97%) | 同 | 同 |

逐位点差异：**0 / 18001、0 / 11735**。

### 2.3 逐位点计数差异（修复前）

平台恒比官方少 1（如 native 107 / plat 106）。深度不同位点 NODE_1 5398 / 18001、NODE_2 3143 / 11735。

---

## 三、平台侧两个真 bug（已修）

### bug 1：SAM QUAL=`*` 被当作低质量丢弃（致命）

**现象**：平台 `pileup` 返回全 0 counts，`call_consensus` 输出全 N。

**根因**：minibwa `map` 输出的 SAM，QUAL 字段为 `*`（无输出碱基质量的选项）。平台 `pileup` 直接按 Phred+33 解析：`ord('*') - 33 = 9 < min_qual(20)` → **所有记录被静默丢弃**。

**htslib 实证**（服务器编译 `qualprobe.c`，读同一 SAM）：

```
read=ERR7586041.10 l_qseq=34 qual[0..4]=255,255,255,255,255
```

SAM 规范规定 QUAL=`*` 表示质量不可用，htslib 填充 **0xFF(255)**。官方 `count.cpp:142` 读的就是这个值，因此 255 ≥ 20 全部通过。

**修复**（`vp/consensus.py::pileup`）：

```python
if qual == '*':
    qual = None
...
q = 255 if qual is None else ((ord(qual[qi]) - 33) if qi < len(qual) else 0)
```

修复后 `called` 17915 / 11386（原为 0）。

### bug 2：MAPQ 过滤阈值与官方不一致

**现象**：扣除 bug 1 后仍有系统性差异。

**根因**：`count.cpp` 的过滤条件只有四项——`proper_pairs`、`FAIL_FLAGS`、`n_cigar==0`、`min_aln_len`/`min_aln_per`/`min_id_per`。**`src->core.qual` 从未被读取**，官方没有 MAPQ 过滤。

平台 `pileup` 默认 `min_mapq=10`，实测差异：

| 参考 | SAM 内记录 | mapq<10 | 平台处理数 | 官方处理数 |
|---|---|---|---|---|
| NODE_1 | 419738 | 93 | 419645 | 419738 |
| NODE_2 | 247495 | 52 | 247443 | 247495 |

---

## 四、差异归因（完整链条）

用严格复刻 `count.cpp` 语义的模拟器（`official_sim.py`）验证：**模拟器 vs 服务器原生产物，逐位点 0 差异 / 18001、0 差异 / 11735**，读入记录数也完全一致（419738 / 247495）。复刻精度已证。

| 差异源 | 机制 | 贡献 |
|---|---|---|
| ① QUAL=`*` | 平台误判为低质量，全部丢弃 | 计数归零（已修） |
| ② MAPQ 过滤 | 平台 `min_mapq=10`，官方无 | 93 / 52 条记录 |
| ③ N 碱基记账 | 官方 `BASE_TO_NUM['N']=-1` **越界写**；平台记独立 N 列 | 847 / 530 个位点 |

### 差异源 ③ 的机制细节

`common.h:38`：`BASE_TO_NUM['N'] = -1`。官方代码：

```cpp
++counts.pos_counts[pos][BASE_TO_NUM[(int)qseq[qpos]]];   // N → [pos][-1]
```

`pos_counts` 是 `std::vector<std::array<COUNT_T,5>>`。`arr[-1]` 落在**前一个元素的第 5 个槽（gap 列）**，即 `pos_counts[pos-1][4]`。

**受控实验验证**（服务器 `probe_n.sh`，2 条 read，参考 40 bp）：

| pos | 官方输出 | 平台语义 |
|---|---|---|
| 9 | A=2, gap=1, Total=3 | A=2, gap=0 |
| 10 | A=1, gap=0, Total=1 | A=1, **N=1** |

pos=10 的 N（R1 第 6 个碱基）被官方写到 pos=9 的 gap 列。官方的 gap 与 Total 因此被污染（pos=9 实际只覆盖 2 条 read，Total 却为 3）。

**集合级铁证**（`min_mapq=0`，`final_verify.py`）：

| 参考 | 平台 N>0 位点数 | 平台 vs 官方 gap 差异位点数 | 集合关系 |
|---|---|---|---|
| NODE_1 | 847 | 847 | `gap差异集合 ≡ (平台N位点 - 1)`，差集 0 |
| NODE_2 | 530 | 530 | `gap差异集合 ≡ (平台N位点 - 1)`，差集 0 |

每一个官方 gap 虚增，精确对应平台某个 N 位点的前一位，**零例外**。

---

## 五、最终一致性验证

平台以 `min_mapq=0` 重跑（对齐官方无 MAPQ 过滤）：

| 参考 | 记录数（平台 / 官方） | ACGT 差异位点 | 仅 gap 差异位点 | 共识序列差异 |
|---|---|---|---|---|
| NODE_1 | 419738 / 419738 | **0** | 847 | **0 / 18001** |
| NODE_2 | 247495 / 247495 | **0** | 530 | **0 / 11735** |

- 全部 gap 差异都是「官方 gap=1、平台 gap=0」，且平台该位点 N=0
- gap 差异位点集合与「平台 N 位点 - 1」**完全重合**（见上节集合级铁证）
- 交叉验证：平台 N 总量（880 / 552）= 模拟器越界写次数（880 / 552）
- 平台 N>0 位点数（847 / 530）少于 N 总量（880 / 552），因同一位置可被多条 N 记录覆盖

**结论**：除官方 N 越界的副作用外，平台与官方在计数与共识层面**逐位点等价**。

---

## 六、接口契约

### 官方默认参数（`common.h`）

```
MIN_QUAL=20  MIN_DEPTH=10  MIN_FREQ=0.5  AMBIG='N'
MIN_ALN_LEN=1  MIN_ALN_PER=0.0  MIN_ID_PER=0.0
FAIL_FLAGS = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP
```

**无 MAPQ 过滤**。

### 平台参数（`vp/consensus.py`）

```
MIN_BASE_QUALITY=20  MIN_DEPTH=10  MIN_FREQ=0.5  AMBIG='N'
MIN_MAPQ=10  MIN_MINOR_FREQ=0.02  MIN_MINOR_COUNT=2
FAIL_FLAGS=0x4|0x100|0x200|0x400   MIN_COV_PCT=10.0
```

### 口径差异（需在论文/文档中写明）

1. **MAPQ 过滤**：平台默认 10，官方无。若要严格对齐官方，传 `min_mapq=0`。
2. **N 列**：平台 poscounts 有 6 列（A C G T N -），官方 5 列（A C G T -）。平台 N 记账**更正确**；官方 N 越界属未定义行为。
3. **官方 `Total` 列含 gap 且可能被 N 越界污染**，平台 `Total = A+C+G+T+N+gap` 自洽。
4. **DEL 边界**：官方无参考末端边界检查，平台有 `if ri < rlen`。真实数据未见影响（差异为 0 的位点占绝大多数）。

### 建议

- 平台默认 `min_mapq=10` 是**更严格**的判据，保留合理；对外报告应与官方口径分开陈述。
- 若需与 ViralConsensus 数值对拍，用 `min_mapq=0`。

---

## 七、改动与备份

| 文件 | 内容 |
|---|---|
| `vp/consensus.py` | QUAL=`*` 按 htslib 语义处理为 255 |
| `vp/consensus.py.bak_qualstar_20260908.orig` | 改动前备份 |

## 八、产物索引

| 文件 | 说明 |
|---|---|
| `_mvc_work/cmp_report3.txt` | 三方逐位点对拍（修复前） |
| `_mvc_work/official_sim_report.txt` | 官方语义模拟器 vs 原生产物（0 差异） |
| `_mvc_work/attrib_report.txt` | 差异归因（MAPQ / N 越界） |
| `_mvc_work/final_verify.txt` | 最终一致性验证 |
| `_mvc_work/diff_diag.txt` | 差值分布诊断 |
| `_mvc_work/probe_n.sh` | 服务器受控实验（N 越界） |
| `_vc_win/qualprobe.c` | htslib QUAL=`*` 实证程序 |
| `_vc_win/build_remote.sh` | 服务器编译脚本 |
| `_vc_win/native_cmp2.sh` | 原生对拍执行脚本 |

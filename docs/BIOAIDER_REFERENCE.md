# BioAider 功能调研与借鉴清单

> 调研对象：BioAider V1.423（Windows 桌面版）+ **V1.727（Linux 版源码解析，
> 见文末"四、V1.727 源码解析"）**。本文档提炼其功能面，对照本平台给出
> 借鉴结论与优先级。反编译参考源码归档于 `docs/bioaider_ref_src/`（仅供
> 算法阅读，不直接复制，版权与引用要求见其 README）。

## 一、BioAider 功能面（来自手册目录 V1.423）

**3.1 SeqTools（序列工具）**
- 3.1.1 Seqformat Convertor：多格式互转（fas/nexus/phylip/meg/paml 等）
- 3.1.2 SeqVary：比对序列差异展示
- 3.1.3 SequenceID Rename：批量重命名序列 ID（映射表）
- 3.1.4 Split Sequence Fragment：序列/比对拆分
- 3.1.5 Combine Gene (Tandem Gene)：多基因串联 + 分区文件（partition）
- 3.1.6 Visual Gene Extractor：按坐标/注释提取基因区段
- 3.1.7 Fast Annotation：参考注释快速映射到序列
- 3.1.9 Correction ambiguous bases：模糊碱基（IUPAC）校正

**3.2 Similar Analysis（相似性）**
- 3.2.1 Sequence Identity Matrix（同 SDT 口径）
- 3.2.2 / 3.2.3 Remove High-Similar / Delete Low-Similar：参考集去冗余/
  筛减（按相似度阈值增删序列）
- 3.2.4 Repeat Fragment Search：重复片段搜索

**3.3 Mutation Tools（突变）**
- 3.3.1 Mutation Analysis：相对参考的突变谱分析
- 3.3.2 Site Counter：位点计数/频率统计
- 3.3.3 Site Screen：位点筛选
- （含密码子级同义/非同义统计 → dN/dS 方向）

**3.4 Drawing（作图）**
- 3.4.1 基因突变棒棒糖图（lollipop chart）
- 3.4.2 常用统计图（频率分布等）

**4. 插件体系**：mafft/muscle/clustal-omega/iqtree/mrbayes/hmmer 等外部
工具以插件方式注册（plug_set.init：名称/状态/描述/路径），可"Add"安装。

## 二、与本平台对照

| BioAider 功能 | 本平台现状 | 结论 |
|---|---|---|
| 下载/溯源/组装/宿主预测/RdRP/报告 | BioAider **没有** | 本平台独有，不受影响 |
| Identity Matrix / 序列差异展示 | SDT 口径矩阵 + MSA SNP 查看器 | 已覆盖，形态更现代 |
| 外部工具管理 | platform.json + 自动探测（22 工具全内置） | 已覆盖（我们还全预装） |
| 格式转换 | 仅 FASTQ↔FASTA | **缺口**：缺 nexus/phylip 等系统发育格式 |
| ID 批量重命名 | 无 | **缺口**（小工具） |
| 基因串联+分区文件 | 无 | **缺口**（分段基因组建树刚需） |
| 模糊碱基校正 | 无 | **缺口**（小工具，序列清洗） |
| 参考集去冗余/筛减 | 仅建树内部抽样（blast/macro/genus/lineage） | **半缺口**：缺独立"按 identity 阈值筛减参考集"工具 |
| 突变谱 + 棒棒糖图 + 位点统计 | 无（路线图"变异分析"项） | **缺口**（BioAider 的呈现形态可参考） |
| 序列级小工具"即传即跑" | 专项分析页已是此形态 | 架构一致，新工具往专项分析分区挂即可 |

## 三、借鉴结论（按优先级）

1. **突变谱 + 棒棒糖图**（并入路线图"变异分析"）：对齐参考 → 突变表
   （位点/参考碱基/替代/频率）→ 交互棒棒糖图。植物病毒场景做成
   "分离物 vs 参考基因组"口径；密码子级同义/非同义一并输出（dN/dS 前置）。
2. **基因串联 + 分区文件**：多分段基因组（双生病毒/呼肠孤/番茄斑萎病毒属）
   多基因座建树刚需；输出 fas + NEXUS/RAxML 分区文件，直接喂 IQ-TREE `-p`。
3. **参考集去冗余**：identity 阈值一键筛减（Remove High-Similar）——
   我们 ⑦ 建树的参考池/NCBI 集合常含高相似冗余，先去冗余再建树更干净。
4. **序列工具小分区**（专项分析新增"🧬 序列工具"组）：ID 批量重命名
   （含映射表导出）、模糊碱基校正、格式转换扩展（输出 nexus/phylip）。
   三个都是百行级小工具，共用 iter_fasta 底座。
5. **不必借鉴**：Qt 桌面架构（Web 已更优）、外部工具插件化安装
   （我们全内置且自动探测，比它的 Not-installed→Add 模式更省心）。

## 四、架构层面的两点观察

- BioAider 把"序列级小工具"与"流程类分析"分层清晰；本平台的
  专项分析页已有同等分层，新小工具照该分区扩展即可，不必新开页面。
- 其配置（plug_set.init 等）是纯文本 TSV/init，我们 platform.json +
  设置中心的方案信息量更大（默认参数/语言/路径一体），维持现状。

## 四、V1.727 Linux 版源码解析（2026-09-05 深夜）

**解析链路**（比 Windows 版好解析得多，已全链路打通）：
PyInstaller ELF → pyinstxtractor-ng（PyPI 安装；自动处理加密——PYZ 为
AES-CTR 加密，密钥 `pyimod00_crypto_key` 模块内明文可见）→ 3.7 pyc →
uncompyle6（需 xdis==6.x，7.x 的 load_module 返回 8 元组会崩；且不能在
提取目录里运行 python——3.7 的 inspect.pyc 会污染 stdlib 导入）。
反编译源码归档：`docs/bioaider_ref_src/`（12 个核心模块，6778 行）。

**关键算法细节（实现借鉴要点）**：

1. **突变分析**（mutation_analysis + mua_fuc）：
   - 参考序列 = FASTA 第一条；其余序列逐条与之做逐位点比对
     （multiprocessing Pool 分块流式读，piecewise_read 控内存）；
   - 每序列输出 [(参考碱基, 位点, 替代碱基)]；两套计数：全部突变 vs
     "精选"（排除 gap 与简并碱基，identify_degenerate_bases 判 IUPAC）；
   - 输出：逐序列突变日志 TSV、位点汇总 TSV（位点/突变/序列数）、
     过滤版汇总、位点频率分布图；密码子模式输出同义/非同义统计
     （带阈值筛选 sy_nonsy_select_value）+ codon_plot；
   - 棒棒糖图数据函数 lollipop_data 也在 mua_fuc（matplotlib 绘制）。
2. **分区模型文件**（partition_make，注意：它只生成分区文件，不含串联）：
   - IQ-TREE：`#nexus; begin sets; charset <名> = <起>-<止>; ...
     charpartition mymodels = <模型>:<名>, ...; end;`
   - MrBayes：charset + lset applyto（model_resolution 把 GTR+G 等拆成
     nucmodel/ntype/rates/ngamma 等参数）+ partition 声明 +
     prset ratepr=variable + unlink statefreq/revmat/shape/pinvar/tratio。
   - 支持核酸/蛋白两种数据类型。
3. **格式转换**（public_functions）：fas→nex/phy/pml 及反向，各 ~30 行
   （nex 带 `begin data; dimensions nchar= ntax=; charstatelabels` 头）。
4. **k-mer 矩阵**（clustering_tools.get_kmer）：k-mer 计数表 + linkage
   聚类（与我们的"参考集去冗余"可合并设计）。
5. **NCBI 批量取序列**：其实现远薄于我们的 ncbi_download（无会话翻页/
   白名单），不必借鉴。
6. 反编译伪影提示：uncompyle6 会把切片 `seq[site:site+1]` 错写成
   `seq[site[:site + 1]]`——阅读时注意甄别，实现时按语义写。

**对借鉴优先级的修正**：之前把"基因串联+分区"列为缺口——V1.727 里
partition_make 只生成"分区文件"（不含串联）；串联（Combine Gene）逻辑在
主程序内未解出，但串联本身是平凡操作（多比对按坐标拼接 + 记录各区段
边界）。实现我们的"基因串联+分区"时：串联自己做，分区文件输出直接按
上面第 2 条的两种格式写。

## 五、模块展示结构对比（2026-09-05 深夜，main.pyc 反编译）

### BioAider 的展示结构（main.py 1117 行，菜单树全解）

```
主窗口 = 中央文件浏览器（QFileSystemModel 目录树 + 信息面板）
        + 底部工具条（引用/帮助/示例/序列查看器拖拽钮）
        + 菜单栏（全部功能入口，两级菜单树）：
  Toolbox        SeqTools        格式转换/SeqVary/重命名/拆分/串联/基因提取/
                                 GB解析/快注释/模糊碱基/顺序提取
                 Similar Analysis  identity矩阵/去高相似/删低相似/重复片段/蛋白定位
                 Download Tools    NCBI 批量取序列
                 Conversion Tools  日期转换/文件合并
  Align          mafft/muscle/clustal/blast + HMMER(hmmbuild/hmmscan/hmmsearch)
  Evolution      Molecular Variation（突变分析/位点计数/位点筛选）
                 Phylogeny（modelfinder/fasttree/iqtree/mrbayes/分区文件）
  Clustering     k-mer 矩阵/层级聚类
  Drawing        棒棒糖图/小提琴图
  Setting        10 种主题/字体/检查更新/插件管理
  Shortcut       外部软件注册
  About
+ "常用工具"快捷菜单：快注释/突变分析/位点计数/矩阵/去高相似/blast/mafft/iqtree/mrbayes
每个功能 = 独立弹窗：拖放输入(Drop_QPushButton) + 参数表格(TableWithCopy)
         + QThread 运行(进度 Signal) + 结果目录打开
```

### 两种范式对比

| 维度 | BioAider（工具导向） | 本平台（数据流导向） |
|---|---|---|
| 组织单位 | 功能=菜单项→独立弹窗 | 页面=样品/结果实体 |
| 主界面 | 文件浏览器（拖文件找工具） | 工作流页面（样品从建到报告） |
| 功能可见性 | 菜单树一眼看全，但两层深 | 专项分析长滚动，功能多时发现性差 |
| 高频入口 | "常用工具"快捷菜单 | 无（总览只有工作流入口） |
| 工具闭环 | 单窗内：参数→运行→进度→开结果 | 工具页运行后看页面底部任务区 |
| 文件输入 | 拖拽到窗口 | 逐个 📁 浏览选择 |
| 主题 | 10 套即时切换 | 单主题 |
| 序列速览 | 拖拽 fasta 即弹查看器 | 无轻量查看器（只有 MSA/树查看器） |
| 流程/项目/队列 | **无** | 样品+13 阶段+队列+项目分组（我们独有） |

结论：BioAider 是"瑞士军刀"（找工具→跑文件），本平台是"流水线"
（样品进报告出）。两种范式服务不同时刻，**不冲突**——我们的专项分析页
本质上就是它的 Toolbox，展示方式可以借它四点长处。

### 借鉴落地项（展示结构，按建议顺序）

1. **专项分析页加"功能目录/锚点侧栏"**：左侧粘性目录列出组与工具
   （锚点跳转），解决长滚动发现性问题——等价它的菜单树，Web 化表达。
2. **常用工具快捷组**：工具卡支持 ⭐ 收藏（localStorage），总览页
   "分析工作流"上方出现"常用工具"条——等价它的 common_tool。
3. **输入框支持拖拽文件**：专项分析所有文件输入框接受从资源管理器
   拖入（Web drop 事件取 path 填入输入框）+ 数据下载页整批拖入文件
   补路径——等价它的 Drop_QPushButton 体验。
4. **工具卡内嵌运行状态**：专项分析每个工具卡内嵌"运行中/进度条/日志
   尾部 + 取消"（把分析流程页阶段卡的 stage-log 机制推广到工具卡），
   免去跑起来后滚到页面底部看任务。
5. **深色主题**：CSS 变量已就绪，加一套 dark 配色 + 设置页切换。
6. **轻量序列查看器**：选/拖一个 fasta → 分页表格展示 ID/长度/预览
   序列 + GC/简并碱基统计——补齐它的 Seq_Viewer。
（不建议借鉴：菜单树本体——页面+锚点目录在 Web 上更优；文件浏览器
主界面——我们以样品/结果为实体更有序。）

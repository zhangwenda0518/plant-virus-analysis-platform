# 植物病毒分析平台（Windows 本地版）

基于 **kunpeng**（超低内存宏基因组分类器）的本地植物病毒诊断与深度分析平台。
鼠标点击即可完成：宿主去除 → 病毒筛查 → 组装 → ORF 预测 → ORF 功能注释 → 进化树/SDT → 引物设计 → 基因组图 → 可视化报告。

## 功能总览

| 阶段 | 工具 | 输出 |
|------|------|------|
| ⓪ Fastp 质控（可选） | fastp | clean reads、质控报告 |
| ⓪b FASTQ→FASTA 转换（可选，推荐） | seqkit fq2fa | conv_R1/R2.fa.gz（分类提速） |
| ① 宿主去除 | kunpeng（宿主库） | kept_R1/R2.fastq.gz、宿主占比统计 |
| ② 病毒筛查 | kunpeng（病毒库） | kreport2、种/属/科汇总表 |
| ③ 组装与分类 | SPAdes（metaviral/rna/meta）+ kunpeng + BLAST | contigs、病毒 contig 注释表 |
| ④ 宿主预测 | ICTV 宿主概率级联 | host_prediction.tsv、桑基/旭日图 |
| ⑥ ORF 预测 | orfipy（getORF 等价）+ pyrodigal + pyrodigal-rv | pep/nt/bed/gff/faa/ffn |
| ⑥b ORF 功能注释 | DIAMOND/MMseqs2/blastp × RefSeq 病毒蛋白库 + 类别投票 + ICTV 映射 | orf_annotation.tsv、功能类别/科分布、注释 GFF3 |
| ⑦ 进化分析 | MAFFT + trimAl 清剪 + FastTree / IQ-TREE（UFBoot+SH-aLRT 双支持值）；SDT 口径 identity 矩阵；可选追加 NCBI 下载参考 | tree.nwk、树图、sdt_matrix.csv + 热图 |
| ⑧ 引物设计 | primer3（保守区 / 全长两种模式） | primers.tsv |
| ⑨ 基因组图 | gbdraw 首选（SVG 圈图 + 线图）；缺则用 dna_features_viewer（DFV，按功能类别着色） | 09_genome_plots/*.svg |
| ⑩ 可视化报告 | plotly + pycirclize + matplotlib + gbdraw/DFV | report.html |

**进度与预估**：任务卡实时显示全局进度条、已运行时长与**预计剩余时间**
（平台按每阶段历史耗时自学习，越用越准）；kunpeng 分类以 chunk 中间盘
增长近似上报进度，SPAdes 按 k-mer 阶段上报。每阶段开始/结束在日志中
输出**资源预估**（线程 / 内存 / 磁盘）与**实际耗时**，分类阶段结束额外
输出"预估 vs 实测"磁盘核对行。

## 界面（按分析工作流布局）

顶部导航即工作流顺序：**总览 · 分析流程 · 结果中心 · 专项分析 · 数据下载 ·
公共数据检索 · 数据库 · 公共病毒组 · LOGAN 溯源 · 提交准备 · 设置 · 手册**
（🌐 一键切换中/英文；📂 打开平台目录）。

- **总览**：数据库状态卡片 + 分析工作流入口（按 建库→分析→结果→专项→溯源 排序）
  + 最近样品
- **分析流程**（原"分析管道"）：新建样品 + 按功能模块分组的 13 步级联流程（详见下文）
- **结果中心**（原"结果"）：已分析样品（报告/目录/文件/删除）+ **专项结果入口**
  （共线性比较图、LOGAN 溯源报告在此聚合）+ SDT / MSA / 进化树交互查看器
- **专项分析**（原"工具箱"）：不分样品流程的独立分析，分两区——
  🧬 单步工具（质控 / 病毒鉴定提取 / 组装 / contig 分类与深度分析四件套）·
  🧬 参考数据与比较分析（NCBI 参考下载、同属共线性比较）
- **数据与库**（原"数据库"）：Taxonomy / 宿主库构建（病毒库预置）
- **设置**：界面语言、全局线程数、NCBI 邮箱、**17 项分析默认参数**，保存即生效；
  分析流程页各步骤卡片初值自动带出设置页默认，仍可逐步覆盖

## 设置中心 · 中英双语 · 任务结果预览

- **中英双语**：导航栏 🌐 一键切换（或到「设置」页切换），前端界面与后端阶段名/
  卡片摘要同步切换；偏好写入 platform.json，重开页面保持。
- **存储与磁盘**（设置 → 存储与磁盘）：盘剩余容量条 + 各目录占用明细
  （databases 按子库拆分，一眼看到哪个库占了多少）；目录统计为后台任务 +
  10 分钟缓存，不拖慢页面。平台盘剩余 < 20GB 时**全站顶部出现黄色水位横幅**。
- **并发闸门**：重任务（kunpeng 分类 / SPAdes / DIAMOND 注释 / 建库 / 建树 /
  SDT）同时最多 2 个，轻任务（下载 / 转换 / 绘图 / 检索）最多 4 个；超出自动
  排队并在任务卡显示「排队中 · 前面还有 N 个任务」。可在 platform.json
  `defaults` 设 `max_heavy_tasks` / `max_light_tasks`（1-8 / 1-16）调整。
- **默认参数**：「设置 → 分析默认参数」集中管理（confidence/subsample/组装模式/
  组装输入/内存/最小 contig/最小 ORF/参考数/建树工具/抽样策略/引物模式/特异性
  检查/出图上限/绘图引擎/强制重跑等）；管道页与 CLI 未显式传参时自动采用。
- **输入输出开放**：文件浏览对话框记忆上次目录；所有输入框均可手填相对/绝对路径
  （建库/分析输入支持平台外绝对路径读取，产物一律限定平台目录内）。
- **任务结果预览**：每个任务（样品分析 / 工具 / 下载 / 比较 / LOGAN）结束后，
  任务卡自动展开**结果预览面板**——样品任务显示各阶段关键数字 + 「打开报告」；
  工具任务显示统计值 + 产物文件清单（可直接下载）；同时弹出自消失 toast 提示，
  页面在后台时发系统通知。
- **阶段产物在线预览**：分析管道 9 个阶段卡片（fastp 质控报告 / ①宿主 / ②筛查 /
  ③组装 / ④宿主桑基图 / ⑥b ORF 注释表 / ⑧引物表 / ⑩主报告）
  带「📊 查看」按钮——表格就地弹窗浏览（含下载），HTML 报告新窗口打开。
- **样品管理**：「结果」页可删除样品（含全部产物；有任务运行时拒绝删除）。
- **公共数据下载**（导航"数据下载"）：粘贴 SRR/ERR/DRR（ENA 直下 FASTQ.GZ）、
  CRR（NGDC/GSA）或数据文件 URL → aria2c 多线程批量下载（断点续传、字节数
  +gzip 完整性校验、NGDC 低并发防封控）→ 一键转入分析流程（自动配对 R1/R2
  建样品）；独立进度窗口可弹出挂机观察。GSA 的 `.sra` 文件由**内置 sracha**
  （纯 Rust 引擎，比 fasterq-dump 快 5-13 倍，Windows 版已随平台编译好）
  自动转换 FASTQ.GZ，无需安装 sra-tools。
- **项目管理与批处理队列**：分析流程页"批量导入"卡（TSV：样品名/R1/R2 +
  项目名）批量建样并入队；队列顺序执行（同时只跑一个样品），实时进度面板；
  样品列表按项目标签筛选。

## 快速开始（图形界面）

1. 双击 **`启动平台.bat`** → 浏览器自动打开平台页面（默认 `http://127.0.0.1:8900`；
   若默认端口被系统占用/保留会自动改用其他端口，以黑窗口里打印的地址为准）
2. 进入 **数据库** 页，按顺序完成（仅首次）：
   - ① 下载 NCBI Taxonomy（~57MB）
   - ② 宿主库：选择任意物种的基因组 FASTA（NCBI 下载、自行组装均可），
     填入宿主物种的 **NCBI TaxID**（如 4081 烟草、112863 枸杞），点击构建。
     TaxID 自动校验：不存在/已被合并（如旧编号并回新编号）会在建库前直接提示，
     不会白跑。序列名自动规范化：任意格式的 FASTA 头（含空格描述、特殊字符、
     重复/空 ID、旧 taxid 标签）均可，平台统一清理并注入标签。
   - 病毒库无需构建：平台已预置（`databases/virus_db`）
3. 打开 **分析管道**（主页「总览」→「分析管道」，或导航栏直达）：
   - **新建样品**：填样品名、选 R1/R2（点击 📁 在平台目录内浏览；单端数据 R2 留空），
     可勾选"截取子样本"先小规模验证；
   - 管道按**功能模块分组**显示（PhyloSuite 风格）：🧹测序数据预处理
     （⓪Fastp质控·可选 → ⓪b FASTQ→FASTA转换·可选 → ①宿主去除）→
     🦠病毒鉴定（②筛查与提取）→ 🧬病毒组装（③组装·分类·病毒contigs提取）→
     🧲宿主预测（④ICTV宿主判定）→
     🔬下游分析（⑥ORF → ⑥b ORF功能注释 → ⑦进化树与SDT → ⑧引物设计 → ⑨基因组图gbdraw）→
     📄报告（⑩可视化），
     **每步的输出自动作为下一步输入**，卡片上标注了摘要（宿主占比/病毒种数/contigs/
     ORF数等）与产物文件（⑩报告卡片带「查看结果」按钮）；
   - **可选步骤自动降级**：fastp（质控）、seqkit（fq2fa 转换）未安装时对应卡片
     灰显"可跳过"、全流程自动绕过；gbdraw 与 dna_features_viewer（DFV）都
     未安装时 ⑨ 基因组图才灰显（安装命令见 requirements.txt 注释）；
   - 每个模块卡片内直接暴露该步参数（③组装模式/组装输入/内存/最小contig、
     ⑥最小ORF、⑦参考数/建树工具、⑧引物模式/特异性检查、⑨出图上限/自备
     fasta+gff/gb），页面右上网格为全局参数（线程数/分类置信度/chunk目录/强制重跑）；
   - 操作：**▶ 运行此步**（单跑该模块）/ **⏩ 依次运行到此**（从第一个未完成步骤连跑）/
     **↻ 重跑此步**（配合"强制重跑"忽略断点）/ 顶部 **▶ 依次运行剩余步骤**；
     上游完成后下游自动解锁，任务结束页面自动刷新状态。
   - 组装说明：metaviral 模式低深度样品可能拼不出 contigs，平台自动改用 rna 模式重试；
     单端数据 metaviral/meta 自动转 rna（SPAdes 宏基因组模式不使用单端 reads）。
4. 完成后到 **分析结果** 页打开报告

> 小提示：首次跑大样品前，可把「子采样」设为 100000（10 万对 reads）先快速验证全流程。

## 命令行（CLI）用法

```bat
python main.py tools                                   :: 查看工具探测状态
python main.py selfcheck                               :: 环境自检（模块/工具/数据库/磁盘）
python main.py init-taxonomy                           :: 下载 NCBI taxonomy
python main.py build-host-db --genome host-db\genome.fa --taxid 4081
python main.py build-virus-db --fasta virus-db\final.cluster.ref.fasta --info virus-db\final.cluster.ref_info.tsv
python main.py analyze --r1 fastq\NX-5_S2_L001_R1_001.fastq.gz --r2 fastq\NX-5_S2_L001_R2_001.fastq.gz --sample NX-5
python main.py analyze --r1 ... --r2 ... --sample T1 --stages virus,orf --subsample 200000
python main.py report --sample NX-5                    :: 重新生成可视化报告
python main.py orfa --sample NX-5                      :: 只跑 ⑥b ORF 功能注释
```

分析阶段可自由组合：`subsample,fastp,fq2fa,host,virus,assembly,hostana,orf,orfa,phylo,primer,gbdraw,report`。
每阶段有断点标记（`.done` 文件），中断后重跑自动跳过已完成阶段（`--force` 强制重跑）。
常用开关：`--no-fq2fa` 关闭预转换；`--gbdraw-max 20` 调整出图条数；
`--gbdraw-fasta x.fa --gbdraw-ann x.gff`（或 `--gbdraw-ann x.gb`）用自备文件出图。

## 公共数据检索（meta_search，源自 MMPV-RNA public_metadata_pipeline）

按物种检索 NCBI SRA + CNCB GSA 公共测序数据，提取统一元数据，勾选 Run
一键转入「数据下载」→ 分析流程，形成"检索 → 下载 → 分析"闭环。

```bat
python main.py meta-search --species "Lycium chinense"               :: SRA+GSA 双引擎检索（详细模式）
python main.py meta-search --species "Lycium chinense" --db sra --source TRANSCRIPTOMIC --no-detailed
python main.py meta-info --runs meta_search\Lycium_chinense\search\sra.list --fill-date
python main.py meta-info --runs SRR24100141                          :: 单个 Run 直接提取
python main.py meta-plot --input meta_search\...\SRA_GSA_Merged_Final.csv   :: SCI 级元数据可视化
python main.py host-genome --species "Lycium barbarum" --include-organelles :: 宿主基因组下载（供建宿主库）
```

- 产物在 `meta_search/<物种>/`：`search/SRA_GSA_Merged_Final.csv`（检索表 +
  download_links.csv 直链）、`info/Global_Unified_Metadata_Core14|Full.csv/tsv`
  （14 列核心表 + 34 列全维表）、`plot/`（时间/机构/组织/地理出版级图组）。
- 详细模式解析 Tissue/Location/发育阶段；可选 `--deepseek-api` 做 AI 元数据
  清洗（Web 端未配 Key 时自动用纯规则模式）。
- 下载环节复用平台「数据下载」（ENA/NGDC 直下 + aria2c + fasterq-dump），
  不引入原管线的 prefetch 下载链路；建宿主库仍用平台 `build-host-db`
  （kunpeng），`host-genome` 只负责把宿主基因组 FASTA 下载到位。

## NCBI 提交准备（submissions，源自 MMPV-RNA submission_gui / virome_submission_pipeline）

以 `unified_metadata.csv` 一张表驱动 GenBank + BioSample 提交准备：在线编辑
（双击单元格、占位符橙色高亮、Ctrl/Shift 多选删行、点表头列排序、每列填充
进度）、批量填充/快速填充、必填字段校验，一键生成 source.src / miuvig.tsv /
assembly.tsv / BioSample 注册表 / template.sbt / 交互式提交报告（report.html），
最终走 NCBI BankIt 网页向导或 Sequin 桌面程序。原 GUI 其余功能面同样齐备：
内置示例（6 行 Demo / 公共数据样例 / 自测样例）、项目另存为、提交文件预览
在线编辑（.sqn 二进制锁定）、BioSample 导出可选含占位符行、物种名 NCBI
Taxonomy 在线校验、序列 FASTA 提取（sequences.fsa）与整包 zip 下载。

```bat
python main.py submit-list                                    :: 列出提交项目
python main.py submit-init --name nx6 --taxonomy taxonomy.tsv --authors "Zhang, Wenda" --title "..."
python main.py submit-init --name demo --demo                 :: 示例数据建表（快速上手）
python main.py submit-fill --name nx6 --column bioproject --value PRJNA123456
python main.py submit-validate --name nx6                     :: 必填字段校验
python main.py submit-export --name nx6                       :: 生成全部提交产物
python main.py submit-sbt --name nx6 --last Zhang --first Wenda --affil "Ningxia University" --city Yinchuan --country China --email me@x.com
```

- 产物在 `submissions/<项目名>/`：unified_metadata.csv（中枢表，Web 端在线
  编辑）、source.src（GenBank source modifiers，含 source_individual/ 按病毒
  拆分）、biosample_template.tsv（BioSample 批量注册）、miuvig.tsv/assembly.tsv
  （结构化注释）、authorset/template.sbt（作者 ASN.1）、report.html（可编辑
  提交报告）、validation_report.txt。
- Windows 适配说明：不引入原管线的 suvtk/tbl2asn（Linux 依赖）；特征表
  (.tbl) 与 Sequin 包构建不在本模块范围——用 BankIt 网页向导提交
  source.src + 序列即可；若后续需要 .sqn，可在装有 suvtk 的环境补跑。

## 病毒参考库（virus_ref / acvirus_db）

```bat
python main.py ref-status                              :: 查看参考库版本与统计
python main.py analyze ... --tree-sampling lineage     :: ⑦建树参考按分类抽样
```

- `databases/virus_ref/`：非冗余植物病毒参考库（8,465 条 98% ANI 聚类代表 +
  5,773 条完整基因组子集 + ~199K 全量 + DATA_VERSION 版本锚点），源自
  plant_virus_db_pipeline。④⑦步分析直接取用：
  ④ 宿主预测用其 ICTV 谱系 + 已知宿主交叉验证（预测宿主 ∉ 已知宿主类别
  → `host_check=WARN`，报告标黄）；⑦ 参考池优先 RefSeq/完整基因组。
  `ref_meta.tsv` 为规范化元数据缓存（首次使用自动构建：Segment 归一 +
  本地 taxonomy 谱系补全，ICTV 种覆盖 10%→100%）。
- `databases/acvirus_db/`：全病毒参考库（来自 246 服务器，20,178 条完整
  基因组 + ICTV 全谱系 taxa.txt，不止植物）。⑦步自动把病毒 contigs 再对比
  该库，非植物病毒/植物库缺代表时也能选到近缘参考；每组输出 `refs.tsv`
  标注参考来源（plant_ref/acvirus）与谱系。
- `databases/ictv_db/`：ICTV VMR 参考库（vp/ictv_db.py）——官方
  [VMR 当前版 xlsx](https://ictv.global/vmr/current?fid=15873) 解析出的
  全病毒分类元数据（MSL41: 22,785 条 accession / 4,068 属 / 393 科）+
  按需下载的参考序列缓存。定位：谱系比 acvirus 旧表更新（MSL 换版的
  种改名/科重分类即时生效），序列层"本地优先、缺了才下"：

```bat
python main.py ictv-status                              :: 库状态（MSL 版本/覆盖）
python main.py ictv-update                              :: 在线下载 VMR → 解析 taxa.txt
python main.py ictv-update --xlsx VMR_MSL42.xlsx        :: 手动放入的 xlsx 解析
python main.py ictv-refs --genus Tobamovirus            :: 选参预览（本地优先）
python main.py ictv-refs --genus Nepovirus --download   :: 缺的 accession 按需下载
```

  分层：①解析层 VMR xlsx → `taxa.txt`（沿用 acvirus 谱系列 + 病毒名/
  基因组完整性/Baltimore/宿主组扩展列，多 accession 分段行拆分）+
  `DATA_VERSION` 锚点；②选参层按属/科/种过滤，本地已有（acvirus_db
  fasta 命中 > gb 缓存）优先、Complete genome 优先；③下载层缺的
  accession 经 NCBI efetch 落 `gb_cache/<acc>.gb`（幂等缓存）并汇总
  `gb_refs.fa`；④兜底：⑦步只读接入（谱系并入分类索引、缓存序列并入
  参考池），不在分析中联网，未下载时自动回退 acvirus_db/植物库。
  下载的 .gb 保留特征注释，可同时供同属共线性比较（synteny）取用。
- ⑦步 `--tree-sampling`：`blast`=按比对 hits（默认）；`macro`=同科建树
  （目标属 + 同科各属 3 条背景）；`genus`=属级树；`lineage`=种级树
  （移植 246 服务器 acvirus_tree_pro 的层级抽样策略）。
- 参考库不做在线自动更新：需要更新时，把新版 4 件套（两个 FASTA +
  两个 Info.tsv + DATA_VERSION）从 pipeline 服务器/246 拷入
  `databases/virus_ref/` 整体替换即可，平台按 DATA_VERSION 自动重建
  元数据缓存与 BLAST 库。

## 宿主预测模块（④，ICTV 级联）

吸收自 MMPV-RNA virome_discovery_pipeline 的 C9 方法，对 ③组装 产出的
病毒 contigs 判定感染宿主类别：

1. contig → **kunpeng 对 contigs 的分类判定**（kunpeng_taxid，kraken2 同构输出
   的 C 行）→ 病毒库 info.tsv 同 taxid 的 ICTV 分类（种/属/科）；
   LCA 判到上级节点（属/科级）时经 nodes/names.dmp 谱系解析还原完整分类
   （学名优先，带跨样品缓存）；BLAST top hit 仅作最后回退；
2. 级联查宿主概率表（`databases/host_prob/{species,genus,family,order}_
   host_probability.tsv`，共 5.9 万条目，源自 ICTV 分类 × 宿主记录交叉统计），
   种→属→科→目逐级回退 + 列错位容错，取层级最深、置信度最高者；
3. 交叉证据：info.tsv 的 NCBI 宿主元数据（accession 官方宿主记录）经 NCBI
   taxonomy 归类到宿主类别（带缓存）；
4. 决策：两者一致 = Agree；不一致时 NCBI 元数据（一手证据）优先；仅其一
   用其一；全无 = Unknown。

输出（`08_host_analysis/`）：host_prediction.tsv（逐 contig 明细）、
host_summary.tsv、按宿主拆分的 `{类别}.classified.fasta`、
sankey_host.html（病毒科→宿主类别桑基图）、sunburst_host.html（科→属→种
旭日图），并自动嵌入 ⑧报告。CLI：`python main.py host-analysis --sample 样品名`。

## 每样品输出（results/<样品>/）

```
00_prep/          fastp_R1/R2.fastq.gz, conv_R1/R2.fa.gz（fq2fa 产物，
                  供 ① 分类直接复用）, fastp_report.html, input.json
01_host_removal/  kept_R1/R2.fastq.gz, stats.json, host.kreport2
02_virus_screen/  virus_summary.tsv（种/属/科+谱系）, summary.json,
                  viral_R1/R2.fastq.gz（提取的病毒 reads，等价 KrakenTools
                  extract_kraken_reads：默认全部 C 行，可按 taxid 子树提取）
03_assembly/      spades/, contigs.filtered.fasta, virus_contigs.tsv, contig_blast.tsv
04_orf/           orfipy_pep/nt/bed, pyrodigal.faa/ffn/gff, pyrodigal_rv.*
04b_orf_annot/    orf_annotation.tsv（逐 ORF 功能注释）, orf_function_summary.tsv,
                  orf_family_summary.tsv, contig_function_profile.tsv,
                  orf_annotation.gff3（含 product/category/evidence 属性）,
                  genome_diagrams/<contig>.svg（功能着色示意图）, summary.json
05_phylo/<组>/    aln.fasta, tree.nwk(+tree.png), sdt_matrix.csv, sdt_input.fas（可拖入 SDT）
06_primer/        primers.tsv（引物序列/位置/Tm/GC/产物大小/二聚体警示）
07_report/        report.html（汇总所有图表，双击浏览器打开）
08_host_analysis/ host_prediction.tsv, host_summary.tsv, sankey/sunburst_host.html
09_genome_plots/  <contig>.circular.svg / .linear.svg（gbdraw，嵌入报告）
logs/             任务日志（含每阶段资源预估/实际耗时/磁盘核对）
```

**组装输入**（分析页下拉框）：默认「病毒 reads」——用阶段②提取的病毒 reads
组装（宿主/杂菌污染最少、最聚焦），病毒 reads 不足 500 对时自动回退去宿主
reads / 原始 reads 并在日志说明；也可手动固定用去宿主或原始 reads。

## 平台目录说明

```
vp/           核心 pipeline 包（GUI 与 CLI 共用）
app.py        Web GUI 服务（仅监听本机）
webapp/       页面模板与静态资源（app.css / app.js / i18n.js）
main.py       CLI 入口
启动平台.bat   双击启动 GUI
环境自检.bat   环境自检（python main.py selfcheck + 内存检查）
scripts/      辅助脚本（build_virus_db.py / package.py / VirusPlatform.spec）
bin/          单文件外部工具（kunpeng / seqkit / crabz / clustalw2 /
              muscle / aria2c / sracha）
tools/        工具套件（Blast / mafft-win / FastTree / iQtree / trimAl /
              Gblocks / SDTv1.3 / diamond / mmseqs / fastp / SPAdes 安装包）
open-virome/  Open-Virome 前端构建源（app.py _FRONTEND_BUILD 引用，勿改名/搬动）
databases/    kunpeng 库（host_db / virus_db）与 taxonomy；examples/ 内置
              各工具示例数据（✨示例按钮共用：2 条病毒基因组、6 条同属
              近缘集、保守区近缘集、示例树、示例 GenBank 与共线性 .gb）
results/      样品结果
logan/        LOGAN 溯源查询任务（<查询名>/ 含查询 FASTA、导入的结果表
              与 trace_report.html 溯源报告）
tasks/, logs/ GUI 任务状态与运行日志
tool_runs/    工具箱独立运行产物（规范化结构，见下节"输出目录管理"）
tests/        自检、造数与集成测试脚本（_it_annotate / _it_compare /
              _it_platform / _it_msa / _it_synteny / _it_phylo /
              _it_submit / _smoke_phylo；make_examples.py 生成内置示例）
docs/         开发文档（DEVELOPMENT_NOTES.md：踩坑记录与目录规范）
platform.json 工具路径、界面语言与分析默认参数配置（可手动改）
```

工具路径在 `vp/config.py` 自动探测（先 `bin/`、`tools/`，兼容旧根目录
布局），`platform.json` 的 `tools` 段可手动覆盖，改动后立即生效。

## 输出目录管理（tool_runs 规范化）

`tool_runs/` 采用「固定目录 + 动态运行区」双层结构：

```
tool_runs/
├─ _archive/   固定：历史运行归档（_archive/<工具>/<运行>/），长期留存
├─ _tmp/       固定：冒烟测试与一次性实验（随时可整体清空）
├─ _scripts/   固定：运维/补丁脚本
├─ _reports/   固定：跨运行汇总产物（统计表、汇总图）
└─ <工具>_<YYYYMMDD_HHMMSS>/   动态活动区（扁平命名，新运行都在这层）
```

- `_` 前缀 = 固定目录，不出现在 GUI「运行目录」列表，也永不被
  archive/clean 当作运行处理；无时间戳命名的目录一律视为杂项。
- 管理动作（含 `--dry-run` 预演，删除不可逆、务必先预演）：

```bat
python main.py tool-runs status                        :: 活动区/固定区统计
python main.py tool-runs organize [--dry-run]          :: 杂项自动归位
python main.py tool-runs archive --before 20260901     :: 归档旧运行
python main.py tool-runs clean --active-days 30 --archive-days 90 --dry-run
```

建议节奏：每月 `archive --older-days 30`；磁盘紧张时 `clean --active-days
30 --archive-days 90`（活动区留 30 天、归档区留 90 天）。设置自定义输出
根时，固定目录会在新根的 `tool_runs/` 下自动预建。

## 软件与数据库分离（打包 / 分发）

平台 60+GB 里约 47GB 是 `databases/`（kunpeng 库、taxonomy、hmm 等）。
为了便于分发，**软件包默认不包含数据库**：

```bat
python scripts/package.py                       :: 软件包（bin/tools/前端/exe，~1GB）
python scripts/package.py --with-db             :: 软件包 + 预置最小病毒库/宿主库（开箱即用）
python scripts/package.py --db-only             :: 对已有发布目录补置数据库（供重打库）
```

- 软件包产物在 `dist/VirusPlatform/`（含「数据库对接说明.txt」），可整体拷到
  任意电脑运行；`platform.json` 为干净配置，不含本机绝对路径，启动自动探测。
- **数据库对接（二选一）**：
  - 方式一：`python main.py db-migrate --to D:\库目录`
    （把本机 databases/host-db/virus-db 复制到目标盘并自动写 database_root；
     `--mode move` 复制校验后删源；`--dry-run` 只预检；`--check` 只校验目标已有库）
  - 方式二：库已拷到目标目录 → 启动平台 →「设置 → 数据库目录」填路径 → 应用。
    保存时会即时反馈宿主库 / 病毒库 / Taxonomy 是否就绪，缺哪个明示哪个。
- 迁移安全：先 robocopy（多线程+断点续传）→ 校验文件数与总字节一致 → **通过后才
  切换配置**；任一环节失败配置不动、源目录原样，可断点重跑。目标目录必须在平台
  目录之外，且目标盘剩余需 ≥ 源体积 + 10GB。
- 清理开发机：迁移确认正常后，原 `databases/` 可自行删除释放空间。

## 依赖安装（首次）

```bat
python -m pip install -r requirements.txt
```

外部工具已内置在平台目录（kunpeng / crabz / seqkit / Blast / mafft-win /
FastTree / iQtree(v2+v3) / trimAl / Gblocks / SDTv1.3 / diamond / mmseqs），
SPAdes 需已安装
（`SPAdes-Windows-4.3.0-dev-Setup.exe`）且 `spades.bat` 在 PATH；
基因组图首选 gbdraw：`pip install git+https://github.com/satoshikawato/gbdraw.git`
（打包版把 gbdraw.exe 放入平台根目录或 PATH 亦可）；未装时自动回退
纯 Python 引擎 dna_features_viewer（随 requirements.txt 安装）。
CLI/卡片可用 `--plot-engine auto|gbdraw|dfv` 强制指定引擎。

## ORF 功能注释模块（⑥b）

对 ⑥ORF 的病毒蛋白（pyrodigal 基因模型优先，回退 pyrodigal_rv / orfipy）
做蛋白级功能注释，回答"每个 ORF 是什么蛋白、属于哪个科属、基因组什么类型"：

1. 蛋白搜索：DIAMOND（首选）→ MMseqs2 → BLAST+ blastp 自动降级，对照
   **NCBI RefSeq 病毒蛋白库**（release viral 全量，首次运行自动从 NCBI
   官方域名下载 ~107MB 并建库，之后全缓存）；
2. 功能类别：**E-value 加权类别投票**（吸收自 LanderDC/annotation_benchmark
   基准结论：序列同源搜索是病毒蛋白注释最准确的主策略——BLASTp 78.4% >
   结构折叠法 72.2% > 蛋白语言模型 63.6%；有效类别权重须 ≥ 最佳
   "假想蛋白"类别的 50% 才判有效）；
3. 科/属/宿主：命中物种经 NCBI taxonomy 谱系解析到属/科，再查 ICTV 科级
   对照表（265 科，源自 LazypipeX）得**宿主类别**与**基因组类型**
   （ssRNA(+)/dsDNA/…），与 ④宿主预测 交叉印证；
4. 输出：orf_annotation.tsv（逐 ORF：产物/物种/identity/覆盖度/E-value/
   科属/宿主/基因组类型/功能类别/证据来源/hmm_hits/cdd_hits）、类别与科分布
   汇总、逐 contig 功能画像、**orf_annotation.gff3**（pyrodigal GFF 追加
   product/organism/category/evidence 属性）、**genome_diagrams/**（每条
   病毒 contig 的线性示意图：ORF 按功能类别着色的方向箭头 + HMM/CDD
   结构域窄条 + 比例尺，littlegenomes 风格适配版，零依赖 SVG）；类别/科
   分布图与 Top 注释表自动嵌入 ⑩报告。

### 层2 HMM / 结构域（序列法兜底，远缘 ORF 补注释）

pyhmmer（Windows 原生 wheel，内存恒定 ~50MB 流式扫描）扫三库，命中过滤用
viral_fams 口径**双门槛**：域级 i-Evalue ≤ 1e-3 **且** HMM 模型覆盖率 ≥ 0.5：

| 库 | profiles | 功能/分类归属 |
|---|---|---|
| VOGDB r236 `vog.hmm` | 49,116 | annotations 共识描述 + 功能类别码 + lca 谱系尾（科/亚科） |
| RVDB v32.0 FAM | 13,679 | 预计算 `rvdb_annotations.tsv.gz`（VOG 同构：LCA 谱系→科/属 + 关键词共识描述 + 类别） |
| vFam-B 2014 | 5,585 | 注释文件 FAMILIES/GENERA 主科/主属 + 代表产物名 |

序列层已注释的 ORF（evidence=seq）保留层1 结果、HMM 命中记入 hmm_hits 列；
序列层未命中的 ORF 由 HMM 兜底（evidence=hmm，产物/类别/科来自 VOG/RVDB
元数据）。另含 **CDD 结构域层**：mmseqs2 搜索 NCBI Cdd（替代 RPS-BLAST，
同 Cenote-Taker3 思路）——把 mmseqs 格式 CDD 库放 `databases/cdd/cdd_db`
（或 platform.json databases.cdd 指定前缀）即自动启用，命中经
viral_cdds_and_pfams_191028.txt（1,580 条精选病毒域列表）标记病毒相关性，
evidence=cdd。首次扫描前平台自动 hmmpress 预压（VOG 约 25-40 分钟，一次性，
期间该库显示"预压进行中"自动跳过）。

CLI：`python main.py orfa --sample 样品名`。
参考库可手动换：替换 `databases/viral_prot/viral_prot.faa` 后删除
`organism_tax.tsv`、`db_info.json` 重跑即自动重建索引与搜索库。

## 基因组图模块（⑨，gbdraw / dna_features_viewer 双引擎）

对 ③组装 的病毒 contigs 逐条出 **圈图 + 线图（SVG）**。两引擎：

- **gbdraw**（默认首选，圈图更精美）：注释自动取 ⑥ORF 的 pyrodigal
  GFF3；安装 `pip install git+https://github.com/satoshikawato/gbdraw.git`。
- **dna_features_viewer（DFV）**：gbdraw 不可用时自动顶上（纯 Python，
  `pip install dna_features_viewer` 即可）；注释优先取 ⑥b 的
  `orf_annotation.gff3`（带 product/category，按结构蛋白/聚合酶等功能
  类别着色），回退 ⑥ 的 pyrodigal GFF3；支持 GenBank 输入。

引擎选择：卡片参数或 CLI `--plot-engine auto|gbdraw|dfv`（默认 auto）。
也可在卡片参数（或 CLI `--gbdraw-fasta/ann`）改用**自备 FASTA + GFF3，
或直接 GenBank `.gb/.gbk` 文件**。出图上限默认 12 条（按长度取最长）。
SVG 内嵌 ⑩报告，结果文件在 `09_genome_plots/`。

## LOGAN 溯源模块（独立页面，导航栏直达，两种提交方式）

回答"这条病毒序列还出现在哪些公开数据里"：把病毒 contig（自动切成
≤2.5kb 查询片段，Logan-Search 单条上限）提交到
[Logan-Search](https://logan-search.org/dashboard)（IndexThePlanet 计划，
对整个 NCBI SRA 全量组装后建立的 k-mer 索引，覆盖 ~2340 万公开样本；
k=31，返回每个样本的共享 k-mer 比例与 ANI 估计）。

**方式一 · 一键批量（推荐）**：任务卡片填通知邮箱（逗号分隔多邮箱自动
轮换防限额）、选 Groups → 点「开始批量提交」。平台以子进程调用内置的
`vp/logan_submit.py`（Selenium 驱动本机 Edge/Chrome，默认 headless）逐条
提交全部未导入片段，自动轮询结果、下载结果表、**导入并生成溯源报告**。
进度为**实时推送**（SSE）：提交脚本在每个关键节点输出结构化进度事件
（正在提交 / 已提交 session / 等待服务器结果（含已等秒数与 HTTP 状态，
每 30s 刷新）/ 下载中 / 已完成 n 条），页面经 EventSource 订阅即时刷新
进度条与日志，不依赖定时轮询；进度条按"已完成条数"驱动、永不回退。
提交侧等待期从整段盲睡改为每 30s 轻量探测，**结果提前就绪即提前下载**
（总等待窗口不变）。自带断点续跑与补漏：个别片段超时未落定，**再点一次
批量提交即只补漏未导入片段**；邮箱池轮换、每 10 条重启浏览器防 session
失效。依赖 selenium（`pip install selenium`，已列入 requirements.txt），
未安装时该面板灰显、手动方式不受影响。
CLI：`python main.py logan-batch --name 查询名 --email a@qq.com[,b@qq.com]
[--group Fast_No_human] [--show-browser]`。

**方式二 · 手动半自动**：复制片段序列 → 自己浏览器打开 Logan-Search 提交
（平台不发任何外部请求）→ 下载结果表（CSV/TSV）→ 回到页面点对应片段
「导入结果」→ 自动聚合出报告。

两种方式产物一致：物种分布条形图、k-mer×ANI 散点、样本类型分布、
逐片段样本清单（可复制 Run 列表），来源为平台样品时自动附
**④ICTV 宿主预测 × LOGAN 实测物种交叉对照**。
入口：导航栏「LOGAN 溯源」；任务存于 `logan/<查询名>/`
（batch_out/ 为批量下载的原始结果表）。其他 CLI：`logan-create /
logan-import / logan-jobs`。
引用：Chikhi et al. 2025, bioRxiv 10.1101/2024.07.30.605881。

## 进度条 · 预计剩余时间 · 资源预估日志

- 任务卡实时显示**全局进度**（按阶段耗时加权，而非步数均分）、已运行
  时长与**预计剩余时间**。剩余时间来自每阶段耗时自学习模型
  （`logs/stage_perf.json`，随使用越来越准；首次运行使用缺省粗估）。
- kunpeng 分类无原生进度输出，平台以 chunk 中间盘增长近似上报；
  SPAdes 解析 spades.log 按 k-mer 阶段上报；过滤/提取/转换按记录数上报。
- 每阶段开始/结束在日志输出 `📊 资源预估`（线程 / 内存 / 磁盘）与实际
  耗时行；kunpeng 分类结束输出 `📊 磁盘核对`（预估 vs 实测 chunk 占用）。

## crabz 加速（可选）

平台根目录放入 `crabz.exe`（[crabz releases](https://github.com/sstadick/crabz/releases)，
已内置）后，**所有 .gz 读写自动改走 crabz 多线程管道**（kept/viral reads、
子样本、fq2fa 产物等），实测写出较 Python gzip(级别9) **快约 10 倍**
（195MB 模拟 FASTQ：21.5s → 2.1s），压缩率相当；未放 crabz 自动回退
Python gzip，无功能差异。

## 已知事项与设计说明

- **中文路径兼容**：SPAdes 与 BLAST(LMDB) 不支持含中文的路径。平台自动处理：
  SPAdes 经 `%TEMP%\vp_spades`（纯 ASCII）中转并把结果拷回；BLAST 库建在 `%TEMP%\vp_blast`。
  其余工具（kunpeng/mafft/FastTree/IQ-TREE）原生支持中文路径。
- **宿主库构建（大基因组内存问题已修复）**：kunpeng convert 阶段按 60 条/批
  整批载入序列且无字节上限，多条大染色体（如玉米 chr01 176MB）同批会触发
  确定性的 ~20.8GB 巨型分配而失败（与机器内存、线程数无关，1.8GB 基因组
  689 条序列必现）。平台在注入 taxid 时自动把每条序列切成 ≤1MB 片段
  （`>{id}.p00001|kraken:taxid|N`），且相邻片段重叠 34bp（k-1），跨切口
  k-mer 全部保留——全量实测 k-mer 数与整条建库完全一致（286,327,790，
  kunpeng estimate 口径一字不差），2GB 内存即可完成 1.8GB 基因组建库
  （约 1.5 分钟，库文件 ~1.6GB）。建库为替换式：重试会自动清理旧 library，
  不会累积重复数据。
- **分类树来源**：汇总表/图表的科属种层级直接来自 kunpeng kreport（轻量解析），
  不加载全量 NCBI taxonomy，避免内存压力。
- **SDT 分析**：平台用 MAFFT 比对 + 成对 gap 删除口径重算全长 identity 矩阵
  （与 SDT 算法一致），输出 CSV/热图；同时生成 `sdt_input.fas`，
  可在「分析结果」页点「启动 SDT」后拖入 SDT v1.3 GUI 交互查看。
- **断点续跑**：每阶段 `.done` 标记；SPAdes 已有组装结果时重跑直接复用。
- **数据安全**：GUI 仅监听 127.0.0.1，所有文件读写限制在平台目录内，
  下载仅限 NCBI 官方域名白名单。

## 常见问题

**Q: 页面报"浏览失败: TypeError: Failed to fetch"或点击无反应？**
这是浏览器连不上平台服务（服务端经测试一切正常），即**服务进程已退出或无响应**：
1. 黑色控制台窗口被关闭了——关闭它就等于退出平台。重新双击
   `VirusPlatform.exe`（或 `启动平台.bat`），浏览器会自动打开新页面；
2. 双击了多次平台，浏览器停在已关闭实例的旧页面——关掉旧标签页即可；
3. 大任务把内存/CPU 占满导致服务暂时无响应——等任务结束（页面顶部出现
   红色"连接已断开"横幅时说明服务失联，横幅消失即恢复）。
新版页面会在失联时自动显示红色诊断横幅，恢复后自动消失。

**Q: 输入文件支持哪些格式？支持压缩吗？**
- 测序数据：FASTQ / FASTQ.gz（`.fastq` `.fq` `.fastq.gz` `.fq.gz`），双端选
  R1/R2，单端只选 R1；
- 基因组/参考：FASTA / FASTA.gz（`.fa` `.fna` `.fasta` `.fas` 及 `.gz`），
  建库（宿主/病毒）与分析输入均支持 gzip 压缩；
- `.zip` / `.rar` / `.tar` 等归档请先解压出里面的数据文件再选择；
- 文件浏览对话框只显示上述数据文件类型，浏览范围限平台目录内。

**Q: 分析时提示"宿主库不可用"？**
先到「数据库构建」页构建宿主库（需要 TaxID），或取消勾选"①宿主去除"阶段。

**Q: SPAdes 报错 67 / non-ASCII？**
旧版本问题，现已自动中转。若仍出现，检查 `%TEMP%` 路径是否含中文。

**Q: 想用 IQ-TREE 更严谨建树？**
分析页把「建树工具」切到 IQ-TREE（自动优先用 v3，缺失回退 v2；
模型自动选择 `-m MFP` + UFBoot 1000 + SH-aLRT 1000 双支持值，
最优模型与对数似然写入 05_phylo summary.json）。

**Q: 进化树想加入更多近缘参考（如整属全基因组）？**
工具箱「⑤ NCBI 参考序列下载」：输入 Entrez 检索式（如
`Tobamovirus[ORGN] AND complete genome[TITL]`）→ 搜索预览 → 下载为
参考集合（会话式批量下载、accession 去重、断点续传、含宿主/分离物
等元数据 TSV）；再到分析管道「⑦进化分析」的"NCBI 参考集合"填集合名
即可把这些参考追加进每组比对。CLI 等价：
`python main.py ncbi-dl "<检索式>" -n <集合名>`（`ncbi-list` 查看已下载），
分析时 `--ncbi-refs <集合名>`。比对后自动经 trimAl（automated1）清剪，
过度修剪时自动回退原比对。

**Q: 想在网页上直接看比对差异（不用 SDT）？**
「结果」页新增 **MSA 查看（SNP-only 变异热图）**：选样品与 05_phylo 分组，
只显示比对中的变异位点——ACGT 彩色字符热图（A 绿 / C 蓝 / G 橙 / T 红）、
共识行、每列变异度柱，悬浮显示位点计数，大比对自动分页；优先展示建树
所用的清剪后比对（aln.trim）。参考 PhyloSuite MSA Viewer 的 SNP-only
设计用平台自有前端实现。

**Q: 引物对宿主特异性检查很慢？**
首次需对 1.8GB 宿主基因组建 BLAST 库（10-30 分钟），此后复用。

**Q: 想直接比较同属病毒的基因组结构（共线性/基因排列）？**
比较基因组组按分析目的组织为五步：**参考序列获取 → 序列比对（MAFFT+trimAl，
可编辑查看器）→ 进化树构建（科/属级，全基因组/CDS/PEP）→ SDT 同一性（属级）
→ 同属共线性（LoVis4u）**。

「参考序列获取」卡整科整属拿序列：**ICTV 界→门→纲→目→科→属→种递进级联
下拉**（选得越深集合越聚焦；计数=accession 条数）（databases/ictv_db
谱系，本地序列优先、缺的 NCBI 补齐；预览离线可用）或 accession 列表 /
Entrez 检索式 / 本机 .gb 导入。每次下载同时产出 GenBank 集合（保留 CDS
注释，供建树与共线性）与 FASTA 参考集（样品流程建树追加参考，集合名即
「NCBI 参考集合」名）。下载可设<b>每属/每种上限</b>（防大属淹没整科集合）；
集合可一键<b>🧬 提取 CDS/PEP</b>（genome / CDS / PEP 分类分目录 + 按基因
拆分，直接送序列比对或基因建树）。

「同属共线性比较」对集合一键出 **LoVis4u 官方出版级 PDF**（平台已修补其
0.2.0 的 Windows 兼容 bug 并指定内置 MMseqs2；按蛋白组相似度排序、可变
基因高亮），另附基因家族清单与成对共享基因比例表；LoVis4u 不可用时自动
回退平台内置绘图。多聚蛋白基因组（如马铃薯 Y 病毒属）自动改用
mat_peptide 成熟肽段做基因。结果存
`databases/gb_collections/<集合名>/compare/`。

「进化树构建（科/属级）」：选集合 → MAFFT + trimAl → NJ / FastTree /
IQ-TREE（含双支持值），页内树查看与 SNP-only MSA 查看一体；也可直接
输入 FASTA 建树。历史项目/测试样品已收进结果中心「归档项目」，各模块
下拉只显示当前样品。
CLI 全家桶：`gb-dl`（检索式/accession 下载）、`gb-import`（本机 .gb
导入）、`gb-list` / `gb-check`（集合巡检：记录/CDS 数/警告）、
`python main.py compare -n tobamo [--lovis4u] [--style lovis|category]
[--min-ident 0.3 --min-cov 0.5 --order ...]`。
注意：比较需要记录带 CDS 注释；FASTA 无注释记录会被跳过并在日志提示。

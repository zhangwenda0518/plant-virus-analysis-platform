# BioAider 反编译参考源码

来源：BioAider V1.727 Linux 版（PyInstaller ELF，经 pyinstxtractor-ng 解包 +
uncompyle6 反编译；PYZ 为 AES-CTR 加密，密钥嵌于 pyimod00_crypto_key 模块）。

用途：仅供本平台实现借鉴功能时**阅读算法逻辑**（参考 `docs/BIOAIDER_REFERENCE.md`）。

⚠️ **安全声明**：本目录是第三方反编译代码的**只读参考件**——平台任何模块
都不导入、不执行这些文件；其中可能包含 BioAider 原作者自身的编程习惯
（如 `shell=True` 启动外部程序），**切勿**把此类写法带回本平台（本平台
命令执行一律参数列表 + shell=False，路径一律 check_path/safe_open）。

这些代码版权归 BioAider 作者所有，本平台**不直接复制其代码**，仅参考算法思路
自行实现；如发布涉及借鉴的功能，按其手册要求引用 BioAider 论文/仓库。

核心文件导读：
- mutation_tools_mua_fuc.py     突变分析核心（single_mutation_nt/codon、频率图、棒棒糖数据）
- mutation_tools_mutation_analysis.py  突变分析主流程（参考=首序列、位点计数、gap/简并过滤）
- mutation_tools_site_counter.py / site_screen.py  位点计数与筛选
- draw_tools_lollipop_draw_run.py  （空壳，实际绘图在 mua_fuc.nuc_aa_plot/codon_plot/lollipop_data）
- phylogeny_tools_partition_make.py  分区模型文件生成（IQ-TREE #nexus charpartition / MrBayes lset+prset）
- custom_libraries_public_functions.py  格式转换 fas↔nex/phy/pml、fasta 读取、突变索引等公共函数
- clustering_tools_get_kmer.py  k-mer 计数矩阵
- download_tools_NCBIseqbatch_get.py  NCBI 批量取序列
- custom_libraries_codon_trans.py  密码子表（3027 行，全翻译表数据）

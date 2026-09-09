# -*- coding: utf-8 -*-
"""vp.public_meta — 公共数据检索与元数据模块。

源自 MMPV-RNA public_metadata_pipeline 的功能子集：

  - search_engine.py   SRA + GSA 双引擎物种检索（合并去重 → 检索结果 CSV）
  - info_engine.py     Run 元数据深度提取 → Core14 / Full 统一元数据表
  - landscape_plot.py  检索/元数据 SCI 级可视化
  - host_genome.py     宿主参考基因组下载（datasets CLI / E-utilities 回退）

下载环节复用平台既有「数据下载」模块（ENA/NGDC 直下 + aria2c + fasterq-dump），
不引入原管线的 gsa_sra.down.py / sra2fastx.py。建宿主库请用平台
build-host-db（kunpeng），不引入原管线的 kraken2/bowtie2 建库脚本。

CLI：python main.py meta-… / host-genome；Web：导航「公共数据检索」页。
"""

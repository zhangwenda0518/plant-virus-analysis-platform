# -*- coding: utf-8 -*-
"""vp.ncbi_submit — NCBI 序列提交准备模块。

源自 MMPV-RNA virome_submission_pipeline / submission_gui 的功能子集：
unified_metadata.csv 一张表驱动 GenBank + BioSample 提交准备。

  - unified_metadata.py  字段契约 / 初表生成 / source.src / miuvig / assembly
  - store.py             提交项目存储：编辑 / 批量填充 / 校验 / BioSample TSV
                         / template.sbt 生成 / 产物导出
  - report_html.py       交互式可编辑提交报告（report.html）

CLI：python main.py submit-… ；Web：导航「提交准备」页。
"""

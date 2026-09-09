# -*- coding: utf-8 -*-
"""Web 层（Flask 路由与任务编排）。

拆分自 app.py 单体，目标结构见 docs/DEVELOPMENT_NOTES.md。

约定（重要）：
  - 所有 blueprint 一律 `from vp.web.state import ...` 取共享状态；
    禁止 `import app`——那会形成 app → blueprint → app 的循环导入。
  - 业务逻辑仍写在 vp/ 下（pipeline / kunpeng / phylo 等），本包只做
    「HTTP 边界 + 任务工厂」的胶水。
  - 每拆一批后用 `python tests/_route_inventory.py` 验证 179 条路由零差异。
"""

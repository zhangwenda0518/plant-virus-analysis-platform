# 平台优化评估报告

> 评估日期：2026-09-07 · 评估范围：全平台代码 + 运行态数据 + 开发笔记
> 结论先行：**功能覆盖已经很完整，短板集中在工程健康度与科研产出闭环**。
> 优先做「并发闸门 / 内存泄漏 / 轮询改 SSE / 磁盘水位」四项，均为小改动大收益。
>
> **已二次核验（2026-09-07 晚）**：文中每条论断都回到源码复核过。核验结果与修正记录见文末
> 「附录 · 核验记录」，其中 2 处初版表述有误已订正，3 处做了精确化。

---

## 一、现状体检

| 维度 | 实测数据 | 判断 |
|---|---|---|
| 代码规模 | `vp/` 34 模块 1.8 万行 + `app.py` 5871 行 + `main.py` 926 行 | 后端分层尚可，Web 层过胖 |
| 路由数量 | `app.py` 内 159 个 `@app.route` | 单文件承载全部路由，定位困难 |
| 前端 | 20 个模板 7504 行 + `app.js` 3316 行 + `i18n.js` 731 行 | 模板偏大（`tools.html` 2060 行） |
| 磁盘占用 | 总计 **64 GB**，其中 `databases/` **47 GB** | 缺少清理与水位预警 |
| 运行产物 | 4 个样品、103 个 `tool_runs/`、100 个 `tasks/*.json`、84 个日志 | 只 trim 了 tasks（100 上限），其余无策略 |
| 测试 | `tests/` 11 个 `_it_*.py` 手工脚本，无 pytest / 无 CI | 159 个路由靠手点，回归靠人肉 |
| 语法检查 | 全量 `py_compile` 通过 | 基线健康 |
| 安全 | `check_path` 拒绝 `..`、写限定平台内、仅监听 127.0.0.1 | 做得扎实，无需大改 |

### 已经做得很好的部分（不要动）

- **路径安全模型**：`vp/utils.py::check_path` 拒绝 `..` 段、写模式限定平台根内，这是正确的设计。**关键的事实是：`app.py` 里裸 `open()` 数量为 0**——所有接收用户输入的 Web 层都走了安全校验，安全边界是守住的。
  - 需要澄清的一点：`utils.py` 的注释写「`open()` 仅出现在以下两个函数内」，但**这个约定在 `vp/` 内部实际被破了 33 处**（assembly 5、hmm_annot 6、public_data 5、utils 5、config 2、contig_annot 2、ictv_db 2、phylo 2、viz 1、genome_diag 1、lovis4u_run 1、sdt_exact 1、universal_ref 1）。这些多为模块内部自己 `os.path.join` 出来的路径，不接受用户输入，风险有限，但失去了 `safe_open` 的 gzip 自动处理与编码统一（这是 CRLF 坑的防线）。建议新代码仍走 `safe_open`，存量可不动。
- **外部工具降级链**：DIAMOND → MMseqs2 → blastp、gbdraw → DFV、aria2c → 内置 HTTP → sracha，缺件不阻断流程。
- **断点续跑**：每阶段 `.done` 标记 + `logs/stage_perf.json` 耗时自学习 + 资源预估日志，长任务的工程化到位。
- **SSE 基建**：`/api/task/<tid>/stream` 已经有差量推送（内容变了才 `yield`），只是目前只有 LOGAN 面板在用——这是后面优化轮询的现成轮子。
- **开发笔记**：`docs/DEVELOPMENT_NOTES.md` 40+ 章节记录踩坑与决策，可维护性远好于同类项目。

---

## 二、P0 · 立即做（小改动，大收益）

> **实施状态（2026-09-07 晚）**：第 1~3 项**已实施并验证通过**（见每节末尾的
> 「已实施」说明），第 4 项（磁盘水位）未做。

### 1. 全局并发闸门 —— 目前没有任何并发上限 ✅ 已实施

**问题**：`TaskManager`（`app.py:167`）用一个 `threading.Lock` 保护字典，但**没有 Semaphore 限制同时运行的任务数**。样品队列 `SampleQueue` 只管样品流程顺序执行；工具箱（`/api/tool/run`）、数据下载、LOGAN 批量、ICTV 下载这些入口都是点一次起一个线程，可以无限叠加。

`platform.json` 里 `threads: 19`，意味着三个工具任务同时跑就能申请 57 线程 + 三份 SPAdes 内存。README 常见问题里那条「大任务把内存/CPU 占满导致服务无响应」正是这个症状——FAQ 把它当成现象解释，其实是可以根治的。

**改法**（约 30 行）：
```python
# TaskManager.__init__
self.slots = threading.Semaphore(2)        # 重任务并发上限
self.light_slots = threading.Semaphore(6)  # 轻量任务（下载/查询）

# start() 增加 weight 参数：'heavy' | 'light'
def start(self, name, fn, log_file=None, weight='heavy'):
    ...
    def _run():
        slot = self.slots if weight == 'heavy' else self.light_slots
        with slot:
            ...原逻辑...
```
再在 `api_tool_run` 里按工具类型标注 weight（SPAdes 组装 / kunpeng 分类 / ORF 注释 / 建树 = heavy；下载 / 检索 / 转换 = light），并在任务卡上显示「排队中 · 前面还有 N 个」。

**配套**：把 `threads` 默认从 19 降到「物理核数 - 2」，并在设置页显示探测到的可用内存，让默认值不至于一上来就打满机器。

**已实施**（`app.py`）：
- 新增 `_slot_limits()`：从 `platform.json` 的 `defaults.max_heavy_tasks` /
  `max_light_tasks` 读上限（缺省 2 / 4，钳制在 1-8 / 1-16），未配置时用缺省值。
- `TaskManager.__init__` 建 `heavy_slots` / `light_slots` 两个 `Semaphore`。
- `start(name, fn, log_file=None, weight='heavy')`——**默认 heavy**，宁可多排队也不放过不认识的任务。
- `_run()` 先用 0.5s 轮询 `slot.acquire(timeout=0.5)` 排队（期间响应取消，不会死锁），拿到名额才置 `running` 并把 `started` 改为真正开跑的时刻。
- 排队状态**对外呈现为 `running`**（`snapshot()` 里映射），这样前端取消按钮、进度条、日志滚动全部零改动可用；排队提示走 `msg`（"排队中 · 前面还有 N 个任务"），前端 `app.js:1112` 会渲染。
- `cancel()` 支持 `queued` 状态；`_recover_interrupted()` 把重启遗留的 `queued` 一并标记失败。
- 工具分级：`LIGHT_TOOLS = {'convert', 'genoplot', 'primer'}`，其余工具 + 样品分析 + 建库 + 建树 + SDT 全为 heavy；ICTV/NCBI/GenBank 下载、CDS 提取、检索绘图、AI 补全、LOGAN 批量、Taxonomy 下载 标为 light。

---

### 2. 任务列表内存泄漏 —— `self.tasks` 只增不减 ✅ 已实施

**问题**：`TaskManager.start()` 每次 `self.tasks[tid] = rec` 并 `self.order.insert(0, tid)`，但**没有任何地方从内存里移除任务**。`_recover_interrupted()` 只 trim 磁盘上的 `tasks/*.json`（保留 100 个），内存字典则是服务跑多久就涨多久。

每个任务持有 `deque(maxlen=500)` 的完整日志 + `result_preview` + `threading.Event`。跑一整天几十个任务后，这是实打实的常驻内存增长，而且 `/api/tasks` 每 2.5 秒会把**全部**任务序列化一遍发出去。

**改法**：
- 内存里只保留最近 200 个任务的完整记录，更早的降级为「仅状态摘要」（清掉 `log` / `result` / `thread` 引用）。
- 或者更简单：任务完成后 5 分钟，把 `rec['log']` 换成空 deque、`rec['thread'] = None`，磁盘日志已持久化，不影响回溯。

**已实施**：新增 `TaskManager._trim_memory()`，在 `_run()` 的 `finally` 里调用（每个任务结束时触发一次稳态回收）：
- 超过 `_KEEP_FULL`(60) 的终态任务 → 清空 `log`、`result`、`result_preview`、`thread`、`procs`（磁盘 `tasks/*.json` 仍在，可回溯）。
- 超过 `_KEEP_TOTAL`(300) 的 → 直接从 `self.tasks` / `self.order` 移除。
- running / queued 的任务一律跳过，不影响正在跑的。

---

### 3. 前端轮询：2.5 秒全量拉 `/api/tasks` ✅ 部分实施

**问题**：`webapp/static/app.js:1026` 是 `pollTimer = setInterval(refreshTasks, 2500)`，每次拉 `tm.list_all()`——所有任务、每个 30 行日志、每个结果预览。任务累积后这个 payload 会持续变大，而且是**无意义的全量**：绝大多数任务处于终态，根本不会变。

现在已经有更好的轮子没用上：`/api/task/<tid>/stream` 的 SSE 实现带差量判断（`data != last` 才推），LOGAN 面板用的就是它。

**改法**（分两步，可拆开做）：
1. **立竿见影**：轮询只拉 running 的任务。`/api/tasks` 增加 `?active=1` 参数，`list_all()` 里过滤掉非 running 的；终态任务由页面在状态变化时增量拉取一次即可。
2. **彻底方案**：把 `refreshTasks` 换成 SSE 订阅「任务总线」事件（新建 / 状态变更 / 完成），浏览器端只在事件到达时更新对应卡片。保留 2.5 秒轮询作为 SSE 断连时的降级——现有的 `setConnBanner` 失联横幅机制要保留。

**已实施（第 1 步 + 自适应降频，未做 SSE 改造）**：
- `list_all(limit=120)`：运行中/排队中的任务给 30 行日志，**终态任务只给最后 5 行**（日志已落盘）；结果预览一律保留（完成后的结果面板依赖它）。列表长度封顶 120。
- `webapp/static/app.js`：`setInterval(refreshTasks, 2500)` 改为 `setTimeout` 递归自适应——有任务在跑 2.5s，**全部空闲降到 8s**；连接异常时（`connDown`）保持 2.5s 快速重试，避免失联横幅要等满 3 个空闲周期才出现。
- 递归调度会丢掉 `setInterval` 的「异常免疫」，所以 `schedulePoll` 里用 try/catch 兜住 `refreshTasks`，否则一次异常就永久断链。

**未做**：SSE 改造（涉及前端任务卡渲染逻辑，改动面大，建议单独一轮做）。

---

### 4. 磁盘水位与清理策略

**问题**：平台 64 GB，`databases/` 独占 47 GB（rvdb_db 11G、hmm 11G、host_db 5.9G、refvirus_db 5.5G、virus_ref 4.5G）。`tool_runs/` 103 个目录、`logs/` 84 个文件、`results/` 2.9 GB 全部无清理策略。`selfcheck` 里虽然检查了磁盘剩余并提示「建议保留 >10GB」，但**只在手动跑 CLI 时才看得到**。

**改法**：
- 服务启动时 + 每隔 30 分钟检查平台所在盘剩余空间，低于 20 GB 时首页顶部显示黄色水位横幅（复用现有的横幅组件）。
- 「结果中心 / 设置」加一个**存储面板**：按目录显示占用 TOP 列表（`databases/` 各库、`results/` 各样品、`tool_runs/`），支持勾选清理 `tool_runs` 与失败任务的中间产物。**`databases/` 默认只读展示，删除需二次确认**——重建一次宿主库代价很高。
- `tool_runs/` 默认保留最近 30 个，超出部分在下次启动时提示（不自动删，避免误伤正在查看的结果）。

---

## 三、P1 · 重点排期（改动大，但决定平台上限）

### 5. `app.py` 按业务组拆 Blueprint

`app.py` 5871 行、159 个路由。你已经在 `NAV_GROUPS`（`app.py:94`）定义好了 7 个业务组：**数据资源 / 样本处理 / 病毒识别和分类分析 / 病毒注释分析 / 比较基因组分析 / 结果中心 / 溯源与提交**——这就是现成的拆分边界。

建议结构：
```
app.py                    # 仅 Flask 装配 + 错误处理器 + 启动（目标 <300 行）
webapp/blueprints/
  resource.py    meta / virome / build
  sample.py      download / samples / hostremoval
  virus.py       鉴定 / 组装 / 再鉴定
  annotate.py    orf / annotation / genome / primer
  compare.py     序列获取 / 比对 / 建树 / SDT / 共线性
  result.py      results / msa / tree / sdt
  trace.py       logan / submit
  api_tasks.py   任务引擎相关 API（TaskManager 单独成模块）
```
迁移顺序建议从 `trace.py`（LOGAN + submit，耦合最少）开始，抽出后立刻跑一遍烟测；`app.py` 里的业务逻辑（比如 `_tool_job_*` 系列，20+ 个任务函数）应该下沉到 `vp/jobs/`——开发笔记第二节自己写了「勿在 app.py 写业务逻辑」，现在这条被破了。

**注意**：`scripts/package.py` 与 `VirusPlatform.spec` 需要同步加 hiddenimports，否则打包后 Blueprint 导入失败（笔记第 17 条踩过同类坑）。

### 6. 多样品比较分析 —— 目前完全缺失

这是**科研价值最高的一块空白**。需要精确说明：单样品内部的丰度分析是有的（`vp/viz.py::fig_abundance_bar` 丰度柱状图、科/属聚合图），但**跨样品的横向比较完全没有**——全代码库搜不到 `cross_sample` / 多样品 / 样品比较 的入口，后端也只有 `api_samples`（列表）、`api_pipeline`（单样品）这类单样品 API，**没有任何跨样品聚合接口**，`results.html` 里也没有汇总/对比/多选相关的功能。

现在平台是「一个样品一条流水线跑到底，出一份报告」，但博士论文需要的是**多样品横向结论**：

- 多样品病毒谱汇总表（样品 × 病毒种，读数/丰度矩阵）
- 样品间病毒组成 Venn / Upset 图、β 多样性（Bray-Curtis）聚类热图
- 按产区 / 年份 / 组织部位分组的差异分析（元数据已经在 `meta_search` 的 Core14 表里了）
- 同属病毒跨样品的 contig 集合比较

建议新建 `vp/cross_sample.py` + 「结果中心 → 多样品比较」页：勾选若干已分析样品 → 读各样品 `02_virus_screen/virus_summary.tsv` 与 `08_host_analysis/host_prediction.tsv` → 出矩阵与图。**数据源都是现成的，主要工作量在聚合与可视化**，`vp/viz.py`（1291 行）里的绘图函数可直接复用。

### 7. 运行清单 manifest —— 补齐可复现性

**实测**：`results/NX-6/00_prep/input.json` 只有三个字段：
```json
{"sample": "NX-6", "r1": "fastq/NX-5_S2_L001_R1_001.fastq.gz", "r2": "fastq/NX-5_S2_L001_R2_001.fastq.gz"}
```
没有参数、没有软件版本、没有数据库版本。全库搜索 `*manifest*` / `*params*` / `run.json` 均无命中——**确实不存在运行清单**。

（补充澄清：平台本身是有版本号的，`vp/__init__.py:4` 定义了 `__version__ = "1.0.0"`，只是**这个值从未被写入任何产物**。`vp/config.py` 里确实没有 VERSION 常量。）

论文方法章节和审稿人追问「参考库版本 / 软件版本 / 参数」时，目前只能翻 `logs/` 人肉拼。

建议每次 `run_analysis` 在样品目录写一份 `run_manifest.json`：
```json
{
  "sample": "NX-6",
  "started": "2026-09-07T19:00:00",
  "platform_version": "1.4.0",
  "stages": {"host": {"tool": "kunpeng", "version": "0.7.12", "params": {...}}},
  "databases": {"virus_db": {"path": "...", "built": "...", "n_seqs": 8465},
                "virus_ref": {"DATA_VERSION": "..."}},
  "params": {"confidence": 0.0, "assembly_mode": "metaviral", ...}
}
```
配套做一个「一键导出方法章节」：把 manifest 渲染成 Markdown 段落（软件 + 版本 + 参数 + 数据库版本），直接贴进论文。对你当前写论文这件事，投入产出比很高。

---

## 四、P2 · 顺手清理

### 8. 冗余静态资源约 100 MB
`webapp/static/virome/static/js/` 里堆了三代构建产物：`main.0f86404a.js`、`main.2fb51f4d.js`、`main.37e52565.js`，每个 6.2 MB，外加三个 16 MB 的 `.map`（合计约 68 MB，算上 chunk 约 100 MB）。`index.html` 只会引用其中一个。清理到只留当前版本即可（先确认 `asset-manifest.json` 指向哪个）。

### 9. 测试脚本 pytest 化 + 路由冒烟

现状要说准确：`tests/` 下有一个 pytest 风格的 `test_logan_trace.py`（含 7 个 `def test_`），说明**写法是有的**；但另外 11 个 `_it_*.py` 是独立脚本靠手工跑，且**全项目没有 `pytest.ini` / `conftest.py` / `pyproject.toml` / 任何 CI 配置**，所以测试实际上不会被自动执行。

建议：
```bash
pip install pytest
```
- 把 `_it_platform.py` 里的 Flask test client 冒烟扩成**全路由自动巡检**：遍历 `app.url_map` 对所有 GET 路由发请求，断言非 500。159 个路由一次跑完，比手点可靠得多。
- `_it_annotate / _it_compare / _it_msa / _it_phylo / _it_submit / _it_synteny` 加 `@pytest.mark.slow`，CI 里可选跳过。
- 至少加一个 `test_paths.py` 覆盖 `check_path` 的路径穿越用例（`../`、绝对窗外路径、`..%2f`），这是安全底线，必须有自动化兜底。

### 10. 全局异常兜底
错误处理器（`app.py:426-438`）覆盖了 `ValueError` / `FileNotFoundError` / `PermissionError` / `HTTPException`，但**没有 generic `Exception` handler**。未捕获异常会以 HTML 500 页返回给 fetch 调用方，前端 `res.json()` 直接抛错，用户看到的是「无反应」而不是错误原因。补一个：
```python
@app.errorhandler(Exception)
def _unhandled(e):
    app.logger.exception(e)
    return jsonify({'error': f'{type(e).__name__}: {e}'}), 500
```
另外 `app.py` 里共 125 处 `except`，其中 **15 处是 `except Exception: pass`**（`TaskManager._persist`、`_log`、`_recover_interrupted` 等）。这些静默吞异常的地方建议至少改成 `app.logger.debug(...)`，否则线上排查时无从下手。

### 11. `_pick_port` 的 TOCTOU 竞态
`app.py::_pick_port()` 逐个 `bind` 探测后**关闭 socket**，再由 `app.run()` 重新 bind，中间存在被其他进程抢占的窗口（虽然概率低，但一旦发生就是启动失败）。改法：探测成功后把 socket 保持打开传给 `app.run(...,  )`——或者更省事，直接 `s.bind()` 拿到端口后不关，用 `SO_REUSEADDR` 交给 Flask 时传入已绑定的 socket。

### 12. `_recover_interrupted` 的双次排序
`app.py:180-184` 先 `files.sort(reverse=True)` 再按 mtime 重排，第一次排序是多余的。小问题，顺手删。

---

## 五、P3 · 暂缓（收益不足以支撑当前成本）

- **i18n 重构**：`i18n.js` 731 行 + 模板里中英双写，加字段要改两处，长期是负担。但当前功能迭代更重要，等字段稳定后再统一收口成单一键值源 + 缺失键检测。
- **`dist/` 产物外置**：4.3 GB 打包产物混在源目录里（还含运行时 `results/logs/tasks`）。清理价值有，但优先级低于上面各项；至少把 `dist/VirusPlatform/results|logs|tasks` 加进 `.gitignore` 概念（或清理脚本）。
- **Flask 内置服务器换 waitress**：`app.run(threaded=True)` 在本机单用户场景够用。如果后面真要多人/多任务并发，再换不迟。

---

## 六、建议执行顺序

| 批次 | 内容 | 预期耗时 | 风险 | 状态 |
|---|---|---|---|---|
| 第 1 批 | 并发闸门 + 任务内存回收 + 轮询瘦身 | 半天 | 低，改动局部 | **已完成** |
| 第 2 批 | 磁盘水位横幅 + 存储面板 | 半天 | 低 | 待做 |
| 第 3 批 | 全局异常兜底 + 冗余资源清理 + 路由冒烟测试 | 半天 | 低 | 待做 |
| 第 4 批 | 多样品比较分析（新模块） | 2-3 天 | 中，需设计交互 |
| 第 5 批 | run manifest + 方法章节导出 | 1 天 | 低 |
| 第 6 批 | Blueprint 拆分（trace → compare → annotate → …） | 3-5 天 | 中，每步需回归 |

第 1~3 批做完，平台的「跑久了变卡、跑满了失联」这类体验问题基本消除；第 4~5 批直接服务于论文产出；第 6 批是为长期维护，可以等阶段性交付后再做。

---

## 七、值得考虑的两个产品向问题

1. **样品 → 项目（Project）的一等公民化**：现在 `project` 只是样品列表的筛选标签。建议把「项目」提升为实体：一个项目 = 一批样品 + 共享元数据（产区/年份/组织）+ 一个汇总视图。这样第 6 条的多样品比较才有天然的作用域，也贴合「宁夏枸杞不同产区病毒组」这类真实研究组织方式。

2. **分析参数模板**：17 项默认参数现在是全局的。建议支持「参数方案」保存/另存为（如「快速验证」「正式分析」），建样品时直接套用。批量跑几十个样品时，避免每次都逐项确认。

---

## 附录 · 核验记录（2026-09-07 晚二次复核）

逐条回到源码验证，结果如下。**❌ = 初版表述有误已订正，⚠️ = 表述不够精确已细化，✅ = 核实无误。**

### 已订正（2 处）

| # | 初版说法 | 实际情况 | 已改为 |
|---|---|---|---|
| ❌1 | 「`vp/config.py` 里连 `VERSION` 常量都没有」（暗示平台无版本号） | `vp/__init__.py:4` 有 `__version__ = "1.0.0"`，只是从未写入产物 | 第七节：明确「平台有版本号但未落盘」，config.py 无 VERSION 仍属实 |
| ❌2 | 「`open()` 只允许出现在 safe_open/check_path 两个函数内，这个约定非常正确」 | 该约定写在 `utils.py` 注释里，但 `vp/` 内部实际有 **33 处裸 `open()`**（13 个模块） | 第一节：改为「app.py 裸 open = 0（安全边界守住），但 vp/ 内部 33 处偏离约定」，并列出分布 |

### 已精确化（3 处）

| # | 初版说法 | 精确说法 |
|---|---|---|
| ⚠️1 | 「多样品比较完全缺失」 | 单样品丰度图是有的（`viz.py::fig_abundance_bar`）；缺的是**跨样品**聚合——无跨样品 API、`results.html` 无汇总/对比入口 |
| ⚠️2 | 「测试无 pytest」 | `tests/test_logan_trace.py` 有 7 个 `def test_`（写法是 pytest 风）；真正缺的是 `pytest.ini`/`conftest.py`/CI，**测试不会被自动执行** |
| ⚠️3 | 「125 处 except 中不少是 pass」 | 实为 **15 处**（占 12%），`_persist`/`_log`/`_recover_interrupted` 等 |

### 核实无误（12 处）

| 论断 | 验证方式与结果 |
|---|---|
| 无全局并发上限 | `grep Semaphore` 无命中；`TaskManager` 只有 `threading.Lock`。4 处 `status=='running'` 检查都是**同一样品互斥**，非全局限流 ✅ |
| `self.tasks` 只增不减 | `TaskManager` 全部方法 = `__init__/_recover_interrupted/start/_persist/snapshot/list_all/cancel`，**无 cleanup**；磁盘 trim 100 个但内存不清理 ✅ |
| 2.5 秒全量轮询 | `app.js:1069 refreshTasks` 直接 `fetch('/api/tasks')` 无参数；无 active 过滤 ✅ |
| SSE 只有 LOGAN 用 | `EventSource` 仅出现在 `logan.html:302`；app.py:5219 注释自述「LOGAN 批量面板用」✅ |
| 无 run manifest | `input.json` 实测仅 `{sample, r1, r2}`；全库 `*manifest*`/`*params*`/`run.json` 零命中 ✅ |
| 159 个路由 | `grep -c '^@app.route'` = 159，去重后仍 159（无重复装饰器）✅ |
| 错误处理器缺兜底 | 仅 `ValueError`/`FileNotFoundError`/`PermissionError`/`HTTPException` 四个，无 generic `Exception` ✅ |
| virome 冗余 100MB | `index.html` 只引用 `main.37e52565.js`，另两个 6.2MB 的 main + 三个 16MB `.map` 均为残留 ✅ |
| 多样品聚合缺失 | 后端仅 `api_samples`(列表)/`api_pipeline`(单样品)，`results.html` 无汇总/对比/多选 ✅ |
| 安全边界 | `app.py` 裸 `open()` = **0**（用户输入入口全部走校验）✅ |
| 磁盘 64GB / databases 47GB | `du` 实测 ✅ |
| 语法基线 | 全量 `py_compile` 通过 ✅ |

### 核验方法备注

- 并发：`grep -n "Semaphore\|maxsize\|Queue(" app.py vp/*.py` → 仅 `sample_queue = SampleQueue()`
- 队列：`SampleQueue._ensure_worker` 确认单 worker（一个 `self.thread`）顺序执行
- 裸 open：`grep -c "[^_a-zA-Z.]open("` 逐模块统计
- 引用关系：`grep -o "main\.[a-f0-9]*\.js" index.html | sort -u`

---

## 附录二 · 第 1 批实施与验证结果（2026-09-07 晚）

改动文件：`app.py`（TaskManager + 任务启动点）、`webapp/static/app.js`（轮询）、
新增 `tests/_it_concurrency.py`（验证脚本）。

### 验证结果

`tests/_it_concurrency.py` 全绿（**ALL PASS**）：

| 验证项 | 结果 |
|---|---|
| heavy 并发峰值 ≤ 2（起了 5 个 heavy 任务） | 峰值 2 ✅ |
| 排队任务对外 `status=running`（前端取消按钮可用） | ✅ |
| 排队 `msg` 含"排队中 · 等待空闲名额" | ✅ |
| 排队中 `cancel()` 返回 True 且状态 → `cancelled` | ✅ |
| 名额正确归还、无死锁（剩余任务全部完成） | 完成 4/4 ✅ |
| light 并发峰值 ≤ 4，且与 heavy 独立 | 峰值 4 ✅ |
| 运行中任务日志 30 行 / 终态任务日志 5 行 | ✅ |
| 路由冒烟：18 个页面 + 6 个核心 API | 全部 200 ✅ |
| `app.js` 语法（`node --check`） | ✅ |

### 需要知道的行为变化

1. **任务启动后不再立即执行**，可能先排队。任务卡会显示"排队中 · 前面还有 N 个任务"。
   这是预期行为——目的是避免多个吃满线程的任务同时跑把机器打挂。
2. **"已运行时长"从真正开跑算起**，不含排队时间。
3. **已完成任务的展开日志默认只回传最后 5 行**（原 30 行）。要看完整日志需打开
   任务对应的日志文件，或让前端在展开时单独请求 `/api/task/<tid>`（未做）。
4. **空闲时轮询从 2.5s 降到 8s**。有任务在跑时仍是 2.5s，观感无变化。
5. 服务重启后，遗留的 `queued` 任务与 `running` 一样被标记为失败（不会幽灵排队）。

### 可调参数

在 `platform.json` 的 `defaults` 段加（需重启生效）：
```json
"max_heavy_tasks": 2,
"max_light_tasks": 4
```
heavy 建议值：1（机器弱/内存小）～ 3（32GB+ 且主要跑小样品）。

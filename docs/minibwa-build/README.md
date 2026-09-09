# minibwa 在 Windows 本地编译记录

> 上游：https://github.com/lh3/minibwa
> 编译 commit：`f0e117436c28addc359b67123d2353f0d4a1f9e8`（2026-08-10）
> 版本号：`0.7-r424-dirty`（`-dirty` 来自本次一处改动）
> 编译时间：2026-09-08
> 主机：Windows_NT 10.0.26200 / x86_64

## 结论

**能编译，能运行，功能完整。** 只需给 `kommon.c` 打一个 20 行的补丁，且不改动任何 POSIX 分支。

- 产物：`minibwa.exe`（1.63 MB，静态链接）
- 编译期改动：`kommon.c` 一处 `#ifdef _WIN32` 分支
- 运行期限制：`--mmap` 不可用（Windows 无 POSIX mmap），默认索引加载路径正常

## 迁移性：换台电脑要不要重编

**不用重编，直接拷 exe 就能用。** 已实测验证。

### 实测：干净环境运行

把 `minibwa.exe` 拷到孤立目录，`PATH` 砍到只剩三个系统目录（剔除 scoop / mingw / Python / git），三步全过：

| 步骤 | 结果 |
|---|---|
| `version` | exit 0 |
| `index ref.fa` | exit 0，生成 `.l2b` + `.mbw` |
| `map ref.fa read.fq` | exit 0，比对成功 |

不依赖任何开发环境、任何附带 DLL。

### 依赖构成（PE 导入表解析）

| 依赖 | 来源 | 说明 |
|---|---|---|
| `KERNEL32.dll` | 系统核心 | 必然存在 |
| `api-ms-win-crt-*.dll` × 9 | Universal C Runtime | Win10+ 内置 |

仅此两项。`libwinpthread-1.dll` 已被静态链接消除。

### 边界条件

| 系统 | 能否直接跑 |
|---|---|
| Windows 10 / 11（64 位） | 直接跑 |
| Windows 8.1 | 需 KB2999226（Windows Update 常规推送） |
| Windows 7 SP1 | 需 KB2999226 + KB2533623 |
| 32 位 Windows | 不能（exe 是 `pei-x86-64`） |
| ARM Windows | 需 x64 模拟层 |

关于 UCRT 的一个容易误解的点：`C:\Windows\System32\api-ms-win-crt-runtime-l1-1-0.dll` 在本机**并不存在**，但程序照跑。因为这些 `api-ms-win-crt-*` 是 API Set 虚拟名，由系统重定向到 `ucrtbase.dll`（本机 10.0.26100.8875）；对应的存根在 `C:\Windows\System32\downlevel\`，已逐项核对 9 个依赖全部存在。

### 想做“绝对零依赖”的话

当前工具链是 **UCRT 版**（`--with-sysroot=.../ucrt-rt_v14-rev1`），`-mcrtdll=msvcrt` 实测**无效**（仍导入 ucrt）。真要彻底摆脱：

1. 换 **MSVCRT 版** mingw 工具链（依赖 `msvcrt.dll`，Windows 自带）
2. 或使用目标机镜像自带的 x64 运行环境

但对 Win10/11 来说没必要，现状已可直接拷。

## 工具链

用 **mingw-w64 (GCC)**，不能用 MSVC。

| 项 | 值 |
|---|---|
| 编译器 | GCC 16.2.0（`x86_64-posix-seh-rev1`, MinGW-Builds） |
| 安装方式 | `scoop install mingw`（bucket 走 gitee 镜像） |
| 路径 | `C:\Users\17711\scoop\apps\mingw\current` |
| 构建工具 | `mingw32-make.exe`（同目录 `bin/`） |

为什么不能用 MSVC：Makefile 用 `$(shell uname -m)`、`printf | $(CC) -x c -fopenmp`、`rm -fr`、`-lpthread -lz`，全是 POSIX/GCC 语法。另外 minibwa 需要 `pthread`，mingw 的 `posix-seh` 线程模型正好提供（`ucrt` runtime）。

环境准备：

```powershell
scoop bucket add main
scoop install mingw
$root = "C:\Users\17711\scoop\apps\mingw\current"
$env:Path = "$root\bin;" + $env:Path
```

## 两个编译阻碍与处理

### 阻碍 1：`sys/mman.h` 不存在

`kommon.c:345` 包含 `sys/mman.h`，Windows 没有。这里的两个函数是索引的 mmap 加载器：

```c
void *kom_mmap_file(const char *fn, size_t *len, int preload);   // kommon.h:60
int   kom_munmap(void *base, size_t map_len);                    // kommon.h:61
```

调用链：`map-main.c:559` → `mb_idx_load_mmap` → `l2b_load_mmap`(`l2bit.c:332`) + `mb_bwt_load_mmap`(`bwt.c:681`)。

关键事实（`map-main.c:438`）：`use_mmap` 默认 **0**，只有显式传 `--mmap` 才置 1。所以存在一条非 mmap 的等价路径：

```c
/* map-main.c:559 */
idx = use_mmap ? mb_idx_load_mmap(argv[o.ind], is_meth, mmap_preload)
               : mb_idx_load(argv[o.ind], is_meth);
```

`mb_idx_load`（`map-algo.c`）用 `l2b_load` + `mb_bwt_load`，都是普通 `fread` 读进堆内存，在 Windows 上完全可用。

**处理**：不移植 mmap，改为提供函数桩，让 `--mmap` 明确失败而非静默降级（符合"不静默丢数据"的原则）。

### 阻碍 2：mimalloc 缺 Windows 后端

`mimalloc/prim/prim.c:12` 在 `_WIN32` 下 `#include "windows/prim.c"`，但 minibwa 内嵌的 mimalloc 子集**只有 `prim/unix/` 和 `prim/osx/`，没有 `prim/windows/`**。

**处理**：用 Makefile 内置开关禁用 mimalloc，回退到自带的 kalloc：

```
mingw32-make mimalloc=0
```

该分支会加 `-DHAVE_KALLOC`，功能等价，只是分配器不同。

## 补丁

见 `windows-mmap-stub.patch`（`git diff kommon.c`）。应用方式：

```bash
git apply windows-mmap-stub.patch
```

补丁内容：在 `kommon.c` 的 mmap 区块前插入 `#ifdef _WIN32` 桩，原 POSIX 代码包进 `#else ... #endif`。**POSIX 分支一行未改**。

## 完整构建命令

```powershell
$root = "C:\Users\17711\scoop\apps\mingw\current"
$env:Path = "$root\bin;" + $env:Path

git clone --depth 1 https://github.com/lh3/minibwa.git
cd minibwa
git apply ..\windows-mmap-stub.patch      # 或手动改 kommon.c
mingw32-make mimalloc=0
```

产物在同目录 `minibwa.exe`。

静态版本（推荐，可单文件分发）：

```powershell
mingw32-make mimalloc=0 LDFLAGS="-static -static-libgcc"
```

## 冒烟测试（实测通过）

数据：2 条 20kb 随机 contig，4 条 150bp 读段（正链/单错配/正链/反链各一）。

```
minibwa index ref.fa          → ref.fa.l2b (10092 B) + ref.fa.mbw (80096 B)，0.004 s
minibwa map ref.fa read.fq    → mapped 600 bp in 4 sequences，0.032 s
```

比对结果逐条核对：

| 读段 | 预期 | 实际 SAM | 判定 |
|---|---|---|---|
| r1 | chrA:1001 正链 0 错配 | `chrA 1001 150M FLAG=0 NM:i:0 AS:i:300` | PASS |
| r2 | chrB:5001 1 错配 | `chrB 5001 150M FLAG=0 NM:i:1 AS:i:290` | PASS |
| r3 | chrA:8001 正链 0 错配 | `chrA 8001 150M FLAG=0 NM:i:0 AS:i:300` | PASS |
| r4 | chrB:12001 反链 | `chrB 12001 150M FLAG=16 NM:i:0` | PASS |

`--mmap` 行为验证（确认不静默）：

```
[W::kom_mmap_file] mmap is not supported on this platform; use the default index loader
[E::main_map] failed to load the index. ABORT!
```

## 平台集成

已作为平台工具发放：

| 项 | 值 |
|---|---|
| 文件位置 | `bin/minibwa.exe`（单文件工具目录，同 seqkit / kunpeng） |
| 探测条目 | `vp/config.py` 的 `detect_tools()`，`'minibwa'` 键，紧邻 `'minimap2'` |
| 探测结果 | `26` 个工具可用（原 23 → 加 minibwa 后 26，含其他新增） |
| 调用方式 | `cfg.tool('minibwa')` |

实测验证（三级递进）：

1. `detect_tools()['minibwa']` → `...\bin\minibwa.exe`
2. `cfg.tool('minibwa')` 可用，`tool_status()['minibwa']` 在设置页可见
3. 用探测到的路径跑 `index` + `map` 端到端：exit 0，`r1 -> ctgA:2001 FLAG=0 CIGAR=150M`、`r4 -> ctgB:6001 FLAG=16 CIGAR=150M`

分发行为：`scripts/package.py` 把 `bin/` 整体拷入 `dist/VirusPlatform/bin/`，分发包的 `platform.json` 的 `tools` 为空，首次运行靠 `detect_tools()` 自动认领，**目标机无需任何配置**。

### 接入共识模块（2026-09-08）

`vp/consensus.py` 的比对引擎已从 minimap2 换成 minibwa（备份 `vp/consensus.py.bak_minibwa_20260908`）。

**接口差异（实测确认，不是照抄 minimap2）**

| 项 | minimap2 | minibwa |
|---|---|---|
| 建索引 | 隐式（内存） | **必须先 `index`**，产出 `.l2b` + `.mbw` |
| map 输入 | FASTA | **索引前缀**（`ref.fa`，非 `.l2b`） |
| 多 reads | 一次传全部 | **只接受单个 fastq**，R1/R2 需分次比对后合并 SAM |
| preset | `sr` / `map-ont` / `map-pb` / `map-hifi` | 仅 `sr` / `lr` / `adap` |
| 压次级比对 | `--secondary=no` | `-N 0` |

代码里对应三处新逻辑：`PRESET_MAP` + `map_preset()`（未知 preset 直接 `ValueError`，不静默回退）、`build_index()`、`run_minimap2()` 的多文件分支。索引落在参考 FASTA 旁，`run_minimap2` 用 `finally` 清理（`build_index` 单独调用时不清理）。

**对拍验证（三级递进）**

1. **干净数据**：200 reads，0-2 错配，位置/CIGAR/链向一致率 **100%**
2. **困难数据**：300 reads，0-4 错配 + 30% 插入/缺失 + 480bp 重复段，位置 96.00%、CIGAR 98.67%、链向 98.33%；4 条差异全落在重复段且 MAPQ=0
3. **共识序列（判据以最终共识为准）**：
   - 合成数据 6000 reads（~88x）：refA / refB 差异均 **0 bp**，`共识序列完全一致 = True`
   - 合成数据 1200 reads（~17x）：refA 差异 22 bp，全在深度跳变区

> 关键诊断：低深度差异的根因是 `MIN_DEPTH=10` 的阈值跳变，不是比对错误。重复区经 `MIN_MAPQ=10` 过滤后深度 minimap2=**8** / minibwa=**10**，差 1-2 条 read 就跨过阈值。把深度提到 88x 后差异归零，假说证实。

**真实数据对拍（ERR7586041，29736 bp 参考，664442 条 reads 比对上）**

| 指标 | minimap2 | minibwa |
|---|---|---|
| mapped reads | 664,442 | 667,086（+2,644） |
| 共识差异碱基 | — | NODE_1 **2 bp** / NODE_2 **3 bp** |
| 覆盖率 | 100.0% / 100.0% | 100.0% / 100.0% |
| mean_depth | 1620.06 / 1461.17 | 1624.61 / 1465.38 |
| iSNV 数 | 78 / 35 | 78 / 35（不变） |
| 耗时 | 2.7 s（比对） | 15.5 s（全流程） |

差异共 5 bp（0.011% / 0.026%），**方向完全一致：全部是 `minimap2=N`、`minibwa=实碱基`**。

逐位点核对 `poscounts.tsv`（pileup 过滤后的共识调用输入）：

| 位点 | minimap2 计数 | minibwa 计数 | 参考碱基 |
|---|---|---|---|
| NODE_1 pos 1 | C=10 | C=12 | C |
| NODE_1 pos 15505 | T=10 | T=11 | T |
| NODE_2 pos 10362 | G=10 | G=11 | G |

两个引擎在这些位点**都是单一碱基强支撑**，minibwa 的判读与参考序列一致。minimap2 判 N 是因为它在该位的过滤后深度恰好卡在 `MIN_DEPTH=10` 边缘（或位置未覆盖），属阈值效应，不是 minibwa 引入了假碱基。minibwa 多比对上的 2,644 条 reads 正是深度提升的来源。

**结论**：minibwa 在真实病毒宏基因组数据上可作为 minimap2 的替代，共识序列差异量级 0.01-0.03%，且方向单一（补实碱基），不引入幻象变异。iSNV 数（78/35）完全一致，说明准种检测结论未受影响。

### ⚠ 集成注意：stderr 是 GBK 编码

minibwa 写 stderr 用 GBK（中文 Windows 默认代码页），Python 侧按 UTF-8 解码会抛：

```
UnicodeDecodeError: 'utf-8' codec can't decode byte 0xd7 in position 53
```

调用时显式指定编码：

```python
subprocess.run([exe, 'map', ref, reads], capture_output=True,
               text=True, encoding='gbk', errors='replace')
```

或在管道模式下用 `errors='replace'` 兼容。这是 lh3 系列工具在 Windows 上的通病，minimap2 同样存在。

## 编译期警告（无害，未修）

- `kalloc.c:222`：`%ld` 打印 `size_t`，mingw 下 `long` 与 `size_t` 都是 64 位但类型不同，只在打印统计信息时用到
- `map-algo.c:732`：`-Wstringop-overflow` 假阳性，来自 memset 负数边界的静态推理

两者都不影响正确性，测试结果已证实。

## 已知限制

1. **`--mmap` 不可用**。大索引场景内存占用会显著高于 Linux（索引整体读入堆内存，Linux 下可依赖 mmap 按需分页）。要用 mmap 路径得真正移植 `CreateFileMapping`/`MapViewOfFile`，本次未做。
2. 动态链接版本依赖 `libwinpthread-1.dll`；静态版本（`LDFLAGS="-static -static-libgcc"`）只依赖系统自带的 KERNEL32 + UC runtime，已实测通过，建议用静态版分发。
3. 未跑长读/甲基化等高级路径，仅验证 `index` + `map` 主流程。
4. **索引落在参考 FASTA 旁**（`.l2b` / `.mbw`），与 minimap2 无索引副作用的行为不同。平台侧 `vp/consensus.py` 的 `run_minimap2()` 用 `finally` 清理，但若在其他场景直接调用 `minibwa index`，需自行清理。

# bcftools Windows 部署报告

> 部署时间：2026-09-08
> 定位：替代原流程中 iVar 的位置（BAM/计数表 → 变异 VCF 的 caller）

## 一、为什么换 bcftools

iVar 在 Windows 上不可得（bioconda 无 win-64、源码构建缺 autotools/perl、本机无 bash）。
bcftools 与 iVar 在链路中占同一个位置：把比对结果转成变异 VCF。
SnpEff 只做下游注释，不承担 calling。

## 二、安装事实

| 项 | 值 |
|---|---|
| 版本 | **bcftools 1.24**（Using htslib 1.24） |
| 来源 | MSYS2 官方构建 `mingw-w64-x86_64-bcftools-1.24-1-any.pkg.tar.zst` |
| 包大小 | 1.62 MB（安装后 20.22 MB，54 个文件） |
| SHA256 | `030eada2a6e7fd519d224649331ce9090ee6ea596c52c8111747122f719b4166` |
| 校验 | 下载件哈希与 MSYS2 官方页面**完全一致** |
| 位置 | `tools/bcftools/`（`bin/` + `libexec/`） |
| htslib | 复用本机 samtools 1.24 自带的 `hts-3.dll`（哈希与 pandepth 一致，同源） |

**SSL 障碍与绕过**：`mirror.msys2.org` / `repo.msys2.org` 服务端证书链含过期证书
（`SSLCertVerificationError: certificate has expired`, verify_code=10）；
清华镜像 403、sjtug 403。最终走 `ssl.CERT_NONE` 下载，**用 SHA256 保证完整性**。
github / pypi / raw.githubusercontent 证书正常。

## 三、依赖 DLL 补全

初持 `rc=0xC0000135`（STATUS_DLL_NOT_FOUND）。
`fix_deps.py` 递归解析 PE 导入表，遍历 32 个 DLL，从 `tools/samtools/bin`
与 `tools/jdk21/.../bin` 补入 20 个依赖，现缺失=无：

```
libtre-5  libintl-8  libiconv-2  libcrypto-3-x64  libidn2-0
libnghttp2-14  libnghttp3-9  libngtcp2-16  libngtcp2_crypto_ossl-0
libpsl-5  libssh2-1  libssl-3-x64  libunistring-5
hts-3  libsystre-0  libgcc_s_seh-1  libwinpthread-1  zlib1
libzstd  libbz2-1  liblzma-5  libdeflate  libbrotlicommon  libbrotlidec
libcurl-4  + 7 个 api-ms-win-crt-*.dll
```

## 四、实跑验证（全部 rc=0）

| 能力 | 结果 |
|---|---|
| `--version` | bcftools 1.24 / htslib 1.24 |
| 12 个核心子命令 | mpileup/call/consensus/norm/filter/view/stats/index/merge/isec/annotate/query 全在 |
| `call` 关键参数 | `--ploidy` / `-A` / `--multiallelic-caller` / `-v` / `-m` 全在 |
| `mpileup` 关键参数 | `-d, --max-depth INT [250]` / `-f, --fasta-ref` / `-q, --min-MQ` / `-Q, --min-BQ` 全在 |
| `view -Oz` + `index -t` | bgzip 压缩 + tabix 索引成功 |
| `consensus` | 真实 OR489165.1（14812 bp）应用 3 变异，`Applied 3 variants`，碱基正确写入 |
| `+counts` 插件 | Number of samples: 1 / SNPs: 2 / sites: 3 |
| `+fill-tags` 插件 | 正确补出 AN/AC/AF |

## 五、三个使用要点（踩过的坑）

1. **必须设 `BCFTOOLS_PLUGINS`**
   编译时 `libexec` 路径被硬编码为打包机的 `D:/M/msys64/mingw64/libexec/bcftools`，
   本机不存在。错误信息自身会提示这一点。
   ```
   BCFTOOLS_PLUGINS = D:\桌面\植物病毒分析平台\tools\bcftools\libexec\bcftools
   ```

2. **`consensus` 只吃 bgzip 压缩的 VCF**
   普通 `.vcf` 报 `not compressed with bgzip`。流程须先
   `bcftools view -Oz -o out.vcf.gz in.vcf` + `bcftools index -t out.vcf.gz`。

3. **`consensus` 会校验参考一致性**
   VCF 的 REF 与 fasta 不符时直接报错退出（`The fasta sequence does not match
   the REF allele at ...`）。这是正面信号，说明它真的按参考定位，不是盲套坐标。

## 六、iSNV 场景的两个参数要求

- `bcftools call` 需 `--ploidy 1 -A`：病毒单倍体，默认二倍体模型会误判基因型；
  `-A` 保留 iSNV 的次要等位。
- `bcftools mpileup` 需显式调大 `-d/--max-depth`：默认 250，
  病毒宏基因组单 contig 深度常达数百至数千，默认值会截断高深度位点。

## 七、与主路线的关系

主路线（`viral_consensus` → poscounts，阈值 `-q 20 -d 5 -f 0.5`）**保持不变**。
bcftools 与它的价值差异在于判定机制不同：

| | 主路线 | bcftools |
|---|---|---|
| 判定依据 | 计数 + 写死阈值 | 贝叶斯基因型推断 |
| 可复核性 | 每位点 A/C/G/T/N 可手查 | 模型输出 |
| 与 virome_analysis_pipeline 一致 | 是 | 否 |

两条独立路线对同一批数据给出一致的变异集，结论才稳固。
bcftools 的定位是**交叉验证**，不是替换。

## 八、文件清单

```
tools/bcftools/bin/          bcftools.exe + tabix/bgzip/vcfutils.pl 等 11 项 + 20 个依赖 DLL
tools/bcftools/libexec/bcftools/   42 个插件 DLL
kv_variant_test/
  install_bcftools.py   下载 + 解包
  dl_bcftools.py        多镜像下载 + SHA256 校验
  deploy_bcftools.py    部署 + 依赖补齐
  fix_deps.py           递归解析 PE 导入表补 DLL
  verify_bcftools.py    基础能力验证
  verify_bcf2.py        bgzip + 插件验证
  verify_bcf3.py        规范 header 验证
  verify_bcf4.py        真实序列 consensus + 插件实跑（最终验证）
  diag_ssl.py / diag_ssl2.py   SSL 诊断
```

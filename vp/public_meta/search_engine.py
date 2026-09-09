#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
🧬 GSA & SRA Global Species-Targeted Retrieval Engine (AI 极致提纯终极版)
核心修复：
1. GSA Excel 解析现在会正确抓取并拼接 'Age unit'。
2. 重写了 AI 提示词，强制剔除 Tissue 中的环境/状态修饰语 (如 under stress, young)。
3. 强化了 Location 的三级自动推理补全 (如 China:Yinchuan -> China, Ningxia, Yinchuan_AI)。
4. 清洗了各种形式的 missing/not collected 等无效信息。

（移植自 MMPV-RNA public_metadata_pipeline/gsa_sra.search.py，纯 requests 实现；
 文件写入统一改走平台 vp.utils.safe_open，sleep 抖动改用 secrets）
"""

import sys
if sys.platform == 'win32' and sys.stdout is not None:
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

import os
import re
import time
import json
import secrets
import argparse
import requests
import pandas as pd
from io import StringIO
from tqdm import tqdm
from urllib.parse import quote
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from openai import OpenAI
import warnings

warnings.filterwarnings("ignore")

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    from vp.utils import safe_open
else:
    from ..utils import safe_open

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Connection": "keep-alive"
}

def get_retry_session(retries=3, backoff_factor=1.5):
    session = requests.Session()
    retry_strategy = Retry(
        total=retries, backoff_factor=backoff_factor,
        status_forcelist=[403, 429, 500, 502, 503, 504]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy, pool_connections=20, pool_maxsize=20)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update(HEADERS)
    return session

# ==========================================
# 🟢 SRA 检索引擎
# ==========================================
class SRAEngine:
    def __init__(self, query, source, out_dir, detailed=False, ncbi_api=None):
        self.query = query
        self.source = source
        self.out_dir = out_dir
        self.detailed = detailed
        self.ncbi_api = ncbi_api
        self.session = get_retry_session()

    def fetch_runinfo(self):
        print("\n" + "="*65)
        mode_str = "详细模式 (XML解析)" if self.detailed else "极速模式 (基础信息)"
        print(f"🟢 启动 SRA 检索引擎 [{mode_str}]")
        print("="*65)

        term = f'"{self.query}"[Organism]'
        if self.source:
            if self.source.upper() == 'TRANSCRIPTOMIC':
                term += ' AND "biomol rna"[Properties]'
            elif self.source.upper() == 'GENOMIC':
                term += ' AND "biomol dna"[Properties]'
            else:
                term += f' AND "{self.source}"[Properties]'

        print(f"🧩 构建 SRA 检索逻辑: {term}")

        try:
            esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {"db": "sra", "term": term, "usehistory": "y", "retmode": "json"}
            if self.ncbi_api: params["api_key"] = self.ncbi_api
            res = self.session.get(esearch_url, params=params, timeout=30).json()
            count = int(res.get('esearchresult', {}).get('count', 0))
            if count == 0:
                print("⚠️ 未在 SRA 找到相关数据。")
                return pd.DataFrame()

            print(f"🎯 锁定 {count} 个 Run。正在下载 RunInfo 基础表...")
            webenv = res['esearchresult']['webenv']
            query_key = res['esearchresult']['querykey']

            efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params_csv = {"db": "sra", "query_key": query_key, "WebEnv": webenv, "rettype": "runinfo", "retmode": "text"}
            if self.ncbi_api: params_csv["api_key"] = self.ncbi_api
            res_csv = self.session.get(efetch_url, params=params_csv, timeout=60)
            df_sra = pd.read_csv(StringIO(res_csv.text))

            if self.detailed:
                print("🧬 [详细模式] 正在批量拉取 SRA XML 获取深层生物学特征...")
                params_xml = {"db": "sra", "query_key": query_key, "WebEnv": webenv, "rettype": "xml", "retmode": "text"}
                if self.ncbi_api: params_xml["api_key"] = self.ncbi_api
                res_xml = self.session.get(efetch_url, params=params_xml, timeout=120)
                xml_text = res_xml.text

                meta_dict = {}
                for pkg_match in re.finditer(r'<EXPERIMENT_PACKAGE>([\s\S]*?)</EXPERIMENT_PACKAGE>', xml_text):
                    pkg = pkg_match.group(1)
                    runs = re.findall(r'<RUN[^>]*accession="([E|S|D]RR\d+)"', pkg)

                    attrs = re.findall(r'<SAMPLE_ATTRIBUTE>\s*<TAG>([\s\S]*?)</TAG>\s*<VALUE>([\s\S]*?)</VALUE>', pkg, re.I)
                    attr_dict = {t.strip().lower(): v.strip() for t, v in attrs}

                    tissue = attr_dict.get('tissue', attr_dict.get('cell type', attr_dict.get('tissue type', pd.NA)))
                    loc = attr_dict.get('geo_loc_name', attr_dict.get('country', attr_dict.get('geographic location', pd.NA)))

                    age_parts = []
                    for k in ['age', 'dev_stage', 'development stage', 'growth stage']:
                        if k in attr_dict and attr_dict[k]:
                            age_parts.append(attr_dict[k])
                    stage = " | ".join(age_parts) if age_parts else pd.NA

                    for r in runs:
                        meta_dict[r] = {'Tissue': tissue, 'Age_GrowthStage': stage, 'Location': loc}

                df_sra['Tissue'] = df_sra['Run'].map(lambda x: meta_dict.get(x, {}).get('Tissue', pd.NA))
                df_sra['Age_GrowthStage'] = df_sra['Run'].map(lambda x: meta_dict.get(x, {}).get('Age_GrowthStage', pd.NA))
                df_sra['Location'] = df_sra['Run'].map(lambda x: meta_dict.get(x, {}).get('Location', pd.NA))

            return df_sra

        except Exception as e:
            print(f"❌ SRA 获取失败: {e}")
            return pd.DataFrame()

# ==========================================
# 🔵 GSA 检索引擎
# ==========================================
class GSAEngine:
    def __init__(self, query, source, out_dir, detailed=False, max_workers=5, progress_cb=None):
        self.query = str(query).strip()
        self.source_filter = str(source).strip() if source else None
        self.center_filter = "NGDC"
        self.detailed = detailed
        self.base_dir = os.path.join(out_dir, "GSA_Results")
        self.d_web = os.path.join(self.base_dir, "0_web_cache")
        self.d_xls = os.path.join(self.base_dir, "1_xls_cache")
        os.makedirs(self.d_web, exist_ok=True)
        if self.detailed:
            os.makedirs(self.d_xls, exist_ok=True)
        self.session = get_retry_session()
        self.max_workers = max_workers
        self.progress_cb = progress_cb  # callback(n_done, n_total, msg)

    def get_accession_list(self):
        url = "https://ngdc.cncb.ac.cn/gsa/search/getAccessionList"
        term = f'&quot;{self.query}&quot;[organism]'
        if self.source_filter: term = f'({term} AND &quot;{self.source_filter}&quot;[source])'
        term = f'({term} AND &quot;{self.center_filter}&quot;[center])'

        clean_term = term.replace('&quot;', '"')
        print(f"🧩 构建 GSA 检索逻辑: {clean_term}")

        payload = f"searchField=&searchTerm={quote(term)}&totalDatas=99999"
        headers = self.session.headers.copy()
        headers["Content-Type"] = "application/x-www-form-urlencoded"

        try:
            res = self.session.post(url, data=payload, headers=headers, timeout=45)
            return sorted(list(set(re.findall(r'(CR[RPAX]\d{6,})', res.text))))
        except Exception as e:
            return []

    def parse_all_from_html(self, html, acc):
        feat = {
            "ReleaseDate": pd.NA, "Organization": pd.NA, "CRA": pd.NA, "CRX": pd.NA,
            "PRJ": pd.NA, "SAMC": pd.NA, "Platform": pd.NA, "TaxID": pd.NA, "ScientificName": pd.NA,
            "LibraryStrategy": "Not_Provided", "LibrarySource": "Not_Provided", "RunRecords": []
        }
        if not html: return feat

        def extract_text(pat):
            m = re.search(pat, html, re.I)
            return re.sub(r'<[^>]+>', '', m.group(1)).strip() if m else pd.NA

        feat["ReleaseDate"] = extract_text(r'<th[^>]*>(?:发布日期|Release date)</th>\s*<td[^>]*>([\s\S]*?)</td>')
        feat["Organization"] = extract_text(r'<th[^>]*>(?:所属单位|Organization)</th>\s*<td[^>]*>([\s\S]*?)</td>')
        feat["ScientificName"] = extract_text(r'wwwtax\.cgi\?id=\d+"[^>]*>([\s\S]*?)</a>')
        feat["CRX"] = extract_text(r'<th[^>]*>(?:实验编号|Accession)</th>\s*<td[^>]*>.*?(CRX\d+)[\s\S]*?</td>')
        feat["PRJ"] = extract_text(r'<th[^>]*>(?:项目编号|BioProject)</th>\s*<td[^>]*>.*?(PRJ[A-Z]*\d+)[\s\S]*?</td>')
        feat["SAMC"] = extract_text(r'<th[^>]*>(?:样本编号|BioSample)</th>\s*<td[^>]*>.*?(SAM[A-Z]*\d+)[\s\S]*?</td>')
        feat["Platform"] = extract_text(r'<th[^>]*>(?:测序平台|Platform)</th>\s*<td[^>]*>([\s\S]*?)</td>')

        cra_m = re.search(r'(CRA\d{6,})', html)
        if cra_m: feat["CRA"] = cra_m.group(1)

        lib_table_m = re.search(r'<th[^>]*>(?:建库信息|Library)</th>\s*<td[^>]*>[\s\S]*?<table[^>]*>([\s\S]*?)</table>', html, re.I)
        if lib_table_m:
            tds = re.findall(r'<td[^>]*>([\s\S]*?)</td>', lib_table_m.group(1), re.I)
            if len(tds) >= 4:
                feat["LibraryStrategy"] = re.sub(r'<[^>]+>', '', tds[2]).strip()
                feat["LibrarySource"] = re.sub(r'<[^>]+>', '', tds[3]).strip()

        file_block_m = re.findall(r'<a href="[^"]*browse/(CRA\d+)/(CRR\d+)"[^>]*>[\s\S]*?</a>\s*</td>\s*<td[^>]*>([\s\S]*?)</td>\s*<td[^>]*>([\s\S]*?)</td>', html, re.I)
        if file_block_m:
            feat["CRA"] = file_block_m[0][0]
            for match in file_block_m:
                if re.findall(r'\.(?:fq|fastq|bam|sra)', match[2], re.I): feat["RunRecords"].append({"Run": match[1]})
        else:
            file_m = re.search(rf'<a[^>]*>{acc}</a>[\s\S]*?</td>\s*<td[^>]*>([\s\S]*?)</td>\s*<td[^>]*>([\s\S]*?)</td>', html, re.I)
            if file_m and re.findall(rf'\.(?:fq|fastq|bam|sra)', file_m.group(1), re.I): feat["RunRecords"].append({"Run": acc})

        return feat

    def fetch_excel_attributes(self, cra_list):
        meta_dict = {}
        print("\n📊 [详细模式] 正在向 CNCB 提交 Excel 级联解析请求...")
        for cra in tqdm(cra_list, desc="解析 Excel 附件"):
            if pd.isna(cra): continue
            xls_path = os.path.join(self.d_xls, str(cra) + '.xlsx')

            if not os.path.exists(xls_path):
                try:
                    rx = self.session.post("https://ngdc.cncb.ac.cn/gsa/file/exportExcelFile", data={"type": 3, "dlAcession": cra}, timeout=30)
                    if len(rx.content) > 1000:
                        with safe_open(xls_path, 'wb') as f: f.write(rx.content)
                except: pass

            if os.path.exists(xls_path):
                try:
                    xls = pd.ExcelFile(xls_path, engine='openpyxl')
                    if 'Sample' in xls.sheet_names:
                        df_samp = pd.read_excel(xls, sheet_name='Sample')
                        cols = df_samp.columns.str.lower().str.strip()
                        df_samp.columns = cols

                        acc_col = 'accession' if 'accession' in cols else 'biosample accession'
                        tissue_col = 'tissue' if 'tissue' in cols else 'organism part'

                        for _, row in df_samp.iterrows():
                            sam = row.get(acc_col)
                            if pd.isna(sam): continue

                            # 【核心修复】智能拼接 age 和 age unit
                            age_val = row.get('age')
                            age_unit = row.get('age unit')
                            stage_val = row.get('dev stage')

                            age_str = ""
                            if pd.notna(age_val):
                                age_str = str(age_val).strip()
                                if pd.notna(age_unit):
                                    age_str += f" {str(age_unit).strip()}"

                            stage_parts = []
                            if age_str: stage_parts.append(age_str)
                            if pd.notna(stage_val): stage_parts.append(str(stage_val).strip())

                            meta_dict[sam] = {
                                'Tissue': row.get(tissue_col, pd.NA),
                                'Age_GrowthStage': " | ".join(stage_parts) if stage_parts else pd.NA,
                                'Location': row.get('geographic location', pd.NA)
                            }
                except: pass
        return meta_dict

    def fetch_download_urls(self):
        """从已缓存的 Run Excel 表提取下载链接和 MD5，返回 {Run: {url_1, md5_1, url_2, md5_2}}"""
        url_dict = {}
        if not os.path.isdir(self.d_xls): return url_dict
        xls_files = [f for f in os.listdir(self.d_xls) if f.endswith('.xlsx')]
        if not xls_files: return url_dict
        print(f"\n📎 [详细模式] 正在从 {len(xls_files)} 个 Run Excel 提取下载链接...")
        for fname in tqdm(xls_files, desc="提取下载链接"):
            xls_path = os.path.join(self.d_xls, fname)
            try:
                xls = pd.ExcelFile(xls_path, engine='openpyxl')
                if 'Run' not in xls.sheet_names: continue
                df_run = pd.read_excel(xls, sheet_name='Run')
                for _, row in df_run.iterrows():
                    acc = row.get('Accession')
                    if pd.isna(acc): continue
                    acc = str(acc).strip()
                    entry = {}
                    for col in df_run.columns:
                        cl = col.strip().lower()
                        if 'download read file1' in cl or 'download  read file1' in cl:
                            v = row.get(col)
                            entry['url_1'] = str(v).strip() if pd.notna(v) and str(v).strip() else None
                        elif 'read file1 md5' in cl:
                            v = row.get(col)
                            entry['md5_1'] = str(v).strip() if pd.notna(v) and str(v).strip() else None
                        elif 'download read file2' in cl or 'download  read file2' in cl:
                            v = row.get(col)
                            entry['url_2'] = str(v).strip() if pd.notna(v) and str(v).strip() else None
                        elif 'read file2 md5' in cl:
                            v = row.get(col)
                            entry['md5_2'] = str(v).strip() if pd.notna(v) and str(v).strip() else None
                    if entry.get('url_1'):
                        url_dict[acc] = entry
            except: pass
        print(f"  ✓ 提取到 {len(url_dict)} 条下载链接")
        return url_dict

    def _fetch_one_acc(self, acc):
        """Fetch and parse a single accession page. Thread-safe."""
        web_cache_f = os.path.join(self.d_web, str(acc) + '.html')
        try:
            if not os.path.exists(web_cache_f):
                time.sleep(0.2 + secrets.randbelow(401) / 1000)
                # Each thread uses its own session for safety
                sess = get_retry_session()
                res_web = sess.get("https://ngdc.cncb.ac.cn/gsa/search?searchTerm=" + str(acc), timeout=30)
                # Atomic-ish write
                tmp = web_cache_f + '.tmp'
                with safe_open(tmp, 'wt') as f:
                    f.write(res_web.text)
                os.replace(tmp, web_cache_f)
            with open(web_cache_f, 'r', encoding='utf-8') as f:
                html = f.read()

            wf = self.parse_all_from_html(html, acc)
            if self.query.lower() not in str(wf.get("ScientificName", "")).lower():
                return None

            records = []
            cra_id = wf.get("CRA", "UNKNOWN_CRA")
            for r in wf.get("RunRecords", [{"Run": acc}]):
                records.append({
                    "Run": r["Run"], "ReleaseDate": wf.get("ReleaseDate", pd.NA),
                    "LibraryStrategy": wf.get("LibraryStrategy"), "LibrarySource": wf.get("LibrarySource"),
                    "BioProject": wf.get("PRJ", pd.NA), "BioSample": wf.get("SAMC", pd.NA),
                    "Platform": wf.get("Platform", pd.NA), "CenterName": wf.get("Organization", pd.NA),
                    "ScientificName": wf.get("ScientificName", pd.NA),
                })
            return (cra_id, records)
        except Exception as e:
            return None

    def fetch_gsa(self):
        from concurrent.futures import ThreadPoolExecutor, as_completed
        print("\n" + "="*65)
        mode_str = "详细模式 (含Excel解析)" if self.detailed else "极速模式 (基础信息)"
        print(f"🔵 启动 GSA 检索引擎 [{mode_str}] [{self.max_workers} workers]")
        print("="*65)

        acc_list = self.get_accession_list()
        if not acc_list:
            return pd.DataFrame()

        total = len(acc_list)
        print(f"🎯 锁定 {total} 个 accession，并行抓取中...")
        if self.progress_cb:
            self.progress_cb(0, total, f"GSA: 0/{total}")

        all_records = []
        unique_cras = set()
        n_done = 0
        last_report = 0

        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futures = {pool.submit(self._fetch_one_acc, acc): acc for acc in acc_list}
            for fut in as_completed(futures):
                n_done += 1
                result = fut.result()
                if result:
                    cra_id, records = result
                    unique_cras.add(cra_id)
                    all_records.extend(records)
                # Report progress every 5% or every 5 items
                if self.progress_cb and (n_done - last_report >= max(1, total // 20)):
                    self.progress_cb(n_done, total, f"GSA: {n_done}/{total} ({len(all_records)} runs)")
                    last_report = n_done

        if self.progress_cb:
            self.progress_cb(total, total, f"GSA: {total}/{total} ({len(all_records)} runs)")

        df_gsa = pd.DataFrame(all_records)

        if self.detailed and not df_gsa.empty and unique_cras:
            excel_meta = self.fetch_excel_attributes(unique_cras)
            df_gsa['Tissue'] = df_gsa['BioSample'].map(lambda x: excel_meta.get(x, {}).get('Tissue', pd.NA))
            df_gsa['Age_GrowthStage'] = df_gsa['BioSample'].map(lambda x: excel_meta.get(x, {}).get('Age_GrowthStage', pd.NA))
            df_gsa['Location'] = df_gsa['BioSample'].map(lambda x: excel_meta.get(x, {}).get('Location', pd.NA))

        return df_gsa

# ==========================================
# 🤖 AI 智能清洗与规范化引擎 (全新提纯规则)
# ==========================================
def run_ai_sanitizer(df, api_key, api_base, model):
    if df.empty or not api_key: return df
    print("\n🤖 启动 AI 洗髓引擎 (深度规范 Tissue, Age_GrowthStage, Location)...")

    client = OpenAI(api_key=api_key, base_url=api_base)
    cols = ['Location', 'Tissue', 'Age_GrowthStage']
    for c in cols:
        if c not in df.columns: df[c] = pd.NA

    prompt = """你是一个极其严苛的生命科学数据清理程序。请严格按以下规则处理输入的JSON，返回清洗后的JSON：
1. **全局去杂**：遇到 'not applicable', 'not collected', 'missing', 'N/A', 'nan', 或只包含无意义编号(如 'tissue9')，一律清空替换为 'Not_Provided'。
2. **Tissue (核心剥离)**：强制剔除所有的状态、环境、时间、物种等修饰语（例如：将 "flower under drought stress in the stage 2", "leaf of Wolfberry", "young leaves" 统统剥离）。你【只能】输出最纯粹的单数英文核心组织词，例如：leaf, root, stem, flower, fruit, seed, anther, stamen, pistil, whole plant 等。绝对不要保留多余信息！
3. **Age_GrowthStage (凝练)**：剥离多余符号和无效重复，合并年龄和时期信息。例如将 "3 | Young Fruit" 规范为 "3 years, young fruit stage"；将 "not applicable | missing" 清除为 "Not_Provided"。
4. **Location (三级正交化)**：必须凭借你的地理知识，强制补全为【国家, 省/州, 市/县_AI】的标准格式（如 "China:Yinchuan" 或 "China:yinchuan" 必须补全为 "China, Ningxia, Yinchuan_AI"；"China:Qaidam Basin" 补全为 "China, Qinghai, Qaidam Basin_AI"）。实在缺失的层级用 Unknown 补全。
仅返回合法 JSON，勿带 ```json 标记。绝对禁止凭空捏造数据！"""

    for idx, row in tqdm(df.iterrows(), total=len(df), desc="AI Sanitization"):
        if pd.isna(row.get('Location')) and pd.isna(row.get('Tissue')) and pd.isna(row.get('Age_GrowthStage')): continue
        target = {k: str(row.get(k, "Not_Provided")) for k in cols}

        try:
            res = client.chat.completions.create(
                model=model,
                messages=[{"role": "system", "content": prompt}, {"role": "user", "content": json.dumps(target, ensure_ascii=False)}],
                temperature=0.1
            )
            clean_str = re.sub(r"^```(?:json)?\s*|\s*```$", "", res.choices[0].message.content.strip(), flags=re.IGNORECASE)
            clean_data = json.loads(clean_str)

            for k in cols:
                if k in clean_data:
                    val = str(clean_data[k]).strip()
                    if val.lower() in ["", "none", "nan", "not_provided", "unknown"]:
                        df.at[idx, k] = pd.NA
                    else:
                        df.at[idx, k] = val
        except: pass
    return df

# ==========================================
# 📊 FileSize_GB 提取
# ==========================================
def _add_filesize_sra(df):
    """SRA: map size_MB → GB."""
    if df is None or df.empty:
        return
    if "size_MB" in df.columns:
        def _to_gb(x):
            try:
                return f"{float(x)/1024:.1f}"
            except (ValueError, TypeError):
                return ""
        df["FileSize_GB"] = df["size_MB"].apply(_to_gb)


def _add_filesize_gsa(df, out_dir):
    """GSA: parse bytes from cached Excel Run sheet → GB."""
    if df is None or df.empty or "Run" not in df.columns:
        return
    xls_dir = os.path.join(out_dir, "GSA_Results", "1_xls_cache")
    df["FileSize_GB"] = ""
    if not os.path.isdir(xls_dir):
        return
    run_sizes = {}
    for fname in os.listdir(xls_dir):
        if not fname.endswith(".xlsx"):
            continue
        xls_path = os.path.join(xls_dir, fname)
        try:
            xls = pd.ExcelFile(xls_path)
            if "Run" not in xls.sheet_names:
                continue
            df_run = pd.read_excel(xls, sheet_name="Run")
            for _, rrow in df_run.iterrows():
                acc = str(rrow.get("Accession", "")).strip()
                if not acc:
                    continue
                total_bytes = 0; file_count = 0
                for col in df_run.columns:
                    cl = col.strip().lower()
                    if "read filename" in cl and "md5" not in cl:
                        val = str(rrow[col]).strip() if pd.notna(rrow[col]) else ""
                        m = re.search(r'\((\d+)\s*bytes?\)', val, re.I)
                        if m:
                            total_bytes += int(m.group(1)); file_count += 1
                if file_count > 0:
                    total_gb = total_bytes / 1073741824
                    run_sizes[acc] = f"{total_gb:.1f}"
        except Exception:
            pass
    for idx, row in df.iterrows():
        run = str(row["Run"]).strip()
        if run in run_sizes:
            df.at[idx, "FileSize_GB"] = run_sizes[run]


# ==========================================
# 🤖 AI Summary (本地)
# ==========================================
def _generate_ai_summary(df, api_key, model, out_dir):
    """Generate SCI writing summary via DeepSeek API."""
    if df.empty or not api_key:
        return
    from collections import Counter

    def _safe_col(name):
        c = name if name in df.columns else None
        if c:
            series = df[c].astype(str).str.strip()
            series = series[~series.isin(["", "NA", "N/A", "Not_Provided", "nan", "None"])]
            return Counter(series)
        return Counter()

    db_counts = _safe_col("Database")
    species_counts = _safe_col("ScientificName")
    tissue_counts = _safe_col("Tissue")
    location_counts = _safe_col("Location")
    center_counts = _safe_col("CenterName")
    age_counts = _safe_col("Age_GrowthStage")
    total = len(df)

    # Data volume
    vol_str = "N/A"
    if "FileSize_GB" in df.columns:
        try:
            gb_vals = pd.to_numeric(df["FileSize_GB"], errors="coerce").dropna()
            if len(gb_vals) > 0:
                vol_str = f"{gb_vals.sum():.1f} GB total, avg {gb_vals.mean():.1f} GB/run"
        except Exception:
            pass

    # Collection years
    years = set()
    for c in ["ReleaseDate", "CollectionDate"]:
        if c in df.columns:
            for v in df[c]:
                parts = str(v).strip().replace("/","-").split("-")
                if parts[0].isdigit() and len(parts[0]) == 4:
                    years.add(parts[0])

    stats_text = f"""Total records: {total}
Data volume: {vol_str}
Database: {', '.join(f'{k}({v})' for k,v in db_counts.most_common())}
Species: {', '.join(f'{k}({v})' for k,v in species_counts.most_common())}
Tissues: {', '.join(f'{k}({v})' for k,v in tissue_counts.most_common())}
Locations (top10): {', '.join(f'{k}({v})' for k,v in location_counts.most_common(10))}
Institutions (top10): {', '.join(f'{k}({v})' for k,v in center_counts.most_common(10))}
Growth stages: {', '.join(f'{k}({v})' for k,v in age_counts.most_common(5))}
Years: {', '.join(sorted(years)) if years else 'N/A'}"""

    prompt = f"""You are a scientific writer. Based on these metadata statistics, write a concise paragraph (150-250 words) for the Data Collection section of a virome/metagenomics paper.
{stats_text}
Requirements: formal scientific English, past tense. Include total runs, database split, data volume, species, tissues, locations, time span, key institutions, sequencing type.
Output ONLY the paragraph, no markdown."""

    try:
        client = OpenAI(api_key=api_key, base_url="https://api.deepseek.com")
        kwargs = {"model": model, "messages": [
            {"role": "system", "content": "You are a scientific writer. Output only the paragraph."},
            {"role": "user", "content": prompt}
        ]}
        kwargs["temperature"] = 0.3
        response = client.chat.completions.create(**kwargs)
        summary = response.choices[0].message.content or ""
        print(f"\n{'='*60}")
        print("📝 AI Summary (Data Collection)")
        print(f"{'='*60}")
        print(summary)
        # Save to file
        summary_file = os.path.join(out_dir, "AI_Summary.txt")
        with safe_open(summary_file, 'wt') as f:
            f.write(summary)
        print(f"\n📁 Summary saved: {summary_file}")
    except Exception as e:
        print(f"\n⚠️ AI Summary failed: {e}")


# ==========================================
# 🌍 数据大一统合并
# ==========================================
def merge_results(df_sra, df_gsa, out_dir, detailed, api_key, api_base, model, do_summary=False):
    if not df_sra.empty:
        df_sra['Database'] = 'SRA'
        _add_filesize_sra(df_sra)
    else: df_sra = pd.DataFrame()

    if not df_gsa.empty:
        df_gsa['Database'] = 'GSA'
        _add_filesize_gsa(df_gsa, out_dir)
    else: df_gsa = pd.DataFrame()

    df_merged = pd.concat([df_sra, df_gsa], ignore_index=True)
    if df_merged.empty: return

    if detailed and api_key:
        df_merged = run_ai_sanitizer(df_merged, api_key, api_base, model)

    target_cols = ['Database', 'Run', 'BioProject', 'BioSample', 'ScientificName']
    if detailed:
        target_cols.extend(['Tissue', 'Age_GrowthStage', 'Location'])
    target_cols.extend(['LibraryStrategy', 'LibrarySource', 'Platform', 'CenterName', 'ReleaseDate', 'FileSize_GB'])

    final_cols = [c for c in target_cols if c in df_merged.columns]

    df_final = df_merged[final_cols].drop_duplicates(subset=['Run'])

    out_file = os.path.join(out_dir, "SRA_GSA_Merged_Final.csv")
    df_final.to_csv(out_file, index=False, encoding='utf-8-sig')

    print(f"\n🎉 大一统合并成功！共汇总 {len(df_final)} 条 Run 记录。")
    print(f"📁 最终输出路径: {out_file}")

def _default_outdir(query):
    """默认输出根目录：meta_search/<物种slug>/search"""
    slug = re.sub(r'[^\w\-]+', '_', str(query or 'query')).strip('_') or 'query'
    try:
        from vp.config import DIRS
        root = DIRS['meta_search']
    except Exception:
        root = os.path.join(os.getcwd(), 'meta_search')
    return os.path.join(root, slug, 'search')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SRA + GSA 双引擎物种检索与提纯工具 (双模式版)")
    parser.add_argument("-q", "--query", required=True, help="输入物种拉丁名 (如 'Lycium barbarum')")
    parser.add_argument("-s", "--source", default="TRANSCRIPTOMIC", help="限制测序类型 (默认: TRANSCRIPTOMIC，'All' 则不限)")
    parser.add_argument("-o", "--outdir", default=None, help="输出根目录（默认 meta_search/<物种>/search）")

    parser.add_argument("--db", default="both", choices=["sra", "gsa", "both"],
                        help="目标数据库: sra (仅NCBI), gsa (仅CNCB), both (默认, 两者)")
    parser.add_argument("--no-detailed", action="store_true", help="关闭详细模式 (默认开启)")
    parser.add_argument("--workers", type=int, default=5,
                        help="GSA 并发爬取线程数 (1-20, 默认 5; 越高越快但易触发限流)")
    parser.add_argument("--summary", action="store_true", help="搜索完成后生成 AI 统计段落 (需 --deepseek-api)")
    parser.add_argument("--ncbi-api", help="NCBI E-utilities API Key (提升速率)")
    parser.add_argument("--deepseek-api", help="DeepSeek API Key (AI 智能清洗及AI统计)")
    parser.add_argument("--deepseek-model", default="deepseek-v4-flash",
                        choices=["deepseek-chat", "deepseek-reasoner", "deepseek-v4-flash", "deepseek-v4-pro"],
                        help="DeepSeek 模型")

    args = parser.parse_args()
    args.detailed = not args.no_detailed  # default True, --no-detailed to disable
    # Normalize "All" / "" source to None (no filter)
    if args.source and args.source.strip().lower() in ("all", ""):
        args.source = None
    if not args.outdir:
        args.outdir = _default_outdir(args.query)
    os.makedirs(args.outdir, exist_ok=True)

    df_sra = pd.DataFrame()
    df_gsa = pd.DataFrame()

    if args.db in ("sra", "both"):
        sra_engine = SRAEngine(args.query, args.source, args.outdir, detailed=args.detailed, ncbi_api=args.ncbi_api)
        df_sra = sra_engine.fetch_runinfo()

    if args.db in ("gsa", "both"):
        gsa_engine = GSAEngine(args.query, args.source, args.outdir, detailed=args.detailed,
                               max_workers=max(1, min(20, args.workers)))
        df_gsa = gsa_engine.fetch_gsa()

    merge_results(df_sra, df_gsa, args.outdir, args.detailed, args.deepseek_api, "https://api.deepseek.com", args.deepseek_model)

    # AI Summary
    if args.summary and args.deepseek_api:
        df_final = pd.read_csv(os.path.join(args.outdir, "SRA_GSA_Merged_Final.csv"))
        _generate_ai_summary(df_final, args.deepseek_api, args.deepseek_model, args.outdir)

    # 生成下载链接文件
    import csv as csv_mod
    links = []
    # SRA: RunInfo download_path 列
    if not df_sra.empty and 'download_path' in df_sra.columns:
        for _, row in df_sra.iterrows():
            url = row.get('download_path')
            if pd.notna(url) and str(url).strip():
                links.append({'Run': row['Run'], 'Database': 'SRA',
                              'url_1': str(url).strip(), 'md5_1': '', 'url_2': '', 'md5_2': ''})
    # GSA: Excel Run 表 (遍历所有已缓存 Excel)
    if not df_gsa.empty and gsa_engine.detailed:
        gsa_urls = gsa_engine.fetch_download_urls()
        for _, row in df_gsa.iterrows():
            info = gsa_urls.get(row['Run'], {})
            if info:
                links.append({'Run': row['Run'], 'Database': 'GSA',
                              'url_1': info.get('url_1', ''), 'md5_1': info.get('md5_1', ''),
                              'url_2': info.get('url_2', ''), 'md5_2': info.get('md5_2', '')})
    if links:
        link_file = os.path.join(args.outdir, "download_links.csv")
        with safe_open(link_file, 'wt') as f:
            w = csv_mod.DictWriter(f, fieldnames=['Run', 'Database', 'url_1', 'md5_1', 'url_2', 'md5_2'])
            w.writeheader()
            w.writerows(links)
        print(f"📎 下载链接已保存: {link_file} ({len(links)} 条)")

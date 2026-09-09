#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
🌍 Global Metadata Unification Engine (SRA + GSA 大一统智能解析引擎 - v12.11 终极完备版)
包含：13大核心字段提取、终极 AI 防篡改指令、Location 强正交化、DeepSeek V4 兼容、GSA 极速免疫防踢爬虫、文献精准溯源。

（移植自 MMPV-RNA public_metadata_pipeline/gsa_sra.info.py；Windows 适配：
 SRA XML 下载由 EDirect 管道改为 E-utilities HTTP；文件写入统一走
 vp.utils.safe_open；HTTP 统一走 requests.Session 实例；datavzrd 缺失时
 干净跳过；sleep 抖动改用 secrets）
"""

import os
import sys
if sys.platform == 'win32' and sys.stdout is not None:
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

import subprocess
import json
import time
import re
import shlex
import argparse
import warnings
import datetime
import secrets
import requests
import builtins
import concurrent.futures
import pandas as pd
import xmltodict
from openai import OpenAI
from tqdm import tqdm
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    from vp.utils import safe_open
else:
    from ..utils import safe_open

# 屏蔽警告
warnings.filterwarnings("ignore", category=UserWarning, module="scipy")
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")
warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")

# 拦截 print，防进度条撕裂
def safe_print(*args, **kwargs):
    kwargs.pop('end', None)
    msg = " ".join(map(str, args))
    tqdm.write(msg)
builtins.print = safe_print

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36"
}

# =========================================================================
# 🧠 终极核心数据结构与全量 AI 提示词库
# =========================================================================

# SRA 专属的 18 个提取维度
CORE_18_SRA = ['query_id', 'Run', 'ReleaseDate', 'CollectionDate', 'Location', 'Source', 'Tissue', 'Age_GrowthStage', 'LibrarySource', 'SRAStudy', 'BioProject', 'ProjectID', 'Sample', 'BioSample', 'ScientificName', 'TaxID', 'SampleName', 'CenterName', 'PMID']
DESC_COLS = ['Study_title', 'Study_abstract']
GSA_AI_10 = ['Run', 'CollectionDate', 'Location', 'Source', 'Tissue', 'Age_GrowthStage', 'ScientificName', 'LibrarySource', 'CenterName', 'BioProject']

# 大一统输出 14 列
FINAL_14 = ['Run', 'ReleaseDate', 'CollectionDate', 'Location', 'Source', 'Tissue', 'Age_GrowthStage', 'ScientificName', 'TaxID', 'LibrarySource', 'CenterName', 'BioProject', 'BioSample', 'PMID']

# ----------------- 【SRA 专用推断提示词（全域信息综合推断版）】 -----------------
PROMPT_SRA_INFER = """你是一个专业的生物信息学数据分析专家。我将提供一段来自NCBI SRA数据库的元数据记录（植物测序实验）。
【核心推断法则 - 全域信息洞察】：
请仔细阅读并基于**所有记录信息**（包括 TITLE, EXPERIMENT_DESIGN, SAMPLE_ATTRIBUTES, ABSTRACT 等全文各个角落）进行综合全局推断，**绝对不要仅仅局限于提取 STUDY_ABSTRACT。** 对于提交极不规范的样本，生境、组织和发育阶段无处不在。

【数据清洗绝对法则】：
1. 遇到 'not collected', 'missing', 'N/A', 'none', 'unknown' 等无意义占位符，直接判定为缺失，必须输出 'Not_Provided'！
2. Tissue 必须使用全小写、单数形式（如 leaves 改为 leaf, Fruits 改为 fruit）。

**【Location 地理溯源最高法则】（格式统一为：国家, 省/州, 市/县_AI）：**
在推断植物的真实生长/采集地点时，必须严格遵循以下优先级寻找线索：
1. **最高优先级（原生定位）**：直接寻找 `geo_loc_name` 或明确的地点标签。
2. **次高优先级（材料提供者）**：寻找 `biomaterial_provider`。如果样本是由特定机构提供的（如 Bonn University Botanical Garden），请推断其所在的国家、省份和城市。
3. **最低兜底（提交机构）**：如果上述均无，且全文无采样地描述，才允许使用 `CenterName` 推测国家、省份和城市。
注意：输出格式必须严格为【国家, 省/州, 市/县_AI】（缺失层级用 Unknown 补全，例如 'Germany, Unknown, Bonn_AI'）。

**必须要输出的字段解释（全部提取，缺失填 Not_Provided）：**
- Run: SRA Run编号 (如 SRR...)
- ReleaseDate: 数据最早入库/公开的确认时间
- CollectionDate: 样本真实的物理采集时间
- Location: 采集地点 (严格遵循三级地理格式，带 _AI)
- Source: 样本来源 (野生 wild、栽培 cultivated 等)
- Tissue: 组织/部位 (全小写单数，如 leaf, root 等)
- Age_GrowthStage: 生长年限或发育阶段 (如 mature, seedling 等)
- LibrarySource: 文库测序类型 (如 TRANSCRIPTOMIC, GENOMIC 等)
- SRAStudy: SRP编号
- BioProject: PRJNA编号
- ProjectID: 项目相关其他ID
- Sample: SRS编号
- BioSample: SAMN编号
- ScientificName: 物种科学名
- TaxID: 物种Taxonomy号
- SampleName: 提交者定义的样本名
- CenterName: 测序/提交中心机构名称
- PMID: PubMed论文编号

**必须输出包含这 18 个键名的合法 JSON 字符串（勿带```json标记）：**
{
  "Run": "...", "ReleaseDate": "...", "CollectionDate": "...", "Location": "...", "Source": "...", 
  "Tissue": "...", "Age_GrowthStage": "...", "LibrarySource": "...", "SRAStudy": "...", 
  "BioProject": "...", "ProjectID": "...", "Sample": "...", "BioSample": "...", 
  "ScientificName": "...", "TaxID": "...", "SampleName": "...", "CenterName": "...", "PMID": "..."
}"""

# ----------------- 【SRA 仲裁提示词的权威修正】 -----------------
PROMPT_SRA_ARBITRATE = """你现在是数据质量最终仲裁官。对比以下同一实验的两份数据：
【数据 A】：本地规则精准提取结果 (绝不含幻觉)。
【数据 B】：AI 全局信息的语义推断结果。

【仲裁绝对指令】：
0. **【清洗过滤】**：遇到 'not collected', 'missing', 'N/A' 等无意义词，直接清除并输出 'Not_Provided'。
1. **【精确 ID 类】**（Run, SRAStudy, BioProject, Sample, BioSample, CenterName, LibrarySource, ProjectID, SampleName, ScientificName, TaxID, ReleaseDate, PMID）：必须采用【数据 A】中的值。如果 A 缺失，才采用 B。
2. **【生物特征类】**（Tissue, Age_GrowthStage, Source, CollectionDate）：**【致命警告】如果【数据 A】中已有明确的值（如 wild, stage 4, 具体日期等），必须 100% 信任 A，绝对不许被 B 篡改或覆盖！** 只有在 A 彻底为空或 Not_Provided 时，才允许采用 B 挖掘出的数据。注意：采用 A 的值时，必须将 Tissue 强行转为小写单数。
3. **【Location 仲裁】**：必须强制统一为【国家, 省/州, 市/县_AI】三级格式（缺失层级用 Unknown 补全）。如果 B 基于 `biomaterial_provider` 推断出了真实地点，无脑信任 B 的推断值。

**必须严格输出如下结构包含 18 个字段的合法 JSON 格式数据（勿带```json标记）：**
{
  "Run": "...", "ReleaseDate": "...", "CollectionDate": "...", "Location": "...", "Source": "...", 
  "Tissue": "...", "Age_GrowthStage": "...", "LibrarySource": "...", "SRAStudy": "...", 
  "BioProject": "...", "ProjectID": "...", "Sample": "...", "BioSample": "...", 
  "ScientificName": "...", "TaxID": "...", "SampleName": "...", "CenterName": "...", "PMID": "..."
}"""

# ----------------- 【GSA 专用增强提示词】 -----------------
PROMPT_GSA_ENHANCE = """你是一个专业的生物信息学数据分析专家。我将提供一条来自 GSA 的原生测序元数据（JSON格式）。
请提纯出最终标准的【10个核心字段】：
1. Run, BioProject, CollectionDate, Source, Age_GrowthStage 等结合语义提炼。
2. 遇到 'not collected', 'missing', 'N/A' 填 'Not_Provided'。
3. Tissue 强制转为全小写单数。
4. ScientificName 综合提炼，去除部位等杂质。
5. LibrarySource 提取极简的测序类型（如果原句是几百字的小作文，强制浓缩为一个单词如 TRANSCRIPTOMIC）。
6. 【Location 推断法则】：必须统一格式为【国家, 省/州, 市/县_AI】（如 'China, Ningxia, Yinchuan_AI'）。
   - 优先查看 Biomaterial provider，如果有，根据提供者推断地点。
   - 如果没有提供者，再根据 CenterName / Organization 推断地点。
   - 实在推断不出省份或城市的，用 Unknown 补全。

【必须严格输出这 10 个键名的合法 JSON 字符串（绝不要使用 ```json 标记）】：
{
  "Run": "...", "CollectionDate": "...", "Location": "...", "Source": "...", "Tissue": "...",
  "Age_GrowthStage": "...", "ScientificName": "...", "LibrarySource": "...", "CenterName": "...", "BioProject": "..."
}"""

# ----------------- 【大一统终极正交化防污染审查提示词】 -----------------
PROMPT_DATA_SANITIZER = """你是一个极其严谨的数据清理程序，你的任务是【去除错位污染】和【浓缩冗余文本】，绝对不是【凭空捏造数据】。
我将给你一条 JSON 数据。请严格按以下规则处理：

1. **全面净化**：遇到 'not collected', 'missing', 'N/A' 等无意义占位符，直接替换为 'Not_Provided'。
2. **ScientificName**：如果里面混入了部位（如 leaf）或来源（如 wild），将其剔除，只保留纯物种名。
3. **Tissue**：强制小写单数！如果原值为空，但 Source 中误填了部位，将部位移到这里。否则保持原样。
4. **Source**：【最高警告】绝不可凭空脑补！如果传入的是 "wild", "cultivated" 等，必须 100% 保持原样！如果被填成了物种名，直接清空，绝对不许瞎编默认值！
5. **Location**：必须强制转为严谨的三级层级结构【国家, 省/州, 市/县_AI】。缺失的层级用 Unknown 补全。
6. **LibrarySource**：【强制浓缩特权】如果传入的是很长的实验操作说明小作文，请强制将其浓缩为一个英文单词（如 TRANSCRIPTOMIC, GENOMIC）。
7. **CollectionDate, Age_GrowthStage, CenterName, BioProject**：必须 100% 保持输入时的原样，一个字都不许改！

【核心纪律】：只能做减法（去污染）、移动（归位）、格式化（转单数/重组地址）和浓缩，绝对禁止做加法凭空捏造数据！
直接输出 9 个键的合法 JSON，勿带 ```json 标记：
{
  "CollectionDate": "...", "Location": "...", "Source": "...", "Tissue": "...",
  "Age_GrowthStage": "...", "ScientificName": "...", "LibrarySource": "...", "CenterName": "...", "BioProject": "..."
}"""

# =========================================================================
# 🛠️ 通用工具与 AI 配置引擎
# =========================================================================

def build_api_kwargs(model, sys_prompt, user_content):
    kwargs = {"model": model, "messages": [{"role": "system", "content": sys_prompt}, {"role": "user", "content": user_content}]}
    if "v4" in model.lower():
        kwargs["stream"] = False
        if "pro" in model.lower():
            kwargs["reasoning_effort"] = "high"
            kwargs["extra_body"] = {"thinking": {"type": "enabled"}}
        else:
            kwargs["temperature"] = 0.1
    else:
        kwargs["temperature"] = 0.1
    return kwargs

def clean_ai_json(text):
    text = re.sub(r'<think>.*?</think>', '', text, flags=re.DOTALL).strip()
    return re.sub(r'^```(?:json)?\s*|\s*```$', '', text, flags=re.IGNORECASE)

def get_retry_session(retries=5, backoff_factor=2.0):
    session = requests.Session()
    retry_strategy = Retry(total=retries, backoff_factor=backoff_factor, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["HEAD", "GET", "POST", "OPTIONS"])
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session

def apply_fill_date(df, fill_date=False, is_gsa=False):
    if not fill_date: return df
    missing_vals = ['nan', 'none', 'null', 'not_provided', 'not provided', 'missing', 'unknown', '', '<na>']
    if 'CollectionDate' not in df.columns: df['CollectionDate'] = pd.NA
    df['CollectionDate'] = df['CollectionDate'].astype(object)
    mask = df['CollectionDate'].astype(str).str.strip().str.lower().isin(missing_vals) | df['CollectionDate'].isna()
    fallback_col = 'LoadDate' if is_gsa and 'LoadDate' in df.columns else ('ReleaseDate' if 'ReleaseDate' in df.columns else None)
    if fallback_col:
        valid_mask = ~df[fallback_col].astype(str).str.strip().str.lower().isin(missing_vals)
        df.loc[mask, 'CollectionDate'] = df.loc[mask, fallback_col].where(valid_mask)
    return df

def save_dual_format(df, base_path):
    df = df.replace(["Not_Provided", "not_provided", "None", ""], pd.NA)
    df = df.loc[:, ~df.columns.duplicated()]
    csv_p = base_path + '.csv'
    tsv_p = base_path + '.tsv'
    with safe_open(csv_p, 'wt') as f:
        f.write('\ufeff')
        df.to_csv(f, index=False, lineterminator='\n')
    with safe_open(tsv_p, 'wt') as f:
        f.write('\ufeff')
        df.to_csv(f, index=False, sep='\t', lineterminator='\n')
    return csv_p

def robust_eutils_get(url, params):
    if os.environ.get("NCBI_API_KEY"): params["api_key"] = os.environ["NCBI_API_KEY"]
    session = get_retry_session()
    try:
        time.sleep(0.2 + secrets.randbelow(201) / 1000)
        res = session.get(url, params=params, headers=HEADERS, timeout=20)
        if res.status_code == 200: return res
    except: pass
    return None

TAXONOMY_CACHE = {}
def resolve_taxonomy(name):
    if pd.isna(name) or str(name).strip().lower() in ['unknown', 'not_provided', 'nan', 'none', '', '<na>']: return pd.NA, name
    name_str = str(name).strip()
    if name_str in TAXONOMY_CACHE: return TAXONOMY_CACHE[name_str]
    s_res = robust_eutils_get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", {"db": "taxonomy", "term": name_str, "retmode": "json"})
    if s_res and (id_list := s_res.json().get('esearchresult', {}).get('idlist', [])):
        taxid = id_list[0]
        f_res = robust_eutils_get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", {"db": "taxonomy", "id": taxid, "retmode": "xml"})
        if f_res and (sci_match := re.search(r'<ScientificName>(.*?)</ScientificName>', f_res.text)):
            TAXONOMY_CACHE[name_str] = (taxid, sci_match.group(1))
            return taxid, sci_match.group(1)
    TAXONOMY_CACHE[name_str] = (pd.NA, name_str)
    return pd.NA, name_str

# =========================================================================
# 🕵️‍♂️ BioProject 溯源引擎 (彻底砸碎连接池，防死锁极速版)
# =========================================================================
def extract_year(date_str: str) -> int:
    m = re.search(r'\b(19\d{2}|20\d{2})\b', str(date_str) or "")
    return int(m.group(1)) if m else 9999

class NCBIAPIClient:
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    def __init__(self, api_key=None):
        self.api_key = api_key
        self.headers = HEADERS.copy()
        self.headers['Connection'] = 'close'

    def _req(self, endpoint, params):
        if self.api_key: params["api_key"] = self.api_key
        params["retmode"] = "json"
        url = self.BASE_URL + '/' + str(endpoint)
        for _ in range(3):
            try:
                time.sleep(0.4)
                r = get_retry_session().get(url, params=params, headers=self.headers, timeout=15)
                if r.status_code == 200: return r.json()
                if r.status_code == 429: time.sleep(2)
            except: time.sleep(1)
        return {}

    def get_bioproject_summary(self, bp_id):
        nid = bp_id.replace("PRJNA", "").replace("PRJEB", "").replace("PRJ", "")
        res = self._req("esummary.fcgi", {"db": "bioproject", "id": nid})
        if "result" in res and nid in res["result"]:
            d = res["result"][nid]
            return {"bioproject_id": bp_id, "title": d.get("project_title", ""), "submission_date": d.get("registration_date", ""), "submitter": d.get("submitter_org", "")}
        return {"bioproject_id": bp_id}

    def get_linked_sra(self, bp_id):
        nid = bp_id.replace("PRJNA", "").replace("PRJEB", "").replace("PRJ", "")
        res = self._req("elink.fcgi", {"dbfrom": "bioproject", "db": "sra", "id": nid})
        try:
            return [l for ls in res.get("linksets", []) for ldb in ls.get("linksetdbs", []) if ldb.get("linkname") == "bioproject_sra" for l in ldb.get("links", [])]
        except: return []

    def get_sra_linked_pubmed(self, sra_ids):
        if not sra_ids: return []
        pmids = []
        for i in range(0, len(sra_ids), 200):
            batch = sra_ids[i:i+200]
            res = self._req("elink.fcgi", {"dbfrom": "sra", "db": "pubmed", "id": ",".join(map(str, batch))})
            try:
                pmids.extend([l for ls in res.get("linksets", []) for ldb in ls.get("linksetdbs", []) if ldb.get("linkname") == "sra_pubmed" for l in ldb.get("links", [])])
            except: pass
        return list(set(pmids))

    def search_exact_id(self, bp_id):
        res = self._req("esearch.fcgi", {"db": "pubmed", "term": bp_id, "retmax": 20})
        return self.get_pubmed_details(res.get("esearchresult", {}).get("idlist", []))

    def get_pubmed_details(self, pmids):
        if not pmids: return []
        articles = []
        for i in range(0, len(pmids), 200):
            batch = pmids[i:i+200]
            res = self._req("esummary.fcgi", {"db": "pubmed", "id": ",".join(map(str, batch))})
            for pmid in batch:
                pmid_s = str(pmid)
                if pmid_s in res.get("result", {}):
                    d = res["result"][pmid_s]
                    doi = next((aid.get("value") for aid in d.get("articleids", []) if aid.get("idtype") == "doi"), "")
                    articles.append({
                        "pmid": pmid_s, "doi": doi, "title": d.get("title", ""),
                        "journal": d.get("fulljournalname", ""), "pubdate": d.get("pubdate", ""),
                        "authors": [a.get("name", "") for a in d.get("authors", [])]
                    })
        return articles

class GSAClientP:
    def __init__(self):
        self.headers = HEADERS.copy()
        self.headers['Connection'] = 'close'

    def get_bioproject_summary(self, bp_id):
        url = "https://ngdc.cncb.ac.cn/gwh/api/public/bioProject/" + str(bp_id)
        for _ in range(3):
            try:
                r = get_retry_session().get(url, headers=self.headers, timeout=15)
                if r.status_code == 200:
                    data = r.json()
                    if data.get("message") == "SUCCESS":
                        return {"bioproject_id": bp_id, "title": data.get("title", ""), "submission_date": data.get("releaseTime", ""), "linked_pmids":[str(pub.get("pubmedId")) for pub in data.get("listPublication", []) if pub.get("pubmedId")]}
                    return {"bioproject_id": bp_id}
            except: time.sleep(1)
        return {"bioproject_id": bp_id}

class EuropePMCClient:
    def __init__(self):
        self.headers = HEADERS.copy()
        self.headers['Connection'] = 'close'

    def search_by_bioproject(self, bp_id):
        for _ in range(2):
            try:
                r = get_retry_session().get("https://www.ebi.ac.uk/europepmc/webservices/rest/search", params={"query": f'"{bp_id}"', "format": "json", "pageSize": 50, "resultType": "core"}, headers=self.headers, timeout=15)
                if r.status_code == 200:
                    return [{"pmid": str(item.get("pmid", "")), "doi": item.get("doi", ""), "title": item.get("title", ""), "journal": item.get("journalTitle", ""), "pubdate": item.get("firstPublicationDate", "")} for item in r.json().get("resultList", {}).get("result",[])]
            except: time.sleep(1)
        return []

def rule_based_primary_expert(articles, bp_info):
    if not articles: return None
    bp_year = extract_year(bp_info.get("submission_date", ""))
    scored = []
    for art in articles:
        if any(kw in art.get("title", "").lower() for kw in ["review", "database", "meta-analysis"]): continue
        diff = extract_year(art.get("pubdate", "")) - bp_year
        if -2 <= diff <= 3: scored.append({"score": 10-abs(diff), "article": art})
    if not scored: return {"pmid": str(min(articles, key=lambda x: extract_year(x.get("pubdate",""))).get("pmid", ""))}
    scored.sort(key=lambda x: -x["score"])
    return {"pmid": str(scored[0]["article"].get("pmid", ""))}

class BioProjectTracer:
    def __init__(self, ncbi_api, openai_client, ai_model, use_scholar=False):
        self.ncbi = NCBIAPIClient(ncbi_api)
        self.gsa = GSAClientP()
        self.epmc = EuropePMCClient()
        self.api_client, self.model = openai_client, ai_model

    def trace(self, bp_id: str) -> dict:
        is_gsa = bp_id.startswith("PRJCA")
        bp_info = self.gsa.get_bioproject_summary(bp_id) if is_gsa else self.ncbi.get_bioproject_summary(bp_id)
        candidates = []
        if is_gsa and bp_info.get("linked_pmids"): 
            candidates.extend(self.ncbi.get_pubmed_details(bp_info["linked_pmids"]))
        if not is_gsa:
            sra_ids = self.ncbi.get_linked_sra(bp_id)
            if sra_ids:
                pmids = self.ncbi.get_sra_linked_pubmed(sra_ids)
                if pmids: candidates.extend(self.ncbi.get_pubmed_details(pmids))
            pubmed_exact = self.ncbi.search_exact_id(bp_id)
            if pubmed_exact: candidates.extend(pubmed_exact)
        try:
            epmc_res = self.epmc.search_by_bioproject(bp_id)
            if epmc_res: candidates.extend(epmc_res)
        except: pass
        
        seen, unique_candidates = set(), []
        for c in candidates:
            uid = str(c.get("pmid")) if c.get("pmid") else c.get("doi")
            if uid and uid not in seen and uid.lower() != "none":
                seen.add(uid)
                unique_candidates.append(c)
        candidates = unique_candidates
        if not candidates: return {}

        if self.api_client:
            prompt = """识别该BioProject首发文献。严格返回JSON格式 {"1": {"pmid": "...", "doi": "...", "category": "PRIMARY"}}。"""
            usr_content = f"项目编号: {bp_id} | 提交日期: {bp_info.get('submission_date', 'Unknown')}\n"
            for i, c in enumerate(candidates, 1):
                usr_content += f"[{i}] PMID:{c.get('pmid')} | DOI:{c.get('doi')} | Title:{c.get('title')}\n"
            try:
                kwargs = build_api_kwargs(self.model, prompt, usr_content)
                kwargs['timeout'] = 30
                res = self.api_client.chat.completions.create(**kwargs)
                parsed = json.loads(clean_ai_json(res.choices[0].message.content))
                for v in parsed.values():
                    if v.get("category") == "PRIMARY":
                        for cand in candidates:
                            if (v.get('pmid') and str(cand.get('pmid')) == str(v['pmid'])) or (v.get('doi') and cand.get('doi') == v['doi']):
                                return cand
            except: pass
        
        r_cls = rule_based_primary_expert(candidates, bp_info)
        if r_cls and r_cls.get("pmid"):
            for cand in candidates:
                if str(cand.get('pmid')) == str(r_cls.get('pmid')): return cand
        return candidates[0] if candidates else {}

# =========================================================================
# 🟢 SRA 高并发引擎
# =========================================================================
class SRAPipeline:
    def __init__(self, mode, api_client, ai_model, out_dir, fill_date):
        self.mode, self.api_client, self.ai_model, self.fill_date = mode, api_client, ai_model, fill_date
        self.base_dir = os.path.join(out_dir, "SRA_Results")
        self.d_xml, self.d_full, self.d_local, self.d_api, self.d_arb = [os.path.join(self.base_dir, x) for x in ["1_raw_xml", "2_full_json", "3_1_local_parsed", "3_2_api_inferred", "3_3_ai_arbitrated"]]
        for d in [self.d_xml, self.d_full, self.d_local, self.d_api, self.d_arb]: os.makedirs(d, exist_ok=True)

    def download_xml(self, srr):
        xp = os.path.join(self.d_xml, srr + '.xml')
        if os.path.exists(xp): return open(xp, 'r', encoding='utf-8').read(), "[✓XML缓存]"
        try:
            # E-utilities 直接拉取 EXPERIMENT_PACKAGE XML（Windows 无 EDirect）
            res = robust_eutils_get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                                    {"db": "sra", "id": srr, "retmode": "xml"})
            if res is not None and "<EXPERIMENT_PACKAGE" in res.text:
                with safe_open(xp, 'wt') as f: f.write(res.text)
                return res.text, "[✓XML下载]"
        except: return None, "[❌XML下载失败]"
        return None, "[❌XML无效]"

    def xml_to_json(self, srr, xml):
        jp = os.path.join(self.d_full, srr + '_full.json')
        if os.path.exists(jp): return json.load(open(jp, 'r', encoding='utf-8')), "[✓JSON缓存]"
        try:
            fl = ('EXPERIMENT_PACKAGE', 'SAMPLE_ATTRIBUTE', 'EXPERIMENT_ATTRIBUTE', 'RUN', 'SRAFile', 'Alternatives', 'EXTERNAL_ID', 'IDENTIFIERS', 'STUDY_LINK')
            d = xmltodict.parse(xml, process_namespaces=False, force_list=fl, dict_constructor=dict)
            with safe_open(jp, 'wt') as f: json.dump(d, f, indent=4, ensure_ascii=False)
            return d, "[✓JSON转换]"
        except: return None, "[❌JSON转换失败]"

    def local_parse(self, srr, data_dict):
        jp = os.path.join(self.d_local, srr + '_local.json')
        if os.path.exists(jp): return json.load(open(jp, 'r', encoding='utf-8')), "[✓Local缓存]"
        
        def get_org_loc(d):
            if isinstance(d, dict):
                if 'Organization' in d and isinstance(d['Organization'], dict):
                    addr = d['Organization'].get('Address', {})
                    if isinstance(addr, list) and addr: addr = addr[0]
                    if isinstance(addr, dict) and any(k in addr for k in ['Country', 'Sub', 'City']):
                        c, s, ct = addr.get('Country', ''), addr.get('Sub', ''), addr.get('City', '')
                        parts = [str(c).strip()] if c else []
                        sub_ct = " ".join([str(x).strip() for x in [s, ct] if x])
                        if sub_ct: parts.append(sub_ct)
                        return ", ".join(parts) if parts else None
                for v in d.values():
                    if res := get_org_loc(v): return res
            elif isinstance(d, list):
                for item in d:
                    if res := get_org_loc(item): return res
            return None
        
        global_org_loc = get_org_loc(data_dict)
        pkgs = data_dict.get('EXPERIMENT_PACKAGE_SET', data_dict).get('EXPERIMENT_PACKAGE', [])
        if isinstance(pkgs, dict): pkgs = [pkgs]
        records = []
        try:
            for pkg in pkgs:
                rec = {"query_id": srr, "CenterName": pkg.get('SUBMISSION', {}).get('@center_name', "Not_Provided"), "PMID": "Not_Provided"}
                st = pkg.get('STUDY', {})
                rec["SRAStudy"], rec["Study_title"], rec["Study_abstract"] = st.get('@accession', "Not_Provided"), st.get('DESCRIPTOR', {}).get('STUDY_TITLE', "Not_Provided"), st.get('DESCRIPTOR', {}).get('STUDY_ABSTRACT', "Not_Provided")
                
                links = st.get('STUDY_LINKS', {}).get('STUDY_LINK', [])
                if isinstance(links, dict): links = [links]
                for slink in links:
                    if isinstance(xref := slink.get('XREF_LINK', {}), dict) and xref.get('DB', '').lower() == 'pubmed':
                        rec["PMID"] = xref.get('ID', "Not_Provided")

                for ib in st.get('IDENTIFIERS', []):
                    for eid in ib.get('EXTERNAL_ID', []):
                        if eid.get('@namespace') == 'BioProject': rec["BioProject"] = eid.get('#text', "Not_Provided")
                
                exp, samp = pkg.get('EXPERIMENT', {}), pkg.get('SAMPLE', {})
                rec["LibrarySource"] = exp.get('DESIGN', {}).get('LIBRARY_DESCRIPTOR', {}).get('LIBRARY_SOURCE', "Not_Provided")
                rec["Sample"], rec["ScientificName"], rec["SampleName"] = samp.get('@accession', "Not_Provided"), samp.get('SAMPLE_NAME', {}).get('SCIENTIFIC_NAME', "Not_Provided"), samp.get('@alias', "Not_Provided")
                rec["TaxID"] = samp.get('SAMPLE_NAME', {}).get('TAXON_ID', "Not_Provided")

                for attr in samp.get('SAMPLE_ATTRIBUTES', {}).get('SAMPLE_ATTRIBUTE', []):
                    if attr.get('TAG'): 
                        val = str(attr.get('VALUE', "Not_Provided")).strip()
                        if val.lower() in ['not collected', 'missing', 'n/a', 'none', 'unknown']: val = "Not_Provided"
                        rec[attr.get('TAG').lower()] = val

                geo_loc = rec.get('geo_loc_name', "Not_Provided")
                if geo_loc != "Not_Provided":
                    rec['Location'] = geo_loc
                elif global_org_loc:
                    rec['Location'] = global_org_loc
                else:
                    rec['Location'] = "Not_Provided"

                for run in pkg.get('RUN_SET', {}).get('RUN', []):
                    rr = rec.copy()
                    rr["Run"] = run.get('@accession', "Not_Provided")
                    file_dates = []
                    for sf in run.get('SRAFiles', {}).get('SRAFile', []):
                        if c := sf.get('@cluster', '').lower().strip():
                            for k, v in sf.items():
                                if k.startswith('@') and k != '@cluster': rr[f"{c}_{k.replace('@', '')}"] = v
                                if k == '@date': file_dates.append(v)
                            for alt in sf.get('Alternatives', []):
                                for k, v in alt.items():
                                    if k == '@date': file_dates.append(v)
                    if run.get('@published'): file_dates.append(run.get('@published'))
                    
                    valid_dates = sorted([str(d).split(' ')[0] for d in file_dates if isinstance(d, str) and len(str(d).strip()) >= 4 and str(d).lower() != "true"])
                    rr["ReleaseDate"] = valid_dates[0] if valid_dates else "Not_Provided"
                    rr["LoadDate"] = "true" 
                    
                    c_date = rr.get('collection_date', rr.get('collection date', 'Not_Provided'))
                    if str(c_date).lower() in ['not_provided', '', 'nan', 'none', '<na>', 'unknown', 'not collected', 'missing']:
                        if valid_dates: rr['collection_date'] = valid_dates[0]
                    records.append(rr)
            with safe_open(jp, 'wt') as f: json.dump(records, f, indent=4, ensure_ascii=False)
            return records, "[✓Local解析]"
        except: return None, "[❌Local解析失败]"

    def api_infer(self, srr, xml):
        jp = os.path.join(self.d_api, srr + '_api.json')
        if os.path.exists(jp): return json.load(open(jp, 'r', encoding='utf-8')), "[✓API缓存]"
        try:
            kwargs = build_api_kwargs(self.ai_model, PROMPT_SRA_INFER, f"【XML 全文】:\n{xml[:20000]}")
            res = self.api_client.chat.completions.create(**kwargs)
            rd = json.loads(clean_ai_json(res.choices[0].message.content))
            rd['query_id'] = srr
            with safe_open(jp, 'wt') as f: json.dump(rd, f, indent=4, ensure_ascii=False)
            return rd, "[✓API推断]"
        except: return {"query_id": srr}, "[❌API推断失败]"

    def api_arbitrate(self, srr, ld, ad):
        jp = os.path.join(self.d_arb, srr + '_arb.json')
        if os.path.exists(jp): return json.load(open(jp, 'r', encoding='utf-8')), "[✓Arb缓存]"
        sl = {k:v for k,v in ld[0].items() if "url" not in k and "md5" not in k and "size" not in k and "abstract" not in k and "title" not in k}
        try:
            kwargs = build_api_kwargs(self.ai_model, PROMPT_SRA_ARBITRATE, f"【数据A】:\n{json.dumps(sl)}\n\n【数据B】:\n{json.dumps(ad)}")
            res = self.api_client.chat.completions.create(**kwargs)
            fd = json.loads(clean_ai_json(res.choices[0].message.content))
            fd['query_id'] = srr
            with safe_open(jp, 'wt') as f: json.dump(fd, f, indent=4, ensure_ascii=False)
            return fd, "[✓Arb仲裁]"
        except: return ad, "[❌Arb仲裁失败]"

    def process_all(self, sra_list, threads):
        all_local, all_arb = [], []
        def _task(srr):
            logs = []
            xml, l1 = self.download_xml(srr); logs.append(l1)
            if not xml: return {"logs": "".join(logs)}
            d_dict, l2 = self.xml_to_json(srr, xml); logs.append(l2)
            if not d_dict: return {"logs": "".join(logs)}
            l_res, l31 = self.local_parse(srr, d_dict); logs.append(l31)
            if self.mode in ['api', 'both'] and l_res:
                a_res, l32 = self.api_infer(srr, xml); logs.append(l32)
                if a_res:
                    arb_res, l33 = self.api_arbitrate(srr, l_res, a_res); logs.append(l33)
                    return {"local": l_res, "arb": arb_res, "logs": "".join(logs)}
            return {"local": l_res, "logs": "".join(logs)}

        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            fut2srr = {executor.submit(_task, s): s for s in sra_list}
            for future in tqdm(concurrent.futures.as_completed(fut2srr), total=len(sra_list), desc="🟢 SRA Processing"):
                srr = fut2srr[future]
                try:
                    res = future.result()
                    tqdm.write(f"  [SRA] {srr.ljust(12)} -> {res.get('logs', '')}")
                    if res.get("local"): all_local.extend(res["local"])
                    if res.get("arb"): all_arb.append(res["arb"])
                except Exception as e: tqdm.write(f"  ❌ [SRA] {srr} 崩溃: {e}")

        if not all_arb and not all_local: return pd.DataFrame()
        df_local = pd.DataFrame(all_local)
        df_arb = pd.DataFrame(all_arb).drop(columns=['Run', 'run'], errors='ignore') if all_arb else pd.DataFrame()
        tech_cols = ['query_id', 'Run', 'PMID', 'TaxID'] + [c for c in df_local.columns if 'url' in c or 'md5' in c or 'size' in c] + DESC_COLS
        df_tech = df_local[list(set(tech_cols).intersection(df_local.columns))]
        df_merged = pd.merge(df_tech, df_arb, on='query_id', how='left') if not df_arb.empty else df_local
        
        final_cols =[c for c in CORE_18_SRA if c in df_merged.columns] +[c for c in DESC_COLS if c in df_merged.columns]
        df_merged = df_merged[final_cols]

        # SRA 漏掉的时间填补函数，在此处补上
        df_merged = apply_fill_date(df_merged, self.fill_date, is_gsa=False)

        save_dual_format(df_merged, os.path.join(self.base_dir, "SRA_Ultimate_Merged"))
        return df_merged


class GSAPipeline:
    def __init__(self, mode, api_client, ai_model, out_dir, fill_date):
        self.mode, self.api_client, self.ai_model, self.fill_date = mode, api_client, ai_model, fill_date
        self.base_dir = os.path.join(out_dir, "GSA_Results")
        # 建立物理存储矩阵
        self.d_web = os.path.join(self.base_dir, "0_web_cache")
        self.d_csv = os.path.join(self.base_dir, "1_csv")
        self.d_xlsx = os.path.join(self.base_dir, "2_xlsx")
        self.d_api = os.path.join(self.base_dir, "3_api_enhanced")
        for d in [self.d_web, self.d_csv, self.d_xlsx, self.d_api]:
            os.makedirs(d, exist_ok=True)
        self.headers = HEADERS.copy()
        self.headers['Connection'] = 'close'
        # 实例级缓存，用于跨 Accession 保持机构地理位置一致性
        self.org_loc_cache: dict = {}

    def sanitize_for_json(self, data_dict):
        """安全清理字典中的 pd.NA/NaN/NaT、Numpy 类型以及 Timestamp 时间格式，彻底杜绝 TypeError 和 ValueError"""
        clean_dict = {}
        for k, v in data_dict.items():
            # 1. 优先拦截复合数据类型 (List, Tuple, Set, Numpy Array)
            if isinstance(v, (list, tuple, set)) or hasattr(v, 'tolist'):
                v_iterable = v.tolist() if hasattr(v, 'tolist') else list(v)
                clean_list = []
                for x in v_iterable:
                    if x is pd.NA or x is None:
                        clean_list.append(None)
                    elif isinstance(x, (float, pd.Timestamp, datetime.datetime)) and pd.isna(x):
                        clean_list.append(None)
                    elif isinstance(x, (pd.Timestamp, datetime.date, datetime.datetime)):
                        clean_list.append(str(x))
                    else:
                        clean_list.append(x)
                clean_dict[k] = clean_list

            # 2. 拦截嵌套字典，递归清洗
            elif isinstance(v, dict):
                clean_dict[k] = self.sanitize_for_json(v)

            # 3. 处理绝对空值 (pd.NA, None)
            elif v is pd.NA or v is None:
                clean_dict[k] = None

            # 4. 处理标量缺失值 (NaN, NaT)
            elif isinstance(v, (float, pd.Timestamp, datetime.datetime)) and pd.isna(v):
                clean_dict[k] = None

            # 5. 处理正常的时间对象强转
            elif isinstance(v, (pd.Timestamp, datetime.date, datetime.datetime)):
                clean_dict[k] = str(v)

            # 6. 处理 Numpy 标量对象 (如 np.int64, np.float32)
            elif hasattr(v, 'item') and callable(getattr(v, 'item')):
                try:
                    val = v.item()
                    clean_dict[k] = str(val) if isinstance(val, (datetime.date, datetime.datetime)) else val
                except Exception:
                    clean_dict[k] = str(v)

            # 7. 终极兜底机制验证
            else:
                try:
                    json.dumps(v)
                    clean_dict[k] = v
                except TypeError:
                    clean_dict[k] = str(v)

        return clean_dict

    def parse_all_from_html(self, html, acc):
        """[一站式解析] 深度提取网页所有元数据，锚定表头防止 example 干扰"""
        feat = {
            "SubmissionDate": pd.NA, "ReleaseDate": pd.NA, "Organization": pd.NA,
            "CRA": pd.NA, "CRX": pd.NA, "PRJ": pd.NA, "SAMC": pd.NA,
            "Title": pd.NA, "TaxID": pd.NA, "ScientificName": pd.NA, "Platform": pd.NA,
            "FileNames": [], "FileSizes_Bytes": [],
            "LibraryStrategy": "Not_Provided", "LibrarySource": "Not_Provided",
            "LibrarySelection": "Not_Provided", "LibraryLayout": "Not_Provided",
            "LibraryConstruction": ""
        }
        if not html: return feat
        h = html.replace('\r', '').replace('\n', ' ')

        # 1. 基础日期与单位提取
        m_date = re.search(r'<th[^>]*>(?:提交日期|Submission Date|Submit Date)</th>\s*<td[^>]*>(.*?)</td>', h, re.I)
        if m_date: feat["SubmissionDate"] = re.sub(r'<[^>]+>', '', m_date.group(1)).strip()
        m_rel = re.search(r'<th[^>]*>(?:发布日期|Release Date)</th>\s*<td[^>]*>(.*?)</td>', h, re.I)
        if m_rel: feat["ReleaseDate"] = re.sub(r'<[^>]+>', '', m_rel.group(1)).strip()
        m_org = re.search(r'<th[^>]*>(?:所属单位|Organization)</th>\s*<td[^>]*>(.*?)</td>', h, re.I)
        if m_org: feat["Organization"] = re.sub(r'<[^>]+>', '', m_org.group(1)).strip()

        # 2. 精确 ID 提取
        m_cra = re.search(r'<th[^>]*>(?:GSA编号|GSA accession)</th>\s*<td[^>]*>.*?>(CRA\d+)</a>', h, re.I)
        if m_cra: feat["CRA"] = m_cra.group(1)
        m_crx = re.search(r'<th[^>]*>(?:实验编号|Experiment accession)</th>\s*<td[^>]*>.*?>(CRX\d+)</a>', h, re.I)
        if m_crx: feat["CRX"] = m_crx.group(1)
        m_prj = re.search(r'<th[^>]*>(?:项目编号|BioProject accession)</th>\s*<td[^>]*>.*?>(PRJ[A-Z0-9]+)</a>', h, re.I)
        if m_prj: feat["PRJ"] = m_prj.group(1)
        m_sam = re.search(r'<th[^>]*>(?:样本编号|BioSample accession)</th>\s*<td[^>]*>.*?>(SAM[A-Z0-9]+)</a>', h, re.I)
        if m_sam: feat["SAMC"] = m_sam.group(1)

        m_tit = re.search(r'<th[^>]*>(?:标题|Title)</th>\s*<td[^>]*>(.*?)</td>', h, re.I)
        if m_tit: feat["Title"] = re.sub(r'<[^>]+>', '', m_tit.group(1)).strip()
        m_tax = re.search(r'wwwtax\.cgi\?id=(\d+)', h, re.I)
        if m_tax: feat["TaxID"] = m_tax.group(1)
        m_sci = re.search(r'Taxonomy/Browser/wwwtax\.cgi\?id=\d+"[^>]*>([^<]+)</a>', h, re.I)
        if m_sci: feat["ScientificName"] = m_sci.group(1)
        m_plat = re.search(r'<th[^>]*>(?:测序平台|Platform)</th>\s*<td[^>]*>\s*(.*?)\s*</td>', h, re.I)
        if m_plat: feat["Platform"] = re.sub(r'<[^>]+>', '', m_plat.group(1)).strip()

        # 3. 文库策略与"小作文"
        lib_m = re.search(r'<th[^>]*>(?:建库信息|Library Information)</th>.*?<th[^>]*>(?:文库布局|Library Layout)</th>.*?</tr>\s*<tr>\s*<td[^>]*>(.*?)</td>\s*<td[^>]*>(.*?)</td>\s*<td[^>]*>\s*(.*?)\s*</td>\s*<td[^>]*>\s*(.*?)\s*</td>\s*<td[^>]*>\s*(.*?)\s*</td>\s*<td[^>]*>\s*(.*?)\s*</td>', h, re.I)
        if lib_m:
            feat["LibraryConstruction"] = re.sub(r'<[^>]+>', '', lib_m.group(2)).strip()
            feat["LibraryStrategy"], feat["LibrarySource"], feat["LibrarySelection"], feat["LibraryLayout"] = [x.strip() for x in lib_m.groups()[2:]]

        # 4. 文件路径与字节反演
        file_m = re.search(rf'<a[^>]*>{acc}</a>.*?</td>\s*<td[^>]*>(.*?)</td>\s*<td[^>]*>(.*?)</td>', h, re.I)
        if file_m:
            feat["FileNames"] = re.findall(rf'({acc}[^<]*\.(?:fq|fastq|bam|sra)(?:\.gz)?)', file_m.group(1))
            sizes_str = re.findall(r'([\d,]+\.\d+|\d+)', re.sub(r'<br\s*/?>', ' ', file_m.group(2)))
            feat["FileSizes_Bytes"] = [str(int(float(s.replace(',', '')) * 1024 * 1024)) for s in sizes_str]
        return feat

    def build_fallback_df(self, acc, wf):
        """[镜像重建] 1:1 还原 GSA 标准 CSV 格式"""
        cra_id = wf.get("CRA") if pd.notna(wf.get("CRA")) else "UNKNOWN_CRA"
        f_names = wf.get("FileNames", [])
        f_sizes = wf.get("FileSizes_Bytes", [])
        d_paths = [f"ftp://download.big.ac.cn/gsa/{cra_id}/{acc}/{fn}" for fn in f_names]
        return pd.DataFrame([{
            "Run": acc, "Run_Accession": acc, "query_id": acc, "Center": "NGDC",
            "ReleaseDate": wf.get("ReleaseDate", pd.NA), "LoadDate": wf.get("SubmissionDate", pd.NA),
            "FileName": "|".join(f_names) if f_names else pd.NA,
            "FileSize": "|".join(f_sizes) if f_sizes else pd.NA,
            "Download_path": "|".join(d_paths) if d_paths else pd.NA,
            "Experiment": wf.get("CRX", pd.NA), "Title": wf.get("Title", pd.NA),
            "LibraryStrategy": wf.get("LibraryStrategy"), "LibrarySource": wf.get("LibrarySource"),
            "LibrarySelection": wf.get("LibrarySelection"), "LibraryLayout": wf.get("LibraryLayout"),
            "Platform": wf.get("Platform", pd.NA), "BioProject": wf.get("PRJ", pd.NA),
            "BioSample": wf.get("SAMC", pd.NA), "TaxID": wf.get("TaxID", pd.NA),
            "ScientificName": wf.get("ScientificName", pd.NA), "Submission": cra_id,
            "Organization": wf.get("Organization", pd.NA)
        }])

    def process_all(self, gsa_list):
        all_dfs = []
        http = get_retry_session()
        for acc in tqdm(gsa_list, desc="🔵 GSA Processing"):
            logs, df_main, web_feat = [], pd.DataFrame(), {}
            csv_f = os.path.join(self.d_csv, acc + '.csv')
            web_cache_f = os.path.join(self.d_web, acc + '.json')

            try:
                # --- 步骤 0: 网页特征缓存 ---
                if os.path.exists(web_cache_f):
                    with open(web_cache_f, 'r', encoding='utf-8') as f: web_feat = json.load(f)
                    logs.append("[✓网页缓存]")
                else:
                    res_web = http.get("https://ngdc.cncb.ac.cn/gsa/search?searchTerm=" + str(acc), headers=self.headers, timeout=15)
                    web_feat = self.parse_all_from_html(res_web.text, acc)
                    with safe_open(web_cache_f, 'wt') as f:
                        json.dump(self.sanitize_for_json(web_feat), f, indent=4, ensure_ascii=False)
                    logs.append("[✓网页提取]")

                # --- 步骤 1: CSV 下载/重建 ---
                api_success = False
                if os.path.exists(csv_f) and os.path.getsize(csv_f) > 500:
                    df_main = pd.read_csv(csv_f); logs.append("[✓CSV缓存]"); api_success = True
                else:
                    sterm = web_feat.get("CRA") if pd.notna(web_feat.get("CRA")) else acc
                    api_path = "getRunInfoByCra" if "CRA" in str(sterm) else "getRunInfo"
                    try:
                        api_url = "https://ngdc.cncb.ac.cn/gsa/search/" + str(api_path)
                        res = http.post(api_url, data={"searchTerm": sterm, "totalDatas": 9999, "downLoadCount": 9999}, headers=self.headers, timeout=20)
                        text = res.content.decode('utf-8', errors='ignore')
                        if "<html" in text[:100].lower() or len(text.strip().splitlines()) <= 1: raise ValueError("API_FAIL")
                        if acc.startswith(('CRR', 'CRX')):
                            lines = text.splitlines(); text = "\n".join([lines[0]] + [l for l in lines[1:] if acc in l]) + "\n"
                        with safe_open(csv_f, 'wt') as f: f.write(text)
                        df_main = pd.read_csv(csv_f); logs.append("[✓CSV下载]"); api_success = True
                    except: pass

                if not api_success or df_main.empty:
                    logs.append("[⚠️容灾重建]"); df_main = self.build_fallback_df(acc, web_feat)
                    with safe_open(csv_f, 'wt') as f:
                        df_main.to_csv(f, index=False, lineterminator='\n')

                # --- 步骤 2: XLSX 三表级联融合 ---
                cra = web_feat.get("CRA")
                final_df = df_main.copy()
                if pd.notna(cra) and cra not in ["Not_Provided", None, ""]:
                    xlsx_f = os.path.join(self.d_xlsx, str(cra) + '.xlsx')

                    if not os.path.exists(xlsx_f):
                        try:
                            rx = http.post("https://ngdc.cncb.ac.cn/gsa/file/exportExcelFile", data={"type": 3, "dlAcession": cra}, headers=self.headers, timeout=30)
                            if len(rx.content) > 2000:
                                with safe_open(xlsx_f, 'wb') as f: f.write(rx.content)
                                logs.append("[✓XLS下载]")
                        except Exception as e: pass

                    if os.path.exists(xlsx_f):
                        try:
                            xls = pd.ExcelFile(xlsx_f, engine='openpyxl')
                            df_run_s = pd.read_excel(xls, sheet_name='Run')
                            df_exp_s = pd.read_excel(xls, sheet_name='Experiment')
                            df_samp_s = pd.read_excel(xls, sheet_name='Sample')

                            for d in [df_run_s, df_exp_s, df_samp_s]:
                                d.columns = d.columns.str.strip()
                            df_run_s = df_run_s.loc[:, ~df_run_s.columns.duplicated()]
                            df_exp_s = df_exp_s.loc[:, ~df_exp_s.columns.duplicated()]
                            df_samp_s = df_samp_s.loc[:, ~df_samp_s.columns.duplicated()]

                            def get_c(df, target):
                                for c in df.columns:
                                    if c.lower() == target.lower(): return c
                                return target

                            run_exp_col = get_c(df_run_s, 'Experiment accession')
                            run_acc_col = get_c(df_run_s, 'Accession')
                            exp_acc_col = get_c(df_exp_s, 'Accession')
                            exp_sam_col = get_c(df_exp_s, 'BioSample accession')
                            sam_acc_col = get_c(df_samp_s, 'Accession')

                            m1 = pd.merge(df_run_s, df_exp_s, left_on=run_exp_col, right_on=exp_acc_col, suffixes=('', '_Exp'))
                            m2 = pd.merge(m1, df_samp_s, left_on=exp_sam_col, right_on=sam_acc_col, suffixes=('', '_Sam'))

                            acc_merged = m2[m2[run_acc_col] == acc].copy()
                            if not acc_merged.empty:
                                merged_name = acc + '_merged.csv'
                                merged_f = os.path.join(self.d_xlsx, merged_name)
                                with safe_open(merged_f, 'wt') as f:
                                    acc_merged.to_csv(f, index=False, lineterminator='\n')

                                df_main_renamed = df_main.rename(columns={'Run accession': 'Run_Accession', 'Run': 'Run_Accession'})
                                df_main_renamed = df_main_renamed.loc[:, ~df_main_renamed.columns.duplicated()]
                                acc_merged = acc_merged.loc[:, ~acc_merged.columns.duplicated()]

                                final_df = pd.merge(df_main_renamed, acc_merged, left_on='Run_Accession', right_on=run_acc_col, how='left', suffixes=('', '_XLS'))

                                # 【防御性加固】：防止 clean_v 内部触发 ValueError
                                invalid_vals = ['not collected', 'missing', 'none', 'na', 'nan', 'unknown']
                                def clean_v(v):
                                    # 1. 拦截数组对象，防止 pd.isna 崩溃
                                    if isinstance(v, (list, tuple, set)) or hasattr(v, 'tolist'):
                                        return str(v)
                                    # 2. 安全判断标量缺失值
                                    if isinstance(v, (float, pd.Timestamp, datetime.datetime)) and pd.isna(v):
                                        return pd.NA
                                    if v is pd.NA or v is None:
                                        return pd.NA

                                    s = str(v).strip()
                                    if s.lower() in invalid_vals or s == "": return pd.NA
                                    return s

                                cul_col = get_c(final_df, 'Cultivar')
                                if cul_col in final_df.columns:
                                    final_df['Source'] = final_df[cul_col].apply(clean_v)

                                tis_col = get_c(final_df, 'Tissue')
                                if tis_col in final_df.columns:
                                    final_df['Tissue'] = final_df[tis_col].apply(clean_v)

                                age_col = get_c(final_df, 'Age')
                                unit_col = get_c(final_df, 'Age unit')
                                dev_col = get_c(final_df, 'Dev stage')

                                age_vals = []
                                for _, r in final_df.iterrows():
                                    a = clean_v(r.get(age_col, pd.NA))
                                    u = clean_v(r.get(unit_col, pd.NA))
                                    dev = clean_v(r.get(dev_col, pd.NA))
                                    res = []
                                    if pd.notna(a): res.append(f"{a} {u}".strip() if pd.notna(u) else a)
                                    if pd.notna(dev): res.append(dev)
                                    age_vals.append(" | ".join(res) if res else pd.NA)
                                final_df['Age_GrowthStage'] = age_vals

                                col_date = get_c(final_df, 'Collection date')
                                if col_date in final_df.columns:
                                    final_df['CollectionDate'] = final_df[col_date]

                                logs.append("[✓XLS级联融合]")
                        except Exception as e:
                            logs.append(f"[⚠️XLS融合失败:{str(e)[:30]}]")

                # --- 步骤 3: AI 深度推断 ---
                if self.mode in ['api', 'both'] and self.api_client:
                    for idx, row in final_df.iterrows():
                        q_id = row.get('Run_Accession', acc)
                        org = row.get('Organization', web_feat.get("Organization"))
                        api_json_name = str(q_id) + '_api.json'
                        api_json_f = os.path.join(self.d_api, api_json_name)
                        ai_res = {}

                        if os.path.exists(api_json_f):
                            try:
                                with open(api_json_f, 'r', encoding='utf-8') as f: ai_res = json.load(f)
                                logs.append("[✓API缓存]")
                                if pd.notna(org) and ai_res.get('Location') and "Unknown" not in ai_res['Location']:
                                    self.org_loc_cache[org] = ai_res['Location']
                            except: pass

                        if not ai_res:
                            raw_context = row.dropna().to_dict()
                            raw_context["Tissue"] = row.get('Tissue_Standard', pd.NA)
                            raw_context["Age_GrowthStage"] = row.get('Age_GrowthStage_Standard', pd.NA)
                            raw_context["Source"] = row.get('Source_Standard', pd.NA)
                            raw_context["Library_Construction_Detail"] = web_feat.get("LibraryConstruction", "")

                            safe_context = self.sanitize_for_json(raw_context)

                            try:
                                kwargs = build_api_kwargs(self.ai_model, PROMPT_GSA_ENHANCE, json.dumps(safe_context, ensure_ascii=False))
                                res = self.api_client.chat.completions.create(**kwargs)
                                ai_res = json.loads(clean_ai_json(res.choices[0].message.content))

                                if pd.notna(org) and ai_res.get('Location') and "Unknown" not in ai_res['Location']:
                                    self.org_loc_cache[org] = ai_res['Location']

                                with safe_open(api_json_f, 'wt') as f:
                                    json.dump(self.sanitize_for_json(ai_res), f, indent=4, ensure_ascii=False)
                                logs.append("[✓API推断]")
                            except Exception as e:
                                error_msg = str(e).replace('\n', ' ')
                                logs.append(f"[⚠️API报错:{type(e).__name__} | {error_msg[:50]}]")
                                ai_res = {}

                        if pd.notna(org) and org in self.org_loc_cache:
                            if ai_res.get('Location') != self.org_loc_cache[org]:
                                ai_res['Location'] = self.org_loc_cache[org]
                                logs.append("[✓位置同步]")

                        for col in GSA_AI_10:
                            val = ai_res.get(col)
                            if pd.notna(val) and str(val).strip().lower() not in ["not_provided", "unknown", "none", "nan", ""]:
                                final_df.at[idx, col] = val

                # --- 步骤 4: 最终字段严格 28 列对齐 ---
                final_df['query_id'] = acc
                final_df['Run_Accession'] = final_df.get('Run_Accession', acc)
                final_df['Run'] = final_df['Run_Accession']
                final_df['Center'] = "NGDC"
                final_df['CenterName'] = web_feat.get("Organization", "NGDC")
                final_df['Organization'] = web_feat.get("Organization", pd.NA)
                final_df['ReleaseDate'] = web_feat.get("ReleaseDate", pd.NA)
                final_df['LoadDate'] = web_feat.get("SubmissionDate", pd.NA)
                final_df['Submission'] = cra

                final_df['Experiment'] = final_df.get('Experiment accession', final_df.get('Experiment', pd.NA))
                final_df['BioProject'] = final_df.get('BioProject accession', final_df.get('BioProject', pd.NA))
                final_df['BioSample'] = final_df.get('BioSample accession', final_df.get('BioSample', pd.NA))

                target_cols = [
                    'Run', 'Run_Accession', 'query_id', 'Center', 'ReleaseDate', 'LoadDate',
                    'FileName', 'FileSize', 'Download_path', 'Experiment', 'Title',
                    'LibraryStrategy', 'LibrarySource', 'LibrarySelection', 'LibraryLayout',
                    'Platform', 'BioProject', 'BioSample', 'TaxID', 'ScientificName',
                    'Submission', 'Organization', 'Location', 'CenterName', 'CollectionDate',
                    'Tissue', 'Source', 'Age_GrowthStage'
                ]
                for c in target_cols:
                    if c not in final_df.columns: final_df[c] = pd.NA

                all_dfs.append(final_df[target_cols])
                tqdm.write(f"  [GSA] {acc.ljust(12)} -> {''.join(dict.fromkeys(logs))}")

            except Exception as e:
                tqdm.write(f"  [GSA] {acc.ljust(12)} -> [❌失败: {str(e)}]")

        if not all_dfs: return pd.DataFrame()
        df_all = pd.concat(all_dfs, ignore_index=True)
        df_all = apply_fill_date(df_all, self.fill_date, is_gsa=True)
        save_dual_format(df_all, os.path.join(self.base_dir, "GSA_Ultimate_Merged"))
        return df_all

def get_lists_by_regex(input_list):
    sra, gsa, unknown = [], [], []
    s_pat = re.compile(r'^[EDS]R[PRXS]\d+|^PRJ[END]\d+|^SAM[END]\d+|^GS[EM]\d+')
    g_pat = re.compile(r'^CR[PRX]\d+|^CRA\d+|^PRJCA?\d+|^SAMC[A-Z]?\d+')
    for acc in input_list:
        acc = acc.strip().upper()
        if s_pat.match(acc): sra.append(acc)
        elif g_pat.match(acc): gsa.append(acc)
        else: unknown.append(acc)
    return sra, gsa, unknown

def generate_global_datavzrd(file_path, out_dir):
    """
    通用报告生成器：自动识别 CSV/TSV，精确匹配列名，带完整错误输出。
    datavzrd 未安装时干净跳过（不影响元数据输出）。
    """
    import shutil
    if not shutil.which('datavzrd'):
        print("ℹ️ 未安装 datavzrd，跳过交互式报告生成（不影响元数据输出）")
        return None

    abs_path = os.path.abspath(file_path).replace('\\', '/')
    file_name = os.path.basename(file_path)

    # 自动判定分隔符
    sep = "\\t" if file_path.endswith('.tsv') else ","

    yaml_name = f"Config_{file_name}.yaml"
    yaml_path = os.path.join(out_dir, yaml_name)
    report_dir = os.path.join(out_dir, f"Report_{file_name.split('.')[0]}")

    # 注意：Python f-string 中 {{value}} 会被转义为 {value}，这正是 datavzrd 需要的格式
    yaml_content = f"""name: "Metadata_Report"
datasets:
  main_table:
    path: "{abs_path}"
    separator: "{sep}"
webview-controls: true
views:
  main_table:
    dataset: main_table
    desc: "Generated Metadata Report for {file_name}"
    render-table:
      columns:
        Run:
          link-to-url:
            "NCBI SRA":
              url: "https://www.ncbi.nlm.nih.gov/sra/?term={{value}}"
            "CNCB GSA":
              url: "https://ngdc.cncb.ac.cn/gsa/search?searchTerm={{value}}"
        BioProject:
          link-to-url:
            "NCBI BioProject":
              url: "https://www.ncbi.nlm.nih.gov/bioproject/?term={{value}}"
            "CNCB BioProject":
              url: "https://ngdc.cncb.ac.cn/bioproject/browse/{{value}}"
        ScientificName:
          link-to-url:
            "NCBI Taxonomy":
              url: "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name={{value}}"
        TaxID:
          link-to-url:
            "Taxonomy ID":
              url: "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={{value}}"
        PMID:
          link-to-url:
            "PubMed":
              url: "https://pubmed.ncbi.nlm.nih.gov/{{value}}/"
            "DOI":
              url: "https://doi.org/{{value}}"
"""
    with safe_open(yaml_path, 'wt') as f:
        f.write(yaml_content)

    try:
        print(f"📊 正在为 {file_name} 生成交互式报告...")
        result = subprocess.run(['datavzrd', yaml_path, '--output', report_dir],
                               capture_output=True, text=True)
        # 完整的状态判定与错误打印，杜绝静默失败
        if result.returncode == 0:
            print(f"✅ 报告已就绪: \033[94m{os.path.join(report_dir, 'index.html')}\033[0m")
            return os.path.join(report_dir, 'index.html')
        else:
            print(f"❌ datavzrd 生成失败！")
            print(f"【详细错误原因】:\n{result.stderr}")
    except Exception as e:
        print(f"❌ datavzrd 遇到致命异常: {e}")
    return None

def standardize_metadata_df(df, api_client, ai_model):
    if df is None: return pd.DataFrame()
    if df.empty: return df
    
    df = df.copy()
    null_words = {'not collected', 'missing', 'none', 'unknown', 'na', 'n/a', 'not applicable', 'not provided', '-', 'null', ''}
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].apply(lambda x: pd.NA if str(x).strip().lower().replace('_', ' ') in null_words else x)
    
    if 'Tissue' in df.columns:
        df['Tissue'] = df['Tissue'].apply(lambda x: str(x).strip().lower() if pd.notna(x) else x)
        mapping = {'leaves': 'leaf', 'fruits': 'fruit', 'roots': 'root', 'stems': 'stem', 'flowers': 'flower', 'seeds': 'seed', 'myceliums': 'mycelium'}
        df['Tissue'] = df['Tissue'].replace(mapping)
        
    if 'Location' in df.columns and api_client:
        unique_locs = df['Location'].dropna().unique().tolist()
        to_fix = [loc for loc in unique_locs if str(loc).strip().lower().replace('_', ' ') not in null_words]
        
        if to_fix:
            mapped_locs = {}
            for i in range(0, len(to_fix), 50):
                batch = to_fix[i:i+50]
                prompt_loc = """你是一个专业的地理位置数据清洗引擎。请将以下地点统一定制化为【国家, 省/州, 市/县_AI】格式。
1. 若只有国家和城市，凭常识补全省份（如 "China:Yinchuan_AI" -> "China, Ningxia, Yinchuan_AI"）。
2. 若缺失城市，用 Unknown 代替。
3. 必须使用逗号加空格分隔层级，且结尾必须带 "_AI"。
只返回 JSON 字典：{"原值": "新值"}"""
                try:
                    kwargs = build_api_kwargs(ai_model, prompt_loc, json.dumps(batch, ensure_ascii=False))
                    res = api_client.chat.completions.create(**kwargs)
                    batch_res = json.loads(clean_ai_json(res.choices[0].message.content))
                    mapped_locs.update(batch_res)
                except: pass
            if mapped_locs:
                df['Location'] = df['Location'].map(mapped_locs).fillna(df['Location'])
    return df

def merge_global_results(df_sra, df_gsa, out_dir, fill_date, mode, api_client, ai_model, ncbi_api, use_scholar):
    if df_sra is None: df_sra = pd.DataFrame()
    if df_gsa is None: df_gsa = pd.DataFrame()

    print("\n" + "="*65)
    print("🌍 正在执行 Global Unification (大一统合并与语义净化正交化)...")
    
    print("🧹 正在对合并前数据进行深度清洗 (清理无效占位符 / 规范部位 / AI地理位置规范化)...")
    df_sra = standardize_metadata_df(df_sra, api_client, ai_model)
    df_gsa = standardize_metadata_df(df_gsa, api_client, ai_model)
    
    if not df_gsa.empty and 'Run_Accession' in df_gsa.columns: 
        df_gsa = df_gsa.rename(columns={'Run_Accession': 'Run'})

    if not df_sra.empty: df_sra = df_sra.loc[:, ~df_sra.columns.duplicated()]
    if not df_gsa.empty: df_gsa = df_gsa.loc[:, ~df_gsa.columns.duplicated()]

    df_global = pd.concat([df_sra, df_gsa], ignore_index=True)
    if df_global.empty: return

    # 🚀 AI 9列正交化自检防火墙 (专治 GSA 的长文本小作文和物种名混入杂质)
    if api_client:
        print("🤖 正在启动 AI 9列正交化自检防火墙 (隔离语义污染并浓缩长文本)...")
        cols_to_check = ['CollectionDate', 'Location', 'Source', 'Tissue', 'Age_GrowthStage', 'ScientificName', 'LibrarySource', 'CenterName', 'BioProject']
        
        for idx, row in tqdm(df_global.iterrows(), total=len(df_global), desc="Data Sanitization"):
            target_dict = {k: row.get(k, "Not_Provided") for k in cols_to_check if k in df_global.columns}
            try:
                kwargs = build_api_kwargs(ai_model, PROMPT_DATA_SANITIZER, json.dumps(target_dict, ensure_ascii=False))
                res = api_client.chat.completions.create(**kwargs)
                clean_data = json.loads(clean_ai_json(res.choices[0].message.content))
                
                for k in cols_to_check:
                    if k in clean_data and k in df_global.columns:
                        val = str(clean_data[k]).strip()
                        # 忽略大小写，防止 AI 输出小写的 not_provided 覆盖掉已填补的时间
                        if val.lower() not in ["", "none", "nan", "not_provided"]:
                            df_global.at[idx, k] = val
            except: pass
    # 强制将 TaxID 转为字符串并去除 .0 后缀
    if 'TaxID' in df_global.columns:
        df_global['TaxID'] = df_global['TaxID'].astype(str).str.replace(r'\.0$', '', regex=True).replace(['nan', '<NA>', 'None', ''], pd.NA)

    if 'TaxID' not in df_global.columns: df_global['TaxID'] = pd.NA
    df_global['TaxID'] = df_global['TaxID'].astype(object)

    print("🧬 正在连接 NCBI Taxonomy 引擎查缺补漏...")
    missing_tax_mask = df_global['TaxID'].isna() | df_global['TaxID'].astype(str).str.strip().str.lower().isin(['not_provided', 'nan', 'none', 'unknown', '', 'pd.na', '<na>'])
    unique_names = df_global.loc[missing_tax_mask, 'ScientificName'].dropna().unique().tolist()
    if unique_names:
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            tax_res = list(tqdm(executor.map(resolve_taxonomy, unique_names), total=len(unique_names), desc="Taxonomy API Check"))
        taxid_map = {orig: res[0] for orig, res in zip(unique_names, tax_res)}
        name_map = {orig: res[1] for orig, res in zip(unique_names, tax_res)}
        df_global.loc[missing_tax_mask, 'TaxID'] = df_global.loc[missing_tax_mask, 'ScientificName'].map(taxid_map)
        df_global.loc[missing_tax_mask, 'ScientificName'] = df_global.loc[missing_tax_mask, 'ScientificName'].map(name_map).fillna(df_global.loc[missing_tax_mask, 'ScientificName'])

    if 'PMID' not in df_global.columns: df_global['PMID'] = pd.NA
    missing_pmid_mask = df_global['PMID'].isna() | df_global['PMID'].astype(str).str.strip().str.lower().isin(['not_provided', '', 'nan', 'unknown', 'none', '<na>'])
    orphan_bps = df_global.loc[missing_pmid_mask, 'BioProject'].dropna().unique().tolist()

    if orphan_bps and api_client:
        bp_res_dir = os.path.join(out_dir, "BioProject_Results")
        bp_cache_dir = os.path.join(bp_res_dir, "1_cache")
        os.makedirs(bp_cache_dir, exist_ok=True)
        
        print("\n" + "="*65 + f"\n📖 发现 {len(orphan_bps)} 个未查明文献的项目，启动深度溯源引擎...")
        tracer = BioProjectTracer(ncbi_api, api_client, ai_model, use_scholar)
        pmid_mapping, trace_details = {}, []
        success, failed = 0, 0
        
        for bp in tqdm(orphan_bps, desc="Literature Tracing"):
            cache_name = str(bp) + '_trace.json'
            cache_file = os.path.join(bp_cache_dir, cache_name)
            if os.path.exists(cache_file):
                with open(cache_file, 'r', encoding='utf-8') as f: pub_info = json.load(f)
                art_id = pub_info.get('pmid') or pub_info.get('doi')
                pmid_mapping[bp] = art_id
                tqdm.write(f"  [Tracer] {bp.ljust(12)} -> [✓本地缓存命中] => {art_id}")
                success += 1
                continue
                
            try:
                pub_info = tracer.trace(bp)
                if pub_info and (pub_info.get('pmid') or pub_info.get('doi')):
                    art_id = pub_info.get('pmid') or pub_info.get('doi')
                    pmid_mapping[bp] = art_id
                    
                    detail = {"BioProject": bp, "Status": "Success", "Article_ID": art_id}
                    detail.update(pub_info)
                    trace_details.append(detail)
                    success += 1
                    
                    with safe_open(cache_file, 'wt') as f: json.dump(pub_info, f, ensure_ascii=False)
                    if 'category' in str(pub_info):
                        tqdm.write(f"  [Tracer] {bp.ljust(12)} -> [✓API精准判定] => {art_id}")
                    else:
                        tqdm.write(f"  [Tracer] {bp.ljust(12)} -> [✓规则引擎兜底] => {art_id}")
                else:
                    trace_details.append({"BioProject": bp, "Status": "Not_Found"})
                    failed += 1
                    tqdm.write(f"  [Tracer] {bp.ljust(12)} -> [❌全网未找到]")
            except Exception as e:
                trace_details.append({"BioProject": bp, "Status": "Error", "Error": str(e)})
                failed += 1
                tqdm.write(f"  [Tracer] {bp.ljust(12)} -> [❌崩溃: {e}]")

        if trace_details:
            trace_f = os.path.join(bp_res_dir, "BioProject_Trace_Details.csv")
            with safe_open(trace_f, 'wt') as f:
                pd.DataFrame(trace_details).to_csv(f, index=False, encoding='utf-8-sig', lineterminator='\n')
        df_global.loc[missing_pmid_mask, 'PMID'] = df_global.loc[missing_pmid_mask, 'BioProject'].map(pmid_mapping).fillna(df_global.loc[missing_pmid_mask, 'PMID'])
        print(f"📊 溯源战报：成功追溯 {success} 项。详情见 BioProject_Results 目录。")

    for col in FINAL_14:
        if col not in df_global.columns: df_global[col] = pd.NA
    df_final_14 = df_global[FINAL_14].copy().replace(["Not_Provided", "not_provided", "None", ""], pd.NA)
    save_dual_format(df_final_14, os.path.join(out_dir, "Global_Unified_Metadata_Core14"))
    full_f = save_dual_format(df_global, os.path.join(out_dir, "Global_Unified_Metadata_Full"))
    generate_global_datavzrd(full_f, out_dir)
    
    print(f"🎉 完美收工！数据已纯净化。")
    print(f"👉 【全维溯源表】: Global_Unified_Metadata_Full.csv")
    print(f"👉 【核心标准表】: Global_Unified_Metadata_Core14.tsv")

def main():
    parser = argparse.ArgumentParser(description="🌍 Global Metadata Engine v12.11 终极完备版")
    parser.add_argument("-i", "--input", required=True, help="混合编号文件（或直接给单个 Run 编号）")
    parser.add_argument("-o", "--outdir", dest="outdir", default=None,
                        help="输出目录（默认 meta_search/info/）")
    parser.add_argument("-m", "--mode", choices=['local', 'api', 'both'], default='both')
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--kimi-api", help="Kimi API Key")
    parser.add_argument("--deepseek-api", help="DeepSeek API Key")
    parser.add_argument("--deepseek-model", default="deepseek-v4-flash", choices=["deepseek-v4-flash", "deepseek-v4-pro"], help="指定 DeepSeek 模型版本 (v4-flash 快速 / v4-pro 深度思考)")
    parser.add_argument("--ncbi-api", help="NCBI API Key")
    parser.add_argument("--use-scholar", action="store_true", help="允许调用 Google Scholar")
    parser.add_argument("--fill-date", action="store_true", help="时间兜底")
    args = parser.parse_args()

    if not args.outdir:
        try:
            from vp.config import DIRS
            args.outdir = os.path.join(DIRS['meta_search'], 'info')
        except Exception:
            args.outdir = './Global_Metadata_Results'

    if args.ncbi_api: os.environ["NCBI_API_KEY"] = args.ncbi_api
    # 环境变量回退 (export DEEPSEEK_API_KEY / NCBI_API_KEY)
    if not args.deepseek_api and os.environ.get("DEEPSEEK_API_KEY"):
        args.deepseek_api = os.environ["DEEPSEEK_API_KEY"]
    if not args.ncbi_api and os.environ.get("NCBI_API_KEY"):
        args.ncbi_api = os.environ["NCBI_API_KEY"]
        os.environ["NCBI_API_KEY"] = args.ncbi_api
    client, model = None, None
    if args.mode in ['api', 'both']:
        if args.kimi_api:
            client, model = OpenAI(api_key=args.kimi_api, base_url="https://api.moonshot.cn/v1"), "moonshot-v1-32k"
            print("🤖 挂载引擎: Kimi")
        elif args.deepseek_api:
            client, model = OpenAI(api_key=args.deepseek_api, base_url="https://api.deepseek.com"), args.deepseek_model
            print(f"🤖 挂载引擎: DeepSeek ({model})")
        else: parser.error("API模式需提供 Key (--deepseek-api 或 export DEEPSEEK_API_KEY)")

    os.makedirs(args.outdir, exist_ok=True)
    input_data = [line.strip() for line in open(args.input, 'r')] if os.path.isfile(args.input) else [args.input]
    input_list = list(dict.fromkeys([x for x in input_data if x]))
    sra_list, gsa_list, _ = get_lists_by_regex(input_list)
    print(f"\n🚀 路由分流 | SRA: {len(sra_list)} | GSA: {len(gsa_list)}")
    
    df_sra, df_gsa = pd.DataFrame(), pd.DataFrame()
    if sra_list: df_sra = SRAPipeline(args.mode, client, model, args.outdir, args.fill_date).process_all(sra_list, args.threads)
    if gsa_list: df_gsa = GSAPipeline(args.mode, client, model, args.outdir, args.fill_date).process_all(gsa_list)
    merge_global_results(df_sra, df_gsa, args.outdir, args.fill_date, args.mode, client, model, args.ncbi_api, args.use_scholar)

if __name__ == "__main__":
    main()

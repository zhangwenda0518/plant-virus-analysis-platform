# -*- coding: utf-8 -*-
"""
NCBI Taxonomy 准备与解析。
- prepare_taxonomy(): 下载 new_taxdump.tar.gz 并解压 nodes/names/merged 到平台 taxonomy 目录
- TaxonomyGraph: taxid -> 名称/rank/祖先链，供 kreport 解析与病毒注释使用
下载仅允许 https + NCBI 官方域名白名单，解析 IP 阻断私网/环回/链路本地地址，禁用重定向。
"""
import os
import socket
import ipaddress
import urllib.request
from urllib.parse import urlparse

from .config import DIRS
from .utils import check_path, safe_open

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
REQUIRED_MEMBERS = ('nodes.dmp', 'names.dmp', 'merged.dmp')
ALLOWED_HOSTS = {'ftp.ncbi.nlm.nih.gov', 'www.ncbi.nlm.nih.gov'}


# ------------------------------------------------------------------
# SSRF 防护：URL 白名单 + IP 边界校验 + 禁止重定向
# ------------------------------------------------------------------
class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None      # 一律不跟随重定向


def validate_ncbi_url(url):
    u = urlparse(str(url))
    if u.scheme != 'https':
        raise ValueError(f"仅允许 https 协议: {url}")
    if u.hostname not in ALLOWED_HOSTS:
        raise ValueError(f"下载域名不在白名单({sorted(ALLOWED_HOSTS)}): {u.hostname}")
    infos = socket.getaddrinfo(u.hostname, 443, proto=socket.IPPROTO_TCP)
    for info in infos:
        ip = ipaddress.ip_address(info[4][0])
        if (ip.is_private or ip.is_loopback or ip.is_link_local
                or ip.is_reserved or ip.is_multicast or ip.is_unspecified):
            raise ValueError(f"目标域名解析到受限地址 {ip}，已阻断")
    return str(url)


# ------------------------------------------------------------------
# taxonomy 准备
# ------------------------------------------------------------------
def taxonomy_ready(tax_dir=None):
    tax_dir = check_path(tax_dir or DIRS['taxonomy'], must_exist=False, in_platform=True)
    nodes = check_path(os.path.join(tax_dir, 'nodes.dmp'), must_exist=False, in_platform=True)
    names = check_path(os.path.join(tax_dir, 'names.dmp'), must_exist=False, in_platform=True)
    return os.path.isfile(nodes) and os.path.isfile(names)


def taxid_lookup(taxid, tax_dir=None):
    """校验 taxid 在 NCBI 分类中是否存在/已合并。

    返回:
      ('ok', None)            - nodes.dmp 中存在
      ('merged', new_taxid)   - 已合并到新 taxid（建议改用新值）
      ('missing', None)       - 不存在（无效或本地 taxonomy 过旧）
    """
    t = int(taxid)
    tax_dir = check_path(tax_dir or DIRS['taxonomy'],
                         must_exist=True, in_platform=True)
    nodes = check_path(os.path.join(tax_dir, 'nodes.dmp'),
                       must_exist=True, in_platform=True)
    with safe_open(nodes) as f:
        for line in f:
            first = line.split('\t', 1)[0]
            if first.isdigit() and int(first) == t:
                return 'ok', None
    merged = check_path(os.path.join(tax_dir, 'merged.dmp'),
                        must_exist=False, in_platform=True)
    if os.path.isfile(merged):
        with safe_open(merged) as f:
            for line in f:
                parts = line.split('\t|\t')
                if len(parts) >= 2 and parts[0].strip().isdigit() \
                        and int(parts[0].strip()) == t:
                    # merged.dmp 字段尾带 "\t|"，取首个制表符前的数字
                    new = parts[1].split('\t')[0].strip()
                    return 'merged', (int(new) if new.isdigit() else new)
    return 'missing', None


def prepare_taxonomy(logger=None, force=False):
    """确保 taxonomy 就绪；未就绪则下载解压。返回 taxonomy 目录。"""
    tax_dir = check_path(DIRS['taxonomy'], must_exist=False, in_platform=True)
    os.makedirs(tax_dir, exist_ok=True)

    if taxonomy_ready(tax_dir) and not force:
        if logger:
            logger.log(f"Taxonomy 已就绪: {tax_dir}")
        return tax_dir

    tar_name = 'new_taxdump.tar.gz'
    tar_path = check_path(os.path.join(DIRS['databases'], tar_name),
                          must_exist=False, in_platform=True)
    os.makedirs(DIRS['databases'], exist_ok=True)

    if not os.path.isfile(tar_path) or force:
        if logger:
            logger.log(f"下载 NCBI new_taxdump (~57MB): {TAXDUMP_URL}")
        _download(TAXDUMP_URL, tar_path, logger)

    if logger:
        logger.log("解压 nodes.dmp / names.dmp / merged.dmp ...")
    import tarfile
    with tarfile.open(tar_path, 'r:gz') as tar:
        for member in tar.getmembers():
            base = os.path.basename(member.name)
            if base in REQUIRED_MEMBERS and member.isfile():
                member.name = base          # 重写为纯文件名，消除任何目录成分
                tar.extract(member, tax_dir, filter='data')   # data 过滤器阻断穿越/链接

    if not taxonomy_ready(tax_dir):
        raise RuntimeError("taxonomy 解压后校验失败（缺少 nodes.dmp/names.dmp），请检查下载文件")
    if logger:
        logger.log(f"Taxonomy 准备完成: {tax_dir}")
    return tax_dir


def _download(url, dest, logger=None, chunk=1024 * 512):
    """带进度回调的安全下载（URL 白名单 + IP 校验 + 禁重定向）。"""
    url = validate_ncbi_url(url)
    tmp = check_path(str(dest) + '.part', must_exist=False, in_platform=True)
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (vp-platform)'})
    opener = urllib.request.build_opener(_NoRedirect)
    with opener.open(req, timeout=120) as resp, safe_open(tmp, 'wb') as out:
        total = int(resp.headers.get('Content-Length', 0))
        done, last_pct = 0, -1
        while True:
            block = resp.read(chunk)
            if not block:
                break
            out.write(block)
            done += len(block)
            if total:
                pct = done * 100 // total
                if pct != last_pct and pct % 5 == 0:
                    last_pct = pct
                    if logger:
                        logger.log(f"  下载进度 {pct}% ({done // 1048576}/{total // 1048576} MB)")
    final = check_path(dest, must_exist=False, in_platform=True)
    os.replace(tmp, final)


# ------------------------------------------------------------------
# Taxonomy 解析
# ------------------------------------------------------------------
class TaxonomyGraph:
    """轻量 NCBI taxonomy：parent/rank/name/merged，支持祖先链查询。"""

    def __init__(self, tax_dir=None):
        self.tax_dir = check_path(tax_dir or DIRS['taxonomy'], must_exist=True)
        self.parent = {}       # taxid -> parent_taxid
        self.rank = {}         # taxid -> rank
        self.names = {}        # taxid -> scientific name
        self.merged = {}       # old_taxid -> new_taxid
        self._loaded = False

    def load(self):
        if self._loaded:
            return
        nodes = check_path(os.path.join(self.tax_dir, 'nodes.dmp'), must_exist=True)
        names = check_path(os.path.join(self.tax_dir, 'names.dmp'), must_exist=True)
        with safe_open(nodes) as f:
            for line in f:
                parts = [x.strip() for x in line.split('|')]
                if len(parts) >= 3:
                    t, p, r = int(parts[0]), int(parts[1]), parts[2]
                    self.parent[t] = p
                    self.rank[t] = r
        with safe_open(names) as f:
            for line in f:
                parts = [x.strip() for x in line.split('|')]
                if len(parts) >= 4 and parts[3] in ('scientific name', ''):
                    self.names[int(parts[0])] = parts[1]
        merged = os.path.join(self.tax_dir, 'merged.dmp')
        if os.path.isfile(merged):
            with safe_open(merged) as f:
                for line in f:
                    parts = [x.strip() for x in line.split('|')]
                    if len(parts) >= 2:
                        try:
                            self.merged[int(parts[0])] = int(parts[1])
                        except ValueError:
                            pass
        self._loaded = True

    def resolve(self, taxid):
        taxid = int(taxid)
        return self.merged.get(taxid, taxid)

    def name(self, taxid):
        self.load()
        return self.names.get(self.resolve(taxid), f"taxid:{taxid}")

    def rank_of(self, taxid):
        self.load()
        return self.rank.get(self.resolve(taxid), 'no rank')

    def lineage(self, taxid):
        """从自身到根的 [(taxid, rank, name), ...]；环路与超长保护。"""
        self.load()
        chain, seen = [], set()
        t = self.resolve(taxid)
        for _ in range(64):
            if t in seen or t not in self.parent:
                break
            seen.add(t)
            chain.append((t, self.rank.get(t, 'no rank'), self.names.get(t, f'taxid:{t}')))
            if t == 1:
                break
            t = self.parent[t]
        return chain

    def lineage_names(self, taxid):
        """从根到自身的名称列表（GUI 显示用）。"""
        return [name for _, _, name in reversed(self.lineage(taxid))]

    def get_ancestor_at_rank(self, taxid, rank):
        """返回指定 rank 的祖先 taxid（如 'family' / 'genus' / 'species'），无则 None。"""
        for t, r, _ in self.lineage(taxid):
            if r == rank:
                return t
        return None

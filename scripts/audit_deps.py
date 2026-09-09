# -*- coding: utf-8 -*-
"""依赖审计：扫描全部源码的 import，找出 requirements.txt 缺失的第三方依赖。

跑法：  C:\\Python312\\python.exe scripts\\audit_deps.py
输出：  已安装 / 缺失 / 未被 requirements 覆盖 三张表，退出码 1 表示有缺失。

打包前必跑：PyInstaller 不会自动带依赖，缺一个就是运行期 ModuleNotFoundError。
"""
import ast
import os
import sys
import importlib.util

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 标准库（sys.stdlib_module_names 在 3.10+ 可用，这里补一份兜底）
try:
    STDLIB = set(sys.stdlib_module_names)
except AttributeError:                                   # pragma: no cover
    STDLIB = set()
STDLIB |= {'os', 'sys', 're', 'json', 'time', 'shutil', 'glob', 'csv',
           'math', 'random', 'uuid', 'threading', 'subprocess', 'collections',
           'argparse', 'dataclasses', 'gzip', 'io', 'itertools', 'functools',
           'pathlib', 'types', 'typing', 'warnings', 'unicodedata', 'hashlib',
           'tempfile', 'traceback', 'string', 'struct', 'copy', 'enum',
           'multiprocessing', 'socket', 'webbrowser', 'email', 'urllib',
           'http', 'xml', 'sqlite3', 'zipfile', 'tarfile', 'base64', 'errno',
           'platform', 'stat', 'textwrap', 'difflib', 'bisect', 'heapq',
           'operator', 'contextlib', 'inspect', 'abc', 'numbers', 'decimal',
           'datetime', 'calendar', 'locale', 'gettext', 'logging'}

# 本地模块（平台自己的包）
LOCAL = {'vp', 'app', 'main', 'scripts', 'tests', 'known_virus_suite'}


def _local_submodules():
    """扫描平台本地包目录下的模块名（如 known_virus_suite/kv_*.py）。

    同包内既用 `from known_virus_suite import x` 也用裸 `import kv_common`，
    后者的顶层名不在 LOCAL 里会误报为缺失第三方包。
    """
    names = set()
    for pkg in ('known_virus_suite', 'vp'):
        d = os.path.join(ROOT, pkg)
        if not os.path.isdir(d):
            continue
        for fn in os.listdir(d):
            if fn.endswith('.py'):
                names.add(fn[:-3])
    return names


LOCAL_SUBMODULES = None  # 延迟到 ROOT 定义之后初始化

# 导入名 → pip 包名（不一致的才需要列）
PIP_NAME = {
    'flask': 'Flask', 'werkzeug': 'Werkzeug', 'jinja2': 'Jinja2',
    'plotly': 'plotly', 'matplotlib': 'matplotlib', 'numpy': 'numpy',
    'pandas': 'pandas', 'Bio': 'biopython', 'dna_features_viewer':
        'dna_features_viewer', 'pyrodigal': 'pyrodigal', 'pyhmmer': 'pyhmmer',
    'primer3': 'primer3-py', 'selenium': 'selenium', 'requests': 'requests',
    'openpyxl': 'openpyxl', 'PIL': 'Pillow', 'lxml': 'lxml',
    'sklearn': 'scikit-learn', 'scipy': 'SciPy', 'yaml': 'PyYAML',
    'psutil': 'psutil', 'tqdm': 'tqdm', 'xgboost': 'xgboost',
    'openai': 'openai', 'docx': 'python-docx', 'pptx': 'python-pptx',
    'pypdf': 'pypdf', 'fitz': 'PyMuPDF', 'cv2': 'opencv-python',
    'streamlit': 'streamlit', 'pycirclize': 'pycirclize',
}

# 已知可选依赖：缺失时功能降级，不阻断运行
OPTIONAL = {'selenium', 'gbdraw', 'dna_features_viewer', 'pyhmmer',
            'primer3', 'lovis4u', 'openai', 'psutil', 'pycirclize',
            'distinctipy', 'taxburst'}


def iter_py_files():
    # docs/bioaider_ref_src 是 BioAider 的参考源码（他项目 GUI，仅供对照），
    # 不是平台运行依赖；tests/ 的造数脚本同理不参与打包；
    # git-repo/ vendor/ _archive/ 是第三方镜像与归档，也不算平台依赖
    # （曾因此把 git-repo/MultiVirusConsensus 的 pysam 误报为必需依赖）。
    skip = {'dist', 'tools', 'open-virome', 'node_modules', '.git',
            '__pycache__', 'build', 'host-db', 'virus-db', 'databases',
            'results', 'downloads', 'logs', 'tasks', 'tool_runs',
            'meta_search', 'logan', 'submissions', 'fastq', 'uploads',
            'bin', 'logan', 'docs', 'tests',
            'git-repo', 'vendor', '_archive', '.pytest_cache'}
    for cur, dirs, files in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in skip and not d.startswith('.')]
        for fn in files:
            if fn.endswith('.py'):
                yield os.path.join(cur, fn)


def collect_imports():
    """返回 {顶层模块名: [引用它的文件相对路径]}，并跳过 try/except ImportError 块。"""
    found = {}
    for path in iter_py_files():
        try:
            with open(path, encoding='utf-8', errors='replace') as f:
                tree = ast.parse(f.read(), filename=path)
        except (SyntaxError, UnicodeDecodeError, ValueError):
            continue
        rel = os.path.relpath(path, ROOT).replace(os.sep, '/')
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    top = a.name.split('.')[0]
                    found.setdefault(top, []).append(rel)
            elif isinstance(node, ast.ImportFrom):
                if node.level:            # 相对导入（from .xxx）→ 本地
                    continue
                if node.module:
                    top = node.module.split('.')[0]
                    found.setdefault(top, []).append(rel)
    return found


def analyze():
    """完整分析：返回 {third, missing, uncovered, installed, hard_missing}。"""
    found = collect_imports()
    local_sub = _local_submodules()   # 本地包内的子模块（kv_*.py 等）
    third = {m: fs for m, fs in found.items()
             if m not in STDLIB and m not in LOCAL
             and m not in local_sub and not m.startswith('_')}

    req_path = os.path.join(ROOT, 'requirements.txt')
    req = set()
    if os.path.isfile(req_path):
        with open(req_path, encoding='utf-8') as f:
            for line in f:
                line = line.split('#', 1)[0].strip()
                if not line:
                    continue
                name = line.split('==')[0].split('>=')[0].split('[')[0]
                req.add(name.strip().lower().replace('-', '_'))

    installed, missing = {}, []
    for m in sorted(third):
        ok = importlib.util.find_spec(m) is not None
        installed[m] = ok
        if not ok:
            missing.append(m)

    uncovered = []
    for m in sorted(third):
        pip = PIP_NAME.get(m, m)
        if pip.lower().replace('-', '_') not in req:
            uncovered.append((m, pip))

    hard_missing = [m for m in missing if m not in OPTIONAL]
    return {'third': third, 'missing': missing, 'uncovered': uncovered,
            'installed': installed, 'hard_missing': hard_missing,
            'req': req,
            'all_ready': not hard_missing and not missing}


def missing_required():
    """返回缺失的必需依赖名列表（空 = 全部就绪）。供 selfcheck / CI 用。"""
    return analyze()['hard_missing']


def main():
    res = analyze()
    third, missing = res['third'], res['missing']
    uncovered, installed = res['uncovered'], res['installed']
    hard_missing = res['hard_missing']
    req = res['req']

    print('扫描到第三方依赖 %d 个（引用文件已去重统计）' % len(third))
    print('\n── 未安装 ──')
    if not missing:
        print('  （无）')
    for m in missing:
        tag = '可选' if m in OPTIONAL else '必需'
        print('  ✘ %-24s [%s]  pip install %s   ← %s'
              % (m, tag, PIP_NAME.get(m, m), ', '.join(sorted(set(third[m]))[:3])))

    print('\n── requirements.txt 未覆盖 ──')
    if not uncovered:
        print('  （无）')
    for m, pip in uncovered:
        tag = '可选' if m in OPTIONAL else '必需'
        print('  △ %-24s [%s]  建议加入: %s' % (m, tag, pip))

    print('\n── 已安装且已覆盖 ──')
    for m in sorted(third):
        pip = PIP_NAME.get(m, m)
        if installed[m] and pip.lower().replace('-', '_') in req:
            print('  ✔ %s' % m)

    hard_missing = [m for m in missing if m not in OPTIONAL]
    print('\n' + '=' * 46)
    if hard_missing:
        print('缺失必需依赖 %d 个: %s' % (len(hard_missing), ', '.join(hard_missing)))
        print('安装: pip install ' + ' '.join(PIP_NAME.get(m, m)
                                              for m in hard_missing))
        return 1
    if missing:
        print('仅缺可选依赖（对应功能会降级）: %s' % ', '.join(missing))
    else:
        print('全部依赖已安装')
    return 0


if __name__ == '__main__':
    sys.exit(main())

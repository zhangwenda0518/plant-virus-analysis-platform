# -*- coding: utf-8 -*-
"""对 VOG HMM 库重建 hmmpress 索引（.h3f/.h3i/.h3m/.h3p）。

背景：databases/hmm/vogdb/vog_all.hmm 原有的 .h3* 索引残缺（仅 8KB），
导致 VOG 扫描每次都走慢速的逐条解析路径。本脚本用 pyhmmer 的 hmmpress
重新生成正确的二进制索引，让 hmmscan 吃到加速红利。

用法（用平台实际运行的系统 Python）：
    C:\\Python312\\python.exe scripts/press_vogdb.py
"""
import os
import sys
import time

import pyhmmer
import pyhmmer.plan7
import pyhmmer.easel

_REAL = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'databases', 'hmm', 'vogdb', 'vog_all.hmm')


def _ascii_path(path):
    """pyhmmer 底层 C 库在 Windows 打不开非 ASCII 路径；经 %TEMP% 的
    ASCII junction（vp_ascii_root → 平台根）绕行。与 vp/hmm_annot.py
    的 _ascii_path 同思路。"""
    if str(path).isascii():
        return path
    import subprocess
    temp_root = os.path.join(os.environ.get('TEMP', ''), 'vp_ascii_root')
    if not os.path.isdir(temp_root):
        platform_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        subprocess.run(['cmd', '/c', 'mklink', '/J', temp_root, platform_root],
                       capture_output=True)
    if os.path.isdir(temp_root):
        rel = os.path.relpath(path, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        cand = os.path.join(temp_root, rel)
        if str(cand).isascii():
            return cand
    return path


def main():
    VOG_HMM = _ascii_path(_REAL)
    if not os.path.isfile(VOG_HMM):
        print(f'✘ 未找到 {_REAL}')
        return 1
    print(f'源库(真实): {_REAL} ({os.path.getsize(_REAL)/1024**3:.2f} GB)')
    print(f'源库(ASCII): {VOG_HMM}')

    t0 = time.time()
    alpha = pyhmmer.easel.Alphabet.amino()
    hmms = []
    n = 0
    # 流式读取，避免一次性把 4.4G 文本全部读进内存（仍会生成全部 HMM 对象）
    with pyhmmer.plan7.HMMFile(VOG_HMM) as hf:
        for hmm in hf:
            hmms.append(hmm)
            n += 1
    print(f'读取 {n:,} 个 HMM，耗时 {time.time()-t0:.1f}s')

    t1 = time.time()
    pyhmmer.hmmpress(hmms, VOG_HMM)   # 生成 vog_all.hmm.h3f/.h3i/.h3m/.h3p
    print(f'hmmpress 完成，耗时 {time.time()-t1:.1f}s')

    # 验证生成的索引
    print('\n── 生成的索引文件 ──')
    ok = True
    for ext in ('.h3f', '.h3i', '.h3m', '.h3p'):
        p = VOG_HMM + ext
        if os.path.isfile(p):
            sz = os.path.getsize(p)
            print(f'  ✔ {os.path.basename(p):24s} {sz/1024**2:8.1f} MB')
            if ext in ('.h3m', '.h3p') and sz < 1024 * 1024:
                print(f'    ⚠ 该文件偏小，可能未完整生成')
                ok = False
        else:
            print(f'  ✘ 缺失 {os.path.basename(p)}')
            ok = False
    print(f'\n{"✔ 索引重建成功" if ok else "✘ 索引可能不完整，请检查"}')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

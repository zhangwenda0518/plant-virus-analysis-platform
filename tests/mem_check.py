# -*- coding: utf-8 -*-
"""环境内存自检：验证当前环境进程可用的内存上限。

用法: 双击 环境自检.bat（或 python tests/mem_check.py）
正常桌面环境应能分配数 GB；若只能分配约 2GB，说明程序运行在受限沙箱中。
"""
import ctypes
import sys

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')


class MEMORYSTATUSEX(ctypes.Structure):
    _fields_ = [('dwLength', ctypes.c_ulong),
                ('dwMemoryLoad', ctypes.c_ulong),
                ('ullTotalPhys', ctypes.c_ulonglong),
                ('ullAvailPhys', ctypes.c_ulonglong),
                ('ullTotalPageFile', ctypes.c_ulonglong),
                ('ullAvailPageFile', ctypes.c_ulonglong),
                ('ullTotalVirtual', ctypes.c_ulonglong),
                ('ullAvailVirtual', ctypes.c_ulonglong),
                ('ullAvailExtendedVirtual', ctypes.c_ulonglong)]


def main():
    stat = MEMORYSTATUSEX()
    stat.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
    gb = lambda v: v / 2 ** 30
    print('=' * 52)
    print('  植物病毒分析平台 - 环境内存自检')
    print('=' * 52)
    print(f'  物理内存总量        : {gb(stat.ullTotalPhys):6.1f} GB')
    print(f'  物理内存当前可用    : {gb(stat.ullAvailPhys):6.1f} GB')
    print(f'  提交上限(RAM+页面)  : {gb(stat.ullTotalPageFile):6.1f} GB')
    print(f'  进程虚拟地址可用    : {gb(stat.ullAvailVirtual):6.1f} GB')
    print('-' * 52)

    # 实际分配测试：逐步申请 512MB，目标 8GB，随后立即释放
    blocks = []
    chunk = 512 * 2 ** 20
    target = 8 * 2 ** 30
    got = 0
    try:
        while got < target:
            blocks.append(bytearray(chunk))     # 触碰内存确保真实分配
            blocks[-1][0:4096] = b'\x01' * 4096
            got += chunk
            print(f'\r  正在测试分配 ... {gb(got):4.1f} GB', end='', flush=True)
    except MemoryError:
        pass
    finally:
        blocks = None
    print()
    print('-' * 52)
    if got >= 4 * 2 ** 30:
        print(f'  [PASS] 本进程实际可分配 {gb(got):.1f} GB 内存。')
        print('         宿主库建库 / SPAdes 组装在此环境可正常运行。')
    elif got >= 3 * 2 ** 30:
        print(f'  [OK]   本进程实际可分配 {gb(got):.1f} GB 内存。')
        print('         建库/组装可用；建议同时只跑一个重内存任务。')
    else:
        print(f'  [WARN] 本进程仅能分配 {gb(got):.1f} GB！')
        print('         若这是直接双击运行的结果，请确认：')
        print('         - 没有在 AI 助手/沙箱/容器的终端里运行')
        print('         - 其他大型程序未占满内存（关闭后重试）')
    print('=' * 52)


if __name__ == '__main__':
    main()
    try:
        input('\n按回车键退出...')
    except EOFError:
        pass

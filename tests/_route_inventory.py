# -*- coding: utf-8 -*-
"""路由清单基线：app.py 拆分（单体化治理）前后的行为等价性守卫。

背景：app.py 有 180 个路由。拆成 vp/web/ 下的多个 Blueprint 后，
最大的风险是「漏注册某个 blueprint」「某条路由的方法变了」
「某条路由被静默改名」。tests/_it_platform.py 只覆盖 21 个页面 +
少数 API，不足以守住全部路由，本脚本补齐这一层。

用法：
    python tests/_route_inventory.py --save      # 拆分前：写入基线
    python tests/_route_inventory.py             # 拆分后：比对，有差异则退出码 1

比对口径：
  - rule + methods  ：严格比对，任何增删改都算失败
  - endpoint 名      ：Blueprint 化后必然变化（page_tools -> tools.page_tools），
                       故只作为「重命名报告」输出，不判失败

只读，不启动分析任务。
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

BASELINE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        '_route_baseline.json')


def collect():
    """收集全部路由：rule -> sorted(methods)，另记 endpoint 便于报告。"""
    import app as appmod
    rules = {}
    endpoints = {}
    for r in appmod.app.url_map.iter_rules():
        if r.endpoint == 'static':
            continue
        methods = sorted(m for m in r.methods if m not in ('HEAD', 'OPTIONS'))
        rules.setdefault(r.rule, set()).update(methods)
        endpoints[r.rule] = r.endpoint
    return ({k: sorted(v) for k, v in sorted(rules.items())},
            {k: endpoints[k] for k in sorted(endpoints)})


def main():
    rules, endpoints = collect()
    print(f'路由数: {len(rules)}')

    if '--save' in sys.argv:
        with open(BASELINE, 'w', encoding='utf-8', newline='\n') as fh:
            json.dump({'rules': rules, 'endpoints': endpoints}, fh,
                      ensure_ascii=False, indent=1, sort_keys=True)
        print(f'已写入基线: {os.path.relpath(BASELINE)}')
        return 0

    if not os.path.isfile(BASELINE):
        print(f'缺少基线文件 {BASELINE}，先运行 --save')
        return 2

    with open(BASELINE, 'r', encoding='utf-8') as fh:
        base = json.load(fh)
    base_rules = base['rules']
    base_eps = base.get('endpoints', {})

    missing = sorted(set(base_rules) - set(rules))
    added = sorted(set(rules) - set(base_rules))
    changed = []
    for rule in sorted(set(base_rules) & set(rules)):
        if base_rules[rule] != rules[rule]:
            changed.append((rule, base_rules[rule], rules[rule]))

    renamed = [(r, base_eps.get(r), endpoints.get(r))
               for r in sorted(set(base_rules) & set(rules))
               if base_eps.get(r) != endpoints.get(r)]

    ok = True
    if missing:
        ok = False
        print(f'\n✗ 丢失路由 {len(missing)} 条:')
        for r in missing:
            print(f'    {r}  {base_rules[r]}')
    if added:
        ok = False
        print(f'\n✗ 新增路由 {len(added)} 条（拆分不应新增）:')
        for r in added:
            print(f'    {r}  {rules[r]}')
    if changed:
        ok = False
        print(f'\n✗ 方法变化 {len(changed)} 条:')
        for r, b, n in changed:
            print(f'    {r}: {b} -> {n}')

    if renamed:
        print(f'\n· endpoint 重命名 {len(renamed)} 条（Blueprint 化预期如此，'
              f'只要模板未用 url_for 即无影响）:')
        for r, b, n in renamed[:15]:
            print(f'    {r}: {b} -> {n}')
        if len(renamed) > 15:
            print(f'    ... 另有 {len(renamed) - 15} 条')

    print('\n' + ('✔ 路由清单与基线一致' if ok else '✗ 路由清单存在差异'))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())

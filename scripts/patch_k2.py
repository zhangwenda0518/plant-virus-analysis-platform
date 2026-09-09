# -*- coding: utf-8 -*-
"""Kraken2 库转换端点 + UI 按钮。

【已执行完毕，存档备查】2026-09-07 核验：
补丁已落入 webapp/templates/build.html 内联 <script>（buildUniversalDb /
convertKraken2 两个函数在那里定义，不在 app.js），后端 /api/convert_kraken2
与 /api/build_universal_db 均已存在。
注意：本脚本第二段 patch 的目标文件写的是 build.html 而非 app.js（原意可能是
app.js），但结果上函数可用，不再追改。再次运行本脚本会因 assert count==1 失败。
"""
import io

def patch(path, pairs):
    s = io.open(path, encoding='utf-8').read()
    for old, new in pairs:
        assert s.count(old) == 1, f'{path}: count={s.count(old)} for {old[:70]!r}'
        s = s.replace(old, new, 1)
    io.open(path, 'w', encoding='utf-8').write(s)
    print(path, 'OK')

# ---------- app.py：端点 + dbs 状态 ----------
patch('app.py', [
    ("""@app.route('/api/build_universal_db', methods=['POST'])""",
     """@app.route('/api/convert_kraken2', methods=['POST'])
def api_convert_kraken2():
    \"\"\"Kraken2 库包/目录 → kunpeng 分片库（kunpeng hashshard，方式 C）。

    body: {tar: "databases/k2_viral_20260626.tar.gz"（平台内路径）,
           name: "k2viral"（目标库目录名 databases/<name>）,
           hash_capacity: "1G"}
    \"\"\"
    body = request.get_json(force=True) or {}
    tar = (body.get('tar') or '').strip()
    name = (body.get('name') or 'k2viral').strip()
    if not tar:
        abort(400, '缺少 Kraken2 库包路径')
    check_path(tar, must_exist=True)
    if not re.fullmatch(r'[A-Za-z0-9_\\-]+', name):
        abort(400, '库名仅限字母数字-_')
    db_dir = os.path.join(DIRS['databases'], name)

    def job(log, prog, cancel):
        from vp.kunpeng import convert_kraken2
        logger = TaskLogger(callback=log)
        prog('convert', 0.2, '解包 + hashshard 转换（见日志）')
        res = convert_kraken2(tar, db_dir,
                              hash_capacity=body.get('hash_capacity') or '1G',
                              logger=logger)
        prog('done', 1.0, '完成')
        logger.close()
        return {'db_dir': res}

    tid = tm.start(cfg.tr(f'Kraken2 库转换 {name}', f'Kraken2 convert {name}'),
                   job)
    return jsonify({'task': tid})


@app.route('/api/build_universal_db', methods=['POST'])"""),
    ("""    # 通用参考库（refvirus / rvdb）
    for key, sub in (('refvirus', 'refvirus_db'), ('rvdb', 'rvdb_db')):""",
     """    # 通用参考库（refvirus / rvdb / 其它 kunpeng 库目录自动探测）
    for key, sub in (('refvirus', 'refvirus_db'), ('rvdb', 'rvdb_db'),
                     ('k2viral', 'k2viral_db')):"""),
])
print('app.py k2 endpoints OK')

# ---------- build.html：Kraken2 转换按钮 ----------
patch('webapp/templates/build.html', [
    ("""        <div style="display:flex;gap:10px;margin-top:10px">
          <button class="btn" onclick="buildUniversalDb('refvirus', this)"
                  title="NCBI RefSeq Viral（推荐先建这个，速度快）">
            🧬 构建 RefSeq 病毒库</button>
          <button class="btn" onclick="buildUniversalDb('rvdb', this)"
                  title="RVDB C-RVDB（大库，耗时长）">
            🧬 构建 RVDB 库</button>
        </div>""",
     """        <div style="display:flex;gap:10px;margin-top:10px;flex-wrap:wrap">
          <button class="btn" onclick="buildUniversalDb('refvirus', this)"
                  title="NCBI RefSeq Viral（推荐先建这个，速度快）">
            🧬 构建 RefSeq 病毒库</button>
          <button class="btn" onclick="buildUniversalDb('rvdb', this)"
                  title="RVDB C-RVDB（大库，耗时长）">
            🧬 构建 RVDB 库</button>
        </div>
        <div style="display:flex;gap:10px;margin-top:10px;flex-wrap:wrap;align-items:end">
          <div class="gp-item">
            <label>Kraken2 库包路径（.tar.gz，官方预构建库）</label>
            <div class="filerow">
              <input id="k2Tar" type="text" value="databases/k2_viral_20260626.tar.gz">
              <button class="btn small" onclick="browse('k2Tar')">📁</button>
            </div>
          </div>
          <div class="gp-item">
            <label>hash capacity</label>
            <input id="k2Cap" type="text" value="1G" style="width:80px">
          </div>
          <button class="btn" onclick="convertKraken2(this)"
                  title="kunpeng hashshard：Kraken2 索引一次性转为 kunpeng 分片库">
            🔄 从 Kraken2 库转换</button>
        </div>
        <div id="k2Stat" class="dbbadge" style="margin-top:6px"></div>"""),
    ("""    if ($('rvdbStat')) $('rvdbStat').innerHTML =
      badge(!!(d.rvdb && d.rvdb.ready), 'RVDB 库');""",
     """    if ($('rvdbStat')) $('rvdbStat').innerHTML =
      badge(!!(d.rvdb && d.rvdb.ready), 'RVDB 库');
    if ($('k2Stat')) $('k2Stat').innerHTML =
      badge(!!(d.k2viral && d.k2viral.ready), 'Kraken2 转换库 (k2viral)');"""),
    ("""/* 通用病毒库构建 */
async function buildUniversalDb(source, btn) {""",
     """/* Kraken2 库转换（kunpeng hashshard） */
async function convertKraken2(btn) {
  const tar = val('k2Tar');
  if (!tar) { alert('请填写 Kraken2 库包路径'); return; }
  try {
    const r = await fetch('/api/convert_kraken2', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ tar, name: 'k2viral_db',
                             hash_capacity: val('k2Cap') || '1G' })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    taskLogOpen.add(d.task);
    startPolling();
    watchTaskBtn(d.task, btn, '⏳ 转换中…', () => loadDbs());
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

/* 通用病毒库构建 */
async function buildUniversalDb(source, btn) {"""),
])
print('K2 PATCH ALL OK')

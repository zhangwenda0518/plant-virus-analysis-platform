/* 植物病毒分析平台 GUI 逻辑（原生 JS，无依赖；多语言文案来自 i18n.js） */
'use strict';

const $ = (id) => document.getElementById(id);

function esc(s) {
  return String(s ?? '').replace(/[&<>\"']/g,
    c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}

// 输入框取值（trim）；meta/download/submit 等页的历史页面依赖此全局函数
function val(id) {
  const el = $(id);
  return el && el.value != null ? String(el.value).trim() : '';
}

// 样品创建统一入口（P2-12）：pipeline 页与 samples 页共用，避免重复 fetch 样板。
// payload 字段与 /api/pipeline/create 契约一致：sample/r1/r2/project/subsample
async function apiCreateSample(payload) {
  try {
    const r = await fetch('/api/pipeline/create', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(payload)});
    if (!r.ok) {
      let msg = t('c.fail', '创建失败');
      try { msg = (await r.json()).error || msg; } catch (e) {}
      return { ok: false, error: msg };
    }
    return { ok: true, data: await r.json() };
  } catch (e) {
    return { ok: false, error: t('c.connFail', '无法连接平台服务') + ': ' + e, conn: true };
  }
}

// ---------------- 按钮防重复点击 / 任务看护 ----------------
// lockBtn(btn, '⏳ 检索中…') → 返回解锁函数；按钮置灰防二次点击。
function lockBtn(btn, busyText) {
  if (!btn) return () => {};
  if (btn.dataset.origHtml == null) btn.dataset.origHtml = btn.innerHTML;
  btn.disabled = true;
  btn.classList.add('btn-busy');
  btn.innerHTML = busyText || '⏳ 运行中…';
  let done = false;
  return function unlock() {
    if (done) return;
    done = true;
    btn.disabled = false;
    btn.classList.remove('btn-busy');
    btn.innerHTML = btn.dataset.origHtml;
  };
}

// 瞬时请求包装：请求期间锁定按钮（fetch 顺序操作也适用）
async function withBtn(btn, fn, busyText) {
  if (btn && btn.disabled) return undefined;   // 已在执行，忽略本次点击
  const unlock = lockBtn(btn, busyText);
  try {
    return await fn(unlock);
  } finally {
    unlock();
  }
}

// 任务看护：锁定按钮直到后台任务进入终态（done/failed/cancelled）。
// 返回终态任务快照；onEnd(t) 在解锁前回调（刷新列表等）。
// opts.cancelable: 锁定期间在按钮旁显示「⏹ 停止」，终态后自动移除。
async function watchTaskBtn(tid, btn, busyText, onEnd, opts) {
  const unlock = lockBtn(btn, busyText);
  let stopBtn = null;
  if (opts && opts.cancelable) {
    stopBtn = document.createElement('button');
    stopBtn.className = 'btn small danger';
    stopBtn.textContent = '⏹ 停止';
    stopBtn.onclick = async () => {
      if (!confirm(t('tk.confirmStop', '停止该任务？'))) return;
      stopBtn.disabled = true; stopBtn.textContent = '停止中…';
      try { await fetch('/api/task/' + tid + '/cancel', {method: 'POST'}); } catch (e) {}
      stopBtn.disabled = false; stopBtn.textContent = '⏹ 停止';
    };
    btn.insertAdjacentElement('afterend', stopBtn);
  }
  const finish = (snap) => {
    if (stopBtn) stopBtn.remove();
    unlock();
    if (onEnd) onEnd(snap);
    return snap;
  };
  let misses = 0;              // 任务快照不可达计数（404/网络）
  for (;;) {
    await new Promise(r => setTimeout(r, 3000));
    let snap = null, httpOk = true;
    try {
      const r = await fetch('/api/task/' + tid);
      httpOk = r.ok;
      snap = await r.json();
    } catch (e) { httpOk = false; }
    if (!httpOk) {
      // 服务重启后内存任务丢失（404）或网络故障：最多重试 10 次（30s）后
      // 解锁按钮，避免永久卡在「运行中」。
      if (++misses >= 10) return finish(null);
      continue;
    }
    misses = 0;
    if (snap && ['done', 'failed', 'cancelled'].includes(snap.status)) {
      return finish(snap);
    }
  }
}

function fmtSize(n) {
  if (n == null) return '';
  if (n > 1e9) return (n / 1e9).toFixed(2) + ' GB';
  if (n > 1e6) return (n / 1e6).toFixed(1) + ' MB';
  if (n > 1e3) return (n / 1e3).toFixed(0) + ' KB';
  return n + ' B';
}

// 秒 → 人读时长
function fmtDur(sec) {
  if (sec == null || sec < 0) return '';
  const s = Math.round(sec);
  if (s >= 3600) return `${Math.floor(s / 3600)} h ${String(Math.floor((s % 3600) / 60)).padStart(2, '0')} min`;
  if (s >= 60) return `${Math.floor(s / 60)} min ${String(s % 60).padStart(2, '0')} s`;
  return `${s} s`;
}

// ---------------- Toast 通知 ----------------
function toast(title, body, opts) {
  let box = $('toastBox');
  if (!box) {
    box = document.createElement('div');
    box.id = 'toastBox';
    document.body.appendChild(box);
  }
  const el = document.createElement('div');
  el.className = 'toast' + ((opts && opts.kind) ? ' ' + opts.kind : '');
  el.innerHTML = `<div class="t-title">${esc(title)}</div>` +
    (body ? `<div class="t-body">${esc(body)}</div>` : '') +
    ((opts && opts.actions)
      ? `<div class="t-actions">${opts.actions.map((a, i) =>
          `<button class="btn small ${a.primary ? 'primary' : ''}" data-ta="${i}">${esc(a.label)}</button>`).join('')}</div>`
      : '');
  if (opts && opts.actions) {
    el.querySelectorAll('[data-ta]').forEach(b => {
      b.onclick = () => { opts.actions[+b.dataset.ta].onClick && opts.actions[+b.dataset.ta].onClick(); el.remove(); };
    });
  }
  box.appendChild(el);
  const ttl = (opts && opts.ttl) != null ? opts.ttl : 8000;
  if (ttl > 0) setTimeout(() => el.remove(), ttl);
  return el;
}

/* 后台完成时发系统通知（页面不可见才打扰；需用户授权） */
function sysNotify(title, body) {
  try {
    if (document.visibilityState === 'visible' || !window.Notification) return;
    if (Notification.permission === 'granted') new Notification(title, {body});
    else if (Notification.permission === 'default') Notification.requestPermission();
  } catch (e) { /* 通知不可用忽略 */ }
}

// ---------------- 服务断连检测 ----------------
let connFails = 0;
function setConnBanner(down) {
  let b = $('conn-banner');
  if (down) {
    if (!b) {
      b = document.createElement('div');
      b.id = 'conn-banner';
      b.style.cssText = 'position:sticky;top:0;z-index:99;background:#c62828;color:#fff;'
        + 'padding:10px 16px;font-size:14px;line-height:1.6;';
      document.body.prepend(b);
    }
    b.innerHTML = '<b>⚠ ' + t('conn.lost', '与平台服务的连接已断开') + '</b><br>'
      + t('conn.lostHint',
          '1. 黑色控制台窗口是否仍在运行（关闭该窗口 = 退出平台）；<br>'
        + '2. 若已关闭：重新双击 VirusPlatform.exe（或 启动平台.bat）；<br>'
        + '3. 大任务占满内存/CPU 时服务可能暂时无响应，任务结束后本提示自动消失。');
    startConnProbe();
  } else if (b) {
    b.remove();
  }
}

// 断连横幅自愈探针：每 5s 探测一次，服务恢复即自动摘掉横幅
let _connProbe = null;
function startConnProbe() {
  if (_connProbe) return;
  _connProbe = setInterval(async () => {
    try {
      const r = await fetch('/api/tasks');
      if (r.ok) { clearInterval(_connProbe); _connProbe = null; setConnBanner(false); }
    } catch (e) { /* 仍未恢复，继续探测 */ }
  }, 5000);
}

// 日志框跟随尾部（运行中任务日志持续追加时自动滚到底）
function autoscrollLogs() {
  document.querySelectorAll('pre[data-autoscroll="1"]').forEach(el => {
    el.scrollTop = el.scrollHeight;
  });
}

// ---------------- 粘贴序列输入 ----------------
/* pasteSeq(inputId, ext)：弹出粘贴对话框，确认后写入 uploads/ 并回填路径。
   ext: '.fasta' | '.fastq' | '.tsv' */
function pasteSeq(inputId, ext) {
  ext = ext || '.fasta';
  var mask = document.createElement('div');
  mask.className = 'dlgmask';
  mask.style.zIndex = 300;
  var isFq = ext === '.fastq' || ext === '.fq';
  var fmtHint = isFq ? 'FASTQ（@ 开头）' : 'FASTA（> 开头）';
  mask.innerHTML = '<div class="dlg" style="width:640px">' +
    '<div class="dlghead"><b>📋 粘贴' + fmtHint + ' 序列</b>' +
    '<button class="btn small" onclick="this.closest(\'.dlgmask\').remove()">✕</button></div>' +
    '<textarea id="pasteTA" rows="14" style="width:calc(100% - 32px);margin:10px 16px;' +
    'font-family:monospace;font-size:11.5px;border:1px solid #ddd;border-radius:4px;' +
    'padding:8px;resize:vertical" placeholder="粘贴 ' + fmtHint + ' 序列文本…"></textarea>' +
    '<div style="padding:0 16px 12px;text-align:right">' +
    '<button class="btn small" onclick="this.closest(\'.dlgmask\').remove()">取消</button> ' +
    '<button class="btn small primary" id="pasteOk">✓ 确认粘贴</button></div></div>';
  document.body.appendChild(mask);
  mask.querySelector('#pasteTA').focus();
  mask.querySelector('#pasteOk').addEventListener('click', function() {
    var text = mask.querySelector('#pasteTA').value.trim();
    if (!text) { alert('请粘贴序列内容'); return; }
    var btn = this;
    btn.disabled = true; btn.textContent = '写入中…';
    fetch('/api/paste_input', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ text: text, ext: ext })
    }).then(function(r) { return r.json(); }).then(function(d) {
      if (d.path) {
        var inp = $(inputId);
        if (inp) inp.value = d.path;
        mask.remove();
      } else { alert(d.error || '写入失败'); btn.disabled = false; btn.textContent = '✓ 确认粘贴'; }
    }).catch(function(e) { alert('无法连接: ' + e); btn.disabled = false; btn.textContent = '✓ 确认粘贴'; });
  });
}

// ---------------- 内置示例数据 ----------------
// 平台自带示例（两条植物病毒完整基因组，~6-7kb），各工具「✨ 示例」按钮共用
const EXAMPLE_FASTA = 'databases/examples/example_viral_contigs.fasta';
// 6 条同属近缘基因组（3 参考 + 3 受控突变衍生株）：MSA / 结构比较 / SDT 等多序列示例
const EXAMPLE_SET_FASTA = 'databases/examples/example_virus_set.fasta';
// 5 条同种近缘序列（≥90% 一致）：保守区引物设计示例（需全表保守区段）
const EXAMPLE_CONSERVED_FASTA = 'databases/examples/example_conserved_set.fasta';
// 示例树（上集建树产物）：进化树查看器「✨ 示例」
const EXAMPLE_TREE_NWK = 'databases/examples/example_tree.nwk';
// 示例 GenBank（含 CDS 注释）：基因组图谱 GenBank 模式「✨ 示例」
const EXAMPLE_GENBANK_GB = 'databases/examples/example_genome.gb';
// 共线性比较离线示例（3 条同属小基因组 .gb，逗号分隔供 s_files 导入）
const EXAMPLE_SYNTENY_GBS = 'databases/examples/example_synteny_A.gb,databases/examples/example_synteny_B.gb,databases/examples/example_synteny_C.gb';
function fillExample(inputId, path) {
  var inp = $(inputId);
  if (!inp) return;
  inp.value = path || EXAMPLE_FASTA;
  if (typeof toast === 'function') {
    toast('已填入内置示例', inp.value, { ttl: 3000 });
  }
}

/* 「参考序列获取」卡「✨ 示例」：预填递进选择（界→科→属）+ 集合名 */
function seqPrepExample() {
  var name = $('s_name');
  if (name) name.value = 'dianthovirus_ictv';
  ['Realm', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'].forEach(c => {
    delete document.body.dataset['ictv_' + c];
  });
  document.body.dataset['ictv_Realm'] = 'Riboviria';
  document.body.dataset['ictv_Family'] = 'Tombusviridae';
  document.body.dataset['ictv_Genus'] = 'Dianthovirus';
  ictvCascadeRefetch().then(function () {
    if (typeof toast === 'function')
      toast('已填入示例：界→科→属 递进选择', 'Dianthovirus（6 条），点「🔍 预览」离线查看', { ttl: 3500 });
  });
}

const ICTV_RANK_MAIN = ['Realm', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'];
let ictvRankMeta = [];      // [{col, zh}]（含 Subfamily 等附加级）

/* 级联选参：任意一级变更 → 清空更深等级 → 重新拉取各级选项 */
async function ictvCascadeRefetch() {
  const box = $('spCascade');
  if (!box) return;
  const qs = ICTV_RANK_MAIN.filter(c => _v('sp_c_' + c))
    .map(c => `${c}=${encodeURIComponent(_v('sp_c_' + c))}`).join('&');
  try {
    const r = await fetch('/api/ictv/cascade' + (qs ? '?' + qs : ''));
    if (!r.ok) { box.innerHTML = '<p class="hint" style="color:#b91c1c">' + esc((await r.json()).error || '加载失败') + '</p>'; return; }
    const d = await r.json();
    ictvRankMeta = d.ranks.map(x => ({ col: x.col, zh: x.zh }));
    const deepestSel = ICTV_RANK_MAIN.filter(c => _v('sp_c_' + c)).pop();
    $('spScope').textContent = deepestSel
      ? `当前最深选择：${deepestSel} = ${_v('sp_c_' + deepestSel)}`
      : '（尚未选择，先从「界」开始逐级下钻）';
    box.innerHTML = ICTV_RANK_MAIN.map(col => {
      const meta = d.ranks.find(x => x.col === col);
      const opts = (meta && meta.options) || [];
      return `<div class="gp-item"><label>${esc(col === 'Realm' ? '界 Realm' : meta ? meta.zh + ' ' + col : col)}</label>` +
        `<select id="sp_c_${col}" onchange="ictvCascadeChanged('${col}')">` +
        `<option value="">（全部）</option>` +
        opts.slice(0, 600).map(o => `<option value="${esc(o.name)}">${esc(o.name)}（${o.n}）</option>`).join('') +
        `</select></div>`;
    }).join('');
    // 回填已选值（后端按已选过滤返回，值仍在选项里）
    ICTV_RANK_MAIN.forEach(c => {
      const sel = $('sp_c_' + c);
      const v = document.body.dataset['ictv_' + c] || '';
      if (sel && v && [...sel.options].some(o => o.value === v)) sel.value = v;
    });
  } catch (e) { box.innerHTML = '<p class="hint" style="color:#b91c1c">无法连接: ' + esc(e) + '</p>'; }
}

function ictvCascadeChanged(col) {
  // 该级变更 → 更深等级全部清空（dataset 记录当前选择）
  const idx = ICTV_RANK_MAIN.indexOf(col);
  ICTV_RANK_MAIN.slice(idx).forEach(c => { delete document.body.dataset['ictv_' + c]; });
  const v = _v('sp_c_' + col);
  if (v) document.body.dataset['ictv_' + col] = v;
  ictvCascadeRefetch();
}

function ictvSelectedLevels() {
  const out = {};
  ICTV_RANK_MAIN.forEach(c => {
    const v = _v('sp_c_' + c) || document.body.dataset['ictv_' + c] || '';
    if (v) out[c] = v;
  });
  return out;
}

async function ictvPreview() {
  const levels = ictvSelectedLevels();
  if (!Object.keys(levels).length) { alert('请至少选择一个分类等级'); return; }
  const box = $('spPreview');
  box.innerHTML = '<p class="hint">选参中…</p>';
  try {
    const r = await fetch('/api/ictv/preview', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(Object.assign({}, levels, {
        genome: _v('sp_genome') || 'complete',
        limit: +_v('sp_limit') || 200,
        per_genus: +(_v('sp_per_genus') || 0),
        per_species: +(_v('sp_per_species') || 0) }))});
    if (!r.ok) { box.innerHTML = '<p class="hint" style="color:#b91c1c">选参失败: ' + esc((await r.json()).error || '') + '</p>'; return; }
    const d = await r.json();
    box.innerHTML = `<p class="hint"><b>${esc(d.scope)}</b> 命中 <b>${d.total}</b> 条 · 抽样 <b>${esc(d.sampling)}</b>`
      + `（${d.total > d.n_shown ? '按本地优先显示前 ' + d.n_shown : '全部显示'}）</p>`
      + '<table class="table" style="width:100%;border-collapse:collapse">'
      + '<tr><th>Accession</th><th>来源</th><th>完整度</th><th>物种</th><th>基因组</th><th>宿主</th></tr>'
      + d.rows.map(x => `<tr><td>${esc(x.acc)}</td>`
        + `<td>${x.source === 'ncbi' ? '🌐 ncbi（需下载）' : '💾 ' + esc(x.source)}</td>`
        + `<td>${esc(x.coverage || '—')}</td><td>${esc(x.species)}</td>`
        + `<td>${esc(x.genome || '—')}</td><td>${esc(x.host || '—')}</td></tr>`).join('')
      + '</table>';
  } catch (e) { box.innerHTML = '<p class="hint" style="color:#b91c1c">无法连接: ' + esc(e) + '</p>'; }
}

async function ictvDownload(btn) {
  const levels = ictvSelectedLevels();
  const coll = _v('s_name');
  if (!Object.keys(levels).length) { alert('请至少选择一个分类等级'); return; }
  if (!coll) { alert('请填集合名'); return; }
  try {
    const r = await fetch('/api/ictv/download', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(Object.assign({}, levels, {
        coll, genome: _v('sp_genome') || 'complete',
        limit: +_v('sp_limit') || 200,
        per_genus: +(_v('sp_per_genus') || 0),
        per_species: +(_v('sp_per_species') || 0) }))});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    taskLogOpen.add(d.task);
    startPolling();
    _watchTask(d.task, () => { loadGbCollections(); loadNcbiCollections(); loadTbColls(); },
               btn, '⏳ ICTV 序列下载中…');
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

/* ================= 序列比对模块：运行 / 彩色查看器 / 编辑保存 ================= */
let alData = null, alPageNo = 0, alEditing = false;

const AL_NT_COLORS = { A: '#2e7d32', T: '#c62828', U: '#c62828', G: '#e65100',
                       C: '#1565c0', N: '#78909c', '-': '#eceff1', '?': '#eceff1' };
const AL_AA_COLORS = (() => {
  const groups = {
    hydrophobic: 'AVLIMFWC', positive: 'KRH', negative: 'DE',
    polar: 'STNQ', special: 'PG', stop: '*' };
  const colors = { hydrophobic: '#90caf9', positive: '#ef9a9a',
                   negative: '#a5d6a7', polar: '#ffe082', special: '#ce93d8',
                   stop: '#b0bec5' };
  const map = { '-': '#eceff1', X: '#cfd8dc', U: '#cfd8dc', B: '#cfd8dc',
                Z: '#cfd8dc', '*': '#b0bec5' };
  for (const g in groups)
    for (const ch of groups[g]) map[ch] = colors[g];
  return map;
})();
function alColor(ch, type) {
  const m = type === 'aa' ? AL_AA_COLORS : AL_NT_COLORS;
  return (m[ch] || '#cfd8dc');
}

async function alignRun(btn) {
  const seqs = _v('al_fa');
  if (!seqs) { alert('请选择或粘贴 FASTA（≥2 条序列）'); return; }
  try {
    const r = await fetch('/api/tool/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ tool: 'align',
                             params: { seqs,
                                       strategy: _v('al_strategy') || 'auto',
                                       trimal: _v('al_trimal') || 'automated1',
                                       max_n: +_v('al_maxn') || 100 } })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    taskLogOpen.add(d.task);
    startPolling();
    _watchTask(d.task, () => {
      // 完成后自动打开清剪结果（无清剪则原比对）
      const res = (TMlastResult(d.task) || {});
      const p = res.aln_used === 'aln.trim.fasta'
        ? `tool_runs/${d.run}/aln.trim.fasta` : `tool_runs/${d.run}/aln.fasta`;
      $('alPath').value = p;
      alignLoad();
    }, btn, '⏳ 比对中…');
  } catch (e) { alert('无法连接平台服务: ' + e); }
}
function TMlastResult(tid) {
  // 任务面板轮询的数据不在 JS 侧持久化：退化为直接请求比对文件即可
  return null;
}

async function alignLoad() {
  const p = _v('alPath');
  if (!p) { alert('请填比对 FASTA 路径'); return; }
  const box = $('alBox');
  box.innerHTML = '<p class="hint">加载中…</p>';
  try {
    const r = await fetch('/api/align/file?path=' + encodeURIComponent(p));
    if (!r.ok) { box.innerHTML = '<p class="hint" style="color:#b91c1c">' + esc((await r.json()).error || '加载失败') + '</p>'; return; }
    alData = await r.json();
    alPageNo = 0;
    alignRender();
  } catch (e) { box.innerHTML = '<p class="hint" style="color:#b91c1c">无法连接: ' + esc(e) + '</p>'; }
}

function alignRender() {
  const box = $('alBox');
  if (!alData) { box.innerHTML = '<p class="hint">暂无数据</p>'; $('alNav').style.display = 'none'; return; }
  if (!alData.aligned)
    $('alMeta').textContent = '⚠ 序列长度不一致（未比对或 FASTA）：彩色视图按最长序列展示，建议先运行比对';
  else
    $('alMeta').textContent = '';
  const cols = +($('alCols') && $('alCols').value) || 100;
  const total = Math.max(1, Math.ceil(alData.cols / cols));
  alPageNo = Math.min(Math.max(0, alPageNo), total - 1);
  const start = alPageNo * cols;
  const shown = alData.seqs.map(s => (s + '-'.repeat(alData.cols)).slice(start, start + cols));
  const nameW = Math.max(...alData.names.map(n => n.length), 8);
  let html = '<div style="font-family:Consolas,monospace;font-size:12.5px;white-space:pre;line-height:1.55">';
  html += '<div style="color:#64748b">' + ' '.repeat(nameW) +
    ' ' + String(start + 1).padStart(cols / 2) + '</div>';
  alData.names.forEach((name, i) => {
    let line = '<span style="color:#334155;font-weight:600">' +
      esc(name.padEnd(nameW)) + '</span> ';
    for (const ch of shown[i])
      line += `<span style="color:${alColor(ch, alData.type)};font-weight:700">${ch === '-' ? '-' : esc(ch)}</span>`;
    html += '<div>' + line + '</div>';
  });
  html += '</div>';
  box.innerHTML = html;
  $('alNav').style.display = '';
  $('alPageInfo').textContent = `列 ${start + 1}–${Math.min(start + cols, alData.cols)} / ${alData.cols}`;
  // 图例
  const key = alData.type === 'aa'
    ? ['AVLIMFWC|疏水', 'KRH|正电', 'DE|负电', 'STNQ|极性', 'PG|特殊', '-|gap']
    : ['A', 'T/U', 'G', 'C', 'N', '-'].map(x => x + '|' + ({ A: 'A', 'T/U': 'T/U', G: 'G', C: 'C', N: 'N', '-': 'gap' }[x]));
  $('alLegend').innerHTML = key.map(kv => {
    const [chars, label] = kv.split('|');
    const c = alColor(chars[0], alData.type);
    return `<span style="margin-right:12px;font-size:12px"><span style="display:inline-block;width:12px;height:12px;background:${c};vertical-align:-2px;margin-right:3px;border-radius:2px"></span>${esc(label)}</span>`;
  }).join('');
  const dl = $('alDl');
  dl.style.display = '';
  dl.href = '/tool_runs/' + alData.path.replace(/^.*tool_runs[\/\\]/, '');
  dl.download = alData.path.split(/[\\/]/).pop();
  dl.textContent = t('tk.dlRef', '⬇ 下载') + ' ' + alData.path.split(/[\\/]/).pop();
}

function alignPage(d) {
  if (!alData) return;
  const cols = +($('alCols') && $('alCols').value) || 100;
  const total = Math.max(1, Math.ceil(alData.cols / cols));
  alPageNo = Math.min(total - 1, Math.max(0, alPageNo + d));
  alignRender();
}

function alignEditToggle() {
  alEditing = !alEditing;
  const wrap = $('alEditWrap');
  wrap.style.display = alEditing ? '' : 'none';
  $('alEditBtn').textContent = alEditing ? '👁 退出编辑' : '✏️ 编辑模式';
  if (alEditing && alData) {
    // FASTA 化（70 列换行）
    let out = '';
    alData.names.forEach((n, i) => {
      out += '>' + n + '\n';
      const s = alData.seqs[i];
      for (let j = 0; j < s.length; j += 70) out += s.slice(j, j + 70) + '\n';
    });
    $('alEdit').value = out;
    $('alEditMsg').textContent = '';
  }
}

async function alignSave(btn) {
  const content = $('alEdit').value;
  if (!content.trim()) { alert('编辑区为空'); return; }
  btn.disabled = true;
  try {
    const r = await fetch('/api/align/save', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ path: alData.path, content }) });
    const d = await r.json();
    if (!r.ok) { $('alEditMsg').textContent = '❌ ' + (d.error || '保存失败'); return; }
    $('alEditMsg').textContent = `✅ 已保存 ${d.n} 条${d.aligned ? '（等长比对）' : '（⚠ 长度不一致——未比对状态）'}：${d.saved}`;
    // 打开保存后的副本继续查看
    $('alPath').value = d.saved;
    alignLoad();
  } catch (e) { $('alEditMsg').textContent = '无法连接: ' + e; }
  btn.disabled = false;
}

function alignSend(target) {
  if (!alData) { alert('请先打开要比对的文件'); return; }
  const p = alData.path;
  if (target === 'treebuild') {
    $('qt_fa').value = p;
    location.hash = '#t-treebuild';
  } else {
    $('sd_fa').value = p;
    $('sd_aligned').checked = alData.aligned;
    location.hash = '#t-sdt';
  }
  if (typeof renderModuleTree === 'function') renderModuleTree();
}

function tbMolChanged() {
  const mol = ($('tbMolecule') && $('tbMolecule').value) || 'genome';
  const gene = $('tbGene');
  if (gene) gene.disabled = mol === 'genome';
}

/* ================= 进化树构建：集合下拉 + 建树 ================= */
async function loadTbColls() {
  const sel = $('tbColl');
  if (!sel) return;
  try {
    const cols = await (await fetch('/api/gb/collections')).json();
    sel.innerHTML = cols.length
      ? cols.map(c => `<option value="${esc(c.name)}">${esc(c.name)}（${c.n_records} 条记录）</option>`).join('')
      : '<option value="">（暂无集合——先到「参考序列获取」下载）</option>';
  } catch (e) { sel.innerHTML = '<option value="">（无法连接）</option>'; }
}

async function tbBuild(btn) {
  const coll = _v('tbColl');
  if (!coll) { alert('请先选择 GenBank 集合（到「参考序列获取」卡下载）'); return; }
  await gbBuildTree(coll, btn);
  loadGbCollections();
}

/* 通用文本复制（引物序列等小段文本） */
async function copyText(text, tag) {
  try {
    await navigator.clipboard.writeText(text);
  } catch (e) {
    const ta = document.createElement('textarea');
    ta.value = text; document.body.appendChild(ta);
    ta.select(); document.execCommand('copy'); ta.remove();
  }
  if (typeof toast === 'function') toast('已复制', tag || '', { ttl: 2500 });
}

/* Metabuli 同款示例病毒（TMV / PVY / CMV / PSTVd / Mix）：
   来源 <DEMO-IP>/metabuli examples，已固化到 databases/examples/ */
const EXAMPLE_METABULI = {
  tmv: 'databases/examples/example_tmv.fasta',
  pvy: 'databases/examples/example_pvy.fasta',
  cmv: 'databases/examples/example_cmv.fasta',
  pstvd: 'databases/examples/example_pstvd.fasta',
  mix: 'databases/examples/example_mix.fasta'
};
function vEx(key, inputId) {
  const path = EXAMPLE_METABULI[key];
  if (!path) return;
  const inp = $(inputId || 'c_fa');
  if (inp) inp.value = path;
  if (key === 'pstvd') {                 // 359nt 类病毒：低于默认 500bp 过滤线
    const ml = $('c_minlen');
    if (ml && +ml.value > 300) { ml.value = 200; }
  }
  const names = { tmv: 'TMV 烟草花叶病毒', pvy: 'PVY 马铃薯Y病毒',
                  cmv: 'CMV 黄瓜花叶病毒(RNA1-3)', pstvd: 'PSTVd 马铃薯纺锤块茎类病毒',
                  mix: 'Mix 四病毒混合(6条)' };
  if (typeof toast === 'function')
    toast('已填入示例病毒', names[key] + ' · ' + path, { ttl: 3500 });
}

// ---------------- 表格分页（每页 20 行） ----------------
const PER_PAGE = 20;
// fnExpr 里可用 {p} 占位目标页码，如 "metaGo" 或 "dlFileGo('abc',{p})"
function pagerHtml(page, totalPages, fnExpr) {
  totalPages = Math.max(1, Math.floor(totalPages || 1));
  if (totalPages <= 1) return '';
  page = Math.min(Math.max(1, page || 1), totalPages);
  const call = p => fnExpr.includes('{p}')
    ? fnExpr.replace('{p}', p) : `${fnExpr}(${p})`;
  const nums = [];
  for (let n = 1; n <= totalPages; n++) {
    if (n === 1 || n === totalPages || Math.abs(n - page) <= 2) nums.push(n);
  }
  let html = `<div class="pager" style="display:flex;gap:4px;align-items:center;flex-wrap:wrap;margin:8px 0">`
    + `<button class="btn small" ${page <= 1 ? 'disabled' : ''} onclick="${call(page - 1)}">◀</button>`;
  let prev = 0;
  for (const n of nums) {
    if (n - prev > 1) html += '<span class="hint">…</span>';
    html += `<button class="btn small${n === page ? ' primary' : ''}" onclick="${call(n)}">${n}</button>`;
    prev = n;
  }
  html += `<button class="btn small" ${page >= totalPages ? 'disabled' : ''} onclick="${call(page + 1)}">▶</button>`
    + `<span class="hint" style="margin-left:4px">第 ${page} / ${totalPages} 页</span></div>`;
  return html;
}

// ---------------- 模块历史运行（折叠 · 上下游衔接 · 删除） ----------------
/* 每个工具卡/独立模块页的持久历史：数据来自 /api/tool/runs（tool_runs/ 目录，
   重启不丢）。条目可展开文件下载、一键复用为下游模块输入（跨页暂存回填）、
   可删除。默认折叠，避免过多。 */
const RH_PAGES = {
  f_r1: '/tools#t-fastp', f_r2: '/tools#t-fastp',
  hr_r1: '/hostremoval', hr_r2: '/hostremoval',
  i_input: '/tools#t-identify', i_input2: '/tools#t-identify',
  a_r1: '/tools#t-assemble', a_r2: '/tools#t-assemble',
  c_fa: '/tools#t-contigs',
  hp_tsv: '/hostpredict', hp_fa: '/hostpredict',
  of_fa: '/orf', oa_fa: '/annotation', oa_run: '/annotation',
  gp_fa: '/genome', gp_ann: '/genome', pr_fa: '/primer',
  al_fa: '/tools#t-align', qt_fa: '/tools#t-treebuild',
  tv_file: '/tools#t-treebuild', sd_fa: '/tools#t-sdt'
};

function rhPick(files, sub, excl) {
  return ((files || []).filter(f => (f.path || '').includes(sub)
                            && !(excl && (f.path || '').includes(excl)))[0] || {}).path;
}
function rhFqPair(files) {
  const fqs = (files || []).map(f => f.path).filter(p => /\.fastq\.gz$/.test(p));
  if (!fqs.length) return null;
  const num = p => { const b = p.replace(/\.fastq\.gz$/, '').split(/[\/]/).pop();
                     return b.endsWith('_1') || b.endsWith('.1') || b.endsWith('_R1') ? 1
                          : b.endsWith('_2') || b.endsWith('.2') || b.endsWith('_R2') ? 2 : 0; };
  const r1 = fqs.find(p => num(p) === 1), r2 = fqs.find(p => num(p) === 2);
  if (r1 && r2) return [r1, r2];
  return [fqs[0], null];                       // 单端
}
const RH_DEFS = [
  { box: 'rh-convert', prefixes: ['convert_'],
    chains: files => {
      const out = [];
      const fa = rhPick(files, '.fa.gz') || rhPick(files, '.fasta.gz');
      const fq = rhFqPair(files);
      if (fq) out.push(
        { label: '→ 质控', page: RH_PAGES.f_r1,
          fills: { f_r1: fq[0], ...(fq[1] ? { f_r2: fq[1] } : {}) } },
        { label: '→ 宿主去除', page: RH_PAGES.hr_r1,
          fills: { hr_r1: fq[0], ...(fq[1] ? { hr_r2: fq[1] } : {}) } },
        { label: '→ 病毒识别', page: RH_PAGES.i_input,
          fills: { i_input: fq[0], ...(fq[1] ? { i_input2: fq[1] } : {}) } });
      if (fa) out.push(
        { label: '→ 识别(FASTA)', page: RH_PAGES.i_input,
          fills: { i_input: fa, i_type: 'fasta' } },
        { label: '→ ORF', page: RH_PAGES.of_fa, fills: { of_fa: fa } },
        { label: '→ 比对', page: RH_PAGES.al_fa, fills: { al_fa: fa } },
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: fa } },
        { label: '→ SDT', page: RH_PAGES.sd_fa, fills: { sd_fa: fa } });
      return out;
    } },
  { box: 'rh-fastp', prefixes: ['fastp_'],
    chains: files => {
      const p = f => rhPick(files, f);
      const r1 = p('fastp_R1'), r2 = p('fastp_R2');
      if (!r1) return [];
      return [
        { label: '→ 宿主去除', page: RH_PAGES.hr_r1,
          fills: { hr_r1: r1, hr_r2: r2 || '' } },
        { label: '→ 病毒识别', page: RH_PAGES.i_input,
          fills: { i_input: r1, i_input2: r2 || '' } },
        { label: '→ 组装', page: RH_PAGES.a_r1,
          fills: { a_r1: r1, a_r2: r2 || '' } }];
    } },
  { box: 'rh-hostremoval', prefixes: ['hostremoval_'],
    chains: files => {
      const r1 = rhPick(files, 'kept_R1'), r2 = rhPick(files, 'kept_R2');
      if (!r1) return [];
      return [
        { label: '→ 病毒识别', page: RH_PAGES.i_input,
          fills: { i_input: r1, i_input2: r2 || '' } },
        { label: '→ 组装', page: RH_PAGES.a_r1,
          fills: { a_r1: r1, a_r2: r2 || '' } },
        { label: '→ 质控', page: RH_PAGES.f_r1,
          fills: { f_r1: r1, f_r2: r2 || '' } }];
    } },
  { box: 'rh-identify', prefixes: ['identify_'],
    chains: files => {
      const fa = rhPick(files, 'viral_sequences');
      if (!fa) return [];
      return [
        { label: '→ 组装', page: RH_PAGES.a_r1, fills: { a_r1: fa } },
        { label: '→ ORF', page: RH_PAGES.of_fa, fills: { of_fa: fa } },
        { label: '→ 比对', page: RH_PAGES.al_fa, fills: { al_fa: fa } },
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: fa } },
        { label: '→ SDT', page: RH_PAGES.sd_fa, fills: { sd_fa: fa } }];
    } },
  { box: 'rh-assemble', prefixes: ['assemble_'],
    chains: files => {
      const fa = rhPick(files, 'contigs.filtered.fasta');
      if (!fa) return [];
      return [
        { label: '→ contigs 分类', page: RH_PAGES.c_fa, fills: { c_fa: fa } },
        { label: '→ ORF', page: RH_PAGES.of_fa, fills: { of_fa: fa } },
        { label: '→ 比对', page: RH_PAGES.al_fa, fills: { al_fa: fa } },
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: fa } },
        { label: '→ SDT', page: RH_PAGES.sd_fa, fills: { sd_fa: fa } },
        { label: '→ 图谱', page: RH_PAGES.gp_fa, fills: { gp_fa: fa } },
        { label: '→ 引物', page: RH_PAGES.pr_fa, fills: { pr_fa: fa } }];
    } },
  { box: 'rh-contigs', prefixes: ['contigs_'],
    chains: files => {
      const fa = rhPick(files, 'viral_contigs.fasta');
      const tsv = rhPick(files, 'virus_classification.tsv');
      const out = [];
      if (fa) out.push(
        { label: '→ ORF', page: RH_PAGES.of_fa, fills: { of_fa: fa } },
        { label: '→ 功能注释', page: RH_PAGES.oa_fa, fills: { oa_fa: fa } },
        { label: '→ 比对', page: RH_PAGES.al_fa, fills: { al_fa: fa } },
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: fa } },
        { label: '→ SDT', page: RH_PAGES.sd_fa, fills: { sd_fa: fa } },
        { label: '→ 图谱', page: RH_PAGES.gp_fa, fills: { gp_fa: fa } },
        { label: '→ 引物', page: RH_PAGES.pr_fa, fills: { pr_fa: fa } });
      if (fa && tsv) out.push(
        { label: '→ 宿主预测', page: RH_PAGES.hp_tsv,
          fills: { hp_tsv: tsv, hp_fa: fa } });
      return out;
    } },
  { box: 'rh-orf', prefixes: ['orf_'],
    chains: (files, run) => {
      const out = [];
      const faa = rhPick(files, 'pyrodigal.faa');
      const nt = rhPick(files, 'pyrodigal.ffn');
      if (faa) out.push(
        { label: '→ 比对(蛋白)', page: RH_PAGES.al_fa, fills: { al_fa: faa } },
        { label: '→ SDT(AA)', page: RH_PAGES.sd_fa, fills: { sd_fa: faa } });
      if (nt) out.push(
        { label: '→ SDT(NT)', page: RH_PAGES.sd_fa, fills: { sd_fa: nt } });
      out.push({ label: '→ 功能注释(该运行)', page: '/annotation',
                 fills: { oa_run: run.name } });
      return out;
    } },
  { box: 'rh-orfa', prefixes: ['orfa_'],
    chains: (files, run) => {
      const gff = rhPick(files, 'orf_annotation.gff3');
      const fa = rhPick(files, 'viral_contigs.fasta');
      const out = [];
      if (gff && fa) out.push(
        { label: '→ 基因组图谱(注释)', page: '/genome',
          fills: { gp_ann: gff, gp_fa: fa } });
      const realRun = (run.name || '').replace(/^orfa_/, 'orf_');
      out.push({ label: '→ 功能注释(orf 运行)', page: '/annotation',
                 fills: { oa_run: realRun } });
      return out;
    } },
  { box: 'rh-genoplot', prefixes: ['genoplot_'], chains: () => [] },
  { box: 'rh-primer', prefixes: ['primer_'], chains: () => [] },
  { box: 'rh-align', prefixes: ['align_'],
    chains: files => {
      const aln = rhPick(files, 'aln.trim.fasta') || rhPick(files, 'aln.fasta');
      if (!aln) return [];
      return [
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: aln } },
        { label: '→ SDT(已比对)', page: RH_PAGES.sd_fa,
          fills: { sd_fa: aln, sd_aligned: true } }];
    } },
  { box: 'rh-treebuild', prefixes: ['structcmp_', 'quicktree_'],
    chains: files => {
      const out = [];
      const aln = rhPick(files, 'aln.fasta');
      const nwk = rhPick(files, 'nj.nwk') || rhPick(files, 'tree.nwk');
      if (aln) out.push(
        { label: '→ 建树', page: RH_PAGES.qt_fa, fills: { qt_fa: aln } },
        { label: '→ SDT(已比对)', page: RH_PAGES.sd_fa,
          fills: { sd_fa: aln, sd_aligned: true } });
      if (nwk) out.push(
        { label: '→ 树查看', page: RH_PAGES.tv_file, fills: { tv_file: nwk } });
      return out;
    } },
  { box: 'rh-sdt', prefixes: ['sdt_', 'identity_'],
    chains: files => {
      const out = [];
      const csv = rhPick(files, 'sdt_matrix.csv');
      if (csv) out.push({ label: '⬇ 矩阵 CSV', page: '', fills: {}, dl: csv });
      return out;
    } },
  { box: 'rh-hostpredict', prefixes: ['contigs_'],
    chains: files => {
      const tsv = rhPick(files, 'virus_classification.tsv');
      const fa = rhPick(files, 'viral_contigs.fasta');
      if (!tsv) return [];
      return [{ label: '→ 载入分类表', page: '/hostpredict',
                fills: { hp_tsv: tsv, ...(fa ? { hp_fa: fa } : {}) } }];
    } }
];

let RH_DATA = null;
const RH_OPEN = new Set();

async function rhRefresh() {
  try {
    const r = await fetch('/api/tool/runs');
    RH_DATA = await r.json();
  } catch (e) { RH_DATA = []; }
  for (const def of RH_DEFS) {
    const box = $(def.box);
    if (!box) continue;
    const runs = RH_DATA.filter(x => def.prefixes.some(p => x.name.startsWith(p)));
    const n = runs.length;
    const head = box.querySelector('.rh-head');
    if (head) head.querySelector('.rh-count').textContent = n;
    const body = box.querySelector('.rh-body');
    if (!body) continue;
    if (!n) { body.innerHTML = '<p class="hint">（暂无历史运行）</p>'; continue; }
    body.innerHTML = runs.map(run => {
      const chains = (def.chains || (() => []))(run.files || [], run) || [];
      const chainBtns = chains.map((c, i) => c.dl
        ? `<a class="btn small" href="/tool_runs/${encodeURIComponent(run.name)}/${c.dl}" download>${esc(c.label)}</a> `
        : `<button class="btn small" onclick="rhGo(this, '${esc(run.name)}', '${esc(def.box)}', ${i})">${esc(c.label)}</button> `).join('');
      const open = RH_OPEN.has(run.name);
      // 关键产物：只展示主结果/报告/序列（primary/green），其余收进“全部文件”
      const keyed = (run.files || []).map(f => {
        const base = String(f.path).split('/').pop();
        const { label, color } = dlLabel(base, f.path);
        return { f, base, label, color };
      });
      const keys = keyed.filter(x => x.color === 'primary' || x.color === 'green');
      const rest = keyed.filter(x => x.color !== 'primary' && x.color !== 'green');
      const keyHtml = keys.length ? keys.map(x => {
        const cls = x.color ? ` dlbtn ${x.color}` : ' dlbtn';
        const low = x.base.toLowerCase();
        // 基因组图 (SVG)：直接内联预览
        if (low.endsWith('.svg')) {
          const normPath = String(x.f.path).replace(/\\/g, '/');
          const segs = normPath.split('/').map(encodeURIComponent).join('/');
          return `<div style="background:#fff;border-radius:8px;margin:6px 0;padding:8px"><img src="/tool_runs/${encodeURIComponent(run.name)}/${segs}" alt="${esc(x.base)}" style="display:block;width:100%;height:220px;object-fit:contain;border:1px solid var(--line-100);border-radius:8px;background:#fff" loading="lazy"><a class="${cls.trim()}" href="/tool_runs/${encodeURIComponent(run.name)}/${encodeURIComponent(normPath)}" title="${esc(x.f.path)}【${esc(x.f.size)}】">⬇ 下载 ${esc(x.label)}</a></div>`;
        }
        return `<a class="${cls.trim()}" href="/tool_runs/${encodeURIComponent(run.name)}/${encodeURIComponent(String(x.f.path).replace(/\\/g, '/'))}" title="${esc(x.f.path)}【${esc(x.f.size)}】">⬇ ${esc(x.label)}</a>`;
      }).join('') : '<p class="hint">（暂无关键产物）</p>';
      const restHtml = rest.length ? `<details class="rpt-sec" style="margin:8px 0 0"><summary>全部文件（${rest.length}）</summary><div class="rpt-sec-body" style="display:flex;flex-wrap:wrap;gap:4px">${rest.map(x => {
        const cls = x.color ? ` dlbtn ${x.color}` : ' dlbtn';
        return `<a class="${cls.trim()}" href="/tool_runs/${encodeURIComponent(run.name)}/${esc(x.f.path)}" title="${esc(x.f.path)}【${esc(x.f.size)}】">⬇ ${esc(x.label)}</a>`;
      }).join('')}</div></details>` : '';
      const fileHtml = open ? `<div class="dlbtn-row">${keyHtml}</div>${restHtml}` : '';
      const rptBtn = (run.name.startsWith('identify_') || run.name.startsWith('contigs_'))
        ? `<button class="btn small" onclick="rhOpenReport('${esc(run.name)}')" title="查看该运行的分类报告（桑基图/分类表/明细）">📊 报告</button> `
        : '';
      return `<div style="border-bottom:1px dashed var(--line-100);padding:5px 0">
        <div style="display:flex;gap:6px;align-items:center;flex-wrap:wrap">
          <b class="mono" style="font-size:12px;cursor:pointer" onclick="rhFiles('${esc(run.name)}')" title="点击展开/收起文件列表">${esc(run.name)}</b>
          <span class="hint">${(run.files || []).length} 文件</span>
          ${chainBtns}
          ${rptBtn}
          <button class="btn small" onclick="rhFiles('${esc(run.name)}')">📁 文件</button>
          <button class="btn small danger" onclick="rhDel(this, '${esc(run.name)}')" title="删除该运行目录">🗑</button>
        </div>
        ${open ? (run.name.startsWith('contigs_') ? `<div id="rhcontig-${esc(run.name)}" style="margin-top:8px"><p class="hint">加载病毒序列分类明细…</p></div>` : '') + `<div style="margin-top:4px" class="dlbtn-row">${fileHtml}</div>` : ''}
      </div>`;
    }).join('');
  }
}

function rhToggle(boxId) {
  const body = $(boxId).querySelector('.rh-body');
  body.style.display = body.style.display === 'none' ? '' : 'none';
}
function rhFiles(name) {
  if (RH_OPEN.has(name)) RH_OPEN.delete(name); else RH_OPEN.add(name);
  rhRefresh();
  if (RH_OPEN.has(name) && name.startsWith('contigs_')) rhLoadContigDetail(name);
}

/* 历史运行点开 → 展示该 run 的病毒序列分类（contig 明细）关键信息。
   读取 virus_classification.tsv，展示每条病毒 contig 的分类/长度/得分。 */
async function rhLoadContigDetail(run) {
  const box = document.getElementById('rhcontig-' + run);
  if (!box) {   // rhRefresh 异步渲染可能未完成，稍后重试
    setTimeout(() => rhLoadContigDetail(run), 300);
    return;
  }
  try {
    const r = await fetch(`/tool_runs/${encodeURIComponent(run)}/virus_classification.tsv`);
    if (!r.ok) { box.innerHTML = '<p class="hint">（无病毒分类明细）</p>'; return; }
    const text = await r.text();
    const lines = text.split(/\r?\n/).filter(Boolean);
    if (lines.length < 2) { box.innerHTML = '<p class="hint">（无病毒分类明细）</p>'; return; }
    const head = lines[0].split('\t');
    const idx = k => head.indexOf(k);
    const ic = idx('contig'), it = idx('taxon'), il = idx('length'), isc = idx('score'), ih = idx('host');
    const rows = lines.slice(1).map(l => l.split('\t'));
    const show = rows.slice(0, 12);
    const n = rows.length;
    const trs = show.map(row => {
      const c = ic >= 0 ? row[ic] : '';
      const t = it >= 0 ? row[it] : '';
      const l = il >= 0 ? row[il] : '';
      const s = isc >= 0 ? row[isc] : '';
      const h = ih >= 0 ? row[ih] : '';
      return `<tr><td class="mono" style="font-size:11px">${esc(c)}</td><td>${esc(t)}</td><td class="num">${esc(l)}</td><td class="num">${esc(s)}</td><td>${esc(h)}</td></tr>`;
    }).join('');
    box.innerHTML = `<div style="font-size:12.5px;font-weight:700;margin-bottom:6px;color:#1a5276">病毒序列分类（contig 明细 · ${n} 条）</div>
      <div style="max-height:260px;overflow:auto;border:1px solid #eee;border-radius:6px">
      <table class="tb"><thead><tr><th>Contig</th><th>分类</th><th class="num">长度</th><th class="num">得分</th><th>宿主</th></tr></thead>
      <tbody>${trs}</tbody></table></div>
      ${n > show.length ? `<p class="hint" style="margin-top:4px">共 ${n} 条，仅显示前 ${show.length} 条（点“📊 报告”看完整）</p>` : ''}`;
  } catch (e) { box.innerHTML = '<p class="hint">明细加载失败</p>'; }
}
async function rhDel(btn, name) {
  if (!confirm(`删除运行 ${name}（目录与全部产物，不可恢复）？`)) return;
  try {
    const r = await fetch(`/api/tool/runs/${encodeURIComponent(name)}/delete`,
                          { method: 'POST' });
    if (!r.ok) { alert((await r.json()).error || '删除失败'); return; }
    RH_OPEN.delete(name);
    rhRefresh();
  } catch (e) { alert('无法连接: ' + e); }
}
function rhGo(btn, runName, boxId, idx) {
  const def = RH_DEFS.find(d => d.box === boxId);
  const run = (RH_DATA || []).find(x => x.name === runName);
  if (!def || !run) return;
  const chain = (def.chains(run.files || [], run) || [])[idx];
  if (!chain) return;
  // 暂存全部回填值 → 跳目标页 → rhApply 落位（同页目标立即填充）
  let fills = Object.assign({}, chain.fills);
  if (chain.page) {
    const samePage = !chain.page.startsWith('/') ||
      location.pathname === chain.page.split('#')[0];
    const missing = Object.keys(fills).filter(k => !$(k));
    if (samePage && !missing.length) {
      rhApplyFills(fills);
      toast('已回填下游输入', Object.keys(fills).join(', '), { ttl: 3000 });
      if (chain.page.includes('#')) location.hash = chain.page.split('#')[1];
      return;
    }
    try { sessionStorage.setItem('vp_runfill', JSON.stringify(fills)); } catch (e) {}
    location.href = chain.page;
  } else {
    rhApplyFills(fills);
  }
}
function rhApplyFills(fills) {
  for (const [id, val] of Object.entries(fills)) {
    const el = $(id);
    if (!el) continue;
    if (el.type === 'checkbox') el.checked = !!val;
    else if (el.tagName === 'SELECT') {
      if ([...el.options].some(o => o.value === val)) el.value = val;
    } else el.value = val;
  }
}
function rhApplyStashed() {
  let fills = null;
  try { fills = JSON.parse(sessionStorage.getItem('vp_runfill') || 'null'); } catch (e) {}
  if (!fills) return;
  const applied = [];
  for (const id of Object.keys(fills)) {
    if ($(id)) { rhApplyFills({ [id]: fills[id] }); applied.push(id); delete fills[id]; }
  }
  if (applied.length) {
    if (Object.keys(fills).length)
      try { sessionStorage.setItem('vp_runfill', JSON.stringify(fills)); } catch (e) {}
    else sessionStorage.removeItem('vp_runfill');
    toast('已回填历史运行的衔接输入', applied.join(', '), { ttl: 4000 });
  }
}
window.__rhOnTasks = tasks => {
  // 任何任务刚结束（running → 终态）时刷新历史区
  const done = new Set(tasks.filter(t => ['done', 'failed', 'cancelled'].includes(t.status))
                            .map(t => t.id));
  const newly = [...done].filter(id => !RH_SEEN.has(id));
  if (RH_SEEN.size && newly.length) rhRefresh();
  RH_SEEN.clear(); tasks.forEach(t => RH_SEEN.add(t.id));
};
const RH_SEEN = new Set();
function rhInitContainers() {
  for (const def of RH_DEFS) {
    const box = $(def.box);
    if (!box || box.dataset.init) continue;
    box.dataset.init = '1';
    box.innerHTML = `<div class="rh-head" style="display:flex;gap:8px;align-items:center;margin-top:8px">
        <button class="btn small" onclick="rhToggle('${def.box}')">🕘 历史运行（<span class="rh-count">…</span>）</button>
        <span class="hint">重启不丢 · 可复用为下游输入 · 可删除（默认折叠）</span>
      </div>
      <div class="rh-body" style="display:none;margin-top:4px"></div>`;
  }
}
document.addEventListener('DOMContentLoaded', () => {
  rhInitContainers();
  rhRefresh();
  setTimeout(rhApplyStashed, 350);   // 等各页自身 init（下拉等）先就绪
});

// ---------------- 跨模块文件交接（下载页 → 其它模块预填） ----------------
// 把文件参数暂存到 sessionStorage，再跳到目标模块页；目标页用 takePrefill 取出并预填输入框。
function sendTo(page, params) {
  try {
    sessionStorage.setItem('vp_prefill',
      JSON.stringify(Object.assign({ page: page, ts: Date.now() }, params)));
  } catch (e) {}
  if (page === 'fastp' || page === 'identify') {
    location.href = '/tools?g=sample#t-' + page;
  } else {
    location.href = '/' + page;
  }
}

// 取出本页的预填参数（page 不匹配则不动；取出即销毁，避免二次误填）
function takePrefill(page) {
  try {
    const d = JSON.parse(sessionStorage.getItem('vp_prefill') || 'null');
    if (!d || d.page !== page) return null;
    sessionStorage.removeItem('vp_prefill');
    return d;
  } catch (e) { return null; }
}

// ---------------- 任务轮询 + 结果预览 ----------------
let pollTimer = null;
const taskLogOpen = new Set();   // 展开日志的任务 id
const notifiedTasks = new Set(); // 已弹过结束通知的任务
let prevRunning = new Set();     // 上次轮询时运行中的任务

// 自适应轮询：有任务在跑时 2.5s（保持进度实时），全部空闲时降到 8s。
// 任务列表每 2.5s 全量回传，长时间挂着页面会产生大量无谓请求，
// 空闲降频可在不影响体验的前提下显著减轻服务端与浏览器负担。
const POLL_BUSY_MS = 2500;
const POLL_IDLE_MS = 8000;
let lastTaskList = [];
let connDown = false;

function startPolling() {
  if (pollTimer) clearTimeout(pollTimer);
  schedulePoll(0);
}

function schedulePoll(delay) {
  if (pollTimer) clearTimeout(pollTimer);
  pollTimer = setTimeout(async () => {
    // 递归调度不同于 setInterval：refreshTasks 抛异常会中断整条轮询链，
    // 这里必须兜住（原实现用 setInterval 不会有这个问题）
    try {
      await refreshTasks();
    } catch (e) {
      console.error('[poll] refreshTasks failed:', e);
    }
    const busy = lastTaskList.some(t => t.status === 'running') || connDown;
    schedulePoll(busy ? POLL_BUSY_MS : POLL_IDLE_MS);
  }, delay);
}

function taskStatusText(s) {
  return {running: t('c.st.running'), done: t('c.st.done'),
          failed: t('c.st.failed'), cancelled: t('c.st.cancelled')}[s] || s;
}

/* 任务结果预览面板（后端 result 字段为规范化预览结构） */
function taskResultHtml(res) {
  if (!res || !res.kind) return '';
  const chips = [];
  let links = [];
  if (res.kind === 'sample') {
    (res.stages || []).forEach(st => {
      if (st.status === 'done' && st.summary) {
        chips.push(`<span class="tr-chip" title="${esc(st.name)}">${esc(st.name)}：${esc(st.summary)}</span>`);
      }
    });
    if (res.report) {
      links.push(`<a class="btn small primary" href="${esc(res.report)}" target="_blank">${t('c.openReport')}</a>`);
      links.push(`<a class="btn small" href="/results" target="_blank">${t('nav.results')}</a>`);
    }
  } else if (res.kind === 'tool') {
    for (const [k, v] of Object.entries(res.stats || {})) {
      if (v === null || v === '') continue;
      const s = String(v);
      const isPath = /[\\/]/.test(s);
      const label = isPath ? s.split(/[\\/]/).pop() : s;
      chips.push(`<span class="tr-chip" title="${esc(k)}=${esc(s)}">${esc(k)}=${esc(label)}</span>`);
    }
    (res.files || []).forEach(f => {
      links.push(`<a class="btn small" href="/tool_runs/${esc(f.path)}" download title="${esc(f.path)}">📄 ${esc(f.path.split('/').pop())} <span class="hint">(${esc(f.size)})</span></a>`);
    });
  }
  if (!chips.length && !links.length) return '';
  return `<div class="task-result">
    <div class="tr-title">${t('c.resultTitle')}</div>
    ${chips.length ? `<div class="tr-stats">${chips.slice(0, 14).join('')}</div>` : ''}
    ${links.length ? `<div class="tr-links">${links.join('')}</div>` : ''}
  </div>`;
}

async function refreshTasks() {
  let tasks;
  try {
    const r = await fetch('/api/tasks');
    tasks = await r.json();
    lastTaskList = tasks || [];
    connFails = 0;
    connDown = false;
    setConnBanner(false);
    if (window.__rhOnTasks) { try { window.__rhOnTasks(tasks); } catch (e2) {} }
  } catch (e) {
    // 连接异常时保持短间隔重试，避免失联横幅要等满 3 次空闲周期才出现
    connDown = true;
    if (++connFails >= 3) setConnBanner(true);
    return;
  }
  try {
  const box = $('tasks');
  // 页面没有「任务」汇总区（如工具页已移除底部任务面板）时跳过面板渲染；
  // 连接横幅、结束通知、工具卡内日志注入（injectStageLogs）仍照常执行。
  if (!box) {
    /* 无面板：仅保留通知与卡内日志注入 */
  } else if (!tasks.length) {
    box.innerHTML = `<p class="hint">${t('c.noTasks')}</p>`;
  } else {
  box.innerHTML = tasks.map(task => {
    const pct = Math.round((task.pct || 0) * 100);
    const stat = taskStatusText(task.status);
    const logHtml = taskLogOpen.has(task.id)
      ? `<pre data-autoscroll="${task.status === 'running' ? 1 : 0}">${esc((task.log || []).join('\n'))}</pre>` : '';
    const err = task.error ? `<div class="err">${esc(task.error)}</div>` : '';
    let timing = '';
    if (task.status === 'running') {
      const el = fmtDur(Date.now() / 1000 - (task.started || Date.now() / 1000));
      timing = `<span class="hint" style="margin-left:10px">${t('c.elapsed')} ${el}` +
        (task.eta != null ? ` · ${t('c.eta')} ~${fmtDur(task.eta)}` : '') + '</span>';
    } else if (task.finished && task.started) {
      timing = `<span class="hint" style="margin-left:10px">${t('c.used')} ${fmtDur(task.finished - task.started)}</span>`;
    }
    const preview = task.status === 'done' ? taskResultHtml(task.result) : '';
    return `<div class="task ${task.status}">
      <div style="display:flex;justify-content:space-between;align-items:center">
        <span class="tname">${esc(task.name)}</span>
        <span><span class="tstat ${task.status}">${stat}</span>${timing}
        ${task.status === 'running' ? `<button class="btn small danger" onclick="cancelTask('${task.id}')">${t('c.cancel')}</button>` : ''}
        </span>
      </div>
      <div class="tmsg">${esc(task.msg || '')}</div>
      <div class="bar"><div style="width:${pct}%"></div></div>
      ${err}${preview}
      <button class="btn small" onclick="toggleLog('${task.id}')">${taskLogOpen.has(task.id) ? t('c.collapseLog') : t('c.expandLog')}</button>
      ${logHtml}
    </div>`;
  }).join('');
  }

  // 任务结束 → toast/系统通知 + 刷新管道视图
  const nowRunning = new Set(tasks.filter(x => x.status === 'running').map(x => x.id));
  for (const task of tasks) {
    if (task.status !== 'running' && prevRunning.has(task.id)
        && !notifiedTasks.has(task.id)) {
      notifiedTasks.add(task.id);
      if (task.status === 'done') {
        const hasReport = task.result && task.result.kind === 'sample' && task.result.report;
        toast(`${t('c.resultNotify')} · ${task.name}`,
              task.result && task.result.kind === 'sample' ? '' : (task.msg || ''),
              {actions: hasReport ? [
                {label: t('c.openReport'), primary: true,
                 onClick: () => window.open(task.result.report, '_blank')},
                {label: t('c.close'), onClick: () => {}},
              ] : [{label: t('c.close'), onClick: () => {}}]});
        sysNotify(t('c.resultNotify'), task.name);
      } else if (task.status === 'failed') {
        toast(`${t('c.failedNotify')} · ${task.name}`, task.error || '',
              {kind: 'failed', ttl: 12000, actions: [
                {label: t('c.close'), onClick: () => {}}]});
        sysNotify(t('c.failedNotify'), task.name);
      }
      loadSamples();
      if (curSample) selectSample(curSample);
      if (typeof loadToolRuns === 'function') loadToolRuns();
      break;   // 一次轮询只弹一条，避免轰炸
    }
  }
  prevRunning = nowRunning;
  } catch (e) { console.warn('任务面板渲染失败:', e); }
  injectStageLogs(tasks);
  autoscrollLogs();
}

// ---------------- 阶段卡片内嵌运行日志 ----------------
// 任务名按服务端语言生成；两种语言的标签都参与匹配
const STAGE_LABELS = {
  subsample:['预处理(子采样)', 'Prep (subsample)'],
  fastp:    ['⓪ Fastp 质控', '⓪ Fastp QC'],
  fq2fa:    ['⓪b 序列转换(FASTQ→FASTA)', '⓪b FASTQ→FASTA'],
  host:     ['① 宿主去除', '① Host removal'],
  virus:    ['② 病毒筛查与提取', '② Virus screening'],
  assembly: ['③ 组装·分类·提取', '③ Assembly'],
  verify:   ['③b 候选序列验证', '③b Candidate verify'],
  hostana:  ['④ 宿主预测(ICTV)', '④ Host prediction (ICTV)'],
  orf:      ['⑥ ORF 预测', '⑥ ORF prediction'],
  orfa:     ['⑥b ORF 功能注释', '⑥b ORF annotation'],
  phylo:    ['⑦ 进化树与 SDT', '⑦ Phylogeny'],
  primer:   ['⑧ 引物设计', '⑧ Primer design'],
  gbdraw:   ['⑨ 基因组图', '⑨ Genome plots'],
  report:   ['⑩ 可视化报告', '⑩ Visual report'],
};

const TOOL_LABELS = {
  fastp:    ['质控预处理', 'QC preprocess'],
  hostremoval: ['宿主去除与序列提取', 'Host removal'],
  hostpredict: ['宿主预测', 'Host prediction'],
  orf:      ['ORF 预测', 'ORF predict'],
  orfa:     ['功能注释', 'ORF annotate'],
  genoplot: ['基因组图谱', 'Genome plots'],
  primer:   ['引物设计', 'Primer design'],
  identify: ['病毒鉴定', 'Virus identify'],
  assemble: ['病毒组装', 'Virus assembly'],
  contigs:  ['contig分类', 'Contig classify'],
  verify:   ['候选序列验证', 'Candidate verify'],
  consensus: ['共识序列与变异', 'Consensus & variants'],
  ncbi:     ['NCBI下载', 'NCBI download'],
  gbdown:   ['GenBank下载', 'GenBank download'],
  gbimport: ['GenBank导入', 'GenBank import'],
  synteny:  ['同属比较', 'Synteny'],
};

// 数据库构建页：建库任务名 → 卡片内嵌日志容器（buildrun-<key>）。
// 两种语言的标签都参与匹配；refvirus/rvdb 用各自独有词区分。
const BUILD_TASK_LABELS = {
  taxonomy: ['Taxonomy', 'taxonomy'],
  host:     ['宿主库构建', 'Host DB build'],
  refvirus: ['RefSeq Viral'],
  rvdb:     ['RVDB 库构建', 'C-RVDB'],
  k2:       ['Kraken2 库转换', 'Kraken2 convert'],
};

function injectStageLogs(tasks) {
  snapshotLogScroll();
  document.querySelectorAll('.stage-log').forEach(el => { el.innerHTML = ''; });
  document.querySelectorAll('.toolrun').forEach(el => { el.innerHTML = ''; });
  document.querySelectorAll('.buildrun').forEach(el => { el.innerHTML = ''; });
  for (const task of tasks) {
    if (!['running', 'done', 'failed'].includes(task.status)) continue;
    // 样品管道任务：按后端上报的当前阶段（task.stage）精确归属——
    // 日志只进"正在跑的这一张"阶段卡，不再整条链每张卡都重复一份。
    // stage 尚未上报（刚启动）时退化为名称匹配，且仅当名称只命中
    // 一个阶段时注入，避免多阶段链一次性铺满所有卡。
    // 另按样品名过滤：其它样品的管道任务不出现在当前样品的卡上。
    const nmStage = task.name || '';
    const sampleOk = !curSample || nmStage.includes(curSample);
    if (sampleOk) {
      let hitKeys = [];
      if (task.stage && STAGE_LABELS[task.stage]) {
        hitKeys = [task.stage];
      } else if (!task.stage) {
        hitKeys = Object.entries(STAGE_LABELS)
          .filter(([, labels]) => labels.some(label => nmStage.includes(label)))
          .map(([key]) => key);
        if (hitKeys.length > 1) hitKeys = [];   // 多阶段链未上报 → 等首个进度
      }
      for (const key of hitKeys) {
        const box = $('log-' + key);
        if (box) box.innerHTML = stageLogHtml(task);
      }
    }
    // 专项分析工具任务（按工具标签匹配 → 卡内嵌 toolrun 容器）
    for (const [key, labels] of Object.entries(TOOL_LABELS)) {
      const nm = task.name || '';
      if (labels.some(label => nm.includes(label))
          && (nm.includes('工具·') || nm.includes('Tool·')
              || nm.includes('NCBI') || nm.includes('GenBank')
              || nm.includes('同属') || nm.includes('Synteny'))) {
        const boxId = 'toolrun-' + key;
        const box = $(boxId);
        if (box) box.innerHTML = toolRunHtml(task, boxId);
      }
    }
    // 建库任务 → 数据库构建页对应卡片内嵌日志（.buildrun 容器，按任务名匹配）
    for (const [key, labels] of Object.entries(BUILD_TASK_LABELS)) {
      const nm = task.name || '';
      if (labels.some(label => nm.includes(label))) {
        const box = $('buildrun-' + key);
        if (box) box.innerHTML = stageLogHtml(task);
      }
    }
  }
  autoscrollLogs();
}

// 输出区标签状态（boxId → 'result' | 'log' | 'files'）。
// 输出区已改为纵向平铺（日志/结果/文件同时展示），不再渲染标签页按钮；
// 这里保留状态与 setOutTab 仅为向后兼容（防御性，不再被触发）。
const OUT_TAB_STATE = {};
const OUT_TASKS = {};

/* 工具产物文件名 → 语义化下载名 + 按钮配色。
   用于把原始英文文件名（run.log / contigs.filtered.fasta 等）翻译成
   「运行日志 / 组装结果」这类可读名称，并按主结果/辅助产物配色。 */
const DL_LABEL_SPECS = [
  // [文件名包含, 语义化名, 按钮色]（先匹配先得，结果类放前）
  // —— 组装 / 分类主结果 ——
  ['contigs.filtered.fasta', '组装结果 (contigs)', 'primary'],
  ['transcripts.fasta', '组装结果 (transcripts)', 'primary'],
  ['viral_contigs.fasta', '病毒序列 (viral contigs)', 'primary'],
  ['virus_classification.tsv', '病毒分类结果', 'primary'],
  ['viral_ids.tsv', '分类序列 ID', 'primary'],
  ['viral_sequences.', '病毒序列', 'primary'],
  ['contigs.fasta', '组装结果 (contigs)', 'primary'],
  ['scaffolds.fasta', '组装结果 (scaffolds)', 'primary'],
  ['pyrodigal.faa', 'ORF 蛋白', 'primary'],
  ['pyrodigal.ffn', 'ORF 核苷酸', 'primary'],
  ['host_prediction.tsv', '宿主预测结果', 'primary'],
  // —— 报告 / 汇总 ——
  ['report.html', '分析报告', 'green'],
  ['kreport', '分类报告 (kreport)', 'green'],
  ['keep_R1.fastq.gz', '保留 reads R1', 'green'],
  ['keep_R2.fastq.gz', '保留 reads R2', 'green'],
  ['kept_R1.fastq.gz', '保留 reads R1', 'green'],
  ['kept_R2.fastq.gz', '保留 reads R2', 'green'],
  ['taxburst', 'Krona 旭日图', 'green'],
  // —— 图 / 结构文件 ——
  ['.svg', '基因组图 (SVG)', 'green'],
  ['.gfa', '组装图 (GFA)', 'amber'],
  ['before_rr.fasta', '纠错前序列', 'amber'],
  ['fastp.', '质控报告', 'amber'],
  // —— 日志 / 配置 ——
  ['run.log', '运行日志', 'gray'],
  ['spades.log', 'SPAdes 日志', 'gray'],
  ['warnings.log', '警告日志', 'gray'],
  ['spades.sh', 'SPAdes 脚本', 'gray'],
  ['spades.yaml', 'SPAdes 配置', 'gray'],
  ['run_spades.sh', 'SPAdes 脚本', 'gray'],
  ['run_spades.yaml', 'SPAdes 配置', 'gray'],
  ['params.txt', '参数表', 'gray'],
  ['input_dataset.yaml', '输入数据描述', 'gray'],
  ['dataset.info', '数据描述', 'gray'],
  ['.json', '结果数据 (JSON)', 'gray'],
  ['summary.', '阶段汇总', 'gray']
];
function dlLabel(base, path) {
  const low = String(base || path || '').toLowerCase();
  for (const [key, label, color] of DL_LABEL_SPECS) {
    if (low.includes(key)) return { label, color };
  }
  return { label: base, color: '' };   // 未知 → 原始文件名
}

function outFilesHtml(task) {
  const files = (task.result && task.result.files) || [];
  if (!files.length) {
    return `<p class="hint" style="margin:6px 0">${
      task.status === 'running' ? t('c.filesRunning') : t('c.noFiles')}</p>`;
  }
  // 链接必须含运行名（/tool_runs/<run>/<file>），task.result.run 由后端注入；
  // download 属性带运行名前缀：不同运行的同名产物（viral_ids.tsv 等）
  // 下载到本地不再互相覆盖。
  const run = (task.result && task.result.run) || '';
  const btns = files.map(f => {
    const normPath = String(f.path).replace(/\\/g, '/');
    const segs = normPath.split('/').map(encodeURIComponent).join('/');
    const base = String(f.path).split(/[\\/]/).pop();
    const dl = run ? `${run}_${base}` : base;
    const href = run ? `/tool_runs/${encodeURIComponent(run)}/${segs}` : '#';
    const { label, color } = dlLabel(base, f.path);
    const cls = color ? ` dlbtn ${color}` : ' dlbtn';
    const low = base.toLowerCase();
    // 基因组图 (SVG)：直接内联预览，不再只给下载按钮。
    if (low.endsWith('.svg') && run) {
      const img = `<img src="/tool_runs/${encodeURIComponent(run)}/${segs}" alt="${esc(base)}" style="display:block;width:100%;height:220px;object-fit:contain;border:1px solid var(--line-100);border-radius:8px;margin:8px 0;background:#fff" loading="lazy">`;
      return `<div class="svg-prev" style="border-radius:8px;padding:8px;background:#fff">${img}<a class="${cls.trim()}" href="${href}" download="${esc(dl)}" title="${esc(f.path)}【${esc(f.size)}】">⬇ 下载 ${esc(label)}</a></div>`;
    }
    return `<a class="${cls.trim()}" href="${href}" download="${esc(dl)}" title="${esc(f.path)}【${esc(f.size)}】">⬇ ${esc(label)}</a>`;
  }).join('');
  return `<div class="dlbtn-row">${btns}</div>`;
}

/* 工具运行的关键数字 chips（仅标量；路径取文件名；run 名已在状态徽标里，不重复） */
const TOOL_STAT_LABELS = {
  n_classified: ['已分类序列', 'Classified seqs'],
  n_extracted:  ['提取病毒序列', 'Extracted viral seqs'],
  n_contigs:    ['过滤后 contigs', 'Contigs kept'],
  n_viral:      ['病毒 contigs', 'Viral contigs'],
  total_bp:     ['总长度(bp)', 'Total bp'],
};

function toolStatChips(res) {
  if (!res || res.kind !== 'tool') return '';
  const zh = (typeof VP_LANG === 'undefined' || VP_LANG !== 'en');
  const chips = [];
  for (const [k, v] of Object.entries(res.stats || {})) {
    if (v === null || v === '' || k === 'run') continue;
    const lab = (TOOL_STAT_LABELS[k] || [k, k])[zh ? 0 : 1];
    const s = String(v);
    const shown = /[\\/]/.test(s) ? s.split(/[\\/]/).pop() : s;
    chips.push(`<span class="tr-chip" title="${esc(k)}=${esc(s)}">${esc(lab)}=${esc(shown)}</span>`);
  }
  return chips.length ? `<div class="tr-stats">${chips.slice(0, 14).join('')}</div>` : '';
}

function toolRunHtml(task, boxId) {
  boxId = boxId || '';
  if (boxId) { OUT_TASKS[boxId] = task; if (!OUT_TAB_STATE[boxId]) OUT_TAB_STATE[boxId] = 'result'; }
  const pct = Math.round((task.pct || 0) * 100);
  const stat = taskStatusText(task.status);
  const err = task.error ? `<div class="err" style="margin:6px 0 0">${esc(task.error)}</div>` : '';
  const chips = task.status === 'done' ? toolStatChips(task.result) : '';
  const lines = esc((task.log || []).join('\n'));
  // 纵向平铺：状态/进度 → 关键数字 → 运行日志 → Downloads（下载）。
  // 产物文件统一在 Downloads 一处列出（分类报告面板内不再重复）。
  return `<div style="margin:10px 0 2px">
    <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:4px">
      <span class="sbadge ${task.status}">${stat} · ${esc(task.name)}</span>
      ${task.status === 'running' ? `<button class="btn small danger" onclick="cancelTask('${task.id}')">${t('c.cancel')}</button>` : ''}
    </div>
    <div class="bar"><i style="width:${pct}%"></i></div>
    ${err}${chips}
    <div style="font-size:12.5px;font-weight:700;margin:10px 0 4px">${t('c.runlog')}</div>
    <pre class="logbox" style="max-height:180px;overflow-y:auto" data-autoscroll="${task.status === 'running' ? 1 : 0}">${lines || t('c.noLog')}</pre>
    <div style="font-size:12.5px;font-weight:700;margin:10px 0 4px">${t('c.downloads')}</div>
    ${outFilesHtml(task)}
  </div>`;
}

function setOutTab(boxId, tab) {
  OUT_TAB_STATE[boxId] = tab;
  const task = OUT_TASKS[boxId];
  const box = $(boxId);
  if (task && box) box.innerHTML = toolRunHtml(task, boxId);
}
document.addEventListener('click', ev => {
  const btn = ev.target.closest?.('.out-tab');
  if (!btn) return;
  const holder = btn.closest('.toolrun');
  if (holder) setOutTab(holder.id, btn.dataset.tab);
});

function stageLogHtml(task) {
  const stat = taskStatusText(task.status);
  const lines = esc((task.log || []).join('\n'));
  const closed = logClosedSet().has(task.id);
  return `
    <div style="display:flex;justify-content:space-between;align-items:center;margin:8px 0 4px">
      <span class="sbadge ${task.status}">${stat} · ${t('c.liveLog')}</span>
      <span>
        <button class="btn small" onclick="logToggle('${task.id}')" title="折叠/展开日志">${closed ? '▸' : '▾'} ${t('c.logBtn', '日志')}</button>
        <button class="btn small" onclick="copyStageLog(this)">${t('c.copyLog')}</button>
        ${task.status === 'running' ? `<button class="btn small danger" onclick="cancelTask('${task.id}')">${t('c.cancel')}</button>` : ''}
      </span>
    </div>
    <pre class="logbox" style="margin:0;${closed ? 'display:none' : ''}" data-tid="${esc(task.id)}" data-autoscroll="${task.status === 'running' ? 1 : 0}">${lines || t('c.noLog')}</pre>`;
}

function copyStageLog(btn) {
  const pre = btn.closest('.stage-log')?.querySelector('pre.logbox');
  if (!pre) return;
  const text = pre.textContent;
  (navigator.clipboard?.writeText(text) ?? Promise.reject())
    .then(() => { btn.textContent = t('c.copied'); setTimeout(() => btn.textContent = t('c.copyLog'), 1500); })
    .catch(() => {
      const ta = document.createElement('textarea');
      ta.value = text; document.body.appendChild(ta); ta.select();
      document.execCommand('copy'); ta.remove();
      btn.textContent = t('c.copied'); setTimeout(() => btn.textContent = t('c.copyLog'), 1500);
    });
}

/* 日志滚动策略：重绘前记录每个日志框 scrollTop 与是否贴底；重绘后——
   运行中且原本贴底 → 跟随新日志；其余（含已结束任务、用户上翻阅读中）
   → 精确恢复原位置，不再被 2.5s 轮询重绘拉回顶部。 */
const _logScrollMemo = new WeakMap();
function snapshotLogScroll() {
  document.querySelectorAll('.toolrun, .stage-log, .buildrun').forEach(box => {
    const st = [...box.querySelectorAll('pre.logbox')].map(pre => ({
      top: pre.scrollTop,
      atBottom: pre.scrollHeight - pre.scrollTop - pre.clientHeight < 16
    }));
    _logScrollMemo.set(box, st);
  });
}
function restoreLogScroll() {
  document.querySelectorAll('.toolrun, .stage-log, .buildrun').forEach(box => {
    const memo = _logScrollMemo.get(box) || [];
    [...box.querySelectorAll('pre.logbox')].forEach((pre, i) => {
      const st = memo[i];
      if (pre.dataset.autoscroll === '1' && (!st || st.atBottom)) {
        pre.scrollTop = pre.scrollHeight;
      } else if (st) {
        pre.scrollTop = st.top;
      }
    });
  });
}
function autoscrollLogs() { restoreLogScroll(); }

function toggleLog(tid) {
  taskLogOpen.has(tid) ? taskLogOpen.delete(tid) : taskLogOpen.add(tid);
  refreshTasks();
}
function logClosedSet() {
  try { return new Set(JSON.parse(sessionStorage.getItem('vp_logclosed') || '[]')); }
  catch (e) { return new Set(); }
}
function logToggle(tid) {
  const s = logClosedSet();
  s.has(tid) ? s.delete(tid) : s.add(tid);
  try { sessionStorage.setItem('vp_logclosed', JSON.stringify([...s])); } catch (e) {}
  refreshTasks();
}

async function cancelTask(tid) {
  await fetch('/api/task/' + tid + '/cancel', {method: 'POST'});
  refreshTasks();
}

// ---------------- 库状态 ----------------
async function loadDbs() {
  try {
    const r = await fetch('/api/dbs');
    const d = await r.json();
    const badge = (ok, label) =>
      `<span class="badge ${ok ? 'ok' : 'no'}">${label} ${ok ? '✔' : '✘'}</span>`;
    const html = badge(d.taxonomy.ready, 'NCBI Taxonomy')
               + badge(d.host.ready, t('hm.db.host'))
               + badge(d.virus.ready, t('hm.db.virus'));
    const put = (id, text) => { const el = $(id); if (el) el.innerHTML = text; };
    put('dbStatus', html);
    put('taxStat', badge(d.taxonomy.ready, 'Taxonomy'));
    put('hostStat', badge(d.host.ready, t('hm.db.host')));
    put('virusStat', badge(d.virus.ready, t('hm.db.virus')));
  } catch (e) {}
}

async function loadTools() {
  const box = $('toolList');
  if (!box) return;
  try {
    const r = await fetch('/api/tools');
    const d = await r.json();
    box.innerHTML = Object.entries(d.tools).map(([k, v]) =>
      `<div class="trow"><span>${esc(k)}</span>
       <span class="${v ? 'ok' : 'no'}">${v ? '✔ ' + esc(v.split('\\').pop()) : '✘ ' + t('bd.notFound', '未找到')}</span></div>`
    ).join('');
  } catch (e) {}
}

// ---------------- 文件浏览（记忆上次目录） ----------------
let browseTarget = null;
let browseMode = 'file';   // 'file' 选数据文件 | 'dir' 选目录
let browseCwd = '.';

function browse(inputId) {
  browseMode = 'file';
  browseTarget = inputId;
  $('dlgMask').style.display = 'flex';
  const ttl = $('dlgTitle');
  if (ttl) ttl.textContent = t('c.chooseFile');
  const last = localStorage.getItem('vp_browse_cwd');
  loadBrowse(last || '.');
}

function browseDir(inputId) {
  browseMode = 'dir';
  browseTarget = inputId;
  $('dlgMask').style.display = 'flex';
  const ttl = $('dlgTitle');
  if (ttl) ttl.textContent = t('c.chooseDir', '选择目录');
  const last = localStorage.getItem('vp_browse_cwd');
  loadBrowse(last || '.');
}

function pickDir() {
  if (browseCwd === '此电脑') { alert('请先进入某个目录'); return; }
  if (browseTarget && $(browseTarget)) $(browseTarget).value = browseCwd;
  closeDlg();
}

function _joinBrowse(cwd, name) {
  /* 拼接子路径：
     - '此电脑' 层级的子项就是盘符绝对路径（如 D:\），直接返回；
     - 盘根目录（D:\ 带尾反斜杠）直接拼名字；
     - 其余用 / 连接（服务端 normpath 兼容正斜杠）。 */
  if (cwd === '此电脑') return name;
  if (/[\\/]$/.test(cwd)) return cwd + name;
  return cwd + '/' + name;
}

async function loadBrowse(path) {
  try {
    // 选目录模式带 all=1：服务端不限文件扩展名，方便查看目录里有什么
    const r = await fetch('/api/browse?path=' + encodeURIComponent(path) +
                          (browseMode === 'dir' ? '&all=1' : ''));
    if (!r.ok) { alert('无法打开目录'); return; }
    const d = await r.json();
    browseCwd = d.cwd;
    try { localStorage.setItem('vp_browse_cwd', d.cwd === '此电脑' ? '' : d.cwd); } catch (e) {}
    $('dlgPath').textContent = '📁 ' + d.cwd;
    const pick = $('dlgPickDir');
    if (pick) pick.style.display = browseMode === 'dir' ? '' : 'none';
    const items = [];
    if (d.parent !== '') {
      items.push(`<div class="fitem dir" onclick="loadBrowse(${JSON.stringify(d.parent).replace(/"/g, '&quot;')})">⬆ ..</div>`);
    } else if (d.cwd !== '此电脑') {
      items.push(`<div class="fitem dir" onclick="loadBrowse('此电脑')">🖥 此电脑</div>`);
    }
    for (const dir of d.dirs) {
      const target = _joinBrowse(d.cwd, dir);
      items.push(`<div class="fitem dir" onclick="loadBrowse(${JSON.stringify(target).replace(/"/g, '&quot;')})">📂 ${esc(dir)}${browseMode === 'dir' ? ` <span class="hint">${t('c.dirHint', '双击进入 · 或到上级选此目录')}</span>` : ''}</div>`);
    }
    if (browseMode === 'file') {
      for (const f of d.files) {
        const full = _joinBrowse(d.cwd, f.name);
        items.push(`<div class="fitem" onclick="pickFile(${JSON.stringify(full).replace(/"/g, '&quot;')})">
          🧬 ${esc(f.name)} <span class="hint">${esc(f.size)}</span></div>`);
      }
    } else if (d.files.length) {
      // 目录模式：文件仅作内容预览（不可选，选中的是目录本身）
      for (const f of d.files) {
        items.push(`<div class="fitem" style="cursor:default">
          📄 ${esc(f.name)} <span class="hint">${esc(f.size)}</span></div>`);
      }
    }
    if (browseMode === 'dir' && !d.dirs.length && !d.files.length) {
      items.push(`<p class="hint">${t('c.emptyDir', '（空目录）')}</p>`);
    }
    $('dlgList').innerHTML = items.join('') || `<p class="hint">${t('c.noMatchFile')}</p>`;
  } catch (e) {
    setConnBanner(true);
  }
}

function pickFile(fullPath) {
  if (browseTarget && $(browseTarget)) $(browseTarget).value = fullPath;
  closeDlg();
}

function closeDlg() { $('dlgMask').style.display = 'none'; }

// ---------------- 管道（模块化逐级运行） ----------------
let curSample = null;
let curStages = [];          // 当前样品的管道状态（runUpTo 计算用）
let SERVER_DEFAULTS = {};    // 设置页的分析默认参数（渲染初值用）

async function loadDefaults() {
  try {
    const d = await (await fetch('/api/settings')).json();
    SERVER_DEFAULTS = d.defaults || {};
  } catch (e) { SERVER_DEFAULTS = {}; }
}

function defVal(key, fallback) {
  const v = SERVER_DEFAULTS[key];
  return (v === undefined || v === null || v === '') ? fallback : v;
}

async function loadSamples() {
  try {
    const r = await fetch('/api/samples');
    const list = await r.json();
    const box = $('samples');
    if (!box) return;
    if (!list.length) {
      box.innerHTML = `<p class="hint">${t('pp.noSamples')}</p>`;
      return;
    }
    // 项目筛选下拉（去重）
    const sel = $('projFilter');
    if (sel) {
      const projs = [...new Set(list.map(s => s.project).filter(Boolean))];
      const cur = sel.value;
      sel.innerHTML = `<option value="">${t('pp.allProjects', '全部项目')}</option>` +
        projs.map(p => `<option value="${esc(p)}">${esc(p)}</option>`).join('');
      if (projs.includes(cur)) sel.value = cur;
    }
    const filt = sel ? sel.value : '';
    const shown = filt ? list.filter(s => s.project === filt) : list;
    box.innerHTML = shown.length ? shown.map(s => {
      const pct = Math.round(s.done / (s.total || 7) * 100);
      return `<div class="sample-item ${s.name === curSample ? 'active' : ''}"
            onclick="selectSample('${esc(s.name)}')">
        <button class="btn small si-del" title="${t('rs.confirmClear', '清除结果（保留样品与输入信息）')}"
                onclick="event.stopPropagation();clearSample('${esc(s.name)}', this)">🧹</button>
        <button class="btn small danger si-del" title="${t('rs.del', '删除样品及其全部文件')}"
                onclick="event.stopPropagation();deleteSample('${esc(s.name)}', this)">🗑</button>
        <div class="nm"><span>${esc(s.name)}</span><span class="pct">${pct}%</span></div>
        ${s.project ? `<div class="hint" style="margin:2px 0 0">🏷 ${esc(s.project)}</div>` : ''}
        <div class="pr-bar"><div style="width:${pct}%"></div></div>
        <div class="pr">${t('pp.progress')} ${s.done}/${s.total} ${t('pp.steps')}</div></div>`;
    }).join('') : `<p class="hint">${t('pp.noSamples')}</p>`;
  } catch (e) { setConnBanner(true); }
}

// ---------------- 批量导入 + 批处理队列 ----------------
async function batchCreate(btn) {
  return withBtn(btn, async () => {
    const errBox = $('bErr');
    errBox.textContent = '';
    const project = ($('batchProject')?.value || '').trim();
    const lines = ($('batchText')?.value || '').split(/\r?\n/)
      .map(l => l.trim()).filter(l => l && !l.startsWith('#'));
    if (!lines.length) { errBox.textContent = t('pp.needRows', '请至少粘贴一行样品信息'); return; }
    const created = [];
    for (const line of lines) {
      const parts = line.split(/\t| {2,}|,/).map(s => s.trim());
      const [sample, r1, r2] = parts;
      if (!sample || !r1) { errBox.textContent = `${t('pp.badRow', '格式错误（需要 样品名+R1）')}: ${line}`; return; }
      const res = await apiCreateSample({ sample, r1, r2: r2 || null, project });
      if (!res.ok) {
        if (res.conn) setConnBanner(true);
        errBox.textContent = `${sample}: ${res.error}`; return;
      }
      created.push(sample);
    }
    $('batchText').value = '';
    toast(t('pp.batchCreate'), `${created.length} ${t('dl.samples')}`, {ttl: 4000});
    loadSamples();
    if ($('batchEnq')?.checked && created.length) {
      await fetch('/api/queue/add', {
        method: 'POST', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({ samples: created, project })});
      loadQueue();
    }
  }, '⏳ 批量创建中…');
}

const Q_ST = { queued: '排队中', running: '运行中', done: '已完成',
               failed: '失败', cancelled: '已取消' };
async function loadQueue() {
  const box = $('queueBox');
  if (!box) return;
  try {
    const q = await (await fetch('/api/queue')).json();
    if (!q.items || !q.items.length) {
      box.innerHTML = `<p class="hint">${t('pp.queueEmpty')}</p>`;
      return;
    }
    box.innerHTML = '<table class="tbl dl-mini"><tr>' +
      `<th>${t('pp.qSample')}</th><th>${t('pp.qStatus')}</th><th>${t('pp.qOps')}</th></tr>` +
      q.items.map(it => {
        const st = Q_ST[it.status] || it.status;
        const statCls = it.status === 'done' ? 'done'
          : it.status === 'failed' ? 'failed'
          : it.status === 'running' ? 'running' : 'queued';
        return `<tr><td><b>${esc(it.sample)}</b>` +
          (it.project ? ` <span class="hint" style="margin:0">🏷 ${esc(it.project)}</span>` : '') +
          `</td><td><span class="tstat ${statCls}">${st}</span>` +
          (it.error ? `<div class="hint" style="color:#b03a2e;margin:2px 0 0">${esc(it.error)}</div>` : '') +
          `</td><td>${it.status === 'queued'
            ? `<button class="btn small danger" onclick="queueRemove('${esc(it.id)}')">✖</button>` : ''}</td></tr>`;
      }).join('') + '</table>';
  } catch (e) { /* 静默，下一次轮询重试 */ }
}

async function queueRemove(id) {
  await fetch(`/api/queue/${encodeURIComponent(id)}/remove`, {method: 'POST'});
  loadQueue();
}

async function queueClear() {
  await fetch('/api/queue/clear_finished', {method: 'POST'});
  loadQueue();
}

async function selectSample(name) {
  curSample = name;
  try { localStorage.setItem('vp_last_sample', name); } catch (e) {}
  loadSamples();
  try {
    const r = await fetch('/api/pipeline/' + encodeURIComponent(name) +
                          '?lang=' + encodeURIComponent(VP_LANG));
    if (!r.ok) { alert(t('c.loadFail', '读取失败')); return; }
    renderPipe(await r.json());
  } catch (e) { setConnBanner(true); }
}

/* 管道 DAG 视图：③组装后的下游分支互不依赖（依赖关系见 vp/pipeline.py STAGE_REGISTRY），
   并排展示为紧凑节点，点击节点展开完整卡片（参数/产物/日志） */
const PIPE_FANOUT = ['hostana', 'orf', 'orfa', 'phylo', 'primer', 'gbdraw'];

/* 分析模板：一键按预设组合运行（stages=null 表示全部可用阶段；
   exclude 表示「全部可用阶段里剔除这些」） */
const PIPE_TEMPLATES = [
  { id: 'fast', label: '⚡ 快速筛查',
    tip: '质控 → 宿主去除 → 病毒筛查 → 报告（最快出结果）',
    stages: ['fastp', 'host', 'virus', 'report'] },
  { id: 'std',  label: '🎯 标准分析',
    tip: '质控 → 转换 → 宿主去除 → 病毒筛查 → 组装 → 验证 → ORF → 进化树 → 报告',
    stages: ['fastp', 'fq2fa', 'host', 'virus', 'assembly', 'verify', 'orf', 'phylo', 'report'] },
  { id: 'full', label: '🔬 完整注释',
    tip: '全部可用阶段（含宿主预测 / 功能注释 / 引物设计 / 基因组图），不含子采样',
    stages: null, exclude: ['subsample'] },
];

function runTemplate(sample, id) {
  const tpl = PIPE_TEMPLATES.find(x => x.id === id);
  if (!tpl) return;
  let usable = curStages.filter(s => s.status !== 'unavailable').map(s => s.stage);
  if (tpl.exclude) usable = usable.filter(k => !tpl.exclude.includes(k));
  const stages = tpl.stages ? tpl.stages.filter(k => usable.includes(k)) : usable;
  if (!stages.length) { alert(t('pp.tplNone', '模板所需阶段当前不可用')); return; }
  runStages(sample, stages);
}

let curBranch = '';   // 当前展开的分支节点（fan-out 详情卡）
let cardIdx = 0;      // 阶段卡片序号（renderPipe 每次渲染前重置）

function toggleBranch(stage) {
  curBranch = curBranch === stage ? '' : stage;
  document.querySelectorAll('#pipe .pipe-branch-details > .stage-card').forEach(el => {
    el.style.display = el.dataset.stage === curBranch ? '' : 'none';
  });
  document.querySelectorAll('#pipe .pipe-node').forEach(el =>
    el.classList.toggle('on', el.dataset.stage === curBranch));
}

function renderPipe(d) {
  const stMap = {done: t('pp.done'), ready: t('pp.ready'),
                 blocked: t('pp.blocked'), unavailable: t('pp.unavailable'),
                 running: t('pp.running'), skipped: t('pp.skipped')};
  curStages = d.stages;
  cardIdx = 0;                     // 每次渲染重置序号
  const byKey = {}; d.stages.forEach(s => byKey[s.stage] = s);
  const saved = collectParams();   // 重渲染前保存参数值
  const gmap = {};
  (d.groups || []).forEach(([g, ks]) => ks.forEach(k => gmap[k] = g));

  const card = (s, hidden) => {
    const i = cardIdx++;
    const canRun = s.status === 'ready' || s.status === 'done'
                || s.status === 'skipped';
    const isFresh = s.status === 'ready' || s.status === 'skipped';
    const btns = canRun ? `
      <button class="btn small ${isFresh ? 'primary' : ''}"
              onclick="runStages('${esc(d.sample)}', ['${s.stage}'])">
        ${s.status === 'done' ? t('pp.rerun') : t('pp.runThis')}</button>
      <button class="btn small" onclick="runUpTo('${esc(d.sample)}', '${s.stage}')">
        ${t('pp.runTo')}</button>` :
      (s.status === 'unavailable' ?
        `<span class="hint" style="margin:0">${t('pp.installHint')}</span>` : '');
    const viewBtn = s.view
      ? ` <button class="btn small" onclick="previewStageFile('${esc(d.sample)}','${esc(s.dir)}','${esc(s.view.split('/').pop())}','${esc(s.stage)}')">📊 ${t('c.view')}</button>` : '';
    const files = (s.outputs || []).map(o =>
      `<span>📄 ${esc(o.name.split('/').pop())} · ${fmtSize(o.size)}</span>`).join('');
    const params = renderStageParams(s.stage, saved);
    return `<div class="stage-card ${s.status}" data-stage="${s.stage}"${hidden ? ' style="display:none"' : ''}>
      <div class="stage-head"><span class="no">${i + 1}</span>
        <span class="nm">${esc(s.name)}</span>
        <span class="sbadge ${s.status}">${stMap[s.status] || s.status}</span></div>
      <div class="stage-sum">${esc(s.summary || t('pp.notRun'))}</div>
      ${files ? `<div class="stage-files">${files}</div>` : ''}
      ${params}
      <div class="stage-btns">${btns}${viewBtn}</div>
      <div class="stage-log" id="log-${s.stage}"></div></div>`;
  };
  const arrow = () => `<div class="pipe-arrow">↓</div>`;

  // 依赖分层：预处理链 → ②病毒筛查 → ③组装 → 下游并行分支 → ⑩报告；
  // 未来新增的未知阶段兜底追加在报告前
  const pre = ((d.groups || [])[0] || [null, []])[1].filter(k => byKey[k]).map(k => byKey[k]);
  const fanout = PIPE_FANOUT.map(k => byKey[k]).filter(Boolean);
  const known = new Set([...pre.map(s => s.stage), 'virus', 'assembly', 'report', ...PIPE_FANOUT]);
  const extras = d.stages.filter(s => !known.has(s.stage));

  let html = `
    <div style="display:flex;justify-content:space-between;align-items:center;
                flex-wrap:wrap;gap:8px;margin-bottom:8px">
      <div><b style="font-size:17px;color:var(--green-900)">🧪 ${esc(d.sample)}</b>
        <div class="hint" style="margin:4px 0 0">${t('pp.input')}${esc(d.r1 || t('pp.none'))}
          ${d.r2 ? ' ＋ ' + esc(d.r2) : ` ${t('pp.single')}`}</div></div>
      <button class="btn primary" onclick="runStages('${esc(d.sample)}', null)">
        ${t('pp.runAll')}</button>
    </div>
    <div class="pipe-templates">
      <span class="hint" style="margin:0">${t('pp.tplHint', '分析模板：')}</span>
      ${PIPE_TEMPLATES.map(tpl =>
        `<button class="btn small" title="${esc(tpl.tip)}"
                 onclick="runTemplate('${esc(d.sample)}', '${tpl.id}')">${esc(tpl.label)}</button>`).join('')}
    </div>`;

  if (pre.length) {
    html += `<div class="group-title">${esc(gmap[pre[0].stage] || '')}</div>` +
      pre.map(s => card(s)).join(arrow()) + arrow();
  }
  if (byKey.virus) {
    html += `<div class="group-title">${esc(gmap.virus || '')}</div>` +
      card(byKey.virus) + arrow();
  }
  if (byKey.assembly) {
    html += `<div class="group-title">${esc(gmap.assembly || '')}</div>` + card(byKey.assembly);
  }
  if (fanout.length) {
    html += `<div class="group-title">${esc(gmap.orf || '')} ·
        ${t('pp.branchHint', '并行分支（互不依赖，点节点展开参数与日志）')}</div>
      <div class="pipe-branch-row${fanout.length > 1 ? ' fanned' : ''}">
        ${fanout.map(s => `
          <div class="pipe-node ${s.status}" data-stage="${s.stage}"
               onclick="toggleBranch('${s.stage}')">
            <b>${esc(s.name)}</b>
            <span class="sbadge ${s.status}">${stMap[s.status] || s.status}</span>
            ${(s.status === 'ready' || s.status === 'done' ||
               s.status === 'skipped') ? `
              <button class="btn small" title="${t('pp.runThis')}"
                onclick="event.stopPropagation();if(curBranch!=='${s.stage}')toggleBranch('${s.stage}');
                         runStages('${esc(d.sample)}', ['${s.stage}'])">▶</button>` : ''}
          </div>`).join('')}
      </div>
      <div class="pipe-branch-details">
        ${fanout.map(s => card(s, true)).join('')}
      </div>` + arrow();
  }
  if (byKey.report) {
    html += `<div class="group-title">${esc(gmap.report || '')}</div>` + card(byKey.report);
  }
  extras.forEach(s => { html += arrow() + card(s); });
  $('pipe').innerHTML = html;
}

// 各阶段专属参数（渲染进对应卡片；label 为 i18n 键 pp.<key>）
const STAGE_PARAMS = {
  subsample: [
    { id: 'subsample', label: 'pp.subPairs', type: 'number', def: 100000 },
  ],
  fastp: [
    { id: 'fastp_dedup', label: 'pp.fastp_dedup', type: 'chk', def: false },
  ],
  fq2fa: [
    { id: 'do_fq2fa', label: 'pp.do_fq2fa', type: 'chk', def: true },
  ],
  host: [
    { id: 'db_host', label: 'pp.db_host', type: 'dir', def: '' },
  ],
  virus: [
    { id: 'db_virus', label: 'pp.db_virus', type: 'dir', def: '' },
  ],
  assembly: [
    { id: 'assembly_mode', label: 'pp.assembly_mode', type: 'select', def: 'metaviral',
      opts: [['metaviral', 'metaviral'], ['rna', 'rna'], ['meta', 'meta'],
             ['isolate', 'isolate']] },
    { id: 'assembly_input', label: 'pp.assembly_input', type: 'select', def: 'virus',
      opts: [['virus', '② virus'], ['kept', '① host-free'], ['raw', 'raw']] },
    { id: 'memory', label: 'pp.memory', type: 'number', def: 64 },
    { id: 'min_contig_len', label: 'pp.min_contig_len', type: 'number', def: 200 },
  ],
  orf: [
    { id: 'min_orf_aa', label: 'pp.min_orf_aa', type: 'number', def: 100 },
  ],
  verify: [
    { id: 'verify_host', label: 'pp.verify_host', type: 'select', def: 'all',
      opts: [['all', 'all'], ['plants', 'plants'], ['fungi', 'fungi'],
             ['bacteria', 'bacteria'], ['invertebrates', 'invertebrates'],
             ['vertebrates', 'vertebrates'], ['algae', 'algae'],
             ['protozoa', 'protozoa'], ['archaea', 'archaea']] },
    { id: 'verify_methods', label: 'pp.verify_methods', type: 'multichk',
      def: ['blastx', 'cdd'],
      opts: [['blastx', 'blastx'], ['cdd', 'CDD']] },
    { id: 'verify_combine', label: 'pp.verify_combine', type: 'select',
      def: 'union',
      opts: [['union', 'union'], ['intersect', 'intersect']] },
  ],
  phylo: [
    { id: 'do_trim', label: 'pp.do_trim', type: 'chk', def: true },
    { id: 'top_n_refs', label: 'pp.top_n_refs', type: 'number', def: 10 },
    { id: 'tree_tool', label: 'pp.tree_tool', type: 'select', def: 'fasttree',
      opts: [['fasttree', 'FastTree'], ['nj', 'NJ（快速）'], ['iqtree', 'IQ-TREE']] },
    { id: 'tree_sampling', label: 'pp.tree_sampling', type: 'select', def: 'blast',
      opts: [['blast', 'blast'], ['macro', 'macro'], ['genus', 'genus'],
             ['lineage', 'lineage']] },
    { id: 'ncbi_refs', label: 'pp.ncbi_refs', type: 'text', def: '' },
  ],
  primer: [
    { id: 'primer_mode', label: 'pp.primer_mode', type: 'select', def: 'conserved',
      opts: [['conserved', 'conserved'], ['plain', 'plain']] },
    { id: 'specificity', label: 'pp.specificity', type: 'chk', def: false },
  ],
  gbdraw: [
    { id: 'gbdraw_max', label: 'pp.gbdraw_max', type: 'number', def: 12 },
    { id: 'plot_engine', label: 'pp.plot_engine', type: 'select', def: 'auto',
      opts: [['auto', 'auto'], ['gbdraw', 'gbdraw'], ['dfv', 'DFV']] },
    { id: 'gbdraw_fasta', label: 'pp.gbdraw_fasta', type: 'file', def: '' },
    { id: 'gbdraw_ann', label: 'pp.gbdraw_ann', type: 'file', def: '' },
  ],
  report: [],
};

function renderStageParams(stage, saved) {
  const defs = STAGE_PARAMS[stage] || [];
  if (!defs.length) return '';
  const items = defs.map(p => {
    // 初值优先级：页面已填 > 设置页默认 > 代码缺省
    const dv = defVal(p.id, p.def);
    const cur = saved[p.id] !== undefined && saved[p.id] !== null ? saved[p.id] : dv;
    if (p.type === 'select') {
      const opts = p.opts.map(([v, txt]) =>
        `<option value="${v}" ${String(cur) === String(v) ? 'selected' : ''}>${txt}</option>`).join('');
      return `<div><label>${t(p.label)}</label><select id="${p.id}">${opts}</select></div>`;
    }
    if (p.type === 'dir') {
      const ph = { db_host: t('pp.dbHostPh'), db_virus: t('pp.dbVirusPh') }[p.id]
                 || t('pp.defPh');
      return `<div style="grid-column:1/-1"><label>${t(p.label)}</label>
        <div class="filerow">
          <input id="${p.id}" type="text" value="${cur}" placeholder="${ph}">
          <button class="btn small" onclick="browseDir('${p.id}')">📁</button>
        </div></div>`;
    }
    if (p.type === 'file') {
      const ph = { gbdraw_fasta: t('pp.gFastaPh'),
                   gbdraw_ann: t('pp.gAnnPh') }[p.id] || '';
      return `<div style="grid-column:1/-1"><label>${t(p.label)}</label>
        <div class="filerow">
          <input id="${p.id}" type="text" value="${cur}" placeholder="${ph}">
          <button class="btn small" onclick="browse('${p.id}')">📁</button>
        </div></div>`;
    }
    if (p.type === 'chk') {
      const on = saved[p.id] !== undefined ? saved[p.id] : !!dv;
      return `<label class="chk" style="align-self:end;margin-bottom:8px">
        <input type="checkbox" id="${p.id}" ${on ? 'checked' : ''}> ${t(p.label)}</label>`;
    }
    if (p.type === 'multichk') {
      // 多选（如验证证据方法）：name=p.id，值取勾选项；全不勾时 collectParams 回传 undefined，
      // 服务端回退默认值。
      const vals = Array.isArray(cur) ? cur.map(String)
                 : String(cur || '').split(/[,;]/).filter(Boolean);
      const boxes = p.opts.map(([v, txt]) =>
        `<label class="chk" style="margin-right:14px">
          <input type="checkbox" name="${p.id}" value="${v}"
            ${vals.includes(String(v)) ? 'checked' : ''}> ${txt}</label>`).join('');
      return `<div style="grid-column:1/-1"><label>${t(p.label)}</label>
        <div>${boxes}</div></div>`;
    }
    if (p.type === 'text') {
      return `<div style="grid-column:1/-1"><label>${t(p.label)}</label>
        <input id="${p.id}" type="text" value="${cur}" placeholder="${t('pp.ncbi_ph')}"></div>`;
    }
    return `<div><label>${t(p.label)}</label>
      <input id="${p.id}" type="number" value="${cur}"></div>`;
  }).join('');
  return `<div class="stage-params">${items}</div>`;
}

function collectParams() {
  return {
    threads: +($('threads')?.value) || undefined,
    confidence: +($('confidence')?.value) || 0,
    chunk_dir: ($('chunk_dir')?.value || '').trim() || null,
    subsample: +($('subsample')?.value) || undefined,
    db_host: ($('db_host')?.value || '').trim() || null,
    db_virus: ($('db_virus')?.value || '').trim() || null,
    assembly_mode: $('assembly_mode')?.value,
    assembly_input: $('assembly_input')?.value,
    memory: +($('memory')?.value) || undefined,
    min_contig_len: +($('min_contig_len')?.value) || undefined,
    min_orf_aa: +($('min_orf_aa')?.value) || undefined,
    do_trim: $('do_trim') === null ? undefined : !!($('do_trim')?.checked),
    top_n_refs: +($('top_n_refs')?.value) || undefined,
    tree_tool: $('tree_tool')?.value,
    tree_sampling: $('tree_sampling')?.value,
    ncbi_refs: ($('ncbi_refs')?.value || '').trim() || null,
    primer_mode: $('primer_mode')?.value,
    specificity: !!($('specificity')?.checked),
    fastp_dedup: !!($('fastp_dedup')?.checked),
    do_fq2fa: $('do_fq2fa') === null ? undefined : !!($('do_fq2fa')?.checked),
    gbdraw_max: +($('gbdraw_max')?.value) || undefined,
    plot_engine: $('plot_engine')?.value,
    gbdraw_fasta: ($('gbdraw_fasta')?.value || '').trim() || null,
    gbdraw_ann: ($('gbdraw_ann')?.value || '').trim() || null,
    force: !!($('force')?.checked),
    verify_host: $('verify_host')?.value,
    verify_methods: (() => {
      const els = document.querySelectorAll('input[name="verify_methods"]:checked');
      return els.length ? [...els].map(e => e.value) : undefined;
    })(),
    verify_combine: $('verify_combine')?.value,
  };
}

async function runStages(sample, stages) {
  const body = Object.fromEntries(
    Object.entries(collectParams()).filter(([_k, v]) => v !== undefined));
  if (stages) body.stages = stages;
  if (body.force && !confirm(t('pp.forceConfirm',
      '强制重跑会忽略断点标记，所选步骤将重新执行。继续？'))) return;
  try {
    const r = await fetch('/api/pipeline/' + encodeURIComponent(sample) + '/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body)});
    if (!r.ok) { alert(t('c.startFail', '启动失败') + ': ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    taskLogOpen.add(d.task);
    startPolling();
    const tb = $('tasks');
    if (tb) tb.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  } catch (e) { setConnBanner(true); alert(t('c.connFail', '无法连接平台服务') + ': ' + e); }
}

function runUpTo(sample, stage) {
  const effective = curStages.filter(s => s.status !== 'unavailable');
  const firstUndone = effective.findIndex(s => s.status !== 'done');
  const effOrder = effective.map(s => s.stage);
  const start = firstUndone < 0 ? effOrder.indexOf(stage) : firstUndone;
  const end = effOrder.indexOf(stage);
  const stages = effOrder.slice(Math.min(start, end), end + 1);
  if (!stages.length) { alert(t('pp.alreadyDone', '该步骤已完成（可用「重跑此步」单独重跑）')); return; }
  runStages(sample, stages);
}

async function createSample(btn) {
  return withBtn(btn, async () => {
    const r1 = $('r1').value.trim();
    const errBox = $('cErr');
    errBox.textContent = '';
    if (!r1) { errBox.textContent = t('pp.needR1', '请选择 R1 FASTQ 文件'); return; }
    const body = {
      sample: $('sample').value.trim(),
      r1, r2: $('r2').value.trim() || null,
      project: $('sampleProject')?.value.trim() || null,
      subsample: ($('subsample_on')?.checked)
        ? (+($('subsample')?.value) || defVal('subsample', 100000)) : 0,
    };
    const res = await apiCreateSample(body);
    if (!res.ok) {
      if (res.conn) setConnBanner(true);
      errBox.textContent = res.error; return;
    }
    $('r1').value = ''; $('r2').value = ''; $('sample').value = '';
    selectSample(res.data.sample);
    loadSamples();
  }, '⏳ 创建中…');
}

// ---------------- 建库 ----------------
async function buildTaxonomy(btn) {
  let r;
  try {
    r = await fetch('/api/build_taxonomy', {method: 'POST'});
  } catch (e) { alert(t('c.startFail', '启动失败') + ': ' + e); return; }
  if (!r.ok) { alert(t('c.startFail', '启动失败') + ': ' + ((await r.json().catch(() => ({}))).error || r.status)); return; }
  const d = await r.json();
  taskLogOpen.add(d.task);
  startPolling();
  watchTaskBtn(d.task, btn, '⏳ 构建中…', () => loadDbs());
}

async function buildHostDb(btn) {
  const taxid = +$('hostTaxid').value;
  if (!taxid) { alert(t('bd.needTaxid', '请填写宿主物种的 NCBI TaxID')); return; }
  const body = {
    genome: $('hostGenome').value.trim(),
    taxid,
    hash_capacity: $('hostCap').value || '256M',
    threads: +$('hostThreads').value || undefined,
    rebuild: $('hostRebuild').checked,
    clean_mid: $('hostCleanMid').checked,
  };
  let r;
  try {
    r = await fetch('/api/build_host_db', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(body)});
  } catch (e) { alert(t('c.startFail', '启动失败') + ': ' + e); return; }
  if (!r.ok) { alert(t('c.startFail', '启动失败') + ': ' + ((await r.json().catch(() => ({}))).error || r.status)); return; }
  const d = await r.json();
  taskLogOpen.add(d.task);
  startPolling();
  watchTaskBtn(d.task, btn, '⏳ 建库中…', () => loadDbs());
}

async function buildVirusDb(btn) {
  const body = {
    fasta: $('virusFasta').value.trim(),
    info: $('virusInfo').value.trim(),
    hash_capacity: $('virusCap').value || '64M',
    threads: +$('virusThreads').value || undefined,
    rebuild: $('virusRebuild').checked,
    clean_mid: !!($('virusCleanMid') && $('virusCleanMid').checked),
  };
  if (!body.fasta || !body.info) { alert(t('bd.needFastaInfo', '请填写病毒参考 FASTA 与 info 表路径')); return; }
  let r;
  try {
    r = await fetch('/api/build_virus_db', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(body)});
  } catch (e) { alert(t('c.startFail', '启动失败') + ': ' + e); return; }
  if (!r.ok) { alert(t('c.startFail', '启动失败') + ': ' + ((await r.json().catch(() => ({}))).error || r.status)); return; }
  const d = await r.json();
  taskLogOpen.add(d.task);
  startPolling();
  watchTaskBtn(d.task, btn, '⏳ 建库中…', () => loadDbs());
}

// ---------------- 结果页 ----------------
async function openDir(sample) {
  await fetch('/api/open_report_dir/' + encodeURIComponent(sample));
}

async function deleteSample(sample, btn) {
  const msg = t('pv.confirmDel', '删除样品「{name}」将连同其全部结果文件与输入档案一并移除，不可恢复。')
    .replace('{name}', sample);
  if (!confirm(msg)) return;
  btn.disabled = true;
  try {
    const r = await fetch(`/api/samples/${encodeURIComponent(sample)}/delete`,
                          {method: 'POST'});
    if (!r.ok) {
      alert(t('rs.delFail', '删除失败') + ': ' + ((await r.json()).error || ''));
      btn.disabled = false;
      return;
    }
    toast(t('rs.del', '🗑 删除') + ' · ' + sample, '', {ttl: 3000});
    btn.closest('tr')?.remove();
    btn.closest('.sample-item')?.remove();
    /* 流程页：删除的是当前选中样品时清除选中态并刷新列表 */
    if (typeof curSample !== 'undefined' && curSample === sample) {
      curSample = null;
      try { localStorage.removeItem('vp_last_sample'); } catch (e) {}
    }
    loadSamples();   /* 无样品列表容器时内部自动跳过 */
  } catch (e) { btn.disabled = false; setConnBanner(true); }
}

/* 清除结果 ≠ 删除样品：只删分析产物，保留样品与输入档案（00_prep/input.json），
   之后重跑管道即从第一步重新分析 */
async function clearSample(sample, btn) {
  const msg = t('rs.confirmClear',
    '清除样品「{name}」的全部分析结果文件？（保留样品与输入信息，可重新分析）')
    .replace('{name}', sample);
  if (!confirm(msg)) return;
  btn.disabled = true;
  try {
    const r = await fetch(`/api/samples/${encodeURIComponent(sample)}/clear`,
                          {method: 'POST'});
    if (!r.ok) {
      alert(t('rs.clearFail', '清除失败') + ': ' + ((await r.json()).error || ''));
      btn.disabled = false;
      return;
    }
    toast(t('rs.cleared', '🧹 结果已清除') + ' · ' + sample, '', {ttl: 3000});
    btn.disabled = false;
    /* 流程页：刷新列表并重载当前管道视图；结果页：整页刷新完成度 */
    if ($('pipe')) {
      loadSamples();
      if (typeof curSample !== 'undefined' && curSample === sample) selectSample(sample);
    } else {
      location.reload();
    }
  } catch (e) { btn.disabled = false; setConnBanner(true); }
}

async function toggleFiles(sample, btn) {
  const row = btn.closest('tr').nextElementSibling;
  const box = row.querySelector('.filelist');
  if (row.style.display === 'none') {
    row.style.display = '';
    if (!box.dataset.loaded) {
      const r = await fetch('/api/sample_files/' + encodeURIComponent(sample));
      const d = await r.json();
      box.innerHTML = (d.files || []).map(f =>
        `<a href="/api/samples/${encodeURIComponent(sample)}/${f.path}?dl=1" download>${esc(f.path)} (${esc(f.size)})</a>`
      ).join('') || `<span class="hint">${t('rs.noFiles')}</span>`;
      box.dataset.loaded = '1';
    }
  } else {
    row.style.display = 'none';
  }
}

// ---------------- 主页总览 ----------------
async function loadHome() {
  try {
    const [r, rkv] = await Promise.all([
      fetch('/api/dbs'),
      fetch('/api/kv_index_list').catch(() => null),
    ]);
    const dbs = await r.json();
    // 病毒鉴定库（salmon/minibwa 比对索引）不在 /api/dbs 体系内，单独探
    let kvLibs = [];
    try { kvLibs = ((await rkv.json()).libs) || []; } catch (e) {}
    // 四张库卡均指向 /build（库的查看与构建入口都在那里，与分析流程页无关）
    const defs = [
      ['0️⃣', t('hm.db.tax'), 'taxonomy', t('hm.db.tax.d'), '/build'],
      ['1️⃣', t('hm.db.host'), 'host', t('hm.db.host.d'), '/build'],
      ['🦠', t('hm.db.virus'), 'virus', t('hm.db.virus.d'), '/build'],
    ];
    let html = defs.map(([icon, name, key, desc, href]) => {
      const st = (dbs[key] || {}).ready;
      const label = st ? t('hm.ready')
        : (key === 'virus' ? t('hm.virusNotReady') : t('hm.notReady'));
      return `<a class="module-card" href="${href}">
        <span class="mc-state ${st ? 'st-ok' : 'st-no'}">${label}</span>
        <div class="mc-icon">${icon}</div>
        <div class="mc-name">${name}</div>
        <div class="mc-desc">${desc}</div></a>`;
    }).join('');
    // 第 4 张：病毒鉴定库（有库才称就绪）
    const kvSt = kvLibs.length > 0;
    const kvEng = (kvLibs[0] && kvLibs[0].minibwa && kvLibs[0].salmon) ? 'minibwa+salmon'
      : (kvLibs[0] && kvLibs[0].salmon ? 'salmon'
      : (kvLibs[0] && kvLibs[0].minibwa ? 'minibwa' : ''));
    const kvDesc = t('hm.db.kvidx.d') + (kvSt && kvEng ? ` · ${kvEng}` : '');
    html += `<a class="module-card" href="/build">
      <span class="mc-state ${kvSt ? 'st-ok' : 'st-no'}">
        ${kvSt ? t('hm.db.kvidxN') + ' (' + kvLibs.length + ')' : t('hm.db.kvidxNone')}</span>
      <div class="mc-icon">🔍</div>
      <div class="mc-name">${t('hm.db.kvidx')}</div>
      <div class="mc-desc">${kvDesc}</div></a>`;
    $('dbCards').innerHTML = html;
  } catch (e) { setConnBanner(true); }
  try {
    const r = await fetch('/api/samples');
    const list = (await r.json()).slice(0, 6);
    $('homeSamples').innerHTML = list.length ? list.map(s => {
      const pct = Math.round(s.done / (s.total || 7) * 100);
      return `<div class="module-card home-sample"
        onclick="location.href='/pipeline?sample=${encodeURIComponent(s.name)}'">
        <div class="mc-name">🧪 ${esc(s.name)}</div>
        <div class="pr-bar" style="height:6px;background:var(--line-100);border-radius:4px;overflow:hidden;margin:8px 0 4px">
          <div style="height:100%;width:${pct}%;background:linear-gradient(90deg,var(--green-600),#43b47c)"></div></div>
        <div class="mc-desc">${t('hm.progress')} ${s.done}/${s.total} ${t('pp.steps')}（${pct}%）</div></div>`;
    }).join('') : `<p class="hint">${t('hm.noSamples')}</p>`;
  } catch (e) { setConnBanner(true); }
}

// ---------------- 阶段结果预览（表格/JSON 就地弹窗，HTML 新窗口打开） ----------------
function previewStageFile(sample, stageDir, filename, stageKey) {
  // HTML 报告类（fastp/桑基/旭日/主报告）直接整页打开，体验更好
  if (filename.endsWith('.html')) {
    window.open(`/api/samples/${encodeURIComponent(sample)}/${encodeURIComponent(stageDir)}/${encodeURIComponent(filename)}`, '_blank');
    return;
  }
  const mask = document.createElement('div');
  mask.className = 'dlgmask';
  mask.style.zIndex = 120;
  mask.innerHTML = `<div class="dlg" style="width:900px">
    <div class="dlghead"><b>📄 ${esc(filename)}</b>
      <span>
        <button class="btn small" onclick="window.open('/api/samples/${encodeURIComponent(sample)}/${encodeURIComponent(stageDir)}/${encodeURIComponent(filename)}?dl=1','_blank')">⬇ 下载</button>
        <button class="btn small" onclick="this.closest('.dlgmask').remove()">✕</button>
      </span></div>
    <div class="pv-body" style="padding:6px 18px 16px;overflow:auto;max-height:70vh">
      <p class="hint">${t('c.loading')}</p></div></div>`;
  document.body.appendChild(mask);
  mask.addEventListener('click', e => { if (e.target === mask) mask.remove(); });
  const body = mask.querySelector('.pv-body');
  fetch(`/api/samples/${encodeURIComponent(sample)}/${encodeURIComponent(stageDir)}/${encodeURIComponent(filename)}`)
    .then(r => { if (!r.ok) throw new Error('HTTP ' + r.status); return r.text(); })
    .then(text => {
      if (filename.endsWith('.json')) {
        let obj;
        try { obj = JSON.parse(text); } catch (e) { obj = null; }
        body.innerHTML = obj
          ? `<pre style="background:#0d1f16;color:#b8e6c9;padding:12px;border-radius:8px;font-size:11.5px;white-space:pre-wrap">${esc(JSON.stringify(obj, null, 2))}</pre>`
          : `<pre>${esc(text.slice(0, 200000))}</pre>`;
        return;
      }
      // TSV/CSV → 表格（每页 20 行分页，全量行数保留）
      const rows = text.split(/\r?\n/).filter(l => l.trim());
      if (!rows.length) { body.innerHTML = `<p class="hint">${t('c.noData')}</p>`; return; }
      const sep = filename.endsWith('.csv') ? ',' : '\t';
      pvCells = rows.map(l => l.split(sep));
      pvPage = 1;
      renderPvTable(body);
    })
    .catch(e => { body.innerHTML = `<p class="err">${t('c.loadFail', '读取失败')}: ${esc(e)}</p>`; });
}

// 阶段文件预览表格分页 + 排序
let pvCells = null, pvPage = 1, pvSortI = -1, pvSortDir = 1;
function pvGo(p) {
  pvPage = p;
  const body = document.querySelector('.pv-body');
  if (body) renderPvTable(body);
}
function pvSortBy(i) {
  if (pvSortI === i) pvSortDir = -pvSortDir;
  else { pvSortI = i; pvSortDir = 1; }
  const body = document.querySelector('.pv-body');
  if (body) renderPvTable(body);
}
function renderPvTable(body) {
  if (!pvCells || !pvCells.length) return;
  const ncol = Math.max(...pvCells.map(c => c.length));
  let dataRows = pvCells.slice(1);
  if (pvSortI >= 0) {
    dataRows.sort((a, b) => {
      const va = a[pvSortI] || '', vb = b[pvSortI] || '';
      const na = parseFloat(va), nb = parseFloat(vb);
      const c = (!isNaN(na) && !isNaN(nb)) ? na - nb
        : String(va).localeCompare(String(vb), 'zh');
      return c * pvSortDir;
    });
  }
  const pages = Math.max(1, Math.ceil(dataRows.length / PER_PAGE));
  pvPage = Math.min(Math.max(1, pvPage), pages);
  const rows = dataRows.slice((pvPage - 1) * PER_PAGE, pvPage * PER_PAGE);
  let html = '<table class="tbl"><thead><tr>';
  for (let i = 0; i < ncol; i++) {
    const arrow = pvSortI === i ? (pvSortDir === 1 ? ' ▲' : ' ▼') : '';
    html += `<th class="sortable" onclick="pvSortBy(${i})">${esc((pvCells[0][i] || '').slice(0, 40))}${arrow}</th>`;
  }
  html += '</tr></thead><tbody>';
  for (const cells of rows) {
    html += '<tr>';
    for (let i = 0; i < ncol; i++) {
      const v = (cells[i] || '');
      html += `<td title="${esc(v.slice(0, 300))}">${esc(v.length > 60 ? v.slice(0, 57) + '…' : v)}</td>`;
    }
    html += '</tr>';
  }
  html += '</tbody></table>';
  html += `<p class="hint">共 ${dataRows.toLocaleString()} 行（不含表头）· 每页 ${PER_PAGE} 行</p>`
    + pagerHtml(pvPage, pages, 'pvGo');
  body.innerHTML = html;
}

// ---------------- 工具注册表 + 常用收藏（跨页共享） ----------------
const TOOL_REGISTRY = [
  { id: 't-convert',  href: '/tools#t-convert',  ic: '🔄', zh: '⓪ 格式转换',         en: '⓪ Format convert' },
  { id: 't-fastp',    href: '/tools#t-fastp',    ic: '🧹', zh: '① 序列质控',         en: '① Sequence QC' },
  { id: 't-identify', href: '/tools#t-identify', ic: '🦠', zh: '② 病毒识别和分类',   en: '② Virus classify' },
  { id: 't-assemble', href: '/tools#t-assemble', ic: '🧩', zh: '③ 病毒组装',         en: '③ Assembly' },
  { id: 't-contigs',  href: '/tools#t-contigs',  ic: '🔎', zh: '④ 病毒 contig 深度分析', en: '④ Contig deep-dive' },
  { id: 't-ncbi',     href: '/tools#t-ncbi',     ic: '⬇', zh: '⑤ 参考序列下载',    en: '⑤ NCBI references' },
  { id: 't-synteny',  href: '/tools#t-synteny',  ic: '🧬', zh: '⑥ 同属共线性比较',   en: '⑥ Synteny' },
  { id: 't-msa',      href: '/tools#t-msa',      ic: '🔤', zh: '⑦ 多序列比对/MSA 查看', en: '⑦ MSA viewer' },
  { id: 't-tree',     href: '/tools#t-tree',     ic: '🌳', zh: '⑧ 进化树查看器',     en: '⑧ Tree viewer' },
  { id: 't-sdt',      href: '/tools#t-sdt',      ic: '📐', zh: '⑨ SDT 分析和绘制',   en: '⑨ SDT matrix' },
  { id: 't-seqview',  href: '/results#seqview',  ic: '📄', zh: '序列查看器',         en: 'Sequence viewer' },
  { id: 'dl',         href: '/download',          ic: '📥', zh: '公共数据下载',       en: 'Public data downloads' },
];

function getFavs() {
  try { return JSON.parse(localStorage.getItem('vp_fav_tools') || '[]'); }
  catch (e) { return []; }
}

function toggleFav(id) {
  const favs = getFavs();
  const i = favs.indexOf(id);
  if (i >= 0) favs.splice(i, 1); else favs.push(id);
  try { localStorage.setItem('vp_fav_tools', JSON.stringify(favs)); } catch (e) {}
  document.querySelectorAll(`.starbtn[data-tool="${id}"]`)
    .forEach(b => b.classList.toggle('on', favs.includes(id)));
  if (typeof renderFavStrip === 'function') renderFavStrip();
  return favs.includes(id);
}

function toolName(tool) {
  return VP_LANG === 'en' ? tool.en : tool.zh;
}

function renderFavStrip() {
  const box = $('favStrip');
  if (!box) return;
  const favs = getFavs();
  const section = $('favSection');
  if (!favs.length) {
    box.innerHTML = '';
    if (section) section.style.display = 'none';
    return;
  }
  if (section) section.style.display = '';
  box.innerHTML = favs.map(id => {
    const tool = TOOL_REGISTRY.find(x => x.id === id);
    return tool ? `<a class="fav-chip" href="${tool.href}"><span class="ic">${tool.ic}</span>${esc(toolName(tool))}</a>` : '';
  }).join('');
}

// ---------------- 深色主题 ----------------
function applyTheme(theme) {
  if (theme === 'dark') document.documentElement.setAttribute('data-theme', 'dark');
  else document.documentElement.removeAttribute('data-theme');
  try { localStorage.setItem('vp_theme', theme); } catch (e) {}
}

// ---------------- 拖拽文件 → 上传 → 填路径 ----------------
function wireDropzones() {
  document.querySelectorAll('.filerow').forEach(row => {
    if (row.dataset.dropWired) return;
    const input = row.querySelector('input[type=text]');
    if (!input) return;
    row.classList.add('drop-able');
    row.dataset.dropWired = '1';
    row.addEventListener('dragover', e => {
      if ([...(e.dataTransfer.types || [])].includes('Files')) {
        e.preventDefault();
        row.classList.add('dragover');
      }
    });
    row.addEventListener('dragleave', () => row.classList.remove('dragover'));
    row.addEventListener('drop', async e => {
      row.classList.remove('dragover');
      const file = [...(e.dataTransfer.files || [])][0];
      if (!file) return;                       // 非文件拖放（如文本）不接管
      e.preventDefault();
      const old = input.placeholder;
      input.placeholder = t('dz.uploading', '⬆ 上传中…') + ' ' + file.name;
      try {
        const fd = new FormData();
        fd.append('file', file);
        const r = await fetch('/api/upload', { method: 'POST', body: fd });
        const d = await r.json();
        if (!r.ok) throw new Error(d.error || 'upload failed');
        input.value = d.path;
        input.placeholder = old;
        toast(t('dz.ok', '已导入上传文件'), `${file.name} → ${d.path}`, {ttl: 4000});
      } catch (err) {
        input.placeholder = old;
        toast(t('dz.fail', '上传失败'), String(err), {kind: 'failed', ttl: 6000});
      }
    });
  });
}

// ---------------- 序列查看器（专项分析页） ----------------
let svPath = '', svPage = 0;
async function seqviewLoad(page) {
  const p = ($('svPath')?.value || '').trim();
  const box = $('svBox');
  if (!p) { box.innerHTML = `<p class="hint">${t('sv.needPath', '请先选择或拖入 FASTA 文件')}</p>`; return; }
  if (p !== svPath) svPage = 0;
  svPage = Math.max(0, page || 0);
  box.innerHTML = `<p class="hint">${t('c.loading')}</p>`;
  try {
    const r = await fetch(`/api/seqview?path=${encodeURIComponent(p)}&page=${svPage}`);
    const d = await r.json();
    if (!r.ok) { box.innerHTML = `<p class="hint" style="color:#b91c1c">${esc(d.error || '读取失败')}</p>`; return; }
    svPath = p;
    $('svMeta').textContent =
      `${d.total} ${t('sv.records', '条')} · ${fmtSize(d.total_bp)}` +
      (d.truncated ? ` · ${t('sv.capped', '仅统计前 20000 条')}` : '');
    if (!d.rows || !d.rows.length) {
      box.innerHTML = `<p class="hint">${t('c.noData')}</p>`;
      $('svPager').style.display = 'none';
      return;
    }
    box.innerHTML = '<table class="tbl"><tr>' +
      `<th>ID</th><th>${t('sv.len', '长度')}</th><th>GC%</th><th>${t('sv.deg', '简并%')}</th><th>${t('sv.preview', '序列预览')} (300bp)</th></tr>` +
      d.rows.map(x => `<tr><td class="mono" title="${esc(x.id)}">${esc(x.id.length > 42 ? x.id.slice(0, 40) + '…' : x.id)}</td>` +
        `<td class="mono">${x.len.toLocaleString()}</td><td class="mono">${x.gc}</td><td class="mono">${x.deg}</td>` +
        `<td><span class="seqview-preview" title="${esc(x.preview || '')}">${esc(x.preview || '')}</span></td></tr>`).join('') +
      '</table>';
    $('svPager').style.display = 'flex';
    $('svPageInfo').textContent = `${t('pp.progress').split(' ')[0] || ''} ${svPage + 1} / ${d.pages}`;
  } catch (e) { box.innerHTML = `<p class="hint" style="color:#b91c1c">${t('c.loadFail', '读取失败')}: ${esc(e)}</p>`; }
}

// ---------------- 初始化 ----------------
window.addEventListener('DOMContentLoaded', () => {
  const q = new URLSearchParams(location.search);
  // 建库页
  if ($('hostTaxid') !== null) loadDbs();
  // 主页
  if ($('dbCards')) loadHome();
  // 分析管道页
  if ($('samples')) {
    loadSamples();
    loadQueue();
    setInterval(loadQueue, 4000);
    loadDefaults().then(() => {
      // 全局参数初值应用设置页默认
      const cf = $('confidence');
      if (cf && !cf.value && SERVER_DEFAULTS.confidence) cf.value = SERVER_DEFAULTS.confidence;
      if (q.get('sample')) selectSample(q.get('sample'));
    });
    if (q.get('focus')) setTimeout(() => {
      const el = $('sample'); if (el) { el.focus(); el.scrollIntoView({ behavior: 'smooth' }); }
    }, 200);
    const sub = $('subsample_on');
    if (sub) sub.addEventListener('change', () => {
      $('sub_row').style.display = sub.checked ? '' : 'none';
    });
    if (location.hash.startsWith('#grp-')) {
      setTimeout(() => {
        const g = $(location.hash.slice(1));
        if (g) g.scrollIntoView({ behavior: 'smooth', block: 'start' });
      }, 800);
    }
  }
  if ($('tasks')) startPolling();
  // 全页通用：拖拽接线 / 收藏星 / 常用条
  wireDropzones();
  wireStars();
  renderFavStrip();
});

// ---------------- 收藏星渲染（专项分析各工具卡标题） ----------------
function wireStars() {
  document.querySelectorAll('section.card > h2[id]').forEach(h2 => {
    const id = h2.id;
    if (!id.startsWith('t-') || h2.querySelector('.starbtn')) return;
    const btn = document.createElement('button');
    btn.className = 'starbtn';
    btn.type = 'button';
    btn.dataset.tool = id;
    btn.title = t('fav.toggle', '加入/移出常用工具');
    btn.textContent = '☆';
    btn.classList.toggle('on', getFavs().includes(id));
    btn.onclick = () => {
      const on = toggleFav(id);
      btn.textContent = on ? '★' : '☆';
      toast(t('fav.strip', '常用工具'),
            on ? t('fav.added', '★ 已加入常用（总览页可见）')
               : t('fav.removed', '已移出常用'), {ttl: 2000});
    };
    h2.appendChild(btn);
  });
}

// ---------------- 比较基因组：⑤ 参考序列下载 / ⑥ GenBank 集合与共线性 ----------------
/* 原 tools.html 内联脚本迁移至此；输入取值用 _v（避免与页面内联 val 重名）。 */
function _v(id) { return ($(id)?.value || '').trim(); }

function _watchTask(tid, onDone, btn, busyText) {
  const unlock = lockBtn(btn, busyText || '⏳ 运行中…');
  taskLogOpen.add(tid);
  startPolling();
  const timer = setInterval(async () => {
    try {
      const snap = await (await fetch('/api/task/' + tid)).json();
      if (['done', 'failed', 'cancelled'].includes(snap.status)) {
        clearInterval(timer);
        unlock();
        if (onDone) onDone(snap);
      }
    } catch (e) { /* 轮询继续 */ }
  }, 5000);
}

async function ncbiSearch() {
  const term = _v('n_term'), db = $('n_db').value;
  if (!term) { alert('请输入检索式'); return; }
  const box = $('ncbiPreview');
  box.innerHTML = '<p class="hint">检索中…</p>';
  try {
    const r = await fetch('/api/ncbi/search', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ term, db, limit: 50 })});
    if (!r.ok) { box.innerHTML = '<p class="hint" style="color:#b91c1c">检索失败: ' + esc((await r.json()).error || '') + '</p>'; return; }
    const d = await r.json();
    const head = '<tr><th>Accession</th><th>描述</th><th>物种</th><th>长度</th><th>更新</th></tr>';
    box.innerHTML = `<p class="hint">命中 <b>${d.count}</b> 条（预览前 ${d.rows.length} 条）</p>` +
      '<table class="table" style="width:100%;border-collapse:collapse">' + head +
      d.rows.map(x => `<tr><td>${esc(x.acc)}</td>` +
        `<td title="${esc(x.title)}">${esc((x.title || '').slice(0, 60))}</td>` +
        `<td>${esc(x.organism)}</td><td>${esc(x.length)}</td><td>${esc(x.updated)}</td></tr>`).join('') +
      '</table>';
  } catch (e) { box.innerHTML = '<p class="hint" style="color:#b91c1c">无法连接: ' + esc(e) + '</p>'; }
}

async function ncbiDownload(btn) {
  const term = _v('n_term'), name = _v('n_name');
  if (!term || !name) { alert('检索式与集合名均为必填'); return; }
  const unlock = lockBtn(btn, '⏳ 下载中…');
  try {
    const r = await fetch('/api/ncbi/download', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ term, name, db: $('n_db').value,
                             max: +_v('n_max') || 100 })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    taskLogOpen.add(d.task);
    startPolling();
    const timer = setInterval(async () => {
      try {
        const snap = await (await fetch('/api/task/' + d.task)).json();
        if (['done', 'failed', 'cancelled'].includes(snap.status)) {
          clearInterval(timer); unlock(); loadNcbiCollections();
        }
      } catch (e) { /* 轮询继续 */ }
    }, 5000);
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

async function loadNcbiCollections() {
  const box = $('ncbiCollections');
  if (!box) return;
  try {
    const cols = await (await fetch('/api/ncbi/collections')).json();
    box.innerHTML = cols.length
      ? '<table class="table" style="width:100%;border-collapse:collapse"><tr><th>集合</th><th>序列数</th><th>库</th><th>检索式</th><th>日期</th></tr>' +
        cols.map(c => `<tr><td><b>${esc(c.name)}</b></td><td>${c.n_seqs}</td><td>${esc(c.db)}</td>` +
          `<td title="${esc(c.query)}">${esc((c.query || '').slice(0, 50))}</td><td>${esc(c.date)}</td></tr>`).join('') + '</table>' +
        '<p class="hint">分析管道 → 「NCBI 参考集合」填集合名即可把这些参考追加进进化树比对。</p>'
      : '<p class="hint">（暂无集合，先搜索并下载）</p>';
  } catch (e) { box.innerHTML = '<p class="hint">加载失败</p>'; }
}

async function gbDownload(btn) {
  const name = _v('s_name'), term = _v('s_term'), accs = _v('s_acc');
  if (!name || (!term && !accs) || (term && accs)) {
    alert('集合名必填；检索式与 accession 二选一'); return;
  }
  try {
    const r = await fetch('/api/gb/download', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name, term, accessions: accs,
                             max: +_v('s_max') || 50 })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    _watchTask(d.task, loadGbCollections, btn, '⏳ 下载中…');
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

async function gbImport(btn) {
  const name = _v('s_name'), files = _v('s_files');
  if (!name || !files) { alert('集合名与本机 .gb 路径必填'); return; }
  try {
    const r = await fetch('/api/gb/import', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name,
                             files: files.split(',').map(s => s.trim()).filter(Boolean) })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    _watchTask(d.task, loadGbCollections, btn, '⏳ 导入中…');
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

async function runCompareFor(name, btn) {
  try {
    const r = await fetch('/api/compare/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ collection: name,
                             min_ident: parseFloat(_v('s_min_ident')) || 0.30,
                             min_cov: parseFloat(_v('s_min_cov')) || 0.50 })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    _watchTask(d.task, () => { loadGbCollections(); loadSyntenyResult(name); }, btn, '⏳ 比较中…');
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

/* 页内结果展示：读 /api/compare/preview 把共线性图 + 成对表 + 家族表内嵌到卡内 */
async function loadSyntenyResult(name) {
  const box = $('syntenyResult');
  if (!box) return;
  if (!name) { box.innerHTML = '<p class="hint">选择集合点「▶ 比较」查看结果</p>'; return; }
  try {
    const r = await fetch('/api/compare/preview?name=' + encodeURIComponent(name));
    if (!r.ok) { box.innerHTML = '<p class="hint">（暂无可展示的结果，先运行比较）</p>'; return; }
    const d = await r.json();
    box.innerHTML = '<div style="max-height:600px;overflow:auto;border:1px solid #eee;border-radius:8px;padding:12px;background:#fff">' + (d.html || '<p class="hint">无内容</p>') + '</div>';
  } catch (e) { box.innerHTML = '<p class="hint">结果加载失败: ' + e + '</p>'; }
}

/* ================= CDS/PEP 提取（PhyloSuite 布局） ================= */
async function gbExtract(name, btn) {
  try {
    const r = await fetch('/api/gb/extract', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    taskLogOpen.add(d.task);
    startPolling();
    _watchTask(d.task, () => {
      loadGbCollections();
      const sel = $('exColl');
      if (sel) sel.value = name;
      gbExtractFilesLoad();
    }, btn, '⏳ 提取中…');
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

async function gbExtractFilesLoad() {
  const box = $('exFiles');
  if (!box) return;
  const name = ($('exColl') && $('exColl').value) || '';
  if (!name) { box.innerHTML = '<p class="hint">（先选择集合；未提取的集合点上方「🧬 提取」）</p>'; return; }
  try {
    const r = await fetch('/api/gb/extract_files?name=' + encodeURIComponent(name));
    if (!r.ok) {
      box.innerHTML = '<p class="hint">' + esc((await r.json()).error || '未提取') + '</p>';
      return;
    }
    const items = await r.json();
    const label = k => ({ 'genome.fa': '🧬 genome.fa（全基因组）',
                          'CDS.fa': '🧬 CDS.fa（全部 CDS 核酸）',
                          'PEP.fa': '🧬 PEP.fa（全部蛋白）',
                          'genes.tsv': '📋 genes.tsv（基因×基因组摘要）',
                          'cds_gene': '📕', 'pep_gene': '📘' }[k] || k);
    box.innerHTML = '<table class="table" style="width:100%;border-collapse:collapse">' +
      '<tr><th>产物</th><th style="width:260px">操作</th></tr>' +
      items.map(x => `<tr><td>${label(x.kind).startsWith('📕') || label(x.kind).startsWith('📘')
        ? `${x.kind === 'cds_gene' ? '📕 CDS/' : '📘 PEP/'}<b>${esc(x.path.split('/').pop())}</b>（按基因）`
        : label(x.kind)}</td>` +
        `<td>` +
        `<a class="btn small" href="${esc(x.path)}" download>⬇ 下载</a> ` +
        `<button class="btn small" onclick="exSend('${esc(x.path)}', 'align')">→ 送比对</button> ` +
        `<button class="btn small" onclick="exSend('${esc(x.path)}', 'treebuild')">→ 送建树</button>` +
        `</td></tr>`).join('') + '</table>';
  } catch (e) { box.innerHTML = '<p class="hint">无法连接: ' + esc(e) + '</p>'; }
}

function exSend(path, target) {
  if (target === 'align') {
    $('al_fa').value = path;
    location.hash = '#t-align';
  } else {
    $('qt_fa').value = path;
    location.hash = '#t-treebuild';
  }
  if (typeof renderModuleTree === 'function') renderModuleTree();
}

/* 集合列表「🌳 建树」按钮：切换到进化树构建卡并选中该集合 */
async function tbBuildFor(name, btn) {
  const sel = $('tbColl');
  if (sel) {
    await loadTbColls();
    sel.value = name;
  }
  location.hash = '#t-treebuild';
  if (typeof renderModuleTree === 'function') renderModuleTree();
  await tbBuild(btn);
}

/* 显式巡检：重解析全部 .gb，重建清单与警告（列表页只读已生成的清单） */
async function gbInspect(name, btn) {
  try {
    const r = await fetch('/api/gb/inspect', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    _watchTask(d.task, loadGbCollections, btn, '⏳ 巡检中…');
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

/* 集合建树：全基因组 MAFFT 比对 → FastTree/IQ-TREE；
   产物在结果中心 MSA / 进化树 / SDT 查看器以「🧬 集合」样品展示 */
async function gbBuildTree(name, btn) {
  const treeTool = ($('s_tree_tool')?.value) || 'fasttree';
  try {
    const r = await fetch('/api/gb/phylo', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name, tree_tool: treeTool,
                             molecule: $('tbMolecule')?.value || 'genome',
                             gene: ($('tbGene') && $('tbGene').value || '').trim() })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    _watchTask(d.task, loadGbCollections, btn, '⏳ 建树中…');
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

/* 「▶ 比较」实际使用的阈值/风格来自上方表单——在集合表上方实时提示 */
function updateGbParamHint() {
  const el = $('gbParamHint');
  if (!el) return;
  const mi = _v('s_min_ident') || '0.30', mc = _v('s_min_cov') || '0.50';
  el.innerHTML = `「▶ 比较」以 <b>LoVis4u</b> 出图（不可用时自动回退内置绘图）；` +
    `蛋白聚类参数：identity ≥ ${esc(mi)} · 双侧覆盖度 ≥ ${esc(mc)}。` +
    `建树请用「进化树构建（科/属级）」卡。`;
}

async function loadGbCollections() {
  /* 集合列表双渲染：
     #gbCollections（同属共线性比较卡）= ▶ 比较 / 🌳 建树 + 结果链接；
     #gbCollectionsMgmt（GenBank 集合管理卡）= 🔍 巡检。 */
  const cmpBox = $('gbCollections'), mgmtBox = $('gbCollectionsMgmt');
  if (!cmpBox && !mgmtBox) return;
  updateGbParamHint();
  try {
    const cols = await (await fetch('/api/gb/collections')).json();
    const gbRow = c => {
      const cds = (c.records || []).reduce((s, r) => s + (+r.cds || 0), 0);
      const src = c.source === 'query' ? esc((c.term || '').slice(0, 36))
        : (c.source === 'accessions' ? esc(c.accessions) + ' 个 accession' : '本机导入');
      const cmp = `/compare/${encodeURIComponent(c.name)}`;
      const results = [
        c.has_compare ? `<a href="${cmp}/" target="_blank">📊 共线性图</a>` : '',
        c.has_lovis4u ? `<a href="${cmp}/lovis4u.pdf" target="_blank">📄 PDF</a>` : '',
        c.has_phylo ? `<a href="/results#msa" title="结果中心查看比对/树/SDT">🌳 MSA·树</a>` : '',
      ].filter(Boolean).join(' ') || '—';
      const warnN = (c.warnings || []).length;
      const warnTip = warnN ? esc(c.warnings.join(' | ')) : '';
      return `<tr><td><b>${esc(c.name)}</b>${warnN ? ` <span class="hint" style="color:#b45309" title="${warnTip}">⚠${warnN}</span>` : ''}</td>` +
        `<td>${c.n_records}</td><td>${cds}</td>` +
        `<td>${src}</td><td>${esc(c.date)}</td><td>${results}</td>` +
        `<td style="white-space:nowrap">[[ACTIONS]]</td></tr>`;
    };
    const mkTable = rows => rows.length
      ? '<table class="table" style="width:100%;border-collapse:collapse">' +
        '<tr><th>集合</th><th>记录数</th><th>CDS 数</th><th>来源</th><th>日期</th><th>结果</th><th></th></tr>' +
        rows.join('') + '</table>'
      : '<p class="hint">（暂无集合）</p>';
    if (cmpBox) {
      cmpBox.innerHTML = mkTable(cols.map(c => gbRow(c).replace('[[ACTIONS]]',
        `<button class="btn small primary" data-name="${esc(c.name)}" onclick="runCompareFor(this.dataset.name, this)">▶ 比较</button> ` +
        `<button class="btn small" data-name="${esc(c.name)}" onclick="tbBuildFor(this.dataset.name, this)" title="到「进化树构建」卡用该集合建树">🌳 建树</button>`)));
    }
    if (mgmtBox) {
      mgmtBox.innerHTML = mkTable(cols.map(c => gbRow(c).replace('[[ACTIONS]]',
        `<button class="btn small primary" data-name="${esc(c.name)}" onclick="gbExtract(this.dataset.name, this)" title="提取 genome / CDS / PEP（分类分目录 + 按基因拆分）">🧬 提取</button> ` +
        `<button class="btn small" data-name="${esc(c.name)}" onclick="gbInspect(this.dataset.name, this)" title="重解析全部 .gb，重建清单与警告">🔍 巡检</button>`)));
    }
    const exSel = $('exColl');
    if (exSel) {
      const cur = exSel.value;
      exSel.innerHTML = '<option value="">（选择集合查看提取产物）</option>' +
        cols.map(c => `<option value="${esc(c.name)}">${esc(c.name)}${c.has_extract ? '（已提取）' : ''}</option>`).join('');
      if (cur && cols.some(c => c.name === cur)) exSel.value = cur;
    }
    // 共线性比较卡的集合选择下拉（#synColl）
    const synSel = $('synColl');
    if (synSel) {
      const curSyn = synSel.value;
      synSel.innerHTML = '<option value="">（选择集合）</option>' +
        cols.map(c => `<option value="${esc(c.name)}">${esc(c.name)}（${c.n_records} 条）</option>`).join('');
      if (curSyn && cols.some(c => c.name === curSyn)) synSel.value = curSyn;
      else if (cols.length === 1) synSel.value = cols[0].name;
    }
    const wbox = $('gbWarnings');
    if (wbox) {
      const warned = cols.filter(c => (c.warnings || []).length);
      wbox.innerHTML = warned.map(c =>
        `<details><summary class="hint" style="color:#b45309;cursor:pointer">⚠ ${esc(c.name)}：${c.warnings.length} 条巡检警告</summary>` +
        '<ul class="hint" style="margin:4px 0 8px">' +
        c.warnings.map(w => `<li>${esc(w)}</li>`).join('') + '</ul></details>').join('');
    }
  } catch (e) {
    if (cmpBox) cmpBox.innerHTML = '<p class="hint">加载失败</p>';
    if (mgmtBox) mgmtBox.innerHTML = '<p class="hint">加载失败</p>';
  }
}

['s_min_ident', 's_min_cov', 's_style', 's_lovis4u'].forEach(id => {
  const el = $(id);
  if (el) el.addEventListener('change', updateGbParamHint);
});
loadNcbiCollections();
loadGbCollections();
ictvCascadeRefetch();
loadTbColls();

/* ================= 比较基因组独立工作区：⑦ MSA / ⑧ 进化树 / ⑨ SDT ================= */
/* 三张卡片各自的输入/参数/输出都独立于结果中心；渲染器为共享实现。 */

/* ---- 样品/集合列表（MSA / 树 / SDT 卡与 results 页共用，缓存一次） ---- */
let sampleListCache = null;
async function fetchSampleList() {
  if (!sampleListCache) {
    try {
      sampleListCache = await (await fetch('/api/msa/samples')).json();
    } catch (e) { sampleListCache = []; }
  }
  return sampleListCache;
}
function sampleLabel(x) { return x.label || x.sample; }

/* ---- 共享渲染：SNP-only 热图 ---- */
const SNP_COLORS = { A: '#4daf4a', C: '#377eb8', G: '#ffb200',
                     T: '#e41a1c', U: '#e41a1c', '-': '#e8edf2' };

function snpTableHtml(d, cols, pageNo, showCons, showDiv) {
  /* d = /api/msa/data 或 /api/tool/msa_data 的返回。返回 {html, info, meta}。 */
  const nChunks = Math.max(1, Math.ceil(d.n_snp / cols));
  const start = pageNo * cols, end = Math.min(d.n_snp, start + cols);
  const pos = d.positions.slice(start, end);
  const tips = pos.map((p, ci) => {
    const counts = {};
    for (const r of d.rows) counts[r[ci]] = (counts[r[ci]] || 0) + 1;
    return `${t('c.snpPos')} ${p} · ` + Object.entries(counts)
      .sort((a, b) => b[1] - a[1]).map(([c, n]) => `${c}:${n}`).join(' ');
  });
  let html = '<table style="border-collapse:collapse;font-family:Consolas,monospace;font-size:11px">';
  if (showDiv) {
    html += '<tr><td style="position:sticky;left:0;background:#fff;font-size:9px;color:#94a3b8;padding-right:6px;text-align:right">变异度</td>';
    for (let ci = start; ci < end; ci++) {
      const h = Math.round((d.diversity[ci] || 0) * 26);
      html += `<td style="width:16px;height:28px;vertical-align:bottom;padding:0">` +
        `<div style="height:${h}px;background:#dc2626;opacity:.75" title="${esc(tips[ci - start])}"></div></td>`;
    }
    html += '</tr>';
  }
  if (showCons) {
    html += '<tr><td style="position:sticky;left:0;background:#fff;font-size:9px;color:#94a3b8;padding-right:6px;text-align:right">' + t('c.consensus') + '</td>';
    for (let ci = start; ci < end; ci++) {
      const c = d.consensus[ci];
      html += `<td style="width:16px;text-align:center;background:#f1f5f9;font-weight:bold" ` +
        `title="${esc(tips[ci - start])}">${esc(c)}</td>`;
    }
    html += '</tr>';
  }
  for (let ri = 0; ri < d.names.length; ri++) {
    const nm = d.names[ri];
    const label = nm.length > 26 ? nm.slice(0, 25) + '…' : nm;
    html += `<tr><td style="position:sticky;left:0;background:#fff;font-size:10px;max-width:180px;` +
      `overflow:hidden;text-overflow:ellipsis;white-space:nowrap;padding-right:6px" title="${esc(nm)}">${esc(label)}</td>`;
    const row = d.rows[ri];
    for (let ci = start; ci < end; ci++) {
      const c = row[ci];
      const bg = SNP_COLORS[c] || '#984ea3';
      const fg = c === '-' ? '#cbd5e1' : '#fff';
      html += `<td style="width:16px;text-align:center;background:${bg};color:${fg}" ` +
        `title="${esc(nm)} · ${esc(tips[ci - start])}">${c === '-' ? '·' : esc(c)}</td>`;
    }
    html += '</tr>';
  }
  html += '</table>';
  return { html,
    info: d.n_snp ? ' ' + t('c.pageOf').replace('{p}', pageNo + 1).replace('{n}', nChunks).replace('{a}', pos[0]).replace('{b}', pos[pos.length - 1]) + ' ' : '',
    meta: `序列 ${d.n_seq} 条${d.seqs_truncated ? '（超出上限仅显示前 80）' : ''}` +
      ` · ${t('c.vsFull').replace('{l}', d.aln_len).replace('{n}', d.n_snp)}` };
}

/* ---- 共享渲染：identity 热图（plotly）与矩阵表 ---- */
function identityHeatmap(divId, names, matrix, title, zminManual) {
  const n = names.length;
  const lowest = Math.min(...matrix.flat());
  const zmin = Math.max(0, zminManual != null && zminManual !== '' ?
                        +zminManual : lowest - 5);
  const box = $(divId);
  Plotly.newPlot(divId, [{
    type: 'heatmap', z: matrix, x: names, y: names,
    zmin, zmax: 100,
    colorscale: [[0, '#b2182b'], [0.5, '#f7f7f7'], [1, '#1a9850']],
    colorbar: { title: { text: 'identity %' }, thickness: 12 },
    hovertemplate: '%{y} vs %{x}<br>%{z:.1f}%<extra></extra>',
  }], {
    title, width: Math.min(Math.max(box.clientWidth, 560), 1100),
    height: Math.min(Math.max(n * 24 + 160, 480), 900),
    margin: { l: 150, b: 120, t: 60, r: 20 },
    xaxis: { tickangle: -45, automargin: true },
    yaxis: { automargin: true, autorange: 'reversed' },
  }, { responsive: false });
  return zmin;
}

function identityMatrixTable(names, matrix) {
  const head = '<tr><th></th>' + names.map(n => `<th title="${esc(n)}">${esc(n.slice(0, 14))}</th>`).join('') + '</tr>';
  return '<tbody>' + head + names.map((n, i) =>
    `<tr><th style="text-align:left;white-space:nowrap" title="${esc(n)}">${esc(n.slice(0, 22))}</th>` +
    matrix[i].map(v => `<td>${v == null ? '' : (+v).toFixed(1)}</td>`).join('') + '</tr>').join('') + '</tbody>';
}

/* ---- 共享渲染：Archaeopteryx.js 树查看器 ----
   面板自带：矩形/环形/无根布局、系统发育图/聚类图切换、支持值与支持圆点、
   梯化排序、中点重根、字号/节点/枝宽滑块、双搜索框、子树上钻、
   Download 菜单（PNG/SVG/PDF/phyloXML/Newick/Nexus）。 */

let _tvViewer = null;

function _treeInternalLabelsAllNumeric(newick) {
  /* FastTree tree.nwk 与 IQ-TREE treefile/contree 把支持值写成内部节点名
     （如 )0.753: 或 )95.2:）；全部为数字时按支持值解析，否则保留为内部名。 */
  const labels = [];
  const re = /\)[ \t]*([^,():;\[\]\s]*)[ \t]*:/g;
  let m;
  while ((m = re.exec(newick)) !== null) {
    const lab = (m[1] || '').trim();
    if (lab) labels.push(lab);
  }
  return labels.length > 0 && labels.every(l => /^[0-9.]+$/.test(l));
}

function renderTreeTo(box, newick, opts) {
  /* opts: {layout}：'rectangular' | 'circular' | 'unrooted'。返回 tips 数。 */
  if (_tvViewer) { try { _tvViewer.destroy(); } catch (e) { /* 重建即可 */ } _tvViewer = null; }
  box.innerHTML = '';
  const asConf = _treeInternalLabelsAllNumeric(newick);
  /* 同一容器重新 launch 即为官方的换树方式；launch 会整体重建面板。 */
  _tvViewer = archaeopteryx.launchArchaeopteryx(box, 'tree.nwk', newick, {
    layout: opts.layout || 'rectangular',
    nhConfidenceValuesAsInternalNames: asConf,
    nhConfidenceValuesInBrackets: !asConf,
    enableDownloads: true,
    enableAccessToDatabases: false,
    pngExportScale: 4,
  });
  /* 面板中支持值复选框默认不勾选；树带支持值时自动勾上。 */
  setTimeout(() => {
    const cb = box.querySelector('#conf_cb');
    if (asConf && cb && !cb.checked) cb.click();
  }, 60);
  let tips = 0;
  try {
    tips = forester.getAllExternalNodes(archaeopteryx.parseTree('tree.nwk', newick)).length;
  } catch (e) { tips = 0; }
  return tips;
}

/* ---- ⑦ MSA 卡片 ---- */
let msData = null, msPageNo = 0, msRun = null;

async function msaRunAlign(btn) {
  const seqs = _v('ms_fa');
  if (!seqs) { alert('请选择或粘贴 FASTA（≥2 条序列）'); return; }
  try {
    const r = await fetch('/api/tool/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ tool: 'structcmp',
                             params: { seqs, max_n: +_v('ms_maxn') || 30 } })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    msRun = d.run;
    _watchTask(d.task, () => msaToolLoadRun(), btn, '⏳ 比对中…');
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

async function msaToolLoadRun() {
  try {
    const r = await fetch(`/api/tool/msa_data?run=${encodeURIComponent(msRun)}`);
    if (!r.ok) { $('msBox').innerHTML = '<p class="hint" style="color:#b91c1c">' + esc((await r.json()).error || '加载失败') + '</p>'; return; }
    msData = await r.json();
    msPageNo = 0;
    $('msDl').innerHTML = `<a class="btn small" href="/tool_runs/${encodeURIComponent(msRun)}/aln.fasta" download>⬇ 下载比对 aln.fasta（运行 ${esc(msRun)}）</a>`;
    msaToolRender();
  } catch (e) { $('msBox').innerHTML = '<p class="hint" style="color:#b91c1c">无法连接: ' + esc(e) + '</p>'; }
}

async function msaToolLoadGroups() {
  const s = _v('msSample'), sel = $('msGroup');
  if (!s) { sel.innerHTML = '<option value="">（先选样品）</option>'; return; }
  const list = await fetchSampleList();
  const item = list.find(x => x.sample === s);
  sel.innerHTML = (item?.groups || []).map(g =>
    `<option value="${esc(g.group)}">${esc(g.group)}${g.aln === 'trim' ? '（清剪后）' : ''}</option>`).join('')
    || '<option value="">（该样品没有比对）</option>';
  if (sel.value) msaToolView();
}

async function msaToolView() {
  const s = _v('msSample'), g = _v('msGroup');
  const box = $('msBox');
  if (!s || !g) return;
  box.innerHTML = '<p class="hint">加载中…</p>';
  try {
    const r = await fetch(`/api/msa/data?sample=${encodeURIComponent(s)}&group=${encodeURIComponent(g)}`);
    if (!r.ok) { box.innerHTML = '<p class="hint" style="color:#b91c1c">' + esc((await r.json()).error || '加载失败') + '</p>'; return; }
    msData = await r.json(); msRun = null; msPageNo = 0;
    $('msDl').innerHTML = '';
    msaToolRender();
  } catch (e) { box.innerHTML = '<p class="hint" style="color:#b91c1c">无法连接: ' + esc(e) + '</p>'; }
}

function msaToolRender() {
  const box = $('msBox');
  if (!msData) { box.innerHTML = `<p class="hint">${t('c.noData2')}</p>`; $('msNav').style.display = 'none'; return; }
  if (!msData.n_snp) {
    box.innerHTML = `<p class="hint">${t('c.noSnp')}</p>`;
    $('msNav').style.display = 'none'; $('msMeta').textContent = ''; return;
  }
  const out = snpTableHtml(msData, +_v('msCols') || 100, msPageNo,
                           $('msCons').checked, $('msDiv').checked);
  box.innerHTML = out.html;
  $('msNav').style.display = '';
  $('msPageInfo').textContent = out.info;
  $('msMeta').textContent = out.meta;
}

function msaToolPage(d) {
  if (!msData) return;
  const nChunks = Math.max(1, Math.ceil(msData.n_snp / (+_v('msCols') || 100)));
  msPageNo = Math.min(nChunks - 1, Math.max(0, msPageNo + d));
  msaToolRender();
}

async function msaToolInit() {
  const sel = $('msSample');
  if (!sel) return;
  const list = await fetchSampleList();
  sel.innerHTML = list.length
    ? list.map(x => `<option value="${esc(x.sample)}">${esc(sampleLabel(x))}</option>`).join('')
    : '<option value="">（暂无比对结果）</option>';
  if (sel.value) msaToolLoadGroups();
}

/* ---- ⑧ 进化树查看器卡片 ---- */
let tvData = null, tvSource = '';

async function treeToolLoadFile() {
  const p = _v('tv_file');
  if (!p) { alert('请选择 Newick 树文件'); return; }
  try {
    const r = await fetch(`/api/tree/file?path=${encodeURIComponent(p)}`);
    if (!r.ok) { $('tvMeta').textContent = esc((await r.json()).error || '加载失败'); return; }
    const d = await r.json();
    tvData = { newick: d.newick, meta: { file: d.file, tool: '' } };
    tvSource = 'file';
    treeToolRender();
  } catch (e) { $('tvMeta').textContent = '无法连接: ' + esc(e); }
}

async function treeToolLoadGroups() {
  const s = _v('tvSample'), sel = $('tvGroup');
  if (!s) { sel.innerHTML = '<option value="">（先选样品）</option>'; return; }
  const list = await fetchSampleList();
  const item = list.find(x => x.sample === s);
  const groups = (item?.groups || []).filter(g => g.trees && g.trees.length);
  sel.innerHTML = groups.map(g => `<option value="${esc(g.group)}">${esc(g.group)}</option>`).join('')
    || '<option value="">（该样品没有树文件）</option>';
  if (sel.value) treeToolLoadFiles();
}

function treeToolLoadFiles() {
  const sel = $('tvFileSel');
  fetchSampleList().then(list => {
    const item = list.find(x => x.sample === _v('tvSample'));
    const g = (item?.groups || []).find(x => x.group === _v('tvGroup'));
    sel.innerHTML = ((g && g.trees) || []).map(tr =>
      `<option value="${esc(tr.file)}">${esc(tr.label)}（${esc(tr.file)}）</option>`).join('')
      || '<option value="">（无树文件）</option>';
    if (sel.value) treeToolView();
  });
}

async function treeToolView() {
  const s = _v('tvSample'), g = _v('tvGroup'), f = _v('tvFileSel');
  if (!s || !g || !f) return;
  $('tvMeta').textContent = '加载中…';
  try {
    const r = await fetch(`/api/tree/data?sample=${encodeURIComponent(s)}&group=${encodeURIComponent(g)}&file=${encodeURIComponent(f)}`);
    if (!r.ok) { tvData = null; $('tvMeta').textContent = esc((await r.json()).error || '加载失败'); return; }
    const d = await r.json();
    tvData = { newick: d.newick, meta: d };
    tvSource = 'sample';
    treeToolRender();
  } catch (e) { $('tvMeta').textContent = '无法连接: ' + esc(e); }
}

function treeToolRender() {
  const box = $('tvBox');
  if (!tvData) { box.innerHTML = `<p class="hint">${t('c.noData2')}</p>`; $('tvMeta').textContent = ''; return; }
  let tips = 0;
  try {
    tips = renderTreeTo(box, tvData.newick, { layout: 'rectangular' });
  } catch (e) {
    box.innerHTML = '<p class="hint" style="color:#b91c1c">树解析失败: ' + esc(e) + '</p>';
    return;
  }
  const m = tvData.meta || {};
  $('tvMeta').textContent = [
    `${tips} 个序列`, m.tool, m.model ? '模型 ' + m.model : '',
    m.logl ? 'logL ' + m.logl : '', m.file,
    '查看器面板：布局切换 / 显示与支持值 / 缩放 / 搜索 / Download 导出',
  ].filter(Boolean).join(' · ');
}

async function treeToolInit() {
  const sel = $('tvSample');
  if (!sel) return;
  const list = await fetchSampleList();
  const withTrees = list.map(x => ({
    sample: x.sample, label: sampleLabel(x),
    groups: (x.groups || []).filter(g => g.trees && g.trees.length)
  })).filter(x => x.groups.length);
  sel.innerHTML = withTrees.length
    ? withTrees.map(x => `<option value="${esc(x.sample)}">${esc(x.label)}</option>`).join('')
    : '<option value="">（暂无树文件，先建树）</option>';
  if (sel.value) treeToolLoadGroups();
}

/* ---- ⑧① 快速建树（FASTA → MAFFT → NJ / FastTree） ---- */
let qtRun = null;

async function quickTreeRun(btn) {
  const seqs = _v('qt_fa');
  if (!seqs) { alert('请选择或粘贴 FASTA（≥2 条序列）'); return; }
  const method = _v('qtMethod') || 'nj';
  try {
    const r = await fetch('/api/tool/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ tool: 'quicktree',
                             params: { seqs, method,
                                       max_n: +_v('qt_maxn') || 100 } })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    qtRun = d.run;
    $('tvMeta').textContent =
      `⏳ MAFFT 比对 + ${method === 'fasttree' ? 'FastTree' : 'NJ'} 建树中…`;
    taskLogOpen.add(d.task);
    startPolling();
    _watchTask(d.task, () => quickTreeLoad(), btn, '⏳ 建树中…');
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

async function quickTreeLoad() {
  try {
    const r = await fetch(`/api/tool/quicktree_data?run=${encodeURIComponent(qtRun)}`);
    if (!r.ok) { $('tvMeta').textContent = esc((await r.json()).error || '加载失败'); return; }
    const d = await r.json();
    tvData = { newick: d.newick, meta: { file: d.file, tool: d.tool } };
    tvSource = 'quicktree';
    treeToolRender();
    const dl = $('qtDl');
    dl.style.display = '';
    dl.href = `/tool_runs/${encodeURIComponent(qtRun)}/${encodeURIComponent(d.file)}`;
    dl.download = d.file;
  } catch (e) { $('tvMeta').textContent = '无法连接: ' + esc(e); }
}

/* ---- ⑨ SDT 卡片（精确引擎：逐对 MAFFT + Get_Similarity，替代 SDT exe） ---- */
let sdData = null, sdRun = null;

async function sdtRunMatrix(btn) {
  const seqs = _v('sd_fa');
  if (!seqs) { alert('请选择或粘贴 FASTA（≥2 条序列）'); return; }
  const mode = _v('sd_mode') || 'nt';          // nt | aa | ntaa
  const ntaa = mode === 'ntaa';
  try {
    const r = await fetch('/api/tool/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(ntaa ? {
        /* NT+AA 同一性表模式（BioAider 口径） */
        tool: 'identity',
        threads: +_v('sd_threads') || null,
        params: { nt_seqs: seqs, aa_seqs: _v('sd_aa'),
                  max_n: +_v('sd_maxn') || 30,
                  aligned: !!$('sd_aligned').checked,
                  palette: _v('sd_palette') || 'sdt' } } : {
        /* NT / AA · SDT 精确矩阵模式（序列类型自动判别兜底） */
        tool: 'sdt',
        threads: +_v('sd_threads') || null,
        params: { seqs, max_n: +_v('sd_maxn') || 30,
                  seqtype: mode,
                  orient: !!$('sd_orient').checked,
                  aligned: !!$('sd_aligned').checked,
                  palette: _v('sd_palette') || 'sdt' } })});
    if (!r.ok) { alert('启动失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    sdRun = d.run;
    sdData = null;
    $('idtResult').style.display = 'none';
    $('sdExact').innerHTML = '';
    $('sdMeta').textContent = ntaa
      ? '⏳ NT+AA 同一性计算中（NT 与 AA 各做逐对比对）…'
      : (mode === 'aa'
        ? '⏳ AA 蛋白同一性逐对比对中…'
        : '⏳ SDT 逐对精确比对中（耗时与序列对数成正比，可用线程越多越快）…');
    taskLogOpen.add(d.task);
    startPolling();
    _watchTask(d.task, () => ntaa ? identityLoad(sdRun) : sdtExactLoad(),
               btn, '⏳ 矩阵分析中…');
  } catch (e) { alert('无法连接平台服务: ' + e); }
}

async function identityLoad(run) {
  try {
    const r = await fetch(`/tool_runs/${encodeURIComponent(run)}/identity.json`);
    if (!r.ok) { alert('加载失败: ' + ((await r.json()).error || '')); return; }
    const d = await r.json();
    const base = `/tool_runs/${encodeURIComponent(run)}`;
    $('idtHeat').innerHTML = d.aa
      ? `<a href="${base}/identity_composite.png" target="_blank"><img src="${base}/identity_composite.png"
           style="max-width:100%;max-height:720px;border:1px solid #dfe5ec;border-radius:8px"></a>`
      : '<p class="hint">有效 AA 序列不足 2 条，未出复合热图（仅 NT 矩阵）</p>';
    const rows = (d.pairs || []).map((p, i) =>
      `<tr><td>${i + 1}</td><td style="font-family:monospace" title="${esc(p.a)}">${esc(p.a)}</td>` +
      `<td style="font-family:monospace" title="${esc(p.b)}">${esc(p.b)}</td>` +
      `<td>${p.nt == null ? '—' : p.nt.toFixed(2)}</td>` +
      `<td>${p.aa == null ? '—' : p.aa.toFixed(2)}</td></tr>`).join('');
    $('idtTable').innerHTML =
      '<tbody><tr><th>#</th><th>序列 A</th><th>序列 B</th><th>NT Identity (%)</th><th>AA Identity (%)</th></tr>' +
      rows + '</tbody>';
    $('idtDl').innerHTML =
      `<a class="btn small" href="${base}/identity_table.csv?dl=1">⬇ 逐对同一性表 CSV</a>` +
      ` <a class="btn small" href="${base}/nt_matrix.csv?dl=1">⬇ NT 矩阵 CSV</a>` +
      (d.aa ? ` <a class="btn small" href="${base}/aa_matrix.csv?dl=1">⬇ AA 矩阵 CSV</a>` +
        ` <a class="btn small" href="${base}/identity_composite.png?dl=1">⬇ 复合热图 PNG</a>` +
        ` <a class="btn small" href="${base}/identity_composite.pdf?dl=1">⬇ 复合热图 PDF（矢量）</a>` : '') +
      ` <a class="btn small" href="${base}/nt_heatmap.png?dl=1">⬇ NT 热图 PNG</a>`;
    $('idtResult').style.display = '';
    $('idtResult').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  } catch (e) { alert('无法连接: ' + e); }
}

async function sdtExactLoad() {
  try {
    const r = await fetch(`/tool_runs/${encodeURIComponent(sdRun)}/sdt_matrix.json`);
    if (!r.ok) { $('sdMeta').textContent = '加载失败: ' + (await r.json()).error; return; }
    sdData = await r.json();
    const base = `/tool_runs/${encodeURIComponent(sdRun)}`;
    $('sdExact').innerHTML = `
      <h4 style="font-size:14px;color:#1a5276;margin:8px 0 8px">SDT 热图（聚类排序 · 三角）</h4>
      <a href="${base}/sdt_heatmap.png" target="_blank"><img src="${base}/sdt_heatmap.png"
         style="max-width:100%;max-height:720px;border:1px solid #dfe5ec;border-radius:8px"></a>
      <h4 style="font-size:14px;color:#1a5276;margin:16px 0 8px">Identity 分布</h4>
      <a href="${base}/sdt_distribution.png" target="_blank"><img src="${base}/sdt_distribution.png"
         style="max-width:520px;width:100%;border:1px solid #dfe5ec;border-radius:8px"></a>`;
    $('sdDl').innerHTML =
      `<a class="btn small" href="${base}/sdt_matrix.csv?dl=1">⬇ 矩阵 CSV</a>` +
      ` <a class="btn small" href="${base}/sdt_heatmap.png?dl=1">⬇ 热图 PNG</a>` +
      ` <a class="btn small" href="${base}/sdt_heatmap.pdf?dl=1">⬇ 热图 PDF（矢量）</a>` +
      ` <a class="btn small" href="${base}/sdt_distribution.png?dl=1">⬇ 分布图 PNG</a>` +
      ` <a class="btn small" href="${base}/sdt_distribution.pdf?dl=1">⬇ 分布图 PDF（矢量）</a>`;
    const names = sdData.names, m = sdData.matrix;
    $('sdTable').innerHTML = identityMatrixTable(names, m);
    const st = (sdData.seqtype || 'nt').toUpperCase();
    $('sdMeta').textContent =
      `${names.length} 条序列（${st}） · ${sdData.pairs} 对独立比对 · SDT v1.3 口径` +
      (sdData.aligned ? ' · 已比对直算' : '') +
      (sdData.palette === 'sdt' ? ' · SDT 经典色阶' : ` · ${sdData.palette} 色阶`);
    $('sdExact').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  } catch (e) { $('sdMeta').textContent = '无法连接: ' + esc(e); }
}

async function sdtToolLoadGroups() {
  const s = _v('sdSample'), sel = $('sdGroup');
  if (!s) { sel.innerHTML = '<option value="">（先选样品）</option>'; return; }
  const list = await fetchSampleList();
  const item = list.find(x => x.sample === s);
  sel.innerHTML = (item?.groups || []).map(g =>
    `<option value="${esc(g.group)}">${esc(g.group)}</option>`).join('')
    || '<option value="">（该样品没有矩阵/比对）</option>';
  if (sel.value) sdtToolView();
}

async function sdtToolView() {
  const s = _v('sdSample'), g = _v('sdGroup');
  if (!s || !g) return;
  $('sdMeta').textContent = '计算/加载矩阵中…';
  try {
    const r = await fetch(`/api/sdt/data?sample=${encodeURIComponent(s)}&group=${encodeURIComponent(g)}`);
    if (!r.ok) { sdData = null; $('sdMeta').textContent = esc((await r.json()).error || '加载失败'); return; }
    sdData = await r.json();
    sdRun = null;
    $('sdDl').innerHTML = '';
    sdtToolRender();
  } catch (e) { $('sdMeta').textContent = '无法连接: ' + esc(e); }
}

function sdtToolRender() {
  const box = $('sdHeat'), meta = $('sdMeta');
  if (!sdData || !sdData.names?.length) {
    box.style.display = 'none'; meta.textContent = ''; return;
  }
  const zmin = identityHeatmap('sdHeat', sdData.names, sdData.matrix,
    `SDT 全长成对 identity 矩阵（${esc(sdRun || (_v('sdSample') + ' · ' + _v('sdGroup')))}）`,
    _v('sdZmin'));
  box.style.display = '';
  $('sdTable').innerHTML = identityMatrixTable(sdData.names, sdData.matrix);
  meta.textContent =
    `${sdData.names.length} 条序列 · 色标 ${zmin.toFixed(0)}–100%` +
    (sdData.source ? ` · 数据来源 ${sdData.source}` : '');
}

async function sdtToolInit() {
  const sel = $('sdSample');
  if (!sel) return;
  const list = await fetchSampleList();
  sel.innerHTML = list.length
    ? list.map(x => `<option value="${esc(x.sample)}">${esc(sampleLabel(x))}</option>`).join('')
    : '<option value="">（暂无矩阵/比对，先跑流程或集合建树）</option>';
  if (sel.value) sdtToolLoadGroups();
}

msaToolInit();
treeToolInit();
sdtToolInit();

// ---------------- 磁盘水位与存储面板 ----------------
// 磁盘剩余每次实时取（毫秒级）；目录占用要遍历 6 万+ 文件，后端异步统计，
// 前端轮询等待结果，避免首页被 47GB 的扫描拖住。
async function checkStorageBanner() {
  try {
    const r = await fetch('/api/storage');
    if (!r.ok) return;
    const d = await r.json();
    const low = (d.disks || []).filter(x => x.warn);
    let b = $('storage-banner');
    if (!low.length) { if (b) b.remove(); return; }
    if (!b) {
      b = document.createElement('div');
      b.id = 'storage-banner';
      b.style.cssText = 'position:sticky;top:0;z-index:98;background:#f9a825;'
        + 'color:#3e2723;padding:10px 16px;font-size:14px;line-height:1.6;';
      const cb = $('conn-banner');
      if (cb && cb.parentNode) cb.parentNode.insertBefore(b, cb.nextSibling);
      else document.body.prepend(b);
    }
    b.innerHTML = '<b>⚠ ' + t('st.lowTitle', '磁盘空间不足') + '</b><br>'
      + low.map(x => esc(x.label) + t('st.disk', '盘') + '（' + esc(x.path) + '）'
        + t('st.onlyLeft', ' 仅剩 ') + esc(x.free_h) + '，'
        + t('st.suggest', '建议保留 ') + d.threshold_gb + 'GB '
        + t('st.above', '以上')).join('<br>')
      + '<br>' + t('st.lowHint',
          '可到「设置 → 存储与磁盘」查看占用明细；或把数据库目录改到空间更大的盘。');
  } catch (e) { /* 磁盘检查失败不应影响正常使用 */ }
}

function _storageRow(x, indent) {
  const pad = indent ? 'padding-left:' + (indent * 16) + 'px;' : '';
  const kindTxt = ({ database: t('st.kDb', '数据库'),
                     output: t('st.kOut', '产物'),
                     input: t('st.kIn', '输入'),
                     build: t('st.kBuild', '打包产物') })[x.kind] || x.kind;
  return '<tr>'
    + '<td style="' + pad + '">' + (indent ? '└ ' : '') + esc(x.name) + '</td>'
    + '<td>' + esc(kindTxt) + '</td>'
    + '<td style="text-align:right">' + esc(x.size_h) + '</td>'
    + '<td style="text-align:right">' + (x.files || 0).toLocaleString() + '</td>'
    + '</tr>';
}

async function loadStoragePanel(times) {
  const box = $('storagePanel');
  if (!box) return;
  let d;
  try {
    const r = await fetch('/api/storage?detail=1');
    if (!r.ok) {
      box.innerHTML = '<p class="hint">' + t('st.loadFail', '加载失败') + '</p>';
      return;
    }
    d = await r.json();
  } catch (e) {
    box.innerHTML = '<p class="hint">' + t('st.loadFail', '加载失败') + '</p>';
    return;
  }

  const disks = (d.disks || []).map(x => {
    const pct = x.used_pct || 0;
    const col = x.warn ? '#c62828' : (pct > 85 ? '#f9a825' : '#2e7d32');
    return '<div style="margin-bottom:8px">'
      + '<div style="display:flex;justify-content:space-between;font-size:13px">'
      + '<span>' + esc(x.label) + t('st.disk', '盘') + ' · '
      + '<span style="color:#666">' + esc(x.path) + '</span></span>'
      + '<span>' + esc(x.free_h) + ' ' + t('st.free', '可用') + ' / '
      + esc(x.total_h) + '</span></div>'
      + '<div style="height:8px;background:#eee;border-radius:4px;overflow:hidden">'
      + '<div style="height:100%;width:' + pct + '%;background:' + col + '"></div>'
      + '</div></div>';
  }).join('');

  const items = d.scan.items || [];
  let rows = '';
  for (const it of items) {
    rows += _storageRow(it, 0);
    for (const c of (it.children || [])) rows += _storageRow(c, 1);
  }
  const table = items.length
    ? '<table class="tbl" style="width:100%;font-size:13px;margin-top:6px">'
      + '<thead><tr><th>' + t('st.colName', '目录') + '</th><th>'
      + t('st.colKind', '类别') + '</th><th style="text-align:right">'
      + t('st.colSize', '占用') + '</th><th style="text-align:right">'
      + t('st.colFiles', '文件数') + '</th></tr></thead>'
      + '<tbody>' + rows + '</tbody></table>'
    : '';

  const status = d.scan.ready
    ? t('st.at', '统计于 ') + esc(d.scan.scanned_at || '')
    : (d.scan.scanning
        ? t('st.scanning', '正在统计（数据库较大，约需几十秒）…')
        : t('st.waiting', '等待统计…'));

  box.innerHTML = disks
    + '<div style="display:flex;align-items:center;gap:8px;margin:10px 0 4px">'
    + '<b style="font-size:13px">' + t('st.breakdown', '目录占用') + '</b>'
    + '<span class="hint">' + status + '</span>'
    + '<button class="btn small" onclick="rescanStorage()" style="margin-left:auto">'
    + t('st.rescan', '重新统计') + '</button></div>'
    + (d.scan.error
        ? '<p class="hint" style="color:#c62828">' + esc(d.scan.error) + '</p>' : '')
    + table
    + '<p class="hint" style="margin-top:8px">' + t('st.tip',
        '提示：dist（打包产物）、tool_runs（工具运行）、logs（日志）可安全清理；'
        + 'databases 下的库请不要手删——重建宿主库代价很高。'
        + '空间紧张时建议把「数据库目录」改到更大的盘。') + '</p>';

  if (!d.scan.ready && times > 0) {
    setTimeout(function () { loadStoragePanel(times - 1); }, 3000);
  }
}

async function rescanStorage() {
  try { await fetch('/api/storage/rescan', { method: 'POST' }); } catch (e) {}
  loadStoragePanel(40);
}

checkStorageBanner();
loadStoragePanel(0);

// ---------------- 顶部导航任务徽章 ----------------
// 「任务中心」链接旁显示运行中任务数；每 15s 静默刷新，
// 仅当页面渲染了徽章（topnav）时启用。
(function navTasksBadgeLoop() {
  const badge = document.getElementById('navTasksBadge');
  if (!badge) return;
  const tick = async () => {
    try {
      const list = await (await fetch('/api/tasks?logs=0')).json();
      const arr = Array.isArray(list) ? list : (list.active || []);
      const n = arr.filter(x => x.status === 'running').length;
      badge.textContent = n > 99 ? '99+' : (n || '');
      badge.style.display = n ? 'inline-block' : 'none';
      badge.title = n ? ('运行中任务 ' + n + ' 个') : '';
    } catch (e) { /* 服务不可达时静默 */ }
    setTimeout(tick, 15000);
  };
  tick();
})();

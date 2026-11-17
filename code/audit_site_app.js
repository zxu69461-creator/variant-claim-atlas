// varclaims reading check, site version: access code -> items and saved progress from /api; every answer is kept in this browser
// and saved to the server (debounced); export / import of a progress file; submit locks the answers.
let ITEMS = [];
let CODE = '', SLOT = 0;
const state = {answers:{}, submitted:null, cur:0, view:'login'};
const $ = id => document.getElementById(id);
const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
const LSK = () => 'vc_audit_site_v1_' + SLOT;

function lsLoad(){ try{ const s = JSON.parse(localStorage.getItem(LSK()) || '{}'); return s; }catch(e){ return {}; } }
function lsSave(){ try{ localStorage.setItem(LSK(), JSON.stringify({answers:state.answers, submitted:state.submitted, cur:state.cur})); }catch(e){} }
function progressDoc(){ return {app:'varclaims-audit', set_version:'v1', slot:SLOT, answers:state.answers, cur:state.cur}; }

function complete(it){ const a = state.answers[it.id]; if (!a) return false;
  if (it.kind === 'pair') return !!a.v; return it.rsids.every(r => a.v && a.v[r]); }
function nDone(){ return ITEMS.filter(complete).length; }

function setSave(txt, off){ const s = $('saveS'); s.textContent = txt; s.classList.toggle('off', !!off); }
async function api(path, method, body){
  const r = await fetch('/api/' + path, {method, headers:{'content-type':'application/json', 'x-audit-code':CODE}, body: body == null ? undefined : JSON.stringify(body)});
  const j = await r.json().catch(() => ({}));
  if (!r.ok) { const e = new Error(j.error || ('http_' + r.status)); e.status = r.status; e.body = j; throw e; }
  return j;
}
let timer = null, inflight = false, again = false;
function push(){
  lsSave();
  if (state.submitted) return;
  setSave('正在保存…');
  clearTimeout(timer); timer = setTimeout(flush, 700);
}
async function flush(){
  if (inflight){ again = true; return; }
  inflight = true;
  try{ const j = await api('progress', 'PUT', progressDoc()); setSave('已保存到服务器 ' + new Date(j.saved_at).toLocaleTimeString('zh-CN')); }
  catch(e){ if (e.status === 409){ state.submitted = e.body.submitted_at; lsSave(); renderProg(); setSave('已提交，答案已锁定'); }
            else setSave('服务器保存失败，已保存在本机，稍后自动重试', true); setTimeout(() => { if (!state.submitted) flush(); }, 15000); }
  inflight = false;
  if (again){ again = false; flush(); }
}

function mark(text, rsids, outcome){
  let h = esc(text);
  if (outcome){ const o = esc(outcome).replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
    if (o.length >= 3) h = h.replace(new RegExp('(' + o + ')', 'gi'), '<mark class="oc">$1</mark>'); }
  for (const r of rsids){ const n = r.slice(2);
    h = h.replace(new RegExp('\\b(rs\\s?' + n + ')\\b', 'gi'), '<mark class="rs">$1</mark>'); }
  return h;
}

function show(v){ for (const k of ['Login','Intro','Item','Grid','End']) $('v' + k).hidden = (k.toLowerCase() !== v);
  $('hdr').hidden = (v === 'login'); state.view = v; window.scrollTo(0, 0); }

function renderProg(){ const d = nDone(), n = ITEMS.length || 1;
  $('progT').textContent = `已完成 ${d} / ${ITEMS.length}` + (state.submitted ? ' · 已提交' : '');
  $('progI').style.width = (100 * d / n) + '%';
  $('introProg').textContent = d ? `已完成 ${d} 题，上次停在第 ${state.cur + 1} 题` : '';
}

const PAIR_OPTS = [['yes','有关联','摘要报告存在关联'],['no','报告无关联','摘要说没有关联或不显著'],
                   ['none','无结论','提到了，但没下结论'],['notstudied','不相关','摘要没研究这对关系']];

function renderItem(){
  const it = ITEMS[state.cur], a = state.answers[it.id] || {}, locked = !!state.submitted;
  const part = it.kind === 'pair' ? '第一部分 · 判断一对关系' : '第二部分 · 查漏';
  let q, body;
  if (it.kind === 'pair'){
    q = `<p class="q">摘要对 <span class="rsid">${esc(it.rsid)}</span> 与「${esc(it.outcome)}」的结论是？</p>`;
    body = `<div class="choices" role="group" aria-label="选项">` + PAIR_OPTS.map(([v,l,s]) =>
      `<button type="button" data-v="${v}" class="${a.v === v ? 'on' : ''}" ${locked ? 'disabled' : ''}>${l}<small>${s}</small></button>`).join('') + `</div>`;
  } else {
    q = `<p class="q">下列变异中，摘要报告与某个疾病或性状存在有统计学支持的关联的是？</p>`;
    body = `<div class="box stack">` + it.rsids.map(r => { const x = (a.v || {})[r];
      return `<div class="rsrow"><span class="rsid">${esc(r)}</span><span class="seg">
        <button type="button" class="${x === 'yes' ? 'on' : ''}" data-r="${r}" data-x="yes" ${locked ? 'disabled' : ''}>有</button>
        <button type="button" class="n ${x === 'no' ? 'on' : ''}" data-r="${r}" data-x="no" ${locked ? 'disabled' : ''}>无</button></span></div>`; }).join('') + `</div>`;
  }
  const rs = it.kind === 'pair' ? [it.rsid] : it.rsids;
  $('vItem').innerHTML = `
    <div class="meta">${part} · ${it.id} · 第 ${state.cur + 1} / ${ITEMS.length} 题</div>
    ${q}
    <article class="box paper"><h3>${mark(it.title, rs, it.outcome)}</h3><div class="abs">${mark(it.abstract || '(无摘要)', rs, it.outcome)}</div></article>
    ${body}
    <textarea id="note" placeholder="备注（可不填）" ${locked ? 'disabled' : ''}>${esc(a.note || '')}</textarea>
    <div class="nav"><button type="button" id="bPrev" ${state.cur === 0 ? 'disabled' : ''}>上一题</button>
      <button type="button" id="bHome">说明</button>
      <button type="button" class="pri" id="bNext">${state.cur === ITEMS.length - 1 ? '完成' : '下一题'}</button></div>`;
  $('vItem').querySelectorAll('.choices button').forEach(b => b.onclick = () => {
    setAns(it.id, {v: b.dataset.v}); renderItem(); setTimeout(() => go(1), 260); });
  $('vItem').querySelectorAll('.seg button').forEach(b => b.onclick = () => {
    const v = Object.assign({}, (state.answers[it.id] || {}).v || {}); v[b.dataset.r] = b.dataset.x;
    setAns(it.id, {v}); renderItem(); });
  $('note').onchange = e => setAns(it.id, {note: e.target.value});
  $('bPrev').onclick = () => go(-1); $('bNext').onclick = () => go(1); $('bHome').onclick = () => { renderProg(); show('intro'); };
}

function setAns(id, patch){ if (state.submitted) return;
  state.answers[id] = Object.assign({}, state.answers[id] || {}, patch, {at: new Date().toISOString()});
  renderProg(); push(); }

function go(d){ const n = state.cur + d;
  if (n >= ITEMS.length){ renderEnd(); show('end'); return; }
  state.cur = Math.max(0, n); lsSave(); renderItem(); show('item'); }

function renderGrid(){ $('gridBox').innerHTML = ITEMS.map((it, i) =>
  `<button type="button" data-i="${i}" class="${complete(it) ? 'done' : ''} ${i === state.cur ? 'cur' : ''}">${it.id}</button>`).join('');
  $('gridBox').querySelectorAll('button').forEach(b => b.onclick = () => { state.cur = +b.dataset.i; renderItem(); show('item'); }); }

function renderEnd(){ const d = nDone(), n = ITEMS.length, open = ITEMS.filter(i => !complete(i));
  let h = `<h1>${state.submitted ? '已提交，谢谢' : '全部题目'}</h1>`;
  if (state.submitted) h += `<p class="lead">提交时间 ${esc(new Date(state.submitted).toLocaleString('zh-CN'))}。答案已锁定，可以回看。</p>`;
  else if (open.length) h += `<p>还有 ${open.length} 题未完成：${open.slice(0, 12).map(i => i.id).join('、')}${open.length > 12 ? ' 等' : ''}。</p>
      <div class="row"><button class="pri" type="button" id="bOpen">去第一道未完成的题</button></div>`;
  else h += `<p>${d} / ${n} 题都已完成。提交后答案锁定，结果会进入分析。</p>
      <div class="row" id="subRow"><button class="pri" type="button" id="bSubmit">提交</button></div>`;
  h += `<div class="row"><button type="button" id="bEndHome">回到说明</button></div>`;
  $('vEnd').innerHTML = h;
  if ($('bOpen')) $('bOpen').onclick = firstOpen;
  if ($('bSubmit')) $('bSubmit').onclick = () => { $('subRow').innerHTML =
      `<span>确认提交？提交后不能再修改。</span><button class="pri" type="button" id="bYes">确认提交</button><button type="button" id="bNo">再看看</button>`;
    $('bNo').onclick = renderEnd; $('bYes').onclick = submit; };
  $('bEndHome').onclick = () => { renderProg(); show('intro'); };
}

async function submit(){
  clearTimeout(timer);
  try{ const j = await api('submit', 'POST', progressDoc()); state.submitted = j.submitted_at; }
  catch(e){ if (e.status === 409) state.submitted = e.body.submitted_at;
    else { $('vEnd').insertAdjacentHTML('beforeend', '<p class="notice">提交没有成功，请检查网络后再点一次，或在说明页导出进度文件发给我。</p>'); return; } }
  lsSave(); renderProg(); setSave('已提交，答案已锁定'); renderEnd(); }

function firstOpen(){ const i = ITEMS.findIndex(x => !complete(x)); state.cur = i < 0 ? 0 : i; renderItem(); show('item'); }

$('bStart').onclick = () => { if (nDone() === ITEMS.length){ renderEnd(); show('end'); } else { renderItem(); show('item'); } };
$('bRules').onclick = () => { renderProg(); show('intro'); $('rulesBody').scrollIntoView(); };
$('bGrid').onclick = () => { renderGrid(); show('grid'); };
$('bGridBack').onclick = () => { if (state.view === 'grid'){ renderItem(); show('item'); } };
$('bFirstOpen').onclick = firstOpen;

$('bExport').onclick = async () => {
  const txt = JSON.stringify({form: 'varclaims_reader_check_v1', exported: new Date().toISOString(), answers: state.answers, submitted: state.submitted});
  try{ const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([txt], {type:'application/json'}));
    a.download = '变异判读核对_进度.json'; document.body.appendChild(a); a.click(); a.remove(); $('ioMsg').textContent = '已导出。'; return; }catch(e){}
  const t = $('exportTxt'); t.hidden = false; t.value = txt; t.select();
  try{ await navigator.clipboard.writeText(txt); $('ioMsg').textContent = '无法直接保存文件，进度已复制到剪贴板（也显示在上面）。'; }
  catch(e){ $('ioMsg').textContent = '无法直接保存文件，请复制上面的文字保存。'; }
};
$('fImport').onchange = async e => { const f = e.target.files[0]; if (!f) return;
  try{ const j = JSON.parse(await f.text()); if (j.form !== 'varclaims_reader_check_v1') throw 0;
    let k = 0; for (const [id, a] of Object.entries(j.answers || {})){
      if (!ITEMS.some(i => i.id === id)) continue;
      const cur = state.answers[id]; if (!cur || (a.at || '') > (cur.at || '')){ state.answers[id] = a; k++; } }
    renderProg(); push(); $('ioMsg').textContent = `已导入 ${k} 题。`; }
  catch(err){ $('ioMsg').textContent = '这个文件不是本问卷导出的进度文件。'; } };

async function login(code, quiet){
  CODE = String(code || '').replace(/[\s-]+/g, '').toUpperCase();
  $('loginMsg').textContent = '正在登录…';
  let j;
  try{ j = await api('login', 'POST', {code: CODE, name: $('nameIn').value}); }
  catch(e){ $('loginMsg').textContent = e.status === 401 ? '访问码不对，请再核对一下。' : '连不上服务器，请稍后再试。'; return; }
  ITEMS = j.items; SLOT = j.slot; $('loginMsg').textContent = '';
  try{ localStorage.setItem('vc_audit_site_code', CODE); }catch(e){}
  const local = lsLoad(), remote = (j.progress && j.progress.answers) || {};
  state.answers = Object.assign({}, remote);
  let newer = false;
  for (const [id, a] of Object.entries(local.answers || {})){ const r = remote[id];
    if (!r || (a.at || '') > (r.at || '')){ state.answers[id] = a; newer = true; } }
  state.cur = Number.isInteger(local.cur) ? local.cur : ((j.progress && j.progress.cur) || 0);
  state.submitted = j.submitted_at || null;
  $('who').textContent = (j.name ? j.name + ' · ' : '') + '读者 ' + j.label;
  lsSave(); renderProg(); show('intro');
  setSave(state.submitted ? '已提交，答案已锁定' : (j.updated_at ? '已连接服务器，自动保存' : '已连接服务器'));
  if (newer && !state.submitted) push();
}
$('bLogin').onclick = () => login($('codeIn').value);
$('codeIn').onkeydown = e => { if (e.key === 'Enter') login($('codeIn').value); };
$('bLogout').onclick = () => { try{ localStorage.removeItem('vc_audit_site_code'); }catch(e){} location.reload(); };
(() => { let c = ''; try{ c = localStorage.getItem('vc_audit_site_code') || ''; }catch(e){}
  show('login'); if (c){ $('codeIn').value = c; login(c, true); } })();

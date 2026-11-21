/* Variant–Trait Claim Atlas: a static web tool over data/*.json (built by scripts/34_lookup_data.py) and atlas/Data_S3 (ClinVar table).
   Views (URL hash): #lookup (default; #rs123, #gene=SYMBOL, #trait=Label, #q=text), #check (#check=rs123|trait), #batch, #clinvar, #about.
   Everything runs in the browser; no data leave it. */
(function () {
  'use strict';
  const $ = (s) => document.querySelector(s);
  const out = $('#out');
  const esc = (s) => String(s == null ? '' : s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const cap = (s) => (s ? s[0].toUpperCase() + s.slice(1) : s);
  const pct = (x) => (x == null || isNaN(x) ? '—' : (100 * x).toFixed(x < 0.1 ? 1 : 0) + '%');
  const num = (x) => (x == null ? '—' : Number(x).toLocaleString('en-US'));
  const DV = '?v=1.1.0';                                   // data version: bump when data/ changes
  let V, T, META, byRs, byGene, CV = null, shards = {};

  const obj = (t) => t.rows.map((r) => Object.fromEntries(t.cols.map((c, i) => [c, r[i]])));
  const cls = (x) => (x.reported < 5 || x.testable === 0 ? '' : x.supported >= 0.5 ? 'mostly supported' : x.supported >= 0.05 ? 'partly supported' : 'rarely supported');

  async function load() {
    const [v, t, m] = await Promise.all(['variants', 'traits', 'meta'].map((n) => fetch('data/' + n + '.json' + DV).then((r) => r.json())));
    V = obj(v); T = obj(t); META = m;
    byRs = new Map(V.map((x) => [x.rsid, x]));
    byGene = new Map();
    V.forEach((x) => { if (x.gene) { const g = x.gene.toUpperCase(); (byGene.get(g) || byGene.set(g, []).get(g)).push(x); } });
    $('#traitlist').innerHTML = T.slice(0, 3000).map((t) => `<option value="${esc(cap(t.trait))}">`).join('');
    $('#intro').textContent = `The 201,720 disease or trait claims of the atlas include ${num(META.claims)} that report an association or its absence, from ${num(META.papers)} PubMed papers on ${num(META.variants)} variants; each is set against genome-wide significant associations in the GWAS Catalog (release of ${META.gwas_catalog}).`;
    $('#rel').textContent = `Sources: PubMed ${META.pubmed}; GWAS Catalog ${META.gwas_catalog}; Experimental Factor Ontology ${META.efo}; ClinVar ${META.clinvar}. Built ${META.built}.`;
    route();
  }
  async function claims(rs) {
    const s = rs % META.shards;
    if (!shards[s]) shards[s] = obj(await fetch('data/claims/' + String(s).padStart(2, '0') + '.json' + DV).then((r) => r.json()));
    return shards[s].filter((c) => c.rsid === rs);
  }
  function download(name, rows) {
    const q = (v) => { const s = v == null ? '' : String(v); return /[",\n\t]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; };
    const blob = new Blob([rows.map((r) => r.map(q).join(',')).join('\n') + '\n'], { type: 'text/csv' });
    const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500);
  }

  function meter(sup, nul) {
    if (sup == null) return '';
    return `<div class="meter" role="img" aria-label="supported ${pct(sup)}, matched null ${pct(nul)}"><i style="width:${Math.min(100, 100 * sup)}%"></i>${nul != null ? `<u style="left:${Math.min(100, 100 * nul)}%"></u>` : ''}</div>
      <p class="muted" style="font-size:12.5px">Bar: share of testable reported associations with genome-wide support. Tick: matched null.</p>`;
  }
  function kpis(x, extra) {
    return `<div class="kpis">
      <div class="kpi"><b>${num(x.papers)}</b><span>papers</span></div>
      <div class="kpi"><b>${num(x.reported)}</b><span>reported associations</span></div>
      <div class="kpi"><b>${num(x.absent)}</b><span>reports of no association</span></div>
      <div class="kpi"><b>${num(x.testable)}</b><span>testable reported associations</span></div>
      <div class="kpi s"><b>${pct(x.supported)}</b><span>with genome-wide support</span></div>
      <div class="kpi n"><b>${pct(x.null)}</b><span>matched null</span></div>${extra || ''}</div>`;
  }
  const supCell = (c) => (c.status !== 'reported' ? '<span class="nt">—</span>' : !c.testable ? '<span class="nt">not testable</span>' : c.supported ? `<span class="yes">yes${c.earlier ? ', earlier' : ''}</span>` : '<span class="no">no</span>');
  const claimRow = (c) => `<tr><td>${c.year}</td><td><a href="https://pubmed.ncbi.nlm.nih.gov/${c.pmid}/" target="_blank" rel="noopener">${c.pmid}</a></td><td class="t">${esc(cap(c.trait))}</td><td>${c.status === 'reported' ? 'associated' : 'no association'}</td><td class="num">${num(c.n)}</td><td>${esc(c.design)}</td><td>${esc(c.population)}</td><td>${supCell(c)}</td></tr>`;
  const claimHead = '<tr><th>Year</th><th>PMID</th><th>Trait</th><th>Claim</th><th class="num">n</th><th>Design</th><th>Population</th><th>Genome-wide support</th></tr>';
  const claimCsv = (cl) => [['rsid', 'pmid', 'year', 'trait', 'claim', 'stated_n', 'design', 'population', 'testable', 'genome_wide_support', 'support_earlier_than_claim', 'first_genome_wide_year']]
    .concat(cl.map((c) => ['rs' + c.rsid, c.pmid, c.year, c.trait, c.status === 'reported' ? 'associated' : 'no association', c.n, c.design, c.population, c.testable, c.status === 'reported' ? c.supported : '', c.status === 'reported' ? c.earlier : '', c.first_gw_year]));

  // ---------------------------------------------------------------- look up
  async function showVariant(rs) {
    const x = byRs.get(rs);
    if (!x) { out.innerHTML = `<p>rs${rs} is not in the atlas. The atlas holds variants named by rsID in a PubMed title or abstract with a disease or trait claim.</p>`; return; }
    const extra = `<div class="kpi"><b>${num(x.median_n)}</b><span>median stated sample size</span></div>
      <div class="kpi"><b>${x.first_gw_year || '—'}</b><span>first genome-wide report for a claimed trait</span></div>
      ${x.clinvar ? `<div class="kpi"><b style="font-size:16px">${esc(x.clinvar)}</b><span>ClinVar assertion</span></div>` : ''}`;
    const sizes = META.size_classes.map((lab, i) => {
      const [k, s] = x.size_classes[i]; const sh = k ? s / k : null;
      return `<div class="sz"><span>${lab} participants: <b>${k ? pct(sh) : '—'}</b> <span class="muted">(${s} of ${k})</span></span>
        <div class="meter" style="height:8px"><i style="width:${sh == null ? 0 : 100 * sh}%"></i></div></div>`;
    }).join('');
    const tt = x.top_traits.map(([t, k, tt, s]) => `<tr><td class="t"><a href="#check=rs${rs}|${encodeURIComponent(t)}">${esc(cap(t))}</a></td><td class="num">${k}</td><td class="num">${tt}</td><td class="num">${tt ? pct(s / tt) : '—'}</td></tr>`).join('');
    const c = cls(x);
    out.innerHTML = `<h2>rs${rs}${x.gene ? ` <span class="muted" style="font-weight:400">· <i>${esc(x.gene)}</i></span>` : ''}${c ? ` <span class="muted" style="font-weight:400;font-size:14px">· ${c}</span>` : ''}</h2>
      ${kpis(x, extra)}${meter(x.supported, x.null)}
      <h3>Support by stated sample size</h3><div class="sizes">${sizes}</div>
      <h3>Most-claimed traits <span style="text-transform:none;letter-spacing:0;font-weight:400">(click to check the claim)</span></h3><div class="tw"><table><thead><tr><th>Trait</th><th class="num">Reported</th><th class="num">Testable</th><th class="num">Supported</th></tr></thead><tbody>${tt}</tbody></table></div>
      <div class="row"><h3>Every claim</h3><div style="display:flex;gap:8px;flex-wrap:wrap"><input id="cf" placeholder="Filter by trait, year, design" style="width:220px"><button class="sec" id="dl">Download CSV</button></div></div>
      <div class="tw" id="cl"><p class="muted">Loading claims…</p></div>`;
    const cl = await claims(rs);
    const draw = (f) => {
      const s = f.toLowerCase(); const L = s ? cl.filter((c) => (c.trait + ' ' + c.year + ' ' + c.design + ' ' + c.population).toLowerCase().includes(s)) : cl;
      $('#cl').innerHTML = `<table><thead>${claimHead}</thead><tbody>${L.map(claimRow).join('')}</tbody></table><p class="muted" style="font-size:12.5px">${L.length} of ${cl.length} claims</p>`;
    };
    draw('');
    $('#cf').addEventListener('input', (e) => draw(e.target.value));
    $('#dl').addEventListener('click', () => download(`rs${rs}_claims.csv`, claimCsv(cl)));
  }
  function showTrait(name) {
    const t = T.find((r) => r.trait === name) || T.find((r) => r.trait.toLowerCase() === name.toLowerCase());
    if (!t) return showSearch(name);
    const rows = t.top_variants.map(([rs, g, k, tt, s]) => `<tr><td><a href="#rs${rs}">rs${rs}</a></td><td><i>${esc(g)}</i></td><td class="num">${k}</td><td class="num">${tt}</td><td class="num">${tt ? pct(s / tt) : '—'}</td><td><a href="#check=rs${rs}|${encodeURIComponent(t.trait)}">check</a></td></tr>`).join('');
    out.innerHTML = `<h2>${esc(cap(t.trait))} <span class="muted" style="font-weight:400;font-size:14px">· ${esc(t.area)}</span></h2>${kpis(t)}${meter(t.supported, t.null)}
      <h3>Most-claimed variants</h3><div class="tw"><table><thead><tr><th>Variant</th><th>Gene</th><th class="num">Reported</th><th class="num">Testable</th><th class="num">Supported</th><th></th></tr></thead><tbody>${rows}</tbody></table></div>`;
  }
  function showGene(g) {
    const L = (byGene.get(g.toUpperCase()) || []).slice().sort((a, b) => b.reported - a.reported);
    if (!L.length) return showSearch(g);
    out.innerHTML = `<div class="row"><h2><i>${esc(g.toUpperCase())}</i> <span class="muted" style="font-weight:400;font-size:14px">· ${L.length} variant${L.length > 1 ? 's' : ''} (gene from dbSNP or the GWAS Catalog)</span></h2><button class="sec" id="dl">Download CSV</button></div>
      <div class="tw"><table><thead><tr><th>Variant</th><th class="num">Papers</th><th class="num">Reported</th><th class="num">Testable</th><th class="num">Supported</th><th class="num">Matched null</th><th>Class</th><th>ClinVar</th></tr></thead><tbody>${
      L.slice(0, 300).map((x) => `<tr><td><a href="#rs${x.rsid}">rs${x.rsid}</a></td><td class="num">${num(x.papers)}</td><td class="num">${num(x.reported)}</td><td class="num">${num(x.testable)}</td><td class="num">${pct(x.supported)}</td><td class="num">${pct(x.null)}</td><td>${cls(x)}</td><td>${esc(x.clinvar)}</td></tr>`).join('')}</tbody></table></div>`;
    $('#dl').addEventListener('click', () => download(`${g.toUpperCase()}_variants.csv`, batchCsv(L.map((x) => ({ q: 'rs' + x.rsid, x })))));
  }
  function showSearch(q) {
    const s = q.toLowerCase();
    const hits = T.filter((t) => t.trait.toLowerCase().includes(s)).slice(0, 40);
    out.innerHTML = hits.length
      ? `<h2>Traits matching “${esc(q)}”</h2><div class="list">${hits.map((t) => `<a href="#trait=${encodeURIComponent(t.trait)}">${esc(cap(t.trait))} <span class="muted">· ${num(t.reported)} reported associations</span></a>`).join('')}</div>`
      : `<p>Nothing matches “${esc(q)}”. Try an rsID (rs1801133), a gene symbol (FTO) or part of a disease name (asthma).</p>`;
  }
  function home() {
    out.innerHTML = `<h2>What this tool does</h2>
      <p>Before you cite, curate or follow up a published claim that a variant is associated with a disease or trait, see how often claims about that variant and trait are backed by genome-wide evidence, in studies of which size.</p>
      <div class="grid3">
        <a class="tool" href="#rs1800795"><b>Look up</b><span>A variant, a gene or a trait: support of its claims against the matched null, by sample size, with every claim linked to PubMed.</span></a>
        <a class="tool" href="#check=rs1800795|type 2 diabetes"><b>Check a claim</b><span>One variant and one trait: does genome-wide evidence support the association, and since when?</span></a>
        <a class="tool" href="#batch"><b>Batch check</b><span>Paste or upload a list of variants (or variant–trait pairs) and download the calibration as CSV.</span></a>
        <a class="tool" href="#clinvar"><b>ClinVar risk variants</b><span>The ${num(239)} variants that ClinVar lists as risk factors, associations or protective variants, with their genome-wide support.</span></a>
      </div>`;
  }

  // ---------------------------------------------------------------- check a claim
  async function showCheck(arg) {
    const [a, b] = (arg || '').split('|');
    out.innerHTML = `<h2>Check a claim</h2>
      <p class="muted">Enter a variant and the disease or trait it is claimed to affect.</p>
      <form id="cf2"><div class="fld" style="flex:0 1 180px"><label for="crs">rsID</label><input id="crs" placeholder="rs1800795" value="${esc(a || '')}"></div>
      <div class="fld"><label for="ctr">Disease or trait (part of the name is enough)</label><input id="ctr" list="traitlist" placeholder="type 2 diabetes" value="${esc(b ? decodeURIComponent(b) : '')}"></div>
      <button type="submit">Check</button></form><div id="cres"></div>`;
    $('#cf2').addEventListener('submit', (e) => { e.preventDefault(); const r = $('#crs').value.trim().toLowerCase().replace(/^(\d)/, 'rs$1'); const t = $('#ctr').value.trim(); if (r) location.hash = 'check=' + r + '|' + encodeURIComponent(t); });
    if (!a) return;
    const m = a.match(/^rs(\d+)$/i); const res = $('#cres');
    if (!m) { res.innerHTML = '<p>Enter an rsID such as rs1800795.</p>'; return; }
    const rs = +m[1]; const x = byRs.get(rs);
    if (!x) { res.innerHTML = `<div class="verdict nt"><b>rs${rs} has no disease or trait claim in the atlas.</b><span>The atlas holds variants named by rsID in a PubMed title or abstract.</span></div>`; return; }
    const cl = await claims(rs); const tq = (b ? decodeURIComponent(b) : '').toLowerCase().trim();
    const M = tq ? cl.filter((c) => c.trait.toLowerCase().includes(tq)) : cl;
    if (!M.length) {
      res.innerHTML = `<div class="verdict nt"><b>No claim on rs${rs} names a trait matching “${esc(tq)}”.</b><span>Traits claimed for this variant:</span></div>
        <div class="list">${x.top_traits.map(([t]) => `<a href="#check=rs${rs}|${encodeURIComponent(t)}">${esc(cap(t))}</a>`).join('')}</div>`;
      return;
    }
    const R = M.filter((c) => c.status === 'reported'), Tt = R.filter((c) => c.testable), S = Tt.filter((c) => c.supported), E = S.filter((c) => c.earlier);
    const A = M.filter((c) => c.status !== 'reported'); const yrs = S.map((c) => c.first_gw_year).filter(Boolean);
    const traits = [...new Set(M.map((c) => cap(c.trait)))];
    const ns = R.map((c) => c.n).filter((n) => n != null).sort((p, q) => p - q); const med = ns.length ? ns[Math.floor(ns.length / 2)] : null;
    let kind, head, body;
    if (!Tt.length) { kind = 'nt'; head = 'Not testable against genome-wide evidence'; body = 'The GWAS Catalog has no association for this trait or a narrower term, so the claims cannot be calibrated.'; }
    else if (S.length / Tt.length >= 0.5) { kind = ''; head = `Supported genome-wide: ${S.length} of ${Tt.length} testable reported associations`; body = `A genome-wide significant association of this variant or a proxy with the trait was first reported in ${Math.min(...yrs)}; ${E.length} of the supported claims came after that report.`; }
    else if (S.length) { kind = ''; head = `Partly supported: ${S.length} of ${Tt.length} testable reported associations`; body = `Genome-wide support exists for some of the matched traits (first reported in ${Math.min(...yrs)}).`; }
    else { kind = 'no'; head = `No genome-wide support: 0 of ${Tt.length} testable reported associations`; body = `For comparison, ${pct(x.null)} of frequency-matched random variants would be supported for this variant's claimed traits. Absence of support is not evidence of no effect.`; }
    res.innerHTML = `<div class="verdict ${kind}"><b>${esc(head)}</b><span>${esc(body)}</span>
        <span class="muted" style="font-size:13px">rs${rs}${x.gene ? ' (' + esc(x.gene) + ')' : ''}; matched traits: ${esc(traits.slice(0, 6).join('; '))}${traits.length > 6 ? '; …' : ''}. ${R.length} reported associations (median stated sample size ${num(med)}), ${A.length} reports of no association.${x.clinvar ? ' ClinVar: ' + esc(x.clinvar) + '.' : ''}</span></div>
      <div class="row"><h3>Matching claims</h3><button class="sec" id="dl">Download CSV</button></div>
      <div class="tw"><table><thead>${claimHead}</thead><tbody>${M.map(claimRow).join('')}</tbody></table></div>
      <p class="muted" style="font-size:12.5px"><a href="#rs${rs}">All claims on rs${rs}</a></p>`;
    $('#dl').addEventListener('click', () => download(`rs${rs}_${tq.replace(/\W+/g, '_') || 'all'}_claims.csv`, claimCsv(M)));
  }

  // ---------------------------------------------------------------- batch
  const batchCsv = (L) => [['query', 'rsid', 'gene', 'trait_filter', 'papers', 'reported_associations', 'reports_of_no_association', 'testable', 'supported_share', 'matched_null', 'class', 'median_stated_n', 'first_genome_wide_year', 'clinvar', 'in_atlas']]
    .concat(L.map(({ q, x, t, r }) => r ? [q, 'rs' + r.rs, x && x.gene, t, '', r.R, r.A, r.T, r.T ? (r.S / r.T).toFixed(3) : '', x && x.null, '', '', r.y || '', x && x.clinvar, 'yes']
      : x ? [q, 'rs' + x.rsid, x.gene, '', x.papers, x.reported, x.absent, x.testable, x.supported, x.null, cls(x), x.median_n, x.first_gw_year, x.clinvar, 'yes'] : [q, '', '', t || '', '', '', '', '', '', '', '', '', '', '', 'no']));
  function showBatch() {
    out.innerHTML = `<h2>Batch check</h2>
      <p class="muted">Paste rsIDs (any separator), or one variant–trait pair per line (<code>rs1800795, type 2 diabetes</code>), or upload a text or CSV file. Up to 10,000 lines; the file is read in your browser and not sent anywhere.</p>
      <textarea id="bt" placeholder="rs9939609&#10;rs1800795&#10;rs738409, fatty liver disease"></textarea>
      <div class="row"><div style="display:flex;gap:8px;flex-wrap:wrap"><button id="brun">Check</button><button class="sec" id="bex">Load example</button><label class="sec" style="display:inline-flex;align-items:center;gap:6px;cursor:pointer"><input type="file" id="bfile" accept=".txt,.csv,.tsv" style="width:auto;border:0;padding:0"> </label></div><button class="sec" id="bdl" disabled>Download CSV</button></div>
      <div id="bres"></div>`;
    let last = [];
    $('#bex').addEventListener('click', () => { $('#bt').value = 'rs9939609\nrs1800795\nrs738409, fatty liver disease\nrs1801133, coronary artery\nrs7903146\nrs2910164\nrs671, alcohol'; });
    $('#bfile').addEventListener('change', async (e) => { const f = e.target.files[0]; if (f) $('#bt').value = await f.text(); });
    $('#bdl').addEventListener('click', () => download('variant_claim_atlas_batch.csv', batchCsv(last)));
    $('#brun').addEventListener('click', async () => {
      const lines = $('#bt').value.split(/\r?\n/).map((l) => l.trim()).filter(Boolean).slice(0, 10000); const L = [];
      for (const l of lines) {
        const pair = l.match(/^(rs\d+)\s*[,;\t]\s*(\S.*)$/i);
        if (pair) L.push({ q: l, rs: +pair[1].slice(2), t: pair[2].replace(/^["']|["']$/g, '').trim() });
        else for (const m of l.matchAll(/rs(\d+)/gi)) L.push({ q: 'rs' + m[1], rs: +m[1] });
      }
      const res = $('#bres'); res.innerHTML = `<p class="muted">Checking ${L.length} entries…</p>`;
      last = [];
      for (const e of L) {
        const x = byRs.get(e.rs);
        if (e.t && x) {
          const cl = (await claims(e.rs)).filter((c) => c.trait.toLowerCase().includes(e.t.toLowerCase()));
          const R = cl.filter((c) => c.status === 'reported'), Tt = R.filter((c) => c.testable), S = Tt.filter((c) => c.supported);
          const yrs = S.map((c) => c.first_gw_year).filter(Boolean);
          last.push({ q: e.q, x, t: e.t, r: { rs: e.rs, R: R.length, A: cl.length - R.length, T: Tt.length, S: S.length, y: yrs.length ? Math.min(...yrs) : null } });
        } else last.push({ q: e.q, x, t: e.t });
      }
      const found = last.filter((r) => r.x).length;
      res.innerHTML = `<p>${found} of ${last.length} entries are in the atlas.</p><div class="tw"><table><thead><tr><th>Query</th><th>Gene</th><th class="num">Reported</th><th class="num">Testable</th><th class="num">Supported</th><th class="num">Matched null</th><th>Class or verdict</th><th>ClinVar</th></tr></thead><tbody>${
        last.map(({ q, x, r }) => !x ? `<tr><td>${esc(q)}</td><td colspan="7" class="nt">not in the atlas</td></tr>` : r
          ? `<tr><td><a href="#check=rs${r.rs}|${encodeURIComponent(q.split(/[,;\t]/).slice(1).join(',').trim())}">${esc(q)}</a></td><td><i>${esc(x.gene)}</i></td><td class="num">${r.R}</td><td class="num">${r.T}</td><td class="num">${r.T ? pct(r.S / r.T) : '—'}</td><td class="num">${pct(x.null)}</td><td>${!r.R ? 'no claim on this trait' : !r.T ? 'not testable' : r.S / r.T >= 0.5 ? 'supported' : r.S ? 'partly supported' : 'not supported'}</td><td>${esc(x.clinvar)}</td></tr>`
          : `<tr><td><a href="#rs${x.rsid}">${esc(q)}</a></td><td><i>${esc(x.gene)}</i></td><td class="num">${num(x.reported)}</td><td class="num">${num(x.testable)}</td><td class="num">${pct(x.supported)}</td><td class="num">${pct(x.null)}</td><td>${cls(x) || 'fewer than 5 reported'}</td><td>${esc(x.clinvar)}</td></tr>`).join('')}</tbody></table></div>`;
      $('#bdl').disabled = !last.length;
    });
  }

  // ---------------------------------------------------------------- ClinVar
  async function showClinvar() {
    out.innerHTML = '<p class="muted">Loading the ClinVar table…</p>';
    if (!CV) {
      const txt = await fetch('atlas/Data_S3_clinvar_risk_variants.tsv' + DV).then((r) => r.text());
      const [h, ...rows] = txt.trim().split('\n').map((l) => l.split('\t'));
      CV = rows.map((r) => Object.fromEntries(h.map((k, i) => [k, r[i]])));
    }
    const classes = ['all', 'mostly (>=50%)', 'partly (5-50%)', 'rarely (<5%)', 'fewer than 5 reported associations'];
    out.innerHTML = `<div class="row"><h2>ClinVar risk variants</h2><button class="sec" id="dl">Download table</button></div>
      <p class="muted">Variants listed in ClinVar (release of ${esc(META.clinvar)}) as risk factors, associations or protective variants that have reported associations in the held-out papers of the study, with the share supported genome-wide and whether genome-wide evidence supports the ClinVar condition itself (Data S3 of the paper).</p>
      <form id="cvf"><div class="fld"><label for="cvq">Filter by rsID, gene or condition</label><input id="cvq" placeholder="IL6, PADI4, rheumatoid"></div>
      <div class="fld" style="flex:0 1 260px"><label for="cvc">Support class</label><select id="cvc">${classes.map((c) => `<option>${c}</option>`).join('')}</select></div></form>
      <div class="tw" id="cvt"></div>`;
    let key = 'reported_associations', dir = -1;
    const draw = () => {
      const s = $('#cvq').value.toLowerCase(), c = $('#cvc').value;
      const L = CV.filter((r) => (c === 'all' || r.support_class === c)
        && (!s || (r.rsid + ' ' + r.gene_clinvar + ' ' + r.clinvar_conditions).toLowerCase().includes(s)))
        .sort((a, b) => { const p = +a[key], q = +b[key]; return (isNaN(p) ? String(a[key]).localeCompare(b[key]) : p - q) * dir; });
      $('#cvt').innerHTML = `<table><thead><tr><th data-k="rsid">Variant</th><th data-k="gene_clinvar">Gene</th><th>ClinVar assertion</th><th class="num" data-k="reported_associations">Reported</th><th class="num" data-k="share_supported">Supported</th><th class="num" data-k="matched_null">Matched null</th><th>Class</th><th>Support for the ClinVar condition</th></tr></thead><tbody>${
        L.map((r) => `<tr><td><a href="#${r.rsid}">${esc(r.rsid)}</a></td><td><i>${esc(r.gene_clinvar)}</i></td><td class="t" title="${esc(r.clinvar_conditions)}">${esc(r.clinvar_significance)}</td><td class="num">${esc(r.reported_associations)}</td><td class="num">${pct(+r.share_supported)}</td><td class="num">${pct(+r.matched_null)}</td><td>${esc(r.support_class)}</td><td>${r.clinvar_condition_in_ontology === 'yes' ? (r.genome_wide_support_for_clinvar_condition === 'yes' ? '<span class="yes">yes</span>' : '<span class="no">no</span>') : '<span class="nt">condition not mapped</span>'}</td></tr>`).join('')}</tbody></table><p class="muted" style="font-size:12.5px">${L.length} of ${CV.length} variants. Click a column heading to sort.</p>`;
      document.querySelectorAll('#cvt th[data-k]').forEach((th) => th.addEventListener('click', () => { if (key === th.dataset.k) dir = -dir; else { key = th.dataset.k; dir = -1; } draw(); }));
    };
    $('#cvq').addEventListener('input', draw); $('#cvc').addEventListener('change', draw); draw();
    $('#dl').addEventListener('click', () => { const h = Object.keys(CV[0]); download('clinvar_risk_variants.csv', [h].concat(CV.map((r) => h.map((k) => r[k])))); });
  }

  // ---------------------------------------------------------------- about
  function showAbout() {
    out.innerHTML = `<h2>About and download</h2>
      <p>The atlas records what each PubMed abstract that names a dbSNP identifier claimed about a variant and a disease or trait, links the trait to the Experimental Factor Ontology, and compares each claim with genome-wide significant associations of the same variant or a proxy in the NHGRI-EBI GWAS Catalog, against frequency-matched random variants. Claims were read by an open-weight language model run locally; each claim quotes its sentence from the abstract and was checked on a blinded sample.</p>
      <h3>Download</h3>
      <ul style="margin:0;padding-left:20px;display:flex;flex-direction:column;gap:4px">
        <li><a href="atlas/Data_S1_variant_trait_claim_atlas.tsv.gz">Full atlas of claim records</a> (tab-separated, compressed, 11 MB)</li>
        <li><a href="atlas/Data_S3_clinvar_risk_variants.tsv">ClinVar risk variants</a> (tab-separated)</li>
        <li><a href="https://github.com/zxu69461-creator/variant-claim-atlas/tree/main/code">Analysis code and the dated analysis plan</a></li>
        <li>Archived release with all files: <a href="https://doi.org/10.5281/zenodo.24314560">doi:10.5281/zenodo.24314560</a></li>
      </ul>
      <h3>Programmatic access</h3>
      <p>The tool reads static JSON files that can be fetched directly: <code>data/variants.json</code> (one row per variant), <code>data/traits.json</code> (one row per trait) and <code>data/claims/NN.json</code> (claim records, file <code>NN</code> = rsID number modulo 64, two digits). Each file holds <code>cols</code> and <code>rows</code>.</p>
      <h3>Cite</h3>
      <p>Hu Y, Xu Z, Jiang S, Huang L, Pan M, Luo S, Zhou H, Guo X. A computational atlas calibrating 201,720 published variant–trait claims against genome-wide evidence and ClinVar assertions. Manuscript under review, 2026. Data and code: doi:10.5281/zenodo.24314560.</p>`;
  }

  // ---------------------------------------------------------------- routing
  function tab(t) {
    document.querySelectorAll('nav.tabs a').forEach((a) => { if (a.dataset.tab === t) a.setAttribute('aria-current', 'page'); else a.removeAttribute('aria-current'); });
    $('#searchcard').hidden = t !== 'lookup';
  }
  function route() {
    const h = decodeURIComponent(location.hash.slice(1)); let m;
    if (h === 'check' || h.startsWith('check=')) { tab('check'); return showCheck(h.slice(6)); }
    if (h === 'batch') { tab('batch'); return showBatch(); }
    if (h === 'clinvar') { tab('clinvar'); return showClinvar(); }
    if (h === 'about') { tab('about'); return showAbout(); }
    tab('lookup');
    if (!h || h === 'lookup') return home();
    if ((m = h.match(/^rs(\d+)$/i))) return showVariant(+m[1]);
    if ((m = h.match(/^gene=(.+)$/))) return showGene(m[1]);
    if ((m = h.match(/^trait=(.+)$/))) return showTrait(m[1]);
    if ((m = h.match(/^q=(.+)$/))) return interpret(m[1]);
    return home();
  }
  function interpret(q) {
    q = q.trim();
    if (/^rs\d+$/i.test(q)) location.hash = q.toLowerCase();
    else if (byGene.has(q.toUpperCase())) location.hash = 'gene=' + q.toUpperCase();
    else { const t = T.find((r) => r.trait.toLowerCase() === q.toLowerCase()); location.hash = t ? 'trait=' + encodeURIComponent(t.trait) : 'q=' + encodeURIComponent(q); if (!t) showSearch(q); }
  }
  $('#f').addEventListener('submit', (e) => { e.preventDefault(); const q = $('#q').value; if (q.trim()) interpret(q); });
  window.addEventListener('hashchange', route);
  load().catch((e) => { out.innerHTML = '<p>The atlas data could not be loaded: ' + esc(e.message) + '</p>'; });
})();

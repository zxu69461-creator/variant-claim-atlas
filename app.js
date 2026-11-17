/* Variant–Trait Claim Atlas: static lookup over data/*.json (built by scripts/34_lookup_data.py). */
(function () {
  'use strict';
  const $ = (s) => document.querySelector(s);
  const out = $('#out');
  const esc = (s) => String(s == null ? '' : s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const cap = (s) => (s ? s[0].toUpperCase() + s.slice(1) : s);
  const pct = (x) => (x == null ? '—' : (100 * x).toFixed(x < 0.1 ? 1 : 0) + '%');
  const num = (x) => (x == null ? '—' : x.toLocaleString('en-US'));
  let V, T, META, byRs, byGene, shards = {};

  const DV = '?v=1.0.1';                                   // data version: bump when data/ changes
  const obj = (t) => t.rows.map((r) => Object.fromEntries(t.cols.map((c, i) => [c, r[i]])));
  async function load() {
    const [v, t, m] = await Promise.all(['variants', 'traits', 'meta'].map((n) => fetch('data/' + n + '.json' + DV).then((r) => r.json())));
    V = obj(v); T = obj(t); META = m;
    byRs = new Map(V.map((x) => [x.rsid, x]));
    byGene = new Map();
    V.forEach((x) => { if (x.gene) { const g = x.gene.toUpperCase(); (byGene.get(g) || byGene.set(g, []).get(g)).push(x); } });
    $('#intro').textContent = `${num(META.claims)} published claims from ${num(META.papers)} PubMed papers on ${num(META.variants)} variants, each set against genome-wide significant associations in the GWAS Catalog (release of ${META.gwas_catalog}).`;
    $('#rel').textContent = `Sources: PubMed ${META.pubmed}; GWAS Catalog ${META.gwas_catalog}; Experimental Factor Ontology ${META.efo}; ClinVar ${META.clinvar}. Built ${META.built}.`;
    route();
  }
  async function claims(rs) {
    const s = rs % META.shards;
    if (!shards[s]) shards[s] = obj(await fetch('data/claims/' + String(s).padStart(2, '0') + '.json' + DV).then((r) => r.json()));
    return shards[s].filter((c) => c.rsid === rs);
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
    const tt = x.top_traits.map(([t, k, tt, s]) => `<tr><td class="t"><a href="#trait=${encodeURIComponent(t)}">${esc(cap(t))}</a></td><td class="num">${k}</td><td class="num">${tt}</td><td class="num">${tt ? pct(s / tt) : '—'}</td></tr>`).join('');
    out.innerHTML = `<h2>rs${rs}${x.gene ? ` <span class="muted" style="font-weight:400">· <i>${esc(x.gene)}</i></span>` : ''}</h2>
      ${kpis(x, extra)}${meter(x.supported, x.null)}
      <h3>Support by stated sample size</h3><div class="sizes">${sizes}</div>
      <h3>Most-claimed traits</h3><div class="tw"><table><thead><tr><th>Trait</th><th class="num">Reported</th><th class="num">Testable</th><th class="num">Supported</th></tr></thead><tbody>${tt}</tbody></table></div>
      <h3>Every claim</h3><div class="tw" id="cl"><p class="muted">Loading claims…</p></div>`;
    const cl = await claims(rs);
    const st = (c) => (!c.testable ? '<span class="nt">not testable</span>' : c.supported ? `<span class="yes">yes${c.earlier ? ', earlier' : ''}</span>` : '<span class="no">no</span>');
    $('#cl').innerHTML = `<table><thead><tr><th>Year</th><th>PMID</th><th>Trait</th><th>Claim</th><th class="num">n</th><th>Design</th><th>Population</th><th>Genome-wide support</th></tr></thead><tbody>${
      cl.map((c) => `<tr><td>${c.year}</td><td><a href="https://pubmed.ncbi.nlm.nih.gov/${c.pmid}/" target="_blank" rel="noopener">${c.pmid}</a></td><td class="t">${esc(cap(c.trait))}</td><td>${c.status === 'reported' ? 'associated' : 'no association'}</td><td class="num">${num(c.n)}</td><td>${esc(c.design)}</td><td>${esc(c.population)}</td><td>${c.status === 'reported' ? st(c) : '<span class="nt">—</span>'}</td></tr>`).join('')}</tbody></table>`;
  }

  function showTrait(name) {
    const t = T.find((r) => r.trait === name);
    if (!t) return showSearch(name);
    const rows = t.top_variants.map(([rs, g, k, tt, s]) => `<tr><td><a href="#rs${rs}">rs${rs}</a></td><td><i>${esc(g)}</i></td><td class="num">${k}</td><td class="num">${tt}</td><td class="num">${tt ? pct(s / tt) : '—'}</td></tr>`).join('');
    out.innerHTML = `<h2>${esc(cap(t.trait))} <span class="muted" style="font-weight:400;font-size:14px">· ${esc(t.area)}</span></h2>${kpis(t)}${meter(t.supported, t.null)}
      <h3>Most-claimed variants</h3><div class="tw"><table><thead><tr><th>Variant</th><th>Gene</th><th class="num">Reported</th><th class="num">Testable</th><th class="num">Supported</th></tr></thead><tbody>${rows}</tbody></table></div>`;
  }

  function showGene(g) {
    const L = (byGene.get(g.toUpperCase()) || []).slice().sort((a, b) => b.reported - a.reported);
    if (!L.length) return showSearch(g);
    out.innerHTML = `<h2><i>${esc(g.toUpperCase())}</i> <span class="muted" style="font-weight:400;font-size:14px">· ${L.length} variant${L.length > 1 ? 's' : ''} (gene from dbSNP or the GWAS Catalog)</span></h2>
      <div class="tw"><table><thead><tr><th>Variant</th><th class="num">Papers</th><th class="num">Reported</th><th class="num">Testable</th><th class="num">Supported</th><th class="num">Matched null</th><th>ClinVar</th></tr></thead><tbody>${
      L.slice(0, 200).map((x) => `<tr><td><a href="#rs${x.rsid}">rs${x.rsid}</a></td><td class="num">${num(x.papers)}</td><td class="num">${num(x.reported)}</td><td class="num">${num(x.testable)}</td><td class="num">${pct(x.supported)}</td><td class="num">${pct(x.null)}</td><td>${esc(x.clinvar)}</td></tr>`).join('')}</tbody></table></div>`;
  }

  function showSearch(q) {
    const s = q.toLowerCase();
    const hits = T.filter((t) => t.trait.toLowerCase().includes(s)).slice(0, 40);
    out.innerHTML = hits.length
      ? `<h2>Traits matching “${esc(q)}”</h2><div class="list">${hits.map((t) => `<a href="#trait=${encodeURIComponent(t.trait)}">${esc(cap(t.trait))} <span class="muted">· ${num(t.reported)} reported associations</span></a>`).join('')}</div>`
      : `<p>Nothing matches “${esc(q)}”. Try an rsID (rs1801133), a gene symbol (FTO) or part of a disease name (asthma).</p>`;
  }

  function route() {
    const h = decodeURIComponent(location.hash.slice(1));
    if (!h) { out.innerHTML = '<p class="muted">Search an rsID, a gene or a trait, or pick an example above.</p>'; return; }
    let m;
    if ((m = h.match(/^rs(\d+)$/i))) return showVariant(+m[1]);
    if ((m = h.match(/^gene=(.+)$/))) return showGene(m[1]);
    if ((m = h.match(/^trait=(.+)$/))) return showTrait(m[1]);
    if ((m = h.match(/^q=(.+)$/))) return interpret(m[1]);
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

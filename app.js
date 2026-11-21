/* Variant–Trait Claim Atlas: a static web tool over data/*.json (built by code/34_lookup_data.py and code/35_lookup_extras.py) and
   atlas/Data_S3 (ClinVar table).
   Views (URL hash): #lookup (default; #rs123, #gene=SYMBOL, #trait=Label, #q=text), #check (#check=rs123|trait), #paper (#paper=PMID),
   #batch, #clinvar, #help, #about.
   Everything runs in the browser. Pasted lists, uploaded files and pasted text never leave it; screening a paper by PMID sends only
   the PMID to Europe PMC to retrieve its title and abstract. */
(function () {
  'use strict';
  const APP = '1.2.2', RELEASED = '21 November 2026';
  const DV = '?v=' + APP;                                  // data version: bump when data/ changes
  const $ = (s) => document.querySelector(s);
  const out = $('#out');
  const esc = (s) => String(s == null ? '' : s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const cap = (s) => (s ? s[0].toUpperCase() + s.slice(1) : s);
  const pct = (x) => (x == null || isNaN(x) ? '—' : (100 * x).toFixed(x < 0.1 && x > 0 ? 1 : 0) + '%');
  const num = (x) => (x == null ? '—' : Number(x).toLocaleString('en-US'));
  const reEsc = (s) => s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  const ext = (u, t) => `<a href="${u}" target="_blank" rel="noopener">${t}</a>`;
  let V, T, META, byRs, byGene, byTrait, CV = null, LEX = null, PAP = null, shards = {};

  const obj = (t) => t.rows.map((r) => Object.fromEntries(t.cols.map((c, i) => [c, r[i]])));
  const getJSON = (p) => fetch(p + DV).then((r) => { if (!r.ok) throw new Error(p + ': HTTP ' + r.status); return r.json(); });
  const cls = (x) => (x.reported < 5 || x.testable === 0 ? '' : x.supported >= 0.5 ? 'mostly supported' : x.supported >= 0.05 ? 'partly supported' : 'rarely supported');

  async function load() {
    const [v, t, m] = await Promise.all(['variants', 'traits', 'meta'].map((n) => getJSON('data/' + n + '.json')));
    V = obj(v); T = obj(t); META = m;
    byRs = new Map(V.map((x) => [x.rsid, x]));
    byTrait = new Map(T.map((x) => [x.trait, x]));
    byGene = new Map();
    V.forEach((x) => { if (x.gene) { const g = x.gene.toUpperCase(); (byGene.get(g) || byGene.set(g, []).get(g)).push(x); } });
    $('#traitlist').innerHTML = T.slice(0, 3000).map((r) => `<option value="${esc(cap(r.trait))}">`).join('');
    $('#intro').textContent = `The 201,720 disease or trait claims of the atlas include ${num(META.claims)} that report an association or its absence, from ${num(META.papers)} PubMed papers on ${num(META.variants)} variants; each is set against genome-wide significant associations in the GWAS Catalog (release of ${META.gwas_catalog}).`;
    $('#rel').textContent = `Sources: PubMed ${META.pubmed}; GWAS Catalog ${META.gwas_catalog}; Experimental Factor Ontology ${META.efo}; ClinVar ${META.clinvar}. Data built ${META.built}.`;
    $('#ver').innerHTML = `Version ${APP} (${RELEASED}). No cookies, no analytics, no registration; what you paste or upload stays in your browser. Questions and problem reports: <a href="https://github.com/zxu69461-creator/variant-claim-atlas/issues">GitHub issues</a>.`;
    route();
  }
  async function claims(rs) {
    const s = rs % META.shards;
    if (!shards[s]) shards[s] = obj(await getJSON('data/claims/' + String(s).padStart(2, '0') + '.json'));
    return shards[s].filter((c) => c.rsid === rs);
  }
  async function lexicon() {
    if (!LEX) {
      const L = obj(await getJSON('data/lexicon.json'));
      LEX = { byTrait: new Map(L.map((r) => [r.trait, r])), byPhrase: new Map() };
      L.forEach((r) => r.phrases.forEach((p) => { if (!LEX.byPhrase.has(p)) LEX.byPhrase.set(p, []); LEX.byPhrase.get(p).push(r.trait); }));
    }
    return LEX;
  }
  async function papers() {
    if (!PAP) PAP = new Map(obj(await getJSON('data/papers.json')).map((r) => [r.pmid, r]));
    return PAP;
  }
  function download(name, rows) {
    const q = (v) => { const s = v == null ? '' : String(v); return /[",\n\t]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; };
    const blob = new Blob([rows.map((r) => r.map(q).join(',')).join('\n') + '\n'], { type: 'text/csv' });
    const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500);
  }

  // ---------------------------------------------------------------- shared pieces
  // Traits named by a query: trait names that contain it, and traits for which abstracts used the query as a phrase (T2D, BMI, …).
  async function resolveTraits(q) {
    const s = q.toLowerCase().trim(); const set = new Set();
    if (!s) return set;
    // names match at the start of a word; short queries (abbreviations) only as a whole word, so that CAD does not match 'academic'
    const re = new RegExp('(?:^|[^a-z0-9])' + reEsc(s) + (s.length <= 4 ? '(?![a-z0-9])' : ''));
    T.forEach((r) => { if (re.test(r.trait.toLowerCase())) set.add(r.trait); });
    const L = await lexicon();
    (L.byPhrase.get(s) || []).forEach((t) => set.add(t));
    return set;
  }
  // The verdict for a set of claims on one variant (the rules are those of the Help page).
  function verdict(M) {
    const R = M.filter((c) => c.status === 'reported'), Tt = R.filter((c) => c.testable), S = Tt.filter((c) => c.supported), E = S.filter((c) => c.earlier);
    const yrs = S.map((c) => c.first_gw_year).filter(Boolean); const y = yrs.length ? Math.min(...yrs) : null;
    const kind = !R.length ? 'none' : !Tt.length ? 'nt' : S.length / Tt.length >= 0.5 ? 'sup' : S.length ? 'part' : 'no';
    const label = { none: 'no reported association', nt: 'not testable', sup: 'supported', part: 'partly supported', no: 'not supported' }[kind];
    return { R, Tt, S, E, A: M.length - R.length, y, kind, label };
  }
  const vcls = (k) => (k === 'sup' || k === 'part' ? 'yes' : k === 'no' ? 'no' : 'nt');
  function meter(sup, nul) {
    if (sup == null) return '';
    return `<div class="meter" role="img" aria-label="supported ${pct(sup)}, matched null ${pct(nul)}"><i style="width:${Math.min(100, 100 * sup)}%"></i>${nul != null ? `<u style="left:${Math.min(100, 100 * nul)}%"></u>` : ''}</div>
      <p class="muted" style="font-size:12.5px">Bar: share of testable reported associations with genome-wide support. Tick: matched null (frequency-matched random variants, same traits).</p>`;
  }
  const mini = (s) => (s == null ? '' : `<span class="mini" aria-hidden="true"><i style="width:${Math.min(100, 100 * s)}%"></i></span>`);
  function kpis(x, extra) {
    return `<div class="kpis">
      <div class="kpi"><b>${num(x.papers)}</b><span>papers</span></div>
      <div class="kpi"><b>${num(x.reported)}</b><span>reported associations</span></div>
      <div class="kpi"><b>${num(x.absent)}</b><span>reports of no association</span></div>
      <div class="kpi"><b>${num(x.testable)}</b><span>testable reported associations</span></div>
      <div class="kpi s"><b>${pct(x.supported)}</b><span>with genome-wide support</span></div>
      <div class="kpi n"><b>${pct(x.null)}</b><span>matched null</span></div>${extra || ''}</div>`;
  }
  const rsLinks = (rs) => `<div class="links"><span>Elsewhere:</span>${ext(`https://www.ncbi.nlm.nih.gov/snp/rs${rs}`, 'dbSNP')}${ext(`https://www.ebi.ac.uk/gwas/variants/rs${rs}`, 'GWAS Catalog')}${ext(`https://www.ncbi.nlm.nih.gov/clinvar/?term=rs${rs}`, 'ClinVar')}${ext(`https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?text=rs${rs}`, 'LitVar 2.0')}${ext(`https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=rs${rs}`, 'Ensembl')}${ext(`https://platform.opentargets.org/search?q=rs${rs}`, 'Open Targets')}</div>`;
  const supCell = (c) => (c.status !== 'reported' ? '<span class="nt">—</span>' : !c.testable ? '<span class="nt">not testable</span>' : c.supported ? `<span class="yes">yes${c.earlier ? ', earlier' : ''}</span>` : '<span class="no">no</span>');
  const claimRow = (c) => `<tr><td>${c.year}</td><td><a href="https://pubmed.ncbi.nlm.nih.gov/${c.pmid}/" target="_blank" rel="noopener">${c.pmid}</a></td><td class="t">${esc(cap(c.trait))}</td><td>${c.status === 'reported' ? 'associated' : 'no association'}</td><td class="num">${num(c.n)}</td><td>${esc(c.design)}</td><td>${esc(c.population)}</td><td>${supCell(c)}</td></tr>`;
  const claimHead = '<tr><th>Year</th><th>PMID</th><th>Trait</th><th>Claim</th><th class="num">n</th><th>Design</th><th>Population</th><th>Genome-wide support</th></tr>';
  const claimCsv = (cl) => [['rsid', 'pmid', 'year', 'trait', 'claim', 'stated_n', 'design', 'population', 'testable', 'genome_wide_support', 'support_earlier_than_claim', 'first_genome_wide_year']]
    .concat(cl.map((c) => ['rs' + c.rsid, c.pmid, c.year, c.trait, c.status === 'reported' ? 'associated' : 'no association', c.n, c.design, c.population, c.testable, c.status === 'reported' ? c.supported : '', c.status === 'reported' ? c.earlier : '', c.first_gw_year]));
  const notInAtlas = (rs) => `<div class="verdict nt"><b>rs${rs} has no disease or trait claim in the atlas.</b><span>The atlas holds variants named by rsID in a PubMed title or abstract (PubMed ${esc(META.pubmed)}) with a claim about a disease or trait. Genome-wide associations of this variant, if any, are listed in the GWAS Catalog.</span></div>${rsLinks(rs)}`;

  // Claims per year: supported, not supported and not testable reported associations, and reports of no association, with the year of
  // the first genome-wide report for a claimed trait.
  function timeline(cl, gw) {
    const ys = cl.map((c) => c.year); const y0 = Math.min(...ys, gw || 9999), y1 = Math.max(...ys);
    const n = y1 - y0 + 1, W = 760, H = 170, L = 34, R = 8, Tp = 10, B = 24, bw = (W - L - R) / n;
    const by = {}; for (let y = y0; y <= y1; y++) by[y] = [0, 0, 0, 0];
    cl.forEach((c) => { const k = c.status !== 'reported' ? 3 : !c.testable ? 2 : c.supported ? 0 : 1; by[c.year][k]++; });
    const mx = Math.max(1, ...Object.values(by).map((v) => v[0] + v[1] + v[2] + v[3]));
    const sy = (v) => ((H - Tp - B) * v) / mx;
    let g = '';
    for (let y = y0; y <= y1; y++) {
      const [s, x, q, a] = by[y]; let base = H - B; const xx = L + (y - y0) * bw + 0.12 * bw, w = Math.max(1, 0.76 * bw);
      [[s, 's'], [x, 'x'], [q, 'q']].forEach(([v, c]) => { if (v) { const h = sy(v); base -= h; g += `<rect class="${c}" x="${xx.toFixed(1)}" y="${base.toFixed(1)}" width="${w.toFixed(1)}" height="${h.toFixed(1)}"><title>${y}: ${v}</title></rect>`; } });
      if (a) { const h = sy(a); base -= h; g += `<rect class="a" x="${(xx + 0.5).toFixed(1)}" y="${(base + 0.5).toFixed(1)}" width="${Math.max(0.5, w - 1).toFixed(1)}" height="${Math.max(0.5, h - 1).toFixed(1)}"><title>${y}: ${a} reports of no association</title></rect>`; }
    }
    const step = n > 30 ? 10 : n > 12 ? 5 : n > 6 ? 2 : 1; let ax = '';
    for (let y = Math.ceil(y0 / step) * step; y <= y1; y += step) ax += `<text x="${(L + (y - y0 + 0.5) * bw).toFixed(1)}" y="${H - 8}" text-anchor="middle">${y}</text>`;
    const yt = [0, Math.round(mx / 2), mx].filter((v, i, a) => a.indexOf(v) === i).map((v) => `<text x="${L - 6}" y="${(H - B - sy(v) + 4).toFixed(1)}" text-anchor="end">${v}</text>`).join('');
    const gl = gw ? `<line class="gw" x1="${(L + (gw - y0) * bw).toFixed(1)}" x2="${(L + (gw - y0) * bw).toFixed(1)}" y1="${Tp}" y2="${H - B}"/><text x="${(L + (gw - y0) * bw + ((gw - y0) / n > 0.6 ? -4 : 4)).toFixed(1)}" y="${Tp + 9}" text-anchor="${(gw - y0) / n > 0.6 ? 'end' : 'start'}" style="fill:var(--first)">first genome-wide report, ${gw}</text>` : '';
    return `<svg class="tl" viewBox="0 0 ${W} ${H}" role="img" aria-label="Claims per year from ${y0} to ${y1}"><line class="ax" x1="${L}" x2="${W - R}" y1="${H - B}" y2="${H - B}"/>${g}${gl}${ax}${yt}</svg>
      <div class="key"><span><i style="background:var(--sup)"></i>supported</span><span><i style="background:var(--warn)"></i>not supported</span><span><i style="background:var(--null);opacity:.55"></i>not testable</span><span><i style="border:1px solid var(--muted)"></i>no association reported</span>${gw ? '<span><i style="background:var(--first);height:2px;vertical-align:3px"></i>first genome-wide report for a claimed trait</span>' : ''}</div>`;
  }

  // ---------------------------------------------------------------- look up
  async function showVariant(rs) {
    const x = byRs.get(rs);
    if (!x) { out.innerHTML = notInAtlas(rs); return; }
    const extra = `<div class="kpi"><b>${num(x.median_n)}</b><span>median stated sample size</span></div>
      <div class="kpi"><b>${x.first_gw_year || '—'}</b><span>first genome-wide report for a claimed trait</span></div>
      ${x.clinvar ? `<div class="kpi"><b style="font-size:16px">${esc(x.clinvar)}</b><span>ClinVar assertion</span></div>` : ''}`;
    const sizes = META.size_classes.map((lab, i) => {
      const [k, s] = x.size_classes[i]; const sh = k ? s / k : null;
      return `<div class="sz"><span>${lab} participants: <b>${k ? pct(sh) : '—'}</b> <span class="muted">(${s} of ${k})</span></span>
        <div class="meter" style="height:8px"><i style="width:${sh == null ? 0 : 100 * sh}%"></i></div></div>`;
    }).join('');
    const tt = x.top_traits.map(([t, k, tt, s]) => `<tr><td class="t"><a href="#check=rs${rs}|${encodeURIComponent(t)}">${esc(cap(t))}</a></td><td class="num">${k}</td><td class="num">${tt}</td><td class="num">${tt ? mini(s / tt) + pct(s / tt) : '—'}</td></tr>`).join('');
    const c = cls(x);
    out.innerHTML = `<h2>rs${rs}${x.gene ? ` <span class="muted" style="font-weight:400">· <a href="#gene=${encodeURIComponent(x.gene)}"><i>${esc(x.gene)}</i></a></span>` : ''}${c ? ` <span class="muted" style="font-weight:400;font-size:14px">· ${c}</span>` : ''}</h2>
      ${rsLinks(rs)}
      ${kpis(x, extra)}${meter(x.supported, x.null)}
      <h3>Claims per year</h3><div id="tl"><p class="muted">Loading claims…</p></div>
      <h3>Support by stated sample size</h3><div class="sizes">${sizes}</div>
      <h3>Most-claimed traits <span style="text-transform:none;letter-spacing:0;font-weight:400">(click to check the claim)</span></h3><div class="tw"><table><thead><tr><th>Trait</th><th class="num">Reported</th><th class="num">Testable</th><th class="num">Supported</th></tr></thead><tbody>${tt}</tbody></table></div>
      <div class="row"><h3>Every claim</h3><div style="display:flex;gap:8px;flex-wrap:wrap"><input id="cf" placeholder="Filter by trait, year, design" aria-label="Filter claims" style="width:220px"><button class="sec" id="dl">Download CSV</button></div></div>
      <div class="tw" id="cl"><p class="muted">Loading claims…</p></div>`;
    const cl = await claims(rs);
    $('#tl').innerHTML = timeline(cl, x.first_gw_year);
    const draw = (f) => {
      const s = f.toLowerCase(); const L = s ? cl.filter((c) => (c.trait + ' ' + c.year + ' ' + c.design + ' ' + c.population).toLowerCase().includes(s)) : cl;
      $('#cl').innerHTML = `<table><thead>${claimHead}</thead><tbody>${L.map(claimRow).join('')}</tbody></table><p class="muted" style="font-size:12.5px">${L.length} of ${cl.length} claims</p>`;
    };
    draw('');
    $('#cf').addEventListener('input', (e) => draw(e.target.value));
    $('#dl').addEventListener('click', () => download(`rs${rs}_claims.csv`, claimCsv(cl)));
  }
  async function showTrait(name) {
    const t = byTrait.get(name) || T.find((r) => r.trait.toLowerCase() === name.toLowerCase());
    if (!t) return showSearch(name);
    const rows = t.top_variants.map(([rs, g, k, tt, s]) => `<tr><td><a href="#rs${rs}">rs${rs}</a></td><td>${g ? `<a href="#gene=${encodeURIComponent(g)}"><i>${esc(g)}</i></a>` : ''}</td><td class="num">${k}</td><td class="num">${tt}</td><td class="num">${tt ? mini(s / tt) + pct(s / tt) : '—'}</td><td><a href="#check=rs${rs}|${encodeURIComponent(t.trait)}">check</a></td></tr>`).join('');
    out.innerHTML = `<h2>${esc(cap(t.trait))}</h2><div id="tlinks"></div>${kpis(t)}${meter(t.supported, t.null)}
      <h3>Most-claimed variants</h3><div class="tw"><table><thead><tr><th>Variant</th><th>Gene</th><th class="num">Reported</th><th class="num">Testable</th><th class="num">Supported</th><th></th></tr></thead><tbody>${rows}</tbody></table></div>`;
    const L = await lexicon(); const e = (L.byTrait.get(t.trait) || {}).efo;
    if (e) $('#tlinks').innerHTML = `<div class="links"><span>Ontology term ${esc(e.replace('_', ':'))}:</span>${ext(`https://www.ebi.ac.uk/ols4/search?q=${encodeURIComponent(e)}`, 'Ontology Lookup Service')}${ext(`https://www.ebi.ac.uk/gwas/efotraits/${encodeURIComponent(e)}`, 'GWAS Catalog trait page')}</div>`;
  }
  function showGene(g) {
    const G = g.toUpperCase(); const L = (byGene.get(G) || []).slice().sort((a, b) => b.reported - a.reported);
    if (!L.length) return showSearch(g);
    out.innerHTML = `<div class="row"><h2><i>${esc(G)}</i> <span class="muted" style="font-weight:400;font-size:14px">· ${L.length} variant${L.length > 1 ? 's' : ''} (gene from dbSNP or the GWAS Catalog)</span></h2><button class="sec" id="dl">Download CSV</button></div>
      <div class="links"><span>Elsewhere:</span>${ext(`https://www.ncbi.nlm.nih.gov/gene/?term=${encodeURIComponent(G)}%5Bsym%5D+AND+human%5Borgn%5D`, 'NCBI Gene')}${ext(`https://www.ebi.ac.uk/gwas/genes/${encodeURIComponent(G)}`, 'GWAS Catalog gene page')}${ext(`https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?text=${encodeURIComponent(G)}`, 'LitVar 2.0')}</div>
      <div class="tw"><table><thead><tr><th>Variant</th><th class="num">Papers</th><th class="num">Reported</th><th class="num">Testable</th><th class="num">Supported</th><th class="num">Matched null</th><th>Class</th><th>ClinVar</th></tr></thead><tbody>${
      L.slice(0, 300).map((x) => `<tr><td><a href="#rs${x.rsid}">rs${x.rsid}</a></td><td class="num">${num(x.papers)}</td><td class="num">${num(x.reported)}</td><td class="num">${num(x.testable)}</td><td class="num">${x.supported == null ? '—' : mini(x.supported) + pct(x.supported)}</td><td class="num">${pct(x.null)}</td><td>${cls(x)}</td><td>${esc(x.clinvar)}</td></tr>`).join('')}</tbody></table></div>
      ${L.length > 300 ? `<p class="muted" style="font-size:12.5px">The 300 most-claimed of ${L.length} variants are shown; the CSV holds all.</p>` : ''}`;
    $('#dl').addEventListener('click', () => download(`${G}_variants.csv`, batchCsv(L.map((x) => ({ q: 'rs' + x.rsid, x })))));
  }
  async function showSearch(q) {
    const s = q.toLowerCase();
    const set = await resolveTraits(q);
    const hits = T.filter((t) => set.has(t.trait)).slice(0, 40);
    out.innerHTML = hits.length
      ? `<h2>Traits matching “${esc(q)}”</h2><div class="list">${hits.map((t) => `<a href="#trait=${encodeURIComponent(t.trait)}">${esc(cap(t.trait))} <span class="muted">· ${num(t.reported)} reported associations</span></a>`).join('')}</div>`
      : `<p>Nothing matches “${esc(s)}”. Try an rsID (rs1801133), a gene symbol (FTO) or part of a disease name (asthma).</p>`;
  }
  function home() {
    out.innerHTML = `<h2>Check a published variant–trait claim against genome-wide evidence</h2>
      <p>Before you cite, curate, review or follow up a published claim that a variant is associated with a disease or trait, see how often claims about that variant and trait are backed by genome-wide significant associations, compared with frequency-matched random variants, and in studies of which size.</p>
      <div class="kpis">
        <div class="kpi"><b>${num(META.claims)}</b><span>claims of an association or its absence</span></div>
        <div class="kpi"><b>${num(META.papers)}</b><span>PubMed papers</span></div>
        <div class="kpi"><b>${num(META.variants)}</b><span>variants</span></div>
        <div class="kpi"><b>${num(META.traits)}</b><span>ontology-linked traits</span></div>
      </div>
      <div class="grid3">
        <a class="tool" href="#rs1800795"><b>Look up</b><span>A variant, a gene or a trait: support of its claims against the matched null, by sample size and year, with every claim linked to PubMed.</span></a>
        <a class="tool" href="#check=rs1800795|type 2 diabetes"><b>Check a claim</b><span>One variant and one trait (names and abbreviations such as T2D): does genome-wide evidence support the association, and since when?</span></a>
        <a class="tool" href="#paper=42749416"><b>Screen a paper</b><span>Enter a PubMed ID, or paste an abstract or a manuscript: every variant–trait pair named in it is checked.</span></a>
        <a class="tool" href="#batch"><b>Batch check</b><span>Paste or upload up to 10,000 variants or variant–trait pairs and download the calibration as CSV.</span></a>
        <a class="tool" href="#clinvar"><b>ClinVar risk variants</b><span>The 239 variants that ClinVar lists as risk factors, associations or protective variants, with their genome-wide support.</span></a>
        <a class="tool" href="#help"><b>Help and programmatic access</b><span>A tutorial with worked examples, how to read every number, and Python and R clients for the data files.</span></a>
      </div>
      <div class="chips" aria-label="Examples"><span class="muted" style="font-size:13px;align-self:center">Try:</span><a href="#rs7903146">rs7903146</a><a href="#check=rs7903146|T2D">rs7903146 and T2D</a><a href="#check=rs1800795|type 2 diabetes">rs1800795 and type 2 diabetes</a><a href="#gene=MTHFR">MTHFR</a><a href="#trait=Rheumatoid arthritis">rheumatoid arthritis</a><a href="#paper=20967696">PMID 20967696</a></div>`;
  }

  // ---------------------------------------------------------------- check a claim
  async function showCheck(arg) {
    const [a, b] = (arg || '').split('|'); const bq = b ? decodeURIComponent(b) : '';
    out.innerHTML = `<h2>Check a claim</h2>
      <p class="muted">Enter a variant and the disease or trait it is claimed to affect. Part of the name, the full name or a common abbreviation (T2D, CAD, BMI) is enough.</p>
      <form id="cf2"><div class="fld" style="flex:0 1 180px"><label for="crs">rsID</label><input id="crs" placeholder="rs1800795" value="${esc(a || '')}"></div>
      <div class="fld"><label for="ctr">Disease or trait</label><input id="ctr" list="traitlist" placeholder="type 2 diabetes" value="${esc(bq)}"></div>
      <button type="submit">Check</button></form>
      <div class="chips"><span class="muted" style="font-size:13px;align-self:center">Examples:</span><a href="#check=rs7903146|type 2 diabetes">rs7903146, type 2 diabetes</a><a href="#check=rs1800795|type 2 diabetes">rs1800795, type 2 diabetes</a><a href="#check=rs1801133|CAD">rs1801133, CAD</a><a href="#check=rs738409|fatty liver disease">rs738409, fatty liver disease</a></div>
      <div id="cres"></div>`;
    $('#cf2').addEventListener('submit', (e) => { e.preventDefault(); const r = $('#crs').value.trim().toLowerCase().replace(/\s+/g, '').replace(/^(\d)/, 'rs$1'); const t = $('#ctr').value.trim(); if (r) location.hash = 'check=' + r + '|' + encodeURIComponent(t); });
    if (!a) return;
    const m = a.match(/^rs(\d+)$/i); const res = $('#cres');
    if (!m) { res.innerHTML = '<p>Enter an rsID such as rs1800795.</p>'; return; }
    const rs = +m[1]; const x = byRs.get(rs);
    if (!x) { res.innerHTML = notInAtlas(rs); return; }
    res.innerHTML = '<p class="muted">Checking…</p>';
    const cl = await claims(rs); const tq = bq.trim();
    const set = tq ? await resolveTraits(tq) : null;
    const M = set ? cl.filter((c) => set.has(c.trait)) : cl;
    if (!M.length) {
      res.innerHTML = `<div class="verdict nt"><b>No claim on rs${rs} names a trait matching “${esc(tq)}”.</b><span>Traits claimed for this variant:</span></div>
        <div class="list">${x.top_traits.map(([t]) => `<a href="#check=rs${rs}|${encodeURIComponent(t)}">${esc(cap(t))}</a>`).join('')}</div>${rsLinks(rs)}`;
      return;
    }
    const v = verdict(M); const traits = [...new Set(M.map((c) => cap(c.trait)))];
    const ns = v.R.map((c) => c.n).filter((n) => n != null).sort((p, q) => p - q); const med = ns.length ? ns[Math.floor(ns.length / 2)] : null;
    let kind, head, body;
    if (v.kind === 'none') { kind = 'nt'; head = 'Only reports of no association'; body = `The ${v.A} matching claims report no association, so there is no reported association to calibrate.`; }
    else if (v.kind === 'nt') { kind = 'nt'; head = 'Not testable against genome-wide evidence'; body = 'The GWAS Catalog has no association for this trait or a narrower term, so the claims cannot be calibrated.'; }
    else if (v.kind === 'sup') { kind = ''; head = `Supported genome-wide: ${v.S.length} of ${v.Tt.length} testable reported associations`; body = `A genome-wide significant association of this variant or a proxy with the trait was first reported in ${v.y}; ${v.E.length} of the supported claims came after that report.`; }
    else if (v.kind === 'part') { kind = ''; head = `Partly supported: ${v.S.length} of ${v.Tt.length} testable reported associations`; body = `Genome-wide support exists for some of the matched traits (first reported in ${v.y}).`; }
    else { kind = 'no'; head = `No genome-wide support: 0 of ${v.Tt.length} testable reported associations`; body = `For comparison, ${pct(x.null)} of frequency-matched random variants would be supported for this variant's claimed traits. Absence of support is not evidence of no effect.`; }
    res.innerHTML = `<div class="verdict ${kind}"><b>${esc(head)}</b><span>${esc(body)}</span>
        <span class="muted" style="font-size:13px">rs${rs}${x.gene ? ' (' + esc(x.gene) + ')' : ''}; matched traits: ${esc(traits.slice(0, 6).join('; '))}${traits.length > 6 ? '; …' : ''}. ${v.R.length} reported associations (median stated sample size ${num(med)}), ${v.A} reports of no association.${x.clinvar ? ' ClinVar: ' + esc(x.clinvar) + '.' : ''}</span></div>
      ${timeline(M, v.y)}
      <div class="row"><h3>Matching claims</h3><button class="sec" id="dl">Download CSV</button></div>
      <div class="tw"><table><thead>${claimHead}</thead><tbody>${M.map(claimRow).join('')}</tbody></table></div>
      <p class="muted" style="font-size:12.5px"><a href="#rs${rs}">All claims on rs${rs}</a> · <a href="#help">How the verdict is reached</a></p>${rsLinks(rs)}`;
    $('#dl').addEventListener('click', () => download(`rs${rs}_${tq.replace(/\W+/g, '_') || 'all'}_claims.csv`, claimCsv(M)));
  }

  // ---------------------------------------------------------------- screen a paper or a text
  const RS_RE = /\brs\s?(\d{2,10})\b/gi;
  // A phrase as a regular expression: words joined by spaces, hyphens or apostrophes; abbreviations only where written in capitals.
  function phraseRe(p) {
    const abbr = !p.includes(' ') && p.length <= 6 && /[a-z]/.test(p);
    const body = p.split(/[\s'-]+/).filter(Boolean).map(reEsc).join("[\\s'\\u2019-]+");
    // the leading group takes the character before the phrase (no look-behind, for older browsers); spans() strips it
    return abbr ? new RegExp('(^|[^A-Za-z0-9])(' + reEsc(p.toUpperCase()) + 's?)(?![A-Za-z0-9])') : new RegExp('(^|[^A-Za-z0-9])(' + body + ')(?![A-Za-z0-9])', 'i');
  }
  function spans(r, text) {
    return [...text.matchAll(new RegExp(r.source, r.flags + 'g'))].map((m) => [m.index + m[1].length, m.index + m[0].length]);
  }
  async function screenText(text) {
    const L = await lexicon();
    const ids = [...new Set([...text.matchAll(RS_RE)].map((m) => +m[1]))];
    const cand = [], CL = new Map();
    for (const rs of ids) {
      if (!byRs.has(rs)) continue;
      const cl = await claims(rs); CL.set(rs, cl);
      const groups = new Map(); cl.forEach((c) => { const k = c.trait.toLowerCase(); if (!groups.has(k)) groups.set(k, []); groups.get(k).push(c.trait); });
      for (const [k, names] of groups) {
        // only traits linked to the ontology can be calibrated; a phrase that is a gene symbol (IL6) does not name a trait
        const linked = names.some((t) => (L.byTrait.get(t) || {}).efo);
        const ph = [...new Set(names.flatMap((t) => (L.byTrait.get(t) || {}).phrases || [k]))].sort((p, q) => q.length - p.length);
        for (const p of ph) {
          if (!p.includes(' ') && byGene.has(p.toUpperCase().replace(/-/g, ''))) continue;
          const sp = spans(phraseRe(p), text);
          if (sp.length) { cand.push({ rs, k, t: names[0], phrase: p, sp, linked }); break; }
        }
      }
    }
    // a trait named only inside the name of another matched trait (diabetes mellitus in type 1 diabetes mellitus) is not counted
    const all = cand.flatMap((c) => c.sp);
    const inside = ([s, e]) => all.some(([a, b]) => a <= s && b >= e && b - a > e - s);
    const kept = cand.filter((c) => c.linked && c.sp.some((x) => !inside(x)));
    // sentences, to flag pairs named in the same sentence
    const ends = [...text.matchAll(/[.!?;](?=\s+[A-Z0-9(]|\s*$)/g)].map((m) => m.index);
    const sent = (i) => ends.filter((e) => e < i).length;
    const rsSent = new Map(); for (const m of text.matchAll(RS_RE)) { const r = +m[1]; if (!rsSent.has(r)) rsSent.set(r, new Set()); rsSent.get(r).add(sent(m.index)); }
    const rows = [], marks = [];
    for (const rs of ids) {
      const x = byRs.get(rs);
      if (!x) { rows.push({ rs, x: null }); continue; }
      const K = kept.filter((c) => c.rs === rs).map((c) => ({ ...c, near: c.sp.some(([a]) => rsSent.get(rs).has(sent(a))) })).sort((a, b) => b.near - a.near);
      K.forEach((c) => { marks.push(c.phrase); rows.push({ rs, x, t: c.t, phrase: c.phrase, near: c.near, v: verdict(CL.get(rs).filter((cl) => cl.trait.toLowerCase() === c.k)) }); });
      if (!K.length) rows.push({ rs, x, t: null });
    }
    return { ids, rows, marks };
  }
  function highlight(text, marks) {
    const sp = [];
    for (const m of text.matchAll(RS_RE)) sp.push([m.index, m.index + m[0].length, 'rs']);
    [...new Set(marks)].forEach((p) => spans(phraseRe(p), text).forEach(([s, e]) => sp.push([s, e, ''])));
    sp.sort((a, b) => a[0] - b[0] || b[1] - a[1]);
    let o = '', i = 0;
    for (const [s, e, c] of sp) { if (s < i) continue; o += esc(text.slice(i, s)) + `<mark${c ? ' class="rs"' : ''}>${esc(text.slice(s, e))}</mark>`; i = e; }
    return o + esc(text.slice(i));
  }
  const screenCsv = (rows) => [['rsid', 'gene', 'in_atlas', 'trait_named_in_text', 'phrase_found', 'same_sentence', 'verdict', 'reported_associations', 'testable', 'supported', 'first_genome_wide_year', 'variant_supported_share', 'variant_matched_null', 'clinvar']]
    .concat(rows.map((r) => ['rs' + r.rs, r.x && r.x.gene, r.x ? 'yes' : 'no', r.t || '', r.phrase || '', r.v ? (r.near ? 'yes' : 'no') : '', r.v ? r.v.label : r.x ? 'no trait of this variant named in the text' : '', r.v ? r.v.R.length : '', r.v ? r.v.Tt.length : '', r.v ? r.v.S.length : '', r.v ? r.v.y || '' : '', r.x ? r.x.supported : '', r.x ? r.x.null : '', r.x ? r.x.clinvar : '']));
  function screenTable(rows) {
    return `<div class="tw"><table><thead><tr><th>Variant</th><th>Gene</th><th>Trait named in the text</th><th>Same sentence</th><th>Verdict</th><th class="num">Supported of testable</th><th>First genome-wide report</th><th>Variant overall</th></tr></thead><tbody>${rows.map((r) => !r.x
      ? `<tr><td>rs${r.rs}</td><td colspan="7" class="nt">no claim in the atlas · ${ext(`https://www.ebi.ac.uk/gwas/variants/rs${r.rs}`, 'GWAS Catalog')} · ${ext(`https://www.ncbi.nlm.nih.gov/snp/rs${r.rs}`, 'dbSNP')}</td></tr>`
      : !r.t ? `<tr><td><a href="#rs${r.rs}">rs${r.rs}</a></td><td><i>${esc(r.x.gene)}</i></td><td class="nt" colspan="5">no trait claimed for this variant in the atlas is named in the text</td><td>${mini(r.x.supported)}${pct(r.x.supported)} supported${cls(r.x) ? ' (' + cls(r.x) + ')' : ''}</td></tr>`
      : `<tr><td><a href="#rs${r.rs}">rs${r.rs}</a></td><td><i>${esc(r.x.gene)}</i></td><td class="t"><a href="#check=rs${r.rs}|${encodeURIComponent(r.t)}">${esc(cap(r.t))}</a> <span class="muted">(“${esc(r.phrase)}”)</span></td><td>${r.near ? 'yes' : '—'}</td><td class="${vcls(r.v.kind)}">${r.v.label}</td><td class="num">${r.v.Tt.length ? `${r.v.S.length} of ${r.v.Tt.length}` : '—'}</td><td>${r.v.y || '—'}</td><td>${mini(r.x.supported)}${pct(r.x.supported)} supported</td></tr>`).join('')}</tbody></table></div>`;
  }
  async function renderScreen(box, text, inAtlas) {
    box.innerHTML = '<p class="muted">Screening…</p>';
    const { ids, rows, marks } = await screenText(text);
    if (!ids.length) { box.innerHTML = '<div class="verdict nt"><b>No rsID found in the text.</b><span>The screen looks for dbSNP identifiers such as rs1800795 or rs 1800795.</span></div>'; return; }
    const pairs = rows.filter((r) => r.v), sup = pairs.filter((r) => r.v.kind === 'sup' || r.v.kind === 'part').length, no = pairs.filter((r) => r.v.kind === 'no').length;
    box.innerHTML = `<div class="verdict ${no && !sup ? 'no' : pairs.length ? '' : 'nt'}"><b>${ids.length} variant${ids.length > 1 ? 's' : ''} named; ${pairs.length} variant–trait pair${pairs.length === 1 ? '' : 's'} matched to the atlas${pairs.length ? `: ${sup} with genome-wide support, ${no} without` : ''}.</b>
        <span>Each pair is a variant named in the text and a disease or trait that the text names and that published claims about this variant concern (pairs named in the same sentence first); the verdict summarizes those published claims${inAtlas ? '' : ' (the screened text itself is not part of the atlas)'}.</span></div>
      <div class="row"><h3>Variant–trait pairs</h3><button class="sec" id="sdl">Download CSV</button></div>${screenTable(rows)}
      <h3>Text with variants and matched traits marked</h3><p class="abs">${highlight(text, marks)}</p>`;
    $('#sdl').addEventListener('click', () => download('variant_claim_atlas_screen.csv', screenCsv(rows)));
  }
  async function showPaper(arg) {
    const pm = (arg || '').replace(/\D/g, '');
    out.innerHTML = `<h2>Screen a paper or a text</h2>
      <p class="muted">Enter a PubMed ID to screen its title and abstract, or paste any text (an abstract, a manuscript, a report). Every rsID in it is found, the diseases and traits it names are matched to the published claims about that variant, and each variant–trait pair receives the verdict of <a href="#check">Check a claim</a>.</p>
      <form id="pf"><div class="fld" style="flex:0 1 220px"><label for="pmid">PubMed ID</label><input id="pmid" inputmode="numeric" placeholder="20967696" value="${esc(pm)}"></div><button type="submit">Screen paper</button></form>
      <div class="chips"><span class="muted" style="font-size:13px;align-self:center">Examples:</span><a href="#paper=42749416">PMID 42749416 (2026, not in the atlas)</a><a href="#paper=20967696">PMID 20967696 (2010)</a><a href="#paper=32685554">PMID 32685554 (2020)</a></div>
      <details id="pdt"><summary style="cursor:pointer;color:var(--sup);font-size:14px">Or paste text (an abstract or a manuscript)</summary>
      <div style="display:flex;flex-direction:column;gap:8px;margin-top:8px"><textarea id="ptx" aria-label="Text to screen" placeholder="Paste an abstract or a manuscript. Text is screened in your browser and not sent anywhere."></textarea>
      <div style="display:flex;gap:8px;flex-wrap:wrap"><button id="prun" type="button">Screen text</button><button class="sec" id="pex" type="button">Load example text</button></div></div></details>
      <div id="pres" style="display:flex;flex-direction:column;gap:12px"></div>`;
    $('#pf').addEventListener('submit', (e) => { e.preventDefault(); const v = $('#pmid').value.replace(/\D/g, ''); if (v) location.hash = 'paper=' + v; });
    $('#pex').addEventListener('click', () => { $('#ptx').value = 'We genotyped TCF7L2 rs7903146 and IL6 rs1800795 in 812 patients with type 2 diabetes (T2D) and 790 controls. The T allele of rs7903146 was associated with T2D (OR 1.41). The rs1800795 C allele was associated with higher T2D risk and with obesity. MTHFR rs1801133 was not associated with coronary artery disease.'; });
    $('#prun').addEventListener('click', () => { const t = $('#ptx').value.trim(); if (t) renderScreen($('#pres'), t, false); });
    if (!pm) return;
    const box = $('#pres'); box.innerHTML = '<p class="muted">Retrieving the abstract from Europe PMC…</p>';
    const P = await papers(); const rec = P.get(+pm);
    let title = '', text = '', meta = '';
    try {
      const r = await fetch(`https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:${pm}%20AND%20SRC:MED&resultType=core&format=json`).then((r) => r.json());
      const a = r.resultList && r.resultList.result && r.resultList.result[0];
      if (a) {
        const plain = (h) => { const d = document.createElement('div'); d.innerHTML = h || ''; return d.textContent; };
        title = plain(a.title); text = plain(a.abstractText);
        meta = [a.authorString && a.authorString.split(',').slice(0, 3).join(',') + (a.authorString.split(',').length > 3 ? ' et al.' : ''), a.journalInfo && a.journalInfo.journal && a.journalInfo.journal.title, a.pubYear].filter(Boolean).join(' · ');
      }
    } catch (e) { /* offline or blocked: fall back to the atlas records */ }
    let atlas = '', acl = [];
    if (rec) {
      for (const rs of rec.rsids) (await claims(rs)).forEach((c) => { if (c.pmid === +pm) acl.push(c); });
      atlas = `<div class="row"><h3>Claims read from this paper in the atlas</h3><button class="sec" id="adl">Download CSV</button></div>
        <div class="tw"><table><thead><tr><th>Variant</th>${claimHead.slice(4)}</thead><tbody>${acl.map((c) => `<tr><td><a href="#rs${c.rsid}">rs${c.rsid}</a></td>` + claimRow(c).slice(4)).join('')}</tbody></table></div>
        <p class="muted" style="font-size:12.5px">Genome-wide support in this table is for each claim as read from the paper; “earlier” means a genome-wide report preceded the paper.</p>`;
    }
    const head = `<div class="verdict nt" style="border-left-color:var(--sup)"><b>${title ? esc(title) : 'PMID ' + esc(pm)}</b><span class="muted" style="font-size:13px">${meta ? esc(meta) + ' · ' : ''}${ext(`https://pubmed.ncbi.nlm.nih.gov/${pm}/`, 'PubMed ' + esc(pm))}${rec ? ` · in the atlas (${rec.rsids.length} variant${rec.rsids.length > 1 ? 's' : ''})` : ' · not in the atlas'}</span></div>`;
    if (!text && !rec) { box.innerHTML = head + '<p>No abstract could be retrieved for this PubMed ID (Europe PMC). Paste the text instead.</p>'; return; }
    box.innerHTML = head + atlas + '<h3>Screen of the title and abstract</h3><div id="ps" style="display:flex;flex-direction:column;gap:12px"></div>';
    if (rec) $('#adl').addEventListener('click', () => download(`PMID${pm}_claims.csv`, claimCsv(acl)));
    if (text) await renderScreen($('#ps'), (title ? title + ' ' : '') + text, !!rec);
  }

  // ---------------------------------------------------------------- batch
  const batchCsv = (L) => [['query', 'rsid', 'gene', 'trait_query', 'matched_traits', 'verdict', 'papers', 'reported_associations', 'reports_of_no_association', 'testable', 'supported', 'supported_share', 'matched_null', 'class', 'median_stated_n', 'first_genome_wide_year', 'clinvar', 'in_atlas']]
    .concat(L.map(({ q, x, t, r }) => r ? [q, 'rs' + r.rs, x && x.gene, t, r.traits, r.label, '', r.R, r.A, r.T, r.S, r.T ? (r.S / r.T).toFixed(3) : '', x && x.null, '', '', r.y || '', x && x.clinvar, 'yes']
      : x ? [q, 'rs' + x.rsid, x.gene, '', '', '', x.papers, x.reported, x.absent, x.testable, '', x.supported, x.null, cls(x), x.median_n, x.first_gw_year, x.clinvar, 'yes'] : [q, '', '', t || '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'no']));
  function showBatch() {
    out.innerHTML = `<h2>Batch check</h2>
      <p class="muted">Paste rsIDs (any separator), or one variant–trait pair per line (<code>rs1800795, type 2 diabetes</code>; comma, semicolon, colon or tab), or upload a text or CSV file. Up to 10,000 lines; the list is read in your browser and not sent anywhere. <a href="examples/batch_example.csv">Example file</a>.</p>
      <textarea id="bt" aria-label="Variants or variant–trait pairs" placeholder="rs9939609&#10;rs1800795&#10;rs738409, fatty liver disease"></textarea>
      <div class="row"><div style="display:flex;gap:8px;flex-wrap:wrap;align-items:center"><button id="brun">Check</button><button class="sec" id="bex">Load example</button><label style="display:inline-flex;align-items:center;gap:6px;cursor:pointer">Upload <input type="file" id="bfile" accept=".txt,.csv,.tsv" style="width:auto;border:0;padding:0"></label></div><button class="sec" id="bdl" disabled>Download CSV</button></div>
      <div id="bres"></div>`;
    let last = [];
    $('#bex').addEventListener('click', () => { $('#bt').value = 'rs9939609\nrs1800795\nrs738409, fatty liver disease\nrs1801133, CAD\nrs7903146, T2D\nrs2910164\nrs671, alcohol'; });
    $('#bfile').addEventListener('change', async (e) => { const f = e.target.files[0]; if (f) $('#bt').value = await f.text(); });
    $('#bdl').addEventListener('click', () => download('variant_claim_atlas_batch.csv', batchCsv(last)));
    $('#brun').addEventListener('click', async () => {
      const lines = $('#bt').value.split(/\r?\n/).map((l) => l.trim()).filter(Boolean).slice(0, 10000); const L = [];
      for (const l of lines) {
        const pair = l.match(/^(rs\s?\d+)\s*[,;:\t]\s*(\S.*)$/i);
        if (pair) L.push({ q: l, rs: +pair[1].replace(/\D/g, ''), t: pair[2].replace(/^["']|["']$/g, '').trim() });
        else for (const m of l.matchAll(RS_RE)) L.push({ q: 'rs' + m[1], rs: +m[1] });
      }
      const res = $('#bres'); res.innerHTML = `<p class="muted">Checking ${L.length} entries…</p>`;
      last = []; const tcache = new Map();
      for (let i = 0; i < L.length; i++) {
        const e = L[i], x = byRs.get(e.rs);
        if (e.t && x) {
          if (!tcache.has(e.t.toLowerCase())) tcache.set(e.t.toLowerCase(), await resolveTraits(e.t));
          const set = tcache.get(e.t.toLowerCase());
          const cl = (await claims(e.rs)).filter((c) => set.has(c.trait)); const v = verdict(cl);
          last.push({ q: e.q, x, t: e.t, r: { rs: e.rs, R: v.R.length, A: v.A, T: v.Tt.length, S: v.S.length, y: v.y, label: cl.length ? v.label : 'no claim on this trait', kind: cl.length ? v.kind : 'none', traits: [...new Set(cl.map((c) => c.trait))].join('; ') } });
        } else last.push({ q: e.q, x, t: e.t });
        if (i % 250 === 249) { res.innerHTML = `<p class="muted">Checked ${i + 1} of ${L.length} entries…</p>`; await new Promise((r) => setTimeout(r)); }
      }
      const found = last.filter((r) => r.x).length;
      res.innerHTML = `<p>${found} of ${last.length} entries are in the atlas.${last.length > 2000 ? ' The first 2,000 rows are shown; the CSV holds all.' : ''}</p><div class="tw"><table><thead><tr><th>Query</th><th>Gene</th><th class="num">Reported</th><th class="num">Testable</th><th class="num">Supported</th><th class="num">Matched null</th><th>Class or verdict</th><th>ClinVar</th></tr></thead><tbody>${
        last.slice(0, 2000).map(({ q, x, r }) => !x ? `<tr><td>${esc(q)}</td><td colspan="7" class="nt">not in the atlas</td></tr>` : r
          ? `<tr><td><a href="#check=rs${r.rs}|${encodeURIComponent(q.split(/[,;:\t]/).slice(1).join(',').trim())}">${esc(q)}</a></td><td><i>${esc(x.gene)}</i></td><td class="num">${r.R}</td><td class="num">${r.T}</td><td class="num">${r.T ? mini(r.S / r.T) + pct(r.S / r.T) : '—'}</td><td class="num">${pct(x.null)}</td><td class="${vcls(r.kind)}">${r.label}</td><td>${esc(x.clinvar)}</td></tr>`
          : `<tr><td><a href="#rs${x.rsid}">${esc(q)}</a></td><td><i>${esc(x.gene)}</i></td><td class="num">${num(x.reported)}</td><td class="num">${num(x.testable)}</td><td class="num">${x.supported == null ? '—' : mini(x.supported) + pct(x.supported)}</td><td class="num">${pct(x.null)}</td><td>${cls(x) || 'fewer than 5 reported'}</td><td>${esc(x.clinvar)}</td></tr>`).join('')}</tbody></table></div>`;
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
        L.map((r) => `<tr><td><a href="#${r.rsid}">${esc(r.rsid)}</a></td><td><i>${esc(r.gene_clinvar)}</i></td><td class="t" title="${esc(r.clinvar_conditions)}">${esc(r.clinvar_significance)}</td><td class="num">${esc(r.reported_associations)}</td><td class="num">${mini(+r.share_supported)}${pct(+r.share_supported)}</td><td class="num">${pct(+r.matched_null)}</td><td>${esc(r.support_class)}</td><td>${r.clinvar_condition_in_ontology === 'yes' ? (r.genome_wide_support_for_clinvar_condition === 'yes' ? '<span class="yes">yes</span>' : '<span class="no">no</span>') : '<span class="nt">condition not mapped</span>'}</td></tr>`).join('')}</tbody></table><p class="muted" style="font-size:12.5px">${L.length} of ${CV.length} variants. Click a column heading to sort; hover over an assertion for the ClinVar conditions.</p>`;
      document.querySelectorAll('#cvt th[data-k]').forEach((th) => th.addEventListener('click', () => { if (key === th.dataset.k) dir = -dir; else { key = th.dataset.k; dir = -1; } draw(); }));
    };
    $('#cvq').addEventListener('input', draw); $('#cvc').addEventListener('change', draw); draw();
    $('#dl').addEventListener('click', () => { const h = Object.keys(CV[0]); download('clinvar_risk_variants.csv', [h].concat(CV.map((r) => h.map((k) => r[k])))); });
  }

  // ---------------------------------------------------------------- help
  function showHelp() {
    const B = 'https://zxu69461-creator.github.io/variant-claim-atlas/';
    out.innerHTML = `<div class="help" style="display:flex;flex-direction:column;gap:12px">
      <h2>Help</h2>
      <p>The tool answers one question at the point of use: <i>how often is a published claim like this one backed by genome-wide evidence?</i> It holds the disease and trait claims made about named variants (rsIDs) in PubMed titles and abstracts, and sets each claim against genome-wide significant associations in the NHGRI-EBI GWAS Catalog, with frequency-matched random variants as the chance expectation. All views run in your browser and return their results at once.</p>
      <h3>Quick start</h3>
      <ol>
        <li><b>Look up</b> a variant, a gene or a trait in the search box (<a href="#rs7903146">rs7903146</a>, <a href="#gene=IL6">IL6</a>, <a href="#trait=Asthma">asthma</a>). The variant page shows the support of its claims against the matched null, by stated sample size and by year, the most-claimed traits and every claim with its PubMed link.</li>
        <li><b>Check a claim</b>: enter an rsID and a trait (<a href="#check=rs7903146|T2D">rs7903146 and T2D</a>). The verdict states how many testable reported associations are supported genome-wide and since when.</li>
        <li><b>Screen a paper</b>: enter a PubMed ID (<a href="#paper=42749416">42749416</a>, a 2026 paper that is not in the atlas) or paste a text. Each variant–trait pair named in it is checked.</li>
        <li><b>Batch check</b> a list of up to 10,000 variants or variant–trait pairs (<a href="examples/batch_example.csv">example file</a>) and download the table.</li>
        <li><b>ClinVar risk variants</b>: filter the 239 variants that ClinVar lists as risk factors, associations or protective variants by support class, gene or condition.</li>
      </ol>
      <h3>Accepted input</h3>
      <dl>
        <dt>Variants</dt><dd>dbSNP identifiers: <code>rs1800795</code>, <code>RS1800795</code>, <code>rs 1800795</code> or the number alone in Check a claim.</dd>
        <dt>Diseases and traits</dt><dd>Part of the ontology name (<code>diabetes</code>), the full name, or a phrase that abstracts used for it, including abbreviations (<code>T2D</code>, <code>CAD</code>, <code>BMI</code>). Matching is case-insensitive.</dd>
        <dt>Batch lists</dt><dd>rsIDs separated by anything, or one <code>rsID, trait</code> pair per line (comma, semicolon, colon or tab); text, CSV or TSV files.</dd>
        <dt>Papers and texts</dt><dd>A PubMed ID (title and abstract are retrieved from Europe PMC), or any pasted text. Abbreviations in a text are matched only where written in capitals.</dd>
      </dl>
      <h3>How to read the results</h3>
      <dl>
        <dt>Reported association</dt><dd>A claim that the variant is associated with the disease or trait; a <i>report of no association</i> states its absence.</dd>
        <dt>Testable</dt><dd>The trait, or a narrower ontology term, has at least one association in the GWAS Catalog, so the claim can be compared with genome-wide evidence.</dd>
        <dt>Supported</dt><dd>The variant, or a proxy in linkage disequilibrium (r² ≥ 0.8 in the 1000 Genomes European or East Asian reference), has an association at P &lt; 5 × 10⁻⁸ with the same trait or a broader or narrower ontology term, in a GWAS Catalog paper other than the claim's own. Support does not require the same allele or direction of effect.</dd>
        <dt>Earlier</dt><dd>That genome-wide report was published before the claim, so the claim restated a known association.</dd>
        <dt>Matched null</dt><dd>The share of frequency-matched random variants that would be supported for the same traits: the support expected by chance. A support well above the matched null indicates a real genome-wide signal.</dd>
        <dt>Variant class</dt><dd>For variants with at least five reported associations: <i>mostly supported</i> (≥ 50% of testable ones), <i>partly supported</i> (5–50%), <i>rarely supported</i> (&lt; 5%).</dd>
        <dt>Verdict for a claim</dt><dd><i>Supported</i> when at least half of the testable reported associations on that variant and trait are supported; <i>partly supported</i> when some are; <i>not supported</i> when none is; <i>not testable</i> when no matching claim is testable.</dd>
      </dl>
      <h3>Worked examples (live output)</h3>
      <ul>
        <li><a href="#check=rs7903146|type 2 diabetes">rs7903146 (<i>TCF7L2</i>) and type 2 diabetes</a>: supported; the genome-wide association was first reported in 2007 and nearly all later claims restate it.</li>
        <li><a href="#check=rs1800795|type 2 diabetes">rs1800795 (<i>IL6</i>) and type 2 diabetes</a>: none of the testable reported associations is supported genome-wide, although the variant has been claimed in hundreds of papers.</li>
        <li><a href="#check=rs738409|fatty liver disease">rs738409 (<i>PNPLA3</i>) and fatty liver disease</a>: supported, with every testable reported association backed genome-wide.</li>
        <li><a href="#clinvar">ClinVar risk variants</a> filtered to <i>rarely (&lt;5%)</i>: risk variants whose published associations are almost never supported genome-wide.</li>
        <li><a href="#paper=42749416">A 2026 paper</a> screened from its abstract, and <a href="#paper=20967696">a 2010 paper</a> whose claims are in the atlas.</li>
      </ul>
      <h3>Programmatic access</h3>
      <p>The data behind the tool are static files under <code>${B}</code>, free to fetch without a key. Each JSON file holds <code>cols</code> (column names) and <code>rows</code> (one array per row).</p>
      <div class="tw"><table><thead><tr><th>File</th><th>Content</th></tr></thead><tbody>
        <tr><td><code>data/meta.json</code></td><td class="t">releases of the sources and counts</td></tr>
        <tr><td><code>data/variants.json</code></td><td class="t">one row per variant: gene, papers, reported, absent, testable, supported share, matched null, median stated n, support by sample-size class, top traits, ClinVar assertion, first genome-wide year</td></tr>
        <tr><td><code>data/traits.json</code></td><td class="t">one row per trait: area, counts, supported share, matched null, most-claimed variants</td></tr>
        <tr><td><code>data/claims/NN.json</code></td><td class="t">every claim record; <code>NN</code> = rsID number modulo 64, two digits (rs7903146 → <code>data/claims/50.json</code>)</td></tr>
        <tr><td><code>data/lexicon.json</code></td><td class="t">per trait: ontology identifier and the phrases abstracts used for it</td></tr>
        <tr><td><code>data/papers.json</code></td><td class="t">per PubMed ID: year and the variants with a claim</td></tr>
        <tr><td><code>atlas/Data_S1_variant_trait_claim_atlas.tsv.gz</code></td><td class="t">the full atlas (tab-separated)</td></tr>
        <tr><td><code>atlas/Data_S3_clinvar_risk_variants.tsv</code></td><td class="t">the ClinVar table</td></tr>
      </tbody></table></div>
      <p>Ready-to-run clients reproduce the views and the verdict: <a href="examples/vca_client.py">Python</a> (standard library only) and <a href="examples/vca_client.R">R</a> (package jsonlite). For example:</p>
      <pre>python vca_client.py check rs7903146 "type 2 diabetes"
python vca_client.py batch my_variants.txt -o results.csv
Rscript vca_client.R check rs1800795 "type 2 diabetes"</pre>
      <h3>Frequently asked questions</h3>
      <dl>
        <dt>My variant is not in the atlas</dt><dd>The atlas holds variants named by rsID in a PubMed title or abstract with a disease or trait claim (PubMed 2026 baseline). Papers naming a variant only in the full text, by HGVS name or by a legacy name (such as −174G&gt;C) are not included; the screen of a paper or text still finds rsIDs written in it.</dd>
        <dt>Does “not supported” mean the claim is false?</dt><dd>No. It means that no genome-wide significant association of the variant or a proxy with the trait has been reported. Effects too small, too rare or too population-specific for genome-wide studies can remain undetected.</dd>
        <dt>Does “supported” mean the variant is causal?</dt><dd>No. Support means that the locus is associated genome-wide with the trait; the claimed variant may tag a causal variant nearby.</dd>
        <dt>How were claims read?</dt><dd>By an open-weight language model run locally, which quotes the sentence of each claim; a blinded reading of held-out papers found 92% precision for reported associations. See the paper and the analysis code.</dd>
        <dt>How does this differ from LitVar 2.0, the GWAS Catalog or ClinVar?</dt><dd>${ext('https://www.ncbi.nlm.nih.gov/research/litvar2/', 'LitVar 2.0')} finds the papers that mention a variant; the ${ext('https://www.ebi.ac.uk/gwas/', 'GWAS Catalog')} lists genome-wide associations; ${ext('https://www.ncbi.nlm.nih.gov/clinvar/', 'ClinVar')} lists clinical assertions. This tool reads what each paper claimed and sets every claim against genome-wide evidence and a matched null, so the reliability of a claim can be judged; every variant page links to these resources.</dd>
      </dl>
      <h3>Privacy, browsers and maintenance</h3>
      <ul>
        <li>No registration, no cookies, no analytics or tracking. Lists, files and texts that you paste or upload are processed in your browser and never sent anywhere. Screening a paper by PubMed ID sends only that ID to Europe PMC. The host, GitHub Pages, keeps standard server logs.</li>
        <li>Tested in current Chrome, Safari and Firefox on desktop and on phone-sized screens.</li>
        <li>The authors maintain the tool for at least five years. Each data update is a new numbered version archived at Zenodo (${ext('https://doi.org/10.5281/zenodo.24314560', 'doi:10.5281/zenodo.24314560')}), where every earlier version stays available. Report problems at ${ext('https://github.com/zxu69461-creator/variant-claim-atlas/issues', 'GitHub issues')}.</li>
      </ul>
      <h3>Version history</h3>
      <ul>
        <li><b>1.2.0–1.2.2</b> (21 November 2026): screening of a paper or a text; trait abbreviations and phrases in Check a claim and Batch check; claims per year on variant pages; links to dbSNP, the GWAS Catalog, ClinVar, LitVar 2.0, Ensembl, Open Targets and the ontology; this help page; Python and R clients.</li>
        <li><b>1.1.0–1.1.1</b> (21 November 2026): Check a claim, Batch check with CSV export, ClinVar risk variants, filters and downloads.</li>
        <li><b>1.0.0–1.0.2</b> (17 November 2026): look-up of variants, genes and traits; atlas data and analysis code.</li>
      </ul>
    </div>`;
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
        <li>Clients for programmatic access: <a href="examples/vca_client.py">Python</a>, <a href="examples/vca_client.R">R</a>; data files listed under <a href="#help">Help</a></li>
        <li>Archived release with all files: <a href="https://doi.org/10.5281/zenodo.24314560">doi:10.5281/zenodo.24314560</a></li>
      </ul>
      <h3>Sources</h3>
      <p>PubMed (${esc(META.pubmed)}) and PubTator3; NHGRI-EBI GWAS Catalog (release of ${esc(META.gwas_catalog)}); 1000 Genomes phase 3; Experimental Factor Ontology ${esc(META.efo)}, Mondo and MeSH; ClinVar (release of ${esc(META.clinvar)}); dbSNP gene names.</p>
      <h3>Cite</h3>
      <p>Hu Y, Xu Z, Jiang S, Huang L, Pan M, Luo S, Zhou H, Guo X. Variant–Trait Claim Atlas: checking 201,720 published variant–trait claims against genome-wide evidence and ClinVar. Manuscript under review, 2026. Data and code: doi:10.5281/zenodo.24314560.</p>
      <h3>Licence</h3>
      <p>Code under the MIT licence; atlas and tool data under CC BY 4.0. Free for any use with attribution. The atlas records what each abstract claimed; the abstracts themselves are not redistributed.</p>`;
  }

  // ---------------------------------------------------------------- routing
  function tab(t) {
    document.querySelectorAll('nav.tabs a').forEach((a) => { if (a.dataset.tab === t) a.setAttribute('aria-current', 'page'); else a.removeAttribute('aria-current'); });
    $('#searchcard').hidden = t !== 'lookup';
  }
  function route() {
    const h = decodeURIComponent(location.hash.slice(1)); let m;
    window.scrollTo(0, 0);
    if (h === 'check' || h.startsWith('check=')) { tab('check'); return showCheck(h.slice(6)); }
    if (h === 'paper' || h.startsWith('paper=')) { tab('paper'); return showPaper(h.slice(6)); }
    if (h === 'batch') { tab('batch'); return showBatch(); }
    if (h === 'clinvar') { tab('clinvar'); return showClinvar(); }
    if (h === 'help') { tab('help'); return showHelp(); }
    if (h === 'about') { tab('about'); return showAbout(); }
    tab('lookup');
    if (!h || h === 'lookup') return home();
    if ((m = h.match(/^rs(\d+)$/i))) return showVariant(+m[1]);
    if ((m = h.match(/^gene=(.+)$/))) return showGene(m[1]);
    if ((m = h.match(/^trait=(.+)$/))) return showTrait(m[1]);
    if ((m = h.match(/^q=(.+)$/))) return interpret(m[1]);
    return home();
  }
  async function interpret(q) {
    q = q.trim();
    if (/^rs\s?\d+$/i.test(q)) { location.hash = q.toLowerCase().replace(/\s/g, ''); return; }
    if (byGene.has(q.toUpperCase())) { location.hash = 'gene=' + q.toUpperCase(); return; }
    // a trait whose name equals the query is opened directly only when it is also the most-claimed match (T2D opens the list, led by
    // type 2 diabetes mellitus, not the small unlinked trait called 'T2D'); a single match is opened directly
    const set = await resolveTraits(q); const t = T.find((r) => r.trait.toLowerCase() === q.toLowerCase());
    if (t) set.add(t.trait);
    const best = T.find((r) => set.has(r.trait));
    if (best && (set.size === 1 || best === t)) { location.hash = 'trait=' + encodeURIComponent(best.trait); return; }
    if (location.hash !== '#q=' + encodeURIComponent(q)) history.replaceState(null, '', '#q=' + encodeURIComponent(q));
    showSearch(q);
  }
  const safe = (f) => (...a) => Promise.resolve(f(...a)).catch((e) => { out.innerHTML = `<p>Something went wrong: ${esc(e.message)}. Reload the page; if it persists, please report it at <a href="https://github.com/zxu69461-creator/variant-claim-atlas/issues">GitHub issues</a>.</p>`; });
  $('#f').addEventListener('submit', (e) => { e.preventDefault(); const q = $('#q').value; if (q.trim()) safe(interpret)(q); });
  window.addEventListener('hashchange', safe(route));
  load().catch((e) => { out.innerHTML = '<p>The atlas data could not be loaded: ' + esc(e.message) + '. Check the connection and reload.</p>'; });
})();

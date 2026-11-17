"""Claim-trait mapping (charter amendment 2): PubTator3 disease id (MeSH or OMIM) -> EFO term, by route efo_xref, mondo_xref, then the
nearest mapped MeSH ancestor (tree level >= 2). Ontology files only; the GWAS Catalog is read only to report testability afterwards."""
import sys, gzip, collections
import xml.etree.ElementTree as ET
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import pandas as pd, duckdb
from rdata.proj import Project
P = Project(__file__); R = P.root / 'raw'
MESH = '/Volumes/ResearchData/Claude/store/raw/NLM_MeSH_2026/desc2026.gz'

def mondo_xrefs(obo):
    """MONDO term -> list of (xref id, is_equivalence); obsolete terms skipped."""
    out, cur, obs = collections.defaultdict(list), None, False
    for line in open(obo):
        line = line.rstrip('\n')
        if line.startswith('['): cur, obs = None, False; continue
        if line.startswith('id: MONDO:'): cur = line[4:].replace(':', '_')
        elif cur and line == 'is_obsolete: true': obs = True; out.pop(cur, None)
        elif cur and not obs and line.startswith('xref: '):
            x = line[6:].split(' ')[0]
            out[cur].append((x, 'MONDO:equivalentTo' in line))
    return out

def mesh_trees(path):
    ui2tn, tn2ui = {}, {}
    for _, el in ET.iterparse(gzip.open(path), events=('end',)):
        if el.tag != 'DescriptorRecord': continue
        ui = el.findtext('DescriptorUI'); tns = [t.text for t in el.iter('TreeNumber')]
        ui2tn[ui] = tns
        for t in tns: tn2ui[t] = ui
        el.clear()
    return ui2tn, tn2ui

with P.run(__file__, params={'routes': 'efo_xref > mondo_xref > mesh_ancestor (level >= 2)', 'mondo_xref': 'equivalence marks preferred'}) as run:
    run.input('MONDO_202611', path=R / 'mondo.obo', url='https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.obo')
    run.input('NLM_MeSH_2026', path=MESH)
    run.input('pubtator_disease_ids', path=P.root / 'work' / 'disease_ids.tsv')
    con = duckdb.connect(str(P.root / 'db' / 'varclaims.duckdb'), read_only=True)
    efo_terms = set(con.sql('select efo_id from efo_terms').df().efo_id)
    efo_mesh = con.sql('select * from efo_mesh').df()
    gw = con.sql("""select distinct regexp_extract(efo_uri, '[^/]+$') as efo_id from gwas_catalog_rsid_trait""").df()
    anc = con.sql('select * from efo_ancestors').df(); con.close()

    direct = collections.defaultdict(set)
    for e, m in zip(efo_mesh.efo_id, efo_mesh.mesh_id): direct[m].add(e)
    mx = mondo_xrefs(R / 'mondo.obo')
    via_mondo = collections.defaultdict(lambda: [set(), set()])  # source id -> [equivalence targets, any targets]
    for mo, xs in mx.items():
        efo_of = {mo} if mo in efo_terms else {x.replace(':', '_') for x, _ in xs if x.startswith('EFO:')} & efo_terms
        if not efo_of: continue
        for x, eq in xs:
            if x.startswith(('MESH:', 'OMIM:')):
                via_mondo[x][0 if eq else 1] |= efo_of
    def map12(i):
        if direct.get(i): return direct[i], 'efo_xref'
        eq, anyx = via_mondo.get(i, (set(), set()))
        if eq: return eq, 'mondo_xref'
        if anyx: return anyx, 'mondo_xref'
        return set(), None

    ui2tn, tn2ui = mesh_trees(MESH)
    ids = pd.read_csv(P.root / 'work' / 'disease_ids.tsv', sep='\t')
    rows = []
    for i, n in zip(ids.d, ids.n):
        terms, route = map12(i); via = None
        if not terms and i.startswith('MESH:'):
            tns = ui2tn.get(i[5:], [])
            for up in range(1, 20):
                lvl = {'.'.join(t.split('.')[:-up]) for t in tns if t.count('.') >= up}
                lvl = {t for t in lvl if '.' in t}  # first tree level (e.g. C14) not used
                if not lvl: break
                hit = {}
                for t in lvl:
                    a = 'MESH:' + tn2ui[t]; s, _ = map12(a)
                    if s: hit[a] = s
                if hit:
                    terms = set().union(*hit.values()); route = 'mesh_ancestor'; via = '|'.join(sorted(hit)); break
        for e in sorted(terms) or [None]:
            rows.append((i, n, e, route, via))
    m = pd.DataFrame(rows, columns=['pubtator_id', 'n_mentions', 'efo_id', 'route', 'via_mesh'])

    # testability (report only): the term or a descendant has a GWAS Catalog association
    gws = set(gw.efo_id.dropna()) - {""}; testable = set(gws) | set(anc[anc.efo_id.isin(gws)].ancestor)
    m['testable'] = m.efo_id.isin(testable)
    run.table('claim_trait_map', m, note='PubTator3 disease id -> EFO term by charter amendment 2 routes (efo_xref, mondo_xref, mesh_ancestor)')
    per = m.groupby('pubtator_id').agg(n=('n_mentions', 'first'), route=('route', 'first'), testable=('testable', 'max'))
    per['route'] = per.route.fillna('unmapped')
    s = per.groupby('route').agg(ids=('n', 'size'), mentions=('n', 'sum'))
    s['mention_share'] = s.mentions / per.n.sum()
    print(s.to_string(float_format=lambda x: f'{x:.3f}'))
    print('testable share of mentions:', f'{per.n[per.testable].sum() / per.n.sum():.3f}')
    print('largest unmapped:'); print(per[per.route == 'unmapped'].sort_values('n', ascending=False).head(12).to_string())

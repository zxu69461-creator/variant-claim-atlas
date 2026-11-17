"""EFO ontology tables for trait matching (charter amendment 1): terms, is_a edges, transitive ancestors, and MeSH cross-references.
EFO ids are written as in GWAS Catalog URIs' last path segment (e.g. EFO_0000400, MONDO_0005148)."""
import sys, collections
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import pandas as pd
from rdata.proj import Project
P = Project(__file__); R = P.root / 'raw'
with P.run(__file__, params={'file': 'efo.obo (EBISPOT latest, 2026-11-04)'}) as run:
    run.input('EFO_202611', path=R / 'efo.obo', url='https://github.com/EBISPOT/efo/releases/latest/download/efo.obo')
    terms, isa, mesh, cur, obs = {}, [], [], None, False
    norm = lambda x: x[4:] if x.startswith('efo:') else x.replace(':', '_')  # EFO's own terms are written efo:EFO_0000400
    for line in open(R / 'efo.obo'):
        line = line.rstrip('\n')
        if line.startswith('['):
            cur = None; obs = False
            continue
        if line.startswith('id: '): cur = norm(line[4:])
        elif cur is None: continue
        elif line.startswith('name: '): terms[cur] = line[6:]
        elif line == 'is_obsolete: true': obs = True; terms.pop(cur, None)
        elif line.startswith('is_a: ') and not obs: isa.append((cur, norm(line[6:].split(' ')[0])))
        elif line.startswith('xref: MESH:') and not obs: mesh.append((cur, line[6:].split(' ')[0]))
    par = collections.defaultdict(set)
    for c, p in isa: par[c].add(p)
    memo = {}
    def anc(t):
        if t in memo: return memo[t]
        memo[t] = set()
        s = set()
        for p in par.get(t, ()): s |= {p} | anc(p)
        memo[t] = s; return s
    anc_rows = [(t, a) for t in terms for a in anc(t)]
    run.table('efo_terms', pd.DataFrame(list(terms.items()), columns=['efo_id', 'name']), note='EFO terms (non-obsolete)')
    run.table('efo_isa', pd.DataFrame(isa, columns=['child', 'parent']), note='EFO is_a edges')
    run.table('efo_ancestors', pd.DataFrame(anc_rows, columns=['efo_id', 'ancestor']), note='EFO transitive is_a ancestors (excluding self)')
    run.table('efo_mesh', pd.DataFrame(mesh, columns=['efo_id', 'mesh_id']).drop_duplicates(), note='EFO -> MeSH cross-references (MESH:Dxxxxxx)')
    print(len(terms), 'terms;', len(isa), 'is_a;', len(anc_rows), 'ancestor rows;', len(mesh), 'MeSH xrefs')

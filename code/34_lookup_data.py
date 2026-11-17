"""Data files for the public lookup page of the atlas (2026-11-17). Every disease or trait claim of both splits that reports an
association or its absence, with its genome-wide comparison where the claim is testable. Output lookup_site/data/:
  variants.json  one row per variant: summary of its claims (support, matched null, sample-size classes, top traits, ClinVar)
  traits.json    one row per linked trait: summary and the most-claimed variants
  claims/NN.json the claim records, sharded by rsID modulo 64
  meta.json      releases and counts
Numbers are descriptive (both splits together); the confirmatory estimates of the paper use held-out papers only."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent)); import figdata as F
OUT = F.ROOT / 'lookup_site' / 'data'; (OUT / 'claims').mkdir(parents=True, exist_ok=True)
NSHARD = 64
SIZE = [0, 500, 2000, 10000, np.inf]

with F.con() as c:
    d = c.sql("""select v.pmid, v.rsid, v.efo, v.outcome, v.outcome_type, v.status, v.pub_year, v.n, v.design, v.pop_group, v.dev,
                 coalesce(m1.testable, m2.testable) testable, coalesce(m1.replicated, m2.replicated) replicated,
                 coalesce(m1.rep_prior, m2.rep_prior) rep_prior, coalesce(m1.first_hit_year, m2.first_hit_year) first_hit_year,
                 coalesce(m1.null_share, m2.null_share) null_share
                 from variant_claims v
                 left join variant_claim_match m1 using (pmid, rsid, trait_key, pop_group)
                 left join variant_claim_match_heldout m2 using (pmid, rsid, trait_key, pop_group)
                 where v.outcome_type <> 'drug_response' and v.status in ('positive', 'null')""").df()
d['testable'] = d.testable.fillna(False).astype(bool)
d['replicated'] = d.replicated.fillna(False).astype(bool) & d.testable
d['rep_prior'] = d.rep_prior.fillna(False).astype(bool) & d.testable
d['trait'] = [F.trait_label(e, o[:80] if isinstance(o, str) else '') for e, o in zip(d.efo, d.outcome)]
d['area'] = F.area_group(F.area_of(d.efo))
d['n'] = pd.to_numeric(d.n, errors='coerce')
genes = F.genes()
d['gene'] = d.rsid.map(genes).fillna('')
cnt = d[d.status == 'positive'].groupby('rsid').size()
gl = F.gene_label(cnt[cnt >= 5].index.tolist())                 # dbSNP gene for the much-claimed variants (cached)
d['gene'] = [gl.get(r) or g for r, g in zip(d.rsid, d.gene)]

cv = pd.read_csv(F.ROOT / 'raw' / 'variant_summary.txt.gz', sep='\t', usecols=['ClinicalSignificance', 'RS# (dbSNP)', 'Assembly'], dtype=str)
cv = cv[cv.Assembly == 'GRCh38'].copy(); cv['rsid'] = pd.to_numeric(cv['RS# (dbSNP)'], errors='coerce'); cv = cv[cv.rsid > 0]
sig = cv.ClinicalSignificance.fillna('').str.lower().str.replace('association not found', '', regex=False)
cv['a'] = np.where(sig.str.contains('risk factor'), 'risk factor', np.where(sig.str.contains('association'), 'association',
                   np.where(sig.str.contains('protective'), 'protective', '')))
CV = cv[cv.a != ''].sort_values('a', key=lambda s: s.map({'risk factor': 0, 'association': 1, 'protective': 2})).drop_duplicates('rsid')
CV = dict(zip(CV.rsid.astype('int64'), CV.a))

r1 = lambda x: None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 3)   # noqa: E731

# ---- per variant
rows = []
for rs, x in d.groupby('rsid'):
    p = x[x.status == 'positive']; t = p[p.testable]
    sz = []
    for lo, hi in zip(SIZE[:-1], SIZE[1:]):
        s = t[(t.n >= lo) & (t.n < hi)]
        sz.append([len(s), int(s.replicated.sum())])
    top = p.groupby('trait').agg(k=('pmid', 'size'), s=('replicated', 'sum'), tt=('testable', 'sum')).sort_values('k', ascending=False).head(6)
    rows.append([int(rs), x.gene.iloc[0], int(x.pmid.nunique()), len(p), int((x.status == 'null').sum()), len(t),
                 r1(t.replicated.mean()) if len(t) else None, r1(t.null_share.mean()) if len(t) else None,
                 None if p.n.isna().all() else int(np.nanmedian(p.n)), sz,
                 [[k, int(v.k), int(v.tt), int(v.s)] for k, v in top.iterrows()],
                 CV.get(int(rs), ''), None if x.first_hit_year.isna().all() else int(x.first_hit_year.min())])
cols = ['rsid', 'gene', 'papers', 'reported', 'absent', 'testable', 'supported', 'null', 'median_n', 'size_classes', 'top_traits', 'clinvar', 'first_gw_year']
(OUT / 'variants.json').write_text(json.dumps({'cols': cols, 'rows': rows}, separators=(',', ':')))

# ---- per trait
tr = []
for k, x in d.groupby('trait'):
    if not k:
        continue
    p = x[x.status == 'positive']; t = p[p.testable]
    top = p.groupby('rsid').agg(k=('pmid', 'size'), s=('replicated', 'sum'), tt=('testable', 'sum'), g=('gene', 'first')).sort_values('k', ascending=False).head(8)
    tr.append([k, x.area.iloc[0], int(x.pmid.nunique()), len(p), int((x.status == 'null').sum()), len(t),
               r1(t.replicated.mean()) if len(t) else None, r1(t.null_share.mean()) if len(t) else None,
               [[int(rs), v.g, int(v.k), int(v.tt), int(v.s)] for rs, v in top.iterrows()]])
tcols = ['trait', 'area', 'papers', 'reported', 'absent', 'testable', 'supported', 'null', 'top_variants']
(OUT / 'traits.json').write_text(json.dumps({'cols': tcols, 'rows': sorted(tr, key=lambda r: -r[3])}, separators=(',', ':')))

# ---- claims, sharded
ccols = ['pmid', 'rsid', 'trait', 'status', 'year', 'n', 'design', 'population', 'testable', 'supported', 'earlier', 'first_gw_year']
pop = {'EUR': 'European', 'EAS': 'East Asian', 'other': 'Other or mixed', 'not_stated': 'Not stated'}
d['shard'] = d.rsid % NSHARD
for s, x in d.groupby('shard'):
    out = [[int(r.pmid), int(r.rsid), r.trait, 'reported' if r.status == 'positive' else 'absent', int(r.pub_year),
            None if np.isnan(r.n) else int(r.n), r.design, pop.get(r.pop_group, r.pop_group), bool(r.testable), bool(r.replicated),
            bool(r.rep_prior), None if pd.isna(r.first_hit_year) else int(r.first_hit_year)] for r in x.sort_values(['rsid', 'pub_year']).itertuples()]
    (OUT / 'claims' / f'{s:02d}.json').write_text(json.dumps({'cols': ccols, 'rows': out}, separators=(',', ':')))

meta = dict(claims=len(d), papers=int(d.pmid.nunique()), variants=len(rows), traits=len(tr), shards=NSHARD, size_classes=['<500', '500–1,999', '2,000–9,999', '≥10,000'],
            gwas_catalog='29 October 2026', clinvar='16 November 2026', pubmed='2026 baseline', efo='v3.95.0', built='2026-11-17')
(OUT / 'meta.json').write_text(json.dumps(meta, indent=1))
print(meta, sum(f.stat().st_size for f in OUT.rglob('*.json')) // 1e6, 'MB')

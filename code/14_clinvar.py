"""Additional analysis (2026-11-16, after the AJHG desk rejection; not part of the dated plan, reported as an additional analysis).
Definitions written here before ClinVar was read.

Question: which literature claims reach clinical variant curation, and how well are they supported genome-wide?
ClinVar (variant_summary.txt.gz, NCBI release of 2026-11-16) is the archive clinical laboratories and curators consult. Besides
Mendelian pathogenicity it holds non-Mendelian assertions ("risk factor", "association", "protective"), which are mostly drawn from
published association studies.

Definitions
- ClinVar rows on GRCh38 with a dbSNP identifier (RS# > 0). An rsID is a ClinVar risk variant when any of its rows has a clinical
  significance containing "risk factor", "association" (not "association not found") or "protective". The gene is ClinVar's
  GeneSymbol of that row (first symbol when several); the phenotypes are the PhenotypeList of those rows.
- Claims: the held-out analysis set of the primary estimates (11_confirm.load(): held-out papers, testable, disease/trait/severity
  claims, reported associations and reports of no association), reported associations only.
- Per group (claims on ClinVar risk variants; claims on all other variants): claims, papers, share with genome-wide support, mean
  matched-null share, and the excess (support minus null) with a paper-cluster bootstrap 95% CI (1,000 resamples, seed 20261108).
- Per ClinVar risk variant claimed in the held-out papers: claims, papers, distinct traits, share supported, mean null share, and
  whether the GWAS Catalog holds any genome-wide significant association (P <= 5e-8) at that rsID for any trait.
Outputs: tables additional_clinvar_groups, additional_clinvar_variants.
"""
import sys
from pathlib import Path
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import numpy as np, pandas as pd, duckdb
from rdata.proj import Project
HERE = Path(__file__).resolve().parent
P = Project(__file__)

src = (HERE / '11_confirm.py').read_text().split('\nwith P.run(')[0]
ns = {'__file__': str(HERE / '11_confirm.py'), '__name__': 'confirm'}
argv, sys.argv = sys.argv, ['11_confirm.py', 'heldout']
exec(compile(src, str(HERE / '11_confirm.py'), 'exec'), ns)
sys.argv = argv

CV = P.root / 'raw' / 'variant_summary.txt.gz'
SEED, B = 20261108, 1000


def excess(x):
    g = x.groupby('pmid').agg(w=('rsid', 'size'), r=('replicated', 'sum'), u=('null_share', 'sum'))
    W, R, U = g.w.values.astype(float), g.r.values.astype(float), g.u.values
    rng = np.random.default_rng(SEED)
    m = rng.multinomial(len(g), np.full(len(g), 1 / len(g)), size=B).astype(float)
    bs = (m @ R - m @ U) / (m @ W)
    return (R.sum() - U.sum()) / W.sum(), np.quantile(bs, .025), np.quantile(bs, .975)


with P.run(__file__, params={'clinvar': 'variant_summary 2026-11-16', 'assembly': 'GRCh38', 'split': 'heldout',
                             'risk_terms': 'risk factor | association (not "association not found") | protective',
                             'bootstrap': B, 'seed': SEED}) as run:
    run.input('clinvar_variant_summary_202611', path=str(CV), note='NCBI ClinVar tab-delimited variant summary')
    cv = pd.read_csv(CV, sep='\t', usecols=['GeneSymbol', 'ClinicalSignificance', 'RS# (dbSNP)', 'PhenotypeList', 'Assembly',
                                            'ReviewStatus'], dtype=str, compression='gzip')
    cv = cv[(cv.Assembly == 'GRCh38')].copy()
    cv['rsid'] = pd.to_numeric(cv['RS# (dbSNP)'], errors='coerce')
    cv = cv[cv.rsid > 0]
    sig = cv.ClinicalSignificance.fillna('').str.lower().str.replace('association not found', '', regex=False)
    cv['risk'] = sig.str.contains('risk factor') | sig.str.contains('association') | sig.str.contains('protective')
    rv = cv[cv.risk].copy()
    rv['gene'] = rv.GeneSymbol.fillna('').str.split(r'[;|,]').str[0]
    rvar = rv.groupby('rsid').agg(gene=('gene', 'first'), significance=('ClinicalSignificance', lambda s: '; '.join(sorted(set(s)))[:200]),
                                  phenotypes=('PhenotypeList', lambda s: '; '.join(sorted(set('|'.join(s.fillna('')).split('|')) - {''}))[:300]),
                                  records=('risk', 'size')).reset_index()
    rvar['rsid'] = rvar.rsid.astype('int64')
    print(len(cv.rsid.unique()), 'ClinVar rsIDs on GRCh38;', len(rvar), 'risk / association / protective rsIDs')

    d = ns['load']()
    pos = d[d.status == 'positive'].copy()
    pos['clinvar_risk'] = pos.rsid.isin(set(rvar.rsid))
    rows = []
    for name, x in (('ClinVar risk variants', pos[pos.clinvar_risk]), ('Other variants', pos[~pos.clinvar_risk])):
        e, lo, hi = excess(x)
        rows.append(dict(group=name, variants=x.rsid.nunique(), claims=len(x), papers=x.pmid.nunique(), supported=x.replicated.mean(),
                         null=x.null_share.fillna(0).mean(), excess=e, ci_low=lo, ci_high=hi))
    G = pd.DataFrame(rows)
    run.table('additional_clinvar_groups', G, note='held-out reported associations on ClinVar risk/association/protective variants vs other variants: support, matched null, excess (additional analysis, 2026-11-16)')

    con = duckdb.connect(str(P.root / 'db' / 'varclaims.duckdb'), read_only=True)
    gw = set(con.sql("select distinct rsid from gwas_catalog_rsid_trait where pvalue_mlog >= 7.30103 and rsid is not null").df().rsid)
    con.close()
    x = pos[pos.clinvar_risk]
    V = x.groupby('rsid').agg(claims=('pmid', 'size'), papers=('pmid', 'nunique'), traits=('trait_key', 'nunique'),
                              supported=('replicated', 'mean'), null=('null_share', 'mean')).reset_index()
    V = V.merge(rvar, on='rsid', how='left')
    V['gwas_any'] = V.rsid.isin(gw)
    V = V.sort_values('claims', ascending=False)
    run.table('additional_clinvar_variants', V, note='ClinVar risk/association/protective variants claimed in held-out papers: claims, support, null, any genome-wide association (additional analysis, 2026-11-16)')
    pd.set_option('display.width', 220)
    print(G.to_string(index=False, float_format=lambda v: f'{v:.4f}'))
    print(V.gwas_any.value_counts().to_string())
    print(V.head(30)[['rsid', 'gene', 'claims', 'papers', 'traits', 'supported', 'null', 'gwas_any', 'significance']].to_string(index=False, float_format=lambda v: f'{v:.3f}'))

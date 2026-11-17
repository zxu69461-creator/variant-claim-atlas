"""Additional analyses on the held-out papers (2026-11-10; not part of the dated plan, reported as additional analyses).
Both definitions were written here before either was run.

A. Claims from genome-wide studies and meta-analyses excluded. The five primary estimates are re-estimated, with the estimators of
   11_confirm.py unchanged, on held-out claims whose design is neither 'GWAS' nor 'meta-analysis' (the two design classes of the
   secondary breakdown), so that the estimates rest on case-control, cohort, family-based and other single-sample studies.
   -> table additional_sensitivity (same columns as confirm_sensitivity; set 'x_no_gwas_meta').

B. Agreement of the allele. Among held-out reported associations on disease risk or a quantitative trait that have genome-wide
   support and state an allele (A, C, G or T) with a direction (increased or decreased), the claim is compared with genome-wide
   significant GWAS Catalog associations of the same rsID (no proxies) for the same ontology term (equal terms only, so that the
   direction refers to the same outcome), from papers other than the claim's own, that give an allele (A, C, G or T) and a direction
   of effect (the catalog's "unit increase" / "unit decrease"; otherwise an odds ratio above or below 1). Alleles are mapped to the
   two alleles of the variant in the 1000 Genomes reference, allowing the complementary strand; A/T and C/G variants, alleles that
   match neither reference allele and catalog associations without a direction are left out. A claim agrees with a catalog
   association when both name the same allele with the same direction or opposite alleles with opposite directions; a claim compared
   with several catalog associations agrees when most of them do (ties: the association with the smallest P value decides).
   -> table additional_allele_agreement (group, claims, papers, agree, share, ci_low, ci_high; paper-cluster bootstrap, 1,000 resamples).
"""
import sys
from pathlib import Path
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import numpy as np, pandas as pd, duckdb
from rdata.proj import Project
HERE = Path(__file__).resolve().parent
P = Project(__file__)

# the estimators of the confirmatory script, loaded unchanged (everything above its run block)
src = (HERE / '11_confirm.py').read_text().split('\nwith P.run(')[0]
ns = {'__file__': str(HERE / '11_confirm.py'), '__name__': 'confirm'}
argv, sys.argv = sys.argv, ['11_confirm.py', 'heldout']
exec(compile(src, str(HERE / '11_confirm.py'), 'exec'), ns)
sys.argv = argv
sys.path.insert(0, str(HERE)); import figdata as F

COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def to_ref(allele, a1, a2):
    """0 / 1 = the allele is the first / second reference allele (either strand); None = it matches neither."""
    if allele in (a1, COMP[a1]): return 0
    if allele in (a2, COMP[a2]): return 1
    return None


with P.run(__file__, params={'split': 'heldout', 'excluded_designs': 'GWAS, meta-analysis', 'bootstrap': 1000, 'seed': 20261108}) as run:
    # ---- A
    d = ns['load']()
    sub = d[~d.design.isin(['GWAS', 'meta-analysis'])]
    fam = ns['family'](sub).assign(split='heldout', set='x_no_gwas_meta')
    run.table('additional_sensitivity', fam, note='held-out: primary estimates C1-C5 without claims from genome-wide studies and meta-analyses (additional analysis, 2026-11-10)')
    pd.set_option('display.width', 200)
    print(len(d), 'claims in the analysis set;', len(sub), 'after excluding genome-wide studies and meta-analyses')
    print(fam[['test', 'estimate', 'ci_low', 'ci_high', 'p', 'n', 'n_papers', 'extra']].to_string(index=False, float_format=lambda x: f'{x:.4g}'))

    # ---- B
    con = duckdb.connect(str(P.root / 'db' / 'varclaims.duckdb'), read_only=True)
    con.execute(f"set file_search_path = '{P.root}'")
    cl = con.sql("""select v.pmid, v.rsid, v.trait_key, v.pop_group, v.outcome_type, v.risk_allele, v.direction, v.efo
                    from variant_claims v join variant_claim_match_heldout m using (pmid, rsid, trait_key, pop_group)
                    where v.outcome_type in ('disease', 'trait') and v.status = 'positive' and m.testable and m.replicated
                      and v.risk_allele in ('A', 'C', 'G', 'T') and v.direction in ('increased', 'decreased')""").df()
    g = con.sql("""select try_cast(snp_id_current as bigint) rsid, pubmedid, try_cast(pvalue_mlog as double) mlog,
                          regexp_extract(strongest_snp_risk_allele, '-([ACGT])$', 1) g_allele, try_cast(or_or_beta as double) eff,
                          "95pct_ci_text" ci, mapped_trait_uri
                   from gwas_catalog_assoc where try_cast(pvalue_mlog as double) >= 7.30103 and mapped_trait_uri is not null
                     and try_cast(snp_id_current as bigint) is not null""").df()
    con.close()
    g = g[g.g_allele.isin(list('ACGT'))].copy()
    ci = g.ci.fillna('')
    g['g_sign'] = np.where(ci.str.contains('increase'), 1, np.where(ci.str.contains('decrease'), -1,
                  np.where(g.eff > 1, 1, np.where((g.eff > 0) & (g.eff < 1), -1, 0))))
    g = g[g.g_sign != 0]
    g = g.assign(gefo=g.mapped_trait_uri.str.split(', ')).explode('gefo')
    g['gefo'] = g.gefo.str.extract(r'([^/]+)$')[0]
    cl['claim_id'] = np.arange(len(cl))
    ce = cl[['claim_id', 'pmid', 'rsid', 'efo']].explode('efo').rename(columns={'efo': 'gefo'})
    j = ce.merge(g[['rsid', 'gefo', 'pubmedid', 'mlog', 'g_allele', 'g_sign']], on=['rsid', 'gefo'])
    j = j[j.pubmedid != j.pmid.astype(str)].drop_duplicates(['claim_id', 'pubmedid', 'mlog', 'g_allele', 'g_sign'])
    fr = pd.concat([pd.read_csv(P.root / 'work' / 'null' / f'{p}_freq.frq.gz', sep=r'\s+', usecols=['SNP', 'A1', 'A2']) for p in ('eur', 'eas')])
    fr = fr[fr.SNP.str.match(r'^rs\d+$') & fr.A1.isin(list('ACGT')) & fr.A2.isin(list('ACGT'))].drop_duplicates('SNP')
    fr['rsid'] = fr.SNP.str[2:].astype('int64')
    fr = fr[fr.A1.map(COMP) != fr.A2]                                         # A/T and C/G variants left out
    j = j.merge(fr[['rsid', 'A1', 'A2']], on='rsid').merge(cl[['claim_id', 'risk_allele', 'direction', 'outcome_type']], on='claim_id')
    j['c_idx'] = [to_ref(a, x, y) for a, x, y in zip(j.risk_allele, j.A1, j.A2)]
    j['g_idx'] = [to_ref(a, x, y) for a, x, y in zip(j.g_allele, j.A1, j.A2)]
    j = j[j.c_idx.notna() & j.g_idx.notna()].copy()
    j['c_sign'] = np.where(j.direction == 'increased', 1, -1)
    j['agree'] = (j.c_idx == j.g_idx) == (j.c_sign == j.g_sign)
    j = j.sort_values(['claim_id', 'mlog'], ascending=[True, False])
    per = j.groupby('claim_id').agg(pmid=('pmid', 'first'), outcome_type=('outcome_type', 'first'), k=('agree', 'size'), a=('agree', 'sum'),
                                    top=('agree', 'first')).reset_index()
    per['agree'] = np.where(2 * per.a > per.k, True, np.where(2 * per.a < per.k, False, per.top))
    rows = []
    for name, x in (('Disease risk or quantitative trait', per), ('Disease risk', per[per.outcome_type == 'disease']),
                    ('Quantitative trait', per[per.outcome_type == 'trait'])):
        est, lo, hi = F.boot_share(x, 'agree')
        rows.append(dict(group=name, claims=len(x), papers=x.pmid.nunique(), agree=int(x.agree.sum()), share=est, ci_low=lo, ci_high=hi))
    res = pd.DataFrame(rows)
    run.table('additional_allele_agreement', res,
              note='held-out supported reported associations stating an allele and direction: agreement with the GWAS Catalog allele and direction at the same rsID and ontology term (additional analysis, 2026-11-10)')
    print(len(cl), 'supported claims with an allele and direction;', len(per), 'comparable with a catalog association')
    print(res.to_string(index=False, float_format=lambda x: f'{x:.4f}'))

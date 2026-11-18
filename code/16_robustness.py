"""Additional analyses on the held-out papers (2026-11-17, after an external pre-submission review; not part of the dated plan,
reported as additional analyses). Definitions written here before any of them was run.

A. Narrower trait matching. The matching rule of the primary analysis accepts a GWAS Catalog term equal to the claim term or any
   ancestor or descendant of it (disease / measurement / phenotype branches, roots excluded). Two narrower rules are applied to the
   same analysis set (held-out, testable, the claims of 11_confirm.load()), with the matched null recomputed under the same rule:
   (exact) the two terms are equal; (one step) equal, or one is the direct parent of the other in EFO (efo_isa). Support, earlier
   support, first hit year and null share are recomputed with everything else unchanged (variant or LD proxy, P <= 5e-8, a paper
   other than the claim's own, reference population by claim population). The five primary estimates are re-estimated with the
   estimators of 11_confirm.py.   -> table additional_trait_rule (columns of confirm_sensitivity; set 'exact' / 'one_step').

B. ClinVar conditions. For the ClinVar risk variants of Figure 4 (held-out, >= 5 reported associations, classes mostly >= 50%,
   partly 5-50%, rarely < 5% supported): (i) whether the variant or a European LD proxy has any genome-wide significant association
   (P <= 5e-8) with any trait in the GWAS Catalog; (ii) whether it has one for a condition that ClinVar itself lists for the variant
   on its risk-factor / association / protective records (PhenotypeIDS MONDO and HP identifiers that are EFO terms in the three
   branches; matched with the rule of the primary analysis). Variants whose conditions carry no such identifier are counted apart.
   -> table additional_clinvar_condition (one row per variant).

C. Allele and direction, restated from 13_additional.py part B (already computed 2026-11-10, non-palindromic variants, exact rsID,
   exact term): the agreement is tested against the 50% expected if allele and direction were unrelated (binomial test on claims).
"""
import sys
from pathlib import Path
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import numpy as np, pandas as pd, duckdb
from scipy import stats
from rdata.proj import Project
HERE = Path(__file__).resolve().parent
P = Project(__file__); N = P.root / 'work' / 'null'

src = (HERE / '11_confirm.py').read_text().split('\nwith P.run(')[0]
ns = {'__file__': str(HERE / '11_confirm.py'), '__name__': 'confirm'}
argv, sys.argv = sys.argv, ['11_confirm.py', 'heldout']
exec(compile(src, str(HERE / '11_confirm.py'), 'exec'), ns)
sys.argv = argv
m8 = (HERE / '08_match.py').read_text()
exec(compile(m8.split('\nwith P.run(')[0].split('P = Project(__file__)')[1].split('\n', 1)[1], '08_match', 'exec'), ns)   # read_ld, read_frq, mbin, BINS
ROOTS = ('MONDO_0000001', 'EFO_0001444', 'EFO_0000651')
CV = P.root / 'raw' / 'variant_summary.txt.gz'


def rematch(con, d, rule):
    """Support and null share for the claims in d under a narrower trait rule ('exact' or 'one_step')."""
    c = con.sql("""select pmid, rsid, trait_key, pop_group, efo, pub_year from variant_claims where outcome_type <> 'drug_response' and not dev""").df()
    c = d[['pmid', 'rsid', 'trait_key', 'pop_group']].drop_duplicates().merge(c, on=['pmid', 'rsid', 'trait_key', 'pop_group'])
    c['claim_id'] = np.arange(len(c))
    con.register('c', c)
    con.execute("""create or replace temp table branch as select efo_id from efo_ancestors where ancestor in ('MONDO_0000001','EFO_0001444','EFO_0000651')
                   union select efo_id from efo_terms where efo_id like 'MONDO_%'""")
    con.execute(f"delete from branch where efo_id in {ROOTS}")
    con.execute("create or replace temp table ct as select distinct claim_id, unnest(efo) as efo from c")
    con.execute("delete from ct where efo not in (select efo_id from branch)")
    con.execute("""create or replace temp table g as select rsid, pubmedid, try_cast(left(date, 4) as int) as gyear, pvalue_mlog,
                   regexp_extract(efo_uri, '[^/]+$') as gefo from gwas_catalog_rsid_trait where efo_uri is not null""")
    con.execute("delete from g where gefo not in (select efo_id from branch)")
    one = "" if rule == 'exact' else """
        union select i.child, i.parent from efo_isa i join e on i.child = e.efo join t on i.parent = t.gefo
        union select i.parent, i.child from efo_isa i join e on i.parent = e.efo join t on i.child = t.gefo"""
    con.execute(f"""create or replace temp table tm as with e as (select distinct efo from ct), t as (select distinct gefo from g)
                    select e.efo, t.gefo from e join t on e.efo = t.gefo {one}""")
    con.execute("""create or replace temp table cp as select claim_id, rsid as prsid from c
        union select c.claim_id, l.proxy_rsid from c join ld_proxies l on l.query_rsid = c.rsid
              and l.pop = case when c.pop_group = 'EAS' then 'eas' else 'eur' end""")
    hits = con.sql("""select cp.claim_id, min(g.gyear) as first_hit_year, bool_or(g.gyear < c.pub_year) as rep_prior
        from cp join c on c.claim_id = cp.claim_id join ct on ct.claim_id = cp.claim_id
        join tm on tm.efo = ct.efo join g on g.rsid = cp.prsid and g.gefo = tm.gefo
        where g.pvalue_mlog >= 7.30103 and g.pubmedid <> cast(c.pmid as varchar) group by 1""").df()
    nl = pd.concat([ns['read_ld'](N / 'null_eur_proxies.ld.gz', 'eur'), ns['read_ld'](N / 'null_eas_proxies.ld.gz', 'eas')])
    nulls = pd.read_csv(N / 'null_rsids.txt', header=None, names=['r']).r.str[2:].astype('int64')
    frq = ns['read_frq'](N / 'eur_freq.frq.gz')
    nb = pd.DataFrame({'null_rsid': nulls}).merge(frq, left_on='null_rsid', right_on='rsid', how='left'); nb['bin'] = ns['mbin'](nb.maf)
    con.register('nb', nb[['null_rsid', 'bin']]); con.register('nl', nl)
    con.execute("""create or replace temp table np as select pop, null_rsid, null_rsid as prsid from nb, (values ('eur'), ('eas')) v(pop)
                   union select pop, query_rsid, proxy_rsid from nl""")
    con.execute("create or replace temp table nh as select distinct np.pop, np.null_rsid, g.gefo from np join g on g.rsid = np.prsid where g.pvalue_mlog >= 7.30103")
    c['bin'] = ns['mbin'](c.rsid.map(frq.set_index('rsid').maf)); c['ref_pop'] = np.where(c.pop_group == 'EAS', 'eas', 'eur')
    con.register('c2', c[['claim_id', 'bin', 'ref_pop']])
    con.execute("create or replace temp table sig as select claim_id, string_agg(efo, '|' order by efo) as tset from ct group by 1")
    ns_ = con.sql("""
        with k as (select distinct s.tset, c2.bin, c2.ref_pop from sig s join c2 using (claim_id)),
        kt as (select k.*, unnest(string_split(k.tset, '|')) as efo from k),
        km as (select distinct kt.tset, kt.bin, kt.ref_pop, tm.gefo from kt join tm on tm.efo = kt.efo),
        hit as (select km.tset, km.bin, km.ref_pop, count(distinct nh.null_rsid) as h from km
                join nh on nh.gefo = km.gefo and nh.pop = km.ref_pop join nb on nb.null_rsid = nh.null_rsid and nb.bin = km.bin group by 1, 2, 3),
        den as (select bin, count(*) as d from nb group by 1)
        select s.claim_id, coalesce(hit.h, 0) / den.d as null_share from sig s join c2 using (claim_id)
        left join hit on hit.tset = s.tset and hit.bin = c2.bin and hit.ref_pop = c2.ref_pop join den on den.bin = c2.bin""").df()
    m = c[['claim_id', 'pmid', 'rsid', 'trait_key', 'pop_group']].merge(hits, on='claim_id', how='left').merge(ns_, on='claim_id', how='left')
    m['replicated'] = m.first_hit_year.notna(); m['rep_prior'] = m.rep_prior.fillna(False).astype(bool)
    m['null_share'] = m.null_share.fillna(0.0)
    out = d.drop(columns=['replicated', 'rep_prior', 'first_hit_year', 'null_share']).merge(
        m.drop(columns='claim_id'), on=['pmid', 'rsid', 'trait_key', 'pop_group'], how='left')
    out['replicated'] = out.replicated.fillna(False).astype(bool); out['rep_prior'] = out.rep_prior.fillna(False).astype(bool)
    out['null_share'] = out.null_share.fillna(0.0)
    return out


with P.run(__file__, params={'split': 'heldout', 'rules': 'exact, one_step (efo_isa)', 'clinvar': 'variant_summary 2026-11-16'}) as run:
    run.input('clinvar_variant_summary_202611', path=str(CV), note='NCBI ClinVar tab-delimited variant summary')
    d = ns['load']()
    con = duckdb.connect(str(P.root / 'db' / 'varclaims.duckdb'), read_only=True)
    con.execute(f"set file_search_path = '{P.root}'")
    fam = []
    for rule in ('exact', 'one_step'):
        x = rematch(con, d, rule)
        f = ns['family'](x).assign(split='heldout', set=rule)
        f['support'] = x[x.status == 'positive'].replicated.mean(); f['null'] = x[x.status == 'positive'].null_share.mean()
        fam.append(f)
        print(rule, f"support {f.support[0]:.4f} null {f.null[0]:.4f}")
        print(f[['test', 'estimate', 'ci_low', 'ci_high', 'p', 'n']].to_string(index=False, float_format=lambda v: f'{v:.4g}'))

    # ---- B
    cv = pd.read_csv(CV, sep='\t', usecols=['ClinicalSignificance', 'RS# (dbSNP)', 'Assembly', 'PhenotypeIDS', 'PhenotypeList'], dtype=str)
    cv = cv[cv.Assembly == 'GRCh38'].copy(); cv['rsid'] = pd.to_numeric(cv['RS# (dbSNP)'], errors='coerce'); cv = cv[cv.rsid > 0]
    sig = cv.ClinicalSignificance.fillna('').str.lower().str.replace('association not found', '', regex=False)
    cv = cv[sig.str.contains('risk factor') | sig.str.contains('association') | sig.str.contains('protective')]
    ids = cv.assign(t=cv.PhenotypeIDS.fillna('').str.split(r'[|,;]')).explode('t')
    ids['efo'] = ids.t.str.extract(r'(MONDO:\d{7}|HP:\d{7})$')[0].str.replace(':', '_')
    ids = ids.dropna(subset=['efo'])[['rsid', 'efo']].drop_duplicates(); ids['rsid'] = ids.rsid.astype('int64')
    pos = d[d.status == 'positive']
    V = pos.groupby('rsid').agg(n=('pmid', 'size'), sup=('replicated', 'mean')).reset_index()
    V = V[V.rsid.isin(set(cv.rsid.astype('int64')))].copy()                 # every claimed ClinVar risk variant (2026-11-17: all, not only >= 5)
    V['cls'] = np.where(V.n < 5, 'under5', np.where(V.sup >= 0.5, 'mostly', np.where(V.sup >= 0.05, 'partly', 'rarely')))
    con.register('V', V[['rsid']]); con.register('ids', ids[ids.rsid.isin(set(V.rsid))])
    con.execute("""create or replace temp table vp as select rsid, rsid prsid from V
                   union select V.rsid, l.proxy_rsid from V join ld_proxies l on l.query_rsid = V.rsid""")   # 1000 Genomes EUR and EAS proxies
    any_ = set(con.sql("""select distinct vp.rsid from vp join gwas_catalog_rsid_trait g on g.rsid = vp.prsid where g.pvalue_mlog >= 7.30103""").df().rsid)
    con.execute("create or replace temp table ib as select * from ids where efo in (select efo_id from branch)")
    with_term = set(con.sql('select distinct rsid from ib').df().rsid)
    cond = set(con.sql("""with gg as (select rsid, regexp_extract(efo_uri, '[^/]+$') gefo from gwas_catalog_rsid_trait
                                      where pvalue_mlog >= 7.30103 and efo_uri is not null),
            tmm as (select e.efo, e.efo gefo from (select distinct efo from ib) e
                    union select a.efo_id, a.ancestor from efo_ancestors a where a.efo_id in (select efo from ib)
                    union select a.ancestor, a.efo_id from efo_ancestors a where a.ancestor in (select efo from ib))
        select distinct ib.rsid from ib join tmm on tmm.efo = ib.efo join vp on vp.rsid = ib.rsid
        join gg on gg.rsid = vp.prsid and gg.gefo = tmm.gefo where tmm.gefo not in ('MONDO_0000001','EFO_0001444','EFO_0000651')""").df().rsid)
    con.close()
    run.table('additional_trait_rule', pd.concat(fam), note='held-out: primary estimates C1-C5 under exact and one-step EFO trait matching, null recomputed (additional analysis, 2026-11-17)')
    pl = cv.assign(rsid=cv.rsid.astype('int64')).groupby('rsid').PhenotypeList.agg(lambda x: '|'.join(x.fillna('')))
    ids_any = cv.assign(rsid=cv.rsid.astype('int64')).groupby('rsid').PhenotypeIDS.agg(lambda x: '|'.join(x.fillna('')))
    def why(r):
        if r in with_term: return ''
        txt = pl.get(r, '').lower(); idt = ids_any.get(r, '')
        conds = [c for c in txt.split('|') if c and c not in ('not provided', 'not specified', 'see cases')]
        if not conds: return 'no named condition (not provided / not specified)'
        if 'MONDO:' in idt or 'HP:' in idt: return 'condition term outside the disease, measurement and phenotype branches'
        return 'condition without a Mondo or HPO identifier'
    V['no_term_reason'] = [why(r) for r in V.rsid]
    V['gwas_any'] = V.rsid.isin(any_); V['cond_term'] = V.rsid.isin(with_term); V['cond_support'] = V.rsid.isin(cond)
    run.table('additional_clinvar_condition', V, note='held-out ClinVar risk variants (>=5 reported associations): any genome-wide association, ClinVar condition with an EFO/MONDO/HP term, genome-wide support for that condition (additional analysis, 2026-11-17)')
    print(V.groupby('cls').agg(variants=('rsid', 'size'), gwas_any=('gwas_any', 'sum'), cond_term=('cond_term', 'sum'),
                               cond_support=('cond_support', 'sum')).to_string())

    # ---- C
    with duckdb.connect(str(P.root / 'db' / 'varclaims.duckdb'), read_only=True) as c2:
        a = c2.sql("select * from additional_allele_agreement").df()
    a['p_vs_half'] = [stats.binomtest(int(k), int(n), 0.5).pvalue for k, n in zip(a.agree, a.claims)]
    print(a.to_string(index=False))

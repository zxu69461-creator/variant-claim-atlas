"""Arbiter A (charter amendments 1-3): each variant claim against the GWAS Catalog. A claim is replicated when its rsID or an LD proxy
(r2 >= 0.8, 1000G EAS for EAS claims, EUR otherwise) has a P < 5e-8 association for a matching EFO term (equal, or one an ancestor of the
other, both inside the disease / measurement / phenotype branches, roots excluded) from a paper other than the claim's own. Prior support
= a GWAS paper dated before the claim's publication year; later = same year or after. Testable = a claim term, or a descendant, has any
GWAS Catalog association; testable_prior = such an association was published before the claim year. Null (amendment 3): share of
20,000 random GWAS Catalog rsIDs in the claim variant's 1000G EUR MAF bin whose own variant or proxy hits a matching term."""
import sys
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import pandas as pd, numpy as np, duckdb
from rdata.proj import Project
P = Project(__file__); N = P.root / 'work' / 'null'
SPLIT = sys.argv[1] if len(sys.argv) > 1 else 'dev'   # 'dev' (development) or 'heldout' (amendment 4, after the reader check)
assert SPLIT in ('dev', 'heldout')
OUT = 'variant_claim_match' if SPLIT == 'dev' else 'variant_claim_match_heldout'
ROOTS = ('MONDO_0000001', 'EFO_0001444', 'EFO_0000651')
BINS = [0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.5000001]

def read_ld(f, pop):
    d = pd.read_csv(f, sep=r'\s+', usecols=['SNP_A', 'SNP_B'])
    d = d[d.SNP_A.str.match(r'^rs\d+$') & d.SNP_B.str.match(r'^rs\d+$')]
    return pd.DataFrame({'pop': pop, 'query_rsid': d.SNP_A.str[2:].astype('int64'), 'proxy_rsid': d.SNP_B.str[2:].astype('int64')})

def read_frq(f):
    d = pd.read_csv(f, sep=r'\s+', usecols=['SNP', 'MAF'])
    d = d[d.SNP.str.match(r'^rs\d+$')]
    return pd.DataFrame({'rsid': d.SNP.str[2:].astype('int64'), 'maf': pd.to_numeric(d.MAF, errors='coerce')}).drop_duplicates('rsid')

def mbin(m):
    m = np.asarray(m, float)
    b = np.digitize(m, BINS).astype(object)
    b[np.isnan(m)] = 'absent'; b[(~np.isnan(m)) & (m < 0.01)] = 'lt0.01'
    return b.astype(str)

with P.run(__file__, params={'p_threshold': 5e-8, 'r2_min': 0.8, 'null_n': 20000, 'null_seed': 20261107,
                             'maf_bins': '0.01,0.05,0.10,0.20,0.30,0.40,0.50', 'split': SPLIT}) as run:
    for f in ('null_eur_proxies.ld.gz', 'null_eas_proxies.ld.gz', 'eur_freq.frq.gz'): run.input(f, path=N / f)
    con = duckdb.connect(str(P.root / 'db' / 'varclaims.duckdb'), read_only=True)
    c = con.sql("""select row_number() over (order by pmid, rsid, trait_key, pop_group) as claim_id, pmid, rsid, trait_key, pop_group,
                   efo, pub_year from variant_claims where outcome_type <> 'drug_response' and """ + ('dev' if SPLIT == 'dev' else 'not dev')).df()
    con.register('c', c)
    con.execute("""create temp table branch as select efo_id from efo_ancestors where ancestor in ('MONDO_0000001','EFO_0001444','EFO_0000651')
                   union select efo_id from efo_terms where efo_id like 'MONDO_%'""")
    con.execute(f"delete from branch where efo_id in {ROOTS}")
    con.execute("""create temp table ct as select distinct claim_id, unnest(efo) as efo from c""")
    con.execute("delete from ct where efo not in (select efo_id from branch)")
    con.execute("""create temp table g as select rsid, pubmedid, try_cast(left(date, 4) as int) as gyear, pvalue_mlog,
                   regexp_extract(efo_uri, '[^/]+$') as gefo from gwas_catalog_rsid_trait where efo_uri is not null""")
    con.execute("delete from g where gefo not in (select efo_id from branch)")
    # term pairs: claim term e matches GWAS term t when equal or either is an ancestor of the other
    con.execute("""create temp table tm as
        with e as (select distinct efo from ct), t as (select distinct gefo from g)
        select e.efo, t.gefo from e join t on e.efo = t.gefo
        union select a.efo_id, a.ancestor from efo_ancestors a join e on a.efo_id = e.efo join t on a.ancestor = t.gefo
        union select a.ancestor, a.efo_id from efo_ancestors a join e on a.ancestor = e.efo join t on a.efo_id = t.gefo""")
    # testable: a claim term or a descendant has any GWAS association (any P); earliest year
    con.execute("""create temp table tt as
        with t as (select gefo, min(gyear) y from g group by 1)
        select ct.claim_id, min(t.y) as first_gwas_year from ct join (select efo, efo as d from (select distinct efo from ct)
             union select a.ancestor, a.efo_id from efo_ancestors a where a.ancestor in (select efo from ct)) dd on dd.efo = ct.efo
        join t on t.gefo = dd.d group by 1""")
    # proxies for claim variants: self plus LD proxies in the claim's reference population
    con.execute("""create temp table cp as select claim_id, rsid as prsid from c
        union select c.claim_id, l.proxy_rsid from c join ld_proxies l on l.query_rsid = c.rsid
              and l.pop = case when c.pop_group = 'EAS' then 'eas' else 'eur' end""")
    hits = con.sql("""select cp.claim_id, min(g.gyear) as first_hit_year, count(distinct g.pubmedid) as n_hit_papers,
            bool_or(g.gyear < c.pub_year) as rep_prior, bool_or(g.gyear >= c.pub_year) as rep_later
        from cp join c on c.claim_id = cp.claim_id join ct on ct.claim_id = cp.claim_id
        join tm on tm.efo = ct.efo join g on g.rsid = cp.prsid and g.gefo = tm.gefo
        where g.pvalue_mlog >= 7.30103 and g.pubmedid <> cast(c.pmid as varchar) group by 1""").df()
    tt = con.sql('select * from tt').df(); nterm = con.sql('select claim_id, count(*) n_branch_terms from ct group by 1').df()
    nprox = con.sql('select claim_id, count(*) - 1 as n_proxies from cp group by 1').df()

    # null sample
    nl = pd.concat([read_ld(N / 'null_eur_proxies.ld.gz', 'eur'), read_ld(N / 'null_eas_proxies.ld.gz', 'eas')])
    nulls = pd.read_csv(N / 'null_rsids.txt', header=None, names=['r']).r.str[2:].astype('int64')
    frq = read_frq(N / 'eur_freq.frq.gz')
    nb = pd.DataFrame({'null_rsid': nulls}).merge(frq, left_on='null_rsid', right_on='rsid', how='left')
    nb['bin'] = mbin(nb.maf)
    con.register('nb', nb[['null_rsid', 'bin']])
    con.register('nl', nl)
    con.execute("""create temp table np as select pop, null_rsid, null_rsid as prsid from nb, (values ('eur'), ('eas')) v(pop)
                   union select pop, query_rsid, proxy_rsid from nl""")
    con.execute("""create temp table nh as select distinct np.pop, np.null_rsid, g.gefo from np join g on g.rsid = np.prsid where g.pvalue_mlog >= 7.30103""")
    c['bin'] = mbin(c.rsid.map(frq.set_index('rsid').maf)); c['ref_pop'] = np.where(c.pop_group == 'EAS', 'eas', 'eur')
    con.register('c2', c[['claim_id', 'bin', 'ref_pop']])
    # null share per claim: computed per distinct (term set, bin, reference population)
    con.execute("""create temp table sig as select claim_id, string_agg(efo, '|' order by efo) as tset from ct group by 1""")
    nullshare = con.sql("""
        with k as (select distinct s.tset, c2.bin, c2.ref_pop from sig s join c2 using (claim_id)),
        kt as (select k.*, unnest(string_split(k.tset, '|')) as efo from k),
        km as (select distinct kt.tset, kt.bin, kt.ref_pop, tm.gefo from kt join tm on tm.efo = kt.efo),
        hit as (select km.tset, km.bin, km.ref_pop, count(distinct nh.null_rsid) as h from km
                join nh on nh.gefo = km.gefo and nh.pop = km.ref_pop join nb on nb.null_rsid = nh.null_rsid and nb.bin = km.bin group by 1, 2, 3),
        den as (select bin, count(*) as d from nb group by 1)
        select s.claim_id, coalesce(hit.h, 0) / den.d as null_share, den.d as null_n from sig s join c2 using (claim_id)
        left join hit on hit.tset = s.tset and hit.bin = c2.bin and hit.ref_pop = c2.ref_pop join den on den.bin = c2.bin""").df()
    con.close()

    m = c[['claim_id', 'pmid', 'rsid', 'trait_key', 'pop_group', 'bin']].rename(columns={'bin': 'maf_bin'})
    m['maf_eur'] = m.rsid.map(frq.set_index('rsid').maf)
    for d in (nterm, nprox, tt, hits, nullshare): m = m.merge(d, on='claim_id', how='left')
    m['n_branch_terms'] = m.n_branch_terms.fillna(0).astype(int)
    m['testable'] = m.first_gwas_year.notna()
    py = c.set_index('claim_id').pub_year
    m['testable_prior'] = m.first_gwas_year < m.claim_id.map(py)
    m['replicated'] = m.first_hit_year.notna()
    for col in ('rep_prior', 'rep_later'): m[col] = m[col].fillna(False).astype(bool)
    run.table(OUT, m.drop(columns=['claim_id']),
              note=f'{SPLIT} split: GWAS Catalog arbiter per variant claim (amendments 1-4): testability, replication any/prior/later, matched null share')
    print(len(m), 'claims;', f"testable {m.testable.mean():.3f}; replicated among testable {m[m.testable].replicated.mean():.3f}")

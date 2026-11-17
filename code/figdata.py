"""Shared data access for figures and the manuscript builder. `claims(split)` returns the analysis set of charter amendment 4.2 for one
split ('dev' or 'heldout'), with trait area (GWAS Catalog parent category reached through EFO ancestry) and the mapped gene of the
variant; `objects()` returns every claim record (both splits) for descriptive panels, which use no arbiter result."""
from functools import lru_cache
from pathlib import Path
import duckdb, numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / 'db' / 'varclaims.duckdb'
MATCH = {'dev': 'variant_claim_match', 'heldout': 'variant_claim_match_heldout'}
# GWAS Catalog parent categories, most specific first; a claim takes the first category any of its terms falls under
AREAS = [('Cancer', 'Cancer'), ('Cardiovascular disease', 'Cardiovascular disease'), ('Metabolic disorder', 'Metabolic disease'),
         ('Immune system disorder', 'Immune disease'), ('Neurological disorder', 'Neurological disease'),
         ('Digestive system disorder', 'Digestive disease'), ('Lipid or lipoprotein measurement', 'Lipid measurement'),
         ('Cardiovascular measurement', 'Cardiovascular measurement'), ('Body measurement', 'Body measurement'),
         ('Hematological measurement', 'Blood measurement'), ('Inflammatory measurement', 'Inflammatory marker'),
         ('Other disease', 'Other disease'), ('Other measurement', 'Other measurement'), ('Other trait', 'Other trait'),
         ('Biological process', 'Other trait')]


def con():
    c = duckdb.connect(str(DB), read_only=True)
    c.execute(f"set file_search_path = '{ROOT}'")        # Parquet-backed views use project-relative paths
    return c


@lru_cache(None)
def genes():
    """rsID -> most frequent mapped gene in the GWAS Catalog (ties broken alphabetically)."""
    with con() as c:
        return c.sql("""select rsid, arg_max(g, n) gene from (
            select try_cast(snp_id_current as bigint) rsid, trim(split_part(split_part(mapped_gene, ' - ', 1), ',', 1)) g, count(*) n
            from gwas_catalog_assoc where snp_id_current is not null and mapped_gene is not null and mapped_gene <> '' group by 1, 2)
            where rsid is not null group by 1""").df().set_index('rsid').gene


@lru_cache(None)
def _areas():
    """EFO term -> trait area. Terms the GWAS Catalog uses take its own parent category; other terms take the most specific category
    among their EFO ancestors that the GWAS Catalog uses."""
    with con() as c:
        par = c.sql("""select regexp_extract(efo_uri, '[^/]+$') efo, parent_term, count(*) n from gwas_efo_parent
                       where efo_uri is not null group by 1, 2""").df()
        anc = c.sql("select efo_id, ancestor from efo_ancestors").df()
    rank = {t: i for i, (t, _) in enumerate(AREAS)}
    label = dict(AREAS)
    par = par[par.parent_term.isin(rank)].assign(r=lambda z: z.parent_term.map(rank))
    own = par.sort_values(['efo', 'r']).drop_duplicates('efo')          # a term with several categories keeps the most specific
    omap = dict(zip(own.efo, own.parent_term))
    a = anc[anc.ancestor.isin(omap)].assign(term=lambda z: z.ancestor.map(omap))
    a['r'] = a.term.map(rank)
    inh = a.sort_values('r').drop_duplicates('efo_id')
    m = dict(zip(inh.efo_id, inh.term)); m.update(omap)
    return {k: label[v] for k, v in m.items()}, {label[t]: rank[t] for t, _ in AREAS}


def gene_label(rsids):
    """rsID -> first gene listed by dbSNP (NCBI esummary), cached in work/fig/dbsnp_genes.tsv; GWAS Catalog mapped gene as fallback."""
    import json, subprocess
    f = ROOT / 'work' / 'fig' / 'dbsnp_genes.tsv'
    cache = pd.read_csv(f, sep='\t', dtype={'rsid': 'int64', 'genes': str}).fillna('') if f.exists() else pd.DataFrame(columns=['rsid', 'genes'])
    miss = [int(r) for r in rsids if int(r) not in set(cache.rsid)]
    for i in range(0, len(miss), 200):
        ids = ','.join(map(str, miss[i:i + 200]))
        out = subprocess.run(['ncurl', '-s', f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=snp&id={ids}&retmode=json'],
                             capture_output=True, text=True, check=True).stdout
        j = json.loads(out)['result']
        add = pd.DataFrame([(int(u), ','.join(g['name'] for g in j[u].get('genes', []))) for u in j['uids']], columns=['rsid', 'genes'])
        if add.empty:
            raise RuntimeError('dbSNP returned no records')
        cache = pd.concat([cache, add]); cache.to_csv(f, sep='\t', index=False)
    m = dict(zip(cache.rsid, cache.genes.str.split(',').str[0]))
    g = genes()
    return {int(r): (m.get(int(r)) or g.get(int(r)) or '') for r in rsids}


def area_of(efo_lists):
    m, rank = _areas()
    out = []
    for terms in efo_lists:
        hits = [m[t] for t in (terms if terms is not None else []) if t in m]
        out.append(min(hits, key=rank.get) if hits else 'Unmapped')
    return out


@lru_cache(None)
def claims(split='dev'):
    with con() as c:
        d = c.sql(f"""select v.pmid, v.rsid, v.trait_key, v.pop_group, v.status, v.pub_year, v.n, v.design, v.retracted, v.pop_inferred,
                    v.link_route, v.map_route, v.outcome, v.outcome_type, v.efo, v.country, m.testable, m.replicated, m.rep_prior, m.rep_later,
                    m.first_hit_year, m.first_gwas_year, m.null_share, m.maf_eur
                    from variant_claims v join {MATCH[split]} m using (pmid, rsid, trait_key, pop_group)
                    where v.outcome_type <> 'drug_response' and v.status in ('positive','null') and m.testable""").df()
    for col in ('replicated', 'rep_prior', 'rep_later', 'pop_inferred'):
        d[col] = d[col].astype(bool)
    d['n'] = pd.to_numeric(d.n, errors='coerce')
    d['area'] = area_of(d.efo)
    d['gene'] = d.rsid.map(genes())
    d['w'] = 1.0 / d.groupby('pmid').pmid.transform('size')
    return d


@lru_cache(None)
def objects():
    with con() as c:
        return c.sql("""select pmid, rsid, trait_key, trait_id, pop_group, status, outcome_type, pub_year, country, design, n, dev
                        from variant_claims""").df()


# trait areas as drawn: the quantitative categories are one area; one hue per area in every figure (semantic table "area")
QUANT = {'Body measurement', 'Cardiovascular measurement', 'Blood measurement', 'Lipid measurement', 'Inflammatory marker',
         'Other measurement', 'Other trait'}
AREA_ORDER = ['Cancer', 'Cardiovascular disease', 'Metabolic disease', 'Neurological disease', 'Immune disease', 'Digestive disease',
              'Quantitative traits', 'Other disease', 'Not mapped']
AREA_SHORT = {'Cardiovascular disease': 'Cardiovascular', 'Metabolic disease': 'Metabolic', 'Neurological disease': 'Neurological',
              'Immune disease': 'Immune', 'Digestive disease': 'Digestive', 'Quantitative traits': 'Quantitative', 'Other disease': 'Other disease'}


def area_group(area):
    """The drawn trait area of each claim: quantitative categories merged, claims without an ontology area 'Not mapped'."""
    a = pd.Series(area)
    return np.where(a.isin(QUANT), 'Quantitative traits', np.where(a == 'Unmapped', 'Not mapped', a))


@lru_cache(None)
def claim_areas():
    """Every disease and trait claim (both splits) with its drawn trait area: pmid, rsid, pub_year, o (outcome, lower case), trait_id, grp."""
    with con() as c:
        C = c.sql("""select pmid, rsid, efo, pub_year, lower(outcome) o, trait_id from variant_claims
                     where outcome_type <> 'drug_response'""").df()
    C['grp'] = area_group(area_of(C.efo))
    assert set(C.grp) <= set(AREA_ORDER), set(C.grp) - set(AREA_ORDER)
    return C.drop(columns='efo')


@lru_cache(None)
def variant_areas():
    """(W, top): W = papers per variant and trait area (rows rsid, columns AREA_ORDER); top = for each variant the area claimed
    in most papers among the mapped areas ('Not mapped' only where no claim of the variant has an area)."""
    C = claim_areas()
    W = C.drop_duplicates(['pmid', 'rsid', 'grp']).groupby(['rsid', 'grp']).size().unstack(fill_value=0).reindex(columns=AREA_ORDER, fill_value=0)
    M = W.drop(columns='Not mapped')
    return W, M.idxmax(axis=1).where(M.sum(axis=1) > 0, 'Not mapped')


def boot_share(df, col, B=1000, seed=20261108, weight=None):
    """Share of `col` with a PMID-cluster bootstrap 95% interval."""
    w = np.ones(len(df)) if weight is None else df[weight].values
    g = pd.DataFrame({'pmid': df.pmid.values, 'w': w, 'y': w * df[col].astype(float).values}).groupby('pmid').sum()
    W, Y = g.w.values, g.y.values
    if len(g) == 0:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    m = rng.multinomial(len(g), np.full(len(g), 1 / len(g)), size=B)
    bs = (m @ Y) / (m @ W)
    return Y.sum() / W.sum(), np.quantile(bs, 0.025), np.quantile(bs, 0.975)


def check(cond, statement):
    """Titles and text state findings; when the data do not support a statement on this split, say so loudly."""
    print(('TITLE OK      ' if cond else 'TITLE CHECK FAILED  ') + statement)


@lru_cache(None)
def efo_names():
    with con() as c:
        return c.sql("select efo_id, name from efo_terms").df().set_index('efo_id').name.to_dict()


def trait_label(efo_list, fallback=''):
    """Readable trait name: the first ontology term of the claim, first letter upper case."""
    n = efo_names()
    for t in (efo_list if efo_list is not None else []):
        if t in n:
            s = n[t]; return s[0].upper() + s[1:]
    return fallback


@lru_cache(None)
def positions():
    """Genome position of every claimed variant (rsid, chr, bp, src): dbSNP GRCh38 where fetched (scripts/40_positions.py), else the
    GWAS Catalog position (GRCh38), else the position in the 1000 Genomes reference used for linkage disequilibrium (GRCh37; the
    difference is far below the resolution of a genome-wide plot)."""
    parts = []
    f = ROOT / 'work' / 'fig' / 'dbsnp_pos.tsv'
    if f.exists():
        p = pd.read_csv(f, sep='\t', dtype=str).fillna('')
        p = p[(p.chr != '') & (p.bp != '')]
        parts.append(pd.DataFrame({'rsid': p.rsid.astype('int64'), 'chr': p.chr, 'bp': p.bp.astype('int64'), 'src': 'dbsnp'}))
    with con() as c:
        g = c.sql("""select try_cast(snp_id_current as bigint) rsid, any_value(chr_id) chr, any_value(try_cast(chr_pos as bigint)) bp
                     from gwas_catalog_assoc where try_cast(chr_pos as bigint) is not null and try_cast(snp_id_current as bigint) is not null
                     group by 1""").df()
    parts.append(g.assign(src='gwas_catalog'))
    for pop in ('eur', 'eas'):
        l = pd.read_csv(ROOT / 'work' / 'ld' / f'{pop}_proxies.ld.gz', sep=r'\s+', usecols=['CHR_A', 'BP_A', 'SNP_A'], dtype=str).drop_duplicates('SNP_A')
        l = l[l.SNP_A.str.startswith('rs')]
        parts.append(pd.DataFrame({'rsid': l.SNP_A.str[2:].astype('int64'), 'chr': l.CHR_A, 'bp': l.BP_A.astype('int64'), 'src': '1000g'}))
    return pd.concat(parts).drop_duplicates('rsid', keep='first').reset_index(drop=True)


@lru_cache(None)
def variants():
    """One row per claimed variant (disease and trait claims): papers, claims and distinct traits claimed, genome position, and the
    number of distinct traits with a genome-wide significant association recorded for the same rsID in the GWAS Catalog (any trait;
    a description of the two records, not the claim-level comparison of the later figures)."""
    o = objects(); o = o[o.outcome_type != 'drug_response']
    v = o.groupby('rsid').agg(papers=('pmid', 'nunique'), claims=('pmid', 'size'), traits=('trait_id', 'nunique')).reset_index()
    with con() as c:
        g = c.sql("""select rsid, count(distinct coalesce(efo_uri, mapped_trait)) gw_traits, count(distinct pubmedid) gw_studies
                     from gwas_catalog_rsid_trait where pvalue_mlog >= 7.30103 and rsid is not null group by 1""").df()
    v = v.merge(positions(), how='left', on='rsid').merge(g, how='left', on='rsid')
    v[['gw_traits', 'gw_studies']] = v[['gw_traits', 'gw_studies']].fillna(0).astype(int)
    f = ROOT / 'work' / 'fig' / 'dbsnp_pos.tsv'                  # first gene listed by dbSNP, for labels
    if f.exists():
        p = pd.read_csv(f, sep='\t', dtype=str).fillna('')
        v = v.merge(pd.DataFrame({'rsid': p.rsid.astype('int64'), 'gene': p.gene}), how='left', on='rsid')
    else:
        v['gene'] = ''
    return v

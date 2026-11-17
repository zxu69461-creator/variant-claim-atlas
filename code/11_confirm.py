"""Confirmatory estimates (charter amendment 4.3-4.5) for one split: `python 11_confirm.py heldout` (confirmatory) or `dev` (code check
on the development split). C1 replication above the matched null (PMID cluster bootstrap); C2-C5 logistic GEE with PMID clusters;
Holm over C1-C5; secondary breakdowns; sensitivity sets (a)-(e)."""
import sys, warnings
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import numpy as np, pandas as pd, duckdb
import statsmodels.api as sm
from scipy import stats
from rdata.proj import Project
warnings.filterwarnings('ignore')
P = Project(__file__)
SPLIT = sys.argv[1] if len(sys.argv) > 1 else 'dev'; assert SPLIT in ('dev', 'heldout')
MT = 'variant_claim_match' if SPLIT == 'dev' else 'variant_claim_match_heldout'
SEED, B = 20261108, 2000

def load():
    con = duckdb.connect(str(P.root / 'db' / 'varclaims.duckdb'), read_only=True)
    d = con.sql(f"""select v.pmid, v.rsid, v.trait_key, v.pop_group, v.status, v.pub_year, v.n, v.design, v.retracted, v.pop_inferred,
                    v.link_route, v.map_route, m.testable, m.replicated, m.rep_prior, m.first_hit_year, m.null_share
                    from variant_claims v join {MT} m using (pmid, rsid, trait_key, pop_group)
                    where v.outcome_type <> 'drug_response' and v.status in ('positive','null') and m.testable""").df(); con.close()
    d['w'] = 1.0 / d.groupby('pmid').pmid.transform('size')
    for c in ('replicated', 'rep_prior'): d[c] = d[c].astype(bool)
    d['n'] = pd.to_numeric(d.n, errors='coerce'); d['pop_inferred'] = d.pop_inferred.astype(bool)
    return d

def c1(d, weighted=False):
    x = d[d.status == 'positive']; w = x.w.values if weighted else np.ones(len(x))
    g = pd.DataFrame({'pmid': x.pmid.values, 'w': w, 'r': w * x.replicated.values, 'u': w * x.null_share.fillna(0).values}).groupby('pmid').sum()
    W, R, U = g.w.values, g.r.values, g.u.values
    est = (R.sum() - U.sum()) / W.sum(); rng = np.random.default_rng(SEED); bs = []
    for k in range(0, B, 100):
        m = rng.multinomial(len(g), np.full(len(g), 1 / len(g)), size=min(100, B - k)).astype(np.float32)
        bs.append((m @ R - m @ U) / (m @ W))
    bs = np.concatenate(bs); se = bs.std(ddof=1)
    return dict(test='C1', estimand='replicated share minus matched null (positive claims)', estimate=est,
                ci_low=np.quantile(bs, .025), ci_high=np.quantile(bs, .975), p=2 * stats.norm.sf(abs(est / se)), n=len(x), n_papers=len(g),
                extra=f'replicated {x.replicated.mean() if not weighted else np.average(x.replicated, weights=w):.4f}; null {x.null_share.mean():.4f}')

def gee(x, y, xvar, scale, test, estimand, weighted=False):
    x = x.dropna(subset=[xvar]).sort_values('pmid')
    if x[y].nunique() < 2 or len(x) < 30: return dict(test=test, estimand=estimand, n=len(x), extra='too few')
    X = sm.add_constant(x[xvar].values / scale)
    mod = sm.GEE(x[y].astype(float).values, X, groups=x.pmid.values, family=sm.families.Binomial(),
                 cov_struct=sm.cov_struct.Exchangeable(), weights=x.w.values if weighted else None)
    r = mod.fit(); b, s = r.params[1], r.bse[1]
    return dict(test=test, estimand=estimand, estimate=float(np.exp(b)), ci_low=float(np.exp(b - 1.959964 * s)),
                ci_high=float(np.exp(b + 1.959964 * s)), p=float(r.pvalues[1]), n=len(x), n_papers=x.pmid.nunique(), extra=f'mean outcome {x[y].mean():.4f}')

def family(d, weighted=False):
    pos, nul = d[d.status == 'positive'], d[d.status == 'null']
    out = [c1(d, weighted)]
    t = pos[pos.n > 0].assign(lg=lambda z: np.log10(z.n))
    out.append(gee(t, 'replicated', 'lg', 1, 'C2', 'OR of replication per 10-fold stated sample size (positive claims)', weighted))
    t = pos[pos.replicated & (pos.pub_year >= 2005)]
    out.append(gee(t, 'rep_prior', 'pub_year', 5, 'C3', 'OR per 5 years that a replicated positive claim restates an earlier GWAS hit', weighted))
    t = pos[~pos.rep_prior & pos.pub_year.between(2005, 2020)].copy()
    t['y5'] = t.replicated & (t.first_hit_year - t.pub_year <= 5)
    out.append(gee(t, 'y5', 'pub_year', 5, 'C4', 'OR per 5 years of GWAS support within 5 years (new positive claims, 2005-2020)', weighted))
    t = nul[nul.pub_year >= 2005]
    out.append(gee(t, 'rep_prior', 'pub_year', 5, 'C5', 'OR per 5 years that a null claim concerns an earlier genome-wide hit for its trait', weighted))
    f = pd.DataFrame(out)
    if f.p.notna().all():
        o = np.argsort(f.p.values); adj = np.empty(len(f)); run_max = 0
        for i, j in enumerate(o): run_max = max(run_max, min(1, (len(f) - i) * f.p.values[j])); adj[j] = run_max
        f['p_holm'] = adj
    return f

def secondary(d):
    pos = d[d.status == 'positive'].copy()
    pos['size_bin'] = pd.cut(pos.n, [0, 500, 2000, 10000, np.inf], right=False, labels=['<500', '500-1,999', '2,000-9,999', '>=10,000']).astype(str)
    pos['era'] = pd.cut(pos.pub_year, [0, 2010, 2015, 2020, 3000], right=False, labels=['<2010', '2010-14', '2015-19', '2020+']).astype(str)
    pos['design_g'] = pos.design.where(pos.design.isin(['family', 'case-control', 'cohort', 'meta-analysis', 'GWAS']), 'other')
    rows = []
    for by in ('design_g', 'pop_group', 'size_bin', 'era'):
        for k, g in pos.groupby(by):
            r = g.replicated.mean(); se = np.sqrt(r * (1 - r) / len(g))
            rows.append(dict(by=by.removesuffix('_g'), level=k, n=len(g), share=len(g) / len(pos), replicated=r,
                             ci_low=max(0, r - 1.96 * se), ci_high=min(1, r + 1.96 * se), null=g.null_share.mean()))
    t = d[~d.rep_prior & d.pub_year.between(2005, 2020)].copy()
    t['y5'] = t.replicated & (t.first_hit_year - t.pub_year <= 5)
    t['era'] = pd.cut(t.pub_year, [2005, 2010, 2015, 2018, 2021], right=False, labels=['2005-09', '2010-14', '2015-17', '2018-20']).astype(str)
    for (st, e), g in t.groupby(['status', 'era']):
        rows.append(dict(by=f'new_{st}_5y', level=e, n=len(g), share=np.nan, replicated=g.y5.mean(), ci_low=np.nan, ci_high=np.nan, null=g.null_share.mean()))
    return pd.DataFrame(rows)

with P.run(__file__, params={'split': SPLIT, 'seed': SEED, 'bootstrap': B, 'charter': 'amendment 4'}) as run:
    d = load()
    prim = family(d).assign(split=SPLIT, set='primary')
    sens = []
    for name, sub, wt in [('a_no_link_d', d[d.link_route != 'd_single'], False), ('b_no_mesh_ancestor', d[d.map_route != 'mesh_ancestor'], False),
                          ('c_population_stated', d[~d.pop_inferred], False), ('d_no_retracted', d[d.retracted.fillna(0) == 0], False),
                          ('e_paper_weighted', d, True)]:
        sens.append(family(sub, wt).assign(split=SPLIT, set=name))
    sec = secondary(d).assign(split=SPLIT)
    sfx = '' if SPLIT == 'heldout' else '_dev'
    run.table('confirm_primary' + sfx, prim, note=f'{SPLIT}: amendment 4.3 primary family C1-C5, Holm')
    run.table('confirm_sensitivity' + sfx, pd.concat(sens), note=f'{SPLIT}: amendment 4.5 sensitivity sets')
    run.table('confirm_secondary' + sfx, sec, note=f'{SPLIT}: amendment 4.4 breakdowns')
    pd.set_option('display.width', 200)
    print(prim[['test', 'estimate', 'ci_low', 'ci_high', 'p', 'p_holm', 'n', 'extra']].to_string(index=False, float_format=lambda x: f'{x:.4g}'))
    print(pd.concat(sens)[['set', 'test', 'estimate', 'ci_low', 'ci_high', 'p']].to_string(index=False, float_format=lambda x: f'{x:.3g}'))

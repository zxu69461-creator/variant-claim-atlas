"""Manuscript numbers: computes every value printed in the text from the project database for one split and writes them to
manuscript/values.json (held-out, the confirmatory manuscript) or values_dev.json. Every number in the text comes from here; nothing is
typed by hand.
  ajhg: stops after the values. The AJHG manuscript (manuscript/template_ajhg.md) is assembled by the shared toolchain: `ms build`,
        `ms check`, `ms pack` with manuscript/ms.toml (references in manuscript/refs.tsv, gene symbols by scripts/ajhg_genes.py).
  base: also fills manuscript/template_en.md, the earlier Nature-style draft, into manuscript/varclaims_<split>.md (+ .docx); `dev`
        writes a preview from the development split (same text, development numbers, reader-check slots marked).
Usage: python 30_build_manuscript.py [dev|heldout] [base|ajhg]"""
import json, re, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent)); import figdata as F
SPLIT = sys.argv[1] if len(sys.argv) > 1 else 'dev'
JOURNAL = sys.argv[2] if len(sys.argv) > 2 else 'base'      # 'base' (Nature-style draft) or 'ajhg'
sfx = '' if SPLIT == 'heldout' else '_dev'
ROOT = F.ROOT; M = ROOT / 'manuscript'
V = {}


def pct(x, d=1):
    return f'{x * 100:.{d}f}%'


def fp(p):
    """P value in journal style."""
    if p < 1e-300: return 'P < 10^−300^'
    if p < 0.001:
        m, e = f'{p:.1e}'.split('e'); return f'P = {m} × 10^{int(e)}^'.replace('-', '−')
    return f'P = {p:.3f}' if p < 0.01 else f'P = {p:.2f}'


def n(x):
    return f'{int(x):,}'


with F.con() as c:
    q = lambda s: c.sql(s).df()
    o = F.objects()
    V['n_abstracts'] = n(q("select count(*) k from read_parquet('work/reader/items_all.parquet')").k[0])
    V['n_papers'] = n(o.pmid.nunique()); V['n_claims'] = n(len(o))
    nd = o[o.outcome_type != 'drug_response']
    V['n_nondrug'] = n(len(nd)); V['n_drug'] = n((o.outcome_type == 'drug_response').sum())
    V['n_rsids'] = n(nd.rsid.nunique()); V['n_traitids'] = n(nd.trait_id.nunique())
    V['n_linked'] = n(q("select count(*) k from variant_claims where outcome_type <> 'drug_response' and len(efo) > 0").k[0])
    V['pct_linked'] = pct(q("select avg((len(efo) > 0)::int) k from variant_claims where outcome_type <> 'drug_response'").k[0], 0)
    V['n_journals'] = n(o.drop_duplicates('pmid').shape[0] and q("select count(distinct journal) k from variant_claims").k[0])
    V['n_countries'] = n(o[o.country != 'unknown'].country.nunique())
    st = nd.status.value_counts(normalize=True)
    V['pct_pos_all'], V['pct_null_all'] = pct(st['positive'], 0), pct(st['null'], 0)
    yr = nd.groupby('pub_year').size()
    V['n_pre2005'] = n(nd.pub_year.lt(2005).sum())
    peak = yr.loc[2005:2024].idxmax(); V['peak_year'] = str(peak); V['peak_claims'] = n(yr[peak])
    pc = o.drop_duplicates('pmid').country.value_counts(normalize=True)
    V['china_share'] = pct(pc.get('China', 0), 0); V['us_share'] = pct(pc.get('United States', 0), 0)
    V['n_dev_papers'] = n(o[o.dev].pmid.nunique()); V['n_ho_papers'] = n(o[~o.dev].pmid.nunique())
    V['n_dev_claims'] = n((nd.dev).sum()); V['n_ho_claims'] = n((~nd.dev).sum())
    med = nd[(nd.n > 0)].n.median(); V['median_n'] = n(med)
    V['pct_with_n'] = pct((nd.n > 0).mean(), 0)
    top = nd.groupby('rsid').pmid.nunique().sort_values(ascending=False).head(3)
    gl = F.gene_label(top.index)
    V['top_variants'] = ', '.join(f'{gl[r]} rs{r} ({n(k)} papers)' for r, k in top.items())

    # analysis set on this split
    mt = F.MATCH[SPLIT]
    base = q(f"""select v.status, m.testable from variant_claims v join {mt} m using (pmid, rsid, trait_key, pop_group)
                 where v.outcome_type <> 'drug_response' and v.status in ('positive','null')""")
    V['n_split_claims'] = n(len(base)); V['pct_testable'] = pct(base.testable.mean(), 0)
    d = F.claims(SPLIT); P = d[d.status == 'positive']; N = d[d.status == 'null']
    V['n_pos'] = n(len(P)); V['n_pos_papers'] = n(P.pmid.nunique()); V['n_null'] = n(len(N)); V['n_null_papers'] = n(N.pmid.nunique())
    V['split_name'] = 'confirmatory' if SPLIT == 'heldout' else 'development'

    pr = q(f"select * from confirm_primary{sfx}").set_index('test')
    se = q(f"select * from confirm_sensitivity{sfx}")
    sec = q(f"select * from confirm_secondary{sfx}")
    r = pr.loc['C1']
    rep, nul = [float(x.split()[1]) for x in r.extra.split(';')]
    V['c1_rep'] = pct(rep); V['c1_null'] = pct(nul); V['c1_diff'] = f'{r.estimate * 100:.1f}'
    V['c1_ci'] = f'{r.ci_low * 100:.1f}–{r.ci_high * 100:.1f}'; V['c1_p'] = fp(r.p_holm); V['c1_ratio'] = f'{rep / nul:.0f}'
    V['c1_times'] = 'more than ten times' if rep / nul >= 10 else f'about {rep / nul:.0f} times'
    F.check(0.17 <= rep <= 0.25, 'text: one in five reported associations supported')
    for t in ('C2', 'C3', 'C4', 'C5'):
        r = pr.loc[t]; k = t.lower()
        V[f'{k}_or'] = f'{r.estimate:.2f}'; V[f'{k}_ci'] = f'{r.ci_low:.2f}–{r.ci_high:.2f}'; V[f'{k}_p'] = fp(r.p_holm)
        V[f'{k}_n'] = n(r.n); V[f'{k}_papers'] = n(r.n_papers)
        V[f'{k}_mean'] = pct(float(r.extra.split()[-1]))
    for t in ('C1', 'C2', 'C3', 'C4', 'C5'):
        x = se[se.test == t].estimate * (100 if t == 'C1' else 1)
        V[f'{t.lower()}_sens'] = (f'{x.min():.1f}–{x.max():.1f} percentage points' if t == 'C1' else f'{x.min():.2f}–{x.max():.2f}')
    V['sens_all_dir'] = 'in the same direction in all five sensitivity analyses'
    ok = all(((se[se.test == t].estimate > (0 if t == 'C1' else 1)) == (pr.loc[t].estimate > (0 if t == 'C1' else 1))).all() for t in pr.index)
    F.check(ok, 'text: all sensitivity estimates in the primary direction')

    def sv(by, level, col='replicated'):
        x = sec[(sec.by == by) & (sec.level == level)]; return pct(float(x[col].iloc[0]), 0)
    for lv, k in [('GWAS', 'gwas'), ('meta-analysis', 'meta'), ('cohort', 'cohort'), ('case-control', 'cc'), ('family', 'family')]:
        V[f'des_{k}'] = sv('design', lv)
    for lv, k in [('<500', 'lt500'), ('500-1,999', 'lt2k'), ('2,000-9,999', 'lt10k'), ('>=10,000', 'ge10k')]:
        V[f'size_{k}'] = sv('size_bin', lv)
    sz = sec[sec.by == 'size_bin'].set_index('level').replicated
    V['size_ratio'] = f"{sz['>=10,000'] / sz['<500']:.1f} times"
    popby = [b for b in sec.by.unique() if b.startswith('po')][0]
    for lv, k in [('EUR', 'eur'), ('EAS', 'eas'), ('other', 'oth'), ('not_stated', 'ns')]:
        V[f'pop_{k}'] = sv(popby, lv)
    for era, k in [('2005-09', '0509'), ('2010-14', '1014'), ('2015-17', '1517'), ('2018-20', '1820')]:
        V[f'y5p_{k}'] = sv('new_positive_5y', era); V[f'y5n_{k}'] = sv('new_null_5y', era)
        V[f'y5null_{k}'] = sv('new_positive_5y', era, 'null')
    # support timing and examples
    V['pct_earlier'] = pct((P.replicated & P.rep_prior).mean(), 0); V['pct_later'] = pct((P.replicated & ~P.rep_prior).mean(), 0)
    rp = P[P.replicated]
    V['share_prior_2020'] = pct(rp[rp.pub_year >= 2020].rep_prior.mean(), 0); V['share_prior_lt2010'] = pct(rp[rp.pub_year < 2010].rep_prior.mean(), 0)
    nn = N[N.pub_year >= 2005]
    V['null_prior_lt2010'] = pct(nn[nn.pub_year < 2010].rep_prior.mean(), 1); V['null_prior_2020'] = pct(nn[nn.pub_year >= 2020].rep_prior.mean(), 1)
    V['pct_null_rep'] = pct(N.replicated.mean(), 0)
    ar = []
    QUANT = {'Body measurement', 'Cardiovascular measurement', 'Blood measurement', 'Lipid measurement', 'Inflammatory marker',
             'Other measurement', 'Other trait'}
    Pg = P.assign(grp=lambda z: z.area.where(~z.area.isin(QUANT), 'Quantitative traits'))
    for g, x in Pg[Pg.grp != 'Unmapped'].groupby('grp'):
        if len(x) >= 150: ar.append((g, x.replicated.mean(), x.null_share.mean()))
    ar.sort(key=lambda z: -z[1])
    V['area_top'] = f'{ar[0][0].lower()} ({pct(ar[0][1], 0)})'; V['area_low'] = f'{ar[-1][0].lower()} ({pct(ar[-1][1], 0)})'
    V['area_null_max'] = pct(max(a[2] for a in ar), 1)
    ex = d[d.replicated].groupby(['rsid', 'trait_key']).agg(k=('pmid', 'nunique'), fh=('first_hit_year', 'min'), y0=('pub_year', 'min'),
                                                            efo=('efo', 'first')).sort_values('k', ascending=False).reset_index().drop_duplicates('rsid')
    ex = ex.head(3); gl = F.gene_label(ex.rsid)
    V['examples'] = '; '.join(f'{gl[r.rsid]} rs{r.rsid} with {F.trait_label(r.efo).lower()} ({n(r.k)} papers; first genome-wide report {int(r.fh)})'
                              for r in ex.itertuples())

    # additional values for the AJHG version
    def share(series, key): return pct(float((series == key).mean()), 0)
    V['out_disease'], V['out_trait'], V['out_sev'] = (share(nd.outcome_type, k) for k in ('disease', 'trait', 'severity'))
    for k, lab in [('case-control', 'cc'), ('cohort', 'cohort'), ('meta-analysis', 'meta'), ('GWAS', 'gwas'), ('family', 'family')]:
        V[f'des_share_{lab}'] = share(nd.design, k)
    rt = q("select link_route, map_route from variant_claims where outcome_type <> 'drug_response'")
    for k, lab in [('a_name', 'name'), ('b_tokens', 'tokens'), ('c_efo_label', 'efo'), ('d_single', 'single')]:
        V[f'link_{lab}'] = pct(float((rt.link_route == k).mean()), 1)
    for k, lab in [('efo_xref', 'xref'), ('mondo_xref', 'mondo'), ('mesh_ancestor', 'anc')]:
        V[f'map_{lab}'] = pct(float((rt.map_route == k).mean()), 1)
    V['gwas_assoc'] = n(q("select count(*) k from gwas_catalog_assoc").k[0])
    g = q("select count(*) a, count(distinct rsid) r, count(distinct pubmedid) p from gwas_catalog_rsid_trait")
    V['gwas_rsid_assoc'], V['gwas_rsids'], V['gwas_papers'] = n(g.a[0]), n(g.r[0]), n(g.p[0])
    ld = q("select pop, count(*) k from ld_proxies group by 1").set_index('pop').k
    V['ld_eur'], V['ld_eas'] = n(ld['eur']), n(ld['eas'])
    era = lambda y: pd.cut(y, [0, 2010, 2015, 2020, 3000], right=False, labels=['<2010', '2010-14', '2015-19', '2020+'])
    med = nd[nd.n > 0].groupby(era(nd[nd.n > 0].pub_year), observed=True).n.median()
    V['medn_lt2010'], V['medn_2020'] = n(med['<2010']), n(med['2020+'])
    rp2 = P[P.replicated].assign(e=lambda z: era(z.pub_year))
    rr = rp2.groupby('e', observed=True).rep_prior.mean()
    for k, lab in [('<2010', 'lt2010'), ('2010-14', '1014'), ('2015-19', '1519'), ('2020+', '2020')]: V[f'rest_{lab}'] = pct(rr[k], 0)
    nn2 = N[N.pub_year >= 2005].assign(e=lambda z: era(z.pub_year)).groupby('e', observed=True).rep_prior.mean()
    V['null_prior_1014'], V['null_prior_1519'] = pct(nn2['2010-14'], 1), pct(nn2['2015-19'], 1)
    Ps = P[P.n > 0].assign(e=lambda z: era(z.pub_year), sb=lambda z: pd.cut(z.n, [0, 500, 2000, 10000, 1e9], right=False, labels=[0, 1, 2, 3]))
    gr = Ps.pivot_table(index='e', columns='sb', values='replicated', aggfunc='mean', observed=True)
    V['grad_2020_lt500'], V['grad_2020_ge10k'] = pct(gr.loc['2020+', 0], 0), pct(gr.loc['2020+', 3], 0)
    for t, lab in [('disease', 'disease'), ('trait', 'trait'), ('severity', 'sev')]: V[f'sup_{lab}'] = pct(P[P.outcome_type == t].replicated.mean(), 0)
    mb = pd.cut(P.maf_eur, [-1, 0.01, 0.40, 0.51], right=False, labels=['low', 'mid', 'high'])
    V['maf_low'], V['maf_high'] = pct(P[mb == 'low'].replicated.mean(), 0), pct(P[mb == 'high'].replicated.mean(), 0)
    import statsmodels.api as sm
    def adj_or(x, y, scale_year=5):
        x = x.sort_values('pmid'); X = sm.add_constant(np.column_stack([x.pub_year / scale_year, np.log10(x.n)]))
        r = sm.GEE(x[y].astype(float).values, X, groups=x.pmid.values, family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()).fit()
        b, s_ = r.params[1], r.bse[1]; return f'{np.exp(b):.2f}', f'{np.exp(b - 1.959964 * s_):.2f}–{np.exp(b + 1.959964 * s_):.2f}'
    t4 = P[~P.rep_prior & P.pub_year.between(2005, 2020) & (P.n > 0)].assign(y5=lambda z: z.replicated & (z.first_hit_year - z.pub_year <= 5))
    V['c4adj_or'], V['c4adj_ci'] = adj_or(t4, 'y5')
    t3 = P[P.replicated & (P.pub_year >= 2005) & (P.n > 0)]
    V['c3adj_or'], V['c3adj_ci'] = adj_or(t3, 'rep_prior')
    rcr = q("select * from reader_check_result").set_index('measure') if SPLIT == 'heldout' else None
    if rcr is not None:
        for m_, k in [('pair_error', 'rc_err'), ('no_conclusion', 'rc_nc'), ('recall', 'rc_recall'), ('recall_specificity', 'rc_recspec')]:
            x = rcr.loc[m_]; V[k] = f'{x.estimate * 100:.0f}% ({int(x.k)} of {int(x.n)})'
    # Table 1: composition and support by era (composition: all papers; support: this split)
    E = nd.assign(e=era(nd.pub_year))
    comp = E.groupby('e', observed=True).agg(papers=('pmid', 'nunique'), claims=('pmid', 'size'), eas=('pop_group', lambda s: (s == 'EAS').mean()),
                                              medn=('n', lambda s: s[s > 0].median()), pos=('status', lambda s: (s == 'positive').mean()))
    sp = P.assign(e=era(P.pub_year)).groupby('e', observed=True).agg(k=('replicated', 'size'), rep=('replicated', 'mean'), nul=('null_share', 'mean'), prior=('rep_prior', 'mean'))
    lab = {'<2010': 'Before 2010', '2010-14': '2010–2014', '2015-19': '2015–2019', '2020+': '2020 onward'}
    t1 = ['| Era | Papers | Claims | East Asian samples | Median stated sample size | Reporting an association | Testable reported associations (held-out) | With genome-wide support | Matched null | Support earlier than the claim |',
          '|---|---|---|---|---|---|---|---|---|---|']
    for e_ in comp.index:
        c_, s_ = comp.loc[e_], sp.loc[e_]
        t1.append(f'| {lab[e_]} | {n(c_.papers)} | {n(c_.claims)} | {pct(c_.eas, 0)} | {n(c_.medn)} | {pct(c_.pos, 0)} | {n(s_.k)} | {pct(s_.rep, 1)} | {pct(s_.nul, 1)} | {pct(s_.prior, 1)} |')
    V['n_nondrug_papers'] = n(nd.pmid.nunique())
    assert int(comp.papers.sum()) == nd.pmid.nunique() and int(comp.claims.sum()) == len(nd)
    V['TABLE1'] = '\n'.join(t1) + (f"\n\nPapers, claims, populations and sample sizes are from all {V['n_nondrug_papers']} papers with disease or trait claims "
                                   f"({V['n_nondrug']} claims; claims on drug response are not included); the last four columns are from held-out papers. "
                                   'Shares of claims; support and matched null among testable reported associations.')
    for k_, lab_ in [('<2010', 'lt2010'), ('2010-14', '1014'), ('2015-19', '1519'), ('2020+', '2020')]: V[f'eas_{lab_}'] = pct(comp.eas[k_], 0)
    assert comp.eas[['2010-14', '2015-19', '2020+']].min() > 0.25          # text: more than a quarter of claims in every era since 2010
    # Table 2: most frequently claimed variants among held-out testable reported associations
    T2 = P.groupby('rsid').agg(papers=('pmid', 'nunique'), traits=('trait_id', 'nunique') if 'trait_id' in P else ('trait_key', 'nunique'),
                               rep=('replicated', 'mean'), nul=('null_share', 'mean'), fh=('first_hit_year', 'min')).sort_values('papers', ascending=False).head(12)
    gl2 = F.gene_label(T2.index)
    t2 = ['| Variant | Gene | Papers | Distinct traits | With genome-wide support | Matched null | First genome-wide report for a claimed trait |', '|---|---|---|---|---|---|---|']
    for r_, x in T2.iterrows():
        t2.append(f"| rs{r_} | *{gl2[r_]}* | {n(x.papers)} | {n(x.traits)} | {pct(x.rep, 0)} | {pct(x.nul, 1)} | {'—' if pd.isna(x.fh) else int(x.fh)} |")
    V['TABLE2'] = '\n'.join(t2) + '\n\nHeld-out testable reported associations; gene from dbSNP; support and matched null as in Figure 2; the last column is the earliest genome-wide significant report among the traits claimed for the variant (—, none).'
    # Table 2 text: two groups among the most-claimed variants; development estimates for comparison
    hi = T2[T2.rep >= 0.5]; lo = T2[T2.rep <= 0.10]
    fmt = lambda ix: ', '.join(f'*{gl2[r_]}* rs{r_}' for r_ in ix[:-1]) + (' and ' if len(ix) > 1 else '') + f'*{gl2[ix[-1]]}* rs{ix[-1]}'
    V['t2_n'] = str(len(T2)); V['t2_high'] = fmt(list(hi.index)); V['t2_low'] = fmt(list(lo.index))
    V['t2_high_range'] = f'{hi.rep.min() * 100:.0f}–{hi.rep.max() * 100:.0f}%'; V['t2_low_range'] = f'{lo.rep.min() * 100:.0f}–{lo.rep.max() * 100:.0f}%'
    V['t2_high_years'] = f'{int(hi.fh.min())} and {int(hi.fh.max())}'; V['t2_low_traits'] = f'{int(lo.traits.min())} and {int(lo.traits.max())}'
    F.check(len(hi) >= 2 and len(lo) >= 2 and len(hi) + len(lo) >= len(T2) - 2, 'text: most-claimed variants separate into two groups')
    dv = q("select * from confirm_primary_dev").set_index('test')
    V['dev_c1'] = f"{dv.loc['C1'].estimate * 100:.1f}"
    for t_ in ('C2', 'C3', 'C4', 'C5'): V[f'dev_{t_.lower()}'] = f"{dv.loc[t_].estimate:.2f}"
    # reader check
    rc = q("select * from reader_check_result").set_index('measure') if 'reader_check_result' in set(q('show tables').name) else None
    if rc is not None and SPLIT == 'heldout':
        for m, k in [('positive_precision', 'rc_pos'), ('null_precision', 'rc_null'), ('pair_error', 'rc_err'), ('recall', 'rc_recall')]:
            x = rc.loc[m]; V[k] = f'{x.estimate * 100:.0f}% (95% CI {x.ci_low * 100:.0f}–{x.ci_high * 100:.0f}%; {int(x.k)} of {int(x.n)})'

    # ---- numbers quoted by the legends and text of the figures redrawn on 2026-11-10 (scripts 41-45); all papers unless stated
    from scipy.stats import spearmanr
    VV = F.variants(); CA = F.claim_areas()
    V['fig1_rho'] = f"{spearmanr(VV.papers, VV.gw_traits)[0]:.2f}"
    V['n_var_gw'] = n((VV.gw_traits > 0).sum()); V['pct_var_gw'] = pct((VV.gw_traits > 0).mean(), 0)
    V['n_var_placed'] = n((VV.chr.astype(str).isin([str(i) for i in range(1, 23)] + ['X']) & VV.bp.notna()).sum())
    cl = np.sort(VV.claims.values)[::-1]; V['top1_share'] = pct(cl[:int(0.01 * len(cl))].sum() / cl.sum(), 0)
    ac = CA.grp.value_counts()
    V['pct_area'] = pct(1 - ac['Not mapped'] / len(CA), 0)
    named = sorted([a for a in F.AREA_ORDER if a not in ('Not mapped', 'Other disease')], key=lambda a: -ac[a])
    lab = [f"{a.lower()} ({n(ac[a])})" for a in named]
    V['area_counts'] = ', '.join(lab[:-1]) + ' and ' + lab[-1]
    fto = VV[VV.rsid == 9939609].iloc[0]
    assert fto.papers == VV.papers.max()
    V['fto_papers'], V['fto_traits'], V['fto_gw'] = n(fto.papers), n(fto.traits), n(fto.gw_traits)
    sp = q("""select v.pub_year, m.first_hit_year, m.n_hit_papers from variant_claims v join variant_claim_match_heldout m
              using (pmid, rsid, trait_key, pop_group) where v.pmid = 18426866 and v.rsid = 9939609 and v.outcome = 'BMI'""").iloc[0]
    V['fto_bmi_studies'], V['fto_bmi_first'] = n(sp.n_hit_papers), str(int(sp.first_hit_year))
    fy = CA[CA.rsid == 9939609].drop_duplicates('pmid').groupby('pub_year').size().reindex(range(2008, 2026), fill_value=0)
    V['fto_min_per_year'] = n(fy.min())                     # papers per year claiming rs9939609, 2008-2025
    # Figure 2: timing shares to one decimal (they add to the supported share); trait area by allele-frequency class
    V['pct_earlier1'] = pct((P.replicated & P.rep_prior).mean(), 1); V['pct_later1'] = pct((P.replicated & ~P.rep_prior).mean(), 1)
    Pm = P.assign(grp=F.area_group(P.area), bin=pd.cut(P.maf_eur, [-1, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.51], right=False).astype(str))
    gm = Pm[Pm.grp != 'Not mapped'].groupby(['grp', 'bin']).agg(s=('replicated', 'mean'), nu=('null_share', 'mean'), k=('replicated', 'size'))
    gm = gm[gm.k >= 30]
    assert (gm.s > gm.nu).all()
    V['n_area_maf_cells'] = str(len(gm)); V['area_maf_min'] = pct(gm.s.min(), 0); V['area_maf_max'] = pct(gm.s.max(), 0)
    # Figure 3: the size gradient within trait areas, designs and populations (classes with 30 or more claims)
    Ps = P[P.n > 0].assign(grp=F.area_group(P[P.n > 0].area), sb=pd.cut(P[P.n > 0].n, [0, 500, 2000, 10000, 1e9], right=False, labels=[0, 1, 2, 3]).astype(int))
    for col, keep in (('grp', [a for a in F.AREA_ORDER if a != 'Not mapped']), ('design', ['family', 'case-control', 'cohort', 'meta-analysis', 'GWAS']),
                      ('pop_group', ['EUR', 'other', 'EAS', 'not_stated'])):
        gs_ = Ps[Ps[col].isin(keep)].groupby([col, 'sb']).replicated.agg(['mean', 'size'])
        gs_ = gs_[gs_['size'] >= 30]['mean'].unstack()
        assert all(r.dropna().iloc[-1] > r.dropna().iloc[0] for _, r in gs_.iterrows()), col
    # Figure 4: restatement by three-year period of publication
    rp4 = P[P.replicated]
    V['rest_0507'] = pct(rp4[rp4.pub_year.between(2005, 2007)].rep_prior.mean(), 0)
    V['rest_2325'] = pct(rp4[rp4.pub_year.between(2023, 2025)].rep_prior.mean(), 0)

    # ---- values added with the pre-submission revision of 2026-11-10
    V['pct_pos_with_n'] = pct((P.n > 0).mean(), 0)
    assert int((P.n > 0).sum()) == int(V['c2_n'].replace(',', ''))
    V['pct_null_earlier'] = pct(N.rep_prior.mean(), 1); V['pct_null_later'] = pct((N.replicated & ~N.rep_prior).mean(), 1)
    V['pct_allele'] = pct(q("select avg((risk_allele in ('A','C','G','T'))::int) k from variant_claims where outcome_type <> 'drug_response'").k[0], 0)
    if SPLIT == 'heldout':
        ra = q("select reader_status, human, count(*) k from reader_check_answers where part = 'pair' group by all")
        rk = lambda s_, h_: int(ra[(ra.reader_status == s_) & (ra.human == h_)].k.sum())
        V['rc_null_rest'], V['rc_null_nc'], V['rc_null_yes'] = str(rk('null', 'none') + rk('null', 'yes') + rk('null', 'notstudied')), str(rk('null', 'none')), str(rk('null', 'yes'))
        assert rk('null', 'no') + int(V['rc_null_rest']) == 40
        # the design label against the sample: associations reported by genome-wide studies, by stated sample size
        Pg_ = P[(P.n > 0) & (P.design == 'GWAS')]
        gg = Pg_.groupby(pd.cut(Pg_.n, [0, 500, 2000, 10000, 1e9], right=False, labels=[0, 1, 2, 3]), observed=True).replicated.agg(['mean', 'size'])
        assert gg['size'].min() >= 30
        V['gwas_lt500'], V['gwas_ge10k'] = pct(gg.loc[0, 'mean'], 0), pct(gg.loc[3, 'mean'], 0)
        # additional analysis: the primary estimates without claims from genome-wide studies and meta-analyses (scripts/13_additional.py)
        xs = q("select * from additional_sensitivity").set_index('test')
        r = xs.loc['C1']; rep_, nul_ = [float(x.split()[1]) for x in r.extra.split(';')]
        V['x_c1_rep'], V['x_c1_null'] = pct(rep_), pct(nul_)
        V['x_c1_diff'], V['x_c1_ci'] = f'{r.estimate * 100:.1f}', f'{r.ci_low * 100:.1f}–{r.ci_high * 100:.1f}'
        for t in ('C2', 'C3', 'C4', 'C5'):
            r = xs.loc[t]; V[f'x_{t.lower()}_or'] = f'{r.estimate:.2f}'; V[f'x_{t.lower()}_ci'] = f'{r.ci_low:.2f}–{r.ci_high:.2f}'
        for t in xs.index:                                          # same side as the primary estimate, interval excluding no effect
            ref = 0 if t == 'C1' else 1
            assert (xs.loc[t].ci_low > ref) if pr.loc[t].estimate > ref else (xs.loc[t].ci_high < ref), t
        V['x_share'] = pct(d.design.isin(['GWAS', 'meta-analysis']).mean(), 0)
    # software versions, read from the environment that ran the analysis
    import platform, importlib.metadata as md
    V['ver_python'] = platform.python_version()
    for pk in ('pandas', 'numpy', 'scipy', 'statsmodels', 'matplotlib', 'pyvista'): V[f'ver_{pk}'] = md.version(pk)

# values for the manuscript toolchain: `ms build` fills manuscript/template_ajhg.md from this file (manuscript/ms.toml)
vf = M / ('values.json' if SPLIT == 'heldout' else f'values_{SPLIT}.json')
vf.write_text(json.dumps(V, ensure_ascii=False, indent=1, sort_keys=True) + '\n')
if JOURNAL == 'ajhg':   # since 2026-11-09 the AJHG version is assembled by ms (references from refs.tsv, gene symbols by scripts/ajhg_genes.py)
    print('written', vf, f'({len(V)} values) | build the AJHG manuscript with: ms build --lang all   (in {ROOT})')
    sys.exit(0)


# references (PubMed-verified 2026-11-09); numbered by first citation in the text
REFS = {
 'ioannidis2001': 'Ioannidis, J. P., Ntzani, E. E., Trikalinos, T. A. & Contopoulos-Ioannidis, D. G. Replication validity of genetic association studies. *Nat. Genet.* **29**, 306–309 (2001).',
 'hirschhorn2002': 'Hirschhorn, J. N., Lohmueller, K., Byrne, E. & Hirschhorn, K. A comprehensive review of genetic association studies. *Genet. Med.* **4**, 45–61 (2002).',
 'colhoun2003': 'Colhoun, H. M., McKeigue, P. M. & Davey Smith, G. Problems of reporting genetic associations with complex outcomes. *Lancet* **361**, 865–872 (2003).',
 'ioannidis2005': 'Ioannidis, J. P. A. Why most published research findings are false. *PLoS Med.* **2**, e124 (2005).',
 'ioannidis2006': 'Ioannidis, J. P., Trikalinos, T. A. & Khoury, M. J. Implications of small effect sizes of individual genetic variants on the design and interpretation of genetic association studies of complex diseases. *Am. J. Epidemiol.* **164**, 609–614 (2006).',
 'visscher2017': 'Visscher, P. M. et al. 10 years of GWAS discovery: biology, function, and translation. *Am. J. Hum. Genet.* **101**, 5–22 (2017).',
 'abdellaoui2023': 'Abdellaoui, A., Yengo, L., Verweij, K. J. H. & Visscher, P. M. 15 years of GWAS discovery: realizing the promise. *Am. J. Hum. Genet.* **110**, 179–194 (2023).',
 'siontis2010': 'Siontis, K. C., Patsopoulos, N. A. & Ioannidis, J. P. Replication of past candidate loci for common diseases and phenotypes in 100 genome-wide association studies. *Eur. J. Hum. Genet.* **18**, 832–837 (2010).',
 'border2019': 'Border, R. et al. No support for historical candidate gene or candidate gene-by-interaction hypotheses for major depression across multiple large samples. *Am. J. Psychiatry* **176**, 376–387 (2019).',
 'duncan2011': 'Duncan, L. E. & Keller, M. C. A critical review of the first 10 years of candidate gene-by-environment interaction research in psychiatry. *Am. J. Psychiatry* **168**, 1041–1049 (2011).',
 'grant2006': 'Grant, S. F. et al. Variant of transcription factor 7-like 2 (TCF7L2) gene confers risk of type 2 diabetes. *Nat. Genet.* **38**, 320–323 (2006).',
 'klein2005': 'Klein, R. J. et al. Complement factor H polymorphism in age-related macular degeneration. *Science* **308**, 385–389 (2005).',
 'romeo2008': 'Romeo, S. et al. Genetic variation in PNPLA3 confers susceptibility to nonalcoholic fatty liver disease. *Nat. Genet.* **40**, 1461–1465 (2008).',
 'sollis2023': 'Sollis, E. et al. The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource. *Nucleic Acids Res.* **51**, D977–D985 (2023).',
 'malone2010': 'Malone, J. et al. Modeling sample variables with an Experimental Factor Ontology. *Bioinformatics* **26**, 1112–1118 (2010).',
 'wei2024': 'Wei, C.-H. et al. PubTator 3.0: an AI-powered literature resource for unlocking biomedical knowledge. *Nucleic Acids Res.* **52**, W540–W546 (2024).',
 '1000g2015': '1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature* **526**, 68–74 (2015).',
 'chang2015': 'Chang, C. C. et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience* **4**, 7 (2015).',
}


def cite(text):
    order = []
    def rep(m):
        keys = [k.strip().lstrip('@') for k in m.group(1).split(';')]
        for k in keys:
            if k not in REFS: raise KeyError(f'unknown reference {k}')
            if k not in order: order.append(k)
        nums = sorted(order.index(k) + 1 for k in keys)
        return '^' + ','.join(map(str, nums)) + '^'
    body = re.sub(r'\[(@[^\]]+)\]', rep, text)
    return body.replace('{{REFERENCES}}', '\n'.join(f'{i}. {REFS[k]}' for i, k in enumerate(order, 1)))

tpl = cite((M / 'template_en.md').read_text())
miss = sorted(set(re.findall(r'\{\{(\w+)\}\}', tpl)) - set(V))
out = re.sub(r'\{\{(\w+)\}\}', lambda m: V.get(m.group(1), f'⟦{m.group(1)}⟧'), tpl)
dst = M / f'varclaims_{SPLIT}.md'
if SPLIT != 'heldout':
    out = ('> Development-split preview: numbers are from the development papers and the reading check is pending; the confirmatory '
           'manuscript is built from the held-out papers.\n\n') + out
dst.write_text(out)
print('written', dst, '| unfilled:', ', '.join(miss) or 'none')
subprocess.run(['pandoc', str(dst), '-o', str(dst.with_suffix('.docx'))], check=True)

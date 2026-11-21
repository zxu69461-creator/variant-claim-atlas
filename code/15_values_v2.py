"""Numbers for the v2 manuscript (2026-11-17, figures/final_2026-11-17), added to manuscript/values.json after
30_build_manuscript.py (run that first, then this). Same definitions as the figure scripts scripts/v2/fig1_v3.py, fig3_v3.py and
fig4_v3.py and the additional ClinVar analysis scripts/14_clinvar.py."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); import figdata as F
VJ = F.ROOT / 'manuscript' / 'values.json'
V = json.loads(VJ.read_text())
pct = lambda x, d=0: f'{x * 100:.{d}f}%'          # noqa: E731
num = lambda x: f'{int(x):,}'                     # noqa: E731

# ---- ClinVar groups (additional analysis)
with F.con() as c:
    G = c.sql('select * from additional_clinvar_groups').df().set_index('group')
    CV = c.sql('select * from additional_clinvar_variants').df()
g = G.loc['ClinVar risk variants']; o = G.loc['Other variants']
V.update(cv_vars=num(g.variants), cv_claims=num(g.claims), cv_papers=num(g.papers), cv_sup=pct(g.supported, 1), cv_null=pct(g.null, 1),
         cv_excess=f'{g.excess * 100:.1f}', cv_ci=f'{g.ci_low * 100:.1f}–{g.ci_high * 100:.1f}', oth_sup=pct(o.supported, 1), oth_null=pct(o.null, 1))
c10 = CV[CV.claims >= 10]
V.update(cv_ge10=num(len(c10)), cv_ge10_mostly=num((c10.supported >= 0.5).sum()), cv_ge10_rarely=num((c10.supported < 0.05).sum()))

# ---- Figure 4 classes (held-out, variants with >= 5 reported associations), as in fig4_v3.py
cv = pd.read_csv(F.ROOT / 'raw' / 'variant_summary.txt.gz', sep='\t', usecols=['ClinicalSignificance', 'RS# (dbSNP)', 'Assembly', 'ReviewStatus'], dtype=str)
cv = cv[cv.Assembly == 'GRCh38'].copy(); cv['rsid'] = pd.to_numeric(cv['RS# (dbSNP)'], errors='coerce'); cv = cv[cv.rsid > 0]
sig = cv.ClinicalSignificance.fillna('').str.lower().str.replace('association not found', '', regex=False)
cv['rf'], cv['as_'], cv['pr'] = sig.str.contains('risk factor'), sig.str.contains('association'), sig.str.contains('protective')
cv = cv[cv.rf | cv.as_ | cv.pr]
RANK = {'no assertion criteria provided': 0, 'criteria provided, single submitter': 1, 'criteria provided, conflicting classifications': 1,
        'criteria provided, multiple submitters, no conflicts': 2, 'reviewed by expert panel': 2, 'practice guideline': 2}
cv['rv'] = cv.ReviewStatus.map(RANK).fillna(0)
R = cv.groupby('rsid').agg(rf=('rf', 'any'), as_=('as_', 'any'), rv=('rv', 'max')).reset_index(); R['rsid'] = R.rsid.astype('int64')
R['assert'] = np.where(R.rf, 'Risk factor', np.where(R.as_, 'Association', 'Protective'))


def split(name):
    d = F.claims(name); pos = d[d.status == 'positive']
    Vv = pos.groupby('rsid').agg(n=('pmid', 'size'), sup=('replicated', 'mean')).reset_index().merge(R, on='rsid', how='left')
    Vv['cv'] = Vv['assert'].notna(); Vv = Vv[Vv.n >= 5].copy()
    Vv['cls'] = np.where(Vv.sup >= 0.5, 'mostly', np.where(Vv.sup >= 0.05, 'partly', 'rarely'))
    return pos, Vv


posH, VH = split('heldout'); posD, VD = split('dev')
rk = VH[VH.cv]; ot = VH[~VH.cv]; rf = VH[VH['assert'] == 'Risk factor']
V.update(f4_n_risk=num(len(rk)), f4_n_rf=num(len(rf)), f4_n_other=num(len(ot)),
         f4_rf_mostly=pct((rf.cls == 'mostly').mean()), f4_rf_rarely=pct((rf.cls == 'rarely').mean()),
         f4_risk_mostly=pct((rk.cls == 'mostly').mean()), f4_risk_rarely=pct((rk.cls == 'rarely').mean()),
         f4_oth_mostly=pct((ot.cls == 'mostly').mean()), f4_oth_rarely=pct((ot.cls == 'rarely').mean()),
         f4_rf_claims_mostly=pct(rf.n[rf.cls == 'mostly'].sum() / rf.n.sum()))
med = [rk[rk.rv == k].sup.median() for k in range(3)]; nn = [int((rk.rv == k).sum()) for k in range(3)]
V.update(f4_rv_med0=pct(med[0]), f4_rv_med1=pct(med[1]), f4_rv_med2=pct(med[2]), f4_rv_n0=str(nn[0]), f4_rv_n1=str(nn[1]), f4_rv_n2=str(nn[2]),
         f4_rv_p=f'P = {mannwhitneyu(rk[rk.rv == 0].sup, rk[rk.rv == 2].sup).pvalue:.2g}')
V.update(f4_e_p_ho=f'{mannwhitneyu(rk.sup, ot.sup).pvalue:.1e}'.replace('e-0', ' × 10^−').replace('e-', ' × 10^−') + '^',
         f4_e_p_dev=f'{mannwhitneyu(VD[VD.cv].sup, VD[~VD.cv].sup).pvalue:.2f}', f4_n_risk_dev=num(VD.cv.sum()))
rows = posH.merge(rk[['rsid', 'cls']], on='rsid')
pk = {c: int(rows[rows.cls == c].groupby('pub_year').size().reindex(range(2003, 2026), fill_value=0).rolling(3, center=True).mean().idxmax()) for c in ('mostly', 'rarely')}
V.update(f4_peak_mostly=str(pk['mostly']), f4_peak_rarely=str(pk['rarely']))

# ---- Figure 1 (both splits, descriptive), as in fig1_v3.py
Vt = F.variants(); dd = pd.concat([F.claims('dev'), F.claims('heldout')]); p = dd[dd.status == 'positive']
S = p.groupby('rsid').agg(n_t=('pmid', 'size'), sup=('replicated', 'mean')).reset_index(); Vt = Vt.merge(S, on='rsid', how='left')
x = Vt[Vt.n_t >= 10]
V.update(f1_hl=num(((Vt.n_t >= 5) & (Vt.sup > 0.5)).sum()), f1_n10=num(len(x)), f1_lt5=num((x.sup < 0.05).sum()), f1_gt80=num((x.sup > 0.8).sum()),
         f1_lt5_pct=pct((x.sup < 0.05).mean()), f1_gt80_pct=pct((x.sup > 0.8).mean()))

# ---- Figure 3B (held-out pairs), as in fig3_v3.py
d = F.claims('heldout')
pairs = (d[d.replicated].groupby(['rsid', 'trait_key']).agg(n=('pmid', 'nunique'), fh=('first_hit_year', 'min')).sort_values('n', ascending=False)
         .reset_index().drop_duplicates('rsid').head(22))
nb = na = 0
for pr in pairs.itertuples():
    q = d[(d.rsid == pr.rsid) & (d.trait_key == pr.trait_key)].drop_duplicates('pmid'); b = (q.pub_year < pr.fh).sum(); nb += b; na += len(q) - b
V.update(f3_pairs=str(len(pairs)), f3_after=num(na), f3_before=num(nb))

# ---- additional analyses of 2026-11-17 (scripts/16_robustness.py): narrower trait rules, ClinVar conditions, allele agreement
with F.con() as c:
    TR = c.sql('select * from additional_trait_rule').df()
    CC = c.sql('select * from additional_clinvar_condition').df()
    AA = c.sql('select * from additional_allele_agreement').df().set_index('group')
for rule in ('exact', 'one_step'):
    t = TR[TR.set == rule].set_index('test')
    V[f'tr_{rule}_sup'] = pct(t.support.iloc[0], 1); V[f'tr_{rule}_null'] = pct(t.null.iloc[0], 1)
    V[f'tr_{rule}_c1'] = f"{t.loc['C1'].estimate * 100:.1f}"
    V[f'tr_{rule}_c1ci'] = f"{t.loc['C1'].ci_low * 100:.1f}–{t.loc['C1'].ci_high * 100:.1f}"
    for k in ('C2', 'C3', 'C4', 'C5'):
        V[f'tr_{rule}_{k.lower()}'] = f'{t.loc[k].estimate:.2f}'
        V[f'tr_{rule}_{k.lower()}ci'] = f'{t.loc[k].ci_low:.2f}–{t.loc[k].ci_high:.2f}'
r_ = CC[CC.cls == 'rarely']; m_ = CC[CC.cls == 'mostly']
t_ = CC[CC.cond_term]
V.update(cc_all=num(len(CC)), cc_all_term=num(len(t_)), cc_all_condsup=num(t_.cond_support.sum()), cc_all_condsup_pct=pct(t_.cond_support.mean()),
         cc_rare_noterm=num(len(r_) - r_.cond_term.sum()))
V.update(cc_rare=num(len(r_)), cc_rare_nogw=num((~r_.gwas_any).sum()), cc_rare_term=num(r_.cond_term.sum()), cc_rare_condsup=num(r_.cond_support.sum()),
         cc_most=num(len(m_)), cc_most_term=num(m_.cond_term.sum()), cc_most_condsup=num(m_.cond_support.sum()),
         cc_most_condsup_pct=pct(m_.cond_support.sum() / m_.cond_term.sum()))
a = AA.loc['Disease risk or quantitative trait']
from scipy.stats import binomtest
V.update(al_claims=num(a.claims), al_share=pct(a.share, 1), al_ci=f'{a.ci_low * 100:.1f}–{a.ci_high * 100:.1f}',
         al_p=f'{binomtest(int(a.agree), int(a.claims), 0.5).pvalue:.0e}'.replace('e-', ' × 10^−') + '^')

# ---- Figure 2D median and the pyCirclize version (previously added by hand; now computed here)
_p = F.claims('heldout'); _p = _p[(_p.status == 'positive') & (_p.n > 0)]
from importlib.metadata import version as _ver
V.update(f2_median=f'{int(_p.n.median()):,}', ver_pycirclize=_ver('pycirclize'))

# ---- the web tool's claim checks shown in Figure 5 (lookup_site/data/claims, both splits)
def _chk(rs, t):
    sh = json.loads((F.ROOT / 'lookup_site' / 'data' / 'claims' / f'{rs % 64:02d}.json').read_text())
    C = [dict(zip(sh['cols'], r)) for r in sh['rows'] if r[1] == rs and t in r[2].lower()]
    R = [c for c in C if c['status'] == 'reported']; Tt = [c for c in R if c['testable']]; S = [c for c in Tt if c['supported']]
    return len(R), len(Tt), len(S), len(C) - len(R), min([c['first_gw_year'] for c in S if c['first_gw_year']] or [0])
r_ = _chk(1800795, 'type 2 diabetes'); t_ = _chk(7903146, 'type 2 diabetes')
V.update(chk_il6_rep=str(r_[0]), chk_il6_test=str(r_[1]), chk_il6_sup=str(r_[2]), chk_il6_abs=str(r_[3]),
         chk_tcf_rep=str(t_[0]), chk_tcf_test=str(t_[1]), chk_tcf_sup=str(t_[2]), chk_tcf_first=str(t_[4]))

V['rc_pos_short'] = V['rc_pos'].split(' ')[0]                     # reading-check precision alone, for the abstract

# ---- language-model outputs that reached the token limit or could not be parsed (work/reader/reader_all_v2.jsonl)
_n = _t = _pf = _ca = _ct = 0
for _l in open(F.ROOT / 'work' / 'reader' / 'reader_all_v2.jsonl'):
    _r = json.loads(_l); _n += 1
    if _r['claims'] is None:
        _pf += 1; continue
    _ca += len(_r['claims'])
    if _r.get('truncated'):
        _t += 1; _ct += len(_r['claims'])
V.update(trunc_n=num(_t), trunc_pct=pct(_t / _n, 1), trunc_entries_pct=pct(_ct / _ca, 1), parse_fail=num(_pf), parse_fail_pct=pct(_pf / _n, 2))

# ---- the lookup page's record for rs1800795 (lookup_site/data/variants.json, both splits)
_lv = json.loads((F.ROOT / 'lookup_site' / 'data' / 'variants.json').read_text())
_row = dict(zip(_lv['cols'], next(r for r in _lv['rows'] if r[0] == 1800795)))
V.update(lk_il6_papers=num(_row['papers']), lk_il6_reported=num(_row['reported']), lk_il6_testable=num(_row['testable']),
         lk_il6_sup=pct(_row['supported'], 1), lk_il6_null=pct(_row['null'], 1))

V['TABLE2'] = V['TABLE2'].replace('support and matched null as in Figure 2', 'support and matched null as in Figure 1')
VJ.write_text(json.dumps(V, ensure_ascii=False, indent=1))
print({k: V[k] for k in V if k.startswith(('cv_', 'f4_', 'f1_', 'f3_', 'tr_', 'cc_', 'al_'))})

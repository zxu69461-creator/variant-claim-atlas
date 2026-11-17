"""Figure 1 v3 (2026-11-17): rebuilt at the proportions and density of the user's picked references after the user judged v2 not a
reproduction. No panel titles (statements go to the legend). Descriptive panels over all papers of both splits.
a  highlighted Manhattan skyline (PMC12900644_F1): papers claiming each variant; chromosomes alternate maroon / blue; variants whose
   reported associations are mostly supported genome-wide (>= 5 testable, > 50%) in green; named variants as diamonds / open circles.
b  every variant with >= 10 reported associations, square scatter with marginal distributions (PMC11090817_F2 a).
c  six named variants (three mostly supported, three almost never) x eight trait areas (rows shared across facets, area colour bar on the
   left as the reference's pathway bar) x publication year; dot size = papers in that area and year, filled green when most of those
   associations are supported, tan otherwise (PMC12900648_F7 d).
d  claims per year by trait area; e the most-claimed variant (FTO rs9939609) through time (kept from the submission, restyled).
Output figures/v3/fig1."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
from matplotlib.collections import LineCollection as LC
import figkit as fk
sys.path.insert(0, str(Path(__file__).resolve().parents[1])); import figdata as F
sys.path.insert(0, str(Path(__file__).resolve().parent)); import roles as RL
sys.path.insert(0, '/Volumes/ResearchData/Claude/figtaste')
from layout.canvas import figure, axes, ftext

MAR, NAV = RL.G2, RL.G3                                  # PMC12900644_F1 chromosome alternation
GRN, GRE = RL.SUP, RL.dark(RL.SUP)                                  # PMC12900644_F1 highlight / outline
TAN, TAND = RL.NOT, RL.NOT_D                                 # PMC11090817_F2 second colour
INK, GD, GL, GRID = '#1a1a1a', '#555555', '#c8c8c8', '#e8e8e8'
FS, FT, FL = 6.5, 6.0, 9.5
e = fk.enc('varclaims'); AREAS = F.AREA_ORDER; C = {a: e.color('area', a) for a in AREAS}
OUT = RL.OUT

V = F.variants()
dd = pd.concat([F.claims('dev'), F.claims('heldout')]); pos = dd[dd.status == 'positive'].copy(); pos['grp'] = F.area_group(pos.area)
S = pos.groupby('rsid').agg(n_t=('pmid', 'size'), sup=('replicated', 'mean')).reset_index()
V = V.merge(S, on='rsid', how='left'); V['hl'] = (V.n_t >= 5) & (V.sup > 0.5)
CL = F.claim_areas()
GENE = {9939609: 'FTO', 738409: 'PNPLA3', 7903146: 'TCF7L2', 429358: 'APOE', 671: 'ALDH2', 12979860: 'IFNL4',
        1801133: 'MTHFR', 1800629: 'TNF', 1800795: 'IL6', 2228570: 'VDR', 2910164: 'MIR146A', 4680: 'COMT', 6265: 'BDNF'}
NUP = [9939609, 738409, 7903146, 429358, 671, 12979860]; NDN = [1801133, 1800629, 1800795, 2228570, 2910164, 4680, 6265]
fig = figure(180, 176)


def letter(x, y, s):
    ftext(fig, x, y, s, fontsize=FL, fontweight='bold', va='top')


def style(ax):
    for s_ in ('top', 'right'): ax.spines[s_].set_visible(False)
    ax.tick_params(labelsize=FT, length=2, width=0.5, pad=1.5)
    for s_ in ('left', 'bottom'): ax.spines[s_].set_linewidth(0.6)


# ------------------------------------------------------------------ a  Manhattan skyline
letter(0, 0.5, 'a')
ax = axes(fig, 8, 2, 171, 50)
P = V[V.chr.notna()].copy(); P['x'], spans = fk.hero.genome(P.chr, P.bp, sex=True); P = P[P.x.notna()]
base, hl = P[~P.hl], P[P.hl]
ob = base.chr.map(lambda s: 22 if s == 'X' else int(s) - 1).values % 2 == 1
ax.add_collection(LC([((x, 0), (x, h)) for x, h in zip(base.x, base.papers)], colors=np.where(ob, NAV, MAR), linewidths=0.5, rasterized=True))
ax.add_collection(LC([((x, 0), (x, h)) for x, h in zip(hl.x, hl.papers)], colors=GRN, linewidths=0.7, rasterized=True, zorder=3))
top = P.papers.max()
for i, (lb, x0, x1) in enumerate(spans):
    ax.plot([x0 + 0.002, x1 - 0.002], [-8, -8], color=GD, lw=0.5)
    ax.text((x0 + x1) / 2, -36 - (13 if lb in ('19', '21') else 0), lb, ha='center', va='center', fontsize=FT)
for v in (200, 400, 600):
    ax.plot([0, 1], [v, v], color=GRID, lw=0.3, zorder=0)
for v in (0, 200, 400, 600, 800):
    ax.text(-0.005, v, f'{v}', ha='right', va='center', fontsize=FT)
ax.plot([-0.002, -0.002], [0, 800], color=INK, lw=0.6)
ax.text(-0.035, 400, 'Papers claiming the variant', rotation=90, ha='center', va='center', fontsize=FS)
for rs in NUP + NDN:
    r = P[P.rsid == rs]
    if r.empty: continue
    r = r.iloc[0]; up = rs in NUP; left = rs in (2910164, 7903146)
    ax.scatter([r.x], [r.papers], marker='D' if up else 'o', s=11 if up else 8, color=GRN if up else 'white', edgecolor=GRE if up else INK, lw=0.5, zorder=5)
    ax.text(r.x + (-0.004 if left else 0.004), r.papers + 8, GENE[rs], ha='right' if left else 'left', va='bottom', fontsize=FT, fontstyle='italic', zorder=6)
ax.plot([0.015, 0.025], [790, 790], color=GRN, lw=2.2, solid_capstyle='butt'); ax.text(0.03, 790, f'reported associations mostly supported genome-wide ({int(P.hl.sum()):,} variants)', fontsize=FT, va='center')
ax.plot([0.015, 0.0195], [735, 735], color=MAR, lw=2.2, solid_capstyle='butt'); ax.plot([0.0205, 0.025], [735, 735], color=NAV, lw=2.2, solid_capstyle='butt')
ax.text(0.03, 735, f'other claimed variants, chromosomes alternating ({len(P) - int(P.hl.sum()):,})', fontsize=FT, va='center')
ax.scatter([0.0185], [680], marker='D', s=11, color=GRN, edgecolor=GRE, lw=0.5); ax.scatter([0.032], [680], marker='o', s=8, color='white', edgecolor=INK, lw=0.5)
ax.text(0.04, 680, 'named: mostly / rarely supported', fontsize=FT, va='center')
n_claims, n_papers = len(CL), CL.pmid.nunique()
ax.text(0.47, 790, f'{n_claims:,}', fontsize=9, fontweight='bold', ha='left', va='center')
ax.text(0.47, 735, f'disease and trait claims on {len(V):,} variants\nin {n_papers:,} papers', fontsize=FT, color=GD, va='top', linespacing=1.05)
ax.set_xlim(-0.045, 1.005); ax.set_ylim(-62, 820); ax.axis('off')

# ------------------------------------------------------------------ b  square scatter with marginals
letter(0, 55, 'b')
x = V[V.n_t >= 10].copy(); lx = np.log10(x.n_t); hi = (x.sup > 0.5).values
axm = axes(fig, 11, 68, 44, 44); axt = axes(fig, 11, 58.5, 44, 8.5, sharex=axm); axr = axes(fig, 56.5, 68, 8.5, 44, sharey=axm)
rng = np.random.default_rng(1)
jx = lx + rng.uniform(-0.012, 0.012, len(x)); jy = np.clip(x.sup + rng.uniform(-0.01, 0.01, len(x)), -0.01, 1.01)
axm.scatter(jx[~hi], jy[~hi], s=6, color=TAN, edgecolor=TAND, lw=0.15, rasterized=True)
axm.scatter(jx[hi], jy[hi], s=6, color=GRN, edgecolor=GRE, lw=0.15, rasterized=True)
nl = pos.null_share.mean(); axm.axhline(nl, color=INK, lw=0.5, ls=(0, (2, 1.5)))
for rs, tx, ty in [(9939609, 2.45, 0.68), (1800795, 2.45, 0.30), (738409, 2.15, 0.80), (1801133, 2.2, 0.14)]:
    r = x[x.rsid == rs].iloc[0]; up = r.sup > 0.5
    axm.scatter([np.log10(r.n_t)], [r.sup], s=11, marker='D' if up else 'o', color=GRN if up else 'white', edgecolor=GRE if up else INK, lw=0.5, zorder=5)
axm.set_xticks(np.log10([10, 30, 100, 300])); axm.set_xticklabels(['10', '30', '100', '300']); axm.set_xlim(0.97, 2.78)
axm.set_yticks([0, .25, .5, .75, 1]); axm.set_yticklabels(['0', '25', '50', '75', '100']); axm.set_ylim(-0.03, 1.03)
axm.set_xlabel('Reported associations of the variant', fontsize=FS, labelpad=1); axm.set_ylabel('% supported genome-wide', fontsize=FS, labelpad=1); style(axm)
xb = np.linspace(1, 2.78, 24); yb = np.linspace(0, 1, 21)
for msk, c_ in ((~hi, TAN), (hi, GRN)):
    axt.hist(lx[msk], bins=xb, color=c_, alpha=0.7, lw=0); axr.hist(x.sup[msk], bins=yb, orientation='horizontal', color=c_, alpha=0.7, lw=0)
axt.axis('off'); axr.axis('off')
c0 = np.histogram(x.sup, bins=yb)[0]; axr.set_xlim(0, c0.max() * 1.6)
axr.text(c0[0] * 1.05, 0.025, f'{int((x.sup < 0.05).sum()):,}', fontsize=FT, color=TAND, va='center')
axr.text(c0[-1] * 1.15, 0.975, f'{int((x.sup > 0.8).sum()):,}', fontsize=FT, color=GRE, va='center')
axt.text(2.76, axt.get_ylim()[1] * 0.95, f'{len(x):,} variants\ndashed line: matched null', fontsize=FT - 0.5, color=GD, ha='right', va='top', linespacing=1.1)

# ------------------------------------------------------------------ c  six variants x trait areas x year, shared rows
letter(70, 55, 'c')
SIX = [9939609, 738409, 7903146, 1801133, 1800795, 1800629]
AR = [a for a in AREAS if a != 'Not mapped']; ym = {a: i for i, a in enumerate(AR[::-1])}
X0, FW, GAP = 98, 12.2, 1.2
bar = axes(fig, X0 - 2.2, 64, 1.4, 48); bar.set_ylim(-0.5, len(AR) - 0.5); bar.axis('off')
for a in AR: bar.barh(ym[a], 1, height=1, color=C[a], lw=0)
for a in AR: ftext(fig, X0 - 3, 64 + (len(AR) - 0.5 - ym[a]) / len(AR) * 48, F.AREA_SHORT.get(a, a).replace(' disease', ''), fontsize=FT, ha='right', va='center')
for j, rs in enumerate(SIX):
    q = pos[(pos.rsid == rs) & pos.grp.isin(AR)]
    g = q.groupby(['grp', 'pub_year']).agg(n=('pmid', 'nunique'), s=('replicated', 'mean')).reset_index()
    fa = axes(fig, X0 + j * (FW + GAP), 64, FW, 48)
    for yv in range(len(AR)): fa.axhline(yv, color=GRID, lw=0.3, zorder=0)
    for xv in (2010, 2020): fa.axvline(xv, color=GRID, lw=0.3, zorder=0)
    fa.axvline(2004.5, color=INK, lw=0.5, ls=(0, (2.5, 2)))
    sup = g.s > 0.5
    fa.scatter(g.pub_year, g.grp.map(ym), s=4 + 7 * g.n, c=np.where(sup, GRN, TAN), edgecolors=np.where(sup, GRE, TAND), lw=0.3, alpha=0.9, zorder=3)
    fa.set_xlim(2004, 2026.5); fa.set_ylim(-0.6, len(AR) - 0.4); fa.set_yticks([]); fa.set_xticks([2010, 2020]); fa.set_xticklabels(["'10", "'20"])
    for s_ in ('top', 'right', 'left'): fa.spines[s_].set_visible(False)
    fa.tick_params(labelsize=FT - 0.5, length=1.5, pad=1)
    allq = pos[pos.rsid == rs]
    fa.text(0.5, 1.085, GENE[rs], transform=fa.transAxes, ha='center', fontsize=FS, fontstyle='italic', fontweight='bold')
    fa.text(0.5, 1.02, f'{allq.replicated.mean() * 100:.0f}%', transform=fa.transAxes, ha='center', fontsize=FT, color=GRE if allq.replicated.mean() > 0.5 else TAND)
ftext(fig, X0 + 3 * (FW + GAP) - GAP / 2, 116.6, 'Publication year', fontsize=FS, ha='center')
ftext(fig, X0 + 1.5 * (FW + GAP) - GAP / 2, 57.2, 'mostly supported', fontsize=FT, ha='center', color=GRE)
ftext(fig, X0 + 4.5 * (FW + GAP) - GAP / 2, 57.2, 'almost never supported', fontsize=FT, ha='center', color=TAND)
lg = axes(fig, 72, 92, 12, 20); lg.axis('off'); lg.set_xlim(0, 1); lg.set_ylim(0, 1)
lg.text(0, 0.98, 'papers per\narea-year', fontsize=FT - 0.5, va='top', linespacing=1.0)
for k, n_ in enumerate((1, 3, 6)):
    lg.scatter([0.15], [0.62 - k * 0.2], s=4 + 7 * n_, color=GD); lg.text(0.38, 0.62 - k * 0.2, str(n_), fontsize=FT - 0.5, va='center')

# ------------------------------------------------------------------ d  claims per year by area (kept)
letter(0, 121, 'd')
ax = axes(fig, 10, 126, 66, 44)
YA = CL.groupby(['pub_year', 'grp']).size().unstack(fill_value=0).reindex(columns=AREAS, fill_value=0).reindex(range(2000, 2026), fill_value=0)
ax.stackplot(YA.index, *[YA[a] / 1e3 for a in AREAS], colors=[C[a] for a in AREAS], lw=0.3, edgecolor='white')
tot = YA.sum(axis=1)
cum = np.cumsum([YA.loc[2025, a] for a in AREAS]) / 1e3
ys = [cum[i] - YA.loc[2025, a] / 2e3 for i, a in enumerate(AREAS)]
ys2 = list(ys)
for i in range(1, len(ys2)): ys2[i] = max(ys2[i], ys2[i - 1] + 0.72)
for a, y0, y1 in zip(AREAS, ys, ys2):
    ax.text(2026.6, y1, F.AREA_SHORT.get(a, a).replace(' disease', ''), fontsize=FT - 0.5, va='center', color=INK)
    ax.plot([2025.1, 2026.4], [y0, y1], color=GL, lw=0.4)
ax.annotate(f'{tot.max() / 1e3:.1f} thousand in 2015', (2015, tot.max() / 1e3), (2004, tot.max() / 1e3 * 0.98), fontsize=FT, ha='left', va='center',
            arrowprops=dict(arrowstyle='-', lw=0.4, color=GD))
ax.set_xlim(2000, 2025); ax.set_xticks([2000, 2010, 2020]); ax.set_ylim(0, tot.max() / 1e3 * 1.06)
ax.set_ylabel('Claims per year (thousands)', fontsize=FS, labelpad=1); style(ax)

# ------------------------------------------------------------------ e  FTO through time (kept)
letter(98, 121, 'e')
FTO = CL[CL.rsid == 9939609].copy()
merge = {'bmi': 'Body mass index', 'body mass index': 'Body mass index', 'obesity': 'Obesity', 'type 2 diabetes': 'Type 2 diabetes',
         'type 2 diabetes mellitus': 'Type 2 diabetes', 'overweight': 'Overweight', 'waist circumference': 'Waist circumference', 'metabolic syndrome': 'Metabolic syndrome'}
FTO['t'] = FTO.o.map(merge).fillna('other')
years = np.arange(2007, 2026)
per = lambda q: q.drop_duplicates('pmid').groupby('pub_year').size().reindex(years, fill_value=0).values   # noqa: E731
ax = axes(fig, 108, 126, 30, 44)
go = [a for a in AREAS if (FTO.grp == a).any()]
ax.stackplot(years, *[per(FTO[FTO.grp == a]) for a in go], colors=[C[a] for a in go], lw=0.3, edgecolor='white')
ax.set_xlim(2007, 2025); ax.set_xticks([2010, 2020]); ax.set_ylabel('FTO rs9939609: papers per year', fontsize=FS, labelpad=1); style(ax)
ax.text(0.0, 1.02, f'{FTO.pmid.nunique()} papers, {FTO.trait_id.nunique()} traits', transform=ax.transAxes, fontsize=FT, va='bottom')
named = ['Obesity', 'Body mass index', 'Type 2 diabetes', 'Waist circumference', 'Metabolic syndrome', 'Overweight']
SHORT = {'Waist circumference': 'Waist circumf.', 'Metabolic syndrome': 'Metabolic synd.'}
for k, t in enumerate(named):
    q = FTO[FTO.t == t]; i, j = divmod(k, 2)
    mm = axes(fig, 143 + j * 18.5, 128 + i * 15.0, 15.5, 10)
    mm.fill_between(years, 0, per(q), color=GD, lw=0); mm.set_xlim(2007, 2025); mm.set_ylim(0, 20); mm.set_yticks([0, 20] if j == 0 else [])
    mm.set_xticks([2010, 2020] if i == 2 else []); mm.set_xticklabels(["'10", "'20"] if i == 2 else []); style(mm); mm.tick_params(labelsize=FT - 1)
    mm.text(2007.3, 23.5, SHORT.get(t, t), fontsize=FT - 1, va='bottom')
    mm.text(2024.8, 19, f'{q.pmid.nunique()}', fontsize=FT - 1, va='top', ha='right', color=GD)


def blacken(fig):
    """User 2026-11-17: every non-essential text in black; only white labels on filled marks stay white."""
    import matplotlib.text, matplotlib.colors as _mc
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text() and _mc.to_rgb(t.get_color()) != (1.0, 1.0, 1.0):
            t.set_color('black')


blacken(fig)
import os as _os
if _os.environ.get('VC_AUDIT'):
    sys.path.insert(0, str(Path(__file__).resolve().parent)); import vcaudit; vcaudit.audit(fig, 'fig1')
fig.savefig(OUT / 'fig1.png', dpi=600, facecolor='white'); fig.savefig(OUT / 'fig1.pdf', facecolor='white'); print('wrote', OUT / 'fig1')

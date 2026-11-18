"""Figure 2 v3 (2026-11-17): rebuilt at the proportions and density of the user's main reference PMC10335933_F3 (Nature Genetics 2023,
Fig. 3), after the user judged v2 not a reproduction ("完全做不到复刻参考图，增加信息密度和各部分元素的比例").
Like the reference: no panel titles (statements go to the legend), panel boxes in the reference's proportions (a 1/3 + b 2/3 on top;
a large circos c with d and e stacked beside it), curves with an inset pie, ~100 sorted stacked bars with a group dot strip, a genome
circos with heat, bar and line tracks around a grey legend disc, a histogram with an inset pie, and group curves.
Data: held-out reported associations with a stated sample size (50,312) for a, b, d, e; c uses every claimed variant with a genome
position (both splits; descriptive) in 5-Mb windows.
a  % supported across 40 equal-count size bins: all supported, of which restating an earlier genome-wide report, and came first;
   matched null dashed; inset pie of size classes.
b  100 equal-count size bins: stacked bars restated / came first; dot strip = the most common study population of the bin.
c  circos, 5-Mb windows: outer ideogram; heat = claims per window; bars = % of reported associations supported; bars = % from studies
   of 2,000 or more; inner line = GWAS Catalog associations (P <= 5e-8) per window.
d  stated sample sizes; inset pie of study populations (colours shared with the dot strip in b and the curves in e).
e  % supported against sample size per study population.
Output figures/v3/fig2."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib, matplotlib.colors as mc, matplotlib.ticker
sys.path.insert(0, str(Path(__file__).resolve().parents[1])); import figdata as F
sys.path.insert(0, str(Path(__file__).resolve().parent)); import roles as RL
sys.path.insert(0, '/Volumes/ResearchData/Claude/figtaste')
from layout.canvas import figure, axes, ftext
from pycirclize import Circos
matplotlib.rcParams['savefig.bbox'] = 'standard'   # pycirclize sets 'tight', which re-sizes the page

# palette sampled from PMC10335933_F3
TEAL, LBLUE, DTEAL = RL.SUP, RL.REST, RL.dark(RL.SUP)
PALE, SLATE = RL.light(RL.SUP, 0.7), RL.dark(RL.SUP, 0.3)   # size-class ramp of the 'supported' role
RED, GREEN = RL.FIRST, RL.REST
INK, GD, GRID = '#1a1a1a', '#555555', '#e6e6e6'
FS, FT, FL = 6.5, 6.0, 9.5
POP = {'EUR': 'European', 'EAS': 'East Asian', 'other': 'Other or mixed', 'not_stated': 'Not stated'}
PCOL = {'European': RL.G1, 'East Asian': RL.G2, 'Other or mixed': RL.G3, 'Not stated': RL.G4}
OUT = RL.OUT

d = F.claims('heldout'); P = d[(d.status == 'positive') & (d.n > 0)].copy()
P['lg'] = np.log10(P.n); P['pop'] = P.pop_group.map(POP)
P['first'] = P.replicated & ~P.rep_prior; P['rest'] = P.replicated & P.rep_prior
with F.con() as c:
    C2 = c.sql("select * from confirm_primary where test = 'C2'").df().iloc[0]
F.check(len(P) == 50312, 'held-out reported associations with a stated sample size: 50,312')


def letter(x, y, s):
    ftext(fig, x, y, s, fontsize=FL, fontweight='bold', va='top')


def style(ax):
    for s_ in ('top', 'right'): ax.spines[s_].set_visible(False)
    ax.tick_params(labelsize=FT, length=2, width=0.5, pad=1.5)
    for s_ in ('left', 'bottom'): ax.spines[s_].set_linewidth(0.6)


def qbins(x, k):
    x = x.sort_values('lg').copy(); x['b'] = pd.qcut(x.lg.rank(method='first'), k, labels=False); return x


fig = figure(180, 150)

# ------------------------------------------------------------------ a  three curves + inset pie (reference a)
letter(0, 0.5, RL.PL('a'))
ax = axes(fig, 9, 3, 50, 44)
Q = qbins(P, 40)
g = Q.groupby('b').agg(lg=('lg', 'median'), s=('replicated', 'mean'), r=('rest', 'mean'), f=('first', 'mean'), nl=('null_share', 'mean'), m=('pmid', 'size'))
for col_, c_ in (('s', TEAL), ('r', LBLUE), ('f', RED)):
    se = np.sqrt(g[col_] * (1 - g[col_]) / g.m) * 196
    ax.fill_between(g.lg, g[col_] * 100 - se, g[col_] * 100 + se, color=c_, alpha=0.22, lw=0)
ax.plot(g.lg, g.s * 100, color=TEAL, lw=1.4)
ax.plot(g.lg, g.r * 100, color=LBLUE, lw=1.4)
ax.plot(g.lg, g.f * 100, color=RED, lw=1.4)
ax.plot(g.lg, g.nl * 100, color='#8a8a8a', lw=0.9, ls=(0, (2.5, 1.5)))
ax.set_xticks([2, 3, 4, 5]); ax.set_xticklabels(['100', '1k', '10k', '100k']); ax.set_xlim(1.75, 5.2)
ax.set_ylim(0, 52); ax.set_yticks([0, 10, 20, 30, 40, 50])
ax.set_xlabel('Stated sample size', fontsize=FS, labelpad=1); ax.set_ylabel('% of reported associations', fontsize=FS, labelpad=1)
style(ax)
for k, (lab, col) in enumerate([('supported genome-wide', TEAL), ('restated an earlier report', LBLUE), ('came before the report', RED), ('matched null', '#8a8a8a')]):
    ax.plot([1.85, 2.05], [49 - k * 4.0] * 2, color=col, lw=1.4 if k < 3 else 0.9, ls='-' if k < 3 else (0, (2.5, 1.5)))
    ax.text(2.1, 49 - k * 4.0, lab, fontsize=FT, va='center')
ax.text(1.85, 33.5, f'OR {C2.estimate:.2f} per 10-fold\nlarger sample', fontsize=FT, color=DTEAL, fontweight='bold', va='top', linespacing=1.05)
cls = pd.cut(P.n, [0, 500, 2000, 10000, np.inf], labels=['<500', '0.5–2k', '2–10k', '≥10k'], right=False).value_counts().sort_index()
pa = ax.inset_axes([0.62, 0.08, 0.2, 0.2])
pa.pie(cls.values, colors=[PALE, LBLUE, TEAL, SLATE], startangle=90, counterclock=False, wedgeprops=dict(lw=0.6, edgecolor='white'))
pa.set_aspect('equal')
for k, (lab, v) in enumerate(cls.items()):
    ang = np.deg2rad(90 - (cls.values[:k].sum() + v / 2) / cls.sum() * 360)
    pa.text(1.35, 0.85 - k * 0.55, f'{lab} {v / cls.sum() * 100:.0f}%', ha='left', va='center', fontsize=FT - 0.5)
    pa.add_patch(matplotlib.patches.Rectangle((1.1, 0.72 - k * 0.55), 0.18, 0.26, color=[PALE, LBLUE, TEAL, SLATE][k], lw=0, clip_on=False))

# ------------------------------------------------------------------ b  100 sorted stacked bars + population dot strip (reference b)
letter(64, 0.5, RL.PL('b'))
ax = axes(fig, 74, 3, 104, 41)
Q = qbins(P, 100)
B = Q.groupby('b').agg(r=('rest', 'mean'), f=('first', 'mean'), med=('n', 'median'), pop=('pop', lambda s: s.mode().iat[0]))
xb = np.arange(100)
ax.bar(xb, B.r * 100, width=0.86, color=TEAL, lw=0)
ax.bar(xb, B.f * 100, bottom=B.r * 100, width=0.86, color=RED, lw=0)
ax.set_xlim(-0.8, 99.8); ax.set_ylim(0, 62); ax.set_yticks([0, 20, 40, 60]); ax.set_xticks([])
ax.set_ylabel('% supported genome-wide', fontsize=FS, labelpad=1); style(ax)
ds = axes(fig, 74, 44.6, 104, 2.2); ds.set_xlim(-0.8, 99.8); ds.set_ylim(0, 1); ds.axis('off')
ds.scatter(xb, np.full(100, 0.5), s=5, c=[PCOL[p] for p in B['pop']], lw=0)
for i in (0, 33, 66, 99):
    ds.text(i, -0.9, f'{int(B.med.iloc[i]):,}', ha='center' if 0 < i < 99 else ('left' if i == 0 else 'right'), va='top', fontsize=FT)
ftext(fig, 126, 51.2, '100 equal-count bins of claims by stated sample size (bin median); dot: most common study population',
      fontsize=FS, ha='center', va='top')
for k, (lab, col) in enumerate([('restated an earlier genome-wide report', TEAL), ('came before the report', RED)]):
    ax.add_patch(matplotlib.patches.Rectangle((2, 56 - k * 5.2), 1.6, 3.2, color=col, lw=0)); ax.text(4.6, 57.6 - k * 5.2, lab, fontsize=FT, va='center')
for k, lab in enumerate(PCOL):
    ax.scatter([44], [57.6 - k * 4.2], s=8, color=PCOL[lab], lw=0); ax.text(45.5, 57.6 - k * 4.2, lab, fontsize=FT, va='center')
pc = P['pop'].value_counts().reindex(list(PCOL))
pb = ax.inset_axes([0.66, 0.50, 0.20, 0.50])
pb.pie(pc.values, colors=[PCOL[k] for k in pc.index], startangle=90, counterclock=False, wedgeprops=dict(lw=0.6, edgecolor='white')); pb.set_aspect('equal')
for k, (lab, v) in enumerate(pc.items()):
    ang = np.deg2rad(90 - (pc.values[:k].sum() + v / 2) / pc.sum() * 360)
    if v / pc.sum() < 0.12: continue
    pb.text(0.62 * np.cos(ang), 0.62 * np.sin(ang), f'{v / pc.sum() * 100:.0f}%', ha='center', va='center', fontsize=FT - 0.5, color='white')

# ------------------------------------------------------------------ c  genome circos (reference c)
letter(0, 55, RL.PL('c'))
V = F.variants()
dd = pd.concat([F.claims('dev'), F.claims('heldout')]); pp = dd[dd.status == 'positive'].copy()
S = pp.groupby('rsid').agg(n_t=('pmid', 'size'), sup=('replicated', 'sum'), big=('n', lambda s: (pd.to_numeric(s, errors='coerce') >= 2000).sum()),
                           nn=('n', lambda s: pd.to_numeric(s, errors='coerce').notna().sum())).reset_index()
V = V[V.chr.notna() & (V.chr != 'X')].merge(S, on='rsid', how='left').fillna({'n_t': 0, 'sup': 0, 'big': 0, 'nn': 0})
with F.con() as cc:
    G = cc.sql("select chr_id chr, try_cast(chr_pos as bigint) bp from gwas_catalog_rsid_trait where pvalue_mlog >= 7.30103").df()
G = G[G.chr.isin([str(i) for i in range(1, 23)])].dropna()
W = 5_000_000
V['w'] = (V.bp // W).astype(int); G['w'] = (G.bp // W).astype(int)
CHR = [str(i) for i in range(1, 23)]
L = {c: int(max(V[V.chr == c].bp.max(), G[G.chr == c].bp.max())) + 1 for c in CHR}
cax = axes(fig, 4, 63, 84, 84, projection='polar')
circ = Circos({c: L[c] for c in CHR}, space=1.6, start=4, end=356)
heat = RL.seq(RL.SUP, 0.92)
for sec in circ.sectors:
    c = sec.name; v = V[V.chr == c]; gg = G[G.chr == c]
    nw = L[c] // W + 1; x = np.minimum(np.arange(nw) * W + W / 2, L[c] - W * 0.46)
    agg = v.groupby('w').agg(cl=('claims', 'sum'), nt=('n_t', 'sum'), sup=('sup', 'sum'), big=('big', 'sum'), nn=('nn', 'sum')).reindex(range(nw), fill_value=0)
    t0 = sec.add_track((96, 100)); t0.axis(fc=RL.G2 if int(c) % 2 else RL.G3, ec='none')
    sec.text(c, r=105, size=FT, orientation='horizontal')
    t1 = sec.add_track((89, 95))
    hv = np.minimum(np.log10(agg.cl.values + 1), 3.6); t1.heatmap(hv.reshape(1, -1), vmin=0, vmax=3.6, cmap=heat, rect_kws=dict(ec='none'))
    t2 = sec.add_track((74, 88)); t2.axis(fc=RL.light(RL.SUP, 0.93), ec='none')
    sh = np.where(agg.nt > 5, agg.sup / agg.nt.replace(0, np.nan), 0); t2.bar(x, np.nan_to_num(sh) * 100, width=W * 0.9, vmin=0, vmax=100, color=TEAL, lw=0)
    t3 = sec.add_track((60, 73))
    bs = np.where(agg.nn > 5, agg.big / agg.nn.replace(0, np.nan), 0); t3.bar(x, np.nan_to_num(bs) * 100, width=W * 0.9, vmin=0, vmax=100, color=RL.G2, lw=0)
    t4 = sec.add_track((48, 59))
    gw = gg.groupby('w').size().reindex(range(nw), fill_value=0).values
    t4.fill_between(x, np.minimum(np.log10(gw + 1), 3.6), y2=0, vmin=0, vmax=3.6, color=GREEN, lw=0.4, alpha=0.9)
circ.plotfig(ax=cax)
fig.set_size_inches(180 / 25.4, 150 / 25.4)
cax.add_artist(matplotlib.patches.Circle((0, 0), 0.45, transform=cax.transData._b, color='#ececec', zorder=0)) if False else None
ins = axes(fig, 28.5, 88.5, 35, 34); ins.axis('off'); ins.set_xlim(0, 1); ins.set_ylim(0, 1)
ins.add_patch(matplotlib.patches.Circle((0.5, 0.5), 0.5, color='#ececec', lw=0))
for k, (lab, col) in enumerate([('Claims per 5 Mb', None), ('% supported', TEAL), ('% from studies ≥ 2,000', RL.G2), ('GWAS associations', GREEN)]):
    y = 0.74 - k * 0.15
    if col is None:
        for j, cc_ in enumerate(heat(np.linspace(0, 1, 5))): ins.add_patch(matplotlib.patches.Rectangle((0.62 + j * 0.05, y - 0.035), 0.05, 0.07, color=cc_, lw=0))
    else:
        ins.add_patch(matplotlib.patches.Rectangle((0.62, y - 0.035), 0.1, 0.07, color=col, lw=0))
    ins.text(0.58, y, lab, ha='right', va='center', fontsize=FT)

# ------------------------------------------------------------------ d  sample-size histogram + population pie (reference d)
letter(98, 55, RL.PL('d'))
ax = axes(fig, 109, 59, 69, 34)
bins = np.linspace(1, 6, 46)
for lab in PCOL:
    pass
h, e_ = np.histogram(P.lg, bins=bins)
ax.bar(e_[:-1], h, width=np.diff(e_) * 0.9, align='edge', color=RL.G3, lw=0)
ax.set_xticks([1, 2, 3, 4, 5, 6]); ax.set_xticklabels(['10', '100', '1k', '10k', '100k', '1M']); ax.set_xlim(1, 6)
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
ax.set_xlabel('Stated sample size', fontsize=FS, labelpad=1); ax.set_ylabel('Reported associations', fontsize=FS, labelpad=1); style(ax)
md = np.log10(P.n.median()); ax.axvline(md, color=INK, lw=0.6, ls=(0, (2, 1.5)))
ax.set_ylim(0, h.max() * 1.22); ax.text(md + 0.06, h.max() * 1.18, f'median {int(P.n.median()):,}', fontsize=FT, va='top')
cls2 = pd.cut(P.n, [0, 500, 2000, 10000, np.inf], labels=['<500', '0.5–2k', '2–10k', '≥10k'], right=False)
dc = P.design.where(P.design.isin(['case-control', 'cohort', 'meta-analysis', 'GWAS']), 'other').value_counts().reindex(['case-control', 'cohort', 'meta-analysis', 'GWAS', 'other'])
DC = {'case-control': RL.G4, 'cohort': RL.G3, 'meta-analysis': RL.G2, 'GWAS': RL.G1, 'other': '#ECE8E2'}
pd_ = ax.inset_axes([0.47, 0.42, 0.24, 0.58])
pd_.pie(dc.values, colors=[DC[k] for k in dc.index], startangle=90, counterclock=False, wedgeprops=dict(lw=0.6, edgecolor='white')); pd_.set_aspect('equal')
for k, (lab, v) in enumerate(dc.items()):
    ax.add_patch(matplotlib.patches.Rectangle((4.62, h.max() * (0.93 - k * 0.12)), 0.1, h.max() * 0.07, color=DC[lab], lw=0, clip_on=False))
    ax.text(4.76, h.max() * (0.965 - k * 0.12), f'{lab} {v / dc.sum() * 100:.0f}%', fontsize=FT - 0.5, va='center')

# ------------------------------------------------------------------ e  curves per population (reference e)
letter(98, 97, RL.PL('e'))
ax = axes(fig, 109, 101, 69, 38)
for lab, col in PCOL.items():
    x = qbins(P[P['pop'] == lab], 8)
    gg = x.groupby('b').agg(lg=('lg', 'median'), s=('replicated', 'mean'), m=('pmid', 'size'))
    se = np.sqrt(gg.s * (1 - gg.s) / gg.m) * 196
    ax.fill_between(gg.lg, gg.s * 100 - se, gg.s * 100 + se, color=col, alpha=0.18, lw=0)
    ax.plot(gg.lg, gg.s * 100, color=col, lw=1.3)
    ax.scatter(gg.lg, gg.s * 100, s=14, color=col, lw=0)
ax.set_xticks([2, 3, 4, 5]); ax.set_xticklabels(['100', '1k', '10k', '100k']); ax.set_xlim(1.6, 5.3); ax.set_ylim(0, 52)
ax.set_xlabel('Stated sample size', fontsize=FS, labelpad=1); ax.set_ylabel('% supported genome-wide', fontsize=FS, labelpad=1); style(ax)
for k, (lab, col) in enumerate(PCOL.items()):
    ax.plot([1.75, 1.95], [48 - k * 4.5] * 2, color=col, lw=1.3); ax.text(2.0, 48 - k * 4.5, lab, fontsize=FT, va='center')


def blacken(fig):
    """User 2026-11-17: every non-essential text in black; only white labels on filled marks stay white."""
    import matplotlib.text, matplotlib.colors as _mc
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text() and _mc.to_rgb(t.get_color()) != (1.0, 1.0, 1.0):
            t.set_color('black')


blacken(fig)
import os as _os
if _os.environ.get('VC_AUDIT'):
    sys.path.insert(0, str(Path(__file__).resolve().parent)); import vcaudit; vcaudit.audit(fig, 'fig2')
fig.set_size_inches(180 / 25.4, 150 / 25.4); fig.savefig(OUT / 'fig2.png', dpi=600, facecolor='white', bbox_inches=None); fig.savefig(OUT / 'fig2.pdf', facecolor='white', bbox_inches=None); print('wrote', OUT / 'fig2')

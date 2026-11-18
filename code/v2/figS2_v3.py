"""Supplementary Figure 2 v3 (2026-11-17): the AJHG figure 2 (now S2) redrawn to the rules of the v3 main figures: lower-case letters
without title lines (statements in the legend), colour roles of roles.py (palette SF), all non-essential text black, no text on
filled badges, text-on-mark audit (VC_AUDIT=1). Output figures/colour/<PAL>/figS2.
Figure 2 (design record figures/design.toml). Genome-wide support of reported associations against the matched expectation.
a, of 100 reported associations and of 100 frequency-matched random variants, the number with genome-wide support for the claimed trait
(vermilion: the genome-wide report came after the claim; blue: it came first), then the same count for disease, quantitative-trait and
severity claims; rows share one scale, ten blocks of ten. b, observed share and matched null by minor-allele frequency (band, 95%
bootstrap interval over papers). c, the count per 100 in each trait area (the area's colour of Figure 1; tick, matched null; 95%
bootstrap interval beside the count). d, the share in each trait area at each allele frequency, rows as in c.
Usage: python 42_fig2_support.py [dev|heldout]."""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
from matplotlib.colors import to_hex
from matplotlib.patches import Rectangle
import figkit as fk
sys.path.insert(0, str(Path(__file__).resolve().parents[1])); import figdata as F
sys.path.insert(0, str(Path(__file__).resolve().parent)); import roles as RL
SPLIT = sys.argv[1] if len(sys.argv) > 1 else 'heldout'
e = fk.enc('varclaims')
VER, GD, GR, GL, WH = RL.SUP, fk.c.grey_dark, fk.c.grey, fk.c.grey_light, fk.c.white
GN = fk.c['grey@s0.2']                                         # notes and secondary text
FB, FT, FA = fk.pt('base'), fk.pt('tick'), fk.pt('annot')
d = F.claims(SPLIT); P = d[d.status == 'positive'].copy()
P['grp'] = F.area_group(P.area)
LATER, EARLIER = RL.FIRST, RL.REST                      # semantic table "support": Later (the claim came first), Earlier


def per100(x):
    """(later-only, earlier) per 100 claims; the supported total is rounded first, then split."""
    tot = int(round(100 * x.replicated.mean())); early = min(int(round(100 * (x.replicated & x.rep_prior).mean())), tot)
    return tot - early, early


fig = fk.figure('double', height=155, profile=fk.target())
gs = fig.add_gridspec(2, 2, height_ratios=[0.6, 1], width_ratios=[1.5, 1])

# ------------------------------------------------------------------ a  of 100 reported associations, of 100 matched random variants
ax = fig.add_subplot(gs[0, 0]); fk.panel(ax, RL.PL('a'), '')
k = per100(P); nul = int(round(100 * P.null_share.mean()))
F.check(sum(k) >= 5 * nul and k[1] > k[0], 'a: support far above the matched expectation; most of it preceded the claim')
rows = [dict(label='Reported associations', parts=[(k[0], LATER), (k[1], EARLIER)], note=f'{len(P):,} claims', hero=True, color='ink'),
        dict(label='Random variants of the same allele frequency, same traits', parts=[(nul, 'grey')], value=f'{nul} in 100', note='by chance')]
lab = {'disease': 'Disease', 'trait': 'Quantitative trait', 'severity': 'Severity or progression'}
for i, t in enumerate(('disease', 'trait', 'severity')):
    x = P[P.outcome_type == t]; kk = per100(x)
    rows.append(dict(label=lab[t], parts=[(kk[0], LATER), (kk[1], EARLIER)], note=f'{len(x):,} claims', **({'head': 'By kind of outcome'} if i == 0 else {})))
    F.check(F.boot_share(x, 'replicated')[1] > x.null_share.mean(), f'a: {lab[t]} above the matched null')
fk.forms.meter(ax, rows, value_share=0.33)
y1, h1 = ax._figkit_meter[0]                                   # the two tones: a key above the rows, the counts inside the first bar
ax.set_ylim(ax.get_ylim()[0], 0.5)
ax.text(0, 0.2, 'Genome-wide report published', ha='left', va='center', fontsize=FT, color=GD)
for x0, col, txt in ((53.0, LATER, 'after the claim'), (83.0, EARLIER, 'before the claim')):
    ax.add_patch(Rectangle((x0, 0.08), 2.4, 0.24, fc=col, ec='none', gid='figkit-unit'))
    ax.text(x0 + 3.6, 0.2, txt, ha='left', va='center', fontsize=FT, color=GD)
ax.text(k[0] / 2, y1 + h1 / 2, f'{k[0]}', ha='center', va='center', fontsize=FT, fontweight='bold', color=fk.c.ink, zorder=4)
ax.text(k[0] + 2 + k[1] / 2, y1 + h1 / 2, f'{k[1]}', ha='center', va='center', fontsize=FT, fontweight='bold', color=fk.c.ink, zorder=4)

# ------------------------------------------------------------------ b  minor-allele frequency
ax = fig.add_subplot(gs[0, 1]); fk.panel(ax, RL.PL('b'), '')
order = ['lt0.01', '1', '2', '3', '4', '5', '6']
names = {'lt0.01': '<1', '1': '1–5', '2': '5–10', '3': '10–20', '4': '20–30', '5': '30–40', '6': '40–50'}
P['bin'] = pd.cut(P.maf_eur, [-1, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.51], right=False, labels=order).astype(str)
rws = []
for b in order:
    x = P[P.bin == b]
    if len(x) < 30:
        continue
    es, l, h = F.boot_share(x, 'replicated'); rws.append((names[b], es, l, h, x.null_share.mean(), len(x)))
xx = np.arange(len(rws))
ax.fill_between(xx, 0, [r[4] * 100 for r in rws], color=GL, lw=0, zorder=1)
ax.plot(xx, [r[4] * 100 for r in rws], color=GR, lw=0.9, zorder=2)
ax.fill_between(xx, [r[2] * 100 for r in rws], [r[3] * 100 for r in rws], color=RL.light(RL.SUP, 0.6), lw=0, zorder=3)
ax.plot(xx, [r[1] * 100 for r in rws], color=VER, lw=1.6, zorder=4); ax.scatter(xx, [r[1] * 100 for r in rws], s=16, color=VER, lw=0, zorder=5)
ax.text(-0.35, rws[0][2] * 100 - 2.2, 'reported associations', ha='left', va='center', fontsize=FB, color=VER)
ax.text(xx[-1] - 0.3, rws[-1][4] * 100 + 1.6, 'matched null', ha='right', va='bottom', fontsize=FB, color=GN)
ax.vlines(xx[-1] + 0.35, rws[-1][4] * 100, rws[-1][1] * 100, color=GD, lw=0.8)
ax.hlines([rws[-1][4] * 100, rws[-1][1] * 100], xx[-1] + 0.25, xx[-1] + 0.35, color=GD, lw=0.8)
ax.text(xx[-1] + 0.5, (rws[-1][1] + rws[-1][4]) * 50, f'+{(rws[-1][1] - rws[-1][4]) * 100:.0f}\npoints', ha='left', va='center', fontsize=FT, linespacing=1.0)
F.check(all(r[2] > r[4] for r in rws), 'b: observed above the matched null in every frequency bin')
ax.set_xticks(xx); ax.set_xticklabels([r[0] for r in rws]); ax.set_xlim(-0.75, len(rws) + 1.25)
ax.set_ylim(0, 35); ax.set_yticks([0, 10, 20, 30])
fk.chart(ax, grid='y', unit='% with genome-wide support', inside=True)
ax.set_xlabel('Minor-allele frequency of the claimed variant (%)', fontsize=FB, color=GD)
fk.pin_inside(ax)

# ------------------------------------------------------------------ c  every trait area, in its colour: count per 100, matched null, interval
ax = fig.add_subplot(gs[1, 0]); fk.panel(ax, RL.PL('c'), '')
R = []
for g, x in P[P.grp != 'Not mapped'].groupby('grp'):
    est, lo, hi = F.boot_share(x, 'replicated')
    R.append(dict(grp=g, n=len(x), est=est, lo=lo, hi=hi, null=x.null_share.mean()))
R = pd.DataFrame(R).sort_values('est', ascending=False).reset_index(drop=True)
F.check((R.lo > R.null).all(), 'c: support above the matched null in every trait area')
fk.forms.meter(ax, [dict(label=r.grp, parts=[(round(100 * r.est), e.color('area', r.grp))],
                         ref=(100 * r.null, None),
                         note=f'95% CI {100 * r.lo:.0f}–{100 * r.hi:.0f}; {r.n:,} claims') for i, r in enumerate(R.itertuples())], value_share=0.33)
axm = ax

# ------------------------------------------------------------------ d  trait area by allele frequency, rows as in c
ax = fig.add_subplot(gs[1, 1]); fk.panel(ax, RL.PL('d'), '')
T = P.pivot_table(index='grp', columns='bin', values='replicated', aggfunc='mean').reindex(index=R.grp, columns=order) * 100
Nn = P.pivot_table(index='grp', columns='bin', values='replicated', aggfunc='size').reindex(index=R.grp, columns=order)
Nu = P.pivot_table(index='grp', columns='bin', values='null_share', aggfunc='mean').reindex(index=R.grp, columns=order) * 100
ok = Nn >= 30
F.check(bool((T[ok] > Nu[ok]).all().all()), f'd: observed above the matched null in all {int(ok.sum().sum())} area-by-frequency cells with 30 or more claims')
cm = RL.seq(RL.SUP, 0.92); vmin, vmax = float(np.nanmin(T[ok].values)), float(np.nanmax(T[ok].values))   # the whole ramp over the observed range
for (y0, h0), g in zip(axm._figkit_meter, R.grp):
    for j, b in enumerate(order):
        if not ok.loc[g, b]:
            continue
        v = T.loc[g, b]; col = to_hex(cm(0.04 + 0.94 * (v - vmin) / (vmax - vmin)))
        ax.add_patch(Rectangle((j + 0.04, y0 - 0.2), 0.92, 0.86, fc=col, ec='none', gid='figkit-unit'))
        ax.text(j + 0.5, y0 + 0.23, f'{v:.0f}', ha='center', va='center', fontsize=FT, color=fk.on(col))
ax.set_xlim(0, len(order)); ax.set_ylim(*axm.get_ylim())
ax.set_xticks(np.arange(len(order)) + 0.5); ax.set_xticklabels([names[b] for b in order]); ax.set_yticks([])
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_xlabel('Minor-allele frequency (%)', fontsize=FB, color=GD)
ax._figkit_kind = 'matrix'; ax._figkit_ramp = True


def blacken(fig):
    """User 2026-11-17: every non-essential text in black; only white labels on filled marks stay white."""
    import matplotlib.text, matplotlib.colors as _mc
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text() and _mc.to_rgb(t.get_color()) != (1.0, 1.0, 1.0):
            t.set_color('black')


blacken(fig)
if os.environ.get('VC_AUDIT'):
    import vcaudit; vcaudit.audit(fig, 'figS2')
saved = fk.save(fig, RL.OUT / 'figS2', strict=False); print(saved.report)

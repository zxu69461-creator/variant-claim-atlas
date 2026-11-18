"""Supplementary Figure 3 v3 (2026-11-17): the AJHG figure 3 (now S3) redrawn to the rules of the v3 main figures: lower-case letters
without title lines (statements in the legend), colour roles of roles.py (palette SF), all non-essential text black, no text on
filled badges, text-on-mark audit (VC_AUDIT=1). Output figures/colour/<PAL>/figS3.
Figure 3 (design record figures/design.toml). Study size and genome-wide support. a, share of reported associations with genome-wide
support by stated sample size (thick line and band, all claims with a 95% bootstrap interval over papers; thin lines, each trait area in
its colour of Figure 1; grey floor, matched null), with the distribution of supported and unsupported claims along the same axis on the
top edge; the odds ratio is from a logistic model with papers as clusters. b-d, share supported in every stratum, by size class
(cells with 30 or more claims) and over all sizes (bar; pale stub, matched null): b, study design; c, study population; d, publication
era. One colour scale for all cells and bars. Usage: python 43_fig3_size.py [dev|heldout]."""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
import figkit as fk
sys.path.insert(0, str(Path(__file__).resolve().parents[1])); import figdata as F
sys.path.insert(0, str(Path(__file__).resolve().parent)); import roles as RL
SPLIT = sys.argv[1] if len(sys.argv) > 1 else 'heldout'
VER, GD, GR, GL, WH = RL.SUP, fk.c.grey_dark, fk.c.grey, fk.c.grey_light, fk.c.white
GN = fk.c['grey@s0.2']                                         # notes and secondary text
FB, FT = fk.pt('base'), fk.pt('tick')
e = fk.enc('varclaims')
d = F.claims(SPLIT); P = d[d.status == 'positive'].copy()
P['grp'] = F.area_group(P.area)
with F.con() as c:
    C = c.sql(f"select * from confirm_primary{'' if SPLIT == 'heldout' else '_dev'}").df().set_index('test')

fig = fk.figure('double', height=135, profile=fk.target())
gs = fig.add_gridspec(3, 2, width_ratios=[1.45, 1], height_ratios=[7.9, 4.3, 5.2])

# ------------------------------------------------------------------ a  hero: the share supported rises between the claims themselves
ax = fig.add_subplot(gs[:, 0]); fk.panel(ax, RL.PL('a'), '')
S = P[P.n > 0].copy(); S['lg'] = np.log10(S.n)
S['bin'] = pd.qcut(S.lg, 20, labels=False, duplicates='drop')
rows = []
for b, x in S.groupby('bin'):
    es, lo, hi = F.boot_share(x, 'replicated'); rows.append((x.lg.median(), es, lo, hi, x.null_share.mean(), len(x)))
B = pd.DataFrame(rows, columns=['x', 'est', 'lo', 'hi', 'null', 'n'])
ax.fill_between(B.x, 0, B.null * 100, color=GL, lw=0, zorder=1)
ax.plot(B.x, B.null * 100, color=GN, lw=0.9, zorder=2)
# every trait area in its colour: the same gradient, area by area (four size classes, each at its median sample size)
S['sb'] = pd.cut(S.n, [0, 500, 2000, 10000, 1e9], right=False, labels=[0, 1, 2, 3]).astype(int)
areas = [a for a in F.AREA_ORDER if a != 'Not mapped']
rise = []
for a in areas:
    g = S[S.grp == a].groupby('sb').agg(x=('lg', 'median'), y=('replicated', 'mean'))
    ax.plot(g.x, g.y * 100, color=e.color('area', a), lw=0.9, marker='o', ms=2.2, zorder=3, solid_capstyle='round')
    rise.append(bool((np.diff(g.y.values) > -0.005).all()))
F.check(all(rise), 'a: the share supported rises with sample size in every trait area')
ax.fill_between(B.x, B.lo * 100, B.hi * 100, color=RL.light(RL.SUP, 0.6), lw=0, zorder=4, alpha=0.9)
ax.plot(B.x, B.est * 100, color=VER, lw=2.0, zorder=5)
ax.scatter(B.x, B.est * 100, s=12, color=VER, lw=0, zorder=6)
r = C.loc['C2']
F.check(r.ci_low > 1, 'a: support rises with stated sample size')
fk.key(ax, 1.72, 70.0, f'{r.estimate:.2f}', f'odds ratio of genome-wide support per 10-fold larger sample\n(95% CI {r.ci_low:.2f}–{r.ci_high:.2f}); '
       f'{int(r.n):,} reported associations, {int(r.n_papers):,} papers', color=VER)
for k_, a in enumerate(areas):                                 # the key of the thin lines, where the data leave room
    cx, cy = 1.74 + (k_ // 4) * 0.98, 62.0 - (k_ % 4) * 3.6
    ax.plot([cx, cx + 0.13], [cy, cy], color=e.color('area', a), lw=1.6, solid_capstyle='butt')
    ax.text(cx + 0.18, cy, F.AREA_SHORT.get(a, a), ha='left', va='center', fontsize=FT, color=GD)
ax.text(B.x.iloc[-1] + 0.06, B.est.iloc[-1] * 100, 'all reported\nassociations', ha='left', va='center', fontsize=FB, color=VER, linespacing=1.05)
ax.text(B.x.iloc[-1] + 0.06, B.null.iloc[-1] * 100 + 1.2, 'matched null', ha='left', va='center', fontsize=FB, color=GN)
ax.set_xticks([2, 3, 4, 5]); ax.set_xticklabels(['100', '1,000', '10,000', '100,000'])
ax.set_xlim(1.62, 5.9); ax.set_ylim(0, 78); ax.set_yticks([0, 10, 20, 30, 40, 50, 60])
fk.chart(ax, grid='y', unit='% with genome-wide support')
ax.set_xlabel('Stated sample size of the claiming study (log scale)', fontsize=FB, color=GD)
n_sup, n_not = int(S.replicated.sum()), int((~S.replicated).sum())           # where the claims are: the distribution on the top edge
mg = fk.marginal(ax, [S.lg[S.replicated].clip(1.7, 5.75), S.lg[~S.replicated].clip(1.7, 5.75)], ['vermilion', 'grey_light'], bins=np.arange(1.7, 5.85, 0.1),
            size=0.2, gap=0.085, labels=[f'{n_sup:,} supported claims', f'{n_not:,} claims without support'])
for p_ in mg.patches[:len(np.arange(1.7, 5.85, 0.1)) - 1]:
    p_.set_facecolor(RL.SUP)
fk.pin_inside(ax)

# ------------------------------------------------------------------ b-d  every stratum: share supported by size class, and over all sizes
from matplotlib.colors import to_hex
from matplotlib.patches import Rectangle
Q = P[P.n > 0].copy()
Q['sb'] = pd.cut(Q.n, [0, 500, 2000, 10000, 1e9], right=False, labels=[0, 1, 2, 3]).astype(int)
Q['era'] = pd.cut(Q.pub_year, [0, 2010, 2015, 2020, 3000], right=False, labels=['Before 2010', '2010–14', '2015–19', '2020 on']).astype(str)
SIZES = ['<0.5', '0.5–2', '2–10', '≥10']                       # stated sample size, thousands


def block(col, keys, names):
    """rows of one stratifying variable: (name, share by size class with 30 or more claims, share over all sizes, matched null)"""
    out = []
    for k_, nm in zip(keys, names):
        x = Q[Q[col] == k_]; allx = P[P[col] == k_] if col != 'era' else x
        g = x.groupby('sb').replicated.agg(['mean', 'size']).reindex(range(4))
        out.append((nm, [m * 100 if n_ >= 30 else np.nan for m, n_ in zip(g['mean'], g['size'])], allx.replicated.mean() * 100, allx.null_share.mean() * 100))
    return out


des = sorted(block('design', ['family', 'case-control', 'cohort', 'meta-analysis', 'GWAS'], ['Family', 'Case–control', 'Cohort', 'Meta-analysis', 'Genome-wide']),
             key=lambda t: -t[2])
pop = block('pop_group', ['EUR', 'other', 'EAS', 'not_stated'], [e.kw('population', k_)['label'] for k_ in ['EUR', 'other', 'EAS', 'not_stated']])
era = block('era', ['Before 2010', '2010–14', '2015–19', '2020 on'], ['Before 2010', '2010–14', '2015–19', '2020 on'])
F.check({des[0][0], des[1][0]} == {'Genome-wide', 'Meta-analysis'}, 'b: genome-wide and pooled designs highest')
rising = lambda rows: all(np.nanmax(np.where(np.isnan(r[1]), -1, r[1])[-2:]) > r[1][0] and np.all(np.diff([v for v in r[1] if not np.isnan(v)]) > -0.5) for r in rows)  # noqa: E731
F.check(rising(des), 'b: within every design the share supported rises with sample size')
F.check(rising(pop), 'c: the size gradient holds in every study population')
F.check(rising(era), 'd: the size gradient holds in every era')
cells = np.array([v for rows in (des, pop, era) for r in rows for v in r[1] + [r[2]] if not np.isnan(v)])
vmin, vmax = cells.min(), cells.max()
cm = RL.seq(RL.SUP, 0.92)
tone = lambda v: to_hex(cm(0.04 + 0.94 * (v - vmin) / (vmax - vmin)))      # noqa: E731  one scale for every cell and bar of b-d
LW, BW, BMAX = 3.7, 3.0, 40.0                                   # label column and bar column in cell widths; bar length at 40%


def strata(ax, rows, header=False):
    n = len(rows)
    for i_, (nm, vals, tot, nul) in enumerate(rows):
        y = n - 1 - i_
        ax.text(LW - 0.2, y + 0.5, nm, ha='right', va='center', fontsize=FT, color=GD)
        for j_, v in enumerate(vals):
            if np.isnan(v):
                continue
            ax.add_patch(Rectangle((LW + j_ + 0.03, y + 0.06), 0.94, 0.88, fc=tone(v), ec='none', gid='figkit-unit'))
            ax.text(LW + j_ + 0.5, y + 0.5, f'{v:.0f}', ha='center', va='center', fontsize=FT, color=fk.on(tone(v)))
        x0 = LW + 4.4
        ax.add_patch(Rectangle((x0, y + 0.2), tot / BMAX * BW, 0.6, fc=tone(tot), ec='none', gid='figkit-unit'))
        ax.add_patch(Rectangle((x0, y + 0.2), nul / BMAX * BW, 0.6, fc=GL, ec='none', gid='figkit-unit'))
        ax.text(x0 + tot / BMAX * BW + 0.12, y + 0.5, f'{tot:.0f}%', ha='left', va='center', fontsize=FT, color=GD, fontweight='bold')
    top = n
    if header:
        for j_, lab in enumerate(SIZES):
            ax.text(LW + j_ + 0.5, n + 0.16, lab, ha='center', va='bottom', fontsize=fk.pt('annot'), color=GD, linespacing=0.95)
        ax.text(LW + 0.03, n + 0.72, 'Stated sample size (thousands)', ha='left', va='bottom', fontsize=fk.pt('annot'), color=GN)
        ax.text(LW + 4.4, n + 0.16, 'All sizes', ha='left', va='bottom', fontsize=fk.pt('annot'), color=GD)
        top = n + 1.35
    ax.set_xlim(0, LW + 4.4 + BW + 1.25); ax.set_ylim(-0.05, top); ax.set_axis_off()
    ax._figkit_kind = 'matrix'; ax._figkit_ramp = True


ax = fig.add_subplot(gs[0, 1]); fk.panel(ax, RL.PL('b'), '')
strata(ax, des, header=True)
ax.text(LW + 4.4 + des[-1][3] / BMAX * BW + 0.12, 0.5 - 0.62, 'pale stub: matched null', ha='left', va='center', fontsize=fk.pt('annot'), color=GN)
ax.set_ylim(-0.75, ax.get_ylim()[1])
ax = fig.add_subplot(gs[1, 1]); fk.panel(ax, RL.PL('c'), '')
strata(ax, pop)
ax = fig.add_subplot(gs[2, 1]); fk.panel(ax, RL.PL('d'), '')
strata(ax, era)
ax.text(LW + 2.0, -0.5, '% with genome-wide support', ha='center', va='center', fontsize=FT, color=GN)
ax.set_ylim(-0.95, ax.get_ylim()[1])


def blacken(fig):
    """User 2026-11-17: every non-essential text in black; only white labels on filled marks stay white."""
    import matplotlib.text, matplotlib.colors as _mc
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text() and _mc.to_rgb(t.get_color()) != (1.0, 1.0, 1.0):
            t.set_color('black')


blacken(fig)
if os.environ.get('VC_AUDIT'):
    import vcaudit; vcaudit.audit(fig, 'figS3')
saved = fk.save(fig, RL.OUT / 'figS3', strict=False); print(saved.report)

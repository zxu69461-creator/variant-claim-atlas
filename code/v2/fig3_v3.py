"""Figure 3 v2 (2026-11-17; was Figure 4 of the AJHG submission). Panel a kept; b-d redrawn (user: "fig3修改bcd"; no reference fitted, forms chosen
by Claude): b, the 20 most-claimed supported variant-trait pairs as barcodes, one tick per paper; c, five-year yield by single year;
d, reports of no association at established variants, single years and three-year means. Original record: Figure 4 (design record figures/design.toml). From anticipating genome-wide reports to restating them. a, reported associations with
genome-wide support, by the number of years between the first genome-wide report of the association and the claim; one thin sheet per
three-year period of publication, each scaled to 100% of that period's supported claims (vermilion: the claim preceded or coincided
with the genome-wide report; blue: the claim restated it; both deepen with height on one scale); the odds ratio per five years is from
a logistic model with papers as clusters. b, the most frequently claimed variant-trait pairs, papers per year (vermilion before, blue
after the first genome-wide report, a blue line). c, new reported
associations (no earlier genome-wide support) that gained it within five years, by era, with authors' reports of no association and the
matched null. d, authors' reports of no association that concern a variant already genome-wide significant for the trait, by year.
Usage: python 44_fig4_time.py [dev|heldout]."""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
import figkit as fk
sys.path.insert(0, str(Path(__file__).resolve().parents[1])); import figdata as F
sys.path.insert(0, str(Path(__file__).resolve().parent)); import roles as RL
SPLIT = sys.argv[1] if len(sys.argv) > 1 else 'heldout'
VER, BLU, GD, GR, GL = RL.FIRST, RL.REST, fk.c.grey_dark, fk.c.grey, fk.c.grey_light
GN = fk.c['grey@s0.2']
FB, FT = fk.pt('base'), fk.pt('tick')
e = fk.enc('varclaims')
d = F.claims(SPLIT); d['grp'] = F.area_group(d.area); P = d[d.status == 'positive'].copy(); N = d[d.status == 'null'].copy()
with F.con() as c:
    C = c.sql(f"select * from confirm_primary{'' if SPLIT == 'heldout' else '_dev'}").df().set_index('test')

fig = fk.figure('double', height=225, profile=fk.target())
gs = fig.add_gridspec(3, 2, width_ratios=[1, 1], height_ratios=[1.42, 1.25, 0.85])

# ------------------------------------------------------------------ a  hero: the ridge of supported claims crosses the genome-wide report
ax = fig.add_subplot(gs[0, :]); fk.panel(ax, 'a', 'Supported claims moved from anticipating genome-wide reports to restating them')
R = P[P.replicated & P.pub_year.between(2005, 2025)].copy()
R['lag'] = (R.pub_year - R.first_hit_year).clip(-18, 18).astype(int)
periods = [(2005, 2007), (2008, 2010), (2011, 2013), (2014, 2016), (2017, 2019), (2020, 2022), (2023, 2025)]
lags = np.arange(-18, 19)
S = fk.solid.Scene(ax, angle=34, depth=0.62)
WB = dict(fc=fk.c.white, ec='none', pad=1.2, alpha=0.85)       # wall text on a white patch, readable over the grid and the slabs
XW, ZS, TH, DY = 16.0, 0.115, 0.07, 0.86                       # x extent; height per percentage point; sheet thickness; depth step
xs = (lags - lags.min()) / (lags.max() - lags.min()) * XW
half = (xs[1] - xs[0]) / 2
x0, x1, ny = xs[0] - half - 0.15, xs[-1] + half + 0.15, len(periods)
S.floor(x0, x1, -0.2, ny * DY + 0.1)
S.wall(x0, x1, ny * DY + 0.1, [t * ZS for t in (0, 10, 20)], y0=-0.2, labels=['0', '10%', '20%'])
xz = (0 - lags.min()) / (lags.max() - lags.min()) * XW - half      # the year of the first genome-wide report
S.line([(xz, -0.2, 0), (xz, ny * DY + 0.1, 0)], color='grey_dark', lw=0.7)
S.line([(xz, ny * DY + 0.1, 0), (xz, ny * DY + 0.1, 24 * ZS)], color='grey_dark', lw=0.7)
shares = []
WARM, COOL = RL.seq(RL.FIRST), RL.seq(RL.REST)     # a ramp of each role colour: the claim came first / the report came first
allpct = [R[R.pub_year.between(a, b)].groupby('lag').size().reindex(lags, fill_value=0).values / len(R[R.pub_year.between(a, b)]) * 100 for a, b in periods]
ZTOP = max(p.max() for p in allpct) * ZS                       # one height scale for the colour of every sheet
for j, (a, b) in enumerate(periods[::-1]):                     # the latest period in front, where the ridge is lowest
    x = R[R.pub_year.between(a, b)]
    pct = x.groupby('lag').size().reindex(lags, fill_value=0).values / len(x) * 100
    y = j * DY
    for part, cm in ((lags <= 0, WARM), (lags >= 1, COOL)):
        nz = np.where(part & (pct > 0))[0]                      # draw each sheet only where it has claims
        if len(nz) == 0:
            continue
        sel = np.arange(nz[0], nz[-1] + 1)
        px = np.repeat(np.concatenate([xs[sel] - half, [xs[sel][-1] + half]]), 2)[1:-1]
        S.slab(px, np.repeat(pct[sel] * ZS, 2), y, TH, None, gradient=cm, zmax=ZTOP, gamma=0.5, lo=0.14)
    k0 = int(np.where(lags == 0)[0][0])                         # the seam: the year of the first genome-wide report on every sheet
    S.line([(xz + 2 * half, y, 0), (xz + 2 * half, y, max(pct[k0], pct[k0 + 1]) * ZS)], color='ink', lw=0.9, key=(-y, 1.0), z=2)
    rest = 100 * x.rep_prior.mean(); shares.append(rest)
    S.text(x1, y + TH / 2, 0, f'{a}–{str(b)[2:]}', ha='left', va='center', dx=3, fontsize=FB, color='grey_dark')
    S.text(x1, y + TH / 2, 0, f'{rest:.0f}% restate', ha='left', va='center', dx=34, fontsize=FB, color='blue@s0.2', fontweight='bold')
for t in (-15, -10, -5, 0, 5, 10, 15):
    u = (t - lags.min()) / (lags.max() - lags.min()) * XW
    S.text(u, -0.2, 0, f'{t:+d}' if t else '0', ha='center', va='top', dy=-3, fontsize=FT)
S.text(xz, ny * DY + 0.1, 24 * ZS, 'year of the first genome-wide report', ha='right', va='center', dx=-16, fontsize=FB, color='grey_dark', bbox=WB)
S.text(xs[1], ny * DY + 0.1, 11.5 * ZS, 'the claim came first', ha='left', va='center', fontsize=FB, color='vermilion', fontweight='bold', bbox=WB)
S.text(xs[24], ny * DY + 0.1, 17 * ZS, 'the claim restates the report', ha='left', va='center', fontsize=FB, color='blue@s0.2', fontweight='bold', bbox=WB)
S.text((x0 + x1) / 2, -0.2, 0, 'Years from the first genome-wide report of the association to the claim', ha='center', va='top', dy=-13,
       fontsize=FB, color='grey_dark')
r = C.loc['C3']; F.check(r.ci_low > 1 and shares[0] > shares[-1], 'a: restatement share rises with publication year')
_unused = (f'OR {r.estimate:.2f} ({r.ci_low:.2f}–{r.ci_high:.2f}) per 5 years\nthat a supported claim restates\nan earlier report; '
       f'{int(r.n):,} claims')
ax.text(0.985, 0.585, f'OR {r.estimate:.2f} ({r.ci_low:.2f}–{r.ci_high:.2f}) per 5 years that a\nsupported claim restates an earlier report\n({int(r.n):,} claims)', transform=ax.transAxes, ha='right', va='center', fontsize=FB, linespacing=1.15)
S.render(pad=(0.02, 0.02, 0.02, 0.02))
# the two ends of the story as key numbers: the share of supported claims that restated an earlier genome-wide report
fk.key(ax, 0.985, 0.90, f'{shares[-1]:.0f}%', f'restating\nin {periods[0][0]}–{str(periods[0][1])[2:]}', ha='right', transform=ax.transAxes, color=fk.c['blue@s0.2'])
fk.key(ax, 0.985, 0.72, f'{shares[0]:.0f}%', f'in {periods[-1][0]}–{str(periods[-1][1])[2:]}', ha='right', transform=ax.transAxes, color=fk.c['blue@s0.2'])

# ------------------------------------------------------------------ b  barcodes: every paper on the 20 most-claimed supported pairs
ax = fig.add_subplot(gs[1, :]); fk.panel(ax, 'b', '')
pairs = (d[d.replicated].groupby(['rsid', 'trait_key']).agg(n=('pmid', 'nunique'), fh=('first_hit_year', 'min'),
         tn=('efo', lambda s: F.trait_label(s.iloc[0]))).sort_values('n', ascending=False).reset_index()
         .drop_duplicates('rsid').head(22).sort_values('fh', ascending=False).reset_index(drop=True))
gl = F.gene_label(pairs.rsid)
gl.update({17782313: 'MC4R', 2046210: 'CCDC170', 35705950: 'MUC5B', 3803662: 'TOX3', 6983267: 'CCAT2'})   # the locus names used in the GWAS literature
rng = np.random.default_rng(3)
nb = na = 0
for k, pr in enumerate(pairs.itertuples()):
    x = d[(d.rsid == pr.rsid) & (d.trait_key == pr.trait_key)].drop_duplicates('pmid')
    xs_ = x.pub_year.values + rng.uniform(-0.42, 0.42, len(x))
    bef = x.pub_year.values < pr.fh; nb += bef.sum(); na += (~bef).sum()
    ax.vlines(xs_[~bef], k - 0.42, k + 0.42, color=BLU, lw=0.8, alpha=0.8)
    ax.vlines(xs_[bef], k - 0.42, k + 0.42, color=VER, lw=1.0)
    ax.scatter([pr.fh - 0.5], [k], marker='D', s=9, color='white', edgecolor=fk.c.ink, lw=0.6, zorder=6)
    tr = ax.text(2003.8, k, f'rs{pr.rsid}', ha='right', va='center', fontsize=FT)
    ax.annotate(gl[pr.rsid], xy=(0, 0.5), xycoords=tr, xytext=(-2.5, 0), textcoords='offset points', ha='right', va='center',
                fontsize=FT, fontstyle='italic')
    ax.text(2025.9, k, pr.tn, ha='left', va='center', fontsize=FT, color=GD)
    ax.text(2034.6, k, f'{pr.n}', ha='right', va='center', fontsize=FT, color=GD)
ax.set_xlim(1999.6, 2025.7); ax.set_ylim(-0.7, len(pairs) + 0.9); ax.set_yticks([]); ax.set_xticks([2005, 2010, 2015, 2020, 2025])
for s_ in ('top', 'right', 'left'): ax.spines[s_].set_visible(False)
ax.spines['bottom'].set_bounds(2004.0, 2025.7)
ax.text(2004.0, len(pairs) - 0.05, f'one tick per paper: {na:,} after the report (teal), {nb:,} before it (orange); open diamond, first genome-wide report',
        ha='left', va='bottom', fontsize=FT, color=GD)
ax.set_xlabel('Publication year')

# ------------------------------------------------------------------ c  five-year yield of new claims, single years
ax = fig.add_subplot(gs[2, 0]); fk.panel(ax, 'c', '')
def yield5(x):
    x = x[~x.rep_prior].copy(); x['y5'] = x.replicated & (x.first_hit_year - x.pub_year <= 5); return x
rows = []
for y in range(2005, 2021):
    xp = yield5(P[P.pub_year == y]); xn = yield5(N[N.pub_year == y])
    es, lo, hi = F.boot_share(xp, 'y5'); rows.append((y, es, lo, hi, xn.y5.mean(), xp.null_share.mean(), len(xp)))
T = pd.DataFrame(rows, columns=['y', 'pos', 'lo', 'hi', 'neg', 'null', 'n'])
ax2 = ax.twinx(); ax2.bar(T.y, T.n, width=0.8, color=fk.c['grey@t0.8'], lw=0, zorder=0); ax2.set_ylim(0, T.n.max() * 2.4)
ax2.set_yticks([]); ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)
ax.text(2004.6, 20.4, f'grey bars: new claims per year (up to {int(T.n.max()):,})', fontsize=FT, va='center')
ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
ax.fill_between(T.y, T.lo * 100, T.hi * 100, color=RL.light(RL.FIRST, 0.6), lw=0, zorder=1)
Tn = [F.boot_share(yield5(N[N.pub_year == y]), 'y5') for y in T.y]
ax.fill_between(T.y, [t[1] * 100 for t in Tn], [t[2] * 100 for t in Tn], color=fk.c['grey@t0.6'], lw=0, zorder=1)
ax.scatter(T.y, T.pos * 100, s=14, color=VER, zorder=3)
sm = T.pos.rolling(3, center=True, min_periods=2).mean()
ax.plot(T.y, sm * 100, color=VER, lw=1.2)
ax.plot(T.y, T.neg.rolling(3, center=True, min_periods=2).mean() * 100, color=GD, lw=1.0)
ax.plot(T.y, T.null * 100, color=GR, lw=0.8, ls=(0, (3, 2)))
ax.text(2020.8, 6.3, 'reported\nassociations', ha='left', va='center', fontsize=FT, color=VER, linespacing=1.0)
ax.text(2020.8, 3.6, 'no-association\nreports', ha='left', va='center', fontsize=FT, color=GD, linespacing=1.0)
ax.text(2020.8, 1.4, 'matched\nnull', ha='left', va='center', fontsize=FT, color=fk.c['grey@s0.2'], linespacing=1.0)
r = C.loc['C4']; F.check(r.ci_high < 1, 'c: five-year yield of new claims falls')
ax.text(2004.6, 18.8, 'lines: three-year means', ha='left', fontsize=FT, color=GD)
ax.text(2013.4, 18.8, f'OR {r.estimate:.2f} ({r.ci_low:.2f}–{r.ci_high:.2f}) per 5 years', fontsize=FB, color=fk.c.ink, va='center')
ax.set_xlim(2004.3, 2025.4); ax.set_ylim(0, 23.5); ax.set_yticks([0, 5, 10, 15]); ax.set_xticks([2005, 2010, 2015, 2020])
fk.chart(ax, grid='y', unit='% supported within 5 years')
fk.pin_inside(ax)

# ------------------------------------------------------------------ d  reports of no association at established variants
ax = fig.add_subplot(gs[2, 1]); fk.panel(ax, 'd', '')
Nn = N[N.pub_year.between(2005, 2025)]
rows = [(y,) + F.boot_share(x, 'rep_prior') + (len(x),) for y, x in Nn.groupby('pub_year') if len(x) >= 50]
Z = pd.DataFrame(rows, columns=['y', 'est', 'lo', 'hi', 'n'])
ax.bar(Z.y, Z.est * 100, width=0.86, color=fk.c['grey@t0.4'], lw=0, zorder=2)
ax.vlines(Z.y, Z.lo * 100, Z.hi * 100, color=fk.c['grey@s0.35'], lw=0.6, zorder=3)

Z3 = Z.set_index('y').est.rolling(3, center=True, min_periods=2).mean()
ax.plot(Z3.index, Z3.values * 100, color=GD, lw=1.6, zorder=3)
r = C.loc['C5']; F.check(r.ci_low > 1, 'd: null claims at established loci increase')
ax.text(2005.0, 23.0, f'OR {r.estimate:.2f} ({r.ci_low:.2f}–{r.ci_high:.2f}) per 5 years', fontsize=FB, color=fk.c.ink, va='center')
ax.text(2005.0, 21.2, f'{int(r.n):,} reports of no association;\nbars, single years (95% CI); line, three-year mean', ha='left', va='top', fontsize=FB,
        color=GD, linespacing=1.15)
ax.set_xlim(2004.4, 2025.6); ax.set_ylim(0, 27.5); ax.set_yticks([0, 5, 10, 15]); ax.set_xticks(range(2005, 2026, 5))
fk.chart(ax, grid='y', unit='% concerning a variant already genome-wide significant')
fk.pin_inside(ax)                                              # in-plot notes stay out of the layout computation


def blacken(fig):
    """User 2026-11-17: every non-essential text in black; only white labels on filled marks stay white."""
    import matplotlib.text, matplotlib.colors as _mc
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text() and _mc.to_rgb(t.get_color()) != (1.0, 1.0, 1.0):
            t.set_color('black')


blacken(fig)
import os as _os
if _os.environ.get('VC_AUDIT'):
    sys.path.insert(0, str(Path(__file__).resolve().parent)); import vcaudit; vcaudit.audit(fig, 'fig3')
saved = fk.save(fig, RL.OUT / 'fig3', strict=False); print(saved.report)

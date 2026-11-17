"""Figure 4 v3 (2026-11-17): ClinVar risk variants hold both histories. Rebuilt at the proportions and density of the user's main
reference PMC11106717_F5 after the user judged v2 not a reproduction. Like the reference: no panel titles; a framed overview with large
pies under name chips and a side legend (A); a column of small multiples with paired bars (B); a schematic beside a box plot with
significance brackets (C); two tiled tables, one per dataset, with rounded vertical dataset chips and header bars (D); stacked density
panels with dashed medians and significance chips (E); a 2 x 2 grid of density panels (F). The two datasets are the two splits of the
corpus (held-out papers, development papers), so D-F show each result twice.
Additional analysis of 2026-11-16 (scripts/14_clinvar.py; definitions fixed before ClinVar was read): ClinVar variant_summary of
2026-11-16, GRCh38 rows with an rsID; risk variant = clinical significance with "risk factor", "association" (not "association not
found") or "protective". Variant classes on variants with >= 5 reported associations: mostly (>= 50% supported), partly (5-50%),
rarely (< 5%). Output figures/v3/fig4."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib, matplotlib.patches as mp
from scipy.stats import gaussian_kde, mannwhitneyu
sys.path.insert(0, str(Path(__file__).resolve().parents[1])); import figdata as F
sys.path.insert(0, str(Path(__file__).resolve().parent)); import roles as RL
sys.path.insert(0, '/Volumes/ResearchData/Claude/figtaste')
from layout.canvas import figure, axes, ftext
matplotlib.rcParams['savefig.bbox'] = 'standard'

# palette sampled from PMC11106717_F5
ORG, GRN, BLU = RL.SUP, RL.PART, RL.NOT_D           # classes: mostly / partly / rarely supported
CVR, OTH = RL.CLIN, '#8b919c'                           # ClinVar risk variants / not in ClinVar
HO, DEV = '#c33c44', '#2c2e5c'                            # dataset chips: held-out / development (reference: FAFB red, hemibrain navy)
PINK = [RL.light(RL.CLIN, t) for t in (0.92, 0.75, 0.58, 0.38, 0.15)]
INK, GD, GL, GRID = '#1a1a1a', '#555555', '#c8c8c8', '#ececec'
FS, FT, FL = 6.5, 6.0, 9.5
CLS = ['mostly', 'partly', 'rarely']; CCOL = {'mostly': ORG, 'partly': GRN, 'rarely': BLU}
CLAB = {'mostly': '≥ 50%', 'partly': '5–50%', 'rarely': '< 5%'}
GROUPS = ['Risk factor', 'Association', 'Protective', 'Not in ClinVar']
OUT = RL.OUT

# ---------------------------------------------------------------- data
cv = pd.read_csv(F.ROOT / 'raw' / 'variant_summary.txt.gz', sep='\t', usecols=['ClinicalSignificance', 'RS# (dbSNP)', 'Assembly', 'ReviewStatus'], dtype=str)
cv = cv[cv.Assembly == 'GRCh38'].copy(); cv['rsid'] = pd.to_numeric(cv['RS# (dbSNP)'], errors='coerce'); cv = cv[cv.rsid > 0]
sig = cv.ClinicalSignificance.fillna('').str.lower().str.replace('association not found', '', regex=False)
cv['rf'], cv['as_'], cv['pr'] = sig.str.contains('risk factor'), sig.str.contains('association'), sig.str.contains('protective')
cv = cv[cv.rf | cv.as_ | cv.pr]
RANK = {'no assertion criteria provided': 0, 'criteria provided, single submitter': 1, 'criteria provided, conflicting classifications': 1,
        'criteria provided, multiple submitters, no conflicts': 2, 'reviewed by expert panel': 2, 'practice guideline': 2}
cv['rv'] = cv.ReviewStatus.map(RANK).fillna(0)
R = cv.groupby('rsid').agg(rf=('rf', 'any'), as_=('as_', 'any'), pr=('pr', 'any'), rv=('rv', 'max')).reset_index(); R['rsid'] = R.rsid.astype('int64')
R['assert'] = np.where(R.rf, 'Risk factor', np.where(R.as_, 'Association', 'Protective'))


def split(name):
    d = F.claims(name); d['grp'] = F.area_group(d.area)
    pos = d[d.status == 'positive'].merge(R[['rsid', 'assert', 'rv']], on='rsid', how='left')
    pos['cv'] = pos['assert'].notna()
    V = pos.groupby('rsid').agg(n=('pmid', 'size'), sup=('replicated', 'mean'), grp=('grp', lambda s: s.mode().iat[0])).reset_index()
    V = V.merge(R[['rsid', 'assert', 'rv']], on='rsid', how='left'); V['cv'] = V['assert'].notna()
    V = V[V.n >= 5].copy()
    V['cls'] = np.where(V.sup >= 0.5, 'mostly', np.where(V.sup >= 0.05, 'partly', 'rarely'))
    V['g'] = V['assert'].fillna('Not in ClinVar')
    return pos, V


posH, VH = split('heldout'); posD, VD = split('dev')
F.check(VH.cv.sum() > 100, f'{int(VH.cv.sum())} held-out and {int(VD.cv.sum())} development ClinVar risk variants with >= 5 reported associations')

fig = figure(180, 167)


def letter(x, y, s):
    ftext(fig, x, y, s, fontsize=FL, fontweight='bold', va='top')


def style(ax):
    for s_ in ('top', 'right'): ax.spines[s_].set_visible(False)
    ax.tick_params(labelsize=FT, length=2, width=0.5, pad=1.5)
    for s_ in ('left', 'bottom'): ax.spines[s_].set_linewidth(0.6)


def chip(ax, x, y, text, col, size=FT, ha='center', rot=0, pad=0.3, tr=None):
    """Plain black label (user 2026-11-17: no capsule backgrounds behind text)."""
    return ax.text(x, y, text, ha=ha, va='center', fontsize=size, color=INK, fontweight='bold', rotation=rot, transform=tr or ax.transAxes)


def bkde(v, xs, bw):
    v = np.asarray(v, float); k = gaussian_kde(np.concatenate([v, -v, 2 - v]), bw_method=bw)   # reflected at 0 and 1
    y = k(xs) * 3; return y / y.max()


# ------------------------------------------------------------------ A  framed overview
letter(0, 0.5, 'a')
fr = axes(fig, 3, 4, 128, 55); fr.set_xlim(0, 1); fr.set_ylim(0, 1); fr.set_xticks([]); fr.set_yticks([])
for s_ in fr.spines.values(): s_.set_visible(True); s_.set_linewidth(0.6); s_.set_color(INK)
fr.text(0.015, 0.95, 'variant class', fontsize=FS, fontweight='bold', va='top')
fr.text(0.015, 0.88, '(share of its reported\nassociations supported\ngenome-wide)', fontsize=FT - 0.5, color=GD, va='top', linespacing=1.05)
for k, c in enumerate(CLS):
    fr.add_patch(mp.Rectangle((0.018, 0.62 - k * 0.09), 0.026, 0.06, color=CCOL[c], lw=0))
    fr.text(0.052, 0.65 - k * 0.09, CLAB[c] + (' supported' if c != 'partly' else ''), fontsize=FS, va='center')
fr.text(0.015, 0.30, 'ClinVar assertion', fontsize=FS, fontweight='bold', va='top')
fr.text(0.015, 0.23, 'held-out papers;\nvariants with ≥ 5\nreported associations', fontsize=FT - 0.5, color=GD, va='top', linespacing=1.05)
for j, gname in enumerate(GROUPS):
    x0 = 0.215 + j * 0.198
    q = VH[VH.g == gname]; cnt = q.cls.value_counts().reindex(CLS, fill_value=0)
    if j == 3: fr.plot([x0 - 0.012] * 2, [0.05, 0.95], color=GL, lw=0.7, ls=(0, (3, 2)))
    chip(fr, x0 + 0.09, 0.90, gname, CVR if j < 3 else OTH, size=FS)
    pa = fr.inset_axes([x0 - 0.005, 0.11, 0.19, 0.72])
    pa.pie(cnt.values, colors=[CCOL[c] for c in CLS], startangle=90, counterclock=False, wedgeprops=dict(lw=0.7, edgecolor='white'))
    pa.set_aspect('equal')
    for k, (c, v) in enumerate(cnt.items()):
        ang = np.deg2rad(90 - (cnt.values[:k].sum() + v / 2) / cnt.sum() * 360)
        pa.text(0.55 * np.cos(ang), 0.55 * np.sin(ang), (f'{v / cnt.sum() * 100:.0f}%\n({v})' if v / cnt.sum() >= 0.16 else f'{v / cnt.sum() * 100:.0f}%'), ha='center', va='center', fontsize=FT - 1, linespacing=1.0, clip_on=False,
                color='white', fontweight='bold')
    gl = F.gene_label(q.sort_values('n', ascending=False).rsid.head(40).tolist())
    for k, c in enumerate(('mostly', 'rarely')):
        top = q[q.cls == c].sort_values('n', ascending=False).head(3)
        fr.text(x0 + 0.09, 0.085 - k * 0.055, ', '.join(dict.fromkeys(g for g in (gl.get(r, '') for r in top.rsid) if g)), ha='center',
                fontsize=FT - 0.5, fontstyle='italic', color=CCOL[c])
    fr.text(x0 + 0.09, 0.825, f'n = {len(q):,}', ha='center', fontsize=FT, color=GD)

# ------------------------------------------------------------------ B  small multiples by trait area
letter(134, 0.5, 'b')
AR = [a for a in F.AREA_ORDER if a not in ('Not mapped',)]
for i, a in enumerate(AR):
    sb = axes(fig, 137, 4.5 + i * 6.6, 25, 5.9)
    for k, c in enumerate(CLS):
        for m, (cvf, alpha) in enumerate(((True, 1.0), (False, 0.45))):
            q = VH[(VH.grp == a) & (VH.cv == cvf)]
            v = (q.cls == c).mean() * 100 if len(q) else np.nan
            sb.bar(k * 2.6 + m, v, width=0.9, color=CCOL[c], alpha=alpha, lw=0)
    sb.set_xlim(-0.7, 7.5); sb.set_ylim(0, 100); sb.set_yticks([0, 100] if i == len(AR) - 1 else []); sb.set_xticks([])
    style(sb); sb.tick_params(labelsize=FT - 1, length=1.5)
    sb.text(7.6, 70, F.AREA_SHORT.get(a, a).replace(' disease', ''), fontsize=FT - 0.5, ha='left', va='center')
bl = axes(fig, 137, 58.6, 40, 3); bl.axis('off'); bl.set_xlim(0, 1); bl.set_ylim(0, 1)
bl.add_patch(mp.Rectangle((0.0, 0.25), 0.06, 0.5, color=GD, lw=0)); bl.text(0.08, 0.5, 'ClinVar risk', fontsize=FT - 0.5, va='center')
bl.add_patch(mp.Rectangle((0.48, 0.25), 0.06, 0.5, color=GD, alpha=0.45, lw=0)); bl.text(0.56, 0.5, 'not in ClinVar', fontsize=FT - 0.5, va='center')
ftext(fig, 137, 2.2, '% of variants in each class', fontsize=FT - 0.5, color=GD)

# ------------------------------------------------------------------ C  schematic + box plot with brackets
letter(0, 63, 'c')
sc = axes(fig, 2, 67, 40, 37); sc.set_xlim(0, 1); sc.set_ylim(0, 1); sc.axis('off')
sc.plot([0.05, 0.95], [0.86, 0.86], color=GD, lw=1.4, solid_capstyle='round')
for xx in np.linspace(0.08, 0.92, 22): sc.plot([xx, xx], [0.83, 0.89], color=GL, lw=0.5)
sc.scatter([0.5], [0.86], s=26, color=CVR, zorder=3, lw=0); sc.text(0.5, 0.95, 'rs1800795', ha='center', fontsize=FT, color=CVR, fontweight='bold')
for x0, title, sub, col in ((0.01, 'literature', '225 papers,\n138 traits', '#8e7bb5'), (0.345, 'ClinVar', 'risk factor', CVR), (0.68, 'GWAS', 'none at\nP < 5e-8', DEV)):
    sc.add_patch(mp.FancyBboxPatch((x0, 0.30), 0.30, 0.32, boxstyle='round,pad=0.01,rounding_size=0.04', fc='white', ec=col, lw=0.9))
    sc.text(x0 + 0.15, 0.53, title, ha='center', fontsize=FT - 0.5, fontweight='bold', color=col)
    sc.text(x0 + 0.15, 0.40, sub, ha='center', va='center', fontsize=FT - 1.5, color=GD, linespacing=1.0)
    sc.annotate('', (x0 + 0.15, 0.63), (0.5, 0.82), arrowprops=dict(arrowstyle='-|>', lw=0.5, color=GD, mutation_scale=5))
sc.text(0.5, 0.17, 'one IL6 promoter variant:\nlisted as a risk factor,\n3% of its associations supported', ha='center', va='center', fontsize=FT - 1, color=GD, linespacing=1.05)
bx = axes(fig, 52, 70, 18, 30)
q = VH[VH.cv]; data = [q[q.rv == k].sup.values * 100 for k in range(3)]
bp = bx.boxplot(data, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=0.8), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
cols = [PINK[1], PINK[2], PINK[3]]
for b_, c_ in zip(bp['boxes'], cols): b_.set(facecolor=c_, edgecolor=CVR, lw=0.7)
rng = np.random.default_rng(2)
for k, dd in enumerate(data):
    bx.scatter(k + 1 + rng.uniform(-0.2, 0.2, len(dd)), dd, s=2, color=CVR, alpha=0.7, lw=0, zorder=3)
    bx.text(k + 1, -14, f'n={len(dd)}', ha='center', fontsize=FT - 1, color=GD)
for (i, j, y) in ((0, 1, 108), (0, 2, 121)):
    p = mannwhitneyu(data[i], data[j]).pvalue; st = '****' if p < 1e-4 else '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 0.05 else 'n.s.'
    bx.plot([i + 1, i + 1, j + 1, j + 1], [y - 3, y, y, y - 3], color=INK, lw=0.5); bx.text((i + j) / 2 + 1, y + 1, st, ha='center', fontsize=FT - 0.5)
bx.set_xticks([1, 2, 3]); bx.set_xticklabels(['none', 'single', 'multi'], fontsize=FT - 1); bx.set_ylim(-20, 130); bx.set_yticks([0, 50, 100])
bx.set_ylabel('% supported (per variant)', fontsize=FT, labelpad=1); bx.set_xlabel('review criteria of the assertion', fontsize=FT - 0.5, labelpad=6); style(bx)
bx.spines['left'].set_bounds(0, 100)

# ------------------------------------------------------------------ D  two tiled tables, one per dataset
letter(73, 63, 'd')
for b_i, (V, name, col) in enumerate(((VH, 'held-out', HO), (VD, 'development', DEV))):
    X0 = 84 + b_i * 44
    T = pd.DataFrame(index=GROUPS[:3], columns=CLS, dtype=float)
    for gname in GROUPS[:3]:
        q = V[V.g == gname]; w = q.n; tot = w.sum()
        for c in CLS: T.loc[gname, c] = w[q.cls == c].sum() / tot * 100
    allq = V[V.cv]
    hb = axes(fig, X0, 69, 33, 3.6); hb.axis('off'); hb.set_xlim(0, 3); hb.set_ylim(0, 1)
    for k, c in enumerate(CLS):
        sh = allq.n[allq.cls == c].sum() / allq.n.sum()
        hb.add_patch(mp.Rectangle((k + 0.02, 0), 0.96, 1, color=['#e3e3e3', '#cdcdcd', '#b3b3b3'][k], lw=0))
        hb.text(k + 0.5, 0.5, f'{sh * 100:.0f}%', ha='center', va='center', fontsize=FS)
    tb = axes(fig, X0, 73.4, 33, 21); tb.set_xlim(0, 3); tb.set_ylim(0, 3); tb.invert_yaxis(); tb.axis('off')
    for i_, gname in enumerate(GROUPS[:3]):
        for k, c in enumerate(CLS):
            v = T.loc[gname, c]; fc = PINK[min(int(v / 15), 4)]
            tb.add_patch(mp.Rectangle((k, i_), 1, 1, fc=fc, ec='#a0a0a0', lw=0.4))
            tb.text(k + 0.5, i_ + 0.5, f'{v:.1f}%', ha='center', va='center', fontsize=FS)
        if b_i == 1: tb.text(3.08, i_ + 0.5, gname, ha='left', va='center', fontsize=FS)
    for k, c in enumerate(CLS):
        tb.text(k + 0.75, 3.1, CLAB[c], ha='right', va='top', rotation=40, fontsize=FT)
    vc = axes(fig, X0 - 7.5, 70, 5, 24); vc.axis('off')
    chip(vc, 0.5, 0.5, name, col, size=FS, rot=90, pad=0.45)
ftext(fig, 84, 103.5, 'claims on variants of each class, by ClinVar assertion; header: all risk variants', fontsize=FT - 0.5, color=GD)

# ------------------------------------------------------------------ E  stacked densities of the variant-level share
letter(0, 108, 'e')
xs = np.linspace(0, 1, 300)
for r_, (V, name, col) in enumerate(((VH, 'held-out', HO), (VD, 'development', DEV))):
    ea = axes(fig, 16, 113 + r_ * 25, 60, 19)
    for lab, msk, c_ in (('ClinVar risk', V.cv, CVR), ('not in ClinVar', ~V.cv, OTH)):
        yk = bkde(V[msk].sup, xs, 0.08); ea.fill_between(xs * 100, 0, yk, color=c_, alpha=0.18, lw=0); ea.plot(xs * 100, yk, color=c_, lw=1.2)
        h_, b_ = np.histogram(V[msk].sup, bins=np.linspace(0, 1, 21)); ea.bar(b_[:-1] * 100 + (1.2 if c_ == CVR else 2.6), h_ / h_.max(), width=1.3, color=c_, alpha=0.55, lw=0, align='edge')
        ea.axvline(np.median(V[msk].sup) * 100, color=c_, lw=0.7, ls=(0, (3, 2)))
    p = mannwhitneyu(V[V.cv].sup, V[~V.cv].sup).pvalue
    ea.text(78, 0.82, '****' if p < 1e-4 else f'P = {p:.2g}', fontsize=FS, ha='center')
    ea.add_patch(mp.Rectangle((70, 0.62), 5, 0.12, color=CVR, lw=0)); ea.text(76, 0.68, 'vs', fontsize=FT - 1, va='center'); ea.add_patch(mp.Rectangle((81, 0.62), 5, 0.12, color=OTH, lw=0))
    ea.set_xlim(0, 100); ea.set_ylim(0, 1.08); ea.set_yticks([0, 0.5, 1]); style(ea)
    if r_ == 1: ea.set_xlabel("% of the variant's reported associations supported", fontsize=FT, labelpad=1)
    else: ea.set_xticklabels([])
    ce = axes(fig, 2, 113 + r_ * 25, 5, 19); ce.axis('off'); chip(ce, 0.5, 0.5, name, col, size=FT, rot=90, pad=0.4)
ftext(fig, 10.2, 125, 'max-normalised density', fontsize=FT, rotation=90, va='center', ha='center')
lg = axes(fig, 78, 113, 18, 10); lg.axis('off'); lg.set_xlim(0, 1); lg.set_ylim(0, 1)
for k, (lab, c_) in enumerate((('ClinVar risk\nvariants', CVR), ('not in\nClinVar', OTH))):
    lg.plot([0, 0.18], [0.82 - k * 0.5] * 2, color=c_, lw=1.3); lg.text(0.24, 0.82 - k * 0.5, lab, fontsize=FT - 0.5, va='center', linespacing=1.0)
ftext(fig, 78, 126, 'dashed: medians', fontsize=FT - 1, color=GD)

# ------------------------------------------------------------------ F  2 x 2 densities of publication year
letter(98, 108, 'f')
ys = np.linspace(2000, 2025, 300)
for r_, (pos, V, name, col) in enumerate(((posH, VH, 'held-out', HO), (posD, VD, 'development', DEV))):
    for c_i, gname in enumerate(('Risk factor', 'Association')):
        fa = axes(fig, 108 + c_i * 37, 113 + r_ * 25, 31, 19)
        rows = pos.merge(V[V.g == gname][['rsid', 'cls']], on='rsid')
        for c in ('mostly', 'rarely'):
            yy = rows[rows.cls == c].pub_year.values.astype(float)
            if len(yy) < 10: continue
            k = gaussian_kde(yy, bw_method=0.25 if len(yy) > 200 else 0.45)(ys); fa.fill_between(ys, 0, k / k.max(), color=CCOL[c], alpha=0.2, lw=0); fa.plot(ys, k / k.max(), color=CCOL[c], lw=1.2)
            hy = np.bincount(yy.astype(int) - 2000, minlength=26)[:26]; fa.bar(np.arange(2000, 2026) + (-0.2 if c == 'mostly' else 0.2), hy / hy.max() * 0.35, width=0.4, color=CCOL[c], alpha=0.6, lw=0)
            fa.axvline(np.median(yy), color=CCOL[c], lw=0.7, ls=(0, (3, 2)))
        p = mannwhitneyu(rows[rows.cls == 'mostly'].pub_year, rows[rows.cls == 'rarely'].pub_year).pvalue
        fa.text(2002, 0.9, '****' if p < 1e-4 else '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 0.05 else 'n.s.', fontsize=FS)
        fa.set_xlim(2000, 2025); fa.set_ylim(0, 1.1); fa.set_yticks([0, 0.5, 1]); fa.set_xticks([2005, 2015]); style(fa)
        if r_ == 0: fa.set_title(gname, fontsize=FS, pad=2, color=CVR)
        if r_ == 0: fa.set_xticklabels([])
        if c_i == 1: fa.set_yticklabels([])
    cf = axes(fig, 99.5, 113 + r_ * 25, 5, 19); cf.axis('off'); chip(cf, 0.5, 0.5, name, col, size=FT, rot=90, pad=0.4)
ftext(fig, 142, 165.0, 'publication year of the claims', fontsize=FT, ha='center', va='bottom')
lk = axes(fig, 104, 105.5, 60, 2.5); lk.axis('off'); lk.set_xlim(0, 1); lk.set_ylim(0, 1)
for k, c in enumerate(('mostly', 'rarely')):
    lk.plot([k * 0.52, k * 0.52 + 0.06], [0.5, 0.5], color=CCOL[c], lw=1.3); lk.text(k * 0.52 + 0.08, 0.5, c + ' supported variants', fontsize=FT - 0.5, va='center')


def blacken(fig):
    """User 2026-11-17: every non-essential text in black; only white labels on filled marks stay white."""
    import matplotlib.text, matplotlib.colors as _mc
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text() and _mc.to_rgb(t.get_color()) != (1.0, 1.0, 1.0):
            t.set_color('black')


blacken(fig)
import os as _os
if _os.environ.get('VC_AUDIT'):
    sys.path.insert(0, str(Path(__file__).resolve().parent)); import vcaudit; vcaudit.audit(fig, 'fig4')
fig.savefig(OUT / 'fig4.png', dpi=600, facecolor='white'); fig.savefig(OUT / 'fig4.pdf', facecolor='white'); print('wrote', OUT / 'fig4')

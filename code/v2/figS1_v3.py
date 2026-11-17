"""Supplementary Figure 1 v3 (2026-11-17): the AJHG figure S1 redrawn to the rules of the v3 main figures: lower-case letters
without title lines (statements in the legend), colour roles of roles.py (palette SF), all non-essential text black, no text on
filled badges, text-on-mark audit (VC_AUDIT=1). Output figures/colour/<PAL>/figS1.
Supplementary Figure 1. The five primary estimates under the prespecified sensitivity analyses: all claims; without outcomes linked
only because the abstract tagged a single disease; without traits placed through a broader ontology ancestor; populations stated in the
abstract; without retracted papers; one weight per paper; and, as an additional analysis added after the plan (open circles, below the
dotted line), without claims from genome-wide studies and meta-analyses. One column per estimate, one row per analysis; the band is the
95% interval of the primary estimate, the dashed line the value of no effect. Stage 1 draws the journal-neutral master (figures/master);
`python3 -m figkit fit ajhg` re-renders it for the journal. Usage: python 45_figS1_sensitivity.py [dev|heldout]."""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
import figkit as fk
sys.path.insert(0, str(Path(__file__).resolve().parents[1])); import figdata as F
sys.path.insert(0, str(Path(__file__).resolve().parent)); import roles as RL
SPLIT = sys.argv[1] if len(sys.argv) > 1 else 'heldout'
sfx = '' if SPLIT == 'heldout' else '_dev'
e = fk.enc('varclaims')
VER, GD, GR, GL, WH = RL.SUP, fk.c.grey_dark, fk.c.grey, fk.c.grey_light, fk.c.white
GN = fk.c['grey@s0.2']
FB, FT = fk.pt('base'), fk.pt('tick')
with F.con() as c:
    pr = c.sql(f"select * from confirm_primary{sfx}").df().assign(set='primary')
    se = c.sql(f"select * from confirm_sensitivity{sfx}").df()
    ad = c.sql("select * from additional_sensitivity").df() if SPLIT == 'heldout' else se.iloc[:0]   # scripts/13_additional.py
    tr = c.sql("select * from additional_trait_rule").df() if SPLIT == 'heldout' else se.iloc[:0]   # scripts/16_robustness.py
ad = pd.concat([ad, tr[[c_ for c_ in ad.columns if c_ in tr.columns]]])
D = pd.concat([pr, se, ad])
names = {'primary': 'All claims', 'a_no_link_d': 'Outcome linked by name or ontology', 'b_no_mesh_ancestor': 'Trait mapped without broadening',
         'c_population_stated': 'Population stated in the abstract', 'd_no_retracted': 'Retracted papers removed', 'e_paper_weighted': 'One weight per paper'}
if len(ad):
    names['x_no_gwas_meta'] = 'Genome-wide studies and\nmeta-analyses excluded'
    names['exact'] = 'Same ontology term'
    names['one_step'] = 'Ontology terms one step apart'
NADD = 3 if len(ad) else 0                                        # rows of additional analyses, drawn last with open circles
tests = {'C1': ('a', '', 'Excess support over the\nmatched null (points)', 0, 100, [0, 10, 20]),
         'C2': ('b', '', 'Study size\n(OR per 10-fold sample size)', 1, 1, [1, 1.5, 2]),
         'C3': ('c', '', 'Restatement\n(OR per 5 years)', 1, 1, [1, 1.5, 2, 2.5]),
         'C4': ('d', '', 'Five-year yield of new claims\n(OR per 5 years)', 1, 1, [0.7, 0.8, 0.9, 1]),
         'C5': ('e', '', 'Discordant reports\n(OR per 5 years)', 1, 1, [1, 1.2, 1.4])}

fig = fk.figure('double', height=100, kind='supp', profile=fk.target())
gs = fig.add_gridspec(1, 5)
yy = np.arange(len(names))[::-1]
for k, (t, (let, title, xl, ref, mult, ticks)) in enumerate(tests.items()):
    ax = fig.add_subplot(gs[0, k]); fk.panel(ax, let, title)
    x = D[D.test == t].set_index('set').reindex(list(names))
    est, lo, hi = x.estimate.values * mult, x.ci_low.values * mult, x.ci_high.values * mult
    F.check(bool(((lo > ref) if est[0] > ref else (hi < ref)).all()), f'{let}: every sensitivity estimate on the side of the primary estimate')
    ax.axvspan(lo[0], hi[0], color=RL.light(RL.SUP, 0.82), lw=0, zorder=1)          # the primary estimate's interval, through all rows
    ax.axvline(ref, color=GR, lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.hlines(yy, lo, hi, color=[VER] + [GD] * (len(yy) - 1), lw=1.2, zorder=3)
    m = len(yy) - NADD                                          # planned rows: filled circles; additional rows: open circles
    ax.scatter(est[:m], yy[:m], s=[30] + [20] * (m - 1), color=[VER] + [GD] * (m - 1), lw=0.6, edgecolor=WH, zorder=4)
    if NADD:
        ax.scatter(est[m:], yy[m:], s=20, facecolor=WH, edgecolor=GD, lw=1.0, zorder=4)
        ax.axhline(yy[m] + 0.5, color=GR, lw=0.7, ls=(0, (1, 2)), zorder=2)
    if ref == 1:
        ax.set_xscale('log')
    span = (min(ticks), max(max(ticks), hi.max() * 1.02)) if est[0] > ref else (min(min(ticks), lo.min() * 0.985), max(ticks))
    ax.set_xlim(span[0] - (0.6 if ref == 0 else 0), span[1] * (1.22 if ref == 1 and est[0] > ref else 1.0) + (5.5 if ref == 0 else 0))
    for yi, v, h_, l_ in zip(yy, est, hi, lo):                  # the estimate printed beside its interval, away from the reference
        right = est[0] > ref
        ax.annotate(f'{v:.1f}' if ref == 0 else f'{v:.2f}', (max(h_, hi[0]) if right else min(l_, lo[0]), yi), (3 if right else -3, 0), textcoords='offset points',
                    ha='left' if right else 'right', va='center', fontsize=FT, color=VER if yi == yy[0] else GD,
                    fontweight='bold' if yi == yy[0] else 'normal')
    ax.set_xticks(ticks); ax.set_xticklabels([f'{v:g}' for v in ticks]); ax.minorticks_off()
    ax.set_yticks(yy); ax.set_yticklabels(list(names.values()) if k == 0 else [], fontsize=FT)
    ax.tick_params(axis='y', length=0)
    if k == 0:
        ax.get_yticklabels()[0].set_color(VER); ax.get_yticklabels()[0].set_fontweight('bold')
    ax.set_ylim(-0.6, len(names) - 0.4)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel(xl, fontsize=FT, color=GD, linespacing=1.1)
    ax._figkit_stat = True                                      # estimates and intervals are printed in the panel


def blacken(fig):
    """User 2026-11-17: every non-essential text in black; only white labels on filled marks stay white."""
    import matplotlib.text, matplotlib.colors as _mc
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text() and _mc.to_rgb(t.get_color()) != (1.0, 1.0, 1.0):
            t.set_color('black')


blacken(fig)
if os.environ.get('VC_AUDIT'):
    import vcaudit; vcaudit.audit(fig, 'figS1')
saved = fk.save(fig, RL.OUT / 'figS1', strict=False); print(saved.report)

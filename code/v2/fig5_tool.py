"""Figure 5 (CSBJ revision, 2026-11-21): the web tool (version 1.2), as screenshots of the live site (light theme, 2x pixel density,
1,048-pixel window) captured from https://zxu69461-creator.github.io/variant-claim-atlas/ on 2026-11-21: (a) look up rs1800795;
(b) check a claim, rs1800795 and type 2 diabetes; (c) screen a paper, PMID 42749416; (d) batch check of the example list; (e) ClinVar
risk variants filtered to the rarely supported class; (f) programmatic access on the help page.
Screenshots in figures/tool/v12/. Output RL.OUT/fig5 (FIG_OUT=figures/csbj, PANEL_PAREN=1)."""
import os, sys
from pathlib import Path
import matplotlib, numpy as np
from PIL import Image
import figkit as fk
sys.path.insert(0, str(Path(__file__).resolve().parent)); import roles as RL
matplotlib.rcParams['savefig.bbox'] = 'standard'
D = Path(__file__).resolve().parents[2] / 'figures' / 'tool' / 'v12'
names = ['a_lookup.png', 'b_check.png', 'c_paper.png', 'd_batch.png', 'e_clinvar.png', 'f_help.png']
H = 800                                                     # crop every screenshot to the same CSS height (2x pixels)
ims = []
for n in names:
    im = Image.open(D / n).convert('RGB'); w, h = im.size
    im = im.crop((0, 0, w, min(h, 2 * H)))
    if im.size[1] < 2 * H:
        bg = Image.new('RGB', (w, 2 * H), 'white'); bg.paste(im, (0, 0)); im = bg
    ims.append(np.asarray(im))
fig = fk.figure('double', height=225, profile=fk.target())
W, Hf = fig.get_size_inches()
gap, top, lab = 0.012, 0.03, 0.03
pw = (1 - 3 * gap) / 2
ph = pw * W * ims[0].shape[0] / ims[0].shape[1] / Hf
for k, im in enumerate(ims):
    c, r = k % 2, k // 2
    x0 = gap + c * (pw + gap); y0 = 1 - top - (r + 1) * ph - r * lab
    ax = fig.add_axes([x0, y0, pw, ph]); ax.imshow(im, interpolation='lanczos'); ax.set_axis_off()
    for s in ('left', 'right', 'top', 'bottom'):
        ax.spines[s].set_visible(False)
    ax.add_patch(matplotlib.patches.Rectangle((0, 0), im.shape[1] - 1, im.shape[0] - 1, fill=False, ec=fk.c['grey@t0.4'], lw=0.5))
    fig.text(x0, y0 + ph + 0.006, RL.PL('abcdef'[k]), fontsize=9, fontweight='bold', ha='left', va='bottom')
saved = fk.save(fig, RL.OUT / 'fig5', strict=False); print(saved.report)

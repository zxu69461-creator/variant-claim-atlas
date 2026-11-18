"""Supplementary Figure S1 (CSBJ version, 2026-11-17): experimental and technical design, as the journal asks for at the start of
Materials and Methods. A left-to-right flow from the PubMed baseline to the atlas, with the split into development and held-out papers
and the comparisons made on held-out papers. Numbers from manuscript/values.json. Plain outlined boxes, black text, no fills behind text.
Output RL.OUT/figS0_design (FIG_OUT=figures/csbj)."""
import json, os, sys
from pathlib import Path
import matplotlib
from matplotlib.patches import FancyArrowPatch, Rectangle
import figkit as fk
sys.path.insert(0, str(Path(__file__).resolve().parents[1])); import figdata as F
sys.path.insert(0, str(Path(__file__).resolve().parent)); import roles as RL
matplotlib.rcParams['savefig.bbox'] = 'standard'
V = json.loads((F.ROOT / 'manuscript' / 'values.json').read_text())
FB, FT = fk.pt('base'), fk.pt('tick')

fig = fk.figure('double', height=82, kind='supp', profile=fk.target())
ax = fig.add_axes([0.01, 0.02, 0.98, 0.96]); ax.set_xlim(0, 180); ax.set_ylim(0, 82); ax.set_axis_off()


def box(x, y, w, h, head, body, edge=None, lw=0.8):
    ax.add_patch(Rectangle((x, y), w, h, fc='white', ec=edge or RL.G2, lw=lw, zorder=2))
    ax.text(x + 2, y + h - 2.2, head, ha='left', va='top', fontsize=FB, fontweight='bold', zorder=3)
    ax.text(x + 2, y + h - 7.0, body, ha='left', va='top', fontsize=FT, linespacing=1.25, zorder=3)


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=7, lw=0.8, color=RL.G1, zorder=1))


# row 1: from abstracts to claim records
H1 = 21
box(2, 59, 40, H1, 'Literature corpus', f"PubMed 2026 baseline: titles\nand abstracts naming a dbSNP\nidentifier ({V['n_abstracts']})")
box(48, 59, 40, H1, 'Claim extraction', 'local language model\n(Qwen3.5-9B); each result quotes\nits sentence; only variants the\nabstract names')
box(94, 59, 40, H1, 'Claim records', f"{V['n_claims']} records from\n{V['n_papers']} papers; one per paper,\nvariant, trait and population")
box(140, 59, 38, H1, 'Trait linking', f"outcome linked to EFO\n({V['pct_linked']} of disease or trait\nclaims); population, sample size")
for a, b in ((42, 48), (88, 94), (134, 140)):
    arrow(a, 69.5, b, 69.5)

# row 2: the split, made by PMID before any comparison
ax.text(2, 54, 'split by PMID before any claim was\ncompared with genome-wide evidence', ha='left', va='center', fontsize=FT, linespacing=1.15)
box(30, 31, 64, 18, 'Development papers (1/8)', f"{V['n_dev_papers']} papers; exploration of the rules\nand estimators, then a dated analysis plan;\nresults reported alongside as a check")
box(104, 31, 74, 18, 'Held-out papers (7/8)', f"{V['n_ho_papers']} papers; blinded reading check (150 items),\nthen five primary estimates (Holm), sensitivity analyses,\nnarrower trait matching and ClinVar (16 November 2026)", edge=RL.SUP, lw=1.0)
arrow(150, 59, 62, 49)
arrow(165, 59, 165, 49)

# row 3: the reference both splits are compared with
box(30, 7, 148, 18, 'Genome-wide comparison (the same rules for every claim)', 'GWAS Catalog (29 October 2026): the same variant or an LD proxy (r$^2$ $\\geq$ 0.8, 1000 Genomes EUR or EAS) at P < 5 $\\times$ 10$^{-8}$,\nfor the same or a broader or narrower EFO term, in a paper other than the claim\'s own;\nmatched null from 20,000 random GWAS Catalog variants in the same allele-frequency bin', edge=RL.SUP, lw=1.0)
arrow(62, 25, 62, 31); arrow(141, 25, 141, 31)

# outputs
ax.text(2, 2.5, 'Outputs:', fontsize=FB, fontweight='bold', va='center')
ax.text(19, 2.5, 'the atlas (Data S1), the analysis code and dated plan (Data S2), the ClinVar table (Data S3) and the online lookup page', fontsize=FT, va='center')


def blacken(fig):
    import matplotlib.text, matplotlib.colors as _mc
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text() and _mc.to_rgb(t.get_color()) != (1.0, 1.0, 1.0):
            t.set_color('black')


blacken(fig)
if os.environ.get('VC_AUDIT'):
    import vcaudit; vcaudit.audit(fig, 'figS0_design')
saved = fk.save(fig, RL.OUT / 'figS0_design', strict=False); print(saved.report)

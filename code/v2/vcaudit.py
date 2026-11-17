"""Text audit for the varclaims v3 figures (figtaste layout audit): text on text and text on marks (lines, points, bar and band edges).
Call audit(fig, name) after the figure is drawn; writes figures/v3/<name>.audit.svg and prints every issue."""
import sys, json
from pathlib import Path
import matplotlib
sys.path.insert(0, '/Volumes/ResearchData/Claude/figtaste')
from layout.audit import mpl_boxes, check_texts, check_marks

OUT = Path('/Volumes/ResearchData/Claude/varclaims/figures/v3')


def audit(fig, name, min_px=6):
    matplotlib.rcParams['svg.fonttype'] = 'none'
    svg = OUT / f'{name}.audit.svg'
    fig.savefig(svg, facecolor='white')
    texts, axes = mpl_boxes(fig)
    W, H = fig.get_size_inches() * 25.4
    iss = check_texts(texts, {}, (0, 0, W, H))
    iss2, _ = check_marks(svg, texts, min_px=min_px)
    allx = iss + iss2
    (OUT / f'{name}.audit.json').write_text(json.dumps(allx, indent=0))
    print(f'AUDIT {name}: {len(allx)} issues'); [print('  ', x) for x in allx]
    return allx

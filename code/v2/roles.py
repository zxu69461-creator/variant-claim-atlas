"""Colour roles for the varclaims figures (palettes.json; system chosen with env PAL, output folder figures/colour/<PAL>)."""
import json, os
from pathlib import Path
import numpy as np
from matplotlib.colors import to_rgb, to_hex, LinearSegmentedColormap

HERE = Path(__file__).resolve().parent
D = json.loads((HERE / 'palettes.json').read_text())
S = os.environ.get('PAL', D['default']); R = D['systems'][S]


def mix(a, b, t):
    a, b = np.array(to_rgb(a)), np.array(to_rgb(b)); return to_hex((1 - t) * a + t * b)


def dark(c, t=0.35): return mix(c, '#000000', t)
def light(c, t=0.6): return mix(c, '#FFFFFF', t)
def seq(c, lo=0.9): return LinearSegmentedColormap.from_list('s', [light(c, lo), c, dark(c, 0.45)])


SUP, NOT, PART, REST, FIRST, CLIN = (R[k] for k in ('cSup', 'cNot', 'cPart', 'cRest', 'cFirst', 'cClin'))
NOT_D = dark(NOT, 0.45)                                   # the 'not supported' role as a line or outline
G1, G2, G3, G4 = '#4A4440', '#7D746B', '#A69C91', '#CFC7BD'  # warm greys for context variables (populations, designs, chromosomes)
NULL = '#8A8A8A'
OUT = Path(os.environ['FIG_OUT']) if os.environ.get('FIG_OUT') else Path('/Volumes/ResearchData/Claude/varclaims/figures/colour') / S
OUT.mkdir(parents=True, exist_ok=True)


def PL(s):
    """Panel letter as the target journal writes it: '(a)' when PANEL_PAREN=1 (CSBJ), else 'a'."""
    return f'({s})' if os.environ.get('PANEL_PAREN') else s

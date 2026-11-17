"""Journal-neutral text convention for the v3 manuscript (2026-11-17): human gene symbols in italics everywhere after the title page
except the reference list. No MIM numbers (those were AJHG-specific, scripts/ajhg_genes.py)."""
import re

GENES = ['TCF7L2', 'FTO', 'COMT', 'TNF', 'PNPLA3', 'MTHFR', 'ALDH2', 'IL6', 'VDR', 'MIR146A', 'MIR196A2', 'APOE', 'ESR1', 'SOD2',
         'IL1B', 'NOS3', 'BDNF', 'FGFR2', 'CFH']


def apply(text, lang):
    if lang != 'en' or '## Abstract' not in text:
        return text
    head, body = text.split('## Abstract', 1)
    refs = ''
    if '## References' in body:
        body, rest = body.split('## References', 1)
        nxt = re.search(r'\n## ', rest)
        refs, after = (rest[:nxt.start()], rest[nxt.start():]) if nxt else (rest, '')
        body_after = after
    else:
        body_after = ''
    def it(t):
        for g in GENES:
            t = re.sub(rf'(?<![*\w/-]){g}(?![*\w-])', f'*{g}*', t)
        return t
    out = head + '## Abstract' + it(body)
    if refs or body_after:
        out += '## References' + refs + it(body_after)
    return out

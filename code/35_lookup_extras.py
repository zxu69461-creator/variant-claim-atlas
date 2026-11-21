"""Extra data files for the web tool (2026-11-21), next to those of 34_lookup_data.py, for screening a paper or a pasted text and for
links to the ontology. Output lookup_site/data/:
  lexicon.json  one row per linked trait: its ontology identifier and the disease or trait phrases that abstracts used for it
                (the wording read from the abstracts, lower case; a phrase is kept for a trait when it was read at least twice and
                at least 80% of its uses map to that trait), so that a trait named in new text can be matched to the atlas
  papers.json   one row per paper: PMID, year and the variants with a claim in it
Trait labels are built exactly as in 34_lookup_data.py, so rows join on the label."""
import json, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent)); import figdata as F
OUT = F.ROOT / 'lookup_site' / 'data'

with F.con() as c:
    d = c.sql("""select v.pmid, v.rsid, v.efo, v.outcome, v.pub_year from variant_claims v
                 where v.outcome_type <> 'drug_response' and v.status in ('positive', 'null')""").df()
d['trait'] = [F.trait_label(e, o[:80] if isinstance(o, str) else '') for e, o in zip(d.efo, d.outcome)]
d['efo1'] = [next((t for t in (e if e is not None else []) if t in F.efo_names()), '') for e in d.efo]

norm = lambda s: re.sub(r'\s+', ' ', re.sub(r'[^\w\s\'-]', ' ', str(s).lower())).strip()   # noqa: E731
d['phrase'] = d.outcome.map(norm)
d = d[d.trait != '']
use = d.groupby(['phrase', 'trait']).size().rename('k').reset_index()
tot = use.groupby('phrase').k.transform('sum')
keep = use[(use.k >= 2) & (use.k / tot >= 0.8) & (use.phrase.str.len() >= 3) & ~use.phrase.str.fullmatch(r'[\d\s\'-]*')]
GENERIC = set('disease diseases disorder disorders risk susceptibility level levels serum plasma concentration concentrations '
               'syndrome type with and the of in to for onset early late development presence severity measurement trait'.split())
words = lambda s: {w for w in re.split(r"[\s'-]+", s) if len(w) >= 4 and w not in GENERIC}   # noqa: E731
abbr = lambda p: ' ' not in p and len(p) <= 6 and re.search('[a-z]', p) is not None             # noqa: E731
keep = keep.assign(gk=keep.groupby('phrase').k.transform('sum'))
lex = []
for t, x in d.groupby('trait'):
    lab = norm(t); lw = words(lab)
    k = keep[keep.trait == t].sort_values('k', ascending=False)
    # a phrase is kept when it shares a content word with the trait name, or when it is an abbreviation read at least three times
    # (abbreviations are matched in new text only where written in capitals)
    ph = [p for p, kk in zip(k.phrase, k.k) if p != lab and ((lw & words(p)) or (abbr(p) and kk >= 3))]
    ph = [lab] + ph[:24]
    lex.append([t, x.efo1.mode().iloc[0] if (x.efo1 != '').any() else '', ph])
(OUT / 'lexicon.json').write_text(json.dumps({'cols': ['trait', 'efo', 'phrases'], 'rows': lex}, separators=(',', ':'), allow_nan=False))

pp = [[int(p), int(x.pub_year.iloc[0]), sorted({int(r) for r in x.rsid})] for p, x in d.groupby('pmid')]
(OUT / 'papers.json').write_text(json.dumps({'cols': ['pmid', 'year', 'rsids'], 'rows': pp}, separators=(',', ':'), allow_nan=False))
print(len(lex), 'traits,', sum(len(r[2]) for r in lex), 'phrases;', len(pp), 'papers;',
      (OUT / 'lexicon.json').stat().st_size // 1000, 'kB,', (OUT / 'papers.json').stat().st_size // 1000, 'kB')

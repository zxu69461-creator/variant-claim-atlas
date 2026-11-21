#!/usr/bin/env python3
"""Client for the Variant–Trait Claim Atlas (https://zxu69461-creator.github.io/variant-claim-atlas/).

Reads the static data files of the web tool and reproduces its views and its verdict. Standard library only (Python 3.8+).

  python vca_client.py variant rs7903146                   # summary of the claims on a variant
  python vca_client.py check rs7903146 "type 2 diabetes"   # verdict for one variant-trait claim (names or abbreviations: T2D)
  python vca_client.py batch my_list.txt -o results.csv    # rsIDs, or "rsID, trait" pairs, one per line
  python vca_client.py claims rs1800795 -o claims.csv      # every claim on a variant

Use --base to read a local copy of the repository instead of the web site (python vca_client.py --base ./variant-claim-atlas ...).
As a module: c = Atlas(); c.check(7903146, "T2D")."""
import argparse, csv, json, re, sys, urllib.request
from pathlib import Path

BASE = 'https://zxu69461-creator.github.io/variant-claim-atlas/'


class Atlas:
    def __init__(self, base=BASE):
        self.base = base if base.endswith('/') or not base.startswith('http') else base + '/'
        self._cache = {}
        self.meta = self._get('data/meta.json')

    # ---- data access
    def _get(self, path):
        if path not in self._cache:
            if self.base.startswith('http'):
                with urllib.request.urlopen(self.base + path, timeout=120) as r:
                    self._cache[path] = json.load(r)
            else:
                self._cache[path] = json.loads((Path(self.base) / path).read_text())
        return self._cache[path]

    def _table(self, path):
        t = self._get(path)
        return [dict(zip(t['cols'], r)) for r in t['rows']]

    def variants(self):
        if 'V' not in self._cache:
            self._cache['V'] = {v['rsid']: v for v in self._table('data/variants.json')}
        return self._cache['V']

    def claims(self, rs):
        rs = _rs(rs)
        return [c for c in self._table(f"data/claims/{rs % self.meta['shards']:02d}.json") if c['rsid'] == rs]

    def traits_for(self, query):
        """Trait names that contain the query, and traits for which abstracts used the query as a phrase (as on the web site)."""
        q = query.lower().strip()
        if 'T' not in self._cache:
            self._cache['T'] = [t['trait'] for t in self._table('data/traits.json')]
            ph = {}
            for r in self._table('data/lexicon.json'):
                for p in r['phrases']:
                    ph.setdefault(p, []).append(r['trait'])
            self._cache['P'] = ph
        # names match at the start of a word; short queries (abbreviations) only as a whole word
        rx = re.compile(r'(?<![a-z0-9])' + re.escape(q) + (r'(?![a-z0-9])' if len(q) <= 4 else ''))
        return {t for t in self._cache['T'] if rx.search(t.lower())} | set(self._cache['P'].get(q, []))

    # ---- views
    def variant(self, rs):
        return self.variants().get(_rs(rs))

    def check(self, rs, trait=None):
        """Verdict for the claims on variant rs that concern the trait (all claims when trait is None)."""
        cl = self.claims(rs)
        if trait:
            keep = self.traits_for(trait)
            cl = [c for c in cl if c['trait'] in keep]
        rep = [c for c in cl if c['status'] == 'reported']
        tst = [c for c in rep if c['testable']]
        sup = [c for c in tst if c['supported']]
        yrs = [c['first_gw_year'] for c in sup if c['first_gw_year']]
        if not cl:
            v = 'no claim on this trait'
        elif not rep:
            v = 'no reported association'
        elif not tst:
            v = 'not testable'
        elif len(sup) / len(tst) >= 0.5:
            v = 'supported'
        elif sup:
            v = 'partly supported'
        else:
            v = 'not supported'
        return dict(rsid=f'rs{_rs(rs)}', trait_query=trait or '', matched_traits='; '.join(sorted({c['trait'] for c in cl})), verdict=v,
                    reported=len(rep), absent=len(cl) - len(rep), testable=len(tst), supported=len(sup),
                    first_genome_wide_year=min(yrs) if yrs else '', earlier=sum(c['earlier'] for c in sup))

    def batch(self, lines):
        out = []
        for line in lines:
            line = line.strip()
            m = re.match(r'^(rs\s?\d+)\s*[,;:\t]\s*(\S.*)$', line, re.I)
            items = [(m.group(1), m.group(2).strip('"\' '))] if m else [(f'rs{x}', None) for x in re.findall(r'rs\s?(\d+)', line, re.I)]
            for rs, trait in items:
                x = self.variant(rs)
                if x is None:
                    out.append(dict(rsid=f'rs{_rs(rs)}', trait_query=trait or '', verdict='not in the atlas'))
                elif trait:
                    out.append(dict(gene=x['gene'], **self.check(rs, trait)))
                else:
                    out.append(dict(rsid=f'rs{x["rsid"]}', gene=x['gene'], verdict='', reported=x['reported'], testable=x['testable'],
                                    supported_share=x['supported'], matched_null=x['null'], clinvar=x['clinvar'],
                                    first_genome_wide_year=x['first_gw_year'] or ''))
        return out


def _rs(x):
    return int(re.sub(r'\D', '', str(x)))


def _write(rows, path):
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    f = open(path, 'w', newline='') if path else sys.stdout
    w = csv.DictWriter(f, cols); w.writeheader(); w.writerows(rows)
    if path:
        f.close(); print(f'{len(rows)} rows written to {path}', file=sys.stderr)


def main():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument('--base', default=BASE, help='web site or local folder holding data/ (default: the web site)')
    s = a.add_subparsers(dest='cmd', required=True)
    p = s.add_parser('variant'); p.add_argument('rsid')
    p = s.add_parser('check'); p.add_argument('rsid'); p.add_argument('trait', nargs='?')
    p = s.add_parser('claims'); p.add_argument('rsid'); p.add_argument('-o', '--out')
    p = s.add_parser('batch'); p.add_argument('file'); p.add_argument('-o', '--out')
    o = a.parse_args(); c = Atlas(o.base)
    if o.cmd == 'variant':
        print(json.dumps(c.variant(o.rsid), indent=1))
    elif o.cmd == 'check':
        print(json.dumps(c.check(o.rsid, o.trait), indent=1))
    elif o.cmd == 'claims':
        _write(c.claims(o.rsid), o.out)
    else:
        _write(c.batch(Path(o.file).read_text().splitlines()), o.out)


if __name__ == '__main__':
    main()

"""Genome position (GRCh38) and first listed gene of every claimed variant, from dbSNP (NCBI esummary), for the genome-wide overview of
Figure 1. Resumable; writes work/fig/dbsnp_pos.tsv (rsid, chr, bp, gene, snp_class). usage: python3 scripts/40_positions.py"""
import json, subprocess, sys, time
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent)); import figdata as F
OUT = F.ROOT / 'work' / 'fig' / 'dbsnp_pos.tsv'
ids = sorted(set(int(r) for r in F.objects().query("outcome_type != 'drug_response'").rsid))
done = pd.read_csv(OUT, sep='\t', dtype=str) if OUT.exists() else pd.DataFrame(columns=['rsid', 'chr', 'bp', 'gene', 'snp_class'])
have = set(done.rsid.astype('int64')) if len(done) else set()
todo = [i for i in ids if i not in have]
# variants without a position from the GWAS Catalog or the 1000 Genomes reference first, then the most-claimed ones
fb = F.variants().set_index('rsid')
todo.sort(key=lambda i: (fb.src.get(i) in ('gwas_catalog', '1000g'), -int(fb.papers.get(i, 0))))
print(f'{len(ids)} claimed variants, {len(have)} cached, {len(todo)} to fetch', flush=True)
rows, fails = [], 0
for k in range(0, len(todo), 200):
    chunk = todo[k:k + 200]
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=snp&retmode=json&id=' + ','.join(map(str, chunk))
    j = None
    for attempt in range(4):
        r = subprocess.run(['ncurl', '-s', '-m', '90', url], capture_output=True, text=True)
        try:
            j = json.loads(r.stdout)['result']; break
        except Exception:
            time.sleep(2 + 3 * attempt)
    if j is None:
        fails += 1; print('failed chunk', k, flush=True); continue
    for u in j.get('uids', []):
        d = j[u]
        cp = (d.get('chrpos') or '').split(':')
        rows.append((int(u), cp[0] if len(cp) == 2 else '', cp[1] if len(cp) == 2 else '',
                     (d.get('genes') or [{}])[0].get('name', '') if d.get('genes') else '', d.get('snp_class', '')))
    got = {r[0] for r in rows}
    rows += [(i, '', '', '', 'not_found') for i in chunk if i not in got and str(i) not in j.get('uids', [])]   # merged or withdrawn ids
    if (k // 200) % 10 == 9 or k + 200 >= len(todo):
        done = pd.concat([done, pd.DataFrame(rows, columns=done.columns).astype(str)]).drop_duplicates('rsid', keep='last')
        done.to_csv(OUT, sep='\t', index=False); rows = []
        print(f'{len(done)} / {len(ids)}', flush=True)
    time.sleep(0.34)
if len(done) == 0:
    raise SystemExit('dbSNP returned no records')
print('done; failed chunks:', fails)

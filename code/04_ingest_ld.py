"""Ingest LD proxies (plink 1.9 --r2, r2 >= 0.8, 1 Mb window, 1000G EUR and EAS MAGMA g1000 bfiles; M5 job scripts/m5/run_04_ld.sh) for
every rsID named in a candidate abstract. One row per (query rsID, proxy rsID, population); self-pairs kept (r2 = 1)."""
import sys, gzip
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import pandas as pd
from rdata.proj import Project
P = Project(__file__); W = P.root / 'work' / 'ld'
with P.run(__file__, params={'r2_min': 0.8, 'window_kb': 1000, 'reference': '1000G phase 3 (MAGMA g1000 EUR/EAS)'}) as run:
    out = []
    for pop in ('eur', 'eas'):
        f = W / f'{pop}_proxies.ld.gz'; run.input(f'ld_{pop}', path=f)
        d = pd.read_csv(f, sep=r'\s+', usecols=['CHR_A', 'BP_A', 'SNP_A', 'BP_B', 'SNP_B', 'R2'])
        d = d[d.SNP_A.str.match(r'^rs\d+$') & d.SNP_B.str.match(r'^rs\d+$')]
        out.append(pd.DataFrame({'pop': pop, 'chr': d.CHR_A.astype(int), 'query_rsid': d.SNP_A.str[2:].astype('int64'), 'query_bp': d.BP_A,
                                 'proxy_rsid': d.SNP_B.str[2:].astype('int64'), 'proxy_bp': d.BP_B, 'r2': d.R2}))
    ld = pd.concat(out, ignore_index=True)
    run.parquet('ld_proxies', ld, note='LD proxies r2>=0.8 within 1 Mb for candidate-abstract rsIDs, 1000G EUR and EAS, hg19 positions')
    print(ld.groupby('pop').agg(pairs=('r2', 'size'), queries=('query_rsid', 'nunique')))

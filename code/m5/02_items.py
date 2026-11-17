"""Reader items: rsID candidate abstracts + the diseases PubTator3 tags in each (M5 job).
Usage: 02_items.py CAND_DIR DISEASE2PUBTATOR3_GZ OUT_PARQUET
Disease ids are PubTator3 MeSH (or OMIM) concept ids; records with no tagged disease are kept (traits such as BMI or lipid levels are often
untagged) and flagged. Writes the candidate columns plus diseases (ids) and disease_names (first mention)."""
import sys, gzip, collections
import duckdb
cand_dir, d2p, out_p = sys.argv[1:4]
con = duckdb.connect()
cand = con.sql(f"select * from read_parquet('{cand_dir}/cand_*.parquet', union_by_name=true) where pmid is not null").df()
want = set(cand.pmid.astype(str))
dis = collections.defaultdict(dict)
with gzip.open(d2p, 'rt') as fh:
    for line in fh:
        p = line.split('\t', 1)[0]
        if p not in want: continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 4 or f[2] in ('-', ''): continue
        for cid in f[2].split(';'):
            dis[int(p)].setdefault(cid, f[3].split('|')[0])
cand['diseases'] = cand.pmid.map(lambda p: sorted(dis.get(p, {})))
cand['disease_names'] = cand.pmid.map(lambda p: [dis[p][k] for k in sorted(dis.get(p, {}))])
print(f'candidates {len(cand)}, with tagged disease {(cand.diseases.map(len) > 0).sum()}')
cand.to_parquet(out_p, index=False)

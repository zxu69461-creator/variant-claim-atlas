"""Reader check sample (charter amendment 4.1): 80 positive and 40 null held-out records with a linked trait (one per abstract), plus
30 further held-out abstracts with 1-6 listed rsIDs for recall. Reader-facing items carry only title, abstract, rsID and the reader's
outcome text; the key (reader status) stays in work/audit/key_v1.tsv."""
import sys, json, hashlib
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import pandas as pd, numpy as np, duckdb
from rdata.proj import Project
P = Project(__file__); A = P.root / 'work' / 'audit'; SEED = 20261108
with P.run(__file__, params={'seed': SEED, 'n_pos': 80, 'n_null': 40, 'n_recall': 30}) as run:
    con = duckdb.connect(str(P.root / 'db' / 'varclaims.duckdb'), read_only=True)
    rec = con.sql("""select pmid, rsid, trait_key, pop_group, outcome, status from variant_claims
                     where not dev and outcome_type <> 'drug_response' and trait_id is not null and status in ('positive','null')""").df()
    anypos = con.sql("select distinct pmid, rsid from variant_claims where not dev and status = 'positive'").df(); con.close()
    it = duckdb.sql(f"select pmid, title, abstract, rsids from '{P.root / 'work/reader/items_all.parquet'}'").df().set_index('pmid')
    rng = np.random.default_rng(SEED)
    rec = rec.sample(frac=1, random_state=SEED).drop_duplicates('pmid')  # one record per abstract, uniform within status
    pos = rec[rec.status == 'positive'].head(80); nul = rec[~rec.pmid.isin(pos.pmid) & (rec.status == 'null')].head(40)
    pairs = pd.concat([pos, nul]).sample(frac=1, random_state=SEED + 1)
    used = set(pairs.pmid)
    ho = [p for p in it.index if hashlib.md5(str(p).encode()).hexdigest()[-1] not in '01' and p not in used and 1 <= len(it.rsids[p]) <= 6]
    rpm = list(rng.choice(sorted(ho), 30, replace=False))
    items, key = [], []
    for k, r in enumerate(pairs.itertuples(), 1):
        iid = f'P{k:03d}'; m = it.loc[r.pmid]
        items.append(dict(id=iid, kind='pair', title=m.title, abstract=m.abstract, rsid=f'rs{r.rsid}', outcome=r.outcome))
        key.append(dict(id=iid, pmid=r.pmid, rsid=r.rsid, trait_key=r.trait_key, pop_group=r.pop_group, reader_status=r.status))
    ap = set(zip(anypos.pmid, anypos.rsid))
    for k, p in enumerate(rpm, 1):
        iid = f'R{k:03d}'; m = it.loc[p]; rs = sorted(set(int(x) for x in m.rsids))
        items.append(dict(id=iid, kind='recall', title=m.title, abstract=m.abstract, rsids=[f'rs{x}' for x in rs]))
        for x in rs: key.append(dict(id=iid, pmid=p, rsid=x, reader_status='positive' if (p, x) in ap else 'none'))
    json.dump(items, open(A / 'items_v1.json', 'w'), ensure_ascii=False)
    pd.DataFrame(key).to_csv(A / 'key_v1.tsv', sep='\t', index=False)
    run.table('reader_audit_sample', pd.DataFrame(key), note='reader check sample and reader key (amendment 4.1); reader-facing items in work/audit/items_v1.json')
    print(len(items), 'items;', pairs.status.value_counts().to_dict(), '; recall rsIDs', sum(len(i.get('rsids', [])) for i in items))

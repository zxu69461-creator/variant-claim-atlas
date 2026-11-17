"""Claim records (charter amendment 3) from the local reader pass: one row per (PMID, rsID, linked trait, population group), with the
outcome-linking route, the EFO terms of the linked trait (amendment 2 map), population group and paper metadata. Rebuilt as the reader
pass grows; the development split is the `dev` flag."""
import sys, re, json, hashlib, collections
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import pandas as pd, duckdb
from rdata.proj import Project
P = Project(__file__); W = P.root / 'work' / 'reader'

EAS = re.compile(r'(?i)chinese|\bhan\b|japanese|korean|taiwanese|vietnamese|\bthai\b|malay|mongolian|east asian')
EUR_NAT = ('european|caucasian|white|british|english|scottish|welsh|irish|french|german|italian|spanish|portuguese|dutch|belgian|swiss|'
           'austrian|danish|swedish|norwegian|finnish|icelandic|polish|czech|slovak|hungarian|romanian|bulgarian|greek|croatian|serbian|'
           'slovenian|estonian|latvian|lithuanian|russian|ukrainian|belarusian|maltese|cypriot|luxembourg')
EUR = re.compile(rf'(?i)\b({EUR_NAT})')
EAS_C = {'China', 'Japan', 'South Korea', 'Taiwan', 'Hong Kong'}
EUR_C = {'United Kingdom', 'Germany', 'Spain', 'Italy', 'Poland', 'France', 'Netherlands', 'Belgium', 'Switzerland', 'Austria', 'Denmark',
         'Sweden', 'Norway', 'Finland', 'Iceland', 'Ireland', 'Portugal', 'Greece', 'Czech Republic', 'Czechia', 'Slovakia', 'Hungary',
         'Romania', 'Bulgaria', 'Croatia', 'Serbia', 'Slovenia', 'Estonia', 'Latvia', 'Lithuania', 'Russia', 'Ukraine', 'Belarus', 'Malta',
         'Cyprus', 'Luxembourg', 'Bosnia and Herzegovina', 'North Macedonia', 'Montenegro', 'Albania'}

def pop_group(text, country):
    t = (text or '').strip()
    if t and t.lower() not in ('not_stated', 'not stated', 'null', 'none'):
        e, u = bool(EAS.search(t)), bool(EUR.search(t))
        return ('EAS' if e and not u else 'EUR' if u and not e else 'other'), False
    return ('EAS' if country in EAS_C else 'EUR' if country in EUR_C else 'not_stated'), True

TOK = re.compile(r'[a-z0-9]{3,}')
def link(outcome, ids, names, efo_lab):
    o = (outcome or '').lower().strip()
    if not o: return None, None
    cand = [(n.lower(), i) for i, n in zip(ids, names) if n]
    hit = [(len(n), i) for n, i in cand if n == o or n in o or o in n]
    if hit: return max(hit)[1], 'a_name'
    ot = set(TOK.findall(o))
    best = max(((len(ot & set(TOK.findall(n))) / len(ot | set(TOK.findall(n))), i) for n, i in cand if ot), default=(0, None))
    if best[0] >= 0.5: return best[1], 'b_tokens'
    if o in efo_lab: return 'EFO=' + efo_lab[o], 'c_efo_label'
    if len(set(ids)) == 1: return ids[0], 'd_single'
    return None, None

def efo_labels(obo):
    lab, cur = {}, None
    norm = lambda x: x[4:] if x.startswith('efo:') else x.replace(':', '_')
    for line in open(obo):
        if line.startswith('['): cur = None
        elif line.startswith('id: '): cur = norm(line[4:].strip())
        elif cur and (line.startswith('name: ') or (line.startswith('synonym: "') and ' EXACT' in line)):
            s = line[6:].strip() if line.startswith('name: ') else line.split('"')[1]
            lab.setdefault(s.lower(), cur)
    return lab

with P.run(__file__, params={'reader': 'Qwen3.5-9B-MLX-4bit, prompt v2', 'rules': 'charter amendment 3'}) as run:
    run.input('reader_all_v2', path=W / 'reader_all_v2.jsonl'); run.input('items_all', path=W / 'items_all.parquet')
    run.input('EFO_202611', path=P.root / 'raw' / 'efo.obo')
    con = duckdb.connect(str(P.root / 'db' / 'varclaims.duckdb'), read_only=True)
    tmap = con.sql('select pubtator_id, efo_id, route from claim_trait_map where efo_id is not null').df()
    terms = set(con.sql('select efo_id from efo_terms').df().efo_id); con.close()
    id2efo = tmap.groupby('pubtator_id').efo_id.agg(lambda s: sorted(set(s))).to_dict()
    id2route = tmap.groupby('pubtator_id').route.first().to_dict()
    lab = {k: v for k, v in efo_labels(P.root / 'raw' / 'efo.obo').items() if v in terms}
    it = duckdb.sql(f"select pmid, pub_year, journal, issn, country, ptypes, retracted, rsids, diseases, disease_names from '{W / 'items_all.parquet'}'"
                    ).df().set_index('pmid')
    rows, stat = [], collections.Counter()
    for l in open(W / 'reader_all_v2.jsonl'):
        r = json.loads(l); stat['abstracts'] += 1
        if r['claims'] is None: stat['parse_fail'] += 1; continue
        m = it.loc[r['pmid']]; listed = set(int(x) for x in m.rsids)
        ids, names = list(m.diseases), list(m.disease_names)
        for c in r['claims']:
            if not isinstance(c, dict): stat['not_object'] += 1; continue
            rs = re.fullmatch(r'(?i)rs ?(\d+)', str(c.get('rsid') or '').strip())
            if not rs or int(rs.group(1)) not in listed: stat['rsid_not_listed'] += 1; continue
            ot = c.get('outcome_type')
            if ot not in ('disease', 'trait', 'severity', 'drug_response'): stat['outcome_type_other'] += 1; continue
            tid, lroute = link(c.get('outcome'), ids, names, lab)
            if tid is None: efo, mroute = [], None
            elif tid.startswith('EFO='): efo, mroute = [tid[4:]], 'efo_label'
            else: efo, mroute = id2efo.get(tid, []), id2route.get(tid)
            pg, inferred = pop_group(c.get('population'), m.country)
            rows.append(dict(pmid=r['pmid'], rsid=int(rs.group(1)), outcome=c.get('outcome'), outcome_type=ot, trait_id=tid,
                             link_route=lroute, map_route=mroute, efo=efo, pop_group=pg, pop_inferred=inferred, population=c.get('population'),
                             associated=c.get('associated'), direction=c.get('direction'), risk_allele=c.get('risk_allele'),
                             effect=c.get('effect'), p_value=c.get('p_value'), n_cases=c.get('n_cases'), n_controls=c.get('n_controls'),
                             n_total=c.get('n_total'), design=c.get('design'), quote=c.get('quote'), truncated=r['truncated']))
            stat['entries'] += 1
    e = pd.DataFrame(rows)
    for col in ('n_cases', 'n_controls', 'n_total'): e[col] = pd.to_numeric(e[col], errors='coerce')
    e['trait_key'] = e.trait_id.fillna('outcome:' + e.outcome.fillna('').str.lower().str.strip())
    g = e.groupby(['pmid', 'rsid', 'trait_key', 'pop_group'], sort=False)
    rec = g.agg(outcome=('outcome', 'first'), outcome_type=('outcome_type', 'first'), trait_id=('trait_id', 'first'),
                link_route=('link_route', 'first'), map_route=('map_route', 'first'), efo=('efo', 'first'), pop_inferred=('pop_inferred', 'all'),
                population=('population', 'first'), n_entries=('associated', 'size'),
                any_yes=('associated', lambda s: (s == 'yes').any()), all_no=('associated', lambda s: (s == 'no').all()),
                direction=('direction', 'first'), risk_allele=('risk_allele', 'first'), effect=('effect', 'first'), p_value=('p_value', 'first'),
                n_cases=('n_cases', 'max'), n_controls=('n_controls', 'max'), n_total=('n_total', 'max'), design=('design', 'first'),
                quote=('quote', 'first')).reset_index()
    rec['status'] = rec.any_yes.map({True: 'positive'}).fillna(rec.all_no.map({True: 'null', False: 'unclear'}))
    rec = rec.drop(columns=['any_yes', 'all_no'])
    rec['n'] = rec.n_total.fillna(rec.n_cases.fillna(0) + rec.n_controls.fillna(0)).replace(0, pd.NA)
    meta = it[['pub_year', 'journal', 'issn', 'country', 'ptypes', 'retracted']]
    rec = rec.join(meta, on='pmid')
    rec['dev'] = rec.pmid.map(lambda p: hashlib.md5(str(p).encode()).hexdigest()[-1] in '01')
    run.parquet('variant_claims', rec, note='variant-trait claim records (charter amendment 3) from the local reader pass; dev flag = development split')
    print(dict(stat)); print(len(rec), 'records;', rec.dev.sum(), 'development')
    d = rec[rec.dev & (rec.outcome_type != 'drug_response')]
    print(d.status.value_counts().to_dict()); print(d.link_route.value_counts(dropna=False).to_dict())
    print('mapped to EFO:', f'{(d.efo.str.len() > 0).mean():.3f}'); print(d.pop_group.value_counts().to_dict())

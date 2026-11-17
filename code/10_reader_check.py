"""Reader check (charter amendment 4.1): score the first author's answers against the reader key. Input: work/audit/responses_v1.json
(the questionnaire's `answers` documents, {item id: {v, note, at}}). Writes reader_check_result and work/audit/gate.json; exit 0 when
the gate passes (positive precision >= 0.85 and null precision >= 0.75), exit 3 when it fails, exit 2 when answers are incomplete."""
import sys, json, math
sys.path.insert(0, '/Volumes/ResearchData/Claude/rdata')
import pandas as pd
from rdata.proj import Project
P = Project(__file__); A = P.root / 'work' / 'audit'

def wilson(k, n, z=1.959964):
    if n == 0: return (float('nan'),) * 2
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d

with P.run(__file__, params={'gate_pos_precision': 0.85, 'gate_null_precision': 0.75}) as run:
    run.input('responses_v1', path=A / 'responses_v1.json'); run.input('key_v1', path=A / 'key_v1.tsv')
    ans = json.load(open(A / 'responses_v1.json')); key = pd.read_csv(A / 'key_v1.tsv', sep='\t', keep_default_na=False, dtype={'reader_status': str})
    pairs = key[key.id.str.startswith('P')].copy(); rec = key[key.id.str.startswith('R')].copy()
    pairs['human'] = pairs.id.map(lambda i: (ans.get(i) or {}).get('v'))
    rec['human'] = [((ans.get(i) or {}).get('v') or {}).get(f'rs{r}') for i, r in zip(rec.id, rec.rsid)]
    missing = int(pairs.human.isna().sum() + rec.human.isna().sum())
    rows = []
    def add(measure, k, n, note=''):
        lo, hi = wilson(k, n); rows.append(dict(measure=measure, k=int(k), n=int(n), estimate=k / n if n else float('nan'), ci_low=lo, ci_high=hi, note=note))
    pp = pairs[pairs.reader_status == 'positive']; pn = pairs[pairs.reader_status == 'null']
    add('positive_precision', (pp.human == 'yes').sum(), pp.human.notna().sum(), 'reader positive judged associated')
    add('null_precision', (pn.human == 'no').sum(), pn.human.notna().sum(), 'reader null judged not associated')
    add('pair_error', (pairs.human == 'notstudied').sum(), pairs.human.notna().sum(), 'pair not studied in the abstract')
    add('no_conclusion', (pairs.human == 'none').sum(), pairs.human.notna().sum(), 'pair studied, no conclusion')
    hy = rec[rec.human == 'yes']
    add('recall', (hy.reader_status == 'positive').sum(), len(hy), 'human-yes rsIDs with a reader-positive record')
    hn = rec[rec.human == 'no']
    add('recall_specificity', (hn.reader_status == 'none').sum(), len(hn), 'human-no rsIDs without a reader-positive record')
    res = pd.DataFrame(rows)
    pos, nul = res.set_index('measure').estimate[['positive_precision', 'null_precision']]
    gate = 'incomplete' if missing else ('pass' if pos >= 0.85 and nul >= 0.75 else 'fail')
    run.table('reader_check_result', res, note=f'reader check (amendment 4.1); gate {gate}; missing answers {missing}')
    run.table('reader_check_answers', pd.concat([pairs.assign(part='pair'), rec.assign(part='recall')]), note='per-item reader key and human answer')
    json.dump(dict(gate=gate, positive_precision=pos, null_precision=nul, missing=missing), open(A / 'gate.json', 'w'))
    print(res.to_string(index=False, float_format=lambda x: f'{x:.3f}')); print('GATE', gate)
sys.exit({'pass': 0, 'fail': 3, 'incomplete': 2}[gate])

"""Local reader model (MLX on M5) for variant-disease association claims. Usage: 05_llm_read.py MODEL ITEMS_PARQUET OUT_JSONL [BS]
ITEMS: pmid, title, abstract, rsids (ints), disease_names (PubTator3 first mentions). One prompt per abstract; the model may only use the
listed rsIDs. Output: one JSON line per pmid with the parsed claim list (or the raw text when parsing fails). Resumes.
Adapted from depclaims/scripts/m5/05_llm_read.py (same model, decoding and batching)."""
import sys, os, json, re, time
import pandas as pd
from mlx_lm import load, generate, batch_generate
from mlx_lm.sample_utils import make_sampler
model_id, items_p, out_p = sys.argv[1:4]; BS = int(sys.argv[4]) if len(sys.argv) > 4 else 8
SYS = ("You extract genetic association results from biomedical abstracts. Report only what the authors themselves tested and found, as "
       "stated in the title or abstract. Answer with JSON only.")
TPL = """Title: {title}
Abstract: {abstract}

Candidate variants: {rsids}
Diseases or traits tagged in the text: {diseases}

List every result in which the authors tested one candidate variant for association with a disease, disease risk, or a human trait
(including quantitative traits such as BMI or lipid levels, and disease severity or progression), in human participants. Include results
reported as not associated. Do not count variants only genotyped for description, associations with drug response or adverse drug
reactions (record those with outcome_type "drug_response"), functional or cell-line experiments, or results the authors only cite from
earlier work. If a result is given separately for several populations or subgroups, give one entry per population.
For each result copy, word for word from the title or abstract, the sentence that states it ("quote").
Return JSON: {{"claims": [{{"rsid": "rs...", "outcome": "disease or trait as written", "outcome_type": "disease|trait|severity|drug_response",
"associated": "yes|no|unclear", "risk_allele": "A|C|G|T|null", "direction": "increased|decreased|null", "effect": "e.g. OR 1.45 (1.12-1.88) or null",
"p_value": "as written or null", "population": "e.g. Han Chinese, European, mixed, not_stated", "n_cases": int or null, "n_controls": int or null,
"n_total": int or null, "design": "case-control|cohort|GWAS|family|meta-analysis|other", "quote": "..."}}]}}.
Write the JSON compactly on one line per claim, without indentation. Return {{"claims": []}} if there is none."""
done = set()
if os.path.exists(out_p):
    for l in open(out_p):
        try: done.add(json.loads(l)['pmid'])
        except Exception: pass
it = pd.read_parquet(items_p); it = it[~it.pmid.isin(done)]
model, tok = load(model_id); samp = make_sampler(temp=0.0)
def prompt_of(r):
    msg = [{'role': 'system', 'content': SYS}, {'role': 'user', 'content': TPL.format(title=r.title, abstract=r.abstract or '',
           rsids=', '.join(f'rs{x}' for x in r.rsids), diseases=', '.join(r.disease_names) or 'none tagged')}]
    try: return tok.apply_chat_template(msg, add_generation_prompt=True, tokenize=False, enable_thinking=False)
    except TypeError: return tok.apply_chat_template(msg, add_generation_prompt=True, tokenize=False)
def parse(txt):
    txt = re.sub(r'(?s)<think>.*?</think>', '', txt).strip(); m = re.search(r'(?s)\{.*\}', txt)
    try: return (json.loads(m.group(0))['claims'] if m else None), txt
    except Exception: pass
    part = []  # truncated output: keep every complete claim object before the cut
    for o in re.findall(r'\{[^{}]*\}', txt):
        try: part.append(json.loads(o))
        except Exception: pass
    return (part or None), txt
t0 = time.time(); rows = list(it.itertuples())
with open(out_p, 'a') as out:
    for k in range(0, len(rows), BS):
        chunk = rows[k:k + BS]
        toks = [tok.encode(prompt_of(r), add_special_tokens=False) for r in chunk]
        res = batch_generate(model, tok, toks, max_tokens=2500, sampler=samp).texts if BS > 1 else [generate(model, tok, prompt=prompt_of(chunk[0]), max_tokens=2500, sampler=samp)]
        for r, txt in zip(chunk, res):
            claims, txt = parse(txt)
            out.write(json.dumps({'pmid': int(r.pmid), 'claims': claims, 'raw': None if claims is not None else txt[:4000], 'truncated': claims is not None and not txt.rstrip().endswith('}')}) + '\n')
        out.flush(); print(k + len(chunk), f'{(time.time() - t0) / (k + len(chunk)):.2f}s/item', flush=True)
print('done', len(it), f'{time.time() - t0:.0f}s')

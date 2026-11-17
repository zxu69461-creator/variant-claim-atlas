"""AJHG supplemental package in manuscript/ajhg/: Data S1 (claim atlas, all records, with the genome-wide comparison for testable claims in
both splits), Data S2 (analysis code, a README and the dated analysis plan), and the Markdown source of the supplemental PDF
(Supplemental Methods with the analysis plan and the extraction instructions, Figure S1 from figures/ajhg/figS1, Tables S1-S2)."""
import importlib.util, re, sys, zipfile
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent)); import figdata as F
O = F.ROOT / 'manuscript' / 'ajhg'
with F.con() as c:
    m = pd.concat([c.sql(f"select *, '{s}' as split_matched from {t}").df() for s, t in (('development', 'variant_claim_match'), ('held-out', 'variant_claim_match_heldout'))])
    v = c.sql("""select pmid, rsid, trait_key, pop_group, outcome, outcome_type, trait_id, link_route, map_route, array_to_string(efo, '|') as efo_terms,
                 pop_inferred, population, status, direction, risk_allele, effect, p_value, n_cases, n_controls, n_total, n as stated_n, design,
                 pub_year, journal, country, retracted, case when dev then 'development' else 'held-out' end as split from variant_claims""").df()
    rc = c.sql("select measure, k, n, estimate, ci_low, ci_high from reader_check_result").df()
    pr = c.sql("select test, estimate, ci_low, ci_high, p_holm, n, n_papers from confirm_primary").df().assign(set='primary')
    se = c.sql("select set, test, estimate, ci_low, ci_high, p, n, n_papers from confirm_sensitivity").df()
    ad = c.sql("select set, test, estimate, ci_low, ci_high, p, n, n_papers from additional_sensitivity").df()
    ndev, nho = c.sql("select count(distinct pmid) filter (where dev), count(distinct pmid) filter (where not dev) from variant_claims").fetchone()
    dv = c.sql("select test, extra from confirm_primary_dev").df().set_index('test').extra
cols = ['pmid', 'rsid', 'trait_key', 'pop_group', 'testable', 'replicated', 'rep_prior', 'first_hit_year', 'null_share', 'maf_eur', 'n_proxies']
a = v.merge(m[cols], on=['pmid', 'rsid', 'trait_key', 'pop_group'], how='left')
a = a.rename(columns={'replicated': 'genome_wide_support', 'rep_prior': 'support_before_claim', 'first_hit_year': 'first_matching_gwas_year',
                      'null_share': 'matched_null_share', 'maf_eur': 'maf_1000g_eur', 'n_proxies': 'ld_proxies'})
a.insert(1, 'rsid_label', 'rs' + a.rsid.astype(str))
a.to_csv(O / 'Data_S1_variant_trait_claim_atlas.tsv.gz', sep='\t', index=False, compression='gzip')
# Data S2: the code that produced the atlas, the estimates and the figures, in order of use; superseded figure scripts, the project
# template and the earlier manuscript builder are left out
SKIP = {'00_template.py', '21_fig1_overview.py', '22_fig2_support.py', '23_fig3_size.py', '24_fig4_time.py', '25_figS1_sensitivity.py', '31_build_docx.py'}
plan = (F.ROOT / 'CHARTER.md').read_text()
HOUSE = r' The \w+ directory is never touched\.'               # one sentence about local file storage, not about the analysis
assert len(re.findall(HOUSE, plan)) == 1 and plan.count('Reader: the first author (one reader).') == 1
plan = re.sub(HOUSE, '', plan).replace('Reader: the first author (one reader).', 'Reader: the first author (one reader) [Z.X.; the order of authors was set later].')
plan = ('> Analysis plan of the study as written, with the dates of its amendments (4, 7 and 8 November 2026). The text is unchanged except '
        'that one sentence about local file storage was removed and one note in square brackets was added. "Reader" is the language model '
        'in amendments 1-3 and section 4.1 names the human reader; "arbiter" is the source of genome-wide evidence. The comparison with '
        'FinnGen considered in amendment 1 (arbiter B) was not part of the confirmatory plan (amendment 4) and was not carried out.\n\n') + plan
README = """Analysis code for "A literature-wide atlas of genetic association claims calibrated against genome-wide evidence"

ANALYSIS_PLAN.md   the analysis plan with its dated amendments; amendment 4 (8 November 2026) is the confirmatory plan
m5/                candidate abstracts from the PubMed baseline (01, 02, pmparse), linkage-disequilibrium proxies and the matched-null
                   sample (run_04_ld, run_06_null), and the language-model pass over the abstracts (05_llm_read; run_05)
03-06              GWAS Catalog, linkage disequilibrium, ontology tables and the mapping of disease identifiers to ontology terms
07_claims          claim records from the model output
08_match           each claim against the GWAS Catalog: testability, support at any time and earlier, matched null
09, 10, 12, 12b    the blinded reading check: sample, reading page, collection of answers and scoring
11_confirm         the five primary estimates, secondary breakdowns and the five planned sensitivity analyses
13_additional      additional analyses added after the plan
30, 32             numbers quoted in the text; supplemental tables and Data S1
40-45, figdata     variant positions and the figures (Figures 1-4 and S1)
ajhg_genes         gene symbols and MIM numbers in the text
run_chain.sh       the order in which the confirmatory steps were run

The scripts read a local project database built by steps 03-08 from the public inputs listed in the article; paths at the top of
each script point to that database and to a shared helper for recording inputs and outputs, and have to be adapted before re-use.
"""
with zipfile.ZipFile(O / 'Data_S2_analysis_code.zip', 'w', zipfile.ZIP_DEFLATED) as z:
    for f in sorted((F.ROOT / 'scripts').rglob('*')):
        if f.is_file() and f.suffix in ('.py', '.sh', '.html', '.js', '.tsv') and '__pycache__' not in f.parts and f.name not in SKIP:
            z.write(f, 'code/' + str(f.relative_to(F.ROOT / 'scripts')))
    z.writestr('code/ANALYSIS_PLAN.md', plan)
    z.writestr('code/README.txt', README)
names = {'positive_precision': 'Reported associations judged associated (precision)', 'null_precision': 'Reported absences judged not associated (precision)',
         'pair_error': 'Pair not studied in the abstract', 'no_conclusion': 'Pair studied without a conclusion', 'recall': 'Recall of reported associations',
         'recall_specificity': 'Variants without a reported association correctly absent'}
t1 = ['| Measure | Count | Total | Estimate (95% CI) |', '|---|---|---|---|'] + [
    f"| {names[r.measure]} | {r.k} | {r.n} | {r.estimate * 100:.1f}% ({r.ci_low * 100:.1f}–{r.ci_high * 100:.1f}) |" for r in rc.itertuples()]
setn = {'primary': 'All claims', 'a_no_link_d': 'Outcome linked by name or ontology', 'b_no_mesh_ancestor': 'Trait mapped without broadening',
        'c_population_stated': 'Population stated in the abstract', 'd_no_retracted': 'Retracted papers removed', 'e_paper_weighted': 'One weight per paper',
        'x_no_gwas_meta': 'Genome-wide studies and meta-analyses excluded^a^'}
testn = {'C1': 'Support above matched null (percentage points)', 'C2': 'OR per 10-fold sample size', 'C3': 'Restatement, OR per 5 years',
         'C4': 'Five-year yield of new claims, OR per 5 years', 'C5': 'Discordant reports at established loci, OR per 5 years'}
allr = pd.concat([pr, se, ad])
t2 = ['| Estimate | Analysis set | Estimate (95% CI) | Claims | Papers |', '|---|---|---|---|---|']
for t in testn:
    for s_ in setn:
        r = allr[(allr.test == t) & (allr.set == s_)].iloc[0]; k = 100 if t == 'C1' else 1
        t2.append(f"| {testn[t]} | {setn[s_]} | {r.estimate * k:.2f} ({r.ci_low * k:.2f}–{r.ci_high * k:.2f}) | {int(r.n):,} | {int(r.n_papers):,} |")
# the extraction instructions, word for word from the script that ran the model
rd = (F.ROOT / 'scripts' / 'm5' / '05_llm_read.py').read_text()
SYS = ''.join(re.findall(r'"([^"]*)"', re.search(r'SYS = \((.*?)\)\n', rd, re.S).group(1)))
TPL = re.search(r'TPL = """(.*?)"""', rd, re.S).group(1).replace('{{', '{').replace('}}', '}')
assert SYS.startswith('You extract genetic association results') and TPL.startswith('Title: {title}') and TPL.rstrip().endswith('if there is none.')
HARD = ('Title: {title}', 'Abstract: {abstract}', 'Candidate variants: {rsids}')       # lines of the instruction that end with a line break
quote = lambda t: '\n'.join(('> ' + l + ('\\' if l in HARD else '')) if l.strip() else '>' for l in t.strip().split('\n'))
d1 = [float(x.split()[1]) for x in dv['C1'].split(';')]         # development papers: supported share and matched null
md = f"""# Supplemental information

## Supplemental methods

### Analysis plan

The analysis followed a written plan with dated amendments, which is reproduced in full in Data S2 and summarized in this section.

*Rules fixed before any claim was compared with genome-wide evidence (4 and 7 November 2026).* The unit of analysis is the claim record
(paper, variant, linked trait and population group). The rules for linking the outcome text to a disease identifier and to ontology
terms, for grouping populations, for matching a claim trait to a GWAS Catalog trait, for testability, for support at any time and
earlier than the claim, and for the matched null are those given in Material and methods. They were written before the claims of any
paper, development papers included, were compared with the GWAS Catalog, and the ontology mapping uses only ontology files, so it could
not be adjusted to the comparison.

*Development and held-out papers.* Papers were assigned by a fixed pseudorandom rule that depends on the PMID alone (given in Data S2),
which places one-eighth of the papers ({ndev:,}) in the development set, the only papers used for exploration, and the other {nho:,}
in the held-out set.

*Confirmatory plan (8 November 2026).* The plan was written after the language model had read all abstracts and after exploration in the
development papers, and before any held-out claim was compared with the GWAS Catalog; it specified the six points below.

1. Reading check. From held-out records with a linked trait, 80 reported associations and 40 reported absences were drawn at random (at
   most one per abstract; seed 20261108), with 30 further held-out abstracts naming one to six variants for recall. The reader saw the
   title, the abstract, the variant and the outcome text of the record, and never the model's classification, direction or population
   or any genome-wide result. The held-out claims were to be compared with genome-wide evidence only if the precision of reported
   associations was at least 85% and that of reported absences at least 75%; otherwise the instructions were to be revised and all
   abstracts read again. Recall was to be reported without a threshold.
2. Analysis set. Held-out claim records on a disease, a quantitative trait or disease severity that report an association or its
   absence and are testable; papers are the clustering unit.
3. Five primary estimates, two-sided and Holm-adjusted at α = 0.05, with the direction expected from the development papers: (i) the
   share of reported associations with genome-wide support minus the mean matched null, by bootstrap of papers (development papers:
   {d1[0] * 100:.0f}% against {d1[1] * 100:.0f}%); (ii) the odds ratio of support per ten-fold larger stated sample size (expected to
   rise); (iii) among supported associations published from 2005, the odds ratio per five years that the support was earlier than the
   claim (expected to rise); (iv) among new reported associations published in 2005–2020, the odds ratio per five years of support first
   published within five years of the claim (expected to fall); (v) among reported absences published from 2005, the odds ratio per five
   years that the variant was already genome-wide significant for the trait (expected to rise). Estimates (ii)–(v) use logistic
   generalized estimating equations with an exchangeable working correlation.
4. Secondary breakdowns, reported with intervals and without tests: support by design, population group, stated sample size (fewer than
   500, 500–1,999, 2,000–9,999 and 10,000 or more participants) and era.
5. Five sensitivity analyses, each repeating the five primary estimates: outcomes linked only through the single disease of an abstract
   excluded; traits mapped through a broader MeSH ancestor excluded; populations stated in the abstract only; retracted papers excluded;
   one weight per paper.
6. Effect estimates are not pooled across studies, and development and held-out results are reported separately.

Analyses added after the plan are listed as additional analyses in Material and methods; they include the trait areas, the breakdowns
within trait areas, designs, populations and eras, the adjustment for stated sample size and the exclusion of claims from genome-wide
studies and meta-analyses. The reading check was done by one author (Z.X.).

### Extraction instructions

Every abstract was presented to the language model with the system message and the instruction below; the fields in braces were filled
with the title, the abstract, the dbSNP identifiers found in them and the disease names tagged by PubTator3 ("none tagged" when there
were none). The text is reproduced as used.

System message:

{quote(SYS)}

Instruction:

{quote(TPL)}

The output was parsed as JSON; when the output was cut off at the token limit, every complete entry before the cut was kept. An entry
was used only when its variant is one of the candidate variants listed for the abstract.

**Figure S1. The primary estimates under the sensitivity analyses**

![](../../figures/ajhg/figS1_sensitivity.png)

Each column shows one primary estimate for all held-out claims (red), for the five prespecified sensitivity analyses (filled circles:
outcomes linked only through the single disease of an abstract excluded, traits mapped through a broader MeSH ancestor excluded,
populations stated in the abstract only, retracted papers excluded, and one weight per paper) and, below the dotted line, for an
additional analysis added after the plan, in which claims from genome-wide studies and meta-analyses were excluded (open circles).
Points are estimates and lines 95% confidence intervals, with the estimate printed beside each; the band is the 95% confidence interval
of the primary estimate and the dashed line the value of no effect.

**Table S1. Blinded reading check of the extracted claims in held-out papers**

{chr(10).join(t1)}

Eighty reported associations and 40 reported absences were drawn at random (at most one per abstract) and read by one author (Z.X.)
without the model's answer; recall was assessed in 30 further abstracts. Wilson 95% confidence intervals.

**Table S2. Primary estimates in held-out papers under each sensitivity analysis**

{chr(10).join(t2)}

Odds ratios from logistic generalized estimating equations with papers as clusters; the first estimate from 2,000 bootstrap resamples of
papers. ^a^Additional analysis added after the plan; the other five sensitivity analyses were prespecified.

**Data S1. Atlas of variant–trait claim records** (tab-separated, compressed; one row per paper, variant, linked trait and population
group, with the genome-wide comparison for testable claims).

**Data S2. Analysis code**, with the dated analysis plan.
"""
(O / 'supplemental_information.md').write_text(md)
print(len(a), 'atlas rows;', round((O / 'Data_S1_variant_trait_claim_atlas.tsv.gz').stat().st_size / 1e6, 1), 'MB; code zip',
      round((O / 'Data_S2_analysis_code.zip').stat().st_size / 1e3), 'kB')

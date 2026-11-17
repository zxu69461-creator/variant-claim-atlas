"""Data files for the v3 manuscript (2026-11-17): Data S1 (unchanged atlas, copied from the AJHG package build), Data S2 (analysis code
with the curated analysis plan, now including the ClinVar and robustness analyses and the v3 figure scripts) and Data S3 (ClinVar risk
variants claimed in held-out papers, per variant). Output manuscript/v3/data/."""
import shutil, sys, zipfile
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent)); import figdata as F
O = F.ROOT / 'manuscript' / 'v3' / 'data'; O.mkdir(parents=True, exist_ok=True)
OLD = F.ROOT / 'manuscript' / 'ajhg'
shutil.copy2(OLD / 'Data_S1_variant_trait_claim_atlas.tsv.gz', O / 'Data_S1_variant_trait_claim_atlas.tsv.gz')

# ---- Data S3
with F.con() as c:
    V = c.sql('select * from additional_clinvar_variants').df()
    C = c.sql('select rsid, cls, gwas_any as gwas_any_proxy, cond_term, cond_support from additional_clinvar_condition').df()
V = V.merge(C, on='rsid', how='left')
cls = {'mostly': 'mostly (>=50%)', 'partly': 'partly (5-50%)', 'rarely': 'rarely (<5%)'}
S3 = pd.DataFrame({
    'rsid': 'rs' + V.rsid.astype(str), 'gene_clinvar': V.gene, 'clinvar_significance': V.significance, 'clinvar_conditions': V.phenotypes,
    'reported_associations': V.claims, 'papers': V.papers, 'distinct_traits': V.traits,
    'share_supported': V.supported.round(4), 'matched_null': V['null'].round(4),
    'support_class': V.cls.map(cls).fillna('fewer than 5 reported associations'),
    'genome_wide_association_any_trait_variant': V.gwas_any,
    'genome_wide_association_any_trait_variant_or_proxy': V.gwas_any_proxy.map({True: 'yes', False: 'no'}).fillna(''),
    'clinvar_condition_in_ontology': V.cond_term.map({True: 'yes', False: 'no'}).fillna(''),
    'genome_wide_support_for_clinvar_condition': V.cond_support.map({True: 'yes', False: 'no'}).fillna('')})
S3['genome_wide_association_any_trait_variant'] = S3.genome_wide_association_any_trait_variant.map({True: 'yes', False: 'no'})
S3 = S3.sort_values('reported_associations', ascending=False)
S3.to_csv(O / 'Data_S3_clinvar_risk_variants.tsv', sep='\t', index=False)

# ---- Data S2
with zipfile.ZipFile(OLD / 'Data_S2_analysis_code.zip') as z:
    plan = z.read('code/ANALYSIS_PLAN.md').decode()
SKIP = {'00_template.py', '21_fig1_overview.py', '22_fig2_support.py', '23_fig3_size.py', '24_fig4_time.py', '25_figS1_sensitivity.py',
        '31_build_docx.py', 'ajhg_genes.py', 'pilot_f1.py', 'fig1_v2.py', 'fig2_v2.py', 'fig3_v2.py', 'fig4_v2.py', 'lay_colour.py', 'colour_page.json'}
README = """Analysis code for "Calibrating 212,710 published variant-trait claims against genome-wide evidence"

ANALYSIS_PLAN.md   the analysis plan with its dated amendments; amendment 4 (8 November 2026) is the confirmatory plan
m5/                candidate abstracts from the PubMed baseline (01, 02, pmparse), linkage-disequilibrium proxies and the matched-null
                   sample (run_04_ld, run_06_null), and the language-model pass over the abstracts (05_llm_read; run_05)
03-06              GWAS Catalog, linkage disequilibrium, ontology tables and the mapping of disease identifiers to ontology terms
07_claims          claim records from the model output
08_match           each claim against the GWAS Catalog: testability, support at any time and earlier, matched null
09, 10, 12, 12b    the blinded reading check: sample, reading page, collection of answers and scoring
11_confirm         the five primary estimates, secondary breakdowns and the five planned sensitivity analyses
13_additional      additional analyses: genome-wide studies and meta-analyses excluded; agreement of allele and direction
14_clinvar         additional analysis: ClinVar risk-factor, association and protective variants (Data S3)
16_robustness      additional analyses: narrower trait matching; ClinVar conditions
30, 15, 32, 33     numbers quoted in the text; supplemental tables and Data S1-S3
34_lookup_data     data files of the online lookup page
40, figdata, v2/   variant positions and the figures (v2/fig1-4_v3, v2/figS1-S3_v3; colour roles in v2/roles.py and palettes.json)
run_chain.sh       the order in which the confirmatory steps were run

The scripts read a local project database built by steps 03-08 from the public inputs listed in the article; paths at the top of
each script point to that database and to a shared helper for recording inputs and outputs, and have to be adapted before re-use.
"""
with zipfile.ZipFile(O / 'Data_S2_analysis_code.zip', 'w', zipfile.ZIP_DEFLATED) as z:
    for f in sorted((F.ROOT / 'scripts').rglob('*')):
        if (f.is_file() and f.suffix in ('.py', '.sh', '.html', '.js', '.tsv', '.json') and '__pycache__' not in f.parts and f.name not in SKIP
                and not f.name.startswith('4') or f.name == '40_positions.py') and f.is_file() and '__pycache__' not in f.parts \
                and f.suffix in ('.py', '.sh', '.html', '.js', '.tsv', '.json') and f.name not in SKIP:
            z.write(f, 'code/' + str(f.relative_to(F.ROOT / 'scripts')))
    z.writestr('code/ANALYSIS_PLAN.md', plan)
    z.writestr('code/README.txt', README)
print(len(S3), 'ClinVar variants in Data S3;', (O / 'Data_S2_analysis_code.zip').stat().st_size // 1000, 'kB code')

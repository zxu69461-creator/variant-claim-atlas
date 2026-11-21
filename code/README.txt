Analysis code for "Variant-Trait Claim Atlas: checking 201,720 published variant-trait claims against genome-wide evidence and ClinVar"

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
34_lookup_data     data files of the web tool (variants, traits, claims)
35_lookup_extras   data files of the web tool for screening papers and texts (trait phrases, papers)
40, figdata, v2/   variant positions and the figures (v2/fig1-4_v3, v2/fig5_tool, v2/figS0_design, v2/figS1-S3_v3; colour roles in v2/roles.py and palettes.json)
run_chain.sh       the order in which the confirmatory steps were run

The scripts read a local project database built by steps 03-08 from the public inputs listed in the article; paths at the top of
each script point to that database and to a shared helper for recording inputs and outputs, and have to be adapted before re-use.

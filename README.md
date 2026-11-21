# Variant–Trait Claim Atlas

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24314560.svg)](https://doi.org/10.5281/zenodo.24314560)

Published claims that a named genetic variant is associated with a disease or trait (201,720 disease or trait claims from 71,577 PubMed abstracts), each set
against genome-wide significant associations in the NHGRI-EBI GWAS Catalog, with frequency-matched random variants as the chance
expectation.

**Online tool:** https://zxu69461-creator.github.io/variant-claim-atlas/ (runs in the browser; free, no installation, registration,
cookies or tracking; tested in current Chrome, Safari and Firefox on desktop and phone-sized screens)

- **Look up** a variant (rs1800795), a gene (MTHFR) or a disease or trait (asthma): support of its published claims against the matched
  null, by stated sample size and by year, with links to dbSNP, the GWAS Catalog, ClinVar, LitVar 2.0, Ensembl and Open Targets, and
  every claim linked to PubMed and downloadable as CSV.
- **Check a claim**: one variant and one trait (rs7903146 and T2D; names, parts of names and common abbreviations) give a verdict
  (supported, partly supported, not supported or not testable), the year of the first genome-wide report and the matching claims.
- **Screen a paper**: a PubMed ID (title and abstract from Europe PMC) or any pasted text; every rsID in it is found, the diseases and
  traits it names are matched to the published claims about that variant, and each variant–trait pair is checked.
- **Batch check**: paste or upload up to 10,000 rsIDs or variant–trait pairs (`examples/batch_example.csv`) and download the results.
- **ClinVar risk variants**: the 239 variants that ClinVar lists as risk factors, associations or protective variants, with their
  genome-wide support and whether genome-wide evidence supports the listed condition.
- **Help**: tutorial with worked examples, definitions of every number, the verdict rules, FAQ and version history.
- **Programmatic access**: the tool reads static files (`data/*.json`, listed under Help) that can be fetched directly; clients in
  Python (`examples/vca_client.py`, standard library only) and R (`examples/vca_client.R`, jsonlite) reproduce the views and the verdict:

```
python examples/vca_client.py check rs7903146 "type 2 diabetes"
python examples/vca_client.py batch my_variants.txt -o results.csv
Rscript examples/vca_client.R check rs1800795 "type 2 diabetes"
```

Lists, files and texts pasted or uploaded into the tool are processed in the browser and never sent anywhere; screening a paper by PubMed
ID sends only that ID to Europe PMC. The tool is maintained by the authors for at least five years; every data update is released as a
new version archived at Zenodo, where earlier versions stay available. Problems and questions: GitHub issues.

## Contents

| Path | What it is |
|---|---|
| `index.html`, `app.js`, `data/` | the web tool and its data (per variant, per trait, every claim, trait phrases, papers) |
| `examples/` | Python and R clients and an example batch file |
| `atlas/Data_S1_variant_trait_claim_atlas.tsv.gz` | the full atlas: one row per paper, variant, linked trait and population group, with the genome-wide comparison for testable claims |
| `atlas/Data_S3_clinvar_risk_variants.tsv` | variants listed in ClinVar as risk factors, associations or protective variants, with their genome-wide comparison |
| `code/` | analysis code and the dated analysis plan (`code/README.txt` lists the steps) |

## Definitions

A reported association is **supported** when the claimed variant, or a proxy in linkage disequilibrium (r² ≥ 0.8 in the 1000 Genomes
European or East Asian reference), has an association at P < 5 × 10⁻⁸ with the same trait or a broader or narrower Experimental Factor
Ontology term, in a GWAS Catalog paper other than the claim's own. The **matched null** is the share of frequency-matched random variants
that would be supported for the same traits. A claim is **testable** when its trait has any association in the GWAS Catalog. Support
does not require the same allele or direction. Claims were read by a language model run locally (Qwen3.5-9B) and checked on a blinded
sample of held-out papers.

Inputs: PubMed 2026 baseline, PubTator3, GWAS Catalog (release of 29 October 2026), 1000 Genomes phase 3, EFO v3.95.0, Mondo, MeSH 2026,
ClinVar (release of 16 November 2026).

## Citation

Hu Y, Xu Z, Jiang S, Huang L, Pan M, Luo S, Zhou H, Guo X. Variant–Trait Claim Atlas: checking 201,720 published variant–trait claims
against genome-wide evidence and ClinVar. Manuscript under review, 2026. Archived release: https://doi.org/10.5281/zenodo.24314560

## Licences

Code: MIT (`LICENSE`). Atlas data and lookup data: CC BY 4.0 (`LICENSE-DATA`). The atlas records what each abstract claimed; the
abstracts themselves are not redistributed.

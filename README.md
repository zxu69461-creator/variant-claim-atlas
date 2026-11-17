# Variant–Trait Claim Atlas

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24314560.svg)](https://doi.org/10.5281/zenodo.24314560)

Published claims that a named genetic variant is associated with a disease or trait (201,720 claims from 71,577 PubMed abstracts), each set
against genome-wide significant associations in the NHGRI-EBI GWAS Catalog, with frequency-matched random variants as the chance
expectation.

**Online lookup:** https://zxu69461-creator.github.io/variant-claim-atlas/ — search an rsID (rs1800795), a gene symbol (MTHFR) or a disease
or trait (asthma) to see how often its published claims have genome-wide support, by study size, with every claim linked to PubMed.

## Contents

| Path | What it is |
|---|---|
| `index.html`, `app.js`, `data/` | the static lookup page and its data (per variant, per trait, and every claim) |
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

Hu Y, Xu Z, Jiang S, Huang L, Pan M, Luo S, Zhou H, Guo X. Calibrating 212,710 published variant–trait claims against genome-wide
evidence exposes a split extending into ClinVar. Manuscript submitted, 2026. Archived release: https://doi.org/10.5281/zenodo.24314560

## Licences

Code: MIT (`LICENSE`). Atlas data and lookup data: CC BY 4.0 (`LICENSE-DATA`). The atlas records what each abstract claimed; the
abstracts themselves are not redistributed.

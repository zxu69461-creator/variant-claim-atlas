> Analysis plan of the study as written, with the dates of its amendments (4, 7 and 8 November 2026). The text is unchanged except that one sentence about local file storage was removed and one note in square brackets was added. "Reader" is the language model in amendments 1-3 and section 4.1 names the human reader; "arbiter" is the source of genome-wide evidence. The comparison with FinnGen considered in amendment 1 (arbiter B) was not part of the confirmatory plan (amendment 4) and was not carried out.

# varclaims — charter (2026-11-04, written before any claim was read or matched)

## Object
Every published claim that a named dbSNP variant (rsID) is associated with a disease or trait, extracted from PubMed titles and abstracts
(2026 baseline, all years), each linked to the association of the same variant (or an LD proxy from the 1000 Genomes reference) with the
matched endpoint in FinnGen R13, plus paper-level metadata (year, journal, country, study population, sample size where stated,
retraction notices).

## Why it can conflict with current understanding
The candidate-variant literature is widely held to have failed to replicate (Siontis 2010 Eur J Hum Genet, 291 associations before 2008),
and GWAS is held to have replaced it. Whether candidate-variant claims kept being published after 2010, how many there are, which of them a
modern biobank supports, and whether support depends on era, population, sample size or venue, has not been measured across the literature.

## Ruler (fixed before matching)
- Claim unit: (PMID, rsID, disease/trait, direction where stated, population), from the abstract; one record per distinct tuple.
- Disease/trait mapped to a FinnGen R13 endpoint by an explicit, pre-written mapping (PubTator3 MeSH disease id -> FinnGen endpoint);
  unmapped claims are reported, never counted as unsupported.
- Support: the variant (or best LD proxy, r2 >= 0.8 in the matching 1000G population) has P < 5e-8 for the endpoint (strict) or nominal
  P < 0.05 with the claimed direction where stated (lenient); each claim is read against a matched null (frequency-matched random variants,
  same endpoint), never pooled across claims.
- Claims from non-European populations are reported separately.

## Rules
Exploration on a development split (md5(str(pmid)) last hex digit in {0, 1}) is free; confirmatory numbers come only after a dated plan
is committed. No change of ruler after seeing results; additions are labelled post hoc. No meta-analytic framing, no pooling of effect
estimates across studies.

## Amendment 1 (2026-11-04, before any claim was matched to an arbiter)
FinnGen R13 full summary statistics are local for 32 endpoints only, so a second arbiter is added and the FinnGen rule is narrowed.
- Arbiter A, GWAS Catalog (release 2026-10-29, all ancestries): a claim is replicated if the same rsID, or a proxy with r2 >= 0.8 in the
  1000 Genomes super-population matching the claim's population (EUR or EAS reference; other populations use EUR and are reported
  separately), has a GWAS Catalog association with P < 5e-8 for the claim's trait, from a study other than the claim's own paper.
  Trait match: the claim outcome mapped to EFO (PubTator3 MeSH id -> EFO MeSH cross-reference; claims without a MeSH-mapped EFO term are
  reported as unmapped); match if the GWAS EFO term equals the claim term, or one is an ancestor of the other within the EFO disease or
  measurement branch.
  Testability: the claim's EFO term (or a descendant) has at least one GWAS Catalog study; claims on untestable traits are reported
  separately, never counted as unreplicated.
  Era: support is read both from all studies and from studies published before the claim (prior support) and after it (later support).
- Arbiter B, FinnGen R13: only for claims whose trait maps to a locally available full endpoint; lookup of the rsID (or proxy) with the
  claimed direction where stated; used for powered nulls (P >= 0.05 with the endpoint's case count reported).
- Matched null for every claim: frequency-matched random rsIDs from the same candidate gene region are not available for most claims, so
  the null is the share of all GWAS Catalog-tested rsIDs with an association for the same EFO term, at the same minor-allele-frequency bin
  (1000G), computed per claim and averaged over claims.

## Amendment 2 (2026-11-04, before any claim was matched to an arbiter)
The claim-trait mapping in amendment 1 is widened to recover diseases that PubTator3 tags with MeSH or OMIM ids that EFO does not
cross-reference. Steps are applied in order and the first that yields a term is used; every claim keeps a mapping-route label.
1. Direct: MeSH id -> EFO term through the EFO MeSH cross-reference (route `efo_xref`).
2. MONDO: MeSH or OMIM id -> MONDO term through a MONDO cross-reference marked as equivalence (`MONDO:equivalentTo`) when such marks
   exist, otherwise any MONDO cross-reference to that id; the MONDO term is used as is when EFO contains it, otherwise through its EFO
   cross-reference (route `mondo_xref`).
3. Nearest mapped MeSH ancestor: walk up the MeSH 2026 tree from every tree number of the descriptor, one level at a time, and take the
   closest level at which any ancestor maps by step 1 or 2 (all ancestors at that level are kept). Ancestors at the first tree level (for
   example C14, Cardiovascular Diseases) are not used (route `mesh_ancestor`).
Claims mapped by step 3 are broader than the outcome the authors named; they are reported with the rest and again separately, and the
primary replication estimate is also given restricted to steps 1 and 2. A descriptor or OMIM id still unmapped after step 3 stays
unmapped. Mapping uses only ontology files (EFO, MONDO, MeSH 2026), never the GWAS Catalog associations, so it cannot be tuned to the
matching result.

## Amendment 3 (2026-11-07, after the reader pass of the development split, before any claim was matched to an arbiter)
Operational rules that amendments 1 and 2 left open. Only parse statistics of the reader output have been looked at.
- Claim records: reader entries whose rsID is one of the rsIDs listed for that abstract and whose outcome_type is disease, trait or
  severity (drug_response entries are kept as a separate set). One record per (PMID, rsID, linked trait, population group); a record
  is positive when any of its entries says associated = yes, null when all say no, otherwise unclear.
- Outcome linking, in order: (a) a PubTator3 disease id of the same abstract whose tagged name equals the outcome text, or is contained
  in it or contains it (case-insensitive; the longest such name wins); (b) a PubTator3 disease id whose name shares >= 50% of the
  outcome's word tokens (Jaccard, tokens of >= 3 letters); (c) an EFO term whose label or exact synonym equals the outcome text
  (case-insensitive); (d) the only PubTator3 disease id of the abstract when it has exactly one. The id is then mapped to EFO by
  amendment 2. The route is kept on every record and the estimate is also given without route (d).
- Population group: EAS when the population text names Chinese, Han, Japanese, Korean, Taiwanese, Vietnamese, Thai, Malay, Mongolian
  or East Asian participants; EUR when it names European, Caucasian or white participants or a European nationality; otherwise
  `other`. When the population is not stated, the first author's country decides (China, Japan, South Korea, Taiwan, Hong Kong ->
  EAS; European countries -> EUR; others -> `not_stated`), and the record is labelled inferred. `other` and `not_stated` use the EUR
  LD reference and are reported separately.
- Trait match within the disease or measurement branch: both terms descend from EFO disease (EFO_0000408), measurement
  (EFO_0001444) or phenotype (EFO_0000651) (MONDO terms count as disease), and one equals or is an ancestor of the other; the
  branch roots themselves never match.
- Matched null: 20,000 rsIDs drawn at random (seed 20261107) from the current rsIDs that appear in the GWAS Catalog, with LD proxies and
  1000G EUR minor-allele frequency computed exactly as for claim variants. For each claim the null is the share of null rsIDs in the
  claim variant's MAF bin (0.01-0.05, 0.05-0.10, 0.10-0.20, 0.20-0.30, 0.30-0.40, 0.40-0.50; below 0.01 or absent from 1000G form
  their own bin) whose own variant or proxy has a GWAS Catalog hit matching the claim's trait. Null rsIDs are themselves hits for some
  trait, so this null is higher than one drawn from random genomic variants (conservative).
- Clarification (2026-11-07, before matching): EFO_0000408 is obsolete in the EFO release used; the disease root is MONDO_0000001
  (disease), which is used in its place. Roots never matched: MONDO_0000001, EFO_0001444, EFO_0000651.

## Amendment 4 (2026-11-08) — confirmatory plan, written before any held-out claim was matched to an arbiter
State when written: the reader pass is complete (82,013 abstracts); claim records exist for all abstracts; only development records
(md5 last hex digit 0 or 1) have been matched (scripts/08). Directions below are those seen on the development split; the held-out
split (all other records) is analysed once, by scripts/10 and 11 as committed with this amendment, after the reader check passes.

### 4.1 Reader check (gate before matching the held-out split)
- Sample (seed 20261108), from held-out records with outcome_type disease, trait or severity and a linked trait: 80 positive and 40 null
  records, at most one per abstract; plus 30 further held-out abstracts with 1-6 listed rsIDs for recall. Records are drawn uniformly
  within status; the reported precision is per status (no reweighting needed).
- Reader-facing material: title, abstract, the rsID and the reader's outcome text for the record (to define the item); never the
  reader's status, direction, population or any arbiter result. Recall items show the abstract and its listed rsIDs only.
- Questions: (Q1) does the abstract report that this variant is associated with this disease or trait — associated / reported as not
  associated / no conclusion for this pair / the pair is not studied in the abstract (extraction error); (Q2) for recall abstracts, for
  each listed rsID, does the abstract report a statistically supported association of that variant with any disease or trait (yes / no).
- Measures: positive precision = share of reader-positive records judged associated; null precision = share of reader-null records
  judged not associated; pair error = share judged not studied; recall = share of rsIDs judged yes that have a reader-positive record
  in that abstract. Wilson 95% intervals.
- Gate: the held-out split is matched only if positive precision >= 0.85 and null precision >= 0.75 (point estimates). If either
  fails, no held-out matching is run; the prompt is revised under a new amendment and the full pass re-read. Recall is reported, not
  gated (missed claims shrink the object; they do not enter the replication estimates).
- Reader: the first author (one reader) [Z.X.; the order of authors was set later]. The check measures reading, not replication; its sample is excluded from nothing.

### 4.2 Analysis set
Held-out records, outcome_type not drug_response, status positive or null, testable (amendment 1). Matching is scripts/08 unchanged
except for the split filter. Era by publication year; "new" = no prior support (rep_prior false). Clustering unit: PMID.

### 4.3 Confirmatory estimates (primary family, Holm-adjusted over C1-C5, two-sided alpha 0.05)
- C1 Replication above chance: share of positive claims replicated (any time) minus mean matched null share; PMID cluster bootstrap
  (2,000 resamples, seed 20261108). Expected from development: 22% vs 2%.
- C2 Sample size: logistic GEE (exchangeable, PMID clusters) of replicated on log10(stated n) among positive claims with n; odds ratio
  per 10-fold n. Expected: increasing.
- C3 Restatement: among replicated positive claims published 2005 or later, logistic GEE of rep_prior on publication year (per 5
  years). Expected: increasing.
- C4 New-claim yield, fixed window: among new positive claims published 2005-2020, replicated with first GWAS hit within 5 years of
  publication; logistic GEE on publication year (per 5 years). Expected: decreasing. Also reported: the same share for new null
  claims and the matched null, per era (2005-09, 2010-14, 2015-17, 2018-20).
- C5 Null claims contradicting earlier GWAS evidence: among authors' null claims published 2005 or later, logistic GEE of rep_prior on
  publication year (per 5 years). Expected: increasing.

### 4.4 Secondary (reported with intervals, no hypothesis test)
Replication by design (family, case-control, cohort, meta-analysis, GWAS), by population group (EUR, EAS, other, not stated), by
sample-size bin (<500, 500-1,999, 2,000-9,999, >=10,000) and by era; for each: claim share, replicated share, matched null.

### 4.5 Pre-specified sensitivity analyses (C1-C5 repeated)
(a) without outcome-link route (d); (b) without mapping route mesh_ancestor; (c) only populations stated in the abstract (not inferred);
(d) without retracted papers; (e) one weight per paper (each PMID's claims weighted 1/k).

### 4.6 Not done
No pooling of effect estimates across studies, no meta-analytic summary of the claims. Development and held-out results are reported
separately; held-out estimates are the confirmatory numbers.

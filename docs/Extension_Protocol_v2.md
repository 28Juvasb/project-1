# Extension Protocol v2.0 — WP1 to WP3

**Written:** 2026-10-12. **Status:** DRAFT for OSF registration. Register this file on OSF Registries (Open-Ended Registration) **before** any expression data from GSE50397 or GSE76895 are opened, and before WP2 screening starts. Until it is registered, only file checksums and sample metadata may be inspected.
**Relation to the primary analysis:** the primary analysis (`docs/Methodology_Protocol_v1.md`, tags `protocol-v1` to `analysis-complete-v1`) is complete and locked. This extension adds new, separately registered work and changes none of its numbers. Anything not stated here follows the primary protocol and CLAUDE.md.
**Decision authority:** Juvas (project lead). Supervisor review: Dr. Thirumurugan.
**Disclosure to put in the OSF record:** the primary protocol was locked in git on 2026-10-05 (tag `protocol-v1`) and deposited on OSF later; only this extension is registered before its data are opened.

## Aim and questions
Aim: quantify how duplicated public islet datasets and selection outside cross-validation inflate reported accuracy of T2D gene classifiers.
- Q1 (WP1): which public human islet T2D series on Affymetrix arrays share donors or arrays?
- Q2 (WP2): how often do published studies that use these series treat overlapping series as independent validation or replicate, or merge them, and how often is gene selection inside cross-validation?
- Q3 (WP3): for the same discovery data, how much larger is the gap between reported and independent performance for the common-practice pipeline than for our nested pipeline?

## WP1 — Overlap map
- **Series:** GSE41762, GSE38642, GSE50397, GSE20966, GSE76894, GSE76895, GSE25724 (and any other public human islet T2D series found by a GEO search done on the registration date; list every one searched in the log). Include a series only if raw CEL files are public.
- **Primary criterion:** identical decompressed-CEL MD5 between series (an identical raw file is the same array). **Secondary:** covariate matching (sex + diabetes group + age + BMI, exact; near-match ±1 year and ±0.1 BMI reported separately, as in `scripts/02_donor_overlap.R`); series-header provenance. No expression values are read for this step.
- **Pilot result disclosed before registration (2026-10-13):** a metadata-and-checksum pilot found 67 of 89 GSE50397 arrays byte-identical to GSE41762 arrays (56 also in GSE38642), with deposited ages shifted by 1–2 years for most shared donors, which is why exact covariate matching is only secondary. GSE54279 (128 Uppsala islets, HbA1c only) is added for the checksum comparison only. No expression values of GSE50397, GSE54279 or GSE76895 were read. GSE50397 has no diabetes-status field and is not used in WP3, so its labels are not inferred.
- **Outputs:** pairwise table (n shared arrays, n covariate-matched donors, n label conflicts), Figure 1.
- **Contact:** the submitters of GSE41762 and GSE38642 are asked to confirm shared donors and T2D definitions. Their answers are reported as received.

## WP2 — Literature audit
- **Searches (run once, dates logged):** Europe PMC and PubMed. Population: studies of human islet or pancreas expression in T2D that name at least one WP1 series by accession. Query: each accession as free text OR ACCESSION_ID, combined with OR; time window 2012-01-01 to the search date.
- **Screening:** title/abstract then full text, by two independent reviewers (project lead plus one labmate); disagreements resolved by Dr. Thirumurugan. Exclusions: reviews, conference abstracts, non-human data, no full text available (counted and listed).
- **Extraction per paper:** series used and their roles (discovery, validation, replicate, merged); whether any used pair is an overlap pair from WP1; whether the paper states the overlap; feature selection inside or outside cross-validation (yes / no / not reported); reported performance (AUC, CI, validation type); code or data availability; PROBAST+AI development-analysis signalling questions (item list fixed before extraction).
- **Primary outcome:** proportion of included papers that use an overlap pair as independent validation or replicate, or merge it. **Secondary:** proportion with selection inside cross-validation; proportion reporting an AUC CI; median reported AUC by validation type.
- **Agreement:** Cohen's kappa for each extraction item on all papers. Report counts with exact binomial 95% CIs. No inferential tests across subgroups beyond these.
- **Tone rule:** results describe practices, never authors.

## WP3 — Leakage demonstration
- **Discovery data:** GSE41762 (as in the primary analysis, 77 samples, post-QC data from stage 1).
- **Common-practice pipeline (fixed here, rebuilt from the audited papers' methods before looking at any WP3 result):** limma on the full discovery cohort, keep senescence-union genes with nominal p < 0.05, LASSO (alpha 1, lambda.1se) on the full cohort, logistic model on selected genes, apparent AUC on the same cohort, then 'validation' on GSE38642 (post-QC, normalised independently). Variants used by audited papers may be added as pre-listed pipelines in a dated addendum before running.
- **Our pipeline:** the frozen nested pipeline and frozen final model from the primary analysis. No refit.
- **External test:** GSE76895 (IMIDIA pancreatectomy cohort; groups ND and T2D only, IGT and type 3c excluded; the groups and n are confirmed from the manifest before scoring). Each frozen model (ours; each common-practice pipeline's model) is scored once; a lock file is written as in `scripts/09_external_validation.R`.
- **Metrics:** AUC with 95% CI (stratified bootstrap, 2000). **Primary contrast:** (reported/apparent AUC minus independent AUC) for the common pipeline versus ours, with a bootstrap 95% CI for the difference. **Secondary:** AUC on GSE38642 for the common pipeline (invalid validation) versus GSE76895 (independent).
- **Caveat stated up front:** GSE76895 islets were laser-capture microdissected from pancreatectomy patients, so tissue and isolation differ from the discovery islets.
- **Rules:** no tuning on any test set; no pipeline change after seeing WP3 results; all deviations logged in `logs/deviations.md` with date, reason and original value.

## Deferred
Molecular docking is not part of this study (no candidate passed protocol §5.4). Ligand files remain in `ligands/`.

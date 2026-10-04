# Methodology Protocol v1 — Cross-Tissue Senescence Hub Genes in Type 2 Diabetes

**Version:** 1.0, locked 2026-10-05 (git tag `protocol-v1`)
**Decision authority:** Juvas (project lead). Supervisor: Dr. Thirumurugan.
**Status of data at lock time (stated plainly):** all GEO files had been downloaded and the *sample metadata* (counts, group labels, covariates) inspected. **No expression values had been normalised, plotted or analysed.** Every analytical choice below was fixed before seeing any expression result.

This is the execution spec referenced by `CLAUDE.md` and `IMPLEMENTATION_PLAN.md`. Section numbers are fixed; scripts cite them. Any change after lock goes to `logs/deviations.md` with date, reason and original value (CLAUDE.md rule 7). Decisions and their rationale are summarised in `docs/pre_analysis_plan.md`.

---

## 1. Design and cohorts

Two parallel arms (islet, adipose). Each arm: discovery → nested-CV-validated hub-gene panel → frozen-model external validation. Then cross-tissue comparison, enrichment, single-cell specificity (islet only) and docking.

| Role | Accession | Tissue | Platform | Groups used | n |
|---|---|---|---|---|---|
| Islet discovery | GSE41762 | Whole islets | GPL6244 (HuGene 1.0 ST) | 57 ND / 20 T2D | 77 |
| Islet sensitivity (direction only) | GSE38642, non-overlapping donors only | Whole islets | GPL6244 | determined in Task 2 | ~12 |
| Islet external validation | GSE20966 | LCM beta-cells | GPL1352 (U133 X3P) | 10 ND / 10 T2D | 20 |
| Islet single-cell | GSE81608 | Islet cells | GPL16791 (SMARTer) | 651 ND / 949 T2D cells | 1600 cells |
| Adipose discovery | GSE27951 (GPL570 subseries) | Subcutaneous | GPL570 | 11 ND / 12 T2D (10 IGT excluded) | 23 |
| Adipose external validation | GSE71416 | Omental, morbidly obese | GPL570 | 6 ND / 14 T2D | 20 |
| Adipose sanity (direction only) | GSE16415 | Visceral | GPL2986 | 5 ND / 5 T2D | 10 |

**Group definitions.** GSE27951: the `clinical status` field (NGT → ND, DM → T2D, IGT excluded). Three samples have title/field discordance (GSM691135, GSM691139, GSM691140); the field is used, and a sensitivity rerun excludes GSM691135 and GSM691140 (GSM691139 is IGT by field and already excluded). GSE71416: sample titles (no disease field in GEO). All other series: the GEO status field.

**Fixed rules** (from CLAUDE.md): discovery cohorts never merged; validation sets scored once with a frozen model; no batch correction across cohorts; nested-CV AUC is the headline and is always reported beside the apparent AUC; all AUCs carry a 95% CI (`pROC::ci.auc`, 2000 stratified bootstrap replicates); `set.seed(42)`.

---

## 2. Preprocessing

### 2.1 Normalisation
Each GSE is normalised **separately**, from raw files where available.
- GPL6244 (GSE41762, GSE38642): `oligo::rma()` on CEL files, transcript-cluster level (`target = "core"`). GSE41762: all 77 CELs normalised together; the submitters' donor 1–48 / 49–80 split is kept as the covariate `submitter_norm_batch`. Before use, the batch split parsed from titles is verified against the GEO series description; if it cannot be verified, that is logged and the covariate is still used as parsed.
- GPL570 (GSE27951, GSE71416): `affy::rma()` on CEL files. GSE27951: only the 23 included samples (11 ND + 12 T2D) are normalised; the excluded IGT samples are not used anywhere.
- GPL1352 (GSE20966): `affy::rma()` on CEL files with `u133x3pcdf`.
- GPL2986 (GSE16415): submitter matrix (`GSE16415_all_data_no_controls.txt.gz`). If values are not already log-scale (max > 50), apply log2(x + 1); then `limma::normalizeBetweenArrays(method = "quantile")`.
- No ComBat or any cross-cohort batch correction.

### 2.2 Probe-to-gene collapse
- Annotation: `hugene10sttranscriptcluster.db`, `hgu133plus2.db`, `u133x3p.db`; GPL2986 from the GPL annotation table.
- Drop control probes (Affymetrix `AFFX-` prefixes; HuGene transcript clusters with category ≠ main) and probes mapping to no or several symbols.
- Collapse to HGNC symbol with `WGCNA::collapseRows(method = "MaxMean")`.
- **Expression filter (per dataset):** a gene is "expressed" if its log2 expression exceeds the dataset's 25th percentile of all values in ≥ 20% of samples. The expressed set defines the enrichment background (§5.2) and the hypergeometric universe (§5.1). Filtering for modelling happens inside folds (§3).

### 2.3 Quality control
- RLE and NUSE (`affyPLM::fitPLM` for affy; `oligo::fitProbeLevelModel` for oligo). **Exclude arrays with median NUSE > 1.05.** No NUSE for GSE16415 (no probe-level data); RLE-type boxplots only.
- PCA (top 5000 genes by variance) and hierarchical clustering (1 − Pearson, average linkage). **Samples > 3 SD from the mean on PC1 or PC2 are flagged, not excluded.** Plots go to the project lead, who decides; any exclusion is logged with its reason.
- Every exclusion goes to `logs/excluded_samples.csv` (gsm, gse, reason, date).

### 2.4 Sex check and covariates
- Expression-based sex: XIST vs mean(RPS4Y1, DDX3Y, KDM5D, UTY where present), classified by k-means (k = 2).
- Where sex is recorded (GSE41762, GSE38642, GSE20966 if available), mismatches are reported and flagged. They are not auto-excluded; the project lead decides.
- **GSE27951 has no sex field in GEO.** Expression-inferred sex is used as its sex covariate.
- **limma model, both arms:** `~ status + age + sex`; GSE41762 also includes `submitter_norm_batch`. Pre-specified sensitivity analysis: add BMI.

---

## 3. Nested cross-validation engine (stage 3)

Implemented once as `nested_cv(expr, pheno, geneset, params)`, called for both arms.

**Fold structure (stratified by status):**
| Arm | Outer | Inner (LASSO, SVM-RFE, ridge) |
|---|---|---|
| Islet (n=77, 20 T2D) | 10-fold × 5 repeats = 50 folds | 10-fold |
| Adipose (n=23, 12 T2D) | 5-fold × 10 repeats = 50 folds | 3-fold |

**Inside each outer fold, using the training split only:**
1. **Screen:** limma (§2.4 covariates) on training samples; keep senescence-union genes (§3.1) with nominal p < 0.05. **Fallback:** if fewer than 10 genes pass, take the top 30 senescence-union genes by p. The fallback count is reported.
2. **LASSO:** `glmnet`, alpha = 1, lambda = `lambda.1se` from inner CV; selected = nonzero coefficients.
3. **SVM-RFE:** `caret::rfe` with a linear SVM, inner CV, subset sizes 1–min(30, p); selected = the optimal subset.
4. **Random Forest + Boruta:** `Boruta(maxRuns = 500, ntree = 1000)`; selected = "Confirmed" (Tentative is not selected).
5. **features_k** = intersection of LASSO, SVM-RFE and Boruta. **Fallback 1:** genes selected by ≥ 2 of 3. **Fallback 2:** top 5 LASSO genes by abs(coefficient); if LASSO selected nothing, top 5 screen genes by p. The fallback level used in each fold is recorded.
6. **Classifier:** ridge logistic regression (`glmnet`, alpha = 0, lambda by inner CV) on features_k, with features z-scored using training-fold means and SDs only.
7. Predict on the held-out fold (scaled with training parameters). Store out-of-fold probabilities.

**Outputs:** pooled out-of-fold AUC with 95% CI (per repeat, then averaged across repeats for the point estimate; CI from the pooled predictions via 2000 stratified bootstraps). Also: per-gene selection frequency across the 50 folds (`results/stage3/stability_selection.csv`), the fallback frequency table, and per-fold feature lists.

**Mandatory self-tests (assertions; failure stops the run):**
1. No sample index appears in both train and test in any fold.
2. The matrix passed to limma in each fold has exactly the training-fold number of columns.
3. **Permutation test:** 20 independent label shuffles, each running the full pipeline with 1 repeat. **Pass:** mean AUC in [0.45, 0.55] and 95th percentile < 0.65. Failure means leakage; stop and report.

### 3.1 Senescence gene universe
Union (not intersection) of SenMayo, CellAge and GenAge-human, mapped to current HGNC symbols (`limma::alias2SymbolTable` with `org.Hs.eg.db`); unmapped symbols are reported and dropped. SenMayo: the 125-gene list from Saul et al. 2022 Supplementary Data 1 if retrievable; otherwise the MSigDB `SAUL_SEN_MAYO` 124-gene copy, with the discrepancy named.

### 3.2 Final panel and apparent AUC
The full pipeline (steps 1–6) runs once on the complete discovery cohort to give the **final panel**. The apparent AUC of this model on its own training data is reported **only beside** the nested-CV AUC, labelled "optimistic (apparent)".

### 3.3 Descriptive DEG table
Full-cohort limma, reported for description only and **never fed into selection**: BH FDR < 0.05 and abs(log2FC) ≥ 0.263 (1.2-fold).

---

## 4. Co-expression network analysis

### 4.1 WGCNA
**Role (fixed):** descriptive and corroborative. WGCNA does **not** choose the final panel (that is §3.2, which the nested CV measures). The paper reports which panel genes are hubs of trait-associated modules.
- Input: the full normalised discovery matrix (not senescence-filtered), top 5000 genes by MAD.
- Network: signed hybrid; `blockwiseModules` with `minModuleSize = 30`, `mergeCutHeight = 0.25`, `deepSplit = 2`, `TOMType = "signed"`, `corType = "bicor"` (`maxPOutliers = 0.05`).
- Soft threshold: the lowest beta in 1–20 reaching scale-free R² ≥ 0.85. If none does, report the curve and **stop**; never pick a beta silently.
- Module–trait: eigengene vs status, age, sex and BMI (Pearson). Modules kept if abs(r) > 0.3 and p < 0.05 vs status. Modules also correlated (p < 0.05) with age, sex or BMI are flagged as possibly confounded.
- Hubs: kME > 0.8 **and** abs(GS) > 0.2, intersected with senescence-union genes that are DEGs in §3.3.
- **Adipose stability (n=23):** 100 bootstrap resamples, each re-running module detection; preservation of each reference module via `modulePreservation` Zsummary (> 2 = preserved). A module is retained only if preserved in **> 70%** of resamples. If none clears 70%, the adipose WGCNA is labelled **exploratory**. The threshold is never lowered.

---

## 5. Integration

### 5.1 Cross-tissue comparison
- Intersect the two final panels by HGNC symbol → shared, islet-specific and adipose-specific.
- **Primary test:** permutation null. Repeat 1000 times: shuffle labels within each discovery cohort independently, rerun the screen + LASSO selection (§3 steps 1–2, full cohort), and record the overlap size. Empirical p = (count(null ≥ observed) + 1) / 1001. Boruta and SVM-RFE are excluded from the permutation loop for compute reasons; this simplification is fixed here, not a deviation.
- **Secondary test:** hypergeometric, universe = senescence-union genes expressed (§2.2) in both discovery datasets.
- **Verdict rule:** the overlap is reported as a finding only if the primary p < 0.05. Otherwise the text says the overlap is consistent with chance.
- UpSet plot.

### 5.2 Enrichment
`clusterProfiler` GO (BP, MF, CC) and KEGG, run separately on the shared, islet-specific and adipose-specific sets. Universe = expressed genes (§2.2) in the relevant tissue (union of both tissues for the shared set). BH-adjusted p < 0.05; minGSSize 10. A set with < 10 genes is declared underpowered and annotated narratively (gene-level function, senescence evidence) instead of tested.

### 5.3 Single-cell specificity (GSE81608)
- QC: MAD-based (3 MADs) on log genes detected, log library size and mitochondrial fraction.
- Cell types: the authors' labels if present in the metadata; otherwise Leiden clustering and annotation by INS, GCG, SST, PPY, PRSS1 and COL1A1.
- For each islet panel gene: detection fraction and mean expression per cell type.
- T2D vs ND: **pseudobulk** (sum counts per donor × cell type; donors with < 10 cells of a type dropped), limma-voom or edgeR. **No per-cell tests.**
- "Beta-cell restricted" = significant (FDR < 0.05) in beta-cell pseudobulk and not significant in alpha or delta. With ~6 T2D donors the tests are low-powered, so absence of a significant difference is reported as inconclusive, not as evidence of no effect.

### 5.4 Docking target selection rule
1. Candidates: shared panel genes; if none, the islet final panel.
2. Each candidate must pass all of:
   - (a) nested-CV selection frequency ≥ 50%
   - (b) same direction of effect in external validation (GSE20966 for islet; GSE71416 for adipose-specific)
   - (c) a usable structure (§6.1)
3. Rank by selection frequency. Dock the top 1, or the top 2 if tied.
4. If no candidate passes, docking is reported as not performed under the rule. The rule is not relaxed.

---

## 6. Molecular docking

### 6.1 Structure
Prefer an experimental PDB structure (resolution < 2.5 Å, no mutation within 8 Å of the site, ≥ 90% sequence coverage of the domain). Otherwise use AlphaFold DB, recording per-residue pLDDT; no binding site may be defined in a region with pLDDT < 70. Prepare with PDBFixer (add missing atoms, protonate at pH 7.4) and meeko/OpenBabel to PDBQT.

### 6.2 Binding site
Use an annotated or co-crystal ligand site if one exists. Otherwise predict with fpocket and P2Rank, and proceed only where the top pockets agree (centroid distance < 4 Å). Grid box coordinates and size are recorded explicitly (box = ligand extent + 8 Å, minimum 20 Å per side).

### 6.3 Docking runs
AutoDock Vina 1.2, exhaustiveness 32 (64 for rutin), 20 poses, 3 runs with seeds 42, 43 and 44. Ligands: naringenin (CID 932), embelin (CID 3218), rutin (CID 5280805) and LSC (structure from the lab; docked only once it is supplied). Ligands are protonated at pH 7.4 (OpenBabel) and prepared with meeko.

### 6.4 Controls (mandatory)
- **Positive:** redock the co-crystallised or known ligand; accepted if RMSD < 2 Å.
- **Negative:** 40 decoys per phytocompound, property-matched DUD-E style from ZINC (MW ± 25 Da, logP ± 1, HBD ± 1, HBA ± 2, rotatable bonds ± 2), docked under identical settings.
- Each phytocompound's best ΔG is reported as a **percentile of its decoy distribution**. A ΔG without its percentile is not a result.

### 6.5 Output table
Ligand, best ΔG, mean ± SD ΔG across seeds, decoy percentile, contact residues (≤ 4 Å), H-bonds, and pose convergence (pairwise RMSD of the top poses across seeds).

---

## 7. External validation
- Models frozen after §3.2: the panel genes, ridge coefficients, and the training means and SDs.
- **Islet → GSE20966** (cross-platform): match by HGNC symbol and z-score each gene within GSE20966. Panel genes absent from U133 X3P are reported. The model is refit on the available genes using the **discovery data only** (no validation labels), and that refit is reported.
- **Adipose → GSE71416** (same platform): z-score within dataset, apply the frozen model.
- AUC with 95% CI (2000 bootstraps) and ROC with confidence bands. **One pass, no tuning.**
- **Direction only, no AUC:** GSE38642 unique donors (islet) and GSE16415 (adipose). Report per-gene direction concordance, with a sign test.
- Pre-stated expectation: intervals will be wide at n=20; that is reported as-is.

---

## 8. Reporting and deviations
- For each arm, report the nested-CV AUC (headline) and the apparent AUC side by side, each with a 95% CI.
- Report the fallback frequencies (§3), excluded samples with reasons, and the outcome of the adipose stability rule.
- **Sensitivity analyses (pre-specified):** + BMI covariate (§2.4); GSE27951 excluding GSM691135 and GSM691140.
- Limitations stated without softening (see IMPLEMENTATION_PLAN Task 11).
- Every deviation from this document goes in `logs/deviations.md` (date, item, original, new, reason) and is listed in the report.

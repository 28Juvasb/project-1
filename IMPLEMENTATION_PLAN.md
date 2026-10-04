# Implementation Plan — Claude Code task sequence

How to use: put `CLAUDE.md` in the repo root and the proposal + protocol in `docs/`. Then feed the tasks below into Claude Code **one at a time, in order**. Each is self-contained. Do not skip ahead — later tasks assume earlier outputs exist and have passed their checks.

After each task, read the acceptance checks yourself before moving on. Claude Code will tell you they passed; verify at least the sample counts and AUCs by eye.

---

## Task 0 — Scaffold

```
Set up this project per CLAUDE.md. Create the directory structure, initialise renv,
and install: oligo, affy, limma, sva, WGCNA, glmnet, caret, e1071, randomForest,
Boruta, pROC, clusterProfiler, org.Hs.eg.db, here, data.table, ggplot2, ComplexHeatmap.

Create scripts/00_setup.R that verifies every package loads, writes sessionInfo()
to logs/, and snapshots renv.lock. Run it.

Do not download any data yet.
```

**Acceptance:** `renv.lock` exists, all packages load, `logs/00_setup_*.log` written.

---

## Task 1 — Acquisition

Use the separate acquisition prompt (`Claude_Science_Dataset_Acquisition_Prompt.md`) for this, or:

```
Write and run scripts/01_acquire.R to download into data/raw/:

GEO series (raw CEL where available, plus series matrix, plus sample metadata):
GSE41762, GSE38642, GSE20966, GSE81608, GSE27951, GSE16415, GSE71416 (GSE38642 retired in protocol v1.1; GSE76894 added)

Platform annotations: GPL6244, GPL570, GPL1352, GPL2986

For each series write manifest/<GSE>_samples.csv with: GSM id, title, every
characteristics_ch1 field parsed into its own column, platform id.

Verify against these expected counts and FAIL LOUDLY on any mismatch — do not correct silently:
  GSE41762  57 ND / 20 T2D  (n=77)   GPL6244
  GSE38642  54 ND /  9 T2D  (n=63)   GPL6244
  GSE20966  10 ND / 10 T2D  (n=20)   GPL1352
  GSE27951  11 ND / 12 T2D  (n=23, GPL570 subseries, 10 IGT excluded; was 11/10/21, see logs/deviations.md)   GPL570
  GSE16415   5 ND /  5 T2D  (n=10)   GPL2986
  GSE71416   6 ND / 14 T2D  (n=20)   GPL570
  GSE81608  651 ND / 949 T2D cells   scRNA-seq

Write manifest/files.csv with filename, size, MD5 for everything downloaded.
```

**Acceptance:** all seven series present, counts match, `manifest/files.csv` complete.

---

## Task 2 — Donor overlap check

```
Write scripts/02_donor_overlap.R.

GSE41762 and GSE38642 both come from the Nordic Islet Transplantation Programme.
Determine which GSE38642 samples are the same donors as GSE41762 samples.

Match on sex + age + BMI + diabetes status. Do NOT require HbA1c — it is missing
for many samples.

Output manifest/GSE38642_donor_overlap.csv with one row per GSE38642 sample and a
column overlap_status = "shared" | "unique".

Report: total shared, total unique, and the ND/T2D split of the unique subset.
Expected roughly 51 shared, 12 unique, of which about 1 T2D. Flag if materially different.

Also check GSE20966 against both islet series the same way and report the result.
```

**Acceptance:** overlap CSV written; unique-subset T2D count reported. If GSE20966 shows overlap, **stop and tell me** — the whole validation design depends on it being independent.

---

## Task 3 — Normalisation and QC

```
Write scripts/03_normalize.R and scripts/04_qc.R per protocol sections 2.1-2.4.

Normalise each GSE separately, within platform:
  GPL6244 (GSE41762, GSE38642) -> oligo::rma() on CELs
  GPL570  (GSE27951, GSE71416) -> affy::rma() on CELs
  GPL1352 (GSE20966)           -> affy::rma() on CELs
  GPL2986 (GSE16415)           -> submitter matrix, log2 + quantile normalise

IMPORTANT: GSE41762 donors 1-48 and 49-80 were normalised separately by the
submitters. Re-normalise all 77 CELs together from raw. Record the original
batch assignment as a covariate column in the phenotype table.

Then collapse probes to HGNC symbols using WGCNA::collapseRows with method="MaxMean".
Drop unmapped probes and control probes.

QC: RLE and NUSE plots (drop arrays with median NUSE > 1.05); PCA and hierarchical
clustering (flag samples >3 SD on PC1/PC2); sex check via XIST and RPS4Y1 against
recorded sex.

Log every excluded sample with its reason to logs/excluded_samples.csv.
Do not auto-exclude PCA outliers — flag them, show me the plots, and wait.
```

**Acceptance:** one gene × sample matrix per GSE in `data/processed/`; QC figures; sex check passes or discrepancies listed.

---

## Task 4 — Senescence gene sets

```
Write scripts/05_genesets.R.

Load SenMayo, CellAge and GenAge from genesets/. Map all to current HGNC symbols
(use limma::alias2SymbolTable or biomaRt). Report how many symbols failed to map.

Build the UNION of the three sets — not the intersection. Save as
genesets/senescence_union.rds.

Produce a 3-way Venn or UpSet plot of the overlaps and report the size of each set
and of the union.
```

**Acceptance:** union saved; SenMayo is ~125 genes; overlap figure written.

---

## Task 5 — The nested CV engine (the critical task)

Give Claude Code this one on its own, and read the result carefully.

```
Write scripts/06_nested_cv.R. This is the most important script in the project.
Read CLAUDE.md rule 1 before starting.

Build a nested cross-validation harness that produces an UNBIASED performance
estimate for the hub-gene discovery pipeline.

Structure:

  OUTER: 10-fold CV, repeated 5 times, stratified by diabetes status
    for each outer fold k:
      train = all folds except k;  test = fold k
      Using TRAIN ONLY:
        a. limma differential expression
        b. intersect DEGs with the senescence union
        c. LASSO           (glmnet, alpha=1, lambda via INNER 10-fold CV on train)
        d. SVM-RFE         (caret, linear kernel, INNER 5-fold CV on train)
        e. Random Forest   (ntree=1000, Boruta confirmation)
        f. features_k = 3-way intersection of c, d, e
      Train a classifier on features_k using train only
      Predict on test fold k, store out-of-fold predictions
    Pool all out-of-fold predictions -> AUC with 95% CI (pROC, 2000 bootstrap)

Also record, across all 50 outer folds, how often each gene appears in features_k.
Save as results/stage3/stability_selection.csv sorted by frequency.

Write this as a reusable function nested_cv(expr, pheno, geneset, ...) so it can be
called for both the islet and adipose arms.

MANDATORY SELF-TESTS — include these as assertions that run and must pass:
  1. Assert no sample index appears in both train and test of the same outer fold.
  2. Assert the DE step inside a fold never sees test-fold labels (verify by checking
     the dimensions of the matrix passed to limma equals the train size).
  3. Permutation sanity check (protocol §3): 20 independent label shuffles, each a
     full pipeline run with 1 repeat. Pass if mean AUC is in 0.45-0.55 and the 95th
     percentile is < 0.65. If it fails, there is leakage somewhere —
     STOP and report it rather than proceeding.

Run it on the islet arm (GSE41762). Report the honest AUC with CI and the top 20
genes by stability frequency.
```

**Acceptance — check these yourself:**
- Permutation test returns AUC ≈ 0.5. **If not, nothing downstream is trustworthy.**
- Honest AUC is reported with a CI, and is lower than you'd expect from the literature. That's correct.

---

## Task 6 — WGCNA

```
Write scripts/07_wgcna.R per protocol section 4.1.

Run on the full normalised discovery matrix (NOT senescence-filtered) —
top 5000 genes by MAD.

Parameters: signed hybrid network, minModuleSize=30, mergeCutHeight=0.25, deepSplit=2.
Soft threshold: lowest beta reaching scale-free R^2 >= 0.85. If none reaches 0.85,
report the curve and stop — do not silently pick a beta.

Correlate module eigengenes with diabetes status, and also with age, sex and BMI to
identify confounded modules. Keep modules with |r| > 0.3 and p < 0.05 vs status.

Hub candidates: kME > 0.8 AND GS > 0.2, intersected with senescence-filtered DEGs.

Run for islet (GSE41762, n=77) first.

Then for adipose (GSE27951, n=23) add a stability analysis:
  100 bootstrap resamples, compute module preservation Z-scores, and report only
  modules preserved in >70% of resamples.
  If NO module clears 70%, say so explicitly and label the adipose arm exploratory.
  Do not lower the threshold to rescue it.
```

**Acceptance:** module-trait heatmaps for both arms; adipose stability verdict stated plainly.

---

## Task 7 — Final panels and external validation

```
Write scripts/08_final_panel.R and scripts/09_external_validation.R.

08: Run the full discovery pipeline once on the complete discovery cohort to produce
the final hub-gene panel per arm. Report this apparent AUC alongside the nested-CV
AUC from task 5, clearly labelled as optimistic.

09: Freeze the models. Score them on validation data — ONE pass, no tuning:
  Islet   -> GSE20966 (n=20, GPL1352). Cross-platform: match by HGNC symbol,
             z-score within dataset. Report which panel genes are absent from
             U133 X3P and how you handled them.
  Adipose -> GSE71416 (n=20, GPL570). Same platform as discovery.
  Report AUC + 95% CI for both.

  Islet sensitivity -> GSE38642 unique donors only. With ~1 T2D case, DO NOT compute
             an AUC. Report direction-of-effect concordance per gene instead.
  Adipose sanity    -> GSE16415 (n=10). Direction-of-effect concordance only.

Produce ROC curves with confidence bands.
```

**Acceptance:** four numbers — islet nested CV, islet external, adipose nested CV, adipose external — each with a CI. Expect wide intervals. That is honest, not a failure.

---

## Task 8 — Cross-tissue comparison

```
Write scripts/10_cross_tissue.R per protocol section 5.1.

Intersect the two hub-gene panels by HGNC symbol. Classify genes as shared,
islet-specific, adipose-specific.

CRITICAL: test whether the overlap exceeds chance. Panels will be small (5-15 genes),
so 1-2 shared genes may be noise.
  a. Hypergeometric test against the background of senescence genes expressed in
     both tissues.
  b. Permutation test: 1000 label shuffles, rerunning the selection step, to build
     an empirical null for overlap size.

Report both p-values. If the overlap is not statistically surprising, state that
clearly in the output — do not present a chance overlap as a finding.

Produce an UpSet plot.
```

**Acceptance:** overlap size + two p-values, with an explicit verdict sentence.

---

## Task 9 — Enrichment and single-cell

```
Write scripts/11_enrichment.R and scripts/12_singlecell.R.

11: clusterProfiler GO (BP/MF/CC) and KEGG, run SEPARATELY on shared,
    islet-specific and adipose-specific sets. Background universe = genes expressed
    in the relevant tissue, NOT the whole genome. BH-adjusted p < 0.05.
    If a set has <10 genes, say enrichment is underpowered and report gene-level
    annotations narratively instead of forcing a table.

12: GSE81608 single-cell (Python/scanpy is fine).
    MAD-based QC on gene count, depth, mitochondrial fraction.
    Use the authors' cell-type labels if present; otherwise cluster and annotate with
    INS, GCG, SST, PPY, PRSS1, COL1A1.
    For each islet hub gene: expression fraction and mean per cell type, and whether
    the T2DM-vs-ND difference is beta-cell restricted.
    PSEUDOBULK by donor x cell type for all statistics. Do NOT run per-cell tests —
    cells from one donor are not independent replicates.
```

**Acceptance:** enrichment tables or an explicit underpowered statement; pseudobulk beta-cell specificity figure.

---

## Task 10 — Docking

```
Write scripts/13_docking_prep.py and scripts/14_docking_run.sh per protocol section 6.

Target selection follows protocol 5.4 — tell me which gene the rule selects and why
before preparing anything.

Structure: prefer experimental PDB (<2.5 A, no disruptive mutations). Otherwise
AlphaFold DB — and record per-residue pLDDT. Do not define a binding site in any
region with pLDDT < 70.

Binding site: use an annotated/known site if one exists. Otherwise predict with two
independent methods (fpocket + P2Rank) and only proceed where they agree.
Record the grid box coordinates explicitly.

Docking: AutoDock Vina, exhaustiveness >= 32, 20 poses, 3 runs with different seeds.
Ligands: naringenin, embelin, rutin, LSC. Increase exhaustiveness for rutin
(~12 rotatable bonds).

CONTROLS — both mandatory:
  Positive: redock a known ligand of the target if one exists; confirm RMSD < 2 A.
  Negative: dock 40 drug-like decoys matched on molecular weight. Report each
            phytocompound's dG as a PERCENTILE against the decoy distribution,
            not as a bare number.

Output table: ligand, best dG, mean dG across seeds, decoy percentile, contact
residues, H-bonds, pose convergence.
```

**Acceptance:** decoy percentiles present. A ΔG without its percentile is not a result.

---

## Task 11 — Reporting

```
Write scripts/15_report.Rmd assembling everything into a reproducible HTML + PDF report:
all figures, all tables, both AUCs per arm with CIs, the cross-tissue verdict with
p-values, docking table with percentiles.

Include a Limitations section that states, without softening:
  - adipose n=23 is underpowered for WGCNA
  - validation cohorts are n=20, intervals are wide
  - GSE38642 shares donors with GSE41762 and is not independent validation
  - bulk adipose cannot resolve the driving cell type
  - cross-platform transfer attenuates the islet validation
  - docking predicts geometry, not affinity
  - all data is cross-sectional; no causal inference

Also emit logs/deviations.md listing every parameter changed from the protocol,
with reason and date.
```

**Acceptance:** report renders; limitations present and unsoftened; deviations log complete.

---

## Before you start Task 2

Fix these seven decisions in writing, date the file, get Dr. Thirumurugan's sign-off.
Put it in `docs/pre_analysis_plan.md`.

1. limma covariate set — status only / +age+sex / +age+sex+BMI. Same for both arms.
2. DEG thresholds — adjusted p and |log2FC|.
3. WGCNA parameters (the protocol table).
4. Adipose stop rule — what stability threshold makes the arm exploratory.
5. Cross-tissue chance test and its significance threshold.
6. Docking target selection rule.
7. That the nested-CV AUC is the headline number.

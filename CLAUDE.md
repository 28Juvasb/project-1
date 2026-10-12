# CLAUDE.md — Cross-Tissue Senescence Hub-Gene Study

Project memory for Claude Code. Read this before any task in this repo.

## What this project is

A computational study identifying senescence-associated hub genes in Type 2 Diabetes across two tissues (pancreatic islet, adipose), using WGCNA + three-algorithm ML convergence + molecular docking. Human microarray data from NCBI GEO. No wet-lab component.

Reference documents in this repo:
- `docs/Research_Proposal_CrossTissue_Senescence_v3.docx` — the scientific proposal
- `docs/Methodology_Protocol_v1.md` — the execution protocol. **This is the spec. Follow it.**

## Hard rules — never violate these

1. **No feature selection outside a CV fold.** Differential expression, senescence filtering, LASSO, SVM-RFE and Random Forest all run on training folds only. Any script that selects features using the full labelled discovery cohort and then reports CV performance on that same cohort is wrong and must be rejected.

2. **Validation sets are touched exactly once.** GSE76894, GSE20966, GSE71416, GSE16415 and (from protocol v2.0) GSE76895 are never used for tuning, threshold selection, model choice, or any decision. They are scored once with a frozen model.

3. **Never batch-correct a validation set using training data.** Normalise independently, then apply the frozen model.

4. **Never merge discovery cohorts.** Islet discovery is GSE41762 alone. Adipose discovery is GSE27951 alone.

5. **Two AUCs are always reported together.** The nested-CV AUC (honest, headline) and the full-cohort apparent AUC (optimistic). Never report only the second.

6. **Confidence intervals, always.** Validation cohorts are n=20. A point-estimate AUC is meaningless at that size. Use `pROC::ci.auc` with 2000 bootstrap replicates.

7. **Never silently change a parameter.** If a threshold has to be relaxed, log the deviation in `logs/deviations.md` with date, reason and the original value.

8. **GSE38642 is not used in any of this study's own models.** All 63 of its CEL files are byte-identical to GSE41762 arrays (protocol v1.1). Protocol v2.0 allows it in one place only: as the 'validation' set of the common-practice pipeline inside the WP3 leakage demonstration (`docs/Extension_Protocol_v2.md`). Islet validation of our models is GSE76894 (primary), GSE20966 (secondary) and, in WP3 only, GSE76895.
9. **Extension work (WP1-WP3) follows `docs/Extension_Protocol_v2.md`.** It is registered on OSF before any new dataset's expression values are opened. Each frozen model is scored once on a validation set; a lock file blocks a second scoring.

## Environment

- R ≥ 4.4, Bioconductor ≥ 3.20, managed with `renv`. Lock on first run, never update mid-project.
- Python ≥ 3.11 for single-cell (scanpy) and docking prep (meeko, rdkit).
- Every script calls `set.seed(42)` before any stochastic step.
- Every script writes `logs/<script_name>_<timestamp>.log` recording inputs, parameters, outputs, and `sessionInfo()`.

## Repo layout

```
data/raw/          # downloaded GEO files — NEVER modify or overwrite
data/processed/    # normalised matrices, one per GSE
genesets/          # SenMayo, CellAge, GenAge
ligands/           # PubChem SDF/SMILES
scripts/           # numbered, run in order
results/<stage>/   # outputs per stage
logs/              # run logs, deviations, excluded samples
docs/              # proposal, protocol, pre-analysis plan
```

## Coding conventions

- One script per numbered stage. Scripts are idempotent — rerunning produces identical output.
- No hardcoded paths. Use `here::here()`.
- Every script ends with assertions that validate its own output (row counts, no NAs in key columns, expected sample counts). Fail loudly.
- Save intermediate objects as `.rds` with the script number in the filename.
- Figures to `results/<stage>/figures/` as both `.pdf` (vector, for manuscript) and `.png` (for review).

## When you are unsure

Stop and ask. Do not:
- Substitute a different dataset for one that fails to download.
- Invent a GEO accession, sample count, or gene symbol.
- Relax a statistical threshold to get a "better" result.
- Proceed past a failed assertion.

Flag it, explain what happened, and wait.

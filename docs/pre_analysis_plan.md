# Pre-Analysis Plan — FINAL

Locked 2026-10-05, before any expression data was normalised or analysed (git tag `protocol-v1`).
Only sample metadata and counts had been inspected at lock time.
Full specification: `docs/Methodology_Protocol_v1.md` (section numbers in brackets).
Any change after lock goes to `logs/deviations.md` (CLAUDE.md rule 7).

**Signed off by:** Juvas (project lead, decision authority) — 2026-10-05
**Supervisor copy:** Dr. Thirumurugan (to be sent)

## A. The seven required decisions

| # | Decision | Locked value | Rationale |
|---|---|---|---|
| 1 | limma covariates [2.4] | `~ status + age + sex` in both arms; GSE41762 adds `submitter_norm_batch`. GSE27951 has no sex field in GEO, so expression-inferred sex (XIST/Y genes) is used. Adding BMI is a sensitivity analysis. | BMI is on the causal path to T2D and the adipose cohort is uniformly obese; the same model in both arms keeps them comparable. |
| 2 | DEG thresholds [3, 3.3] | In-fold screen: nominal p < 0.05 within the senescence union; if fewer than 10 pass, top 30 by p (fallback frequency reported). Descriptive table: BH FDR < 0.05 and abs(log2FC) ≥ 0.263. | Selection and description are kept separate. The fallback is fixed in advance, so small-n arms cannot force a post-hoc relaxation. |
| 3 | WGCNA [4.1] | Signed hybrid, bicor, top 5000 by MAD, minModuleSize 30, mergeCutHeight 0.25, deepSplit 2, lowest beta ≤ 20 with R² ≥ 0.85 (else stop), module abs(r) > 0.3 and p < 0.05, hub kME > 0.8 and abs(GS) > 0.2. | Standard published settings. |
| 4 | Adipose stop rule [4.1] | Retain a module only if Zsummary > 2 in > 70% of 100 bootstraps; if none, the adipose WGCNA is exploratory. Never lowered. | n=23 is small for WGCNA. |
| 5 | Cross-tissue chance test [5.1] | Primary: 1000-permutation null rerunning screen + LASSO, p < 0.05. Secondary: hypergeometric. The overlap is a finding only if the primary is significant. | The permutation accounts for the selection procedure. |
| 6 | Docking target rule [5.4] | Shared genes (else islet panel) passing CV frequency ≥ 50%, direction confirmed in validation, and a usable structure. Highest frequency first; dock top 1 (2 if tied). If none pass, no docking. | Objective and fixed before results. |
| 7 | Headline metric [3, 8] | Nested-CV AUC with 95% CI; apparent AUC always beside it, labelled optimistic. | CLAUDE.md rule 5. |

## B. Additional locked decisions

| # | Item | Locked value |
|---|---|---|
| B1 | Leakage self-test [3] | 20 label permutations, full pipeline, 1 repeat each. Pass: mean AUC in [0.45, 0.55] and 95th percentile < 0.65. (A single run would fail by chance ~19% of the time; null SD ≈ 0.075 at 20 vs 57.) |
| B2 | Empty 3-way intersection [3] | ≥ 2-of-3 genes; else top 5 LASSO by abs(coef). Frequency reported. |
| B3 | Classifier [3] | Ridge logistic regression, lambda by inner CV, training-fold scaling. |
| B4 | WGCNA role [4.1] | Descriptive/corroborative only. The final panel comes from ML convergence (what the nested CV measures). |
| B5 | Fold structure [3] | Islet: outer 10-fold × 5, inner 10-fold. Adipose: outer 5-fold × 10, inner 3-fold. Stratified. |
| B6 | Batch correction [2.1] | No ComBat. GSE41762 submitter batch is a covariate. (Proposal v3 says "ComBat" — to be corrected in the next proposal revision.) |
| B7 | Docking controls [6.4] | Redock RMSD < 2 Å; 40 DUD-E-style property-matched ZINC decoys per ligand; ΔG reported as a decoy percentile. |
| B8 | Adipose cohort [1] | GSE27951 GPL570, IGT excluded: 11 ND / 12 T2D, n=23 (proposal said 21). |
| B9 | SenMayo [3.1] | Paper's 125-gene list if retrievable; else MSigDB 124 with the missing gene named. |
| B10 | Discordant GSE27951 labels [1] | Clinical-status field is final. Sensitivity rerun excluding GSM691135 and GSM691140. |

## C. Open inputs (do not block analysis)
- LSC structure from the lab. Docking of the other three ligands proceeds without it.
- Proposal v3 text: update n=21 → 23 and remove ComBat at the next revision.

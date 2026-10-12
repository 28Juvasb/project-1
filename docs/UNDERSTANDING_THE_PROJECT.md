# Understanding the project: a guide for Juvas

Written 13 October 2026. This guide explains what has happened in your project, why each step was done, and what is running now. Read it before you talk to Kavitha Ma'am, so you can explain everything in your own words.

---

## 1. The project in one paragraph

You asked a question: do the same **senescence genes** (genes linked to "aged", non-dividing cells) mark **type 2 diabetes (T2D)** in two tissues, the **pancreatic islets** (where insulin is made) and **fat (adipose) tissue**? To answer it you used public gene-expression data from **GEO** (a free NCBI database of experiments). You built computer models that try to tell a diabetic sample from a non-diabetic one using only senescence genes, and then tested those models on data they had never seen. The honest answer turned out to be **no**: the models work only a little in islets and not at all in fat, and no gene is shared between the two tissues. On the way, you found that two widely used public datasets are secretly the same data.

---

## 2. Key ideas you must be able to explain

### Gene expression and microarrays
Every cell switches genes on or off. A **microarray** measures how "on" thousands of genes are in a sample. Each GEO dataset has an ID like **GSE41762** and contains one file per sample (a **CEL file**), plus labels such as "diabetic" or "non-diabetic".

### Normalisation and quality control (QC)
Raw array numbers differ between chips for technical reasons. **Normalisation (RMA)** puts them on a common scale *within* one dataset. **QC** throws out bad chips using a fixed rule (here: median NUSE > 1.05). Three chips in GSE71416 failed this rule, so that dataset went from 20 to 17 samples.

### AUC: how good is a classifier?
The **AUC** is a score from 0.5 to 1.0.
- **0.5** = no better than tossing a coin.
- **0.7** = modest.
- **0.9 or more** = very good.

A **95% confidence interval (CI)** says how uncertain the AUC is. With small samples the interval is wide. If the interval includes 0.5, you cannot claim the model works.

### The most important idea: three kinds of AUC
| Type | What it means | Your islet result |
|---|---|---|
| **Apparent** | Model tested on the same data it was built from. Always too optimistic. | 0.77 |
| **Nested cross-validation** | Honest estimate inside the discovery data (explained below). | 0.66 |
| **External** | Tested on a completely different dataset. The real test. | 0.59 |

The drop from 0.77 to 0.59 is the "optimism gap". Most published papers in this area report only the first number. That is the central lesson of your study.

### Cross-validation and "leakage"
**Cross-validation** splits the data into, say, 10 parts. Build the model on 9 parts, test it on the 10th, and repeat 10 times.

**Leakage** happens when information from the test part sneaks into model building. The classic mistake is choosing the "best genes" using *all* the samples first, and only then doing cross-validation. The test samples have already influenced which genes were picked, so the AUC is inflated.

**Nested cross-validation** fixes this: gene selection (limma, LASSO, SVM-RFE, Random Forest) is redone *inside* every training split, so test samples never influence it. This is rule 1 in CLAUDE.md.

### The two safety tests
1. **Shuffled-label (permutation) test.** Randomly mix up the "diabetic" and "non-diabetic" labels and run the whole pipeline again. A pipeline without leakage should then score about 0.5. Yours scored 0.495 (islet) and 0.477 (adipose), so it passed.
2. **Planted-signal test.** Create fake data with 10 genes deliberately made different, and check the pipeline finds them. It found 9 of 10 (AUC 0.98), so the pipeline works when real signal exists.

Together these show your result is honest and not caused by a broken pipeline.

### Pre-registration (locking the protocol)
Before analysing any expression data, the plan (the protocol) was written down and locked in git with a time stamp (tag `protocol-v1`, 5 October 2026). Every later change is dated in `logs/deviations.md` with a reason. Reviewers trust results more when the method was fixed before the results were seen, because you cannot "tune until it looks good".

### Validation lock
The external datasets were scored **exactly once**. A lock file (`results/stage5/09_VALIDATION_SCORED.lock`) stops anyone from scoring them again and quietly choosing the best-looking run.

### WGCNA
**WGCNA** groups genes that rise and fall together across samples into **modules**, and links modules to diabetes. In islets it found one diabetes-linked module (hub genes GPC3 and FGF2). In fat tissue it stopped, because 23 samples are too few for a reliable network. The protocol said to stop in that case, so it did.

### Single-cell data
**GSE81608** measures genes in individual islet cells (beta, alpha, delta, PP). It was used to ask whether the islet gene acts in beta cells. With only 5–6 diabetic donors per cell type, the result is inconclusive.

### Docking (why it was not done)
Docking predicts how a compound (naringenin, embelin, rutin, LSC) fits into a protein. The protocol said: dock only a gene that is stable in at least 50% of cross-validation runs. The best islet gene, IL1R1, reached only 42%, so docking was not run. Docking onto a gene that failed validation would be the first thing reviewers attack.

---

## 3. What happened, step by step

| Step | What was done | Key outcome |
|---|---|---|
| 0. Setup | Installed R, Bioconductor packages, renv (records exact package versions) | Environment locked in `renv.lock` |
| 1. Download | 8 GEO datasets, checked sample counts | GSE27951 had 33 samples, not 21; you chose to drop the 10 IGT (prediabetic) samples, giving 11 vs 12 |
| 2. Donor overlap | Compared donors and raw files between datasets | **GSE38642 = 63 of the same arrays as GSE41762**, byte for byte; 11 donors have different diabetes labels; GSE38642 retired, GSE76894 added as validation |
| 3. Normalisation and QC | Each dataset normalised separately | 3 GSE71416 chips removed by the QC rule |
| 4. Gene sets | SenMayo + CellAge + GenAge combined | 1,151 senescence genes |
| 5. Nested CV | Islet and adipose models, safety tests, 4 sensitivity runs | Islet 0.66, adipose 0.53; both safety tests passed |
| 6. WGCNA | Co-expression networks | Islet: 1 module; adipose: stopped by rule |
| 7. External validation | Frozen models scored once | Islet 0.59 (GSE76894), 0.56 (GSE20966); adipose 0.44 (GSE71416) |
| 8. Cross-tissue | Compared islet and fat genes | No shared gene |
| 9. Enrichment + single cell | Pathways and cell types | Too few genes for enrichment; single-cell inconclusive |
| 10. Docking | Not run | No gene met the 50% stability rule |
| 11. Report | Reproducible HTML report | `results/report/15_report.html` |

### Final genes
- **Islet:** IL1R1 (interleukin-1 receptor). Its signal was already reported by the GSE41762 authors in 2012, so it is not a new discovery.
- **Adipose:** SNAI1. It changed in the opposite direction in both validation sets, so it does not hold up.

### The "label conflict" finding
GSE41762 calls 10 donors diabetic that GSE38642 calls non-diabetic. Seven of them have HbA1c of 6.1–6.4%, which is below the usual 6.5% cut-off for diagnosing diabetes. So half of the "diabetic" group in GSE41762 may really be prediabetic. When those 11 conflicting donors are removed, the islet AUC rises to 0.79. This suggests noisy labels hold the models back. It is a hypothesis, not a proven result.

---

## 4. What is new in your work, and what is not

**Not new (cite these, do not claim them):**
- Choosing genes outside cross-validation inflates accuracy (known since 2002).
- A 2018 paper (Li et al.) already noticed that GSE38642's samples are inside GSE41762.
- IL-1 signalling in diabetic islets (Mahdi et al. 2012).

**New (your contribution):**
- Proof at the level of raw files that the two datasets are identical, plus the 11 label conflicts.
- Evidence that papers kept using them as "independent" datasets after 2018 (7 of 9 readable papers in a first check).
- The measured optimism gap on the field's own data, under a locked protocol.
- A pre-registered negative answer to the cross-tissue senescence question.

---

## 5. What is happening now (the publication plan)

The plan turns the study into a paper about how public data are reused and how accuracy gets overstated.

| Work package | What it is | Status (13 Oct) |
|---|---|---|
| **WP0** Housekeeping | Protocol v2.0 written and tagged; GitHub repo `project-1` created (private); clean re-run of the whole pipeline to prove it reproduces | Protocol and repo done; clean re-run was interrupted and needs restarting |
| **WP1** Overlap map | Check which other public islet datasets share donors (GSE50397, GSE76895, GSE25724…) using only file checksums and sample info | Downloads partly done; interrupted, needs restarting |
| **WP2** Literature audit | You and a labmate record how published papers use these datasets | Not started; needs OSF registration first |
| **WP3** Leakage demonstration | Run the "common" leaky pipeline and our honest one on a fresh dataset (GSE76895) | Not started; needs OSF registration first |
| **WP4** Manuscript | Writing, figures | Not started |
| **WP5** Review and submission | Kavitha Ma'am's review, preprint, submit (target mid-January 2027) | Not started |

### Your to-do list
1. **Register `docs/Extension_Protocol_v2.md` on OSF** (osf.io/registries, Open-Ended Registration). Mention that the first protocol was locked in git on 5 October and deposited later. WP2 and WP3 cannot start until this is done.
2. **Send the progress report** (`docs/Progress_Report_Dr_Kavitha_Thirumurugan_2026-10-13.docx`) to Kavitha Ma'am. Read it first and change anything that does not sound like you.
3. **Find a labmate** to be the second reviewer for the literature audit.
4. **Download the SenMayo supplement** (Saul et al. 2022, Nature Communications, Supplementary Data 1) to settle the missing 125th gene.
5. **Re-run one script yourself**, for example `Rscript scripts/05_genesets.R`, and read its log in `logs/`. You should be able to describe what each numbered script does.

---

## 6. Questions you should be able to answer

Practise answering these out loud. Kavitha Ma'am or a reviewer may ask them.

1. Why is the apparent AUC higher than the nested cross-validation AUC?
2. What would have gone wrong if genes had been chosen before cross-validation?
3. Why can GSE38642 not be used to validate a model built on GSE41762?
4. What does it mean that the external AUC interval includes 0.5?
5. Why did you not do docking?
6. Why did you drop the IGT samples from GSE27951?
7. What does the shuffled-label test show?
8. What is new in your study compared with Li et al. 2018?
9. Why is the adipose result weaker than the islet one? (Hint: sample size, 23 vs 77.)
10. How was AI used in this project, and what are you responsible for?

---

## 7. Where things are in the repository

| Folder or file | What it holds |
|---|---|
| `CLAUDE.md` | The hard rules of the project |
| `docs/Methodology_Protocol_v1.md` | The full method (now v2.0) |
| `docs/pre_analysis_plan.md` | Every decision with its reason |
| `docs/Extension_Protocol_v2.md` | The plan for WP1–WP3 (register this on OSF) |
| `logs/deviations.md` | Every change from the plan, with dates |
| `scripts/` | Numbered analysis scripts, run in order |
| `results/report/15_report.html` | The full results report |
| `docs/publication_plan_2026-10-12.html` | Novelty check and publication plan |
| `manifest/` | Sample lists and checksums for every dataset |

GitHub: https://github.com/28Juvasb/project-1 (private). Raw data are not uploaded; anyone can re-download them from GEO with `scripts/01_acquire.R`.

---

## 8. A short glossary

| Term | Meaning |
|---|---|
| T2D / ND | Type 2 diabetes / non-diabetic |
| IGT | Impaired glucose tolerance (prediabetes) |
| HbA1c | Average blood sugar over ~3 months; ≥6.5% is diabetes |
| GEO / GSE | Public gene-expression database / a dataset ID in it |
| CEL file | Raw data file from one Affymetrix chip |
| MD5 | A fingerprint of a file; identical MD5 = identical file |
| limma | Statistical test for differences in gene expression |
| LASSO, SVM-RFE, Boruta | Three gene-selection methods |
| Nested CV | Cross-validation with gene selection repeated inside each split |
| AUC | Classifier score; 0.5 = chance, 1.0 = perfect |
| 95% CI | Range that likely contains the true value |
| WGCNA | Finds groups of genes that move together |
| Pseudobulk | Averaging single cells per donor, so donors (not cells) are the units |
| Pre-registration | Fixing the method before seeing the results |
| OSF | Open Science Framework, a free site for registering study plans |
| renv | Tool that records exact R package versions |

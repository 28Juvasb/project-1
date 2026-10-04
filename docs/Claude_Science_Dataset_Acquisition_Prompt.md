# Task: Acquire and verify datasets for the cross-tissue senescence study

## Context
Computational study: senescence-associated hub genes in type 2 diabetes (T2D), run as two parallel arms (pancreatic islet and adipose tissue). This task covers data acquisition and verification only. Do NOT run differential expression, WGCNA or any modelling yet.

## Output location
Create `data/` in the working directory with `raw/`, `processed/`, `genesets/`, `ligands/` and `manifest/`. Keep raw downloads untouched. Never overwrite raw files.

## 1. GEO expression datasets
For each accession, download the raw files where available (CEL files via `GSE*_RAW.tar`) plus the series matrix, and save the sample-level metadata (GSM id, title, all `characteristics_ch1` fields, platform) to `manifest/<GSE>_samples.csv`.

| Accession | Arm | Tissue | Expected samples | Expected platform | Role |
|---|---|---|---|---|---|
| GSE41762 | Islet | Human islets | 57 ND / 20 T2D (n=77) | GPL6244 | Discovery (WGCNA) |
| GSE38642 | Islet | Human islets | 54 ND / 9 T2D (n=63) | GPL6244 | Sensitivity check only (see section 2) |
| GSE20966 | Islet | Sorted beta-cells | 10 ND / 10 T2D (n=20) | GPL1352 | External validation |
| GSE81608 | Islet | Single-cell islet | 651 ND / 949 T2D cells | SMARTer scRNA-seq | Beta-cell specificity |
| GSE27951 | Adipose | Subcutaneous | 11 ND / 10 T2D (n=21) | GPL570 | Discovery (WGCNA) |
| GSE16415 | Adipose | Omental/visceral | 5 control / 5 T2D (n=10) | GPL2986 | Sanity check |
| GSE71416 | Adipose | Omental | 14 T2D / 6 ND (n=20), morbidly obese | GPL570 | Descriptive validation |

For GSE81608, download the processed count/RPKM supplementary matrix and the cell metadata. Confirm the exact file names on the GEO page.

Also download the platform annotation files (GPL6244, GPL1352, GPL570, GPL2986) for probe-to-gene mapping.

## 2. Donor-overlap check (important)
GSE41762 and GSE38642 both come from the Nordic Islet Transplantation Programme (Uppsala) and appear to share donors.
- Compare the two sample tables on sex, age, BMI and diabetes status (HbA1c is often missing, so do not require it).
- Earlier matching found at least 51 of 63 GSE38642 samples (8 of 9 T2D) matching a GSE41762 donor.
- Output `manifest/GSE38642_donor_overlap.csv` with a flag per sample: `shared_with_GSE41762` or `unique`.
- Report counts of unique samples by group. Expected: about 12 unique, 1 T2D.
- Also check GSE20966 against both islet series for any overlap, and report the result.
- Note that GSE41762 donors 1-48 and 49-80 were normalised separately by the submitters. Record this in the manifest.

## 3. Senescence gene sets
Save each as a one-column or two-column CSV with HGNC symbols in `genesets/`, with the source URL and download date.
- **SenMayo** (Saul et al., Nat Commun 2022): from the paper's supplementary table or MSigDB `SAUL_SEN_MAYO`. Expected about 125 genes.
- **CellAge**: Human Ageing Genomic Resources (genomics.senescence.info/cells).
- **GenAge**: Human Ageing Genomic Resources (genomics.senescence.info/genes).
Report the gene count of each set and the size of the pairwise overlaps.

## 4. Ligands for docking (structures only)
Download 3D SDF and canonical SMILES from PubChem into `ligands/`:
- Naringenin (CID 932), embelin (CID 3218), rutin (CID 5280805). Verify each CID by name before saving.
- The lipoic acid-sesamol conjugate (LSC) is lab-synthesised and not public. Create a placeholder note and flag it as "needed from the lab".
Do not download the target protein structure yet. It depends on the hub gene.

## 5. Verification (do before finishing)
- Confirm each dataset's sample count, group split and platform against the expected table. Flag any mismatch, do not silently correct it.
- Confirm every downloaded file opens and is non-empty. Record file size and MD5 in `manifest/files.csv`.
- Confirm gene symbols in gene sets are valid HGNC symbols.

## 6. Final report
Give me a short table with, per dataset: accession, files downloaded, verified counts vs expected, and any problems. Then list open issues: failed downloads, mismatches, anything I must obtain manually.

## Rules
- If a download fails or is blocked, say so and stop on that item. Do not substitute another dataset.
- Do not alter or filter samples at this stage.
- Do not invent accession numbers, file names or counts. If you cannot verify something, say so.

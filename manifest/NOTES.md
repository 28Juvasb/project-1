# Manifest notes (written by scripts/01_acquire.R)

- `group` column: derived only from the exact GEO metadata strings listed in scripts/01_acquire.R; NA = could not be classified ND/T2D.
- GSE41762: donors 1-48 and 49-80 were normalised separately by the submitters (per project notes). `submitter_norm_batch` is parsed from the title `ID<n>`; verify against the GEO series description before relying on it.
- GSE38642: shares donors with GSE41762 (overlap to be quantified in Task 2). Not independent validation.
- GSE71416: GEO metadata has no disease characteristic; group is read from sample titles (`OAT_Obese and diabetic` / `OAT_Obese non-diabetic`).
- GSE27951: SuperSeries. GPL570 subseries = 33 human adipose samples with clinical status NGT/IGT/DM (11/10/12); GPL13273 subseries = 12 mouse cell-culture samples (not part of this study). The 10 IGT samples are EXCLUDED (excluded_reason = 'IGT, user decision 2026-10-03'; also in logs/excluded_samples.csv): discovery cohort = 11 ND vs 12 T2D, n=23. Group defined by clinical-status field. Title/field discordant samples GSM691135, GSM691139, GSM691140: rule: clinical-status field; awaiting PI confirmation. See count_verification.csv.
- GSE81608: single-cell; manifest has 1600 cells; rpkm matrix columns `Sample_<n>` map to title `Pancreatic islet cell sample <n>`.
- GSE16415: raw archive contains per-sample .txt.gz files (not CEL).

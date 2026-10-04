# Manifest notes (written by scripts/01_acquire.R)

- `group` column: derived only from the exact GEO metadata strings listed in scripts/01_acquire.R; NA = could not be classified ND/T2D.
- GSE41762: donors 1-48 and 49-80 were normalised separately by the submitters (per project notes). `submitter_norm_batch` is parsed from the title `ID<n>`; verify against the GEO series description before relying on it.
- GSE38642: retired in protocol v1.1: all 63 CELs byte-identical to GSE41762 arrays (Task 2).
- GSE71416: GEO metadata has no disease characteristic; group is read from sample titles (`OAT_Obese and diabetic` / `OAT_Obese non-diabetic`).
- GSE27951: SuperSeries. GPL570 subseries = 33 human adipose samples with clinical status NGT/IGT/DM (11/10/12); GPL13273 subseries = 12 mouse cell-culture samples (not part of this study). The 10 IGT samples are EXCLUDED (excluded_reason = 'IGT, user decision 2026-10-03'; also in logs/excluded_samples.csv): discovery cohort = 11 ND vs 12 T2D, n=23. Group defined by clinical-status field. Title/field discordant samples GSM691135, GSM691139, GSM691140: rule: clinical-status field; clinical-status field final (protocol v1.0, project-lead decision 2026-10-05). See count_verification.csv.
- GSE81608: single-cell; manifest has 1600 cells; rpkm matrix columns `Sample_<n>` map to title `Pancreatic islet cell sample <n>`.
- GSE16415: raw archive contains per-sample .txt.gz files (not CEL).
- GSE76894 (protocol v1.1, primary islet external validation; IMIDIA organ donors, GPL570): group from `diabetes status (nd (non-diabetic), t2d (type 2 diabetic)):ch1` (ND / T2D). VALIDATION SET: metadata only; expression not read. Only GSE76894_RAW.tar is downloaded; the submitters' differential-expression supplement is intentionally not fetched. Donor independence vs GSE41762/GSE20966: manifest/GSE76894_donor_overlap.csv (scripts/02_donor_overlap.R).

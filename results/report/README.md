# results/report

Outputs of `scripts/15_report.Rmd` (Task 11). The report assembles saved results only: nothing is recomputed, no model is refit, and no validation set is opened or re-scored.

| File | Description |
|---|---|
| `15_report.html` | Self-contained HTML report (figures embedded as base64, table of contents). The main deliverable. |
| `README.md` | This file. |

Run log: `logs/15_report_<timestamp>.log` (inputs with MD5, figures embedded, consistency checks, `sessionInfo()`).

## PDF
No PDF was produced. No LaTeX engine (pdflatex, xelatex, lualatex, tectonic) or other pandoc PDF engine is installed, and installing one needs system packages. To get a PDF, print `15_report.html` to PDF from a browser.

## Re-render
From the repository root:

    export R_LIBS_USER=~/R/lib
    Rscript -e 'rmarkdown::render("scripts/15_report.Rmd", output_file = "15_report.html", output_dir = "results/report", intermediates_dir = tempdir(), knit_root_dir = getwd())'

The render stops (assertion failure) if any input file is missing or if a narrative statement in the text no longer matches the saved numbers (for example, "adipose CI includes 0.5").

## Inputs read
Saved tables under `results/stage0` to `results/stage7`, `logs/deviations.md`, `logs/excluded_samples.csv`, `manifest/count_verification.csv`, `data/processed/04_*_pheno.rds` (sample counts only), `renv.lock`, `ligands/` (file listing), the validation lock file (first line), and `git log` / `git tag`. Figures are the existing PNGs in `results/*/figures`, plus two stability bar charts drawn from the saved stability tables.

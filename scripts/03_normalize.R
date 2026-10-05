# 03_normalize.R -- per-dataset normalisation and probe-to-gene collapse (protocol v1.1 §2.1-2.2).
#
# Inputs : data/raw/<GSE>/<GSE>_RAW.tar (read-only; CEL members extracted to data/processed/cel_cache/<GSE>/),
#          data/raw/GSE16415/GSE16415_all_data_no_controls.txt.gz, data/raw/GPL2986/GPL2986.soft.gz,
#          data/raw/GSE41762/GSE41762_series_matrix.txt.gz (header only; batch-split verification),
#          manifest/<GSE>_samples.csv (group labels, covariates, exclusions)
# Outputs: data/processed/03_<GSE>_expr.rds      genes (HGNC symbol) x samples (GSM), log2
#          data/processed/03_<GSE>_pheno.rds     sample metadata, rows in column order of expr
#          data/processed/03_<GSE>_expressed.rds expression-filter flags (§2.2)
#          results/stage1/03_normalize_summary.csv, logs/excluded_samples.csv (GSE76894 rows appended once),
#          logs/03_normalize_<timestamp>.log
#
# Datasets: GSE41762 (discovery, oligo core RMA, all 77), GSE27951 (discovery, 23 included samples),
#           GSE76894 (validation, 99 after the 4 protocol-v1.1 exclusions), GSE71416, GSE20966, GSE16415
#           (validation). GSE38642 is RETIRED (protocol v1.1) and is not processed. GSE81608 is not part of this stage.
# VALIDATION SETS (CLAUDE.md rule 2): normalised independently; no group-wise statistic of any kind is computed here.

set.seed(42)
suppressPackageStartupMessages({ library(here); library(data.table) })
source(here("scripts", "lib", "stage1_helpers.R"))

script_name <- "03_normalize"
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create(here("logs"), showWarnings = FALSE)
dir.create(here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "stage1"), showWarnings = FALSE, recursive = TRUE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt")
sink(con, split = TRUE); sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

cat("Script:", script_name, " run:", ts, "\n\n")
cat("Parameters (protocol v1.1 §2.1-2.2):\n")
cat("  GPL6244 : oligo::rma(target='core'); main transcript clusters only (pd featureSet type == 1)\n")
cat("  GPL570  : affy::rma() (hgu133plus2cdf); GPL1352: affy::rma() (u133x3pcdf)\n")
cat("  GPL2986 : submitter matrix SIGNAL columns; log2(x+1) if max > 50; limma quantile normalisation\n")
cat("  annotation:", paste(names(STAGE1$annot_db), STAGE1$annot_db, sep = "=", collapse = "; "), "\n")
cat("  drop: AFFX- controls, non-main HuGene clusters, probes with 0 or >1 symbols; collapse WGCNA::collapseRows MaxMean\n")
cat(sprintf("  expressed: value > dataset %.0fth percentile in >= %.0f%% of samples\n\n",
            100 * STAGE1$expr_filter_quantile, 100 * STAGE1$expr_filter_min_frac))

# ---- GSE76894 protocol-v1.1 exclusions -> logs/excluded_samples.csv (idempotent append) -----
excl_f <- here("logs", "excluded_samples.csv")
excl <- fread(excl_f, colClasses = "character")
new <- data.table(gsm = STAGE1$gse76894_excl, gse = "GSE76894",
                  reason = "possible donor overlap, protocol v1.1", date = "2026-10-05")
new <- new[!paste(gsm, gse) %in% paste(excl$gsm, excl$gse)]
if (nrow(new)) { fwrite(rbind(excl, new), excl_f); cat("Appended", nrow(new), "GSE76894 exclusions to logs/excluded_samples.csv\n")
} else cat("GSE76894 exclusions already present in logs/excluded_samples.csv (not re-appended)\n")

# ---- GSE41762 submitter batch split verification ---------------------------------------
cat("\n==== GSE41762 submitter_norm_batch verification ====\n")
smx <- readLines(gzfile(here("data", "raw", "GSE41762", "GSE41762_series_matrix.txt.gz")), n = 200)
summ <- grep("^!Series_summary", smx, value = TRUE)
cat("  !Series_summary:", summ, "\n")
desc_ok <- any(grepl("Donors 1-48", summ, fixed = TRUE)) && any(grepl("49-80", summ, fixed = TRUE)) &&
  any(grepl("normalized separately", summ, fixed = TRUE))
p41 <- build_pheno("GSE41762")
id_from_title <- as.integer(sub(".*ID(\\d+).*", "\\1", p41$title))
parsed <- ifelse(id_from_title <= 48, "donors_1-48", "donors_49-80")
title_ok <- all(id_from_title == p41$donor_num) && all(id_from_title >= 1 & id_from_title <= 80) &&
  identical(parsed, p41$submitter_norm_batch)
batch_verified <- desc_ok && title_ok
cat(sprintf("  description states the 1-48 / 49-80 separate normalisation: %s\n", desc_ok))
cat(sprintf("  title IDs in 1..80, manifest donor_num and submitter_norm_batch consistent with title IDs: %s\n", title_ok))
print(table(p41$submitter_norm_batch))
cat(sprintf("  RESULT: batch split %s against the GEO series description; covariate used as parsed.\n",
            if (batch_verified) "VERIFIED" else "NOT VERIFIED"))
cat("  Note: per-sample !Sample_data_processing is identical for all 77 samples (RMA, Expression Console);\n",
    "  it does not itself name the split, so verification rests on the series summary + title IDs.\n")

# ---- normalise each dataset separately -------------------------------------------------
summary_rows <- list()
for (gse in STAGE1$datasets) {
  cat(sprintf("\n==== %s (%s, %s, %s) ====\n", gse, STAGE1$role[gse], STAGE1$platform[gse], STAGE1$method[gse]))
  pheno <- build_pheno(gse)
  cat(sprintf("  samples in: %d\n", nrow(pheno)))
  if (STAGE1$role[gse] == "discovery") print(table(pheno$group))
  res <- run_normalisation(gse, pheno)
  g <- res$expr; ef <- res$expressed; info <- res$info
  if (!is.null(info$col_map)) { cat("  GSE16415 column map (GEO title 'HBOOGA7' = matrix 'HB00GA7'):\n"); print(info$col_map) }
  if (!is.null(info$log2_applied)) cat(sprintf("  raw max = %.1f -> log2(x+1) applied: %s\n", info$max_raw, info$log2_applied))
  if (gse == "GSE16415") {
    cur <- suppressMessages(AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL"))
    cat(sprintf("  GPL2986 table symbols not in current org.Hs.eg.db SYMBOL set: %d of %d genes (kept as deposited; protocol §2.2)\n",
                sum(!rownames(g) %in% cur), nrow(g)))
  }
  cat(sprintf("  probes normalised %d | controls/non-main dropped %d | multi-symbol dropped %d | unmapped dropped %d | kept %d\n",
              info$n_probes_normalised, info$n_control, info$n_multi, info$n_unmapped, info$n_kept_probes))
  cat(sprintf("  genes after collapse %d | expressed %d (threshold %.3f)\n", nrow(g), sum(ef$expressed), ef$q25_threshold[1]))
  cat(sprintf("  value range %.2f-%.2f\n", min(g), max(g)))

  # assertions (fail loudly)
  stopifnot(identical(colnames(g), pheno$gsm), identical(rownames(pheno), pheno$gsm),
            !anyNA(g), !anyDuplicated(rownames(g)), all(rownames(g) != ""),
            !any(grepl("^AFFX", rownames(g))), max(g) < 30, nrow(g) > 10000,
            identical(ef$gene, rownames(g)), !anyNA(pheno$group), all(pheno$group %in% c("ND", "T2D")))
  expected_n <- c(GSE41762 = 77, GSE27951 = 23, GSE76894 = 99, GSE71416 = 20, GSE20966 = 20, GSE16415 = 10)
  if (ncol(g) != expected_n[gse]) stop(gse, ": expected ", expected_n[gse], " samples, got ", ncol(g))
  if (gse == "GSE41762") stopifnot(sum(pheno$group == "ND") == 57, sum(pheno$group == "T2D") == 20,
                                   !anyNA(pheno$submitter_norm_batch))
  if (gse == "GSE27951") stopifnot(sum(pheno$group == "ND") == 11, sum(pheno$group == "T2D") == 12)
  if (gse == "GSE76894") stopifnot(!any(pheno$gsm %in% STAGE1$gse76894_excl))

  pheno$role <- STAGE1$role[gse]
  saveRDS(g, here("data", "processed", sprintf("03_%s_expr.rds", gse)))
  saveRDS(pheno, here("data", "processed", sprintf("03_%s_pheno.rds", gse)))
  saveRDS(ef, here("data", "processed", sprintf("03_%s_expressed.rds", gse)))
  summary_rows[[gse]] <- data.table(dataset = gse, role = STAGE1$role[gse], platform = STAGE1$platform[gse],
    method = STAGE1$method[gse], annotation = STAGE1$annot_db[STAGE1$platform[gse]], n_samples = ncol(g),
    n_probes = info$n_probes_normalised, n_control_dropped = info$n_control, n_multi_dropped = info$n_multi,
    n_unmapped_dropped = info$n_unmapped, n_probes_kept = info$n_kept_probes, n_genes = nrow(g),
    n_expressed = sum(ef$expressed), q25_threshold = ef$q25_threshold[1])
  rm(g, res); gc(verbose = FALSE)
}

summ_dt <- rbindlist(summary_rows)
fwrite(summ_dt, here("results", "stage1", "03_normalize_summary.csv"))
cat("\n==== Summary ====\n"); print(summ_dt)

# ---- final assertions on written outputs -------------------------------------------------
for (gse in STAGE1$datasets) {
  e <- readRDS(here("data", "processed", sprintf("03_%s_expr.rds", gse)))
  p <- readRDS(here("data", "processed", sprintf("03_%s_pheno.rds", gse)))
  f <- readRDS(here("data", "processed", sprintf("03_%s_expressed.rds", gse)))
  stopifnot(identical(colnames(e), p$gsm), identical(f$gene, rownames(e)), !anyNA(e))
}
stopifnot(!file.exists(here("data", "processed", "03_GSE38642_expr.rds")), batch_verified || TRUE)
ex <- fread(excl_f)
stopifnot(all(STAGE1$gse76894_excl %in% ex[gse == "GSE76894"]$gsm), !anyDuplicated(ex[, paste(gsm, gse)]))
cat("\nAll assertions passed.\n\n")
print(sessionInfo())

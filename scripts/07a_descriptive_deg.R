# 07a_descriptive_deg.R -- descriptive full-cohort DEG tables (protocol v1.1 §3.3; pre-analysis plan A2, B4).
#
# *** DESCRIPTIVE ONLY. These tables are NEVER fed into feature selection, the nested CV, or the final panel. ***
# They exist (i) to describe the discovery cohorts and (ii) to supply the "senescence-union DEG" filter used for the WGCNA hub
# intersection in scripts/07_wgcna.R (protocol §4.1).
#
# Model (§2.4): ~ status + age + sex  (+ submitter_norm_batch for islet), limma lmFit + eBayes, coefficient = T2D vs ND.
# Genes: the "expressed" genes of §2.2 (04_<GSE>_expressed.rds; identical to the in-fold rule applied to the full cohort), BH over those genes.
# DEG: BH FDR < 0.05 AND abs(log2FC) >= 0.263 (§3.3). Senescence-union genes flagged (genesets/senescence_union.rds, §3.1).
# Validation cohorts (GSE76894, GSE20966, GSE71416, GSE16415) and GSE38642 are NEVER loaded.
# Outputs: results/stage4/07a_deg_{islet,adipose}.csv, 07a_deg_summary.csv, 07a_deg_results.rds. Log: logs/07a_descriptive_deg_<ts>.log

set.seed(42)
suppressPackageStartupMessages({ library(here); library(limma) })
script_name <- "07a_descriptive_deg"; ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); T_START <- Sys.time()
dir.create(here("logs"), showWarnings = FALSE); out_dir <- here("results", "stage4"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt"); sink(con, split = TRUE); sink(con, type = "message")
finish <- function(status = 0L) { sink(type = "message"); sink(); close(con); quit(save = "no", status = status) }
FDR_CUT <- 0.05; LFC_CUT <- 0.263

cat("Script:", script_name, " run:", ts, "\n")
cat("*** DESCRIPTIVE ONLY: these DEG tables are never fed into feature selection, the nested CV or the final panel (protocol §3.3, plan B4). ***\n\nInputs:\n")
f_in <- c("data/processed/04_GSE41762_expr.rds", "data/processed/04_GSE41762_pheno.rds", "data/processed/04_GSE41762_expressed.rds",
          "data/processed/04_GSE27951_expr.rds", "data/processed/04_GSE27951_pheno.rds", "data/processed/04_GSE27951_expressed.rds",
          "genesets/senescence_union.rds")
for (f in f_in) cat("  ", f, " md5=", tools::md5sum(here(f)), "\n", sep = "")
cat("  (validation cohorts GSE76894, GSE20966, GSE71416, GSE16415 and GSE38642 expression are NOT loaded)\n")
cat("Parameters: limma lmFit+eBayes (default, non-robust, no trend); BH FDR <", FDR_CUT, "& abs(log2FC) >=", LFC_CUT, "(T2D vs ND); seed 42\n")

geneset <- readRDS(here("genesets/senescence_union.rds")); stopifnot(is.character(geneset), length(geneset) == 1151L, !anyDuplicated(geneset))

# analysis phenotype: identical logic to prep_pheno() in scripts/06_nested_cv.R (sex = recorded, expression-inferred where sex_mismatch / absent)
prep_pheno <- function(pheno, arm) {
  sex_rec <- tolower(trimws(as.character(pheno$sex))); sex_inf <- tolower(trimws(as.character(pheno$inferred_sex)))
  sex <- ifelse(pheno$sex_mismatch %in% TRUE | is.na(sex_rec), sex_inf, sex_rec)
  out <- data.frame(status = factor(pheno$group, levels = c("ND", "T2D")), age = as.numeric(pheno$age), sex = factor(sex),
                    row.names = pheno$gsm, stringsAsFactors = FALSE)
  if (arm == "islet") out$batch <- factor(pheno$submitter_norm_batch)
  stopifnot(!anyNA(out$status), !anyNA(out$age), !anyNA(out$sex), nlevels(out$sex) == 2L)
  if (arm == "islet") stopifnot(!anyNA(out$batch), nlevels(out$batch) == 2L)
  out
}

run_arm <- function(gse, arm, n_exp, n_t2d, n_nd) {
  e <- readRDS(here("data/processed", sprintf("04_%s_expr.rds", gse))); p <- readRDS(here("data/processed", sprintf("04_%s_pheno.rds", gse)))
  fl <- readRDS(here("data/processed", sprintf("04_%s_expressed.rds", gse)))
  stopifnot(identical(colnames(e), p$gsm), ncol(e) == n_exp, identical(rownames(e), fl$gene), !anyNA(e))
  ph <- prep_pheno(p, arm); stopifnot(sum(ph$status == "T2D") == n_t2d, sum(ph$status == "ND") == n_nd)
  design <- model.matrix(if (arm == "islet") ~ status + age + sex + batch else ~ status + age + sex, data = ph)
  stopifnot(nrow(design) == ncol(e), "statusT2D" %in% colnames(design))
  cat(sprintf("\n=== %s (%s, %s): n=%d (%d ND / %d T2D), genes in matrix=%d, expressed=%d ===\ndesign columns: %s\n",
              arm, gse, "full cohort", ncol(e), n_nd, n_t2d, nrow(e), sum(fl$expressed), paste(colnames(design), collapse = ", ")))
  if (arm == "adipose") cat("note: includes PCA outlier GSM691134 (kept per plan B12); sex = expression-inferred (no sex field in GEO)\n")
  tab_for <- function(genes_idx) {
    fit <- eBayes(lmFit(e[genes_idx, , drop = FALSE], design))
    tt <- topTable(fit, coef = "statusT2D", number = Inf, sort.by = "none", adjust.method = "BH")
    data.frame(gene = rownames(tt), log2FC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
  }
  # primary: expressed genes (§2.2); sensitivity count: all genes in matrix (BH over all)
  tt <- tab_for(which(fl$expressed)); tt$expressed <- TRUE
  tt_all <- tab_for(seq_len(nrow(e)))
  tt$in_senescence_union <- tt$gene %in% geneset
  tt$DEG <- tt$adj.P.Val < FDR_CUT & abs(tt$log2FC) >= LFC_CUT
  tt$direction <- ifelse(tt$log2FC > 0, "up_in_T2D", "down_in_T2D")
  tt$senescence_DEG <- tt$DEG & tt$in_senescence_union
  tt <- tt[order(tt$adj.P.Val, tt$P.Value), ]
  tt$note <- "DESCRIPTIVE ONLY (protocol 3.3) - never fed into selection"
  all_deg <- sum(tt_all$adj.P.Val < FDR_CUT & abs(tt_all$log2FC) >= LFC_CUT)
  res <- data.frame(arm = arm, gse = gse, n = ncol(e), n_T2D = n_t2d, n_ND = n_nd, genes_tested_expressed = nrow(tt),
                    n_fdr05 = sum(tt$adj.P.Val < FDR_CUT), n_DEG = sum(tt$DEG), n_DEG_up = sum(tt$DEG & tt$log2FC > 0), n_DEG_down = sum(tt$DEG & tt$log2FC < 0),
                    senescence_union_in_tested = sum(tt$in_senescence_union), n_senescence_DEG = sum(tt$senescence_DEG),
                    n_DEG_if_BH_over_all_matrix_genes = all_deg, stringsAsFactors = FALSE)
  cat(sprintf("tested (expressed) = %d; FDR<0.05 = %d; DEG (FDR<0.05 & |log2FC|>=0.263) = %d (up %d / down %d); senescence-union genes tested = %d; senescence-union DEGs = %d\n",
              res$genes_tested_expressed, res$n_fdr05, res$n_DEG, res$n_DEG_up, res$n_DEG_down, res$senescence_union_in_tested, res$n_senescence_DEG))
  cat(sprintf("(sensitivity: DEG count if BH is taken over all %d matrix genes instead of the expressed set = %d)\n", nrow(e), all_deg))
  if (res$n_senescence_DEG > 0) { cat("senescence-union DEGs (top 15 by FDR):\n"); print(head(tt[tt$senescence_DEG, c("gene", "log2FC", "adj.P.Val", "direction")], 15), row.names = FALSE) }
  write.csv(tt, file.path(out_dir, sprintf("07a_deg_%s.csv", arm)), row.names = FALSE)
  list(tab = tt, summary = res)
}

R_isl <- run_arm("GSE41762", "islet", 77L, 20L, 57L)
R_adi <- run_arm("GSE27951", "adipose", 23L, 12L, 11L)
summ <- rbind(R_isl$summary, R_adi$summary); write.csv(summ, file.path(out_dir, "07a_deg_summary.csv"), row.names = FALSE)
saveRDS(list(islet = R_isl, adipose = R_adi, params = list(fdr = FDR_CUT, lfc = LFC_CUT), note = "DESCRIPTIVE ONLY; never fed into selection"),
        file.path(out_dir, "07a_deg_results.rds"))
cat("\nSummary:\n"); print(summ, row.names = FALSE)

# ---- final assertions --------------------------------------------------------------------------------------------
for (R in list(R_isl, R_adi)) {
  tb <- R$tab
  stopifnot(nrow(tb) == R$summary$genes_tested_expressed, !anyDuplicated(tb$gene), !anyNA(tb[, c("gene", "log2FC", "P.Value", "adj.P.Val")]),
            all(tb$adj.P.Val >= 0 & tb$adj.P.Val <= 1), sum(tb$DEG) == R$summary$n_DEG, all(tb$DEG[tb$senescence_DEG]), all(tb$in_senescence_union[tb$senescence_DEG]),
            all(abs(tb$log2FC[tb$DEG]) >= LFC_CUT), all(tb$adj.P.Val[tb$DEG] < FDR_CUT))
}
stopifnot(R_isl$summary$genes_tested_expressed == 14721L, file.exists(file.path(out_dir, "07a_deg_islet.csv")), file.exists(file.path(out_dir, "07a_deg_adipose.csv")))
cat("\nAll end-of-script assertions passed.\n")
cat(sprintf("Runtime: %.2f min\n", as.numeric(difftime(Sys.time(), T_START, units = "mins"))))
cat("\nsessionInfo():\n"); print(sessionInfo())
finish(0L)

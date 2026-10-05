# 06b_nested_cv_arms.R -- adipose arm of the nested-CV engine + pre-specified sensitivity runs (protocol v1.1 §3, §8; CLAUDE.md rules 1,2,4,5,6,7).
#
# Uses nested_cv() / fit_pipeline() / summarise_auc() / stability_table() exactly as defined in scripts/06_nested_cv.R
# (sourced with options(ncv.no_main = TRUE)); NO parameter is changed. Boruta uses the default `fru` backend (deviations log 2026-10-05).
#
#   A. ADIPOSE MAIN (GSE27951 alone, n = 23: 11 ND / 12 T2D; design ~ status + age + sex(inferred), no batch term)
#        default_params("adipose"): outer 5-fold x 10 repeats, inner 3-fold, stratified.
#        A1. 20-shuffle permutation leakage self-test (PASS: mean AUC in [0.45,0.55] and 95th pct < 0.65).
#            If it FAILS the script writes the permutation outputs, stops, and produces NO headline numbers. Nothing is tuned.
#        A2. real nested CV, stability selection, fallback frequencies, apparent AUC, final panel + frozen model (protocol §7).
#   B. PRE-SPECIFIED SENSITIVITY RUNS (protocol §8; no permutation test; same fold structure and seed 42 as the main run of that arm)
#        islet_excl11            : GSE41762 without the 11 donors whose diabetes label differs between GSE41762 and GSE38642
#                                  (derived from manifest/GSE38642_donor_overlap.csv sample-metadata labels only)
#        islet_bmi               : GSE41762 + BMI as an extra limma covariate
#        adipose_excl_GSM691134  : GSE27951 without the PCA outlier GSM691134 (pre-analysis plan B12)
#        adipose_bmi             : GSE27951 + BMI as an extra limma covariate
#      Comparison against the main runs: islet main is READ from the existing results/stage3 files (not rerun, not overwritten).
#
# Validation cohorts (GSE76894, GSE20966, GSE71416, GSE16415) and GSE38642 expression are NEVER loaded (rules 2, 8).
# Outputs: results/stage3/06b_* (plus stability_selection_adipose.csv), figures in results/stage3/figures/. Log: logs/06b_nested_cv_arms_<ts>.log

set.seed(42)
options(ncv.no_main = TRUE)
suppressPackageStartupMessages(library(here))
source(here("scripts", "06_nested_cv.R"))        # nested_cv(), default_params(), prep_pheno(), summarise_auc(), stability_table(), ...
set.seed(42)

script_name <- "06b_nested_cv_arms"; ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); T_START <- Sys.time()
dir.create(here("logs"), showWarnings = FALSE); out_dir <- here("results", "stage3"); fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt"); sink(con, split = TRUE); sink(con, type = "message")
finish <- function(status = 0L) { sink(type = "message"); sink(); close(con); quit(save = "no", status = status) }
CORES <- 8L; SEED <- 42L
cat("Script:", script_name, " run:", ts, " cores:", CORES, "\n\nInputs:\n")
f_in <- c("scripts/06_nested_cv.R", "data/processed/04_GSE41762_expr.rds", "data/processed/04_GSE41762_pheno.rds",
          "data/processed/04_GSE27951_expr.rds", "data/processed/04_GSE27951_pheno.rds", "data/processed/04_GSE27951_expressed.rds",
          "genesets/senescence_union.rds", "manifest/GSE38642_donor_overlap.csv",
          "results/stage3/06_nested_cv_islet_result.rds", "results/stage3/stability_selection.csv", "results/stage3/06_final_panel_islet.csv")
for (f in f_in) cat("  ", f, " md5=", tools::md5sum(here(f)), "\n", sep = "")
cat("  (validation cohorts GSE76894, GSE20966, GSE71416, GSE16415 and GSE38642 expression are NOT loaded; GSE38642 manifest is sample metadata only)\n")

geneset <- readRDS(here("genesets/senescence_union.rds")); stopifnot(is.character(geneset), length(geneset) == 1151L, !anyDuplicated(geneset))
load_arm <- function(gse, arm) {
  e <- readRDS(here("data/processed", sprintf("04_%s_expr.rds", gse))); p <- readRDS(here("data/processed", sprintf("04_%s_pheno.rds", gse)))
  stopifnot(identical(colnames(e), p$gsm)); list(expr = e, pheno = prep_pheno(p, arm), raw = p)
}
add_bmi <- function(ph, raw, design) {              # BMI as an extra continuous limma covariate (protocol §2.4 sensitivity)
  ph$bmi <- as.numeric(raw$bmi[match(rownames(ph), raw$gsm)]); stopifnot(!anyNA(ph$bmi))
  list(pheno = ph, design = as.formula(paste(paste(deparse(design), collapse = ""), "+ bmi"), env = globalenv()))
}

# ---- data and parameters -----------------------------------------------------------------------------------------
ad <- load_arm("GSE27951", "adipose"); PA <- default_params("adipose")
stopifnot(ncol(ad$expr) == 23L, table(ad$pheno$status)[["T2D"]] == 12L, table(ad$pheno$status)[["ND"]] == 11L)
isl <- load_arm("GSE41762", "islet"); PI <- default_params("islet")
stopifnot(ncol(isl$expr) == 77L, table(isl$pheno$status)[["T2D"]] == 20L, table(isl$pheno$status)[["ND"]] == 57L)
cat("\nPre-run checks: adipose GSE27951 n =", ncol(ad$expr), "genes =", nrow(ad$expr), " ND/T2D =", paste(table(ad$pheno$status), collapse = "/"),
    "; sex (inferred, no recorded sex) F/M =", paste(table(ad$pheno$sex), collapse = "/"), "; age range", paste(range(ad$pheno$age), collapse = "-"), "\n")
fl <- readRDS(here("data/processed/04_GSE27951_expressed.rds"))
q <- quantile(ad$expr, PA$expr_q, names = FALSE); mine <- rownames(ad$expr)[rowMeans(ad$expr > q) >= PA$expr_frac]
stopifnot(setequal(mine, fl$gene[fl$expressed]))
cat("In-fold expression filter applied to the full adipose cohort reproduces 04_GSE27951_expressed.rds exactly (", length(mine), " genes); union genes expressed: ",
    sum(geneset %in% mine), "\n", sep = "")
cat("\nParameters (adipose): outer ", PA$outer_k, "-fold x ", PA$outer_reps, " reps; inner ", PA$inner_k, "-fold; design ", deparse(PA$design),
    "; screen p<", PA$screen_p, " (fallback top ", PA$screen_top, " if <", PA$screen_min, "); SVM-RFE sizes 1..", PA$rfe_max,
    "; Boruta maxRuns=", PA$boruta_runs, " ntree=", PA$boruta_ntree, " (fru backend); ridge lambda.min; LASSO lambda.1se; seed ", SEED, "\n", sep = "")
cat("Parameters (islet sensitivity): outer ", PI$outer_k, "-fold x ", PI$outer_reps, " reps; inner ", PI$inner_k, "-fold; design ", deparse(PI$design), "; seed ", SEED, "\n", sep = "")

# ---- A1. adipose permutation leakage test -----------------------------------------------------------------------
cat("\n=== A1. adipose permutation leakage test (20 label shuffles, full pipeline, 1 repeat of 5-fold outer CV; same code path as islet) ===\n")
PP <- PA; PP$outer_reps <- 1L
t_perm <- system.time(perm <- par_map(seq_len(20), function(i) {
  set.seed(4200L + i); ph <- ad$pheno; ph$status <- sample(ph$status)            # independent shuffle of status; covariates untouched
  stopifnot(table(ph$status) == table(ad$pheno$status))
  r <- nested_cv(ad$expr, ph, geneset, PP, seed = SEED + 100L * i, cores = 1L, keep_folds = FALSE)
  c(auc = r$auc_mean, same_as_real = mean(ph$status == ad$pheno$status),
    fb = mean(vapply(r$folds, `[[`, "", "level") != "L0_intersection"), nf = mean(vapply(r$folds, function(f) length(f$features), 0)))
}, CORES))[3]
pm <- do.call(rbind, perm); pa <- pm[, "auc"]
perm_res <- list(auc = pa, mean = mean(pa), sd = sd(pa), q95 = unname(quantile(pa, .95)), min = min(pa), max = max(pa), median = median(pa))
perm_pass <- perm_res$mean >= 0.45 && perm_res$mean <= 0.55 && perm_res$q95 < 0.65
write.csv(data.frame(shuffle = 1:20, auc = round(pa, 4), prop_labels_equal_to_real = round(pm[, "same_as_real"], 3),
                     frac_folds_fallback = round(pm[, "fb"], 2), mean_panel_size = round(pm[, "nf"], 2)),
          file.path(out_dir, "06b_permutation_test_adipose.csv"), row.names = FALSE)
saveRDS(c(perm_res, list(pass = perm_pass, minutes = t_perm / 60)), file.path(out_dir, "06b_permutation_test_adipose.rds"))
cat("20 permutation AUCs:", paste(sprintf("%.3f", sort(pa)), collapse = " "), "\n")
cat(sprintf("mean %.3f  median %.3f  SD %.3f  95th pct %.3f  min %.3f  max %.3f  (%.1f min)\n", perm_res$mean, perm_res$median, perm_res$sd, perm_res$q95,
            perm_res$min, perm_res$max, t_perm / 60))
cat("Criteria (fixed): mean in [0.45,0.55] AND 95th percentile < 0.65 -> ", if (perm_pass) "PASS" else "FAIL", "\n", sep = "")
draw_perm <- function() {
  stripchart(pa, vertical = TRUE, method = "jitter", pch = 19, col = "grey30", ylim = c(0.1, 0.9), ylab = "nested-CV AUC (shuffled labels)",
             main = "Adipose permutation test (20 shuffles)", xlim = c(0.5, 1.5))
  abline(h = c(0.45, 0.55), lty = 3); abline(h = 0.65, lty = 2, col = "red"); points(1, mean(pa), pch = 4, cex = 2, lwd = 2, col = "blue")
}
pdf(file.path(fig_dir, "06b_adipose_permutation.pdf"), 4.5, 4.5); draw_perm(); dev.off()
png(file.path(fig_dir, "06b_adipose_permutation.png"), 900, 900, res = 200); draw_perm(); dev.off()
if (!perm_pass) {
  cat("\n*** ADIPOSE PERMUTATION TEST FAILED. STOPPING. No headline numbers or sensitivity runs are produced. Thresholds and pipeline NOT adjusted. ***\n")
  finish(2L)
}

# ---- generic runner ----------------------------------------------------------------------------------------------
run_one <- function(label, expr, pheno, P, geneset) {
  cat("\n=== ", label, ": n = ", ncol(expr), " (ND/T2D ", paste(table(pheno$status), collapse = "/"), "); design ", deparse(P$design),
      "; outer ", P$outer_k, "x", P$outer_reps, ", inner ", P$inner_k, " ===\n", sep = "")
  t_real <- system.time(real <- nested_cv(expr, pheno, geneset, P, seed = SEED, cores = CORES))[3]
  cat(sprintf("nested CV done in %.1f min; folds = %d\n", t_real / 60, real$n_folds))
  stab <- stability_table(real$folds)
  lv <- vapply(real$folds, `[[`, "", "level")
  fb <- data.frame(level = c("L0_intersection", "L1_>=2of3", "L2_top5_LASSO", "L3_top5_screen"))
  fb$n_folds <- as.integer(table(factor(lv, fb$level))); fb$frequency <- fb$n_folds / length(lv)
  sfb <- data.frame(item = c("screen_fallback_top30 (<10 pass p<0.05)", "LASSO_empty"),
                    n_folds = c(sum(vapply(real$folds, `[[`, NA, "screen_fallback")), sum(vapply(real$folds, function(f) length(f$lasso) == 0, NA))))
  sfb$frequency <- sfb$n_folds / length(lv)
  yst <- setNames(as.character(pheno$status), rownames(pheno))
  pf <- do.call(rbind, lapply(real$folds, function(f) data.frame(
    run = label, rep = f$rep, fold = f$fold, n_train = f$n_train, n_test = length(f$test_ids), n_T2D_test = sum(yst[f$test_ids] == "T2D"),
    n_ND_test = sum(yst[f$test_ids] == "ND"), n_screen = length(f$screen), screen_fallback = f$screen_fallback, n_lasso = length(f$lasso),
    n_svmrfe = length(f$svmrfe), n_boruta = length(f$boruta), boruta_tentative = f$boruta_tentative, level = f$level, panel_size = length(f$features),
    panel = paste(f$features, collapse = ";"), n_warnings = length(f$warnings), stringsAsFactors = FALSE)))
  cat("fallback levels:\n"); print(fb, row.names = FALSE); print(sfb, row.names = FALSE)
  cat(sprintf("panel size per fold: median %g (range %g-%g); screen genes per fold median %g (range %g-%g); Boruta tentative per fold median %g; warnings in %d folds\n",
              median(pf$panel_size), min(pf$panel_size), max(pf$panel_size), median(pf$n_screen), min(pf$n_screen), max(pf$n_screen),
              median(pf$boruta_tentative), sum(pf$n_warnings > 0)))
  cat("T2D per outer test fold (all ", nrow(pf), " folds): ", paste(names(table(pf$n_T2D_test)), table(pf$n_T2D_test), sep = " T2D x", collapse = "; "),
      "; test-fold size range ", min(pf$n_test), "-", max(pf$n_test), "; folds with 0 T2D in test: ", sum(pf$n_T2D_test == 0), "\n", sep = "")
  wmsg <- unique(unlist(lapply(real$folds, `[[`, "warnings"))); if (length(wmsg)) cat("distinct warning messages:\n", paste(" -", wmsg, collapse = "\n"), "\n")
  sm <- summarise_auc(real$oof, P$outer_reps, yst)
  set.seed(SEED); final <- fit_pipeline(expr, pheno, geneset, P)
  p_app <- predict_model(final$model, expr); r_app <- roc_obj(pheno$status, p_app)
  set.seed(SEED); ci_app <- as.numeric(pROC::ci.auc(r_app, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE, progress = "none"))
  stopifnot(abs(as.numeric(pROC::auc(r_app)) - fast_auc(p_app, pheno$status)) < 1e-9)
  cat("final panel (", length(final$features), " genes, ", final$level, "): ", paste(final$features, collapse = ", "), "\n", sep = "")
  cat(sprintf("selected by LASSO %d, SVM-RFE %d, Boruta %d (screen %d genes)\n", length(final$lasso), length(final$svmrfe), length(final$boruta), length(final$screen)))
  cat(sprintf("NESTED-CV AUC (honest): %.3f  [per-repeat AUCs: %s]\n", sm$auc_mean_over_repeats, paste(sprintf("%.3f", sm$auc_per_repeat), collapse = ", ")))
  cat(sprintf("   95%% CI, sample-level stratified bootstrap (HEADLINE; conservative): %.3f - %.3f\n", sm$ci_sample_level[1], sm$ci_sample_level[2]))
  cat(sprintf("   95%% CI, protocol (pROC::ci.auc, 2000 stratified bootstraps on pooled OOF predictions): %.3f - %.3f   [pooled-stack AUC %.3f]\n",
              sm$ci_protocol[1], sm$ci_protocol[2], sm$auc_pooled))
  cat(sprintf("APPARENT AUC, optimistic (apparent): %.3f   95%% CI (pROC, 2000 bootstraps) %.3f - %.3f\n", as.numeric(pROC::auc(r_app)), ci_app[1], ci_app[3]))
  cat("Top 20 genes by stability (final-panel selection frequency over ", real$n_folds, " outer folds):\n", sep = ""); print(head(stab, 20), row.names = FALSE)
  list(label = label, real = real, stab = stab, fb = fb, sfb = sfb, pf = pf, sm = sm, final = final, auc_app = as.numeric(pROC::auc(r_app)), ci_app = ci_app[c(1, 3)],
       r_app = r_app, minutes = t_real / 60, n = ncol(expr), n_t2d = sum(pheno$status == "T2D"), P = P)
}

save_run <- function(R, stab_file) {
  l <- R$label
  write.csv(R$stab[, c("gene", "frequency", "n_folds_final", "freq_screen", "freq_lasso", "freq_svmrfe", "freq_boruta")], stab_file, row.names = FALSE)
  write.csv(rbind(cbind(kind = "fallback_level", R$fb), cbind(kind = "other", setNames(R$sfb, names(R$fb)))),
            file.path(out_dir, sprintf("06b_fallback_frequency_%s.csv", l)), row.names = FALSE)
  saveRDS(lapply(R$real$folds, function(f) f[c("rep", "fold", "train_ids", "test_ids", "n_train", "limma_ncol", "n_expressed", "n_union_expressed", "n_pass_nominal",
                                              "screen_fallback", "screen", "lasso", "svmrfe", "svm_size", "boruta", "boruta_tentative", "boruta_runs", "features", "level", "lambda", "warnings", "elapsed")]),
          file.path(out_dir, sprintf("06b_per_fold_features_%s.rds", l)))
  write.csv(R$real$oof, file.path(out_dir, sprintf("06b_oof_predictions_%s.csv", l)), row.names = FALSE)
  fm <- R$final$model
  saveRDS(list(genes = fm$genes, coef = fm$coef, intercept = fm$intercept, mu = fm$mu, sd = fm$sd, lambda = fm$lambda, fallback_level = R$final$level,
               lasso = R$final$lasso, svmrfe = R$final$svmrfe, boruta = R$final$boruta, screen = R$final$screen, params = R$P[names(R$P) != "design"],
               design = deparse(R$P$design), source = sprintf("%s (04_ post-QC matrix); frozen model for protocol §7 (only the main adipose run is intended for validation)", l)),
          file.path(out_dir, sprintf("06b_final_model_%s.rds", l)))
  write.csv(data.frame(gene = fm$genes, ridge_coef_z = fm$coef, train_mean = fm$mu, train_sd = fm$sd, fallback_level = R$final$level),
            file.path(out_dir, sprintf("06b_final_panel_%s.csv", l)), row.names = FALSE)
  draw <- function() {
    op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1)); on.exit(par(op))
    plot(R$sm$roc, col = "navy", lwd = 2, legacy.axes = TRUE, main = sprintf("%s: out-of-fold vs apparent", l), cex.main = 0.9)
    plot(R$r_app, col = "firebrick", lwd = 2, lty = 2, add = TRUE, legacy.axes = TRUE)
    legend("bottomright", bty = "n", cex = 0.7, col = c("navy", "firebrick"), lwd = 2, lty = c(1, 2),
           legend = c(sprintf("nested CV (honest) %.2f [%.2f-%.2f]", R$sm$auc_mean_over_repeats, R$sm$ci_sample_level[1], R$sm$ci_sample_level[2]),
                      sprintf("optimistic (apparent) %.2f", R$auc_app)))
    top <- head(R$stab, 15); par(mar = c(4.2, 6, 3, 1))
    barplot(rev(top$frequency), names.arg = rev(top$gene), horiz = TRUE, las = 1, xlim = c(0, 1), col = "steelblue", xlab = "final-panel selection frequency",
            main = sprintf("Stability (top 15; %d folds)", R$real$n_folds), cex.names = 0.7, cex.main = 0.9)
    abline(v = 0.5, lty = 3)
  }
  pdf(file.path(fig_dir, sprintf("06b_%s_roc_stability.pdf", l)), 9, 4.8); draw(); dev.off()
  png(file.path(fig_dir, sprintf("06b_%s_roc_stability.png", l)), 1800, 960, res = 200); draw(); dev.off()
}

# ---- A2. adipose main --------------------------------------------------------------------------------------------
cat("\n\n################ A2. ADIPOSE MAIN (GSE27951 alone) ################\n")
RA <- run_one("adipose_main", ad$expr, ad$pheno, PA, geneset)
save_run(RA, file.path(out_dir, "stability_selection_adipose.csv"))

# ---- B. sensitivity runs -----------------------------------------------------------------------------------------
cat("\n\n################ B. PRE-SPECIFIED SENSITIVITY RUNS ################\n")
ov <- read.csv(here("manifest/GSE38642_donor_overlap.csv"), stringsAsFactors = FALSE)
disc <- ov[!is.na(ov$sameID_GSE41762_gsm) & nzchar(ov$sameID_GSE41762_gsm) & ov$group != ov$sameID_GSE41762_group,
           c("gsm", "group", "sameID_GSE41762_gsm", "sameID_GSE41762_group")]
names(disc) <- c("GSE38642_gsm", "GSE38642_group", "GSE41762_gsm", "GSE41762_group")
disc$GSE41762_group_in_analysis <- as.character(isl$pheno[disc$GSE41762_gsm, "status"])
stopifnot(nrow(disc) == 11L, !anyDuplicated(disc$GSE41762_gsm), all(disc$GSE41762_gsm %in% rownames(isl$pheno)),
          identical(disc$GSE41762_group_in_analysis, disc$GSE41762_group))
cat("11 GSE41762 donors with a different diabetes label in GSE38642 (sample-metadata labels only; ", sum(disc$GSE41762_group == "T2D"), " T2D / ",
    sum(disc$GSE41762_group == "ND"), " ND in GSE41762):\n", sep = ""); print(disc, row.names = FALSE)
write.csv(disc, file.path(out_dir, "06b_islet_excluded_discordant_donors.csv"), row.names = FALSE)

keep <- setdiff(colnames(isl$expr), disc$GSE41762_gsm)
ph1 <- droplevels(isl$pheno[keep, , drop = FALSE]); ex1 <- isl$expr[, keep, drop = FALSE]
stopifnot(ncol(ex1) == 66L, nlevels(ph1$batch) == 2L, nlevels(ph1$sex) == 2L)
R1 <- run_one("islet_excl11", ex1, ph1, PI, geneset)

b2 <- add_bmi(isl$pheno, isl$raw, PI$design); PI2 <- PI; PI2$design <- b2$design
R2 <- run_one("islet_bmi", isl$expr, b2$pheno, PI2, geneset)

keep <- setdiff(colnames(ad$expr), "GSM691134"); stopifnot(length(keep) == 22L, "GSM691134" %in% colnames(ad$expr), ad$pheno["GSM691134", "status"] == "ND")
ph3 <- droplevels(ad$pheno[keep, , drop = FALSE]); ex3 <- ad$expr[, keep, drop = FALSE]
R3 <- run_one("adipose_excl_GSM691134", ex3, ph3, PA, geneset)

b4 <- add_bmi(ad$pheno, ad$raw, PA$design); PA4 <- PA; PA4$design <- b4$design
R4 <- run_one("adipose_bmi", ad$expr, b4$pheno, PA4, geneset)

for (R in list(R1, R2, R3, R4)) save_run(R, file.path(out_dir, sprintf("06b_stability_selection_%s.csv", R$label)))

# ---- comparison with main runs -----------------------------------------------------------------------------------
im <- readRDS(here("results/stage3/06_nested_cv_islet_result.rds")); ist <- read.csv(here("results/stage3/stability_selection.csv"))
ipan <- read.csv(here("results/stage3/06_final_panel_islet.csv"))$gene
main <- list(islet = list(auc = im$nested$auc_mean_over_repeats, ci_p = im$nested$ci_protocol, ci_s = im$nested$ci_sample_level, app = im$apparent$auc,
                          panel = ipan, top10 = head(ist$gene, 10)),
             adipose = list(auc = RA$sm$auc_mean_over_repeats, ci_p = RA$sm$ci_protocol, ci_s = RA$sm$ci_sample_level, app = RA$auc_app,
                            panel = RA$final$features, top10 = head(RA$stab$gene, 10)))
runs <- list(RA, R1, R2, R3, R4)
arm_of <- c(adipose_main = "adipose", islet_excl11 = "islet", islet_bmi = "islet", adipose_excl_GSM691134 = "adipose", adipose_bmi = "adipose")
auc_tab <- do.call(rbind, lapply(runs, function(R) data.frame(run = R$label, n = R$n, n_T2D = R$n_t2d, nested_auc = R$sm$auc_mean_over_repeats,
  ci_lo_sample_level_headline = R$sm$ci_sample_level[1], ci_hi_sample_level_headline = R$sm$ci_sample_level[2],
  ci_lo_protocol = R$sm$ci_protocol[1], ci_hi_protocol = R$sm$ci_protocol[2], pooled_stack_auc = R$sm$auc_pooled,
  apparent_auc_optimistic = R$auc_app, apparent_ci_lo = R$ci_app[1], apparent_ci_hi = R$ci_app[2],
  per_repeat_auc_min = min(R$sm$auc_per_repeat), per_repeat_auc_max = max(R$sm$auc_per_repeat),
  final_panel = paste(R$final$features, collapse = ";"), final_level = R$final$level, minutes = R$minutes, stringsAsFactors = FALSE)))
rownames(auc_tab) <- NULL; write.csv(auc_tab, file.path(out_dir, "06b_auc_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(runs, `[[`, "pf")), file.path(out_dir, "06b_per_fold_summary.csv"), row.names = FALSE)
top10_tab <- do.call(rbind, lapply(runs, function(R) cbind(run = R$label, head(R$stab[, c("gene", "frequency", "n_folds_final", "freq_screen", "freq_lasso", "freq_svmrfe", "freq_boruta")], 10))))
write.csv(top10_tab, file.path(out_dir, "06b_top10_stability.csv"), row.names = FALSE)

cmp <- do.call(rbind, lapply(runs[-1], function(R) {
  m <- main[[arm_of[[R$label]]]]; fp <- R$final$features; t10 <- head(R$stab$gene, 10)
  data.frame(run = R$label, vs_main = paste0(arm_of[[R$label]], "_main"), n = R$n, n_T2D = R$n_t2d,
             main_nested_auc = m$auc, main_ci_sample = sprintf("%.3f-%.3f", m$ci_s[1], m$ci_s[2]), main_ci_protocol = sprintf("%.3f-%.3f", m$ci_p[1], m$ci_p[2]), main_apparent = m$app,
             sens_nested_auc = R$sm$auc_mean_over_repeats, sens_ci_sample = sprintf("%.3f-%.3f", R$sm$ci_sample_level[1], R$sm$ci_sample_level[2]),
             sens_ci_protocol = sprintf("%.3f-%.3f", R$sm$ci_protocol[1], R$sm$ci_protocol[2]), sens_apparent = R$auc_app,
             auc_diff_sens_minus_main = R$sm$auc_mean_over_repeats - m$auc,
             main_panel = paste(m$panel, collapse = ";"), sens_panel = paste(fp, collapse = ";"),
             panel_shared = paste(intersect(m$panel, fp), collapse = ";"), panel_jaccard = length(intersect(m$panel, fp)) / length(union(m$panel, fp)),
             top10_overlap_n = length(intersect(m$top10, t10)), top10_shared = paste(intersect(m$top10, t10), collapse = ";"),
             sens_top10 = paste(t10, collapse = ";"), main_top10 = paste(m$top10, collapse = ";"), stringsAsFactors = FALSE)
}))
write.csv(cmp, file.path(out_dir, "06b_sensitivity_comparison.csv"), row.names = FALSE)

cat("\n\n================ RESULTS SUMMARY ================\n")
cat(sprintf("ISLET MAIN (read from results/stage3, not rerun): nested %.3f; sample-level CI %.3f-%.3f; protocol CI %.3f-%.3f; apparent %.3f; panel %s; top10 %s\n",
            main$islet$auc, main$islet$ci_s[1], main$islet$ci_s[2], main$islet$ci_p[1], main$islet$ci_p[2], main$islet$app, paste(main$islet$panel, collapse = ","),
            paste(main$islet$top10, collapse = ",")))
print(auc_tab[, c("run", "n", "n_T2D", "nested_auc", "ci_lo_sample_level_headline", "ci_hi_sample_level_headline", "ci_lo_protocol", "ci_hi_protocol",
                  "apparent_auc_optimistic", "apparent_ci_lo", "apparent_ci_hi", "final_panel")], row.names = FALSE)
cat("\nComparison with main runs:\n"); print(cmp[, c("run", "main_nested_auc", "sens_nested_auc", "auc_diff_sens_minus_main", "main_panel", "sens_panel", "panel_jaccard", "top10_overlap_n")], row.names = FALSE)
cat("\nTop-10 stability genes per run:\n"); for (R in runs) cat(sprintf("  %-24s %s\n", R$label, paste(sprintf("%s(%.2f)", head(R$stab$gene, 10), head(R$stab$frequency, 10)), collapse = ", ")))
cat(sprintf("  %-24s %s\n", "islet_main (reference)", paste(sprintf("%s(%.2f)", ist$gene[1:10], ist$frequency[1:10]), collapse = ", ")))

saveRDS(list(permutation = c(perm_res, list(pass = perm_pass)), auc_summary = auc_tab, comparison = cmp, discordant_donors = disc, params_adipose = PA[names(PA) != "design"],
             params_islet = PI[names(PI) != "design"], runtime_min = setNames(vapply(runs, `[[`, 0, "minutes"), vapply(runs, `[[`, "", "label")),
             permutation_min = t_perm / 60), file.path(out_dir, "06b_nested_cv_arms_result.rds"))

# ---- final assertions --------------------------------------------------------------------------------------------
stopifnot(length(perm) == 20L, perm_pass)
for (R in runs) {
  stopifnot(R$real$n_folds == 50L, nrow(R$real$oof) == R$n * R$P$outer_reps, !anyNA(R$real$oof$prob), sum(R$fb$n_folds) == 50L, nrow(R$stab) > 0,
            all(R$stab$frequency >= 0 & R$stab$frequency <= 1), length(R$final$features) >= 1, all(R$final$features %in% geneset),
            R$sm$auc_mean_over_repeats >= 0 && R$sm$auc_mean_over_repeats <= 1, !anyNA(c(R$sm$ci_sample_level, R$sm$ci_protocol, R$ci_app)),
            all(R$pf$n_train + R$pf$n_test == R$n), all(vapply(R$real$folds, function(f) f$limma_ncol == f$n_train, NA)),
            all(vapply(R$real$folds, function(f) !any(f$train_ids %in% f$test_ids), NA)), nrow(R$pf) == 50L)
}
stopifnot(RA$n == 23L, RA$n_t2d == 12L, ncol(ex1) == 66L, ncol(ex3) == 22L, nrow(cmp) == 4L, !anyNA(cmp$sens_nested_auc),
          file.exists(file.path(out_dir, "stability_selection_adipose.csv")),
          !any(grepl("GSM691134", colnames(ex3))), !any(disc$GSE41762_gsm %in% colnames(ex1)),
          identical(read.csv(here("results/stage3/stability_selection.csv"))$gene, ist$gene))      # islet stability file untouched
cat("\nAll end-of-script assertions passed (5 runs x 50 folds; OOF rows = n x repeats; no NA; fallback tables sum to 50; excluded samples absent; islet stability file untouched).\n")
cat(sprintf("\nRuntime: permutation %.1f min; adipose main %.1f min; islet_excl11 %.1f; islet_bmi %.1f; adipose_excl_GSM691134 %.1f; adipose_bmi %.1f; total %.1f min\n",
            t_perm / 60, RA$minutes, R1$minutes, R2$minutes, R3$minutes, R4$minutes, as.numeric(difftime(Sys.time(), T_START, units = "mins"))))
cat("\nsessionInfo():\n"); print(sessionInfo())
finish(0L)

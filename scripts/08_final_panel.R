# 08_final_panel.R -- assemble the final hub-gene panel results per arm (protocol §3.2, §8; CLAUDE.md rules 2, 5, 6).
#
# DISCOVERY DATA ONLY. This script recomputes nothing: every number is read from the saved Task 5/5b outputs in results/stage3/
# (the full-cohort pipeline, nested CV, permutation tests and the four pre-specified sensitivity runs were already executed there).
# Validation cohorts are never opened here.
#
# Output: results/stage5/08_final_panels.csv (one row per run: main arms first, then sensitivity runs), results/stage5/08_final_panels.md.
# Every nested-CV AUC is shown with BOTH CIs (headline = sample-level stratified bootstrap; protocol = pROC::ci.auc on pooled OOF
# predictions) and always beside the apparent AUC, labelled "optimistic (apparent)" (CLAUDE.md rule 5).
# Log: logs/08_final_panel_<timestamp>.log

set.seed(42)
suppressPackageStartupMessages(library(here))
script_name <- "08_final_panel"; ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create(here("logs"), showWarnings = FALSE); out_dir <- here("results", "stage5"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt"); sink(con, split = TRUE); sink(con, type = "message")
cat("Script:", script_name, " run:", ts, "\nParameters: none (read-only assembly of saved discovery results); seed 42\n\nInputs:\n")
f_in <- c("results/stage3/06_final_panel_islet.csv", "results/stage3/06_auc_summary_islet.csv", "results/stage3/06_fallback_frequency_islet.csv",
          "results/stage3/stability_selection.csv", "results/stage3/06_nested_cv_islet_result.rds", "results/stage3/06_permutation_test_islet.csv",
          "results/stage3/06b_final_panel_adipose_main.csv", "results/stage3/06b_auc_summary.csv", "results/stage3/06b_fallback_frequency_adipose_main.csv",
          "results/stage3/stability_selection_adipose.csv", "results/stage3/06b_permutation_test_adipose.rds", "results/stage3/06b_sensitivity_comparison.csv")
for (r in c("islet_excl11", "islet_bmi", "adipose_excl_GSM691134", "adipose_bmi")) {
  f_in <- c(f_in, sprintf("results/stage3/06b_final_panel_%s.csv", r), sprintf("results/stage3/06b_fallback_frequency_%s.csv", r),
            sprintf("results/stage3/06b_stability_selection_%s.csv", r))
}
for (f in f_in) { stopifnot(file.exists(here(f))); cat("  ", f, " md5=", tools::md5sum(here(f)), "\n", sep = "") }

rd <- function(f) read.csv(here("results/stage3", f), stringsAsFactors = FALSE, check.names = FALSE)
fmt <- function(x, d = 3) sprintf(paste0("%.", d, "f"), x)
ci_txt <- function(lo, hi) sprintf("%s-%s", fmt(lo), fmt(hi))

# ---- per-run assembly ---------------------------------------------------------------------------------------------------
fb_row <- function(fb, lev) fb$frequency[fb$kind == "fallback_level" & fb$level == lev]
fb_oth <- function(fb, key) fb$frequency[fb$kind == "other" & startsWith(fb$level, key)]
stab_txt <- function(stab, genes) paste(sprintf("%s:%s", genes, fmt(stab$frequency[match(genes, stab$gene)], 2)), collapse = ";")

# main arms ------------------------------------------------------------------------------------------------------------------
isl_panel <- rd("06_final_panel_islet.csv"); isl_auc <- rd("06_auc_summary_islet.csv"); isl_fb <- rd("06_fallback_frequency_islet.csv"); isl_st <- rd("stability_selection.csv")
isl_res <- readRDS(here("results/stage3/06_nested_cv_islet_result.rds"))
ada_all <- rd("06b_auc_summary.csv"); ada <- ada_all[ada_all$run == "adipose_main", ]; ada_panel <- rd("06b_final_panel_adipose_main.csv")
ada_fb <- rd("06b_fallback_frequency_adipose_main.csv"); ada_st <- rd("stability_selection_adipose.csv")
ada_perm <- readRDS(here("results/stage3/06b_permutation_test_adipose.rds"))

ih <- isl_auc[startsWith(isl_auc$metric, "nested-CV AUC (honest"), ]; ia <- isl_auc[startsWith(isl_auc$metric, "apparent"), ]
stopifnot(nrow(ih) == 1L, nrow(ia) == 1L, nrow(ada) == 1L,
          abs(ih$auc - isl_res$nested$auc_mean_over_repeats) < 1e-12, abs(ia$auc - isl_res$apparent$auc) < 1e-12)       # csv vs rds agree

main_row <- function(arm, run, n, n_t2d, panel, nested, ci_s, ci_p, app, ci_app, fb, stab, perm) data.frame(
  arm = arm, run = run, role = "main", n = n, n_T2D = n_t2d, final_panel = paste(panel$gene, collapse = ";"), final_panel_size = nrow(panel),
  full_cohort_fallback_level = unique(panel$fallback_level), ridge_coef_z = paste(fmt(panel$ridge_coef_z, 3), collapse = ";"),
  panel_gene_stability_freq = stab_txt(stab, panel$gene),
  nested_cv_auc_headline = nested, nested_ci_lo_sample_level_HEADLINE = ci_s[1], nested_ci_hi_sample_level_HEADLINE = ci_s[2],
  nested_ci_lo_protocol_pROC = ci_p[1], nested_ci_hi_protocol_pROC = ci_p[2],
  apparent_auc_optimistic = app, apparent_ci_lo_pROC = ci_app[1], apparent_ci_hi_pROC = ci_app[2],
  freq_fold_level_L0_intersection = fb_row(fb, "L0_intersection"), freq_fold_level_L1_2of3 = fb_row(fb, "L1_>=2of3"),
  freq_fold_level_L2_top5_LASSO = fb_row(fb, "L2_top5_LASSO"), freq_fold_level_L3_top5_screen = fb_row(fb, "L3_top5_screen"),
  freq_screen_fallback_top30 = fb_oth(fb, "screen_fallback"), freq_LASSO_empty = fb_oth(fb, "LASSO_empty"),
  permutation_mean_auc = perm$mean, permutation_q95_auc = perm$q95, stringsAsFactors = FALSE)

rows <- list(
  main_row("islet", "islet_main (GSE41762)", 77L, 20L, isl_panel, ih$auc, c(ih$ci_lo_sample_level, ih$ci_hi_sample_level), c(ih$ci_lo_protocol, ih$ci_hi_protocol),
           ia$auc, c(ia$ci_lo_protocol, ia$ci_hi_protocol), isl_fb, isl_st, isl_res$permutation),
  main_row("adipose", "adipose_main (GSE27951)", ada$n, ada$n_T2D, ada_panel, ada$nested_auc, c(ada$ci_lo_sample_level_headline, ada$ci_hi_sample_level_headline),
           c(ada$ci_lo_protocol, ada$ci_hi_protocol), ada$apparent_auc_optimistic, c(ada$apparent_ci_lo, ada$apparent_ci_hi), ada_fb, ada_st, ada_perm))
# the 06b run for adipose_main must agree with the final-panel file
stopifnot(identical(ada$final_panel, paste(ada_panel$gene, collapse = ";")), identical(paste(isl_panel$gene, collapse = ";"), "IL1R1"))

# sensitivity runs ------------------------------------------------------------------------------------------------------------
sens <- c(islet_excl11 = "islet", islet_bmi = "islet", adipose_excl_GSM691134 = "adipose", adipose_bmi = "adipose")
for (r in names(sens)) {
  a <- ada_all[ada_all$run == r, ]; pn <- rd(sprintf("06b_final_panel_%s.csv", r)); fb <- rd(sprintf("06b_fallback_frequency_%s.csv", r)); st <- rd(sprintf("06b_stability_selection_%s.csv", r))
  stopifnot(nrow(a) == 1L, identical(a$final_panel, paste(pn$gene, collapse = ";")))
  x <- main_row(sens[[r]], r, a$n, a$n_T2D, pn, a$nested_auc, c(a$ci_lo_sample_level_headline, a$ci_hi_sample_level_headline), c(a$ci_lo_protocol, a$ci_hi_protocol),
                a$apparent_auc_optimistic, c(a$apparent_ci_lo, a$apparent_ci_hi), fb, st, list(mean = NA_real_, q95 = NA_real_))
  x$role <- "sensitivity (pre-specified, protocol §8)"; rows[[length(rows) + 1L]] <- x
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL

# add sensitivity deltas versus the main run of the same arm (read from the saved comparison file)
cmp <- rd("06b_sensitivity_comparison.csv")
tab$sens_nested_auc_minus_main <- NA_real_; tab$sens_panel_jaccard_vs_main <- NA_real_; tab$sens_top10_overlap_n_vs_main <- NA_integer_; tab$sens_panel_shared_with_main <- NA_character_
for (i in seq_len(nrow(cmp))) {
  j <- which(tab$run == cmp$run[i]); stopifnot(length(j) == 1L, abs(tab$nested_cv_auc_headline[j] - cmp$sens_nested_auc[i]) < 1e-9)
  tab$sens_nested_auc_minus_main[j] <- cmp$auc_diff_sens_minus_main[i]; tab$sens_panel_jaccard_vs_main[j] <- cmp$panel_jaccard[i]
  tab$sens_top10_overlap_n_vs_main[j] <- cmp$top10_overlap_n[i]; tab$sens_panel_shared_with_main[j] <- cmp$panel_shared[i]
}
write.csv(tab, file.path(out_dir, "08_final_panels.csv"), row.names = FALSE)

# ---- markdown -----------------------------------------------------------------------------------------------------------------
md <- c("# Final hub-gene panels (discovery cohorts only)", "",
        "Generated by `scripts/08_final_panel.R` from saved Task 5/5b outputs (nothing recomputed; validation cohorts not opened).",
        "Nested-CV AUC is the honest headline. The apparent AUC (full-cohort model scored on its own training data) is **optimistic** and is shown only beside it (CLAUDE.md rule 5).",
        "Headline CI = sample-level stratified bootstrap (2000); the protocol CI = `pROC::ci.auc` on pooled out-of-fold predictions (2000 stratified bootstraps) and is anti-conservative because the repeats are not independent (deviations log 2026-10-05).", "",
        "## Main arms", "",
        "| Arm (cohort) | n (T2D) | Final panel (ridge coef, z-scale) | Full-cohort fallback level | Nested-CV AUC (honest) | 95% CI headline (sample-level) | 95% CI protocol (pROC) | Apparent AUC, optimistic (apparent) [pROC 95% CI] |",
        "|---|---|---|---|---|---|---|---|")
for (i in 1:2) with(tab[i, ], md <<- c(md, sprintf("| %s | %d (%d) | %s (%s) | %s | %s | %s | %s | %s [%s] |", run, n, n_T2D, final_panel, ridge_coef_z, full_cohort_fallback_level,
  fmt(nested_cv_auc_headline), ci_txt(nested_ci_lo_sample_level_HEADLINE, nested_ci_hi_sample_level_HEADLINE), ci_txt(nested_ci_lo_protocol_pROC, nested_ci_hi_protocol_pROC),
  fmt(apparent_auc_optimistic), ci_txt(apparent_ci_lo_pROC, apparent_ci_hi_pROC))))
md <- c(md, "", "## Stability, fallback frequencies and leakage test (main arms; 50 outer folds each)", "",
        "| Arm | Panel-gene stability (freq. in per-fold final panel) | L0 intersection | L1 >=2 of 3 | L2 top-5 LASSO | L3 top-5 screen | Screen fallback (<10 pass) | LASSO empty | Permutation test (20 shuffles): mean / 95th pct AUC |", "|---|---|---|---|---|---|---|---|---|")
for (i in 1:2) with(tab[i, ], md <<- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s / %s (pass) |", arm, panel_gene_stability_freq, fmt(freq_fold_level_L0_intersection, 2),
  fmt(freq_fold_level_L1_2of3, 2), fmt(freq_fold_level_L2_top5_LASSO, 2), fmt(freq_fold_level_L3_top5_screen, 2), fmt(freq_screen_fallback_top30, 2), fmt(freq_LASSO_empty, 2),
  fmt(permutation_mean_auc), fmt(permutation_q95_auc))))
md <- c(md, "", "Top-5 genes by stability, islet: ", paste0("`", paste(sprintf("%s %s", head(isl_st$gene, 5), fmt(head(isl_st$frequency, 5), 2)), collapse = ", "), "`"), "",
        "Top-5 genes by stability, adipose: ", paste0("`", paste(sprintf("%s %s", head(ada_st$gene, 5), fmt(head(ada_st$frequency, 5), 2)), collapse = ", "), "`"), "",
        "## Pre-specified sensitivity runs (protocol §8; descriptive, no decision made from them)", "",
        "| Arm | Run | n (T2D) | Final panel | Nested-CV AUC | CI headline | CI protocol | Apparent AUC, optimistic (apparent) | Nested AUC minus main | Panel genes shared with main | Panel gene stability |", "|---|---|---|---|---|---|---|---|---|---|---|")
for (i in 3:nrow(tab)) with(tab[i, ], md <<- c(md, sprintf("| %s | %s | %d (%d) | %s | %s | %s | %s | %s | %+.3f | %s | %s |", arm, run, n, n_T2D, final_panel, fmt(nested_cv_auc_headline),
  ci_txt(nested_ci_lo_sample_level_HEADLINE, nested_ci_hi_sample_level_HEADLINE), ci_txt(nested_ci_lo_protocol_pROC, nested_ci_hi_protocol_pROC), fmt(apparent_auc_optimistic),
  sens_nested_auc_minus_main, if (nzchar(sens_panel_shared_with_main)) sens_panel_shared_with_main else "none", panel_gene_stability_freq)))
md <- c(md, "", "Reading note (stated plainly): the islet panel is a single gene (IL1R1) and the adipose panel a single gene (SNAI1); both come from the full-cohort pipeline of protocol §3.2. ",
        "The adipose nested-CV CI includes 0.5 under both interval methods, and the adipose apparent AUC (0.871) is far above its nested AUC (0.532). ",
        "Panels are unstable across sensitivity runs (islet main vs excl11: no shared panel gene; adipose main vs excl GSM691134: no shared panel gene).", "")
writeLines(md, file.path(out_dir, "08_final_panels.md"))

cat("\n=== results (08_final_panels.csv) ===\n")
print(tab[, c("arm", "run", "n", "n_T2D", "final_panel", "nested_cv_auc_headline", "nested_ci_lo_sample_level_HEADLINE", "nested_ci_hi_sample_level_HEADLINE",
              "nested_ci_lo_protocol_pROC", "nested_ci_hi_protocol_pROC", "apparent_auc_optimistic")], row.names = FALSE)
cat("\nOutputs: results/stage5/08_final_panels.csv (md5 ", tools::md5sum(file.path(out_dir, "08_final_panels.csv")), "), results/stage5/08_final_panels.md\n", sep = "")

# ---- final assertions -----------------------------------------------------------------------------------------------------------
num_cols <- c("nested_cv_auc_headline", "nested_ci_lo_sample_level_HEADLINE", "nested_ci_hi_sample_level_HEADLINE", "nested_ci_lo_protocol_pROC", "nested_ci_hi_protocol_pROC",
              "apparent_auc_optimistic", "apparent_ci_lo_pROC", "apparent_ci_hi_pROC")
stopifnot(nrow(tab) == 6L, !anyNA(tab[, num_cols]), all(as.matrix(tab[, num_cols]) >= 0 & as.matrix(tab[, num_cols]) <= 1),
          all(tab$nested_ci_lo_sample_level_HEADLINE <= tab$nested_cv_auc_headline + 1e-9 & tab$nested_cv_auc_headline <= tab$nested_ci_hi_sample_level_HEADLINE + 1e-9),
          !anyNA(tab$final_panel_size), tab$n[1] == 77L, tab$n_T2D[1] == 20L, tab$n[2] == 23L, tab$n_T2D[2] == 12L,
          abs(sum(tab[1, grep("^freq_fold_level", names(tab))]) - 1) < 1e-9, abs(sum(tab[2, grep("^freq_fold_level", names(tab))]) - 1) < 1e-9,
          file.exists(file.path(out_dir, "08_final_panels.md")))
cat("\nAll end-of-script assertions passed (6 rows; AUCs and CIs in [0,1] with no NA; nested CI brackets the point estimate; fallback levels sum to 1; counts 77/20 and 23/12).\n")
cat("\nsessionInfo():\n"); print(sessionInfo())
sink(type = "message"); sink(); close(con)

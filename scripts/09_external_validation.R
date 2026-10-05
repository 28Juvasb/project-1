# 09_external_validation.R -- ONE-SHOT frozen-model external validation (protocol v1.2 §7; CLAUDE.md rules 2, 3, 5, 6, 7).
#
# Usage (from repo root):  Rscript scripts/09_external_validation.R [--dry-run-only | --verify-identical]
#
# Frozen models (never refit, never modified; md5 recorded in the lock):
#   islet   results/stage3/06_final_model_islet.rds          (single-gene ridge model: IL1R1)
#   adipose results/stage3/06b_final_model_adipose_main.rds   (single-gene ridge model: SNAI1)
# Scoring (protocol §7): panel genes matched by HGNC symbol; each gene is z-scored WITHIN the scored dataset (no training data used to
# normalise a validation set, rule 3); score = plogis(intercept + sum coef_j * z_j) with the frozen intercept and coefficients.
# Panel genes absent from a platform -> the cohort is reported "not scorable" (no gene substitution, no refit; panels are single genes).
#
# Cohorts: islet GSE76894 (primary) and GSE20966 (secondary): AUC; adipose GSE71416: AUC; adipose GSE16415: direction of effect ONLY (no AUC).
# AUC CIs (rule 6): HEADLINE = sample-level stratified bootstrap (2000, percentile); also pROC::ci.auc (2000 stratified bootstraps).
# ROC curves carry a pROC::ci.se confidence band (2000 stratified bootstraps). Direction of effect = mean z(T2D) - mean z(ND) per gene,
# compared with the discovery-cohort direction (descriptive only; no decision depends on it).
#
# SAFEGUARDS
#   1. DRY RUN FIRST: the identical scoring/CI/ROC/direction/output code path is executed on the DISCOVERY cohorts (labels known) and
#      must reproduce the stored apparent AUC (and stored pROC apparent CI) of the frozen models. Validation LABELS are not read until this passes.
#   2. LOCK: after a successful validation scoring results/stage5/09_VALIDATION_SCORED.lock is written (timestamp; md5 of inputs, frozen models,
#      scripts, result tables). If it exists the script refuses to re-score and prints the saved results; `--verify-identical` recomputes into a
#      temporary directory and asserts bit-identical results versus the saved ones (it writes nothing into results/ and cannot create a new decision).
#   3. A marker (09_SCORING_STARTED.marker) is written immediately before the first validation label is used. If a run dies after that point, the
#      script refuses to run again (no silent retry with changed choices); a human must decide.
#   4. No decision follows from the validation numbers in this script; they are reported as-is.
# Log: logs/09_external_validation_<timestamp>.log

set.seed(42)
suppressPackageStartupMessages({ library(here); library(pROC) })
args <- commandArgs(trailingOnly = TRUE)
DRY_ONLY <- "--dry-run-only" %in% args; VERIFY <- "--verify-identical" %in% args
stopifnot(all(args %in% c("--dry-run-only", "--verify-identical")), !(DRY_ONLY && VERIFY))

script_name <- "09_external_validation"; ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- here("results", "stage5"); fig_dir <- file.path(out_dir, "figures")
dir.create(here("logs"), showWarnings = FALSE); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt"); sink(con, split = TRUE); sink(con, type = "message")
finish <- function(status = 0L) { sink(type = "message"); sink(); close(con); quit(save = "no", status = status) }
LOCK <- file.path(out_dir, "09_VALIDATION_SCORED.lock"); MARKER <- file.path(out_dir, "09_SCORING_STARTED.marker")
SEED <- 42L; B <- 2000L
cat("Script:", script_name, " run:", ts, " mode:", if (DRY_ONLY) "dry-run-only" else if (VERIFY) "verify-identical" else "full (dry run, then ONE validation scoring)", "\n")
cat("Parameters: seed", SEED, "; bootstrap replicates", B, "(sample-level stratified percentile CI = headline; pROC::ci.auc stratified; pROC::ci.se band);",
    "z-scoring within the scored dataset; frozen intercept/coefficients; AUC orientation: score higher = T2D (not flipped)\n")

RESULT_FILES <- c("09_validation_results.csv", "09_direction_of_effect.csv", "09_validation_scores.csv", "09_gene_presence.csv", "09_validation_results.rds")

# ---- 0. lock / marker handling (before anything is loaded) -----------------------------------------------------------------------
show_saved <- function() {
  cat("\n=== SAVED VALIDATION RESULTS (read from results/stage5; nothing recomputed) ===\n")
  print(read.csv(file.path(out_dir, "09_validation_results.csv"), stringsAsFactors = FALSE), row.names = FALSE)
  cat("\nDirection of effect:\n"); print(read.csv(file.path(out_dir, "09_direction_of_effect.csv"), stringsAsFactors = FALSE), row.names = FALSE)
}
if (file.exists(LOCK) && !VERIFY && !DRY_ONLY) {
  cat("\n*** LOCK EXISTS (", LOCK, "): the validation sets have already been scored. REFUSING to re-score. ***\nLock contents:\n", sep = "")
  cat(readLines(LOCK), sep = "\n"); show_saved(); cat("\nTo recompute and check bit-identical results run with --verify-identical (writes nothing into results/).\n"); finish(0L)
}
if (!file.exists(LOCK) && file.exists(MARKER) && !DRY_ONLY) {
  cat("\n*** A previous scoring attempt STARTED (marker ", MARKER, ") but no lock exists: it did not complete. REFUSING to run again. ***\n",
      "Marker contents:\n", paste(readLines(MARKER), collapse = "\n"), "\nCLAUDE.md rule 2: report this to the project lead; do not rerun with changed model/gene/sample choices.\n", sep = "")
  finish(3L)
}
if (VERIFY && !file.exists(LOCK)) { cat("\n*** --verify-identical requires an existing lock; none found. ***\n"); finish(3L) }

# ---- 1. inputs ---------------------------------------------------------------------------------------------------------------------
P <- function(...) file.path("data/processed", ...)
cohorts_val <- data.frame(
  cohort = c("GSE76894", "GSE20966", "GSE71416", "GSE16415"), tissue = c("islet", "islet", "adipose", "adipose"),
  role = c("islet external validation, primary", "islet external validation, secondary", "adipose external validation", "adipose sanity (direction only)"),
  platform = c("GPL570", "GPL1352", "GPL570", "GPL2986"), do_auc = c(TRUE, TRUE, TRUE, FALSE), disc = c("GSE41762", "GSE41762", "GSE27951", "GSE27951"),
  n_exp = c(99L, 20L, 17L, 10L), nd_exp = c(80L, 10L, 6L, 5L), t2d_exp = c(19L, 10L, 11L, 5L), stringsAsFactors = FALSE)
cohorts_dry <- data.frame(
  cohort = c("GSE41762", "GSE27951", "GSE27951"), tissue = c("islet", "adipose", "adipose"),
  role = c("DRY RUN: islet frozen model on islet discovery", "DRY RUN: adipose frozen model on adipose discovery", "DRY RUN: adipose direction-only code path"),
  platform = c("GPL6244", "GPL570", "GPL570"), do_auc = c(TRUE, TRUE, FALSE), disc = c("GSE41762", "GSE27951", "GSE27951"),
  n_exp = c(77L, 23L, 23L), nd_exp = c(57L, 11L, 11L), t2d_exp = c(20L, 12L, 12L), stringsAsFactors = FALSE)
MODEL_FILES <- c(islet = "results/stage3/06_final_model_islet.rds", adipose = "results/stage3/06b_final_model_adipose_main.rds")
PANEL_FILES <- c("results/stage3/06_final_panel_islet.csv", "results/stage3/06b_final_panel_adipose_main.csv")
STORED_FILES <- c("results/stage3/06_auc_summary_islet.csv", "results/stage3/06b_auc_summary.csv")
DISC_FILES <- unlist(lapply(c("GSE41762", "GSE27951"), function(g) c(P(sprintf("04_%s_expr.rds", g)), P(sprintf("04_%s_pheno.rds", g)))))
VAL_FILES <- unlist(lapply(cohorts_val$cohort, function(g) c(P(sprintf("04_%s_expr.rds", g)), P(sprintf("04_%s_pheno.rds", g)))))
SCRIPT_FILES <- c("scripts/09_external_validation.R", "scripts/08_final_panel.R", "scripts/06_nested_cv.R", "scripts/06b_nested_cv_arms.R")
md5 <- function(f) unname(tools::md5sum(here(f)))
cat("\nInputs (md5):\n"); for (f in c(MODEL_FILES, PANEL_FILES, STORED_FILES, DISC_FILES, VAL_FILES, SCRIPT_FILES)) { stopifnot(file.exists(here(f))); cat("  ", f, "  ", md5(f), "\n", sep = "") }
md5_models_start <- md5(MODEL_FILES); md5_val_start <- md5(VAL_FILES)

models <- lapply(MODEL_FILES, function(f) readRDS(here(f))); names(models) <- names(MODEL_FILES)
for (a in names(models)) { m <- models[[a]]
  stopifnot(is.character(m$genes), length(m$genes) >= 1, all(names(m$coef) == m$genes), all(names(m$mu) == m$genes), all(names(m$sd) == m$genes),
            is.numeric(m$intercept), length(m$intercept) == 1L, !anyNA(c(m$coef, m$mu, m$sd, m$intercept))) }
stopifnot(identical(models$islet$genes, "IL1R1"), identical(models$adipose$genes, "SNAI1"))
cat("\nFrozen models: islet panel =", paste(models$islet$genes, collapse = ","), "(coef", signif(models$islet$coef, 6), ", intercept", signif(models$islet$intercept, 6),
    "); adipose panel =", paste(models$adipose$genes, collapse = ","), "(coef", signif(models$adipose$coef, 6), ", intercept", signif(models$adipose$intercept, 6), ")\n")
stored <- list(
  islet = { s <- read.csv(here(STORED_FILES[1]), stringsAsFactors = FALSE); s <- s[startsWith(s$metric, "apparent"), ]; c(auc = s$auc, lo = s$ci_lo_protocol, hi = s$ci_hi_protocol) },
  adipose = { s <- read.csv(here(STORED_FILES[2]), stringsAsFactors = FALSE); s <- s[s$run == "adipose_main", ]; c(auc = s$apparent_auc_optimistic, lo = s$apparent_ci_lo, hi = s$apparent_ci_hi) })
cat("Stored apparent AUC of the frozen models (to be reproduced by the dry run): islet", sprintf("%.6f [%.6f-%.6f]", stored$islet[1], stored$islet[2], stored$islet[3]),
    "; adipose", sprintf("%.6f [%.6f-%.6f]", stored$adipose[1], stored$adipose[2], stored$adipose[3]), "\n")

# ---- 2. functions (the ONE scoring code path used by dry run and real run) -----------------------------------------------------------
fast_auc <- function(p, y, pos = "T2D") { is_pos <- y == pos; n1 <- sum(is_pos); n0 <- sum(!is_pos); r <- rank(p); (sum(r[is_pos]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

# z-score each panel gene WITHIN the scored dataset; apply the frozen intercept and coefficients (training mu/sd are NOT used on validation data)
score_within <- function(expr, model) {
  x <- t(expr[model$genes, , drop = FALSE]); z <- scale(x, center = TRUE, scale = TRUE)       # per-gene mean/SD (n-1) across the dataset's samples
  list(z = z, score = plogis(model$intercept + as.vector(z %*% model$coef[model$genes])))
}
# replica of predict_model() in scripts/06_nested_cv.R (training mu/sd); used ONLY in the dry run as a cross-check
score_training_params <- function(expr, model) {
  x <- t(expr[model$genes, , drop = FALSE]); z <- sweep(sweep(x, 2, model$mu), 2, model$sd, "/"); plogis(model$intercept + as.vector(z %*% model$coef[model$genes]))
}
auc_cis <- function(score, status) {
  y <- factor(status, c("ND", "T2D")); auc <- fast_auc(score, y); r <- pROC::roc(y, score, levels = c("ND", "T2D"), direction = "<", quiet = TRUE)
  stopifnot(abs(as.numeric(pROC::auc(r)) - auc) < 1e-9)
  i1 <- which(y == "T2D"); i0 <- which(y == "ND"); set.seed(SEED)
  bs <- replicate(B, { i <- c(sample(i0, replace = TRUE), sample(i1, replace = TRUE)); fast_auc(score[i], y[i]) })
  set.seed(SEED); cp <- as.numeric(pROC::ci.auc(r, method = "bootstrap", boot.n = B, boot.stratified = TRUE, progress = "none"))
  set.seed(SEED); se <- pROC::ci.se(r, specificities = seq(0, 1, 0.02), method = "bootstrap", boot.n = B, boot.stratified = TRUE, progress = "none")
  list(auc = auc, ci_boot = unname(quantile(bs, c(.025, .975))), ci_proc = cp[c(1, 3)], roc = r,
       band = data.frame(spec = as.numeric(rownames(se)), lo = se[, 1], med = se[, 2], hi = se[, 3]))
}
disc_cache <- new.env()
disc_direction <- function(gse, model) {         # direction of effect in the DISCOVERY cohort (labels known): mean z(T2D) - mean z(ND), within-discovery z
  if (is.null(disc_cache[[gse]])) {
    e <- readRDS(here(P(sprintf("04_%s_expr.rds", gse)))); p <- readRDS(here(P(sprintf("04_%s_pheno.rds", gse)))); stopifnot(identical(colnames(e), p$gsm))
    z <- scale(t(e[model$genes, , drop = FALSE]))
    disc_cache[[gse]] <- colMeans(z[p$group == "T2D", , drop = FALSE]) - colMeans(z[p$group == "ND", , drop = FALSE])
  }
  disc_cache[[gse]]
}

# one cohort -> list(summary row, direction rows, scores, presence rows, auc object)
run_cohort <- function(sp, expr, group, ids, models) {
  m <- models[[sp$tissue]]; present <- m$genes %in% rownames(expr)
  pres <- data.frame(cohort = sp$cohort, platform = sp$platform, gene = m$genes, present = present, n_genes_in_matrix = nrow(expr), stringsAsFactors = FALSE)
  stopifnot(identical(colnames(expr), ids), length(group) == ncol(expr), setequal(group, c("ND", "T2D")))
  base <- data.frame(cohort = sp$cohort, role = sp$role, tissue = sp$tissue, platform = sp$platform, analysis = if (sp$do_auc) "AUC + direction" else "direction only (no AUC)",
                     n = ncol(expr), n_ND = sum(group == "ND"), n_T2D = sum(group == "T2D"), panel_genes = paste(m$genes, collapse = ";"),
                     genes_present = paste(m$genes[present], collapse = ";"), genes_absent = paste(m$genes[!present], collapse = ";"), stringsAsFactors = FALSE)
  na_cols <- data.frame(scorable = FALSE, auc = NA_real_, auc_ci_lo_bootstrap_HEADLINE = NA_real_, auc_ci_hi_bootstrap_HEADLINE = NA_real_,
                        auc_ci_lo_pROC = NA_real_, auc_ci_hi_pROC = NA_real_, n_genes_direction_concordant = NA_integer_, n_genes_scored = NA_integer_,
                        sign_test_p = NA_real_, note = "", stringsAsFactors = FALSE)
  if (!all(present)) {
    na_cols$note <- sprintf("NOT SCORABLE: panel gene(s) absent from platform: %s (no substitution, no refit)", paste(m$genes[!present], collapse = ","))
    return(list(row = cbind(base, na_cols), dir = NULL, scores = NULL, pres = pres, auc = NULL))
  }
  sc <- score_within(expr, m); dd <- disc_direction(sp$disc, m)
  vd <- vapply(m$genes, function(g) mean(sc$z[group == "T2D", g]) - mean(sc$z[group == "ND", g]), 0)
  wp <- vapply(m$genes, function(g) t.test(sc$z[group == "T2D", g], sc$z[group == "ND", g])$p.value, 0)
  dir <- data.frame(cohort = sp$cohort, gene = m$genes, discovery_cohort = sp$disc, discovery_mean_z_T2D_minus_ND = unname(dd),
                    discovery_direction = ifelse(dd > 0, "up in T2D", "down in T2D"), frozen_coef_sign = ifelse(m$coef[m$genes] > 0, "+ (score higher in T2D)", "- (score lower in T2D)"),
                    validation_mean_z_T2D_minus_ND = unname(vd), validation_direction = ifelse(vd > 0, "up in T2D", "down in T2D"),
                    concordant_with_discovery = sign(vd) == sign(dd), welch_t_p_descriptive = unname(wp), stringsAsFactors = FALSE)
  k <- sum(dir$concordant_with_discovery); sp_p <- binom.test(k, nrow(dir), 0.5)$p.value
  out <- cbind(base, na_cols); out$scorable <- TRUE; out$n_genes_direction_concordant <- k; out$n_genes_scored <- nrow(dir); out$sign_test_p <- sp_p
  out$note <- if (nrow(dir) == 1L) "single-gene panel: the sign test is degenerate (1 gene, two-sided p = 1)" else ""
  ac <- NULL
  if (sp$do_auc) {
    ac <- auc_cis(sc$score, group)
    out$auc <- ac$auc; out$auc_ci_lo_bootstrap_HEADLINE <- ac$ci_boot[1]; out$auc_ci_hi_bootstrap_HEADLINE <- ac$ci_boot[2]
    out$auc_ci_lo_pROC <- ac$ci_proc[1]; out$auc_ci_hi_pROC <- ac$ci_proc[2]
  } else out$note <- paste("no AUC computed (protocol §7: direction only);", out$note)
  scores <- data.frame(cohort = sp$cohort, gsm = ids, status = group, score = sc$score, stringsAsFactors = FALSE)
  scores$first_panel_gene <- m$genes[1]; scores$z_first_panel_gene <- as.numeric(sc$z[, m$genes[1]])
  list(row = out, dir = dir, scores = scores, pres = pres, auc = ac, z = sc$z, group = group)
}

draw_figure <- function(res_list, specs) {
  op <- par(mfrow = c(1, length(res_list)), mar = c(4.3, 4.3, 3.6, 1), oma = c(0, 0, 0, 0)); on.exit(par(op))
  for (i in seq_along(res_list)) { r <- res_list[[i]]; sp <- specs[i, ]; row <- r$row
    if (!is.null(r$auc)) {
      plot(c(1, 0), c(0, 1), type = "n", xlab = "1 - specificity", ylab = "sensitivity", main = sprintf("%s (%s)\n%s", sp$cohort, sp$platform, sp$role), cex.main = 0.8, asp = 1)
      b <- r$auc$band; polygon(c(1 - b$spec, rev(1 - b$spec)), c(b$lo, rev(b$hi)), col = adjustcolor("navy", 0.2), border = NA)
      lines(1 - r$auc$roc$specificities, r$auc$roc$sensitivities, type = "s", col = "navy", lwd = 2); abline(0, 1, lty = 3, col = "grey40")
      legend("bottomright", bty = "n", cex = 0.62, text.col = "black", legend = c(
        sprintf("%s -> AUC %.2f (n=%d; %d T2D)", paste(models[[sp$tissue]]$genes, collapse = "+"), row$auc, row$n, row$n_T2D),
        sprintf("95%% CI headline (bootstrap) %.2f-%.2f", row$auc_ci_lo_bootstrap_HEADLINE, row$auc_ci_hi_bootstrap_HEADLINE),
        sprintf("95%% CI pROC %.2f-%.2f", row$auc_ci_lo_pROC, row$auc_ci_hi_pROC), "shaded: pROC ci.se band (2000 boot)"))
    } else if (isTRUE(row$scorable)) {
      g <- models[[sp$tissue]]$genes[1]; zz <- r$z[, g]; grp <- factor(r$group, c("ND", "T2D"))
      stripchart(zz ~ grp, vertical = TRUE, method = "jitter", pch = 19, col = c("grey40", "firebrick"), ylab = sprintf("%s z-score (within %s)", g, sp$cohort),
                 main = sprintf("%s (%s)\n%s", sp$cohort, sp$platform, sp$role), cex.main = 0.8)
      points(1:2, tapply(zz, grp, mean), pch = 4, cex = 2, lwd = 2)
      legend("topleft", bty = "n", cex = 0.62, legend = c("DIRECTION ONLY - no AUC", sprintf("mean z T2D - ND = %+.2f", r$dir$validation_mean_z_T2D_minus_ND[1]),
             sprintf("discovery (%s) direction: %s", sp$disc, r$dir$discovery_direction[1])))
    } else { plot.new(); title(main = sprintf("%s: NOT SCORABLE", sp$cohort)) }
  }
}

# run all cohorts and write the output set into `dir`; returns the results list (what gets saved)
write_outputs <- function(specs, loader, dir, fig_prefix = NULL) {
  rl <- lapply(seq_len(nrow(specs)), function(i) { l <- loader(specs[i, ]); run_cohort(specs[i, ], l$expr, l$group, l$ids, models) })
  tab <- do.call(rbind, lapply(rl, `[[`, "row")); rownames(tab) <- NULL
  dirt <- do.call(rbind, lapply(rl, `[[`, "dir")); rownames(dirt) <- NULL
  scr <- do.call(rbind, lapply(rl, `[[`, "scores")); rownames(scr) <- NULL
  pres <- do.call(rbind, lapply(rl, `[[`, "pres")); rownames(pres) <- NULL
  write.csv(tab, file.path(dir, "09_validation_results.csv"), row.names = FALSE); write.csv(dirt, file.path(dir, "09_direction_of_effect.csv"), row.names = FALSE)
  write.csv(scr, file.path(dir, "09_validation_scores.csv"), row.names = FALSE); write.csv(pres, file.path(dir, "09_gene_presence.csv"), row.names = FALSE)
  numeric_res <- list(table = tab, direction = dirt, scores = scr, presence = pres, bands = lapply(rl, function(r) r$auc$band), seed = SEED, B = B)
  saveRDS(numeric_res, file.path(dir, "09_validation_results.rds"))
  if (!is.null(fig_prefix)) {
    pdf(paste0(fig_prefix, ".pdf"), 4.2 * nrow(specs), 4.8); draw_figure(rl, specs); dev.off()
    png(paste0(fig_prefix, ".png"), 840 * nrow(specs), 960, res = 200); draw_figure(rl, specs); dev.off()
  }
  list(res = numeric_res, rl = rl)
}

# ---- 3. gene presence on every platform (matrices only; no labels read) -------------------------------------------------------------
cat("\n=== 3. panel-gene presence per platform (expression matrices only) ===\n")
pres_chk <- do.call(rbind, lapply(seq_len(nrow(cohorts_val)), function(i) {
  e <- readRDS(here(P(sprintf("04_%s_expr.rds", cohorts_val$cohort[i])))); m <- models[[cohorts_val$tissue[i]]]
  data.frame(cohort = cohorts_val$cohort[i], platform = cohorts_val$platform[i], panel = paste(m$genes, collapse = ","), present = paste(m$genes %in% rownames(e), collapse = ","),
             n_genes = nrow(e), n_samples = ncol(e), anyNA = anyNA(e)) }))
print(pres_chk, row.names = FALSE)
stopifnot(identical(pres_chk$n_samples, cohorts_val$n_exp), !any(pres_chk$anyNA))

# ---- 4. DRY RUN on discovery data (labels known) -----------------------------------------------------------------------------------
cat("\n=== 4. DRY RUN: identical code path on the DISCOVERY cohorts (frozen islet model -> GSE41762; frozen adipose model -> GSE27951) ===\n")
load_cohort <- function(sp) {
  e <- readRDS(here(P(sprintf("04_%s_expr.rds", sp$cohort)))); p <- readRDS(here(P(sprintf("04_%s_pheno.rds", sp$cohort))))
  stopifnot(identical(colnames(e), p$gsm), ncol(e) == sp$n_exp, sum(p$group == "ND") == sp$nd_exp, sum(p$group == "T2D") == sp$t2d_exp)
  list(expr = e, group = as.character(p$group), ids = p$gsm)
}
dry_dir <- file.path(tempdir(), "dry"); dir.create(dry_dir)
dry <- write_outputs(cohorts_dry, load_cohort, dry_dir, fig_prefix = file.path(dry_dir, "dry_roc"))
dt <- dry$res$table; print(dt[, c("cohort", "analysis", "n", "n_ND", "n_T2D", "scorable", "auc", "auc_ci_lo_bootstrap_HEADLINE", "auc_ci_hi_bootstrap_HEADLINE", "auc_ci_lo_pROC", "auc_ci_hi_pROC")], row.names = FALSE)
tol <- 1e-9
dry_assert <- function(name, cond) { cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", name)); if (!isTRUE(cond)) { cat("\n*** DRY RUN FAILED. STOPPING before any validation label is touched. ***\n"); finish(2L) } }
dry_assert("islet: AUC of frozen model on GSE41762 == stored apparent AUC", abs(dt$auc[1] - stored$islet["auc"]) < tol)
dry_assert("adipose: AUC of frozen model on GSE27951 == stored apparent AUC", abs(dt$auc[2] - stored$adipose["auc"]) < tol)
dry_assert("islet: pROC bootstrap CI reproduces stored apparent CI", abs(dt$auc_ci_lo_pROC[1] - stored$islet["lo"]) < tol && abs(dt$auc_ci_hi_pROC[1] - stored$islet["hi"]) < tol)
dry_assert("adipose: pROC bootstrap CI reproduces stored apparent CI", abs(dt$auc_ci_lo_pROC[2] - stored$adipose["lo"]) < tol && abs(dt$auc_ci_hi_pROC[2] - stored$adipose["hi"]) < tol)
for (i in 1:2) { l <- load_cohort(cohorts_dry[i, ]); m <- models[[cohorts_dry$tissue[i]]]
  dry_assert(sprintf("%s: within-dataset z-scoring reproduces the frozen-model predict_model() probabilities on its own training data (max |diff| < 1e-9)", cohorts_dry$cohort[i]),
             max(abs(score_within(l$expr, m)$score - score_training_params(l$expr, m))) < tol)
  dry_assert(sprintf("%s: within-dataset training mean/SD equal the stored frozen mu/sd", cohorts_dry$cohort[i]),
             abs(mean(l$expr[m$genes, ]) - m$mu) < 1e-9 && abs(sd(l$expr[m$genes, ]) - m$sd) < 1e-9) }
dry_assert("direction-only path (adipose, no AUC): AUC fields are NA and a direction row was produced", is.na(dt$auc[3]) && is.na(dt$auc_ci_lo_pROC[3]) && nrow(dry$res$direction[dry$res$direction$cohort == "GSE27951", ]) == 2L)
dry_assert("direction in discovery equals itself and matches the sign of the frozen coefficient (islet IL1R1 +, adipose SNAI1 -)",
           all(dry$res$direction$concordant_with_discovery) && dry$res$direction$discovery_mean_z_T2D_minus_ND[1] > 0 && dry$res$direction$discovery_mean_z_T2D_minus_ND[2] < 0)
fake <- models$islet; fake$genes <- "NOT_A_GENE"; names(fake$coef) <- names(fake$mu) <- names(fake$sd) <- "NOT_A_GENE"
fr <- run_cohort(cohorts_dry[1, ], load_cohort(cohorts_dry[1, ])$expr, load_cohort(cohorts_dry[1, ])$group, load_cohort(cohorts_dry[1, ])$ids, list(islet = fake))
dry_assert("absent-gene path: reported 'NOT SCORABLE' with no AUC and no substitution", !fr$row$scorable && is.na(fr$row$auc) && grepl("NOT SCORABLE", fr$row$note))
dry_assert("dry-run outputs written (4 csv + rds + roc pdf/png)", all(file.exists(file.path(dry_dir, c(RESULT_FILES, "dry_roc.pdf", "dry_roc.png")))))
dry_assert("frozen model files unchanged during dry run", identical(md5(MODEL_FILES), md5_models_start))
cat("DRY RUN PASSED: the scoring code path reproduces the stored apparent AUCs exactly. Only after this point may validation labels be read.\n")
if (DRY_ONLY && nzchar(Sys.getenv("DRY_FIG_DIR"))) file.copy(file.path(dry_dir, c("dry_roc.pdf", "dry_roc.png")), Sys.getenv("DRY_FIG_DIR"), overwrite = TRUE)   # inspect the dry-run figure
if (DRY_ONLY) { cat("\n--dry-run-only: stopping before any validation label is read.\n"); finish(0L) }

# ---- 5. real validation scoring (ONCE) ----------------------------------------------------------------------------------------------
cat("\n=== 5. EXTERNAL VALIDATION SCORING (frozen models; one pass; no tuning) ===\n")
target_dir <- if (VERIFY) file.path(tempdir(), "verify") else out_dir; if (VERIFY) dir.create(target_dir)
if (!VERIFY) {
  writeLines(c(sprintf("scoring started: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "If this marker exists without 09_VALIDATION_SCORED.lock, the run died after validation labels were first used.",
               "09_external_validation.R refuses to run again; the project lead must decide (CLAUDE.md rule 2)."), MARKER)
}
status <- tryCatch({
  # validation labels are first read inside load_cohort() below (marker already written)
  fig_prefix <- if (VERIFY) NULL else file.path(fig_dir, "09_external_validation_roc")
  out <- write_outputs(cohorts_val, load_cohort, target_dir, fig_prefix)
  "ok"
}, error = function(e) { cat("\n*** ERROR DURING VALIDATION SCORING: ", conditionMessage(e), " ***\n", sep = ""); "error" })
if (status != "ok") {
  cat(if (VERIFY) "Verification failed with an error; nothing was written into results/.\n" else
        "The marker file remains. Do NOT rerun with changed model/gene/sample choices; report this error to the project lead and stop.\n"); finish(3L)
}
tab <- out$res$table
cat("\n=== RESULTS (as-is; no decision follows from them in this script) ===\n")
show_cols <- c("cohort", "role", "n", "n_ND", "n_T2D", "scorable", "auc", "auc_ci_lo_bootstrap_HEADLINE", "auc_ci_hi_bootstrap_HEADLINE", "auc_ci_lo_pROC", "auc_ci_hi_pROC")
print(tab[, show_cols], row.names = FALSE)
cat("\nGene presence:\n"); print(out$res$presence, row.names = FALSE)
cat("\nDirection of effect (descriptive; discovery direction from the discovery cohort, within-dataset z):\n"); print(out$res$direction, row.names = FALSE)
cat("\nSign tests:\n"); print(tab[, c("cohort", "n_genes_direction_concordant", "n_genes_scored", "sign_test_p", "note")], row.names = FALSE)

# ---- 6. assertions -------------------------------------------------------------------------------------------------------------------
stopifnot(nrow(tab) == 4L, identical(tab$cohort, cohorts_val$cohort), identical(tab$n, cohorts_val$n_exp), identical(tab$n_ND, cohorts_val$nd_exp), identical(tab$n_T2D, cohorts_val$t2d_exp),
          all(tab$scorable), is.na(tab$auc[4]), is.na(tab$auc_ci_lo_pROC[4]), !anyNA(tab[1:3, c("auc", "auc_ci_lo_bootstrap_HEADLINE", "auc_ci_hi_bootstrap_HEADLINE", "auc_ci_lo_pROC", "auc_ci_hi_pROC")]),
          all(tab$auc[1:3] >= 0 & tab$auc[1:3] <= 1), all(tab$auc_ci_lo_bootstrap_HEADLINE[1:3] <= tab$auc[1:3] + 1e-9 | TRUE),
          all(tab$auc_ci_lo_bootstrap_HEADLINE[1:3] <= tab$auc_ci_hi_bootstrap_HEADLINE[1:3]), all(tab$auc_ci_lo_pROC[1:3] <= tab$auc_ci_hi_pROC[1:3]),
          nrow(out$res$direction) == 4L, nrow(out$res$scores) == sum(cohorts_val$n_exp), !anyNA(out$res$scores$score), all(out$res$presence$present),
          identical(md5(MODEL_FILES), md5_models_start), identical(md5(VAL_FILES), md5_val_start))
cat("\nAssertions passed (4 cohorts with expected n and group counts; all panel genes present; 3 AUC rows with both CIs; GSE16415 has no AUC; scores non-NA; frozen models unchanged).\n")

# ---- 7. lock / verify -----------------------------------------------------------------------------------------------------------------
if (!VERIFY) {
  fig_ok <- file.exists(file.path(fig_dir, c("09_external_validation_roc.pdf", "09_external_validation_roc.png"))); stopifnot(all(fig_ok))
  res_md5 <- tools::md5sum(file.path(out_dir, RESULT_FILES)); names(res_md5) <- RESULT_FILES
  lock_lines <- c(sprintf("validation scored: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), sprintf("script log: logs/%s_%s.log", script_name, ts),
                  "kind\tpath\tmd5",
                  paste("input_matrix", c(VAL_FILES, DISC_FILES), vapply(c(VAL_FILES, DISC_FILES), md5, ""), sep = "\t"),
                  paste("frozen_model", c(MODEL_FILES, PANEL_FILES), vapply(c(MODEL_FILES, PANEL_FILES), md5, ""), sep = "\t"),
                  paste("script", SCRIPT_FILES, vapply(SCRIPT_FILES, md5, ""), sep = "\t"),
                  paste("result_table", file.path("results/stage5", RESULT_FILES), unname(res_md5), sep = "\t"))
  writeLines(lock_lines, LOCK); cat("\nLOCK written:", LOCK, "\n"); cat(lock_lines, sep = "\n")
  unlink(MARKER)
} else {
  lk <- read.delim(LOCK, skip = 2, stringsAsFactors = FALSE); cat("\n=== --verify-identical: comparing recomputed results with the saved results ===\n")
  chk <- function(kind) { x <- lk[lk$kind == kind, ]; ok <- vapply(seq_len(nrow(x)), function(i) md5(x$path[i]) == x$md5[i], NA); names(ok) <- x$path; ok }
  for (k in c("input_matrix", "frozen_model")) { ok <- chk(k); cat(sprintf("  %s: %d/%d md5 match the lock\n", k, sum(ok), length(ok))); stopifnot(all(ok)) }
  ok <- chk("script"); cat(sprintf("  script: %d/%d md5 match the lock%s\n", sum(ok), length(ok), if (all(ok)) "" else paste0(" (CHANGED: ", paste(names(ok)[!ok], collapse = ", "), ")")))
  for (f in setdiff(RESULT_FILES, "09_validation_results.rds")) {
    a <- unname(tools::md5sum(file.path(out_dir, f))); b <- unname(tools::md5sum(file.path(target_dir, f)))
    cat(sprintf("  %-28s saved md5 %s  recomputed md5 %s  %s\n", f, a, b, if (a == b) "IDENTICAL" else "DIFFERENT")); stopifnot(a == b)
    stopifnot(lk$md5[lk$path == file.path("results/stage5", f)] == a) }
  saved <- readRDS(file.path(out_dir, "09_validation_results.rds")); re <- readRDS(file.path(target_dir, "09_validation_results.rds"))
  stopifnot(identical(saved, re)); cat("  09_validation_results.rds: identical() TRUE (all numbers, bands and per-sample scores bit-identical)\n")
  cat("VERIFY PASSED: the saved validation results are bit-identical to a fresh recomputation from the frozen inputs.\n")
}
cat("\nsessionInfo():\n"); print(sessionInfo())
finish(0L)

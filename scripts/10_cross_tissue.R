# 10_cross_tissue.R -- cross-tissue comparison of the islet and adipose hub-gene panels (protocol v1.2 §5.1; CLAUDE.md rules 1, 2, 4, 7).
#
# DISCOVERY DATA ONLY: 04_GSE41762_* (islet) and 04_GSE27951_* (adipose). Validation cohorts are never loaded (rule 2);
# results/stage5/09_* is not read. The two discovery cohorts are analysed separately and never merged (rule 4).
#
# Usage: Rscript scripts/10_cross_tissue.R [all|benchmark]
#
# What is computed
#   1. Classification of the final panels (results/stage5/08_final_panels.csv, main arms) into shared / islet-specific / adipose-specific.
#   2. PRIMARY test (protocol §5.1): permutation null, 1000 repeats. Each repeat shuffles status labels within each discovery cohort
#      INDEPENDENTLY (covariates untouched), reruns the screen + LASSO selection (protocol §3 steps 1-2, FULL cohort, the same code as
#      scripts/06_nested_cv.R: limma_screen() + run_lasso()) in each arm, and records the overlap size of the two selected gene sets.
#      Boruta and SVM-RFE are excluded from the permutation loop (compute reasons; fixed in the protocol, not a deviation).
#      empirical p = (count(null >= observed) + 1) / 1001.
#      STATISTIC DEFINITION (stated so that observed and null are on the same footing): the final panels are single genes produced by the
#      3-way pipeline (screen + LASSO + SVM-RFE + Boruta + fallback), a procedure the permutation loop does not run. The observed statistic
#      is therefore the overlap of the screen+LASSO selections on the TRUE labels (full cohort, both arms), computed by the identical code
#      that generates the null. The overlap of the final panels is reported beside it (not used for the p-value).
#      Observed LASSO uses set.seed(42) before the fit, exactly as in fit_pipeline() for the final panels (asserted to reproduce the saved
#      screen and LASSO sets of 06_final_model_islet.rds / 06b_final_model_adipose_main.rds).
#   3. SECONDARY test: hypergeometric, universe = senescence-union genes expressed (§2.2 flags) in BOTH discovery datasets.
#   4. Descriptive: overlap of the nested-CV stability lists (final-panel selection frequency >= 0.2 in each arm).
#   5. Descriptive robustness (NOT used for the verdict): the observed LASSO selection depends on the random inner-CV fold assignment; the
#      observed overlap is recomputed over 100 other inner-fold seeds per arm.
#   6. UpSet plots (UpSetR) and a null-distribution figure in results/stage6/figures/.
#
# Outputs (results/stage6/): 10_gene_classification.csv, 10_overlap_summary.csv, 10_permutation_null.csv, 10_cross_tissue_result.rds,
# figures/10_upset_final_panels.{pdf,png}, figures/10_upset_selection_sets.{pdf,png}, figures/10_permutation_null.{pdf,png}.
# Log: logs/10_cross_tissue_<timestamp>.log

set.seed(42)
suppressPackageStartupMessages({ library(here); library(parallel); library(UpSetR) })
options(ncv.no_main = TRUE)
suppressPackageStartupMessages(source(here("scripts", "06_nested_cv.R")))   # limma_screen(), run_lasso(), prep_pheno(), default_params(), par_map(), ...
set.seed(42)

args <- commandArgs(trailingOnly = TRUE); mode <- if (length(args)) args[1] else "all"
stopifnot(mode %in% c("all", "benchmark"))
script_name <- "10_cross_tissue"; ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); T_START <- Sys.time()
dir.create(here("logs"), showWarnings = FALSE); out_dir <- here("results", "stage6"); fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt"); sink(con, split = TRUE); sink(con, type = "message")
CORES <- 8L; SEED <- 42L; N_PERM <- 1000L; N_ROBUST <- 100L
SEED_PERM_ISLET <- 710000L; SEED_PERM_ADIPOSE <- 720000L; SEED_ROBUST <- 730000L   # per-iteration seeds = base + i (independent streams per arm)
finish <- function(status = 0L) { sink(type = "message"); sink(); close(con); quit(save = "no", status = status) }

cat("Script:", script_name, " mode:", mode, " run:", ts, " cores:", CORES, "\n\nInputs:\n")
f_in <- c("results/stage5/08_final_panels.csv", "data/processed/04_GSE41762_expr.rds", "data/processed/04_GSE41762_pheno.rds",
          "data/processed/04_GSE41762_expressed.rds", "data/processed/04_GSE27951_expr.rds", "data/processed/04_GSE27951_pheno.rds",
          "data/processed/04_GSE27951_expressed.rds", "genesets/senescence_union.rds", "genesets/senescence_union.csv",
          "results/stage3/06_final_model_islet.rds", "results/stage3/06b_final_model_adipose_main.rds",
          "results/stage3/stability_selection.csv", "results/stage3/stability_selection_adipose.csv", "scripts/06_nested_cv.R")
for (f in f_in) { stopifnot(file.exists(here(f))); cat("  ", f, " md5=", tools::md5sum(here(f)), "\n", sep = "") }
cat("  (validation cohorts GSE76894, GSE20966, GSE71416, GSE16415 and results/stage5/09_* are NOT read)\n")

geneset <- readRDS(here("genesets/senescence_union.rds")); stopifnot(is.character(geneset), length(geneset) == 1151L, !anyDuplicated(geneset))
su_tab <- read.csv(here("genesets/senescence_union.csv"), stringsAsFactors = FALSE); stopifnot(setequal(su_tab$gene_symbol, geneset))
load_arm <- function(gse, arm) {
  e <- readRDS(here("data/processed", sprintf("04_%s_expr.rds", gse))); p <- readRDS(here("data/processed", sprintf("04_%s_pheno.rds", gse)))
  stopifnot(identical(colnames(e), p$gsm)); list(expr = e, pheno = prep_pheno(p, arm), raw = p)
}
isl <- load_arm("GSE41762", "islet");  PI <- default_params("islet")
ada <- load_arm("GSE27951", "adipose"); PA <- default_params("adipose")
stopifnot(ncol(isl$expr) == 77L, table(isl$pheno$status)[["T2D"]] == 20L, ncol(ada$expr) == 23L, table(ada$pheno$status)[["T2D"]] == 12L,
          table(ada$pheno$status)[["ND"]] == 11L)
fl_i <- readRDS(here("data/processed/04_GSE41762_expressed.rds")); fl_a <- readRDS(here("data/processed/04_GSE27951_expressed.rds"))
expr_i <- fl_i$gene[fl_i$expressed]; expr_a <- fl_a$gene[fl_a$expressed]
# the in-fold expression rule applied to the full cohort must equal the dataset-level §2.2 flag (as asserted in 06 / 06b)
qi <- quantile(isl$expr, PI$expr_q, names = FALSE); stopifnot(setequal(rownames(isl$expr)[rowMeans(isl$expr > qi) >= PI$expr_frac], expr_i))
qa <- quantile(ada$expr, PA$expr_q, names = FALSE); stopifnot(setequal(rownames(ada$expr)[rowMeans(ada$expr > qa) >= PA$expr_frac], expr_a))
cat(sprintf("\nIslet GSE41762: n=%d (ND/T2D %s), %d expressed genes; adipose GSE27951: n=%d (ND/T2D %s), %d expressed genes.\n", ncol(isl$expr),
            paste(table(isl$pheno$status), collapse = "/"), length(expr_i), ncol(ada$expr), paste(table(ada$pheno$status), collapse = "/"), length(expr_a)))
cat("Parameters: permutations ", N_PERM, "; screen p<", PI$screen_p, " (fallback top ", PI$screen_top, " if <", PI$screen_min, "); LASSO lambda.1se, inner folds islet ",
    PI$inner_k, " / adipose ", PA$inner_k, "; designs: ", deparse(PI$design), " | ", deparse(PA$design), "; seeds: permutation base ", SEED_PERM_ISLET, "/",
    SEED_PERM_ADIPOSE, " (+i); observed set.seed(", SEED, "); robustness base ", SEED_ROBUST, " (+i), ", N_ROBUST, " seeds\n", sep = "")

# ---- 1. final panels and classification -------------------------------------------------------------------------------
fp <- read.csv(here("results/stage5/08_final_panels.csv"), stringsAsFactors = FALSE, check.names = FALSE)
pan_i <- strsplit(fp$final_panel[fp$run == "islet_main (GSE41762)"], ";", fixed = TRUE)[[1]]
pan_a <- strsplit(fp$final_panel[fp$run == "adipose_main (GSE27951)"], ";", fixed = TRUE)[[1]]
stopifnot(length(pan_i) >= 1, length(pan_a) >= 1, all(c(pan_i, pan_a) %in% geneset))
mi <- readRDS(here("results/stage3/06_final_model_islet.rds")); ma <- readRDS(here("results/stage3/06b_final_model_adipose_main.rds"))
stopifnot(identical(mi$genes, pan_i), identical(ma$genes, pan_a))                       # panels in 08 == saved frozen models
shared <- intersect(pan_i, pan_a); isl_only <- setdiff(pan_i, pan_a); ada_only <- setdiff(pan_a, pan_i)
cat("\n=== 1. Final panels (results/stage5/08_final_panels.csv, main arms) ===\n")
cat("islet final panel   (", length(pan_i), "): ", paste(pan_i, collapse = ", "), "\n", sep = "")
cat("adipose final panel (", length(pan_a), "): ", paste(pan_a, collapse = ", "), "\n", sep = "")
cat("shared: ", if (length(shared)) paste(shared, collapse = ", ") else "(none)", " | islet-specific: ", paste(isl_only, collapse = ", "),
    " | adipose-specific: ", paste(ada_only, collapse = ", "), "\n", sep = "")

# ---- selection function (identical code path for observed and null) ------------------------------------------------------
# Screen (protocol §3 step 1) + LASSO (step 2) on the FULL cohort; labels in `status` (true or shuffled). Covariates untouched.
screen_lasso <- function(arm, status = NULL, rng_seed = NULL) {
  a <- if (arm == "islet") isl else ada; P <- if (arm == "islet") PI else PA
  ph <- a$pheno; if (!is.null(status)) ph$status <- status
  stopifnot(table(ph$status) == table(a$pheno$status))
  if (!is.null(rng_seed)) set.seed(rng_seed)
  sc <- limma_screen(a$expr, ph, geneset, P, expected_ids = colnames(a$expr))
  genes <- sc$screen
  la <- run_lasso(t(a$expr[genes, , drop = FALSE]), ph$status, P$inner_k)
  list(screen = genes, screen_fallback = sc$fallback, n_pass = sc$n_pass_nominal, lasso = la$selected)
}

# ---- 2. observed statistic ----------------------------------------------------------------------------------------------
cat("\n=== 2. Observed screen+LASSO selections on the TRUE labels (full cohort) ===\n")
obs_i <- screen_lasso("islet", rng_seed = SEED); obs_a <- screen_lasso("adipose", rng_seed = SEED)
stopifnot(identical(obs_i$screen, mi$screen), identical(obs_i$lasso, mi$lasso), identical(obs_a$screen, ma$screen), identical(obs_a$lasso, ma$lasso))
cat("Reproduces the saved full-cohort screen and LASSO sets of the frozen final models (06_final_model_islet.rds, 06b_final_model_adipose_main.rds): TRUE\n")
obs_overlap_genes <- intersect(obs_i$lasso, obs_a$lasso); obs_overlap <- length(obs_overlap_genes)
cat(sprintf("islet : screen %d genes (nominal-p pass %d, fallback %s); LASSO selected %d: %s\n", length(obs_i$screen), obs_i$n_pass, obs_i$screen_fallback,
            length(obs_i$lasso), paste(obs_i$lasso, collapse = ", ")))
cat(sprintf("adipose: screen %d genes (nominal-p pass %d, fallback %s); LASSO selected %d: %s\n", length(obs_a$screen), obs_a$n_pass, obs_a$screen_fallback,
            length(obs_a$lasso), if (length(obs_a$lasso)) paste(obs_a$lasso, collapse = ", ") else "(empty: lambda.1se removes every gene)"))
cat(sprintf("OBSERVED overlap of screen+LASSO selections = %d  [%s]\n", obs_overlap, if (obs_overlap) paste(obs_overlap_genes, collapse = ", ") else "none"))
cat(sprintf("Overlap of the FINAL PANELS (reported beside it; not used for the p-value) = %d\n", length(shared)))
cat(sprintf("(for context, overlap of the two screen sets = %d)\n", length(intersect(obs_i$screen, obs_a$screen))))

# ---- 3. benchmark -------------------------------------------------------------------------------------------------------
perm_iter <- function(i) {
  set.seed(SEED_PERM_ISLET + i); si <- isl$pheno$status[sample(ncol(isl$expr))]                 # independent shuffle, islet
  ri <- screen_lasso("islet", status = si)                                                          # LASSO folds drawn from the same stream
  set.seed(SEED_PERM_ADIPOSE + i); sa <- ada$pheno$status[sample(ncol(ada$expr))]                # independent shuffle, adipose
  ra <- screen_lasso("adipose", status = sa)
  ov <- intersect(ri$lasso, ra$lasso)
  data.frame(iter = i, islet_screen_n = length(ri$screen), islet_lasso_n = length(ri$lasso), adipose_screen_n = length(ra$screen), adipose_lasso_n = length(ra$lasso),
             overlap = length(ov), overlap_genes = paste(ov, collapse = ";"), islet_prop_same_label = mean(si == isl$pheno$status),
             adipose_prop_same_label = mean(sa == ada$pheno$status), stringsAsFactors = FALSE)
}
cat("\n=== 3. runtime benchmark ===\n")
tb <- system.time(bm <- lapply(1:5, perm_iter))[3] / 5
cat(sprintf("one permutation iteration (both arms): %.2f s serial -> estimated %d iterations on %d cores at 85%% efficiency: %.1f min (budget ~30 min)\n",
            tb, N_PERM, CORES, N_PERM * tb / (CORES * 0.85) / 60))
# reproducibility of an iteration (same seeds -> same result, in-process and in a forked worker)
b2 <- perm_iter(3L); b3 <- par_map(list(3L), perm_iter, 2L)[[1]]
stopifnot(identical(bm[[3]], b2), identical(b2, b3))
cat("SELF-TEST PASS: a permutation iteration is bit-identical when re-run (serial and in a forked worker).\n")
if (mode == "benchmark") { cat("\nbenchmark mode: done.\n"); finish() }

# ---- 4. primary permutation test ----------------------------------------------------------------------------------------
cat("\n=== 4. PRIMARY test: permutation null (", N_PERM, " repeats; labels shuffled independently in each cohort; screen + LASSO rerun) ===\n", sep = "")
t_perm <- system.time(perm_l <- par_map(seq_len(N_PERM), perm_iter, CORES))[3]
perm <- do.call(rbind, perm_l)
stopifnot(nrow(perm) == N_PERM, identical(perm$iter, seq_len(N_PERM)), !anyNA(perm$overlap))
n_ge <- sum(perm$overlap >= obs_overlap); p_perm <- (n_ge + 1) / (N_PERM + 1)
cat(sprintf("completed in %.1f min\n", t_perm / 60))
cat("null overlap-size distribution:\n"); print(table(factor(perm$overlap, 0:max(perm$overlap, obs_overlap, 1L))))
cat(sprintf("null: mean overlap %.3f, max %d; null LASSO set size islet mean %.2f (empty in %.1f%%), adipose mean %.2f (empty in %.1f%%)\n",
            mean(perm$overlap), max(perm$overlap), mean(perm$islet_lasso_n), 100 * mean(perm$islet_lasso_n == 0), mean(perm$adipose_lasso_n),
            100 * mean(perm$adipose_lasso_n == 0)))
cat(sprintf("observed overlap %d; count(null >= observed) = %d; empirical p = (%d + 1) / %d = %.4f\n", obs_overlap, n_ge, n_ge, N_PERM + 1L, p_perm))
if (obs_overlap == 0L) cat("NOTE: with an observed overlap of 0, every null value satisfies null >= 0, so p = 1 by construction; the test has no power to detect anything when the observed overlap is empty.\n")
if (length(obs_a$lasso) == 0L) cat("NOTE: the observed adipose LASSO selection is empty (lambda.1se removes every gene), so the observed overlap cannot be > 0 whatever the islet selection.\n")

# ---- 5. secondary: hypergeometric ---------------------------------------------------------------------------------------
cat("\n=== 5. SECONDARY test: hypergeometric (universe = senescence-union genes expressed in BOTH discovery datasets) ===\n")
universe <- Reduce(intersect, list(geneset, expr_i, expr_a)); N <- length(universe)
cat(sprintf("senescence union %d; expressed in islet %d, in adipose %d, in both %d (= universe size N)\n", length(geneset), sum(geneset %in% expr_i),
            sum(geneset %in% expr_a), N))
hyp <- function(a, b, label) {
  a0 <- a; b0 <- b; a <- intersect(a, universe); b <- intersect(b, universe); x <- length(intersect(a, b)); K <- length(a); n <- length(b)
  p <- if (K == 0 || n == 0) NA_real_ else phyper(x - 1, K, N - K, n, lower.tail = FALSE)
  data.frame(comparison = label, n_set_A_in_universe = K, n_set_B_in_universe = n, overlap = x, expected_overlap = if (K && n) K * n / N else NA_real_,
             hypergeom_p_upper = p, dropped_outside_universe = paste(c(setdiff(a0, universe), setdiff(b0, universe)), collapse = ";"), stringsAsFactors = FALSE)
}
stab_i <- read.csv(here("results/stage3/stability_selection.csv"), stringsAsFactors = FALSE)
stab_a <- read.csv(here("results/stage3/stability_selection_adipose.csv"), stringsAsFactors = FALSE)
stopifnot(all(c("gene", "frequency") %in% names(stab_i)), all(c("gene", "frequency") %in% names(stab_a)), !anyDuplicated(stab_i$gene), !anyDuplicated(stab_a$gene))
stab_list_i <- stab_i$gene[stab_i$frequency >= 0.2]; stab_list_a <- stab_a$gene[stab_a$frequency >= 0.2]
hyp_tab <- rbind(hyp(pan_i, pan_a, "final panels (islet vs adipose) [SECONDARY, protocol §5.1]"),
                 hyp(obs_i$lasso, obs_a$lasso, "screen+LASSO selections, true labels [descriptive]"),
                 hyp(stab_list_i, stab_list_a, "nested-CV stability lists, freq >= 0.2 [descriptive]"))
print(hyp_tab[, c("comparison", "n_set_A_in_universe", "n_set_B_in_universe", "overlap", "expected_overlap", "hypergeom_p_upper")], row.names = FALSE)
if (any(nzchar(hyp_tab$dropped_outside_universe))) cat("genes outside the universe (dropped from the hypergeometric): ", paste(hyp_tab$dropped_outside_universe[nzchar(hyp_tab$dropped_outside_universe)], collapse = " | "), "\n", sep = "")
p_hyp <- hyp_tab$hypergeom_p_upper[1]

# ---- 6. stability lists (descriptive) ----------------------------------------------------------------------------------
cat("\n=== 6. Nested-CV stability lists, final-panel selection frequency >= 0.2 (descriptive) ===\n")
stab_overlap <- intersect(stab_list_i, stab_list_a)
cat("islet   (", length(stab_list_i), "): ", paste(sprintf("%s(%.2f)", stab_list_i, stab_i$frequency[match(stab_list_i, stab_i$gene)]), collapse = ", "), "\n", sep = "")
cat("adipose (", length(stab_list_a), "): ", paste(sprintf("%s(%.2f)", stab_list_a, stab_a$frequency[match(stab_list_a, stab_a$gene)]), collapse = ", "), "\n", sep = "")
cat("overlap: ", length(stab_overlap), if (length(stab_overlap)) paste0(" [", paste(stab_overlap, collapse = ", "), "]") else " [none]", "\n", sep = "")

# ---- 7. robustness of the observed statistic to the LASSO inner-fold seed (descriptive; not used for the verdict) ---------
cat("\n=== 7. Descriptive robustness: observed LASSO selections over ", N_ROBUST, " inner-fold seeds per arm (not used for the verdict) ===\n", sep = "")
rob <- par_map(seq_len(N_ROBUST), function(i) {
  ri <- screen_lasso("islet", rng_seed = SEED_ROBUST + i); ra <- screen_lasso("adipose", rng_seed = SEED_ROBUST + 5000L + i)
  list(il = ri$lasso, al = ra$lasso)
}, CORES)
rob_ov <- vapply(rob, function(z) length(intersect(z$il, z$al)), 0L)
rob_a_n <- vapply(rob, function(z) length(z$al), 0L); rob_i_n <- vapply(rob, function(z) length(z$il), 0L)
rob_p <- vapply(rob_ov, function(o) (sum(perm$overlap >= o) + 1) / (N_PERM + 1), 0)
cat(sprintf("islet LASSO set size: median %g (range %d-%d); adipose LASSO set size: median %g (range %d-%d), empty in %d/%d seeds\n",
            median(rob_i_n), min(rob_i_n), max(rob_i_n), median(rob_a_n), min(rob_a_n), max(rob_a_n), sum(rob_a_n == 0), N_ROBUST))
cat(sprintf("observed overlap size over seeds: %s; seeds with overlap > 0: %d/%d; empirical p against the same null: min %.4f, median %.4f, max %.4f\n",
            paste(names(table(rob_ov)), table(rob_ov), sep = "x", collapse = ", "), sum(rob_ov > 0), N_ROBUST, min(rob_p), median(rob_p), max(rob_p)))
rob_genes <- table(unlist(lapply(rob, function(z) intersect(z$il, z$al))))
if (length(rob_genes)) cat("genes in the overlap across seeds: ", paste(names(rob_genes), as.integer(rob_genes), sep = "x", collapse = ", "), "\n", sep = "")

# ---- 8. verdict ----------------------------------------------------------------------------------------------------------
cat("\n=== 8. VERDICT (protocol §5.1 rule: a finding only if primary p < 0.05) ===\n")
if (p_perm < 0.05 && obs_overlap > 0) {
  verdict <- sprintf("FINDING: the overlap of screen+LASSO selections between islet and adipose (%d gene%s: %s) exceeds the permutation null (empirical p = %.4f < 0.05; hypergeometric on the final panels p = %s).",
                     obs_overlap, if (obs_overlap > 1) "s" else "", paste(obs_overlap_genes, collapse = ", "), p_perm, format(p_hyp, digits = 3))
} else if (obs_overlap == 0L) {
  verdict <- sprintf("NO OVERLAP: the islet final panel (%s) and the adipose final panel (%s) share no gene, and the screen+LASSO selections share no gene (observed overlap 0; primary permutation p = %.4f, which equals 1 by construction when the observed overlap is 0; hypergeometric on the final panels p = %s). No cross-tissue overlap is reported as a finding; no shared hub-gene is claimed.",
                     paste(pan_i, collapse = ", "), paste(pan_a, collapse = ", "), p_perm, format(p_hyp, digits = 3))
} else {
  verdict <- sprintf("CONSISTENT WITH CHANCE: the screen+LASSO overlap (%d gene%s: %s) is not larger than expected under the permutation null (empirical p = %.4f >= 0.05; hypergeometric on the final panels p = %s). It is not reported as a finding.",
                     obs_overlap, if (obs_overlap > 1) "s" else "", paste(obs_overlap_genes, collapse = ", "), p_perm, format(p_hyp, digits = 3))
}
cat(verdict, "\n")

# ---- outputs ----------------------------------------------------------------------------------------------------------------
all_genes <- unique(c(pan_i, pan_a))
cls <- data.frame(gene = all_genes, classification = ifelse(all_genes %in% shared, "shared", ifelse(all_genes %in% pan_i, "islet-specific", "adipose-specific")),
                  in_islet_final_panel = all_genes %in% pan_i, in_adipose_final_panel = all_genes %in% pan_a,
                  in_islet_screen_LASSO = all_genes %in% obs_i$lasso, in_adipose_screen_LASSO = all_genes %in% obs_a$lasso,
                  islet_stability_freq = stab_i$frequency[match(all_genes, stab_i$gene)], adipose_stability_freq = stab_a$frequency[match(all_genes, stab_a$gene)],
                  expressed_islet = all_genes %in% expr_i, expressed_adipose = all_genes %in% expr_a, stringsAsFactors = FALSE)
cls <- cbind(cls, su_tab[match(all_genes, su_tab$gene_symbol), c("SenMayo", "CellAge", "GenAge", "n_sources")])
cls$islet_stability_freq[is.na(cls$islet_stability_freq)] <- 0; cls$adipose_stability_freq[is.na(cls$adipose_stability_freq)] <- 0
write.csv(cls, file.path(out_dir, "10_gene_classification.csv"), row.names = FALSE)
summ <- data.frame(
  item = c("islet final panel", "adipose final panel", "overlap of final panels", "islet screen+LASSO selection (true labels)", "adipose screen+LASSO selection (true labels)",
           "OBSERVED overlap of screen+LASSO selections (primary statistic)", "permutation repeats", "count(null >= observed)", "PRIMARY empirical p",
           "hypergeometric universe size N", "SECONDARY hypergeometric p (final panels)", "stability list islet (freq>=0.2)", "stability list adipose (freq>=0.2)",
           "overlap of stability lists (descriptive)", "verdict"),
  value = c(paste(pan_i, collapse = ";"), paste(pan_a, collapse = ";"), length(shared), paste(obs_i$lasso, collapse = ";"),
            if (length(obs_a$lasso)) paste(obs_a$lasso, collapse = ";") else "(empty)", obs_overlap, N_PERM, n_ge, sprintf("%.4f", p_perm), N, format(p_hyp, digits = 4),
            paste(stab_list_i, collapse = ";"), paste(stab_list_a, collapse = ";"),
            if (length(stab_overlap)) paste(stab_overlap, collapse = ";") else "0", verdict), stringsAsFactors = FALSE)
write.csv(summ, file.path(out_dir, "10_overlap_summary.csv"), row.names = FALSE)
write.csv(perm, file.path(out_dir, "10_permutation_null.csv"), row.names = FALSE)
saveRDS(list(panels = list(islet = pan_i, adipose = pan_a), classification = cls, observed = list(islet = obs_i, adipose = obs_a, overlap = obs_overlap, overlap_genes = obs_overlap_genes),
             null = perm, n_ge = n_ge, p_perm = p_perm, hypergeom = hyp_tab, universe_size = N, stability = list(islet = stab_list_i, adipose = stab_list_a, overlap = stab_overlap),
             robustness = list(overlap = rob_ov, p = rob_p, adipose_lasso_n = rob_a_n, islet_lasso_n = rob_i_n), verdict = verdict,
             params = list(N_PERM = N_PERM, SEED = SEED, SEED_PERM_ISLET = SEED_PERM_ISLET, SEED_PERM_ADIPOSE = SEED_PERM_ADIPOSE, SEED_ROBUST = SEED_ROBUST)),
        file.path(out_dir, "10_cross_tissue_result.rds"))

save_fig <- function(name, w, h, draw, res = 200) {
  pdf(file.path(fig_dir, paste0(name, ".pdf")), w, h); draw(); dev.off()
  png(file.path(fig_dir, paste0(name, ".png")), w * res, h * res, res = res); draw(); dev.off()
}
# UpSet 1: final panels (the protocol's UpSet)
sets1 <- list(`Islet final panel` = pan_i, `Adipose final panel` = pan_a)
save_fig("10_upset_final_panels", 6, 4.5, function() print(upset(fromList(sets1), nsets = 2, order.by = "freq", text.scale = 1.3, mainbar.y.label = "Genes in intersection",
                                                                  sets.x.label = "Panel size")))
# UpSet 2: selection sets (empty sets cannot be drawn by UpSetR and are named in the log instead)
sets2 <- list(`Islet screen+LASSO` = obs_i$lasso, `Adipose screen+LASSO` = obs_a$lasso, `Islet stability >=0.2` = stab_list_i, `Adipose stability >=0.2` = stab_list_a)
empty2 <- names(sets2)[lengths(sets2) == 0]; sets2 <- sets2[lengths(sets2) > 0]
if (length(empty2)) cat("UpSet (selection sets): empty set(s) not drawn: ", paste(empty2, collapse = ", "), "\n", sep = "")
save_fig("10_upset_selection_sets", 8, 5, function() print(upset(fromList(sets2), nsets = length(sets2), order.by = "freq", text.scale = 1.2, mainbar.y.label = "Genes in intersection",
                                                                  sets.x.label = "Set size")))
save_fig("10_permutation_null", 5.5, 4.5, function() {
  tb <- table(factor(perm$overlap, 0:max(perm$overlap, obs_overlap, 1L)))
  bp <- barplot(as.numeric(tb), names.arg = names(tb), xlab = "Overlap size of screen+LASSO selections", ylab = sprintf("Permutations (of %d)", N_PERM),
                main = "Cross-tissue overlap: permutation null", col = "grey70")
  abline(v = bp[as.character(obs_overlap) == names(tb)], col = "red", lwd = 2, lty = 2)
  legend("topright", bty = "n", legend = sprintf("observed = %d, p = %.3f", obs_overlap, p_perm), col = "red", lwd = 2, lty = 2)
})

# ---- final assertions ----------------------------------------------------------------------------------------------------------
stopifnot(file.exists(file.path(out_dir, c("10_gene_classification.csv", "10_overlap_summary.csv", "10_permutation_null.csv", "10_cross_tissue_result.rds"))),
          file.exists(file.path(fig_dir, c("10_upset_final_panels.pdf", "10_upset_final_panels.png", "10_upset_selection_sets.pdf", "10_upset_selection_sets.png",
                                           "10_permutation_null.pdf", "10_permutation_null.png"))),
          nrow(perm) == N_PERM, !anyNA(perm$overlap), p_perm >= 1 / (N_PERM + 1), p_perm <= 1, nrow(cls) == length(all_genes), !anyNA(cls$classification),
          setequal(cls$gene[cls$classification == "shared"], shared), setequal(cls$gene[cls$classification == "islet-specific"], isl_only),
          setequal(cls$gene[cls$classification == "adipose-specific"], ada_only), N > 0, all(perm$islet_prop_same_label < 1), all(perm$adipose_prop_same_label < 1),
          all(perm$overlap <= pmin(perm$islet_lasso_n, perm$adipose_lasso_n)), obs_overlap <= min(length(obs_i$lasso), length(obs_a$lasso)))
cat("\nAll end-of-script assertions passed.\n")
cat(sprintf("\nRuntime: benchmark %.2f s/iteration; permutation %.1f min; total %.1f min\n", tb, t_perm / 60, as.numeric(difftime(Sys.time(), T_START, units = "mins"))))
cat("\nsessionInfo():\n"); print(sessionInfo())
finish()

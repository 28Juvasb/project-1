# 06_nested_cv.R -- nested cross-validation engine (protocol v1.1 §3; CLAUDE.md rules 1,2,5,6,7).
#
# Defines the reusable nested_cv(expr, pheno, geneset, params) (islet and adipose arms) and, when run as a script
# (Rscript scripts/06_nested_cv.R [all|benchmark|smoke]), executes the ISLET arm (GSE41762 alone, rule 4):
#   0. adipose smoke test (code path only; NO adipose performance is computed-for-report or printed)
#   1. runtime benchmark (one outer fold + one full-cohort fit) and runtime estimates
#   2. synthetic positive control (planted genes must be recovered)
#   3. permutation leakage test (20 shuffles x full pipeline x 1 repeat; PASS: mean AUC in [0.45,0.55], 95th pct < 0.65)
#      -> if it FAILS the script stops and NO headline numbers are produced
#   4. real nested CV (10-fold x 5 repeats), stability selection, fallback table
#   5. apparent (full-cohort) AUC, labelled "optimistic (apparent)", beside the nested AUC; final panel + frozen model
# Validation cohorts (GSE76894, GSE20966, GSE71416, GSE16415) are NEVER loaded here (rule 2).
#
# Design decisions stated explicitly (all fixed before running; none tuned on results):
#  * Expression filter (§2.2) for modelling is recomputed INSIDE each fold from the training samples only (25th percentile of
#    the training matrix values; gene "expressed" if above it in >= 20% of training samples). Nothing from test samples
#    (labels or expression) enters any selection step. The full-cohort application is asserted to reproduce 04_<GSE>_expressed.rds.
#  * limma runs on all in-fold-expressed genes (so variance moderation uses the whole expressed transcriptome); the screen then
#    keeps senescence-union genes with nominal p < 0.05 (fallback: top 30 union genes by p if < 10 pass).
#  * LASSO: cv.glmnet type.measure = "deviance" (glmnet default; "auc" is unusable with < 10 obs per inner fold), stratified
#    inner folds, lambda.1se. Selected = non-zero coefficients.
#  * SVM-RFE: caret::rfe with custom e1071 linear-SVM functions (kernlab not installed; caret's own svmLinear needs kernlab).
#    Genuine SVM-RFE ranking by |w| of the linear SVM (C = 1, features scaled), inner stratified CV, metric = ROC AUC of
#    the SVM decision value, sizes 1..min(30,p), optimal subset = caret::pickSizeBest.
#  * Boruta 10.0.0 default importance source (fru random forest), maxRuns = 500, ntree = 1000, 1 thread per worker; Confirmed only.
#  * Ridge classifier: lambda.min from inner stratified CV (the protocol fixes lambda.1se only for LASSO); features z-scored with
#    training means/SDs; standardize = FALSE in glmnet (already z-scored). A single-feature panel is handled with a zero dummy column.
#  * fallback levels: L0 3-way intersection non-empty; L1 genes selected by >= 2 of 3; L2 top 5 LASSO by |coef|;
#    L3 LASSO empty -> top 5 screen genes by p.
#  * CI (rule 6, protocol §3): pROC::ci.auc, 2000 stratified bootstraps on the pooled out-of-fold predictions (all repeats stacked).
#    Because the same sample is predicted once per repeat, stacking treats repeats as independent; this CI is anti-conservative.
#    A second, sample-level stratified bootstrap (resample donors, average per-repeat AUC) is therefore reported beside it.
#
# Outputs (results/stage3/): see the end of the script. Log: logs/06_nested_cv_<timestamp>.log

set.seed(42)
suppressPackageStartupMessages({
  library(here); library(limma); library(glmnet); library(caret); library(e1071); library(Boruta); library(pROC)
  library(parallel)
})

# ================================================================================================================
# helpers
# ================================================================================================================
fast_auc <- function(p, y, pos = "T2D") {          # Mann-Whitney AUC, ties = 0.5; positive class = T2D
  is_pos <- y == pos; n1 <- sum(is_pos); n0 <- sum(!is_pos)
  r <- rank(p); (sum(r[is_pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
roc_obj <- function(y, p) pROC::roc(response = y, predictor = p, levels = c("ND", "T2D"), direction = "<", quiet = TRUE)

stratified_foldid <- function(y, k) {              # stratified assignment of each sample to 1..k (uses current RNG stream)
  id <- integer(length(y))
  for (cl in levels(y)) { ix <- which(y == cl); id[ix] <- sample(rep_len(seq_len(k), length(ix))) }
  id
}

make_outer_folds <- function(y, k, reps, seed) {
  out <- list()
  for (r in seq_len(reps)) {
    set.seed(seed + r); f <- stratified_foldid(y, k)
    for (j in seq_len(k)) out[[length(out) + 1L]] <- list(rep = r, fold = j, test = which(f == j))
  }
  out
}

# analysis phenotype: status (ND/T2D factor), age, sex (inferred where recorded sex mismatches expression, pre-analysis plan B12),
# batch (GSE41762 only). Row order = pheno row order; rownames = GSM.
prep_pheno <- function(pheno, arm = c("islet", "adipose")) {
  arm <- match.arg(arm)
  sex_rec <- tolower(trimws(as.character(pheno$sex))); sex_inf <- tolower(trimws(as.character(pheno$inferred_sex)))
  mism <- pheno$sex_mismatch %in% TRUE
  sex <- ifelse(mism | is.na(sex_rec), sex_inf, sex_rec)
  out <- data.frame(status = factor(pheno$group, levels = c("ND", "T2D")), age = as.numeric(pheno$age),
                    sex = factor(sex), row.names = pheno$gsm, stringsAsFactors = FALSE)
  if (arm == "islet") out$batch <- factor(pheno$submitter_norm_batch)
  stopifnot(!anyNA(out$status), !anyNA(out$age), !anyNA(out$sex), nlevels(out$sex) == 2L)
  if (arm == "islet") stopifnot(!anyNA(out$batch), nlevels(out$batch) == 2L)
  out
}

default_params <- function(arm = c("islet", "adipose")) {
  arm <- match.arg(arm)
  base <- list(screen_p = 0.05, screen_min = 10L, screen_top = 30L, expr_q = 0.25, expr_frac = 0.20,
               rfe_max = 30L, boruta_runs = 500L, boruta_ntree = 1000L, fallback_top_lasso = 5L)
  if (arm == "islet") c(base, list(outer_k = 10L, outer_reps = 5L, inner_k = 10L, design = ~ status + age + sex + batch, arm = "islet"))
  else c(base, list(outer_k = 5L, outer_reps = 10L, inner_k = 3L, design = ~ status + age + sex, arm = "adipose"))
}

# ================================================================================================================
# pipeline steps (each receives TRAINING data only)
# ================================================================================================================
limma_screen <- function(expr_tr, pheno_tr, geneset, params, expected_ids) {
  # SELF-TEST 2 (per fold): matrix handed to limma has exactly the training columns, in training order
  stopifnot(ncol(expr_tr) == length(expected_ids), identical(colnames(expr_tr), expected_ids),
            nrow(pheno_tr) == ncol(expr_tr), identical(rownames(pheno_tr), colnames(expr_tr)))
  q <- quantile(expr_tr, params$expr_q, names = FALSE)
  expressed <- rowMeans(expr_tr > q) >= params$expr_frac          # §2.2 rule, training samples only
  m <- expr_tr[expressed, , drop = FALSE]
  stopifnot(ncol(m) == length(expected_ids), identical(colnames(m), expected_ids))
  design <- model.matrix(params$design, pheno_tr)
  stopifnot(nrow(design) == ncol(m), qr(design)$rank == ncol(design), "statusT2D" %in% colnames(design))
  fit <- eBayes(lmFit(m, design))
  tt <- topTable(fit, coef = "statusT2D", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  su <- tt[tt$gene %in% geneset, c("gene", "logFC", "P.Value")]
  su <- su[order(su$P.Value), ]
  passed <- su$gene[su$P.Value < params$screen_p]
  fallback <- length(passed) < params$screen_min
  screen <- if (fallback) head(su$gene, params$screen_top) else passed
  list(screen = screen, n_pass_nominal = length(passed), fallback = fallback, n_expressed = sum(expressed),
       n_union_expressed = nrow(su), limma_ncol = ncol(m), screen_table = su, expressed_genes = rownames(m))
}

run_lasso <- function(x, y, inner_k) {
  foldid <- stratified_foldid(y, inner_k)
  cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1, foldid = foldid, type.measure = "deviance")
  cf <- coef(cvfit, s = "lambda.1se"); nm <- rownames(cf)[as.vector(cf != 0)]
  list(selected = setdiff(nm, "(Intercept)"), coef = setNames(as.vector(cf)[-1], rownames(cf)[-1]), lambda = cvfit$lambda.1se)
}

svmFuncs <- list(
  summary = function(data, lev = NULL, model = NULL) caret::twoClassSummary(data, lev, model),
  fit = function(x, y, first, last, ...) e1071::svm(x = as.matrix(x), y = y, kernel = "linear", type = "C-classification",
                                                   cost = 1, scale = TRUE),
  pred = function(object, x) {
    dv <- as.vector(attr(predict(object, as.matrix(x), decision.values = TRUE), "decision.values"))
    lv <- object$levels                                      # e1071: positive decision value favours levels[1]
    out <- data.frame(pred = factor(ifelse(dv > 0, lv[1], lv[2]), levels = lv))
    out[[lv[1]]] <- plogis(dv); out[[lv[2]]] <- 1 - plogis(dv); out
  },
  rank = function(object, x, y) {
    w <- drop(t(object$coefs) %*% object$SV)                 # weight vector in scaled feature space
    v <- data.frame(var = colnames(object$SV), Overall = abs(w), stringsAsFactors = FALSE)
    v[order(v$Overall, decreasing = TRUE), , drop = FALSE]
  },
  selectSize = caret::pickSizeBest, selectVar = caret::pickVars
)

run_svmrfe <- function(xs, y, inner_k, rfe_max) {     # xs: samples x genes with safe column names
  p <- ncol(xs); sizes <- seq_len(min(rfe_max, p))
  idx <- caret::createFolds(y, k = inner_k, returnTrain = TRUE)
  ctrl <- caret::rfeControl(functions = svmFuncs, method = "cv", number = inner_k, index = idx,
                            indexOut = lapply(idx, function(i) setdiff(seq_along(y), i)),
                            verbose = FALSE, returnResamp = "none", allowParallel = FALSE)
  fit <- caret::rfe(x = as.data.frame(xs), y = y, sizes = sizes, rfeControl = ctrl, metric = "ROC")
  list(selected = fit$optVariables, opt_size = fit$optsize)
}

run_boruta <- function(xs, y, runs, ntree) {
  b <- Boruta::Boruta(as.data.frame(xs), y, maxRuns = runs, ntree = ntree, num.threads = 1)
  list(selected = names(b$finalDecision)[b$finalDecision == "Confirmed"],
       n_tentative = sum(b$finalDecision == "Tentative"), n_runs = nrow(b$ImpHistory))
}

# Full pipeline on a training set -> frozen model. y/pheno_tr/expr_tr contain TRAINING samples only.
fit_pipeline <- function(expr_tr, pheno_tr, geneset, params) {
  y <- pheno_tr$status; ids <- colnames(expr_tr)
  sc <- limma_screen(expr_tr, pheno_tr, geneset, params, expected_ids = ids)
  genes <- sc$screen; safe <- paste0("G", seq_along(genes))
  x <- t(expr_tr[genes, , drop = FALSE]); xs <- x; colnames(xs) <- safe
  la <- run_lasso(x, y, params$inner_k)
  sv <- run_svmrfe(xs, y, params$inner_k, params$rfe_max); sv$selected <- genes[match(sv$selected, safe)]
  bo <- run_boruta(xs, y, params$boruta_runs, params$boruta_ntree); bo$selected <- genes[match(bo$selected, safe)]
  votes <- table(factor(c(la$selected, sv$selected, bo$selected), levels = genes))
  inter <- Reduce(intersect, list(la$selected, sv$selected, bo$selected))
  if (length(inter) > 0) { feat <- inter; level <- "L0_intersection" }
  else {
    two <- names(votes)[votes >= 2]
    if (length(two) > 0) { feat <- two; level <- "L1_>=2of3" }
    else if (length(la$selected) > 0) {
      o <- order(abs(la$coef[la$selected]), decreasing = TRUE); feat <- la$selected[o][seq_len(min(params$fallback_top_lasso, length(o)))]
      level <- "L2_top5_LASSO"
    } else { feat <- head(sc$screen_table$gene[sc$screen_table$gene %in% genes], params$fallback_top_lasso); level <- "L3_top5_screen" }
  }
  stopifnot(length(feat) >= 1, all(feat %in% genes))
  # ridge on features, training-only z-scoring
  xf <- x[, feat, drop = FALSE]; mu <- colMeans(xf); sdv <- apply(xf, 2, sd); sdv[sdv == 0] <- 1
  z <- sweep(sweep(xf, 2, mu), 2, sdv, "/"); zz <- if (ncol(z) == 1L) cbind(z, dummy0 = 0) else z
  rg <- cv.glmnet(zz, y, family = "binomial", alpha = 0, standardize = FALSE, foldid = stratified_foldid(y, params$inner_k),
                  type.measure = "deviance")
  cf <- coef(rg, s = "lambda.min"); b0 <- as.numeric(cf[1]); b <- setNames(as.vector(cf)[-1][seq_along(feat)], feat)
  list(features = feat, level = level, lasso = la$selected, svmrfe = sv$selected, svm_size = sv$opt_size, boruta = bo$selected,
       boruta_tentative = bo$n_tentative, boruta_runs = bo$n_runs, screen = genes, screen_info = sc[c("n_pass_nominal", "fallback",
       "n_expressed", "n_union_expressed", "limma_ncol")], model = list(genes = feat, mu = mu, sd = sdv, coef = b, intercept = b0,
       lambda = rg$lambda.min), ridge_check = list(cvfit = rg, zz = zz))
}

predict_model <- function(model, expr_new) {          # frozen model applied with TRAINING mu/sd (z-scoring with training params)
  x <- t(expr_new[model$genes, , drop = FALSE]); z <- sweep(sweep(x, 2, model$mu), 2, model$sd, "/")
  plogis(model$intercept + as.vector(z %*% model$coef[model$genes]))
}

# ================================================================================================================
# one outer fold, and the nested CV driver
# ================================================================================================================
run_fold <- function(task, expr, pheno, geneset, params, base_seed, keep_fit = FALSE) {
  t0 <- proc.time()[3]
  set.seed(base_seed + task$rep * 1000L + task$fold)               # per-fold derived seed -> deterministic under any scheduling
  test <- task$test; train <- setdiff(seq_len(ncol(expr)), test)
  # SELF-TEST 1 (per fold)
  stopifnot(length(intersect(train, test)) == 0L, length(train) + length(test) == ncol(expr))
  tr_ids <- colnames(expr)[train]; te_ids <- colnames(expr)[test]
  warns <- character(0)
  fit <- withCallingHandlers(
    fit_pipeline(expr[, train, drop = FALSE], pheno[train, , drop = FALSE], geneset, params),   # test columns are never passed in
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  stopifnot(fit$screen_info$limma_ncol == length(train))                                       # SELF-TEST 2
  prob <- predict_model(fit$model, expr[, test, drop = FALSE])
  res <- list(rep = task$rep, fold = task$fold, train_ids = tr_ids, test_ids = te_ids, n_train = length(train),
              limma_ncol = fit$screen_info$limma_ncol, limma_saw_test_id = any(te_ids %in% tr_ids),
              n_expressed = fit$screen_info$n_expressed, n_union_expressed = fit$screen_info$n_union_expressed,
              n_pass_nominal = fit$screen_info$n_pass_nominal, screen_fallback = fit$screen_info$fallback,
              screen = fit$screen, lasso = fit$lasso, svmrfe = fit$svmrfe, svm_size = fit$svm_size, boruta = fit$boruta,
              boruta_tentative = fit$boruta_tentative, boruta_runs = fit$boruta_runs, features = fit$features, level = fit$level,
              lambda = fit$model$lambda,
              pred = data.frame(rep = task$rep, fold = task$fold, id = te_ids, truth = as.character(pheno$status[test]), prob = prob,
                                stringsAsFactors = FALSE),
              warnings = unique(warns), elapsed = proc.time()[3] - t0)
  if (keep_fit) res$fit <- fit
  res
}

par_map <- function(X, FUN, cores) {
  if (cores <= 1L) return(lapply(X, FUN))
  r <- mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
  bad <- vapply(r, function(z) inherits(z, "try-error") || is.null(z), logical(1))
  if (any(bad)) stop("worker failure: ", paste(unique(vapply(r[bad], function(z) as.character(z)[1], "")), collapse = " | "))
  r
}

# nested_cv(): expr = genes x samples (colnames = sample IDs), pheno from prep_pheno() (rownames = sample IDs, status factor),
# geneset = character vector, params = default_params(arm). Returns per-fold results, OOF predictions and per-repeat AUCs.
nested_cv <- function(expr, pheno, geneset, params, seed = 42L, cores = 1L, keep_folds = TRUE) {
  stopifnot(is.matrix(expr), identical(colnames(expr), rownames(pheno)), is.factor(pheno$status), identical(levels(pheno$status), c("ND", "T2D")))
  tasks <- make_outer_folds(pheno$status, params$outer_k, params$outer_reps, seed)
  # SELF-TEST 1 (global): within every repeat, test sets partition the samples exactly once; no train/test overlap in any fold
  for (r in seq_len(params$outer_reps)) {
    tt <- unlist(lapply(Filter(function(z) z$rep == r, tasks), `[[`, "test"))
    stopifnot(!anyDuplicated(tt), setequal(tt, seq_len(ncol(expr))))
  }
  res <- par_map(tasks, function(tk) run_fold(tk, expr, pheno, geneset, params, base_seed = seed), cores)
  oof <- do.call(rbind, lapply(res, `[[`, "pred"))
  # global assertions for self-tests 1 and 2
  for (f in res) {
    stopifnot(length(intersect(f$train_ids, f$test_ids)) == 0L, f$limma_ncol == f$n_train, !f$limma_saw_test_id,
              f$n_train == ncol(expr) - length(f$test_ids))
  }
  stopifnot(!anyNA(oof$prob), all(oof$prob >= 0 & oof$prob <= 1), all(table(oof$id) == params$outer_reps))
  auc_rep <- vapply(seq_len(params$outer_reps), function(r) { o <- oof[oof$rep == r, ]; fast_auc(o$prob, factor(o$truth, c("ND", "T2D"))) }, 0)
  out <- list(oof = oof, auc_per_rep = auc_rep, auc_mean = mean(auc_rep), params = params, n_folds = length(res))
  if (keep_folds) out$folds <- res else out$folds <- lapply(res, function(f) f[c("rep", "fold", "level", "screen_fallback", "features", "elapsed")])
  out
}

# AUC summaries with CIs (rule 6)
summarise_auc <- function(oof, outer_reps, ids_status, B = 2000L, seed = 42L) {
  y <- factor(oof$truth, c("ND", "T2D"))
  per_rep <- vapply(seq_len(outer_reps), function(r) { o <- oof$rep == r; fast_auc(oof$prob[o], y[o]) }, 0)
  set.seed(seed); r_all <- roc_obj(y, oof$prob)
  ci_p <- as.numeric(pROC::ci.auc(r_all, method = "bootstrap", boot.n = B, boot.stratified = TRUE, progress = "none"))
  # sample-level stratified bootstrap (resample donors; AUC = mean over repeats on the resampled donor set)
  ids <- names(ids_status); P <- sapply(seq_len(outer_reps), function(r) { o <- oof[oof$rep == r, ]; o$prob[match(ids, o$id)] })
  stopifnot(nrow(P) == length(ids), !anyNA(P))
  i1 <- which(ids_status == "T2D"); i0 <- which(ids_status == "ND"); set.seed(seed + 1L)
  bs <- replicate(B, { i <- c(sample(i0, replace = TRUE), sample(i1, replace = TRUE)); yb <- ids_status[i]
    mean(apply(P[i, , drop = FALSE], 2, fast_auc, y = yb)) })
  list(auc_mean_over_repeats = mean(per_rep), auc_per_repeat = per_rep, auc_pooled = as.numeric(pROC::auc(r_all)),
       ci_protocol = c(lower = ci_p[1], upper = ci_p[3]), ci_sample_level = unname(quantile(bs, c(.025, .975))), roc = r_all)
}

# stability selection from per-fold lists
stability_table <- function(folds) {
  n <- length(folds)
  cnt <- function(field) table(unlist(lapply(folds, function(f) unique(f[[field]]))))
  f_final <- cnt("features"); genes <- names(f_final)
  all_g <- unique(unlist(lapply(folds, function(f) c(f$features, f$lasso, f$svmrfe, f$boruta, f$screen))))
  g <- function(tab) as.numeric(tab[match(all_g, names(tab))]) %|% 0
  `%|%` <- function(a, b) { a[is.na(a)] <- b; a }
  d <- data.frame(gene = all_g, n_folds_final = g(f_final), stringsAsFactors = FALSE)
  d$frequency <- d$n_folds_final / n
  d$freq_screen <- g(cnt("screen")) / n; d$freq_lasso <- g(cnt("lasso")) / n
  d$freq_svmrfe <- g(cnt("svmrfe")) / n; d$freq_boruta <- g(cnt("boruta")) / n
  d <- d[order(-d$n_folds_final, -d$freq_lasso, d$gene), ]; rownames(d) <- NULL
  d[d$n_folds_final > 0 | d$freq_lasso > 0 | d$freq_svmrfe > 0 | d$freq_boruta > 0, ]
}

# ================================================================================================================
# synthetic positive control
# ================================================================================================================
make_synthetic <- function(geneset, pheno, n_genes = 3000L, n_plant = 10L, shift_sd = 1.5, seed = 7L) {
  set.seed(seed); n <- nrow(pheno); gs <- sort(geneset)
  noise_genes <- paste0("SYNTH", seq_len(n_genes - length(gs)))
  nm <- c(gs, noise_genes); mu <- runif(length(nm), 5, 11); sdg <- 0.5
  X <- matrix(rnorm(length(nm) * n, mu, sdg), nrow = length(nm), dimnames = list(nm, rownames(pheno)))
  planted <- sample(gs, n_plant)
  t2d <- pheno$status == "T2D"
  X[planted, t2d] <- X[planted, t2d] + shift_sd * sdg
  ph <- pheno; ph$age <- round(rnorm(n, 55, 10)); ph$sex <- factor(sample(c("female", "male"), n, TRUE))   # covariates independent of status
  if ("batch" %in% names(ph)) ph$batch <- factor(sample(c("A", "B"), n, TRUE))
  list(expr = X, pheno = ph, planted = planted)
}

# ================================================================================================================
# main
# ================================================================================================================
if (sys.nframe() == 0L && !isTRUE(getOption("ncv.no_main"))) {
  args <- commandArgs(trailingOnly = TRUE); mode <- if (length(args)) args[1] else "all"
  stopifnot(mode %in% c("all", "benchmark", "smoke"))
  script_name <- "06_nested_cv"; ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); T_START <- Sys.time()
  dir.create(here("logs"), showWarnings = FALSE); out_dir <- here("results", "stage3"); fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt"); sink(con, split = TRUE); sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)
  CORES <- 8L; SEED <- 42L
  cat("Script:", script_name, " mode:", mode, " run:", ts, " cores:", CORES, "\n\nInputs:\n")
  f_in <- c("data/processed/04_GSE41762_expr.rds", "data/processed/04_GSE41762_pheno.rds", "data/processed/04_GSE41762_expressed.rds",
            "genesets/senescence_union.rds", "data/processed/04_GSE27951_expr.rds", "data/processed/04_GSE27951_pheno.rds")
  for (f in f_in) cat("  ", f, " md5=", tools::md5sum(here(f)), "\n", sep = "")
  cat("  (validation cohorts GSE76894, GSE20966, GSE71416, GSE16415 are NOT loaded)\n")

  geneset <- readRDS(here("genesets/senescence_union.rds")); stopifnot(is.character(geneset), length(geneset) == 1151L, !anyDuplicated(geneset))
  load_arm <- function(gse, arm) {
    e <- readRDS(here("data/processed", sprintf("04_%s_expr.rds", gse))); p <- readRDS(here("data/processed", sprintf("04_%s_pheno.rds", gse)))
    stopifnot(identical(colnames(e), p$gsm)); list(expr = e, pheno = prep_pheno(p, arm), raw = p)
  }
  isl <- load_arm("GSE41762", "islet"); P <- default_params("islet")
  stopifnot(ncol(isl$expr) == 77L, table(isl$pheno$status)[["T2D"]] == 20L, table(isl$pheno$status)[["ND"]] == 57L)
  cat("\nPre-run checks (islet GSE41762): n =", ncol(isl$expr), "genes =", nrow(isl$expr), " ND/T2D =", paste(table(isl$pheno$status), collapse = "/"), "\n")
  cat("sex covariate: recorded, with expression-inferred value for sex_mismatch rows (", paste(isl$raw$gsm[isl$raw$sex_mismatch %in% TRUE], collapse = ","), ")\n", sep = "")
  # expressed-gene filter reproduces the dataset-level flag when applied to the full cohort (in-fold rule == §2.2 rule)
  fl <- readRDS(here("data/processed/04_GSE41762_expressed.rds"))
  q <- quantile(isl$expr, P$expr_q, names = FALSE); mine <- rownames(isl$expr)[rowMeans(isl$expr > q) >= P$expr_frac]
  stopifnot(setequal(mine, fl$gene[fl$expressed]))
  cat("In-fold expression filter applied to the full cohort reproduces 04_GSE41762_expressed.rds exactly (", length(mine), " genes).\n", sep = "")
  cat("senescence union genes: ", length(geneset), "; present in matrix: ", sum(geneset %in% rownames(isl$expr)), "; expressed: ", sum(geneset %in% mine), "\n", sep = "")
  cat("\nParameters (islet): outer ", P$outer_k, "-fold x ", P$outer_reps, " reps; inner ", P$inner_k, "-fold; design ", deparse(P$design),
      "; screen p<", P$screen_p, " (fallback top ", P$screen_top, " if <", P$screen_min, "); SVM-RFE sizes 1..", P$rfe_max,
      "; Boruta maxRuns=", P$boruta_runs, " ntree=", P$boruta_ntree, "; ridge lambda.min; LASSO lambda.1se; seed ", SEED, "\n", sep = "")

  # ---- 0. adipose smoke test (code path only; no performance is reported) --------------------------------------
  cat("\n=== 0. adipose smoke test (code path only; NO adipose performance is reported or saved) ===\n")
  ad <- load_arm("GSE27951", "adipose"); PA <- default_params("adipose"); PA$outer_reps <- 1L
  stopifnot(ncol(ad$expr) == 23L)
  smk <- nested_cv(ad$expr, ad$pheno, geneset, PA, seed = SEED, cores = CORES)
  stopifnot(smk$n_folds == PA$outer_k * PA$outer_reps)
  cat("adipose smoke: ", smk$n_folds, " folds ran (outer ", PA$outer_k, ", inner ", PA$inner_k, "); self-tests 1-2 passed; fallback levels: ",
      paste(names(table(vapply(smk$folds, `[[`, "", "level"))), table(vapply(smk$folds, `[[`, "", "level")), sep = "=", collapse = ", "),
      "; screen fallback in ", sum(vapply(smk$folds, `[[`, NA, "screen_fallback")), " folds. (AUC deliberately not printed.)\n", sep = "")
  rm(smk); gc(verbose = FALSE)
  if (mode == "smoke") { cat("\nsmoke mode: done.\n"); quit(save = "no") }

  # ---- 1. benchmark ------------------------------------------------------------------------------------------------
  cat("\n=== 1. runtime benchmark (islet) ===\n")
  tasks1 <- make_outer_folds(isl$pheno$status, P$outer_k, 1L, SEED)
  tb <- system.time(f1 <- run_fold(tasks1[[1]], isl$expr, isl$pheno, geneset, P, base_seed = SEED))[3]
  tf <- system.time(fullfit <- fit_pipeline(isl$expr, isl$pheno, geneset, P))[3]
  cat(sprintf("single outer fold: %.1f s (screen %d genes, features %d, level %s); full-cohort pipeline: %.1f s\n", tb, length(f1$screen),
              length(f1$features), f1$level, tf))
  eff <- 0.85
  est <- c(real_50_folds = 50 * tb / (CORES * eff), permutation_200_folds = 200 * tb / (CORES * eff), synthetic_20_folds = 20 * tb / (CORES * eff), final_fit = tf)
  cat(sprintf("estimated wall time (8 cores, 85%% parallel efficiency): real 50 folds %.1f min; permutation 20x10 folds %.1f min; synthetic 20 folds %.1f min; final fit %.1f min; TOTAL %.1f min\n",
              est[1] / 60, est[2] / 60, est[3] / 60, est[4] / 60, sum(est) / 60))
  # reproducibility self-test (same seed -> same result, in-process and in a forked worker), and test-label irrelevance
  f1b <- run_fold(tasks1[[1]], isl$expr, isl$pheno, geneset, P, base_seed = SEED)
  f1c <- par_map(list(tasks1[[1]]), function(tk) run_fold(tk, isl$expr, isl$pheno, geneset, P, base_seed = SEED), 2L)[[1]]
  ph_s <- isl$pheno; te <- tasks1[[1]]$test; set.seed(99); ph_s$status[te] <- sample(ph_s$status[te])
  f1d <- run_fold(tasks1[[1]], isl$expr, ph_s, geneset, P, base_seed = SEED)
  stopifnot(identical(f1$features, f1b$features), isTRUE(all.equal(f1$pred$prob, f1b$pred$prob)),
            identical(f1$features, f1c$features), isTRUE(all.equal(f1$pred$prob, f1c$pred$prob)),
            identical(f1$features, f1d$features), isTRUE(all.equal(f1$pred$prob, f1d$pred$prob)))
  cat("SELF-TEST (extra) PASS: a fold re-run (serial, forked) is bit-identical; scrambling TEST-fold labels leaves selected features and predictions unchanged.\n")
  # ridge predict check and svm direction check
  rc <- fullfit$ridge_check; zz <- rc$zz
  stopifnot(isTRUE(all.equal(as.vector(predict(rc$cvfit, zz, s = "lambda.min", type = "response")),
                             as.vector(plogis(fullfit$model$intercept + as.vector(zz[, seq_along(fullfit$model$genes), drop = FALSE] %*% fullfit$model$coef))), tolerance = 1e-8)))
  sv_chk <- svmFuncs$fit(matrix(c(rnorm(20, 0), rnorm(20, 3)), ncol = 1, dimnames = list(NULL, "a")), factor(rep(c("ND", "T2D"), each = 20), c("ND", "T2D")))
  pp <- svmFuncs$pred(sv_chk, matrix(c(0, 3), ncol = 1, dimnames = list(NULL, "a"))); stopifnot(pp$ND[1] > 0.5, pp$T2D[2] > 0.5)
  cat("SELF-TEST (extra) PASS: manual frozen-model prediction equals glmnet predict; SVM decision-value direction (positive = ND) verified.\n")
  if (mode == "benchmark") { cat("\nbenchmark mode: done.\n"); quit(save = "no") }

  # ---- 2. synthetic positive control ------------------------------------------------------------------------------
  cat("\n=== 2. synthetic positive control ===\n")
  cat("Design (fixed before running): 3000 genes (the 1151 senescence-union names + 1849 SYNTH noise genes), 77 samples with the real ND/T2D labels,\n",
      "10 randomly chosen union genes shifted +1.5 SD in T2D, all other genes pure noise (SD 0.5), covariates random and independent of status.\n",
      "Pipeline = production nested_cv, 10-fold x 2 repeats. Pre-stated PASS: nested AUC >= 0.90 AND >= 8 of the top-10 stability genes are planted.\n", sep = "")
  syn <- make_synthetic(geneset, isl$pheno)
  PS <- P; PS$outer_reps <- 2L
  t_syn <- system.time(sy <- nested_cv(syn$expr, syn$pheno, geneset, PS, seed = SEED, cores = CORES))[3]
  sy_auc <- summarise_auc(sy$oof, PS$outer_reps, setNames(as.character(syn$pheno$status), rownames(syn$pheno)))
  sy_stab <- stability_table(sy$folds)
  top10 <- head(sy_stab$gene, 10); n_rec <- sum(top10 %in% syn$planted)
  cat(sprintf("synthetic nested-CV AUC = %.3f (protocol CI %.3f-%.3f; sample-level CI %.3f-%.3f)\n", sy_auc$auc_mean_over_repeats,
              sy_auc$ci_protocol[1], sy_auc$ci_protocol[2], sy_auc$ci_sample_level[1], sy_auc$ci_sample_level[2]))
  cat("planted genes:", paste(sort(syn$planted), collapse = ", "), "\n")
  cat("planted-gene stability frequencies:", paste(sprintf("%s=%.2f", syn$planted, sy_stab$frequency[match(syn$planted, sy_stab$gene)]), collapse = ", "), "\n")
  cat("top-10 stability genes:", paste(sprintf("%s(%.2f)%s", sy_stab$gene[1:10], sy_stab$frequency[1:10], ifelse(top10 %in% syn$planted, "*", "")), collapse = ", "), " (* = planted)\n")
  cat(sprintf("non-planted gene with highest final-panel frequency: %.2f\n", max(sy_stab$frequency[!sy_stab$gene %in% syn$planted])))
  cat("synthetic fallback levels:", paste(names(table(vapply(sy$folds, `[[`, "", "level"))), table(vapply(sy$folds, `[[`, "", "level")), sep = "=", collapse = ", "), sprintf("(%.0f s)\n", t_syn))
  syn_pass <- sy_auc$auc_mean_over_repeats >= 0.90 && n_rec >= 8
  cat("POSITIVE CONTROL:", if (syn_pass) "PASS" else "FAIL", " (", n_rec, "/10 planted in top 10)\n", sep = "")
  saveRDS(list(auc = sy_auc[names(sy_auc) != "roc"], planted = syn$planted, stability = sy_stab, pass = syn_pass), file.path(out_dir, "06_synthetic_control.rds"))
  stopifnot(syn_pass)

  # ---- 3. permutation leakage test --------------------------------------------------------------------------------
  cat("\n=== 3. permutation leakage test (20 label shuffles, full pipeline, 1 repeat of 10-fold outer CV) ===\n")
  PP <- P; PP$outer_reps <- 1L
  t_perm <- system.time(perm <- par_map(seq_len(20), function(i) {
    set.seed(4200L + i); ph <- isl$pheno; ph$status <- sample(ph$status)          # independent shuffle of status; covariates untouched
    stopifnot(table(ph$status) == table(isl$pheno$status))
    r <- nested_cv(isl$expr, ph, geneset, PP, seed = SEED + 100L * i, cores = 1L, keep_folds = FALSE)
    c(auc = r$auc_mean, same_as_real = mean(ph$status == isl$pheno$status),
      fb = mean(vapply(r$folds, `[[`, "", "level") != "L0_intersection"), nf = mean(vapply(r$folds, function(f) length(f$features), 0)))
  }, CORES))[3]
  pm <- do.call(rbind, perm); pa <- pm[, "auc"]
  perm_res <- list(auc = pa, mean = mean(pa), sd = sd(pa), q95 = unname(quantile(pa, .95)), min = min(pa), max = max(pa))
  perm_pass <- perm_res$mean >= 0.45 && perm_res$mean <= 0.55 && perm_res$q95 < 0.65
  write.csv(data.frame(shuffle = 1:20, auc = round(pa, 4), prop_labels_equal_to_real = round(pm[, "same_as_real"], 3),
                       frac_folds_fallback = round(pm[, "fb"], 2), mean_panel_size = round(pm[, "nf"], 2)),
            file.path(out_dir, "06_permutation_test_islet.csv"), row.names = FALSE)
  cat("20 permutation AUCs:", paste(sprintf("%.3f", sort(pa)), collapse = " "), "\n")
  cat(sprintf("mean %.3f  SD %.3f  95th pct %.3f  min %.3f  max %.3f  (%.1f min)\n", perm_res$mean, perm_res$sd, perm_res$q95, perm_res$min, perm_res$max, t_perm / 60))
  cat("Criteria (fixed): mean in [0.45,0.55] AND 95th percentile < 0.65 -> ", if (perm_pass) "PASS" else "FAIL", "\n", sep = "")
  saveRDS(perm_res, file.path(out_dir, "06_permutation_test_islet.rds"))
  if (!perm_pass) {
    cat("\n*** PERMUTATION TEST FAILED: possible leakage. STOPPING. No headline numbers are produced. Thresholds and pipeline NOT adjusted. ***\n")
    quit(save = "no", status = 2)
  }

  # ---- 4. real nested CV --------------------------------------------------------------------------------------------
  cat("\n=== 4. real nested CV (islet GSE41762, 10-fold x 5 repeats) ===\n")
  t_real <- system.time(real <- nested_cv(isl$expr, isl$pheno, geneset, P, seed = SEED, cores = CORES))[3]
  cat(sprintf("done in %.1f min; folds = %d\n", t_real / 60, real$n_folds))
  stab <- stability_table(real$folds)
  lv <- vapply(real$folds, `[[`, "", "level")
  fb <- data.frame(level = c("L0_intersection", "L1_>=2of3", "L2_top5_LASSO", "L3_top5_screen"))
  fb$n_folds <- as.integer(table(factor(lv, fb$level))); fb$frequency <- fb$n_folds / length(lv)
  sfb <- data.frame(item = c("screen_fallback_top30 (<10 pass p<0.05)", "LASSO_empty"),
                    n_folds = c(sum(vapply(real$folds, `[[`, NA, "screen_fallback")), sum(vapply(real$folds, function(f) length(f$lasso) == 0, NA))))
  sfb$frequency <- sfb$n_folds / length(lv)
  psz <- vapply(real$folds, function(f) length(f$features), 0)
  cat("fallback levels:\n"); print(fb, row.names = FALSE); print(sfb, row.names = FALSE)
  cat(sprintf("panel size per fold: median %g (range %g-%g); screen genes per fold median %g; Boruta tentative per fold median %g; warnings in %d folds\n",
              median(psz), min(psz), max(psz), median(vapply(real$folds, function(f) length(f$screen), 0)), median(vapply(real$folds, `[[`, 0, "boruta_tentative")),
              sum(vapply(real$folds, function(f) length(f$warnings) > 0, NA))))
  wmsg <- unique(unlist(lapply(real$folds, `[[`, "warnings"))); if (length(wmsg)) cat("distinct warning messages:\n", paste(" -", wmsg, collapse = "\n"), "\n")
  idst <- setNames(as.character(isl$pheno$status), rownames(isl$pheno))
  sm <- summarise_auc(real$oof, P$outer_reps, idst)

  # ---- 5. apparent AUC, final panel ---------------------------------------------------------------------------------
  cat("\n=== 5. full-cohort pipeline (final panel, protocol §3.2) and apparent AUC ===\n")
  set.seed(SEED); final <- fit_pipeline(isl$expr, isl$pheno, geneset, P)
  p_app <- predict_model(final$model, isl$expr); r_app <- roc_obj(isl$pheno$status, p_app)
  set.seed(SEED); ci_app <- as.numeric(pROC::ci.auc(r_app, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE, progress = "none"))
  stopifnot(abs(as.numeric(pROC::auc(r_app)) - fast_auc(p_app, isl$pheno$status)) < 1e-9)
  cat("final panel (", length(final$features), " genes, ", final$level, "): ", paste(final$features, collapse = ", "), "\n", sep = "")
  cat(sprintf("selected by LASSO %d, SVM-RFE %d, Boruta %d (screen %d genes)\n", length(final$lasso), length(final$svmrfe), length(final$boruta), length(final$screen)))

  cat("\n================ RESULTS (islet, GSE41762 discovery) ================\n")
  cat(sprintf("NESTED-CV AUC (honest, headline): %.3f  [per-repeat AUCs: %s]\n", sm$auc_mean_over_repeats, paste(sprintf("%.3f", sm$auc_per_repeat), collapse = ", ")))
  cat(sprintf("   95%% CI, protocol (pROC::ci.auc, 2000 stratified bootstraps on pooled OOF predictions): %.3f - %.3f   [pooled-stack AUC %.3f]\n",
              sm$ci_protocol[1], sm$ci_protocol[2], sm$auc_pooled))
  cat(sprintf("   95%% CI, sample-level stratified bootstrap (conservative; repeats are not independent): %.3f - %.3f\n", sm$ci_sample_level[1], sm$ci_sample_level[2]))
  cat(sprintf("APPARENT AUC, optimistic (apparent): %.3f   95%% CI (pROC, 2000 bootstraps) %.3f - %.3f   [full-cohort model scored on its own training data]\n",
              as.numeric(pROC::auc(r_app)), ci_app[1], ci_app[3]))
  cat("Top 20 genes by stability (final-panel selection frequency over ", real$n_folds, " outer folds):\n", sep = ""); print(head(stab, 20), row.names = FALSE)

  # ---- outputs ----------------------------------------------------------------------------------------------------------
  write.csv(stab[, c("gene", "frequency", "n_folds_final", "freq_screen", "freq_lasso", "freq_svmrfe", "freq_boruta")], file.path(out_dir, "stability_selection.csv"), row.names = FALSE)
  write.csv(rbind(cbind(kind = "fallback_level", fb), cbind(kind = "other", setNames(sfb, names(fb)))), file.path(out_dir, "06_fallback_frequency_islet.csv"), row.names = FALSE)
  saveRDS(lapply(real$folds, function(f) f[c("rep", "fold", "train_ids", "test_ids", "n_train", "limma_ncol", "n_expressed", "n_union_expressed", "n_pass_nominal",
                                            "screen_fallback", "screen", "lasso", "svmrfe", "svm_size", "boruta", "boruta_tentative", "boruta_runs", "features", "level", "lambda", "warnings", "elapsed")]),
          file.path(out_dir, "06_per_fold_features_islet.rds"))
  write.csv(real$oof, file.path(out_dir, "06_oof_predictions_islet.csv"), row.names = FALSE)
  saveRDS(list(genes = final$model$genes, coef = final$model$coef, intercept = final$model$intercept, mu = final$model$mu, sd = final$model$sd, lambda = final$model$lambda,
               fallback_level = final$level, lasso = final$lasso, svmrfe = final$svmrfe, boruta = final$boruta, screen = final$screen, params = P[names(P) != "design"],
               design = deparse(P$design), source = "GSE41762 full cohort (04_ post-QC matrix); frozen model for protocol §7"),
          file.path(out_dir, "06_final_model_islet.rds"))
  write.csv(data.frame(gene = final$model$genes, ridge_coef_z = final$model$coef, train_mean = final$model$mu, train_sd = final$model$sd, fallback_level = final$level),
            file.path(out_dir, "06_final_panel_islet.csv"), row.names = FALSE)
  auc_tab <- data.frame(metric = c("nested-CV AUC (honest, headline) mean of per-repeat AUCs", "nested-CV AUC pooled stack", "apparent AUC, optimistic (apparent)"),
                        auc = c(sm$auc_mean_over_repeats, sm$auc_pooled, as.numeric(pROC::auc(r_app))),
                        ci_lo_protocol = c(sm$ci_protocol[1], sm$ci_protocol[1], ci_app[1]), ci_hi_protocol = c(sm$ci_protocol[2], sm$ci_protocol[2], ci_app[3]),
                        ci_lo_sample_level = c(sm$ci_sample_level[1], NA, NA), ci_hi_sample_level = c(sm$ci_sample_level[2], NA, NA))
  write.csv(auc_tab, file.path(out_dir, "06_auc_summary_islet.csv"), row.names = FALSE)
  saveRDS(list(nested = sm[names(sm) != "roc"], apparent = list(auc = as.numeric(pROC::auc(r_app)), ci = ci_app[c(1, 3)]), permutation = perm_res,
               synthetic_pass = syn_pass, params = P[names(P) != "design"], runtime_min = c(permutation = t_perm / 60, real = t_real / 60), estimates_min = est / 60),
          file.path(out_dir, "06_nested_cv_islet_result.rds"))
  # figure
  draw <- function() {
    op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1)); on.exit(par(op))
    plot(sm$roc, col = "navy", lwd = 2, legacy.axes = TRUE, main = "Islet discovery: out-of-fold vs apparent", print.auc = FALSE)
    plot(r_app, col = "firebrick", lwd = 2, lty = 2, add = TRUE, legacy.axes = TRUE)
    legend("bottomright", bty = "n", cex = 0.8, col = c("navy", "firebrick"), lwd = 2, lty = c(1, 2),
           legend = c(sprintf("nested CV (honest) %.2f [%.2f-%.2f]", sm$auc_mean_over_repeats, sm$ci_protocol[1], sm$ci_protocol[2]),
                      sprintf("optimistic (apparent) %.2f", as.numeric(pROC::auc(r_app)))))
    stripchart(pa, vertical = TRUE, method = "jitter", pch = 19, col = "grey30", ylim = c(0.2, 0.8), ylab = "AUC", main = "Permutation test (20 shuffles)", xlim = c(0.5, 1.5))
    abline(h = c(0.45, 0.55), lty = 3); abline(h = 0.65, lty = 2, col = "red"); points(1, mean(pa), pch = 4, cex = 2, lwd = 2, col = "blue")
  }
  pdf(file.path(fig_dir, "06_islet_roc_permutation.pdf"), 9, 4.5); draw(); dev.off()
  png(file.path(fig_dir, "06_islet_roc_permutation.png"), 1800, 900, res = 200); draw(); dev.off()

  # ---- final assertions ---------------------------------------------------------------------------------------------
  stopifnot(real$n_folds == 50L, nrow(real$oof) == 77L * 5L, !anyNA(real$oof$prob), sum(fb$n_folds) == 50L, nrow(stab) > 0,
            all(stab$frequency >= 0 & stab$frequency <= 1), identical(order(-stab$n_folds_final, -stab$freq_lasso, stab$gene), seq_len(nrow(stab))),
            length(final$features) >= 1, all(final$features %in% geneset), sm$auc_mean_over_repeats >= 0 && sm$auc_mean_over_repeats <= 1,
            all(vapply(real$folds, function(f) f$limma_ncol == f$n_train && f$n_train == 77 - length(f$test_ids), NA)),
            all(vapply(real$folds, function(f) !any(f$train_ids %in% f$test_ids), NA)), length(perm) == 20L)
  cat("\nAll end-of-script assertions passed (50 folds; 5 predictions per donor; no NA; fallback table sums to 50; self-tests 1-3 passed).\n")
  cat(sprintf("\nRuntime: benchmark fold %.1f s; synthetic %.1f min; permutation %.1f min; real %.1f min; total %.1f min\n", tb, t_syn / 60, t_perm / 60, t_real / 60,
              as.numeric(difftime(Sys.time(), T_START, units = "mins"))))
  cat("\nsessionInfo():\n"); print(sessionInfo())
}

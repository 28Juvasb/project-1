# Diagnostic (NOT the pipeline): swap ONLY the Boruta importance backend. LASSO and SVM-RFE selections are taken from the stored production
# per-fold results (they do not depend on Boruta); Boruta is re-run with ranger (getImpRfZ) on each of the 50 training splits (same screen genes),
# and the full-cohort final fit is re-run with fru / ranger / randomForest. Validation cohorts are not loaded.
options(ncv.no_main = TRUE)
.libPaths(c(path.expand("~/R/lib"), .libPaths()))
library(here); library(parallel); source(here("scripts/06_nested_cv.R"))
out_dir <- here("results/stage3/diagnostics")
geneset <- readRDS(here("genesets/senescence_union.rds"))
e <- readRDS(here("data/processed/04_GSE41762_expr.rds")); p <- readRDS(here("data/processed/04_GSE41762_pheno.rds"))
ph <- prep_pheno(p, "islet"); P <- default_params("islet")
stored <- readRDS(here("results/stage3/06_per_fold_features_islet.rds"))
mk_boruta <- function(imp) function(xs, y, runs, ntree) {
  b <- Boruta::Boruta(as.data.frame(xs), y, maxRuns = runs, ntree = ntree, num.threads = 1, getImp = imp)
  list(selected = names(b$finalDecision)[b$finalDecision == "Confirmed"], n_tentative = sum(b$finalDecision == "Tentative"), n_runs = nrow(b$ImpHistory))
}
final_fit <- function(nm) {
  imp <- switch(nm, fru = Boruta::getImpFruZ, ranger = Boruta::getImpRfZ, rf = Boruta::getImpLegacyRfZ)
  run_boruta <<- mk_boruta(imp)             # local to this forked worker
  set.seed(42); f <- suppressWarnings(fit_pipeline(e, ph, geneset, P))
  list(kind = "final", backend = nm, features = f$features, level = f$level, lasso = f$lasso, svmrfe = f$svmrfe, boruta = f$boruta, coef = f$model$coef)
}
fold_fit <- function(i) {
  st <- stored[[i]]; run_boruta <<- mk_boruta(Boruta::getImpRfZ)
  tr <- match(st$train_ids, colnames(e)); te <- match(st$test_ids, colnames(e)); genes <- st$screen
  set.seed(42L + st$rep * 1000L + st$fold)
  xs <- t(e[genes, tr]); colnames(xs) <- paste0("G", seq_along(genes)); y <- ph$status[tr]
  bo <- suppressWarnings(run_boruta(xs, y, P$boruta_runs, P$boruta_ntree)); bsel <- genes[match(bo$selected, colnames(xs))]
  la <- st$lasso; sv <- st$svmrfe
  votes <- table(factor(c(la, sv, bsel), levels = genes)); inter <- Reduce(intersect, list(la, sv, bsel))
  if (length(inter)) { feat <- inter; lev <- "L0_intersection" } else { two <- names(votes)[votes >= 2]
    if (length(two)) { feat <- two; lev <- "L1_>=2of3" } else { feat <- NULL; lev <- "L2/L3_not_computed" } }
  prob <- rep(NA_real_, length(te))
  if (!is.null(feat)) {
    x <- t(e[genes, tr]); xf <- x[, feat, drop = FALSE]; mu <- colMeans(xf); sdv <- apply(xf, 2, sd); sdv[sdv == 0] <- 1
    z <- sweep(sweep(xf, 2, mu), 2, sdv, "/"); zz <- if (ncol(z) == 1L) cbind(z, dummy0 = 0) else z
    rg <- cv.glmnet(zz, y, family = "binomial", alpha = 0, standardize = FALSE, foldid = stratified_foldid(y, P$inner_k), type.measure = "deviance")
    cf <- coef(rg, s = "lambda.min"); b0 <- as.numeric(cf[1]); b <- setNames(as.vector(cf)[-1][seq_along(feat)], feat)
    xt <- t(e[feat, te, drop = FALSE]); zt <- sweep(sweep(xt, 2, mu), 2, sdv, "/"); prob <- plogis(b0 + as.vector(zt %*% b[feat]))
  }
  list(kind = "fold", i = i, rep = st$rep, fold = st$fold, features = feat, level = lev, boruta = bsel, stored_boruta = st$boruta, stored_features = st$features,
       stored_level = st$level, pred = data.frame(rep = st$rep, id = st$test_ids, truth = as.character(ph$status[te]), prob = prob))
}
jobs <- c(as.list(c("fru", "ranger", "rf")), as.list(seq_along(stored)))
t0 <- Sys.time()
res <- mclapply(jobs, function(j) if (is.character(j)) final_fit(j) else fold_fit(j), mc.cores = 8, mc.preschedule = FALSE, mc.set.seed = FALSE)
cat("elapsed min:", format(as.numeric(difftime(Sys.time(), t0, units = "mins")), digits = 3), "\n")
bad <- vapply(res, function(z) !is.list(z), TRUE); cat("failed jobs:", sum(bad), "\n"); if (any(bad)) print(res[bad][1])
saveRDS(res, file.path(out_dir, "d5_backend_swap.rds"))
fin <- Filter(function(z) z$kind == "final", res); fo <- Filter(function(z) z$kind == "fold", res)
for (f in fin) cat(sprintf("FINAL fit, Boruta backend %s: panel = %s (%s); LASSO %d, SVM-RFE %d, Boruta %d\n", f$backend, paste(f$features, collapse = ","), f$level, length(f$lasso), length(f$svmrfe), length(f$boruta)))
jac <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))
J <- vapply(fo, function(z) jac(z$stored_boruta, z$boruta), 0); Jf <- vapply(fo, function(z) jac(z$stored_features, z$features), 0)
cat(sprintf("Boruta confirmed set, fru(production) vs ranger, 50 folds: Jaccard mean %.3f median %.3f min %.3f; sizes fru median %g ranger median %g\n", mean(J), median(J), min(J),
  median(vapply(fo, function(z) length(z$stored_boruta), 0L)), median(vapply(fo, function(z) length(z$boruta), 0L))))
cat(sprintf("features_k, production vs ranger-swap: identical in %d/50 folds; Jaccard mean %.3f; fallback level identical in %d/50; levels swap: L0 %d vs prod %d\n",
  sum(Jf == 1, na.rm = TRUE), mean(Jf, na.rm = TRUE), sum(vapply(fo, function(z) z$level == z$stored_level, TRUE)),
  sum(vapply(fo, function(z) z$level == "L0_intersection", TRUE)), sum(vapply(fo, function(z) z$stored_level == "L0_intersection", TRUE))))
oof <- do.call(rbind, lapply(fo, `[[`, "pred")); ok <- !anyNA(oof$prob)
if (ok) { auc_rep <- vapply(1:5, function(r) { o <- oof[oof$rep == r, ]; fast_auc(o$prob, factor(o$truth, c("ND", "T2D"))) }, 0)
  cat("ranger-swap nested AUC per repeat:", sprintf("%.3f", auc_rep), " mean", sprintf("%.3f", mean(auc_rep)), "(production 0.659)\n") } else cat("NA predictions present (L2/L3 folds):", sum(vapply(fo, function(z) z$level == "L2/L3_not_computed", TRUE)), "\n")
tab <- sort(table(unlist(lapply(fo, `[[`, "features"))) / 50, decreasing = TRUE); tab0 <- sort(table(unlist(lapply(fo, `[[`, "stored_features"))) / 50, decreasing = TRUE)
cat("top-10 stability, ranger-swap:", paste(sprintf("%s(%.2f)", names(head(tab, 10)), head(tab, 10)), collapse = " "), "\n")
cat("top-10 stability, production :", paste(sprintf("%s(%.2f)", names(head(tab0, 10)), head(tab0, 10)), collapse = " "), "\n")

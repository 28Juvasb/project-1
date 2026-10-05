# Diagnostic (not part of the pipeline): Boruta importance backends compared on the same training splits.
# Backends: fru-Z (production default), fru-Z rerun with a different seed (noise floor), ranger permutation Z (getImpRfZ),
# randomForest Z (getImpLegacyRfZ). maxRuns=500, ntree=1000 for all.
options(ncv.no_main = TRUE)
.libPaths(c(path.expand("~/R/lib"), .libPaths()))
library(here); library(parallel)
source(here("scripts/06_nested_cv.R"))
out_dir <- here("results/stage3/diagnostics"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
geneset <- readRDS(here("genesets/senescence_union.rds"))
e <- readRDS(here("data/processed/04_GSE41762_expr.rds")); p <- readRDS(here("data/processed/04_GSE41762_pheno.rds"))
ph <- prep_pheno(p, "islet"); P <- default_params("islet")
tasks <- make_outer_folds(ph$status, P$outer_k, P$outer_reps, 42L)
stored <- readRDS(here("results/stage3/06_per_fold_features_islet.rds"))
sel <- as.integer(commandArgs(trailingOnly = TRUE)); if (!length(sel)) sel <- c(1L, 11L, 23L, 37L, 46L)
jac <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))
one <- function(i) {
  tk <- tasks[[i]]; test <- tk$test; train <- setdiff(seq_len(ncol(e)), test)
  sc <- limma_screen(e[, train], ph[train, ], geneset, P, colnames(e)[train])
  genes <- sc$screen; xs <- t(e[genes, train]); colnames(xs) <- paste0("G", seq_along(genes)); y <- ph$status[train]
  runb <- function(imp, seed, nt = 1) {
    set.seed(seed); t0 <- proc.time()[3]
    b <- suppressWarnings(Boruta::Boruta(as.data.frame(xs), y, maxRuns = 500, ntree = 1000, getImp = imp, num.threads = nt))
    list(conf = genes[match(names(b$finalDecision)[b$finalDecision == "Confirmed"], colnames(xs))],
         tent = genes[match(names(b$finalDecision)[b$finalDecision == "Tentative"], colnames(xs))], t = proc.time()[3] - t0,
         runs = nrow(b$ImpHistory))
  }
  st <- stored[[i]]
  r <- list(task = i, screen = genes, stored_boruta = st$boruta, lasso = st$lasso, svmrfe = st$svmrfe, stored_features = st$features,
            stored_level = st$level,
            fru_a = runb(Boruta::getImpFruZ, 11), fru_b = runb(Boruta::getImpFruZ, 22),
            ranger = runb(Boruta::getImpRfZ, 11), rf = runb(Boruta::getImpLegacyRfZ, 11))
  r
}
res <- mclapply(sel, one, mc.cores = min(8L, length(sel)), mc.set.seed = FALSE)
stopifnot(!any(vapply(res, function(z) inherits(z, "try-error"), TRUE)))
saveRDS(res, file.path(out_dir, "d2_backend_compare.rds"))
tab <- do.call(rbind, lapply(res, function(r) {
  inter <- function(bo) Reduce(intersect, list(r$lasso, r$svmrfe, bo))
  data.frame(task = r$task, screen = length(r$screen), lasso = length(r$lasso), svm = length(r$svmrfe),
    n_stored = length(r$stored_boruta), n_fruA = length(r$fru_a$conf), n_fruB = length(r$fru_b$conf), n_ranger = length(r$ranger$conf), n_rf = length(r$rf$conf),
    tent_fru = length(r$fru_a$tent), tent_ranger = length(r$ranger$tent), tent_rf = length(r$rf$tent),
    J_stored_fruA = jac(r$stored_boruta, r$fru_a$conf), J_fruA_fruB = jac(r$fru_a$conf, r$fru_b$conf),
    J_fru_ranger = jac(r$fru_a$conf, r$ranger$conf), J_fru_rf = jac(r$fru_a$conf, r$rf$conf), J_ranger_rf = jac(r$ranger$conf, r$rf$conf),
    lasso_in_boruta_fru = length(intersect(r$lasso, r$stored_boruta)), lasso_in_ranger = length(intersect(r$lasso, r$ranger$conf)),
    svm_in_boruta_fru = length(intersect(r$svmrfe, r$stored_boruta)), svm_in_ranger = length(intersect(r$svmrfe, r$ranger$conf)),
    feat_fru = paste(inter(r$stored_boruta), collapse = "/"), feat_ranger = paste(inter(r$ranger$conf), collapse = "/"), feat_rf = paste(inter(r$rf$conf), collapse = "/"),
    sec_fru = round(r$fru_a$t), sec_ranger = round(r$ranger$t), sec_rf = round(r$rf$t))
}))
options(width = 250); print(tab); write.csv(tab, file.path(out_dir, "d2_backend_compare.csv"), row.names = FALSE)

# Diagnostic (not part of the pipeline): where do NaN importances come from in Boruta/fru?
options(ncv.no_main = TRUE)
.libPaths(c(path.expand("~/R/lib"), .libPaths()))
library(here)
source(here("scripts/06_nested_cv.R"))
out_dir <- here("results/stage3/diagnostics"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
geneset <- readRDS(here("genesets/senescence_union.rds"))
e <- readRDS(here("data/processed/04_GSE41762_expr.rds")); p <- readRDS(here("data/processed/04_GSE41762_pheno.rds"))
ph <- prep_pheno(p, "islet"); P <- default_params("islet")
tasks <- make_outer_folds(ph$status, P$outer_k, P$outer_reps, 42L)
stored <- readRDS(here("results/stage3/06_per_fold_features_islet.rds"))
sel <- c(1L, 11L, 23L, 37L, 46L)                      # spread over repeats; none touches validation data
rec <- new.env(); rec$rows <- list()
mk_imp <- function(fold_tag) function(x, y, ntree = 500, num.threads = 0, ...) {
  m <- fru::fru(x, y, trees = ntree, importance = TRUE, threads = num.threads)
  sc <- fru::importance(m, scale = TRUE); rw <- fru::importance(m, scale = FALSE)
  isn <- is.nan(sc) | is.na(sc); sh <- grepl("^shadow", names(sc))
  rec$rows[[length(rec$rows) + 1L]] <- data.frame(tag = fold_tag, n_attr = length(sc), n_real = sum(!sh), n_shadow = sum(sh),
     nan_real = sum(isn & !sh), nan_shadow = sum(isn & sh), raw_zero_real = sum(rw == 0 & !sh), raw_zero_nan = sum(rw == 0 & isn),
     nan_with_nonzero_raw = sum(isn & rw != 0),
     nan_names = paste(names(sc)[isn & !sh], collapse = ","))
  sc
}
res <- list()
for (i in sel) {
  tk <- tasks[[i]]; test <- tk$test; train <- setdiff(seq_len(ncol(e)), test)
  sc <- limma_screen(e[, train], ph[train, ], geneset, P, colnames(e)[train])
  genes <- sc$screen; xs <- t(e[genes, train]); colnames(xs) <- paste0("G", seq_along(genes)); y <- ph$status[train]
  # data checks
  v <- apply(xs, 2, var); 
  set.seed(1000 + i)
  b <- suppressWarnings(Boruta::Boruta(as.data.frame(xs), y, maxRuns = 500, ntree = 1000, num.threads = 1, getImp = mk_imp(i)))
  df <- do.call(rbind, rec$rows); df <- df[df$tag == i, ]
  cat(sprintf("fold task %d (rep %d fold %d): screen=%d genes, min var=%.4g, dup columns=%d | iterations=%d, iterations with any NaN=%d (%.1f%%), NaN real-gene cells=%d, NaN shadow cells=%d, NaN with nonzero raw=%d\n",
      i, tk$rep, tk$fold, length(genes), min(v), sum(duplicated(as.data.frame(xs), MARGIN = 2)), nrow(df), sum(df$nan_real + df$nan_shadow > 0),
      100 * mean(df$nan_real + df$nan_shadow > 0), sum(df$nan_real), sum(df$nan_shadow), sum(df$nan_with_nonzero_raw)))
  cat("   median n_attr per iteration", median(df$n_attr), "; iterations NaN-in-real", sum(df$nan_real > 0), "; median/max NaN per iteration (real)", median(df$nan_real), max(df$nan_real),
      "; median raw-zero real", median(df$raw_zero_real), "\n")
  tb <- table(unlist(strsplit(df$nan_names[df$nan_names != ""], ",")))
  cat("   distinct real genes ever NaN:", length(tb), "; top:", paste(head(sort(tb, TRUE), 5), names(head(sort(tb, TRUE), 5)), collapse = " "), "\n")
  res[[as.character(i)]] <- list(df = df, boruta = names(b$finalDecision)[b$finalDecision == "Confirmed"], nan_tab = tb, genes = genes,
                                 stored = stored[[i]]$boruta, status = table(b$finalDecision))
}
saveRDS(res, file.path(out_dir, "d1_nan_summary.rds"))

# Diagnostic: anatomy of NaN cells from fru scaled importance (one fold, full Boruta run, record raw and scaled values of NaN cells)
options(ncv.no_main = TRUE)
.libPaths(c(path.expand("~/R/lib"), .libPaths()))
library(here); source(here("scripts/06_nested_cv.R"))
out_dir <- here("results/stage3/diagnostics")
geneset <- readRDS(here("genesets/senescence_union.rds"))
e <- readRDS(here("data/processed/04_GSE41762_expr.rds")); p <- readRDS(here("data/processed/04_GSE41762_pheno.rds"))
ph <- prep_pheno(p, "islet"); P <- default_params("islet")
tasks <- make_outer_folds(ph$status, P$outer_k, P$outer_reps, 42L)
i <- 1L; tk <- tasks[[i]]; train <- setdiff(seq_len(ncol(e)), tk$test)
sc <- limma_screen(e[, train], ph[train, ], geneset, P, colnames(e)[train])
genes <- sc$screen; xs <- t(e[genes, train]); colnames(xs) <- paste0("G", seq_along(genes)); y <- ph$status[train]
cells <- list(); it <- 0L; allrows <- list()
imp <- function(x, y, ntree = 500, num.threads = 0, ...) {
  it <<- it + 1L
  m <- fru::fru(x, y, trees = ntree, importance = TRUE, threads = num.threads)
  s <- fru::importance(m, scale = TRUE); r <- fru::importance(m, scale = FALSE)
  bad <- which(is.na(s))
  allrows[[it]] <<- data.frame(it = it, n_attr = length(s), n_nan = length(bad), n_nan_rawNaN = sum(is.na(r[bad])), n_nan_rawZero = sum(r[bad] == 0, na.rm = TRUE),
                               n_nan_rawPos = sum(r[bad] > 0, na.rm = TRUE), n_nan_rawNeg = sum(r[bad] < 0, na.rm = TRUE), n_raw_exact0 = sum(r == 0, na.rm = TRUE),
                               n_raw_NaN = sum(is.na(r)), tries = round(sqrt(ncol(x))))
  if (length(bad)) cells[[length(cells) + 1L]] <<- data.frame(it = it, n_attr = length(s), var = names(s)[bad], raw = r[bad], scaled = s[bad],
       var_x = vapply(names(s)[bad], function(n) var(as.numeric(x[[n]])), 0), n_unique = vapply(names(s)[bad], function(n) length(unique(x[[n]])), 0L))
  s
}
set.seed(1001)
b <- suppressWarnings(Boruta::Boruta(as.data.frame(xs), y, maxRuns = 500, ntree = 1000, num.threads = 1, getImp = imp))
A <- do.call(rbind, allrows); C <- do.call(rbind, cells)
saveRDS(list(A = A, C = C), file.path(out_dir, "d3_nan_anatomy.rds"))
cat("iterations:", nrow(A), " with NaN:", sum(A$n_nan > 0), "\n")
cat("NaN cells:", sum(A$n_nan), " of which raw NaN:", sum(A$n_nan_rawNaN), " raw==0:", sum(A$n_nan_rawZero), " raw>0:", sum(A$n_nan_rawPos), " raw<0:", sum(A$n_nan_rawNeg), "\n")
cat("n_attr in NaN iterations: "); print(summary(A$n_attr[A$n_nan > 0])); cat("n_attr overall: "); print(summary(A$n_attr))
cat("iteration index of NaN events:", A$it[A$n_nan > 0], "\n")
cat("attributes with exactly zero raw importance per iteration (all attrs): "); print(summary(A$n_raw_exact0))
print(head(C, 15)); cat("NaN-cell shadow vs real:", table(grepl("^shadow", C$var)), "\n")
# rate of NaN as function of n_attr
A$bin <- cut(A$n_attr, c(0, 20, 40, 80, 160, 400)); print(aggregate(cbind(nan_iter = n_nan > 0, n_nan = n_nan) ~ bin, A, function(z) c(sum = sum(z), n = length(z))))

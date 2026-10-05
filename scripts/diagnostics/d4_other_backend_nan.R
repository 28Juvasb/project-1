# Diagnostic: do the ranger / randomForest backends emit NaN importances on the same training split (task 1)? Also fru raw (unscaled).
options(ncv.no_main = TRUE)
.libPaths(c(path.expand("~/R/lib"), .libPaths()))
library(here); library(parallel); source(here("scripts/06_nested_cv.R"))
out_dir <- here("results/stage3/diagnostics")
geneset <- readRDS(here("genesets/senescence_union.rds"))
e <- readRDS(here("data/processed/04_GSE41762_expr.rds")); p <- readRDS(here("data/processed/04_GSE41762_pheno.rds"))
ph <- prep_pheno(p, "islet"); P <- default_params("islet")
tasks <- make_outer_folds(ph$status, P$outer_k, P$outer_reps, 42L)
i <- 1L; train <- setdiff(seq_len(ncol(e)), tasks[[i]]$test)
sc <- limma_screen(e[, train], ph[train, ], geneset, P, colnames(e)[train])
genes <- sc$screen; xs <- t(e[genes, train]); colnames(xs) <- paste0("G", seq_along(genes)); y <- ph$status[train]
cnt <- function(f, nm) { env <- new.env(); env$n <- 0L; env$it <- 0L
  g <- function(x, y, ntree = 500, ...) { r <- f(x, y, ntree = ntree, ...); env$it <- env$it + 1L; env$n <- env$n + sum(is.na(r)); r }
  list(g = g, env = env) }
bk <- list(ranger = Boruta::getImpRfZ, rf = Boruta::getImpLegacyRfZ, fruraw = Boruta::getImpFruRaw)
r <- mclapply(names(bk), function(nm) { w <- cnt(bk[[nm]], nm); set.seed(5)
  nw <- 0L; b <- withCallingHandlers(Boruta::Boruta(as.data.frame(xs), y, maxRuns = 500, ntree = 1000, getImp = w$g, num.threads = 1),
        warning = function(x) { nw <<- nw + 1L; invokeRestart("muffleWarning") })
  list(nm = nm, nan_cells = w$env$n, warns = nw, conf = length(which(b$finalDecision == "Confirmed")), tent = sum(b$finalDecision == "Tentative"), runs = nrow(b$ImpHistory)) },
  mc.cores = 3, mc.set.seed = FALSE)
for (z in r) cat(sprintf("%s: NaN cells=%d, Boruta NaN warnings=%d, Confirmed=%d, Tentative=%d, iterations=%d\n", z$nm, z$nan_cells, z$warns, z$conf, z$tent, z$runs))

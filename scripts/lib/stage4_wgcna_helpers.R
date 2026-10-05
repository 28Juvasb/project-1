# stage4_wgcna_helpers.R -- helpers for scripts/07_wgcna.R (protocol v1.1 §4.1). Sourced, not run.

# Save a base-graphics figure as both vector PDF and PNG (CLAUDE.md figure convention).
save_fig <- function(dir, name, fn, width = 8, height = 6, res = 150) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  pdf(file.path(dir, paste0(name, ".pdf")), width = width, height = height); fn(); dev.off()
  png(file.path(dir, paste0(name, ".png")), width = width, height = height, units = "in", res = res); fn(); dev.off()
  invisible(file.path(dir, paste0(name, c(".pdf", ".png"))))
}

# Fixed WGCNA network parameters (protocol §4.1) -- one place so reference run and resamples cannot diverge.
wgcna_params <- function() list(networkType = "signed hybrid", TOMType = "signed", corType = "bicor", maxPOutliers = 0.05,
                                minModuleSize = 30L, mergeCutHeight = 0.25, deepSplit = 2L, maxBlockSize = 6000L)

run_blockwise <- function(datExpr, beta, nThreads = 1L) {
  P <- wgcna_params()
  WGCNA::blockwiseModules(datExpr, power = beta, networkType = P$networkType, TOMType = P$TOMType, corType = P$corType, maxPOutliers = P$maxPOutliers,
                          minModuleSize = P$minModuleSize, mergeCutHeight = P$mergeCutHeight, deepSplit = P$deepSplit, maxBlockSize = P$maxBlockSize,
                          numericLabels = TRUE, saveTOMs = FALSE, nThreads = nThreads, verbose = 0)
}

# ADIPOSE STABILITY (protocol §4.1): B stratified bootstrap resamples of SAMPLES. For each resample: re-run module detection (same parameters, same beta,
# same 5000 genes), then modulePreservation of every reference module (reference = full-cohort modules; test = the resample, with the resample's own
# module labels supplied so the overlap statistics are also computed). A reference module is "preserved" in a resample if Zsummary.pres > zcut.
# Parallelism: one forked worker per resample (WGCNA internal threads = 1 inside workers).
run_stability <- function(datExpr, status, ref_labels, beta, B = 100L, nPerm = 200L, zcut = 2, cores = 8L, seed = 42L) {
  stopifnot(nrow(datExpr) == length(status), length(ref_labels) == ncol(datExpr), is.factor(status))
  P <- wgcna_params()
  one <- function(b) {
    WGCNA::disableWGCNAThreads(); set.seed(seed + b)
    idx <- unlist(lapply(levels(status), function(l) { i <- which(status == l); i[sample.int(length(i), length(i), replace = TRUE)] }))
    datB <- datExpr[idx, , drop = FALSE]; rownames(datB) <- paste0("b", b, "_", seq_along(idx), "_", rownames(datExpr)[idx])
    t0 <- Sys.time()
    r <- tryCatch({
      netB <- run_blockwise(datB, beta, nThreads = 1L)
      mp <- WGCNA::modulePreservation(list(ref = list(data = datExpr), test = list(data = datB)),
                                      multiColor = list(ref = ref_labels, test = netB$colors), referenceNetworks = 1, testNetworks = list(2),
                                      nPermutations = nPerm, networkType = P$networkType, corFnc = "bicor", corOptions = paste0("use = 'p', maxPOutliers = ", P$maxPOutliers),
                                      randomSeed = seed + b, verbose = 0, savePermutedStatistics = FALSE, parallelCalculation = FALSE)
      Z <- mp$preservation$Z[[1]][[2]]
      list(ok = TRUE, Z = data.frame(module = rownames(Z), Zsummary = Z$Zsummary.pres, Zdensity = Z$Zdensity.pres, Zconnectivity = Z$Zconnectivity.pres,
                                     stringsAsFactors = FALSE), n_modules_boot = length(setdiff(unique(netB$colors), 0L)))
    }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
    r$b <- b; r$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs")); r
  }
  res <- parallel::mclapply(seq_len(B), one, mc.cores = cores, mc.preschedule = FALSE)
  ok <- vapply(res, function(r) is.list(r) && isTRUE(r$ok), NA)
  list(results = res, ok = ok)
}

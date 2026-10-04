# 00_setup.R -- verify every package loads, log sessionInfo(), snapshot renv.lock.
# Stage 0. Idempotent: rerunning rewrites the log and refreshes the lockfile.

set.seed(42)
library(here)

script_name <- "00_setup"
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create(here("logs"), showWarnings = FALSE, recursive = TRUE)
log_path <- here("logs", sprintf("%s_%s.log", script_name, ts))
con <- file(log_path, open = "wt")
sink(con, split = TRUE)
sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

cat("Script:", script_name, "\nTimestamp:", ts, "\nInputs: none (package check only)\n")
cat("Parameters: set.seed(42)\n\n")

pkgs <- c(
  "oligo", "affy", "limma", "sva", "WGCNA", "glmnet", "caret", "e1071",
  "randomForest", "Boruta", "pROC", "clusterProfiler", "org.Hs.eg.db", "here",
  "data.table", "ggplot2", "ComplexHeatmap",
  # additions required by later stages
  "GEOquery", "affyPLM", "Biobase", "impute", "preprocessCore", "UpSetR",
  "rmarkdown", "knitr", "enrichplot", "pd.hugene.1.0.st.v1",
  "hugene10sttranscriptcluster.db", "hgu133plus2.db", "hgu133plus2cdf",
  "pd.u133.x3p", "u133x3p.db", "u133x3pcdf"
)

res <- vapply(pkgs, function(p) {
  ok <- suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
  if (ok) as.character(utils::packageVersion(p)) else NA_character_
}, character(1))
status <- data.frame(package = pkgs, version = unname(res), loaded = !is.na(res))
print(status, row.names = FALSE)

stopifnot("R >= 4.4 required" = getRversion() >= "4.4.0")
bioc <- as.character(BiocManager::version())
cat("\nR:", R.version.string, "\nBioconductor:", bioc, "\n")
stopifnot("Bioconductor >= 3.20 required" = package_version(bioc) >= "3.20")

dir.create(here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
renv::settings$snapshot.type("all")
renv::snapshot(type = "all", prompt = FALSE)
cat("\nrenv.lock written:", file.exists(here("renv.lock")), "\n\n")
print(sessionInfo())

# Assertions: fail loudly
missing <- status$package[!status$loaded]
if (length(missing)) stop("Packages failed to load: ", paste(missing, collapse = ", "))
stopifnot(file.exists(here("renv.lock")))
cat("\nALL ASSERTIONS PASSED\n")

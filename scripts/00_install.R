# 00_install.R -- one-time package installation into the renv project library.
# Run once, then 00_setup.R snapshots renv.lock. Never update packages mid-project (CLAUDE.md).
# Usage: Rscript scripts/00_install.R   (from the repo root)

options(repos = c(CRAN = "https://cloud.r-project.org"))
stopifnot(requireNamespace("renv", quietly = TRUE))

cran_pkgs <- c(
  "Boruta", "here", "data.table", "ggplot2", "glmnet", "caret", "e1071", "randomForest",
  "pROC", "WGCNA", "UpSetR", "rmarkdown", "knitr", "pheatmap",
  "ggrepel", "VennDiagram", "readr", "jsonlite", "digest"
)

bioc_pkgs <- c(
  "oligo", "affy", "affyPLM", "limma", "sva", "GEOquery", "Biobase",
  "clusterProfiler", "org.Hs.eg.db", "ComplexHeatmap", "enrichplot",
  "impute", "preprocessCore", "GO.db",
  # platform annotation / probe-level packages (GPL6244, GPL570, GPL1352)
  "pd.hugene.1.0.st.v1", "hugene10sttranscriptcluster.db",
  "hgu133plus2.db", "hgu133plus2cdf",
  "pd.u133.x3p", "u133x3p.db", "u133x3pcdf"
)

renv::install(c(cran_pkgs, paste0("bioc::", bioc_pkgs)), prompt = FALSE)
cat("Install finished.\n")

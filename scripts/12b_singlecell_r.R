# 12b_singlecell_r.R -- R helper for scripts/12_singlecell.py (called via subprocess; do not run stand-alone out of order).
#
# mode "map"   : Entrez gene ID -> HGNC symbol table from org.Hs.eg.db (GSE81608 rows are Entrez IDs, not symbols).
#                Output: results/stage7/12_entrez_symbol_map.csv
# mode "limma" : pseudobulk T2D-vs-ND test per cell type with limma lmFit + eBayes (protocol §5.3; INPUT IS mean log2(RPKM+1), NOT
#                summed counts -- no raw counts exist for GSE81608; deviation logged by 12_singlecell.py).
#                Inputs : results/stage7/12_pseudobulk_meta.csv, results/stage7/12_pseudobulk_<celltype>.csv (genes x pseudobulk samples)
#                Outputs: results/stage7/12_limma_allgenes.csv.gz, results/stage7/12_limma_model_info.csv
# Models per cell type (ND = reference; logFC = T2D - ND):
#   unadjusted: ~ group
#   adjusted  : ~ group + age + sex   (only if residual df >= 8 after adding the covariates and both sexes are present)
#   primary   = adjusted when feasible, else unadjusted (rule fixed before looking at any result).
# Cell types with < 2 donors in either group are not testable. No per-cell tests anywhere.
set.seed(42)
suppressPackageStartupMessages(library(here))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L, args[1] %in% c("map", "limma"))
mode <- args[1]
out_dir <- here("results", "stage7"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (mode == "map") {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  ids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
  sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  tab <- data.frame(entrez = as.integer(ids), symbol = unname(sym), stringsAsFactors = FALSE)
  tab <- tab[!is.na(tab$symbol), ]
  write.csv(tab, file.path(out_dir, "12_entrez_symbol_map.csv"), row.names = FALSE, quote = TRUE)
  cat("org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db")), " entrez->symbol rows:", nrow(tab), "\n")
  quit(save = "no", status = 0L)
}

suppressPackageStartupMessages(library(limma))
cat("limma", as.character(packageVersion("limma")), "; R", R.version.string, "\n")
meta <- read.csv(file.path(out_dir, "12_pseudobulk_meta.csv"), stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("sample", "celltype", "donor", "group", "age", "sex", "n_cells") %in% names(meta)))
MIN_DETECT_DONORS <- 3L   # gene kept in a cell type's fit if pseudobulk value > 0 in >= 3 donors
MIN_RESID_DF_ADJ  <- 8L
all_res <- list(); info <- list()
for (ct in unique(meta$celltype)) {
  m <- meta[meta$celltype == ct, ]
  f <- file.path(out_dir, sprintf("12_pseudobulk_%s.csv", ct))
  mat <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  stopifnot(setequal(colnames(mat), m$sample), !anyNA(mat))
  mat <- mat[, m$sample, drop = FALSE]
  n_nd <- sum(m$group == "ND"); n_t2d <- sum(m$group == "T2D"); n <- nrow(m)
  row <- data.frame(celltype = ct, n_donors = n, n_ND = n_nd, n_T2D = n_t2d, n_genes_tested = NA_integer_,
                    adjusted_feasible = FALSE, resid_df_adjusted = NA_integer_, resid_df_unadjusted = NA_integer_, status = "", stringsAsFactors = FALSE)
  if (min(n_nd, n_t2d) < 2L) {
    row$status <- sprintf("NOT TESTABLE: %d ND / %d T2D donors with >= 10 cells (need >= 2 per group)", n_nd, n_t2d)
    info[[ct]] <- row; cat(ct, ":", row$status, "\n"); next
  }
  keep <- rowSums(mat > 0) >= MIN_DETECT_DONORS
  y <- mat[keep, , drop = FALSE]; row$n_genes_tested <- nrow(y)
  grp <- factor(m$group, levels = c("ND", "T2D")); sex <- factor(m$sex); age <- as.numeric(m$age)
  d0 <- model.matrix(~ grp); row$resid_df_unadjusted <- n - ncol(d0)
  adj_ok <- nlevels(sex) == 2L && (n - 4L) >= MIN_RESID_DF_ADJ
  models <- list(unadjusted = d0)
  if (adj_ok) { models$adjusted <- model.matrix(~ grp + age + sex); row$adjusted_feasible <- TRUE; row$resid_df_adjusted <- n - ncol(models$adjusted) }
  row$status <- "tested"; info[[ct]] <- row
  for (mn in names(models)) {
    des <- models[[mn]]; stopifnot(qr(des)$rank == ncol(des))
    fit <- eBayes(lmFit(y, des))
    tt <- topTable(fit, coef = "grpT2D", number = Inf, sort.by = "none")
    all_res[[paste(ct, mn)]] <- data.frame(celltype = ct, model = mn, gene = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
                                           P.Value = tt$P.Value, adj.P.Val_allgenes = tt$adj.P.Val, sigma_post = sqrt(fit$s2.post),
                                           df_total = fit$df.residual[1] + fit$df.prior, stringsAsFactors = FALSE)
    cat(sprintf("%s [%s]: %d ND / %d T2D donors, %d genes fitted, residual df %d, genes FDR<0.05 (all-gene BH) = %d\n",
                ct, mn, n_nd, n_t2d, nrow(y), fit$df.residual[1], sum(tt$adj.P.Val < 0.05)))
  }
}
res <- do.call(rbind, all_res); info <- do.call(rbind, info)
stopifnot(!is.null(res), nrow(res) > 0, !anyNA(res$P.Value))
gz <- gzfile(file.path(out_dir, "12_limma_allgenes.csv.gz"), "wt"); write.csv(res, gz, row.names = FALSE); close(gz)
write.csv(info, file.path(out_dir, "12_limma_model_info.csv"), row.names = FALSE)
cat("\nsessionInfo():\n"); print(sessionInfo())

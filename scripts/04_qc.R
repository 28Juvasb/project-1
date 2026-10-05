# 04_qc.R -- array quality control, outlier flagging and sex check (protocol v1.1 §2.3-2.4).
#
# Inputs : data/processed/03_<GSE>_{expr,pheno,expressed}.rds (from 03_normalize.R),
#          CEL files via data/processed/cel_cache/<GSE>/ (extracted read-only from data/raw tars),
#          logs/excluded_samples.csv
# Outputs: data/processed/04_<GSE>_expr.rds      post-QC matrix (= 03 matrix unless the NUSE rule excluded arrays,
#                                                in which case the remaining arrays are re-normalised from raw)
#          data/processed/04_<GSE>_pheno.rds     pheno + QC columns + inferred sex (GSE27951: sex covariate = inferred)
#          data/processed/04_<GSE>_expressed.rds expression-filter flags on the post-QC matrix
#          results/stage1/qc_summary.csv, results/stage1/04_qc.rds,
#          results/stage1/figures/04_<GSE>_{rle_nuse,pca,hclust,sexcheck}.{pdf,png},
#          logs/excluded_samples.csv (NUSE exclusions only), logs/04_qc_<timestamp>.log
#
# Rules: the ONLY automatic exclusion is median NUSE > 1.05 (pre-specified, objective).
#        PCA outliers (> 3 SD on PC1 or PC2) and sex mismatches are FLAGGED, never excluded here (project lead decides).
# VALIDATION SETS (GSE76894, GSE71416, GSE20966, GSE16415; CLAUDE.md rule 2): QC is label-blind. Figures are coloured by
#        batch / sex, never by disease group, and no ND-vs-T2D statistic is computed. Discovery sets may be coloured by group.

# oligo::fitProbeLevelModel segfaults inside multithreaded OpenBLAS (dgemm_kernel_ZEN, called from preprocessCore
# PLM threads) on this machine. OpenBLAS reads its thread count at load time, so the script re-executes itself
# with OPENBLAS_NUM_THREADS=1 (numerical results are unaffected; only BLAS threading changes).
if (Sys.getenv("OPENBLAS_NUM_THREADS") != "1") {
  self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(self) != 1) stop("Run 04_qc.R via Rscript with OPENBLAS_NUM_THREADS=1")
  st <- system2(file.path(R.home("bin"), "Rscript"), shQuote(self), env = "OPENBLAS_NUM_THREADS=1")
  quit(save = "no", status = st)
}

set.seed(42)
suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2) })
source(here("scripts", "lib", "stage1_helpers.R"))

script_name <- "04_qc"
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
fig_dir <- here("results", "stage1", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt")
sink(con, split = TRUE); sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

cat("Script:", script_name, " run:", ts, "\n\nInputs:\n")
for (gse in STAGE1$datasets) for (k in c("expr", "pheno")) {
  f <- here("data", "processed", sprintf("03_%s_%s.rds", gse, k))
  cat("  ", sub(paste0(here(), "/"), "", f), " md5=", tools::md5sum(f), "\n", sep = "")
}
cat(sprintf("\nParameters: NUSE exclusion median > %.2f | PCA top %d genes by variance, prcomp centred, unscaled | flag |z| > %d on PC1 or PC2\n",
            STAGE1$nuse_max, STAGE1$pca_top, STAGE1$pca_sd))
cat("  hclust: same top genes, distance 1 - Pearson r, average linkage\n")
cat("  sex: XIST vs mean(", paste(STAGE1$y_genes, collapse = ", "), " where present), k-means k=2 on z-scored (XIST, Y-mean), nstart=25, seed 42;\n",
    "       cluster with higher (Y - XIST) centroid = male\n", sep = "")
cat("  PLM: affyPLM::fitPLM (affy) / oligo::fitProbeLevelModel(target='core') (oligo); GSE16415: RLE from normalised gene matrix, no NUSE\n\n")

save_gg <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
}
save_base <- function(fun, name, w, h) {
  pdf(file.path(fig_dir, paste0(name, ".pdf")), width = w, height = h); fun(); dev.off()
  png(file.path(fig_dir, paste0(name, ".png")), width = w, height = h, units = "in", res = 150); fun(); dev.off()
}
pal <- function(x) {
  lv <- sort(unique(as.character(x)))
  cols <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")[seq_along(lv)], lv)
  cols
}

# colour variable used in figures: group for discovery; batch/sex for validation (never group)
colour_var <- function(gse, p) {
  if (STAGE1$role[gse] == "discovery") return(list(name = "group", values = p$group))
  if (gse == "GSE76894") return(list(name = "batch", values = paste0("batch ", p$batch)))
  if (gse %in% c("GSE71416", "GSE20966")) return(list(name = "recorded sex", values = p$sex))
  list(name = "inferred sex", values = p$inferred_sex)
}

qc_all <- list(); summary_rows <- list(); excl_new <- list()

for (gse in STAGE1$datasets) {
  cat(sprintf("\n==================== %s (%s, %s) ====================\n", gse, STAGE1$role[gse], STAGE1$platform[gse]))
  pheno <- readRDS(here("data", "processed", sprintf("03_%s_pheno.rds", gse)))
  expr  <- readRDS(here("data", "processed", sprintf("03_%s_expr.rds", gse)))
  stopifnot(identical(colnames(expr), pheno$gsm))
  n_in <- nrow(pheno)

  # ---------------- RLE / NUSE (§2.3) ----------------
  if (STAGE1$method[gse] == "matrix") {
    rle_v <- expr - apply(expr, 1, median)
    nuse_v <- NULL
  } else {
    raw <- read_raw(gse, pheno$gsm)
    if (STAGE1$method[gse] == "oligo") {
      plm <- oligo::fitProbeLevelModel(raw, target = "core", verbose = FALSE)
      rle_v <- oligo::RLE(plm, type = "values"); nuse_v <- oligo::NUSE(plm, type = "values")
    } else {
      plm <- affyPLM::fitPLM(raw, verbosity.level = 0,
                             output.param = list(weights = FALSE, residuals = FALSE, varcov = "none", resid.SE = TRUE))
      rle_v <- affyPLM::RLE(plm, type = "values"); nuse_v <- affyPLM::NUSE(plm, type = "values")
    }
    rm(raw, plm); gc(verbose = FALSE)
    colnames(rle_v) <- colnames(nuse_v) <- pheno$gsm
  }
  colnames(rle_v) <- pheno$gsm
  rle_med <- apply(rle_v, 2, median, na.rm = TRUE)
  rle_iqr <- apply(rle_v, 2, IQR, na.rm = TRUE)
  nuse_med <- if (is.null(nuse_v)) setNames(rep(NA_real_, n_in), pheno$gsm) else apply(nuse_v, 2, median, na.rm = TRUE)
  nuse_excl <- !is.na(nuse_med) & nuse_med > STAGE1$nuse_max
  cat(sprintf("  median NUSE range: %s | RLE median range %.3f..%.3f | RLE IQR range %.3f..%.3f\n",
              if (is.null(nuse_v)) "n/a (no probe-level data)" else sprintf("%.4f..%.4f", min(nuse_med), max(nuse_med)),
              min(rle_med), max(rle_med), min(rle_iqr), max(rle_iqr)))
  cat(sprintf("  NUSE exclusions (median > %.2f): %d %s\n", STAGE1$nuse_max, sum(nuse_excl),
              if (any(nuse_excl)) paste0("(", paste(sprintf("%s=%.4f", pheno$gsm[nuse_excl], nuse_med[nuse_excl]), collapse = ", "), ")") else ""))

  cv <- colour_var(gse, within(pheno, inferred_sex <- NA))
  draw_rle_nuse <- function() {
    cvq <- colour_var(gse, within(pheno, inferred_sex <- NA))
    cols <- if (all(is.na(cvq$values))) rep("grey80", n_in) else pal(cvq$values)[as.character(cvq$values)]
    np <- if (is.null(nuse_v)) 1 else 2
    par(mfrow = c(np, 1), mar = c(6, 4, 3, 1), cex.axis = if (n_in > 50) 0.45 else 0.7)
    boxplot(rle_v, outline = FALSE, las = 2, col = cols, main = paste(gse, "RLE", if (is.null(nuse_v)) "(normalised gene matrix)" else "(PLM)"),
            ylab = "RLE"); abline(h = 0, lty = 2)
    if (!is.null(nuse_v)) {
      boxplot(nuse_v, outline = FALSE, las = 2, col = cols, main = paste(gse, "NUSE (PLM); dashed = 1.05 exclusion limit"),
              ylab = "NUSE", ylim = range(c(0.9, 1.1, apply(nuse_v, 2, quantile, c(.25, .75), na.rm = TRUE))))
      abline(h = c(1, STAGE1$nuse_max), lty = c(1, 2), col = c("grey40", "red"))
    }
    if (!all(is.na(cvq$values))) legend("topright", fill = pal(cvq$values), legend = names(pal(cvq$values)),
                                        title = cvq$name, cex = 0.7, bg = "white")
  }
  save_base(draw_rle_nuse, sprintf("04_%s_rle_nuse", gse), w = max(8, n_in * 0.12), h = if (is.null(nuse_v)) 5 else 9)
  rm(rle_v, nuse_v); gc(verbose = FALSE)

  # ---------------- apply NUSE exclusion; re-normalise remaining arrays if needed ----------
  if (any(nuse_excl)) {
    for (s in pheno$gsm[nuse_excl]) excl_new[[length(excl_new) + 1]] <- data.table(gsm = s, gse = gse,
      reason = sprintf("median NUSE %.4f > 1.05 (protocol 2.3)", nuse_med[s]), date = format(Sys.Date()))
    pheno_qc <- pheno[!nuse_excl, ]
    cat(sprintf("  re-normalising %d remaining arrays from raw (protocol 2.1 settings)\n", nrow(pheno_qc)))
    rn <- run_normalisation(gse, pheno_qc)
    expr <- rn$expr; ef <- rn$expressed
  } else {
    pheno_qc <- pheno
    ef <- readRDS(here("data", "processed", sprintf("03_%s_expressed.rds", gse)))
  }
  stopifnot(identical(colnames(expr), pheno_qc$gsm))

  # ---------------- sex check (§2.4) ----------------
  yg <- intersect(STAGE1$y_genes, rownames(expr))
  if (length(yg) == 0) stop(gse, ": no Y gene present; sex check impossible")
  has_xist <- "XIST" %in% rownames(expr)
  ysc <- colMeans(expr[yg, , drop = FALSE])
  xist <- if (has_xist) expr["XIST", ] else setNames(rep(NA_real_, ncol(expr)), colnames(expr))
  set.seed(42)
  if (has_xist) {
    sm <- scale(cbind(XIST = xist, Y = ysc))
    km <- kmeans(sm, centers = 2, nstart = 25)
    male_cl <- which.max(km$centers[, "Y"] - km$centers[, "XIST"])
  } else {
    # DEVIATION (logs/deviations.md): XIST is not annotated on this platform's protocol annotation source
    # (hugene10sttranscriptcluster.db for GPL6244; GPL2986 table for GSE16415). k-means k=2 on the Y-gene mean alone.
    cat("  *** XIST ABSENT from the collapsed matrix (not annotated on this platform's annotation source):\n",
        "      k-means k=2 on the z-scored Y-gene mean ONLY (deviation from §2.4, logged in logs/deviations.md) ***\n")
    sm <- scale(cbind(Y = ysc))
    km <- kmeans(sm, centers = 2, nstart = 25)
    male_cl <- which.max(km$centers[, "Y"])
  }
  inferred <- ifelse(km$cluster == male_cl, "male", "female")
  cat(sprintf("  Y-mean by inferred sex: female max %.2f, male min %.2f (gap %.2f log2)\n",
              max(ysc[inferred == "female"]), min(ysc[inferred == "male"]),
              min(ysc[inferred == "male"]) - max(ysc[inferred == "female"])))
  rec <- pheno_qc$sex
  mismatch <- ifelse(is.na(rec), NA, rec != inferred)
  cat(sprintf("  sex check genes: %s%s (Y genes absent: %s)\n", if (has_xist) "XIST + " else "NO XIST; Y only: ", paste(yg, collapse = ","),
              paste(setdiff(STAGE1$y_genes, yg), collapse = ",") ))
  cat("  k-means centres (raw scale):\n")
  print(round(cbind(XIST = tapply(xist, inferred, mean), Ymean = tapply(ysc, inferred, mean), n = table(inferred)), 3))
  if (all(is.na(rec))) {
    cat("  recorded sex: none in GEO -> inferred sex reported", if (gse == "GSE27951") "and used as the sex covariate (§2.4)" else "only", "\n")
  } else {
    print(table(recorded = rec, inferred = inferred, useNA = "ifany"))
    cat(sprintf("  sex mismatches: %d %s\n", sum(mismatch, na.rm = TRUE),
                if (any(mismatch, na.rm = TRUE)) paste0("(", paste(sprintf("%s rec=%s inf=%s", pheno_qc$gsm[which(mismatch)],
                  rec[which(mismatch)], inferred[which(mismatch)]), collapse = "; "), ")") else ""))
  }
  if (gse == "GSE41762") {
    i30 <- which(pheno_qc$donor_num == 30)
    cat(sprintf("  donor ID30 (%s): recorded GSE41762 sex = %s (GSE38642 deposit says male); expression-inferred = %s; XIST=%.2f Ymean=%.2f\n",
                pheno_qc$gsm[i30], rec[i30], inferred[i30], xist[i30], ysc[i30]))
  }
  pheno_qc$inferred_sex <- inferred
  pheno_qc$sex_mismatch <- mismatch
  pheno_qc$XIST <- xist; pheno_qc$Y_mean <- ysc
  if (gse == "GSE27951") { pheno_qc$sex <- inferred; pheno_qc$sex_source <- "expression-inferred (protocol 2.4; no sex field in GEO)"
  } else if (gse == "GSE16415") { pheno_qc$sex_source <- "none recorded; inferred_sex column is expression-based"
  } else pheno_qc$sex_source <- "recorded (GEO)"

  sd_df <- data.frame(XIST = if (has_xist) xist else rank(ysc, ties.method = "first"), Y = ysc, inferred = inferred,
                      recorded = ifelse(is.na(rec), "not recorded", rec), gsm = pheno_qc$gsm,
                      mismatch = !is.na(mismatch) & mismatch)
  p_sex <- ggplot(sd_df, aes(XIST, Y, colour = recorded, shape = inferred)) + geom_point(size = 2.5) +
    geom_text(data = sd_df[sd_df$mismatch, ], aes(label = gsm), vjust = -0.9, size = 3, show.legend = FALSE) +
    labs(title = sprintf("%s sex check (Y = mean of %s)", gse, paste(yg, collapse = ",")),
         subtitle = sprintf("k-means k=2%s; mismatches labelled (%d)",
                            if (has_xist) "" else " on Y-gene mean only (XIST not annotated on this platform)", sum(sd_df$mismatch)),
         x = if (has_xist) "XIST (log2)" else "sample rank by Y-gene mean (XIST unavailable)", y = "Y-gene mean (log2)", colour = "recorded sex", shape = "inferred sex") + theme_bw()
  save_gg(p_sex, sprintf("04_%s_sexcheck", gse), 7, 5.5)

  # ---------------- PCA + hierarchical clustering (§2.3) ----------------
  ntop <- min(STAGE1$pca_top, nrow(expr))
  top <- order(apply(expr, 1, var), decreasing = TRUE)[seq_len(ntop)]
  xt <- expr[top, , drop = FALSE]
  pc <- prcomp(t(xt), center = TRUE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  z1 <- as.numeric(scale(pc$x[, 1])); z2 <- as.numeric(scale(pc$x[, 2]))
  flag <- abs(z1) > STAGE1$pca_sd | abs(z2) > STAGE1$pca_sd
  cat(sprintf("  PCA on top %d genes: PC1 %.1f%%, PC2 %.1f%% variance; max |z| PC1 %.2f, PC2 %.2f\n",
              ntop, ve[1], ve[2], max(abs(z1)), max(abs(z2))))
  if (nrow(pheno_qc) <= 10) cat(sprintf("  note: with n=%d the largest attainable |z| is (n-1)/sqrt(n) = %.2f, so the 3 SD flag cannot fire\n",
                                        nrow(pheno_qc), (nrow(pheno_qc) - 1) / sqrt(nrow(pheno_qc))))
  cat(sprintf("  PCA-flagged (>3 SD, NOT excluded): %d %s\n", sum(flag),
              if (any(flag)) paste0("(", paste(sprintf("%s PC1z=%.2f PC2z=%.2f", pheno_qc$gsm[flag], z1[flag], z2[flag]), collapse = "; "), ")") else ""))
  pheno_qc$PC1_z <- z1; pheno_qc$PC2_z <- z2; pheno_qc$pca_outlier_flag <- flag

  cv <- colour_var(gse, pheno_qc)
  pdf_df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], col = as.character(cv$values), gsm = pheno_qc$gsm, flag = flag,
                       shp = if (gse == "GSE41762") pheno_qc$submitter_norm_batch else if (gse == "GSE76894") pheno_qc$sex else "all")
  shp_name <- if (gse == "GSE41762") "submitter_norm_batch" else if (gse == "GSE76894") "recorded sex" else NULL
  mu <- colMeans(pc$x[, 1:2]); sdv <- apply(pc$x[, 1:2], 2, sd)
  p_pca <- ggplot(pdf_df, aes(PC1, PC2, colour = col, shape = shp)) + geom_point(size = 2.5) +
    geom_vline(xintercept = mu[1] + c(-1, 1) * STAGE1$pca_sd * sdv[1], linetype = 2, colour = "grey50") +
    geom_hline(yintercept = mu[2] + c(-1, 1) * STAGE1$pca_sd * sdv[2], linetype = 2, colour = "grey50") +
    geom_text(data = pdf_df[pdf_df$flag, ], aes(label = gsm), vjust = -0.9, size = 3, show.legend = FALSE) +
    labs(title = sprintf("%s PCA (top %d genes by variance)%s", gse, ntop,
                         if (STAGE1$role[gse] == "validation") " - validation set, label-blind" else ""),
         subtitle = sprintf("dashed = mean +/- 3 SD; flagged (not excluded): %d", sum(flag)),
         x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), colour = cv$name, shape = shp_name) +
    theme_bw()
  if (is.null(shp_name)) p_pca <- p_pca + guides(shape = "none")
  save_gg(p_pca, sprintf("04_%s_pca", gse), 8, 6)

  hc <- hclust(as.dist(1 - cor(xt, method = "pearson")), method = "average")
  lab_col <- if (all(is.na(cv$values))) setNames(rep("black", ncol(xt)), colnames(xt)) else
    setNames(pal(cv$values)[as.character(cv$values)], colnames(xt))
  lab_col[is.na(lab_col)] <- "black"
  dend <- dendrapply(as.dendrogram(hc), function(n) {
    if (is.leaf(n)) { l <- attr(n, "label"); attr(n, "nodePar") <- list(pch = NA, lab.col = lab_col[l],
                                                                         lab.cex = if (ncol(xt) > 50) 0.45 else 0.7) }
    n })
  draw_hc <- function() {
    par(mar = c(7, 4, 3, 1))
    plot(dend, main = sprintf("%s hierarchical clustering (1 - Pearson, average; top %d genes)", gse, ntop),
         ylab = "1 - r")
    if (!all(is.na(cv$values))) legend("topright", fill = pal(cv$values), legend = names(pal(cv$values)), title = cv$name,
                                       cex = 0.7, bg = "white")
  }
  save_base(draw_hc, sprintf("04_%s_hclust", gse), w = max(8, ncol(xt) * 0.13), h = 6)

  # ---------------- per-sample QC table ----------------
  qs <- data.table(dataset = gse, role = STAGE1$role[gse], gsm = pheno$gsm, group = pheno$group,
                   median_NUSE = round(nuse_med[pheno$gsm], 4), RLE_median = round(rle_med[pheno$gsm], 4),
                   RLE_IQR = round(rle_iqr[pheno$gsm], 4), nuse_excluded = nuse_excl)
  m <- match(qs$gsm, pheno_qc$gsm)
  qs[, `:=`(PC1_z = round(pheno_qc$PC1_z[m], 3), PC2_z = round(pheno_qc$PC2_z[m], 3),
            outlier_flag = pheno_qc$pca_outlier_flag[m],
            recorded_sex = if (gse == "GSE27951") NA_character_ else pheno$sex,
            inferred_sex = pheno_qc$inferred_sex[m], sex_mismatch = pheno_qc$sex_mismatch[m])]
  summary_rows[[gse]] <- qs
  pheno_qc$median_NUSE <- nuse_med[pheno_qc$gsm]; pheno_qc$RLE_IQR <- rle_iqr[pheno_qc$gsm]

  # assertions per dataset
  stopifnot(identical(colnames(expr), pheno_qc$gsm), identical(ef$gene, rownames(expr)), !anyNA(expr),
            !anyNA(pheno_qc$inferred_sex), nrow(qs) == n_in, nrow(pheno_qc) == n_in - sum(nuse_excl),
            all(pheno_qc$inferred_sex %in% c("male", "female")))
  if (gse == "GSE27951") stopifnot(!anyNA(pheno_qc$sex))

  saveRDS(expr, here("data", "processed", sprintf("04_%s_expr.rds", gse)))
  saveRDS(pheno_qc, here("data", "processed", sprintf("04_%s_pheno.rds", gse)))
  saveRDS(ef, here("data", "processed", sprintf("04_%s_expressed.rds", gse)))
  qc_all[[gse]] <- list(nuse_median = nuse_med, rle_median = rle_med, rle_iqr = rle_iqr, nuse_excluded = pheno$gsm[nuse_excl],
                        pca_var_explained = ve[1:5], pca_scores = pc$x[, 1:5], top_genes = rownames(xt),
                        sex_genes = c(if (has_xist) "XIST", yg), kmeans_centers_scaled = km$centers, hclust = hc)
  rm(expr, xt, pc); gc(verbose = FALSE)
}

# ---- write NUSE exclusions (idempotent) --------------------------------------------------
excl_f <- here("logs", "excluded_samples.csv")
ex <- fread(excl_f, colClasses = "character")
if (length(excl_new)) {
  en <- rbindlist(excl_new)
  en <- en[!paste(gsm, gse) %in% paste(ex$gsm, ex$gse)]
  if (nrow(en)) { ex <- rbind(ex, en); fwrite(ex, excl_f) }
  cat("\nNUSE exclusions logged to logs/excluded_samples.csv\n")
} else cat("\nNo array exceeded median NUSE 1.05; nothing appended to logs/excluded_samples.csv\n")

qc <- rbindlist(summary_rows)
fwrite(qc, here("results", "stage1", "qc_summary.csv"))
saveRDS(qc_all, here("results", "stage1", "04_qc.rds"))

cat("\n==== Per-dataset QC overview ====\n")
print(qc[, .(n_in = .N, nuse_excluded = sum(nuse_excluded), pca_flagged = sum(outlier_flag, na.rm = TRUE),
             sex_recorded = sum(!is.na(recorded_sex)), sex_mismatch = sum(sex_mismatch, na.rm = TRUE),
             max_median_NUSE = if (all(is.na(median_NUSE))) NA_real_ else max(median_NUSE, na.rm = TRUE), max_RLE_IQR = max(RLE_IQR)), by = dataset])
cat("\nFlagged samples (PCA > 3 SD or sex mismatch) -- NOT excluded; project lead decides:\n")
print(qc[outlier_flag == TRUE | sex_mismatch == TRUE, .(dataset, gsm, PC1_z, PC2_z, outlier_flag, recorded_sex, inferred_sex, sex_mismatch)])

# ---- final assertions ----------------------------------------------------------------------
exp_n <- c(GSE41762 = 77, GSE27951 = 23, GSE76894 = 99, GSE71416 = 20, GSE20966 = 20, GSE16415 = 10)
stopifnot(all(qc[, .N, by = dataset][, N == exp_n[dataset]]),
          setequal(names(qc), c("dataset", "role", "gsm", "group", "median_NUSE", "RLE_median", "RLE_IQR", "nuse_excluded",
                                "PC1_z", "PC2_z", "outlier_flag", "recorded_sex", "inferred_sex", "sex_mismatch")),
          !anyDuplicated(qc$gsm), all(!is.na(qc[nuse_excluded == FALSE]$PC1_z)),
          all(!is.na(qc[dataset != "GSE16415"]$median_NUSE)))
for (gse in STAGE1$datasets) {
  for (k in c("pca", "hclust", "sexcheck", "rle_nuse")) for (ext in c("pdf", "png"))
    stopifnot(file.exists(file.path(fig_dir, sprintf("04_%s_%s.%s", gse, k, ext))))
  e <- readRDS(here("data", "processed", sprintf("04_%s_expr.rds", gse)))
  p <- readRDS(here("data", "processed", sprintf("04_%s_pheno.rds", gse)))
  stopifnot(identical(colnames(e), p$gsm), ncol(e) == exp_n[gse] - sum(qc[dataset == gse]$nuse_excluded))
}
cat("\nAll assertions passed.\n\n")
print(sessionInfo())

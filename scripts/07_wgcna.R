# 07_wgcna.R -- WGCNA co-expression network analysis, both arms (protocol v1.1 §4.1; pre-analysis plan A3, A4, B4).
#
# ROLE (fixed, §4.1 / plan B4): DESCRIPTIVE and CORROBORATIVE only. Nothing here feeds back into the ML panel or the nested CV.
#
# Per arm (islet GSE41762 n=77; adipose GSE27951 n=23; discovery cohorts, never merged):
#   input   : full normalised matrix, top 5000 genes by MAD (NOT senescence-filtered)
#   network : signed hybrid, bicor (maxPOutliers = 0.05), TOMType signed, blockwiseModules single block (maxBlockSize 6000 >= 5000),
#             minModuleSize 30, mergeCutHeight 0.25, deepSplit 2 (other blockwiseModules arguments left at package defaults,
#             incl. pamRespectsDendro = TRUE, reassignThreshold = 1e-6)
#   beta    : lowest power in 1..20 whose scale-free fit statistic SFT.R.sq >= 0.85.
#             Fit statistic = WGCNA::pickSoftThreshold "SFT.R.sq" = signed R^2 = -sign(slope) * R^2 of log10 p(k) ~ log10 k (10 connectivity bins),
#             i.e. a positive slope can never pass. If NO beta reaches 0.85: the curve is saved and that arm is STOPPED; no beta is picked.
#   traits  : module eigengene vs status (T2D=1), age, sex (male=1), BMI -- Pearson r, Student p. Kept: abs(r) > 0.3 and p < 0.05 vs status.
#             Flag "possibly confounded" if p < 0.05 vs age, sex or BMI. The unassigned set (label 0, grey) is never a module.
#   hubs    : kME (bicor of gene with own-module eigengene) > 0.8 AND abs(GS) > 0.2 (GS = Pearson r gene vs status), intersected with
#             senescence-union genes that are DEGs in the descriptive §3.3 table (results/stage4/07a_deg_<arm>.csv)
#   panels  : nested-CV final-panel genes (results/stage3/06_final_panel_islet.csv, 06b_final_panel_adipose_main.csv) located in modules -- corroboration only
#   adipose : bootstrap stability (100 stratified resamples, module detection re-run, modulePreservation Zsummary > 2; retained only if preserved in
#             > 70% of resamples; else EXPLORATORY). Runs ONLY if the adipose soft-threshold criterion is met.
# Validation cohorts (GSE76894, GSE20966, GSE71416, GSE16415) and GSE38642 are NEVER loaded.
# Outputs: results/stage4/07_*.csv / 07_*.rds, figures in results/stage4/figures/ (pdf + png). Log: logs/07_wgcna_<ts>.log

set.seed(42)
suppressPackageStartupMessages({ library(here); library(WGCNA); library(matrixStats); library(ggplot2) })
source(here("scripts", "lib", "stage4_wgcna_helpers.R"))
script_name <- "07_wgcna"; ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); T_START <- Sys.time()
dir.create(here("logs"), showWarnings = FALSE); out_dir <- here("results", "stage4"); fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt"); sink(con, split = TRUE); sink(con, type = "message")
finish <- function(status = 0L) { sink(type = "message"); sink(); close(con); quit(save = "no", status = status) }
CORES <- 8L; allowWGCNAThreads(CORES)
POWERS <- 1:20; R2_CUT <- 0.85; TOPN <- 5000L; TRAIT_R <- 0.3; ALPHA <- 0.05; KME_CUT <- 0.8; GS_CUT <- 0.2
N_BOOT <- 100L; N_PERM <- 200L; Z_CUT <- 2; FRAC_CUT <- 0.70
WP <- wgcna_params()

cat("Script:", script_name, " run:", ts, " cores:", CORES, " WGCNA", as.character(packageVersion("WGCNA")), "\n")
cat("ROLE: descriptive / corroborative only (protocol 4.1, plan B4). Nothing here feeds the ML panel or the nested CV.\n\nInputs:\n")
f_in <- c("data/processed/04_GSE41762_expr.rds", "data/processed/04_GSE41762_pheno.rds", "data/processed/04_GSE27951_expr.rds", "data/processed/04_GSE27951_pheno.rds",
          "genesets/senescence_union.rds", "results/stage4/07a_deg_islet.csv", "results/stage4/07a_deg_adipose.csv",
          "results/stage3/06_final_panel_islet.csv", "results/stage3/06b_final_panel_adipose_main.csv", "scripts/lib/stage4_wgcna_helpers.R")
for (f in f_in) cat("  ", f, " md5=", tools::md5sum(here(f)), "\n", sep = "")
cat("  (validation cohorts GSE76894, GSE20966, GSE71416, GSE16415 and GSE38642 expression are NOT loaded)\n")
cat(sprintf("\nParameters: top %d genes by MAD; powers %d..%d; R2 cut %.2f (SFT.R.sq, signed R^2); network %s, TOM %s, corType %s maxPOutliers %.2f; minModuleSize %d; mergeCutHeight %.2f; deepSplit %d; maxBlockSize %d\n",
            TOPN, min(POWERS), max(POWERS), R2_CUT, WP$networkType, WP$TOMType, WP$corType, WP$maxPOutliers, WP$minModuleSize, WP$mergeCutHeight, WP$deepSplit, WP$maxBlockSize))
cat(sprintf("Module-trait: Pearson; keep abs(r) > %.1f & p < %.2f vs status; confounding flag p < %.2f vs age/sex/BMI. Hubs: kME > %.1f & abs(GS) > %.1f & senescence-union DEG (07a).\n",
            TRAIT_R, ALPHA, ALPHA, KME_CUT, GS_CUT))
cat(sprintf("Adipose stability (only if beta criterion met): %d stratified bootstrap resamples; nPermutations = %d (script parameter, not a deviation); Zsummary > %g = preserved; retain if preserved in > %.0f%%.\n",
            N_BOOT, N_PERM, Z_CUT, 100 * FRAC_CUT))

geneset <- readRDS(here("genesets/senescence_union.rds")); stopifnot(is.character(geneset), length(geneset) == 1151L)

# analysis phenotype: identical logic to prep_pheno() in scripts/06_nested_cv.R (sex recorded, expression-inferred where mismatch / absent)
prep_pheno <- function(pheno, arm) {
  sex_rec <- tolower(trimws(as.character(pheno$sex))); sex_inf <- tolower(trimws(as.character(pheno$inferred_sex)))
  sex <- ifelse(pheno$sex_mismatch %in% TRUE | is.na(sex_rec), sex_inf, sex_rec)
  out <- data.frame(status = factor(pheno$group, levels = c("ND", "T2D")), age = as.numeric(pheno$age), sex = factor(sex),
                    bmi = as.numeric(pheno$bmi), row.names = pheno$gsm, stringsAsFactors = FALSE)
  stopifnot(!anyNA(out$status), !anyNA(out$age), !anyNA(out$sex), nlevels(out$sex) == 2L)
  out
}

mod_colour <- function(lab) ifelse(lab == 0, "grey", WGCNA::labels2colors(lab))   # numeric labels 1.. map to colours in size order

run_arm <- function(gse, arm, n_exp, n_t2d, n_nd, panel_file) {
  cat(sprintf("\n\n################ %s arm (%s) ################\n", toupper(arm), gse))
  e <- readRDS(here("data/processed", sprintf("04_%s_expr.rds", gse))); p <- readRDS(here("data/processed", sprintf("04_%s_pheno.rds", gse)))
  stopifnot(identical(colnames(e), p$gsm), ncol(e) == n_exp, !anyNA(e))
  ph <- prep_pheno(p, arm); stopifnot(sum(ph$status == "T2D") == n_t2d, sum(ph$status == "ND") == n_nd)
  deg <- read.csv(here("results/stage4", sprintf("07a_deg_%s.csv", arm)), stringsAsFactors = FALSE)
  # ---- input: top 5000 by MAD of the FULL matrix (not senescence-filtered) ----
  mads <- rowMads(e); names(mads) <- rownames(e)
  ord <- order(mads, decreasing = TRUE); mad_rank <- setNames(seq_along(ord), rownames(e)[ord])
  top <- rownames(e)[ord[seq_len(TOPN)]]
  datExpr <- t(e[top, , drop = FALSE]); stopifnot(nrow(datExpr) == n_exp, ncol(datExpr) == TOPN, !anyDuplicated(colnames(datExpr)))
  gsg <- goodSamplesGenes(datExpr, verbose = 0); stopifnot(gsg$allOK)
  cat(sprintf("n=%d (%d ND / %d T2D); matrix genes=%d; WGCNA input = top %d by MAD (MAD range %.3f..%.3f); senescence-union genes in input: %d\n",
              n_exp, n_nd, n_t2d, nrow(e), TOPN, min(mads[top]), max(mads[top]), sum(top %in% geneset)))
  # ---- soft threshold ----
  t0 <- Sys.time()
  sft <- pickSoftThreshold(datExpr, powerVector = POWERS, networkType = WP$networkType, corFnc = bicor,
                           corOptions = list(maxPOutliers = WP$maxPOutliers), verbose = 0)
  fi <- sft$fitIndices; stopifnot(identical(as.integer(fi$Power), as.integer(POWERS)))
  cat(sprintf("pickSoftThreshold: %.1f s\nScale-free fit (SFT.R.sq = signed R^2, 10 bins; mean.k = mean connectivity):\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  print(data.frame(beta = fi$Power, SFT.R.sq = round(fi$SFT.R.sq, 3), slope = round(fi$slope, 3), mean.k = round(fi$mean.k., 1), median.k = round(fi$median.k., 1), max.k = round(fi$max.k., 1)), row.names = FALSE)
  write.csv(data.frame(arm = arm, beta = fi$Power, SFT.R.sq_signed = fi$SFT.R.sq, slope = fi$slope, truncated.R.sq = fi$truncated.R.sq, mean.k = fi$mean.k.,
                       median.k = fi$median.k., max.k = fi$max.k., meets_cut = fi$SFT.R.sq >= R2_CUT), file.path(out_dir, sprintf("07_softthreshold_%s.csv", arm)), row.names = FALSE)
  ok_beta <- which(fi$SFT.R.sq >= R2_CUT)
  beta <- if (length(ok_beta)) fi$Power[min(ok_beta)] else NA_integer_
  # scale-free figure (always saved, also when the arm stops)
  scale_free_fig <- function() {
    par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    plot(fi$Power, fi$SFT.R.sq, type = "n", xlab = "Soft-thresholding power (beta)", ylab = "Scale-free fit, signed R^2 (SFT.R.sq)",
         main = sprintf("%s (%s): scale-free topology fit", toupper(arm), gse), ylim = c(min(0, min(fi$SFT.R.sq)), 1))
    text(fi$Power, fi$SFT.R.sq, labels = fi$Power, col = "red", cex = 0.85); abline(h = R2_CUT, col = "blue", lty = 2)
    legend("bottomright", legend = sprintf("R^2 = %.2f cut", R2_CUT), lty = 2, col = "blue", bty = "n", cex = 0.85)
    if (is.na(beta)) mtext(sprintf("NO beta reaches %.2f (max %.3f at beta=%d): arm STOPPED", R2_CUT, max(fi$SFT.R.sq), fi$Power[which.max(fi$SFT.R.sq)]), side = 3, line = -1.5, col = "red", cex = 0.8)
    else abline(v = beta, col = "darkgreen", lty = 3)
    plot(fi$Power, fi$mean.k., type = "n", xlab = "Soft-thresholding power (beta)", ylab = "Mean connectivity", main = "Mean connectivity", log = "y")
    text(fi$Power, fi$mean.k., labels = fi$Power, col = "red", cex = 0.85)
  }
  save_fig(fig_dir, sprintf("07_scalefree_%s", arm), scale_free_fig, width = 11, height = 5)
  if (is.na(beta)) {
    msg <- sprintf("STOP: no beta in %d..%d reaches scale-free R^2 >= %.2f (signed R^2); best = %.3f at beta = %d. Curve saved; NO beta chosen; no modules, no trait analysis, no stability run.",
                   min(POWERS), max(POWERS), R2_CUT, max(fi$SFT.R.sq), fi$Power[which.max(fi$SFT.R.sq)])
    cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", toupper(arm), " ARM ", msg, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", sep = "")
    return(list(arm = arm, stopped = TRUE, verdict = msg, best_R2 = max(fi$SFT.R.sq), best_beta = fi$Power[which.max(fi$SFT.R.sq)], n = n_exp, fit = fi))
  }
  cat(sprintf("\nChosen beta = %d (lowest power with SFT.R.sq >= %.2f; SFT.R.sq = %.3f, slope = %.2f, mean k = %.1f)\n", beta, R2_CUT, fi$SFT.R.sq[beta], fi$slope[beta], fi$mean.k.[beta]))
  # ---- modules ----
  t0 <- Sys.time(); set.seed(42)
  net <- run_blockwise(datExpr, beta, nThreads = CORES)
  stopifnot(length(net$dendrograms) == 1L, length(net$blockGenes[[1]]) == TOPN)
  lab <- net$colors; names(lab) <- colnames(datExpr); cols <- mod_colour(lab)
  mods <- sort(setdiff(unique(lab), 0L)); sizes <- table(lab)
  cat(sprintf("blockwiseModules: %.1f s; modules (excl. unassigned) = %d; unassigned genes = %d\n", as.numeric(difftime(Sys.time(), t0, units = "secs")), length(mods), sum(lab == 0)))
  print(data.frame(label = mods, colour = mod_colour(mods), n_genes = as.integer(sizes[as.character(mods)])), row.names = FALSE)
  MEs <- net$MEs; stopifnot(identical(rownames(MEs), rownames(datExpr)), all(paste0("ME", mods) %in% colnames(MEs)))
  # dendrogram
  dendro_fig <- function() plotDendroAndColors(net$dendrograms[[1]], cols[net$blockGenes[[1]]], "Module colours", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                                               main = sprintf("%s (%s): gene dendrogram, beta=%d, %d modules", toupper(arm), gse, beta, length(mods)))
  save_fig(fig_dir, sprintf("07_dendrogram_%s", arm), dendro_fig, width = 11, height = 6)
  # ---- module-trait ----
  traits <- cbind(status_T2D = as.numeric(ph$status == "T2D"), age = ph$age, sex_male = as.numeric(ph$sex == "male"), BMI = ph$bmi)
  rownames(traits) <- rownames(ph); stopifnot(identical(rownames(traits), rownames(datExpr)), !anyNA(traits))
  ME_m <- MEs[, paste0("ME", mods), drop = FALSE]
  mt_r <- cor(ME_m, traits, use = "pairwise.complete.obs"); mt_p <- corPvalueStudent(mt_r, n_exp)
  mt <- data.frame(label = mods, colour = mod_colour(mods), n_genes = as.integer(sizes[as.character(mods)]),
                   r_status = mt_r[, "status_T2D"], p_status = mt_p[, "status_T2D"], r_age = mt_r[, "age"], p_age = mt_p[, "age"],
                   r_sex = mt_r[, "sex_male"], p_sex = mt_p[, "sex_male"], r_BMI = mt_r[, "BMI"], p_BMI = mt_p[, "BMI"], stringsAsFactors = FALSE)
  mt$kept <- abs(mt$r_status) > TRAIT_R & mt$p_status < ALPHA
  mt$confounded_with <- apply(mt, 1, function(z) paste(c("age", "sex", "BMI")[c(as.numeric(z["p_age"]) < ALPHA, as.numeric(z["p_sex"]) < ALPHA, as.numeric(z["p_BMI"]) < ALPHA)], collapse = "+"))
  mt$possibly_confounded <- mt$confounded_with != ""
  cat("\nModule-trait correlations (Pearson; eigengene vs trait):\n")
  print(data.frame(mt[, c("label", "colour", "n_genes")], r_status = round(mt$r_status, 3), p_status = signif(mt$p_status, 3), r_age = round(mt$r_age, 2), p_age = signif(mt$p_age, 2),
                   r_sex = round(mt$r_sex, 2), p_sex = signif(mt$p_sex, 2), r_BMI = round(mt$r_BMI, 2), p_BMI = signif(mt$p_BMI, 2), kept = mt$kept, confounded = mt$confounded_with), row.names = FALSE)
  cat(sprintf("Modules kept (abs(r)>%.1f & p<%.2f vs status): %d of %d%s\n", TRAIT_R, ALPHA, sum(mt$kept), nrow(mt), if (any(mt$kept)) paste0(" -> ", paste(sprintf("M%d(%s)", mt$label[mt$kept], mt$colour[mt$kept]), collapse = ", ")) else ""))
  write.csv(mt, file.path(out_dir, sprintf("07_module_trait_%s.csv", arm)), row.names = FALSE)
  # heatmap
  hm <- do.call(rbind, lapply(colnames(traits), function(tn) data.frame(module = sprintf("M%d %s (n=%d)", mt$label, mt$colour, mt$n_genes), trait = c(status_T2D = "Status (T2D vs ND)", age = "Age", sex_male = "Sex (male)", BMI = "BMI")[[tn]],
                                                                       r = mt_r[, tn], p = mt_p[, tn], stringsAsFactors = FALSE)))
  hm$label <- sprintf("%.2f\n(p=%s)", hm$r, ifelse(hm$p < 0.001, formatC(hm$p, format = "e", digits = 1), formatC(hm$p, format = "f", digits = 3)))
  hm$module <- factor(hm$module, levels = rev(sprintf("M%d %s (n=%d)", mt$label, mt$colour, mt$n_genes)))
  hm$trait <- factor(hm$trait, levels = c("Status (T2D vs ND)", "Age", "Sex (male)", "BMI")); hm$sig <- hm$p < ALPHA
  g <- ggplot(hm, aes(trait, module, fill = r)) + geom_tile(colour = "white") + geom_text(aes(label = label, fontface = ifelse(sig, "bold", "plain")), size = 2.6) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = "Pearson r") +
    labs(title = sprintf("%s (%s, n=%d): module-trait relationships", toupper(arm), gse, n_exp), subtitle = sprintf("beta=%d; kept = |r|>0.3 & p<0.05 vs status; bold = p<0.05. Descriptive only.", beta), x = NULL, y = NULL) +
    theme_minimal(base_size = 10) + theme(panel.grid = element_blank())
  hgt <- max(4, 0.55 * length(mods) + 1.8)
  ggsave(file.path(fig_dir, sprintf("07_module_trait_%s.pdf", arm)), g, width = 7.5, height = hgt); ggsave(file.path(fig_dir, sprintf("07_module_trait_%s.png", arm)), g, width = 7.5, height = hgt, dpi = 150)
  # ---- gene-level: kME, GS, flags ----
  kme_all <- bicor(datExpr, MEs, maxPOutliers = WP$maxPOutliers, use = "pairwise.complete.obs")
  own <- paste0("ME", lab); kME <- ifelse(lab == 0, NA_real_, kme_all[cbind(seq_along(lab), match(own, colnames(kme_all)))])
  gs_r <- as.numeric(cor(datExpr, traits[, "status_T2D"], use = "pairwise.complete.obs")); gs_p <- as.numeric(corPvalueStudent(gs_r, n_exp))
  dg <- deg[match(colnames(datExpr), deg$gene), ]   # NA where the gene is not in the expressed-gene DEG table
  gt <- data.frame(gene = colnames(datExpr), module_label = as.integer(lab), module_colour = cols, kME = kME, GS_status = gs_r, GS_p = gs_p, MAD_rank = as.integer(mad_rank[colnames(datExpr)]),
                   in_senescence_union = colnames(datExpr) %in% geneset, in_07a_table = !is.na(dg$gene), DEG_log2FC = dg$log2FC, DEG_adjP = dg$adj.P.Val,
                   DEG = dg$DEG %in% TRUE, senescence_DEG = dg$senescence_DEG %in% TRUE, stringsAsFactors = FALSE)
  gt$module_kept <- gt$module_label %in% mt$label[mt$kept]
  gt$module_possibly_confounded <- gt$module_label %in% mt$label[mt$possibly_confounded]
  gt$kME_GS_pass <- !is.na(gt$kME) & gt$kME > KME_CUT & abs(gt$GS_status) > GS_CUT
  gt$hub <- gt$kME_GS_pass & gt$module_kept & gt$senescence_DEG     # protocol §4.1 hub rule within trait-associated modules
  write.csv(gt, file.path(out_dir, sprintf("07_gene_table_%s.csv", arm)), row.names = FALSE)
  cat(sprintf("\nHubs: genes with kME>%.1f & |GS|>%.1f: all assigned modules = %d; within kept modules = %d; of those, senescence-union DEGs (07a) = %d\n",
              KME_CUT, GS_CUT, sum(gt$kME_GS_pass), sum(gt$kME_GS_pass & gt$module_kept), sum(gt$hub)))
  cat(sprintf("senescence-union DEGs (07a) present in the 5000-gene WGCNA input: %d of %d total senescence DEGs in 07a\n", sum(gt$senescence_DEG), sum(deg$senescence_DEG)))
  hubs <- gt[gt$hub, c("gene", "module_label", "module_colour", "kME", "GS_status", "GS_p", "DEG_log2FC", "DEG_adjP", "module_possibly_confounded")]
  if (nrow(hubs)) { cat("HUB GENES (kept module & kME>0.8 & |GS|>0.2 & senescence-union DEG):\n"); print(hubs, row.names = FALSE) } else cat("HUB GENES: none\n")
  write.csv(hubs, file.path(out_dir, sprintf("07_hubs_%s.csv", arm)), row.names = FALSE)
  # senescence DEGs that fail a hub criterion (transparency)
  sd_tab <- gt[gt$senescence_DEG, c("gene", "module_label", "module_colour", "kME", "GS_status", "module_kept", "kME_GS_pass", "hub")]
  if (nrow(sd_tab)) { cat("\nAll senescence-union DEGs located in the WGCNA input (hub criteria breakdown):\n"); print(sd_tab, row.names = FALSE) }
  # ---- nested-CV panel genes ----
  panel <- read.csv(here(panel_file), stringsAsFactors = FALSE)
  pn <- do.call(rbind, lapply(panel$gene, function(gn) {
    i <- match(gn, gt$gene)
    data.frame(arm = arm, panel_gene = gn, in_matrix = gn %in% rownames(e), MAD_rank_full_matrix = as.integer(mad_rank[gn]), in_top5000 = !is.na(i),
               module_label = if (is.na(i)) NA_integer_ else gt$module_label[i], module_colour = if (is.na(i)) NA_character_ else gt$module_colour[i],
               kME = if (is.na(i)) NA_real_ else gt$kME[i], GS_status = if (is.na(i)) NA_real_ else gt$GS_status[i], GS_p = if (is.na(i)) NA_real_ else gt$GS_p[i],
               module_kept = if (is.na(i)) NA else gt$module_kept[i], module_possibly_confounded = if (is.na(i)) NA else gt$module_possibly_confounded[i],
               kME_GS_pass = if (is.na(i)) NA else gt$kME_GS_pass[i], senescence_DEG_07a = {j <- match(gn, deg$gene); !is.na(j) && deg$senescence_DEG[j]},
               hub_by_protocol_rule = if (is.na(i)) NA else gt$hub[i], stringsAsFactors = FALSE)
  }))
  cat("\nNested-CV final panel genes vs WGCNA modules (corroboration only):\n"); print(pn, row.names = FALSE)
  write.csv(pn, file.path(out_dir, sprintf("07_panel_genes_modules_%s.csv", arm)), row.names = FALSE)
  list(arm = arm, stopped = FALSE, beta = beta, R2 = fi$SFT.R.sq[beta], slope = fi$slope[beta], mean_k = fi$mean.k.[beta], n = n_exp, n_modules = length(mods), n_unassigned = sum(lab == 0),
       mt = mt, gene_table = gt, hubs = hubs, panel = pn, labels = lab, datExpr = datExpr, status = ph$status, fit = fi, kept = mt$label[mt$kept])
}

A_isl <- run_arm("GSE41762", "islet", 77L, 20L, 57L, "results/stage3/06_final_panel_islet.csv")
A_adi <- run_arm("GSE27951", "adipose", 23L, 12L, 11L, "results/stage3/06b_final_panel_adipose_main.csv")

# ---- adipose stability (only when the adipose beta criterion is met) ------------------------------------------------
stab <- NULL
if (!A_adi$stopped) {
  cat(sprintf("\n=== ADIPOSE STABILITY: %d stratified bootstrap resamples, nPermutations=%d, cores=%d ===\n", N_BOOT, N_PERM, CORES))
  t0 <- Sys.time()
  S <- run_stability(A_adi$datExpr, A_adi$status, A_adi$labels, A_adi$beta, B = N_BOOT, nPerm = N_PERM, zcut = Z_CUT, cores = CORES, seed = 42L)
  cat(sprintf("stability runtime: %.1f min; resamples completed: %d / %d\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(S$ok), N_BOOT))
  if (!all(S$ok)) { cat("FAILED resamples:\n"); print(unique(vapply(S$results[!S$ok], function(r) if (is.list(r)) r$msg else as.character(r), ""))) }
  Zl <- lapply(S$results[S$ok], function(r) transform(r$Z, b = r$b)); Zall <- do.call(rbind, Zl); Zall <- Zall[!Zall$module %in% c("gold", "0"), ]
  mlab <- A_adi$mt$label
  frac <- vapply(mlab, function(m) { z <- Zall$Zsummary[Zall$module == as.character(m)]; sum(z > Z_CUT, na.rm = TRUE) / sum(S$ok) }, 0)
  stab <- data.frame(label = mlab, colour = A_adi$mt$colour, n_genes = A_adi$mt$n_genes, resamples = sum(S$ok),
                     n_preserved = round(frac * sum(S$ok)), fraction_preserved = frac, median_Zsummary = vapply(mlab, function(m) median(Zall$Zsummary[Zall$module == as.character(m)], na.rm = TRUE), 0),
                     retained = frac > FRAC_CUT, trait_kept = A_adi$mt$kept)
  cat("\nAdipose preservation (reference = full-cohort modules):\n"); print(transform(stab, fraction_preserved = round(fraction_preserved, 2), median_Zsummary = round(median_Zsummary, 1)), row.names = FALSE)
  write.csv(stab, file.path(out_dir, "07_adipose_preservation.csv"), row.names = FALSE); write.csv(Zall, file.path(out_dir, "07_adipose_preservation_Z_by_resample.csv"), row.names = FALSE)
  verdict_adi <- if (any(stab$retained)) sprintf("Modules preserved in > 70%% of resamples: %s", paste(sprintf("M%d(%s)", stab$label[stab$retained], stab$colour[stab$retained]), collapse = ", ")) else
    "NO module is preserved in > 70% of resamples: the adipose WGCNA is EXPLORATORY."
  cat("\nADIPOSE VERDICT:", verdict_adi, "\n")
  gb <- ggplot(stab, aes(factor(sprintf("M%d %s", label, colour), levels = sprintf("M%d %s", label, colour)), fraction_preserved, fill = retained)) + geom_col() +
    geom_hline(yintercept = FRAC_CUT, linetype = 2) + scale_fill_manual(values = c(`TRUE` = "#2a9d8f", `FALSE` = "#bbbbbb"), name = "retained (>70%)") + scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Adipose (GSE27951) module stability", subtitle = sprintf("Fraction of %d bootstrap resamples with Zsummary > 2; dashed line = 70%% retention threshold", sum(S$ok)), x = NULL, y = "Fraction preserved") + theme_minimal(base_size = 10)
  ggsave(file.path(fig_dir, "07_adipose_preservation.pdf"), gb, width = 7, height = 4.5); ggsave(file.path(fig_dir, "07_adipose_preservation.png"), gb, width = 7, height = 4.5, dpi = 150)
} else {
  verdict_adi <- paste0("ADIPOSE ARM STOPPED at the soft-threshold step (", A_adi$verdict, ") Stability analysis (100 bootstraps) and adipose module-trait analysis NOT run; adipose WGCNA therefore yields no modules and is not reported (protocol §4.1: never pick a beta silently).")
  cat("\n", verdict_adi, "\n", sep = "")
}

# ---- verdict table + saved objects -------------------------------------------------------------------------------
verd <- data.frame(arm = c("islet", "adipose"), stopped = c(A_isl$stopped, A_adi$stopped),
                   beta = c(A_isl$beta %||% NA, if (A_adi$stopped) NA else A_adi$beta), SFT_R2 = c(A_isl$R2 %||% NA, if (A_adi$stopped) NA else A_adi$R2),
                   best_R2_if_stopped = c(if (A_isl$stopped) A_isl$best_R2 else NA, if (A_adi$stopped) A_adi$best_R2 else NA),
                   best_beta_if_stopped = c(if (A_isl$stopped) A_isl$best_beta else NA, if (A_adi$stopped) A_adi$best_beta else NA),
                   n_modules = c(A_isl$n_modules %||% NA, A_adi$n_modules %||% NA),
                   n_kept = c(if (A_isl$stopped) NA else length(A_isl$kept), if (A_adi$stopped) NA else length(A_adi$kept)),
                   n_hubs = c(if (A_isl$stopped) NA else nrow(A_isl$hubs), if (A_adi$stopped) NA else nrow(A_adi$hubs)),
                   verdict = c(if (A_isl$stopped) A_isl$verdict else "network built", verdict_adi), stringsAsFactors = FALSE)
write.csv(verd, file.path(out_dir, "07_verdicts.csv"), row.names = FALSE)
strip <- function(A) if (A$stopped) A else A[setdiff(names(A), "datExpr")]
saveRDS(list(islet = strip(A_isl), adipose = strip(A_adi), stability = stab, verdict = verd,
             params = list(powers = POWERS, R2_cut = R2_CUT, topn = TOPN, kme = KME_CUT, gs = GS_CUT, trait_r = TRAIT_R, alpha = ALPHA, n_boot = N_BOOT, n_perm = N_PERM, z_cut = Z_CUT, frac_cut = FRAC_CUT, wgcna = WP),
             note = "DESCRIPTIVE / CORROBORATIVE ONLY (protocol 4.1); never fed into selection"), file.path(out_dir, "07_wgcna_results.rds"))

# ---- final assertions --------------------------------------------------------------------------------------------
if (!A_isl$stopped) {
  stopifnot(nrow(A_isl$gene_table) == TOPN, !anyNA(A_isl$gene_table[, c("gene", "module_label", "GS_status")]), nrow(A_isl$mt) == A_isl$n_modules, A_isl$beta %in% POWERS,
            all(A_isl$fit$SFT.R.sq[seq_len(A_isl$beta - 1L)] < R2_CUT), A_isl$R2 >= R2_CUT, all(A_isl$mt$n_genes >= WP$minModuleSize),
            sum(A_isl$mt$n_genes) + A_isl$n_unassigned == TOPN, !anyNA(A_isl$mt[, c("r_status", "p_status")]),
            all(A_isl$gene_table$kME[A_isl$gene_table$module_label != 0] <= 1 + 1e-9, na.rm = TRUE), all(!A_isl$gene_table$hub | A_isl$gene_table$senescence_DEG),
            nrow(A_isl$panel) == 1L, file.exists(file.path(fig_dir, "07_module_trait_islet.pdf")), file.exists(file.path(fig_dir, "07_module_trait_islet.png")),
            file.exists(file.path(fig_dir, "07_dendrogram_islet.pdf")), file.exists(file.path(fig_dir, "07_scalefree_islet.png")))
}
stopifnot(file.exists(file.path(fig_dir, "07_scalefree_adipose.pdf")), file.exists(file.path(fig_dir, "07_scalefree_adipose.png")), nrow(verd) == 2L,
          A_adi$stopped || (!is.null(stab) && nrow(stab) == A_adi$n_modules))
if (A_adi$stopped) stopifnot(A_adi$best_R2 < R2_CUT, is.null(stab))
cat("\nAll end-of-script assertions passed.\n")
cat(sprintf("Runtime: %.2f min\n", as.numeric(difftime(Sys.time(), T_START, units = "mins"))))
cat("\nsessionInfo():\n"); print(sessionInfo())
finish(0L)

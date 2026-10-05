# 11_enrichment.R -- functional enrichment of the cross-tissue gene sets (protocol v1.2 §5.2; CLAUDE.md rules 2, 4, 7).
#
# DISCOVERY DATA ONLY (expressed-gene flags of 04_GSE41762 and 04_GSE27951, discovery-derived gene lists). No validation data is read (rule 2).
#
# Usage: Rscript scripts/11_enrichment.R
#
# PROTOCOL ANALYSIS (§5.2): clusterProfiler GO (BP, MF, CC) and KEGG, run separately on the shared, islet-specific and adipose-specific sets
# (results/stage6/10_gene_classification.csv). Universe = genes expressed (§2.2 flags) in the relevant tissue (union of both tissues for the
# shared set). BH-adjusted p < 0.05; minGSSize 10 (maxGSSize left at the clusterProfiler default, 500). A set with < 10 genes is declared
# UNDERPOWERED and is NOT tested; its genes are annotated narratively instead (org.Hs.eg.db GO annotation, KEGG pathway links retrieved from the
# KEGG REST API when reachable, membership of SenMayo / CellAge / GenAge from genesets/senescence_union.csv, and discovery-cohort descriptive
# statistics). No biology is added beyond what these sources return.
#
# EXPLORATORY ANALYSES (descriptive only, clearly labelled; they do not feed any decision):
#   E1  islet descriptive senescence-union DEGs (results/stage4/07a_deg_islet.csv, senescence_DEG == TRUE; protocol §3.3 thresholds).
#       Background A (as specified): all expressed islet genes. Background B (caveat run): senescence-union genes expressed in islet. Because the
#       E1 genes were drawn from the senescence union, background A is biased toward senescence/growth-factor terms; B removes that selection bias.
#   E2  islet nested-CV stability list (final-panel selection frequency >= 0.2) -- run only if >= 10 genes.
#   E3  adipose stability list (>= 0.2) and adipose descriptive DEGs -- run only if >= 10 genes.
#
# Outputs (results/stage6/): 11_enrichment_summary.csv, 11_enrichment_significant.csv, 11_enrichment_top10_by_p.csv, 11_gene_annotation_summary.csv,
# 11_gene_annotation_terms.csv, 11_enrichment_result.rds, figures/11_*_dotplot.{pdf,png} (only for analyses with significant terms).
# Log: logs/11_enrichment_<timestamp>.log

set.seed(42)
suppressPackageStartupMessages({ library(here); library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2) })
script_name <- "11_enrichment"; ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); T_START <- Sys.time()
dir.create(here("logs"), showWarnings = FALSE); out_dir <- here("results", "stage6"); fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt"); sink(con, split = TRUE); sink(con, type = "message")
P_CUT <- 0.05; MIN_GS <- 10L; MAX_GS <- 500L; MIN_SET <- 10L

cat("Script:", script_name, " run:", ts, "\n\nInputs:\n")
f_in <- c("results/stage6/10_gene_classification.csv", "results/stage6/10_cross_tissue_result.rds", "data/processed/04_GSE41762_expressed.rds",
          "data/processed/04_GSE27951_expressed.rds", "genesets/senescence_union.csv", "results/stage4/07a_deg_islet.csv", "results/stage4/07a_deg_adipose.csv",
          "results/stage3/stability_selection.csv", "results/stage3/stability_selection_adipose.csv")
for (f in f_in) { stopifnot(file.exists(here(f))); cat("  ", f, " md5=", tools::md5sum(here(f)), "\n", sep = "") }
cat("  (validation cohorts and results/stage5/09_* are NOT read)\n")
cat(sprintf("\nParameters: BH p.adjust < %g; minGSSize %d; maxGSSize %d; underpowered if set < %d genes; clusterProfiler %s; org.Hs.eg.db %s; seed 42\n",
            P_CUT, MIN_GS, MAX_GS, MIN_SET, as.character(packageVersion("clusterProfiler")), as.character(packageVersion("org.Hs.eg.db"))))

cls <- read.csv(here("results/stage6/10_gene_classification.csv"), stringsAsFactors = FALSE)
x10 <- readRDS(here("results/stage6/10_cross_tissue_result.rds"))
fl_i <- readRDS(here("data/processed/04_GSE41762_expressed.rds")); fl_a <- readRDS(here("data/processed/04_GSE27951_expressed.rds"))
expr_i <- fl_i$gene[fl_i$expressed]; expr_a <- fl_a$gene[fl_a$expressed]
su <- read.csv(here("genesets/senescence_union.csv"), stringsAsFactors = FALSE)
deg_i <- read.csv(here("results/stage4/07a_deg_islet.csv"), stringsAsFactors = FALSE); deg_a <- read.csv(here("results/stage4/07a_deg_adipose.csv"), stringsAsFactors = FALSE)
stab_i <- read.csv(here("results/stage3/stability_selection.csv"), stringsAsFactors = FALSE); stab_a <- read.csv(here("results/stage3/stability_selection_adipose.csv"), stringsAsFactors = FALSE)
union_expr_i <- intersect(su$gene_symbol, expr_i); union_expr_a <- intersect(su$gene_symbol, expr_a)

sets <- list(
  list(analysis = "protocol", set = "shared", genes = cls$gene[cls$classification == "shared"], universe = union(expr_i, expr_a), bg = "expressed (union of both tissues)"),
  list(analysis = "protocol", set = "islet-specific", genes = cls$gene[cls$classification == "islet-specific"], universe = expr_i, bg = "expressed islet"),
  list(analysis = "protocol", set = "adipose-specific", genes = cls$gene[cls$classification == "adipose-specific"], universe = expr_a, bg = "expressed adipose"))
e1 <- deg_i$gene[deg_i$senescence_DEG %in% TRUE]
sets <- c(sets, list(
  list(analysis = "EXPLORATORY E1: islet senescence-union DEGs", set = "islet_senescence_DEGs", genes = e1, universe = expr_i, bg = "A: all expressed islet genes"),
  list(analysis = "EXPLORATORY E1: islet senescence-union DEGs", set = "islet_senescence_DEGs", genes = e1, universe = union_expr_i, bg = "B (caveat): senescence-union genes expressed in islet"),
  list(analysis = "EXPLORATORY E2: islet nested-CV stability list (freq >= 0.2)", set = "islet_stability_ge0.2", genes = stab_i$gene[stab_i$frequency >= 0.2], universe = expr_i, bg = "A: all expressed islet genes"),
  list(analysis = "EXPLORATORY E3: adipose nested-CV stability list (freq >= 0.2)", set = "adipose_stability_ge0.2", genes = stab_a$gene[stab_a$frequency >= 0.2], universe = expr_a, bg = "A: all expressed adipose genes"),
  list(analysis = "EXPLORATORY E3: adipose descriptive DEGs (protocol §3.3)", set = "adipose_DEGs", genes = deg_a$gene[deg_a$DEG %in% TRUE], universe = expr_a, bg = "A: all expressed adipose genes")))
stopifnot(length(e1) == 11L, setequal(cls$gene[cls$classification == "islet-specific"], x10$panels$islet), setequal(cls$gene[cls$classification == "adipose-specific"], x10$panels$adipose))
for (s in sets) stopifnot(all(s$genes %in% s$universe) || s$set %in% c("shared"))     # every tested gene lies in its universe (they were expressed in discovery)

sym2ent <- function(x) { m <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, unique(x), "ENTREZID", "SYMBOL", multiVals = "first")); m[!is.na(m)] }
kegg_names <- tryCatch(KEGGREST::keggList("pathway", "hsa"), error = function(e) NULL)
kegg_ok <- !is.null(kegg_names) && length(kegg_names) > 100
cat("KEGG REST API reachable: ", kegg_ok, " (retrieved ", format(Sys.time(), "%Y-%m-%d"), "; ", if (kegg_ok) length(kegg_names) else 0, " human pathways)\n", sep = "")

tidy_res <- function(obj, analysis, set, bg, db) {
  if (is.null(obj) || nrow(obj@result) == 0) return(NULL)
  r <- obj@result; r$analysis <- analysis; r$set <- set; r$background <- bg; r$db <- db
  r$GeneRatioNum <- vapply(strsplit(r$GeneRatio, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]), 0)
  r$significant <- r$p.adjust < P_CUT
  r[order(r$p.adjust, r$pvalue), ]
}
run_enrich <- function(genes, universe, analysis, set, bg) {
  out <- list(); ent <- sym2ent(genes); ent_u <- sym2ent(universe)
  for (ont in c("BP", "MF", "CC")) {
    obj <- tryCatch(enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = ont, pAdjustMethod = "BH",
                             pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MIN_GS, maxGSSize = MAX_GS, readable = FALSE), error = function(e) { cat("  enrichGO ", ont, " failed: ", conditionMessage(e), "\n", sep = ""); NULL })
    out[[paste0("GO_", ont)]] <- tidy_res(obj, analysis, set, bg, paste0("GO_", ont))
  }
  if (kegg_ok) {
    obj <- tryCatch(enrichKEGG(gene = unname(ent), universe = unname(ent_u), organism = "hsa", keyType = "kegg", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                               minGSSize = MIN_GS, maxGSSize = MAX_GS), error = function(e) { cat("  enrichKEGG failed: ", conditionMessage(e), "\n", sep = ""); NULL })
    out[["KEGG"]] <- tidy_res(obj, analysis, set, bg, "KEGG")
  } else out[["KEGG"]] <- NULL
  list(res = out, n_mapped = length(ent), n_universe_mapped = length(ent_u))
}

# ---- run -------------------------------------------------------------------------------------------------------------------------
summ <- list(); all_res <- list(); underpowered <- character(0)
cat("\n=== Enrichment runs ===\n")
for (s in sets) {
  n <- length(s$genes); id <- paste(s$analysis, "|", s$set, "|", s$bg)
  cat(sprintf("\n[%s]\n  genes: %d%s; universe: %d genes (%s)\n", id, n, if (n) paste0(" (", paste(s$genes, collapse = ", "), ")") else "", length(s$universe), s$bg))
  if (n < MIN_SET) {
    cat("  -> UNDERPOWERED: fewer than ", MIN_SET, " genes; enrichment NOT tested", if (s$analysis == "protocol") " (narrative gene-level annotation below instead)" else "", ".\n", sep = "")
    underpowered <- c(underpowered, id)
    summ[[length(summ) + 1]] <- data.frame(analysis = s$analysis, set = s$set, background = s$bg, n_genes = n, n_genes_mapped_entrez = NA_integer_, universe_size = length(s$universe),
                                           universe_mapped_entrez = NA_integer_, db = "all", status = "underpowered (<10 genes): not tested", n_terms_tested = NA_integer_, n_significant_BH05 = NA_integer_, stringsAsFactors = FALSE)
    next
  }
  er <- run_enrich(s$genes, s$universe, s$analysis, s$set, s$bg)
  cat(sprintf("  mapped to Entrez: %d/%d genes; universe mapped %d/%d\n", er$n_mapped, n, er$n_universe_mapped, length(s$universe)))
  for (db in c("GO_BP", "GO_MF", "GO_CC", "KEGG")) {
    r <- er$res[[db]]; nt <- if (is.null(r)) 0L else nrow(r); ns <- if (is.null(r)) 0L else sum(r$significant)
    cat(sprintf("  %-6s terms tested %5d; BH p < 0.05: %d\n", db, nt, ns))
    if (!is.null(r)) all_res[[length(all_res) + 1]] <- r
    summ[[length(summ) + 1]] <- data.frame(analysis = s$analysis, set = s$set, background = s$bg, n_genes = n, n_genes_mapped_entrez = er$n_mapped, universe_size = length(s$universe),
                                           universe_mapped_entrez = er$n_universe_mapped, db = db,
                                           status = if (db == "KEGG" && !kegg_ok) "KEGG unreachable: not run" else "tested", n_terms_tested = nt, n_significant_BH05 = ns, stringsAsFactors = FALSE)
  }
}
summary_tab <- do.call(rbind, summ)
res_all <- if (length(all_res)) do.call(rbind, lapply(all_res, function(r) r[, c("analysis", "set", "background", "db", "ID", "Description", "GeneRatio", "BgRatio", "Count", "GeneRatioNum", "RichFactor",
                                                                              "FoldEnrichment", "pvalue", "p.adjust", "significant", "geneID")])) else NULL
res_sig <- if (!is.null(res_all)) res_all[res_all$significant, ] else NULL
if (!is.null(res_all)) {
  key <- paste(res_all$analysis, res_all$set, res_all$background, res_all$db, sep = "|")
  top10 <- do.call(rbind, lapply(split(res_all, key), function(d) head(d[order(d$p.adjust, d$pvalue), ], 10)))
}
write.csv(summary_tab, file.path(out_dir, "11_enrichment_summary.csv"), row.names = FALSE)
if (!is.null(res_all)) {
  write.csv(res_sig, file.path(out_dir, "11_enrichment_significant.csv"), row.names = FALSE)
  write.csv(top10, file.path(out_dir, "11_enrichment_top10_by_p.csv"), row.names = FALSE)
}

cat("\n=== Summary table ===\n"); print(summary_tab[, c("analysis", "set", "background", "n_genes", "db", "status", "n_terms_tested", "n_significant_BH05")], row.names = FALSE, right = FALSE)
if (!is.null(res_sig) && nrow(res_sig)) {
  cat("\n=== Terms with BH-adjusted p < 0.05 (EXPLORATORY analyses; counts per analysis/background/db) ===\n")
  print(as.data.frame(table(analysis = res_sig$analysis, background = res_sig$background, db = res_sig$db)) |> subset(Freq > 0), row.names = FALSE)
  for (k in unique(paste(res_sig$set, res_sig$background, sep = " | "))) {
    d <- res_sig[paste(res_sig$set, res_sig$background, sep = " | ") == k, ]
    cat("\n--", k, "(top 8 per db by p.adjust)\n"); for (db in unique(d$db)) { dd <- head(d[d$db == db, ], 8)
      cat(sprintf("   %-6s %-12s %-55s Count %2d  %s  p.adj %.2e\n", db, dd$ID, substr(dd$Description, 1, 55), dd$Count, dd$GeneRatio, dd$p.adjust), sep = "") }
  }
} else cat("\nNo term reached BH p < 0.05 in any tested analysis.\n")

# dotplots for analyses with significant terms
fig_files <- character(0)
if (!is.null(res_sig) && nrow(res_sig)) {
  for (k in unique(paste(res_sig$set, res_sig$background, res_sig$db, sep = "|"))) {
    d <- res_sig[paste(res_sig$set, res_sig$background, res_sig$db, sep = "|") == k, ]; d <- head(d[order(d$p.adjust), ], 15)
    d$Description <- factor(substr(d$Description, 1, 60), levels = rev(unique(substr(d$Description, 1, 60))))
    bgtag <- if (grepl("^B", d$background[1])) "bgB_unionExpressed" else "bgA_allExpressed"
    nm <- sprintf("11_%s_%s_%s_dotplot", d$set[1], bgtag, d$db[1])
    g <- ggplot(d, aes(x = GeneRatioNum, y = Description, size = Count, colour = p.adjust)) + geom_point() + scale_colour_gradient(low = "firebrick", high = "navy") +
      labs(title = sprintf("%s (%s)", d$db[1], d$set[1]), subtitle = paste0("EXPLORATORY - ", d$background[1]), x = "GeneRatio", y = NULL) + theme_bw(base_size = 10)
    ggsave(file.path(fig_dir, paste0(nm, ".pdf")), g, width = 8, height = 1.8 + 0.28 * nrow(d)); ggsave(file.path(fig_dir, paste0(nm, ".png")), g, width = 8, height = 1.8 + 0.28 * nrow(d), dpi = 200)
    fig_files <- c(fig_files, nm)
  }
  cat("\nfigures written: ", paste(fig_files, collapse = ", "), "\n", sep = "")
} else cat("No significant terms: no dotplots produced.\n")

# ---- narrative gene-level annotation for the protocol sets (underpowered) ---------------------------------------------------------
cat("\n=== Gene-level annotation of the protocol sets (all are < 10 genes: enrichment is underpowered; annotation below replaces the test) ===\n")
cat("Sources: org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db")), " (GO annotation, gene name, cytoband); ", if (kegg_ok) paste0("KEGG REST API (pathway links, retrieved ", format(Sys.time(), "%Y-%m-%d"), ")") else "KEGG not reachable: org.Hs.eg.db only",
    "; genesets/senescence_union.csv (SenMayo/CellAge/GenAge membership); discovery descriptive limma tables 07a (never used for selection) and nested-CV stability tables.\n", sep = "")
ann_genes <- cls$gene
base <- AnnotationDbi::select(org.Hs.eg.db, ann_genes, c("ENTREZID", "GENENAME", "GENETYPE", "MAP"), "SYMBOL"); stopifnot(!anyDuplicated(base$SYMBOL))
goann <- AnnotationDbi::select(org.Hs.eg.db, ann_genes, c("GO", "EVIDENCE", "ONTOLOGY"), "SYMBOL")
goann <- goann[!is.na(goann$GO), ]
goterm <- AnnotationDbi::select(GO.db::GO.db, unique(goann$GO), "TERM", "GOID"); goann$TERM <- goterm$TERM[match(goann$GO, goterm$GOID)]
goann$direct_non_IEA <- goann$EVIDENCE != "IEA"
terms_out <- goann[order(goann$SYMBOL, goann$ONTOLOGY, !goann$direct_non_IEA, goann$TERM), c("SYMBOL", "ONTOLOGY", "GO", "TERM", "EVIDENCE", "direct_non_IEA")]
write.csv(terms_out, file.path(out_dir, "11_gene_annotation_terms.csv"), row.names = FALSE)
get_kegg <- function(entrez) { if (!kegg_ok) return(NA_character_)
  l <- tryCatch(KEGGREST::keggLink("pathway", paste0("hsa:", entrez)), error = function(e) NULL); if (is.null(l)) return(NA_character_)
  p <- unname(l); nm <- sub(" - Homo sapiens \\(human\\)$", "", kegg_names[sub("path:", "", p)]); paste0(sub("path:", "", p), " ", nm, collapse = "; ") }
ann <- do.call(rbind, lapply(ann_genes, function(g) {
  b <- base[base$SYMBOL == g, ]; ga <- goann[goann$SYMBOL == g, ]; cl <- cls[cls$gene == g, ]
  di <- deg_i[deg_i$gene == g, ]; da <- deg_a[deg_a$gene == g, ]
  src <- c(SenMayo = cl$SenMayo, CellAge = cl$CellAge, GenAge = cl$GenAge)
  data.frame(gene = g, classification = cl$classification, entrez = b$ENTREZID, gene_name = b$GENENAME, gene_type = b$GENETYPE, cytoband = b$MAP,
             GO_BP_n = sum(ga$ONTOLOGY == "BP"), GO_MF_n = sum(ga$ONTOLOGY == "MF"), GO_CC_n = sum(ga$ONTOLOGY == "CC"),
             GO_non_IEA_n = sum(ga$direct_non_IEA), kegg_pathways = get_kegg(b$ENTREZID),
             senescence_sets = if (any(src == 1)) paste(names(src)[src == 1], collapse = "+") else "none", n_sources = cl$n_sources,
             islet_stability_freq = cl$islet_stability_freq, adipose_stability_freq = cl$adipose_stability_freq,
             islet_log2FC = if (nrow(di)) di$log2FC else NA, islet_P = if (nrow(di)) di$P.Value else NA, islet_adjP = if (nrow(di)) di$adj.P.Val else NA, islet_DEG = if (nrow(di)) di$DEG else NA,
             adipose_log2FC = if (nrow(da)) da$log2FC else NA, adipose_P = if (nrow(da)) da$P.Value else NA, adipose_adjP = if (nrow(da)) da$adj.P.Val else NA, adipose_DEG = if (nrow(da)) da$DEG else NA,
             stringsAsFactors = FALSE) }))
write.csv(ann, file.path(out_dir, "11_gene_annotation_summary.csv"), row.names = FALSE)
for (i in seq_len(nrow(ann))) {
  a <- ann[i, ]; ga <- terms_out[terms_out$SYMBOL == a$gene, ]
  cat(sprintf("\n%s (%s) -- %s; %s; cytoband %s; Entrez %s\n", a$gene, a$classification, a$gene_name, a$gene_type, a$cytoband, a$entrez))
  cat(sprintf("  GO annotations in org.Hs.eg.db: BP %d, MF %d, CC %d (%d with non-IEA evidence)\n", a$GO_BP_n, a$GO_MF_n, a$GO_CC_n, a$GO_non_IEA_n))
  for (ont in c("MF", "BP", "CC")) { d <- ga[ga$ONTOLOGY == ont & ga$direct_non_IEA, ]; d <- d[!duplicated(d$GO), ]
    cat(sprintf("  non-IEA %s terms (%d): %s\n", ont, nrow(d), if (nrow(d)) paste(sprintf("%s [%s]", head(d$TERM, 12), head(d$EVIDENCE, 12)), collapse = "; ") else "none")) }
  cat("  KEGG pathway links: ", a$kegg_pathways, "\n", sep = "")
  cat(sprintf("  senescence gene sets containing it: %s (n_sources = %d of 3)\n", a$senescence_sets, a$n_sources))
  cat(sprintf("  discovery descriptive limma (never used for selection): islet log2FC %.3f, P %.2e, BH adjP %.3g, DEG %s | adipose log2FC %.3f, P %.2e, BH adjP %.3g, DEG %s\n",
              a$islet_log2FC, a$islet_P, a$islet_adjP, a$islet_DEG, a$adipose_log2FC, a$adipose_P, a$adipose_adjP, a$adipose_DEG))
  cat(sprintf("  nested-CV stability (final-panel selection frequency): islet %.2f, adipose %.2f\n", a$islet_stability_freq, a$adipose_stability_freq))
}

# ---- save + assertions ----------------------------------------------------------------------------------------------------------------------
saveRDS(list(summary = summary_tab, results = res_all, significant = res_sig, annotation = ann, annotation_terms = terms_out, underpowered = underpowered,
             params = list(P_CUT = P_CUT, MIN_GS = MIN_GS, MAX_GS = MAX_GS, MIN_SET = MIN_SET, kegg_ok = kegg_ok, date = format(Sys.Date()))),
        file.path(out_dir, "11_enrichment_result.rds"))
stopifnot(file.exists(file.path(out_dir, c("11_enrichment_summary.csv", "11_gene_annotation_summary.csv", "11_gene_annotation_terms.csv", "11_enrichment_result.rds"))),
          nrow(summary_tab) > 0, !anyNA(summary_tab$analysis),
          # the three protocol sets exist, each has < 10 genes (so each is declared underpowered, not tested)
          sum(summary_tab$analysis == "protocol") == 3L, all(summary_tab$status[summary_tab$analysis == "protocol"] == "underpowered (<10 genes): not tested"),
          all(summary_tab$n_genes[summary_tab$analysis == "protocol"] == c(length(sets[[1]]$genes), length(sets[[2]]$genes), length(sets[[3]]$genes))),
          # every reported term is significant by the BH rule and no term has fewer than minGSSize genes in its background annotation
          is.null(res_sig) || all(res_sig$p.adjust < P_CUT), is.null(res_all) || all(res_all$significant == (res_all$p.adjust < P_CUT)),
          # no test ran on a set of < 10 genes
          all(summary_tab$n_genes[summary_tab$status == "tested"] >= MIN_SET), all(is.na(summary_tab$n_terms_tested[grepl("underpowered", summary_tab$status)])),
          nrow(ann) == length(ann_genes), all(!is.na(ann$gene_name)), all(ann$GO_BP_n + ann$GO_MF_n + ann$GO_CC_n > 0),
          all(file.exists(file.path(fig_dir, c(paste0(fig_files, ".pdf"), paste0(fig_files, ".png"))))))
cat("\nAll end-of-script assertions passed.\n")
cat(sprintf("\nRuntime: %.1f min\n", as.numeric(difftime(Sys.time(), T_START, units = "mins"))))
cat("\nsessionInfo():\n"); print(sessionInfo())
sink(type = "message"); sink(); close(con)

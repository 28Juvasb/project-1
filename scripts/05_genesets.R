# 05_genesets.R -- senescence gene universe (protocol v1.1 §3.1, pre-analysis plan B9).
#
# Inputs : genesets/SenMayo.csv, genesets/CellAge.csv, genesets/GenAge_human.csv (read-only)
#          data/processed/03_<GSE>_{expr,expressed}.rds  (gene names / flags only; no group labels, no expression comparisons)
# Outputs: genesets/senescence_union.rds   character vector of unique HGNC symbols (UNION of the three sets)
#          genesets/senescence_union.csv   membership flags per source set + n_sources
#          results/stage0/05_symbol_mapping.csv   every changed / dropped symbol
#          results/stage0/05_set_summary.csv      per-set mapping outcome, sizes
#          results/stage0/05_overlaps.csv         pairwise, 3-way, union
#          results/stage0/05_union_coverage.csv   union genes present / expressed per processed matrix
#          results/stage0/05_genesets.rds         all intermediate objects
#          results/stage0/figures/05_union_upset.{pdf,png}
#          logs/05_genesets_<timestamp>.log
#
# SenMayo (B9): the Saul et al. 2022 paper (PMID 35974106) describes 125 genes; MSigDB SAUL_SEN_MAYO has 124. The paper's
#   Supplementary Data 1 could not be retrieved in this environment (nature.com redirects to a cookie/IDP challenge, PMC page
#   exposed no supplementary link), so the 124-gene MSigDB copy is used and the discrepancy remains UNRESOLVED. No gene is guessed.
#
# Mapping policy: limma::alias2SymbolTable(species = "Hs") (symbol source: org.Hs.eg.db). Symbols that are not current official
#   symbols and map to >1 current symbol (ambiguous alias) are DROPPED, never resolved by guesswork; unmappable symbols are dropped.

set.seed(42)
suppressPackageStartupMessages({
  library(here); library(limma); library(org.Hs.eg.db); library(AnnotationDbi); library(UpSetR)
})

script_name <- "05_genesets"
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- here("results", "stage0"); fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(here("logs"), showWarnings = FALSE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt")
sink(con, split = TRUE); sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

SENMAYO_EXPECTED <- 124L   # MSigDB copy; paper's 125 not retrievable (see header)
DATASETS <- c(GSE41762 = "discovery islet", GSE27951 = "discovery adipose", GSE76894 = "validation islet (primary)",
              GSE71416 = "validation adipose", GSE20966 = "validation islet (secondary)", GSE16415 = "validation adipose")
SETS <- c(SenMayo = "SenMayo.csv", CellAge = "CellAge.csv", GenAge = "GenAge_human.csv")

cat("Script:", script_name, " run:", ts, "\n\nInputs:\n")
for (f in SETS) { p <- here("genesets", f); cat("  ", sub(paste0(here(), "/"), "", p), " md5=", tools::md5sum(p), "\n", sep = "") }
cat("\nParameters: limma::alias2SymbolTable(species='Hs'), org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db")),
    "| limma", as.character(packageVersion("limma")), "| ambiguous aliases dropped | unmappable dropped | union (not intersection)\n")
cat("SenMayo: 125-gene paper list (Saul 2022 Suppl. Data 1) NOT retrievable -> MSigDB 124-gene copy used; discrepancy UNRESOLVED.\n\n")

official <- keys(org.Hs.eg.db, keytype = "SYMBOL")

# ---- 1-3. load and map ----------------------------------------------------------------------------------------
map_set <- function(set, file) {
  d <- read.csv(here("genesets", file), stringsAsFactors = FALSE)
  inp <- unique(trimws(d$gene_symbol))
  stopifnot(!anyNA(inp), all(nzchar(inp)))
  mapped <- alias2SymbolTable(inp, species = "Hs")
  is_off <- inp %in% official
  # ambiguity: non-official input whose alias resolves to >1 current symbol
  nres <- integer(length(inp))
  for (i in which(!is_off)) nres[i] <- length(alias2Symbol(inp[i], species = "Hs"))
  status <- ifelse(is_off & !is.na(mapped) & mapped == inp, "unchanged",
            ifelse(!is_off & nres > 1, "ambiguous",
            ifelse(is.na(mapped), "unmappable", "updated_alias")))
  # official symbol that alias2SymbolTable nonetheless changed (should not occur) -> treat as updated, but flag
  final <- ifelse(status %in% c("unchanged", "updated_alias"), mapped, NA_character_)
  data.frame(set = set, input = inp, mapped = mapped, final = final, status = status,
             n_candidates = nres, stringsAsFactors = FALSE)
}
maps <- do.call(rbind, Map(map_set, names(SETS), SETS))

raw_rows <- c(SenMayo = nrow(read.csv(here("genesets", SETS[["SenMayo"]]))),
              CellAge = nrow(read.csv(here("genesets", SETS[["CellAge"]]))),
              GenAge  = nrow(read.csv(here("genesets", SETS[["GenAge"]]))))

# cross-check updated aliases against supplied Entrez IDs (descriptive sanity check only)
ent <- lapply(c(CellAge = "CellAge.csv", GenAge = "GenAge_human.csv"), function(f) {
  d <- read.csv(here("genesets", f), stringsAsFactors = FALSE); d[!duplicated(d$gene_symbol), c("gene_symbol", "entrez_id")] })
sym2eg <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = official, keytype = "SYMBOL", columns = "ENTREZID"))
maps$entrez_supplied <- NA_character_; maps$entrez_of_mapped <- NA_character_
for (s in names(ent)) {
  i <- maps$set == s
  maps$entrez_supplied[i] <- as.character(ent[[s]]$entrez_id[match(maps$input[i], ent[[s]]$gene_symbol)])
}
maps$entrez_of_mapped <- sym2eg$ENTREZID[match(maps$final, sym2eg$SYMBOL)]
maps$entrez_agrees <- ifelse(is.na(maps$entrez_supplied) | is.na(maps$entrez_of_mapped), NA,
                             maps$entrez_supplied == maps$entrez_of_mapped)

summ <- do.call(rbind, lapply(names(SETS), function(s) {
  m <- maps[maps$set == s, ]
  fin <- unique(m$final[!is.na(m$final)])
  data.frame(set = s, rows_in_file = raw_rows[[s]], input_symbols = nrow(m),
             mapped_unchanged = sum(m$status == "unchanged"), updated_alias = sum(m$status == "updated_alias"),
             ambiguous_dropped = sum(m$status == "ambiguous"), unmappable_dropped = sum(m$status == "unmappable"),
             final_unique_symbols = length(fin),
             collapsed_duplicates = sum(m$status %in% c("unchanged", "updated_alias")) - length(fin),
             stringsAsFactors = FALSE)
}))
cat("Mapping outcome per set:\n"); print(summ, row.names = FALSE)
chg <- maps[maps$status != "unchanged", c("set", "input", "mapped", "status", "n_candidates", "entrez_supplied", "entrez_of_mapped", "entrez_agrees")]
names(chg)[3] <- "mapped_symbol"
write.csv(chg, file.path(out_dir, "05_symbol_mapping.csv"), row.names = FALSE, na = "")
cat("\nChanged or dropped symbols (", nrow(chg), "), written to results/stage0/05_symbol_mapping.csv:\n", sep = ""); print(chg, row.names = FALSE)
write.csv(summ, file.path(out_dir, "05_set_summary.csv"), row.names = FALSE)

# ---- 4. union --------------------------------------------------------------------------------------------------
sets <- lapply(setNames(names(SETS), names(SETS)), function(s) sort(unique(maps$final[maps$set == s & !is.na(maps$final)])))
union_syms <- sort(unique(unlist(sets, use.names = FALSE)))
uni <- data.frame(gene_symbol = union_syms, stringsAsFactors = FALSE)
for (s in names(sets)) uni[[s]] <- as.integer(union_syms %in% sets[[s]])
uni$n_sources <- rowSums(uni[, names(sets)])
saveRDS(union_syms, here("genesets", "senescence_union.rds"))
write.csv(uni, here("genesets", "senescence_union.csv"), row.names = FALSE)

# ---- 5. overlaps + UpSet --------------------------------------------------------------------------------------
pair <- combn(names(sets), 2, simplify = FALSE)
ov <- rbind(
  data.frame(comparison = names(sets), n = lengths(sets)),
  do.call(rbind, lapply(pair, function(p) data.frame(comparison = paste(p, collapse = " & "), n = length(intersect(sets[[p[1]]], sets[[p[2]]]))))),
  data.frame(comparison = "SenMayo & CellAge & GenAge (3-way)", n = length(Reduce(intersect, sets))),
  data.frame(comparison = "Union", n = length(union_syms)),
  data.frame(comparison = paste0("only ", names(sets)), n = sapply(names(sets), function(s) sum(uni[[s]] == 1 & uni$n_sources == 1)))
)
rownames(ov) <- NULL
write.csv(ov, file.path(out_dir, "05_overlaps.csv"), row.names = FALSE)
cat("\nSet sizes and overlaps:\n"); print(ov, row.names = FALSE)
cat("\nGenes by number of sources:\n"); print(table(uni$n_sources))

draw_upset <- function() {
  print(upset(fromList(sets), sets = names(sets), order.by = "freq", keep.order = FALSE,
              mainbar.y.label = "Intersection size", sets.x.label = "Genes per set", text.scale = 1.3))
}
pdf(file.path(fig_dir, "05_union_upset.pdf"), width = 7, height = 5); draw_upset(); invisible(dev.off())
png(file.path(fig_dir, "05_union_upset.png"), width = 7, height = 5, units = "in", res = 150); draw_upset(); invisible(dev.off())

# ---- 6. coverage (gene names / flags only) -------------------------------------------------------------------
cov <- do.call(rbind, lapply(names(DATASETS), function(g) {
  ex <- readRDS(here("data", "processed", sprintf("03_%s_expressed.rds", g)))
  stopifnot(is.data.frame(ex), all(c("gene", "expressed") %in% names(ex)))
  genes <- rownames(readRDS(here("data", "processed", sprintf("03_%s_expr.rds", g))))
  stopifnot(setequal(genes, ex$gene))
  pres <- union_syms %in% genes
  expr <- union_syms %in% ex$gene[ex$expressed %in% TRUE]
  data.frame(dataset = g, role = DATASETS[[g]], matrix_genes = length(genes), union_size = length(union_syms),
             union_present = sum(pres), pct_present = round(100 * mean(pres), 1),
             union_expressed = sum(expr), pct_expressed = round(100 * mean(expr), 1),
             flag_lt80pct_present = mean(pres) < 0.8, stringsAsFactors = FALSE)
}))
write.csv(cov, file.path(out_dir, "05_union_coverage.csv"), row.names = FALSE)
cat("\nUnion coverage per processed matrix (names/flags only):\n"); print(cov, row.names = FALSE)

saveRDS(list(maps = maps, summary = summ, sets = sets, union = union_syms, union_table = uni, overlaps = ov, coverage = cov),
        file.path(out_dir, "05_genesets.rds"))

# ---- 7. assertions ---------------------------------------------------------------------------------------------
cat("\nAssertions:\n")
chk <- function(ok, msg) { cat(if (isTRUE(ok)) "  PASS " else "  FAIL ", msg, "\n", sep = ""); if (!isTRUE(ok)) stop("Assertion failed: ", msg) }
chk(length(union_syms) == length(unique(maps$final[!is.na(maps$final)])), "union size == number of unique mapped symbols across sets")
chk(!anyNA(union_syms) && all(nzchar(union_syms)) && !anyDuplicated(union_syms), "union has no NA / empty / duplicates")
chk(all(union_syms %in% official), "all union symbols are valid org.Hs.eg.db symbols")
chk(length(sets$SenMayo) == SENMAYO_EXPECTED, sprintf("SenMayo count == %d (MSigDB copy; paper's 125 unresolved)", SENMAYO_EXPECTED))
chk(identical(readRDS(here("genesets", "senescence_union.rds")), union_syms), "senescence_union.rds round-trips")
chk(all(uni$n_sources >= 1 & uni$n_sources <= 3) && nrow(uni) == length(union_syms), "union csv: every gene has 1-3 sources")
chk(length(union_syms) >= max(lengths(sets)), "union >= largest set (union, not intersection)")
chk(all(summ$unmappable_dropped + summ$ambiguous_dropped + summ$mapped_unchanged + summ$updated_alias == summ$input_symbols), "mapping outcomes sum to input symbols")
chk(nrow(cov) == length(DATASETS) && !anyNA(cov$union_present), "coverage computed for all 6 datasets")
chk(file.exists(file.path(fig_dir, "05_union_upset.pdf")) && file.exists(file.path(fig_dir, "05_union_upset.png")), "UpSet figures written")

cat("\nsessionInfo():\n"); print(sessionInfo())

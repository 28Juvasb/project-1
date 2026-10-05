# stage1_helpers.R -- shared functions for scripts/03_normalize.R and scripts/04_qc.R
# (protocol §2.1-2.4). Sourced, never run on its own. No stochastic steps here.
#
# Raw archives in data/raw/ are only READ. CEL members are extracted to
# data/processed/cel_cache/<GSE>/ (git-ignored), never inside data/raw/.

suppressPackageStartupMessages({ library(here); library(data.table) })

STAGE1 <- list(
  datasets   = c("GSE41762", "GSE27951", "GSE76894", "GSE71416", "GSE20966", "GSE16415"),
  role       = c(GSE41762 = "discovery", GSE27951 = "discovery", GSE76894 = "validation",
                 GSE71416 = "validation", GSE20966 = "validation", GSE16415 = "validation"),
  platform   = c(GSE41762 = "GPL6244", GSE27951 = "GPL570", GSE76894 = "GPL570",
                 GSE71416 = "GPL570", GSE20966 = "GPL1352", GSE16415 = "GPL2986"),
  method     = c(GSE41762 = "oligo", GSE27951 = "affy", GSE76894 = "affy",
                 GSE71416 = "affy", GSE20966 = "affy", GSE16415 = "matrix"),
  annot_db   = c(GPL6244 = "hugene10sttranscriptcluster.db", GPL570 = "hgu133plus2.db",
                 GPL1352 = "u133x3p.db", GPL2986 = "GPL2986 annotation table (data/raw/GPL2986)"),
  cdf        = c(GSE27951 = "hgu133plus2cdf", GSE76894 = "hgu133plus2cdf", GSE71416 = "hgu133plus2cdf",
                 GSE20966 = "u133x3pcdf"),
  # protocol v1.1 / pre_analysis_plan B11: possible donor overlap, excluded before scoring
  gse76894_excl = c("GSM2040064", "GSM2040077", "GSM2040120", "GSM2040122"),
  expr_filter_quantile = 0.25,   # §2.2
  expr_filter_min_frac = 0.20,   # §2.2
  nuse_max   = 1.05,             # §2.3
  pca_top    = 5000,             # §2.3
  pca_sd     = 3,                # §2.3
  y_genes    = c("RPS4Y1", "DDX3Y", "KDM5D", "UTY")  # §2.4
)

std_sex <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[x %in% c("m", "male")] <- "male"
  out[x %in% c("f", "female")] <- "female"
  out
}
num <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("", "--", "NA", "null"), NA, x)))

# ---- sample selection + phenotype table from manifest --------------------------------
build_pheno <- function(gse) {
  if (gse == "GSE41762") {
    d <- fread(here("manifest", "GSE41762_samples.csv"), colClasses = "character")
    p <- data.frame(gsm = d$geo_accession, gse = gse, title = d$title, group = d$group,
                    age = num(d[["age:ch1"]]), sex = std_sex(d[["Sex:ch1"]]), bmi = num(d[["bmi:ch1"]]),
                    hba1c = num(d[["hba1c:ch1"]]), donor_num = as.integer(d$donor_num),
                    submitter_norm_batch = d$submitter_norm_batch, batch = d$submitter_norm_batch)
  } else if (gse == "GSE27951") {
    d <- fread(here("manifest", "GSE27951_samples_GPL570.csv"), colClasses = "character")
    d <- d[group %in% c("ND", "T2D") & (is.na(excluded_reason) | excluded_reason == "")]
    p <- data.frame(gsm = d$geo_accession, gse = gse, title = d$title, group = d$group,
                    age = num(d[["age:ch1"]]), sex = NA_character_, bmi = num(d[["bmi:ch1"]]),
                    hba1c = num(d[["hba1c:ch1"]]), fasting_glucose = num(d[["fasting glucose:ch1"]]),
                    fasting_insulin = num(d[["fasting insulin:ch1"]]), clinical_status = d[["clinical status:ch1"]],
                    title_status_discordant = as.logical(d$title_status_discordant), batch = NA_character_)
  } else if (gse == "GSE76894") {
    d <- fread(here("manifest", "GSE76894_samples.csv"), colClasses = "character")
    d <- d[!geo_accession %in% STAGE1$gse76894_excl]
    p <- data.frame(gsm = d$geo_accession, gse = gse, title = d$title, group = d$group,
                    age = num(d[["age (years):ch1"]]), sex = std_sex(d[["gender (m/f):ch1"]]),
                    bmi = num(d[["body mass index (bmi):ch1"]]),
                    icu_glycemia = num(d[["icu (intensive care unit) glycemia (mmol/l):ch1"]]),
                    batch = d[["batch:ch1"]])
  } else if (gse == "GSE71416") {
    d <- fread(here("manifest", "GSE71416_samples.csv"), colClasses = "character")
    p <- data.frame(gsm = d$geo_accession, gse = gse, title = d$title, group = d$group,
                    age = NA_real_, sex = std_sex(d[["Sex:ch1"]]), bmi = NA_real_, batch = NA_character_)
  } else if (gse == "GSE20966") {
    d <- fread(here("manifest", "GSE20966_samples.csv"), colClasses = "character")
    p <- data.frame(gsm = d$geo_accession, gse = gse, title = d$title, group = d$group,
                    age = num(d[["age (years):ch1"]]), sex = std_sex(d[["gender:ch1"]]),
                    bmi = num(d[["bmi (body mass index, kg/m2):ch1"]]), sample_id = d[["sample id:ch1"]],
                    batch = NA_character_)
  } else if (gse == "GSE16415") {
    d <- fread(here("manifest", "GSE16415_samples.csv"), colClasses = "character")
    p <- data.frame(gsm = d$geo_accession, gse = gse, title = d$title, group = d$group,
                    age = NA_real_, sex = NA_character_, bmi = NA_real_, batch = NA_character_)
  } else stop("unknown GSE ", gse)
  rownames(p) <- p$gsm
  p
}

# ---- CEL extraction (read-only on data/raw) ------------------------------------------
cel_files <- function(gse, gsms) {
  tarf <- here("data", "raw", gse, paste0(gse, "_RAW.tar"))
  stopifnot(file.exists(tarf))
  members <- untar(tarf, list = TRUE)
  members <- members[grepl("\\.CEL\\.gz$", members, ignore.case = TRUE)]   # CEL only (GSE71416 rma.chp ignored)
  m_gsm <- sub("^(GSM\\d+).*", "\\1", basename(members))
  if (anyDuplicated(m_gsm)) stop(gse, ": more than one CEL per GSM in archive")
  miss <- setdiff(gsms, m_gsm)
  if (length(miss)) stop(gse, ": no CEL in archive for ", paste(miss, collapse = ","))
  sel <- members[match(gsms, m_gsm)]
  outdir <- here("data", "processed", "cel_cache", gse)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  need <- sel[!file.exists(file.path(outdir, sel))]
  if (length(need)) untar(tarf, files = need, exdir = outdir)
  f <- file.path(outdir, sel)
  stopifnot(all(file.exists(f)))
  names(f) <- gsms
  f
}

read_raw <- function(gse, gsms) {
  f <- cel_files(gse, gsms)
  if (STAGE1$method[gse] == "oligo") {
    raw <- suppressMessages(oligo::read.celfiles(f, sampleNames = gsms, verbose = FALSE))
  } else {
    raw <- affy::ReadAffy(filenames = f, sampleNames = gsms, cdfname = sub("cdf$", "", STAGE1$cdf[gse]),
                          verbose = FALSE)
  }
  raw
}

# ---- normalisation §2.1 -----------------------------------------------------------------
normalize_cel <- function(gse, gsms) {
  raw <- read_raw(gse, gsms)
  if (STAGE1$method[gse] == "oligo") {
    eset <- oligo::rma(raw, target = "core")
  } else {
    stopifnot(affy::cdfName(raw) == sub("cdf$", "", STAGE1$cdf[gse]) ||
              affy::cleancdfname(affy::cdfName(raw)) == STAGE1$cdf[gse])
    eset <- affy::rma(raw)
  }
  m <- Biobase::exprs(eset)
  colnames(m) <- gsms
  rm(raw, eset); gc()
  m
}

gse16415_matrix <- function(pheno) {
  f <- here("data", "raw", "GSE16415", "GSE16415_all_data_no_controls.txt.gz")
  d <- fread(f, skip = 2, sep = "\t", header = TRUE, quote = "\"", colClasses = list(character = c("PROBE_ID")))
  sig <- grep("_SIGNAL$", names(d), value = TRUE)          # raw signal columns (not 'Trim Mean Scaled')
  stopifnot(length(sig) == 10)
  hyb <- sub("_.*$", "", sig)
  # GEO title for GSM412611 is 'HBOOGA7' (letter O) while the matrix header is 'HB00GA7' (zero):
  title_key <- gsub("O", "0", toupper(pheno$title))
  idx <- match(title_key, gsub("O", "0", toupper(hyb)))
  if (anyNA(idx)) stop("GSE16415: cannot map titles to matrix columns: ", paste(pheno$title[is.na(idx)], collapse = ","))
  m <- as.matrix(d[, ..sig])[, idx]
  storage.mode(m) <- "double"
  rownames(m) <- d$PROBE_ID
  colnames(m) <- pheno$gsm
  m <- m[!is.na(rownames(m)) & rownames(m) != "", , drop = FALSE]
  attr(m, "col_map") <- data.frame(gsm = pheno$gsm, title = pheno$title, matrix_column = sig[idx])
  m
}

# ---- annotation §2.2 -------------------------------------------------------------------
hugene_main_tc <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), system.file("extdata", "pd.hugene.1.0.st.v1.sqlite",
                                                        package = "pd.hugene.1.0.st.v1"), flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con))
  t <- DBI::dbGetQuery(con, "select distinct core_mps.meta_fsetid as tc, featureSet.type as type
                             from core_mps inner join featureSet on core_mps.fsetid = featureSet.fsetid")
  t <- as.data.table(t)[, .(all_main = all(!is.na(type) & type == 1)), by = tc]
  as.character(t$tc[t$all_main])
}

probe_to_symbol <- function(gse, probes) {
  plat <- STAGE1$platform[gse]
  log <- list(n_probes = length(probes))
  if (plat == "GPL2986") {
    tab <- read_gpl_table(here("data", "raw", "GPL2986", "GPL2986.soft.gz"))
    sym <- tab[["Gene Symbol"]][match(probes, tab$ID)]
    sym[sym %in% c("", "null", "NA")] <- NA
    ctrl <- rep(FALSE, length(probes))
    multi <- !is.na(sym) & grepl("[,;/ ]", sym)
  } else {
    pkg <- STAGE1$annot_db[plat]
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    db <- get(pkg)
    ctrl <- grepl("^AFFX", probes, ignore.case = TRUE)
    if (plat == "GPL6244") {
      main <- hugene_main_tc()
      ctrl <- ctrl | !(probes %in% main)
    }
    mp <- suppressMessages(AnnotationDbi::mapIds(db, keys = probes, column = "SYMBOL", keytype = "PROBEID",
                                                 multiVals = "list"))
    nsym <- vapply(mp, function(v) length(unique(v[!is.na(v)])), integer(1))
    sym <- vapply(mp, function(v) { v <- unique(v[!is.na(v)]); if (length(v) == 1) v else NA_character_ }, character(1))
    multi <- nsym > 1
    sym[nsym == 0] <- NA
  }
  keep <- !ctrl & !multi & !is.na(sym)
  log$n_control <- sum(ctrl)
  log$n_multi <- sum(multi & !ctrl)
  log$n_unmapped <- sum(is.na(sym) & !ctrl & !multi)
  log$n_kept_probes <- sum(keep)
  list(symbol = unname(sym), keep = keep, log = log)
}

read_gpl_table <- function(f) {
  x <- readLines(f, warn = FALSE)            # file is plain text despite .gz name; readLines handles both
  x <- sub("\r$", "", x)
  b <- which(x == "!platform_table_begin"); e <- which(x == "!platform_table_end")
  stopifnot(length(b) == 1, length(e) == 1)
  fread(text = paste(x[(b + 1):(e - 1)], collapse = "\n"), sep = "\t", quote = "", colClasses = "character")
}

collapse_to_genes <- function(gse, m) {
  a <- probe_to_symbol(gse, rownames(m))
  mk <- m[a$keep, , drop = FALSE]
  sk <- a$symbol[a$keep]
  cr <- WGCNA::collapseRows(mk, rowGroup = sk, rowID = rownames(mk), method = "MaxMean")
  g <- cr$datETcollapsed
  stopifnot(!anyDuplicated(rownames(g)))
  a$log$n_genes <- nrow(g)
  attr(g, "collapse_log") <- a$log
  attr(g, "selected_probe") <- setNames(cr$group2row[, "selectedRowID"], cr$group2row[, "group"])
  g
}

# ---- expression filter §2.2 -------------------------------------------------------------
expressed_flags <- function(g) {
  q <- as.numeric(quantile(g, STAGE1$expr_filter_quantile, na.rm = TRUE))
  frac <- rowMeans(g > q)
  data.frame(gene = rownames(g), frac_above_q25 = frac, expressed = frac >= STAGE1$expr_filter_min_frac,
             q25_threshold = q, row.names = NULL)
}

# ---- one complete normalise+collapse run for a dataset ---------------------------------
run_normalisation <- function(gse, pheno) {
  if (STAGE1$method[gse] == "matrix") {
    m <- gse16415_matrix(pheno)
    col_map <- attr(m, "col_map")
    mx <- max(m, na.rm = TRUE)
    logged <- mx > 50
    if (logged) m <- log2(m + 1)
    m <- limma::normalizeBetweenArrays(m, method = "quantile")
    info <- list(max_raw = mx, log2_applied = logged, col_map = col_map)
  } else {
    m <- normalize_cel(gse, pheno$gsm)
    info <- list()
  }
  stopifnot(identical(colnames(m), pheno$gsm), !anyNA(m))
  info$n_probes_normalised <- nrow(m)
  g <- collapse_to_genes(gse, m)
  ef <- expressed_flags(g)
  list(expr = g, expressed = ef, info = c(info, attr(g, "collapse_log")))
}

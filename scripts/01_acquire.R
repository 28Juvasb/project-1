# 01_acquire.R -- download GEO series, platform annotations and sample manifests.
# Raw files in data/raw/ are NEVER overwritten (existing files are skipped).
# Stage 1a: download + manifests (8 series incl. GSE76894, added for protocol v1.1).  Stage 1b (verification) is appended below once
# the metadata field names have been inspected.

set.seed(42)
suppressPackageStartupMessages({ library(here); library(GEOquery); library(Biobase); library(data.table) })
options(timeout = 7200)

script_name <- "01_acquire"
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create(here("logs"), showWarnings = FALSE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt")
sink(con, split = TRUE); sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

series <- c("GSE41762", "GSE38642", "GSE20966", "GSE81608", "GSE27951", "GSE16415", "GSE71416", "GSE76894")
# GSE76894 (protocol v1.1, primary islet external validation): only the CEL archive is fetched.
# The submitters' GSE76894_DifferentialExpression-T2DvsND.txt.gz (an expression-derived result on a
# validation set) is deliberately NOT downloaded (CLAUDE.md rule 2).
supp_filter <- c(GSE76894 = "_RAW\\.tar$")
platforms <- c("GPL6244", "GPL570", "GPL1352", "GPL2986")
raw <- here("data", "raw"); dir.create(raw, showWarnings = FALSE, recursive = TRUE)
dir.create(here("manifest"), showWarnings = FALSE)
cat("Inputs: GEO accessions", paste(series, collapse = ", "), "\nPlatforms:", paste(platforms, collapse = ", "), "\n\n")

parse_pheno <- function(es) {
  p <- as.data.table(pData(es))
  keep <- c("geo_accession", "title", "platform_id")
  ch <- grep("^characteristics_ch1", names(p), value = TRUE)
  parsed <- grep(":ch1$", names(p), value = TRUE)
  p[, c(keep, parsed, ch), with = FALSE]
}

for (g in series) {
  d <- file.path(raw, g); dir.create(d, showWarnings = FALSE)
  cat("==== ", g, "\n")
  # series matrix + metadata (getGEO skips download if the file already exists in destdir)
  # idempotent/offline: if series matrix files are already on disk, parse them locally (getGEO(accession) always lists the FTP dir)
  sm_local <- list.files(d, pattern = "series_matrix\\.txt\\.gz$", full.names = TRUE)
  if (length(sm_local)) {
    gl <- lapply(sort(sm_local), function(f) getGEO(filename = f, GSEMatrix = TRUE, getGPL = FALSE))
  } else gl <- getGEO(g, GSEMatrix = TRUE, destdir = d, getGPL = FALSE)
  for (i in seq_along(gl)) {
    es <- gl[[i]]
    m <- parse_pheno(es)
    out <- here("manifest", if (length(gl) == 1) sprintf("%s_samples.csv", g) else sprintf("%s_samples_%s.csv", g, annotation(es)))
    fwrite(m, out)
    cat(sprintf("  manifest %s: %d samples, platform %s\n", basename(out), nrow(m), paste(unique(m$platform_id), collapse = ",")))
  }
  # supplementary files (RAW.tar / count matrices): skip if directory already has files
  have <- list.files(d, pattern = "RAW|\\.txt\\.gz$|\\.csv", full.names = FALSE)
  have <- setdiff(have, grep("series_matrix", have, value = TRUE))
  if (length(have) == 0) {
    flt <- if (g %in% names(supp_filter)) supp_filter[[g]] else NULL
    r <- tryCatch(getGEOSuppFiles(g, baseDir = raw, makeDirectory = TRUE, fetch_files = TRUE, filter_regex = flt),
                  error = function(e) { cat("  SUPP DOWNLOAD FAILED:", conditionMessage(e), "\n"); NULL })
    if (!is.null(r)) cat("  supplementary:", paste(basename(rownames(r)), collapse = ", "), "\n")
  } else cat("  supplementary already present:", paste(have, collapse = ", "), "\n")
}

for (p in platforms) {
  d <- file.path(raw, p); dir.create(d, showWarnings = FALSE)
  cat("==== ", p, "\n")
  tryCatch({ if (file.exists(file.path(d, paste0(p, ".soft.gz")))) cat("  already present\n") else { getGEO(p, destdir = d); cat("  ok\n") } },
           error = function(e) cat("  PLATFORM DOWNLOAD FAILED:", conditionMessage(e), "\n"))
}
cat("\nStage 1a complete.\n")

# =====================================================================================
# Stage 1b: derive standardised group, verify counts vs protocol, verify files, MD5.
# Nothing is corrected: mismatches are reported in full, outputs are still written, then
# the script stops with a non-zero status.
# =====================================================================================
cat("\n==== Stage 1b: group derivation and verification\n")
problems <- character(0)
rd <- function(f) fread(here("manifest", f), na.strings = c("", "NA"))
wr <- function(x, f) fwrite(x, here("manifest", f))

# Group is derived from exact observed metadata strings only. Any string not in the map
# gives group = NA and is listed in `unclassified`.
classify <- function(x, nd, t2d) ifelse(x %in% nd, "ND", ifelse(x %in% t2d, "T2D", NA_character_))

m <- list(); note <- list()
# GSE41762: status:ch1
x <- rd("GSE41762_samples.csv")
x[, group := classify(`status:ch1`, "Non-diabetic donor", "Diabetic donor")]
x[, donor_num := as.integer(sub(".*ID\\s*([0-9]+).*", "\\1", title))]
x[, submitter_norm_batch := ifelse(is.na(donor_num), NA_character_, ifelse(donor_num <= 48, "donors_1-48", "donors_49-80"))]
m$GSE41762 <- x
# GSE38642: status:ch1
x <- rd("GSE38642_samples.csv")
x[, group := classify(`status:ch1`, "non-diabetic donors", "T2D donors")]
m$GSE38642 <- x
# GSE20966: disease:ch1
x <- rd("GSE20966_samples.csv")
x[, group := classify(`disease:ch1`, "non-diabetic control", "type 2 diabetes")]
m$GSE20966 <- x
# GSE16415: disease state:ch1
x <- rd("GSE16415_samples.csv")
x[, group := classify(`disease state:ch1`, "control", "type-2 diabetes mellitus")]
m$GSE16415 <- x
# GSE71416: no disease field in characteristics (only Sex); status is encoded ONLY in sample titles
x <- rd("GSE71416_samples.csv")
x[, group := ifelse(grepl("^OAT_Obese non-diabetic", title), "ND", ifelse(grepl("^OAT_Obese and diabetic", title), "T2D", NA_character_))]
x[, group_source := "title (no disease characteristic field in GEO metadata)"]
m$GSE71416 <- x
# GSE81608 (single cells): condition:ch1
x <- rd("GSE81608_samples.csv")
x[, group := classify(`condition:ch1`, "non-diabetic", "T2D")]
m$GSE81608 <- x
# GSE27951: GPL570 subseries = human adipose; GPL13273 subseries = mouse adipocyte cell culture (not part of this study)
x <- rd("GSE27951_samples_GPL570.csv")
x[, status_from_title := sub(".*_(NGT|IGT|DM) \\(.*", "\\1", title)]
x[, title_status_discordant := status_from_title != `clinical status:ch1`]
x[, group := classify(`clinical status:ch1`, "NGT", "DM")]
x[`clinical status:ch1` == "IGT", group_note := "IGT (impaired glucose tolerance): cannot be classified ND or T2D; group left NA"]
# DECISION (PI/user, 2026-10-03): adipose discovery cohort = GSE27951 GPL570 with IGT samples EXCLUDED
# -> 11 ND vs 12 T2D, n = 23 (plan originally expected 11 ND / 10 T2D, n = 21; observed 11 NGT / 10 IGT / 12 DM = 33).
# Group is defined by the clinical-status field, not the title. Title/field discordant samples
# (GSM691135, GSM691139, GSM691140): rule = clinical-status field; clinical-status field final
# (protocol v1.0, project-lead decision 2026-10-05).
x[, excluded_reason := NA_character_]
x[`clinical status:ch1` == "IGT", excluded_reason := "IGT, user decision 2026-10-03"]
x[title_status_discordant == TRUE, group_note := paste0(ifelse(is.na(group_note), "", paste0(group_note, "; ")),
   "title says ", status_from_title, " but characteristics say ", `clinical status:ch1`, " (characteristics field used)")]
m$GSE27951_GPL570 <- x
# GSE76894 (v1.1 validation; metadata only -- no expression is read here):
# field `diabetes status (nd (non-diabetic), t2d (type 2 diabetic)):ch1`, observed values exactly "ND" / "T2D".
x <- rd("GSE76894_samples.csv")
x[, group := classify(`diabetes status (nd (non-diabetic), t2d (type 2 diabetic)):ch1`, "ND", "T2D")]
m$GSE76894 <- x
x2 <- rd("GSE27951_samples_GPL13273.csv"); x2[, group := NA_character_]
x2[, group_note := "non-human (NMRI mouse) cell-culture subseries on GPL13273; no ND/T2D concept; not classified"]
wr(x2, "GSE27951_samples_GPL13273.csv")

for (nm in names(m)) {
  f <- if (nm == "GSE27951_GPL570") "GSE27951_samples_GPL570.csv" else sprintf("%s_samples.csv", nm)
  wr(m[[nm]], f)
}
unc <- rbindlist(lapply(names(m), function(nm) {
  d <- copy(m[[nm]]); if (!"excluded_reason" %in% names(d)) d[, excluded_reason := NA_character_]
  d[is.na(group), .(dataset = nm, geo_accession, title, excluded_reason)] }))
cat("\nSamples with group = NA (cannot be classified ND/T2D):", nrow(unc), "\n")
if (nrow(unc)) { print(unc[, .N, by = dataset]); wr(unc, "unclassified_samples.csv") }
disc <- m$GSE27951_GPL570[title_status_discordant == TRUE, .(geo_accession, title, characteristics_status = `clinical status:ch1`)]
cat("GSE27951 samples whose title and clinical-status field disagree:", nrow(disc), "\n"); if (nrow(disc)) print(disc)

# logs/excluded_samples.csv (gsm,gse,reason,date): append IGT exclusions idempotently
exf <- here("logs", "excluded_samples.csv")
exn <- m$GSE27951_GPL570[!is.na(excluded_reason), .(gsm = geo_accession, gse = "GSE27951", reason = excluded_reason, date = "2026-10-03")]
old <- if (file.exists(exf)) fread(exf, colClasses = "character") else data.table(gsm = character(), gse = character(), reason = character(), date = character())
fwrite(unique(rbind(old, exn[, lapply(.SD, as.character)])), exf)

# ---- expected counts (from protocol; DO NOT edit to fit the data)
# GSE27951_GPL570 expectation = 11 ND / 12 T2D / n=23 after excluding the 10 IGT (user decision 2026-10-03; was 11/10/21).
expected <- data.table(
  dataset = c("GSE41762", "GSE38642", "GSE20966", "GSE27951_GPL570", "GSE16415", "GSE71416", "GSE81608", "GSE76894"),
  exp_ND = c(57, 54, 10, 11, 5, 6, 651, 84), exp_T2D = c(20, 9, 10, 12, 5, 14, 949, 19),
  exp_n = c(77, 63, 20, 23, 10, 20, 1600, 103),
  exp_platform = c("GPL6244", "GPL6244", "GPL1352", "GPL570", "GPL2986", "GPL570", "GPL16791", "GPL570"))
# GSE76894 expectation from protocol v1.1 section 1 (84 ND / 19 T2D, n=103, GPL570).
# GSE81608 expected platform is not given in the protocol ("scRNA-seq"); observed platform is recorded but not asserted.
chk <- expected[, {
  d <- m[[dataset]]; if ("excluded_reason" %in% names(d)) d <- d[is.na(excluded_reason)]  # excluded samples not counted
  .(obs_ND = sum(d$group == "ND", na.rm = TRUE), obs_T2D = sum(d$group == "T2D", na.rm = TRUE),
    obs_other_or_NA = sum(is.na(d$group)), obs_n = nrow(d),
    obs_platform = paste(sort(unique(d$platform_id)), collapse = "/"))
}, by = dataset]
chk <- merge(expected, chk, by = "dataset", sort = FALSE)
chk[, ok_ND := obs_ND == exp_ND][, ok_T2D := obs_T2D == exp_T2D][, ok_n := obs_n == exp_n]
chk[, ok_platform := ifelse(dataset == "GSE81608", NA, obs_platform == exp_platform)]
chk[, PASS := ok_ND & ok_T2D & ok_n & (is.na(ok_platform) | ok_platform)]
cat("\n---- Count verification (expected vs observed) ----\n"); options(width = 250); print(as.data.frame(chk[, .(dataset, exp_ND, obs_ND, exp_T2D, obs_T2D, exp_n, obs_n, obs_other_or_NA, exp_platform, obs_platform, PASS)]))
wr(chk, "count_verification.csv")
bad <- chk[PASS == FALSE]
if (nrow(bad)) problems <- c(problems, sprintf("COUNT MISMATCH: %s", paste(bad$dataset, collapse = ", ")))
cat("\nGSE27951 full breakdown (GPL570 human subseries), by clinical status field:\n"); print(table(m$GSE27951_GPL570$`clinical status:ch1`))
cat("GSE27951 by status implied in sample title:\n"); print(table(m$GSE27951_GPL570$status_from_title))
cat("GSE27951 other subseries GPL13273 (mouse): n =", nrow(x2), "\n")

# ---- GSE81608: matrix columns versus manifest
mf <- here("data", "raw", "GSE81608", "GSE81608_human_islets_rpkm.txt.gz")
hdr <- strsplit(readLines(gzfile(mf), n = 1), "\t")[[1]]
nrows_mat <- length(count.fields(gzfile(mf), sep = "\t", quote = "", comment.char = "")) - 1
cat(sprintf("\nGSE81608 rpkm matrix: %d sample columns (+ %s), %d gene rows\n", length(hdr) - 1, hdr[1], nrows_mat))
mapped <- paste0("Pancreatic islet cell sample ", sub("^Sample_", "", hdr[-1]))
cat("  matrix columns matching manifest titles:", sum(mapped %in% m$GSE81608$title), "of", length(hdr) - 1, "\n")
if (length(hdr) - 1 != 1600 || sum(mapped %in% m$GSE81608$title) != 1600)
  problems <- c(problems, "GSE81608: rpkm matrix columns do not map 1:1 to the 1600 manifest cells")

# ---- notes
writeLines(c(
 "# Manifest notes (written by scripts/01_acquire.R)", "",
 "- `group` column: derived only from the exact GEO metadata strings listed in scripts/01_acquire.R; NA = could not be classified ND/T2D.",
 "- GSE41762: donors 1-48 and 49-80 were normalised separately by the submitters (per project notes). `submitter_norm_batch` is parsed from the title `ID<n>`; verify against the GEO series description before relying on it.",
 "- GSE38642: retired in protocol v1.1: all 63 CELs byte-identical to GSE41762 arrays (Task 2).",
 "- GSE71416: GEO metadata has no disease characteristic; group is read from sample titles (`OAT_Obese and diabetic` / `OAT_Obese non-diabetic`).",
 "- GSE27951: SuperSeries. GPL570 subseries = 33 human adipose samples with clinical status NGT/IGT/DM (11/10/12); GPL13273 subseries = 12 mouse cell-culture samples (not part of this study). The 10 IGT samples are EXCLUDED (excluded_reason = 'IGT, user decision 2026-10-03'; also in logs/excluded_samples.csv): discovery cohort = 11 ND vs 12 T2D, n=23. Group defined by clinical-status field. Title/field discordant samples GSM691135, GSM691139, GSM691140: rule: clinical-status field; clinical-status field final (protocol v1.0, project-lead decision 2026-10-05). See count_verification.csv.",
 "- GSE81608: single-cell; manifest has 1600 cells; rpkm matrix columns `Sample_<n>` map to title `Pancreatic islet cell sample <n>`.",
 "- GSE16415: raw archive contains per-sample .txt.gz files (not CEL).",
 "- GSE76894 (protocol v1.1, primary islet external validation; IMIDIA organ donors, GPL570): group from `diabetes status (nd (non-diabetic), t2d (type 2 diabetic)):ch1` (ND / T2D). VALIDATION SET: metadata only; expression not read. Only GSE76894_RAW.tar is downloaded; the submitters' differential-expression supplement is intentionally not fetched. Donor independence vs GSE41762/GSE20966: manifest/GSE76894_donor_overlap.csv (scripts/02_donor_overlap.R).") , here("manifest", "NOTES.md"))

# ---- file verification + files.csv
files <- list.files(raw, recursive = TRUE, full.names = TRUE)
chk_file <- function(f) {
  sz <- file.size(f); ok <- sz > 0; msg <- ""
  if (ok && grepl("\\.gz$", f)) {
    r <- tryCatch({ cn <- gzfile(f, "rb"); n <- 0
      repeat { b <- readBin(cn, "raw", 1e7); if (!length(b)) break; n <- n + length(b) }
      close(cn); n }, error = function(e) -1, warning = function(w) -1)
    ok <- r > 0; msg <- if (ok) sprintf("gz ok, %d bytes uncompressed", r) else "gz UNREADABLE"
  } else if (ok && grepl("\\.tar$", f)) {
    l <- tryCatch(untar(f, list = TRUE), error = function(e) character(0))
    ok <- length(l) > 0; msg <- if (ok) sprintf("tar ok, %d members", length(l)) else "tar UNREADABLE"
  } else if (!ok) msg <- "EMPTY"
  data.table(filename = sub(paste0(raw, "/"), "", f, fixed = TRUE), size_bytes = sz,
             md5 = unname(tools::md5sum(f)), opens_nonempty = ok, check = msg)
}
fl <- rbindlist(lapply(files, chk_file))
wr(fl, "files.csv"); cat("\n---- files.csv ----\n"); print(as.data.frame(fl[, .(filename, size_bytes, opens_nonempty, check)]))
if (any(!fl$opens_nonempty)) problems <- c(problems, "Unreadable/empty raw files present")
raw_tar_gse <- c("GSE41762", "GSE38642", "GSE20966", "GSE27951", "GSE16415", "GSE71416", "GSE76894")
need <- c(paste0(raw_tar_gse, "/", raw_tar_gse, "_RAW.tar"),
          "GSE81608/GSE81608_human_islets_rpkm.txt.gz", paste0(platforms, "/", platforms, ".soft.gz"))
miss <- setdiff(need, fl$filename); if (length(miss)) problems <- c(problems, paste("Missing expected files:", paste(miss, collapse = ", ")))
sm <- fl[grepl("series_matrix", filename)]; cat("series matrix files:", nrow(sm), "\n")
if (nrow(sm) != 9) problems <- c(problems, sprintf("expected 9 series_matrix files (8 series; GSE27951 has 2), found %d", nrow(sm)))
# GSE76894 RAW.tar must hold one CEL per manifest sample (names only; CEL contents are not read here)
t76 <- untar(here("data", "raw", "GSE76894", "GSE76894_RAW.tar"), list = TRUE)
cel76 <- sub("_.*", "", grep("\\.CEL(\\.gz)?$", t76, ignore.case = TRUE, value = TRUE))
cat(sprintf("GSE76894 RAW.tar: %d members, %d CEL; CEL GSMs matching manifest: %d / %d\n",
    length(t76), length(cel76), sum(cel76 %in% m$GSE76894$geo_accession), nrow(m$GSE76894)))
if (length(cel76) != 103 || !setequal(cel76, m$GSE76894$geo_accession))
  problems <- c(problems, "GSE76894: RAW.tar CEL files do not map 1:1 to the 103 manifest samples")

# ---- final assertions / verdict
stopifnot(all(file.exists(here("manifest", c("files.csv", "count_verification.csv", "NOTES.md")))))
stopifnot(all(sapply(m, function(d) !anyDuplicated(d$geo_accession))))
stopifnot(all(sapply(m, function(d) all(c("geo_accession", "title", "platform_id", "group") %in% names(d)))))
cat("\n==== Stage 1b verdict\n")
sessionInfo() |> print()
if (length(problems)) {
  cat("FAILED:\n", paste0(" - ", problems, collapse = "\n"), "\n")
  sink(type = "message"); sink(); close(con)
  stop("Stage 1b verification FAILED (see log and manifest/count_verification.csv): ", paste(problems, collapse = " | "), call. = FALSE)
}
cat("Stage 1b PASSED: all counts and files verified.\n")

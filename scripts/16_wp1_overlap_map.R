# 16_wp1_overlap_map.R -- Extension Protocol v2.0, WP1: overlap map of public human islet T2D series.
#
# Series (all Affymetrix, CEL files public): GSE41762, GSE38642, GSE50397, GSE20966, GSE76894,
#   GSE76895, GSE25724.  Series search and include/exclude decisions: results/stage8/01b_series_search_log.csv
#
# HARD GATE (Extension Protocol v2.0): until the protocol is registered on OSF, GSE50397 and GSE76895 may be
#   inspected ONLY for file checksums and sample metadata.  This script therefore
#     * reads series-matrix files up to the line "!series_matrix_table_begin" and stops (no expression table
#       is parsed, loaded or summarised, for any series);
#     * hashes CEL files as raw byte streams (md5 of the decompressed file); CEL contents are never parsed;
#     * never calls rma()/getGEO()/affy readers.
#   The GSE76894 / GSE20966 validation sets are treated the same way (CLAUDE.md rule 2).
#
# Inputs : data/raw/<GSE>/<GSE>_RAW.tar and <GSE>_series_matrix.txt.gz for the 7 series (read-only);
#          manifest/{GSE41762,GSE38642,GSE20966,GSE76894}_samples.csv (existing);
#          manifest/GSE38642_donor_overlap.csv (existing result of 02_donor_overlap.R, consistency check)
# Outputs: manifest/{GSE50397,GSE76895,GSE25724}_samples.csv, manifest/files_extension.csv,
#          results/stage8/16_cel_md5.csv, 16_overlap_pairwise.csv, 16_overlap_donor_level.csv,
#          16_series_provenance.csv, 16_wp1_overlap_map.rds, figures/16_overlap_map.{pdf,png},
#          logs/16_wp1_overlap_map_<timestamp>.log
#
# MATCHING RULES (as in scripts/02_donor_overlap.R)
#   exact : sex exact & diabetes group (ND/T2D) exact & age exact (integer years) & BMI equal to 0.1 kg/m2
#   near  : sex exact & group exact & age +-1 year & BMI +-0.1   (reported separately; never merged with exact)
#   Missing BMI on either side (e.g. one GSE20966 donor, GSE76894 GSM2040122): that pair is compared on
#   sex + group + age only (conservative: yields more possible matches), as in 02_donor_overlap.R sections 2-3.
#   Group NA on either side (GSE50397 has no diabetes status in GEO; IGT and type 3c in GSE76895 are NA, not
#   guessed): the group condition is dropped for that pair (conservative) and the pair is flagged group_rule.
#   Label conflict (primary)   : arrays with byte-identical CEL in two series whose ND/T2D groups are both known
#                                and differ.
#   Label conflict (covariate) : sample pairs identical on sex+age+BMI (any group) whose known groups differ;
#                                reported separately.

set.seed(42)
suppressPackageStartupMessages({ library(here); library(data.table); library(ggplot2); library(patchwork) })

script_name <- "16_wp1_overlap_map"
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create(here("logs"), showWarnings = FALSE)
out_dir <- here("results", "stage8"); dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt")
sink(con, split = TRUE); sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

series <- c("GSE41762", "GSE38642", "GSE50397", "GSE20966", "GSE76894", "GSE76895", "GSE25724")
n_expected <- c(GSE41762 = 77, GSE38642 = 63, GSE50397 = 89, GSE20966 = 20, GSE76894 = 103, GSE76895 = 103, GSE25724 = 13)
new_series <- c("GSE50397", "GSE76895", "GSE25724")
ser_order <- series   # copy used inside data.table expressions (avoids column/variable name clash)
raw <- here("data", "raw")
tar_f <- function(s) file.path(raw, s, sprintf("%s_RAW.tar", s))
sm_f  <- function(s) file.path(raw, s, sprintf("%s_series_matrix.txt.gz", s))

cat("Script:", script_name, " run:", ts, "\n")
cat("Series (order = A before B in pairwise table):", paste(series, collapse = ", "), "\n")
cat("Parameters:\n  exact rule : sex & group & age exact & BMI to 0.1\n  near rule  : age +-1, BMI +-0.1, sex & group exact\n",
    " missing BMI : pair matched on sex+group+age only;  group NA : group condition dropped (flagged)\n",
    " CEL hash    : md5 of decompressed CEL bytes, streamed from the tar (nothing extracted to disk)\n",
    " HARD GATE   : series-matrix files read to '!series_matrix_table_begin' only; no expression value is read\n\n", sep = "")
stopifnot(file.exists(tar_f(series)), file.exists(sm_f(series)))
cat("Input file checksums:\n"); for (f in c(tar_f(series), sm_f(series))) cat(sprintf("  %s  %s\n", unname(tools::md5sum(f)), sub(paste0(raw, "/"), "", f, fixed = TRUE)))

num <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("", "--", "NA", "N/A"), NA, x)))

# ---- 0. series-matrix HEADER parsing (stops at table begin) --------------------------------
read_hdr <- function(f) {
  cn <- gzfile(f, "rt"); on.exit(close(cn)); out <- character()
  repeat { l <- readLines(cn, n = 1); if (!length(l) || grepl("^!series_matrix_table_begin", l)) break; out <- c(out, l) }
  out
}
split_tab <- function(l) { p <- strsplit(paste0(l, "\tEND"), "\t", fixed = TRUE)[[1]]; gsub('^"|"$', "", p[-length(p)]) }
hdr <- setNames(lapply(series, function(s) read_hdr(sm_f(s))), series)
cat("\nHeader lines read (table not read):", paste(sprintf("%s=%d", series, lengths(hdr)), collapse = " "), "\n")

build_manifest <- function(h) {
  sl <- h[startsWith(h, "!Sample_")]
  f <- lapply(sl, function(l) { p <- split_tab(l); list(key = sub("^!Sample_", "", p[1]), val = p[-1]) })
  keys <- vapply(f, `[[`, "", "key")
  get1 <- function(k) f[[which(keys == k)[1]]]$val
  m <- data.table(geo_accession = get1("geo_accession"), title = get1("title"), platform_id = get1("platform_id"))
  chl <- f[keys == "characteristics_ch1"]
  parsed <- list(); rawc <- list()
  for (i in seq_along(chl)) {
    v <- chl[[i]]$val; pos <- regexpr(": ", v, fixed = TRUE)
    key <- v[pos > 0][1]; key <- substr(key, 1, regexpr(": ", key, fixed = TRUE) - 1)
    parsed[[paste0(key, ":ch1")]] <- ifelse(pos > 0, substring(v, pos + 2), "")
    rawc[[if (i == 1) "characteristics_ch1" else paste0("characteristics_ch1.", i - 1)]] <- v
  }
  cbind(m, as.data.table(parsed), as.data.table(rawc))
}

# ---- 0b. manifests for the three new series ----------------------------------------------------
mk <- setNames(lapply(new_series, function(s) build_manifest(hdr[[s]])), new_series)
# GSE50397: GEO metadata has age, BMI, gender, year of birth, HbA1c but NO diabetes status -> group NA (not guessed)
mk$GSE50397[, `:=`(group = NA_character_, group_status_raw = "",
                   group_note = "no diabetes status field in GEO metadata (only age, bmi, gender, year of birth, hba1c); group left NA, not inferred from HbA1c")]
# GSE76895: ND / T2D from the status field; IGT and T3cD not classified (raw string kept)
sc <- grep("^diabetes status", names(mk$GSE76895), value = TRUE); stopifnot(length(sc) == 1)
mk$GSE76895[, group_status_raw := get(sc)]
mk$GSE76895[, group := fcase(group_status_raw == "ND", "ND", group_status_raw == "T2D", "T2D", default = NA_character_)]
mk$GSE76895[, group_note := fifelse(is.na(group), paste0(group_status_raw, ": not ND or T2D; group left NA (excluded from ND-vs-T2D analyses)"), NA_character_)]
mk$GSE76895[, cohort := "partially pancreatectomised patients (PPP), laser-capture microdissected islets"]
# GSE25724
mk$GSE25724[, group_status_raw := `disease state:ch1`]
mk$GSE25724[, group := fcase(group_status_raw == "non-diabetic", "ND", group_status_raw == "type 2 diabetes", "T2D", default = NA_character_)]
mk$GSE25724[, group_note := NA_character_]
dir.create(here("manifest"), showWarnings = FALSE)
for (s in new_series) fwrite(mk[[s]], here("manifest", sprintf("%s_samples.csv", s)))
cat("\nNew manifests: ", paste(sprintf("%s n=%d", new_series, vapply(mk, nrow, 1L)), collapse = "; "), "\n", sep = "")
cat("GSE76895 diabetes-status counts (from GEO sample metadata):\n"); print(mk$GSE76895[, .N, by = group_status_raw][order(group_status_raw)])
cat("GSE76895 platform:", unique(mk$GSE76895$platform_id), "; GSE25724 group counts:\n"); print(mk$GSE25724[, .N, by = .(group_status_raw, group)])
cat("GSE50397 group status: none in GEO; HbA1c present for", sum(mk$GSE50397$`hba1c:ch1` != ""), "of", nrow(mk$GSE50397), "arrays\n")

# ---- 1. file checksums for new raw files ----------------------------------------------------
nf <- unlist(lapply(new_series, function(s) c(tar_f(s), sm_f(s))))
fe <- data.table(filename = sub(paste0(raw, "/"), "", nf, fixed = TRUE), size_bytes = file.size(nf), md5 = unname(tools::md5sum(nf)))
fe[, opens_nonempty := size_bytes > 0]
fe[, check := vapply(nf, function(f) {
  if (grepl("\\.tar$", f)) { l <- untar(f, list = TRUE); sprintf("tar ok, %d members", length(l)) }
  else { gz <- gzfile(f, "rt"); x <- readLines(gz, n = 1); close(gz); if (length(x)) "gz ok (first line read)" else "gz UNREADABLE" } }, "")]
fwrite(fe, here("manifest", "files_extension.csv")); cat("\nmanifest/files_extension.csv:\n"); print(fe[, .(filename, size_bytes, md5, check)])

# ---- 2. CEL md5 (decompressed bytes), streamed from each tar ---------------------------------------
helper <- tempfile(fileext = ".sh")
writeLines(c("#!/bin/sh", "printf '%s\\t%s\\n' \"$(gzip -dc | md5sum | cut -c1-32)\" \"$TAR_FILENAME\""), helper); Sys.chmod(helper, "755")
cel_md5_stream <- function(tarf) {
  o <- system2("tar", c("-xf", shQuote(tarf), "--wildcards", shQuote("*.[Cc][Ee][Ll].gz"), paste0("--to-command=", shQuote(helper))), stdout = TRUE, stderr = FALSE)
  p <- tstrsplit(o, "\t", fixed = TRUE)
  data.table(file = p[[2]], gsm = sub("[_.].*", "", basename(p[[2]])), cel_md5 = p[[1]])
}
h <- rbindlist(lapply(series, function(s) { x <- cel_md5_stream(tar_f(s)); x[, series := s]; x }))
setcolorder(h, c("series", "gsm", "file", "cel_md5")); setorder(h, series, gsm)
cat("\nCEL md5 computed:", paste(sprintf("%s=%d", series, h[, .N, by = series][match(ser_order, series)]$N), collapse = " "), "\n")
fwrite(h, file.path(out_dir, "16_cel_md5.csv"))
dup_within <- h[, .(n_dup_md5 = sum(duplicated(cel_md5))), by = series]
cat("Duplicated CEL md5 within a series:", paste(sprintf("%s=%d", dup_within$series, dup_within$n_dup_md5), collapse = " "), "\n")

# ---- 3. standardised covariates ---------------------------------------------------------------------
rd <- function(f) fread(here("manifest", f), colClasses = "character", na.strings = c("", "NA"))
nordic <- function(s, d) data.table(series = s, gsm = d$geo_accession, title = d$title, group = d$group,
                                    sex = tolower(d[["Sex:ch1"]]), age = num(d[["age:ch1"]]), bmi = round(num(d[["bmi:ch1"]]), 1),
                                    hba1c = round(num(d[["hba1c:ch1"]]), 1))
d41 <- rd("GSE41762_samples.csv"); d38 <- rd("GSE38642_samples.csv"); d20 <- rd("GSE20966_samples.csv"); d94 <- rd("GSE76894_samples.csv")
d50 <- rd("GSE50397_samples.csv"); d95 <- rd("GSE76895_samples.csv"); d25 <- rd("GSE25724_samples.csv")
sexm <- c(f = "female", m = "male")
cov_all <- rbindlist(list(
  nordic("GSE41762", d41), nordic("GSE38642", d38),
  data.table(series = "GSE50397", gsm = d50$geo_accession, title = d50$title, group = d50$group, sex = tolower(d50[["gender:ch1"]]),
             age = num(d50[["age:ch1"]]), bmi = round(num(d50[["bmi:ch1"]]), 1), hba1c = round(num(d50[["hba1c:ch1"]]), 1)),
  data.table(series = "GSE20966", gsm = d20$geo_accession, title = d20$title, group = d20$group, sex = tolower(d20[["gender:ch1"]]),
             age = num(d20[["age (years):ch1"]]), bmi = round(num(d20[["bmi (body mass index, kg/m2):ch1"]]), 1), hba1c = NA_real_),
  data.table(series = "GSE76894", gsm = d94$geo_accession, title = d94$title, group = d94$group, sex = unname(sexm[tolower(d94[["gender (m/f):ch1"]])]),
             age = num(d94[["age (years):ch1"]]), bmi = round(num(d94[["body mass index (bmi):ch1"]]), 1), hba1c = NA_real_),
  data.table(series = "GSE76895", gsm = d95$geo_accession, title = d95$title, group = d95$group, sex = unname(sexm[tolower(d95[["gender (m/f):ch1"]])]),
             age = num(d95[["age (years):ch1"]]), bmi = round(num(d95[["body mass index (bmi):ch1"]]), 1), hba1c = NA_real_),
  data.table(series = "GSE25724", gsm = d25$geo_accession, title = d25$title, group = d25$group, sex = tolower(d25[["gender:ch1"]]),
             age = num(sub(" yrs", "", d25[["age:ch1"]])), bmi = round(num(d25[["bmi (kg/m2):ch1"]]), 1), hba1c = NA_real_)
))
cov_all[, series := factor(series, levels = ser_order)]
cat("\nCovariate completeness (n, ND, T2D, other/NA group; NA sex/age/bmi/hba1c):\n")
print(cov_all[, .(n = .N, ND = sum(group == "ND", na.rm = TRUE), T2D = sum(group == "T2D", na.rm = TRUE), group_NA = sum(is.na(group)),
                  sex_NA = sum(is.na(sex)), age_NA = sum(is.na(age)), bmi_NA = sum(is.na(bmi)), hba1c_NA = sum(is.na(hba1c)),
                  nonint_age = sum(age != round(age), na.rm = TRUE)), by = series])
cat("(GSE50397 group is NA by design; GSE76895 NA = IGT/T3cD; GSE50397 BMI/age decimals as deposited.)\n")

# ---- 4. pairwise comparison ---------------------------------------------------------------------------
cnt <- function(x) length(unique(x))
pair_tab <- list(); link_tab <- list()
pairs <- t(combn(series, 2))
for (k in seq_len(nrow(pairs))) {
  sA <- pairs[k, 1]; sB <- pairs[k, 2]
  A <- cov_all[series == sA]; B <- cov_all[series == sB]
  X <- merge(A[, .(k = 1L, gsm_A = gsm, group_A = group, sex_A = sex, age_A = age, bmi_A = bmi, hba1c_A = hba1c)],
             B[, .(k = 1L, gsm_B = gsm, group_B = group, sex_B = sex, age_B = age, bmi_B = bmi, hba1c_B = hba1c)],
             by = "k", allow.cartesian = TRUE)[, k := NULL]
  X[, `:=`(sex_eq = sex_A == sex_B, age_d = abs(age_A - age_B), bmi_d = abs(bmi_A - bmi_B),
           grp_known = !is.na(group_A) & !is.na(group_B))]
  X[, grp_eq := fifelse(grp_known, group_A == group_B, NA)]
  bmi_ok0 <- function(tol) is.na(X$bmi_d) | X$bmi_d <= tol + 1e-9
  grp_ok <- !X$grp_known | X$grp_eq %in% TRUE
  X[, exact := sex_eq & age_d <= 1e-9 & bmi_ok0(0) & grp_ok]
  X[, exact_ig := sex_eq & age_d <= 1e-9 & bmi_ok0(0)]
  X[, near := sex_eq & age_d <= 1 + 1e-9 & bmi_ok0(0.1) & grp_ok]
  X[, label_conflict_cov := exact_ig & grp_known & !(grp_eq %in% TRUE)]
  X[, bmi_missing_rule := is.na(bmi_d)]
  cel <- merge(h[series == sA, .(gsm_A = gsm, cel_md5)], h[series == sB, .(gsm_B = gsm, cel_md5)], by = "cel_md5")
  X[, identical_CEL := FALSE]; X[cel, identical_CEL := TRUE, on = .(gsm_A, gsm_B)]
  X[, label_conflict_cel := identical_CEL & grp_known & !(grp_eq %in% TRUE)]
  stopifnot(nrow(X) == nrow(A) * nrow(B), sum(X$identical_CEL) == nrow(cel))
  nAB <- c(nrow(A), nrow(B))
  X[, group_rule := fifelse(grp_known, "group required", "group dropped (NA in one series)")]
  # provenance from series headers
  pv <- function(s) { l <- hdr[[s]]; g <- function(k) gsub('"', "", sub("^[^\t]*\t", "", l[startsWith(l, paste0("!Series_", k, "\t"))])); list(contrib = gsub(",+", " ", g("contributor")), pm = g("pubmed_id"), design = paste(c(g("overall_design"), g("summary")), collapse = " "), title = paste(g("title"), collapse = " ")) }
  pa <- pv(sA); pb <- pv(sB)
  sh_contrib <- intersect(pa$contrib, pb$contrib); sh_pm <- intersect(pa$pm, pb$pm)
  nord <- grepl("Nordic Islet Transplantation", pa$design) && grepl("Nordic Islet Transplantation", pb$design)
  imid <- grepl("IMIDIA", paste(pa$title, pa$design)) && grepl("IMIDIA", paste(pb$title, pb$design))
  nB_cel <- cnt(X[identical_CEL == TRUE]$gsm_B); nA_cel <- cnt(X[identical_CEL == TRUE]$gsm_A)
  notes <- c(
    if (nB_cel > 0) sprintf("%d/%d %s arrays byte-identical to %s arrays", nB_cel, nrow(B), sB, sA),
    if (nB_cel != nA_cel) sprintf("identical CEL: %d %s vs %d %s arrays (duplicated CELs inside a series)", nB_cel, sB, nA_cel, sA),
    if (any(X$grp_known == FALSE)) sprintf("group NA in %s: group condition dropped", paste(c(if (anyNA(A$group)) sA, if (anyNA(B$group)) sB), collapse = "/")),
    if (any(X$bmi_missing_rule)) "BMI missing for some samples: sex+group+age only for those pairs",
    if (length(sh_contrib)) sprintf("shared header contributors: %s", paste(sh_contrib, collapse = "; ")),
    if (length(sh_pm)) sprintf("shared PubMed ids: %s", paste(sh_pm, collapse = ",")),
    if (nord) "both state source = Nordic Islet Transplantation Programme (Uppsala)",
    if (imid) "both IMIDIA consortium biobank (SuperSeries GSE76896)")
  pair_tab[[k]] <- data.table(
    series_A = sA, series_B = sB, n_arrays_A = nrow(A), n_arrays_B = nrow(B),
    n_identical_CEL = nB_cel, n_identical_CEL_A_side = nA_cel,
    n_covariate_matched = cnt(X[exact == TRUE]$gsm_B), n_covariate_matched_A_side = cnt(X[exact == TRUE]$gsm_A),
    n_covariate_matched_ignoring_group_B = cnt(X[exact_ig == TRUE]$gsm_B),
    n_near_matched = cnt(X[near == TRUE & exact == FALSE]$gsm_B),
    n_near_matched_A_side = cnt(X[near == TRUE & exact == FALSE]$gsm_A),
    n_label_conflict = cnt(X[label_conflict_cel == TRUE]$gsm_B),
    n_label_conflict_covariate = cnt(X[label_conflict_cov == TRUE]$gsm_B),
    n_CEL_identical_also_covariate_matched = cnt(X[identical_CEL == TRUE & exact == TRUE]$gsm_B),
    n_CEL_identical_sex_discordant = cnt(X[identical_CEL == TRUE & sex_eq == FALSE]$gsm_B),
    n_CEL_identical_hba1c_discordant = cnt(X[identical_CEL == TRUE & !is.na(hba1c_A) & !is.na(hba1c_B) & abs(hba1c_A - hba1c_B) > 1e-9]$gsm_B),
    shared_header_contributors = paste(sh_contrib, collapse = "; "), shared_pubmed_ids = paste(sh_pm, collapse = ","),
    notes = paste(notes, collapse = " | "))
  L <- X[identical_CEL | exact_ig | near]
  if (nrow(L)) {
    L[, `:=`(series_A = sA, series_B = sB)]
    L[, link := fifelse(identical_CEL, "identical_CEL", fifelse(exact, "covariate_exact", fifelse(exact_ig, "covariate_exact_ignoring_group", "covariate_near")))]
    L[, covariate_agreement := paste0("sex", fifelse(sex_eq, "=", "!="), ";age", fifelse(age_d == 0, "=", "!="), ";bmi",
        fifelse(is.na(bmi_d), "NA", fifelse(bmi_d <= 1e-9, "=", "!=")), ";hba1c", fifelse(is.na(hba1c_A) | is.na(hba1c_B), "NA", fifelse(abs(hba1c_A - hba1c_B) < 1e-9, "=", "!=")),
        ";group", fifelse(grp_known, fifelse(grp_eq, "=", "!="), "NA"))]
    link_tab[[k]] <- L[, .(series_A, gsm_A, series_B, gsm_B, link, identical_CEL, covariate_exact = exact, covariate_exact_ignoring_group = exact_ig,
                            covariate_near = near, group_A, group_B, label_conflict_cel, label_conflict_covariate = label_conflict_cov,
                            sex_A, sex_B, age_A, age_B, bmi_A, bmi_B, hba1c_A, hba1c_B, covariate_agreement, group_rule, bmi_missing_rule)]
  }
}
pw <- rbindlist(pair_tab); dl <- rbindlist(link_tab)
dl[, link := factor(link, levels = c("identical_CEL", "covariate_exact", "covariate_exact_ignoring_group", "covariate_near"))]
setorder(dl, series_A, series_B, link, gsm_A, gsm_B); dl[, link := as.character(link)]
fwrite(pw, file.path(out_dir, "16_overlap_pairwise.csv")); fwrite(dl, file.path(out_dir, "16_overlap_donor_level.csv"))
cat("\n==== Pairwise overlap table (B-side counts unless stated) ====\n")
print(pw[, .(series_A, series_B, nA = n_arrays_A, nB = n_arrays_B, ident_CEL = n_identical_CEL, cov_exact = n_covariate_matched, cov_exact_Aside = n_covariate_matched_A_side,
             cov_ign_grp = n_covariate_matched_ignoring_group_B, near_only = n_near_matched, label_conf_CEL = n_label_conflict, label_conf_cov = n_label_conflict_covariate)], nrows = 30)
cat("\nNotes:\n"); for (i in seq_len(nrow(pw))) if (nzchar(pw$notes[i])) cat(sprintf("  %s vs %s: %s\n", pw$series_A[i], pw$series_B[i], pw$notes[i]))

# ---- 5. consistency check against 02_donor_overlap.R result (GSE41762 vs GSE38642) --------------------------
prev <- fread(here("manifest", "GSE38642_donor_overlap.csv"), colClasses = "character", na.strings = c("", "NA"))
cc <- pw[series_A == "GSE41762" & series_B == "GSE38642"]
mine <- dl[series_A == "GSE41762" & series_B == "GSE38642" & identical_CEL == TRUE]
cmp <- merge(prev[, .(gsm, prev_cel = identical_CEL_GSE41762_gsm)], mine[, .(gsm = gsm_B, new_cel = gsm_A)], by = "gsm", all.x = TRUE)
n_prev_shared <- sum(prev$overlap_status == "shared")
cat("\n==== Consistency check: GSE41762 vs GSE38642 (reference: manifest/GSE38642_donor_overlap.csv, 02_donor_overlap.R) ====\n")
cat(sprintf("  identical CEL: this pipeline %d / 63 (reference %d / 63); per-sample CEL partner identical to reference: %d / 63\n",
    cc$n_identical_CEL, sum(!is.na(prev$identical_CEL_GSE41762_gsm)), sum(cmp$prev_cel == cmp$new_cel, na.rm = TRUE)))
cat(sprintf("  label conflicts (identical CEL, group differs): this pipeline %d (reference 11 GROUP_DIFF: %d)\n", cc$n_label_conflict, sum(grepl("GROUP_DIFF", prev$sameID_covariates))))
cat(sprintf("  covariate-exact matched GSE38642 samples (sex+group+age+BMI): this pipeline %d (reference 'shared' %d)\n", cc$n_covariate_matched, n_prev_shared))
cat(sprintf("  sex-label discordance among identical CELs: %d (reference %d); HbA1c discordance: %d (reference %d)\n",
    cc$n_CEL_identical_sex_discordant, sum(grepl("SEX_DIFF", prev$sameID_covariates)), cc$n_CEL_identical_hba1c_discordant, sum(grepl("HBA1C_DIFF", prev$sameID_covariates))))
cons_ok <- cc$n_identical_CEL == 63 && cc$n_label_conflict == 11 && cc$n_covariate_matched == n_prev_shared &&
           all(cmp$prev_cel == cmp$new_cel) && cc$n_CEL_identical_sex_discordant == sum(grepl("SEX_DIFF", prev$sameID_covariates)) &&
           cc$n_CEL_identical_hba1c_discordant == sum(grepl("HBA1C_DIFF", prev$sameID_covariates))
cat("  RESULT:", if (cons_ok) "REPRODUCED (63/63 identical CELs, 11 label conflicts, 52 covariate-matched)" else "*** NOT REPRODUCED -- FLAG ***", "\n")

# ---- 6. provenance table ---------------------------------------------------------------------------------------
prov <- rbindlist(lapply(series, function(s) {
  l <- hdr[[s]]; g <- function(k, sep = " | ") paste(gsub('"', "", sub("^[^\t]*\t", "", l[startsWith(l, paste0("!Series_", k, "\t"))])), collapse = sep)
  data.table(series = s, title = g("title"), status = g("status"), submission_date = g("submission_date"), last_update = g("last_update_date"),
             contact_name = gsub(",+", " ", g("contact_name")), contact_email = g("contact_email"), contact_institute = g("contact_institute"),
             contact_city = g("contact_city"), contact_country = g("contact_country"),
             n_contributors = length(l[startsWith(l, "!Series_contributor\t")]), contributors = gsub(",+", " ", g("contributor", "; ")),
             pubmed_ids = g("pubmed_id", ","), relation = g("relation"),
             overall_design = substr(g("overall_design"), 1, 600))
}))
fwrite(prov, file.path(out_dir, "16_series_provenance.csv"))
cat("\n==== Provenance (series-matrix headers) ====\n")
for (i in seq_len(nrow(prov))) cat(sprintf("  %-9s | %s | contact: %s, %s | contributors: %d | PMIDs: %s | %s\n", prov$series[i], substr(prov$title[i], 1, 60),
    prov$contact_name[i], prov$contact_institute[i], prov$n_contributors[i], prov$pubmed_ids[i], prov$relation[i]))

# ---- 7. Figure 1 --------------------------------------------------------------------------------------------------
lab <- vapply(series, function(s) sprintf("%s (%d)", s, sum(cov_all$series == s)), "")
ix <- setNames(seq_along(series), series)
mat <- function(valfun) {   # lower triangle only (row series later in the list than column series)
  rbindlist(lapply(seq_len(nrow(pw)), function(i) {
    r <- pw[i]; v <- valfun(r)
    data.table(x = r$series_A, y = r$series_B, n = v, frac = v / min(r$n_arrays_A, r$n_arrays_B))
  }))
}
pa <- mat(function(r) r$n_identical_CEL)
# covariate-exact: number of donors (arrays) in the smaller series with >=1 exact match
pb <- mat(function(r) if (r$n_arrays_B <= r$n_arrays_A) r$n_covariate_matched else r$n_covariate_matched_A_side)
dg <- data.table(x = series, y = series, n = as.integer(table(cov_all$series)[series]), frac = NA_real_)
mkplot <- function(d, title, tag) {
  d <- rbind(d, dg)
  d[, kind := fifelse(x == y, "diag", fifelse(n == 0, "zero", "pos"))]
  d[, `:=`(x = factor(x, series, lab), y = factor(y, rev(series), rev(lab)))]
  d[kind != "pos", frac := NA_real_]
  d[, tcol := fifelse(kind == "pos" & frac >= 0.5, "k", fifelse(kind == "pos", "w", "g"))]
  ggplot(d, aes(x, y)) +
    geom_tile(aes(fill = frac), colour = "white", linewidth = 0.7) +
    geom_text(aes(label = n, colour = tcol, fontface = ifelse(kind == "diag", "italic", "bold")), size = 2.3, show.legend = FALSE) +
    scale_colour_manual(values = c(w = "white", k = "black", g = "grey45")) +
    scale_fill_viridis_c(option = "E", limits = c(0, 1), na.value = "grey94", name = "Fraction of the\nsmaller series",
                         breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
    labs(title = title, tag = tag, x = NULL, y = NULL) + coord_fixed() +
    theme_minimal(base_size = 6.5) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(size = 5.3, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 5.3), plot.title = element_text(size = 6.5, face = "bold", hjust = 0),
          legend.key.height = unit(2.2, "mm"), legend.key.width = unit(8, "mm"),
          legend.title = element_text(size = 5.3), legend.text = element_text(size = 5.3), plot.tag = element_text(face = "bold", size = 8.5),
          plot.margin = margin(1, 1, 1, 1))
}
f1 <- (mkplot(pa, "Arrays with byte-identical CEL files", "a") / mkplot(pb, "Donors matching on sex, group, age, BMI", "b")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", plot.title.position = "plot")
w <- 88 / 25.4; hgt <- 140 / 25.4
ggsave(file.path(out_dir, "figures", "16_overlap_map.pdf"), f1, width = w, height = hgt, device = "pdf")
ggsave(file.path(out_dir, "figures", "16_overlap_map.png"), f1, width = w, height = hgt, dpi = 300)
cat("\nFigure 1: results/stage8/figures/16_overlap_map.{pdf,png} (88 mm wide). Lower triangle; cell number = arrays (a) or donors (b) in the smaller series of the pair; colour = that number as a fraction of the smaller series (cividis, colour-blind safe); grey = 0; diagonal (italic) = series size. Panel b: group condition dropped where group NA (GSE50397; IGT/T3cD in GSE76895).\n")

saveRDS(list(pairwise = pw, donor_level = dl, provenance = prov, cel_md5 = h, covariates = cov_all, consistency_ok = cons_ok), file.path(out_dir, "16_wp1_overlap_map.rds"))

# ---- assertions ----------------------------------------------------------------------------------------------------
stopifnot(
  all(as.integer(table(cov_all$series)[series]) == n_expected[series]),
  all(h[, .N, by = series][match(ser_order, series)]$N == n_expected[series]),
  all(vapply(mk, nrow, 1L) == n_expected[new_series]),
  !anyNA(h$cel_md5), all(nchar(h$cel_md5) == 32),
  all(h$cel_md5 != "d41d8cd98f00b204e9800998ecf8427e"),           # md5 of empty stream = failed gunzip
  sum(mk$GSE76895$group_status_raw == "ND") == 32, sum(mk$GSE76895$group_status_raw == "T2D") == 36,
  sum(mk$GSE76895$group_status_raw == "IGT") == 15, sum(mk$GSE76895$group_status_raw == "T3cD") == 20,
  sum(mk$GSE25724$group == "ND") == 7, sum(mk$GSE25724$group == "T2D") == 6,
  nrow(pw) == choose(length(series), 2), !anyNA(pw$n_identical_CEL), !anyNA(pw$n_covariate_matched),
  !anyDuplicated(h[, .(series, gsm)]), !anyNA(cov_all$sex), !anyNA(cov_all$age),
  nrow(fe) == 6, all(fe$size_bytes > 0), !anyNA(fe$md5),
  cons_ok
)
cat("\nAll assertions passed.\n\n")
print(sessionInfo())

# 02_donor_overlap.R -- identify GSE38642 donors shared with GSE41762 (both Nordic Islet
# Transplantation Programme, Uppsala) and check GSE20966 (validation) against both.
#
# Inputs : manifest/GSE41762_samples.csv, manifest/GSE38642_samples.csv,
#          manifest/GSE20966_samples.csv, data/raw/<GSE>/<GSE>_series_matrix.txt.gz (headers only),
#          data/raw/GSE38642,GSE41762/*_RAW.tar (read-only; CEL md5 identity check, extracted to tempdir)
#          manifest/GSE76894_samples.csv, data/raw/GSE76894/GSE76894_series_matrix.txt.gz (headers only),
#          data/raw/GSE76894,GSE20966/*_RAW.tar (read-only; CEL md5 identity check only)
#          GSE76894 is a VALIDATION set (protocol v1.1): only metadata and CEL file hashes are used;
#          no expression value is parsed, normalised or summarised (CLAUDE.md rule 2).
# Outputs: manifest/GSE38642_donor_overlap.csv, manifest/GSE20966_donor_overlap.csv,
#          manifest/GSE76894_donor_overlap.csv,
#          results/stage2/02_donor_overlap.rds, logs/02_donor_overlap_<timestamp>.log
#
# MATCHING RULE (primary, decides overlap_status):
#   sex exact AND diabetes group (ND/T2D) exact AND age exact (integer years, as deposited)
#   AND BMI equal to 0.1 kg/m2 (both series deposit BMI with one decimal).
#   No tolerance is applied: ages are deposited as integers and BMI to one decimal in
#   all three series, so exact equality is the strictest sensible rule. HbA1c is NOT required;
#   where both samples have HbA1c, agreement (equal to 0.1 %) is reported as supporting evidence.
#   Near-match diagnostics (rule without group; age +-1 / BMI +-0.1) are reported only and
#   never change overlap_status.

set.seed(42)
suppressPackageStartupMessages({ library(here); library(data.table) })

script_name <- "02_donor_overlap"
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create(here("logs"), showWarnings = FALSE)
dir.create(here("results", "stage2"), showWarnings = FALSE, recursive = TRUE)
con <- file(here("logs", sprintf("%s_%s.log", script_name, ts)), open = "wt")
sink(con, split = TRUE); sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

in_files <- c(GSE41762 = "manifest/GSE41762_samples.csv",
              GSE38642 = "manifest/GSE38642_samples.csv",
              GSE20966 = "manifest/GSE20966_samples.csv",
              GSE76894 = "manifest/GSE76894_samples.csv")
out_38642 <- here("manifest", "GSE38642_donor_overlap.csv")
out_20966 <- here("manifest", "GSE20966_donor_overlap.csv")
out_76894 <- here("manifest", "GSE76894_donor_overlap.csv")
out_rds   <- here("results", "stage2", "02_donor_overlap.rds")

cat("Script:", script_name, " run:", ts, "\n")
cat("Inputs:\n"); for (f in in_files) cat("  ", f, " md5=", tools::md5sum(here(f)), "\n", sep = "")
cat("\nParameters:\n",
    "  primary rule : sex exact & group exact & age exact (integer) & BMI exact to 0.1\n",
    "  HbA1c        : not required; agreement (to 0.1) reported where both present\n",
    "  diagnostics  : (a) same rule ignoring group; (b) age +-1, BMI +-0.1, sex & group exact\n\n", sep = "")

num <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("", "--", "NA"), NA, x)))
title_id <- function(x) as.integer(sub(".*ID(\\d+).*", "\\1", x))

std_nordic <- function(f) {
  d <- fread(here(f), colClasses = "character")
  data.table(gsm = d$geo_accession, title_id = title_id(d$title), group = d$group,
             sex = tolower(d[["Sex:ch1"]]), age = num(d[["age:ch1"]]),
             bmi = round(num(d[["bmi:ch1"]]), 1), hba1c = round(num(d[["hba1c:ch1"]]), 1))
}
a <- std_nordic(in_files["GSE41762"])   # discovery
b <- std_nordic(in_files["GSE38642"])   # sensitivity
g0 <- fread(here(in_files["GSE20966"]), colClasses = "character")
g <- data.table(gsm = g0$geo_accession, sample_id = g0[["sample id:ch1"]], group = g0$group,
                sex = tolower(g0[["gender:ch1"]]), age = num(g0[["age (years):ch1"]]),
                bmi = round(num(g0[["bmi (body mass index, kg/m2):ch1"]]), 1))

cat("Covariate formats / completeness:\n")
for (nm in c("GSE41762", "GSE38642")) {
  d <- if (nm == "GSE41762") a else b
  cat(sprintf("  %s: n=%d  ND=%d T2D=%d  NA: sex=%d age=%d bmi=%d hba1c=%d  non-integer ages=%d  bmi decimals>1=%d\n",
      nm, nrow(d), sum(d$group == "ND"), sum(d$group == "T2D"), sum(is.na(d$sex)), sum(is.na(d$age)),
      sum(is.na(d$bmi)), sum(is.na(d$hba1c)), sum(d$age != round(d$age), na.rm = TRUE),
      sum(abs(d$bmi * 10 - round(d$bmi * 10)) > 1e-9, na.rm = TRUE)))
}
cat(sprintf("  GSE20966: n=%d  ND=%d T2D=%d  NA: sex=%d age=%d bmi=%d  HbA1c field: %s\n\n",
    nrow(g), sum(g$group == "ND"), sum(g$group == "T2D"), sum(is.na(g$sex)), sum(is.na(g$age)),
    sum(is.na(g$bmi)), if (any(grepl("hba1c", names(g0), ignore.case = TRUE))) "present" else "ABSENT"))

# ---- matcher -----------------------------------------------------------------
match_one <- function(q, ref, use_group = TRUE, age_tol = 0, bmi_tol = 0, use_bmi = TRUE) {
  ok <- ref$sex == q$sex & abs(ref$age - q$age) <= age_tol + 1e-9
  if (use_group) ok <- ok & ref$group == q$group
  if (use_bmi) ok <- ok & abs(ref$bmi - q$bmi) <= bmi_tol + 1e-9
  which(ok %in% TRUE)
}
collapse <- function(x) if (length(x)) paste(x, collapse = ";") else ""

# ---- 1. GSE38642 vs GSE41762 ---------------------------------------------------
res <- rbindlist(lapply(seq_len(nrow(b)), function(i) {
  q <- b[i]
  m  <- match_one(q, a)
  mg <- match_one(q, a, use_group = FALSE)
  mn <- match_one(q, a, age_tol = 1, bmi_tol = 0.1)
  hb <- if (length(m)) ifelse(is.na(q$hba1c) | is.na(a$hba1c[m]), NA, abs(q$hba1c - a$hba1c[m]) < 1e-9) else logical(0)
  data.table(gsm = q$gsm, title_id = q$title_id, group = q$group, sex = q$sex, age = q$age, bmi = q$bmi,
             hba1c = q$hba1c,
             matched_GSE41762_gsm = collapse(a$gsm[m]),
             matched_GSE41762_title_id = collapse(a$title_id[m]),
             matched_GSE41762_hba1c = collapse(ifelse(is.na(a$hba1c[m]), "NA", a$hba1c[m])),
             n_matches = length(m),
             hba1c_agrees = if (!length(m) || all(is.na(hb))) NA else all(hb, na.rm = TRUE),
             hba1c_agrees_per_match = collapse(ifelse(is.na(hb), "NA", hb)),
             same_title_id_in_match = if (length(m)) q$title_id %in% a$title_id[m] else NA,
             diag_matches_ignoring_group = if (length(mg)) collapse(paste0(a$gsm[mg], "(", a$group[mg], ", hba1c=", ifelse(is.na(a$hba1c[mg]), "NA", a$hba1c[mg]), ", ID", a$title_id[mg], ")")) else "",
             diag_near_matches_age1_bmi0.1 = collapse(a$gsm[mn]),
             overlap_status = ifelse(length(m) > 0, "shared", "unique"))
}))

# reverse ambiguity: a GSE41762 donor matched by >1 GSE38642 sample
rev_tab <- res[n_matches > 0, .(gsm38642 = gsm, gsm41762 = strsplit(matched_GSE41762_gsm, ";")[[1]]), by = gsm][, gsm := NULL]
rev_amb <- rev_tab[, .(n = .N, gsm38642 = paste(gsm38642, collapse = ";")), by = gsm41762][n > 1]
fwd_amb <- res[n_matches > 1]
# duplicated covariate profiles within each series (would make matching non-identifiable)
dup_prof <- function(d) d[, .(n = .N, gsms = paste(gsm, collapse = ";")), by = .(sex, group, age, bmi)][n > 1]

cat("==== GSE38642 vs GSE41762 (primary rule) ====\n")
n_sh <- sum(res$overlap_status == "shared"); n_un <- sum(res$overlap_status == "unique")
cat(sprintf("  shared = %d   unique = %d   (total %d)\n", n_sh, n_un, nrow(res)))
cat("  unique subset by group:\n"); print(res[overlap_status == "unique", .N, by = group])
cat("  shared subset by group:\n"); print(res[overlap_status == "shared", .N, by = group])
cat(sprintf("  ambiguous (GSE38642 sample -> >1 GSE41762 donor): %d\n", nrow(fwd_amb)))
if (nrow(fwd_amb)) print(fwd_amb[, .(gsm, group, sex, age, bmi, hba1c, matched_GSE41762_gsm, matched_GSE41762_hba1c, hba1c_agrees_per_match)])
cat(sprintf("  ambiguous (GSE41762 donor <- >1 GSE38642 sample): %d\n", nrow(rev_amb)))
if (nrow(rev_amb)) print(rev_amb)
cat("  duplicated (sex,group,age,bmi) profiles within GSE41762:\n"); print(dup_prof(a))
cat("  duplicated (sex,group,age,bmi) profiles within GSE38642:\n"); print(dup_prof(b))
hb_cmp <- res[n_matches > 0 & !is.na(hba1c_agrees)]
cat(sprintf("  HbA1c: shared samples with HbA1c comparable in both = %d; agree = %d; disagree = %d; not comparable = %d\n",
    nrow(hb_cmp), sum(hb_cmp$hba1c_agrees), sum(!hb_cmp$hba1c_agrees), n_sh - nrow(hb_cmp)))
if (any(!hb_cmp$hba1c_agrees)) print(hb_cmp[hba1c_agrees == FALSE, .(gsm, hba1c, matched_GSE41762_gsm, matched_GSE41762_hba1c)])
tid_mm <- res[n_matches > 0 & !same_title_id_in_match]
cat("  shared samples whose GSE41762 match has a DIFFERENT title ID:\n"); print(tid_mm[, .(gsm, title_id, matched_GSE41762_gsm, matched_GSE41762_title_id, hba1c, matched_GSE41762_hba1c)])
cross <- res[overlap_status == "unique" & diag_matches_ignoring_group != ""]
cat(sprintf("  NOTE: %d of %d 'unique' samples match a GSE41762 donor exactly on sex+age+BMI but with the OPPOSITE group label\n",
    nrow(cross), n_un))
cat(sprintf("  title 'IDn' of GSE38642 sample equals title ID of a matched GSE41762 donor: %d / %d shared\n",
    sum(res$same_title_id_in_match, na.rm = TRUE), n_sh))
cat("  unique samples (full detail incl. diagnostics):\n")
print(res[overlap_status == "unique", .(gsm, title_id, group, sex, age, bmi, hba1c,
                                        diag_matches_ignoring_group, diag_near_matches_age1_bmi0.1)])
exp_ok <- abs(n_sh - 51) <= 3 && abs(n_un - 12) <= 3 && sum(res[overlap_status == "unique"]$group == "T2D") <= 2
cat(sprintf("  vs expectation (~51 shared / ~12 unique / ~1 T2D unique): %s\n\n",
    if (exp_ok) "CONSISTENT" else "*** MATERIALLY DIFFERENT -- FLAG ***"))

# ---- 1b. supplementary identity evidence (does NOT change overlap_status) -------
# (i) same donor title ID ("human islets, IDn") in GSE41762, with covariate comparison;
# (ii) byte-identity of the raw CEL files (md5 of decompressed CEL), read-only from data/raw tars.
idm <- a[match(b$title_id, a$title_id)]
res[, sameID_GSE41762_gsm := idm$gsm]
res[, sameID_GSE41762_group := idm$group]
res[, sameID_covariates := fifelse(is.na(idm$gsm), NA_character_,
        paste0(ifelse(idm$sex == sex, "sex=", "SEX_DIFF("), ifelse(idm$sex == sex, "", paste0(idm$sex, ")")),
               ifelse(idm$age == age, ";age=", ";AGE_DIFF"), ifelse(abs(idm$bmi - bmi) < 1e-9, ";bmi=", ";BMI_DIFF"),
               ifelse(is.na(idm$hba1c) | is.na(hba1c), ";hba1c_NA", ifelse(abs(idm$hba1c - hba1c) < 1e-9, ";hba1c=", ";HBA1C_DIFF")),
               ifelse(idm$group == group, ";group=", ";GROUP_DIFF")))]
cel_md5 <- function(tarf) {
  tmp <- file.path(tempdir(), basename(tarf)); dir.create(tmp, showWarnings = FALSE)
  members <- untar(tarf, list = TRUE); members <- members[grepl("\\.CEL\\.gz$", members, ignore.case = TRUE)]
  untar(tarf, files = members, exdir = tmp)
  h <- vapply(members, function(m) {
    zc <- gzfile(file.path(tmp, m), "rb"); on.exit(close(zc))
    raw <- readBin(zc, "raw", n = 2e8); out <- tempfile(); writeBin(raw, out)
    md <- unname(tools::md5sum(out)); unlink(out); md
  }, character(1))
  unlink(tmp, recursive = TRUE)
  data.table(gsm = sub("[_.].*", "", basename(members)), cel_md5 = h)
}
h38 <- cel_md5(here("data", "raw", "GSE38642", "GSE38642_RAW.tar"))
h41 <- cel_md5(here("data", "raw", "GSE41762", "GSE41762_RAW.tar"))
res[, identical_CEL_GSE41762_gsm := h41$gsm[match(h38$cel_md5[match(gsm, h38$gsm)], h41$cel_md5)]]
res[, identical_CEL_GSE41762_gsm := fifelse(is.na(identical_CEL_GSE41762_gsm), "", identical_CEL_GSE41762_gsm)]
n_cel <- sum(res$identical_CEL_GSE41762_gsm != "")
cat("==== Supplementary identity evidence (does not change overlap_status) ====\n")
cat(sprintf("  GSE38642 samples whose title ID exists in GSE41762: %d / 63\n", sum(!is.na(res$sameID_GSE41762_gsm))))
cat(sprintf("  ... of those, sex+age+BMI identical: %d; group label discordant: %d\n",
    sum(grepl("^sex=;age=;bmi=", res$sameID_covariates)), sum(grepl("GROUP_DIFF", res$sameID_covariates))))
cat(sprintf("  CEL files byte-identical (decompressed md5) to a GSE41762 CEL: %d / 63 (CEL md5 computed: GSE38642 %d, GSE41762 %d)\n",
    n_cel, nrow(h38), nrow(h41)))
cat(sprintf("  identical CEL belongs to the same-title-ID GSE41762 sample: %d / %d\n",
    sum(res$identical_CEL_GSE41762_gsm == res$sameID_GSE41762_gsm, na.rm = TRUE), n_cel))
cat("  samples where same-ID record differs on any covariate or group:\n")
print(res[!grepl("^sex=;age=;bmi=;(hba1c=|hba1c_NA);group=$", sameID_covariates),
          .(gsm, title_id, group, sex, age, bmi, hba1c, sameID_GSE41762_gsm, sameID_GSE41762_group, sameID_covariates,
            identical_CEL_GSE41762_gsm, overlap_status)])
if (n_cel == 63) cat("\n  *** FLAG: every GSE38642 array is a byte-identical copy of a GSE41762 array. The 'unique' samples\n",
                     "  *** under the covariate rule are the SAME arrays with a different diabetes label (or sex label).\n",
                     "  *** GSE38642 contains no donor or array independent of GSE41762.\n\n")

fwrite(res, out_38642)

# ---- 2. GSE20966 vs GSE41762 and GSE38642 -------------------------------------
cat("==== GSE20966 vs GSE41762 / GSE38642 ====\n")
cat("  GSE20966 covariates available: sex, age (integer), BMI (0.1; one sample NA), group.  HbA1c: not deposited.\n")
nordic <- list(GSE41762 = a, GSE38642 = b)
g_res <- rbindlist(lapply(seq_len(nrow(g)), function(i) {
  q <- g[i]
  out <- data.table(gsm = q$gsm, sample_id = q$sample_id, group = q$group, sex = q$sex, age = q$age, bmi = q$bmi)
  for (nm in names(nordic)) {
    ref <- nordic[[nm]]
    has_bmi <- !is.na(q$bmi)
    m  <- match_one(q, ref, use_bmi = has_bmi)
    mn <- match_one(q, ref, age_tol = 1, bmi_tol = 0.1, use_bmi = has_bmi)
    out[, (paste0("matched_", nm, "_gsm")) := collapse(ref$gsm[m])]
    out[, (paste0("n_matches_", nm)) := length(m)]
    out[, (paste0("diag_near_", nm)) := collapse(ref$gsm[mn])]
  }
  out[, match_rule := if (is.na(q$bmi)) "sex+group+age only (BMI missing in GSE20966)" else "sex+group+age+BMI"]
  out
}))
g_res[, overlap_status := ifelse(n_matches_GSE41762 + n_matches_GSE38642 > 0, "possible_overlap", "no_covariate_match")]
fwrite(g_res, out_20966)
print(g_res[, .(gsm, sample_id, group, sex, age, bmi, n_matches_GSE41762, n_matches_GSE38642,
                diag_near_GSE41762, diag_near_GSE38642, overlap_status)])

# provenance from series-matrix headers
prov_keys <- c("!Series_title", "!Series_overall_design", "!Series_contributor", "!Series_contact_institute",
               "!Series_contact_city", "!Series_contact_country", "!Series_pubmed_id")
prov <- rbindlist(lapply(c("GSE41762", "GSE38642", "GSE20966"), function(s) {
  f <- here("data", "raw", s, sprintf("%s_series_matrix.txt.gz", s))
  l <- readLines(gzfile(f), n = 200)
  l <- l[sub("\t.*", "", l) %in% prov_keys]
  data.table(series = s, field = sub("^!Series_", "", sub("\t.*", "", l)), value = gsub('"', "", sub("^[^\t]*\t", "", l)))
}))
prov_sum <- prov[, .(value = paste(value, collapse = " | ")), by = .(series, field)]
cat("\n  Provenance (series-matrix headers):\n")
for (i in seq_len(nrow(prov_sum))) cat(sprintf("   %-8s %-18s %s\n", prov_sum$series[i], prov_sum$field[i], substr(prov_sum$value[i], 1, 300)))
g_overlap <- any(g_res$overlap_status == "possible_overlap")
if (g_overlap) {
  cat("\n  ******************************************************************\n",
      " *** STOP: GSE20966 shows covariate-level overlap with Nordic series ***\n",
      " ******************************************************************\n")
} else cat("\n  GSE20966: no sample matches any GSE41762/GSE38642 donor on the covariates available.\n")


# ---- 3. GSE76894 (protocol v1.1 primary islet validation) vs GSE41762 and GSE20966 ----
# Metadata + CEL file hashes only. Rule as above (sex, group, age exact; BMI equal to 0.1).
# GSE76894 deposits BMI with up to 4 decimals; it is rounded to 0.1 before comparison (the
# reference series deposit one decimal), and the deposited value is kept in `bmi_deposited`.
# Where BMI is missing on either side (one GSE20966 donor; GSE76894 GSM2040122, deposited "NA"),
# sex+group+age is used for that comparison (as in section 2) -- conservative: more possible matches.
# Status: "possible_overlap" if the primary rule matches any GSE41762/GSE20966 donor OR the CEL is
# byte-identical to any GSE41762/GSE20966 CEL; otherwise "independent". Near-match diagnostics are
# reported only. Nothing is removed here (removal before scoring is Task 7, protocol section 7).
cat("\n==== GSE76894 vs GSE41762 / GSE20966 ====\n")
v0 <- fread(here(in_files["GSE76894"]), colClasses = "character", na.strings = c("", "NA"))
sex_map <- c(f = "female", m = "male")
v <- data.table(gsm = v0$geo_accession, title = v0$title, group = v0$group,
                sex = unname(sex_map[tolower(v0[["gender (m/f):ch1"]])]),
                age = num(v0[["age (years):ch1"]]),
                bmi_deposited = num(v0[["body mass index (bmi):ch1"]]),
                submitter_batch = v0[["batch:ch1"]])
v[, bmi := round(bmi_deposited, 1)]
hb_cols <- grep("hba1c|a1c|glycated", names(v0), ignore.case = TRUE, value = TRUE)
cat(sprintf("  GSE76894: n=%d  ND=%d T2D=%d  NA: sex=%d age=%d bmi=%d  non-integer ages=%d  HbA1c field: %s\n",
    nrow(v), sum(v$group == "ND"), sum(v$group == "T2D"), sum(is.na(v$sex)), sum(is.na(v$age)), sum(is.na(v$bmi)),
    sum(v$age != round(v$age), na.rm = TRUE), if (length(hb_cols)) paste(hb_cols, collapse = ",") else "ABSENT"))
cat("  other GSE76894 characteristics fields:", paste(setdiff(grep(":ch1$", names(v0), value = TRUE),
    c("gender (m/f):ch1", "age (years):ch1", "body mass index (bmi):ch1", "batch:ch1",
      "diabetes status (nd (non-diabetic), t2d (type 2 diabetic)):ch1")), collapse = ", "), "\n")
cat("  sex x group:\n"); print(v[, .N, by = .(group, sex)][order(group, sex)])
cat(sprintf("  age range %s-%s; BMI range %.2f-%.2f\n", min(v$age, na.rm = TRUE), max(v$age, na.rm = TRUE),
    min(v$bmi_deposited, na.rm = TRUE), max(v$bmi_deposited, na.rm = TRUE)))

# matcher: if BMI is missing on either side, that pair is matched on sex+group+age only
match_v <- function(q, ref, use_group = TRUE, age_tol = 0, bmi_tol = 0) {
  ok <- ref$sex == q$sex & abs(ref$age - q$age) <= age_tol + 1e-9
  if (use_group) ok <- ok & ref$group == q$group
  ok <- ok & (is.na(ref$bmi) | is.na(q$bmi) | abs(ref$bmi - q$bmi) <= bmi_tol + 1e-9)
  which(ok %in% TRUE)
}
refs <- list(GSE41762 = a, GSE20966 = g)
v_res <- rbindlist(lapply(seq_len(nrow(v)), function(i) {
  q <- v[i]
  out <- q[, .(gsm, title, group, sex, age, bmi_deposited, bmi, submitter_batch)]
  for (nm in names(refs)) {
    ref <- refs[[nm]]
    m  <- match_v(q, ref)
    mg <- match_v(q, ref, use_group = FALSE)
    mn <- match_v(q, ref, age_tol = 1, bmi_tol = 0.1)
    out[, (paste0("matched_", nm, "_gsm")) := collapse(ref$gsm[m])]
    out[, (paste0("n_matches_", nm)) := length(m)]
    out[, (paste0("diag_ignoring_group_", nm)) := collapse(paste0(ref$gsm[mg], "(", ref$group[mg], ")")[seq_along(mg)])]
    out[, (paste0("diag_near_age1_bmi0.1_", nm)) := collapse(ref$gsm[mn])]
  }
  out
}))

# CEL byte identity (decompressed md5) against GSE41762 and GSE20966 CELs
h76 <- cel_md5(here("data", "raw", "GSE76894", "GSE76894_RAW.tar"))
h20 <- cel_md5(here("data", "raw", "GSE20966", "GSE20966_RAW.tar"))
ref_md5 <- rbind(h41[, .(gsm, cel_md5, series = "GSE41762")], h20[, .(gsm, cel_md5, series = "GSE20966")])
v_res[, cel_md5 := h76$cel_md5[match(gsm, h76$gsm)]]
v_res[, identical_CEL_gsm := ref_md5$gsm[match(cel_md5, ref_md5$cel_md5)]]
v_res[, identical_CEL_gsm := fifelse(is.na(identical_CEL_gsm), "", identical_CEL_gsm)]
v_res[, covariates_complete := !is.na(sex) & !is.na(group) & !is.na(age) & !is.na(bmi)]
v_res[, matched_gsm := {
  x <- c(strsplit(matched_GSE41762_gsm, ";")[[1]], strsplit(matched_GSE20966_gsm, ";")[[1]], identical_CEL_gsm)
  collapse(unique(x[x != ""])) }, by = gsm]
v_res[, overlap_status := ifelse(matched_gsm != "", "possible_overlap", "independent")]
v_res[, match_rule := fifelse(is.na(bmi), "sex+group+age exact (BMI NA in GSE76894); or identical CEL md5",
                              "sex+group+age exact & BMI equal to 0.1 (sex+group+age vs reference donors with BMI NA); or identical CEL md5")]
fwrite(v_res, out_76894)

n_po <- sum(v_res$overlap_status == "possible_overlap")
cat(sprintf("  primary-rule matches: vs GSE41762 = %d samples, vs GSE20966 = %d samples\n",
    sum(v_res$n_matches_GSE41762 > 0), sum(v_res$n_matches_GSE20966 > 0)))
cat(sprintf("  CEL md5 computed: GSE76894 %d, GSE41762 %d, GSE20966 %d; GSE76894 CELs byte-identical to a reference CEL: %d\n",
    nrow(h76), nrow(h41), nrow(h20), sum(v_res$identical_CEL_gsm != "")))
cat(sprintf("  duplicated CEL md5 within GSE76894: %d\n", sum(duplicated(h76$cel_md5))))
cat(sprintf("  samples with incomplete covariates (matched on sex+group+age only): %d (%s)\n", sum(!v_res$covariates_complete),
    collapse(v_res[covariates_complete == FALSE]$gsm)))
cat(sprintf("  diagnostics: same sex+age+BMI but opposite group: vs GSE41762 %d, vs GSE20966 %d\n",
    sum(v_res$diag_ignoring_group_GSE41762 != "" & v_res$n_matches_GSE41762 == 0),
    sum(v_res$diag_ignoring_group_GSE20966 != "" & v_res$n_matches_GSE20966 == 0)))
cat(sprintf("  diagnostics: near matches (age +-1, BMI +-0.1, sex & group exact): vs GSE41762 %d, vs GSE20966 %d samples\n",
    sum(v_res$diag_near_age1_bmi0.1_GSE41762 != ""), sum(v_res$diag_near_age1_bmi0.1_GSE20966 != "")))
diag_rows <- v_res[diag_ignoring_group_GSE41762 != "" | diag_ignoring_group_GSE20966 != "" |
                   diag_near_age1_bmi0.1_GSE41762 != "" | diag_near_age1_bmi0.1_GSE20966 != ""]
if (nrow(diag_rows)) {
  cat("  near/diagnostic matches (do NOT change overlap_status), with reference covariates:\n")
  print(diag_rows[, .(gsm, group, sex, age, bmi_deposited, diag_ignoring_group_GSE41762, diag_near_age1_bmi0.1_GSE41762,
                      diag_ignoring_group_GSE20966, diag_near_age1_bmi0.1_GSE20966)])
}
cat(sprintf("  overlap_status: independent = %d, possible_overlap = %d (total %d)\n",
    sum(v_res$overlap_status == "independent"), n_po, nrow(v_res)))
if (n_po) {
  cat("\n  ******************************************************************\n",
      " *** FLAG: GSE76894 samples with possible donor overlap (NOT removed here; Task 7) ***\n",
      " ******************************************************************\n")
  po <- v_res[overlap_status == "possible_overlap"]
  print(po[, .(gsm, title, group, sex, age, bmi_deposited, matched_GSE41762_gsm, matched_GSE20966_gsm, identical_CEL_gsm, matched_gsm)])
  for (i in seq_len(nrow(po))) {
    ids <- strsplit(po$matched_gsm[i], ";")[[1]]
    cat("   ", po$gsm[i], "matched reference records:\n")
    print(rbind(a[gsm %in% ids, .(series = "GSE41762", gsm, group, sex, age, bmi, hba1c)],
                g[gsm %in% ids, .(series = "GSE20966", gsm, group, sex, age, bmi, hba1c = NA_real_)]))
  }
} else cat("  GSE76894: no sample matches any GSE41762/GSE20966 donor; no identical CEL.\n")

# provenance (GSE76894 series-matrix headers; extended key set)
prov_keys76 <- c(prov_keys, "!Series_summary", "!Series_relation", "!Series_contact_name", "!Series_contact_department",
                 "!Sample_source_name_ch1", "!Sample_title", "!Sample_contact_institute")
f76 <- here("data", "raw", "GSE76894", "GSE76894_series_matrix.txt.gz")
l76 <- readLines(gzfile(f76), n = 300)
l76 <- l76[seq_len(grep("^!series_matrix_table_begin", l76)[1] - 1)]   # header block only; expression table discarded
l76 <- l76[sub("\t.*", "", l76) %in% prov_keys76]
prov76 <- data.table(series = "GSE76894", field = sub("^!Sample_", "sample_", sub("^!Series_", "", sub("\t.*", "", l76))),
                     value = vapply(strsplit(sub("^[^\t]*\t", "", l76), "\t"), function(z) {
                       z <- unique(gsub('"', "", z)); if (length(z) > 3) paste0(paste(z[1:3], collapse = " | "), " | ... (", length(z), " distinct)") else paste(z, collapse = " | ") }, character(1)))
prov76 <- prov76[, .(value = paste(value, collapse = " | ")), by = .(series, field)]
prov76[field == "contributor", value := paste0(value, "  [", lengths(strsplit(value, " \\| ")), " contributors]")]
cat("\n  GSE76894 provenance (series-matrix headers):\n")
for (i in seq_len(nrow(prov76))) cat(sprintf("   %-8s %-18s %s\n", prov76$series[i], prov76$field[i], substr(prov76$value[i], 1, 1500)))
contrib76 <- strsplit(prov76[field == "contributor"]$value, " \\| ")[[1]]
contrib41 <- strsplit(prov_sum[series == "GSE41762" & field == "contributor"]$value, " \\| ")[[1]]
contrib20 <- strsplit(prov_sum[series == "GSE20966" & field == "contributor"]$value, " \\| ")[[1]]
c76 <- sub("  \\[.*", "", contrib76)
shared_contrib <- list(GSE41762 = intersect(c76, contrib41), GSE20966 = intersect(c76, contrib20))
for (nm in names(shared_contrib)) cat(sprintf("  contributors shared with %s: %s\n", nm,
    if (length(shared_contrib[[nm]])) paste(shared_contrib[[nm]], collapse = "; ") else "none"))
cat("  GSE20966 contact institute / design:", prov_sum[series == "GSE20966" & field %in% c("contact_institute", "contact_country")]$value, "\n")
cat("  islet-isolation centre / organ-donor source: NOT named in the GSE76894 series-matrix headers\n",
    "  (headers state: IMIDIA consortium/biobank, organ donors, 'cultured islets', contact SIB Lausanne).\n", sep = "")

saveRDS(list(rule = "sex exact & group exact & age exact & BMI exact (0.1)", gse38642 = res,
             reverse_ambiguous = rev_amb, gse20966 = g_res, provenance = prov_sum,
             gse76894 = v_res, provenance_gse76894 = prov76, gse76894_shared_contributors = shared_contrib),
        out_rds)
cat("\nOutputs:\n  ", out_38642, "\n  ", out_20966, "\n  ", out_76894, "\n  ", out_rds, "\n", sep = "")

# ---- assertions ------------------------------------------------------------------
chk <- fread(out_38642); chk2 <- fread(out_20966); chk3 <- fread(out_76894, na.strings = c("", "NA"))
stopifnot(
  nrow(a) == 77, nrow(b) == 63, nrow(g) == 20,
  !anyNA(a[, .(sex, group, age, bmi)]), !anyNA(b[, .(sex, group, age, bmi)]),
  nrow(chk) == 63,
  !anyNA(chk$overlap_status), all(chk$overlap_status %in% c("shared", "unique")),
  sum(chk$overlap_status == "shared") + sum(chk$overlap_status == "unique") == 63,
  all((chk$n_matches > 0) == (chk$overlap_status == "shared")),
  setequal(chk$gsm, b$gsm), !anyDuplicated(chk$gsm),
  nrow(h38) == 63, nrow(h41) == 77,
  nrow(chk2) == 20, setequal(chk2$gsm, g$gsm), !anyNA(chk2$overlap_status),
  nrow(v) == 103, nrow(chk3) == 103, setequal(chk3$gsm, v$gsm), !anyDuplicated(chk3$gsm),
  !anyNA(chk3$overlap_status), all(chk3$overlap_status %in% c("independent", "possible_overlap")),
  all((chk3$overlap_status == "possible_overlap") == (!is.na(chk3$matched_gsm) & nzchar(chk3$matched_gsm))),
  sum(v$group == "ND") == 84, sum(v$group == "T2D") == 19, !anyNA(v$group),
  nrow(h76) == 103, setequal(h76$gsm, v$gsm), nrow(h20) == 20, !anyNA(chk3$cel_md5)
)
cat("\nAll assertions passed.\n\n")
print(sessionInfo())

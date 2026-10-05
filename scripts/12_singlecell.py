#!/usr/bin/env python
"""12_singlecell.py -- GSE81608 single-cell specificity (protocol §5.3, IMPLEMENTATION_PLAN Task 9, single-cell part).

GSE81608 (Xin et al. 2016, human islets, 1,600 cells, 12 ND + 6 T2D donors) is DISCOVERY-ADJACENT DESCRIPTIVE data. It is NOT a
validation set, is not scored by any frozen model, and nothing here feeds back into feature selection, the nested CV or the final panel.

Run from anywhere:  .venv-sc/bin/python scripts/12_singlecell.py        (R helper scripts/12b_singlecell_r.R is called via subprocess;
                                                                          needs R_LIBS_USER=~/R/lib, set automatically if unset)

Steps (protocol §5.3):
  1. QC: 3-MAD outliers on log(genes detected), log(library size), mitochondrial fraction.
  2. Cell types: the authors' labels (`cell subtype`), never relabelled; marker sanity check (INS GCG SST PPY PRSS1 COL1A1).
  3. Detection fraction and mean log2(RPKM+1) per cell type for the 9 genes of interest; dotplots.
  4. PSEUDOBULK per donor x cell type (>= 10 QC-passed cells), limma lmFit/eBayes T2D vs ND per cell type; BH across the 9 genes.
     NO per-cell tests.

*** DEVIATION from protocol §5.3 ("sum counts ... limma-voom or edgeR"): the GEO deposit has RPKM only (no raw counts). Pseudobulk
    = MEAN of log2(RPKM+1) over a donor's cells of one type, tested with limma lmFit/eBayes on that log-scale matrix. Logged to
    logs/deviations.md. ***
"""
import datetime
import gzip
import hashlib
import importlib.metadata as imd
import os
import platform
import random
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ROOT = Path(__file__).resolve().parents[1]      # repo root (no hardcoded absolute paths)
os.chdir(ROOT)
SCRIPT = "12_singlecell"
TS = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
(ROOT / "logs").mkdir(exist_ok=True)
OUT = ROOT / "results" / "stage7"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
LOGF = open(ROOT / "logs" / f"{SCRIPT}_{TS}.log", "w")


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOGF.write(s + "\n")
    LOGF.flush()


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_r(mode):
    env = dict(os.environ)
    env.setdefault("R_LIBS_USER", str(Path.home() / "R" / "lib"))
    r = subprocess.run(["Rscript", "scripts/12b_singlecell_r.R", mode], cwd=ROOT, env=env, capture_output=True, text=True)
    log(f"--- Rscript scripts/12b_singlecell_r.R {mode} (exit {r.returncode}) ---")
    log(r.stdout.rstrip())
    if r.stderr.strip():
        log("[R stderr, last 25 lines]\n" + "\n".join(r.stderr.rstrip().splitlines()[-25:]))
    if r.returncode != 0:
        raise RuntimeError(f"R helper failed in mode {mode}")


# ------------------------------------------------------------------------------------------------ parameters
F_RPKM = "data/raw/GSE81608/GSE81608_human_islets_rpkm.txt.gz"
F_MAN = "manifest/GSE81608_samples.csv"
N_MADS = 3
MIN_CELLS = 10                       # donors with < 10 cells of a type are dropped from that type's pseudobulk (protocol §5.3)
FDR_CUT = 0.05
GENES = {"IL1R1": "nested-CV panel", "SNAI1": "nested-CV panel",
         "SERPINE2": "top islet stability", "NLK": "top islet stability", "CDKN1A": "top islet stability",
         "PCK1": "top islet stability", "FXN": "top islet stability",
         "GPC3": "WGCNA hub", "FGF2": "WGCNA hub"}
MARKERS = {"INS": "beta", "GCG": "alpha", "SST": "delta", "PPY": "PP", "PRSS1": None, "COL1A1": None}
CT_ORDER = ["beta", "alpha", "delta", "PP"]

log(f"Script: {SCRIPT}  run: {TS}")
log("*** GSE81608 = DISCOVERY-ADJACENT DESCRIPTIVE DATA. NOT a validation set. Nothing here informs feature selection, the nested CV,")
log("    the final panel or any validation scoring (CLAUDE.md rules 1-4). ***")
log("*** DEVIATION (protocol §5.3): no raw counts in GEO -> pseudobulk = mean log2(RPKM+1) per donor x cell type, limma lmFit/eBayes ***")
log(f"Seed: {SEED}   python {platform.python_version()}  platform {platform.platform()}")
for pkg in ["scanpy", "anndata", "pandas", "numpy", "scipy", "statsmodels", "matplotlib"]:
    log(f"  package {pkg} {imd.version(pkg)}")
log("Device: CPU (numpy/pandas/scipy). GPU (GTX 1650 4GB) NOT used: rapids-singlecell not installed/attempted (1,600 cells; no benefit). "
    "No stochastic step in this script.")
log("Inputs:")
for f in [F_RPKM, F_MAN, "scripts/12b_singlecell_r.R"]:
    log(f"  {f}  md5={md5(ROOT / f)}")
log(f"Parameters: QC {N_MADS} raw MADs (scipy median_abs_deviation scale=1, scverse convention): log(genes detected) two-sided, "
    f"log(library size) two-sided, mito fraction upper-only; pseudobulk min cells/donor/type = {MIN_CELLS}; FDR cut = {FDR_CUT} "
    f"(BH across the {len(GENES)} genes per cell type and model)")

# ------------------------------------------------------------------------------------------------ load
man = pd.read_csv(ROOT / F_MAN)
expr = pd.read_csv(ROOT / F_RPKM, sep="\t", index_col=0)
log(f"\nRPKM matrix: {expr.shape[0]} genes x {expr.shape[1]} cells; manifest rows: {len(man)}")
assert expr.shape[1] == 1600 and len(man) == 1600
assert expr.index.is_unique and not expr.isna().any().any() and (expr.values >= 0).all()
man["col"] = man["title"].str.replace("Pancreatic islet cell sample ", "Sample_", regex=False)
assert man["col"].is_unique and set(man["col"]) == set(expr.columns), "manifest titles do not map 1:1 to matrix columns"
man = man.set_index("col").loc[expr.columns].copy()
obs = pd.DataFrame({"gsm": man["geo_accession"], "celltype": man["cell subtype:ch1"], "donor": man["donor id:ch1"], "group": man["group"],
                    "age": man["age:ch1"].astype(float), "sex": man["gender:ch1"], "ethnicity": man["ethnicity:ch1"]}, index=expr.columns)
assert set(obs["group"]) == {"ND", "T2D"} and set(obs["celltype"]) == {"alpha", "beta", "delta", "PP"}
dchk = obs.groupby("donor")[["group", "age", "sex"]].nunique()
assert (dchk == 1).all().all(), "donor attributes inconsistent across cells"

# Entrez -> symbol (rows are Entrez IDs)
run_r("map")
emap = pd.read_csv(OUT / "12_entrez_symbol_map.csv")
assert emap["entrez"].is_unique
sym = pd.Series(expr.index.map(dict(zip(emap["entrez"], emap["symbol"]))), index=expr.index)
log(f"Entrez IDs mapped to a symbol: {sym.notna().sum()} / {len(sym)} (unmapped rows kept for QC, dropped from symbol-level analyses)")
mt_mask = sym.fillna("").str.startswith("MT-").values
log(f"Mitochondrial genes (symbol starts 'MT-') present in matrix: {int(mt_mask.sum())}")
assert mt_mask.sum() >= 10, "too few MT genes to compute mito fraction"
X = expr.values.astype(np.float64)

# ------------------------------------------------------------------------------------------------ 1. QC
log("\n=== 1. QC ===")
log("NOTE: values are RPKM, not counts. 'library size' = total RPKM per cell as given; 'mito fraction' = sum RPKM of MT- genes / total RPKM "
    "(RPKM-based, so it is a fraction of length-normalised signal, not of reads); 'genes detected' = genes with RPKM > 0.")
n_genes = (X > 0).sum(0)
lib = X.sum(0)
assert (lib > 0).all()
mito = X[mt_mask].sum(0) / lib
qc = pd.DataFrame({"n_genes": n_genes, "total_rpkm": lib, "mito_frac": mito,
                   "log_genes": np.log(n_genes), "log_lib": np.log(lib)}, index=expr.columns)


def mad_flag(v, upper_only=False, scale=1.0):
    med = np.median(v)
    mad = stats.median_abs_deviation(v, scale=scale)
    hi = v > med + N_MADS * mad
    lo = v < med - N_MADS * mad
    return (hi if upper_only else (hi | lo)), (med, mad)


f_g, p_g = mad_flag(qc["log_genes"].values)
f_l, p_l = mad_flag(qc["log_lib"].values)
f_m, p_m = mad_flag(qc["mito_frac"].values, upper_only=True)
qc["out_genes"], qc["out_lib"], qc["out_mito"] = f_g, f_l, f_m
qc["qc_fail"] = f_g | f_l | f_m
qc = pd.concat([obs[["gsm", "celltype", "donor", "group"]], qc], axis=1)
log(f"log genes detected: median {p_g[0]:.3f}, MAD {p_g[1]:.3f} -> fails {f_g.sum()}")
log(f"log library size  : median {p_l[0]:.3f}, MAD {p_l[1]:.3f} -> fails {f_l.sum()}")
log(f"mito fraction     : median {p_m[0]:.5f}, MAD {p_m[1]:.5f} -> fails {f_m.sum()} (upper only)")
log(f"cells failing QC (union): {qc['qc_fail'].sum()} / {len(qc)}  ({100 * qc['qc_fail'].mean():.1f}%); kept {int((~qc['qc_fail']).sum())}")
sc_g, _ = mad_flag(qc["log_genes"].values, scale="normal")
sc_l, _ = mad_flag(qc["log_lib"].values, scale="normal")
sc_m, _ = mad_flag(qc["mito_frac"].values, upper_only=True, scale="normal")
log(f"(information only) with normal-scaled MAD (x1.4826) the union would fail {int((sc_g | sc_l | sc_m).sum())} cells")
qc.to_csv(OUT / "12_qc_cells.csv", index=False)

by_donor = qc.groupby(["group", "donor"]).agg(n_cells=("qc_fail", "size"), n_removed=("qc_fail", "sum")).reset_index()
by_donor["n_kept"] = by_donor["n_cells"] - by_donor["n_removed"]
by_donor["pct_removed"] = (100 * by_donor["n_removed"] / by_donor["n_cells"]).round(1)
by_donor.to_csv(OUT / "12_qc_by_donor.csv", index=False)
by_group = qc.groupby("group").agg(n_cells=("qc_fail", "size"), n_removed=("qc_fail", "sum")).reset_index()
by_group["n_kept"] = by_group["n_cells"] - by_group["n_removed"]
by_group["pct_removed"] = (100 * by_group["n_removed"] / by_group["n_cells"]).round(1)
by_ct = qc.groupby("celltype").agg(n_cells=("qc_fail", "size"), n_removed=("qc_fail", "sum")).reset_index()
by_ct["n_kept"] = by_ct["n_cells"] - by_ct["n_removed"]
by_ct["pct_removed"] = (100 * by_ct["n_removed"] / by_ct["n_cells"]).round(1)
by_group.to_csv(OUT / "12_qc_by_group.csv", index=False)
by_ct.to_csv(OUT / "12_qc_by_celltype.csv", index=False)
log("\nCells removed per group:\n" + by_group.to_string(index=False))
log("\nCells removed per cell type:\n" + by_ct.to_string(index=False))
log("\nCells removed per donor:\n" + by_donor.to_string(index=False))

keep = ~qc["qc_fail"].values
obs_k = obs.loc[keep].copy()
X_k = X[:, keep]
logX = np.log2(X_k + 1.0)

# ------------------------------------------------------------------------------------------------ 2. cell types
log("\n=== 2. Cell types (authors' labels; not relabelled) ===")
ct_tab = pd.crosstab([obs_k["group"], obs_k["donor"]], obs_k["celltype"]).reindex(columns=CT_ORDER, fill_value=0)
ct_tab.to_csv(OUT / "12_celltype_by_donor_postQC.csv")
log("Post-QC cells per donor x cell type:\n" + ct_tab.to_string())
ct_tot = pd.crosstab(obs_k["group"], obs_k["celltype"]).reindex(columns=CT_ORDER, fill_value=0)
ct_tot.to_csv(OUT / "12_celltype_by_condition_postQC.csv")
log("\nPost-QC cells per condition x cell type:\n" + ct_tot.to_string())

# symbol-level matrix (drop unmapped; resolve duplicate symbols by keeping the highest-mean Entrez row)
sym_k = sym.values
has_sym = pd.notna(sym_k)
row_mean = logX.mean(1)
df_s = pd.DataFrame({"symbol": sym_k, "mean": row_mean, "i": np.arange(len(sym_k))})[has_sym]
n_dup = int(df_s["symbol"].duplicated(keep=False).sum())
df_s = df_s.sort_values("mean", ascending=False).drop_duplicates("symbol")
log(f"\nSymbol-level matrix: {len(df_s)} unique symbols ({n_dup} rows shared a symbol; kept highest-mean row per symbol)")
idx_by_symbol = dict(zip(df_s["symbol"], df_s["i"]))
sym_index = sorted(idx_by_symbol)
L = logX[[idx_by_symbol[s] for s in sym_index], :]          # symbols x cells (log2(RPKM+1))
sym_pos = {s: i for i, s in enumerate(sym_index)}
missing = [g for g in list(GENES) + list(MARKERS) if g not in sym_pos]
log(f"Genes of interest / markers absent from matrix: {missing if missing else 'none'}")
assert not missing

log("\nMarker sanity check (mean log2(RPKM+1) per authors' label, QC-passed cells):")
rows = []
for g, exp_ct in MARKERS.items():
    r = {"gene": g, "expected_label": exp_ct if exp_ct else "none present (acinar/stellate labels absent)"}
    for ct in CT_ORDER:
        m = (obs_k["celltype"] == ct).values
        r[f"mean_{ct}"] = L[sym_pos[g], m].mean()
        r[f"det_{ct}"] = (L[sym_pos[g], m] > 0).mean()
    r["highest_label"] = max(CT_ORDER, key=lambda c: r[f"mean_{c}"])
    r["as_expected"] = (r["highest_label"] == exp_ct) if exp_ct else "n/a"
    rows.append(r)
mk = pd.DataFrame(rows)
mk.to_csv(OUT / "12_marker_check.csv", index=False)
log(mk.round(3).to_string(index=False))
for g in ["INS", "GCG", "SST", "PPY"]:
    assert bool(mk.loc[mk["gene"] == g, "as_expected"].iloc[0]), f"marker {g} not highest in its expected label"
log("INS/GCG/SST/PPY highest in beta/alpha/delta/PP respectively: CONFIRMED. PRSS1 and COL1A1: the deposit has only the four endocrine "
    "labels (no acinar/stellate), so there is no expected label; values above show where this ambient/contaminant signal sits.")

# ------------------------------------------------------------------------------------------------ 3. gene expression per cell type
log("\n=== 3. Detection fraction and mean log2(RPKM+1) per cell type ===")
rows = []
for g, cat in GENES.items():
    for ct in CT_ORDER:
        for grp in ["all", "ND", "T2D"]:
            m = (obs_k["celltype"] == ct).values
            if grp != "all":
                m &= (obs_k["group"] == grp).values
            v = L[sym_pos[g], m]
            rows.append({"gene": g, "category": cat, "celltype": ct, "group": grp, "n_cells": int(m.sum()),
                         "detection_fraction": float((v > 0).mean()), "mean_log2_rpkm1": float(v.mean())})
expr_tab = pd.DataFrame(rows)
expr_tab.to_csv(OUT / "12_gene_expression_by_celltype.csv", index=False)
wide = expr_tab[expr_tab.group == "all"].pivot(index="gene", columns="celltype", values=["detection_fraction", "mean_log2_rpkm1"])
wide = wide.reindex(list(GENES))
log(wide.round(3).to_string())


def dotplot(df, rows_order, row_label, fname, title):
    genes = list(GENES)
    fig, ax = plt.subplots(figsize=(7.2, 0.45 * len(rows_order) + 2.2))
    vmax = max(df["mean_log2_rpkm1"].max(), 1e-6)
    for yi, rk in enumerate(rows_order):
        for xi, g in enumerate(genes):
            r = df[(df.gene == g) & (df.row == rk)].iloc[0]
            ax.scatter(xi, yi, s=20 + 380 * r["detection_fraction"], c=[r["mean_log2_rpkm1"]], cmap="viridis", vmin=0, vmax=vmax,
                       edgecolor="0.3", linewidth=0.4)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_yticks(range(len(rows_order)))
    ax.set_yticklabels([row_label[r] for r in rows_order])
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(-0.6, len(rows_order) - 0.4)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=10)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, vmax))
    cb = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("mean log2(RPKM+1)")
    for fr in (0.1, 0.5, 1.0):
        ax.scatter([], [], s=20 + 380 * fr, c="0.75", edgecolor="0.3", linewidth=0.4, label=f"{int(fr * 100)}% cells detected")
    ax.legend(loc="upper left", bbox_to_anchor=(1.18, 1), frameon=False, fontsize=7, title="dot size", title_fontsize=7)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{fname}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


d1 = expr_tab[expr_tab.group == "all"].copy()
d1["row"] = d1["celltype"]
lab1 = {c: f"{c} (n={int((obs_k['celltype'] == c).sum())})" for c in CT_ORDER}
dotplot(d1, CT_ORDER, lab1, "12_dotplot_genes_by_celltype", "GSE81608 (descriptive): detection fraction & mean log2(RPKM+1), QC-passed cells")
d2 = expr_tab[expr_tab.group != "all"].copy()
d2["row"] = d2["celltype"] + "|" + d2["group"]
rows2 = [f"{c}|{g}" for c in CT_ORDER for g in ("ND", "T2D")]
lab2 = {f"{c}|{g}": f"{c} {g} (n={int(((obs_k['celltype'] == c) & (obs_k['group'] == g)).sum())})" for c in CT_ORDER for g in ("ND", "T2D")}
dotplot(d2, rows2, lab2, "12_dotplot_genes_by_celltype_condition",
        "GSE81608 (descriptive, cell-level, no tests): by cell type and condition")

# ------------------------------------------------------------------------------------------------ 4. pseudobulk
log("\n=== 4. Pseudobulk (mean log2(RPKM+1) per donor x cell type; donors with >= %d QC-passed cells of that type) ===" % MIN_CELLS)
meta_rows = []
pb = {}
for ct in CT_ORDER:
    cols, names = [], []
    for donor, d in obs_k[obs_k["celltype"] == ct].groupby("donor"):
        if len(d) < MIN_CELLS:
            continue
        mi = np.where((obs_k["celltype"] == ct).values & (obs_k["donor"] == donor).values)[0]
        cols.append(L[:, mi].mean(1))
        slug = donor.replace(" ", "")
        sid = f"{ct}__{slug}"
        names.append(sid)
        meta_rows.append({"sample": sid, "celltype": ct, "donor": donor, "group": d["group"].iloc[0], "age": d["age"].iloc[0],
                          "sex": d["sex"].iloc[0], "n_cells": len(d)})
    if names:
        pb[ct] = pd.DataFrame(np.column_stack(cols), index=sym_index, columns=names)
        pb[ct].index.name = "gene"
        pb[ct].to_csv(OUT / f"12_pseudobulk_{ct}.csv", float_format="%.6g")
    else:
        stale = OUT / f"12_pseudobulk_{ct}.csv"
        if stale.exists():
            stale.unlink()
meta = pd.DataFrame(meta_rows)
meta.to_csv(OUT / "12_pseudobulk_meta.csv", index=False)
dtab = meta.pivot_table(index="celltype", columns="group", values="donor", aggfunc="count", fill_value=0).reindex(CT_ORDER).fillna(0).astype(int)
for c in ("ND", "T2D"):
    if c not in dtab:
        dtab[c] = 0
dtab = dtab[["ND", "T2D"]]
dtab["n_donors"] = dtab.sum(1)
dtab["testable"] = (dtab[["ND", "T2D"]].min(1) >= 2)
n_donor_total = obs.groupby("group")["donor"].nunique()
dtab.to_csv(OUT / "12_donors_per_celltype.csv")
log(f"Donors in study: {n_donor_total.to_dict()}")
log("Donors (>= %d cells) per group per cell type:\n%s" % (MIN_CELLS, dtab.to_string()))
log("\nPseudobulk sample metadata:\n" + meta.to_string(index=False))
log("\nPOWER LIMITATION: only 6 T2D donors exist in the whole dataset (<= 6 per cell type, 10 ND in alpha/beta). Delta cells: no donor has "
    f">= {MIN_CELLS} QC-passed delta cells, PP cells: too few donors per group -> those cell types cannot be tested without relaxing "
    f"the protocol threshold, which is NOT done. Absence of significance is 'inconclusive', not 'no effect'.")
assert (meta["n_cells"] >= MIN_CELLS).all()
for ct, m in pb.items():
    assert not m.isna().any().any()

# ------------------------------------------------------------------------------------------------ limma (R)
log("\n=== 4b. limma lmFit/eBayes on pseudobulk (R) ===")
run_r("limma")
res = pd.read_csv(OUT / "12_limma_allgenes.csv.gz")
info = pd.read_csv(OUT / "12_limma_model_info.csv")
info.to_csv(OUT / "12_limma_model_info.csv", index=False)
log("\nModel info:\n" + info.to_string(index=False))

# overall check
log("\nOverall check (all fitted genes; BH over all genes within cell type x model):")
ov = []
for (ct, mn), d in res.groupby(["celltype", "model"]):
    ov.append({"celltype": ct, "model": mn, "n_genes": len(d), "n_p_lt_0.05": int((d["P.Value"] < 0.05).sum()),
               "frac_p_lt_0.05": round(float((d["P.Value"] < 0.05).mean()), 4),
               "n_FDR_lt_0.05_allgenes": int((d["adj.P.Val_allgenes"] < 0.05).sum()), "median_abs_t": round(float(d["t"].abs().median()), 3)})
ov = pd.DataFrame(ov)
ov.to_csv(OUT / "12_limma_overall_check.csv", index=False)
log(ov.to_string(index=False))
log("(Expected under no effect: ~5% of genes at p < 0.05. With 6 vs 10 donors few or no genes pass genome-wide FDR; this is a calibration check, "
    "not a finding.)")

# ------------------------------------------------------------------------------------------------ panel results
rows = []
for (ct, mn), d in res.groupby(["celltype", "model"]):
    inf = info[info.celltype == ct].iloc[0]
    primary = (mn == "adjusted") or (mn == "unadjusted" and not bool(inf["adjusted_feasible"]))
    sub = d[d["gene"].isin(GENES)].set_index("gene").reindex(list(GENES))
    ok = sub["P.Value"].notna()
    fdr = pd.Series(np.nan, index=sub.index)
    if ok.sum():
        fdr[ok] = multipletests(sub.loc[ok, "P.Value"], method="fdr_bh")[1]
    df_res = inf["n_donors"] - (4 if mn == "adjusted" else 2)
    n1, n2 = int(inf["n_ND"]), int(inf["n_T2D"])
    tq = stats.t.ppf(0.975, df_res) + stats.t.ppf(0.80, df_res)
    for g in GENES:
        r = sub.loc[g]
        rows.append({"celltype": ct, "model": mn, "primary": primary, "gene": g, "category": GENES[g], "n_ND": n1, "n_T2D": n2,
                     "logFC_T2D_vs_ND": r["logFC"], "AveExpr": r["AveExpr"], "t": r["t"], "P.Value": r["P.Value"],
                     "FDR_BH_9genes": fdr[g], "FDR_BH_allgenes": r["adj.P.Val_allgenes"],
                     "approx_min_detectable_abs_logFC_80pct_power": tq * r["sigma_post"] * np.sqrt(1 / n1 + 1 / n2) if ok[g] else np.nan,
                     "note": "" if ok[g] else "gene not fitted (pseudobulk > 0 in < 3 donors)"})
panel = pd.DataFrame(rows)
panel.to_csv(OUT / "12_panel_pseudobulk_results.csv", index=False)
log("\n=== 4c. T2D vs ND pseudobulk results for the 9 genes (logFC = T2D - ND, log2 scale of RPKM+1; BH across the 9 genes) ===")
show = panel[panel.primary].copy()
log("PRIMARY model (adjusted for age + sex where feasible, else unadjusted):")
log(show[["celltype", "model", "gene", "n_ND", "n_T2D", "logFC_T2D_vs_ND", "P.Value", "FDR_BH_9genes", "FDR_BH_allgenes",
          "approx_min_detectable_abs_logFC_80pct_power", "note"]].round(4).to_string(index=False))
log("\nSENSITIVITY (unadjusted model where an adjusted model was primary):")
sens = panel[(~panel.primary)]
log(sens[["celltype", "model", "gene", "logFC_T2D_vs_ND", "P.Value", "FDR_BH_9genes"]].round(4).to_string(index=False))

# ------------------------------------------------------------------------------------------------ restriction verdicts
log("\n=== 4d. Beta-cell restriction verdict (protocol §5.3) ===")
log("Rule: significant (FDR < 0.05) in beta pseudobulk AND not significant in alpha or delta. Absence of significance = inconclusive.")


def sig(ct, g, model="primary"):
    d = panel[(panel.celltype == ct) & (panel.gene == g)]
    d = d[d.primary] if model == "primary" else d[d.model == model]
    if d.empty or d["FDR_BH_9genes"].isna().all():
        return None                      # not testable
    return bool(d["FDR_BH_9genes"].iloc[0] < FDR_CUT)


verd = []
for g in GENES:
    b, a, dl = sig("beta", g), sig("alpha", g), sig("delta", g)
    st = lambda x: "not testable" if x is None else ("significant" if x else "not significant")
    if b is None:
        v = "NOT ASSESSABLE (beta not testable for this gene)"
    elif b:
        others = [n for n, x in (("alpha", a), ("delta", dl)) if x]
        untest = [n for n, x in (("alpha", a), ("delta", dl)) if x is None]
        if others:
            v = "significant in beta but NOT beta-restricted (also significant in " + ", ".join(others) + ")"
        elif untest:
            v = ("significant in beta; not significant in alpha; " + ", ".join(untest) +
                 " NOT TESTABLE -> restriction criterion cannot be fully assessed (not confirmed)")
        else:
            v = "BETA-RESTRICTED (significant in beta, not in alpha/delta; low-powered comparison)"
    else:
        v = "INCONCLUSIVE (not significant in beta pseudobulk; low power, not evidence of no effect)"
    pb_r = panel[(panel.celltype == "beta") & (panel.gene == g) & panel.primary].iloc[0]
    verd.append({"gene": g, "category": GENES[g], "beta_logFC": pb_r["logFC_T2D_vs_ND"], "beta_P": pb_r["P.Value"],
                 "beta_FDR_9genes": pb_r["FDR_BH_9genes"], "beta": st(b), "alpha": st(a), "delta": st(dl), "verdict": v})
verd = pd.DataFrame(verd)
verd.to_csv(OUT / "12_beta_restriction_verdicts.csv", index=False)
log(verd.round(4).to_string(index=False))

# ------------------------------------------------------------------------------------------------ pseudobulk figure
fig, axes = plt.subplots(3, 3, figsize=(10, 9))
rng = np.random.default_rng(SEED)
xs = {("alpha", "ND"): 0, ("alpha", "T2D"): 1, ("beta", "ND"): 2.5, ("beta", "T2D"): 3.5}
for ax, g in zip(axes.ravel(), GENES):
    for (ct, grp), x in xs.items():
        if ct not in pb:
            continue
        ids = meta[(meta.celltype == ct) & (meta.group == grp)]["sample"]
        v = pb[ct].loc[g, ids].values
        ax.scatter(x + rng.uniform(-0.15, 0.15, len(v)), v, s=22, c="tab:blue" if grp == "ND" else "tab:red", edgecolor="k", linewidth=0.3, zorder=3)
        if len(v):
            ax.hlines(v.mean(), x - 0.28, x + 0.28, color="k", linewidth=1.2, zorder=4)
    pr = panel[(panel.gene == g) & panel.primary]
    txt = []
    for ct in ("alpha", "beta"):
        r = pr[pr.celltype == ct]
        if len(r) and pd.notna(r["P.Value"].iloc[0]):
            txt.append(f"{ct}: logFC {r['logFC_T2D_vs_ND'].iloc[0]:+.2f}, FDR {r['FDR_BH_9genes'].iloc[0]:.2f}")
    ax.set_title(f"{g}\n" + " | ".join(txt), fontsize=8)
    ax.set_xticks(list(xs.values()))
    ax.set_xticklabels(["alpha\nND", "alpha\nT2D", "beta\nND", "beta\nT2D"], fontsize=7)
    ax.set_ylabel("pseudobulk mean log2(RPKM+1)", fontsize=7)
    ax.tick_params(labelsize=7)
fig.suptitle("GSE81608 pseudobulk per donor (descriptive; blue = ND, red = T2D; model-adjusted FDR across 9 genes in titles).\n"
             "Delta/PP not testable (< 10 cells/donor). 6 T2D vs 10 ND donors: low power.", fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.94))
for ext in ("pdf", "png"):
    fig.savefig(FIG / f"12_pseudobulk_beta_alpha_genes.{ext}", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------------------------------------ README + deviations
(OUT / "12_README.txt").write_text(
    "results/stage7 -- GSE81608 single-cell (scripts/12_singlecell.py, scripts/12b_singlecell_r.R)\n"
    "GSE81608 is DISCOVERY-ADJACENT DESCRIPTIVE data, not validation. It informs no feature selection, model or panel decision.\n"
    "Input is RPKM (no raw counts available): pseudobulk = mean log2(RPKM+1) per donor x cell type, limma lmFit/eBayes (deviation from\n"
    "protocol §5.3 'sum counts'; see logs/deviations.md). Cell types = authors' labels; delta and PP have no testable donors at >= 10 cells.\n"
    "BH is across the 9 genes of interest within each cell type/model (also provided: BH across all fitted genes). No per-cell tests.\n")

dev_path = ROOT / "logs" / "deviations.md"
dev_txt = dev_path.read_text()
today = datetime.date.today().isoformat()
new_rows = []
if "scripts/12_singlecell.py (RPKM" not in dev_txt:
    new_rows.append(
        f"| {today} | GSE81608 pseudobulk input (Task 9, scripts/12_singlecell.py (RPKM)) | protocol §5.3: pseudobulk = sum of raw counts per donor x cell type; limma-voom or edgeR | "
        "Pseudobulk = mean of log2(RPKM+1) over a donor's cells of one cell type; limma lmFit + eBayes on that log-scale matrix (no voom: no counts) | "
        "GEO deposit GSE81608 contains only the RPKM matrix, no raw counts. RPKM is length- and depth-normalised, so summing is not valid for voom/edgeR; "
        "means of log values are the closest available analogue. Descriptive/discovery-adjacent analysis only; results are less well calibrated than count-based pseudobulk. |")
if "scripts/12_singlecell.py (delta/PP)" not in dev_txt:
    new_rows.append(
        f"| {today} | GSE81608 delta and PP cell types (Task 9, scripts/12_singlecell.py (delta/PP)) | protocol §5.3: restriction = significant in beta and NOT significant in alpha or delta; donors with >= 10 cells kept | "
        "Threshold NOT relaxed (>= 10 cells). Consequence: delta (max 8 cells in any donor) has no testable donors and PP has too few donors per group (< 2 ND) -> delta/PP pseudobulk tests not performed; "
        "the delta arm of the restriction rule cannot be assessed, so no gene is declared fully beta-restricted unless alpha AND delta can be shown non-significant | "
        "Data limitation, not a parameter change. Recorded so it is not mistaken for a null result. |")
if new_rows:
    if not dev_txt.endswith("\n"):
        dev_txt += "\n"
    dev_path.write_text(dev_txt + "\n".join(new_rows) + "\n")
    log(f"\nAppended {len(new_rows)} row(s) to logs/deviations.md")
else:
    log("\nlogs/deviations.md already contains the Task 9 single-cell rows (idempotent re-run; nothing appended)")

# ------------------------------------------------------------------------------------------------ assertions
log("\n=== Final assertions ===")
assert len(qc) == 1600 and qc["qc_fail"].sum() < 0.3 * 1600
assert int(by_group["n_cells"].sum()) == 1600
assert len(obs_k) == int((~qc["qc_fail"]).sum())
assert set(expr_tab["gene"]) == set(GENES) and not expr_tab[["detection_fraction", "mean_log2_rpkm1"]].isna().any().any()
assert len(meta) > 0 and meta["sample"].is_unique and not meta[["donor", "group", "age", "sex", "n_cells"]].isna().any().any()
assert len(panel) == len(GENES) * len(res.groupby(["celltype", "model"]))
assert panel["P.Value"].dropna().between(0, 1).all() and panel["FDR_BH_9genes"].dropna().between(0, 1).all()
assert panel.groupby(["celltype", "model"]).size().eq(len(GENES)).all()
assert set(res["celltype"]) <= set(dtab.index[dtab["testable"]]), "limma ran on a cell type flagged untestable"
assert "beta" in set(res["celltype"]), "beta not tested"
assert len(verd) == len(GENES)
for f in ["12_qc_cells.csv", "12_qc_by_donor.csv", "12_marker_check.csv", "12_gene_expression_by_celltype.csv", "12_donors_per_celltype.csv",
          "12_panel_pseudobulk_results.csv", "12_beta_restriction_verdicts.csv", "12_limma_allgenes.csv.gz",
          "figures/12_dotplot_genes_by_celltype.pdf", "figures/12_dotplot_genes_by_celltype.png",
          "figures/12_pseudobulk_beta_alpha_genes.pdf", "figures/12_pseudobulk_beta_alpha_genes.png"]:
    assert (OUT / f).exists() and (OUT / f).stat().st_size > 0, f"missing output {f}"
log("All assertions passed.")
log(f"Finished {datetime.datetime.now().isoformat(timespec='seconds')}")
LOGF.close()

results/stage7 -- GSE81608 single-cell (scripts/12_singlecell.py, scripts/12b_singlecell_r.R)
GSE81608 is DISCOVERY-ADJACENT DESCRIPTIVE data, not validation. It informs no feature selection, model or panel decision.
Input is RPKM (no raw counts available): pseudobulk = mean log2(RPKM+1) per donor x cell type, limma lmFit/eBayes (deviation from
protocol §5.3 'sum counts'; see logs/deviations.md). Cell types = authors' labels; delta and PP have no testable donors at >= 10 cells.
BH is across the 9 genes of interest within each cell type/model (also provided: BH across all fitted genes). No per-cell tests.

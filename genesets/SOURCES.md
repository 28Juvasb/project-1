# Gene set sources (downloaded 2026-10-03)

| File | Source | Genes (rows / unique symbols) |
|---|---|---|
| SenMayo.csv | MSigDB human SAUL_SEN_MAYO (C2:CGP), GMT from https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=SAUL_SEN_MAYO&fileType=gmt . MSigDB lists source publication PubMed 35974106 (Saul et al., Nat Commun 2022), exact source Supplementary Data 1, HUMAN_GENE_SYMBOL. | 124 / 124 |
| CellAge.csv | CellAge build 3 (cellage3.tsv, file dated 2023-04-22) from https://genomics.senescence.info/cells/cellAge.zip | 949 / 866 (multiple rows per gene = separate evidence entries; not deduplicated here) |
| GenAge_human.csv | GenAge human (genage_human.csv, file dated 2020-02-09 inside zip) from https://genomics.senescence.info/genes/human_genes.zip ; release notes in source_files/release.html state Build 21 (2023-08-28), 307 human genes | 307 / 307 |

Unmodified downloads are in `source_files/`. Symbols are as supplied by the sources; NOT re-mapped to current HGNC (do at Task 4 union step).
Notes: SenMayo MSigDB version has 124 genes (the paper describes a 125-gene panel; not reconciled). No union built here (Task 4).

## Update 2026-10-05 (Task 4)
Symbols were later mapped to current HGNC with `limma::alias2SymbolTable` (scripts/05_genesets.R). Five aliases were updated (PHB->PHB1, ATP5O->ATP5PO, CTGF->CCN2, H2AFX->H2AX, ARNTL->BMAL1); none dropped. The source files in `source_files/` are unchanged. SenMayo remains the 124-gene MSigDB copy: the paper's 125th gene could not be retrieved (nature.com blocked automated access; PMC had no supplementary link). Union: 1151 genes (`senescence_union.rds`).

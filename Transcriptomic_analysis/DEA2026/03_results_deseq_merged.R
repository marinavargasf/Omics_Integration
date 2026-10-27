library("DESeq2")
library("biomaRt")
library(readr)

setwd("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/results")

################################################################################
## 1. Load DESeq2 results (new model, 2026)
################################################################################

res <- readRDS("./Results_rnaseq/res_deseq_2026.rds")
# Coefficient case_controlcase (control = reference): positive log2FC = higher in cases.
# NO sign flip needed (the old deseq2_output.RData used case as reference -- do not load it).

res <- res[order(res$padj), ]
sig2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 0.15)
cat("Significant genes (padj < 0.05, |log2FC| > 0.15):", nrow(sig2), "\n")
stopifnot(nrow(sig2) == 13)   # must match the 02_deseq2_2026.R output

################################################################################
## 2. Annotate with gene symbols (HGNC)
################################################################################
ensembl_ids <- rownames(res)

gene_map_path <- "../data/gene_map_v100_full.rds"

if (file.exists(gene_map_path)) {
  gene_mapping <- readRDS(gene_map_path)
  cat("Loaded cached gene mapping (", nrow(gene_mapping), "rows) from", gene_map_path, "-- skipping Ensembl query\n")
} else {
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")#, version = 100)
  gene_mapping <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id",
      values = ensembl_ids,
      mart = ensembl
  )
  saveRDS(gene_mapping, file = gene_map_path)
  cat("Saved gene mapping (", nrow(gene_mapping), "rows) to", gene_map_path, "\n")
}

annotations <- data.frame(ensembl_gene_id = ensembl_ids, gene_name = NA, stringsAsFactors = FALSE)
match_idx <- match(annotations$ensembl_gene_id, gene_mapping$ensembl_gene_id)
annotations$gene_name <- gene_mapping$hgnc_symbol[match_idx]

stopifnot(nrow(annotations) == nrow(res))

res$gene_name <- annotations$gene_name
# fall back to Ensembl ID where no HGNC symbol was found or it's blank
missing_symbol <- is.na(res$gene_name) | res$gene_name == ""
res$gene_name[missing_symbol] <- rownames(res)[missing_symbol]

res$abslogFC <- abs(res$log2FoldChange)

# order by effect size, collapse duplicate gene symbols (keep largest |log2FC|)
data_ordered <- res[order(res$abslogFC, decreasing = TRUE), ]
data_ordered <- data_ordered[!duplicated(data_ordered$gene_name), ]
data_ordered <- data_ordered[!is.na(data_ordered$gene_name), ]
rownames(data_ordered) <- data_ordered$gene_name
data_ordered <- data_ordered[, !(names(data_ordered) %in% c("gene_name"))]

data_complete <- data_ordered

# significance filter: same thresholds as 02_deseq2_2026.R
data_filtered <- data_complete[!is.na(data_complete$padj) &
                               data_complete$padj < 0.05 &
                               abs(data_complete$log2FoldChange) > 0.15, ]
data_sorted <- data_filtered[order(data_filtered$abslogFC, decreasing = TRUE), ]

saveRDS(data_complete, file = "./Results_rnaseq/dif_exp_genesymbol_complete.rds")
saveRDS(data_sorted, file = "./Results_rnaseq/dif_exp_genesymbol_sig.rds")

write.csv(data_complete, "./Results_rnaseq/dif_exp_genesymbol_complete.csv")
write.csv(data_sorted, "./Results_rnaseq/dif_exp_genesymbol_sig.csv")

cat("Annotated complete results:", nrow(data_complete), "genes\n")
cat("Annotated + filtered (padj<0.05, |log2FC|>0.15):", nrow(data_sorted), "genes",
    "| higher in cases:", sum(data_sorted$log2FoldChange > 0),
    "| lower in cases:", sum(data_sorted$log2FoldChange < 0), "\n")
print(data_sorted[, c("log2FoldChange", "padj")])
library("DESeq2")
library("biomaRt")
library(readr)

setwd("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/results/Results_rnaseq")

################################################################################
## 1. Load DESeq2 object and extract results
################################################################################

load(file = "deseq2_output.RData")
print(resultsNames(dds))

res <- as.data.frame(results(dds, name = "case_controlcontrol"))

# Coefficient is control-vs-case (case is reference level). Flipping sign here
# to match the historical convention: positive log2FC = higher in cases.
# CONFIRM this is still the direction you want before trusting downstream signs.
res$log2FoldChange <- (-1) * res$log2FoldChange
res$stat <- (-1) * res$stat

res <- res[order(res$padj), ]
sig <- subset(res, padj < 0.05)

sig2 <- subset(sig, abs(log2FoldChange)>0.15)

# dir.create("../../Results_2026", recursive = TRUE, showWarnings = FALSE)
# saveRDS(res, file = "../../Results_2026/res_deseq_full.rds")
# saveRDS(sig, file = "../../Results_2026/res_deseq_sig.rds")
cat("Significant genes (padj < 0.05):", nrow(sig), "\n")

################################################################################
## 2. Annotate with gene symbols (HGNC), pinned Ensembl v100
################################################################################

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")#, version = 100)

ensembl_ids <- rownames(res)
gene_mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
)

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

# significance filter: padj + effect size threshold (confirm 0.15 is still correct)
data_filtered <- data_complete[abs(data_complete$log2FoldChange) > 0.15 & data_complete$padj < 0.05, ]
data_sorted <- data_filtered[order(data_filtered$abslogFC, decreasing = TRUE), ]

saveRDS(data_complete, file = "dif_exp_genesymbol_complete.rds")
saveRDS(data_sorted, file = "dif_exp_genesymbol_sig.rds")

write_csv(data_complete, "dif_exp_genesymbol_complete.csv")
write.csv(data_sorted, "dif_exp_genesymbol_sig.csv")

cat("Annotated complete results:", nrow(data_complete), "genes\n")
cat("Annotated + filtered (padj<0.05, |log2FC|>0.15):", nrow(data_sorted), "genes\n")


kk <- readRDS("~/Downloads/dif_exp_genesymbol_complete (2).rds")

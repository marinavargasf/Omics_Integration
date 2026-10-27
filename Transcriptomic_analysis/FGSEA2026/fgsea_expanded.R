# install.packages("BiocManager")
# BiocManager::install(c("fgsea", "msigdbr", "dplyr", "tibble", "readr"))

#Load required libraries,
library(fgsea)
library(data.table)
library(msigdbr)
library(dplyr)
library(tibble)
library(readr)
library(pathMED)
library(ggplot2)

setwd("~/Documents/Omics_Integration/Transcriptomic_analysis")
dir.create("./FGSEA2026/Results", recursive = TRUE, showWarnings = FALSE)

#Load differential expression results,

expressionData <- as.data.frame(readRDS("./DEA2026/results/Results_rnaseq/normalized_expression_genesymbol.rds"))

deseq_results <- readRDS("./DEA2026/results/Results_rnaseq/dif_exp_genesymbol_complete.rds")
results <- deseq_results
results$gene <- rownames(deseq_results)


#Filter NA and sort by decreasing stat,
ranks <- results %>%
  filter(!is.na(stat)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(desc(stat)) %>%
  dplyr::select(gene, stat) %>%
  deframe()

data(genesetsData)
pathways_genesets <- genesetsData[["go_bp"]]

#Run fgsea

set.seed(42)  # For reproducibility
fgseaRes <- fgsea(
  pathways = pathways_genesets,
  minSize = 5,
  stats = ranks,
  #nperm = 10000
)
scores_fake <- fgseaRes
rownames(scores_fake) <- scores_fake$pathway
traslation <- ann2term(scores_fake)
saveRDS(traslation, "./FGSEA2026/Results/pathway_id_term_translation.rds")

fgseaResTidy <- fgseaRes %>%
  arrange(padj) %>%
  as_tibble()

#Process and save results,
fgseaResTidy_traslated <- fgseaResTidy %>%
  left_join(traslation, by = c("pathway" = "ID"))

fgseaResTidy_traslated %>% filter(term == "negative regulation of apoptotic process")

#Filter top adj pval (predefined threshold)
fgseaResFiltered <- fgseaResTidy_traslated %>%
  filter(padj <= 0.05 & abs(NES) >= 1.5) %>%
  as_tibble()

################################################################################
## Collapse redundant pathways (fgsea::collapsePathways)
## Works on the fgsea results + ranking, NOT on GSVA scores.
################################################################################

sig_dt <- fgseaRes[padj <= 0.1 & abs(NES) >= 1.5]
set.seed(42)
collapsed <- collapsePathways(sig_dt[order(pval)], pathways_genesets, ranks)
main_ids <- collapsed$mainPathways

fgseaResFiltered_collapse <- fgseaResFiltered %>%
  filter(pathway %in% main_ids)

cat("Significant pathways:", nrow(fgseaResFiltered),
    "| after collapsing redundancy:", nrow(fgseaResFiltered_collapse), "\n")

# Which main pathway absorbs each redundant one (useful for the supplementary table)
collapse_map <- tibble(pathway = names(collapsed$parentPathways),
                       parent  = unname(collapsed$parentPathways)) %>%
  left_join(traslation, by = c("pathway" = "ID")) %>%
  left_join(traslation %>% dplyr::rename(parent_term = term), by = c("parent" = "ID"))
write_csv(collapse_map, "./FGSEA2026/Results/fgsea_GO_collapse_parent_map.csv")

################################################################################
## Save tables (leadingEdge is a list column -> flatten for CSV)
################################################################################

save_table <- function(df, suffix = "") {
  df <- as.data.frame(df)
  csv_df <- df
  csv_df$leadingEdge <- sapply(csv_df$leadingEdge, paste, collapse = ";")
  write_csv(csv_df, paste0("./FGSEA2026/Results/fgsea_hallmark_results_filtered_GO", suffix, ".csv"))
  saveRDS(df, paste0("./FGSEA2026/Results/fgsea_hallmark_results_filtered_GO", suffix, ".rds"))
}

save_table(fgseaResFiltered)
save_table(fgseaResFiltered_collapse, "_collapse")

#Optional: plot top pathways,
topPathways <- fgseaResTidy %>%
  top_n(n = 10, wt = -padj) %>%
  arrange(desc(NES))

ggplot(topPathways, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "Top Hallmark Pathways by fgsea") +
  theme_minimal()

# ------------------------------------------------------------------------
# Top GO pathways bar plot -- build once, save to PNG, SVG, and PDF
# ------------------------------------------------------------------------

p_top_pathways <- ggplot(fgseaResFiltered, aes(reorder(term, NES), NES)) +
  geom_col(aes(fill = ((padj < 0.05)&(abs(NES)>2)))) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "Top GO Pathways by fgsea") +
  theme_minimal()

png("./FGSEA2026/Results/fgsea_GO_top_pathways.png", width = 1190, height = 900)
print(p_top_pathways)
dev.off()

svg("./FGSEA2026/Results/fgsea_GO_top_pathways.svg", width = 1190/96, height = 900/96)
print(p_top_pathways)
dev.off()

pdf("./FGSEA2026/Results/fgsea_GO_top_pathways.pdf", width = 1190/96, height = 900/96)
print(p_top_pathways)
dev.off()

################################################################################
## GSVA scores
## Scores are computed once for all significant pathways; each gene set's score
## does not depend on the other sets, so the collapsed version is just a subset.
################################################################################

pathways_sig <- pathways_genesets[names(pathways_genesets) %in% fgseaResFiltered$pathway]
scores <- getScores(expressionData, geneSets = pathways_sig, method = "GSVA", cores = 10)
saveRDS(scores, "./FGSEA2026/Results/pathwayScores_samples.rds")

save_scores <- function(res_table, suffix = "") {
  fs <- scores[rownames(scores) %in% res_table$pathway, , drop = FALSE]
  fs <- as.data.frame(fs) %>%
    rownames_to_column(var = "ID") %>%
    left_join(traslation, by = "ID")
  saveRDS(fs, paste0("./FGSEA2026/Results/pathwayScores_samplesFiltered_GO", suffix, ".rds"))

  nes <- res_table %>% dplyr::select(pathway, NES)
  saveRDS(nes, paste0("./FGSEA2026/Results/nes_GO", suffix, ".RDS"))
}

save_scores(fgseaResFiltered)
save_scores(fgseaResFiltered_collapse, "_collapse")

################################################################################
## Bubble plots -- normal and collapsed, each saved to PNG, SVG, and PDF
################################################################################

make_bubble <- function(res_table, suffix = "", title = "FGSEA GO Enrichment") {
  bubble_df <- res_table %>%
    mutate(
      padj = ifelse(padj == 0, 1e-300, padj),  # fix zero-padj BEFORE computing logFDR
      logFDR = -log10(padj)
    ) %>%
    arrange(padj)

  p <- ggplot(bubble_df, aes(x = logFDR, y = reorder(term, logFDR))) +
    geom_point(aes(size = size, color = NES), alpha = 0.8) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    scale_size(range = c(2, 6)) +
    labs(
      x = "-log10(FDR)",
      y = "Pathway",
      color = "NES",
      size = "Gene set size",
      title = title
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, size = 12)
    )

  base <- paste0("./FGSEA2026/Results/fgsea_GO_bubble", suffix)

  png(paste0(base, ".png"), width = 1200, height = 900, res = 150)
  print(p)
  dev.off()

  svg(paste0(base, ".svg"), width = 1200/150, height = 900/150)
  print(p)
  dev.off()

  pdf(paste0(base, ".pdf"), width = 1200/150, height = 900/150)
  print(p)
  dev.off()

  invisible(p)
}

make_bubble(fgseaResFiltered)
make_bubble(fgseaResFiltered_collapse, "_collapse",
            title = "FGSEA GO Enrichment (non-redundant pathways)")
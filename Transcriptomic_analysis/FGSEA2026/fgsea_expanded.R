# install.packages("BiocManager")
# BiocManager::install(c("fgsea", "msigdbr", "dplyr", "tibble", "readr"))

#Load required libraries,
library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)
library(readr)
library(pathMED)
library(ggplot2)

setwd("~/Documents/Omics_Integration/Transcriptomic_analysis")
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
fgseaResTidy_traslated <- fgseaResTidy%>%
  left_join(traslation, by = c("pathway" = "ID"))  # Ajusta "pathway" al nombre correcto

fgseaResTidy_traslated %>% filter(term == "negative regulation of apoptotic process")
#Filter top adj pval
fgseaResFiltered <- fgseaResTidy_traslated %>%
  filter(padj <= 0.1 & abs(NES) >=1.5) %>%
  as_tibble()

data_table <- as.data.frame(fgseaResFiltered)
#Save results to file,
write_csv(data_table, "./FGSEA2026/Results/fgsea_hallmark_results_filtered_GO.csv")
saveRDS(data_table, "./FGSEA2026/Results/fgsea_hallmark_results_filtered_GO.rds")


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

# plotEnrichment(pathway = pathways_genesets[["hsa05012"]],
#                stats = ranks,
#                gseaParam = 1,
#                ticksSize = 0.2)
# 
# plotEnrichment(pathway = msigdbr_list[["SA_PROGRAMMED_CELL_DEATH"]],
#                stats = ranks,
#                gseaParam = 1,
#                ticksSize = 0.2)

################################################################################
################################################################################

pathways_sig <- pathways_genesets[names(pathways_genesets) %in% fgseaResFiltered$pathway]
scores <- getScores(expressionData, geneSets = pathways_sig, method = "GSVA", cores = 10)
filtered_scores <- scores[rownames(scores) %in% fgseaResFiltered$pathway, ]


filtered_scores <- as.data.frame(filtered_scores)

filtered_scores <- filtered_scores %>%
  rownames_to_column(var = "ID") %>%
  left_join(traslation, by = "ID")  # Usa la columna ID como clave



saveRDS(scores, "./FGSEA2026/Results/pathwayScores_samples.rds")
saveRDS(filtered_scores, "./FGSEA2026/Results/pathwayScores_samplesFiltered_GO.rds")

NES <- fgseaResFiltered %>% dplyr::select(pathway,NES)
saveRDS(NES, "./FGSEA2026/Results/nes_GO.RDS")



# Prepare data
bubble_df <- fgseaResFiltered %>%
  mutate(
    padj = ifelse(padj == 0, 1e-300, padj),  # fix zero-padj BEFORE computing logFDR
    logFDR = -log10(padj)
  ) %>%
  arrange(padj)

# Optional: limit number of pathways (VERY recommended)
#bubble_df <- bubble_df %>% head(25)
nrow(fgseaResFiltered)
# ------------------------------------------------------------------------
# Bubble plot -- build once, save to PNG, SVG, and PDF
# ------------------------------------------------------------------------

p_bubble <- ggplot(bubble_df, aes(x = logFDR, y = reorder(term, logFDR))) +
  geom_point(aes(size = size, color = NES), alpha = 0.8) +
  scale_color_gradient2(
    low = "#2166AC", 
    mid = "white", 
    high = "#B2182B", 
    midpoint = 0
  ) +
  scale_size(range = c(2, 6)) +  # slightly smaller bubbles
  labs(
    x = "-log10(FDR)",
    y = "Pathway",
    color = "NES",
    size = "Gene set size",
    title = "FGSEA GO Enrichment"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

# Higher-resolution PNG with bigger dimensions
png("./FGSEA2026/Results/fgsea_GO_bubble.png", width = 1200, height = 900, res = 150)
print(p_bubble)
dev.off()

svg("./FGSEA2026/Results/fgsea_GO_bubble.svg", width = 1200/150, height = 900/150)
print(p_bubble)
dev.off()

pdf("./FGSEA2026/Results/fgsea_GO_bubble.pdf", width = 1200/150, height = 900/150)
print(p_bubble)
dev.off()

# install.packages("BiocManager")
# BiocManager::install(c("fgsea", "msigdbr", "dplyr", "tibble", "readr"))

#Load required libraries,
library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)
library(readr)
library(pathMED)

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

#Load MSigDB Hallmark gene sets (e.g., for Homo sapiens),

# msigdbr_list <- msigdbr(species = "Homo sapiens", collection = "C2") %>%
#   split(x = .$gene_symbol, f = .$gs_name)
# 
# msigdbr_list_remove <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP") %>%
#   split(x = .$gene_symbol, f = .$gs_name)
# 
# # Remove gene sets present in msigdbr_list_remove from msigdbr_list
# msigdbr_filtered <- msigdbr_list[!names(msigdbr_list) %in% names(msigdbr_list_remove)]
# 
# msigdbr_list <- msigdbr_filtered
# length(msigdbr_list)
# msigdbr_collections(db_species = "Hs")

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

fgseaResTidy <- fgseaRes %>%
  arrange(padj) %>%
  as_tibble()

#Process and save results,
fgseaResTidy_traslated <- fgseaResTidy%>%
  left_join(traslation, by = c("pathway" = "ID"))  # Ajusta "pathway" al nombre correcto


#Filter top adj pval
fgseaResFiltered <- fgseaResTidy_traslated %>%
  filter(padj <= 0.10 & abs(NES) >=1.5) %>%
  as_tibble()

data_table <- as.data.frame(fgseaResFiltered)
#Save results to file,
#write_csv(fgseaResTidy_traslated, "fgsea_hallmark_results.csv")
write_csv(data_table, "./FGSEA2026/Results/fgsea_hallmark_results_filtered_GO.csv")
saveRDS(data_table, "./FGSEA2026/Results/fgsea_hallmark_results_filtered_GO.rds")



#Optional: plot top pathways,
library(ggplot2)
topPathways <- fgseaResTidy %>%
  top_n(n = 10, wt = -padj) %>%
  arrange(desc(NES))

ggplot(topPathways, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "Top Hallmark Pathways by fgsea") +
  theme_minimal()

png("./FGSEA2026/Results/fgsea_GO_top_pathways.png", width = 1190, height = 900)

ggplot(fgseaResFiltered, aes(reorder(term, NES), NES)) +
  geom_col(aes(fill = ((padj < 0.05)&(abs(NES)>2)))) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "Top GO Pathways by fgsea") +
  theme_minimal()

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

#BiocManager::install("pathMED")

#data(genesetsData)
#pathways_genesets <- genesetsData[["kegg"]]

# Filter msigdbr_list to keep only those gene sets present in fgseaResFiltered$pathway
# msigdbr_filtered <- msigdbr_list[names(msigdbr_list) %in% fgseaResFiltered$pathway]
# 
# 

pathways_sig <- pathways_genesets[names(pathways_genesets) %in% fgseaResFiltered$pathway]
scores <- getScores(expressionData, geneSets = pathways_sig, method = "GSVA", cores = 10)
# scores <- getScores(expressionData, geneSets=pathways_genesets, method="GSVA", cores = 10)
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
    logFDR = -log10(padj),
    padj = ifelse(padj == 0, 1e-300, padj)  # avoid Inf
  ) %>%
  arrange(padj)

# Optional: limit number of pathways (VERY recommended)
bubble_df <- bubble_df %>% head(25)

# Higher-resolution PNG with bigger dimensions
png("./FGSEA2026/Results/fgsea_GO_bubble.png", width = 1200, height = 900, res = 150)  # double size for clarity

ggplot(bubble_df, aes(x = logFDR, y = reorder(term, logFDR))) +
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

dev.off()


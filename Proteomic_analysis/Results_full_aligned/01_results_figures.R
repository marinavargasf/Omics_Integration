################################################################################
## Publication-quality volcano plot — proteomics (CSF / plasma)
################################################################################
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(scales)

################################################################################
## 0. Set which compartment you're plotting
################################################################################
compartment <- "CSF"     # "PLA" or "CSF"

setwd("~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned/")

# Symbol-translated, composite-IDs-excluded, collision-disambiguated output
# from the UniProt -> HGNC translation script
data_complete <- readRDS(paste0("full_differential_expression_proteins_", compartment, "_SYMBOL.rds"))

################################################################################
## Prepare data
################################################################################
# NOTE: this uses adj.P.Val (FDR-adjusted), not raw P.Value like the earlier
# proteomics scripts in this pipeline. Confirm this is the criterion you want
# before treating these plots as final -- it changes which proteins count as
# significant, not just the label.
p_cutoff      <- 0.05
log2fc_cutoff <- 0.2

volcano_data <- data_complete %>%
  rownames_to_column("protein_id") %>%
  mutate(
    significance = case_when(
      adj.P.Val <= p_cutoff & logFC >=  log2fc_cutoff ~ "Up in cases",
      adj.P.Val <= p_cutoff & logFC <= -log2fc_cutoff ~ "Down in cases",
      TRUE ~ "Not significant"
    ),
    significance = factor(
      significance,
      levels = c("Up in cases", "Down in cases", "Not significant")
    ),
    neg_log10_p = -log10(adj.P.Val)
  )

# Diagnostic: which proteins clear -log10(p) > 10?
volcano_data %>%
  filter(neg_log10_p > 10) %>%
  arrange(adj.P.Val) %>%
  dplyr::select(protein_id, logFC, adj.P.Val, neg_log10_p)

# Avoid Inf values when adj.P.Val == 0
max_y <- max(
  volcano_data$neg_log10_p[is.finite(volcano_data$neg_log10_p)],
  na.rm = TRUE
)
volcano_data$neg_log10_p[is.infinite(volcano_data$neg_log10_p)] <- max_y

# Symmetric x-axis with extra padding so the data sits centered
max_x <- max(abs(volcano_data$logFC), na.rm = TRUE) * 1.6

################################################################################
## Select proteins to label — network-output genes for this compartment,
## plus proteins shared between CSF and plasma (computed dynamically)
################################################################################
genes_pla_network <- c(
  "SF3B4", "DDC", "RNASE3", "CLEC1B", "FGF21",
  "CDKN2D", "STK11", "PRDX1", "ADH4", "CC2D1A"
)
genes_csf_network <- c(
  "PRL", "DDC", "GCG", "TNFRSF8", "CHRDL2",
  "CBLIF", "O15197", "FABP1", "CLEC4G", "EPHA10"
)

dea_csf <- readRDS("full_differential_expression_proteins_CSF_SYMBOL.rds")
dea_pla <- readRDS("full_differential_expression_proteins_PLA_SYMBOL.rds")

sig_csf <- dea_csf %>% filter(adj.P.Val <= p_cutoff, abs(logFC) >= log2fc_cutoff)
sig_pla <- dea_pla %>% filter(adj.P.Val <= p_cutoff, abs(logFC) >= log2fc_cutoff)

genes_shared_csf_pla <- intersect(rownames(sig_csf), rownames(sig_pla))
cat("Proteins shared between CSF and plasma (n =", length(genes_shared_csf_pla), "):\n")
print(genes_shared_csf_pla)

genes_network_this_compartment <- if (compartment == "PLA") genes_pla_network else genes_csf_network

genes_to_label_manual <- unique(c(genes_network_this_compartment, genes_shared_csf_pla))

proteins_to_label <- volcano_data %>%
  filter(
    significance != "Not significant",
    protein_id %in% genes_to_label_manual
  )

missing_genes <- setdiff(genes_to_label_manual, proteins_to_label$protein_id)
if (length(missing_genes) > 0) {
  message("Requested label genes not found (or not significant) in this dataset: ",
          paste(missing_genes, collapse = ", "))
}

################################################################################
## Volcano plot
################################################################################
volcano_plot <- ggplot(
  volcano_data,
  aes(x = logFC, y = neg_log10_p)
) +
  geom_vline(
    xintercept = c(-log2fc_cutoff, log2fc_cutoff),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  geom_hline(
    yintercept = -log10(p_cutoff),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  geom_point(
    data = filter(volcano_data, significance == "Not significant"),
    aes(color = significance, alpha = significance),
    shape = 16,
    size = 1.4
  ) +
  geom_point(
    data = filter(volcano_data, significance != "Not significant"),
    aes(color = significance, alpha = significance),
    shape = 16,
    size = 2.2
  ) +
  scale_alpha_manual(
    values = c(
      "Up in cases"      = 0.85,
      "Down in cases"    = 0.85,
      "Not significant"  = 0.7
    ),
    guide = "none"
  ) +
  geom_text_repel(
    data = proteins_to_label,
    aes(label = protein_id, color = significance),
    size = 3,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.8,
    point.padding = 0.4,
    force = 2,
    force_pull = 0.5,
    min.segment.length = Inf,
    segment.color = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Up in cases"      = "#619CFF",
      "Down in cases"    = "#D9483F",
      "Not significant"  = "grey60"
    )
  ) +
  scale_x_continuous(limits = c(-max_x, max_x), breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    x = expression(log[2]~"(fold change)"),
    y = expression(-log[10]~"(adj. p-value)"),
    color = NULL
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 11),
    legend.key = element_blank(),
    axis.title = element_text(size = 11, color = "grey30"),
    axis.text = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    panel.grid.major = element_line(color = "#E8E4F0", linewidth = 0.35),
    panel.grid.minor = element_line(color = "#F2F0F7", linewidth = 0.25),
    plot.margin = margin(10, 15, 10, 10)
  )

print(volcano_plot)

################################################################################
## Save
################################################################################
ggsave(paste0("volcano_DAPs_case_vs_control_", compartment, ".pdf"), volcano_plot, width = 10, height = 7)
ggsave(paste0("volcano_DAPs_case_vs_control_", compartment, ".png"), volcano_plot, width = 10, height = 7, dpi = 300)
ggsave(paste0("volcano_DAPs_case_vs_control_", compartment, ".svg"), volcano_plot, width = 10, height = 7)

################################################################################
## Sample overlap sanity check (unrelated to the plot itself)
################################################################################
samples_csf <- read.csv(file="~/Documents/Parkinson/proteomic/Targeted olink/proteomics-CSF-PPEA-D03/releases_2023_v4release_1027_proteomics-CSF-PPEA-D03_proteomics_CSF-PPEA-D03_samples.csv")
samples_pla <- read.csv(file="~/Documents/Parkinson/proteomic/Targeted olink/proteomics-PLA-PPEA-D03/releases_2023_v4release_1027_proteomics-PLA-PPEA-D03_proteomics_PLA-PPEA-D03_samples.csv")

length(intersect(samples_csf$sample_id, samples_pla$sample_id))  # should be 0 -- different specimen types
length(intersect(samples_csf$participant_id, samples_pla$participant_id))  # can overlap -- same people, different fluid draws


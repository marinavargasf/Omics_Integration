################################################################################
## 3. Publication-quality volcano plot
################################################################################
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(scales)

setwd("~/Documents/Omics_Integration/Transcriptomic_analysis/DEA2026/results/Results_rnaseq")

data_complete <- readRDS("dif_exp_genesymbol_complete.rds")
data_sorted   <- readRDS("dif_exp_genesymbol_sig.rds")

################################################################################
## Prepare data
################################################################################
volcano_data <- data_complete %>%
  rownames_to_column("gene_name") %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange >  0.15 ~ "Up in cases",
      padj < 0.05 & log2FoldChange < -0.15 ~ "Down in cases",
      TRUE ~ "Not significant"
    ),
    significance = factor(
      significance,
      levels = c("Up in cases", "Down in cases", "Not significant")
    ),
    neg_log10_padj = -log10(padj)
  )

# Avoid Inf values when padj == 0
max_y <- max(
  volcano_data$neg_log10_padj[is.finite(volcano_data$neg_log10_padj)],
  na.rm = TRUE
)
volcano_data$neg_log10_padj[is.infinite(volcano_data$neg_log10_padj)] <- max_y

# Symmetric x-axis with extra padding so the data sits centered rather than edge-to-edge
max_x <- max(abs(volcano_data$log2FoldChange), na.rm = TRUE) * 1.6

################################################################################
## Select genes to label
################################################################################
genes_to_label <- volcano_data %>%
  filter(significance != "Not significant")

################################################################################
## Volcano plot
################################################################################
volcano_plot <- ggplot(
  volcano_data,
  aes(x = log2FoldChange, y = neg_log10_padj)
) +
  # Fold-change thresholds
  geom_vline(
    xintercept = c(-0.15, 0.15),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  # FDR threshold
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  # Non-significant points: plain filled circles, no border
  geom_point(
    data = filter(volcano_data, significance == "Not significant"),
    aes(color = significance, alpha = significance),
    shape = 16,
    size = 1.4
  ) +
  # Significant points: plain filled circles, no border
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
  # Gene labels (bold; no leader lines)
  geom_text_repel(
    data = genes_to_label,
    aes(label = gene_name, color = significance),
    size = 3,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
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
    y = expression(-log[10]~"(adjusted p-value)"),
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
ggsave(
  "volcano_DEGs_case_vs_control.pdf",
  volcano_plot,
  width = 10,
  height = 7
)
ggsave(
  "volcano_DEGs_case_vs_control.png",
  volcano_plot,
  width = 10,
  height = 7,
  dpi = 300
)
ggsave(
  "volcano_DEGs_case_vs_control.svg",
  volcano_plot,
  width = 10,
  height = 7
)
dev.off()
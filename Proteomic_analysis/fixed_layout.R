################################################################################
# Shared qgraph layout and edge helpers for all network figures
# (CSF, PLA, pathways).
#
#   source("~/Documents/Omics_Integration/Proteomic_analysis/fixed_layout.R")
################################################################################

library(qgraph)

clinical_layout_file <- "~/Documents/Omics_Integration/Proteomic_analysis/clinical_layout.rds"

clinical_order <- c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4", "MOCA",
                    "Mobility39", "ADL39", "Emotional39", "Stigma39", "Social39",
                    "Cognition39", "Communication39", "Discomfort39",
                    "Schwad_ADL", "ESS", "UPSIT")

# Run ONCE, on symm_score BEFORE muting clinical-clinical edges
make_clinical_layout <- function(symm_score) {
  lay <- qgraph(symm_score[clinical_order, clinical_order],
                layout = "spring", DoNotPlot = TRUE)$layout
  rownames(lay) <- clinical_order
  saveRDS(lay, clinical_layout_file)
}

# Every figure: clinical nodes at the saved positions (same in all networks),
# molecular nodes on an outer ring, ordered by the direction of their
# clinical partners
network_layout <- function(mat, r_in = 3.2, r_out = 5.5) {
  clin <- readRDS(clinical_layout_file)
  inn  <- match(rownames(clin), colnames(mat))
  stopifnot(!anyNA(inn))
  out  <- setdiff(seq_len(ncol(mat)), inn)

  xy <- sweep(clin, 2, colMeans(clin))              # centre
  xy <- xy / max(sqrt(rowSums(xy^2))) * r_in        # fit in inner circle

  w    <- abs(mat[out, inn, drop = FALSE])          # |t| to each clinical node
  ang  <- atan2(w %*% xy[, 2], w %*% xy[, 1])       # direction of partners
  slot <- min(ang) + 2 * pi * (seq_along(out) - 1) / length(out)

  lay <- matrix(0, ncol(mat), 2)
  lay[inn, ] <- xy
  lay[out[order(ang)], ] <- r_out * cbind(cos(slot), sin(slot))
  lay
}

# ==========================================================
# Edge label / color helpers
# ==========================================================

make_edge_labels <- function(model, threshold = 3) {
  el <- matrix(NA, nrow = nrow(model), ncol = ncol(model))
  el[abs(model) > threshold] <- round(model[abs(model) > threshold], 2)
  el
}

make_edge_colors <- function(el) {
  ec <- matrix("black", nrow = nrow(el), ncol = ncol(el))
  ec[!is.na(el) & el > 0] <- "#B0E57C"
  ec[!is.na(el) & el < 0] <- "#FF7F7F"
  ec
}

make_edge_colors_gradient <- function(model, max_val = NULL, exponent = 0.6,
                                      pos_col = "#B0E57C", neg_col = "#FF7F7F") {
  ec <- matrix(NA, nrow = nrow(model), ncol = ncol(model))
  if (is.null(max_val)) max_val <- max(abs(model), na.rm = TRUE)

  pos_palette <- colorRampPalette(c("#E6F2E6", pos_col))(100)
  neg_palette <- colorRampPalette(c("#F7E3E3", neg_col))(100)

  pos_idx <- !is.na(model) & model > 0
  neg_idx <- !is.na(model) & model < 0

  pos_norm <- (abs(model[pos_idx]) / max_val)^exponent
  neg_norm <- (abs(model[neg_idx]) / max_val)^exponent

  ec[pos_idx] <- pos_palette[pmin(100, pmax(1, round(pos_norm * 100)))]
  ec[neg_idx] <- neg_palette[pmin(100, pmax(1, round(neg_norm * 100)))]
  ec
}
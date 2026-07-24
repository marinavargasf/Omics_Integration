################################################################################
# Translate UniProt IDs -> HGNC gene symbols, ONE mapping built once and
# applied consistently across expression matrices AND DEA results (CSF + PLA)
################################################################################

library(biomaRt)
library(dplyr)

################################################################################
## 1. Read all four objects
################################################################################

expr_pla <- readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_PLA.rds")
expr_csf <- readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_CSF.rds")

dea_csf <- readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results/differential_expression_proteins_CSF.rds")
dea_pla <- readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results/differential_expression_proteins_PLA.rds")

# ------------------------------------------------------------------------------
# Exclude composite/heterodimer IDs (e.g., "P29459_P29460") from the DEA
# results before mapping -- these don't correspond to a single UniProt
# accession. (Expression matrices may also carry these in colnames; drop
# there too for consistency.)
# ------------------------------------------------------------------------------

drop_composite_cols <- function(df) {
  keep <- !grepl("_", colnames(df)[-c(1, 2)])
  df[, c(colnames(df)[1:2], colnames(df)[-c(1, 2)][keep])]
}

n_composite_dea_csf <- sum(grepl("_", rownames(dea_csf)))
n_composite_dea_pla <- sum(grepl("_", rownames(dea_pla)))
cat("Composite IDs excluded from DEA results — CSF:", n_composite_dea_csf,
    "| Plasma:", n_composite_dea_pla, "\n")

dea_csf <- dea_csf[!grepl("_", rownames(dea_csf)), ]
dea_pla <- dea_pla[!grepl("_", rownames(dea_pla)), ]

################################################################################
## 2. Build ONE UniProt -> symbol mapping covering every ID across all four objects
################################################################################

mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

uniprot_ids <- unique(c(
  colnames(expr_pla)[-c(1, 2)],
  colnames(expr_csf)[-c(1, 2)],
  rownames(dea_csf),
  rownames(dea_pla)
))

results <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters    = "uniprotswissprot",
  values     = uniprot_ids,
  mart       = mart
)
results <- as.data.frame(results)

################################################################################
## 3. Diagnostics BEFORE renaming (run once, applies to everything downstream)
################################################################################

results <- results[results$hgnc_symbol != "", ]

dupe_ids <- results$uniprotswissprot[duplicated(results$uniprotswissprot)]
if (length(dupe_ids) > 0) {
  message("UniProt IDs mapping to multiple symbols (kept first only): ",
          paste(unique(dupe_ids), collapse = ", "))
}
results <- results[!duplicated(results$uniprotswissprot), ]

dupe_symbols <- results$hgnc_symbol[duplicated(results$hgnc_symbol)]
if (length(dupe_symbols) > 0) {
  message("Gene symbols claimed by multiple UniProt IDs (collision risk): ")
  print(results[results$hgnc_symbol %in% dupe_symbols, ])
}

unmapped_ids <- setdiff(uniprot_ids, results$uniprotswissprot)
if (length(unmapped_ids) > 0) {
  message("UniProt IDs with no symbol found (left as UniProt ID): ",
          paste(unmapped_ids, collapse = ", "))
}

id_to_symbol <- setNames(results$hgnc_symbol, results$uniprotswissprot)

# cache the mapping so future reruns don't need to hit Ensembl again
saveRDS(id_to_symbol, "~/Documents/Omics_Integration/Proteomic_analysis/Results/uniprot_to_hgnc_mapping.rds")

################################################################################
## 4. Apply the SAME mapping to all four objects
################################################################################

rename_with_symbol <- function(id_vec, mapping) {
  ifelse(id_vec %in% names(mapping), mapping[id_vec], id_vec)
}

# --- expression matrices: rename colnames ---
colnames(expr_pla) <- rename_with_symbol(colnames(expr_pla), id_to_symbol)
colnames(expr_csf) <- rename_with_symbol(colnames(expr_csf), id_to_symbol)

stopifnot(
  "Duplicate column names in PLA expression matrix after renaming!" =
    !any(duplicated(colnames(expr_pla))),
  "Duplicate column names in CSF expression matrix after renaming!" =
    !any(duplicated(colnames(expr_csf)))
)

# --- DEA results: rename rownames, disambiguate collisions if any occur ---
dea_csf$uniprot_id <- rownames(dea_csf)
dea_pla$uniprot_id <- rownames(dea_pla)

rownames(dea_csf) <- rename_with_symbol(rownames(dea_csf), id_to_symbol)
rownames(dea_pla) <- rename_with_symbol(rownames(dea_pla), id_to_symbol)

disambiguate_dupes <- function(df, label) {
  dupes <- rownames(df)[duplicated(rownames(df)) | duplicated(rownames(df), fromLast = TRUE)]
  if (length(dupes) > 0) {
    message(label, ": duplicate symbol after renaming, appending UniProt ID to disambiguate: ",
            paste(unique(dupes), collapse = ", "))
    rownames(df) <- ifelse(rownames(df) %in% dupes,
                           paste0(rownames(df), "_", df$uniprot_id),
                           rownames(df))
  }
  df
}

dea_csf <- disambiguate_dupes(dea_csf, "CSF DEA")
dea_pla <- disambiguate_dupes(dea_pla, "Plasma DEA")

stopifnot(
  "Duplicate rownames in CSF DEA results after renaming!" = !any(duplicated(rownames(dea_csf))),
  "Duplicate rownames in Plasma DEA results after renaming!" = !any(duplicated(rownames(dea_pla)))
)

################################################################################
## 5. Save all four, symbol-translated and mutually consistent
################################################################################

saveRDS(expr_csf, "~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_CSF_SYMBOL.rds")
saveRDS(expr_pla, "~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_PLA_SYMBOL.rds")
saveRDS(dea_csf, "~/Documents/Omics_Integration/Proteomic_analysis/Results/differential_expression_proteins_CSF_SYMBOL.rds")
saveRDS(dea_pla, "~/Documents/Omics_Integration/Proteomic_analysis/Results/differential_expression_proteins_PLA_SYMBOL.rds")


write.csv(expr_csf, "~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_CSF_SYMBOL.csv", row.names = TRUE)
write.csv(expr_pla, "~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_PLA_SYMBOL.csv", row.names = TRUE)
write.csv(dea_csf, "~/Documents/Omics_Integration/Proteomic_analysis/Results/differential_expression_proteins_CSF_SYMBOL.csv", row.names = TRUE)
write.csv(dea_pla, "~/Documents/Omics_Integration/Proteomic_analysis/Results/differential_expression_proteins_PLA_SYMBOL.csv", row.names = TRUE)
################################################################################
## 6. Recompute shared proteins (DEA results, symbol-translated)
################################################################################

shared_proteins <- intersect(rownames(dea_csf), rownames(dea_pla))
cat("Proteins selected in both CSF and plasma (n =", length(shared_proteins), "):\n")
print(shared_proteins)

cat("CSF DEA — final protein count:", nrow(dea_csf), "\n")
cat("Plasma DEA — final protein count:", nrow(dea_pla), "\n")
cat("CSF expression matrix — protein columns:", ncol(expr_csf) - 2, "\n")
cat("Plasma expression matrix — protein columns:", ncol(expr_pla) - 2, "\n")
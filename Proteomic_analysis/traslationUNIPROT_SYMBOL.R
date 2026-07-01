################################################################################
# Translate UniProt IDs -> HGNC gene symbols (PLA + CSF), with collision checks
################################################################################

filtered_proteins_metadata_PLA <- readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_PLA.rds")
filtered_proteins_metadata_CSF <- readRDS("~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_CSF.rds")

proteins_PLA <- colnames(filtered_proteins_metadata_PLA)[-c(1:2)]
proteins_CSF <- colnames(filtered_proteins_metadata_CSF)[-c(1:2)]

dif_exp_proteins <- intersect(proteins_PLA, proteins_CSF)

# ------------------------------------------------------------------------------
# Query Ensembl for gene symbols
# ------------------------------------------------------------------------------

library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

uniprot_ids <- unique(c(proteins_PLA, proteins_CSF))

results <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters    = "uniprotswissprot",
  values     = uniprot_ids,
  mart       = mart
)
results <- as.data.frame(results)

# ------------------------------------------------------------------------------
# Diagnostics BEFORE renaming — this is what was silently breaking things
# ------------------------------------------------------------------------------

# 1. Drop rows with no symbol (empty string), they're not usable for renaming
results <- results[results$hgnc_symbol != "", ]

# 2. UniProt IDs that returned >1 symbol (one-to-many) — getBM keeps only the
#    first match per ID later via setNames, so flag these for manual review
dupe_ids <- results$uniprotswissprot[duplicated(results$uniprotswissprot)]
if (length(dupe_ids) > 0) {
  message("UniProt IDs mapping to multiple symbols (kept first only): ",
          paste(unique(dupe_ids), collapse = ", "))
}
results <- results[!duplicated(results$uniprotswissprot), ]

# 3. Symbols claimed by more than one UniProt ID — THIS is what causes a
#    silent column loss, since two distinct proteins end up with the same
#    column name and intersect()/colnames() treat them as one.
dupe_symbols <- results$hgnc_symbol[duplicated(results$hgnc_symbol)]
if (length(dupe_symbols) > 0) {
  message("Gene symbols claimed by multiple UniProt IDs (collision risk): ")
  print(results[results$hgnc_symbol %in% dupe_symbols, ])
}

# 4. UniProt IDs queried but absent from results entirely (no Ensembl hit)
unmapped_ids <- setdiff(uniprot_ids, results$uniprotswissprot)
if (length(unmapped_ids) > 0) {
  message("UniProt IDs with no symbol found (left as UniProt ID): ",
          paste(unmapped_ids, collapse = ", "))
}

# ------------------------------------------------------------------------------
# Build the mapping and rename
# ------------------------------------------------------------------------------

id_to_symbol <- setNames(results$hgnc_symbol, results$uniprotswissprot)

rename_with_symbol <- function(colnames_vec, mapping) {
  ifelse(colnames_vec %in% names(mapping), mapping[colnames_vec], colnames_vec)
}

colnames(filtered_proteins_metadata_PLA) <- rename_with_symbol(
  colnames(filtered_proteins_metadata_PLA), id_to_symbol
)
colnames(filtered_proteins_metadata_CSF) <- rename_with_symbol(
  colnames(filtered_proteins_metadata_CSF), id_to_symbol
)

# ------------------------------------------------------------------------------
# Post-rename safety check: make sure renaming didn't create duplicate
# column names within either dataframe (would silently merge two proteins)
# ------------------------------------------------------------------------------

stopifnot(
  "Duplicate column names in PLA after renaming!" =
    !any(duplicated(colnames(filtered_proteins_metadata_PLA))),
  "Duplicate column names in CSF after renaming!" =
    !any(duplicated(colnames(filtered_proteins_metadata_CSF)))
)

ncol(filtered_proteins_metadata_CSF)
 
 saveRDS(filtered_proteins_metadata_CSF, "~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_CSF_SYMBOL.rds")
 saveRDS(filtered_proteins_metadata_PLA, "~/Documents/Omics_Integration/Proteomic_analysis/Results/filtered_proteins_metadata_PLA_SYMBOL.rds")

dif_exp_proteins <- intersect(colnames(filtered_proteins_metadata_CSF),
                              colnames(filtered_proteins_metadata_PLA))

################################################################################
# Translate UniProt IDs -> HGNC gene symbols (CSF + PLA)
# Uses the cached mapping (no Ensembl query unless IDs are missing and
# QUERY_MISSING = TRUE). Produces symbol-named:
#   - full DEA tables
#   - DAP tables (same thresholds as the DEA script: |logFC| >= 0.2, adj.P.Val <= 0.05)
#   - filtered_proteins_metadata (MariNET input)
################################################################################

library(dplyr)

base <- "~/Documents/Omics_Integration/Proteomic_analysis/Results_full_aligned"
f    <- function(...) file.path(base, paste0(...))

LFC_THR       <- 0.2
PADJ_THR      <- 0.05
QUERY_MISSING <- FALSE   # TRUE = query Ensembl only for IDs absent from the cached mapping

################################################################################
## 1. Load cached mapping + manual overrides
################################################################################

mapping_file <- f("uniprot_to_hgnc_mapping_full.rds")
stopifnot("Cached mapping not found" = file.exists(mapping_file))
id_to_symbol <- readRDS(mapping_file)

# Manual fixes for IDs Ensembl does not map via uniprotswissprot
manual_map <- c(O15197 = "EPHB6")
id_to_symbol[names(manual_map)] <- manual_map

################################################################################
## 2. Read full DEA results, drop composite IDs
################################################################################

tissues <- c("CSF", "PLA")

dea_full <- lapply(setNames(tissues, tissues), function(t) {
  d <- readRDS(f("full_differential_expression_proteins_", t, ".rds"))
  comp <- grepl("_", rownames(d))
  comp_sig <- comp & abs(d$logFC) >= LFC_THR & d$adj.P.Val <= PADJ_THR
  message(t, ": ", sum(comp), " composite IDs dropped (", sum(comp_sig), " of them significant)")
  d[!comp, ]
})

################################################################################
## 3. Check coverage of the mapping; optionally fill gaps from Ensembl
################################################################################

all_ids <- unique(unlist(lapply(dea_full, rownames)))
missing <- setdiff(all_ids, names(id_to_symbol))

if (length(missing) > 0 && QUERY_MISSING) {
  library(biomaRt)
  mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  extra <- getBM(attributes = c("uniprotswissprot", "hgnc_symbol"),
                 filters = "uniprotswissprot", values = missing, mart = mart)
  extra <- extra[extra$hgnc_symbol != "" & !duplicated(extra$uniprotswissprot), ]
  id_to_symbol[extra$uniprotswissprot] <- extra$hgnc_symbol
  missing <- setdiff(all_ids, names(id_to_symbol))
}

if (length(missing) > 0) {
  message(length(missing), " UniProt IDs without symbol (kept as UniProt ID): ",
          paste(missing, collapse = ", "))
}

# Save the updated mapping (includes manual fixes / newly queried IDs)
saveRDS(id_to_symbol, mapping_file)

################################################################################
## 4. Translation helper (same rule for rows and columns)
################################################################################

translate_ids <- function(ids, mapping, label) {
  new <- ifelse(ids %in% names(mapping), unname(mapping[ids]), ids)
  dup <- new %in% new[duplicated(new)]
  if (any(dup)) {
    message(label, ": duplicate symbols, appending UniProt ID: ",
            paste(unique(new[dup]), collapse = ", "))
    new[dup] <- paste0(new[dup], "_", ids[dup])
  }
  stopifnot(!any(duplicated(new)))
  new
}

################################################################################
## 5. Full DEA + DAP tables with gene symbols
################################################################################

dap <- list()

for (t in tissues) {
  d <- dea_full[[t]]
  d$uniprot_id <- rownames(d)
  rownames(d)  <- translate_ids(rownames(d), id_to_symbol, paste(t, "DEA"))

  # Full table
  saveRDS(d, f("full_differential_expression_proteins_", t, "_SYMBOL.rds"))
  write.csv(data.frame(gene_symbol = rownames(d), d, check.names = FALSE),
            f("full_differential_expression_proteins_", t, "_SYMBOL.csv"), row.names = FALSE)

  # DAPs — same thresholds as the DEA script
  s <- d[abs(d$logFC) >= LFC_THR & d$adj.P.Val <= PADJ_THR, ]
  s <- s[order(s$P.Value), ]
  dap[[t]] <- s

  saveRDS(s, f("differential_expression_proteins_", t, "_SYMBOL.rds"))
  write.csv(data.frame(gene_symbol = rownames(s), s, check.names = FALSE),
            f("differential_expression_proteins_", t, "_SYMBOL.csv"), row.names = FALSE)

  message(t, ": ", nrow(d), " proteins, ", nrow(s), " DAPs (",
          sum(s$logFC > 0), " up / ", sum(s$logFC < 0), " down)")
}

# Sanity checks
stopifnot(rownames(dap$CSF)[dap$CSF$uniprot_id == "P20711"] == "DDC",
          rownames(dap$CSF)[dap$CSF$uniprot_id == "P01236"] == "PRL")

################################################################################
## 6. filtered_proteins_metadata (MariNET input) with gene symbols
################################################################################

id_cols_fixed <- c("participant_id", "visit_month")

for (t in tissues) {
  m <- readRDS(f("filtered_proteins_metadata_", t, ".rds"))
  prot_cols <- setdiff(colnames(m), id_cols_fixed)

  comp <- prot_cols[grepl("_", prot_cols)]
  if (length(comp) > 0) message(t, " metadata: composite columns dropped: ", paste(comp, collapse = ", "))
  m <- m[, !(colnames(m) %in% comp)]
  prot_cols <- setdiff(prot_cols, comp)

  # Must be exactly the DAP set of this run
  if (!setequal(prot_cols, dap[[t]]$uniprot_id)) {
    warning(t, " metadata columns differ from DAP list — regenerate from the same DEA run.\n",
            "  Only in metadata: ", paste(setdiff(prot_cols, dap[[t]]$uniprot_id), collapse = ", "),
            "\n  Only in DAPs: ",    paste(setdiff(dap[[t]]$uniprot_id, prot_cols), collapse = ", "))
  }

  colnames(m)[match(prot_cols, colnames(m))] <- translate_ids(prot_cols, id_to_symbol, paste(t, "metadata"))

  saveRDS(m, f("filtered_proteins_metadata_", t, "_SYMBOL.rds"))
  write.csv(m, f("filtered_proteins_metadata_", t, "_SYMBOL.csv"), row.names = FALSE)
}

################################################################################
## 7. Summary
################################################################################

shared <- intersect(rownames(dap$CSF), rownames(dap$PLA))
cat("\nDAPs shared by CSF and plasma (n =", length(shared), "):\n")
print(shared)
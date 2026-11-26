# resources/

Shared helper code and local reference caches. No participant-level data.

| File | Versioned | Used by | Origin |
|---|---|---|---|
| `network_layout.R` | yes | 12, 14, 15 | Shared qgraph layout and edge-label helpers for the MariNET network figures |
| `clinical_layout.rds` | no | `network_layout.R` | Fixed positions of the 16 clinical nodes, so all network figures share the same layout. Built with `make_clinical_layout()` from the clinical-only MariNET network |
| `gene_map_v100_full.rds` | no | 03, 04 | Ensembl BioMart (`ensembl_gene_id` → `hgnc_symbol`), built by 03 on the first run |
| `exclude_genes_ENSEMBL_v100.rds` | no | 01 | MT-, ribosomal (RPL/RPS) and globin (HB*) genes from Ensembl, built by 01 on the first run |
| `uniprot_to_hgnc_mapping_full.rds` | no | 09 | Ensembl BioMart (`uniprotswissprot` → `hgnc_symbol`) plus the manual fix `O15197 = EPHB6`, built by 09 on the first run |

The `.rds` files are git-ignored (`resources/*.rds`). The three Ensembl caches are created automatically.
Because Ensembl annotations change between releases, a fresh download may give slightly different gene or protein symbols than the published run.

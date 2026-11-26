# Omics_Integration

Code for the article on **longitudinal multi-omics integration in Parkinson's disease using MariNET** (AMP PD).

The pipeline analyses three molecular layers from the [Accelerating Medicines Partnership® Parkinson's Disease (AMP PD)](https://amp-pd.org/) programme (release v4.0):

- whole-blood **RNA-seq** (BioFIND, PDBP and PPMI cohorts)
- **CSF proteomics** (Olink Explore)
- **plasma proteomics** (Olink Explore)

Each layer goes through differential analysis between PD cases and controls. The selected molecular features (DEGs, GO pathway scores and DAPs) are then linked to longitudinal clinical scales (MDS-UPDRS I–IV, MoCA, PDQ-39, Schwab & England ADL, ESS, UPSIT) with [MariNET](https://CRAN.R-project.org/package=MariNET) networks. These networks are built from linear mixed models with a random effect per participant.

> **Data are not included.** AMP PD and PPMI data are distributed under data use agreements. Access is requested through each programme (see [Input data](#input-data)).

---

## Repository structure

```
Omics_Integration/
├── README.md
├── Omics_Integration.Rproj        # open this: sets the working directory to the repo root
├── scripts/                       # analysis pipeline, run in numerical order
│   ├── 00a_setup.R
│   ├── 00b_clinical_metadata.R
│   ├── 01_transcriptomics_preprocessing_QC.R
│   ├── ...
│   └── 16_DDC_UPDRS4_followup.R
├── resources/
│   └── network_layout.R           # shared qgraph layout and edge helpers (sourced by 12, 14, 15)
├── resources/                     # small fixed reference files, no participant data (see resources/README.md)
├── clinical_processing/
│   └── clinical_data.R            # reading and processing of the clinical variables (sourced by 00b, 12–16)
├── data/                          # NOT distributed: raw inputs (see data/README.md)
└── results/                       # created by the pipeline (git-ignored)
```

All paths are relative to the repository root through [`here::here()`](https://here.r-lib.org/), so the code runs on any machine without editing paths.

---

## Pipeline overview

```
        AMP PD v4.0 raw data  →  00b clinical variables (clinical_processing/)
                                      │
        ┌─────────────────────────────┼──────────────────────────────┐
        │ TRANSCRIPTOMICS             │ PROTEOMICS (CSF + plasma)    │
        │ 01 QC, neutrophil imputation│ 07 Bias diagnostics (swamp)  │
        │ 02 DESeq2 (cases vs ctrls)  │ 08 limma DAPs                │
        │ 03 Gene symbol annotation   │ 09 UniProt → HGNC symbols    │
        │ 04 Batch-corrected matrix   │ 10 Volcano plots             │
        │ 05 Volcano plot             │                              │
        │ 06 fgsea GO:BP + GSVA scores│                              │
        └─────────────────────────────┼──────────────────────────────┘
                                      │ 11 Cohort description
                                      ▼
                       MariNET integration with clinical data
          12 DEGs · 13 GO pathways · 14 CSF DAPs · 15 Plasma DAPs
                                      │
                                      ▼
                   16 DDC–UPDRS IV follow-up (logistic GLMM)
```

---

## Scripts

| # | Script | What it does | Main outputs |
|---|--------|--------------|--------------|
| 00a | `00a_setup.R` | Checks packages (MariNET ≥ 2.0.0), creates `results/`, checks that all inputs exist | — |
| 00b | `00b_clinical_metadata.R` | Sources `clinical_processing/clinical_data.R` (reading and processing of the AMP PD clinical variables) and saves participant-level metadata | `results/clinical/metadata.rds` |
| 01 | `01_transcriptomics_preprocessing_QC.R` | Sample QC (unmapped mismatches < 5 %, RIN ≥ 6.5), neutrophil % imputation (31-gene linear model trained on PPMI hematology), case/control definition, gene filtering (CPM > 10 in ≥ 80 % of samples; MT, ribosomal and globin genes removed), removal of rare Plate × race combinations. swamp/PCA diagnostics before and after covariate correction | `results/transcriptomics/data/`, `results/transcriptomics/bias_correction/` |
| 02 | `02_transcriptomics_DESeq2.R` | DESeq2 with TMM-derived size factors. Design: `~ case_control + visit_month + sex + race + age_at_baseline + Neutrophils + Plate + RIN_Value + Normalization_Volume`, with controls as the reference level | `results/transcriptomics/deseq2/res_deseq_2026.rds` |
| 03 | `03_transcriptomics_annotation.R` | Ensembl → HGNC annotation, DEG selection (padj < 0.05, \|log2FC\| > 0.15) | `dif_exp_genesymbol_{complete,sig}.rds` |
| 04 | `04_transcriptomics_normalization.R` | TMM logCPM + `limma::removeBatchEffect` on nuisance covariates (keeps diagnosis and visit). This matrix is the MariNET input and is **never** passed to DESeq2 | `normalized_expression_genesymbol.rds` |
| 05 | `05_transcriptomics_volcano.R` | Volcano plot of DEGs | `volcano_DEGs_case_vs_control.{pdf,png,svg}` |
| 06 | `06_transcriptomics_fgsea_GO.R` | fgsea on the full ranked list (DESeq2 Wald statistic), GO:BP gene sets from `pathMED`, padj ≤ 0.05 and \|NES\| ≥ 1.5, redundancy collapsed with `collapsePathways`, per-sample GSVA scores | `results/pathways/` |
| 07 | `07_proteomics_bias_diagnostics.R` | swamp/PCA diagnostics of the raw Olink NPX matrices (CSF and plasma). Diagnostic only | `results/proteomics/bias_correction/` |
| 08 | `08_proteomics_DAP_limma.R` | limma, cases vs controls (`~0 + condition + visit_month + sex + age_at_baseline`), DAPs: adj. P ≤ 0.05 and \|logFC\| ≥ 0.2, for CSF and for plasma | `results/proteomics/*differential_expression_proteins_{CSF,PLA}.rds` |
| 09 | `09_proteomics_uniprot_to_symbol.R` | UniProt → HGNC translation using a single fixed mapping; composite (heterodimer) IDs are removed | `*_SYMBOL.rds` (MariNET input) |
| 10 | `10_proteomics_volcano.R` | Volcano plots of DAPs (see [run notes](#run-notes)) | `volcano_DAPs_case_vs_control_{CSF,PLA}.*` |
| 11 | `11_cohort_description.R` | Participants, samples, demographics and visits per omic layer, plus cross-layer overlap | console output |
| 12 | `12_MariNET_transcriptomics.R` | MariNET network: DEGs × clinical variables (PD cases), top-10 genes by max \|t\| | `results/transcriptomics/marinet/` |
| 13 | `13_MariNET_pathways_GO.R` | MariNET network: GSVA GO pathway scores × clinical variables | `results/pathways/GO_clinical_score_matrix_MariNET2.csv` |
| 14 | `14_MariNET_proteomics_CSF.R` | MariNET network: CSF DAPs × clinical variables, top-10 DAPs, and DAPs shared with plasma | `results/proteomics/CSF_*` |
| 15 | `15_MariNET_proteomics_PLA.R` | MariNET network: plasma DAPs × clinical variables, top-10 DAPs, and DAPs shared with CSF (drawn with the CSF layout so both figures match) | `results/proteomics/PLA_*` |
| 16 | `16_DDC_UPDRS4_followup.R` | Follow-up of the CSF DDC–UPDRS IV association: UPDRS IV > 0 modelled with a logistic mixed model (`glmer`), adjusted for PD treatment, with a DDC × time interaction | `results/proteomics/DDC_UPDRS4_followup/` |

**MariNET models (12–15).** These use `MariNET::lmm_integration_analysis()` with random effects `(1|sex) + (1|participant_id) + (1|age_at_baseline)`. Each clinical variable is modelled on the molecular features and the other clinical variables. Each molecular feature is modelled on the clinical variables. The edge weight is the mean of both t-values, and only molecular–clinical edges are kept.

---

## Input data

### 1. AMP PD release v4.0 (raw data)

Request access at <https://amp-pd.org/>, download the files below and place them in `data/raw/` with this layout (original AMP PD file names):

```
data/
├── raw/
│   ├── AMP_PD_Data_v3_2022/
│   │   ├── clinical_metadata_v4_2023/          # full folder as downloaded from AMP PD
│   │   │   ├── releases_2023_v4release_1027_amp_pd_participants.csv
│   │   │   ├── releases_2023_v4release_1027_amp_pd_case_control.csv
│   │   │   ├── ...
│   │   │   └── clinical/                       # all clinical tables and their _dictionary.csv files
│   │   │       ├── releases_2023_v4release_1027_clinical_Demographics.csv
│   │   │       ├── releases_2023_v4release_1027_clinical_PD_Medical_History.csv
│   │   │       ├── releases_2023_v4release_1027_clinical_MDS_UPDRS_Part_{I,II,III,IV}.csv
│   │   │       └── ...                         # MOCA, PDQ_39, UPSIT, Epworth, ... (see clinical_processing/clinical_data.R)
│   │   └── transcriptomics_v4/
│   │       ├── reports/
│   │       │   ├── releases_2023_v4release_1027_rnaseq-WB-RWTS_sequencing_metrics_rna_quality_metrics.csv
│   │       │   └── releases_2023_v4release_1027_rnaseq-WB-RWTS_reports_multiqc_multiqc_data_multiqc_star.txt
│   │       └── subread_feature-counts/
│   │           └── releases_2023_v4release_1027_rnaseq-WB-RWTS_subread_feature-counts_matrix.featureCounts.tsv
│   ├── proteomics/
│   │   └── Targeted olink/
│   │       ├── proteomics-CSF-PPEA-D03/
│   │       │   ├── releases_2023_v4release_1027_proteomics-CSF-PPEA-D03_proteomics_CSF-PPEA-D03_samples.csv
│   │       │   ├── releases_2023_v4release_1027_proteomics-CSF-PPEA-D03_ColumnHeaderGlossary.csv
│   │       │   └── protein-expression/
│   │       │       └── ..._olink-explore_protein-expression_CSF-PPEA-D03_matrix_{cardiometabolic,inflammation,neurology,oncology}.csv
│   │       └── proteomics-PLA-PPEA-D03/        # same structure, PLA instead of CSF
│   └── PPMI/
│       └── Blood_Chemistry___Hematology_10Mar2025.csv
```

### 2. PPMI hematology

`Blood_Chemistry___Hematology_10Mar2025.csv` is the *Blood Chemistry & Hematology* table from [PPMI](https://www.ppmi-info.org/) (LONI IDA, downloaded 10 March 2025). Its measured neutrophil percentages train the model that imputes neutrophil % for every RNA-seq sample (script 01).

### 3. Clinical variables

The clinical variables are read and processed by [`clinical_processing/clinical_data.R`](clinical_processing/README.md), directly from the AMP PD clinical tables. `00b_clinical_metadata.R` saves the participant-level metadata, and scripts 12–16 source the same file to build the visit-level clinical table.

### 4. Fixed reference files (`resources/`)

Small lookup tables with no participant data, kept under version control so that the gene and protein identifiers, and the network layouts, match the published results exactly. See [`resources/README.md`](resources/README.md).

---

## Software

Requires R and Bioconductor (exact versions of the published run: see the `sessionInfo` file below). Install the dependencies with:

```r
install.packages("BiocManager")
BiocManager::install(c(
  # Bioconductor
  "DESeq2", "edgeR", "limma", "biomaRt", "Biobase", "BiocParallel", "fgsea", "pathMED",
  # CRAN
  "here", "dplyr", "tidyr", "purrr", "tibble", "readr", "data.table", "ggplot2", "ggrepel",
  "scales", "pheatmap", "rstatix", "ggpubr", "lme4", "qgraph", "ggeffects", "MuMIn",
  "swamp", "msigdbr", "MariNET", "plyr", "table1", "psych", "bootnet"
))
packageVersion("MariNET")  # must be >= 2.0.0
```

`02_transcriptomics_DESeq2.R` writes the full `sessionInfo()` of the DESeq2 run to `results/transcriptomics/deseq2/sessionInfo_deseq2_2026.txt`.

---

## How to reproduce

1. Clone the repository and open `Omics_Integration.Rproj` in RStudio, or start R from the repository root.
2. Place the input data as described in [Input data](#input-data).
3. Run `scripts/00a_setup.R`. It stops and lists anything that is missing.
4. Run the scripts **in numerical order** (00b → 16), each in a fresh R session:

```bash
Rscript scripts/00a_setup.R
Rscript scripts/00b_clinical_metadata.R
for i in 01 02 03 04 05 06 07 08 09; do Rscript scripts/${i}_*.R; done
Rscript scripts/10_proteomics_volcano.R            # compartment <- "CSF"  (edit to "PLA" and run again)
Rscript scripts/11_cohort_description.R
Rscript scripts/12_MariNET_transcriptomics.R
Rscript scripts/13_MariNET_pathways_GO.R
Rscript scripts/14_MariNET_proteomics_CSF.R
Rscript scripts/15_MariNET_proteomics_PLA.R
Rscript scripts/16_DDC_UPDRS4_followup.R
```

### Run notes

- **10 – volcano plots.** The tissue is set at the top of the script (`compartment <- "CSF"`). Run it once with `"CSF"` and once with `"PLA"`.
- **14 → 15 – shared-protein networks.** The plasma script (15) reads `CSF_proteins_clinical_score_matrix.csv` and `CSF_commonPLA_allproteins_score_matrix.csv` from the CSF script (14), and reuses the CSF layout so that both shared-protein figures place the nodes in the same positions. Run 14 before 15.
- **Expected checkpoints.** The scripts contain `stopifnot()` checks against the published run: 5,731 RNA-seq samples after QC (scripts 02 and 04) and 13 DEGs (script 03). If one of them fails, the input data or package versions differ from the original analysis.
- **Runtime.** Script 02 (DESeq2 on about 5,700 samples) is the slowest step. It runs in parallel with 8 workers (`MulticoreParam(workers = 8)`); change this to suit your machine.

---

## Outputs

```
results/
├── clinical/              # participant-level clinical metadata
├── transcriptomics/
│   ├── data/              # QC'd counts and aligned sample metadata (DESeq2 input)
│   ├── bias_correction/   # swamp / PCA diagnostics
│   ├── deseq2/            # DESeq2 results, DEGs, normalized expression matrix, volcano
│   └── marinet/           # DEG–clinical network
├── pathways/              # fgsea results, GSVA scores, GO–clinical network
└── proteomics/            # DAPs (CSF, plasma), symbol tables, volcanoes, networks
    ├── bias_correction/
    └── DDC_UPDRS4_followup/
```

---

## Citation

If you use this code, please cite the article (reference to be added on publication) and the MariNET package:

> MariNET (R package). <https://CRAN.R-project.org/package=MariNET>

## Acknowledgements

Data used in the preparation of this article were obtained from the AMP PD Knowledge Platform (<https://www.amp-pd.org>), including the BioFIND, PDBP and PPMI cohorts. PPMI hematology data were obtained from the PPMI database (<https://www.ppmi-info.org/access-data-specimens/download-data>).

## Contact

Marina Vargas-Fernández — University of Granada / GENYO.

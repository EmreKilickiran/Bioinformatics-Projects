# High-Dimensional Data Analysis and Statistical Modeling of Gene Regulatory Networks in Liver Cancer

A computational pipeline to quantify transcription factor activity from high-dimensional genomic datasets and identify statistically significant regulatory patterns associated with tumor progression in hepatocellular carcinoma (TCGA-LIHC cohort).

## Overview

This project analyzes large-scale biomedical datasets totaling over **9.5 million data points** - an expression matrix of 20,530 genes across 420 samples and a regulatory network of **1.02 million interaction records** - to infer transcription factor activity and its association with clinical biomarkers.

The VIPER algorithm (Virtual Inference of Protein-activity by Enriched Regulon analysis) is used to transform high-dimensional gene expression matrices into interpretable activity scores. Unlike simple mRNA-level analysis, VIPER estimates TF activity by evaluating the coordinated expression of its ChIP-seq-derived target genes, providing a more functional interpretation of regulatory influence.

### Pipeline Architecture

```
┌──────────────────────────────────────────────────────────────┐
│                    INPUT DATA                                 │
│  TCGA-LIHC: 20,530 genes × 420 samples (RNA-Seq)            │
│  ChIP-Atlas: TF binding targets from HepG2 liver cells       │
│  Regulatory network: 1.02M interaction records                │
└─────────────────────────┬────────────────────────────────────┘
                          │
              ┌───────────▼────────────┐
              │  1. Load Expression    │  TCGA-LIHC normalized
              │     Matrix             │  RNA-Seq counts
              │  2. Classify Samples   │  Tumor vs Normal
              │     (TCGA barcodes)    │  (barcode suffix)
              └───────────┬────────────┘
                          │
            ┌─────────────┴─────────────┐
            │                           │
   ┌────────▼─────────┐      ┌─────────▼────────┐
   │  NFATC3           │      │  NFAT5            │
   │                   │      │                   │
   │  ChIP-Atlas       │      │  ChIP-Atlas       │
   │  regulon (HepG2)  │      │  regulon (HepG2)  │
   └────────┬──────────┘      └─────────┬─────────┘
            │                           │
   ┌────────▼──────────┐      ┌─────────▼─────────┐
   │  VIPER Activity   │      │  VIPER Activity    │
   │  Scoring          │      │  Scoring           │
   └────────┬──────────┘      └─────────┬──────────┘
            │                           │
            └─────────────┬─────────────┘
                          │
              ┌───────────▼────────────┐
              │  Correlation Analysis  │
              │                        │
              │  • Spearman (global)   │
              │  • Stratified by       │
              │    Tumor / Normal      │
              │  • AFP biomarker       │
              └───────────┬────────────┘
                          │
              ┌───────────▼────────────┐
              │  Visualization         │
              │                        │
              │  • Scatter plots       │
              │    (Activity vs AFP)   │
              │  • Color: sample class │
              │  • Statistics overlay  │
              └────────────────────────┘
```

## Key Findings

- **NFATC3** showed a threshold effect: above a certain activity level, all samples were tumors with no normal samples present, suggesting tumor-specific regulation patterns.
- **NFAT5** displayed a positive trend between activity and AFP biomarker expression, concentrated in tumor cases, indicating a potential role in AFP upregulation in liver cancer.
- These findings suggest NFAT-family factors as potential computational biomarkers for liver tumor classification.

## Methods

### VIPER Algorithm
VIPER estimates TF activity by evaluating the expression levels of its target genes, rather than the expression of the TF itself. A regulon object is constructed from ChIP-seq-derived binding targets, where each gene is assigned a binding weight (likelihood). VIPER then assesses how consistently the expression of these target genes reflects activation or repression, producing an activity score for each sample.

### Data Sources
- **Gene expression**: TCGA-LIHC (The Cancer Genome Atlas, Liver Hepatocellular Carcinoma), RNA-Seq normalized counts
- **TF binding targets**: [ChIP-Atlas](https://chip-atlas.org) — liver-specific ChIP-seq experiments in HepG2 cells
- **Sample classification**: TCGA barcode convention (`-01` = Primary Tumor, `-11` = Solid Tissue Normal)

## Requirements

| Package    | Source       | Purpose                           |
|------------|-------------|-----------------------------------|
| viper      | Bioconductor| TF activity inference              |
| tidyverse  | CRAN        | Data manipulation                  |
| ggplot2    | CRAN        | Visualization                      |
| openxlsx   | CRAN        | Excel I/O                          |

## Repository Structure

```
.
├── README.md
├── pipeline.R                          # Complete analysis pipeline
└── data/                               # (not tracked)
    ├── LIHC_expression_data_clean.xlsx  # TCGA-LIHC expression matrix
    ├── NFATC3.Liver.tsv                # ChIP-Atlas: NFATC3 binding targets
    └── NFAT5.Liver.tsv                 # ChIP-Atlas: NFAT5 binding targets
```

## Usage

```r
# Run the full pipeline
source("pipeline.R")

# Or from the command line
Rscript pipeline.R
```

### Data Setup

Place input files in the `data/` directory:

1. **`LIHC_expression_data_clean.xlsx`** - TCGA-LIHC RNA-Seq normalized expression matrix (genes as rows, TCGA sample barcodes as columns)
2. **`NFATC3.Liver.tsv`** - ChIP-Atlas export for NFATC3 in liver tissue (HepG2). Download from [ChIP-Atlas CoLo](https://chip-atlas.org/colo): select Antigen = NFATC3, Cell Type = Liver
3. **`NFAT5.Liver.tsv`** - Same as above for NFAT5

## Output Files

| File | Description |
|------|-------------|
| `NFATC3_activity_scores.xlsx` | Per-sample activity scores, AFP expression, and classification |
| `NFAT5_activity_scores.xlsx` | Same for NFAT5 |
| `NFATC3_vs_AFP.pdf` | Scatter plot: NFATC3 activity vs AFP expression |
| `NFAT5_vs_AFP.pdf` | Scatter plot: NFAT5 activity vs AFP expression |
| `NFATC3_correlation_stats.xlsx` | Spearman correlation statistics (global + stratified) |
| `NFAT5_correlation_stats.xlsx` | Same for NFAT5 |
| `combined_correlation_summary.xlsx` | Consolidated statistics for all TFs |

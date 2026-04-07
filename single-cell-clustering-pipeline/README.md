# Unsupervised Clustering and Dimensionality Reduction Pipeline for Single-Cell Transcriptomic Profiles

A complete single-cell RNA-seq analysis pipeline processing raw count matrices for **310,000+ single cells**, identifying cell-type-specific expression patterns across six sequential liver disease progression stages (NORM → HO → STEA → NASH → FIB → CIRR).

## Overview

This self-directed project builds an end-to-end scRNA-seq pipeline from raw counts to biological interpretation. The pipeline defines custom QC thresholds, performs unsupervised clustering to identify 11 cell types across 26 clusters, and designs a custom differential expression ratio score to quantify directional transcriptional shifts across disease stages.

### Pipeline Architecture

```
┌────────────────────────────────────────────────────────────┐
│                  RAW COUNT MATRIX                           │
│           310,000+ cells × 20,000+ genes                   │
└──────────────────────┬─────────────────────────────────────┘
                       │
         ┌─────────────▼──────────────┐
         │  Module 1: Preprocessing   │
         │                            │
         │  • QC: nFeature, nCount,   │
         │    percent.mt (≤15%)       │
         │  • Residual-based doublet  │
         │    filtering (3×MAD)       │
         │  • LogNormalize (SF=10⁴)   │
         │  • HVG: VST, top 3,000    │
         │  • Scaling                 │
         └─────────────┬──────────────┘
                       │
         ┌─────────────▼──────────────┐
         │  Module 2: Dim. Reduction  │
         │                            │
         │  • PCA (50 PCs computed)   │
         │  • Elbow plot → 20 PCs     │
         │  • Resolution sweep        │
         │    (0.1–1.0) → 0.15        │
         │  • SNN clustering → 26     │
         │    clusters                │
         │  • UMAP embedding          │
         └─────────────┬──────────────┘
                       │
         ┌─────────────▼──────────────┐
         │  Module 3: Annotation      │
         │                            │
         │  • SingleR (pseudobulk)    │
         │  • Reference: MacParland   │
         │    liver atlas (GSE115469) │
         │  • 26 clusters → 11 cell   │
         │    types                   │
         └─────────────┬──────────────┘
                       │
         ┌─────────────▼──────────────┐
         │  Module 4: DE Analysis     │
         │                            │
         │  • Per cell type × stage   │
         │    (vs NORM reference)     │
         │  • LR test, FDR < 0.05    │
         │  • Ratio score heatmap:   │
         │    (Up-Down)/(Up+Down)     │
         │  • Marker gene dot plots   │
         └─────────────┬──────────────┘
                       │
         ┌─────────────▼──────────────┐
         │  Module 5: Enrichment      │
         │                            │
         │  • GO Biological Process   │
         │  • KEGG Pathways           │
         │  • Per cell type (A) and   │
         │    per treatment (B)       │
         └────────────────────────────┘
```

## Key Design Decisions

| Decision | Value | Rationale |
|----------|-------|-----------|
| MT% threshold | 15% | Balances cell viability filtering with liver-specific high mitochondrial content |
| Doublet filtering | 3×MAD residual | Data-driven removal of outliers from nCount–nFeature trend |
| HVGs | 3,000 (VST) | Captures sufficient biological variation for liver heterogeneity |
| PCs for clustering | 20 | Determined from elbow plot inflection point |
| Resolution | 0.15 | Selected from sweep (0.1–1.0) for optimal cluster separation quality |
| Clusters → Cell types | 26 → 11 | SingleR pseudobulk against MacParland liver atlas |
| DE reference | NORM | All stages compared to healthy tissue baseline |
| DE ratio score | (Up−Down)/(Up+Down) | Custom metric for directional balance of transcriptional shifts |

## Results

The pipeline revealed that immune cell types (e.g., Mature B Cells, Hepatocyte subtypes) showed progressively stronger transcriptional shifts from healthy tissue toward advanced fibrosis and cirrhosis stages. Marker gene analysis confirmed cluster-specific expression signatures across treatment groups.

## Requirements

| Package | Source | Purpose |
|---------|--------|---------|
| Seurat | CRAN | Core scRNA-seq framework |
| SingleR | Bioconductor | Cell type annotation |
| SingleCellExperiment | Bioconductor | SCE data structures |
| clusterProfiler | Bioconductor | GO/KEGG enrichment |
| org.Hs.eg.db | Bioconductor | Human gene annotations |
| pheatmap | CRAN | Heatmap visualization |
| data.table | CRAN | Fast data I/O |
| tidyverse | CRAN | Data wrangling |
| ggplot2 | CRAN | Visualization |
| future | CRAN | Parallel computing |

## Repository Structure

```
.
├── README.md
├── run_all.R                               # Pipeline entry point
├── R/
│   ├── 00_config.R                         # Shared paths & parameters
│   ├── 01_preprocessing.R                  # QC → Normalization → HVG → Scaling
│   ├── 02_dimensionality_reduction.R       # PCA → Clustering → UMAP
│   ├── 03_cell_type_annotation.R           # SingleR pseudobulk annotation
│   ├── 04_differential_expression.R        # DE analysis + ratio heatmap + dotplots
│   └── 05_functional_enrichment.R          # GO BP + KEGG (cell-type & treatment level)
└── data/                                   # (not tracked)
    ├── raw.rds                             # Raw count matrix
    ├── GSE115469_Data.csv.gz               # MacParland liver atlas expression
    └── GSE115469_CellClusterType.txt.gz    # MacParland cell type labels
```

## Usage

```r
# Full pipeline
Rscript run_all.R

# Individual modules
Rscript run_all.R 1      # Preprocessing only
Rscript run_all.R 4 5    # DE analysis + enrichment

# HPC (SLURM)
srun --cpus-per-task=8 --mem=64G --time=12:00:00 Rscript run_all.R
```

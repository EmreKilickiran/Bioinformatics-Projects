# Computational Integration and High-Dimensional Correlation Analysis for Oral-Gut-Liver Axis Microbiome Dynamics

A four-stage computational pipeline integrating high-dimensional microbiome (16S rRNA) and transcriptomic (RNA-Seq) datasets, investigating the oral-gut-liver axis in chronic liver disease across 20+ experimental groups spanning three tissue sites and six treatment conditions.

## Overview

The underlying system was modeled as a high-dimensional interaction network between microbial community composition and host gene expression, quantified through cross-domain statistical associations. The pipeline processes data from a multi-site mouse study of periodontitis-driven liver fibroinflammation, using Western Style Diet (WSD) with ligature, antibiotic (AB), and proton pump inhibitor (PPI) treatment arms.

### Pipeline Architecture

```
┌──────────────────────────────────────────────────────────────┐
│              INPUT DATA                                      │
│  16S rRNA: Species/Genus/Phylum abundance matrices           │
│  RNA-Seq:  TPM-normalized gene expression (~20,000 genes)    │
│  Metadata: 3 body sites × 6 treatment conditions             │
└─────────────────────┬────────────────────────────────────────┘
                      │
        ┌─────────────┼─────────────┐
        │             │             │
        ▼             ▼             ▼
┌──────────────┐ ┌──────────┐ ┌──────────────┐
│ Module 1     │ │ Module 2 │ │ Module 3     │
│ Relative     │ │ Alpha    │ │ Beta         │
│ Abundance    │ │ Diversity│ │ Diversity    │
│              │ │          │ │              │
│ • phyloseq   │ │ • Shannon│ │ • Jaccard    │
│ • Phylum top │ │ • Chao1  │ │ • Wt UniFrac │
│   10 (≥1%)   │ │ • Simpson│ │ • PCoA       │
│ • Faceted    │ │ • Kruskal│ │ • PERMANOVA  │
│   barplot    │ │   Wallis │ │   (999 perm) │
│              │ │ • Wilcox │ │ • Spider +   │
│              │ │   pairwis│ │   ellipse    │
└──────────────┘ └──────────┘ └──────┬───────┘
                                     │
                      ┌──────────────┴──────────────┐
                      │                             │
                      ▼                             ▼
              ┌──────────────┐              ┌──────────────┐
              │ Module 4     │              │ Module 5     │
              │ Differential │              │ Cross-Omics  │
              │ Abundance    │              │ Correlation  │
              │              │              │              │
              │ • ANCOM-BC2  │              │ • Spearman   │
              │ • BH correct.│              │   (rank-     │
              │ • Heatmap    │              │   based)     │
              │ • LEfSe-style│              │ • 52M+ pairs │
              │   LDA plot   │              │ • FDR < 0.05 │
              │ • Boxplots   │              │ • Per tissue │
              └──────────────┘              │ • Heatmaps   │
                                            └──────────────┘
```

## Computational Methods

### 1. Taxonomic Composition (Module 1)
Constructs a taxonomic abundance matrix from raw 16S sequencing data using phyloseq in R, applying relative abundance normalization and threshold-based filtering to isolate dominant phyla (top 10, ≥1%) across all experimental groups. Visualized as double-strip faceted barplots (Bodysite × Treatment).

### 2. Alpha Diversity (Module 2)
Computes four diversity indices (Observed, Chao1, Shannon, Simpson) with Kruskal-Wallis global tests and pairwise Wilcoxon rank-sum tests (Benjamini-Hochberg corrected). Results are presented as violin+boxplots with significance annotations and p-value heatmaps.

### 3. Beta Diversity (Module 3)
Computes Jaccard dissimilarity and Weighted UniFrac distance matrices, applies PCoA for dimensionality reduction, and quantifies group separation through global and pairwise PERMANOVA tests (999 permutations). Includes spider diagrams with centroids, 95% confidence ellipses, and bodysite-specific swab comparisons.

### 4. Differential Abundance (Module 4)
Identifies differentially abundant species using ANCOM-BC2 with Benjamini-Hochberg correction for multiple testing. Generates row-scaled heatmaps annotated by Treatment and Bodysite, boxplot/violin plots for top species, and LEfSe-style log2 fold change barplots.

### 5. Cross-Omics Correlation (Module 5)
Integrates 16S rRNA genus-level abundances with TPM-normalized RNA-Seq transcriptomic profiles. The pipeline implements rank-based Spearman correlation with manually computed test statistics (t-distribution conversion, two-tailed p-values) and Benjamini-Hochberg FDR correction, processing over **52 million gene–taxa correlation pairs** through vectorized matrix operations on the HPC cluster. Executed independently across three body sites to capture tissue-specific interaction patterns.

## Requirements

| Package     | Purpose                                    |
|-------------|--------------------------------------------|
| phyloseq    | Microbiome data handling                   |
| vegan       | Ecological statistics (PERMANOVA, Jaccard) |
| ape         | PCoA, phylogenetic analysis                |
| ANCOMBC     | Differential abundance (ANCOM-BC2)         |
| pheatmap    | Heatmap visualization                      |
| ggplot2     | Statistical graphics                       |
| ggpubr      | Significance annotations                   |
| patchwork   | Multi-panel plot composition               |
| ggrepel     | Non-overlapping text labels                |
| tidyverse   | Data wrangling                             |
| openxlsx    | Excel I/O                                  |
| readxl      | Excel import                               |
| reshape2    | Data reshaping                             |

## Repository Structure

```
.
├── README.md
├── run_all.R                           # Main entry point
├── R/
│   ├── 00_config.R                     # Shared configuration & utilities
│   ├── 01_relative_abundance.R         # Phylum-level barplots
│   ├── 02_alpha_diversity.R            # Alpha diversity + Wilcoxon tests
│   ├── 03_beta_diversity.R             # Jaccard + UniFrac + PERMANOVA
│   ├── 04_differential_abundance.R     # ANCOM-BC2 analysis
│   └── 05_cross_omics_correlation.R    # 16S × RNA-Seq Spearman pipeline
├── data/                               # (not tracked - add your data here)
│   ├── DataExport.xlsx
│   ├── MappingFile-v12.txt
│   ├── TPM_normalized_data.csv
│   └── merged_metadata_corrected.csv
└── results/                            # (generated - not tracked)
```

## Usage

```r
# Run the full pipeline
Rscript run_all.R

# Run specific modules
Rscript run_all.R 1        # Only relative abundance
Rscript run_all.R 3 5      # Beta diversity + cross-omics correlation

# Or source individual modules in an R session
source("R/00_config.R")
source("R/03_beta_diversity.R")
```

### Data Setup

Place your input files in the `data/` directory:
- `DataExport.xlsx` - 16S abundance tables (sheets: `Abund_Phylum_wide`, `Abund_Species_wide`, `Abund_Genus_wide`, `Alpha_Diversity_wide`, `Beta_Diversity_wunifrac`)
- `MappingFile-v12.txt` - Sample metadata (tab-delimited)
- `TPM_normalized_data.csv` - RNA-Seq TPM expression matrix
- `merged_metadata_corrected.csv` - Merged metadata linking SampleID ↔ AnimalID

## Key Parameters

| Parameter | Value | Module |
|-----------|-------|--------|
| Abundance threshold | ≥ 1% | 1 |
| Top phyla | 10 | 1 |
| Alpha metrics | Observed, Chao1, Shannon, Simpson | 2 |
| P-value correction | Benjamini-Hochberg | 2, 4, 5 |
| PERMANOVA permutations | 999 | 3 |
| Confidence ellipse | 95% | 3 |
| ANCOM-BC2 prevalence cutoff | 10% | 4 |
| Correlation FDR threshold | 0.05 | 5 |
| Correlation |rho| threshold | 0.5 | 5 |
| Top genes | 100 | 5 |

## Context

This pipeline was developed as part of a HiWi research position at the **Schneider Lab, Center for Regenerative Therapies Dresden (CRTD), TU Dresden**, investigating how periodontitis-driven oral dysbiosis reshapes the gut microbiome and accelerates liver fibroinflammation through the oral-gut-liver axis. 

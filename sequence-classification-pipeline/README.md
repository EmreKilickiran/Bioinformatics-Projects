# End-to-End Computational Pipeline for Automated Sequence Classification and Phylogenetic Analysis

A fully automated, multi-stage computational pipeline to process raw high-throughput 16S rRNA sequencing data into classified taxonomic profiles and phylogenetic tree structures. Applied to gut microbiome characterization under **bile acid-associated liver conditions**.

## Overview

This pipeline automates the transformation of raw paired-end Illumina sequencing data (3.6 GB, V3-V4 region of 16S rRNA gene) into analysis-ready biological profiles through the following computational stages:

```
Raw FASTQ → Lane Merging → Primer Removal → DADA2 Denoising → Taxonomic Classification → Filtering → Phylogenetic Tree → Diversity Analysis → Export
```

### Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    RAW SEQUENCING DATA                          │
│              Paired-end Illumina (2 × 300 bp)                   │
│                    ~3.6 GB, 2 lanes                             │
└──────────────────────┬──────────────────────────────────────────┘
                       │
            ┌──────────▼──────────┐
            │  1. Lane Merging    │  Concatenate L1 + L2 per sample
            │  2. QIIME 2 Import  │  Manifest-based paired-end import
            └──────────┬──────────┘
                       │
            ┌──────────▼──────────┐
            │  3. Primer Removal  │  Cutadapt: 341F / 785R
            │     (Cutadapt)      │  Discard untrimmed reads
            └──────────┬──────────┘
                       │
            ┌──────────▼──────────┐
            │  4. DADA2 Denoising │  Truncation: F=282bp, R=245bp
            │                     │  Error correction, merging,
            │                     │  chimera & PhiX removal
            └──────────┬──────────┘
                       │
              ┌────────┴────────┐
              │                 │
    ┌─────────▼────────┐  ┌────▼─────────────────┐
    │  ASV Feature     │  │  Representative      │
    │  Table           │  │  Sequences           │
    └─────────┬────────┘  └────┬─────────────────┘
              │                │
              │     ┌──────────▼──────────────┐
              │     │ 5. Taxonomy Assignment  │  SILVA 138.1 NR99
              │     │    (Naive Bayes)        │  Region-specific classifier
              │     └──────────┬──────────────┘
              │                │
    ┌─────────▼────────────────▼──────────┐
    │  6. Quality Filtering               │
    │     • Min frequency: 1000 reads     │
    │     • Exclude: mitochondria,        │
    │       chloroplast                   │
    └─────────┬───────────────────────────┘
              │
     ┌────────┴────────────┐
     │                     │
┌────▼──────────┐   ┌──────▼─────────────┐
│ 7. Phylogeny  │   │ 8. Diversity       │
│  MAFFT align  │──►│  Alpha: Shannon,   │
│  FastTree ML  │   │    Faith's PD      │
│  Midpoint     │   │  Beta: Bray-Curtis,│
│  rooting      │   │    UniFrac         │
└───────────────┘   └─────┬──────────────┘
                          │
              ┌───────────▼───────────┐
              │  9. Export            │
              │  BIOM + Newick + TSV  │
              │  → R / phyloseq       │
              └───────────────────────┘
```

## Computational Methods

### Data Preprocessing & Error Correction
Quality-informed parameter selection for paired-end read truncation (forward: 282 bp, reverse: 245 bp), determined from per-base quality score analysis. DADA2 error modeling resolves unique sequence variants from 3.6 GB of raw high-throughput sequencing data, with multi-threaded execution on the HPC cluster.

### Machine Learning-Based Classification
Taxonomic assignment against the SILVA 138.1 reference database (>400,000 entries) using a region-specific (341F–785R) Naive Bayes classifier trained on the target amplicon region. Classifier training is handled by `train_classifier.sh`.

### Phylogenetic Reconstruction
Multiple sequence alignment (MAFFT) and maximum-likelihood tree estimation (FastTree), executed with multi-threaded parallelization on the HPC cluster to compute evolutionary distance matrices.

### Statistical Diversity Analysis
Alpha diversity (Shannon, Faith's PD) and beta diversity (Bray-Curtis, UniFrac) with statistical significance testing (Kruskal-Wallis, PERMANOVA).

## Requirements

| Software     | Version  | Purpose                          |
|-------------|----------|----------------------------------|
| QIIME 2     | 2024.5   | Microbiome bioinformatics        |
| Micromamba  | latest   | Environment management           |
| DADA2       | (bundled)| Denoising & error correction     |
| Cutadapt    | (bundled)| Primer removal                   |
| MAFFT       | (bundled)| Multiple sequence alignment      |
| FastTree    | (bundled)| Phylogenetic tree construction   |
| RESCRIPt    | (bundled)| Reference database management    |
| SLURM       | —        | HPC job scheduling               |

## Repository Structure

```
.
├── README.md                  # This file
├── pipeline.sh                # Main analysis pipeline (Stages 1-9)
├── train_classifier.sh        # SILVA 138.1 classifier training (run once)
└── slurm_submit.sh            # SLURM job submission script
```

## Usage

### 1. Train the Classifier (one-time setup)

Request an interactive HPC session with sufficient memory:
```bash
srun --cpus-per-task=8 --mem=640G --time=2-00:00:00 --pty bash
```

Run the classifier training:
```bash
bash train_classifier.sh
```

### 2. Run the Full Pipeline

```bash
# Full pipeline
bash pipeline.sh

# Or run individual stages (1-9)
bash pipeline.sh 3   # Run only DADA2 denoising
bash pipeline.sh 9   # Run only data export
```

### 3. SLURM Batch Submission

```bash
sbatch slurm_submit.sh
```

## Key Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Forward truncation | 282 bp | Median quality > Q30 |
| Reverse truncation | 245 bp | Median quality > Q30 |
| Read overlap | 83 bp | 282 + 245 - 444 = 83 bp |
| Min sample frequency | 1000 reads | Exclude low-depth samples |
| Rarefaction depth | 7402 reads | Minimum across samples |
| Classifier | Naive Bayes | Region-specific (341F-785R) |
| Reference DB | SILVA 138.1 NR99 | >400k entries |

## Output Summary

| Output | Format | Description |
|--------|--------|-------------|
| Feature table | `.qza` / `.biom` | ASV abundance matrix |
| Taxonomy | `.qza` / `.tsv` | SILVA 138.1 classifications |
| Phylogenetic tree | `.qza` / `.nwk` | Rooted tree (MAFFT + FastTree) |
| Diversity metrics | `.qza` / `.qzv` | Alpha & beta diversity results |
| Export bundle | `.biom` + `.nwk` + `.tsv` | Ready for R/phyloseq import |

## Context

This pipeline was developed as part of a HiWi research position at the **Schneider Lab, Center for Regenerative Therapies Dresden (CRTD), TU Dresden**, investigating gut microbiome composition under bile acid-associated liver conditions.

# =============================================================================
# config.R — Shared Configuration and Utility Functions
# =============================================================================
#
# Computational Integration and High-Dimensional Correlation Analysis
# for Oral-Gut-Liver Axis Microbiome Dynamics
#
# This file defines shared paths, parameters, and helper functions used
# across all analysis modules.
#
# Author : Yunus Emre Kılıçkıran
# Lab    : Schneider Lab, CRTD, TU Dresden
# =============================================================================

# --- Required Packages -------------------------------------------------------
required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "pheatmap", "openxlsx",
  "phyloseq", "vegan", "ape", "patchwork", "ggpubr", "ggrepel",
  "reshape2", "tibble", "stringr", "ANCOMBC", "RColorBrewer",
  "tidyverse"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing: ", pkg)
    if (pkg %in% c("phyloseq", "ANCOMBC")) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# --- File Paths --------------------------------------------------------------
# Modify these paths to match your local or HPC environment.

DATA_PATH <- file.path(
  "data", "DataExport.xlsx"
)

METADATA_PATH <- file.path(
  "data", "MappingFile-v12.txt"
)

RNA_PATH <- file.path(
  "data", "TPM_normalized_data.csv"
)

MERGED_METADATA_PATH <- file.path(
  "data", "merged_metadata_corrected.csv"
)

RESULTS_DIR <- "results"

# Create results directory if it doesn't exist
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# --- Treatment Factor Levels -------------------------------------------------
TREATMENT_ORDER <- c(
  "WSD",
  "WSD_Ligatur",
  "WSD_AB",
  "WSD_Ligatur_AB",
  "WSD_PPI",
  "WSD_Ligatur_PPI"
)

BODYSITE_ORDER <- c("swab", "KOT", "liver")

# --- Treatment Color Palette -------------------------------------------------
TREATMENT_COLORS <- c(
  "WSD"              = "#1f77b4",
  "WSD_Ligatur"      = "#e377c2",
  "WSD_AB"           = "#9467bd",
  "WSD_Ligatur_AB"   = "#8c564b",
  "WSD_PPI"          = "#2ca02c",
  "WSD_Ligatur_PPI"  = "#ff7f0e"
)

BODYSITE_COLORS <- c(
  "swab"  = "gray",
  "KOT"   = "brown",
  "liver" = "purple"
)

# --- Heatmap Color Palette ---------------------------------------------------
HEATMAP_PALETTE <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(50)
HEATMAP_BREAKS  <- seq(-1, 1, length.out = 51)

# --- Helper Functions --------------------------------------------------------

#' Load and prepare the OTU abundance matrix from Excel
#'
#' @param sheet_name Sheet name in DataExport.xlsx (e.g., "Abund_Species_wide")
#' @param taxon_col  Name of the taxon column (e.g., "Species", "Phylum", "Genus")
#' @return A numeric matrix with taxa as rows and samples as columns
load_otu_matrix <- function(sheet_name, taxon_col) {
  otu_wide <- read_excel(DATA_PATH, sheet = sheet_name)
  otu_mat  <- as.matrix(otu_wide[, -1])
  rownames(otu_mat) <- otu_wide[[taxon_col]]
  colnames(otu_mat) <- gsub("^Abundance\\.", "", colnames(otu_mat))
  return(otu_mat)
}

#' Load and prepare sample metadata
#'
#' @return A data.frame with SampleID as rownames
load_metadata <- function() {
  meta <- read.delim(METADATA_PATH, header = TRUE, sep = "\t", check.names = FALSE)
  colnames(meta)[1] <- "SampleID"
  rownames(meta) <- meta$SampleID
  meta <- meta[meta$SampleID != "DNAcommunity", ]
  return(meta)
}

#' Load RNA-Seq TPM data
#'
#' @return A numeric matrix with genes as rows and animal IDs as columns
load_rna_matrix <- function() {
  rna <- read.csv(RNA_PATH, header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE, quote = "")
  colnames(rna)[1] <- "Gene_Name"
  rna$Gene_Name    <- make.unique(rna$Gene_Name)
  rownames(rna)    <- rna$Gene_Name
  rna$Gene_Name    <- NULL
  if ("Length_kb" %in% colnames(rna)) rna$Length_kb <- NULL
  colnames(rna) <- gsub("[[:punct:]]", "", colnames(rna))
  return(as.matrix(rna))
}

#' Custom pairwise PERMANOVA
#'
#' Performs pairwise PERMANOVA tests between all group pairs.
#'
#' @param dist_mat  A dist object
#' @param groups    A factor vector of group assignments
#' @param nperm     Number of permutations (default: 999)
#' @return A data.frame with columns: Group1, Group2, R2, p_value
pairwise_permanova <- function(dist_mat, groups, nperm = 999) {
  groups  <- as.factor(groups)
  combs   <- combn(levels(groups), 2, simplify = FALSE)
  results <- data.frame(
    Group1  = character(),
    Group2  = character(),
    R2      = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )

  for (pair in combs) {
    idx      <- groups %in% pair
    sub_dist <- as.dist(as.matrix(dist_mat)[idx, idx])
    sub_grp  <- droplevels(groups[idx])
    ad       <- adonis2(sub_dist ~ sub_grp, permutations = nperm)
    results  <- rbind(results, data.frame(
      Group1  = pair[1],
      Group2  = pair[2],
      R2      = round(ad$R2[1], 4),
      p_value = ad$`Pr(>F)`[1]
    ))
  }
  return(results)
}

cat("Configuration loaded successfully.\n")
cat("Results will be saved to:", RESULTS_DIR, "\n")

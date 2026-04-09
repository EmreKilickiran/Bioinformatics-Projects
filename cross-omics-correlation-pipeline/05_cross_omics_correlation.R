# 
# 05_cross_omics_correlation.R — 16S rRNA × RNA-Seq Correlation Pipeline
#
# Integrates 16S rRNA genus-level abundances with TPM-normalized RNA-Seq
# transcriptomic profiles to identify cross-domain gene–taxa associations.
#
# The pipeline implements rank-based Spearman correlation with manually
# computed test statistics (t-distribution conversion, two-tailed p-values)
# and Benjamini-Hochberg FDR correction, processing over 52 million
# gene–taxa correlation pairs through vectorized matrix operations.
#
# Outputs (per bodysite or bodysite×treatment):
#   - Full correlation matrix (all gene–taxa pairs)
#   - Significant correlations (FDR < 0.05, |rho| > 0.5)
#   - Heatmap of significant correlations
#   - Heatmap of top 100 genes by FDR-ranked correlation strength
#   - Gene lists and correlation matrices (Excel)


source("R/00_config.R")

cat("\n=== Module 5: Cross-Omics Correlation (16S × RNA-Seq) ===\n")

# --- Configuration -----------------------------------------------------------
# Analysis mode: "bodysite" runs per tissue site (swab, KOT, liver)
#                "bodysite_treatment" runs per tissue × treatment combination
ANALYSIS_MODE <- "bodysite_treatment"

# Groups to exclude from analysis (e.g., antibiotic treatment arms)
EXCLUDE_GROUPS <- c("ABNoLigature", "ABPlusLigature")

# Significance thresholds
FDR_THRESHOLD <- 0.05
RHO_THRESHOLD <- 0.5
TOP_N_GENES   <- 100

# --- Load Data ---------------------------------------------------------------
# 16S genus-level abundance matrix
otu <- read_excel(
  file.path(dirname(DATA_PATH), "DataExport.xlsx"),
  sheet = "Abund_Genus_wide"
)
colnames(otu)[1] <- "Taxonomy"
otu <- as.data.frame(otu)
colnames(otu) <- gsub("^Abundance\\.", "", colnames(otu))
rownames(otu) <- otu$Taxonomy
otu$Taxonomy <- NULL

# RNA-Seq TPM-normalized expression matrix
rna <- as.data.frame(load_rna_matrix())

# Merged metadata (contains both SampleID and AnimalID mapping)
meta <- read.csv(MERGED_METADATA_PATH, stringsAsFactors = FALSE)
meta$AnimalID <- as.character(meta$AnimalID)

# Exclude specified groups
if (length(EXCLUDE_GROUPS) > 0) {
  meta <- meta %>% filter(!Groups %in% EXCLUDE_GROUPS)
  cat("Excluded groups:", paste(EXCLUDE_GROUPS, collapse = ", "), "\n")
}

cat("16S taxa:", nrow(otu), "| RNA genes:", nrow(rna), "\n")


# Core Correlation Function

#' Run Spearman correlation between RNA-Seq genes and 16S taxa
#'
#' @param subset_label  Label for this analysis subset (used in filenames)
#' @param meta_subset   Filtered metadata data.frame
#' @param rna_data      Full RNA expression matrix
#' @param otu_data      Full 16S abundance matrix
#' @param output_dir    Base output directory
run_correlation <- function(subset_label, meta_subset, rna_data, otu_data,
                            output_dir) {

  cat("\n--------------------------------------------------\n")
  cat("Running correlation for:", subset_label, "\n")

  # Create output subdirectory
  sub_dir <- file.path(output_dir, subset_label)
  if (!dir.exists(sub_dir)) dir.create(sub_dir, recursive = TRUE)

  # --- Sample Matching (SampleID ↔ AnimalID) --------------------------------
  # 16S samples are identified by SampleID; RNA-Seq by AnimalID.
  # Metadata maps between the two identifiers.
  common_meta <- meta_subset %>%
    filter(AnimalID %in% colnames(rna_data),
           SampleID %in% colnames(otu_data))

  if (nrow(common_meta) < 3) {
    cat("  Skipping: fewer than 3 matched samples.\n")
    return(NULL)
  }

  samples_rna <- common_meta$AnimalID
  samples_otu <- common_meta$SampleID
  names(samples_rna) <- common_meta$SampleID

  rna_sub <- rna_data[, samples_rna, drop = FALSE]
  otu_sub <- otu_data[, samples_otu, drop = FALSE]

  n_samples <- ncol(rna_sub)
  n_pairs   <- nrow(rna_sub) * nrow(otu_sub)

  cat("  Matched samples:", n_samples, "\n")
  cat("  Gene–taxa pairs:", format(n_pairs, big.mark = ","), "\n")

  # --- Rank-Based Spearman Correlation (Vectorized) --------------------------
  # Rank-transform rows, then compute Pearson on ranks (= Spearman)
  rna_ranked <- t(apply(rna_sub, 1, rank, na.last = "keep"))
  otu_ranked <- t(apply(otu_sub, 1, rank, na.last = "keep"))

  cor_matrix <- cor(t(rna_ranked), t(otu_ranked),
                    method = "pearson", use = "pairwise.complete.obs")

  # --- P-value Computation (t-distribution) ----------------------------------
  t_stat <- cor_matrix * sqrt((n_samples - 2) / (1 - cor_matrix^2))
  p_matrix <- 2 * pt(-abs(t_stat), df = n_samples - 2)

  # --- Flatten to Long Format ------------------------------------------------
  results <- data.frame(
    Gene = rep(rownames(rna_sub), each  = nrow(otu_sub)),
    Taxa = rep(rownames(otu_sub), times = nrow(rna_sub)),
    rho  = as.vector(cor_matrix),
    pval = as.vector(p_matrix),
    stringsAsFactors = FALSE
  )

  # Benjamini-Hochberg FDR correction
  results$FDR <- p.adjust(results$pval, method = "fdr")

  # Save all correlations
  write.xlsx(results, file.path(sub_dir, "all_correlations.xlsx"),
             rowNames = FALSE)

  # Filter significant results
  sig_results <- results %>% filter(FDR < FDR_THRESHOLD, abs(rho) > RHO_THRESHOLD)
  write.xlsx(sig_results, file.path(sub_dir, "significant_correlations.xlsx"),
             rowNames = FALSE)

  cat("  Significant pairs:", nrow(sig_results), "\n")

  # --- Heatmap: Significant Correlations -------------------------------------
  if (nrow(sig_results) > 0) {
    cor_mat_sig <- sig_results %>%
      group_by(Taxa, Gene) %>%
      summarise(rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Gene, values_from = rho) %>%
      column_to_rownames("Taxa") %>%
      as.matrix()
    cor_mat_sig[is.na(cor_mat_sig)] <- 0

    pdf(file.path(sub_dir, "heatmap_significant.pdf"), width = 20, height = 12)
    pheatmap(
      cor_mat_sig,
      color = HEATMAP_PALETTE,
      breaks = HEATMAP_BREAKS,
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      clustering_method = "ward.D2",
      fontsize_row = 4,
      fontsize_col = 6,
      main = paste("Significant Gene–Taxa Correlations —", subset_label)
    )
    dev.off()
  }

  # --- Heatmap: Top N Genes by FDR -------------------------------------------
  top_genes <- results %>%
    arrange(FDR, desc(abs(rho))) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    slice(1:TOP_N_GENES) %>%
    pull(Gene)

  top_df <- results %>% filter(Gene %in% top_genes)

  cor_mat_top <- top_df %>%
    group_by(Taxa, Gene) %>%
    summarise(rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Gene, values_from = rho) %>%
    column_to_rownames("Taxa") %>%
    as.matrix()
  cor_mat_top[is.na(cor_mat_top)] <- 0

  # Sort rows by mean absolute correlation
  row_order <- order(rowMeans(abs(cor_mat_top), na.rm = TRUE), decreasing = TRUE)
  cor_mat_top <- cor_mat_top[row_order, , drop = FALSE]

  pdf(file.path(sub_dir, "heatmap_top100.pdf"), width = 20, height = 12)
  pheatmap(
    cor_mat_top,
    color = HEATMAP_PALETTE,
    breaks = HEATMAP_BREAKS,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    fontsize_row = 4,
    fontsize_col = 6,
    main = paste("Top", TOP_N_GENES, "Genes — RNA–Taxa Correlations —",
                 subset_label)
  )
  dev.off()

  # Save gene lists and matrix
  write.xlsx(data.frame(Gene = top_genes),
             file.path(sub_dir, "top100_gene_list.xlsx"), rowNames = FALSE)
  write.xlsx(top_df %>% arrange(FDR, desc(abs(rho))),
             file.path(sub_dir, "top100_pairs.xlsx"), rowNames = FALSE)
  write.xlsx(
    as.data.frame(cor_mat_top) %>% tibble::rownames_to_column("Taxa"),
    file.path(sub_dir, "top100_cor_matrix.xlsx"), rowNames = FALSE
  )

  cat("  Complete:", subset_label, "\n")
}


# Execute Analysis

if (ANALYSIS_MODE == "bodysite") {
  # Run per bodysite (swab, KOT, liver)
  for (bs in unique(meta$Bodysite)) {
    meta_sub <- meta %>% filter(Bodysite == bs)
    run_correlation(bs, meta_sub, rna, otu, RESULTS_DIR)
  }

} else if (ANALYSIS_MODE == "bodysite_treatment") {
  # Run per bodysite × treatment combination
  combos <- meta %>%
    mutate(label = paste0(Bodysite, "_", Treatment)) %>%
    pull(label) %>%
    unique()

  for (combo in combos) {
    parts     <- strsplit(combo, "_")[[1]]
    bodysite  <- parts[1]
    treatment <- paste(parts[-1], collapse = "_")

    meta_sub <- meta %>% filter(Bodysite == bodysite, Treatment == treatment)
    run_correlation(combo, meta_sub, rna, otu, RESULTS_DIR)
  }
}

cat("\n=== Cross-omics correlation analysis complete ===\n")

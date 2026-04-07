# =============================================================================
# High-Dimensional Data Analysis and Statistical Modeling of Gene Regulatory
# Networks in Liver Cancer
# =============================================================================
#
# Description:
#   A computational pipeline to quantify transcription factor (TF) activity
#   from high-dimensional genomic datasets and identify statistically
#   significant regulatory patterns associated with tumor progression in
#   hepatocellular carcinoma (TCGA-LIHC cohort).
#
#
# Pipeline Overview:
#   1. Load TCGA-LIHC gene expression data and sample metadata
#   2. Load ChIP-Atlas regulon data for each TF (ChIP-seq binding targets)
#   3. Construct VIPER-compatible regulon objects
#   4. Compute TF activity scores using VIPER algorithm
#   5. Merge activity scores with biomarker expression and clinical class
#   6. Perform correlation analysis (Spearman) stratified by sample type
#   7. Generate annotated scatter plots (TF activity vs AFP expression)
#   8. Export all results
#
#
# Data Sources:
#   - Gene expression: TCGA-LIHC (The Cancer Genome Atlas, Liver Hepatocellular
#     Carcinoma), RNA-Seq normalized counts
#   - TF binding targets: ChIP-Atlas (https://chip-atlas.org), liver-specific
#     ChIP-seq experiments in HepG2 cells
#
# Author : Yunus Emre Kılıçkıran
# =============================================================================


# =============================================================================
# 1. SETUP & CONFIGURATION
# =============================================================================

# --- Install and load required packages --------------------------------------
required_cran <- c("tidyverse", "openxlsx", "ggplot2")
required_bioc <- c("viper")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

library(tidyverse)
library(openxlsx)
library(ggplot2)
library(viper)

# --- File paths --------------------------------------------------------------
# Modify these paths to match your local environment.

DATA_DIR    <- "data"
RESULTS_DIR <- "results"

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# Input files
EXPRESSION_FILE <- file.path(DATA_DIR, "LIHC_expression_data_clean.xlsx")

REGULON_FILES <- list(
  NFATC3 = file.path(DATA_DIR, "NFATC3.Liver.tsv"),
  NFAT5  = file.path(DATA_DIR, "NFAT5.Liver.tsv")
)

# --- Analysis parameters -----------------------------------------------------
BIOMARKER_GENE <- "AFP"

# TCGA sample type classification:
# Barcodes ending in "-01" = Primary Tumor, "-11" = Solid Tissue Normal
TUMOR_SUFFIX  <- "-01"
NORMAL_SUFFIX <- "-11"

# Plot aesthetics
SAMPLE_COLORS <- c("Normal" = "#4575B4", "Tumor" = "#D73027")


# =============================================================================
# 2. LOAD GENE EXPRESSION DATA
# =============================================================================

cat("=============================================================\n")
cat(" TF Activity Analysis Pipeline — TCGA-LIHC\n")
cat("=============================================================\n\n")

cat("[Step 1] Loading TCGA-LIHC expression matrix...\n")

expr <- read.xlsx(EXPRESSION_FILE, rowNames = TRUE)
expr_mat <- as.matrix(expr)

cat("  Dimensions:", nrow(expr_mat), "genes ×", ncol(expr_mat), "samples\n")
cat("  Total data points:", format(prod(dim(expr_mat)), big.mark = ","), "\n")

# --- Classify samples as Tumor or Normal based on TCGA barcode ---------------
sample_type <- ifelse(
  grepl(paste0(TUMOR_SUFFIX, "$"), colnames(expr_mat)),
  "Tumor", "Normal"
)
names(sample_type) <- colnames(expr_mat)

cat("  Tumor samples: ", sum(sample_type == "Tumor"), "\n")
cat("  Normal samples:", sum(sample_type == "Normal"), "\n")

# --- Extract biomarker expression --------------------------------------------
if (!BIOMARKER_GENE %in% rownames(expr_mat)) {
  stop("Biomarker gene '", BIOMARKER_GENE, "' not found in expression matrix.")
}

biomarker_expr <- expr_mat[BIOMARKER_GENE, ]
cat("  Biomarker:", BIOMARKER_GENE, "expression loaded.\n\n")


# =============================================================================
# 3. REGULON CONSTRUCTION & VIPER ACTIVITY SCORING
# =============================================================================

#' Build a VIPER regulon object from ChIP-Atlas binding data
#'
#' Constructs a regulon in the format required by the VIPER algorithm:
#'   - tfmode:     Named vector of binding weights (Average score from ChIP-seq)
#'   - likelihood: Named vector of absolute binding weights
#'   - synergy:    Zero vector (no synergy modeled)
#'
#' @param regulon_path Path to ChIP-Atlas TSV file
#' @param tf_name      Name of the transcription factor
#' @return A named list with one element (the TF regulon)
build_regulon <- function(regulon_path, tf_name) {
  cat("[Step 2] Building regulon for:", tf_name, "\n")

  regulon_df <- read_tsv(regulon_path, show_col_types = FALSE)

  cat("  Target genes:", nrow(regulon_df), "\n")
  cat("  Source: ChIP-Atlas, HepG2 liver cells\n")

  # Construct VIPER regulon object
  regulon <- list()
  regulon[["tfmode"]]     <- setNames(regulon_df$Average, regulon_df$Protein)
  regulon[["likelihood"]] <- setNames(abs(regulon_df$Average), regulon_df$Protein)
  regulon[["synergy"]]    <- setNames(rep(0, nrow(regulon_df)), regulon_df$Protein)

  regulon_list <- list(regulon)
  names(regulon_list) <- tf_name

  return(regulon_list)
}


#' Compute TF activity scores using VIPER
#'
#' VIPER estimates TF activity by evaluating the coordinated expression
#' of its ChIP-seq-derived target genes across all samples.
#'
#' @param expr_mat     Gene expression matrix (genes × samples)
#' @param regulon_list VIPER regulon object
#' @param tf_name      Name of the transcription factor
#' @return A named numeric vector of activity scores (one per sample)
compute_viper_activity <- function(expr_mat, regulon_list, tf_name) {
  cat("[Step 3] Computing VIPER activity scores for:", tf_name, "\n")

  activity <- viper(expr_mat, regulon_list, method = "none")

  # Extract the TF row and transpose to a named vector
  activity_scores <- as.numeric(activity[tf_name, ])
  names(activity_scores) <- colnames(activity)

  cat("  Activity scores computed for", length(activity_scores), "samples.\n")
  cat("  Range: [", round(min(activity_scores), 2), ",",
      round(max(activity_scores), 2), "]\n\n")

  return(activity_scores)
}


# =============================================================================
# 4. CORRELATION ANALYSIS & VISUALIZATION
# =============================================================================

#' Perform correlation analysis and generate scatter plot
#'
#' Computes Spearman correlation between TF activity and biomarker expression,
#' both globally and stratified by sample type (Tumor/Normal). Generates an
#' annotated scatter plot.
#'
#' @param activity_scores Named vector of TF activity scores
#' @param biomarker_expr  Named vector of biomarker expression values
#' @param sample_type     Named vector of sample classifications
#' @param tf_name         Name of the transcription factor
#' @param biomarker_name  Name of the biomarker gene
#' @param output_dir      Directory for output files
analyze_and_plot <- function(activity_scores, biomarker_expr, sample_type,
                             tf_name, biomarker_name, output_dir) {

  cat("[Step 4] Correlation analysis:", tf_name, "vs", biomarker_name, "\n")

  # --- Assemble analysis data frame ------------------------------------------
  common_samples <- intersect(names(activity_scores), names(biomarker_expr))

  plot_df <- data.frame(
    Sample     = common_samples,
    Activity   = activity_scores[common_samples],
    Biomarker  = biomarker_expr[common_samples],
    SampleType = sample_type[common_samples],
    stringsAsFactors = FALSE
  )

  # --- Spearman Correlation (global and stratified) --------------------------
  cor_global <- cor.test(plot_df$Activity, plot_df$Biomarker, method = "spearman")
  cor_tumor  <- cor.test(
    plot_df$Activity[plot_df$SampleType == "Tumor"],
    plot_df$Biomarker[plot_df$SampleType == "Tumor"],
    method = "spearman"
  )
  cor_normal <- cor.test(
    plot_df$Activity[plot_df$SampleType == "Normal"],
    plot_df$Biomarker[plot_df$SampleType == "Normal"],
    method = "spearman"
  )

  # --- Print statistics ------------------------------------------------------
  cat("  Global:  rho =", round(cor_global$estimate, 4),
      " p =", signif(cor_global$p.value, 3), "\n")
  cat("  Tumor:   rho =", round(cor_tumor$estimate, 4),
      " p =", signif(cor_tumor$p.value, 3), "\n")
  cat("  Normal:  rho =", round(cor_normal$estimate, 4),
      " p =", signif(cor_normal$p.value, 3), "\n")

  # --- Save statistics to Excel ----------------------------------------------
  stats_df <- data.frame(
    TF        = tf_name,
    Biomarker = biomarker_name,
    Subset    = c("Global", "Tumor", "Normal"),
    N         = c(nrow(plot_df),
                  sum(plot_df$SampleType == "Tumor"),
                  sum(plot_df$SampleType == "Normal")),
    Spearman_rho = c(cor_global$estimate, cor_tumor$estimate, cor_normal$estimate),
    p_value      = c(cor_global$p.value, cor_tumor$p.value, cor_normal$p.value)
  )

  write.xlsx(stats_df,
             file.path(output_dir, paste0(tf_name, "_correlation_stats.xlsx")),
             rowNames = FALSE)

  # --- Save activity scores --------------------------------------------------
  write.xlsx(plot_df,
             file.path(output_dir, paste0(tf_name, "_activity_scores.xlsx")),
             rowNames = FALSE)

  # --- Scatter Plot: TF Activity vs Biomarker Expression ---------------------
  subtitle_text <- sprintf(
    "Spearman rho = %.3f (p = %.2e) | N = %d (Tumor: %d, Normal: %d)",
    cor_global$estimate, cor_global$p.value,
    nrow(plot_df),
    sum(plot_df$SampleType == "Tumor"),
    sum(plot_df$SampleType == "Normal")
  )

  p <- ggplot(plot_df, aes(x = Activity, y = Biomarker,
                           color = SampleType, shape = SampleType)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = SAMPLE_COLORS, name = "Sample Type") +
    scale_shape_manual(values = c("Normal" = 18, "Tumor" = 16),
                       name = "Sample Type") +
    labs(
      title    = paste(tf_name, "Activity vs", biomarker_name, "Expression"),
      subtitle = subtitle_text,
      x        = paste(tf_name, "Activity (VIPER Score)"),
      y        = paste(biomarker_name, "Expression")
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
      legend.position = "right"
    )

  ggsave(
    file.path(output_dir, paste0(tf_name, "_vs_", biomarker_name, ".pdf")),
    plot = p, width = 8, height = 6
  )

  cat("  Plot saved.\n\n")

  return(list(plot = p, stats = stats_df, data = plot_df))
}


# =============================================================================
# 5. MAIN EXECUTION — Run for Each Transcription Factor
# =============================================================================

all_results <- list()

for (tf_name in names(REGULON_FILES)) {

  cat("=============================================================\n")
  cat(" Processing:", tf_name, "\n")
  cat("=============================================================\n\n")

  # Build regulon from ChIP-Atlas data
  regulon_list <- build_regulon(REGULON_FILES[[tf_name]], tf_name)

  # Compute VIPER activity scores
  activity_scores <- compute_viper_activity(expr_mat, regulon_list, tf_name)

  # Correlation analysis + visualization
  result <- analyze_and_plot(
    activity_scores  = activity_scores,
    biomarker_expr   = biomarker_expr,
    sample_type      = sample_type,
    tf_name          = tf_name,
    biomarker_name   = BIOMARKER_GENE,
    output_dir       = RESULTS_DIR
  )

  all_results[[tf_name]] <- result
}


# =============================================================================
# 6. COMBINED SUMMARY
# =============================================================================

cat("=============================================================\n")
cat(" Summary\n")
cat("=============================================================\n\n")

# Combine all statistics into one table
summary_stats <- bind_rows(lapply(all_results, function(r) r$stats))
write.xlsx(summary_stats,
           file.path(RESULTS_DIR, "combined_correlation_summary.xlsx"),
           rowNames = FALSE)

cat("Combined statistics saved.\n\n")

# Print summary table
print(summary_stats)

cat("\n=============================================================\n")
cat(" Pipeline complete. Results saved to:", RESULTS_DIR, "\n")
cat("=============================================================\n")

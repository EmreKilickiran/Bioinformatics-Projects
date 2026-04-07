# =============================================================================
# 01_preprocessing.R — Quality Control, Filtering, Normalization & Scaling
# =============================================================================
#
# Processes raw count matrices for 310,000+ single cells through:
#   1. QC metric computation (nFeature, nCount, percent.mt)
#   2. QC visualization (violin plots, feature scatter)
#   3. Custom filtering: percent.mt <= 15%, residual-based doublet filtering
#      at 3×MAD from the nCount–nFeature linear trend
#   4. Log-normalization (scale factor = 10,000)
#   5. Highly variable gene (HVG) identification (VST, top 3,000)
#   6. Scaling (zero-mean, unit-variance on HVGs)
#
# Input:  raw.rds (raw Seurat/count object)
# Output: seurat_scaled.rds (QC-filtered, normalized, HVG-identified, scaled)
#
# Author : Yunus Emre Kılıçkıran
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# --- Configuration -----------------------------------------------------------
source("R/00_config.R")

cat("\n=== Module 1: Preprocessing (QC → Scaling) ===\n")

# --- QC Parameters -----------------------------------------------------------
MT_MAX        <- 15      # Maximum mitochondrial percentage
RESID_SIGMA   <- 3       # Trend residual threshold (MAD units)
NORM_METHOD   <- "LogNormalize"
SCALE_FACTOR  <- 10000
N_HVG         <- 3000    # Number of highly variable genes
HVG_METHOD    <- "vst"
TOP_HVG_LABEL <- 20      # Number of top HVGs to label in plot

# =============================================================================
# STEP 1: Load Raw Data & Initialize Seurat Object
# =============================================================================

cat("[1/6] Loading raw data...\n")

sc_data <- readRDS(RAW_DATA_FILE)
DefaultAssay(sc_data) <- "RNA"
counts <- GetAssayData(sc_data, layer = "counts")
colnames(counts) <- make.unique(colnames(counts))
seu <- CreateSeuratObject(counts = counts)

cat("  Cells:", ncol(seu), "| Genes:", nrow(seu), "\n")

# Mitochondrial gene percentage
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-|^mt-")

# =============================================================================
# STEP 2: QC Visualization
# =============================================================================

cat("[2/6] Generating QC plots...\n")

Idents(seu) <- factor(rep("Sample", ncol(seu)))

# Violin plots
png(file.path(RESULTS_DIR, "01_qc_violin.png"), width = 2400, height = 800, res = 300)
print(VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0))
dev.off()

# Feature scatter (nCount vs nFeature)
png(file.path(RESULTS_DIR, "01_feature_scatter.png"), width = 1200, height = 1000, res = 150)
print(
  FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    ggtitle("nCount_RNA vs nFeature_RNA") +
    theme_minimal()
)
dev.off()

# =============================================================================
# STEP 3: Quality Filtering
# =============================================================================

cat("[3/6] Filtering cells...\n")

# Residual-based filtering: fit linear trend in log10 space,
# remove cells deviating more than 3×MAD from the trend
df_qc <- data.frame(
  log_counts   = log10(seu$nCount_RNA + 1),
  log_features = log10(seu$nFeature_RNA + 1)
)
fit <- lm(log_features ~ log_counts, data = df_qc)
resid_vec <- residuals(fit)
mad_val   <- mad(resid_vec, center = 0)

keep_resid <- abs(resid_vec) <= RESID_SIGMA * mad_val
keep_mt    <- seu$percent.mt <= MT_MAX
keep_cells <- keep_resid & keep_mt

# QC scatter: kept vs removed
png(file.path(RESULTS_DIR, "01_qc_filter_scatter.png"), width = 1600, height = 1400, res = 180)
keep_factor <- ifelse(keep_cells, "Kept", "Removed")
p <- ggplot(data.frame(nCount = seu$nCount_RNA, nFeat = seu$nFeature_RNA,
                       Status = keep_factor),
            aes(x = nCount, y = nFeat, color = Status)) +
  geom_point(size = 0.4, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Kept" = "#2166AC", "Removed" = "#D73027")) +
  labs(x = "nCount_RNA", y = "nFeature_RNA",
       title = sprintf("QC Filtering: %d → %d cells (%.1f%% kept)",
                       ncol(seu), sum(keep_cells),
                       100 * sum(keep_cells) / ncol(seu))) +
  theme_minimal()
print(p)
dev.off()

seu <- subset(seu, cells = colnames(seu)[keep_cells])

# Save filter summary
sink(file.path(RESULTS_DIR, "01_qc_filter_summary.txt"))
cat("Cells before:", length(keep_cells), "\n")
cat("Cells after :", ncol(seu), "\n")
cat("Kept ratio  :", round(ncol(seu) / length(keep_cells), 4), "\n")
cat("MT% cutoff  :", MT_MAX, "\n")
cat("|residual| <=", RESID_SIGMA, "× MAD (MAD =", round(mad_val, 4), ")\n")
sink()

cat("  Cells after filtering:", ncol(seu), "\n")

# =============================================================================
# STEP 4: Log-Normalization
# =============================================================================

cat("[4/6] Normalizing (LogNormalize, scale.factor =", SCALE_FACTOR, ")...\n")

seu <- NormalizeData(seu, normalization.method = NORM_METHOD,
                     scale.factor = SCALE_FACTOR, verbose = FALSE)

# =============================================================================
# STEP 5: Highly Variable Gene Identification
# =============================================================================

cat("[5/6] Identifying", N_HVG, "highly variable genes (VST)...\n")

seu <- FindVariableFeatures(seu, selection.method = HVG_METHOD,
                            nfeatures = N_HVG, verbose = FALSE)

top_genes <- head(VariableFeatures(seu), TOP_HVG_LABEL)

png(file.path(RESULTS_DIR, "01_variable_features_top20.png"),
    width = 1600, height = 1200, res = 200)
plot1 <- VariableFeaturePlot(seu)
print(LabelPoints(plot = plot1, points = top_genes, repel = TRUE))
dev.off()

cat("  Top HVGs:", paste(head(top_genes, 5), collapse = ", "), "...\n")

# =============================================================================
# STEP 6: Scaling
# =============================================================================

cat("[6/6] Scaling HVGs...\n")

seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)

# --- Save preprocessed object ------------------------------------------------
saveRDS(seu, file = file.path(RESULTS_DIR, "seurat_preprocessed.rds"))
cat("Preprocessing complete. Saved:", file.path(RESULTS_DIR, "seurat_preprocessed.rds"), "\n")

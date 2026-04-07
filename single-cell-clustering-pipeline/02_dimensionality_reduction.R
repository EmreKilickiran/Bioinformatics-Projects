# =============================================================================
# 02_dimensionality_reduction.R — PCA, Clustering & UMAP
# =============================================================================
#
# Performs dimensionality reduction and unsupervised clustering:
#   1. PCA on highly variable genes (50 PCs computed, 25 used)
#   2. Elbow plot for PC selection
#   3. SNN graph construction + multi-resolution clustering (0.1–1.0)
#   4. Final clustering at optimal resolution (0.15)
#   5. UMAP embedding (20 PCs)
#
# Input:  seurat_preprocessed.rds
# Output: seurat_clustered.rds (with PCA, UMAP, cluster assignments)
#
# Author : Yunus Emre Kılıçkıran
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("R/00_config.R")

cat("\n=== Module 2: Dimensionality Reduction & Clustering ===\n")

# --- Parameters --------------------------------------------------------------
N_PCS_COMPUTE  <- 50     # Total PCs to compute
N_PCS_USE      <- 20     # PCs used for neighbors/UMAP
FINAL_RESOLUTION <- 0.15 # Selected after resolution sweep

# Resolution sweep range for evaluation
RESOLUTIONS <- c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0)

# =============================================================================
# STEP 1: PCA
# =============================================================================

cat("[1/4] Running PCA (", N_PCS_COMPUTE, "components)...\n")

seu <- readRDS(file.path(RESULTS_DIR, "seurat_preprocessed.rds"))
seu <- RunPCA(seu, features = VariableFeatures(seu),
              npcs = N_PCS_COMPUTE, verbose = FALSE)

# Print top gene loadings for first 5 PCs
print(seu[["pca"]], dims = 1:5, nfeatures = 5)

# PCA heatmaps for first 5 PCs
for (pc in 1:5) {
  png(file.path(RESULTS_DIR, sprintf("02_pca_heatmap_PC%d.png", pc)),
      width = 1600, height = 1200, res = 180)
  print(DimHeatmap(seu, dims = pc, cells = 500, balanced = TRUE, fast = TRUE))
  dev.off()
}

# Elbow plot for PC selection
png(file.path(RESULTS_DIR, "02_elbow_plot.png"), width = 1200, height = 900, res = 150)
print(ElbowPlot(seu, ndims = N_PCS_COMPUTE))
dev.off()

cat("  Elbow plot saved. Inspect to verify PC selection.\n")

# =============================================================================
# STEP 2: Resolution Sweep
# =============================================================================

cat("[2/4] Testing clustering resolutions:", paste(RESOLUTIONS, collapse = ", "), "\n")

seu <- FindNeighbors(seu, dims = 1:25, verbose = FALSE)

for (res in RESOLUTIONS) {
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
  meta_col <- paste0(DefaultAssay(seu), "_snn_res.", res)

  png(file.path(RESULTS_DIR, sprintf("02_clusters_res_%s.png",
                                      gsub("\\.", "_", as.character(res)))),
      width = 1600, height = 1200, res = 180)
  print(DimPlot(seu, reduction = "pca", group.by = meta_col, label = TRUE))
  dev.off()
}

# =============================================================================
# STEP 3: Final Clustering
# =============================================================================

cat("[3/4] Final clustering at resolution =", FINAL_RESOLUTION, "\n")

# Re-run PCA with selected number of PCs, build neighbors, cluster
seu <- RunPCA(seu, features = VariableFeatures(seu),
              npcs = N_PCS_USE, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:N_PCS_USE, verbose = FALSE)
seu <- FindClusters(seu, resolution = FINAL_RESOLUTION, verbose = FALSE)
Idents(seu) <- paste0("RNA_snn_res.", FINAL_RESOLUTION)

n_clusters <- length(unique(Idents(seu)))
cat("  Clusters found:", n_clusters, "\n")

# =============================================================================
# STEP 4: UMAP
# =============================================================================

cat("[4/4] Computing UMAP embedding (", N_PCS_USE, "PCs)...\n")

seu <- RunUMAP(seu, dims = 1:N_PCS_USE, verbose = FALSE)

png(file.path(RESULTS_DIR, "02_umap.png"), width = 1600, height = 1200, res = 180)
print(DimPlot(seu, reduction = "umap", label = TRUE))
dev.off()

# --- Save --------------------------------------------------------------------
saveRDS(seu, file = file.path(RESULTS_DIR, "seurat_clustered.rds"))
cat("Clustering complete. Saved:", file.path(RESULTS_DIR, "seurat_clustered.rds"), "\n")

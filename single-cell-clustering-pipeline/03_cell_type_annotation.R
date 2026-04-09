# 03_cell_type_annotation.R — Automated Cell Type Annotation (SingleR)
#
# Annotates 26 clusters into 11 cell types using SingleR with pseudobulk
# aggregation against the MacParland et al. human liver atlas (GSE115469).
#
# Method:
#   1. Aggregate expression per cluster (pseudobulk)
#   2. Run SingleR against liver-specific reference
#   3. Map cluster IDs → cell type labels
#   4. Generate annotated UMAP
#
# Input:  seurat_clustered.rds + GSE115469 reference files
# Output: seurat_annotated.rds (with cell_type metadata)


suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(BiocParallel)
  library(Matrix)
})

source("R/00_config.R")

cat("\n=== Module 3: Cell Type Annotation (SingleR) ===\n")

# --- Reference Data Paths ----------------------------------------------------
REFERENCE_EXPR <- file.path(DATA_DIR, "GSE115469_Data.csv.gz")
REFERENCE_META <- file.path(DATA_DIR, "GSE115469_CellClusterType.txt.gz")


# STEP 1: Load Clustered Object

cat("[1/4] Loading clustered object...\n")

obj <- readRDS(file.path(RESULTS_DIR, "seurat_clustered.rds"))
DefaultAssay(obj) <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"

if (is.null(Idents(obj)) || length(unique(Idents(obj))) <= 1) {
  if (!is.null(obj$seurat_clusters)) Idents(obj) <- obj$seurat_clusters
}

# Ensure variable features exist
if (length(VariableFeatures(obj)) == 0) {
  obj <- FindVariableFeatures(obj, nfeatures = 2000, selection.method = "vst")
}
var_genes <- VariableFeatures(obj)

cat("  Clusters:", length(unique(Idents(obj))), "\n")


# STEP 2: Build Liver Reference (MacParland et al.)

cat("[2/4] Building liver reference atlas...\n")

expr_ref <- as.matrix(read.csv(REFERENCE_EXPR, row.names = 1, check.names = FALSE))
meta_ref <- read.delim(REFERENCE_META, header = TRUE, stringsAsFactors = FALSE)
rownames(meta_ref) <- meta_ref$CellName
meta_ref <- meta_ref[colnames(expr_ref), , drop = FALSE]

sce_ref <- SingleCellExperiment(
  assays  = list(logcounts = log1p(expr_ref)),
  colData = meta_ref
)

cat("  Reference cells:", ncol(sce_ref), "\n")


# STEP 3: Pseudobulk SingleR

cat("[3/4] Running SingleR (pseudobulk per cluster)...\n")

# Common genes between query and reference
common_genes <- intersect(var_genes, rownames(sce_ref))
if (length(common_genes) < 500) {
  common_genes <- intersect(rownames(obj[[DefaultAssay(obj)]]), rownames(sce_ref))
}
cat("  Common genes:", length(common_genes), "\n")

# Pseudobulk aggregation
avg_list <- AggregateExpression(obj, assays = DefaultAssay(obj),
                                features = common_genes,
                                group.by = "seurat_clusters", slot = "data")
avg_mat  <- as.matrix(avg_list[[1]])

# Subset reference to matching genes
sce_ref_sub <- sce_ref[rownames(sce_ref) %in% rownames(avg_mat), ]
avg_mat_sub <- avg_mat[rownames(sce_ref_sub), , drop = FALSE]

# Run SingleR
sce_pb <- SingleCellExperiment(assays = list(logcounts = avg_mat_sub))

ncores <- max(1, min(4, parallel::detectCores()))
bp <- SnowParam(workers = ncores, type = "SOCK", progressbar = TRUE)

pred <- SingleR(
  test      = sce_pb,
  ref       = sce_ref_sub,
  labels    = sce_ref_sub$CellType,
  de.method = "wilcox",
  fine.tune = TRUE,
  prune     = TRUE,
  BPPARAM   = bp
)


# STEP 4: Apply Annotations

cat("[4/4] Applying cell type annotations...\n")

# Build cluster → cell type mapping
cluster_to_celltype <- data.frame(
  cluster   = colnames(avg_mat_sub),
  cell_type = pred$labels,
  stringsAsFactors = FALSE
)

# Sort by cluster number
ord <- order(as.numeric(gsub("[^0-9]+", "", cluster_to_celltype$cluster)))
cluster_to_celltype <- cluster_to_celltype[ord, ]

cat("\n  Cluster → Cell Type Mapping:\n")
print(cluster_to_celltype, row.names = FALSE)

# Apply to Seurat object
lut <- setNames(cluster_to_celltype$cell_type, cluster_to_celltype$cluster)
cl  <- as.character(Idents(obj))
Idents(obj) <- cl
obj <- RenameIdents(obj, lut)
obj$cell_type <- as.character(Idents(obj))

# Save mapping table
write.table(cluster_to_celltype,
            file.path(RESULTS_DIR, "03_cluster_celltype_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Annotated UMAP
png(file.path(RESULTS_DIR, "03_umap_annotated.png"),
    width = 1800, height = 1200, res = 180)
print(DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
        ggtitle("UMAP — Cell Type Annotation (SingleR)"))
dev.off()

# --- Save --------------------------------------------------------------------
saveRDS(obj, file = file.path(RESULTS_DIR, "seurat_annotated.rds"))
cat("\nAnnotation complete. Saved:", file.path(RESULTS_DIR, "seurat_annotated.rds"), "\n")

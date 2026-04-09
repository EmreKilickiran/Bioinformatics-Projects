# 04_differential_expression.R — DE Analysis Across Disease Stages
#
# Performs differential expression testing for each cell type across six
# sequential liver disease stages (NORM → HO → STEA → NASH → FIB → CIRR),
# using NORM as the reference baseline.
#
# Computes a custom directional expression ratio score:
#   Score = (Up − Down) / (Up + Down)
# to quantify the directional balance of transcriptional shifts per cell type,
# visualized as a hierarchically clustered heatmap.
#
# Also generates per-cluster marker gene dot plots across treatment groups.
#
# Input:  seurat_annotated.rds (with cell_type and Treatment metadata)
# Output: DE gene tables, ratio heatmaps, marker dot plots


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(future)
})

source("R/00_config.R")

cat("\n=== Module 4: Differential Expression Analysis ===\n")

# --- Parameters --------------------------------------------------------------
DISEASE_STAGES <- c("NORM", "HO", "STEA", "NASH", "FIB", "CIRR")
REFERENCE      <- "NORM"
LFC_THRESHOLD  <- 0.1
MIN_PCT        <- 0.05
TEST_METHOD    <- "LR"     # Likelihood ratio test
PADJ_CUTOFF    <- 0.05

# Parallel workers
plan("multicore", workers = min(6, parallel::detectCores()))


# STEP 1: Load Annotated Object

cat("[1/4] Loading annotated object...\n")

obj <- readRDS(file.path(RESULTS_DIR, "seurat_annotated.rds"))
stopifnot("Treatment" %in% colnames(obj@meta.data))
obj$group <- as.character(obj@meta.data$Treatment)
DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else "RNA"

cell_types <- sort(unique(obj$cell_type))
cat("  Cell types:", length(cell_types), "\n")
cat("  Disease stages:", paste(DISEASE_STAGES, collapse = " → "), "\n")


# STEP 2: DE Testing (Each Cell Type × Disease Stage vs NORM)

cat("[2/4] Running DE tests (", TEST_METHOD, ", ref =", REFERENCE, ")...\n")

out_tsv <- file.path(RESULTS_DIR, "04_DE_genes_vsNORM.tsv")
if (file.exists(out_tsv)) file.remove(out_tsv)

total_jobs <- length(cell_types) * (length(DISEASE_STAGES) - 1)
pb <- txtProgressBar(min = 0, max = total_jobs, style = 3)
job_i <- 0

for (ct in cell_types) {
  sub_ct  <- subset(obj, subset = cell_type == ct)
  g_avail <- intersect(DISEASE_STAGES, unique(sub_ct$group))
  if (!REFERENCE %in% g_avail) { job_i <- job_i + length(DISEASE_STAGES) - 1; next }

  for (g in setdiff(g_avail, REFERENCE)) {
    sub_g <- subset(sub_ct, subset = group %in% c(REFERENCE, g))
    Idents(sub_g) <- "group"

    de <- tryCatch(
      FindMarkers(sub_g, ident.1 = g, ident.2 = REFERENCE,
                  logfc.threshold = LFC_THRESHOLD, min.pct = MIN_PCT,
                  test.use = TEST_METHOD, only.pos = FALSE),
      error = function(e) NULL
    )

    if (!is.null(de) && nrow(de) > 0) {
      de$gene       <- rownames(de)
      de$cell_type  <- ct
      de$case_group <- g
      de$ref_group  <- REFERENCE

      fwrite(de, file = out_tsv, sep = "\t", quote = FALSE,
             row.names = FALSE,
             col.names = !file.exists(out_tsv),
             append = TRUE)
    }

    job_i <- job_i + 1
    setTxtProgressBar(pb, job_i)
  }
}
close(pb)

cat("\n  DE results saved:", out_tsv, "\n")


# STEP 3: DE Ratio Heatmap

cat("[3/4] Computing DE ratio heatmap...\n")

dt <- fread(out_tsv)
dt_sig <- dt[p_val_adj < PADJ_CUTOFF]

scores <- dt_sig[, .(
  up   = sum(avg_log2FC > 0),
  down = sum(avg_log2FC < 0)
), by = .(cell_type, case_group, ref_group)]

scores[, score := ifelse((up + down) > 0, (up - down) / (up + down), 0)]

# Build matrix
mat <- dcast(scores, cell_type ~ case_group, value.var = "score", fill = 0)
rn <- mat$cell_type
mat$cell_type <- NULL
mat <- as.data.frame(mat)
rownames(mat) <- rn

# Add NORM column (reference = 0)
mat$NORM <- 0
keep_cols <- intersect(DISEASE_STAGES, colnames(mat))
mat <- mat[, keep_cols, drop = FALSE]

# Heatmap
breaks <- seq(-1, 1, length.out = 101)
cols   <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(length(breaks) - 1)

png(file.path(RESULTS_DIR, "04_DE_ratio_heatmap.png"),
    width = 1800, height = 1200, res = 200)
pheatmap(mat,
         color = cols, breaks = breaks,
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = paste("DE Ratio per Cell Type (vs", REFERENCE, ")\n",
                      "Score = (Up-Down)/(Up+Down) |",
                      paste(DISEASE_STAGES, collapse = " → ")),
         border_color = NA, fontsize_row = 10, fontsize_col = 11)
dev.off()

write.table(as.data.frame(mat),
            file.path(RESULTS_DIR, "04_DE_ratio_matrix.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)

cat("  Heatmap saved.\n")


# STEP 4: Marker Gene Dot Plots by Treatment

cat("[4/4] Generating marker gene dot plots...\n")

# Find top 3 markers per cluster
all_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = 0.25, verbose = FALSE)

top3 <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  ungroup()

fwrite(as.data.table(top3),
       file.path(RESULTS_DIR, "04_markers_top3.tsv"), sep = "\t")

# Generate dot plot for top-1 marker per cluster
clusters <- sort(unique(top3$cluster))
all_pd <- list()

for (cl in clusters) {
  g <- top3 %>% filter(cluster == cl, rank == 1) %>% pull(gene)
  if (length(g) == 0) next

  sub_obj <- subset(obj, idents = cl)
  dp <- DotPlot(sub_obj, features = g, group.by = "Treatment")
  pd <- dp$data
  pd$cluster <- cl
  all_pd[[as.character(cl)]] <- pd
}

if (length(all_pd) > 0) {
  pd_all <- rbindlist(all_pd, use.names = TRUE, fill = TRUE)
  pd_all$features.plot <- factor(pd_all$features.plot,
                                  levels = rev(unique(pd_all$features.plot)))

  p <- ggplot(pd_all, aes(x = id, y = features.plot)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_color_gradient(low = "lightgrey", high = "red") +
    facet_grid(rows = vars(cluster), scales = "free_y",
               space = "free_y", switch = "y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.placement = "outside",
          strip.background = element_rect(fill = "grey90"),
          axis.title.y = element_blank()) +
    labs(title = "Top Marker Gene per Cluster by Treatment",
         x = "Treatment Group", color = "Avg Exp", size = "% Expressing")

  pdf(file.path(RESULTS_DIR, "04_dotplot_markers_by_treatment.pdf"),
      width = 11, height = 8)
  print(p)
  dev.off()
}

cat("Differential expression analysis complete.\n")

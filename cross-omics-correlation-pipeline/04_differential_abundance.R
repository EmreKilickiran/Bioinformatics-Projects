# =============================================================================
# 04_differential_abundance.R — ANCOM-BC2 Differential Abundance Analysis
# =============================================================================
#
# Identifies differentially abundant taxa across treatment conditions using
# ANCOM-BC2 (Analysis of Compositions of Microbiomes with Bias Correction),
# with Benjamini-Hochberg correction for multiple testing.
#
# Outputs:
#   - Significant species table (q < 0.05)
#   - Heatmap of top differentially abundant species
#   - Boxplot and violin plots for top species
#   - LEfSe-style LDA barplot (log2 fold change)
#
# Author : Yunus Emre Kılıçkıran
# =============================================================================

source("R/00_config.R")

cat("\n=== Module 4: Differential Abundance (ANCOM-BC2) ===\n")

# --- Load Data ---------------------------------------------------------------
otu_mat   <- load_otu_matrix("Abund_Species_wide", "Species")
meta_data <- load_metadata()

common_samples <- intersect(colnames(otu_mat), rownames(meta_data))
otu_mat   <- otu_mat[, common_samples]
meta_data <- meta_data[common_samples, ]

# --- Build Phyloseq Object ---------------------------------------------------
tax_tab <- data.frame(Species = rownames(otu_mat))
rownames(tax_tab) <- rownames(otu_mat)

phy <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax_tab)),
  sample_data(meta_data)
)

cat("Phyloseq:", nsamples(phy), "samples,", ntaxa(phy), "taxa\n")

# =============================================================================
# ANCOM-BC2
# =============================================================================

cat("Running ANCOM-BC2 (this may take a few minutes)...\n")

ancom_res <- ancombc2(
  data         = phy,
  tax_level    = "Species",
  fix_formula  = "Treatment",
  rand_formula = NULL,
  p_adj_method = "BH",
  prv_cut      = 0.10,
  lib_cut      = 0,
  neg_lb       = TRUE,
  global       = TRUE,
  group        = "Treatment"
)

# --- Reshape Results ---------------------------------------------------------
res_long <- ancom_res$res %>%
  pivot_longer(-taxon,
               names_to = c("metric", "contrast"),
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = metric, values_from = value, values_fn = list) %>%
  unnest(cols = c(lfc, se, W, p, q, diff))

sig_res <- res_long %>% filter(q < 0.05)

cat("Significant taxa (q < 0.05):", length(unique(sig_res$taxon)), "\n")

write.xlsx(sig_res,
           file.path(RESULTS_DIR, "04_ANCOMBC2_SignificantSpecies.xlsx"))

# =============================================================================
# Heatmap — Top Differentially Abundant Species
# =============================================================================

sig_taxa <- unique(sig_res$taxon)
top_taxa <- head(sig_taxa, 20)

if (length(top_taxa) > 0) {
  mat_sig <- otu_mat[top_taxa, , drop = FALSE]
  mat_sig <- mat_sig / rowSums(mat_sig)

  # Set factor levels for consistent ordering
  meta_data$Treatment <- factor(meta_data$Treatment, levels = TREATMENT_ORDER)
  meta_data$Bodysite  <- factor(meta_data$Bodysite,  levels = BODYSITE_ORDER)
  meta_data$SampleID  <- rownames(meta_data)
  meta_data <- meta_data %>%
    arrange(Bodysite, Treatment, as.numeric(str_extract(SampleID, "\\d+")))

  mat_sig <- mat_sig[, meta_data$SampleID]

  # Annotation
  annotation_df <- meta_data[, c("Treatment", "Bodysite")]
  ann_colors <- list(
    Treatment = TREATMENT_COLORS,
    Bodysite  = BODYSITE_COLORS
  )

  pdf(file.path(RESULTS_DIR, "04_Heatmap_TopDifferentialSpecies.pdf"),
      width = 20, height = 12)
  pheatmap(
    mat_sig,
    scale            = "row",
    annotation_col   = annotation_df,
    annotation_colors = ann_colors,
    cluster_cols     = FALSE,
    fontsize_row     = 10,
    fontsize_col     = 5,
    angle_col        = 45,
    legend_breaks    = c(-1, 0, 1),
    legend_labels    = c("-1", "0", "+1"),
    breaks           = seq(-1, 1, length = 101),
    color            = colorRampPalette(c("blue", "white", "red"))(100),
    main             = paste("Top", length(top_taxa), "Differentially Abundant Species")
  )
  dev.off()
  cat("Heatmap saved.\n")
}

# =============================================================================
# Boxplots & Violin Plots — Top 5 Species
# =============================================================================

otu_wide <- read_excel(DATA_PATH, sheet = "Abund_Species_wide")
otu_long <- otu_wide %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Sample = gsub("^Abundance\\.", "", Sample)) %>%
  left_join(meta_data, by = c("Sample" = "SampleID"))

sig_to_plot <- head(sig_taxa, 5)

if (length(sig_to_plot) > 0) {
  pdf(file.path(RESULTS_DIR, "04_Boxplots_ViolinPlots.pdf"))
  for (sp in sig_to_plot) {
    sp_data <- otu_long %>% filter(Species == sp)

    p1 <- ggplot(sp_data, aes(x = Treatment, y = Abundance, fill = Treatment)) +
      geom_boxplot() +
      geom_jitter(width = 0.2) +
      theme_bw() +
      ggtitle(paste("Boxplot —", sp)) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    p2 <- ggplot(sp_data, aes(x = Treatment, y = Abundance, fill = Treatment)) +
      geom_violin(trim = FALSE) +
      geom_jitter(width = 0.2) +
      theme_bw() +
      ggtitle(paste("Violin —", sp)) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    print(ggarrange(p1, p2, ncol = 2))
  }
  dev.off()
  cat("Boxplots saved.\n")
}

# =============================================================================
# LEfSe-style LDA Barplot (Log2 Fold Change)
# =============================================================================

lda_data <- sig_res %>%
  group_by(taxon) %>%
  summarize(mean_lfc = mean(lfc, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(abs(mean_lfc)))

pdf(file.path(RESULTS_DIR, "04_LEfSe_LDA_Barplot.pdf"), width = 10, height = 8)
print(
  ggplot(lda_data, aes(x = reorder(taxon, mean_lfc), y = mean_lfc,
                       fill = mean_lfc > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
                      labels = c("Depleted", "Enriched"),
                      name = "Direction") +
    labs(y = "Log2 Fold Change (LDA-style)", x = "Species") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7))
)
dev.off()

cat("Differential abundance analysis complete.\n")

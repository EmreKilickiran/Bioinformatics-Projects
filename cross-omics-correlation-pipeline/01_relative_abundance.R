# =============================================================================
# 01_relative_abundance.R — Taxonomic Composition Visualization
# =============================================================================
#
# Constructs a taxonomic abundance matrix from raw 16S sequencing data using
# phyloseq in R, applying relative abundance normalization and threshold-based
# filtering to isolate dominant phyla across all experimental groups.
#
# Output:
#   - Relative abundance barplot (Phylum level, Top 10 >= 1%)
#   - Double-strip faceted by Bodysite and Treatment
#
# Author : Yunus Emre Kılıçkıran
# =============================================================================

source("R/00_config.R")

cat("\n=== Module 1: Relative Abundance (Phylum Level) ===\n")

# --- Load Data ---------------------------------------------------------------
otu_mat   <- load_otu_matrix("Abund_Phylum_wide", "Phylum")
meta_data <- load_metadata()

# Match samples between OTU table and metadata
common_samples <- intersect(colnames(otu_mat), rownames(meta_data))
otu_mat   <- otu_mat[, common_samples]
meta_data <- meta_data[common_samples, ]

# --- Build Phyloseq Object ---------------------------------------------------
tax_tab <- data.frame(Phylum = rownames(otu_mat))
rownames(tax_tab) <- rownames(otu_mat)

physeq <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax_tab)),
  sample_data(meta_data)
)

cat("Phyloseq object:", nsamples(physeq), "samples,", ntaxa(physeq), "taxa\n")

# --- Relative Abundance Transformation --------------------------------------
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))
df <- psmelt(physeq_rel)

# --- Filter Top 10 Phyla (>= 1% total abundance) ----------------------------
ABUNDANCE_THRESHOLD <- 0.01
TOP_N_PHYLA <- 10

phylum_abund <- df %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total))

keep_taxa <- phylum_abund %>%
  filter(total >= ABUNDANCE_THRESHOLD) %>%
  slice_max(order_by = total, n = TOP_N_PHYLA) %>%
  pull(Phylum)

df$Phylum <- ifelse(df$Phylum %in% keep_taxa, df$Phylum, "Other")

# --- Set Factor Levels for Consistent Ordering --------------------------------
df$Treatment <- factor(df$Treatment, levels = TREATMENT_ORDER)
df$Bodysite  <- factor(df$Bodysite, levels = c("KOT", "liver", "swab"))

# Order samples: Bodysite → Treatment → numeric ID
sample_order <- df %>%
  distinct(Sample, Bodysite, Treatment) %>%
  arrange(Bodysite, Treatment) %>%
  pull(Sample)

df$Sample <- factor(df$Sample, levels = sample_order)

# --- Generate Barplot --------------------------------------------------------
p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ Bodysite + Treatment, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(
    title = "Relative Abundance (Phylum Level, Top 10 >= 1%, Others Grouped)",
    y = "Relative Abundance",
    x = "Sample"
  ) +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
    strip.background = element_rect(fill = "grey90"),
    strip.text.x = element_text(size = 8, face = "bold"),
    legend.key.size = unit(0.3, "cm"),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(ncol = 2))

# --- Save --------------------------------------------------------------------
out_file <- file.path(RESULTS_DIR, "01_RelativeAbundance_Phylum_Top10.pdf")
ggsave(out_file, plot = p, width = 20, height = 9)
cat("Saved:", out_file, "\n")

# 
# 02_alpha_diversity.R - Alpha Diversity Analysis
#
# Computes alpha diversity metrics (Observed, Chao1, Shannon, Simpson) and
# performs statistical comparisons using:
#   - Kruskal-Wallis global test
#   - Pairwise Wilcoxon rank-sum tests with Benjamini-Hochberg correction
#
# Outputs:
#   - Violin + boxplots with significance annotations (Treatment & Bodysite)
#   - Pairwise p-value heatmaps for all metric × grouping combinations
#   - CSV tables of pairwise test results


source("R/00_config.R")

cat("\n=== Module 2: Alpha Diversity ===\n")

# --- Load Alpha Diversity Table ----------------------------------------------
alpha_div <- read_excel(DATA_PATH, sheet = "Alpha_Diversity_wide")
colnames(alpha_div) <- gsub("^value\\.", "", colnames(alpha_div))

METRICS <- c("Observed", "Chao1", "Shannon", "Simpson")

cat("Loaded alpha diversity for", nrow(alpha_div), "samples\n")
cat("Metrics:", paste(METRICS, collapse = ", "), "\n")


# Pairwise Wilcoxon Tests

run_pairwise_wilcoxon <- function(df, group_var, metric) {
  test <- pairwise.wilcox.test(df[[metric]], df[[group_var]],
                               p.adjust.method = "BH")
  mat <- as.data.frame(as.table(test$p.value))
  colnames(mat) <- c("Group1", "Group2", "p_value")
  mat$Metric   <- metric
  mat$Grouping <- group_var
  return(mat)
}

# Run tests for both grouping variables
results_bodysite <- do.call(rbind, lapply(
  METRICS, function(m) run_pairwise_wilcoxon(alpha_div, "Bodysite", m)
))
results_treatment <- do.call(rbind, lapply(
  METRICS, function(m) run_pairwise_wilcoxon(alpha_div, "Treatment", m)
))

# Save test results
write.csv(results_bodysite,
          file.path(RESULTS_DIR, "02_PairwiseWilcoxon_Bodysite.csv"),
          row.names = FALSE)
write.csv(results_treatment,
          file.path(RESULTS_DIR, "02_PairwiseWilcoxon_Treatment.csv"),
          row.names = FALSE)

cat("Pairwise Wilcoxon tests complete.\n")


# P-value Heatmaps

plot_pvalue_heatmap <- function(df, metric, grouping) {
  df_m <- df %>% filter(Metric == metric, Grouping == grouping)
  mat  <- acast(df_m, Group1 ~ Group2, value.var = "p_value")

  ggplot(melt(mat, na.rm = TRUE), aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "red", high = "white", na.value = "grey80",
      trans = "log10", name = "p-value"
    ) +
    theme_minimal() +
    labs(
      title = paste(metric, "— Pairwise Wilcoxon —", grouping),
      x = "", y = ""
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Generate 8 heatmaps (4 metrics × 2 groupings)
heatmaps <- list(
  plot_pvalue_heatmap(results_bodysite,  "Observed", "Bodysite"),
  plot_pvalue_heatmap(results_bodysite,  "Chao1",    "Bodysite"),
  plot_pvalue_heatmap(results_bodysite,  "Shannon",  "Bodysite"),
  plot_pvalue_heatmap(results_bodysite,  "Simpson",  "Bodysite"),
  plot_pvalue_heatmap(results_treatment, "Observed", "Treatment"),
  plot_pvalue_heatmap(results_treatment, "Chao1",    "Treatment"),
  plot_pvalue_heatmap(results_treatment, "Shannon",  "Treatment"),
  plot_pvalue_heatmap(results_treatment, "Simpson",  "Treatment")
)

combined_heatmap <- (heatmaps[[1]] | heatmaps[[2]]) /
                    (heatmaps[[3]] | heatmaps[[4]]) /
                    (heatmaps[[5]] | heatmaps[[6]]) /
                    (heatmaps[[7]] | heatmaps[[8]])

ggsave(file.path(RESULTS_DIR, "02_PairwiseWilcoxon_Heatmaps.pdf"),
       plot = combined_heatmap, width = 16, height = 18)


# Violin + Boxplots with Significance Annotations

plot_alpha_diversity <- function(df, group_var, metric) {
  kw   <- kruskal.test(as.formula(paste(metric, "~", group_var)), data = df)
  kw_p <- signif(kw$p.value, 3)

  ggplot(df, aes(x = .data[[group_var]], y = .data[[metric]],
                 fill = .data[[group_var]])) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = combn(unique(df[[group_var]]), 2, simplify = FALSE),
      label = "p.signif",
      hide.ns = TRUE
    ) +
    theme_bw() +
    labs(
      title = paste(metric, "by", group_var, "| Kruskal-Wallis p =", kw_p),
      y = "Alpha Diversity", x = ""
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8)
    )
}

# Treatment panels
treatment_plots <- (
  plot_alpha_diversity(alpha_div, "Treatment", "Observed") |
  plot_alpha_diversity(alpha_div, "Treatment", "Chao1")
) / (
  plot_alpha_diversity(alpha_div, "Treatment", "Shannon") |
  plot_alpha_diversity(alpha_div, "Treatment", "Simpson")
)

ggsave(file.path(RESULTS_DIR, "02_AlphaDiversity_Treatment.pdf"),
       plot = treatment_plots, width = 16, height = 14)

# Bodysite panels
bodysite_plots <- (
  plot_alpha_diversity(alpha_div, "Bodysite", "Observed") |
  plot_alpha_diversity(alpha_div, "Bodysite", "Chao1")
) / (
  plot_alpha_diversity(alpha_div, "Bodysite", "Shannon") |
  plot_alpha_diversity(alpha_div, "Bodysite", "Simpson")
)

ggsave(file.path(RESULTS_DIR, "02_AlphaDiversity_Bodysite.pdf"),
       plot = bodysite_plots, width = 16, height = 14)

cat("Alpha diversity analysis complete.\n")

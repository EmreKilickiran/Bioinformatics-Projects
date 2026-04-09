# 
# 03_beta_diversity.R - Beta Diversity Analysis
#
# Computes Jaccard dissimilarity and Weighted UniFrac distances, applies
# Principal Coordinates Analysis (PCoA) for dimensionality reduction, and
# quantifies statistical differences between treatment groups through global
# and pairwise PERMANOVA tests (999 permutations).
#
# Outputs:
#   - Jaccard PCoA: Global plots with ellipses per bodysite
#   - Jaccard PCoA: Pairwise comparison plots with spider diagrams
#   - Weighted UniFrac PCoA: Treatment and Bodysite colorings
#   - Weighted UniFrac PCoA: Swab-specific pairwise comparisons
#   - PERMANOVA statistics (Excel)


source("R/00_config.R")

cat("\n=== Module 3: Beta Diversity ===\n")

# --- Load Data ---------------------------------------------------------------
otu_mat   <- load_otu_matrix("Abund_Species_wide", "Species")
meta_data <- load_metadata()

common_samples <- intersect(colnames(otu_mat), rownames(meta_data))
otu_mat   <- otu_mat[, common_samples]
meta_data <- meta_data[common_samples, ]

meta_data$SampleID <- rownames(meta_data)

cat("Samples:", ncol(otu_mat), "| Taxa:", nrow(otu_mat), "\n")


# PART A: JACCARD DISSIMILARITY

cat("\n--- Part A: Jaccard Distance ---\n")

dist_jac  <- vegdist(t(otu_mat), method = "jaccard", binary = FALSE)
pcoa_jac  <- ape::pcoa(dist_jac)
pcoa_jac_df <- data.frame(
  pcoa_jac$vectors[, 1:2],
  SampleID = rownames(pcoa_jac$vectors)
)
pcoa_jac_df <- merge(pcoa_jac_df, meta_data[, c("SampleID", "Bodysite", "Treatment")],
                     by = "SampleID")

var_explained_jac <- pcoa_jac$values$Relative_eig[1:2] * 100

# --- Pairwise comparison plots with spider diagrams -------------------------

make_comparison_plot <- function(df, bodysite, treatments, colors, title_prefix,
                                dist_mat, var_exp) {
  sub_df   <- df %>% filter(Bodysite == bodysite, Treatment %in% treatments)
  sub_dist <- as.dist(as.matrix(dist_mat)[sub_df$SampleID, sub_df$SampleID])

  # PERMANOVA
  adonis_res <- adonis2(sub_dist ~ Treatment, data = sub_df, permutations = 999)
  pval <- adonis_res$`Pr(>F)`[1]
  r2   <- adonis_res$R2[1]

  # Centroids for spider diagram
  centroids <- sub_df %>%
    group_by(Treatment) %>%
    summarise(Cx = mean(Axis.1), Cy = mean(Axis.2), .groups = "drop")
  sub_df <- sub_df %>% left_join(centroids, by = "Treatment")

  p_main <- ggplot(sub_df, aes(Axis.1, Axis.2, color = Treatment)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(aes(label = SampleID), size = 2.5,
                    max.overlaps = 30, show.legend = FALSE) +
    geom_point(data = centroids, aes(x = Cx, y = Cy),
               size = 5, shape = 18, inherit.aes = FALSE,
               color = "black") +
    geom_segment(aes(xend = Cx, yend = Cy), linewidth = 0.5, alpha = 0.6) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(
      title = paste0(title_prefix, " (", bodysite, ")"),
      x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PCoA2 (", round(var_exp[2], 1), "%)")
    ) +
    theme(plot.title = element_text(size = 12, face = "bold"))

  # Annotation box
  annot <- paste0("Jaccard  R² = ", round(r2, 3), "  |  p = ", signif(pval, 3))
  p_annot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = annot, size = 3, hjust = 0.5) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill = NA))

  combined <- p_main / p_annot + plot_layout(heights = c(4, 1))

  return(list(
    plot  = combined,
    stats = data.frame(Bodysite = bodysite, Comparison = title_prefix,
                       R2 = round(r2, 4), p_value = signif(pval, 3))
  ))
}

# Define pairwise comparisons
comparisons <- list(
  list(c("WSD", "WSD_Ligatur"), "WSD vs WSD_Ligatur",
       c("WSD" = "#1f77b4", "WSD_Ligatur" = "#e377c2")),
  list(c("WSD_PPI", "WSD_Ligatur_PPI"), "WSD_PPI vs WSD_Ligatur_PPI",
       c("WSD_PPI" = "#2ca02c", "WSD_Ligatur_PPI" = "#ff7f0e")),
  list(c("WSD_AB", "WSD_Ligatur_AB"), "WSD_AB vs WSD_Ligatur_AB",
       c("WSD_AB" = "#9467bd", "WSD_Ligatur_AB" = "#8c564b"))
)

bodysites  <- c("swab", "liver", "KOT")
all_stats  <- data.frame()

pdf(file.path(RESULTS_DIR, "03a_BetaDiversity_Jaccard_Pairwise.pdf"),
    width = 15, height = 5)
for (bod in bodysites) {
  plots <- list()
  for (comp in comparisons) {
    res <- make_comparison_plot(pcoa_jac_df, bod, comp[[1]], comp[[3]], comp[[2]],
                                dist_jac, var_explained_jac)
    plots <- append(plots, list(res$plot))
    all_stats <- rbind(all_stats, res$stats)
  }
  print(plots[[1]] | plots[[2]] | plots[[3]])
}
dev.off()

write.xlsx(all_stats,
           file.path(RESULTS_DIR, "03a_BetaDiversity_Jaccard_PairwiseStats.xlsx"))

# --- Global PCoA with ellipses per bodysite ----------------------------------

make_global_plot <- function(df, bodysite, dist_mat, var_exp) {
  sub_df   <- df %>% filter(Bodysite == bodysite)
  sub_dist <- as.dist(as.matrix(dist_mat)[sub_df$SampleID, sub_df$SampleID])

  global_res <- adonis2(sub_dist ~ Treatment, data = sub_df, permutations = 999)
  global_p  <- global_res$`Pr(>F)`[1]
  global_r2 <- global_res$R2[1]

  pw_res <- pairwise_permanova(sub_dist, sub_df$Treatment, nperm = 999)

  p <- ggplot(sub_df, aes(Axis.1, Axis.2, color = Treatment)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(aes(label = SampleID), size = 2,
                    max.overlaps = 200, show.legend = FALSE) +
    stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.7) +
    theme_bw() +
    labs(
      title    = paste0("Global PCoA — ", bodysite),
      subtitle = paste0("PERMANOVA: R² = ", round(global_r2, 3),
                        " | p = ", signif(global_p, 3)),
      x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PCoA2 (", round(var_exp[2], 1), "%)")
    ) +
    theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 11, face = "italic")
    )

  return(list(plot = p, pairwise = pw_res))
}

global_results <- list()
pdf(file.path(RESULTS_DIR, "03a_BetaDiversity_Jaccard_Global.pdf"),
    width = 6, height = 5)
for (bod in bodysites) {
  res <- make_global_plot(pcoa_jac_df, bod, dist_jac, var_explained_jac)
  print(res$plot)
  global_results[[bod]] <- res$pairwise
}
dev.off()

# Save pairwise PERMANOVA (one sheet per bodysite)
wb <- createWorkbook()
for (bod in names(global_results)) {
  addWorksheet(wb, bod)
  writeData(wb, bod, global_results[[bod]])
}
saveWorkbook(wb, file.path(RESULTS_DIR, "03a_BetaDiversity_Jaccard_GlobalPairwise.xlsx"),
             overwrite = TRUE)

cat("Jaccard beta diversity complete.\n")


# PART B: WEIGHTED UNIFRAC

cat("\n--- Part B: Weighted UniFrac ---\n")

beta_w   <- read_excel(DATA_PATH, sheet = "Beta_Diversity_wunifrac")
dist_mat_w <- acast(beta_w, Var1 ~ Var2, value.var = "value")
dist_mat_w[lower.tri(dist_mat_w)] <- t(dist_mat_w)[lower.tri(dist_mat_w)]
dist_w <- as.dist(dist_mat_w)

pcoa_w  <- ape::pcoa(dist_w)
pcoa_w_df <- data.frame(
  pcoa_w$vectors[, 1:2],
  SampleID = rownames(pcoa_w$vectors)
)
pcoa_w_df <- merge(pcoa_w_df, meta_data[, c("SampleID", "Bodysite", "Treatment")],
                   by = "SampleID")

var_explained_w <- pcoa_w$values$Relative_eig[1:2] * 100

# --- Treatment and Bodysite PCoA plots (2×2 grid) ---------------------------

make_unifrac_plots <- function(df, treatments, title_prefix, var_exp) {
  sub_df <- df %>% filter(Treatment %in% treatments)

  p_treat <- ggplot(sub_df, aes(Axis.1, Axis.2, color = Treatment, shape = Bodysite)) +
    geom_point(size = 2.5, alpha = 0.9) +
    stat_ellipse(aes(color = Treatment), type = "t", linewidth = 0.6) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    labs(
      title = paste0(title_prefix, " (by Treatment)"),
      x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PCoA2 (", round(var_exp[2], 1), "%)")
    ) +
    theme(plot.title = element_text(size = 12, face = "bold"))

  p_body <- ggplot(sub_df, aes(Axis.1, Axis.2, color = Bodysite, shape = Treatment)) +
    geom_point(size = 2.5, alpha = 0.9) +
    stat_ellipse(aes(color = Bodysite), type = "t", linewidth = 0.6) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    labs(
      title = paste0(title_prefix, " (by Bodysite)"),
      x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PCoA2 (", round(var_exp[2], 1), "%)")
    ) +
    theme(plot.title = element_text(size = 12, face = "bold"))

  return(list(treat = p_treat, body = p_body))
}

unifrac_comparisons <- list(
  list(c("WSD", "WSD_Ligatur"), "WSD vs WSD_Ligatur"),
  list(c("WSD_AB", "WSD_Ligatur_AB"), "WSD_AB vs WSD_Ligatur_AB"),
  list(c("WSD_PPI", "WSD_Ligatur_PPI"), "WSD_PPI vs WSD_Ligatur_PPI"),
  list(c("WSD", "WSD_Ligatur", "WSD_PPI", "WSD_Ligatur_PPI"), "All Non-AB Groups")
)

plots_treat <- list()
plots_body  <- list()
for (comp in unifrac_comparisons) {
  res <- make_unifrac_plots(pcoa_w_df, comp[[1]], comp[[2]], var_explained_w)
  plots_treat <- append(plots_treat, list(res$treat))
  plots_body  <- append(plots_body, list(res$body))
}

pdf(file.path(RESULTS_DIR, "03b_BetaDiversity_UniFrac_Treatment.pdf"),
    width = 12, height = 10)
print((plots_treat[[1]] | plots_treat[[2]]) / (plots_treat[[3]] | plots_treat[[4]]))
dev.off()

pdf(file.path(RESULTS_DIR, "03b_BetaDiversity_UniFrac_Bodysite.pdf"),
    width = 12, height = 10)
print((plots_body[[1]] | plots_body[[2]]) / (plots_body[[3]] | plots_body[[4]]))
dev.off()

# --- Weighted UniFrac PERMANOVA statistics -----------------------------------

permanova_results <- data.frame()
for (comp in unifrac_comparisons) {
  sub_df   <- pcoa_w_df %>% filter(Treatment %in% comp[[1]])
  sub_dist <- as.dist(as.matrix(dist_w)[sub_df$SampleID, sub_df$SampleID])
  adonis_res <- adonis2(sub_dist ~ Treatment, data = sub_df, permutations = 999)

  permanova_results <- rbind(permanova_results, data.frame(
    Comparison = comp[[2]],
    R2         = round(adonis_res$R2[1], 4),
    p_value    = adonis_res$`Pr(>F)`[1],
    Metric     = "Weighted UniFrac"
  ))
}

write.xlsx(permanova_results,
           file.path(RESULTS_DIR, "03b_BetaDiversity_UniFrac_PERMANOVA.xlsx"))

# --- Swab-specific comparisons (Weighted UniFrac) ----------------------------

make_swab_unifrac_plot <- function(df, treatments, colors, title_prefix,
                                   dist_mat, var_exp) {
  sub_df   <- df %>% filter(Bodysite == "swab", Treatment %in% treatments)
  sub_dist <- as.dist(as.matrix(dist_mat)[sub_df$SampleID, sub_df$SampleID])

  adonis_res <- adonis2(sub_dist ~ Treatment, data = sub_df, permutations = 999)
  pval <- adonis_res$`Pr(>F)`[1]
  r2   <- adonis_res$R2[1]

  p_main <- ggplot(sub_df, aes(Axis.1, Axis.2, color = Treatment)) +
    geom_point(size = 3, alpha = 0.9) +
    stat_ellipse(type = "t", linewidth = 0.7) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(
      title = paste0(title_prefix, " (Swabs Only)"),
      x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PCoA2 (", round(var_exp[2], 1), "%)")
    ) +
    theme(plot.title = element_text(size = 12, face = "bold"))

  annot <- paste0("Weighted UniFrac  R² = ", round(r2, 3),
                  "  |  p = ", signif(pval, 3))
  p_annot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = annot, size = 4, hjust = 0.5) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill = NA))

  p_main / p_annot + plot_layout(heights = c(4, 1))
}

p_swab1 <- make_swab_unifrac_plot(
  pcoa_w_df, c("WSD", "WSD_Ligatur"),
  c("WSD" = "#1f77b4", "WSD_Ligatur" = "#e377c2"),
  "WSD vs WSD_Ligatur", dist_w, var_explained_w
)
p_swab2 <- make_swab_unifrac_plot(
  pcoa_w_df, c("WSD_PPI", "WSD_Ligatur_PPI"),
  c("WSD_PPI" = "#2ca02c", "WSD_Ligatur_PPI" = "#ff7f0e"),
  "WSD_PPI vs WSD_Ligatur_PPI", dist_w, var_explained_w
)
p_swab3 <- make_swab_unifrac_plot(
  pcoa_w_df, c("WSD_AB", "WSD_Ligatur_AB"),
  c("WSD_AB" = "#9467bd", "WSD_Ligatur_AB" = "#8c564b"),
  "WSD_AB vs WSD_Ligatur_AB", dist_w, var_explained_w
)

pdf(file.path(RESULTS_DIR, "03b_BetaDiversity_UniFrac_Swabs.pdf"),
    width = 15, height = 7)
print(p_swab1 | p_swab2 | p_swab3)
dev.off()

cat("Weighted UniFrac beta diversity complete.\n")

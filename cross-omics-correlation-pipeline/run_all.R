#!/usr/bin/env Rscript
# 
# run_all.R — Execute the Full Analysis Pipeline
#
# Runs all analysis modules in sequence:
#   1. Relative abundance visualization (phylum level)
#   2. Alpha diversity (Shannon, Chao1, Simpson, Observed)
#   3. Beta diversity (Jaccard + Weighted UniFrac + PERMANOVA)
#   4. Differential abundance (ANCOM-BC2)
#   5. Cross-omics correlation (16S × RNA-Seq, 52M+ pairs)
#
# Usage:
#   Rscript run_all.R          # Run all modules
#   Rscript run_all.R 3        # Run only module 3 (beta diversity)
#   Rscript run_all.R 1 2 3    # Run modules 1, 2, and 3


args <- commandArgs(trailingOnly = TRUE)

modules <- list(
  "1" = "R/01_relative_abundance.R",
  "2" = "R/02_alpha_diversity.R",
  "3" = "R/03_beta_diversity.R",
  "4" = "R/04_differential_abundance.R",
  "5" = "R/05_cross_omics_correlation.R"
)

if (length(args) == 0) {
  to_run <- names(modules)
} else {
  to_run <- args
}

cat(" Oral-Gut-Liver Axis Microbiome Analysis Pipeline\n")
cat(" Modules to run:", paste(to_run, collapse = ", "), "\n")

for (mod in to_run) {
  if (mod %in% names(modules)) {
    cat("\n>>> Running Module", mod, ":", modules[[mod]], "\n")
    source(modules[[mod]])
  } else {
    cat("Warning: Unknown module", mod, "(valid: 1-5)\n")
  }
}

cat(" Pipeline complete. Results saved to: results/\n")

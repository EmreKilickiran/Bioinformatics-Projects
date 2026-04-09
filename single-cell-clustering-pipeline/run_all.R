#!/usr/bin/env Rscript

# run_all.R — Execute the Full scRNA-seq Analysis Pipeline
#
# Usage:
#   Rscript run_all.R          # Run all modules
#   Rscript run_all.R 2        # Run only module 2
#   Rscript run_all.R 4 5      # Run modules 4 and 5

args <- commandArgs(trailingOnly = TRUE)

modules <- list(
  "1" = "R/01_preprocessing.R",
  "2" = "R/02_dimensionality_reduction.R",
  "3" = "R/03_cell_type_annotation.R",
  "4" = "R/04_differential_expression.R",
  "5" = "R/05_functional_enrichment.R"
)

to_run <- if (length(args) == 0) names(modules) else args

cat(" scRNA-seq Liver Disease Progression Pipeline\n")
cat(" 310,000+ cells | 6 disease stages | 11 cell types\n")
cat(" Modules:", paste(to_run, collapse = ", "), "\n")

for (mod in to_run) {
  if (mod %in% names(modules)) {
    cat("\n>>> Module", mod, ":", modules[[mod]], "\n")
    source(modules[[mod]])
  } else {
    cat("Warning: Unknown module", mod, "\n")
  }
}

cat("\n=== Pipeline complete ===\n")

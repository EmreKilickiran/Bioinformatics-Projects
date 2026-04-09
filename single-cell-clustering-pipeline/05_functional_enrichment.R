# 05_functional_enrichment.R — GO & KEGG Enrichment Analysis
#
# Performs over-representation analysis (ORA) on differentially expressed genes:
#   - GO Biological Process (BP)
#   - KEGG Pathways
#
# Two analysis modes:
#   A) Per cell type × disease stage (cell-type-specific pathways)
#   B) Per disease stage (treatment-level, all cell types merged)
#
# Input:  04_DE_genes_vsNORM.tsv (from Module 4)
# Output: Dotplots (PDF) + enrichment tables (TSV) per comparison


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
for (pkg in c("clusterProfiler", "org.Hs.eg.db", "enrichplot")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

source("R/00_config.R")

cat("\n=== Module 5: Functional Enrichment (GO BP + KEGG) ===\n")

# --- Parameters --------------------------------------------------------------
PADJ_CUT <- 0.05
LFC_CUT  <- 0.25
MIN_GENES <- 10
SHOW_N    <- 20     # Top categories in dotplot

# --- Helper Functions --------------------------------------------------------

sanitize <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)

save_enrichment <- function(obj, pdf_path, title_txt, tsv_path = NULL,
                            show_n = SHOW_N) {
  ok <- FALSE
  if (!is.null(obj)) {
    df <- tryCatch(as.data.frame(obj), error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) {
      if (!is.null(tsv_path)) fwrite(df, tsv_path, sep = "\t")
      p <- dotplot(obj, showCategory = show_n) + ggtitle(title_txt)
      ggsave(pdf_path, plot = p, width = 7, height = 7)
      ok <- TRUE
    }
  }
  if (!ok && !is.null(tsv_path)) {
    fwrite(data.table(info = "no enriched terms"), tsv_path, sep = "\t")
  }
  invisible(ok)
}

run_GO <- function(genes_sym, title_txt, out_base) {
  if (length(genes_sym) < MIN_GENES) {
    fwrite(data.table(info = paste("too few genes:", length(genes_sym))),
           paste0(out_base, ".tsv"), sep = "\t")
    return(invisible(NULL))
  }
  ego <- tryCatch(
    enrichGO(gene = genes_sym, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
             ont = "BP", pAdjustMethod = "BH", readable = TRUE),
    error = function(e) NULL
  )
  save_enrichment(ego, paste0(out_base, ".pdf"), title_txt,
                  paste0(out_base, ".tsv"))
}

run_KEGG <- function(genes_sym, title_txt, out_base) {
  if (length(genes_sym) < MIN_GENES) {
    fwrite(data.table(info = paste("too few genes:", length(genes_sym))),
           paste0(out_base, ".tsv"), sep = "\t")
    return(invisible(NULL))
  }
  ent <- tryCatch(
    bitr(genes_sym, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  entrez <- if (!is.null(ent) && "ENTREZID" %in% colnames(ent))
    unique(na.omit(ent$ENTREZID)) else character(0)
  if (length(entrez) < 1) {
    fwrite(data.table(info = "no ENTREZ mapped"), paste0(out_base, ".tsv"),
           sep = "\t")
    return(invisible(NULL))
  }
  ek <- tryCatch(
    enrichKEGG(gene = entrez, organism = "hsa", pAdjustMethod = "BH"),
    error = function(e) NULL
  )
  save_enrichment(ek, paste0(out_base, ".pdf"), title_txt,
                  paste0(out_base, ".tsv"))
}

# --- Load DE Results ---------------------------------------------------------
dt <- fread(file.path(RESULTS_DIR, "04_DE_genes_vsNORM.tsv"))
dt$gene <- as.character(dt$gene)
dt <- dt[!is.na(gene) & nzchar(gene)]

cat("  Total DE entries:", nrow(dt), "\n")


# PART A: Per Cell Type × Disease Stage

cat("\n--- Part A: Cell-type-specific enrichment ---\n")

out_dir_ct <- file.path(RESULTS_DIR, "05_enrichment_celltype")
if (!dir.exists(out_dir_ct)) dir.create(out_dir_ct, recursive = TRUE)

combos <- unique(dt[, .(cell_type, case_group, ref_group)])

for (i in seq_len(nrow(combos))) {
  ct <- combos$cell_type[i]
  cg <- combos$case_group[i]
  rg <- combos$ref_group[i]

  sub <- dt[cell_type == ct & case_group == cg & ref_group == rg]
  if (nrow(sub) == 0) next

  sig  <- sub[p_val_adj <= PADJ_CUT & abs(avg_log2FC) >= LFC_CUT]
  up   <- unique(sig[avg_log2FC > 0, gene])
  down <- unique(sig[avg_log2FC < 0, gene])

  message(sprintf("  %s | %s vs %s | UP=%d DOWN=%d", ct, cg, rg,
                  length(up), length(down)))

  od <- file.path(out_dir_ct, sanitize(ct),
                  paste0(sanitize(cg), "_vs_", sanitize(rg)))
  if (!dir.exists(od)) dir.create(od, recursive = TRUE)

  run_GO(up,   sprintf("%s | %s vs %s | GO BP UP", ct, cg, rg),
         file.path(od, "GO_UP"))
  run_GO(down, sprintf("%s | %s vs %s | GO BP DOWN", ct, cg, rg),
         file.path(od, "GO_DOWN"))
  run_KEGG(up,   sprintf("%s | %s vs %s | KEGG UP", ct, cg, rg),
           file.path(od, "KEGG_UP"))
  run_KEGG(down, sprintf("%s | %s vs %s | KEGG DOWN", ct, cg, rg),
           file.path(od, "KEGG_DOWN"))
}


# PART B: Treatment-Level (All Cell Types Merged)

cat("\n--- Part B: Treatment-level enrichment ---\n")

out_dir_treat <- file.path(RESULTS_DIR, "05_enrichment_treatment")
if (!dir.exists(out_dir_treat)) dir.create(out_dir_treat, recursive = TRUE)

treat_combos <- unique(dt[, .(case_group, ref_group)])

for (i in seq_len(nrow(treat_combos))) {
  cg <- treat_combos$case_group[i]
  rg <- treat_combos$ref_group[i]

  sub <- dt[case_group == cg & ref_group == rg]
  if (nrow(sub) == 0) next

  sig  <- sub[p_val_adj <= PADJ_CUT & abs(avg_log2FC) >= LFC_CUT]
  up   <- unique(sig[avg_log2FC > 0, gene])
  down <- unique(sig[avg_log2FC < 0, gene])

  message(sprintf("  ALL | %s vs %s | UP=%d DOWN=%d", cg, rg,
                  length(up), length(down)))

  od <- file.path(out_dir_treat, paste0(sanitize(cg), "_vs_", sanitize(rg)))
  if (!dir.exists(od)) dir.create(od, recursive = TRUE)

  run_GO(up,   sprintf("%s vs %s | GO BP UP", cg, rg),
         file.path(od, "GO_UP"))
  run_GO(down, sprintf("%s vs %s | GO BP DOWN", cg, rg),
         file.path(od, "GO_DOWN"))
  run_KEGG(up,   sprintf("%s vs %s | KEGG UP", cg, rg),
           file.path(od, "KEGG_UP"))
  run_KEGG(down, sprintf("%s vs %s | KEGG DOWN", cg, rg),
           file.path(od, "KEGG_DOWN"))
}

cat("\nFunctional enrichment analysis complete.\n")

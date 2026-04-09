#!/usr/bin/env bash
# 
# End-to-End Computational Pipeline for Automated Sequence Classification
# and Phylogenetic Analysis
# 
#
# Description:
#   A fully automated, multi-stage computational pipeline to process raw
#   high-throughput 16S rRNA sequencing data into classified taxonomic profiles
#   and phylogenetic tree structures. Applied to gut microbiome characterization
#   under bile acid-associated liver conditions.
#
#
# Pipeline Stages:
#   1. Lane merging & data import
#   2. Primer removal (Cutadapt)
#   3. Quality control & denoising (DADA2)
#   4. Taxonomic classification (SILVA 138.1, Naive Bayes)
#   5. Contaminant filtering (mitochondria, chloroplast)
#   6. Phylogenetic tree reconstruction (MAFFT + FastTree)
#   7. Alpha & beta diversity analysis
#   8. Data export for downstream analysis (phyloseq / R)
=

set -euo pipefail

# CONFIGURATION

# --- Paths ---
PROJECT_DIR=".../16s.bile.new"
RAW_DATA_DIR="${PROJECT_DIR}/01.RawData"
RESULTS_DIR="${PROJECT_DIR}/02.Results"
CLASSIFIER_DIR="${PROJECT_DIR}/03.Classifiers"
EXPORT_DIR="${PROJECT_DIR}/04.Export"

# --- Micromamba / QIIME 2 Environment ---
MAMBA_EXE=".../micromamba/micromamba"
QIIME2_ENV=".../.conda/envs/qiime2-2024.5"

# --- Primer Sequences (V3-V4 region) ---
PRIMER_FWD="CCTACGGGNGGCWGCAG"    # 341F
PRIMER_REV="GACTACHVGGGTATCTAATCC" # 785R

# --- DADA2 Truncation Parameters ---
# Determined from per-base quality score analysis:
# Forward reads: quality drops below Q30 at ~282 bp
# Reverse reads: quality drops below Q30 at ~245 bp
# Overlap: 282 + 245 - 444 = 83 bp (sufficient for merging)
TRUNC_LEN_FWD=282
TRUNC_LEN_REV=245
TRIM_LEFT_FWD=0
TRIM_LEFT_REV=0

# --- Rarefaction Depth ---
# Set to the minimum sample frequency after quality filtering
RAREFACTION_DEPTH=7402

# --- Resource Allocation ---
N_THREADS=8

# --- Metadata ---
METADATA_FILE="${PROJECT_DIR}/metadata_qiime2_final_clean.tsv"


# HELPER FUNCTION

# Wrapper to run QIIME 2 commands via Micromamba
qiime2_run() {
    "${MAMBA_EXE}" run -p "${QIIME2_ENV}" "$@"
}

log_step() {
    echo ""
    echo "================================================================="
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] STEP: $1"
    echo "================================================================="
}


# STAGE 1: LANE MERGING & DATA IMPORT----------------------

# Raw data was sequenced across two lanes (L1, L2).
# Merge corresponding lane files per sample before import.

stage_1_import() {
    log_step "1 - Lane Merging & Data Import"

    cd "${RAW_DATA_DIR}"

    # Merge multi-lane FASTQ files (L1 + L2) per sample
    for sample_dir in */; do
        sample="${sample_dir%/}"
        echo "  Merging lanes for sample: ${sample}"
        cat "${sample}"/*L1_1.fq.gz "${sample}"/*L2_1.fq.gz > "${sample}_R1.fastq.gz"
        cat "${sample}"/*L1_2.fq.gz "${sample}"/*L2_2.fq.gz > "${sample}_R2.fastq.gz"
    done

    # Import paired-end sequences into QIIME 2 artifact format
    # Using a manifest file that maps sample IDs to file paths
    cd "${PROJECT_DIR}"
    qiime2_run qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path manifest_final.tsv \
        --output-path "${RESULTS_DIR}/demux.qza" \
        --input-format PairedEndFastqManifestPhred33V2

    # Generate quality summary visualization
    qiime2_run qiime demux summarize \
        --i-data "${RESULTS_DIR}/demux.qza" \
        --o-visualization "${RESULTS_DIR}/demux-summary.qzv"

    echo "  >> Inspect demux-summary.qzv at https://view.qiime2.org"
    echo "     to determine optimal truncation parameters."
}


# STAGE 2: PRIMER REMOVAL (Cutadapt)----------------------

# Remove V3-V4 amplicon primers (341F / 785R).
# Reads without detectable primers are discarded.

stage_2_primer_removal() {
    log_step "2 - Primer Removal (Cutadapt)"

    qiime2_run qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "${RESULTS_DIR}/demux.qza" \
        --p-front-f "${PRIMER_FWD}" \
        --p-front-r "${PRIMER_REV}" \
        --p-match-adapter-wildcards \
        --p-match-read-wildcards \
        --p-discard-untrimmed \
        --o-trimmed-sequences "${RESULTS_DIR}/demux-trimmed.qza" \
        --verbose

    # Visualize post-trimming quality
    qiime2_run qiime demux summarize \
        --i-data "${RESULTS_DIR}/demux-trimmed.qza" \
        --o-visualization "${RESULTS_DIR}/demux-trimmed-summary.qzv"
}


# STAGE 3: DENOISING WITH DADA2----------------------

# DADA2 performs:
#   - Quality-based read truncation
#   - Error model learning & denoising
#   - Paired-end read merging
#   - Chimera removal
#   - PhiX filtering
#
# Outputs:
#   - ASV feature table (sample x ASV counts)
#   - Representative sequences (one per ASV)
#   - Denoising statistics

stage_3_dada2_denoise() {
    log_step "3 - DADA2 Denoising"

    qiime2_run qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "${RESULTS_DIR}/demux-trimmed.qza" \
        --p-trim-left-f ${TRIM_LEFT_FWD} \
        --p-trim-left-r ${TRIM_LEFT_REV} \
        --p-trunc-len-f ${TRUNC_LEN_FWD} \
        --p-trunc-len-r ${TRUNC_LEN_REV} \
        --o-representative-sequences "${RESULTS_DIR}/rep-seqs.qza" \
        --o-table "${RESULTS_DIR}/table.qza" \
        --o-denoising-stats "${RESULTS_DIR}/denoising-stats.qza" \
        --p-n-threads ${N_THREADS}

    # Visualize DADA2 outputs
    qiime2_run qiime metadata tabulate \
        --m-input-file "${RESULTS_DIR}/denoising-stats.qza" \
        --o-visualization "${RESULTS_DIR}/denoising-stats.qzv"

    qiime2_run qiime feature-table tabulate-seqs \
        --i-data "${RESULTS_DIR}/rep-seqs.qza" \
        --o-visualization "${RESULTS_DIR}/rep-seqs.qzv"

    qiime2_run qiime feature-table summarize \
        --i-table "${RESULTS_DIR}/table.qza" \
        --o-visualization "${RESULTS_DIR}/table-summary.qzv" \
        --m-sample-metadata-file "${METADATA_FILE}"
}


# STAGE 4: TAXONOMIC CLASSIFICATION (SILVA 138.1)----------------------

# Taxonomy is assigned using a Naive Bayes classifier trained on the
# SILVA 138.1 SSURef NR99 database, region-specifically extracted for
# the 341F-785R amplicon region (~444 bp).
#
# Classifier training (run once):
#   See train_classifier.sh for the full training workflow.

stage_4_taxonomy() {
    log_step "4 - Taxonomic Classification"

    local CLASSIFIER="${CLASSIFIER_DIR}/silva-138.1-ssu-nr99-341f-785r-classifier.qza"

    qiime2_run qiime feature-classifier classify-sklearn \
        --i-classifier "${CLASSIFIER}" \
        --i-reads "${RESULTS_DIR}/rep-seqs.qza" \
        --o-classification "${RESULTS_DIR}/taxonomy_16S.qza"

    qiime2_run qiime metadata tabulate \
        --m-input-file "${RESULTS_DIR}/taxonomy_16S.qza" \
        --o-visualization "${RESULTS_DIR}/taxonomy_16S.qzv"
}


# STAGE 5: QUALITY FILTERING----------------------

# Remove:
#   (a) Samples with total read count < 1000 (insufficient depth)
#   (b) ASVs classified as mitochondria or chloroplast (non-bacterial)

stage_5_filtering() {
    log_step "5 - Quality Filtering"

    # Filter low-depth samples
    qiime2_run qiime feature-table filter-samples \
        --i-table "${RESULTS_DIR}/table.qza" \
        --p-min-frequency 1000 \
        --o-filtered-table "${RESULTS_DIR}/table-freq-filtered.qza"

    # Remove mitochondrial and chloroplast sequences
    qiime2_run qiime taxa filter-table \
        --i-table "${RESULTS_DIR}/table-freq-filtered.qza" \
        --i-taxonomy "${RESULTS_DIR}/taxonomy_16S.qza" \
        --p-exclude mitochondria,chloroplast \
        --o-filtered-table "${RESULTS_DIR}/table_16S_final.qza" \
        --verbose

    # Visualize final filtered table
    qiime2_run qiime feature-table summarize \
        --i-table "${RESULTS_DIR}/table_16S_final.qza" \
        --m-sample-metadata-file "${METADATA_FILE}" \
        --o-visualization "${RESULTS_DIR}/table_16S_final.qzv"
}


# STAGE 6: TAXONOMIC BARPLOT VISUALIZATION----------------------

stage_6_taxonomy_barplot() {
    log_step "6 - Taxonomic Barplot"

    qiime2_run qiime taxa barplot \
        --i-table "${RESULTS_DIR}/table_16S_final.qza" \
        --i-taxonomy "${RESULTS_DIR}/taxonomy_16S.qza" \
        --m-metadata-file "${METADATA_FILE}" \
        --o-visualization "${RESULTS_DIR}/taxa-bar-plots.qzv"
}

=
# STAGE 7: PHYLOGENETIC TREE RECONSTRUCTION----------------------

# Build a rooted phylogenetic tree required for phylogenetic diversity
# metrics (Faith's PD, UniFrac).
# Pipeline: MAFFT alignment → alignment masking → FastTree (ML) → midpoint rooting

stage_7_phylogeny() {
    log_step "7 - Phylogenetic Tree Reconstruction (MAFFT + FastTree)"

    qiime2_run qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences "${RESULTS_DIR}/rep-seqs.qza" \
        --o-alignment "${RESULTS_DIR}/aligned-rep-seqs.qza" \
        --o-masked-alignment "${RESULTS_DIR}/masked-aligned-rep-seqs.qza" \
        --o-tree "${RESULTS_DIR}/unrooted-tree.qza" \
        --o-rooted-tree "${RESULTS_DIR}/rooted-tree.qza" \
        --p-n-threads ${N_THREADS}
}


# STAGE 8: ALPHA & BETA DIVERSITY ANALYSIS----------------------

stage_8_diversity() {
    log_step "8 - Alpha & Beta Diversity Analysis"

    qiime2_run qiime diversity core-metrics-phylogenetic \
        --i-phylogeny "${RESULTS_DIR}/rooted-tree.qza" \
        --i-table "${RESULTS_DIR}/table_16S_final.qza" \
        --p-sampling-depth ${RAREFACTION_DEPTH} \
        --m-metadata-file "${METADATA_FILE}" \
        --output-dir "${RESULTS_DIR}/diversity-metrics-results"
}


# STAGE 9: EXPORT FOR DOWNSTREAM ANALYSIS (R / phyloseq)----------------------

# Export QIIME 2 artifacts to standard formats for analysis in R:
#   - Rooted tree → Newick (.nwk)
#   - Feature table → BIOM v1.0 (JSON) with taxonomy & metadata
#   - Taxonomy → TSV

stage_9_export() {
    log_step "9 - Data Export for R / phyloseq"

    mkdir -p "${EXPORT_DIR}"

    # Export rooted phylogenetic tree (Newick format)
    qiime2_run qiime tools export \
        --input-path "${RESULTS_DIR}/rooted-tree.qza" \
        --output-path "${EXPORT_DIR}/tree"

    # Export taxonomy table
    qiime2_run qiime tools export \
        --input-path "${RESULTS_DIR}/taxonomy_16S.qza" \
        --output-path "${EXPORT_DIR}/taxonomy"

    # Fix taxonomy header for BIOM compatibility:
    # Replace 'Feature ID' → '#OTU ID', 'Taxon' → 'taxonomy'
    sed -i '1s/Feature ID/#OTU ID/' "${EXPORT_DIR}/taxonomy/taxonomy.tsv"
    sed -i '1s/Taxon/taxonomy/'     "${EXPORT_DIR}/taxonomy/taxonomy.tsv"

    # Export feature table as BIOM
    qiime2_run qiime tools export \
        --input-path "${RESULTS_DIR}/table_16S_final.qza" \
        --output-path "${EXPORT_DIR}/feature_table"

    # Add taxonomy and metadata annotations to BIOM file
    qiime2_run biom add-metadata \
        -i "${EXPORT_DIR}/feature_table/feature-table.biom" \
        -o "${EXPORT_DIR}/table-with-taxonomy.biom" \
        --observation-metadata-fp "${EXPORT_DIR}/taxonomy/taxonomy.tsv" \
        --sc-separated taxonomy

    qiime2_run biom add-metadata \
        -i "${EXPORT_DIR}/table-with-taxonomy.biom" \
        -o "${EXPORT_DIR}/table-with-taxonomy-metadata.biom" \
        -m "${METADATA_FILE}"

    # Convert to JSON-format BIOM (compatible with phyloseq::import_biom)
    qiime2_run biom convert \
        -i "${EXPORT_DIR}/table-with-taxonomy-metadata.biom" \
        -o "${EXPORT_DIR}/16S_bile_acid.biom" \
        --table-type="OTU table" \
        --to-json

    echo "  >> Exported files:"
    echo "     Tree:     ${EXPORT_DIR}/tree/tree.nwk"
    echo "     Taxonomy: ${EXPORT_DIR}/taxonomy/taxonomy.tsv"
    echo "     BIOM:     ${EXPORT_DIR}/16S_bile_acid.biom"
}


# MAIN EXECUTION----------------------

# Run all stages sequentially. Each stage can also be executed independently
# by calling the individual function.
#
# Usage:
#   ./pipeline.sh              # Run full pipeline
#   ./pipeline.sh <stage_num>  # Run a specific stage (1-9)

main() {
    echo "================================================================="
    echo " 16S rRNA Microbiome Analysis Pipeline"
    echo " Project: Bile Acid-Associated Liver Conditions"
    echo " Platform: TU Dresden HPC Cluster"
    echo "================================================================="

    mkdir -p "${RESULTS_DIR}" "${CLASSIFIER_DIR}" "${EXPORT_DIR}"

    if [[ $# -eq 0 ]]; then
        # Run all stages
        stage_1_import
        stage_2_primer_removal
        stage_3_dada2_denoise
        stage_4_taxonomy
        stage_5_filtering
        stage_6_taxonomy_barplot
        stage_7_phylogeny
        stage_8_diversity
        stage_9_export
    else
        # Run specific stage
        case "$1" in
            1) stage_1_import ;;
            2) stage_2_primer_removal ;;
            3) stage_3_dada2_denoise ;;
            4) stage_4_taxonomy ;;
            5) stage_5_filtering ;;
            6) stage_6_taxonomy_barplot ;;
            7) stage_7_phylogeny ;;
            8) stage_8_diversity ;;
            9) stage_9_export ;;
            *) echo "Error: Invalid stage number. Use 1-9." && exit 1 ;;
        esac
    fi

    log_step "PIPELINE COMPLETE"
}

main "$@"

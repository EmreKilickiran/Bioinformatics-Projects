#!/usr/bin/env bash
# =============================================================================
# SILVA 138.1 Naive Bayes Classifier Training
# =============================================================================
#
# This script trains a region-specific (341F/785R) Naive Bayes taxonomic
# classifier using the SILVA 138.1 SSURef NR99 reference database.
#
# This only needs to be run ONCE. The resulting classifier artifact is
# reused across analyses targeting the same amplicon region.
#
# Steps:
#   1. Download SILVA 138.1 SSURef NR99 sequences and taxonomy
#   2. Cull low-quality sequences (ambiguous bases, short reads)
#   3. Dereplicate sequences at 99% identity
#   4. In-silico PCR: extract the 341F-785R amplicon region (250-500 bp)
#   5. Train Naive Bayes classifier on extracted region
#
# Resource Requirements (HPC):
#   - Memory: ~640 GB (for classifier training step)
#   - CPUs: 8
#   - Time: ~24-48 hours
#
# Author : Yunus Emre Kılıçkıran
# Project: Schneider Lab, CRTD, TU Dresden
# =============================================================================

set -euo pipefail

# --- Environment ---
MAMBA_EXE=".../micromamba/micromamba"
QIIME2_ENV=".../.conda/envs/qiime2-2024.5"
OUTPUT_DIR=".../16s.bile.new/03.Classifiers"

# --- Primer Sequences ---
PRIMER_FWD="CCTACGGGNGGCWGCAG"     # 341F
PRIMER_REV="GACTACHVGGGTATCTAATCC"  # 785R
MIN_AMPLICON_LEN=250
MAX_AMPLICON_LEN=500

qiime2_run() {
    "${MAMBA_EXE}" run -p "${QIIME2_ENV}" "$@"
}

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

echo "[1/5] Downloading SILVA 138.1 SSURef NR99..."
qiime2_run qiime rescript get-silva-data \
    --p-version '138.1' \
    --p-target 'SSURef_NR99' \
    --p-include-species-labels \
    --o-silva-sequences silva-138.1-ssu-nr99-seqs.qza \
    --o-silva-taxonomy silva-138.1-ssu-nr99-tax.qza

echo "[2/5] Culling low-quality sequences..."
qiime2_run qiime rescript cull-seqs \
    --i-sequences silva-138.1-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138.1-ssu-nr99-seqs-clean.qza

echo "[3/5] Dereplicating sequences..."
qiime2_run qiime rescript dereplicate \
    --i-sequences silva-138.1-ssu-nr99-seqs-clean.qza \
    --i-taxa silva-138.1-ssu-nr99-tax.qza \
    --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep.qza \
    --o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep.qza

echo "[4/5] Extracting 341F-785R amplicon region..."
qiime2_run qiime feature-classifier extract-reads \
    --i-sequences silva-138.1-ssu-nr99-seqs-derep.qza \
    --p-f-primer "${PRIMER_FWD}" \
    --p-r-primer "${PRIMER_REV}" \
    --p-min-length ${MIN_AMPLICON_LEN} \
    --p-max-length ${MAX_AMPLICON_LEN} \
    --o-reads silva-138.1-ssu-nr99-341f-785r-seqs.qza

echo "[5/5] Training Naive Bayes classifier..."
echo "       (This step requires ~640 GB RAM and may take several hours)"
qiime2_run qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138.1-ssu-nr99-341f-785r-seqs.qza \
    --i-reference-taxonomy silva-138.1-ssu-nr99-tax-derep.qza \
    --o-classifier silva-138.1-ssu-nr99-341f-785r-classifier.qza

echo ""
echo "Classifier trained successfully:"
echo "  ${OUTPUT_DIR}/silva-138.1-ssu-nr99-341f-785r-classifier.qza"

#!/usr/bin/env bash
#SBATCH --job-name=16s-microbiome-pipeline
#SBATCH --output=pipeline_%j.log
#SBATCH --error=pipeline_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --time=1-00:00:00
#SBATCH --partition=haswell

# SLURM Submission Script for the 16S rRNA Microbiome Analysis Pipeline

# Submit with:  sbatch slurm_submit.sh
# Monitor with: squeue -u $USER
# Cancel with:  scancel <job_id>
#
# Note: Classifier training (train_classifier.sh) requires more memory
# (~640 GB) and should be submitted separately or run interactively.

bash pipeline.sh

#!/bin/bash
#SBATCH --job-name=konuseg_chm13
#SBATCH --partition=hi_end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/konuseg_chm13_%j.out
#SBATCH --error=logs/konuseg_chm13_%j.err

# =============================================================================
# CHM13 Convenience Wrapper
#
# Sets CHM13-specific defaults (input paths, sex, RepeatMasker) and delegates
# to the generic run_konuseg_pipeline.sh.
#
# Usage:
#   sbatch scripts/run_chm13_cnacc_job.sh my_run_v1
#
#   SEGMENTER=fused_lasso \
#   SKIP_PREPROCESS=true \
#   RM_ANNOTATION_BED=output/chm13_my_run_fl_v2/repeat_annotated_w500.bed \
#   sbatch scripts/run_chm13_cnacc_job.sh my_run_fl_v4
#
# All environment variables accepted by run_konuseg_pipeline.sh work here.
# This wrapper only provides CHM13-specific defaults for INPUT_KMERS,
# REPEATMASKER, and SEX.
# =============================================================================

# --- CHM13 defaults (override via env vars) ---
export INPUT_KMERS="${INPUT_KMERS:-/home/klea/recomb/gra_bf/output/final_algo/chm13/chm13_ont_cpn_merged_algo_k72_minseg500.bed_noNoiseFiltering_bloom_filter_kmers.bed}"
export REPEATMASKER="${REPEATMASKER:-/home/klea/recomb/gra_bf/data/human/chm13/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out}"
export SEX="${SEX:-XX}"
export SAMPLE="${SAMPLE:-chm13}"

# Delegate to generic pipeline, prefixing run name with "chm13_"
# NOTE: Cannot use dirname "$0" — SLURM copies the script to a temp dir.
exec scripts/run_konuseg_pipeline.sh "chm13_${1:-run_$(date +%Y%m%d_%H%M)}"

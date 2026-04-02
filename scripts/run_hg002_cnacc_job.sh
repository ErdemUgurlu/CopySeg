#!/bin/bash
#SBATCH --job-name=konuseg_hg002
#SBATCH --partition=hi_end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=logs/konuseg_hg002_%j.out
#SBATCH --error=logs/konuseg_hg002_%j.err

# =============================================================================
# HG002 Convenience Wrapper
#
# Sets HG002-specific defaults (input paths, sex, RepeatMasker) and delegates
# to the generic run_konuseg_pipeline.sh.
#
# Usage:
#   sbatch scripts/run_hg002_cnacc_job.sh my_run_v1
#
#   SEGMENTER=fused_lasso \
#   SKIP_PREPROCESS=true \
#   sbatch scripts/run_hg002_cnacc_job.sh my_run_fl_v2
#
# All environment variables accepted by run_konuseg_pipeline.sh work here.
# This wrapper only provides HG002-specific defaults for INPUT_KMERS,
# REPEATMASKER, SEX, and SAMPLE.
# =============================================================================

# --- HG002 defaults (override via env vars) ---
export INPUT_KMERS="${INPUT_KMERS:-/home/klea/recomb/gra_bf/output/final_algo/hg002_genome_chromosomes/human_hg002_trimmed_chromosomes_whole_genome_ont_cpn_merged_algo_k32_seed0_minseg1000.bed_noNoiseFiltering_bloom_filter_kmers.bed}"
export REPEATMASKER="${REPEATMASKER:-/home/klea/recomb/gra_bf/data/human/chm13/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out}"
export SEX="${SEX:-XY}"
export SAMPLE="${SAMPLE:-hg002}"

# Delegate to generic pipeline, prefixing run name with "hg002_"
# NOTE: Cannot use dirname "$0" — SLURM copies the script to a temp dir.
exec scripts/run_konuseg_pipeline.sh "hg002_${1:-run_$(date +%Y%m%d_%H%M)}"

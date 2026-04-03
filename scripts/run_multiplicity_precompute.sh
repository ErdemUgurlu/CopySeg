#!/bin/bash
#SBATCH --job-name=mult_precompute
#SBATCH --partition=hi_end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/mult_precompute_%j.out
#SBATCH --error=logs/mult_precompute_%j.err

set -euo pipefail

# =============================================================================
# KonuSeg — K-mer Multiplicity Weight Pre-computation
#
# One-time computation per (reference, k_target, k_ref) triple.
# Creates Jellyfish DBs if needed, then computes per-position weights.
#
# Usage:
#   # Default: k=32 target, k=72 reference, CHM13
#   sbatch scripts/run_multiplicity_precompute.sh
#
#   # Custom k sizes:
#   K_TARGET=50 K_REF=72 sbatch scripts/run_multiplicity_precompute.sh
#
#   # Specific reference:
#   REF_FASTA=/path/to/ref.fa sbatch scripts/run_multiplicity_precompute.sh
# =============================================================================

K_TARGET="${K_TARGET:-32}"
K_REF="${K_REF:-72}"

# Reference FASTA — auto-detect
REF_FASTA="${REF_FASTA:-}"
if [ -z "${REF_FASTA}" ]; then
    for candidate in \
        "/home/erdem.ugurlu/KonuSeg_Phase2/data/chm13v2.0.fa" \
        "/home/klea/recomb/gra_bf/data/human/chm13/chm13v2.0.fa"; do
        if [ -f "${candidate}" ]; then
            REF_FASTA="${candidate}"
            break
        fi
    done
fi

if [ -z "${REF_FASTA}" ] || [ ! -f "${REF_FASTA}" ]; then
    echo "ERROR: Reference FASTA not found. Set REF_FASTA=/path/to/ref.fa"
    exit 1
fi

# Output directory
OUTDIR="${OUTDIR:-output/mult_precompute_k${K_TARGET}}"
WEIGHT_DIR="${OUTDIR}/weights_k${K_TARGET}_ref_k${K_REF}"

# Jellyfish DB paths
JF_TARGET="${JF_TARGET:-${OUTDIR}/chm13v2_k${K_TARGET}.jf}"
JF_REF="${JF_REF:-}"

# Auto-detect existing k_ref DB
if [ -z "${JF_REF}" ]; then
    for candidate in \
        "output/k${K_REF}_gt_validation/chm13v2_k${K_REF}.jf" \
        "${OUTDIR}/chm13v2_k${K_REF}.jf"; do
        if [ -f "${candidate}" ]; then
            JF_REF="${candidate}"
            break
        fi
    done
fi

# Jellyfish hash size
JF_HASH_SIZE="${JF_HASH_SIZE:-10G}"

echo "=================================================="
echo "KonuSeg — Multiplicity Weight Pre-computation"
echo "=================================================="
echo "  Reference: ${REF_FASTA}"
echo "  k_target:  ${K_TARGET}"
echo "  k_ref:     ${K_REF}"
echo "  Output:    ${WEIGHT_DIR}"
echo "  JF target: ${JF_TARGET}"
echo "  JF ref:    ${JF_REF:-<will be created>}"
echo "  Job ID:    ${SLURM_JOB_ID:-local}"
echo "  Started:   $(date)"
echo "=================================================="
echo ""

# --- Conda environment ---
if command -v conda &>/dev/null; then
    eval "$(conda shell.bash hook)"
elif [ -f "/home/erdem.ugurlu/miniconda3/etc/profile.d/conda.sh" ]; then
    source "/home/erdem.ugurlu/miniconda3/etc/profile.d/conda.sh"
elif [ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]; then
    source "${HOME}/anaconda3/etc/profile.d/conda.sh"
fi
conda activate konuseg

mkdir -p "${OUTDIR}" "${WEIGHT_DIR}" logs

# =============================================================================
# Step 1: Build Jellyfish DB for k_target (if not existing)
# =============================================================================
if [ -f "${JF_TARGET}" ]; then
    echo "[JF] k_target DB exists: ${JF_TARGET}"
else
    echo "[JF] Building k=${K_TARGET} Jellyfish DB..."
    echo "     jellyfish count -m ${K_TARGET} -s ${JF_HASH_SIZE} -C -t 16"
    time jellyfish count \
        -m "${K_TARGET}" \
        -s "${JF_HASH_SIZE}" \
        -C \
        -t 16 \
        "${REF_FASTA}" \
        -o "${JF_TARGET}"
    echo "[JF] Done: ${JF_TARGET} ($(du -h "${JF_TARGET}" | cut -f1))"
fi
echo ""

# =============================================================================
# Step 2: Build Jellyfish DB for k_ref (if not existing)
# =============================================================================
if [ -z "${JF_REF}" ]; then
    JF_REF="${OUTDIR}/chm13v2_k${K_REF}.jf"
fi

if [ -f "${JF_REF}" ]; then
    echo "[JF] k_ref DB exists: ${JF_REF}"
else
    echo "[JF] Building k=${K_REF} Jellyfish DB..."
    echo "     jellyfish count -m ${K_REF} -s ${JF_HASH_SIZE} -C -t 16"
    time jellyfish count \
        -m "${K_REF}" \
        -s "${JF_HASH_SIZE}" \
        -C \
        -t 16 \
        "${REF_FASTA}" \
        -o "${JF_REF}"
    echo "[JF] Done: ${JF_REF} ($(du -h "${JF_REF}" | cut -f1))"
fi
echo ""

# =============================================================================
# Step 3: Compute per-position weights
# =============================================================================
echo "[WEIGHTS] Computing per-position multiplicity weights..."
echo ""

time python3 scripts/compute_multiplicity_weights.py \
    --ref-fasta "${REF_FASTA}" \
    --k-target "${K_TARGET}" \
    --k-ref "${K_REF}" \
    --jf-target "${JF_TARGET}" \
    --jf-ref "${JF_REF}" \
    --output-dir "${WEIGHT_DIR}" \
    --batch-size 5000000

echo ""
echo "=================================================="
echo "COMPLETE"
echo "  Weight dir: ${WEIGHT_DIR}"
echo "  Usage:      WEIGHT_DIR=${WEIGHT_DIR} sbatch scripts/run_chm13_cnacc_job.sh my_run"
echo "  Finished:   $(date)"
echo "=================================================="

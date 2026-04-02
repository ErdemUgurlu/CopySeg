#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# KonuSeg Parameter Sweep — tek değişken, N job
#
# Usage:
#   # CHM13 penalty factor sweep (default):
#   bash scripts/run_parameter_sweep.sh penalty_factor 2,4,6,8 sweep_pf
#
#   # HG002 penalty factor sweep:
#   SAMPLE_WRAPPER=hg002 \
#   bash scripts/run_parameter_sweep.sh penalty_factor 4,6,8,10,12 sweep_pf_k32
#
#   # Generic pipeline (no wrapper):
#   SAMPLE_WRAPPER=generic \
#   INPUT_KMERS=/path/to/kmers.bed SEX=XY SAMPLE=hg002 \
#   bash scripts/run_parameter_sweep.sh penalty_factor 4,8,12 sweep_pf
#
#   # Min segment length sweep:
#   bash scripts/run_parameter_sweep.sh min_segment_length 3000,5000,10000,15000 sweep_msl
#
# SAMPLE_WRAPPER: chm13 (default), hg002, or generic
# All other parameters use defaults below. Modify DEFAULTS section as needed.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 <param_name> <comma_separated_values> <sweep_prefix>"
    echo ""
    echo "Supported param_name:"
    echo "  penalty_factor      — BIC multiplier (default: 4)"
    echo "  min_segment_length  — minimum segment bp (default: 5000)"
    echo "  merge_cn_tolerance  — cross-state merge CN threshold (default: 0.5)"
    echo "  cv_split_threshold  — CV split threshold (default: 0)"
    echo "  min_kmers           — low-coverage filter (default: 0)"
    echo "  quality_threshold   — hard quality filter (default: 0.7)"
    echo "  lowdup_threshold    — Neutral/LowDup boundary CN (default: 1.5)"
    echo ""
    echo "Environment:"
    echo "  SAMPLE_WRAPPER  — chm13 (default), hg002, or generic"
    echo "  RM_ANNOTATION_BED — path to repeat annotation BED"
    echo ""
    echo "Example:"
    echo "  $0 penalty_factor 2,4,6,8 sweep_pf"
    echo "  SAMPLE_WRAPPER=hg002 $0 penalty_factor 4,8,12 sweep_pf_k32"
    exit 1
fi

PARAM_NAME="$1"
VALUES_CSV="$2"
SWEEP_PREFIX="$3"

# ── SAMPLE WRAPPER ──────────────────────────────────────────────────────────
SAMPLE_WRAPPER="${SAMPLE_WRAPPER:-chm13}"

case "${SAMPLE_WRAPPER}" in
    chm13)
        SUBMIT_SCRIPT="scripts/run_chm13_cnacc_job.sh"
        RUN_PREFIX="chm13_"
        ;;
    hg002)
        SUBMIT_SCRIPT="scripts/run_hg002_cnacc_job.sh"
        RUN_PREFIX="hg002_"
        ;;
    generic)
        SUBMIT_SCRIPT="scripts/run_konuseg_pipeline.sh"
        RUN_PREFIX=""
        ;;
    *)
        echo "ERROR: Unknown SAMPLE_WRAPPER '${SAMPLE_WRAPPER}'. Use: chm13, hg002, generic"
        exit 1
        ;;
esac

# ── DEFAULTS (modify these to match your current best config) ────────────────
DEFAULT_PENALTY_FACTOR="${DEFAULT_PENALTY_FACTOR:-4}"
DEFAULT_MIN_SEGMENT_LENGTH="${DEFAULT_MIN_SEGMENT_LENGTH:-5000}"
DEFAULT_MERGE_CN_TOLERANCE="${DEFAULT_MERGE_CN_TOLERANCE:-0.5}"
DEFAULT_CV_SPLIT_THRESHOLD="${DEFAULT_CV_SPLIT_THRESHOLD:-0}"
DEFAULT_MIN_KMERS="${DEFAULT_MIN_KMERS:-0}"
DEFAULT_QUALITY_THRESHOLD="${DEFAULT_QUALITY_THRESHOLD:-0.7}"
DEFAULT_LOWDUP_THRESHOLD="${DEFAULT_LOWDUP_THRESHOLD:-1.5}"

# ── FIXED SETTINGS ───────────────────────────────────────────────────────────
SEGMENTER="${SEGMENTER:-fused_lasso}"
SKIP_PREPROCESS="${SKIP_PREPROCESS:-true}"
SKIP_REPEAT="${SKIP_REPEAT:-true}"
RM_ANNOTATION_BED="${RM_ANNOTATION_BED:-}"
# PREPROCESS_SRC: path to an existing preprocessed BED.
# Pipeline will symlink this into each sweep run's output dir.
PREPROCESS_SRC="${PREPROCESS_SRC:-}"

# ── SWEEP ────────────────────────────────────────────────────────────────────
IFS=',' read -ra VALUES <<< "$VALUES_CSV"

echo "╔══════════════════════════════════════════════════════╗"
echo "║  KonuSeg Parameter Sweep                            ║"
echo "╠══════════════════════════════════════════════════════╣"
echo "║  Sweep variable:  ${PARAM_NAME}"
echo "║  Values:          ${VALUES_CSV}"
echo "║  N jobs:          ${#VALUES[@]}"
echo "║  Prefix:          ${SWEEP_PREFIX}"
echo "║  Sample wrapper:  ${SAMPLE_WRAPPER}"
echo "║  Submit script:   ${SUBMIT_SCRIPT}"
echo "║  RM annotation:   ${RM_ANNOTATION_BED:-[not set]}"
echo "║  Preprocess src:  ${PREPROCESS_SRC:-[not set — each job preprocesses or needs own BED]}"
echo "╚══════════════════════════════════════════════════════╝"
echo ""

JOBS_SUBMITTED=()

for VAL in "${VALUES[@]}"; do
    # Start with defaults
    PF="${DEFAULT_PENALTY_FACTOR}"
    MSL="${DEFAULT_MIN_SEGMENT_LENGTH}"
    MCT="${DEFAULT_MERGE_CN_TOLERANCE}"
    CVS="${DEFAULT_CV_SPLIT_THRESHOLD}"
    MK="${DEFAULT_MIN_KMERS}"
    QT="${DEFAULT_QUALITY_THRESHOLD}"
    LDT="${DEFAULT_LOWDUP_THRESHOLD}"

    # Override the sweep variable
    case "${PARAM_NAME}" in
        penalty_factor)     PF="${VAL}" ;;
        min_segment_length) MSL="${VAL}" ;;
        merge_cn_tolerance) MCT="${VAL}" ;;
        cv_split_threshold) CVS="${VAL}" ;;
        min_kmers)          MK="${VAL}" ;;
        quality_threshold)  QT="${VAL}" ;;
        lowdup_threshold)   LDT="${VAL}" ;;
        *)
            echo "ERROR: Unknown parameter '${PARAM_NAME}'"
            exit 1
            ;;
    esac

    # Sanitize value for run name (replace . with p)
    VAL_SAFE=$(echo "${VAL}" | tr '.' 'p')
    RUN_NAME="${SWEEP_PREFIX}_${PARAM_NAME}_${VAL_SAFE}"

    echo "─── Submitting: ${RUN_PREFIX}${RUN_NAME} (${PARAM_NAME}=${VAL}) ───"

    SEGMENTER="${SEGMENTER}" \
    SKIP_PREPROCESS="${SKIP_PREPROCESS}" \
    SKIP_REPEAT="${SKIP_REPEAT}" \
    RM_ANNOTATION_BED="${RM_ANNOTATION_BED}" \
    PREPROCESS_SRC="${PREPROCESS_SRC}" \
    PENALTY_FACTOR="${PF}" \
    MIN_SEGMENT_LENGTH="${MSL}" \
    MERGE_CN_TOLERANCE="${MCT}" \
    CV_SPLIT_THRESHOLD="${CVS}" \
    MIN_KMERS="${MK}" \
    LOWDUP_THRESHOLD="${LDT}" \
    sbatch "${SUBMIT_SCRIPT}" "${RUN_NAME}"

    JOBS_SUBMITTED+=("${RUN_NAME}")
    echo ""
done

echo "╔══════════════════════════════════════════════════════╗"
echo "║  Sweep submitted: ${#JOBS_SUBMITTED[@]} jobs                        ║"
echo "╠══════════════════════════════════════════════════════╣"
for JOB in "${JOBS_SUBMITTED[@]}"; do
    echo "║  → ${RUN_PREFIX}${JOB}"
done
echo "╚══════════════════════════════════════════════════════╝"
echo ""
echo "Monitor: squeue -u \$USER"
echo "Results: ls output/${RUN_PREFIX}${SWEEP_PREFIX}_*/validation/"
echo ""
echo "After completion, compare GT results:"
echo "  for d in output/${RUN_PREFIX}${SWEEP_PREFIX}_*/validation/; do"
echo '    echo "=== $(basename $(dirname $d)) ===" && head -20 "$d/gt_evaluation_report.md"'
echo "  done"

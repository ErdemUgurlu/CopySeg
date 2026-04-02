#!/bin/bash
#SBATCH --job-name=konuseg
#SBATCH --partition=hi_end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/konuseg_%j.out
#SBATCH --error=logs/konuseg_%j.err

set -euo pipefail
export PYTHONUNBUFFERED=1   # Python print() output flushed immediately

# =============================================================================
# KonuSeg — Generic CN Calling Pipeline
#
# Works with any sample, any technology (ONT/PacBio), any k-mer size.
#
# Required environment variables:
#   INPUT_KMERS       — path to k-mer BED file (4-col: chrom, start, end, count)
#
# Optional environment variables:
#   REPEATMASKER      — RepeatMasker .out file (for RM-guided bio-filter)
#   SEX               — XX (default) or XY (excludes chrX/chrY from median)
#   SEGMENTER         — hmm (default) or fused_lasso
#   KMER_SIZE         — auto-detected from filename, or 72 default
#   GC_BED            — GC content BED for post-segmentation calibration
#   RM_ANNOTATION_BED — pre-computed repeat annotation BED
#   SKIP_PREPROCESS   — true to reuse existing preprocessed windows
#   SKIP_SEGMENTER    — true to reuse existing segments
#   SKIP_REPEAT       — true to reuse existing repeat annotation
#
# Steps:
#   1.   preprocess_kmer_windows.py   — sliding windows → binned windows
#   1.5. compute_repeat_annotation.py — RM .out → masked_fraction+repeat_class
#   2.   Segmentation (HMM or PELT)
#   3.   validate_cn_accuracy.py      — CN accuracy metrics
#   3b.  evaluate_ground_truth.py     — GT locus evaluation (CHM13-specific)
#
# Usage:
#   INPUT_KMERS=/path/to/kmers.bed sbatch scripts/run_konuseg_pipeline.sh my_run
#
#   INPUT_KMERS=/path/to/kmers.bed \
#   REPEATMASKER=/path/to/RepeatMasker.out \
#   SEGMENTER=fused_lasso SEX=XX \
#   sbatch scripts/run_konuseg_pipeline.sh my_run_fl
# =============================================================================

RUN_NAME="${1:-run_$(date +%Y%m%d_%H%M)}"

# --- Input paths ---
if [ -z "${INPUT_KMERS:-}" ]; then
    echo "ERROR: INPUT_KMERS must be set (path to k-mer BED file)"
    echo "  Usage: INPUT_KMERS=/path/to/kmers.bed sbatch scripts/run_konuseg_pipeline.sh my_run"
    exit 1
fi
ONT_KMERS="${INPUT_KMERS}"
REPEATMASKER="${REPEATMASKER:-}"

# --- Sex karyotype ---
SEX="${SEX:-XX}"

# --- Sample name (for GT evaluation) ---
# Selects built-in GT table: chm13, hg002
# Set to empty string to skip GT evaluation entirely.
SAMPLE="${SAMPLE:-}"

# --- K-mer size ---
# Auto-detect from INPUT_KMERS filename (matches _k32_ or _k72_ pattern).
# Override via: KMER_SIZE=32 sbatch ...
if [ -z "${KMER_SIZE:-}" ]; then
    if [[ "${ONT_KMERS}" =~ _k([0-9]+)_ ]]; then
        KMER_SIZE="${BASH_REMATCH[1]}"
    else
        echo "WARNING: Could not auto-detect k-mer size from filename."
        echo "         Defaulting to KMER_SIZE=72. Set KMER_SIZE explicitly if different."
        KMER_SIZE=72
    fi
fi

# --- Data source label ---
if [[ "${ONT_KMERS}" == *"pacbio"* ]]; then
    DATA_SOURCE="PacBio"
else
    DATA_SOURCE="ONT"
fi

# GC content BED — disabled by default (GC bias typically negligible).
GC_BED="${GC_BED:-}"

# RM annotation BED for per-bin threshold (Step 1).
# Auto-detect from current output dir if not set explicitly.
# For reuse across runs, set RM_ANNOTATION_BED explicitly.

# --- Segmenter selection ---
SEGMENTER="${SEGMENTER:-hmm}"

# --- Skip flags ---
SKIP_PREPROCESS="${SKIP_PREPROCESS:-false}"
SKIP_SEGMENTER="${SKIP_SEGMENTER:-false}"
SKIP_HMM="${SKIP_HMM:-${SKIP_SEGMENTER}}"  # backward-compat alias
SKIP_REPEAT="${SKIP_REPEAT:-false}"

# --- Output paths ---
OUTDIR="output/${RUN_NAME}"
WINDOWS_OUT="${OUTDIR}/cn_w500.bed"
REPEAT_OUT="${OUTDIR}/repeat_annotated_w500.bed"
SEGS_OUT="${OUTDIR}/segs_cnacc_w500.bed"
VALID_DIR="${OUTDIR}/validation"

# RM auto-detect from current output dir
if [ -z "${RM_ANNOTATION_BED:-}" ] && [ -f "${REPEAT_OUT}" ]; then
    RM_ANNOTATION_BED="${REPEAT_OUT}"
fi

mkdir -p "${OUTDIR}" "${VALID_DIR}" logs

# --- chrX correction flag ---
CHRX_FLAG=""
if [ "${SEX}" = "XX" ]; then
    CHRX_FLAG="--no-chrx-correction"
fi

echo "=================================================="
echo "KonuSeg — ${DATA_SOURCE} k=${KMER_SIZE} CN Calling Pipeline"
echo "  Input:      ${DATA_SOURCE} k${KMER_SIZE} k-mer BED"
echo "  Segmenter:  ${SEGMENTER}"
echo "  K-mer size: ${KMER_SIZE}"
echo "  Sex:        ${SEX}"
echo "  Sample:     ${SAMPLE:-[not set — GT evaluation will be skipped]}"
echo "  Validation: validate_cn_accuracy.py + evaluate_ground_truth.py"
echo "=================================================="
echo "Run:      ${RUN_NAME}"
echo "Job ID:   ${SLURM_JOB_ID:-local}"
echo "Started:  $(date)"
echo "Host:     $(hostname)"
echo "=================================================="
echo ""

# =============================================================================
# ENVIRONMENT
# =============================================================================

if command -v conda &>/dev/null; then
    eval "$(conda shell.bash hook)"
elif [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
    source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif [ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]; then
    source "${HOME}/anaconda3/etc/profile.d/conda.sh"
else
    echo "ERROR: Cannot locate conda. Ensure conda is on PATH or installed in ~/miniconda3/."
    exit 1
fi

conda activate konuseg

python3 -c "
import sys
segmenter = '${SEGMENTER}'
required = ['numpy', 'pandas']
if segmenter == 'fused_lasso':
    required += ['ruptures']
else:
    required += ['pomegranate', 'torch']
missing = []
for pkg in required:
    try: __import__(pkg)
    except ImportError: missing.append(pkg)
if missing:
    print('ERROR: missing packages:', missing)
    sys.exit(1)
print('Packages OK (' + segmenter + ')')
"

echo ""

# =============================================================================
# PRE-FLIGHT
# =============================================================================

if [ ! -f "${ONT_KMERS}" ]; then
    if [ "${SKIP_PREPROCESS}" = "true" ]; then
        echo "WARNING: K-mer BED not found: ${ONT_KMERS}"
        echo "         Continuing because SKIP_PREPROCESS=true"
    else
        echo "ERROR: K-mer BED not found: ${ONT_KMERS}"
        exit 1
    fi
    echo "K-mer BED:  [not found — skipping line count]"
else
    echo "K-mer BED:  $(wc -l < "${ONT_KMERS}") lines"
fi
echo "GC BED:     ${GC_BED:-[not set — running without GC calibration]}"
echo "RM BED:     ${RM_ANNOTATION_BED:-[not set — Step 1 will use global threshold]}"
echo "RM .out:    ${REPEATMASKER:-[not set — skip repeat annotation]}"
echo ""

# =============================================================================
# STEP 1: Preprocess k-mer data → window BED
# =============================================================================

echo "=================================================="
echo "Step 1: K-mer preprocessing → ${OUTDIR}/cn_w500.bed"
echo "  Sex: ${SEX}"
echo "=================================================="
echo ""

# If PREPROCESS_SRC is set, symlink preprocessed BED into this run's output dir.
# Useful for sweep jobs that share a single preprocessed BED across runs.
if [ -n "${PREPROCESS_SRC:-}" ] && [ -f "${PREPROCESS_SRC}" ] && [ ! -f "${WINDOWS_OUT}" ]; then
    ln -sf "$(cd "$(dirname "${PREPROCESS_SRC}")" && pwd)/$(basename "${PREPROCESS_SRC}")" "${WINDOWS_OUT}"
    echo "LINK: ${WINDOWS_OUT} → ${PREPROCESS_SRC}"
fi

if [ "${SKIP_PREPROCESS}" = "true" ] && [ -f "${WINDOWS_OUT}" ]; then
    echo "SKIP: ${WINDOWS_OUT} exists (SKIP_PREPROCESS=true)"
else
    BIO_FACTOR="${BIO_FACTOR:-150.0}"

    RM_FLAG=""
    if [ -n "${RM_ANNOTATION_BED:-}" ] && [ -f "${RM_ANNOTATION_BED}" ]; then
        RM_FLAG="--rm-annotation-bed ${RM_ANNOTATION_BED}"
        echo "  Biological filter: RM-guided (Sat=30×, default=${BIO_FACTOR}×)"
    else
        echo "  Biological filter: ${BIO_FACTOR}× Gaussian peak (global)"
    fi

    time python3 scripts/preprocess_kmer_windows.py \
        --input                  "${ONT_KMERS}" \
        --output                 "${WINDOWS_OUT}" \
        --window-size            500 \
        --bio-threshold-factor   "${BIO_FACTOR}" \
        --sex                    "${SEX}" \
        ${RM_FLAG}
    echo ""
    echo "Windows written: ${WINDOWS_OUT}"
fi

echo ""
echo "Window count: $(grep -c -v '^#' "${WINDOWS_OUT}")"
echo ""

# =============================================================================
# STEP 1.5: RepeatMasker annotation
# =============================================================================

echo "=================================================="
echo "Step 1.5: RepeatMasker Window Annotation"
echo "=================================================="
echo ""

if [ -z "${REPEATMASKER}" ]; then
    echo "SKIP: REPEATMASKER not set — no repeat annotation"
    REPEAT_OUT=""
elif [ "${SKIP_REPEAT}" = "true" ] && [ -f "${REPEAT_OUT}" ]; then
    echo "SKIP: ${REPEAT_OUT} exists (SKIP_REPEAT=true)"
elif [ ! -f "${REPEATMASKER}" ]; then
    echo "WARNING: RepeatMasker .out not found: ${REPEATMASKER}"
    echo "         Skipping repeat annotation"
    REPEAT_OUT=""
else
    time python3 scripts/compute_repeat_annotation.py \
        --repeatmasker "${REPEATMASKER}" \
        --windows      "${WINDOWS_OUT}" \
        --output       "${REPEAT_OUT}" \
        --window-size  500
    echo ""
    echo "Repeat annotated windows: ${REPEAT_OUT}"
fi
echo ""

# =============================================================================
# STEP 2: Segmentation
# =============================================================================

echo "=================================================="
echo "Step 2: Segmentation — ${SEGMENTER} / cn-accuracy mode"
echo "  Sex: ${SEX} ${CHRX_FLAG:+(${CHRX_FLAG})}"
echo "=================================================="
echo ""

CV_SPLIT_THRESHOLD="${CV_SPLIT_THRESHOLD:-0}"
if [ "${SEGMENTER}" = "fused_lasso" ]; then
    MERGE_CN_TOLERANCE="${MERGE_CN_TOLERANCE:-2.0}"
    LOWDUP_THRESHOLD="${LOWDUP_THRESHOLD:-1.5}"
    MIN_SEGMENT_LENGTH="${MIN_SEGMENT_LENGTH:-3000}"
else
    MERGE_CN_TOLERANCE="${MERGE_CN_TOLERANCE:-0.5}"
    LOWDUP_THRESHOLD="${LOWDUP_THRESHOLD:-1.87}"
    MIN_SEGMENT_LENGTH="${MIN_SEGMENT_LENGTH:-10000}"
fi
# Penalty factor: k-size-aware default via inverse formula.
# Smaller k-mers → higher noise variance → need stronger penalty to avoid
# excessive changepoints and O(n·K) runtime explosion in PELT.
# Formula: pf(k) = 288/k  (calibration: k=72→4.0, k=32→9.0)
# NOTE: k=32 penalty sweep (pf=4,6,8,10,12) showed IDENTICAL GT CN values.
# For k≤50, penalty factor is irrelevant — the problem is k-mer non-uniqueness
# in the bloom filter data, not segmentation parameters. See docs/roadmap.md §5.
if [ -z "${PENALTY_FACTOR:-}" ]; then
    PENALTY_FACTOR=$(python3 -c "print(round(288.0 / ${KMER_SIZE}, 1))")
fi
echo "  Segmenter:          ${SEGMENTER}"
echo "  CV split threshold: ${CV_SPLIT_THRESHOLD}"
echo "  Merge CN tolerance: ${MERGE_CN_TOLERANCE}"
echo "  LowDup threshold:   ${LOWDUP_THRESHOLD}"
echo "  Penalty factor:     ${PENALTY_FACTOR}"
echo "  Min segment len:    ${MIN_SEGMENT_LENGTH}"

if [ "${SKIP_HMM}" = "true" ] && [ -f "${SEGS_OUT}" ]; then
    echo "SKIP: ${SEGS_OUT} exists (SKIP_HMM/SKIP_SEGMENTER=true)"
else
    GC_FLAG=""
    if [ -n "${GC_BED}" ]; then
        GC_FLAG="--gc-content-bed ${GC_BED}"
        echo "  GC calibration: ENABLED"
    else
        echo "  GC calibration: DISABLED"
    fi

    REPEAT_FLAG=""
    if [ -n "${REPEAT_OUT}" ] && [ -f "${REPEAT_OUT}" ]; then
        REPEAT_FLAG="--repeat-bed ${REPEAT_OUT}"
        echo "  Repeat annotation: ENABLED"
    else
        echo "  Repeat annotation: DISABLED"
    fi
    echo ""

    MIN_KMERS="${MIN_KMERS:-30}"
    echo "  Min k-mer coverage filter: ${MIN_KMERS}"

    if [ "${SEGMENTER}" = "fused_lasso" ]; then
        time python3 scripts/segment_cnv_fused_lasso.py \
            --input                "${WINDOWS_OUT}" \
            --output               "${SEGS_OUT}" \
            --extended \
            ${CHRX_FLAG} \
            --cv-split-threshold   "${CV_SPLIT_THRESHOLD}" \
            --merge-cn-tolerance   "${MERGE_CN_TOLERANCE}" \
            --lowdup-threshold     "${LOWDUP_THRESHOLD}" \
            --penalty-factor       "${PENALTY_FACTOR}" \
            --min-segment-length   "${MIN_SEGMENT_LENGTH}" \
            --min-kmers            "${MIN_KMERS}" \
            ${GC_FLAG} \
            ${REPEAT_FLAG}
    else
        time python3 scripts/segment_cnv_hmm_log_7state.py \
            --input                "${WINDOWS_OUT}" \
            --output               "${SEGS_OUT}" \
            --mode                 cn-accuracy \
            --extended \
            ${CHRX_FLAG} \
            --cv-split-threshold   "${CV_SPLIT_THRESHOLD}" \
            --cn-reclassify-threshold "${LOWDUP_THRESHOLD}" \
            --min-kmers            "${MIN_KMERS}" \
            ${GC_FLAG} \
            ${REPEAT_FLAG}
    fi

    echo ""
    echo "Segments written: ${SEGS_OUT}"
fi

echo ""
echo "Segment count: $(grep -c -v '^#' "${SEGS_OUT}")"
echo ""

# =============================================================================
# STEP 3: CN Accuracy Validation
# =============================================================================

echo "=================================================="
echo "Step 3: CN Accuracy Validation"
echo "=================================================="
echo ""

SAMPLE_FLAG=""
if [ -n "${SAMPLE}" ]; then
    SAMPLE_FLAG="--sample ${SAMPLE}"
fi

time python3 scripts/validate_cn_accuracy.py \
    --segments "${SEGS_OUT}" \
    --output   "${VALID_DIR}/cn_accuracy_report.txt" \
    ${SAMPLE_FLAG} \
    2>&1 | tee "${VALID_DIR}/cn_accuracy_report.txt"

echo ""

# =============================================================================
# STEP 3b: Ground Truth Evaluation (CHM13-specific)
# =============================================================================

echo "=================================================="
echo "Step 3b: Ground Truth Locus Evaluation"
echo "=================================================="
echo ""

if [ -z "${SAMPLE}" ]; then
    echo "SKIP: SAMPLE not set — GT evaluation requires a known sample (chm13, hg002)."
    echo "      Set SAMPLE=chm13 or SAMPLE=hg002 to enable GT evaluation."
else
    echo "  Sample: ${SAMPLE}"
    python3 scripts/evaluate_ground_truth.py \
        --segments "${SEGS_OUT}" \
        --output   "${VALID_DIR}/gt_evaluation_report.md" \
        --kmer-size "${KMER_SIZE}" \
        --sample "${SAMPLE}" \
        2>&1 || true   # non-zero exit should not abort the pipeline
fi

echo ""

# =============================================================================
# SUMMARY
# =============================================================================

echo "=================================================="
echo "${DATA_SOURCE} k=${KMER_SIZE} Pipeline Complete!"
echo "=================================================="
echo ""
echo "Run:       ${RUN_NAME}"
echo "Job ID:    ${SLURM_JOB_ID:-local}"
echo "Completed: $(date)"
echo ""
echo "Output files:"
printf "  %-38s %s\n" "Preprocessed windows:" "${WINDOWS_OUT}"
printf "  %-38s %s\n" "CN-accuracy segments:" "${SEGS_OUT}"
printf "  %-38s %s\n" "CN accuracy report:" "${VALID_DIR}/cn_accuracy_report.txt"
printf "  %-38s %s\n" "GT evaluation report:" "${VALID_DIR}/gt_evaluation_report.md"
echo ""
echo "Target metrics (k=${KMER_SIZE}):"
echo "  State concordance:   > 80%"
echo "  GT expected CN values are k-size dependent (n_copies × identity^k)"
echo "  TP53 CN:             ~1.0  (single-copy control, k-size invariant)"
echo "  See gt_evaluation_report.md for per-locus expected CN at k=${KMER_SIZE}"
echo "=================================================="

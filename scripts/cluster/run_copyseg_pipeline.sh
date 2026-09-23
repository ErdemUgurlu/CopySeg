#!/bin/bash
#SBATCH --job-name=copyseg
#SBATCH --partition=hi_end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/copyseg_%j.out
#SBATCH --error=logs/copyseg_%j.err

set -euo pipefail
export PYTHONUNBUFFERED=1   # Python print() output flushed immediately


RUN_NAME="${1:-run_$(date +%Y%m%d_%H%M)}"

if [ -z "${INPUT_KMERS:-}" ]; then
    echo "ERROR: INPUT_KMERS must be set (path to k-mer BED file)"
    echo "  Usage: INPUT_KMERS=/path/to/kmers.bed sbatch scripts/cluster/run_copyseg_pipeline.sh my_run"
    exit 1
fi
ONT_KMERS="${INPUT_KMERS}"
REPEATMASKER="${REPEATMASKER:-}"

SEX="${SEX:-XX}"

SAMPLE="${SAMPLE:-}"

if [ -z "${KMER_SIZE:-}" ]; then
    if [[ "${ONT_KMERS}" =~ _k([0-9]+)_ ]]; then
        KMER_SIZE="${BASH_REMATCH[1]}"
    else
        echo "WARNING: Could not auto-detect k-mer size from filename."
        echo "         Defaulting to KMER_SIZE=72. Set KMER_SIZE explicitly if different."
        KMER_SIZE=72
    fi
fi

if [[ "${ONT_KMERS}" == *"pacbio"* ]]; then
    DATA_SOURCE="PacBio"
else
    DATA_SOURCE="ONT"
fi

GC_BED="${GC_BED:-}"


SEGMENTER="${SEGMENTER:-hmm}"

SKIP_PREPROCESS="${SKIP_PREPROCESS:-false}"
SKIP_SEGMENTER="${SKIP_SEGMENTER:-false}"
SKIP_HMM="${SKIP_HMM:-${SKIP_SEGMENTER}}"  # backward-compat alias
SKIP_REPEAT="${SKIP_REPEAT:-false}"

OUTDIR="output/${RUN_NAME}"
WINDOWS_OUT="${OUTDIR}/cn_w500.bed"
REPEAT_OUT="${OUTDIR}/repeat_annotated_w500.bed"
SEGS_OUT="${OUTDIR}/segs_cnacc_w500.bed"
VALID_DIR="${OUTDIR}/validation"

if [ -z "${RM_ANNOTATION_BED:-}" ] && [ -f "${REPEAT_OUT}" ]; then
    RM_ANNOTATION_BED="${REPEAT_OUT}"
fi

mkdir -p "${OUTDIR}" "${VALID_DIR}" logs

CHRX_FLAG=""
if [ "${SEX}" = "XX" ]; then
    CHRX_FLAG="--no-chrx-correction"
fi

echo "=================================================="
echo "CopySeg — ${DATA_SOURCE} k=${KMER_SIZE} CN Calling Pipeline"
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

conda activate "${CONDA_ENV:-copyseg}"

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
echo "RM mask:   ${RM_MASK_DIR:-[not set — no k-mer down-weighting]}"
echo ""


REF_FAI="${REF_FAI:-}"
REF_FASTA="${REF_FASTA:-}"
if [ -n "${RM_MASK_DIR:-}" ] && [ -d "${RM_MASK_DIR}" ]; then
    echo "=================================================="
    echo "Step 0.5: RM binary mask — pre-computed"
    echo "  Using: ${RM_MASK_DIR}"
    echo "=================================================="
elif [ -n "${RM_MASK_DIR:-}" ] && [ ! -d "${RM_MASK_DIR}" ] && [ -n "${REPEATMASKER}" ] && { [ -n "${REF_FASTA}" ] || [ -n "${REF_FAI}" ]; }; then
    if [ -n "${REF_FASTA}" ]; then
        REF_FLAG="--ref-fasta ${REF_FASTA}"
        REF_DISPLAY="${REF_FASTA}"
    else
        REF_FLAG="--ref-fai ${REF_FAI}"
        REF_DISPLAY="${REF_FAI}"
    fi
    echo "=================================================="
    echo "Step 0.5: Pre-computing RM binary mask"
    echo "  RepeatMasker: ${REPEATMASKER}"
    echo "  Reference:    ${REF_DISPLAY}"
    echo "  Output:       ${RM_MASK_DIR}"
    echo "=================================================="
    time python3 scripts/pipeline/compute_rm_mask.py \
        --repeatmasker "${REPEATMASKER}" \
        ${REF_FLAG} \
        --output-dir   "${RM_MASK_DIR}"
    echo ""
else
    echo "Step 0.5: RM binary mask — SKIP (RM_MASK_DIR not set)"
fi
echo ""


echo "=================================================="
echo "Step 1: K-mer preprocessing → ${OUTDIR}/cn_w500.bed"
echo "  Sex: ${SEX}"
echo "=================================================="
echo ""

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

    WEIGHT_FLAG=""
    if [ -n "${WEIGHT_DIR:-}" ] && [ -d "${WEIGHT_DIR}" ]; then
        WEIGHT_FLAG="--weight-dir ${WEIGHT_DIR}"
        echo "  Multiplicity correction: ${WEIGHT_DIR}"
    fi

    PW_FLAG=""
    if [ "${PER_WINDOW_CORRECT:-false}" = "true" ]; then
        PW_FLAG="--per-window-correct --pw-percentile ${PW_PERCENTILE}"
        echo "  Per-window multiplicity correction: ENABLED (p${PW_PERCENTILE})"
    fi

    RM_MASK_FLAG=""
    if [ -n "${RM_MASK_DIR:-}" ] && [ -d "${RM_MASK_DIR}" ]; then
        RM_MASK_FLAG="--rm-mask-dir ${RM_MASK_DIR}"
        if [ -n "${RM_CLASS_WEIGHTS:-}" ]; then
            RM_MASK_FLAG="${RM_MASK_FLAG} --rm-class-weights ${RM_CLASS_WEIGHTS}"
            echo "  RM per-class weights: ${RM_MASK_DIR}"
            echo "    weights: ${RM_CLASS_WEIGHTS}"
        elif [ "${RM_EXCLUDE:-false}" = "true" ]; then
            RM_MASK_FLAG="${RM_MASK_FLAG} --rm-exclude"
            echo "  RM k-mer EXCLUSION: ${RM_MASK_DIR} (repeat k-mers fully excluded)"
        else
            REPEAT_KMER_WEIGHT="${REPEAT_KMER_WEIGHT:-0.01}"
            RM_MASK_FLAG="${RM_MASK_FLAG} --repeat-kmer-weight ${REPEAT_KMER_WEIGHT}"
            echo "  RM k-mer down-weighting: ${RM_MASK_DIR} (weight=${REPEAT_KMER_WEIGHT})"
        fi
    fi

    time python3 scripts/pipeline/preprocess_kmer_windows.py \
        --input                  "${ONT_KMERS}" \
        --output                 "${WINDOWS_OUT}" \
        --window-size            500 \
        --bio-threshold-factor   "${BIO_FACTOR}" \
        --sex                    "${SEX}" \
        ${RM_FLAG} \
        ${WEIGHT_FLAG} \
        ${PW_FLAG} \
        ${RM_MASK_FLAG}
    echo ""
    echo "Windows written: ${WINDOWS_OUT}"
fi

echo ""
echo "Window count: $(grep -c -v '^#' "${WINDOWS_OUT}")"
echo ""


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
    time python3 scripts/pipeline/compute_repeat_annotation.py \
        --repeatmasker "${REPEATMASKER}" \
        --windows      "${WINDOWS_OUT}" \
        --output       "${REPEAT_OUT}" \
        --window-size  500
    echo ""
    echo "Repeat annotated windows: ${REPEAT_OUT}"
fi
echo ""


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
if [ -z "${PENALTY_FACTOR:-}" ]; then
    PENALTY_FACTOR=$(python3 -c "print(round(288.0 / ${KMER_SIZE}, 1))")
fi
NOISE_VAR_FLOOR="${NOISE_VAR_FLOOR:-0.02}"
WEIGHT_DIR="${WEIGHT_DIR:-}"
PER_WINDOW_CORRECT="${PER_WINDOW_CORRECT:-false}"
PW_PERCENTILE="${PW_PERCENTILE:-25}"
echo "  Segmenter:          ${SEGMENTER}"
echo "  CV split threshold: ${CV_SPLIT_THRESHOLD}"
echo "  Merge CN tolerance: ${MERGE_CN_TOLERANCE}"
echo "  LowDup threshold:   ${LOWDUP_THRESHOLD}"
echo "  Penalty factor:     ${PENALTY_FACTOR}"
echo "  Noise var floor:   ${NOISE_VAR_FLOOR}"
if [ -n "${WEIGHT_DIR}" ]; then
    echo "  Multiplicity wts:  ${WEIGHT_DIR}"
fi
if [ "${PER_WINDOW_CORRECT}" = "true" ]; then
    echo "  Per-window correct: ENABLED (p${PW_PERCENTILE})"
fi
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
        time python3 scripts/pipeline/segment_cnv_fused_lasso.py \
            --input                "${WINDOWS_OUT}" \
            --output               "${SEGS_OUT}" \
            --extended \
            ${CHRX_FLAG} \
            --cv-split-threshold   "${CV_SPLIT_THRESHOLD}" \
            --merge-cn-tolerance   "${MERGE_CN_TOLERANCE}" \
            --lowdup-threshold     "${LOWDUP_THRESHOLD}" \
            --penalty-factor       "${PENALTY_FACTOR}" \
            --noise-var-floor      "${NOISE_VAR_FLOOR}" \
            --min-segment-length   "${MIN_SEGMENT_LENGTH}" \
            --min-kmers            "${MIN_KMERS}" \
            ${GC_FLAG} \
            ${REPEAT_FLAG}
    else
        time python3 scripts/pipeline/segment_cnv_hmm_log_7state.py \
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

GT_CN_COLUMN="cn_median"
REFINED_SEGS_OUT="${OUTDIR}/segs_cnacc_w500_refined.bed"
if [ -n "${RM_MASK_DIR:-}" ] \
   && [ "${KMER_SIZE}" -le 50 ] \
   && [ "${SKIP_PASS2:-false}" != "true" ]; then
    echo "=================================================="
    echo "Step 2.5: Pass-2 Unique-Only CN Refinement"
    echo "=================================================="
    echo "  Trigger: RM_MASK_DIR set, KMER_SIZE=${KMER_SIZE}<=50"
    echo "  Input:   ${SEGS_OUT} + ${WINDOWS_OUT}"
    echo "  Output:  ${REFINED_SEGS_OUT}"
    echo ""

    PASS2_SEX_FLAG=""
    if [ -n "${SEX:-}" ]; then
        PASS2_SEX_FLAG="--sex ${SEX}"
    fi

    if time python3 scripts/pipeline/refine_cn_unique.py \
            --segments "${SEGS_OUT}" \
            --windows  "${WINDOWS_OUT}" \
            --output   "${REFINED_SEGS_OUT}" \
            ${PASS2_SEX_FLAG} \
            --lowdup-threshold "${LOWDUP_THRESHOLD:-1.25}" ; then
        SEGS_OUT="${REFINED_SEGS_OUT}"
        GT_CN_COLUMN="cn_refined"
        echo ""
        echo "Pass-2 refinement applied. Downstream Step 3/3b will score against"
        echo "cn_refined. Pass-1 cn_median preserved in the same file for diff."
    else
        echo "WARNING: Pass-2 refinement failed (exit $?). Falling back to Pass-1"
        echo "         segments for Step 3/3b. Investigate refine_cn_unique.py logs."
    fi
    echo ""
else
    if [ -z "${RM_MASK_DIR:-}" ]; then
        echo "Step 2.5: SKIP (no RM_MASK_DIR — Pass-1 cn_median used)"
    elif [ "${KMER_SIZE}" -gt 50 ]; then
        echo "Step 2.5: SKIP (KMER_SIZE=${KMER_SIZE}>50 — Pass-2 is for small-k)"
    else
        echo "Step 2.5: SKIP (SKIP_PASS2=true)"
    fi
    echo ""
fi


echo "=================================================="
echo "Step 3: CN Accuracy Validation"
echo "=================================================="
echo ""

SAMPLE_FLAG=""
if [ -n "${SAMPLE}" ]; then
    SAMPLE_FLAG="--sample ${SAMPLE}"
fi

time python3 scripts/pipeline/validate_cn_accuracy.py \
    --segments "${SEGS_OUT}" \
    --output   "${VALID_DIR}/cn_accuracy_report.txt" \
    ${SAMPLE_FLAG} \
    2>&1 | tee "${VALID_DIR}/cn_accuracy_report.txt"

echo ""


echo "=================================================="
echo "Step 3b: Ground Truth Locus Evaluation"
echo "=================================================="
echo ""

if [ -z "${SAMPLE}" ]; then
    echo "SKIP: SAMPLE not set — GT evaluation requires a known sample (chm13, hg002)."
    echo "      Set SAMPLE=chm13 or SAMPLE=hg002 to enable GT evaluation."
else
    echo "  Sample: ${SAMPLE}"
    WINDOWS_FLAG=""
    if [ -f "${WINDOWS_OUT}" ] && [ "${GT_CN_COLUMN:-cn_median}" = "cn_median" ]; then
        WINDOWS_FLAG="--windows ${WINDOWS_OUT}"
        echo "  Window-level CN: ENABLED"
    elif [ "${GT_CN_COLUMN:-cn_median}" = "cn_refined" ]; then
        echo "  Window-level CN: DISABLED (Pass-2 mode — verdict uses segment cn_refined)"
    else
        echo "  Window-level CN: DISABLED (no windows file)"
    fi
    echo "  CN column: ${GT_CN_COLUMN:-cn_median}"
    python3 scripts/pipeline/evaluate_ground_truth.py \
        --segments "${SEGS_OUT}" \
        --output   "${VALID_DIR}/gt_evaluation_report.md" \
        --kmer-size "${KMER_SIZE}" \
        --sample "${SAMPLE}" \
        --cn-column "${GT_CN_COLUMN:-cn_median}" \
        ${WINDOWS_FLAG} \
        2>&1 || true   # non-zero exit should not abort the pipeline
fi

echo ""


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

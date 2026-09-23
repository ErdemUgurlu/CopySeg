# CopySeg — Copy-Number Calling and Segmentation from k-mer Counts

CopySeg segments a genome and assigns **copy number (CN)** to each segment directly from
per-position k-mer counts. Output: BED segments with reliable CN values.

This is a **CN-calling** problem, not a segmental-duplication detection problem. SD-overlap
metrics (SEDEF/BISER F1) are used only as cross-tool validation, not as the primary objective.

## Samples & data

| | CHM13 | HG002 |
|---|---|---|
| Type | Hydatidiform mole (self-validation) | Ashkenazi trio child (cross-validation) |
| Sex | XX | XY |
| Reference | T2T-CHM13v2.0 | HG002 maternal assembly (CM039 accessions) |
| Sequencing | ONT + PacBio, k=72 and k=32 | ONT + PacBio, k=72 and k=32 |
| Input | k-mer count BED (4-col), 500 bp windows | same |

`CN=1.0` = diploid-normal baseline (mole = effectively homozygous diploid; both homologs
collapse to one haploid assembly position). Expected CN values are **not** halved for
"haploid correction" — the diploid-peak normalization already accounts for both homologs.

**k-mer size:** k≥72 is the gold standard (`count ∝ CN`). At k≤50, k-mers fall entirely
inside repeat elements, so `count = coverage × CN × multiplicity` and counts inflate.
A two-pass **unique-only refinement** (Pass-2) extends the pipeline down to k=32.

## Pipeline

```
k-mer BED ─▶ preprocess_kmer_windows.py        # 500bp windows, neutral-band norm, RM-weighting
          ─▶ compute_repeat_annotation.py      # per-window repeat class (RepeatMasker)
          ─▶ segment_cnv_fused_lasso.py (PELT) | segment_cnv_hmm_log_7state.py (HMM)
          ─▶ [k≤50] refine_cn_unique.py        # Pass-2: re-estimate CN from unique positions
          ─▶ validate_cn_accuracy.py + evaluate_ground_truth.py
```

7 CN states (cn-accuracy mode): Neutral, LowDup, HighDup, Amp, MedAmp, HighAmp, ExtremeAmp.
Pass-2 auto-triggers when `RM_MASK_DIR` is set **and** `KMER_SIZE ≤ 50`, writing
`segs_cnacc_w500_refined.bed` (adds `cn_refined`, `state_refined`, `refine_method`).

## Running (SLURM cluster)

```bash
# CHM13 (XX) — k=72, PELT
INPUT_KMERS=/path/chm13_ont_k72.bed SEGMENTER=fused_lasso \
  sbatch scripts/cluster/run_chm13_cnacc_job.sh ont_k72

# CHM13 — k=32 + Pass-2 (RM mask triggers unique-only refinement)
INPUT_KMERS=/path/chm13_pb_k32.bed SEGMENTER=fused_lasso \
  RM_MASK_DIR=output/rm_mask_chm13 sbatch scripts/cluster/run_chm13_cnacc_job.sh pb_k32

# HG002 (XY)
INPUT_KMERS=/path/hg002_pb_k72.bed SEGMENTER=fused_lasso \
  sbatch scripts/cluster/run_hg002_cnacc_job.sh pb_k72

# Generic entry point
INPUT_KMERS=/path/kmers.bed SEX=XY SAMPLE=hg002 SEGMENTER=fused_lasso \
  sbatch scripts/cluster/run_copyseg_pipeline.sh my_run
```

Defaults: `PENALTY_FACTOR = 288/k` (k72→4.0, k32→9.0), `NOISE_VAR_FLOOR = 0.02`.

## Key scripts

| Script | Purpose |
|---|---|
| `preprocess_kmer_windows.py` | k-mer BED → 500 bp windows with CN (any sample/tech/k) |
| `segment_cnv_fused_lasso.py` | Weighted PELT segmenter (primary) |
| `segment_cnv_hmm_log_7state.py` | 7-state HMM segmenter (alternative) |
| `refine_cn_unique.py` | Pass-2 unique-only CN refinement (k≤50) |
| `compute_rm_mask.py` | RepeatMasker → per-base mask (.npy) for k-mer down-weighting |
| `compute_repeat_annotation.py` | RepeatMasker → per-window repeat class |
| `evaluate_ground_truth.py` | GT-locus evaluation (k-size aware, per-sample) |
| `validate_cn_accuracy.py` | CN accuracy metrics |
| `run_copyseg_pipeline.sh` | Generic SLURM pipeline; `run_{chm13,hg002}_cnacc_job.sh` wrappers |

## Validation

- **Ground truth:** per-sample GT loci (CHM13: 8 scored, e.g. TP53=1, AMY1=7, TBC1D3=9;
  HG002: 6 scored). Verdicts: PASS ≤35% error, MARGINAL ≤50%, FAIL >50%.
- **Genome-wide repeat ratio:** healthy runs land at ~40–55% non-neutral (human genome is
  ~54% repetitive).
- **Cross-tool:** SEDEF (CHM13) and BISER (HG002) SD-overlap sensitivity; per-locus MEGABLAST
  paralog counts (alignment-based, independent of k-mer methods).

## Repository layout

Large data and outputs are git-ignored (see `.gitignore`).

```
scripts/
  pipeline/        end-to-end CN caller: preprocess → segment → refine → validate
  cluster/         SLURM runners
  chrom_utils.py   shared chromosome-naming utils
tests/             pytest (in-process logic + integration)
```

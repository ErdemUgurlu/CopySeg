#!/usr/bin/env python3
"""
CopySeg Ground Truth Evaluation Script
=======================================
Evaluates CopySeg CN estimates against AGGREGATE K-MER CN ground truth.

Key principle: k=72 k-mer callers measure aggregate depth across all
near-identical paralogs genome-wide. Expected CN = effective aggregate at
k=72, computed as sum(identity_i^72) across all paralog copies.
Formula: expected_agg_CN ≈ n_copies × mean_pairwise_identity^k
For k=72, Jellyfish-validated values are used. For other k, computed dynamically.

GT value revision history:
  2026-03-27: LPA_KIV2 22→10 (k72 formula, ">99%" identity);
              TBC1D3 9→2 (Sudmant et al. 2010, ~97.9% identity);
              SMN1/2 2→3.0 (SERF1A/B/NAIP inversion block cross-talk);
              rDNA → observatory mode (211/5 formula physically ambiguous);
              per-locus tolerance thresholds added.
  2026-03-29: rDNA 42→219 genome-wide aggregate (Nurk et al. 2022; 219 units, not 211);
              rDNA coords updated per T2T-Ref GitHub; per-NOR unit counts added;
              SMN1/2 3.0→2.0 (SERF1A/B/NAIP contribution unconfirmed; Rochette et al. 2001);
              AMY1/AMY2A: GRCh38 coordinate warning added;
              +TP53 (single-copy control), +NOTCH2NL (SD, CN≈4), +5S_rDNA (observatory).
  2026-03-29b: ASSEMBLY-VALIDATED via Jellyfish k=72 on CHM13v2.0:
              AMY1 coords → CHM13v2.0 (103504385-103513421); CN=7 confirmed;
              AMY2A expected_cn 1→4 (AMY1/AMY2B exonic k-mer sharing);
              SMN1/2 CN=2 confirmed, promoted observatory→normal, tol 0.40→0.25;
              LPA 10→9 (assembly-measured); TBC1D3 2→9 (identity ~99.9%+, gene conversion);
              5S_rDNA 30→117 (actual identity 99.87%, not 98%);
              TP53=1, NOTCH2NL=4, rDNA=219 all confirmed.

Usage:
    python3 scripts/evaluate_ground_truth.py \
        --segments output/chm13_my_clean_v11/segs_chm13_cnacc_w500.bed \
        --output   output/chm13_my_clean_v11/validation/gt_evaluation_report.md
"""

import argparse
import sys
import os
import pandas as pd
import numpy as np
from datetime import datetime

# GROUND TRUTH TABLES
# ---------------------------------------------------------------------------
# Tuple: (chr, start, end, gene_region, expected_agg_cn, tolerance,
#         verdict_mode, category, notes, n_copies, identity)
#
#   expected_agg_cn : Assembly-validated CN at the sample's validated k-mer size
#   tolerance       : PASS threshold as fraction (0.20 = ±20%); None if observatory/skip
#   verdict_mode    : "normal"      → apply tolerance → PASS/MARGINAL/FAIL
#                     "observatory" → report CN only, no pass/fail
#                     "skip"        → omit from evaluation
#   n_copies        : physical copy count (None if not applicable)
#   identity        : mean pairwise nucleotide identity (None if not applicable)
#
# For k == validated_k: expected_cn used directly (most accurate).
# For k != validated_k: expected_cn = n_copies × identity^k (computed dynamically).
# Entries with n_copies=None fall back to validated value at any k.
#
# Coordinates: T2T-CHM13v2.0 (chr prefix). Both CHM13 and HG002 use this reference.
# ---------------------------------------------------------------------------

# ═══════════════════════════════════════════════════════════════════════════
# CHM13 — validated at k=72 (Jellyfish on CHM13v2.0 assembly)
# ═══════════════════════════════════════════════════════════════════════════
_CHM13_VALIDATED_K = 72
_CHM13_RDNA_AGG_CN = 219  # Nurk et al. 2022 Science

GROUND_TRUTH_CHM13 = [
    # AMY1 Cluster (Salivary Amylase)
    # CHM13 has exactly 7 functional AMY1 copies; >99.9% identical.
    # Jellyfish k=72 assembly-validated: measured CN=7.0, IQR [7,7].
    ("chr1",  103504385, 103513421, "AMY1_Cluster",    7,   0.20, "normal",
     "Segmental Dup",
     "AMY1 cluster; 7 copies >99.9% identical; k72 agg CN=7; "
     "ASSEMBLY-VALIDATED: Jellyfish k=72 measured 7.0, IQR [7,7]; "
     "CHM13v2.0 coordinates; Source: Yilmaz et al. 2024 Science; Bolognini et al. 2024 Nature",
     7, 0.999),

    # AMY2A (Pancreatic Amylase)
    # Single physical copy but shares conserved exonic k-mers with AMY1/AMY2B.
    # Jellyfish k=72 assembly-validated: measured CN=4.0, IQR [3,4].
    ("chr1",  103466373, 103474647, "AMY2A",           4,   0.25, "observatory",
     "Segmental Dup",
     "AMY2A; 1 physical copy but AMY1/AMY2B conserved exonic k-mer sharing → aggregate CN=4; "
     "ASSEMBLY-VALIDATED: Jellyfish k=72 measured 4.0, IQR [3,4]; "
     "OLD VALUE 1 WAS WRONG (ignored cross-family k-mer sharing); "
     "OBSERVATORY: 8.3kb locus engulfed by AMY1 segment; 500bp windows cannot separate",
     None, None),

    # SMN Family (Survival Motor Neuron)
    # SMN1+SMN2: 2 copies × ~99.9% identity → k72 CN=2.
    # Jellyfish k=72 assembly-validated: measured CN=2.0, IQR [2,2].
    # SERF1A/B/NAIP contribution: NONE (confirmed by assembly measurement).
    # Mean=49.3 is inflated by Alu/LINE repetitive elements in region,
    # but median is unaffected. Promoted from observatory to normal.
    ("chr5",   71381729,  71423141, "SMN1",            2.0, 0.25, "normal",
     "Segmental Dup",
     "SMN1+SMN2: 2 copies; ASSEMBLY-VALIDATED: Jellyfish k=72 measured 2.0, IQR [2,2]; "
     "SERF1A/B/NAIP k-mer contribution = NONE (assembly-confirmed); "
     "mean=49.3 inflated by Alu/LINE repetitive elements but median unaffected; "
     "OLD VALUES WRONG: 3.0 (SERF1A/B/NAIP hypothesis) and observatory mode; "
     "Source: Rochette et al. 2001 Hum Genet; Zwartkruis et al. 2025",
     2, 0.999),
    ("chr5",   70791126,  70837821, "SMN2",            2.0, 0.25, "normal",
     "Segmental Dup",
     "SMN2 coord; same k-mer pool as SMN1; ASSEMBLY-VALIDATED: Jellyfish k=72 measured 2.0, IQR [2,2]; "
     "promoted from observatory to normal; Source: Rochette et al. 2001; Zwartkruis et al. 2025",
     2, 0.999),

    # LPA KIV-2 (Kringle IV Type 2 VNTR)
    # Jellyfish k=72 assembly-validated: measured CN=9.0, IQR [1,23].
    # Wide IQR because GT window covers full LPA gene including flanking unique regions.
    ("chr6",  161783172, 162011762, "LPA_KIV2_Array",  9,  0.35, "normal",
     "VNTR Array",
     "KIV-2 VNTR; 22 units (T2T-resolved); ASSEMBLY-VALIDATED: Jellyfish k=72 measured 9.0, IQR [1,23]; "
     "wide IQR because coordinate range covers full LPA gene — flanking unique regions (CN=1) pull median down; "
     "prev=10 (from formula 22×0.99^72≈10.6); revised to assembly-measured 9",
     22, 0.99),

    # rDNA arrays (45S ribosomal DNA) — OBSERVATORY MODE
    # Formula "211/5=42" physically wrong unless inter-NOR identity <90% at k=72.
    # rDNA is highly conserved; inter-NOR identity likely >95%.
    # True expected CN: 42×(1 + 4×identity^72) ≈ 47-200 (range, identity unknown).
    # → No pass/fail verdict; CN reported for monitoring only.
    # (known_issues.md D1)
    ("chr13",   5817416,   9348041, "rDNA_chr13",    219, None, "observatory",
     "rDNA Array",
     "45S rDNA; NOR chr13; 76-78 units (largest array); "
     "aggregate CN reflects all 219 units across 5 NORs due to inter-NOR k-mer sharing at k=72 "
     "(inter-NOR identity 99.4-99.7%); Source: Nurk et al. 2022 Science; T2T-Ref GitHub",
     None, None),
    ("chr14",   2099537,   2817811, "rDNA_chr14",    219, None, "observatory",
     "rDNA Array",
     "45S rDNA; NOR chr14; 16 units (smallest array); "
     "aggregate CN reflects all 219 units across 5 NORs due to inter-NOR k-mer sharing at k=72; "
     "Source: Nurk et al. 2022 Science; T2T-Ref GitHub",
     None, None),
    ("chr15",   2506442,   4707485, "rDNA_chr15",    219, None, "observatory",
     "rDNA Array",
     "45S rDNA; NOR chr15; 49-50 units; "
     "aggregate CN reflects all 219 units across 5 NORs due to inter-NOR k-mer sharing at k=72; "
     "Source: Nurk et al. 2022 Science; T2T-Ref GitHub",
     None, None),
    ("chr21",   3108298,   5612715, "rDNA_chr21",    219, None, "observatory",
     "rDNA Array",
     "45S rDNA; NOR chr21; 56 units; "
     "aggregate CN reflects all 219 units across 5 NORs due to inter-NOR k-mer sharing at k=72; "
     "Source: Nurk et al. 2022 Science; T2T-Ref GitHub",
     None, None),
    ("chr22",   4793794,   5720650, "rDNA_chr22",    219, None, "observatory",
     "rDNA Array",
     "45S rDNA; NOR chr22; 21 units; "
     "aggregate CN reflects all 219 units across 5 NORs due to inter-NOR k-mer sharing at k=72; "
     "Source: Nurk et al. 2022 Science; T2T-Ref GitHub",
     None, None),

    # TBC1D3 Family
    # 9 copies in CHM13. Jellyfish k=72 assembly-validated: measured CN=9.0, IQR [7,10].
    # This means inter-copy identity is ~99.9%+ (active gene conversion).
    # Sudmant 2010 ~97.9% identity was WRONG for this assembly — contradicted by T2T data.
    ("chr17",  39044723,  39055625, "TBC1D3_Cluster",  9,  0.25, "normal",
     "Segmental Dup",
     "TBC1D3; 9 copies; ASSEMBLY-VALIDATED: Jellyfish k=72 measured 9.0, IQR [7,10]; "
     "inter-copy identity ~99.9%+ (active gene conversion) → physical CN = k-mer CN; "
     "OLD VALUE 2 WAS WRONG: Sudmant 2010 97.9% identity contradicted by T2T assembly; "
     "97.9% was likely protein-level or allele-type metric, not nucleotide pairwise",
     9, 0.999),

    # TP53 — Single-Copy Negative Control (CN = 1.0)
    # No paralogs, no SD overlap; every 72-mer unique in CHM13 genome.
    ("chr17",   7572544,   7591594, "TP53_Control",    1.0, 0.24, "normal",
     "Single-Copy Control",
     "TP53 tumor suppressor; single-copy, no paralogs, no SD overlap; "
     "19 kbp, 13 exons; every 72-mer is unique in CHM13 genome; "
     "CHM13v2.0 NC_060941.1 complement strand; "
     "Source: NCBI RefSeq RS_2025_08 on T2T-CHM13v2.0",
     1, 1.0),

    # NOTCH2NL Family — CN ≈ 4.0
    # 5 paralogs in CHM13; gene-body identity ~99.7%; k72 agg CN = 5 × 0.997^72 ≈ 4.0
    ("chr1",  145265708, 145345897, "NOTCH2NL_Family", 4.0, 0.30, "normal",
     "Segmental Dup",
     "NOTCH2NL gene family; 5 paralogs in CHM13 (NOTCH2, NLA, NLB, NLC, NLR); "
     "gene-body identity ~99.7%; k72 agg CN = 5 × 0.997^72 ≈ 4.0; "
     "coordinates are for NOTCH2NLA (representative copy); "
     "Source: Fiddes et al. 2018 Cell; Vollger et al. 2022 Science; "
     "Eichler lab 2025 preprint bioRxiv 2025.03.14.643395",
     5, 0.997),

    # NBPF1 — CN ≈ 3 at k=72
    # NBPF gene family member 1; DUF1220 domain amplification region.
    # k72 Jellyfish assembly-validated: measured CN=3, IQR [2,10].
    # Wide IQR due to internal DUF1220 tandem repeats (CN~16 in core) vs flanking (CN~2).
    # Region median reflects mixed architecture — not a uniform duplication.
    ("chr1",   16004103,  16075615, "NBPF1",           3,    0.30, "normal",
     "Segmental Dup",
     "NBPF1; neuroblastoma breakpoint family member 1; DUF1220 domain amplification; "
     "k72 agg CN=3; heterogeneous region: core DUF1220 repeats CN~16, flanking CN~2; "
     "wide IQR [2,10] reflects internal structure; window median 3.3; "
     "k-size sensitive: k72=3, k32=10 (more DUF1220 repeats visible at shorter k)",
     None, None),

    # 5S rDNA Array — CN ≈ 117 (Observatory)
    # ~128 units; actual identity 99.87% (not 98% as initially estimated).
    # Jellyfish k=72 assembly-validated: measured CN=117.0, IQR [92,123].
    ("chr1",  227746662, 228024151, "5S_rDNA_chr1",  117,  None, "observatory",
     "Tandem Array",
     "5S rDNA composite array; ~128 units × 2240bp (AluY + 5S rRNA + dinucleotide); "
     "ASSEMBLY-VALIDATED: Jellyfish k=72 measured 117.0, IQR [92,123]; "
     "actual identity 99.87% (128 × 0.9987^72 ≈ 117); "
     "OLD VALUE 30 WAS WRONG (assumed 98% identity — actual is 99.87%); "
     "OBSERVATORY: retained due to wide IQR [92,123]; "
     "Source: Li 2024 Genome Research; Hoyt et al. 2022 Science",
     None, None),

    # TSPY Array — SKIP (chrY absent in CHM13 female)
    ("chrY",   6100000,   7031000, "TSPY_Array",      45, None, "skip",
     "Array (chrY)",
     "TSPY; 45 copies on HG002 Y chr; CHM13 ONT is female (XX) → no chrY → SKIPPED",
     None, None),
]

# ═══════════════════════════════════════════════════════════════════════════
# HG002 — validated at k=32 (Assembly-Validated)
# Reference: T2T-CHM13v2.0 (same coordinates, different copy numbers)
# ═══════════════════════════════════════════════════════════════════════════
_HG002_VALIDATED_K = 32

GROUND_TRUTH_HG002 = [
    # TP53 Control — single-copy, k-invariant
    ("chr17",  7572544,   7591594, "TP53_Control",     1,   0.24, "normal",
     "Single-Copy Control",
     "TP53; single-copy negative control; k-invariant; same as CHM13",
     1, 1.0),

    # SMN1 — 2 copies, same as CHM13
    ("chr5",  71381729,  71423141, "SMN1",             2,   0.25, "normal",
     "Segmental Dup",
     "SMN1; 2 copies, same as CHM13; validated for HG002",
     2, 1.0),

    # SMN2 — 2 copies, same as CHM13
    ("chr5",  70791126,  70837821, "SMN2",             2,   0.25, "normal",
     "Segmental Dup",
     "SMN2; 2 copies, same as CHM13; validated for HG002",
     2, 1.0),

    # AMY1 Cluster — HG002 has 3 copies (CHM13 has 7)
    ("chr1", 103504385, 103513421, "AMY1_Cluster",     3,   0.20, "normal",
     "Segmental Dup",
     "AMY1 cluster; HG002 has 3 copies (vs CHM13=7); >99.9% identical",
     3, 0.999),

    # AMY2A — observatory (k-mer sharing with AMY1/AMY2B)
    ("chr1", 103466373, 103474647, "AMY2A",            2,   0.25, "observatory",
     "Segmental Dup",
     "AMY2A; HG002-specific expected CN=2; observatory due to AMY1 k-mer sharing",
     None, None),

    # NOTCH2NL Family — ~same as CHM13 (5 paralogs), k32 effect
    ("chr1", 145265708, 145345897, "NOTCH2NL_Family",  4.5, 0.30, "normal",
     "Segmental Dup",
     "NOTCH2NL; 5 paralogs ~99.7% identical (≈CHM13); k32: 5×0.997^32≈4.5",
     5, 0.997),

    # LPA KIV-2 Array — HG002-specific copy number
    ("chr6", 161783172, 162011762, "LPA_KIV2_Array",  14,   0.35, "normal",
     "VNTR Array",
     "LPA KIV-2; HG002-specific copy number (CHM13=9 at k72); k32 validated=14",
     None, None),

    # TBC1D3 Cluster — HG002-specific copy number
    ("chr17", 39044723,  39055625, "TBC1D3_Cluster",  15,   0.25, "normal",
     "Segmental Dup",
     "TBC1D3; HG002-specific copy number (CHM13=9 at k72); k32 validated=15",
     None, None),

    # NBPF1 — HG002-specific copy number
    ("chr1",  16004103,  16075615, "NBPF1",            8,   0.30, "normal",
     "Segmental Dup",
     "NBPF1; HG002-specific k32 CN=8 (CHM13 k72=3, k32=10); "
     "DUF1220 domain amplification; heterogeneous internal structure",
     None, None),

    # 5S rDNA — observatory
    ("chr1", 227746662, 228024151, "5S_rDNA_chr1",    84,   None, "observatory",
     "Tandem Array",
     "5S rDNA; HG002 fewer copies than CHM13 (117 at k72); observatory",
     None, None),

    # rDNA — SKIP (HG002 assembly incomplete for NOR loci)
    ("chr13", 5411627,   5578556, "rDNA_chr13",   0, None, "skip", "rDNA Array",
     "rDNA; HG002 assembly incomplete → cannot validate", None, None),
    ("chr14", 1879084,   2069240, "rDNA_chr14",   0, None, "skip", "rDNA Array",
     "rDNA; HG002 assembly incomplete → cannot validate", None, None),
    ("chr15", 2325498,   2537498, "rDNA_chr15",   0, None, "skip", "rDNA Array",
     "rDNA; HG002 assembly incomplete → cannot validate", None, None),
    ("chr21", 3105606,   3341782, "rDNA_chr21",   0, None, "skip", "rDNA Array",
     "rDNA; HG002 assembly incomplete → cannot validate", None, None),
    ("chr22", 4860006,   5063528, "rDNA_chr22",   0, None, "skip", "rDNA Array",
     "rDNA; HG002 assembly incomplete → cannot validate", None, None),
]

# ═══════════════════════════════════════════════════════════════════════════
# Sample registry
# ═══════════════════════════════════════════════════════════════════════════
SAMPLE_CONFIG = {
    'chm13': {
        'gt': GROUND_TRUTH_CHM13,
        'validated_k': _CHM13_VALIDATED_K,
        'rdna_agg_cn': _CHM13_RDNA_AGG_CN,
        'label': 'CHM13 (female, XX)',
    },
    'hg002': {
        'gt': GROUND_TRUTH_HG002,
        'validated_k': _HG002_VALIDATED_K,
        'rdna_agg_cn': None,  # assembly incomplete
        'label': 'HG002 (male, XY)',
    },
}

# Legacy CHR → CM039 mapping for backward compatibility with CM039-format segments
_CHR_TO_CM039 = {
    "chr1":  "CM039011.1", "chr2":  "CM039012.1", "chr3":  "CM039013.1",
    "chr4":  "CM039014.1", "chr5":  "CM039015.1", "chr6":  "CM039016.1",
    "chr7":  "CM039017.1", "chr8":  "CM039018.1", "chr9":  "CM039019.1",
    "chr10": "CM039020.1", "chr11": "CM039021.1", "chr12": "CM039022.1",
    "chr13": "CM039023.1", "chr14": "CM039024.1", "chr15": "CM039025.1",
    "chr16": "CM039026.1", "chr17": "CM039027.1", "chr18": "CM039028.1",
    "chr19": "CM039029.1", "chr20": "CM039030.1", "chr21": "CM039031.1",
    "chr22": "CM039032.1", "chrX":  "CM039033.1",
}

# ---------------------------------------------------------------------------
# LOAD COPYSEG SEGMENTS
# ---------------------------------------------------------------------------
_FULL_SCHEMA = [
    "chrom", "start", "end", "state", "cn_median", "cn_mean",
    "n_windows", "avg_quality", "min_quality", "cn_std",
    "avg_repeats", "avg_entropy", "max_entropy",
    "masked_fraction", "repeat_class",
    # Extra columns from segment_cnv_fused_lasso.py (D7 fix: previously unnamed)
    "gc_bias_factor", "segment_iqr", "boundary_conf",
]

def load_windows(bed_path: str) -> pd.DataFrame:
    """Load preprocessed window BED (8-col: chrom start end cn ...)."""
    names = ["chrom", "start", "end", "cn", "mean_count", "log_ratio",
             "num_kmers", "num_filtered"]
    df = pd.read_csv(bed_path, sep="\t", comment="#", header=None, names=names)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    df["cn"]    = pd.to_numeric(df["cn"], errors="coerce")
    return df


def load_segments(bed_path: str) -> pd.DataFrame:
    """Load CopySeg BED. Detects column count from first data line (D7 fix)."""
    ncols = len(_FULL_SCHEMA)
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            ncols = len(line.split("\t"))
            break
    names = _FULL_SCHEMA[:ncols]
    df = pd.read_csv(bed_path, sep="\t", comment="#", header=None, names=names)
    df["start"]     = df["start"].astype(int)
    df["end"]       = df["end"].astype(int)
    df["cn_median"] = pd.to_numeric(df["cn_median"], errors="coerce")
    return df

# ---------------------------------------------------------------------------
# INTERSECT: ground truth region → overlapping CopySeg segments
# ---------------------------------------------------------------------------
def intersect_region(segs_df: pd.DataFrame, chrom_cm039: str,
                     gt_start: int, gt_end: int) -> pd.DataFrame:
    """Return all CopySeg segments overlapping [gt_start, gt_end]."""
    mask = (
        (segs_df["chrom"] == chrom_cm039) &
        (segs_df["end"]   > gt_start) &
        (segs_df["start"] < gt_end)
    )
    hits = segs_df[mask].copy()
    if hits.empty:
        return hits
    hits["ovlp_start"] = hits["start"].clip(lower=gt_start)
    hits["ovlp_end"]   = hits["end"].clip(upper=gt_end)
    hits["ovlp_len"]   = hits["ovlp_end"] - hits["ovlp_start"]
    return hits

def weighted_cn(hits: pd.DataFrame, outlier_percentile: float = 95.0) -> tuple:
    """Return (length-weighted mean CN, peak CN, total overlap bp, dominant state).

    Outlier filtering: segments with CN > p95 are excluded from the
    length-weighted mean to prevent ultra-conserved spike regions
    (e.g. SMN1/SMN2 4kb ~CN220 spikes) from dominating the estimate.
    peak_cn is always computed from unfiltered hits.
    """
    if hits.empty:
        return (None, None, 0, "NO_COVERAGE")

    total_bp  = hits["ovlp_len"].sum()
    peak_cn   = hits["cn_median"].max()
    dom_state = hits.loc[hits["ovlp_len"].idxmax(), "state"]

    # Robust estimate: exclude CN > p_outlier
    cn_threshold = np.percentile(hits["cn_median"], outlier_percentile)
    robust_hits  = hits[hits["cn_median"] <= cn_threshold]

    if robust_hits.empty or robust_hits["ovlp_len"].sum() == 0:
        # All segments are outliers → fallback to full weighted mean
        w_cn = (hits["cn_median"] * hits["ovlp_len"]).sum() / total_bp
    else:
        robust_bp = robust_hits["ovlp_len"].sum()
        w_cn = (robust_hits["cn_median"] * robust_hits["ovlp_len"]).sum() / robust_bp

    return (round(w_cn, 4), round(peak_cn, 4), int(total_bp), dom_state)

def window_cn(win_df: pd.DataFrame, chrom: str,
              gt_start: int, gt_end: int) -> float:
    """Median CN from windows overlapping [gt_start, gt_end].

    More accurate than segment-level CN because it measures only windows
    within the GT region, not the entire segment (which may extend beyond).
    """
    mask = (
        (win_df["chrom"] == chrom) &
        (win_df["end"]   > gt_start) &
        (win_df["start"] < gt_end)
    )
    hits = win_df[mask]
    if hits.empty:
        return None
    return round(float(hits["cn"].median()), 4)

# ---------------------------------------------------------------------------
# K-MER-SIZE-AWARE EXPECTED CN
# ---------------------------------------------------------------------------
def compute_expected_cn(exp_cn_validated, n_copies, identity, kmer_size,
                        validated_k=72):
    """Compute expected aggregate CN for a given k-mer size.

    For k == validated_k: returns the assembly-validated value (most accurate).
    For k != validated_k: computes n_copies * identity^k if metadata available,
                          otherwise falls back to validated value.
    """
    if kmer_size == validated_k:
        return exp_cn_validated
    if n_copies is not None and identity is not None:
        return round(n_copies * identity ** kmer_size, 1)
    return exp_cn_validated

# ---------------------------------------------------------------------------
# MAIN EVALUATION
# ---------------------------------------------------------------------------
def evaluate(segments_path: str, kmer_size: int = 72,
             sample: str = 'chm13', windows_path: str = None) -> pd.DataFrame:
    segs = load_segments(segments_path)
    win_df = load_windows(windows_path) if windows_path else None
    cfg = SAMPLE_CONFIG[sample]
    gt_table = cfg['gt']
    validated_k = cfg['validated_k']

    # Auto-detect naming convention: segments may use chr prefix or CM039
    seg_chroms = set(segs['chrom'].unique())
    _use_cm039 = any(c.startswith('CM0') for c in seg_chroms)

    # Windows may use different naming than segments — detect separately
    if win_df is not None:
        win_chroms = set(win_df['chrom'].unique())
        _win_use_cm039 = any(c.startswith('CM0') for c in win_chroms)
    else:
        _win_use_cm039 = False

    rows = []
    for (chrom, start, end, gene, exp_cn_validated, tol, vmode, category, notes,
         n_copies, identity) in gt_table:
        exp_cn = compute_expected_cn(exp_cn_validated, n_copies, identity,
                                     kmer_size, validated_k)
        # Match GT chromosome name to segment naming convention
        if _use_cm039:
            chrom_seg = _CHR_TO_CM039.get(chrom, chrom)
        else:
            chrom_seg = chrom
        # Match GT chromosome name to window naming convention
        chrom_win = _CHR_TO_CM039.get(chrom, chrom) if _win_use_cm039 else chrom

        gt_len = end - start
        hits = intersect_region(segs, chrom_seg, start, end)
        est_cn_seg, peak_cn, cov_bp, dom_state = weighted_cn(hits)

        # Window-level CN: more accurate (measures only within GT region)
        win_cn_val = None
        if win_df is not None:
            win_cn_val = window_cn(win_df, chrom_win, start, end)

        # Use window-level CN for verdict when available, segment as fallback
        est_cn = win_cn_val if win_cn_val is not None else est_cn_seg
        cn_source = "window" if win_cn_val is not None else "segment"

        if vmode == "skip":
            rows.append({
                "Region": gene, "Category": category,
                "Chr": chrom, "GT_Start": start, "GT_End": end,
                "GT_Len_kb": round(gt_len / 1e3, 1),
                "Expected_CN": exp_cn, "Tolerance": "—",
                "Estimated_CN": None, "Seg_CN": None, "Peak_CN": None,
                "Abs_Error": None, "Pct_Error": "—",
                "N_Segs": 0, "Cov_bp": 0, "Cov_Pct": 0.0,
                "Dom_State": "—", "Verdict": "SKIPPED",
                "CN_Source": "—", "Notes": notes,
            })
            continue

        if est_cn is None:
            abs_err = None
            rel_err = None
            pct_err_str = "N/A"
            verdict = "NO DATA"
        else:
            abs_err = round(abs(est_cn - exp_cn), 4)
            rel_err = round((est_cn - exp_cn) / exp_cn * 100, 1) if exp_cn > 0 else None
            pct_err_str = f"{rel_err:+.1f}%" if rel_err is not None else "N/A"

            if vmode == "observatory":
                verdict = "OBSERVATORY"
            elif exp_cn == 0:
                verdict = "OK" if est_cn < 0.5 else "FP"
            elif abs_err <= tol * exp_cn:
                verdict = f"PASS (±{tol*100:.0f}%)"
            elif abs_err <= 0.50 * exp_cn:
                verdict = "MARGINAL (<50%)"
            else:
                dir_flag = "OVER" if est_cn > exp_cn else "UNDER"
                verdict = f"FAIL ({dir_flag})"

        n_segs  = len(hits) if not hits.empty else 0
        cov_pct = round(cov_bp / gt_len * 100, 1) if gt_len > 0 else 0
        tol_str = f"±{tol*100:.0f}%" if tol is not None else "obs"

        rows.append({
            "Region":       gene,
            "Category":     category,
            "Chr":          chrom,
            "GT_Start":     start,
            "GT_End":       end,
            "GT_Len_kb":    round(gt_len / 1e3, 1),
            "Expected_CN":  exp_cn,
            "Tolerance":    tol_str,
            "Estimated_CN": est_cn,
            "Seg_CN":       est_cn_seg,
            "Peak_CN":      peak_cn,
            "Abs_Error":    abs_err,
            "Pct_Error":    pct_err_str,
            "N_Segs":       n_segs,
            "Cov_bp":       cov_bp,
            "Cov_Pct":      cov_pct,
            "Dom_State":    dom_state,
            "Verdict":      verdict,
            "CN_Source":    cn_source,
            "Notes":        notes,
        })
    return pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# MARKDOWN REPORT GENERATOR
# ---------------------------------------------------------------------------
def build_markdown(df: pd.DataFrame, segments_path: str, kmer_size: int = 72,
                   sample: str = 'chm13') -> str:
    run_name = os.path.basename(os.path.dirname(segments_path))
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    cfg = SAMPLE_CONFIG[sample]
    validated_k = cfg['validated_k']

    valid = df[df["Verdict"] != "SKIPPED"]
    scored = valid[valid["Estimated_CN"].notna() & (valid["Verdict"] != "OBSERVATORY")]
    obs    = valid[valid["Verdict"] == "OBSERVATORY"]

    n_pass     = scored["Verdict"].str.startswith("PASS").sum()
    n_marginal = (scored["Verdict"] == "MARGINAL (<50%)").sum()
    n_fail     = scored["Verdict"].str.startswith("FAIL").sum()
    n_nodata   = valid["Estimated_CN"].isna().sum()
    n_obs      = len(obs)

    sdnv = scored[scored["Category"].isin(["Segmental Dup", "VNTR Array"])]
    rdna = valid[valid["Category"] == "rDNA Array"]

    lines = []
    def ln(s=""): lines.append(s)

    ln("# CopySeg Ground Truth Evaluation Report")
    ln()
    ln(f"**Run:** `{run_name}`  ")
    ln(f"**Sample:** {cfg['label']}  ")
    ln(f"**Segments file:** `{segments_path}`  ")
    ln(f"**Generated:** {now}  ")
    ln(f"**k-mer length:** k={kmer_size} | **GT validated at:** k={validated_k} | "
       f"**Reference:** T2T-CHM13v2.0")
    ln()
    ln(f"> **Metodoloji:** expected_cn = k{kmer_size} aggregate CN = n_copies × identity^{kmer_size}.")
    ln("> Fiziksel kopya sayısı değil — tüm yüksek-kimlikli paraloglardaki k-merler toplanır.")
    if kmer_size != validated_k:
        ln(f"> **Not:** k={kmer_size} kullanıldı. k={validated_k} referans değerlerinden "
           f"farklı expected CN hesaplanmıştır (n_copies × identity^{kmer_size}).")
    ln()

    ln("## Özet")
    ln()
    ln("| Metrik | Değer |")
    ln("|--------|-------|")
    ln(f"| PASS | **{n_pass}** |")
    ln(f"| MARGINAL (<50%) | {n_marginal} |")
    ln(f"| FAIL | **{n_fail}** |")
    ln(f"| OBSERVATORY | {n_obs} |")
    ln(f"| Veri yok | {n_nodata} |")
    if not sdnv.empty:
        ln(f"| Median abs. hata (Seg.Dup/VNTR) | {sdnv['Abs_Error'].median():.2f} CN |")
    if not rdna.empty and rdna["Estimated_CN"].notna().any():
        ln(f"| Median est. CN (rDNA, obs.) | {rdna['Estimated_CN'].median():.1f} |")
    ln()

    ln("## Bölge Bazlı Sonuçlar")
    ln()
    has_windows = "Seg_CN" in df.columns and df["CN_Source"].eq("window").any()
    if has_windows:
        ln("| Bölge | Kategori | Exp.CN | Tol | Est.CN | Seg.CN | Peak.CN | %Err | N Seg | Cov% | Src | Sonuç |")
        ln("|-------|----------|--------|-----|--------|--------|---------|------|-------|------|-----|-------|")
    else:
        ln("| Bölge | Kategori | Exp.CN | Tol | Est.CN | Peak.CN | %Err | N Seg | Cov% | Sonuç |")
        ln("|-------|----------|--------|-----|--------|---------|------|-------|------|-------|")
    for _, r in df.iterrows():
        if r["Verdict"] == "SKIPPED":
            if has_windows:
                ln(f"| {r['Region']} | {r['Category']} | {r['Expected_CN']} | — | — | — | — | — | — | — | — | SKIPPED |")
            else:
                ln(f"| {r['Region']} | {r['Category']} | {r['Expected_CN']} | — | — | — | — | — | — | SKIPPED |")
            continue
        est  = f"{r['Estimated_CN']:.2f}" if pd.notna(r["Estimated_CN"]) else "—"
        peak = f"{r['Peak_CN']:.2f}"      if pd.notna(r["Peak_CN"])      else "—"
        if has_windows:
            seg_cn = f"{r['Seg_CN']:.2f}" if pd.notna(r.get("Seg_CN")) else "—"
            src = r.get("CN_Source", "seg")[:3]
            ln(f"| {r['Region']} | {r['Category']} | {r['Expected_CN']} | {r['Tolerance']} | "
               f"{est} | {seg_cn} | {peak} | {r['Pct_Error']} | {r['N_Segs']} | {r['Cov_Pct']}% | {src} | {r['Verdict']} |")
        else:
            ln(f"| {r['Region']} | {r['Category']} | {r['Expected_CN']} | {r['Tolerance']} | "
               f"{est} | {peak} | {r['Pct_Error']} | {r['N_Segs']} | {r['Cov_Pct']}% | {r['Verdict']} |")
    ln()

    ln("## Kategori Analizi")
    ln()
    ln("### 1. Segmental Duplikasyonlar & VNTR")
    ln()
    ln(f"_k{kmer_size} aggregate CN = n_copies × identity^{kmer_size} (fiziksel kopya sayısı değil)_")
    ln()
    sdnv2 = df[df["Category"].isin(["Segmental Dup", "VNTR Array"])
               & (df["Verdict"] != "SKIPPED")]
    for _, r in sdnv2.iterrows():
        est = r["Estimated_CN"]
        exp = r["Expected_CN"]
        if pd.isna(est):
            ln(f"**{r['Region']}** ({r['Chr']}, {r['GT_Len_kb']} kb) — VERİ YOK")
            ln()
            continue
        direction = "overestimation" if est > exp else "underestimation"
        ln(f"**{r['Region']}** ({r['Chr']}, {r['GT_Len_kb']} kb)")
        seg_note = ""
        if has_windows and pd.notna(r.get("Seg_CN")):
            seg_note = f" | Seg.CN: {r['Seg_CN']:.2f}"
        ln(f"- Expected k{kmer_size} CN: **{exp}** | Estimated: **{est:.2f}**{seg_note} | Peak: **{r['Peak_CN']:.2f}** | "
           f"Abs.Error: {r['Abs_Error']:.2f} | {r['Pct_Error']} | Tolerance: {r['Tolerance']}")
        ln(f"- Dom.State: `{r['Dom_State']}` | N segs: {r['N_Segs']} | Kapsama: %{r['Cov_Pct']}")
        ln(f"- Verdict: **{r['Verdict']}** ({direction})")
        ln(f"- {r['Notes']}")
        ln()

    rdna_agg = cfg.get('rdna_agg_cn')
    if rdna_agg and not rdna.empty:
        ln("### 2. rDNA Dizileri — OBSERVATORY")
        ln()
        ln(f"_{rdna_agg} total rDNA units across 5 NOR loci._")
        ln(f"_Inter-NOR identity 99.4-99.7% → k-mers aggregate all {rdna_agg} units genome-wide._")
        ln()
        for _, r in rdna.iterrows():
            est = r["Estimated_CN"]
            if pd.isna(est):
                ln(f"**{r['Region']}** — VERİ YOK")
                continue
            pct_of_ref = round(est / rdna_agg * 100, 1)
            ln(f"**{r['Region']}** ({r['Chr']}, {r['GT_Len_kb']} kb): "
               f"est={est:.1f} | peak={r['Peak_CN']:.1f} | {r['Pct_Error']} vs ref={rdna_agg} | "
               f"ref capture: {pct_of_ref}% | Dom.State: `{r['Dom_State']}`")
        ln()
    elif not rdna.empty:
        ln("### 2. rDNA Dizileri — SKIPPED")
        ln()
        ln(f"_rDNA evaluation skipped for {cfg['label']} (assembly incomplete)._")
        ln()

    ln(f"## k{kmer_size} Aggregate CN Referans Tablosu ({cfg['label']})")
    ln()
    ln(f"| Gen Ailesi | Fiziksel CN | Identity | k{validated_k} Agg CN | k{kmer_size} Agg CN | Kaynak |")
    ln("|-----------|-------------|----------|------------|------------|--------|")
    for entry in cfg['gt']:
        gene = entry[3]
        exp_validated = entry[4]
        n_c = entry[9]
        ident = entry[10]
        vmode = entry[6]
        if vmode == "skip":
            continue
        exp_k = compute_expected_cn(exp_validated, n_c, ident, kmer_size, validated_k)
        id_str = f"{ident*100:.1f}%" if ident is not None else "—"
        nc_str = str(n_c) if n_c is not None else "—"
        k_ref = f"k={validated_k} validated"
        ln(f"| {gene} | {nc_str} | {id_str} | {exp_validated} | **{exp_k}** | {k_ref} |")
    ln()

    ln("## Özet: Kategori Bazlı Performans")
    ln()
    ln("| Kategori | N | Median Exp.CN | Median Est.CN | Median Abs.Err | Yön |")
    ln("|----------|---|--------------|--------------|----------------|-----|")
    for cat, grp in df[df["Verdict"] != "SKIPPED"].groupby("Category"):
        vg = grp[grp["Estimated_CN"].notna()]
        if vg.empty:
            continue
        med_exp = vg["Expected_CN"].median()
        med_est = vg["Estimated_CN"].median()
        med_err = vg["Abs_Error"].median() if vg["Abs_Error"].notna().any() else float("nan")
        direction = "↑ over" if med_est > med_exp else "↓ under"
        err_str = f"{med_err:.2f}" if not np.isnan(med_err) else "obs"
        ln(f"| {cat} | {len(vg)} | {med_exp:.1f} | {med_est:.1f} | {err_str} | {direction} |")
    ln()

    ln("---")
    ln(f"*Rapor otomatik üretildi: evaluate_ground_truth.py — {now}*")
    ln(f"*GT revision: 2026-03-29b ASSEMBLY-VALIDATED (Jellyfish k=72 CHM13v2.0): "
       f"AMY1 coords fixed, AMY2A 1→4, SMN 3→2 (normal), LPA 10→9, TBC1D3 2→9, 5S_rDNA 30→117*")

    return "\n".join(lines)

# ---------------------------------------------------------------------------
# ENTRY POINT
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="CopySeg Ground Truth Evaluator")
    parser.add_argument("--segments", required=True, help="CopySeg BED file")
    parser.add_argument("--output", default=None, help="Output Markdown path")
    parser.add_argument("--kmer-size", type=int, default=72,
                        help="K-mer size for expected CN computation (default: 72)")
    parser.add_argument("--sample", choices=list(SAMPLE_CONFIG.keys()), default="chm13",
                        help="Sample name — selects the appropriate GT table "
                             f"(choices: {', '.join(SAMPLE_CONFIG.keys())}; default: chm13)")
    parser.add_argument("--windows", default=None,
                        help="Window-level CN BED (cn_w500.bed). When provided, uses "
                             "window median within GT region for CN estimation.")
    args = parser.parse_args()

    if not os.path.exists(args.segments):
        print(f"ERROR: segments file not found: {args.segments}", file=sys.stderr)
        sys.exit(1)
    if args.windows and not os.path.exists(args.windows):
        print(f"ERROR: windows file not found: {args.windows}", file=sys.stderr)
        sys.exit(1)

    cfg = SAMPLE_CONFIG[args.sample]
    cn_mode = "window-level" if args.windows else "segment-level"
    print(f"[EVAL] Loading segments: {args.segments} "
          f"(k={args.kmer_size}, sample={args.sample}, "
          f"GT validated at k={cfg['validated_k']}, CN mode: {cn_mode})", file=sys.stderr)
    if args.windows:
        print(f"[EVAL] Loading windows: {args.windows}", file=sys.stderr)
    df = evaluate(args.segments, kmer_size=args.kmer_size, sample=args.sample,
                  windows_path=args.windows)

    print("\n" + "="*80)
    print(f"CopySeg Ground Truth Evaluation — SUMMARY (CN: {cn_mode})")
    print("="*80)
    summary_cols = ["Region", "Category", "Expected_CN", "Tolerance",
                    "Estimated_CN", "Peak_CN", "Pct_Error", "Verdict"]
    if args.windows:
        summary_cols.insert(5, "Seg_CN")
        summary_cols.append("CN_Source")
    print(df[summary_cols].to_string(index=False))
    print("="*80)

    md = build_markdown(df, args.segments, kmer_size=args.kmer_size,
                        sample=args.sample)
    if args.output:
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        with open(args.output, "w") as f:
            f.write(md)
        print(f"\n[EVAL] Markdown report → {args.output}", file=sys.stderr)
    else:
        print("\n" + md)

    sys.exit(0)

if __name__ == "__main__":
    main()

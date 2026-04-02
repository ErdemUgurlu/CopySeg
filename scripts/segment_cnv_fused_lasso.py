#!/usr/bin/env python3
"""
KonuSeg Weighted PELT Segmentation  (Fused-Lasso Formulation)
=============================================================

Drop-in alternative to segment_cnv_hmm_log_7state.py and segment_cnv_cbs.py.
Identical I/O format, identical post-processing chain, but the segmentation
core uses a *quality-weighted* PELT cost that is mathematically equivalent to
the weighted Fused Lasso:

    minimize  Σ wᵢ (yᵢ − θᵢ)²  +  λ Σ |θᵢ − θᵢ₊₁|
              i                     i

where
    yᵢ   = log2(CN)   for window i
    θᵢ   = piecewise-constant estimate (= output CN profile)
    wᵢ   = quality weight  (num_kmers / (num_kmers + num_filtered + 1))
    λ    = BIC per-chromosome (2·log n) or user-supplied --penalty

The objective is solved exactly via PELT (Pruned Exact Linear Time) using a
custom O(1) WeightedL2Cost that precomputes three prefix-sum arrays (Σw,
Σwy, Σwy²) so every segment cost query costs O(1).

KEY DIFFERENCES FROM HMM (segment_cnv_hmm_log_7state.py)
─────────────────────────────────────────────────────────
Issue B1  resolved: no emission model → no wrong-variance problem
Issue B2  resolved: no fixed states → CN=9-31 gap eliminated
Issue B3  resolved: parameters learned from data via BIC selection
Issue B4  resolved: no transition matrix → distance issue irrelevant
Issue B5  resolved: penalty controls granularity (not self_transition)
Issue EK2 resolved: no wasted HomDel/HetDel states on haploid CHM13
Quality weight: wᵢ down-weights noisy windows *without* biasing CN toward 1.0
    → eliminates the 20-37 % underestimation seen with QW>0 in HMM

KEY DIFFERENCES FROM CBS (segment_cnv_cbs.py)
──────────────────────────────────────────────
Quality weights enter the *segmentation cost* (CBS always uses QW=0)
cn_median = quality-weighted mean (CBS uses true median)
C3 fix: merge allows adjacent dup states with small CN difference (not only
    identical-state pairs)
All other CBS fixes are carried forward (C5, C7, EK3).

POST-PROCESSING CHAIN (ORDER MATTERS — mirrors cn-accuracy HMM mode)
─────────────────────────────────────────────────────────────────────
a. filter_small_segments         (min segment length per state)
b. filter_low_quality_segments   (quality < threshold → Neutral; CV filter)
c. merge_nearby_dup_segments     (dup segs + short Neutral gap ≤ max_merge_gap)
   [C3 FIX: allow cross-state dup merge when |ΔCN| < merge_cn_tolerance]
d. split_high_cv_segments        (CV > threshold → split at variance minimum)
e. add_repeat_annotation         (masked_fraction, repeat_class)
   [MUST come before GC calibration — bypass conditions use repeat_class]
f. apply_gc_cn_calibration       (post-segmentation GC bias correction)
   [bypass: Satellite/Low_complexity, GC>0.60, masked_frac>0.85 & GC>0.55]
g. reclassify_by_cn_threshold    (state–CN consistency check)

CONFIDENCE METRIC (replaces HMM entropy)
─────────────────────────────────────────
boundary_conf = Welch t-test |t| at segment boundary
    high value → clear CN change → confident breakpoint
segment_iqr   = IQR of raw CN values within segment
    low value  → homogeneous segment → confident flat CN

Usage:
    python3 scripts/segment_cnv_fused_lasso.py \\
        --input  output/chm13_my_clean_v10/chm13_ont_cn_w500.bed \\
        --output output/chm13_my_clean_v11/segs_flasso_w500.bed \\
        --extended \\
        --gc-content-bed output/c5_c5_v1_chrx_gc/gc_content_w500.bed \\
        --repeat-bed output/chm13_my_clean_v10/repeat_annotated_w500.bed \\
        --min-kmers 30

Install deps (if not present):
    pip install ruptures>=1.1.7 scipy>=1.7.0

Output columns (extended, 15 base + optional extras):
    chrom, start, end, state, cn_median, cn_mean, n_windows,
    avg_quality, min_quality, cn_std, avg_repeats,
    avg_entropy (=0), max_entropy (=0),
    masked_fraction, repeat_class,
    [gc_bias_factor],  [segment_iqr],  [boundary_conf]

Author: KonuSeg Team
Date:   March 2026
"""

import argparse
import multiprocessing as _mp
import os
import sys
import time as _time
import warnings
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

try:
    import ruptures as rpt
    from ruptures.base import BaseCost
    RUPTURES_AVAILABLE = True
except ImportError:
    RUPTURES_AVAILABLE = False
    BaseCost = object  # placeholder so class definition doesn't crash

try:
    from scipy.stats import ttest_ind
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

warnings.filterwarnings('ignore')

# =============================================================================
# CONSTANTS
# =============================================================================

EPSILON = 1e-3       # Prevents log2(0); same as HMM/CBS
WINDOW_BP = 500      # Expected window size in bp (CHM13 ONT default)

# CN → state mapping (post-hoc classification; same as HMM/CBS CN_STATE_BOUNDARIES).
# Threshold at 1.25 keeps backward compatibility with evaluate_ground_truth.py.
# C8 note: HMM emission midpoint (Neutral/LowDup) is CN=1.52; threshold here is
# intentionally kept at 1.25 (conservative) and overridable via --lowdup-threshold.
CN_STATE_BOUNDARIES = [
    (0.30,         'HomDel'),
    (0.70,         'HetDel'),
    (1.25,         'Neutral'),
    (3.00,         'LowDup'),
    (6.00,         'HighDup'),
    (12.00,        'Amp'),
    (22.00,        'MedAmp'),
    (50.00,        'HighAmp'),
    (float('inf'), 'ExtremeAmp'),
]

STATE_ORDER = ['HomDel', 'HetDel', 'Neutral', 'LowDup', 'HighDup', 'Amp',
               'MedAmp', 'HighAmp', 'ExtremeAmp']

DUP_STATES = frozenset({'LowDup', 'HighDup', 'Amp', 'MedAmp', 'HighAmp', 'ExtremeAmp'})


# =============================================================================
# WEIGHTED L2 COST FOR RUPTURES PELT
# =============================================================================

class WeightedL2Cost(BaseCost):
    """O(1) weighted L2 cost function for ruptures PELT.

    Minimises Σ wᵢ (yᵢ − ȳ_w)² for each candidate segment, where ȳ_w is the
    quality-weighted mean.  Three prefix-sum arrays (cumW, cumWY, cumWY2) are
    precomputed in fit() so every error(start, end) call is O(1).

    cost(l, r) = Σwᵢyᵢ² − (Σwᵢyᵢ)² / Σwᵢ
               = cumWY2[r]−cumWY2[l] − (cumWY[r]−cumWY[l])² / (cumW[r]−cumW[l])

    This is equivalent to asking: "how far do the weighted observations deviate
    from their weighted mean within the segment?" — the optimal breakpoint is
    where this cost is minimised across all candidate splits.

    Signal shape expected: (n_samples, 2)
        column 0 → log2_cn  (the signal to segment)
        column 1 → quality  (per-window weight, clamped to [0.01, 1.0])
    """

    model = "weighted_l2"
    min_size = 1

    def fit(self, signal):
        """Precompute prefix sums from (n, 2) signal array.

        Stores self.signal for ruptures BaseCost compatibility
        (BaseCost.n_samples property accesses self.signal.shape[0]).
        """
        # ruptures BaseCost requires self.signal to be set
        self.signal = signal

        if signal.ndim == 1:
            # Fallback: unweighted (quality column missing)
            y = signal.astype(np.float64)
            w = np.ones(len(y), dtype=np.float64)
        else:
            y = signal[:, 0].astype(np.float64)
            w = np.clip(signal[:, 1].astype(np.float64), 0.01, 1.0)

        self.n = len(y)
        self.cumW   = np.concatenate(([0.0], np.cumsum(w)))
        self.cumWY  = np.concatenate(([0.0], np.cumsum(w * y)))
        self.cumWY2 = np.concatenate(([0.0], np.cumsum(w * y ** 2)))
        return self

    def error(self, start, end):
        """Weighted L2 cost of segment [start, end). O(1)."""
        W   = self.cumW[end]   - self.cumW[start]
        WY  = self.cumWY[end]  - self.cumWY[start]
        WY2 = self.cumWY2[end] - self.cumWY2[start]
        if W < 1e-12:
            return 0.0
        return float(WY2 - WY * WY / W)

    def sum_of_costs(self, bkps):
        """Total cost across a list of breakpoints (required by ruptures interface)."""
        cost = 0.0
        prev = 0
        for bp in bkps:
            cost += self.error(prev, bp)
            prev = bp
        return cost


# =============================================================================
# STATE ASSIGNMENT
# =============================================================================

def assign_state_from_cn(cn, lowdup_threshold=1.25):
    """Map calibrated CN to the nearest state.

    Uses CN_STATE_BOUNDARIES with configurable Neutral/LowDup boundary.
    Passing lowdup_threshold=1.52 aligns with the HMM log2 emission midpoint
    (C8 fix — optional, configurable at runtime via --lowdup-threshold).
    """
    boundaries = CN_STATE_BOUNDARIES[:]
    # Override Neutral boundary if requested
    if lowdup_threshold != 1.25:
        boundaries = [(lowdup_threshold if st == 'Neutral' else hi, st)
                      for hi, st in boundaries]
    for upper, state in boundaries:
        if cn < upper:
            return state
    return 'ExtremeAmp'


# =============================================================================
# DATA LOADING
# =============================================================================

def load_data(bed_file, min_kmers=0):
    """Load 8-column windows BED and compute log2_cn + quality.

    Quality weight = num_kmers / (num_kmers + num_filtered + 1).
    Quality weight is stored in the DataFrame for use in segmentation cost
    and CN estimation — it does NOT modify the raw CN values.

    min_kmers: windows with num_kmers < min_kmers are set to CN=1.0 (Neutral
    prior) to suppress false HetDel calls in low-coverage / complex regions.
    """
    print(f"[IO] Loading {bed_file}...")

    try:
        df = pd.read_csv(bed_file, sep='\t', header=None, comment='#')

        if len(df.columns) >= 8:
            df.columns = ['chrom', 'start', 'end', 'cn', 'mean_count', 'log_ratio',
                          'num_kmers', 'num_filtered'][:len(df.columns)]
        elif len(df.columns) >= 4:
            df.columns = (['chrom', 'start', 'end', 'cn'] +
                          [f'col{i}' for i in range(4, len(df.columns))])
        else:
            raise ValueError(f"Unexpected column count: {len(df.columns)}")

        print(f"[IO] Loaded {len(df):,} windows")

        # Low-coverage filter
        if min_kmers > 0 and 'num_kmers' in df.columns:
            low_cov = df['num_kmers'] < min_kmers
            n_low = int(low_cov.sum())
            if n_low > 0:
                df.loc[low_cov, 'cn'] = 1.0
                print(f"[IO] Low-coverage filter (num_kmers < {min_kmers}): "
                      f"{n_low:,} windows ({100 * n_low / len(df):.1f}%) set to CN=1.0")

        # Quality weight — stored but does NOT modify cn
        if 'num_kmers' in df.columns and 'num_filtered' in df.columns:
            df['quality'] = (df['num_kmers'] /
                             (df['num_kmers'] + df['num_filtered'] + 1))
        else:
            df['quality'] = 1.0

        # Log2 CN for PELT segmentation — always raw (quality weighting is in cost)
        df['log2_cn'] = np.log2(df['cn'].values + EPSILON)

        # Remove inf/nan
        before = len(df)
        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.dropna(subset=['log2_cn'])
        if before != len(df):
            print(f"[IO] Removed {before - len(df)} invalid rows")

        print(f"[IO] CN range:    {df['cn'].min():.4f} – {df['cn'].max():.4f}")
        print(f"[IO] Log2 range:  {df['log2_cn'].min():.2f} – {df['log2_cn'].max():.2f}")
        print(f"[IO] Median CN:   {df['cn'].median():.4f}")
        print(f"[IO] Quality:     median={df['quality'].median():.3f}, "
              f"low(<0.5)={100*(df['quality']<0.5).mean():.1f}%")
        return df

    except Exception as e:
        print(f"[ERROR] Failed to load data: {e}")
        sys.exit(1)


def load_gc_content(gc_bed_file):
    """Load GC content per window (chrom, start, end, gc_content)."""
    print(f"[GC] Loading GC content from {gc_bed_file}...")
    gc_df = pd.read_csv(gc_bed_file, sep='\t', comment='#', header=None)
    if len(gc_df.columns) < 4:
        raise ValueError(f"GC BED needs >= 4 columns, got {len(gc_df.columns)}")
    gc_df.columns = ['chrom', 'start', 'end', 'gc_content'][:len(gc_df.columns)]
    gc_df = gc_df[['chrom', 'start', 'gc_content']]
    print(f"[GC] Loaded {len(gc_df):,} windows  (gc mean={gc_df['gc_content'].mean():.4f})")
    return gc_df


def load_repeat_annotation(repeat_bed_file):
    """Load per-window repeat annotation (10-col BED from compute_repeat_annotation.py)."""
    print(f"[RM] Loading repeat annotation from {repeat_bed_file}...")
    rm_df = pd.read_csv(
        repeat_bed_file, sep='\t', comment='#', header=None,
        names=['chrom', 'start', 'end', 'cn', 'mean_count', 'log_ratio',
               'num_kmers', 'num_filtered', 'masked_fraction', 'repeat_class'],
        usecols=[0, 1, 8, 9],
    )
    rm_df.columns = ['chrom', 'start', 'masked_fraction', 'repeat_class']
    pct = (rm_df['masked_fraction'] > 0.5).mean() * 100
    print(f"[RM] Loaded {len(rm_df):,} windows  (>50% masked: {pct:.1f}%)")
    return rm_df


# =============================================================================
# SEGMENT CONSTRUCTION
# =============================================================================

def make_segment_from_windows(windows_df, lowdup_threshold=1.25):
    """Build a segment dict from a slice of the windows DataFrame.

    cn_median  = quality-weighted mean of constituent window CN values.
                 This is the PRIMARY CN estimate for this segment.
                 (C2 note: name 'cn_median' retained for evaluate_ground_truth.py
                  backward compatibility; the value IS a weighted mean, not median.)
    cn_mean    = unweighted mean (secondary, for reference).
    cn_std     = quality-weighted std of CN values.
                 Uses the same quality weights as cn_median so that CV = cn_std/cn_median
                 is self-consistent. Without this, low-quality outlier windows inflate
                 cn_std while barely affecting cn_median → spurious high CV → false
                 CV-split and CV-filter triggers.
    avg_entropy / max_entropy = 0.0 (no HMM → kept for downstream compatibility).
    """
    cn_arr = windows_df['cn'].values
    n = len(cn_arr)

    # Quality weights
    if 'quality' in windows_df.columns:
        w = windows_df['quality'].values.clip(0.01, 1.0)
    else:
        w = np.ones(n)
    w_sum = float(w.sum())
    if w_sum < 1e-12:
        w = np.ones(n)
        w_sum = float(n)

    # Quality-weighted CN mean → primary estimate (no CN distortion)
    cn_wm = float(np.dot(w, cn_arr) / w_sum)

    cn_mean = float(np.mean(cn_arr))
    # Quality-weighted std: consistent with cn_median (both use same weights).
    # σ_w = sqrt(Σwᵢ(cnᵢ - cn_wm)² / Σwᵢ)
    # This prevents low-quality outliers from inflating cn_std while barely
    # moving cn_median, which would produce spurious high CV = cn_std/cn_median.
    cn_std = (float(np.sqrt(np.dot(w, (cn_arr - cn_wm) ** 2) / w_sum))
              if n > 1 else 0.0)

    q_arr = windows_df['quality'].values if 'quality' in windows_df.columns else np.ones(n)
    avg_quality = float(np.mean(q_arr))
    min_quality = float(np.min(q_arr))

    avg_repeats = (float(windows_df['num_filtered'].mean())
                   if 'num_filtered' in windows_df.columns else 0.0)

    return {
        'chrom':       str(windows_df['chrom'].iloc[0]),
        'start':       int(windows_df['start'].iloc[0]),
        'end':         int(windows_df['end'].iloc[-1]),
        'state':       assign_state_from_cn(cn_wm, lowdup_threshold),
        'cn_median':   cn_wm,      # quality-weighted mean (primary CN estimate)
        'cn_mean':     cn_mean,    # unweighted mean (secondary)
        'n_windows':   n,
        'avg_quality': avg_quality,
        'min_quality': min_quality,
        'cn_std':      cn_std,
        'avg_repeats': avg_repeats,
        'avg_entropy': 0.0,        # no HMM — kept for downstream compat
        'max_entropy': 0.0,
    }


# =============================================================================
# CORE SEGMENTATION: WEIGHTED PELT
# =============================================================================

def _estimate_noise_var(y):
    """Estimate within-segment noise variance from log2_cn gradient.

    Uses adjacent differences with 4×MAD clipping to exclude CN transition
    artifacts (large jumps at true segment boundaries skew the variance upward).
    Returns σ²_noise ≥ 0.001 (floor prevents degenerate near-zero penalty).

    For independent Gaussian noise:  Var(diff(y)) = 2 · σ²_noise
    so σ²_noise = 0.5 · Var(diff(y)[clipped])

    Why this matters for BIC penalty:
        pen = 2·log(n)·σ²_noise (variance-normalized BIC)

    Without normalization (plain 2·log(n)):
        For CHM13 σ_noise ≈ 0.15 log2, pen ≈ 25 >> E[max noise gain] ≈ 0.28
        → under-segments: misses small LowDup/HighDup transitions.
    With normalization:
        pen ≈ 2·12.6·0.023 ≈ 0.58, still >> E[max noise gain] ≈ 0.28
        → correctly balanced: detects all real CN changes, ignores noise.
    """
    if len(y) < 4:
        return 0.01
    diffs = np.diff(y)
    mad = float(np.median(np.abs(diffs)))
    # 4×MAD clip: retains >99.99% of Gaussian noise while excluding real CN jumps
    clip_thr = max(4.0 * mad, 0.01)
    diffs_noise = diffs[np.abs(diffs) <= clip_thr]
    if len(diffs_noise) < 2:
        return 0.01
    return max(0.5 * float(np.nanvar(diffs_noise)), 0.001)

def run_wpelt_chrom(df_chrom, penalty, min_size_windows, use_weights=True,
                    lowdup_threshold=1.25):
    """Run weighted PELT on a single chromosome.

    Uses WeightedL2Cost (custom ruptures cost) when use_weights=True,
    falling back to standard L2 when use_weights=False or quality column absent.

    Weighted cost objective for segment [l, r):
        Σᵢ∈[l,r) wᵢ (yᵢ − ȳ_w)²
    where ȳ_w = Σwᵢyᵢ / Σwᵢ  (quality-weighted mean in log2 space).

    Breakpoint minimises total weighted cost across all segments.
    BIC penalty (2 log n) controls the segment count without requiring a
    fixed state count or emission model.
    """
    n = len(df_chrom)
    if n == 0:
        return []

    # Short chromosome: return single segment
    if n < max(min_size_windows * 2, 4):
        return [make_segment_from_windows(df_chrom, lowdup_threshold)]

    y = df_chrom['log2_cn'].values

    has_quality = use_weights and 'quality' in df_chrom.columns
    if has_quality:
        w = df_chrom['quality'].values.clip(0.01, 1.0)
        signal = np.column_stack([y, w])   # shape (n, 2)
        cost_obj = WeightedL2Cost()
        algo = rpt.Pelt(custom_cost=cost_obj, min_size=min_size_windows, jump=5)
    else:
        signal = y
        algo = rpt.Pelt(model='l2', min_size=min_size_windows, jump=5)

    try:
        breakpoints = algo.fit_predict(signal, pen=penalty)
    except Exception as e:
        print(f"[WPELT] WARNING: PELT failed ({e}), returning single segment")
        return [make_segment_from_windows(df_chrom, lowdup_threshold)]

    segments = []
    prev = 0
    for bp in breakpoints:
        if bp <= prev:
            continue
        seg = make_segment_from_windows(df_chrom.iloc[prev:bp], lowdup_threshold)
        segments.append(seg)
        prev = bp

    return segments


def run_segmentation(df, penalty=None, min_size=3000, use_weights=True,
                     lowdup_threshold=1.25, penalty_factor=4.0):
    """Run per-chromosome weighted PELT segmentation.

    penalty: None → BIC = penalty_factor*log(n) per chromosome (recommended).
             float → fixed penalty for all chromosomes.
    min_size: minimum segment size in base pairs (converted to windows internally).
    use_weights: if True, use quality-weighted L2 cost (WeightedL2Cost).
                 if False, use standard L2 (equivalent to CBS/ruptures default).
    penalty_factor: BIC multiplier (default 2.0 = standard BIC). Higher → fewer segments.
    """
    window_sizes = (df['end'] - df['start']).values
    window_bp = int(np.median(window_sizes)) if len(window_sizes) > 0 else WINDOW_BP
    if window_bp != WINDOW_BP:
        print(f"[WPELT] Non-standard window size detected: {window_bp} bp (expected {WINDOW_BP})")

    min_size_windows = max(1, min_size // window_bp)
    weight_mode = "weighted L2" if use_weights else "unweighted L2"
    print(f"[WPELT] Algorithm:          {weight_mode} PELT")
    print(f"[WPELT] min_segment_length: {min_size} bp → {min_size_windows} windows")
    print(f"[WPELT] Penalty:            "
          f"{'variance-normalized BIC: ' + str(penalty_factor) + '·log(n)·σ²_noise per chrom' if penalty is None else penalty}")

    chroms = list(df['chrom'].unique())
    total_t0 = _time.time()

    n_workers = min(len(chroms), int(os.environ.get('KONUSEG_WORKERS', 0))
                    or max(1, _mp.cpu_count() - 1))

    # Build per-chromosome work items
    chrom_work = []
    for chrom in chroms:
        df_chrom = df[df['chrom'] == chrom].sort_values('start').reset_index(drop=True)
        n = len(df_chrom)
        if penalty is not None:
            pen = penalty
            nv = None
        else:
            nv = _estimate_noise_var(df_chrom['log2_cn'].values)
            pen = penalty_factor * np.log(max(n, 2)) * nv
        chrom_work.append((chrom, df_chrom, pen, nv, min_size_windows,
                           use_weights, lowdup_threshold))

    if n_workers > 1:
        print(f"[WPELT] Parallel mode: {n_workers} workers for {len(chroms)} chromosomes")
        with _mp.Pool(n_workers) as pool:
            results = pool.map(_segment_one_chrom_worker, chrom_work)
    else:
        print(f"[WPELT] Sequential mode: {len(chroms)} chromosomes")
        results = [_segment_one_chrom_worker(w) for w in chrom_work]

    # Collect results and print per-chromosome stats
    all_segments = []
    for chrom, segs, elapsed, pen, nv in results:
        all_segments.extend(segs)
        if nv is not None:
            pen_str = f"pen={pen:.4f}, noise_var={nv:.6f}"
        else:
            pen_str = f"pen={pen:.4f} (fixed)"
        n_win = sum(1 for w in chrom_work if w[0] == chrom)
        # Find original window count from work item
        for cw in chrom_work:
            if cw[0] == chrom:
                n_win = len(cw[1])
                break
        print(f"[WPELT] {chrom}: {n_win:,} windows → {len(segs):,} segments "
              f"({elapsed:.1f}s, {pen_str})")

    total_elapsed = _time.time() - total_t0
    print(f"[WPELT] Total: {len(df):,} windows → {len(all_segments):,} segments "
          f"({total_elapsed:.1f}s, {n_workers} workers)")
    return all_segments


def _segment_one_chrom_worker(args):
    """Worker function for parallel chromosome segmentation."""
    chrom, df_chrom, pen, noise_var, min_size_windows, use_weights, lowdup_threshold = args
    try:
        t0 = _time.time()
        segs = run_wpelt_chrom(df_chrom, pen, min_size_windows,
                               use_weights=use_weights,
                               lowdup_threshold=lowdup_threshold)
        elapsed = _time.time() - t0
        return (chrom, segs, elapsed, pen, noise_var)
    except Exception as e:
        raise RuntimeError(f"PELT failed on {chrom} ({len(df_chrom)} windows, "
                           f"pen={pen:.4f}): {e}") from e


# =============================================================================
# CONFIDENCE METRICS (replace HMM entropy)
# =============================================================================

def compute_boundary_confidence(segments, df):
    """Welch t-test boundary confidence at each segment boundary.

    For adjacent segments on the same chromosome, compares the last ≤20 CN
    values of the left segment against the first ≤20 CN values of the right
    segment using Welch's t-test.  abs(t) is stored as 'boundary_conf':

        high boundary_conf → strong CN change → confident breakpoint
        low  boundary_conf → weak CN change   → possible noise breakpoint

    This is the analogue of HMM entropy (inverted: high entropy = uncertain
    boundary; high boundary_conf = certain boundary).
    """
    if not segments or not SCIPY_AVAILABLE:
        if not SCIPY_AVAILABLE:
            print("[BCONF] scipy not available — boundary_conf skipped")
        for seg in segments:
            seg['boundary_conf'] = None
        return segments

    # Build per-chromosome (starts, cn_values) sorted index for fast lookup
    idx = {}
    for chrom, gdf in df.groupby('chrom', sort=False):
        gdf_s = gdf.sort_values('start')
        idx[chrom] = (
            gdf_s['start'].values.astype(np.int64),
            gdf_s['cn'].values.astype(np.float64),
        )

    for seg in segments:
        seg['boundary_conf'] = None

    n_computed = 0
    for i in range(len(segments) - 1):
        left  = segments[i]
        right = segments[i + 1]
        if left['chrom'] != right['chrom']:
            continue

        chrom = left['chrom']
        if chrom not in idx:
            continue
        starts_arr, cn_arr = idx[chrom]

        def _get_cn(start, end):
            lo = int(np.searchsorted(starts_arr, start, side='left'))
            hi = int(np.searchsorted(starts_arr, end,   side='left'))
            return cn_arr[lo:hi]

        cn_l = _get_cn(left['start'],  left['end'])
        cn_r = _get_cn(right['start'], right['end'])
        cn_l = cn_l[-min(len(cn_l), 20):]
        cn_r = cn_r[:min(len(cn_r),  20)]

        if len(cn_l) < 2 or len(cn_r) < 2:
            continue

        stat, _ = ttest_ind(cn_l, cn_r, equal_var=False)
        if not np.isnan(stat):
            left['boundary_conf'] = float(abs(stat))
            n_computed += 1

    print(f"[BCONF] boundary_conf computed for {n_computed:,} / "
          f"{max(len(segments)-1, 0):,} adjacent pairs")
    return segments


def compute_segment_iqr(segments, df):
    """Add segment_iqr (IQR of raw CN values) to each segment.

    IQR is a robust spread measure within a segment.
        low  IQR → homogeneous segment → high CN confidence
        high IQR → heterogeneous → possible unsplit CN transition

    Stored as 'segment_iqr' in each segment dict.
    """
    idx = {}
    for chrom, gdf in df.groupby('chrom', sort=False):
        gdf_s = gdf.sort_values('start')
        idx[chrom] = (
            gdf_s['start'].values.astype(np.int64),
            gdf_s['cn'].values.astype(np.float64),
        )

    for seg in segments:
        chrom = seg['chrom']
        entry = idx.get(chrom)
        if entry is None:
            seg['segment_iqr'] = 0.0
            continue
        starts_arr, cn_arr = entry
        lo = int(np.searchsorted(starts_arr, seg['start'], side='left'))
        hi = int(np.searchsorted(starts_arr, seg['end'],   side='left'))
        vals = cn_arr[lo:hi]
        if len(vals) < 2:
            seg['segment_iqr'] = 0.0
        else:
            seg['segment_iqr'] = float(np.percentile(vals, 75) - np.percentile(vals, 25))

    return segments


# =============================================================================
# POST-PROCESSING HELPERS
# =============================================================================

def _pooled_std(s1, n1, s2, n2, m1, m2):
    """Combined std when merging two segment groups (exact parallel-groups formula).

    pooled_var = (within_SS + between_SS) / (n1 + n2 - 1)
    within_SS  = (n1-1)·s1² + (n2-1)·s2²
    between_SS = n1·n2 / (n1+n2) · (m1 − m2)²

    The between_SS term captures the CN divergence between the two groups —
    critical for split_high_cv_segments to detect merged segments that span a
    genuine CN transition (prevents cv=0 from masking the split opportunity).
    """
    n1 = max(int(n1), 1)
    n2 = max(int(n2), 1)
    within  = max(n1 - 1, 0) * s1 ** 2 + max(n2 - 1, 0) * s2 ** 2
    between = (n1 * n2) / (n1 + n2) * (m1 - m2) ** 2
    return float(np.sqrt(max(0.0, (within + between) / (n1 + n2 - 1))))


def _neutral_merge_dict(prev, seg):
    """Build merged Neutral segment dict using _pooled_std for cn_std.

    Applied consistently in: filter_small_segments, filter_low_quality_segments,
    reclassify_by_cn_threshold.  Fixes C7 (cn_std=0.0 hardcoded) and EK3
    (_pooled_std inconsistently applied).
    """
    total_bp = (prev['end'] - prev['start']) + (seg['end'] - seg['start'])
    pw = (prev['end'] - prev['start']) / total_bp
    cw = (seg['end']  - seg['start']) / total_bp
    return {
        'chrom':       prev['chrom'],
        'start':       prev['start'],
        'end':         seg['end'],
        'state':       'Neutral',
        'cn_median':   prev['cn_median']   * pw + seg['cn_median']   * cw,
        'cn_mean':     prev.get('cn_mean', prev['cn_median']) * pw
                       + seg.get('cn_mean', seg['cn_median']) * cw,
        'n_windows':   prev['n_windows'] + seg['n_windows'],
        'avg_quality': prev.get('avg_quality', 1.0) * pw + seg.get('avg_quality', 1.0) * cw,
        'min_quality': min(prev.get('min_quality', 1.0), seg.get('min_quality', 1.0)),
        # C7 / EK3 FIX: use _pooled_std instead of hardcoded 0.0
        'cn_std':      _pooled_std(
                           prev.get('cn_std', 0.0), prev['n_windows'],
                           seg.get('cn_std',  0.0), seg['n_windows'],
                           prev['cn_median'], seg['cn_median'],
                       ),
        'avg_repeats': prev.get('avg_repeats', 0.0) * pw + seg.get('avg_repeats', 0.0) * cw,
        'avg_entropy': 0.0,
        'max_entropy': 0.0,
    }


# =============================================================================
# POST-PROCESSING CHAIN
# =============================================================================

def filter_small_segments(segments, min_lengths):
    """Reclassify dup/del segments shorter than the state-specific threshold to Neutral.

    Adjacent Neutral segments are merged using _neutral_merge_dict (C7/EK3 fix).
    """
    if not min_lengths or not segments:
        return segments

    reclassified = 0
    result = []
    for seg in segments:
        state  = seg['state']
        length = seg['end'] - seg['start']
        min_len = min_lengths.get(state, 0)
        if min_len > 0 and length < min_len:
            seg = dict(seg)
            seg['state'] = 'Neutral'
            reclassified += 1
        result.append(seg)

    if reclassified:
        print(f"[POST] Reclassified {reclassified:,} small segments to Neutral")

    # Merge adjacent Neutral
    merged = [result[0]]
    for seg in result[1:]:
        prev = merged[-1]
        if (prev['state'] == 'Neutral' and seg['state'] == 'Neutral'
                and prev['chrom'] == seg['chrom']):
            merged[-1] = _neutral_merge_dict(prev, seg)
        else:
            merged.append(seg)

    if len(merged) != len(result):
        print(f"[POST] Merged adjacent Neutral: {len(result):,} → {len(merged):,}")
    return merged


def filter_low_quality_segments(segments, threshold, cv_filter_threshold=0.0):
    """Quality + CV filter: reclassify uncertain dup segments to Neutral.

    Two independent criteria (both applied):
      1. avg_quality < threshold        → low k-mer quality.
      2. cn_std/cn_median > cv_filter_threshold  → within-segment CN heterogeneity.

    Neutral merging uses _neutral_merge_dict (C7/EK3 fix).
    """
    if threshold <= 0 and cv_filter_threshold <= 0:
        return segments

    reclassified_q  = 0
    reclassified_cv = 0
    result = []
    for seg in segments:
        if seg['state'] in DUP_STATES:
            q   = seg.get('avg_quality', 1.0)
            cn_m = seg.get('cn_median', 1.0)
            cn_s = seg.get('cn_std',    0.0)
            cv   = cn_s / (cn_m + EPSILON)

            if threshold > 0 and q < threshold:
                seg = dict(seg)
                seg['state'] = 'Neutral'
                reclassified_q += 1
            elif cv_filter_threshold > 0 and cv > cv_filter_threshold:
                seg = dict(seg)
                seg['state'] = 'Neutral'
                reclassified_cv += 1
        result.append(seg)

    if reclassified_q:
        print(f"[POST] Quality filter (thr={threshold:.2f}): "
              f"reclassified {reclassified_q:,} low-quality dup segs to Neutral")
    if reclassified_cv:
        print(f"[POST] CV filter (thr={cv_filter_threshold:.2f}): "
              f"reclassified {reclassified_cv:,} high-CV dup segs to Neutral")

    if not result:
        return result

    # Merge adjacent Neutral using _neutral_merge_dict (C7 / EK3 fix)
    merged = [result[0]]
    for seg in result[1:]:
        prev = merged[-1]
        if (prev['state'] == 'Neutral' and seg['state'] == 'Neutral'
                and prev['chrom'] == seg['chrom']):
            merged[-1] = _neutral_merge_dict(prev, seg)
        else:
            merged.append(seg)

    if len(merged) != len(result):
        print(f"[POST] Merged adjacent Neutral after quality filter: "
              f"{len(result):,} → {len(merged):,}")
    return merged


def merge_nearby_dup_segments(segments, max_gap, merge_cn_tolerance=0.0):
    """Merge dup segments separated by a short Neutral gap.

    Pattern:  DupA | Neutral(short) | DupB  →  Dup(merged)

    C3 FIX: When merge_cn_tolerance > 0, allows DupA and DupB to be different
    states (e.g. HighDup + Amp) if |cn_median_A - cn_median_B| < tolerance.
    This captures biological SD blocks whose CN varies along the block but
    should not be split by a short Neutral gap artefact.

    Default merge_cn_tolerance=0.0 keeps strict same-state matching (CBS
    behaviour) for backward compatibility.  Set to 2.0 to enable cross-state
    merging for adjacent dup states.

    C4 note: Gap windows are intentionally excluded from merged statistics.
    The merged segment represents the dup signal; including CN≈1.0 gap
    windows would bias the dup CN estimate downward.  Coordinates span the
    gap (start=DupA.start, end=DupB.end).

    Uses _pooled_std for cn_std (EK3 consistent).
    """
    if max_gap <= 0 or len(segments) < 3:
        return segments

    total_merged = 0
    changed = True
    while changed:
        changed = False
        result = []
        i = 0
        while i < len(segments):
            can_merge = False
            if (i + 2 < len(segments)
                    and segments[i]['state'] in DUP_STATES
                    and segments[i + 1]['state'] == 'Neutral'
                    and segments[i + 2]['state'] in DUP_STATES
                    and segments[i]['chrom'] == segments[i + 2]['chrom']
                    and (segments[i + 1]['end'] - segments[i + 1]['start']) <= max_gap):

                s1, s2 = segments[i], segments[i + 2]

                if s1['state'] == s2['state']:
                    can_merge = True
                elif merge_cn_tolerance > 0:
                    # C3 FIX: allow cross-state merge when CN difference is small
                    can_merge = abs(s1['cn_median'] - s2['cn_median']) < merge_cn_tolerance

            if can_merge:
                s1, s2 = segments[i], segments[i + 2]
                # Merged state: use state of segment with more windows
                merged_state = s1['state'] if s1['n_windows'] >= s2['n_windows'] else s2['state']

                total_bp = (s1['end'] - s1['start']) + (s2['end'] - s2['start'])
                w1 = (s1['end'] - s1['start']) / total_bp
                w2 = (s2['end'] - s2['start']) / total_bp

                result.append({
                    'chrom':       s1['chrom'],
                    'start':       s1['start'],
                    'end':         s2['end'],
                    'state':       merged_state,
                    'cn_median':   s1['cn_median'] * w1 + s2['cn_median'] * w2,
                    'cn_mean':     s1.get('cn_mean', s1['cn_median']) * w1
                                   + s2.get('cn_mean', s2['cn_median']) * w2,
                    'n_windows':   s1['n_windows'] + s2['n_windows'],
                    'avg_quality': s1.get('avg_quality', 1.0) * w1
                                   + s2.get('avg_quality', 1.0) * w2,
                    'min_quality': min(s1.get('min_quality', 1.0),
                                       s2.get('min_quality', 1.0)),
                    # EK3: _pooled_std used consistently
                    'cn_std':      _pooled_std(
                                       s1.get('cn_std', 0.0), s1['n_windows'],
                                       s2.get('cn_std', 0.0), s2['n_windows'],
                                       s1['cn_median'], s2['cn_median'],
                                   ),
                    'avg_repeats': s1.get('avg_repeats', 0.0) * w1
                                   + s2.get('avg_repeats', 0.0) * w2,
                    'avg_entropy': 0.0,
                    'max_entropy': 0.0,
                })
                total_merged += 1
                changed = True
                i += 3
            else:
                result.append(segments[i])
                i += 1

        segments = result

    if total_merged:
        print(f"[POST] Merged {total_merged:,} nearby dup pairs (gap ≤ {max_gap} bp"
              + (f", CN tol={merge_cn_tolerance:.1f}" if merge_cn_tolerance > 0 else "") + ")")
    return segments


def split_high_cv_segments(segments, df, cv_threshold=0.6, min_length=3000, max_depth=3):
    """Split dup segments with high within-segment CN variability.

    CV = cn_std / cn_median.  If CV > cv_threshold, the segment likely spans
    a CN transition.  Splitting finds the variance-minimising split point and
    recurses up to max_depth times.

    Uses an O(n) prefix-sum approach to find the optimal split index.
    """
    _t0 = _time.time()

    idx = {}
    for chrom, gdf in df.groupby('chrom', sort=False):
        gdf_s = gdf.sort_values('start')
        raw_qws = (gdf_s['quality'].values.clip(0.01, 1.0).astype(np.float64)
                   if 'quality' in gdf_s.columns
                   else np.ones(len(gdf_s), dtype=np.float64))
        idx[chrom] = (
            gdf_s['start'].values.astype(np.int64),
            gdf_s['cn'].values.astype(np.float64),
            raw_qws,
        )

    def _find_best_split(chrom, seg_start, seg_end, gc_factor):
        entry = idx.get(chrom)
        if entry is None:
            return None
        starts, raw_cns, raw_qws = entry
        lo = int(np.searchsorted(starts, seg_start, side='left'))
        hi = int(np.searchsorted(starts, seg_end,   side='left'))
        n  = hi - lo
        if n < 4:
            return None

        win_cns = raw_cns[lo:hi]
        win_qws = raw_qws[lo:hi]
        # Quality-weighted prefix sums for split finding — consistent with
        # _make_sub which uses quality-weighted CN mean.  Without weighting,
        # low-quality outlier windows near a CN boundary skew the optimal
        # split point away from the true transition.
        w = win_qws.clip(0.01, 1.0)
        wsum  = np.cumsum(w)
        psw   = np.cumsum(w * win_cns)
        psw2  = np.cumsum(w * win_cns ** 2)

        best_score = float('inf')
        best_k     = -1
        for k in range(2, n - 2):
            left_end_bp  = int(starts[lo + k])
            if (left_end_bp  - seg_start) < min_length:
                continue
            if (seg_end - left_end_bp) < min_length:
                break

            Wl  = float(wsum[k - 1])
            SWYl  = float(psw[k - 1])
            SWY2l = float(psw2[k - 1])
            varl = (SWY2l - SWYl * SWYl / Wl) if Wl > 1e-12 else 0.0

            Wr    = float(wsum[n - 1]) - Wl
            SWYr  = float(psw[n - 1])  - SWYl
            SWY2r = float(psw2[n - 1]) - SWY2l
            varr  = (SWY2r - SWYr * SWYr / Wr) if Wr > 1e-12 else 0.0

            score = varl + varr   # weighted cost: Σwᵢ(yᵢ−ȳ)² both sides
            if score < best_score:
                best_score = score
                best_k     = k

        if best_k < 0:
            return None

        split_bp = int(starts[lo + best_k])

        def _make_sub(cn_slice, qw_slice, s, e, parent):
            """Build sub-segment with quality-weighted CN statistics.

            Uses the same quality weights as make_segment_from_windows so that
            CV = cn_std / cn_median is self-consistent within each sub-segment.
            avg_quality / min_quality are computed from the actual sub-segment
            windows (not inherited from the parent segment — Sorun #3 fix).
            """
            w = qw_slice.clip(0.01, 1.0)
            w_sum = float(w.sum())
            if w_sum < 1e-12:
                w = np.ones(len(cn_slice), dtype=np.float64)
                w_sum = float(len(cn_slice))
            raw_med = float(np.dot(w, cn_slice) / w_sum)   # quality-weighted mean
            raw_std = (float(np.sqrt(np.dot(w, (cn_slice - raw_med) ** 2) / w_sum))
                       if len(cn_slice) > 1 else 0.0)       # quality-weighted std
            cal_med = raw_med / gc_factor
            cal_std = raw_std / gc_factor
            sub = dict(parent)
            sub['start']       = s
            sub['end']         = e
            sub['cn_median']   = cal_med
            sub['cn_mean']     = float(np.mean(cn_slice)) / gc_factor
            sub['cn_std']      = cal_std
            sub['n_windows']   = len(cn_slice)
            sub['state']       = assign_state_from_cn(cal_med)
            sub['avg_quality'] = float(np.mean(qw_slice))
            sub['min_quality'] = float(np.min(qw_slice))
            sub['avg_entropy'] = 0.0
            sub['max_entropy'] = 0.0
            return sub

        left_cns  = win_cns[:best_k]
        right_cns = win_cns[best_k:]
        left_qws  = win_qws[:best_k]
        right_qws = win_qws[best_k:]
        return (
            _make_sub(left_cns,  left_qws,  seg_start, split_bp, {}),
            _make_sub(right_cns, right_qws, split_bp,  seg_end,  {}),
        )

    def _recursive_split(seg, depth):
        if depth >= max_depth or seg['state'] not in DUP_STATES:
            return [seg]
        cn_m = seg.get('cn_median', 1.0)
        cn_s = seg.get('cn_std',    0.0)
        if cn_m + EPSILON == 0:
            return [seg]
        cv = cn_s / (cn_m + EPSILON)
        if cv <= cv_threshold:
            return [seg]
        gc_factor = seg.get('gc_bias_factor', 1.0)
        result = _find_best_split(seg['chrom'], seg['start'], seg['end'], gc_factor)
        if result is None:
            return [seg]
        left, right = result
        for key, val in seg.items():
            left.setdefault(key, val)
            right.setdefault(key, val)
        left['start'],  left['end']  = seg['start'], result[0]['end']
        right['start'], right['end'] = result[1]['start'], seg['end']
        return _recursive_split(left, depth + 1) + _recursive_split(right, depth + 1)

    result      = []
    n_candidate = 0
    n_extra     = 0
    for seg in segments:
        cn_m = seg.get('cn_median', 1.0)
        cn_s = seg.get('cn_std',    0.0)
        cv   = cn_s / (cn_m + EPSILON)
        if seg['state'] not in DUP_STATES or cv <= cv_threshold:
            result.append(seg)
            continue
        n_candidate += 1
        sub = _recursive_split(seg, 0)
        n_extra += len(sub) - 1
        result.extend(sub)

    elapsed = _time.time() - _t0
    print(f"[CV-SPLIT] {n_candidate:,} high-CV dup segs → {n_extra:,} extra segs ({elapsed:.1f}s)")
    print(f"[CV-SPLIT] Total: {len(segments):,} → {len(result):,}")
    return result


# =============================================================================
# GC CALIBRATION
# =============================================================================

def compute_segment_mean_gc(segments, df, gc_df):
    """Compute mean GC content for each segment using O(S log W) binary search."""
    n_before = len(df)
    df_gc = df.merge(gc_df[['chrom', 'start', 'gc_content']], on=['chrom', 'start'], how='left')
    df_gc['gc_content'] = df_gc['gc_content'].fillna(0.5)
    assert len(df_gc) == n_before, "GC merge changed row count"

    print(f"[GC-CAL] Computing mean GC for {len(segments):,} segments...")

    gc_by_chrom = {}
    for chrom, grp in df_gc.groupby('chrom', sort=False):
        gs = grp.sort_values('start')
        gc_by_chrom[chrom] = (
            gs['start'].values.astype(np.int64),
            gs['gc_content'].values,
        )

    for seg in segments:
        entry = gc_by_chrom.get(seg['chrom'])
        if entry is None:
            seg['mean_gc'] = 0.5
            continue
        starts, gcs = entry
        lo = int(np.searchsorted(starts, seg['start'], side='left'))
        hi = int(np.searchsorted(starts, seg['end'],   side='left'))
        if hi > lo:
            seg['mean_gc'] = float(np.mean(gcs[lo:hi]))
        else:
            idx_c = min(int(np.searchsorted(starts, (seg['start'] + seg['end']) // 2)),
                        len(starts) - 1)
            seg['mean_gc'] = float(gcs[idx_c])

    print(f"[GC-CAL] Done.")
    return segments


def apply_gc_cn_calibration(segments):
    """Post-segmentation GC-based CN calibration.

    Fits a degree-3 polynomial to Neutral segments: observed_cn ~ mean_gc.
    For all non-Satellite, non-high-GC segments: cn_corrected = cn_median / poly(mean_gc).

    Bypass conditions (gc_bias_factor=1.0):
        - repeat_class in {'Satellite', 'Low_complexity'}
        - mean_gc > 0.60
        - masked_fraction > 0.85 AND mean_gc > 0.55
    """
    neutral_segs = [s for s in segments
                    if s['state'] == 'Neutral'
                    and 'mean_gc' in s
                    and 0.1 < s['mean_gc'] < 0.9]

    if len(neutral_segs) < 20:
        print(f"[GC-CAL] WARNING: only {len(neutral_segs)} neutral segs — skipping calibration")
        return segments

    gc_vals = np.array([s['mean_gc']   for s in neutral_segs])
    cn_vals = np.array([s['cn_median'] for s in neutral_segs])

    cn_ok = (cn_vals > 0.5) & (cn_vals < 2.0)
    if cn_ok.sum() < 20:
        print(f"[GC-CAL] WARNING: too few well-behaved neutral segs — skipping calibration")
        return segments

    gc_vals, cn_vals = gc_vals[cn_ok], cn_vals[cn_ok]
    poly_degree = 3
    coeffs = np.polyfit(gc_vals, cn_vals, poly_degree)
    poly   = np.poly1d(coeffs)

    gc_range_pred = poly(np.array([gc_vals.min(), gc_vals.max()]))
    residual_std  = float(np.std(cn_vals - poly(gc_vals)))
    bias_range    = float(gc_range_pred.max() - gc_range_pred.min())
    print(f"[GC-CAL] poly degree={poly_degree}, fitted on {len(gc_vals)} neutral segs")
    print(f"[GC-CAL] GC range: [{gc_vals.min():.3f}, {gc_vals.max():.3f}]  "
          f"bias range: [{gc_range_pred.min():.3f}, {gc_range_pred.max():.3f}]")
    print(f"[GC-CAL] Residual std: {residual_std:.4f}")

    # ONT data has negligible GC bias (typically <10%). Skip calibration
    # to avoid applying a near-identity polynomial that adds noise.
    if bias_range < 0.10:
        print(f"[GC-CAL] GC bias negligible (range={bias_range:.4f}), skipping calibration")
        return segments

    GC_BYPASS = frozenset({'Satellite', 'Low_complexity'})
    n_corrected = n_clamped = n_bypassed = 0

    for seg in segments:
        if 'mean_gc' not in seg:
            continue
        rc = seg.get('repeat_class', 'None')
        mf = seg.get('masked_fraction', 0.0)
        gc = seg.get('mean_gc', 0.5)

        if rc in GC_BYPASS or gc > 0.60 or (mf > 0.85 and gc > 0.55):
            seg['gc_bias_factor'] = 1.0
            n_bypassed += 1
            continue

        gc_factor_raw = float(poly(gc))
        gc_factor = max(0.6, min(1.8, gc_factor_raw))
        if gc_factor != gc_factor_raw:
            n_clamped += 1

        seg['cn_median'] = seg['cn_median'] / gc_factor
        seg['cn_mean']   = seg.get('cn_mean', seg['cn_median']) / gc_factor
        if seg.get('cn_std', 0.0) > 0:
            seg['cn_std'] = seg['cn_std'] / gc_factor
        seg['gc_bias_factor'] = gc_factor
        n_corrected += 1

    print(f"[GC-CAL] Applied to {n_corrected:,} segs ({n_bypassed:,} bypassed: Sat/high-GC)")
    if n_clamped:
        print(f"[GC-CAL] WARNING: {n_clamped} segs had gc_factor clamped to [0.6, 1.8]")
    return segments


# =============================================================================
# STATE RECLASSIFICATION
# =============================================================================

def reclassify_by_cn_threshold(segments, lowdup_threshold=1.25, hetdel_threshold=0.75):
    """Post-calibration CN-based state correction.

    Rules:
      LowDup → Neutral:  cn_median < lowdup_threshold
      HetDel → Neutral:  cn_median > hetdel_threshold  (GC-bias artefact)
      HetDel → LowDup:   cn_median > 1.5
      Neutral → Dup:     cn_median > lowdup_threshold   (HMM inertia artefact;
                          rare in PELT but possible after GC calibration)

    C5 FIX: Neutral merge updates cn_median and cn_std via _neutral_merge_dict
    (CBS already carries this fix forward — replicated here for completeness).
    EK3: _pooled_std called via _neutral_merge_dict.
    """
    recl = defaultdict(int)
    result = []
    for seg in segments:
        cn = seg.get('cn_median', 1.0)
        if seg['state'] == 'LowDup' and cn < lowdup_threshold:
            seg = dict(seg); seg['state'] = 'Neutral'; recl['ld→n'] += 1
        elif seg['state'] == 'HetDel' and cn > 1.5:
            seg = dict(seg); seg['state'] = 'LowDup';  recl['hd→ld'] += 1
        elif seg['state'] == 'HetDel' and cn > hetdel_threshold:
            seg = dict(seg); seg['state'] = 'Neutral'; recl['hd→n'] += 1
        elif seg['state'] == 'Neutral' and cn > lowdup_threshold:
            seg = dict(seg); seg['state'] = assign_state_from_cn(cn, lowdup_threshold)
            recl['n→dup'] += 1
        result.append(seg)

    for label, count in recl.items():
        print(f"[CN-RECL] Reclassified {count:,} segments: {label}")

    if not result:
        return result

    # C5 FIX: Merge adjacent Neutral using _neutral_merge_dict (updates cn_median + cn_std)
    merged = [result[0]]
    for seg in result[1:]:
        prev = merged[-1]
        if (prev['state'] == 'Neutral' and seg['state'] == 'Neutral'
                and prev['chrom'] == seg['chrom']
                and seg['start'] <= prev['end'] + 1):
            merged[-1] = _neutral_merge_dict(prev, seg)
        else:
            merged.append(seg)

    print(f"[CN-RECL] After merge: {len(result) - len(merged):,} adjacent Neutral merged")
    return merged


# =============================================================================
# REPEAT ANNOTATION
# =============================================================================

def compute_segment_repeat_annotation(segments, df, rm_df):
    """Attach per-segment repeat annotation (masked_fraction, repeat_class)."""
    CLASS_PRIORITY = ['Satellite', 'Simple_repeat', 'LINE', 'SINE', 'LTR',
                      'DNA', 'Low_complexity', 'Other', 'None']

    print(f"[RM] Computing repeat annotation for {len(segments):,} segments...")

    rm_by_chrom = {}
    for chrom, grp in rm_df.groupby('chrom'):
        gs = grp.sort_values('start')
        rm_by_chrom[chrom] = (
            gs['start'].values,
            gs['masked_fraction'].values,
            gs['repeat_class'].values,
        )

    for seg in segments:
        entry = rm_by_chrom.get(seg['chrom'])
        if entry is None:
            seg['masked_fraction'] = 0.0
            seg['repeat_class']    = 'None'
            continue

        starts_arr, mf_arr, cls_arr = entry
        lo = int(np.searchsorted(starts_arr, seg['start']))
        hi = int(np.searchsorted(starts_arr, seg['end']))
        if lo >= hi:
            seg['masked_fraction'] = 0.0
            seg['repeat_class']    = 'None'
            continue

        seg['masked_fraction'] = float(np.mean(mf_arr[lo:hi]))
        cls_counts = Counter(cls_arr[lo:hi])
        max_cnt = max(cls_counts.values())
        tied    = [c for c, n in cls_counts.items() if n == max_cnt]
        seg['repeat_class'] = min(
            tied,
            key=lambda c: CLASS_PRIORITY.index(c) if c in CLASS_PRIORITY else len(CLASS_PRIORITY)
        )

    n_masked = sum(1 for s in segments if s.get('masked_fraction', 0) > 0.5)
    print(f"[RM] Done. Segments with >50% masked: {n_masked:,} / {len(segments):,}")
    return segments


# =============================================================================
# OUTPUT
# =============================================================================

def write_output(segments, output_file, extended=True):
    """Write segments to extended 15-column BED file.

    Base columns (7):
        chrom, start, end, state, cn_median, cn_mean, n_windows

    Extended columns (+8):
        avg_quality, min_quality, cn_std, avg_repeats,
        avg_entropy (=0), max_entropy (=0),
        masked_fraction, repeat_class

    Optional extra columns (appended when present):
        gc_bias_factor, segment_iqr, boundary_conf

    Note on cn_median: this field stores the quality-weighted mean CN
    (not the true median). Name retained for backward compatibility with
    evaluate_ground_truth.py (see C2 in known_issues.md).
    """
    print(f"[IO] Writing {len(segments):,} segments to {output_file}...")

    # Check presence across ALL segments (not just first 20) to avoid silently
    # dropping columns when first N segments happen to be Neutral/unannotated.
    # Use short-circuit any() — stops at first hit, so still fast.
    has_repeat = any('repeat_class'   in s for s in segments)
    has_gc_fac = any('gc_bias_factor' in s for s in segments)
    has_iqr    = any('segment_iqr'    in s for s in segments)
    has_bconf  = any('boundary_conf'  in s for s in segments)

    with open(output_file, 'w') as f:
        if extended:
            hdr = ("#chrom\tstart\tend\tstate\tcn_median\tcn_mean\tn_windows\t"
                   "avg_quality\tmin_quality\tcn_std\tavg_repeats\t"
                   "avg_entropy\tmax_entropy")
            if has_repeat:
                hdr += "\tmasked_fraction\trepeat_class"
            if has_gc_fac:
                hdr += "\tgc_bias_factor"
            if has_iqr:
                hdr += "\tsegment_iqr"
            if has_bconf:
                hdr += "\tboundary_conf"
            f.write(hdr + "\n")
        else:
            f.write("#chrom\tstart\tend\tstate\tcn_median\tcn_mean\tn_windows\n")

        for seg in segments:
            line = (f"{seg['chrom']}\t{seg['start']}\t{seg['end']}\t"
                    f"{seg['state']}\t{seg['cn_median']:.4f}\t"
                    f"{seg.get('cn_mean', seg['cn_median']):.4f}\t{seg['n_windows']}")
            if extended:
                line += (f"\t{seg.get('avg_quality', 1.0):.4f}"
                         f"\t{seg.get('min_quality', 1.0):.4f}"
                         f"\t{seg.get('cn_std',      0.0):.4f}"
                         f"\t{seg.get('avg_repeats', 0.0):.2f}"
                         f"\t0.0000"   # avg_entropy (always 0 — no HMM)
                         f"\t0.0000")  # max_entropy (always 0 — no HMM)
                if has_repeat:
                    line += (f"\t{seg.get('masked_fraction', 0.0):.4f}"
                             f"\t{seg.get('repeat_class', 'None')}")
                if has_gc_fac:
                    line += f"\t{seg.get('gc_bias_factor', 1.0):.4f}"
                if has_iqr:
                    line += f"\t{seg.get('segment_iqr', 0.0):.4f}"
                if has_bconf:
                    bc = seg.get('boundary_conf')
                    line += f"\t{bc:.4f}" if bc is not None else "\t"
            f.write(line + "\n")

    print("[IO] Done.")


def print_statistics(segments):
    """Print per-state summary statistics."""
    print("\n" + "=" * 60)
    print("WEIGHTED PELT SEGMENTATION STATISTICS")
    print("=" * 60)

    state_counts = defaultdict(int)
    state_bases  = defaultdict(int)
    state_cn_sum = defaultdict(float)

    for seg in segments:
        state  = seg['state']
        length = seg['end'] - seg['start']
        state_counts[state] += 1
        state_bases[state]  += length
        state_cn_sum[state] += seg['cn_median'] * length

    total_segs  = len(segments)
    total_bases = sum(state_bases.values())

    print(f"\nTotal segments: {total_segs:,}")
    print(f"Total bases:    {total_bases:,}")
    print("\nState Distribution:")
    print("-" * 60)
    print(f"{'State':<12} {'Segments':>10} {'Pct':>8} {'Bases':>15} {'Mean CN':>10}")
    print("-" * 60)
    for state_name in STATE_ORDER:
        count = state_counts.get(state_name, 0)
        bases = state_bases.get(state_name, 0)
        pct   = 100.0 * count / total_segs if total_segs > 0 else 0.0
        mean_cn = state_cn_sum.get(state_name, 0) / bases if bases > 0 else 0.0
        print(f"{state_name:<12} {count:>10,} {pct:>7.1f}% {bases:>15,} {mean_cn:>10.2f}")
    print("-" * 60)


# =============================================================================
# MAIN
# =============================================================================

def main():
    if not RUPTURES_AVAILABLE:
        print("[ERROR] ruptures is not installed. Run: pip install ruptures>=1.1.7")
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description="KonuSeg Weighted PELT (Fused Lasso) CN Segmentation — "
                    "drop-in alternative to HMM and CBS"
    )

    # --- Required ---
    parser.add_argument("--input",  "-i", required=True,
                        help="Input BED from preprocess_kmer_windows.py "
                             "(8-col: chrom start end cn mean_count log_ratio "
                             "num_kmers num_filtered)")
    parser.add_argument("--output", "-o", required=True,
                        help="Output segments BED file")

    # --- Algorithm ---
    parser.add_argument("--no-weights", action="store_true",
                        help="Disable quality weighting in PELT cost. "
                             "Falls back to standard L2 (CBS-equivalent). "
                             "Default: quality weights enabled.")
    parser.add_argument("--penalty", type=float, default=None,
                        help="PELT penalty override. Default: BIC = 2·log(n) "
                             "per chromosome. Smaller → more segments; larger → fewer.")
    parser.add_argument("--penalty-factor", type=float, default=4.0,
                        help="BIC penalty multiplier. Default: 4.0 (2× standard BIC). "
                             "Larger values → fewer, longer segments. "
                             "Standard BIC = 2.0. Recommended: 4-8 for ONT data.")

    # --- Post-processing ---
    parser.add_argument("--min-kmers", type=int, default=0,
                        help="Min num_kmers per window. Windows below threshold "
                             "→ CN=1.0 (Neutral prior). Recommended: 30.")
    parser.add_argument("--min-segment-length", type=int, default=3000,
                        help="Minimum segment length in bp. Default: 3000.")
    parser.add_argument("--quality-threshold", type=float, default=0.7,
                        help="Hard quality filter: dup segments with avg_quality < "
                             "threshold → Neutral. Default: 0.7. Set 0 to disable.")
    parser.add_argument("--cv-filter-threshold", type=float, default=0.0,
                        help="CV filter: dup segments with cn_std/cn_median > "
                             "threshold → Neutral. Default: 0.0 (disabled). "
                             "Useful range: 0.4–0.8.")
    parser.add_argument("--cv-split-threshold", type=float, default=0.6,
                        help="CV split: dup segments with CV > threshold are split "
                             "at the variance-minimising point. Default: 0.6.")
    parser.add_argument("--merge-cn-tolerance", type=float, default=0.0,
                        help="C3 fix: allow cross-state dup merge when "
                             "|cn_median_A - cn_median_B| < tolerance. "
                             "Default: 0.0 (strict same-state). Set to 2.0 to "
                             "enable HighDup+Amp merging for biological SD blocks.")
    parser.add_argument("--max-merge-gap", type=int, default=10000,
                        help="Max Neutral gap length (bp) for merge_nearby_dup_segments. "
                             "Default: 10000 (matches CLAUDE.md v11 target).")
    parser.add_argument("--lowdup-threshold", type=float, default=1.25,
                        help="CN threshold for Neutral/LowDup boundary in "
                             "reclassify_by_cn_threshold. Default: 1.25 (backward "
                             "compat). Set to 1.52 for HMM log2-midpoint consistency (C8).")
    parser.add_argument("--cn-reclassify-threshold", type=float, default=1.25,
                        help="Deprecated alias for --lowdup-threshold. "
                             "Will be removed in future versions.")

    # --- Optional features ---
    parser.add_argument("--extended", action="store_true",
                        help="Write extended 15-col output. Always enabled for "
                             "compatibility with evaluate_ground_truth.py; "
                             "flag kept for CLI parity with HMM.")
    parser.add_argument("--gc-content-bed",
                        help="GC content BED from compute_gc_content.py. "
                             "Enables post-segmentation GC CN calibration.")
    parser.add_argument("--repeat-bed",
                        help="Repeat-annotated BED from compute_repeat_annotation.py. "
                             "Adds masked_fraction + repeat_class. Does NOT change "
                             "breakpoints — metadata only.")
    parser.add_argument("--no-boundary-conf", action="store_true",
                        help="Skip Welch t-test boundary confidence computation "
                             "(saves ~5s on full genome). Default: computed.")
    parser.add_argument("--no-segment-iqr", action="store_true",
                        help="Skip per-segment IQR computation. Default: computed.")

    # --- HMM/CBS CLI compatibility ---
    parser.add_argument("--mode", default=None,
                        choices=['precision', 'sensitive', 'cn-accuracy'],
                        help="Accepted for CLI compatibility. Weighted PELT always "
                             "runs in cn-accuracy equivalent mode.")
    parser.add_argument("--no-chrx-correction", action="store_true",
                        help="Accepted for CLI compatibility (no-op — PELT uses "
                             "preprocess normalization).")
    parser.add_argument("--no-pomegranate", action="store_true",
                        help="Accepted for CLI compatibility (no-op — PELT does "
                             "not use pomegranate).")
    parser.add_argument("--coarse-output",
                        help="Accepted for CLI compatibility (not implemented).")

    args = parser.parse_args()

    # Handle deprecated alias
    lowdup_threshold = args.lowdup_threshold
    if args.cn_reclassify_threshold != 1.25:
        lowdup_threshold = args.cn_reclassify_threshold
        print(f"[INFO] --cn-reclassify-threshold is deprecated; use --lowdup-threshold")

    use_weights = not args.no_weights

    print("=" * 60)
    print("KonuSeg Weighted PELT (Fused Lasso) Segmentation")
    print("=" * 60)
    print(f"Input:               {args.input}")
    print(f"Output:              {args.output}")
    print(f"Quality weights:     {'yes (weighted L2 cost)' if use_weights else 'no (standard L2)'}")
    print(f"Penalty:             {'BIC per chrom (' + str(args.penalty_factor) + '·log n)' if args.penalty is None else args.penalty}")
    print(f"Penalty factor:      {args.penalty_factor}")
    print(f"Min segment length:  {args.min_segment_length} bp")
    print(f"Min kmers:           {args.min_kmers}")
    print(f"Quality threshold:   {args.quality_threshold}")
    print(f"CV filter threshold: {args.cv_filter_threshold}")
    print(f"CV split threshold:  {args.cv_split_threshold}")
    print(f"Merge CN tolerance:  {args.merge_cn_tolerance} CN units")
    print(f"Max merge gap:       {args.max_merge_gap} bp")
    print(f"LowDup threshold:    {lowdup_threshold}")
    print(f"GC calibration:      {'yes' if args.gc_content_bed else 'no'}")
    print(f"Repeat annotation:   {'yes' if args.repeat_bed else 'no'}")
    print()

    # -------------------------------------------------------------------------
    # Step 1: Load data
    # -------------------------------------------------------------------------
    df = load_data(args.input, min_kmers=args.min_kmers)

    gc_df = None
    if args.gc_content_bed:
        print()
        gc_df = load_gc_content(args.gc_content_bed)

    # -------------------------------------------------------------------------
    # Step 2: Weighted PELT segmentation
    # -------------------------------------------------------------------------
    print()
    segments = run_segmentation(
        df,
        penalty=args.penalty,
        min_size=args.min_segment_length,
        use_weights=use_weights,
        lowdup_threshold=lowdup_threshold,
        penalty_factor=args.penalty_factor,
    )

    # -------------------------------------------------------------------------
    # Step 3: Post-processing chain (ORDER MATTERS — mirrors cn-accuracy HMM)
    # -------------------------------------------------------------------------
    print()

    # 3a. Filter small segments
    min_lengths = {
        'HomDel':     args.min_segment_length,
        'HetDel':     args.min_segment_length,
        'LowDup':     args.min_segment_length,
        'HighDup':    args.min_segment_length,
        'Amp':        args.min_segment_length,
        'MedAmp':     args.min_segment_length,
        'HighAmp':    args.min_segment_length,
        'ExtremeAmp': args.min_segment_length,
    }
    segments = filter_small_segments(segments, min_lengths)

    # 3b. Quality + CV filter
    segments = filter_low_quality_segments(
        segments,
        threshold=args.quality_threshold,
        cv_filter_threshold=args.cv_filter_threshold,
    )

    # 3c. Merge nearby dup segments (C3 fix: merge_cn_tolerance)
    segments = merge_nearby_dup_segments(
        segments,
        max_gap=args.max_merge_gap,
        merge_cn_tolerance=args.merge_cn_tolerance,
    )

    # 3d. CV-based splitting (before GC calibration — operates on raw CN)
    if args.cv_split_threshold > 0:
        print()
        segments = split_high_cv_segments(
            segments, df,
            cv_threshold=args.cv_split_threshold,
            min_length=args.min_segment_length,
        )

    # 3d2. Second merge pass: re-merge fragments created by CV-split
    if args.cv_split_threshold > 0 and args.max_merge_gap > 0:
        segments = merge_nearby_dup_segments(
            segments,
            max_gap=args.max_merge_gap,
            merge_cn_tolerance=args.merge_cn_tolerance,
        )

    # 3e. Repeat annotation (BEFORE GC calibration — bypass uses repeat_class)
    # apply_gc_cn_calibration() checks repeat_class for Satellite/Low_complexity
    # bypass; if repeat annotation runs after GC calibration those bypass conditions
    # are never triggered (repeat_class='None' by default) → GC correction applied
    # to centromeric Satellite arrays where polynomial extrapolation is invalid.
    if args.repeat_bed:
        print()
        rm_df = load_repeat_annotation(args.repeat_bed)
        segments = compute_segment_repeat_annotation(segments, df, rm_df)

    # 3f. GC calibration (post-merge, post-split, post-repeat-annotation)
    if gc_df is not None:
        print()
        segments = compute_segment_mean_gc(segments, df, gc_df)
        segments = apply_gc_cn_calibration(segments)

    # 3g. CN-based state reclassification (after GC calibration)
    print()
    segments = reclassify_by_cn_threshold(
        segments,
        lowdup_threshold=lowdup_threshold,
    )

    # -------------------------------------------------------------------------
    # Step 4: Confidence metrics
    # -------------------------------------------------------------------------
    print()
    if not args.no_boundary_conf:
        segments = compute_boundary_confidence(segments, df)
    if not args.no_segment_iqr:
        segments = compute_segment_iqr(segments, df)
        print(f"[IQR] Segment IQR computed for {len(segments):,} segments")

    # -------------------------------------------------------------------------
    # Step 5: Output
    # -------------------------------------------------------------------------
    print_statistics(segments)
    write_output(segments, args.output, extended=True)


if __name__ == '__main__':
    main()

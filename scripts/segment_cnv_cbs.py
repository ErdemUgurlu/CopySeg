#!/usr/bin/env python3
"""
KonuSeg CBS/PELT Segmentation
==============================
PELT-based drop-in alternative to segment_cnv_hmm_log_7state.py.

Key differences from the HMM approach:
  - No emission model → B1 (wrong variances), B2 (state gap CN=9-31), B3 (no EM),
    B4 (index-based transitions), B5 (1 Mb min segment), EK2 (wasted states) all resolved.
  - ruptures PELT (L2 model) detects breakpoints in log2 space — data-adaptive,
    no state count limitation, O(n) per chromosome.
  - Post-hoc CN assignment via CN_STATE_BOUNDARIES (same as HMM).
  - Identical post-processing chain: filter → merge → GC-cal → reclassify → split-CV.
  - Identical CLI interface and extended 15-column BED output.
  - Additional column: boundary_conf (Welch t-test |t| at segment boundary).

Usage:
    python3 segment_cnv_cbs.py --input windows.bed --output segs.bed --extended \\
        --gc-content-bed gc.bed --repeat-bed rm.bed --min-kmers 30

Install ruptures:
    pip install ruptures>=1.1.7

Output columns (extended, 15+2):
    chrom, start, end, state, cn_median, cn_mean, n_windows,
    avg_quality, min_quality, cn_std, avg_repeats,
    avg_entropy (=0 for CBS), max_entropy (=0 for CBS),
    masked_fraction, repeat_class,
    gc_bias_factor (optional, 16th),
    boundary_conf (optional, 17th)
"""

import argparse
import sys
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd

try:
    import ruptures as rpt
except ImportError:
    print("[ERROR] ruptures is not installed. Run: pip install ruptures>=1.1.7")
    sys.exit(1)

try:
    from scipy.stats import ttest_ind
except ImportError:
    print("[ERROR] scipy is not installed. Run: pip install scipy>=1.7.0")
    sys.exit(1)

warnings.filterwarnings('ignore')

# =============================================================================
# CONSTANTS
# =============================================================================

EPSILON = 1e-3       # Prevents log2(0); same as HMM
WINDOW_BP = 500      # Expected window size in bp (CHM13 ONT pipeline default)

# CN → state mapping (identical to HMM CN_STATE_BOUNDARIES)
CN_STATE_BOUNDARIES = [
    (0.30,         'HomDel'),
    (0.70,         'HetDel'),
    (1.25,         'Neutral'),
    (3.00,         'LowDup'),
    (6.00,         'HighDup'),
    (12.00,        'Amp'),
    (float('inf'), 'HighAmp'),
]

STATE_ORDER = ['HomDel', 'HetDel', 'Neutral', 'LowDup', 'HighDup', 'Amp', 'HighAmp']


# =============================================================================
# STATE ASSIGNMENT
# =============================================================================

def assign_state_from_cn(cn):
    """Map a calibrated CN value to the nearest state using CN_STATE_BOUNDARIES."""
    for upper, state in CN_STATE_BOUNDARIES:
        if cn < upper:
            return state
    return 'HighAmp'


# =============================================================================
# DATA LOADING
# =============================================================================

def load_data(bed_file, min_kmers=0):
    """Load 8-column windows BED and compute log2_cn + quality.

    CBS always uses quality_weight=0: quality is computed for post-filtering
    purposes only, never used to modify the CN signal.

    min_kmers: windows with num_kmers < min_kmers are set to cn=1.0 (Neutral prior)
    to suppress false HetDel calls in low-coverage regions.
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

        # Quality score (signal not modified — CBS always QW=0)
        if 'num_kmers' in df.columns and 'num_filtered' in df.columns:
            df['quality'] = (df['num_kmers'] /
                             (df['num_kmers'] + df['num_filtered'] + 1))
        else:
            df['quality'] = 1.0

        # Log2 CN for PELT segmentation (raw, no quality weighting)
        df['log2_cn'] = np.log2(df['cn'].values + EPSILON)

        # Remove inf/nan
        before = len(df)
        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.dropna(subset=['log2_cn'])
        after = len(df)
        if before != after:
            print(f"[IO] Removed {before - after} invalid rows")

        print(f"[IO] CN range: {df['cn'].min():.4f} - {df['cn'].max():.4f}")
        print(f"[IO] Log2 range: {df['log2_cn'].min():.2f} - {df['log2_cn'].max():.2f}")
        print(f"[IO] Median CN: {df['cn'].median():.4f}")

        return df

    except Exception as e:
        print(f"[ERROR] Failed to load data: {e}")
        sys.exit(1)


def load_gc_content(gc_bed_file):
    """Load GC content per window from BED file (chrom, start, end, gc_content)."""
    print(f"[GC] Loading GC content from {gc_bed_file}...")
    gc_df = pd.read_csv(gc_bed_file, sep='\t', comment='#', header=None)
    if len(gc_df.columns) < 4:
        raise ValueError(f"GC BED needs >= 4 columns, got {len(gc_df.columns)}")
    gc_df.columns = ['chrom', 'start', 'end', 'gc_content'][:len(gc_df.columns)]
    gc_df = gc_df[['chrom', 'start', 'gc_content']]
    print(f"[GC] Loaded {len(gc_df):,} windows (gc mean={gc_df['gc_content'].mean():.4f})")
    return gc_df


def load_repeat_annotation(repeat_bed_file):
    """Load per-window repeat annotation from compute_repeat_annotation.py output.

    Input: 10-column BED (chrom, start, end, cn, mean_count, log_ratio,
                           num_kmers, num_filtered, masked_fraction, repeat_class)
    Returns: DataFrame with (chrom, start, masked_fraction, repeat_class)
    """
    print(f"[RM] Loading repeat annotation from {repeat_bed_file}...")
    rm_df = pd.read_csv(
        repeat_bed_file, sep='\t', comment='#', header=None,
        names=['chrom', 'start', 'end', 'cn', 'mean_count', 'log_ratio',
               'num_kmers', 'num_filtered', 'masked_fraction', 'repeat_class'],
        usecols=[0, 1, 8, 9],
    )
    rm_df.columns = ['chrom', 'start', 'masked_fraction', 'repeat_class']
    pct_masked = (rm_df['masked_fraction'] > 0.5).mean() * 100
    print(f"[RM] Loaded {len(rm_df):,} windows (>50% masked: {pct_masked:.1f}%)")
    return rm_df


# =============================================================================
# CBS/PELT CORE SEGMENTATION
# =============================================================================

def make_segment_from_windows(windows_df):
    """Build a segment dict from a slice of the windows DataFrame.

    Args:
        windows_df: DataFrame slice for one segment (sorted by start).
                    Required columns: chrom, start, end, cn.
                    Optional: quality, num_filtered.

    Returns:
        dict with all keys expected by post-processing functions.
        avg_entropy and max_entropy are always 0.0 (CBS has no HMM entropy).
    """
    cn_arr = windows_df['cn'].values
    n = len(cn_arr)

    cn_median = float(np.median(cn_arr))
    cn_mean = float(np.mean(cn_arr))
    cn_std = float(np.std(cn_arr, ddof=1)) if n > 1 else 0.0

    q_arr = windows_df['quality'].values if 'quality' in windows_df.columns else np.ones(n)
    avg_quality = float(np.mean(q_arr))
    min_quality = float(np.min(q_arr))

    avg_repeats = float(windows_df['num_filtered'].mean()) \
        if 'num_filtered' in windows_df.columns else 0.0

    return {
        'chrom':       str(windows_df['chrom'].iloc[0]),
        'start':       int(windows_df['start'].iloc[0]),
        'end':         int(windows_df['end'].iloc[-1]),
        'state':       assign_state_from_cn(cn_median),
        'cn_median':   cn_median,
        'cn_mean':     cn_mean,
        'n_windows':   n,
        'avg_quality': avg_quality,
        'min_quality': min_quality,
        'cn_std':      cn_std,
        'avg_repeats': avg_repeats,
        'avg_entropy': 0.0,   # CBS has no HMM entropy; kept for downstream compatibility
        'max_entropy': 0.0,
    }


def run_pelt_chrom(df_chrom, penalty, min_size_windows):
    """Segment one chromosome's log2_cn signal using PELT (ruptures L2 model).

    PELT (Pruned Exact Linear Time) finds the globally optimal set of breakpoints
    for the L2 (squared error) cost function. Unlike HMM, it:
      - Requires no emission model → solves B1, B2, B3, B4
      - Controls granularity via penalty (not self_transition) → B5 mitigated
      - Has no fixed state count → EK2 resolved

    Args:
        df_chrom: DataFrame for one chromosome, sorted by start, index reset.
        penalty: float — PELT penalty (BIC = 2*log(n) if computed by caller).
        min_size_windows: int — minimum segment length in window units.

    Returns:
        list of segment dicts.
    """
    signal = df_chrom['log2_cn'].values
    n = len(signal)

    if n == 0:
        return []

    # Single-segment fallback when chromosome is too short for any split
    if n < max(min_size_windows * 2, 4):
        return [make_segment_from_windows(df_chrom)]

    # PELT with L2 cost (change in mean), jump=1 for exact breakpoint placement
    algo = rpt.Pelt(model='l2', min_size=min_size_windows, jump=1)
    breakpoints = algo.fit_predict(signal, pen=penalty)
    # breakpoints: sorted list of end-indices (exclusive), last element == n

    segments = []
    prev = 0
    for bp in breakpoints:
        if bp <= prev:
            continue
        seg = make_segment_from_windows(df_chrom.iloc[prev:bp])
        if seg is not None:
            segments.append(seg)
        prev = bp

    return segments


def run_segmentation(df, penalty=None, min_size=3000):
    """Run PELT segmentation on all chromosomes.

    Args:
        df: full windows DataFrame (all chromosomes)
        penalty: float or None. None → BIC = 2*log(n) per chromosome.
        min_size: minimum segment length in bp.

    Returns:
        list of all segment dicts across all chromosomes.
    """
    import time as _t

    # Estimate window_bp defensively (should be 500 for standard CHM13 pipeline)
    window_sizes = (df['end'] - df['start']).values
    window_bp = int(np.median(window_sizes)) if len(window_sizes) > 0 else WINDOW_BP
    if window_bp != WINDOW_BP:
        print(f"[CBS] Non-standard window size detected: {window_bp} bp (expected {WINDOW_BP})")

    min_size_windows = max(1, min_size // window_bp)
    print(f"[CBS] min_segment_length={min_size} bp → min_size_windows={min_size_windows}")
    print(f"[CBS] PELT penalty: {'BIC (auto)' if penalty is None else penalty}")

    chroms = list(df['chrom'].unique())  # preserve input order
    all_segments = []
    total_t0 = _t.time()

    for chrom in chroms:
        t0 = _t.time()
        df_chrom = df[df['chrom'] == chrom].sort_values('start').reset_index(drop=True)
        n = len(df_chrom)

        # Per-chromosome BIC penalty
        pen = penalty if penalty is not None else (2.0 * np.log(max(n, 2)))

        segs = run_pelt_chrom(df_chrom, pen, min_size_windows)
        all_segments.extend(segs)

        elapsed = _t.time() - t0
        print(f"[CBS] {chrom}: {n:,} windows → {len(segs):,} segments ({elapsed:.1f}s)")

    total_elapsed = _t.time() - total_t0
    print(f"[CBS] Total: {len(df):,} windows → {len(all_segments):,} segments "
          f"({total_elapsed:.1f}s)")
    return all_segments


# =============================================================================
# BOUNDARY CONFIDENCE (CBS-specific, replaces HMM entropy)
# =============================================================================

def compute_boundary_confidence(segments, df):
    """Compute Welch's t-test boundary confidence at each segment boundary.

    For each pair of adjacent segments on the same chromosome, compares
    the last min(n_left_windows, 20) CN values of the left segment against
    the first min(n_right_windows, 20) CN values of the right segment.

    Stores abs(t) as 'boundary_conf' on the LEFT segment of each pair.
    First and last segments of each chromosome receive boundary_conf=None.

    Higher boundary_conf = clearer CN change at the boundary = more confident
    the breakpoint is real. This is the CBS analogue of HMM entropy (inverted:
    high HMM entropy means uncertain boundary; high boundary_conf means certain).

    Args:
        segments: list of segment dicts (chrom, start, end).
        df: original windows DataFrame (columns: chrom, start, cn).

    Returns:
        segments (mutated in-place, 'boundary_conf' key added to each segment).
    """
    if not segments:
        return segments

    # Pre-build per-chromosome index for O(log n) window lookup
    idx = {}
    for chrom, gdf in df.groupby('chrom', sort=False):
        gdf_s = gdf.sort_values('start')
        idx[chrom] = (
            gdf_s['start'].values.astype(np.int64),
            gdf_s['cn'].values.astype(np.float64),
        )

    # Default: no confidence (first/last segment of each chrom)
    for seg in segments:
        seg['boundary_conf'] = None

    n_computed = 0
    for i in range(len(segments) - 1):
        left  = segments[i]
        right = segments[i + 1]
        if left['chrom'] != right['chrom']:
            continue  # chromosome boundary — no t-test

        chrom = left['chrom']
        if chrom not in idx:
            continue

        starts_arr, cn_arr = idx[chrom]

        # Extract CN values for left segment
        lo_l = int(np.searchsorted(starts_arr, left['start'], side='left'))
        hi_l = int(np.searchsorted(starts_arr, left['end'],   side='left'))
        cn_left = cn_arr[lo_l:hi_l]

        # Extract CN values for right segment
        lo_r = int(np.searchsorted(starts_arr, right['start'], side='left'))
        hi_r = int(np.searchsorted(starts_arr, right['end'],   side='left'))
        cn_right = cn_arr[lo_r:hi_r]

        # Use tail/head windows only (max 20 each)
        window_n = 20
        cn_left_slice  = cn_left[-min(len(cn_left),  window_n):]
        cn_right_slice = cn_right[:min(len(cn_right), window_n)]

        if len(cn_left_slice) < 2 or len(cn_right_slice) < 2:
            continue

        stat, _ = ttest_ind(cn_left_slice, cn_right_slice, equal_var=False)
        if not np.isnan(stat):
            left['boundary_conf'] = float(abs(stat))
            n_computed += 1

    print(f"[BCONF] Boundary confidence computed for {n_computed:,} / "
          f"{len(segments) - 1:,} adjacent pairs")
    return segments


# =============================================================================
# POST-PROCESSING HELPERS
# =============================================================================

def _pooled_std(s1, n1, s2, n2, m1, m2):
    """Combined standard deviation when merging two segment groups.

    Uses the exact parallel-groups formula:
      pooled_var = (within_SS + between_SS) / (n1 + n2 - 1)
      within_SS  = (n1-1)*s1² + (n2-1)*s2²
      between_SS = n1*n2 / (n1+n2) * (m1 - m2)²

    Critical for cv-split pre-filter: merged segments with different CN values
    must have non-zero cn_std so split_high_cv_segments() can detect them.
    """
    n1, n2 = max(int(n1), 1), max(int(n2), 1)
    within  = max(n1 - 1, 0) * s1 ** 2 + max(n2 - 1, 0) * s2 ** 2
    between = (n1 * n2) / (n1 + n2) * (m1 - m2) ** 2
    return float(np.sqrt(max(0.0, (within + between) / (n1 + n2 - 1))))


# =============================================================================
# POST-PROCESSING CHAIN (verbatim from segment_cnv_hmm_log_7state.py unless noted)
# =============================================================================

def filter_small_segments(segments, min_lengths):
    """Reclassify small dup/del segments as Neutral.

    Segments shorter than the state-specific threshold are converted to Neutral.
    Adjacent Neutral segments are then merged.
    """
    if not min_lengths or not segments:
        return segments

    reclassified = 0
    result = []
    for seg in segments:
        state = seg['state']
        length = seg['end'] - seg['start']
        min_len = min_lengths.get(state, 0)

        if min_len > 0 and length < min_len:
            seg = dict(seg)
            seg['state'] = 'Neutral'
            reclassified += 1

        result.append(seg)

    if reclassified > 0:
        print(f"[POST] Reclassified {reclassified} small segments to Neutral")

    # Merge adjacent Neutral segments
    merged = [result[0]]
    for seg in result[1:]:
        prev = merged[-1]
        if (seg['state'] == 'Neutral' and prev['state'] == 'Neutral' and
                seg['chrom'] == prev['chrom']):
            total_bp = (prev['end'] - prev['start']) + (seg['end'] - seg['start'])
            pw = (prev['end'] - prev['start']) / total_bp
            cw = (seg['end'] - seg['start']) / total_bp
            merged[-1] = {
                'chrom':       prev['chrom'],
                'start':       prev['start'],
                'end':         seg['end'],
                'state':       'Neutral',
                'cn_median':   prev['cn_median'] * pw + seg['cn_median'] * cw,
                'cn_mean':     prev['cn_mean'] * pw + seg['cn_mean'] * cw,
                'n_windows':   prev['n_windows'] + seg['n_windows'],
                'avg_entropy': prev.get('avg_entropy', 0.0) * pw + seg.get('avg_entropy', 0.0) * cw,
                'max_entropy': max(prev.get('max_entropy', 0.0), seg.get('max_entropy', 0.0)),
            }
        else:
            merged.append(seg)

    if len(merged) != len(result):
        print(f"[POST] Merged adjacent Neutral: {len(result)} -> {len(merged)}")

    return merged


def filter_low_quality_segments(segments, threshold, cv_filter_threshold=0.0):
    """Quality + CV filter: reclassify uncertain dup segments to Neutral.

    Two independent criteria (both applied):
      1. avg_quality < threshold  — k-mer quality too low.
      2. cn_std / cn_median > cv_filter_threshold (if > 0.0) — high within-segment
         CV indicates the segment spans a CN transition. This is the CBS analogue
         of the HMM entropy filter (which is unavailable without a posterior).
         Default cv_filter_threshold=0.0 disables this criterion.

    avg_entropy and max_entropy are always 0.0 for CBS segments. They are preserved
    as zeros in all merge operations to maintain compatibility with downstream tools
    (validate_cn_accuracy.py, evaluate_ground_truth.py).
    """
    if threshold <= 0 and cv_filter_threshold <= 0:
        return segments

    dup_states = {'LowDup', 'HighDup', 'Amp', 'HighAmp'}
    reclassified_q = 0
    reclassified_cv = 0
    result = []
    for seg in segments:
        if seg['state'] in dup_states:
            q   = seg.get('avg_quality', 1.0)
            cn_m = seg.get('cn_median', 1.0)
            cn_s = seg.get('cn_std', 0.0)
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

    if reclassified_q > 0:
        print(f"[POST] Quality filter (threshold={threshold:.2f}): "
              f"reclassified {reclassified_q} low-quality dup segments to Neutral")
    if reclassified_cv > 0:
        print(f"[POST] CV filter (cv_threshold={cv_filter_threshold:.2f}): "
              f"reclassified {reclassified_cv} high-CV dup segments to Neutral")

    if not result:
        return result

    # Merge adjacent Neutral segments
    merged = [result[0]]
    for seg in result[1:]:
        prev = merged[-1]
        if (seg['state'] == 'Neutral' and prev['state'] == 'Neutral' and
                seg['chrom'] == prev['chrom']):
            total_bp = (prev['end'] - prev['start']) + (seg['end'] - seg['start'])
            pw = (prev['end'] - prev['start']) / total_bp
            cw = (seg['end'] - seg['start']) / total_bp
            merged[-1] = {
                'chrom':       prev['chrom'],
                'start':       prev['start'],
                'end':         seg['end'],
                'state':       'Neutral',
                'cn_median':   prev['cn_median'] * pw + seg['cn_median'] * cw,
                'cn_mean':     prev['cn_mean'] * pw + seg['cn_mean'] * cw,
                'n_windows':   prev['n_windows'] + seg['n_windows'],
                'avg_quality': prev.get('avg_quality', 1.0) * pw + seg.get('avg_quality', 1.0) * cw,
                'cn_std':      _pooled_std(
                    prev.get('cn_std', 0.0), prev['n_windows'],
                    seg.get('cn_std', 0.0), seg['n_windows'],
                    prev['cn_median'], seg['cn_median'],
                ),
                'avg_entropy': 0.0,
                'max_entropy': 0.0,
            }
        else:
            merged.append(seg)

    if len(merged) != len(result):
        print(f"[POST] Merged adjacent Neutral after quality/CV filter: "
              f"{len(result)} -> {len(merged)}")
    return merged


def merge_nearby_dup_segments(segments, max_gap):
    """Merge same-state dup segments separated by a short Neutral gap.

    Finds pattern: DupX | Neutral(short) | DupX and merges into one DupX segment.
    Runs iteratively until no more merges are possible.
    """
    if max_gap <= 0 or len(segments) < 3:
        return segments

    dup_states = {'LowDup', 'HighDup', 'Amp', 'HighAmp'}
    total_merged = 0

    changed = True
    while changed:
        changed = False
        result = []
        i = 0
        while i < len(segments):
            if (i + 2 < len(segments) and
                    segments[i]['state'] in dup_states and
                    segments[i+1]['state'] == 'Neutral' and
                    segments[i+2]['state'] == segments[i]['state'] and
                    segments[i]['chrom'] == segments[i+2]['chrom'] and
                    (segments[i+1]['end'] - segments[i+1]['start']) <= max_gap):

                s1, s2 = segments[i], segments[i+2]
                total_bp = (s1['end'] - s1['start']) + (s2['end'] - s2['start'])
                w1 = (s1['end'] - s1['start']) / total_bp
                w2 = (s2['end'] - s2['start']) / total_bp

                result.append({
                    'chrom':       s1['chrom'],
                    'start':       s1['start'],
                    'end':         s2['end'],
                    'state':       s1['state'],
                    'cn_median':   s1['cn_median'] * w1 + s2['cn_median'] * w2,
                    'cn_mean':     s1['cn_mean'] * w1 + s2['cn_mean'] * w2,
                    'n_windows':   s1['n_windows'] + s2['n_windows'],
                    'avg_quality': s1.get('avg_quality', 1.0) * w1 + s2.get('avg_quality', 1.0) * w2,
                    'cn_std':      _pooled_std(
                        s1.get('cn_std', 0.0), s1['n_windows'],
                        s2.get('cn_std', 0.0), s2['n_windows'],
                        s1['cn_median'], s2['cn_median'],
                    ),
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

    if total_merged > 0:
        print(f"[POST] Merged {total_merged} nearby dup segment pairs (gap <= {max_gap} bp)")
    return segments


def compute_segment_mean_gc(segments, df, gc_df):
    """Compute mean GC content for each segment from per-window GC data.

    Uses binary search (O(S log W)) instead of per-segment full scan (O(S*W)).
    Mutates segments in-place, adding 'mean_gc' key.
    """
    n_before = len(df)
    df_gc = df.merge(gc_df[['chrom', 'start', 'gc_content']], on=['chrom', 'start'], how='left')
    df_gc['gc_content'] = df_gc['gc_content'].fillna(0.5)
    assert len(df_gc) == n_before, "GC merge changed row count"

    print(f"[GC-CAL] Computing mean GC per segment for {len(segments):,} segments...")

    gc_by_chrom = {}
    for chrom, gdf in df_gc.groupby('chrom', sort=False):
        gdf_s = gdf.sort_values('start')
        gc_by_chrom[chrom] = (
            gdf_s['start'].values.astype(np.int64),
            gdf_s['gc_content'].values,
        )

    for seg in segments:
        chrom, s, e = seg['chrom'], seg['start'], seg['end']
        entry = gc_by_chrom.get(chrom)
        if entry is None:
            seg['mean_gc'] = 0.5
            continue
        starts, gcs = entry
        lo = int(np.searchsorted(starts, s, side='left'))
        hi = int(np.searchsorted(starts, e, side='left'))
        if hi > lo:
            seg['mean_gc'] = float(np.mean(gcs[lo:hi]))
        else:
            idx = min(int(np.searchsorted(starts, (s + e) // 2)), len(starts) - 1)
            seg['mean_gc'] = float(gcs[idx])

    print(f"[GC-CAL] Done.")
    return segments


def apply_gc_cn_calibration(segments):
    """Post-segmentation GC-based CN calibration.

    Strategy:
      1. Use Neutral segments as ground truth reference (expected CN ≈ 1.0).
      2. Fit a polynomial: observed_cn ~ mean_gc on Neutral segments.
      3. For all segments: cn_corrected = cn_median / poly(mean_gc).

    Requires 'mean_gc' key in each segment (added by compute_segment_mean_gc).
    """
    neutral_segs = [s for s in segments
                    if s['state'] == 'Neutral' and 'mean_gc' in s
                    and 0.1 < s['mean_gc'] < 0.9]

    if len(neutral_segs) < 20:
        print(f"[GC-CAL] WARNING: only {len(neutral_segs)} neutral segments — skipping calibration")
        return segments

    gc_vals = np.array([s['mean_gc'] for s in neutral_segs])
    cn_vals = np.array([s['cn_median'] for s in neutral_segs])

    cn_ok = (cn_vals > 0.5) & (cn_vals < 2.0)
    if cn_ok.sum() < 20:
        print(f"[GC-CAL] WARNING: too few well-behaved neutral segments — skipping calibration")
        return segments

    gc_vals, cn_vals = gc_vals[cn_ok], cn_vals[cn_ok]

    poly_degree = 3
    coeffs = np.polyfit(gc_vals, cn_vals, poly_degree)
    poly = np.poly1d(coeffs)

    pred = poly(gc_vals)
    residual_std = float(np.std(cn_vals - pred))
    gc_range_pred = poly(np.array([gc_vals.min(), gc_vals.max()]))
    print(f"[GC-CAL] Polynomial degree={poly_degree}, fitted on {len(gc_vals)} neutral segs")
    print(f"[GC-CAL] GC range: [{gc_vals.min():.3f}, {gc_vals.max():.3f}]  "
          f"CN bias range: [{gc_range_pred.min():.3f}, {gc_range_pred.max():.3f}]")
    print(f"[GC-CAL] Residual std on neutral segs: {residual_std:.4f}")

    GC_BYPASS_CLASSES = frozenset({'Satellite', 'Low_complexity'})

    n_corrected = 0
    n_clamped = 0
    n_bypassed = 0
    for seg in segments:
        if 'mean_gc' not in seg:
            continue

        skip_gc = False
        rc = seg.get('repeat_class', 'None')
        mf = seg.get('masked_fraction', 0.0)
        gc = seg.get('mean_gc', 0.5)

        if rc in GC_BYPASS_CLASSES:
            skip_gc = True
        elif gc > 0.60:
            skip_gc = True
        elif mf > 0.85 and gc > 0.55:
            skip_gc = True

        if skip_gc:
            seg['gc_bias_factor'] = 1.0
            n_bypassed += 1
            continue

        gc_factor_raw = float(poly(seg['mean_gc']))
        gc_factor = max(0.6, min(1.8, gc_factor_raw))
        if gc_factor != gc_factor_raw:
            n_clamped += 1
        seg['cn_median'] = seg['cn_median'] / gc_factor
        seg['cn_mean']   = seg['cn_mean']   / gc_factor
        if seg.get('cn_std', 0.0) > 0:
            seg['cn_std'] = seg['cn_std'] / gc_factor
        seg['gc_bias_factor'] = gc_factor
        n_corrected += 1

    print(f"[GC-CAL] Applied to {n_corrected:,} segments "
          f"({n_bypassed:,} bypassed: Satellite/high-GC)")
    if n_clamped > 0:
        print(f"[GC-CAL] WARNING: {n_clamped} segments had GC correction factor "
              f"clamped to [0.6, 1.8]")
    return segments


def reclassify_by_cn_threshold(segments, lowdup_threshold=1.25, hetdel_threshold=0.75):
    """Post-GC calibration CN-based state correction.

    1. LowDup → Neutral: cn_median < lowdup_threshold (default 1.25).
    2. HetDel → Neutral: cn_median > hetdel_threshold (GC-bias artefact).
    3. HetDel → LowDup: cn_median > 1.5.
    4. Neutral → Dup: cn_median > lowdup_threshold (HMM inertia artefact).
    """
    recl_lowdup = 0
    recl_hetdel_neutral = 0
    recl_hetdel_lowdup = 0
    recl_neutral_dup = 0
    result = []
    for seg in segments:
        cn = seg.get('cn_median', 1.0)
        if seg['state'] == 'LowDup' and cn < lowdup_threshold:
            seg = dict(seg)
            seg['state'] = 'Neutral'
            recl_lowdup += 1
        elif seg['state'] == 'HetDel' and cn > 1.5:
            seg = dict(seg)
            seg['state'] = 'LowDup'
            recl_hetdel_lowdup += 1
        elif seg['state'] == 'HetDel' and cn > hetdel_threshold:
            seg = dict(seg)
            seg['state'] = 'Neutral'
            recl_hetdel_neutral += 1
        elif seg['state'] == 'Neutral' and cn > lowdup_threshold:
            seg = dict(seg)
            seg['state'] = assign_state_from_cn(cn)
            recl_neutral_dup += 1
        result.append(seg)

    if recl_lowdup > 0:
        print(f"[CN-RECL] Reclassified {recl_lowdup:,} LowDup→Neutral "
              f"(cn_median < {lowdup_threshold:.2f})")
    if recl_hetdel_neutral > 0:
        print(f"[CN-RECL] Reclassified {recl_hetdel_neutral:,} HetDel→Neutral "
              f"(cn_median > {hetdel_threshold:.2f}, GC-bias artefact)")
    if recl_hetdel_lowdup > 0:
        print(f"[CN-RECL] Reclassified {recl_hetdel_lowdup:,} HetDel→LowDup "
              f"(cn_median > 1.50)")
    if recl_neutral_dup > 0:
        print(f"[CN-RECL] Reclassified {recl_neutral_dup:,} Neutral→Dup "
              f"(cn_median > {lowdup_threshold:.2f}, HMM inertia artefact)")

    if not result:
        return result

    # Merge adjacent Neutral segments created by reclassification
    merged = [result[0]]
    for seg in result[1:]:
        prev = merged[-1]
        if (prev['state'] == 'Neutral' and seg['state'] == 'Neutral'
                and prev['chrom'] == seg['chrom']
                and seg['start'] <= prev['end'] + 1):
            total_bp = (prev['end'] - prev['start']) + (seg['end'] - seg['start'])
            pw = (prev['end'] - prev['start']) / total_bp
            cw = (seg['end'] - seg['start']) / total_bp
            merged[-1] = {
                'chrom':       prev['chrom'],
                'start':       prev['start'],
                'end':         seg['end'],
                'state':       'Neutral',
                'cn_median':   prev['cn_median'] * pw + seg['cn_median'] * cw,
                'cn_mean':     prev.get('cn_mean', prev['cn_median']) * pw + seg.get('cn_mean', seg['cn_median']) * cw,
                'n_windows':   prev['n_windows'] + seg['n_windows'],
                'avg_quality': prev.get('avg_quality', 1.0) * pw + seg.get('avg_quality', 1.0) * cw,
                'cn_std':      _pooled_std(
                    prev.get('cn_std', 0.0), prev['n_windows'],
                    seg.get('cn_std', 0.0), seg['n_windows'],
                    prev['cn_median'], seg['cn_median'],
                ),
                'avg_entropy': 0.0,
                'max_entropy': 0.0,
            }
        else:
            merged.append(seg)

    print(f"[CN-RECL] After merge: {len(result) - len(merged):,} adjacent Neutral segments merged")
    return merged


def split_high_cv_segments(segments, df, cv_threshold=0.6, min_length=3000, max_depth=3):
    """Split dup/amp segments with high within-segment CN variability.

    High CV (cn_std / cn_median > cv_threshold) signals that the segment spans
    a CN transition. Splitting produces smaller segments with lower internal
    variance and more accurate CN assignment.

    Algorithm:
      For each dup segment with CV > cv_threshold:
        1. Extract window CN values via binary search.
        2. Find split index minimising weighted within-group variance (O(n) prefix sums).
        3. Build two sub-segments; re-assign state from CN_STATE_BOUNDARIES.
        4. Apply parent's gc_bias_factor to calibrate sub-segment CN.
        5. Recurse up to max_depth times.
    """
    import time as _t
    _t0 = _t.time()

    DUP_STATES = {'LowDup', 'HighDup', 'Amp', 'HighAmp'}

    idx = {}
    for chrom, gdf in df.groupby('chrom', sort=False):
        gdf_s = gdf.sort_values('start')
        idx[chrom] = (
            gdf_s['start'].values.astype(np.int64),
            gdf_s['cn'].values.astype(np.float64),
        )

    def _find_best_split(chrom, seg_start, seg_end, gc_factor):
        if chrom not in idx:
            return None
        starts, raw_cns = idx[chrom]

        lo = int(np.searchsorted(starts, seg_start, side='left'))
        hi = int(np.searchsorted(starts, seg_end,   side='left'))
        n  = hi - lo
        if n < 4:
            return None

        win_cns = raw_cns[lo:hi]
        ps    = np.cumsum(win_cns)
        ps_sq = np.cumsum(win_cns ** 2)

        best_score = float('inf')
        best_k     = -1

        for k in range(2, n - 2):
            left_end_bp  = int(starts[lo + k])
            right_end_bp = seg_end

            if (left_end_bp  - seg_start) < min_length:
                continue
            if (right_end_bp - left_end_bp) < min_length:
                break

            sl  = float(ps[k - 1])
            sql = float(ps_sq[k - 1])
            varl = sql / k - (sl / k) ** 2

            sr  = float(ps[n - 1]) - float(ps[k - 1])
            sqr = float(ps_sq[n - 1]) - float(ps_sq[k - 1])
            nr  = n - k
            varr = sqr / nr - (sr / nr) ** 2

            score = k * varl + nr * varr
            if score < best_score:
                best_score = score
                best_k     = k

        if best_k < 0:
            return None

        split_bp = int(starts[lo + best_k])

        def _make_sub(cn_slice, s, e, parent):
            raw_med = float(np.median(cn_slice))
            raw_std = float(np.std(cn_slice)) if len(cn_slice) > 1 else 0.0
            cal_med = raw_med / gc_factor
            cal_std = raw_std / gc_factor
            sub = dict(parent)
            sub['start']      = s
            sub['end']        = e
            sub['cn_median']  = cal_med
            sub['cn_mean']    = float(np.mean(cn_slice)) / gc_factor
            sub['cn_std']     = cal_std
            sub['n_windows']  = len(cn_slice)
            sub['state']      = assign_state_from_cn(cal_med)
            sub['avg_entropy'] = 0.0
            sub['max_entropy'] = 0.0
            return sub

        left_cns  = win_cns[:best_k]
        right_cns = win_cns[best_k:]
        return (
            _make_sub(left_cns,  seg_start, split_bp, {}),
            _make_sub(right_cns, split_bp,  seg_end,  {}),
        )

    def _recursive_split(seg, depth):
        if depth >= max_depth or seg['state'] not in DUP_STATES:
            return [seg]
        cn_m = seg.get('cn_median', 1.0)
        cn_s = seg.get('cn_std',    0.0)
        cv   = cn_s / (cn_m + EPSILON)
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

    elapsed = _t.time() - _t0
    print(f"[CV-SPLIT] {n_candidate:,} high-CV dup segments → "
          f"{n_extra:,} additional segments ({elapsed:.1f}s)")
    print(f"[CV-SPLIT] Total: {len(segments):,} → {len(result):,} segments")
    return result


def compute_segment_repeat_annotation(segments, df, rm_df):
    """Attach per-segment repeat annotation (masked_fraction, repeat_class).

    For each segment, computes:
      - masked_fraction: mean of per-window masked_fraction
      - repeat_class: most common class by window count (priority: Satellite first)
    """
    CLASS_PRIORITY = ['Satellite', 'Simple_repeat', 'LINE', 'SINE', 'LTR',
                      'DNA', 'Low_complexity', 'Other', 'None']

    print(f"[RM] Computing repeat annotation for {len(segments):,} segments...")

    rm_by_chrom = {}
    for chrom, grp in rm_df.groupby('chrom'):
        grp_s = grp.sort_values('start')
        rm_by_chrom[chrom] = (
            grp_s['start'].values,
            grp_s['masked_fraction'].values,
            grp_s['repeat_class'].values,
        )

    for seg in segments:
        chrom = seg['chrom']
        if chrom not in rm_by_chrom:
            seg['masked_fraction'] = 0.0
            seg['repeat_class'] = 'None'
            continue

        starts_arr, mf_arr, cls_arr = rm_by_chrom[chrom]
        lo = int(np.searchsorted(starts_arr, seg['start']))
        hi = int(np.searchsorted(starts_arr, seg['end']))
        if lo >= hi:
            seg['masked_fraction'] = 0.0
            seg['repeat_class'] = 'None'
            continue

        seg['masked_fraction'] = float(np.mean(mf_arr[lo:hi]))

        from collections import Counter
        cls_counts = Counter(cls_arr[lo:hi])
        if cls_counts:
            max_cnt = max(cls_counts.values())
            tied = [c for c, n in cls_counts.items() if n == max_cnt]
            best = min(tied,
                       key=lambda c: CLASS_PRIORITY.index(c)
                       if c in CLASS_PRIORITY else len(CLASS_PRIORITY))
            seg['repeat_class'] = best
        else:
            seg['repeat_class'] = 'None'

    n_masked = sum(1 for s in segments if s.get('masked_fraction', 0) > 0.5)
    print(f"[RM] Done. Segments with >50% masked: {n_masked:,} / {len(segments):,}")
    return segments


# =============================================================================
# OUTPUT
# =============================================================================

def write_output(segments, output_file, extended=False):
    """Write segments to BED file.

    Base format (7 columns):
      chrom, start, end, state, cn_median, cn_mean, n_windows

    Extended format (15 columns):
      + avg_quality, min_quality, cn_std, avg_repeats,
        avg_entropy (=0 for CBS), max_entropy (=0 for CBS),
        masked_fraction, repeat_class

    Optional extra columns (appended when data present + extended=True):
      gc_bias_factor (col 16), boundary_conf (col 17)

    avg_entropy and max_entropy are always 0.0 for CBS segments.
    boundary_conf is None for first/last segments of each chromosome —
    written as empty string to avoid literal 'None' in output.
    """
    print(f"[IO] Writing {len(segments):,} segments to {output_file}...")

    has_repeat = any('repeat_class' in seg for seg in segments[:10])
    has_gc_fac = any('gc_bias_factor' in seg for seg in segments[:10])
    has_bconf  = any('boundary_conf' in seg for seg in segments[:10])

    with open(output_file, 'w') as f:
        if extended:
            header = ("#chrom\tstart\tend\tstate\tcn_median\tcn_mean\tn_windows\t"
                      "avg_quality\tmin_quality\tcn_std\tavg_repeats\t"
                      "avg_entropy\tmax_entropy")
            if has_repeat:
                header += "\tmasked_fraction\trepeat_class"
            if has_gc_fac:
                header += "\tgc_bias_factor"
            if has_bconf:
                header += "\tboundary_conf"
            f.write(header + "\n")
        else:
            f.write("#chrom\tstart\tend\tstate\tcn_median\tcn_mean\tn_windows\n")

        for seg in segments:
            line = (f"{seg['chrom']}\t{seg['start']}\t{seg['end']}\t"
                    f"{seg['state']}\t{seg['cn_median']:.4f}\t"
                    f"{seg['cn_mean']:.4f}\t{seg['n_windows']}")
            if extended:
                line += (f"\t{seg.get('avg_quality', 1.0):.4f}"
                         f"\t{seg.get('min_quality', 1.0):.4f}"
                         f"\t{seg.get('cn_std', 0.0):.4f}"
                         f"\t{seg.get('avg_repeats', 0.0):.2f}"
                         f"\t{seg.get('avg_entropy', 0.0):.4f}"
                         f"\t{seg.get('max_entropy', 0.0):.4f}")
                if has_repeat:
                    line += (f"\t{seg.get('masked_fraction', 0.0):.4f}"
                             f"\t{seg.get('repeat_class', 'None')}")
                if has_gc_fac:
                    line += f"\t{seg.get('gc_bias_factor', 1.0):.4f}"
                if has_bconf:
                    bc = seg.get('boundary_conf')
                    line += f"\t{bc:.4f}" if bc is not None else "\t"
            f.write(line + "\n")

    print(f"[IO] Done.")


def print_statistics(segments):
    """Print summary statistics."""
    print("\n" + "=" * 60)
    print("CBS SEGMENTATION STATISTICS")
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

    total_segments = len(segments)
    total_bases    = sum(state_bases.values())

    print(f"\nTotal segments: {total_segments:,}")
    print(f"Total bases: {total_bases:,}")
    print("\nState Distribution:")
    print("-" * 60)
    print(f"{'State':<12} {'Segments':>10} {'Pct':>8} {'Bases':>15} {'Mean CN':>10}")
    print("-" * 60)

    for state_name in STATE_ORDER:
        count = state_counts.get(state_name, 0)
        bases = state_bases.get(state_name, 0)
        pct   = 100.0 * count / total_segments if total_segments > 0 else 0
        mean_cn = state_cn_sum.get(state_name, 0) / bases if bases > 0 else 0
        print(f"{state_name:<12} {count:>10,} {pct:>7.1f}% {bases:>15,} {mean_cn:>10.2f}")

    print("-" * 60)
    print(f"\nCompression: {total_segments:,} segments from input windows")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="KonuSeg CBS/PELT CN Segmentation — drop-in alternative to HMM"
    )
    parser.add_argument("--input",  "-i", required=True,
                        help="Input BED file from preprocess_chm13_ont.py "
                             "(8-col: chrom start end cn mean_count log_ratio num_kmers num_filtered)")
    parser.add_argument("--output", "-o", required=True,
                        help="Output segments BED file")
    parser.add_argument("--extended", action="store_true",
                        help="Write extended 15-col output (quality, cn_std, "
                             "avg_entropy=0, boundary_conf, optional GC/repeat cols). "
                             "Default: enabled for CBS (mirrors cn-accuracy HMM mode).")
    parser.add_argument("--gc-content-bed",
                        help="GC content BED from compute_gc_content.py "
                             "(chrom, start, end, gc_content). "
                             "Enables post-segmentation GC CN calibration.")
    parser.add_argument("--repeat-bed",
                        help="Repeat-annotated window BED from compute_repeat_annotation.py "
                             "(10-col). Adds masked_fraction and repeat_class to output. "
                             "Does NOT change CBS breakpoints — metadata only.")
    parser.add_argument("--min-kmers", type=int, default=0,
                        help="Minimum num_kmers per window. Windows below threshold "
                             "are set to CN=1.0 (Neutral prior). Recommended: 30.")
    parser.add_argument("--quality-threshold", type=float, default=0.7,
                        help="Hard quality filter: dup segments with avg_quality "
                             "below this value are reclassified to Neutral. "
                             "Default: 0.7. Set to 0 to disable.")
    parser.add_argument("--cv-filter-threshold", type=float, default=0.0,
                        help="CBS-only: CV filter threshold for quality filtering. "
                             "Dup segments with cn_std/cn_median > this value are "
                             "reclassified to Neutral (high within-segment variability). "
                             "Default: 0.0 (disabled). Useful range: 0.4-0.8.")
    parser.add_argument("--cn-reclassify-threshold", type=float, default=1.25,
                        help="Post-calibration CN threshold: LowDup segments with "
                             "cn_median below this value are reclassified to Neutral. "
                             "Default: 1.25. Set to 0 to disable.")
    parser.add_argument("--cv-split-threshold", type=float, default=0.6,
                        help="CV threshold for post-merge segment splitting. "
                             "Dup segments with cn_std/cn_median above this value "
                             "are split at the variance-minimising breakpoint. "
                             "Default: 0.6. Set to 0 to disable.")
    parser.add_argument("--min-segment-length", type=int, default=3000,
                        help="Minimum segment length in bp (PELT min_size). "
                             "Default: 3000. Equivalent to HMM min_segment_lengths.")
    parser.add_argument("--penalty", type=float, default=None,
                        help="PELT penalty override. Default: None = BIC per "
                             "chromosome (2*log(n)). Smaller penalty → more segments; "
                             "larger → fewer segments.")
    # HMM-compatible arguments (accepted, semantics differ or no-op)
    parser.add_argument("--no-chrx-correction", action="store_true",
                        help="Accepted for CLI compatibility. CBS uses preprocess "
                             "normalization for chrX — no additional correction needed.")
    parser.add_argument("--mode", default=None,
                        choices=['precision', 'sensitive', 'cn-accuracy'],
                        help="Accepted for CLI compatibility. CBS always runs in "
                             "cn-accuracy equivalent mode (QW=0, no HMM emission).")
    parser.add_argument("--no-pomegranate", action="store_true",
                        help="Accepted for CLI compatibility (no-op — CBS does not use pomegranate).")
    parser.add_argument("--coarse-output",
                        help="Accepted for CLI compatibility (not implemented in CBS).")

    args = parser.parse_args()

    # Notify about HMM-only arguments
    if args.mode is not None:
        print(f"[INFO] --mode '{args.mode}' is HMM-only; CBS always runs in cn-accuracy equivalent mode")
    if args.no_pomegranate:
        print(f"[INFO] --no-pomegranate is a no-op for CBS (pomegranate not used)")
    if args.coarse_output:
        print(f"[INFO] --coarse-output is not implemented in CBS segmenter")

    print("=" * 60)
    print("KonuSeg CBS/PELT Segmentation")
    print("=" * 60)
    print(f"Input:               {args.input}")
    print(f"Output:              {args.output}")
    print(f"Penalty:             {'BIC (auto)' if args.penalty is None else args.penalty}")
    print(f"Min segment length:  {args.min_segment_length} bp")
    print(f"Quality threshold:   {args.quality_threshold}")
    print(f"CV filter threshold: {args.cv_filter_threshold}")
    print(f"CV split threshold:  {args.cv_split_threshold}")
    print(f"CN reclassify thr:   {args.cn_reclassify_threshold}")
    print(f"Min kmers:           {args.min_kmers}")
    print(f"GC calibration:      {'yes' if args.gc_content_bed else 'no'}")
    print(f"Repeat annotation:   {'yes' if args.repeat_bed else 'no'}")
    print()

    # -------------------------------------------------------------------------
    # Step 1: Load data
    # -------------------------------------------------------------------------
    df = load_data(args.input, min_kmers=args.min_kmers)

    gc_df = None
    if args.gc_content_bed:
        gc_df = load_gc_content(args.gc_content_bed)
        print()

    # -------------------------------------------------------------------------
    # Step 2: CBS/PELT segmentation
    # -------------------------------------------------------------------------
    print()
    segments = run_segmentation(df, penalty=args.penalty,
                                min_size=args.min_segment_length)

    # -------------------------------------------------------------------------
    # Step 3: Post-processing chain (mirrors HMM cn-accuracy mode)
    # -------------------------------------------------------------------------
    print()

    # 3a. Filter small segments
    min_lengths = {
        'LowDup':  args.min_segment_length,
        'HighDup': args.min_segment_length,
        'Amp':     args.min_segment_length,
    }
    segments = filter_small_segments(segments, min_lengths)

    # 3b. Quality + CV filter
    segments = filter_low_quality_segments(
        segments,
        threshold=args.quality_threshold,
        cv_filter_threshold=args.cv_filter_threshold,
    )

    # 3c. Merge nearby dup segments
    segments = merge_nearby_dup_segments(segments, max_gap=10000)

    # 3d. GC calibration (cn-accuracy mode: post-segmentation only)
    if gc_df is not None:
        print()
        segments = compute_segment_mean_gc(segments, df, gc_df)
        segments = apply_gc_cn_calibration(segments)

    # 3e. CN-based state reclassification
    if args.cn_reclassify_threshold > 0:
        print()
        segments = reclassify_by_cn_threshold(
            segments, lowdup_threshold=args.cn_reclassify_threshold)

    # 3f. CV-based segment splitting
    if args.cv_split_threshold > 0:
        print()
        segments = split_high_cv_segments(
            segments, df, cv_threshold=args.cv_split_threshold,
            min_length=args.min_segment_length)

    # 3g. Repeat annotation (metadata only)
    if args.repeat_bed:
        print()
        rm_df = load_repeat_annotation(args.repeat_bed)
        segments = compute_segment_repeat_annotation(segments, df, rm_df)

    # -------------------------------------------------------------------------
    # Step 4: Boundary confidence (CBS-specific; always computed)
    # -------------------------------------------------------------------------
    print()
    segments = compute_boundary_confidence(segments, df)

    # -------------------------------------------------------------------------
    # Step 5: Output
    # -------------------------------------------------------------------------
    print_statistics(segments)
    # CBS always writes extended output by default (mirrors cn-accuracy HMM)
    do_extended = True  # override: extended always on; --extended flag kept for compatibility
    write_output(segments, args.output, extended=do_extended)


if __name__ == '__main__':
    main()

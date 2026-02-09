#!/usr/bin/env python3
"""
KonuSeg 7-State HMM Segmentation (Log-Seg & Linear-Geno)
=========================================================

This script implements a 7-state HMM for CNV segmentation optimized for
collapsed haploid input with high dynamic range (CN 0.1 to 100+).

Strategy:
- Segmentation in Log2 space (compresses dynamic range)
- Genotyping in Linear space (preserves accurate CN values)

States (Haploid baseline, CN=1 is normal):
  0: HomDel   - Complete deletion (CN ~ 0)
  1: HetDel   - Hemizygous/Half copy (CN ~ 0.5, e.g., chrX in male)
  2: Neutral  - Normal haploid (CN ~ 1)
  3: LowDup   - Low duplication (CN ~ 2)
  4: HighDup  - High duplication (CN ~ 4)
  5: Amp      - Amplification (CN ~ 8)
  6: HighAmp  - High amplification (CN ~ 16+)

Author: KonuSeg Team
Date: February 2026

Modes:
  --mode precision  : Best F1 (QW=3.0 + strict filters) for final output
  --mode sensitive  : High-recall candidate generation for ML classifier (Pass 1)
"""

import numpy as np
import pandas as pd
import torch
import argparse
import sys
import warnings
from collections import defaultdict

warnings.filterwarnings('ignore')

# Try to import pomegranate
try:
    from pomegranate.hmm import DenseHMM
    from pomegranate.distributions import Normal
    POMEGRANATE_AVAILABLE = True
except ImportError:
    POMEGRANATE_AVAILABLE = False
    print("[WARNING] pomegranate not available, using simple Gaussian HMM fallback")


# =============================================================================
# 7-STATE MODEL CONFIGURATION (Log2 Space)
# =============================================================================

EPSILON = 0.001  # For log2 transform

# --- MODE PRESETS ---
MODE_CONFIGS = {
    'precision': {
        'states': {
            0: {"name": "HomDel",  "log2_mean": -6.0, "log2_var": 2.0,  "cn_range": "0-0.1"},
            1: {"name": "HetDel",  "log2_mean": -0.6, "log2_var": 0.4,  "cn_range": "0.5-0.85"},
            2: {"name": "Neutral", "log2_mean":  0.0, "log2_var": 0.28, "cn_range": "0.85-1.5"},
            3: {"name": "LowDup",  "log2_mean":  1.2, "log2_var": 0.3,  "cn_range": "1.6-3"},
            4: {"name": "HighDup", "log2_mean":  2.0, "log2_var": 0.5,  "cn_range": "3-6"},
            5: {"name": "Amp",     "log2_mean":  3.0, "log2_var": 0.8,  "cn_range": "6-12"},
            6: {"name": "HighAmp", "log2_mean":  5.0, "log2_var": 1.5,  "cn_range": "12+"},
        },
        'self_transition': 0.99,
        'quality_weight': 3.0,
        'max_merge_gap': 10000,
        'min_segment_lengths': {"LowDup": 20000, "HighDup": 12000, "Amp": 10000},
    },
    'sensitive': {
        'states': {
            0: {"name": "HomDel",  "log2_mean": -6.0, "log2_var": 2.0,  "cn_range": "0-0.1"},
            1: {"name": "HetDel",  "log2_mean": -0.6, "log2_var": 0.4,  "cn_range": "0.5-0.85"},
            2: {"name": "Neutral", "log2_mean":  0.0, "log2_var": 0.15, "cn_range": "0.85-1.3"},
            3: {"name": "LowDup",  "log2_mean":  0.8, "log2_var": 0.5,  "cn_range": "1.3-3"},
            4: {"name": "HighDup", "log2_mean":  2.0, "log2_var": 0.6,  "cn_range": "3-6"},
            5: {"name": "Amp",     "log2_mean":  3.0, "log2_var": 1.0,  "cn_range": "6-12"},
            6: {"name": "HighAmp", "log2_mean":  5.0, "log2_var": 1.5,  "cn_range": "12+"},
        },
        'self_transition': 0.97,
        'quality_weight': 0.0,
        'max_merge_gap': 5000,
        'min_segment_lengths': {"LowDup": 2000, "HighDup": 2000, "Amp": 2000},
    },
}

# Active configuration (set by --mode, defaults to precision)
STATES = MODE_CONFIGS['precision']['states']
SELF_TRANSITION = MODE_CONFIGS['precision']['self_transition']
QUALITY_WEIGHT = MODE_CONFIGS['precision']['quality_weight']
MAX_MERGE_GAP = MODE_CONFIGS['precision']['max_merge_gap']
MIN_SEGMENT_LENGTHS = MODE_CONFIGS['precision']['min_segment_lengths']


# =============================================================================
# DATA LOADING
# =============================================================================

def load_data(bed_file):
    """Load BED file and compute log2 values."""
    print(f"[IO] Loading {bed_file}...")

    try:
        # Try to read with different column counts
        df = pd.read_csv(bed_file, sep='\t', header=None, comment='#')

        # Determine format based on column count
        if len(df.columns) >= 8:
            # Full format: chrom, start, end, CN, mean_count, log_ratio, num_kmers, num_filtered
            df.columns = ['chrom', 'start', 'end', 'cn', 'mean_count', 'log_ratio',
                          'num_kmers', 'num_filtered'][:len(df.columns)]
        elif len(df.columns) >= 4:
            # Minimal format: chrom, start, end, CN
            df.columns = ['chrom', 'start', 'end', 'cn'] + \
                        [f'col{i}' for i in range(4, len(df.columns))]
        else:
            raise ValueError(f"Unexpected column count: {len(df.columns)}")

        print(f"[IO] Loaded {len(df):,} windows")
        print(f"[IO] Columns: {list(df.columns[:4])}")

        # Compute quality score and adjusted CN
        if 'num_kmers' in df.columns and 'num_filtered' in df.columns and QUALITY_WEIGHT > 0:
            df['quality'] = df['num_kmers'] / (df['num_kmers'] + df['num_filtered'] + 1)
            q = df['quality'].values ** QUALITY_WEIGHT
            df['cn_adjusted'] = df['cn'].values * q + 1.0 * (1.0 - q)
            df['log2_cn'] = np.log2(df['cn_adjusted'].values + EPSILON)
            print(f"[IO] Quality-weighted CN enabled (weight={QUALITY_WEIGHT})")
            print(f"[IO] Quality: median={df['quality'].median():.3f}, "
                  f"low(<0.5)={100*(df['quality']<0.5).mean():.1f}%")
        else:
            df['quality'] = 1.0
            df['log2_cn'] = np.log2(df['cn'].values + EPSILON)

        # Remove inf/nan
        before = len(df)
        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.dropna(subset=['log2_cn'])
        after = len(df)

        if before != after:
            print(f"[IO] Removed {before - after} invalid rows")

        # Statistics
        print(f"[IO] CN range: {df['cn'].min():.4f} - {df['cn'].max():.4f}")
        print(f"[IO] Log2 range: {df['log2_cn'].min():.2f} - {df['log2_cn'].max():.2f}")
        print(f"[IO] Median CN: {df['cn'].median():.4f}")

        return df

    except Exception as e:
        print(f"[ERROR] Failed to load data: {e}")
        sys.exit(1)


# =============================================================================
# HMM MODEL BUILDING
# =============================================================================

def build_transition_matrix(n_states=7, self_prob=0.95):
    """
    Build transition matrix with distance-based penalties.
    Transitions to nearby states are more likely than distant ones.
    """
    transitions = np.zeros((n_states, n_states))

    for i in range(n_states):
        transitions[i, i] = self_prob

        # Distribute remaining probability based on distance
        weights = []
        for j in range(n_states):
            if i != j:
                # Inverse square distance weighting
                distance = abs(i - j)
                weight = 1.0 / (distance ** 2)
                weights.append((j, weight))

        # Normalize weights
        total_weight = sum(w for _, w in weights)
        remaining_prob = 1.0 - self_prob

        for j, weight in weights:
            transitions[i, j] = (weight / total_weight) * remaining_prob

    return transitions


def build_7state_hmm_pomegranate(states_override=None):
    """Build 7-state HMM using pomegranate library."""
    states = states_override if states_override else STATES
    n_states = len(states)

    # Create distributions
    distributions = []
    for i in range(n_states):
        state = states[i]
        # Normal distribution with mean and variance
        dist = Normal(
            means=torch.tensor([state["log2_mean"]], dtype=torch.float32),
            covs=torch.tensor([[state["log2_var"]]], dtype=torch.float32)
        )
        distributions.append(dist)

    # Start probabilities (most likely to start in Neutral)
    starts = torch.tensor([0.01, 0.05, 0.85, 0.05, 0.02, 0.01, 0.01], dtype=torch.float32)

    # Transition matrix
    transitions = build_transition_matrix(n_states, SELF_TRANSITION)
    edges = torch.tensor(transitions, dtype=torch.float32)

    # Build model
    model = DenseHMM(
        distributions=distributions,
        edges=edges,
        starts=starts,
        ends=torch.ones(n_states, dtype=torch.float32) * 0.1,
        max_iter=1,
        verbose=False
    )

    return model


def log_gaussian_pdf(x, mean, var):
    """Log probability density of Gaussian."""
    return -0.5 * np.log(2 * np.pi * var) - 0.5 * ((x - mean) ** 2) / var


def viterbi_simple(observations, states, start_probs, trans_probs):
    """
    Simple Viterbi implementation for fallback when pomegranate is unavailable.
    """
    n_obs = len(observations)
    n_states = len(states)

    # Initialize
    V = np.zeros((n_obs, n_states))
    path = np.zeros((n_obs, n_states), dtype=int)

    # First observation
    for s in range(n_states):
        state = states[s]
        V[0, s] = np.log(start_probs[s] + 1e-10) + \
                  log_gaussian_pdf(observations[0], state["log2_mean"], state["log2_var"])

    # Forward pass
    for t in range(1, n_obs):
        for s in range(n_states):
            state = states[s]
            emission = log_gaussian_pdf(observations[t], state["log2_mean"], state["log2_var"])

            max_prob = -np.inf
            max_state = 0

            for prev_s in range(n_states):
                prob = V[t-1, prev_s] + np.log(trans_probs[prev_s, s] + 1e-10)
                if prob > max_prob:
                    max_prob = prob
                    max_state = prev_s

            V[t, s] = max_prob + emission
            path[t, s] = max_state

    # Backtrack
    best_path = np.zeros(n_obs, dtype=int)
    best_path[-1] = np.argmax(V[-1])

    for t in range(n_obs - 2, -1, -1):
        best_path[t] = path[t + 1, best_path[t + 1]]

    return best_path


# =============================================================================
# SEGMENTATION
# =============================================================================

def segment_chromosome(df_chrom, model=None, use_pomegranate=True):
    """
    Segment a single chromosome using HMM.
    Returns list of segments with state assignments.
    """
    if len(df_chrom) == 0:
        return []

    # Sort by position
    df_chrom = df_chrom.sort_values('start').reset_index(drop=True)

    # Get log2 values
    X = df_chrom['log2_cn'].values

    # Run HMM
    if use_pomegranate and model is not None:
        try:
            X_tensor = torch.tensor(X.reshape(1, -1, 1), dtype=torch.float32)
            with torch.no_grad():
                y_pred = model.predict(X_tensor)[0].numpy()
        except Exception as e:
            print(f"    [WARNING] Pomegranate failed: {e}, using fallback")
            use_pomegranate = False

    if not use_pomegranate:
        # Fallback to simple Viterbi
        start_probs = np.array([0.01, 0.05, 0.85, 0.05, 0.02, 0.01, 0.01])
        trans_probs = build_transition_matrix(7, SELF_TRANSITION)
        states_list = [STATES[i] for i in range(7)]
        y_pred = viterbi_simple(X, states_list, start_probs, trans_probs)

    # Merge consecutive windows with same state into segments
    segments = []

    if len(y_pred) == 0:
        return segments

    current_state = int(y_pred[0])
    start_idx = 0

    for i in range(1, len(y_pred)):
        state = int(y_pred[i])

        if state != current_state:
            # Segment ended - calculate statistics
            subset = df_chrom.iloc[start_idx:i]

            # LINEAR-GENO: Use MEDIAN of raw CN values (not log)
            raw_cn_median = subset['cn'].median()
            raw_cn_mean = subset['cn'].mean()

            chrom = df_chrom.iloc[start_idx]['chrom']
            seg_start = int(df_chrom.iloc[start_idx]['start'])
            seg_end = int(df_chrom.iloc[i-1]['end'])

            avg_quality = subset['quality'].mean() if 'quality' in subset.columns else 1.0
            min_quality = subset['quality'].min() if 'quality' in subset.columns else 1.0
            cn_std = subset['cn'].std() if len(subset) > 1 else 0.0
            avg_repeats = subset['num_filtered'].mean() if 'num_filtered' in subset.columns else 0.0

            segments.append({
                'chrom': chrom,
                'start': seg_start,
                'end': seg_end,
                'state': STATES[current_state]["name"],
                'cn_median': raw_cn_median,
                'cn_mean': raw_cn_mean,
                'n_windows': len(subset),
                'avg_quality': avg_quality,
                'min_quality': min_quality,
                'cn_std': cn_std,
                'avg_repeats': avg_repeats,
            })

            current_state = state
            start_idx = i

    # Last segment
    subset = df_chrom.iloc[start_idx:]
    raw_cn_median = subset['cn'].median()
    raw_cn_mean = subset['cn'].mean()
    avg_quality = subset['quality'].mean() if 'quality' in subset.columns else 1.0
    min_quality = subset['quality'].min() if 'quality' in subset.columns else 1.0
    cn_std = subset['cn'].std() if len(subset) > 1 else 0.0
    avg_repeats = subset['num_filtered'].mean() if 'num_filtered' in subset.columns else 0.0

    chrom = df_chrom.iloc[start_idx]['chrom']
    seg_start = int(df_chrom.iloc[start_idx]['start'])
    seg_end = int(df_chrom.iloc[len(df_chrom)-1]['end'])

    segments.append({
        'chrom': chrom,
        'start': seg_start,
        'end': seg_end,
        'state': STATES[current_state]["name"],
        'cn_median': raw_cn_median,
        'cn_mean': raw_cn_mean,
        'n_windows': len(subset),
        'avg_quality': avg_quality,
        'min_quality': min_quality,
        'cn_std': cn_std,
        'avg_repeats': avg_repeats,
    })

    return segments


def run_segmentation(df, use_pomegranate=True):
    """Run segmentation on all chromosomes."""
    all_segments = []

    # Build model if using pomegranate
    model = None
    if use_pomegranate and POMEGRANATE_AVAILABLE:
        print("[HMM] Building 7-State Log-Space HMM (pomegranate)...")
        try:
            model = build_7state_hmm_pomegranate()
        except Exception as e:
            print(f"[WARNING] Failed to build pomegranate model: {e}")
            model = None
            use_pomegranate = False

    if not use_pomegranate or model is None:
        print("[HMM] Using simple Viterbi fallback...")

    # Build chrX-specific model (male haploid: expected CN ~0.7)
    chrx_model = None
    if use_pomegranate and POMEGRANATE_AVAILABLE:
        import copy
        chrx_states = copy.deepcopy(STATES)
        chrx_states[2]['log2_mean'] = -0.5  # Neutral → CN~0.7 for male chrX
        try:
            chrx_model = build_7state_hmm_pomegranate(states_override=chrx_states)
        except Exception:
            chrx_model = None

    # Process each chromosome
    chroms = df['chrom'].unique()
    print(f"[HMM] Processing {len(chroms)} chromosomes...")

    for chrom in chroms:
        df_chrom = df[df['chrom'] == chrom].copy()

        # Use chrX-specific model for CM039033.1
        if chrom == 'CM039033.1' and chrx_model is not None:
            segments = segment_chromosome(df_chrom, chrx_model, True)
        else:
            segments = segment_chromosome(df_chrom, model, use_pomegranate and model is not None)
        all_segments.extend(segments)

        print(f"  {chrom}: {len(df_chrom):,} windows -> {len(segments):,} segments")

    return all_segments


# =============================================================================
# POST-PROCESSING
# =============================================================================

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
            # Check for pattern: DupX | Neutral(<=max_gap) | DupX
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
                    'chrom': s1['chrom'],
                    'start': s1['start'],
                    'end': s2['end'],
                    'state': s1['state'],
                    'cn_median': s1['cn_median'] * w1 + s2['cn_median'] * w2,
                    'cn_mean': s1['cn_mean'] * w1 + s2['cn_mean'] * w2,
                    'n_windows': s1['n_windows'] + s2['n_windows'],
                    'avg_quality': s1.get('avg_quality', 1.0) * w1 + s2.get('avg_quality', 1.0) * w2,
                    'cn_std': 0.0,
                })
                total_merged += 1
                changed = True
                i += 3  # Skip all three segments
            else:
                result.append(segments[i])
                i += 1

        segments = result

    if total_merged > 0:
        print(f"[POST] Merged {total_merged} nearby dup segment pairs (gap <= {max_gap} bp)")
    return segments


def filter_small_segments(segments, min_lengths):
    """Reclassify small dup/del segments as Neutral.

    Segments shorter than the state-specific threshold are converted to Neutral.
    Adjacent Neutral segments are then merged.
    """
    if not min_lengths:
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
            prev_weight = (prev['end'] - prev['start']) / total_bp
            curr_weight = (seg['end'] - seg['start']) / total_bp

            merged[-1] = {
                'chrom': prev['chrom'],
                'start': prev['start'],
                'end': seg['end'],
                'state': 'Neutral',
                'cn_median': prev['cn_median'] * prev_weight + seg['cn_median'] * curr_weight,
                'cn_mean': prev['cn_mean'] * prev_weight + seg['cn_mean'] * curr_weight,
                'n_windows': prev['n_windows'] + seg['n_windows'],
            }
        else:
            merged.append(seg)

    if len(merged) != len(result):
        print(f"[POST] Merged adjacent Neutral: {len(result)} -> {len(merged)}")

    return merged


def filter_segments_adaptive(segments):
    """Quality-scaled segment filtering (A2 + A4 combined).

    A2: Per-state base thresholds (from Round 13) scaled by segment quality.
        high quality (>0.95) → 0.5x threshold, medium → 1.0x, low (<0.80) → 2.0x
    A4: CN consistency check — very noisy dup segments (CV > threshold) reclassified.
    """
    dup_states = {'LowDup', 'HighDup', 'Amp', 'HighAmp'}
    reclassified_size = 0
    reclassified_cv = 0
    result = []

    for seg in segments:
        state = seg['state']

        if state not in dup_states:
            result.append(seg)
            continue

        length = seg['end'] - seg['start']
        quality = seg.get('avg_quality', 1.0)
        cn_std = seg.get('cn_std', 0.0)
        cn_median = seg.get('cn_median', 1.0)

        should_reclassify = False

        # A2: Quality-scaled per-state threshold
        base_len = MIN_SEGMENT_LENGTHS.get(state, 0)
        if base_len > 0:
            if quality > 0.95:
                mult = QUALITY_MULTIPLIERS['high']
            elif quality > 0.80:
                mult = QUALITY_MULTIPLIERS['medium']
            else:
                mult = QUALITY_MULTIPLIERS['low']

            min_len = int(base_len * mult)
            if length < min_len:
                should_reclassify = True
                reclassified_size += 1

        # A4: CN consistency — high CV means unreliable CN estimate
        if not should_reclassify and cn_median > 0.1:
            cv = cn_std / (cn_median + 0.001)
            if cv > CN_CV_THRESHOLD:
                should_reclassify = True
                reclassified_cv += 1

        if should_reclassify:
            seg = dict(seg)
            seg['state'] = 'Neutral'

        result.append(seg)

    total = reclassified_size + reclassified_cv
    if total > 0:
        print(f"[POST] Adaptive filter: reclassified {total} segments to Neutral "
              f"(size={reclassified_size}, cv={reclassified_cv})")

    # Merge adjacent Neutral segments
    merged = [result[0]]
    for seg in result[1:]:
        prev = merged[-1]
        if (seg['state'] == 'Neutral' and prev['state'] == 'Neutral' and
                seg['chrom'] == prev['chrom']):
            total_bp = (prev['end'] - prev['start']) + (seg['end'] - seg['start'])
            prev_weight = (prev['end'] - prev['start']) / total_bp
            curr_weight = (seg['end'] - seg['start']) / total_bp

            merged[-1] = {
                'chrom': prev['chrom'],
                'start': prev['start'],
                'end': seg['end'],
                'state': 'Neutral',
                'cn_median': prev['cn_median'] * prev_weight + seg['cn_median'] * curr_weight,
                'cn_mean': prev['cn_mean'] * prev_weight + seg['cn_mean'] * curr_weight,
                'n_windows': prev['n_windows'] + seg['n_windows'],
                'avg_quality': prev.get('avg_quality', 1.0) * prev_weight + seg.get('avg_quality', 1.0) * curr_weight,
                'cn_std': 0.0,
            }
        else:
            merged.append(seg)

    if len(merged) != len(result):
        print(f"[POST] Merged adjacent Neutral: {len(result)} -> {len(merged)}")

    return merged


# =============================================================================
# OUTPUT AND STATISTICS
# =============================================================================

def write_output(segments, output_file, extended=False):
    """Write segments to BED file.

    If extended=True, include avg_quality and cn_std columns (for ML classifier input).
    """
    print(f"[IO] Writing {len(segments):,} segments to {output_file}...")

    with open(output_file, 'w') as f:
        if extended:
            f.write("#chrom\tstart\tend\tstate\tcn_median\tcn_mean\tn_windows\t"
                    "avg_quality\tmin_quality\tcn_std\tavg_repeats\n")
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
                         f"\t{seg.get('avg_repeats', 0.0):.2f}")
            f.write(line + "\n")


def print_statistics(segments):
    """Print summary statistics."""
    print("\n" + "=" * 60)
    print("SEGMENTATION STATISTICS")
    print("=" * 60)

    # State distribution
    state_counts = defaultdict(int)
    state_bases = defaultdict(int)
    state_cn_sum = defaultdict(float)

    for seg in segments:
        state = seg['state']
        length = seg['end'] - seg['start']

        state_counts[state] += 1
        state_bases[state] += length
        state_cn_sum[state] += seg['cn_median'] * length

    total_segments = len(segments)
    total_bases = sum(state_bases.values())

    print(f"\nTotal segments: {total_segments:,}")
    print(f"Total bases: {total_bases:,}")

    print("\nState Distribution:")
    print("-" * 60)
    print(f"{'State':<12} {'Segments':>10} {'Pct':>8} {'Bases':>15} {'Mean CN':>10}")
    print("-" * 60)

    for state_id in range(7):
        state_name = STATES[state_id]["name"]
        count = state_counts.get(state_name, 0)
        bases = state_bases.get(state_name, 0)

        pct = 100.0 * count / total_segments if total_segments > 0 else 0
        mean_cn = state_cn_sum.get(state_name, 0) / bases if bases > 0 else 0

        print(f"{state_name:<12} {count:>10,} {pct:>7.1f}% {bases:>15,} {mean_cn:>10.2f}")

    print("-" * 60)

    # Compression ratio
    print(f"\nCompression: {total_segments:,} segments from input windows")


# =============================================================================
# MAIN
# =============================================================================

def apply_mode(mode):
    """Apply mode-specific configuration to global variables."""
    global STATES, SELF_TRANSITION, QUALITY_WEIGHT, MAX_MERGE_GAP, MIN_SEGMENT_LENGTHS

    if mode not in MODE_CONFIGS:
        print(f"[ERROR] Unknown mode: {mode}. Use 'precision' or 'sensitive'.")
        sys.exit(1)

    cfg = MODE_CONFIGS[mode]
    STATES = cfg['states']
    SELF_TRANSITION = cfg['self_transition']
    QUALITY_WEIGHT = cfg['quality_weight']
    MAX_MERGE_GAP = cfg['max_merge_gap']
    MIN_SEGMENT_LENGTHS = cfg['min_segment_lengths']


def main():
    parser = argparse.ArgumentParser(
        description="7-State HMM CNV Segmentation (Log-Seg & Linear-Geno)"
    )
    parser.add_argument("--input", "-i", required=True,
                        help="Input BED file from Phase 1 preprocessing")
    parser.add_argument("--output", "-o", required=True,
                        help="Output segments BED file")
    parser.add_argument("--mode", "-m", choices=['precision', 'sensitive'],
                        default='precision',
                        help="precision: best F1 (QW=3.0, strict filters). "
                             "sensitive: high-recall candidates for ML (QW=0, liberal).")
    parser.add_argument("--no-pomegranate", action="store_true",
                        help="Force use of simple Viterbi (skip pomegranate)")

    args = parser.parse_args()

    # Apply mode configuration
    apply_mode(args.mode)

    mode_label = "HIGH-RECALL CANDIDATE GENERATION" if args.mode == 'sensitive' else "PRECISION (Best F1)"
    print("=" * 60)
    print(f"KonuSeg 7-State HMM Segmentation — {mode_label}")
    print("=" * 60)
    print(f"Input:  {args.input}")
    print(f"Output: {args.output}")
    print(f"Mode:   {args.mode}")
    print()
    print("HMM Parameters:")
    print(f"  SELF_TRANSITION: {SELF_TRANSITION}")
    print(f"  LowDup mean: {STATES[3]['log2_mean']} (CN~{2**STATES[3]['log2_mean']:.2f})")
    print(f"  LowDup var: {STATES[3]['log2_var']}")
    print(f"  Neutral var: {STATES[2]['log2_var']}")
    print(f"  QUALITY_WEIGHT: {QUALITY_WEIGHT}")
    print(f"  MIN_SEGMENT_LENGTHS: {MIN_SEGMENT_LENGTHS}")
    print(f"  MAX_MERGE_GAP: {MAX_MERGE_GAP}")
    print()

    # Load data
    df = load_data(args.input)

    # Run segmentation
    use_pom = POMEGRANATE_AVAILABLE and not args.no_pomegranate
    segments = run_segmentation(df, use_pomegranate=use_pom)

    # Post-processing
    if MAX_MERGE_GAP > 0 or MIN_SEGMENT_LENGTHS:
        print(f"\n[POST] Post-processing: merge_gap={MAX_MERGE_GAP}, min_lengths={MIN_SEGMENT_LENGTHS}")
        segments = filter_small_segments(segments, MIN_SEGMENT_LENGTHS)
        segments = merge_nearby_dup_segments(segments, MAX_MERGE_GAP)

    # Print statistics
    print_statistics(segments)

    # Write output — sensitive mode writes extended features for ML
    write_output(segments, args.output, extended=(args.mode == 'sensitive'))

    print("\n" + "=" * 60)
    print("Segmentation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()

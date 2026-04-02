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
import copy

from chrom_utils import is_chrx

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

EPSILON = 0.001          # For log2 transform
CN_FLOOR_FACTOR = 0.25  # Minimum fraction of raw CN retained as quality-weight anchor

# CN midpoint boundaries for state re-assignment after splitting/reclassification.
# Each tuple: (upper_cn_limit, state_name).
# Used only in cn-accuracy mode (reclassify_by_cn_threshold, split_high_cv_segments).
# Boundaries derived from midpoints between adjacent state log2 means:
#   Neutral=1.0, LowDup=2.3, HighDup=4.0, Amp=8.0, MedAmp=12, HighAmp=22, ExtremeAmp=64+
CN_STATE_BOUNDARIES = [
    (1.20,         'Neutral'),       # tight: true single-copy (CN < 1.2)
    (1.87,         'NoisyNeutral'),  # repeat-elevated neutral (CN 1.2-1.87, LINE/SINE k-mer sharing)
    (3.00,         'LowDup'),
    (6.00,         'HighDup'),
    (12.00,        'Amp'),
    (22.00,        'MedAmp'),
    (50.00,        'HighAmp'),
    (float('inf'), 'ExtremeAmp'),
]

# --- MODE PRESETS ---
MODE_CONFIGS = {
    'precision': {
        'states': {
            0: {"name": "HomDel",  "log2_mean": -6.0, "log2_var": 2.0,  "cn_range": "0-0.1"},
            1: {"name": "HetDel",  "log2_mean": -0.6, "log2_var": 0.4,  "cn_range": "0.5-0.85"},
            2: {"name": "Neutral", "log2_mean":  0.0, "log2_var": 0.28, "cn_range": "0.85-1.5"},
            3: {"name": "LowDup",  "log2_mean":  1.2, "log2_var": 0.4,  "cn_range": "1.6-3"},
            4: {"name": "HighDup", "log2_mean":  2.0, "log2_var": 0.6,  "cn_range": "3-6"},
            5: {"name": "Amp",     "log2_mean":  3.0, "log2_var": 0.9,  "cn_range": "6-12"},
            6: {"name": "HighAmp", "log2_mean":  5.0, "log2_var": 1.6,  "cn_range": "12+"},
        },
        'self_transition': 0.995,
        'quality_weight': 3.0,
        'max_merge_gap': 10000,
        'min_segment_lengths': {"LowDup": 3000, "HighDup": 3000, "Amp": 3000},
    },
    'sensitive': {
        'states': {
            0: {"name": "HomDel",  "log2_mean": -6.0, "log2_var": 2.0,  "cn_range": "0-0.1"},
            1: {"name": "HetDel",  "log2_mean": -0.6, "log2_var": 0.4,  "cn_range": "0.5-0.85"},
            2: {"name": "Neutral", "log2_mean":  0.0, "log2_var": 0.15, "cn_range": "0.85-1.3"},
            3: {"name": "LowDup",  "log2_mean":  0.8, "log2_var": 0.6,  "cn_range": "1.3-3"},
            4: {"name": "HighDup", "log2_mean":  2.0, "log2_var": 0.7,  "cn_range": "3-6"},
            5: {"name": "Amp",     "log2_mean":  3.0, "log2_var": 1.1,  "cn_range": "6-12"},
            6: {"name": "HighAmp", "log2_mean":  5.0, "log2_var": 1.6,  "cn_range": "12+"},
        },
        'self_transition': 0.99,
        'quality_weight': 1.0,      # 2.0 suppressed too many true SDs (-19% TP); 1.0 is a better tradeoff
        'max_merge_gap': 5000,      # was 3000 — reduce 500bp-window fragmentation
        'min_segment_lengths': {"LowDup": 1500, "HighDup": 1500, "Amp": 1500},
    },
    # CN-accuracy mode: optimized for accurate copy-number values (CN=0 to CN=150+).
    # Key design decisions:
    #   - States 0,1 repurposed from HomDel/HetDel (unused on haploid CHM13) to
    #     MedAmp/ExtremeAmp — fills the CN=9-31 gap and narrows HighAmp variance.
    #   - QW=0.0: raw CN used directly, avoids 20-37% underestimation at high CN.
    #   - self_transition=0.999: expected segment ~500kb (balance between stability
    #     and capturing 200kb GT loci like AMY1/LPA; was 0.9995 = 1Mb expectation).
    #   - max_merge_gap=10000: merge nearby dup segments up to 10kb apart.
    #   - quality_threshold=0.7: hard quality filter (not CN weighting).
    #   - All dup states in min_segment_lengths: prevents 500bp HighAmp spikes from
    #     polluting length-weighted CN averages in GT loci (was missing HighAmp).
    'cn-accuracy': {
        'states': {
            0: {"name": "MedAmp",       "log2_mean":  3.58, "log2_var": 0.25, "cn_range": "9-22"},
            1: {"name": "ExtremeAmp",   "log2_mean":  6.0,  "log2_var": 0.15, "cn_range": "50+"},
            2: {"name": "Neutral",      "log2_mean":  0.0,  "log2_var": 0.025,"cn_range": "0.85-1.2"},
            3: {"name": "LowDup",       "log2_mean":  1.2,  "log2_var": 0.08, "cn_range": "1.87-3"},
            4: {"name": "HighDup",      "log2_mean":  2.0,  "log2_var": 0.09, "cn_range": "3-6"},
            5: {"name": "Amp",          "log2_mean":  3.0,  "log2_var": 0.09, "cn_range": "6-12"},
            6: {"name": "HighAmp",      "log2_mean":  4.5,  "log2_var": 0.20, "cn_range": "20-50"},
            7: {"name": "NoisyNeutral", "log2_mean":  0.5,  "log2_var": 0.06, "cn_range": "1.2-1.87"},
        },
        'self_transition': 0.999,
        'quality_weight': 0.0,
        'max_merge_gap': 10000,
        'min_segment_lengths': {
            "LowDup": 3000, "HighDup": 3000, "Amp": 3000,
            "MedAmp": 3000, "HighAmp": 3000, "ExtremeAmp": 3000,
        },
        'quality_threshold': 0.7,
    },
}

# Active configuration (set by --mode, defaults to precision)
STATES = MODE_CONFIGS['precision']['states']
SELF_TRANSITION = MODE_CONFIGS['precision']['self_transition']
QUALITY_WEIGHT = MODE_CONFIGS['precision']['quality_weight']
MAX_MERGE_GAP = MODE_CONFIGS['precision']['max_merge_gap']
MIN_SEGMENT_LENGTHS = MODE_CONFIGS['precision']['min_segment_lengths']
QUALITY_THRESHOLD = 0.0  # 0.0 = disabled; >0 = hard quality filter post-segmentation
# Start probabilities for HMM: per-state prior for first observation.
# Precision/sensitive: [HomDel, HetDel, Neutral, LowDup, HighDup, Amp, HighAmp]
# cn-accuracy:         [MedAmp, ExtremeAmp, Neutral, LowDup, HighDup, Amp, HighAmp, NoisyNeutral]
START_PROBS = np.array([0.01, 0.05, 0.85, 0.05, 0.02, 0.01, 0.01])

# Adaptive filtering constants (used by filter_segments_adaptive)
QUALITY_MULTIPLIERS = {'high': 0.5, 'medium': 1.0, 'low': 2.0}
CN_CV_THRESHOLD = 1.5


# =============================================================================
# DATA LOADING
# =============================================================================

def load_data(bed_file, min_kmers=0):
    """Load BED file and compute log2 values.

    min_kmers: windows with num_kmers < min_kmers are set to CN=1.0 (Neutral prior)
    to suppress false HetDel calls in low-coverage regions (e.g. complex/GC-rich).
    Only applied when num_kmers column is present (8-column format).
    """
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

        # Low-coverage filter: windows with very few contributing k-mer windows
        # have unreliable CN estimates → set to Neutral (CN=1.0) to avoid false HetDel
        if min_kmers > 0 and 'num_kmers' in df.columns:
            low_cov = df['num_kmers'] < min_kmers
            n_low = int(low_cov.sum())
            if n_low > 0:
                df.loc[low_cov, 'cn'] = 1.0
                print(f"[IO] Low-coverage filter (num_kmers < {min_kmers}): "
                      f"{n_low:,} windows ({100*n_low/len(df):.1f}%) set to CN=1.0 (Neutral)")

        # Compute quality score and adjusted CN
        if 'num_kmers' in df.columns and 'num_filtered' in df.columns and QUALITY_WEIGHT > 0:
            df['quality'] = df['num_kmers'] / (df['num_kmers'] + df['num_filtered'] + 1)
            q = df['quality'].values ** QUALITY_WEIGHT
            # Dynamic floor: for amplified regions don't pull all the way to CN=1.
            # cn_floor = max(1.0, raw_cn * CN_FLOOR_FACTOR) so that normal regions
            # (CN~1) are unaffected while high-CN windows retain more signal.
            cn_floor = np.maximum(1.0, df['cn'].values * CN_FLOOR_FACTOR)
            df['cn_adjusted'] = df['cn'].values * q + cn_floor * (1.0 - q)
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

def build_transition_matrix(n_states=7, self_prob=0.95, states=None):
    """
    Build transition matrix with distance-based penalties.
    Transitions to nearby states are more likely than distant ones.

    states: optional dict {0: {'log2_mean': ...}, ...} — when provided, distance is
            computed from log2_mean separation (B4 fix: CN-proportional distances).
            When None, falls back to state-index distance.
    """
    transitions = np.zeros((n_states, n_states))

    for i in range(n_states):
        transitions[i, i] = self_prob

        # Distribute remaining probability based on distance
        weights = []
        for j in range(n_states):
            if i != j:
                # B4 fix: use log2_mean distance when states dict is provided.
                # This makes Amp↔HighAmp (log2 dist=1.5) more probable than
                # MedAmp↔ExtremeAmp (log2 dist=2.42), matching biology.
                if states is not None:
                    distance = abs(states[i]['log2_mean'] - states[j]['log2_mean'])
                    distance = max(distance, 0.1)  # prevent div-by-zero for identical means
                else:
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

    # Start probabilities (mode-aware; most likely to start in Neutral)
    starts = torch.tensor(START_PROBS.tolist(), dtype=torch.float32)

    # Transition matrix — pass states for log2_mean-based distance (B4 fix)
    transitions = build_transition_matrix(n_states, SELF_TRANSITION, states=states)
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


def compute_entropy_fallback(observations, states_list, start_probs, trans_probs):
    """Per-observation HMM posterior entropy via vectorized log-space forward-backward.

    Returns entropy in bits per window: 0 = certain state, log2(7)≈2.81 = uniform.
    High-entropy windows sit on state boundaries or carry noisy signal → likely FP.
    Used as fallback when pomegranate's predict_proba is unavailable.
    """
    n_obs, n_states = len(observations), len(states_list)
    if n_obs == 0:
        return np.zeros(0)

    # Log emission matrix: shape (n_obs, n_states)
    log_emission = np.column_stack([
        log_gaussian_pdf(observations, s["log2_mean"], s["log2_var"])
        for s in states_list
    ])

    log_trans = np.log(trans_probs + 1e-300)   # (n_states, n_states)
    log_start = np.log(start_probs + 1e-300)   # (n_states,)

    # ---- Forward pass (log scale) ----
    log_alpha = np.empty((n_obs, n_states), dtype=np.float64)
    log_alpha[0] = log_start + log_emission[0]
    for t in range(1, n_obs):
        tmp = log_alpha[t - 1, :, np.newaxis] + log_trans   # (n_states, n_states)
        m = tmp.max(axis=0)
        log_alpha[t] = np.log(np.sum(np.exp(tmp - m), axis=0)) + m + log_emission[t]

    # ---- Backward pass (log scale) ----
    log_beta = np.zeros((n_obs, n_states), dtype=np.float64)
    for t in range(n_obs - 2, -1, -1):
        tmp = log_trans + log_emission[t + 1] + log_beta[t + 1]   # (n_states, n_states)
        m = tmp.max(axis=1, keepdims=True)
        log_beta[t] = np.log(np.sum(np.exp(tmp - m), axis=1)) + m[:, 0]

    # ---- Posterior probabilities (normalize per time-step) ----
    log_post = log_alpha + log_beta
    log_Z = log_post.max(axis=1, keepdims=True)
    post = np.exp(log_post - log_Z)
    post /= post.sum(axis=1, keepdims=True)
    post = np.clip(post, 1e-10, 1.0)

    # ---- Shannon entropy per observation (bits) ----
    return -(post * np.log2(post)).sum(axis=1)


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
    X_tensor = None
    if use_pomegranate and model is not None:
        try:
            X_tensor = torch.tensor(X.reshape(1, -1, 1), dtype=torch.float32)
            with torch.no_grad():
                y_pred = model.predict(X_tensor)[0].numpy()
        except Exception as e:
            print(f"    [WARNING] Pomegranate failed: {e}, using fallback")
            use_pomegranate = False
            X_tensor = None

    if not use_pomegranate:
        # Fallback to simple Viterbi
        start_probs = START_PROBS
        states_list = [STATES[i] for i in range(len(STATES))]
        trans_probs = build_transition_matrix(len(STATES), SELF_TRANSITION, states=STATES)
        y_pred = viterbi_simple(X, states_list, start_probs, trans_probs)

    # ---- Rec 3: per-window HMM posterior entropy ----
    # High entropy = model was uncertain → boundary noise or low-signal window → FP indicator
    entropy_per_window = None
    if use_pomegranate and model is not None and X_tensor is not None:
        try:
            with torch.no_grad():
                posteriors = model.predict_proba(X_tensor)[0].numpy()  # (n_obs, n_states)
            entropy_per_window = -(posteriors * np.log2(posteriors + 1e-10)).sum(axis=1)
        except Exception:
            pass  # Fall through to forward-backward below
    if entropy_per_window is None:
        try:
            sp = START_PROBS
            sl = [STATES[i] for i in range(len(STATES))]
            tp = build_transition_matrix(len(STATES), SELF_TRANSITION, states=STATES)
            entropy_per_window = compute_entropy_fallback(X, sl, sp, tp)
        except Exception:
            entropy_per_window = np.zeros(len(X))

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

            seg_ent = entropy_per_window[start_idx:i] if entropy_per_window is not None else None
            avg_entropy = float(seg_ent.mean()) if seg_ent is not None and len(seg_ent) > 0 else 0.0
            max_entropy = float(seg_ent.max()) if seg_ent is not None and len(seg_ent) > 0 else 0.0

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
                'avg_entropy': avg_entropy,
                'max_entropy': max_entropy,
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

    seg_ent = entropy_per_window[start_idx:] if entropy_per_window is not None else None
    avg_entropy = float(seg_ent.mean()) if seg_ent is not None and len(seg_ent) > 0 else 0.0
    max_entropy = float(seg_ent.max()) if seg_ent is not None and len(seg_ent) > 0 else 0.0

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
        'avg_entropy': avg_entropy,
        'max_entropy': max_entropy,
    })

    return segments


def aggregate_coarse_windows(df, coarse_size=5000):
    """Aggregate fine windows into coarser bins for multi-resolution HMM (Rec 1).

    Groups consecutive fine-scale windows into coarse_size-bp bins, taking:
      - cn = median of fine-window CNs
      - num_kmers/num_filtered = sum (for accurate quality re-computation)
    Then re-applies the same quality-weighting formula as load_data() so the
    returned DataFrame is directly usable by run_segmentation().
    """
    df = df.copy()
    df['_bin'] = (df['start'] // coarse_size) * coarse_size

    records = []
    for (chrom, bin_start), grp in df.groupby(['chrom', '_bin'], sort=False):
        cn_med = float(grp['cn'].median())
        nk = int(grp['num_kmers'].sum()) if 'num_kmers' in grp.columns else len(grp)
        nf = int(grp['num_filtered'].sum()) if 'num_filtered' in grp.columns else 0
        mc = float(grp['mean_count'].mean()) if 'mean_count' in grp.columns else cn_med
        records.append({
            'chrom': chrom,
            'start': int(bin_start),
            'end': int(bin_start) + coarse_size,
            'cn': cn_med,
            'mean_count': mc,
            'log_ratio': float(np.log2(max(cn_med, 0.01))),
            'num_kmers': nk,
            'num_filtered': nf,
        })

    coarse_df = pd.DataFrame(records)

    # Re-apply quality-weighting (mirrors load_data logic)
    if QUALITY_WEIGHT > 0 and 'num_kmers' in coarse_df.columns:
        coarse_df['quality'] = (coarse_df['num_kmers'] /
                                (coarse_df['num_kmers'] + coarse_df['num_filtered'] + 1))
        q = coarse_df['quality'].values ** QUALITY_WEIGHT
        cn_floor = np.maximum(1.0, coarse_df['cn'].values * CN_FLOOR_FACTOR)
        coarse_df['cn_adjusted'] = coarse_df['cn'].values * q + cn_floor * (1.0 - q)
        coarse_df['log2_cn'] = np.log2(coarse_df['cn_adjusted'].values + EPSILON)
    else:
        coarse_df['quality'] = 1.0
        coarse_df['log2_cn'] = np.log2(coarse_df['cn'].values + EPSILON)

    coarse_df = coarse_df.replace([np.inf, -np.inf], np.nan).dropna(subset=['log2_cn'])
    print(f"[COARSE] {len(df):,} fine windows → {len(coarse_df):,} coarse windows ({coarse_size}bp bins)")
    return coarse_df


def run_segmentation(df, use_pomegranate=True, chrx_correction=True):
    """Run segmentation on all chromosomes."""
    all_segments = []

    # Build model if using pomegranate
    model = None
    if use_pomegranate and POMEGRANATE_AVAILABLE:
        print(f"[HMM] Building {len(STATES)}-State Log-Space HMM (pomegranate)...")
        try:
            model = build_7state_hmm_pomegranate()
        except Exception as e:
            print(f"[WARNING] Failed to build pomegranate model: {e}")
            model = None
            use_pomegranate = False

    if not use_pomegranate or model is None:
        print("[HMM] Using simple Viterbi fallback...")

    # Build chrX-specific model for male samples (HG002: hemizygous X vs female reference).
    # Skip when --no-chrx-correction is set (female samples, e.g. CHM13, where chrX is diploid).
    chrx_model = None
    if chrx_correction and use_pomegranate and POMEGRANATE_AVAILABLE:
        chrx_states = copy.deepcopy(STATES)
        # HetDel: hemizygous deletion (0 functional copies) → CN~0.25 relative to diploid
        chrx_states[1]['log2_mean'] = -2.0
        # Neutral: 1 copy (hemizygous) → CN~0.5 relative to diploid reference
        chrx_states[2]['log2_mean'] = -1.0
        # LowDup: 2 copies on chrX → CN~1.0 (same as autosomal neutral depth)
        chrx_states[3]['log2_mean'] =  0.0
        # HighDup: 4 copies → CN~2.0 relative to diploid
        chrx_states[4]['log2_mean'] =  1.0
        # Amp: 8 copies → CN~4.0 relative to diploid
        chrx_states[5]['log2_mean'] =  2.0
        # HighAmp: 16+ copies → CN~8.0+ → log2(16 * 0.5) = log2(8) = 3.0 → use 4.0 for headroom
        chrx_states[6]['log2_mean'] =  4.0
        # NoisyNeutral: CN~1.41 autosomal → CN~0.71 on chrX (hemizygous)
        if 7 in chrx_states:
            chrx_states[7]['log2_mean'] = -0.5
        try:
            chrx_model = build_7state_hmm_pomegranate(states_override=chrx_states)
        except Exception:
            chrx_model = None

    # Process each chromosome
    chroms = df['chrom'].unique()
    print(f"[HMM] Processing {len(chroms)} chromosomes...")

    import time as _time
    _t0_total = _time.time()

    for i, chrom in enumerate(chroms):
        df_chrom = df[df['chrom'] == chrom].copy()
        _t0 = _time.time()
        print(f"  [{i+1}/{len(chroms)}] {chrom}: {len(df_chrom):,} windows ...", flush=True)

        # Use chrX-specific model (detects chrX regardless of naming convention)
        if is_chrx(chrom) and chrx_model is not None:
            segments = segment_chromosome(df_chrom, chrx_model, True)
        else:
            segments = segment_chromosome(df_chrom, model, use_pomegranate and model is not None)
        all_segments.extend(segments)

        _elapsed = _time.time() - _t0
        _total   = _time.time() - _t0_total
        print(f"  [{i+1}/{len(chroms)}] {chrom}: {len(segments):,} segments "
              f"({_elapsed:.0f}s | total {_total:.0f}s)", flush=True)

    return all_segments


# =============================================================================
# POST-PROCESSING
# =============================================================================

def _pooled_std(s1, n1, s2, n2, m1, m2):
    """Combined standard deviation when merging two segment groups.

    Uses the exact parallel-groups formula:
      pooled_var = (within_SS + between_SS) / (n1 + n2 - 1)
      within_SS  = (n1-1)*s1² + (n2-1)*s2²
      between_SS = n1*n2 / (n1+n2) * (m1 - m2)²

    This preserves variance even when two segments have the same internal std
    but different means (the between-group term captures that divergence).
    Critical for cv-split pre-filter: merged segments with different CN values
    must have non-zero cn_std so split_high_cv_segments() can detect them.
    """
    n1, n2 = max(int(n1), 1), max(int(n2), 1)
    within  = max(n1 - 1, 0) * s1 ** 2 + max(n2 - 1, 0) * s2 ** 2
    between = (n1 * n2) / (n1 + n2) * (m1 - m2) ** 2
    return float(np.sqrt(max(0.0, (within + between) / (n1 + n2 - 1))))


def merge_nearby_dup_segments(segments, max_gap):
    """Merge same-state dup segments separated by a short Neutral gap.

    Finds pattern: DupX | Neutral(short) | DupX and merges into one DupX segment.
    Runs iteratively until no more merges are possible.
    """
    if max_gap <= 0 or len(segments) < 3:
        return segments

    dup_states = {'LowDup', 'HighDup', 'Amp', 'HighAmp', 'MedAmp', 'ExtremeAmp'}
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
                    segments[i+1]['state'] in ('Neutral', 'NoisyNeutral') and
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
                    'cn_std': _pooled_std(
                        s1.get('cn_std', 0.0), s1['n_windows'],
                        s2.get('cn_std', 0.0), s2['n_windows'],
                        s1['cn_median'], s2['cn_median'],
                    ),
                    'avg_entropy': s1.get('avg_entropy', 0.0) * w1 + s2.get('avg_entropy', 0.0) * w2,
                    'max_entropy': max(s1.get('max_entropy', 0.0), s2.get('max_entropy', 0.0)),
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

    # Merge adjacent Neutral/NoisyNeutral segments (same-state only)
    merged = [result[0]]
    for seg in result[1:]:
        prev = merged[-1]
        if (seg['state'] == prev['state'] and
                seg['state'] in ('Neutral', 'NoisyNeutral') and
                seg['chrom'] == prev['chrom']):
            total_bp = (prev['end'] - prev['start']) + (seg['end'] - seg['start'])
            prev_weight = (prev['end'] - prev['start']) / total_bp
            curr_weight = (seg['end'] - seg['start']) / total_bp

            merged[-1] = {
                'chrom': prev['chrom'],
                'start': prev['start'],
                'end': seg['end'],
                'state': prev['state'],
                'cn_median': prev['cn_median'] * prev_weight + seg['cn_median'] * curr_weight,
                'cn_mean': prev['cn_mean'] * prev_weight + seg['cn_mean'] * curr_weight,
                'n_windows': prev['n_windows'] + seg['n_windows'],
                'avg_entropy': prev.get('avg_entropy', 0.0) * prev_weight + seg.get('avg_entropy', 0.0) * curr_weight,
                'max_entropy': max(prev.get('max_entropy', 0.0), seg.get('max_entropy', 0.0)),
            }
        else:
            merged.append(seg)

    if len(merged) != len(result):
        print(f"[POST] Merged adjacent Neutral/NoisyNeutral: {len(result)} -> {len(merged)}")

    return merged


def filter_segments_adaptive(segments):
    """Quality-scaled segment filtering (A2 + A4 combined).

    A2: Per-state base thresholds (from Round 13) scaled by segment quality.
        high quality (>0.95) → 0.5x threshold, medium → 1.0x, low (<0.80) → 2.0x
    A4: CN consistency check — very noisy dup segments (CV > threshold) reclassified.
    """
    dup_states = {'LowDup', 'HighDup', 'Amp', 'HighAmp', 'MedAmp', 'ExtremeAmp'}
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

    # Merge adjacent Neutral/NoisyNeutral segments (same-state only)
    merged = [result[0]]
    for seg in result[1:]:
        prev = merged[-1]
        if (seg['state'] == prev['state'] and
                seg['state'] in ('Neutral', 'NoisyNeutral') and
                seg['chrom'] == prev['chrom']):
            total_bp = (prev['end'] - prev['start']) + (seg['end'] - seg['start'])
            prev_weight = (prev['end'] - prev['start']) / total_bp
            curr_weight = (seg['end'] - seg['start']) / total_bp

            merged[-1] = {
                'chrom': prev['chrom'],
                'start': prev['start'],
                'end': seg['end'],
                'state': prev['state'],
                'cn_median': prev['cn_median'] * prev_weight + seg['cn_median'] * curr_weight,
                'cn_mean': prev['cn_mean'] * prev_weight + seg['cn_mean'] * curr_weight,
                'n_windows': prev['n_windows'] + seg['n_windows'],
                'avg_quality': prev.get('avg_quality', 1.0) * prev_weight + seg.get('avg_quality', 1.0) * curr_weight,
                'cn_std': _pooled_std(
                    prev.get('cn_std', 0.0), prev['n_windows'],
                    seg.get('cn_std', 0.0), seg['n_windows'],
                    prev['cn_median'], seg['cn_median'],
                ),
                'avg_entropy': prev.get('avg_entropy', 0.0) * prev_weight + seg.get('avg_entropy', 0.0) * curr_weight,
                'max_entropy': max(prev.get('max_entropy', 0.0), seg.get('max_entropy', 0.0)),
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

    Extended format includes quality, CN std, entropy, and — if repeat annotation
    was loaded — masked_fraction and repeat_class columns.
    """
    print(f"[IO] Writing {len(segments):,} segments to {output_file}...")

    has_repeat = any('repeat_class' in seg for seg in segments[:10])

    with open(output_file, 'w') as f:
        if extended:
            header = ("#chrom\tstart\tend\tstate\tcn_median\tcn_mean\tn_windows\t"
                      "avg_quality\tmin_quality\tcn_std\tavg_repeats\t"
                      "avg_entropy\tmax_entropy")
            if has_repeat:
                header += "\tmasked_fraction\trepeat_class"
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

    for state_id in range(len(STATES)):
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
# GC CORRECTION
# =============================================================================

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


def apply_gc_correction(df, gc_df, poly_degree=3):
    """Apply GC bias correction to log2_cn values.

    Strategy:
      1. Merge GC content per window into the main DataFrame.
      2. Identify approximately-neutral windows (log2_cn in [-0.75, 0.75]).
      3. Fit a polynomial: log2_cn ~ gc_content on those neutral windows.
      4. Subtract expected log2_cn(gc) from all windows.

    This removes GC-bias artifacts (e.g., inflated CN in GC-rich pericentromeric
    regions) while preserving real CNV signal.

    Note: correction is derived from raw log2(cn) signal but applied to the
    quality-weighted log2_cn used by the HMM.
    """
    print(f"[GC] Applying GC correction (polynomial degree={poly_degree})...")

    # Merge GC content (left join — windows without GC data keep gc=0.5)
    n_before = len(df)
    df = df.merge(gc_df, on=['chrom', 'start'], how='left')
    df['gc_content'] = df['gc_content'].fillna(0.5)
    assert len(df) == n_before, "GC merge changed row count — check window coordinates"

    # Fit on approximately-neutral windows
    neutral_mask = (df['log2_cn'] >= -0.75) & (df['log2_cn'] <= 0.75)
    n_neutral = neutral_mask.sum()
    gc_has_signal = df['gc_content'].std() > 0.01

    if n_neutral < 100 or not gc_has_signal:
        print(f"[GC] WARNING: {n_neutral} neutral windows or no GC variation — skipping correction")
        return df

    print(f"[GC] Fitting polynomial on {n_neutral:,} neutral windows...")
    gc_neutral = df.loc[neutral_mask, 'gc_content'].values
    log2_neutral = df.loc[neutral_mask, 'log2_cn'].values

    coeffs = np.polyfit(gc_neutral, log2_neutral, poly_degree)
    poly = np.poly1d(coeffs)

    # Apply correction
    expected = poly(df['gc_content'].values)
    df['log2_cn'] = df['log2_cn'] - expected

    print(f"[GC] Correction applied. Bias range: [{expected.min():.3f}, {expected.max():.3f}]")
    print(f"[GC] Post-correction log2_cn: mean={df['log2_cn'].mean():.4f}, "
          f"std={df['log2_cn'].std():.4f}")
    return df


# =============================================================================
# REPEAT ANNOTATION — load and attach to segments
# =============================================================================

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
    print(f"[RM] Loaded {len(rm_df):,} windows "
          f"(>50% masked: {pct_masked:.1f}%)")
    return rm_df


def compute_segment_repeat_annotation(segments, df, rm_df):
    """Attach per-segment repeat annotation (masked_fraction, repeat_class).

    For each segment, computes:
      - masked_fraction: mean of per-window masked_fraction across segment windows
      - repeat_class: most common repeat class by window count
        (tie: Satellite > Simple_repeat > LINE > SINE > LTR > DNA > …)

    Uses binary search on sorted start positions (same approach as compute_segment_mean_gc).
    """
    CLASS_PRIORITY = ['Satellite', 'Simple_repeat', 'LINE', 'SINE', 'LTR',
                      'DNA', 'Low_complexity', 'Other', 'None']

    print(f"[RM] Computing repeat annotation for {len(segments):,} segments...")

    # Build per-chrom sorted arrays for binary search
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

        # Dominant class
        from collections import Counter
        cls_counts = Counter(cls_arr[lo:hi])
        if cls_counts:
            # Among tied classes, pick highest priority
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
# POST-SEGMENTATION GC CN CALIBRATION
# =============================================================================

def compute_segment_mean_gc(segments, df, gc_df):
    """Compute mean GC content for each segment from per-window GC data.

    Uses binary search (O(S log W)) instead of per-segment full scan (O(S*W)).
    For 157K segments × 6M windows the naive approach takes ~26 hours;
    this takes seconds.

    Mutates segments in-place, adding 'mean_gc' key.
    """
    # Merge GC into windows once
    n_before = len(df)
    df_gc = df.merge(gc_df[['chrom', 'start', 'gc_content']], on=['chrom', 'start'], how='left')
    df_gc['gc_content'] = df_gc['gc_content'].fillna(0.5)
    assert len(df_gc) == n_before, "GC merge changed row count"

    print(f"[GC-CAL] Computing mean GC per segment for {len(segments):,} segments "
          f"(binary search)...")

    # Build per-chromosome sorted arrays once — O(W log W) total
    gc_by_chrom = {}
    for chrom, gdf in df_gc.groupby('chrom', sort=False):
        gdf_s = gdf.sort_values('start')
        gc_by_chrom[chrom] = (
            gdf_s['start'].values.astype(np.int64),
            gdf_s['gc_content'].values,
        )

    # Per segment: two binary searches — O(S log W) total
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
            # Segment is smaller than one window — use nearest window
            idx = min(int(np.searchsorted(starts, (s + e) // 2)), len(starts) - 1)
            seg['mean_gc'] = float(gcs[idx])

    print(f"[GC-CAL] Done.")
    return segments


def apply_gc_cn_calibration(segments):
    """Post-segmentation GC-based CN calibration.

    Strategy:
      1. Use Neutral segments as ground truth: at neutral CN, the expected value
         should be ~1.0 regardless of GC. Observed deviation reveals GC bias.
      2. Fit a polynomial: observed_cn ~ mean_gc on Neutral segments.
      3. For all segments: cn_corrected = cn_median / poly(mean_gc).
         This divides out the GC-dependent bias factor.

    This approach is safe for high-CN amplifications because:
      - It does NOT modify the HMM signal (no BISER-SD destruction)
      - Division preserves the relative CN magnitude at all scales
      - At CN=150 with GC-bias=1.1: corrected = 150/1.1 = 136 (not 150/1.1^QW)

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

    # Filter outlier neutral segments (cn very far from 1.0 — likely misclassified)
    cn_ok = (cn_vals > 0.5) & (cn_vals < 2.0)
    if cn_ok.sum() < 20:
        print(f"[GC-CAL] WARNING: too few well-behaved neutral segments — skipping calibration")
        return segments

    gc_vals, cn_vals = gc_vals[cn_ok], cn_vals[cn_ok]

    # Fit polynomial: expected_cn(gc) on neutral segments
    poly_degree = 3
    coeffs = np.polyfit(gc_vals, cn_vals, poly_degree)
    poly = np.poly1d(coeffs)

    # Diagnostics
    pred = poly(gc_vals)
    residual_std = float(np.std(cn_vals - pred))
    gc_range_pred = poly(np.array([gc_vals.min(), gc_vals.max()]))
    print(f"[GC-CAL] Polynomial degree={poly_degree}, fitted on {len(gc_vals)} neutral segs")
    print(f"[GC-CAL] GC range: [{gc_vals.min():.3f}, {gc_vals.max():.3f}]  "
          f"CN bias range: [{gc_range_pred.min():.3f}, {gc_range_pred.max():.3f}]")
    print(f"[GC-CAL] Residual std on neutral segs: {residual_std:.4f}")

    # GC BYPASS: Polynomial training range'i dışındaki bölgeler için kalibrasyon atla
    # Satellite (~67% GC), rDNA (67% GC): training range (max ~55% GC) dışında → unreliable
    GC_BYPASS_CLASSES = frozenset({'Satellite', 'Low_complexity'})

    # Apply to all segments
    n_corrected = 0
    n_clamped = 0
    n_bypassed = 0
    for seg in segments:
        if 'mean_gc' not in seg:
            continue

        # Bypass koşulları (sırasıyla kontrol edilir):
        skip_gc = False
        rc = seg.get('repeat_class', 'None')
        mf = seg.get('masked_fraction', 0.0)
        gc = seg.get('mean_gc', 0.5)

        if rc in GC_BYPASS_CLASSES:
            skip_gc = True                       # Satellite → bypass
        elif gc > 0.60:
            skip_gc = True                       # GC > 60%: polynomial training range dışı
        elif mf > 0.85 and gc > 0.55:
            skip_gc = True                       # Yüksek masked + yüksek GC

        if skip_gc:
            seg['gc_bias_factor'] = 1.0          # Düzeltme uygulanmadı → factor=1.0
            n_bypassed += 1
            continue                             # cn_median değiştirilmez

        gc_factor_raw = float(poly(seg['mean_gc']))
        # Clamp correction factor to [0.6, 1.8] — reject extreme corrections
        gc_factor = max(0.6, min(1.8, gc_factor_raw))
        if gc_factor != gc_factor_raw:
            n_clamped += 1
        seg['cn_median'] = seg['cn_median'] / gc_factor
        seg['cn_mean'] = seg['cn_mean'] / gc_factor
        if seg.get('cn_std', 0.0) > 0:
            seg['cn_std'] = seg['cn_std'] / gc_factor
        seg['gc_bias_factor'] = gc_factor
        n_corrected += 1

    print(f"[GC-CAL] Applied to {n_corrected:,} segments "
          f"({n_bypassed:,} bypassed: Satellite/high-GC)")
    if n_clamped > 0:
        print(f"[GC-CAL] WARNING: {n_clamped} segments had GC correction factor "
              f"clamped to [0.6, 1.8] — polynomial extrapolation outside training GC range. "
              f"Consider extending poly_degree or inspecting GC distribution.")
    return segments


def reclassify_by_cn_threshold(segments, lowdup_threshold=1.87, hetdel_threshold=0.75):
    """Post-GC calibration CN-based state correction.

    Uses CN_STATE_BOUNDARIES via assign_state_from_cn() for routing.
    With NoisyNeutral (8-state model), CN 1.20-1.87 maps to NoisyNeutral.

    Rules:
    1. LowDup CN < threshold → assign_state_from_cn (NoisyNeutral if 1.20-1.87, Neutral if <1.20)
    2. HetDel CN > 1.5 → LowDup (rare mis-assignment)
    3. HetDel CN > hetdel_threshold → Neutral (GC-bias artefact)
    4. Neutral CN > 1.20 → assign_state_from_cn (NoisyNeutral if 1.20-1.87, LowDup+ if >1.87)
    5. NoisyNeutral CN < 1.20 → Neutral
    6. NoisyNeutral CN > 1.87 → assign_state_from_cn (LowDup+)
    7. Universal fallback for all other dup states
    """
    recl_lowdup = 0
    recl_hetdel_neutral = 0
    recl_hetdel_lowdup = 0
    recl_neutral_up = 0
    recl_nn_down = 0
    recl_nn_up = 0
    recl_universal = 0
    result = []
    for seg in segments:
        cn = seg.get('cn_median', 1.0)
        if seg['state'] == 'LowDup' and cn < lowdup_threshold:
            seg = dict(seg)
            seg['state'] = assign_state_from_cn(cn)
            recl_lowdup += 1
        elif seg['state'] == 'HetDel' and cn > 1.5:
            seg = dict(seg)
            seg['state'] = assign_state_from_cn(cn)
            recl_hetdel_lowdup += 1
        elif seg['state'] == 'HetDel' and cn > hetdel_threshold:
            seg = dict(seg)
            seg['state'] = 'Neutral'
            recl_hetdel_neutral += 1
        elif seg['state'] == 'Neutral' and cn > 1.20:
            seg = dict(seg)
            seg['state'] = assign_state_from_cn(cn)
            recl_neutral_up += 1
        elif seg['state'] == 'NoisyNeutral' and cn < 1.20:
            seg = dict(seg)
            seg['state'] = 'Neutral'
            recl_nn_down += 1
        elif seg['state'] == 'NoisyNeutral' and cn > 1.87:
            seg = dict(seg)
            seg['state'] = assign_state_from_cn(cn)
            recl_nn_up += 1
        else:
            # Universal fallback: correct state/CN mismatches not covered above.
            if seg['state'] not in ('Neutral', 'NoisyNeutral', 'HetDel', 'HomDel'):
                correct_state = assign_state_from_cn(cn)
                if correct_state != seg['state']:
                    seg = dict(seg)
                    seg['state'] = correct_state
                    recl_universal += 1
        result.append(seg)

    if recl_lowdup > 0:
        print(f"[CN-RECL] Reclassified {recl_lowdup:,} LowDup→lower "
              f"(cn_median < {lowdup_threshold:.2f})")
    if recl_hetdel_neutral > 0:
        print(f"[CN-RECL] Reclassified {recl_hetdel_neutral:,} HetDel→Neutral "
              f"(cn_median > {hetdel_threshold:.2f}, GC-bias artefact)")
    if recl_hetdel_lowdup > 0:
        print(f"[CN-RECL] Reclassified {recl_hetdel_lowdup:,} HetDel→Dup "
              f"(cn_median > 1.50)")
    if recl_neutral_up > 0:
        print(f"[CN-RECL] Reclassified {recl_neutral_up:,} Neutral→higher "
              f"(cn_median > 1.20)")
    if recl_nn_down > 0:
        print(f"[CN-RECL] Reclassified {recl_nn_down:,} NoisyNeutral→Neutral "
              f"(cn_median < 1.20)")
    if recl_nn_up > 0:
        print(f"[CN-RECL] Reclassified {recl_nn_up:,} NoisyNeutral→Dup "
              f"(cn_median > 1.87)")
    if recl_universal > 0:
        print(f"[CN-RECL] Universal reclassify: {recl_universal:,} segments corrected via assign_state_from_cn")

    # Merge adjacent Neutral/NoisyNeutral segments created by reclassification
    if not result:
        return result
    merged = [result[0]]
    for seg in result[1:]:
        prev = merged[-1]
        if (prev['state'] == seg['state']
                and prev['state'] in ('Neutral', 'NoisyNeutral')
                and prev['chrom'] == seg['chrom']
                and seg['start'] <= prev['end'] + 1):
            total_bp = (prev['end'] - prev['start']) + (seg['end'] - seg['start'])
            pw = (prev['end'] - prev['start']) / total_bp
            cw = (seg['end'] - seg['start']) / total_bp
            merged[-1] = {
                'chrom': prev['chrom'],
                'start': prev['start'],
                'end': seg['end'],
                'state': prev['state'],
                'cn_median': prev['cn_median'] * pw + seg['cn_median'] * cw,
                'cn_mean': prev.get('cn_mean', prev['cn_median']) * pw + seg.get('cn_mean', seg['cn_median']) * cw,
                'n_windows': prev['n_windows'] + seg['n_windows'],
                'avg_quality': prev.get('avg_quality', 1.0) * pw + seg.get('avg_quality', 1.0) * cw,
                'cn_std': _pooled_std(
                    prev.get('cn_std', 0.0), prev['n_windows'],
                    seg.get('cn_std', 0.0), seg['n_windows'],
                    prev['cn_median'], seg['cn_median'],
                ),
                'avg_entropy': prev.get('avg_entropy', 0.0) * pw + seg.get('avg_entropy', 0.0) * cw,
                'max_entropy': max(prev.get('max_entropy', 0.0), seg.get('max_entropy', 0.0)),
            }
        else:
            merged.append(seg)

    print(f"[CN-RECL] After merge: {len(result) - len(merged):,} adjacent same-state segments merged")
    return merged


def assign_state_from_cn(cn):
    """Map a calibrated CN value to the nearest HMM state using CN_STATE_BOUNDARIES."""
    for upper, state in CN_STATE_BOUNDARIES:
        if cn < upper:
            return state
    return 'HighAmp'


def split_high_cv_segments(segments, df, cv_threshold=0.6, min_length=3000, max_depth=3):
    """Split dup/amp segments with high within-segment CN variability.

    High CV (cn_std / cn_median > cv_threshold) signals that the segment spans
    a CN transition — e.g. two distinct-CN regions merged across a gap, or a
    centromere boundary.  Splitting produces smaller segments with lower internal
    variance and more accurate CN assignment.

    Algorithm:
      For each dup segment with CV > cv_threshold:
        1. Extract its raw (pre-GC-calibration) window CN values using binary search.
        2. Find the split index that minimises the weighted within-group variance
           (computed in O(n) via prefix sums).
        3. Build two sub-segments; re-assign their state from CN_STATE_BOUNDARIES.
        4. Apply the parent's gc_bias_factor to calibrate sub-segment CN values.
        5. Recurse up to max_depth times.

    Args:
        segments:     List of segment dicts.
        df:           Original windows DataFrame (columns: chrom, start, end, cn).
        cv_threshold: CV above which a segment is a split candidate (default 0.6).
        min_length:   Minimum sub-segment length in bp (default 3000).
        max_depth:    Maximum recursion depth (default 3, prevents over-splitting).

    Returns:
        New list of segments (always >= len(segments)).
    """
    import time as _t
    _t0 = _t.time()

    DUP_STATES = {'LowDup', 'HighDup', 'Amp', 'HighAmp', 'MedAmp', 'ExtremeAmp'}

    # Build per-chromosome index: chrom → (sorted start array, raw cn array)
    idx = {}
    for chrom, gdf in df.groupby('chrom', sort=False):
        gdf_s = gdf.sort_values('start')
        idx[chrom] = (
            gdf_s['start'].values.astype(np.int64),
            gdf_s['cn'].values.astype(np.float64),
        )

    def _find_best_split(chrom, seg_start, seg_end, gc_factor):
        """Return (left_seg_dict, right_seg_dict) or None if no valid split."""
        if chrom not in idx:
            return None
        starts, raw_cns = idx[chrom]

        lo = int(np.searchsorted(starts, seg_start, side='left'))
        hi = int(np.searchsorted(starts, seg_end,   side='left'))
        n  = hi - lo
        if n < 4:   # need ≥2 windows on each side
            return None

        win_cns = raw_cns[lo:hi]

        # Prefix sums for O(n) variance computation
        ps    = np.cumsum(win_cns)           # ps[k]  = sum(win_cns[0:k+1])
        ps_sq = np.cumsum(win_cns ** 2)

        best_score = float('inf')
        best_k     = -1

        for k in range(2, n - 2):
            left_end_bp  = int(starts[lo + k])
            right_end_bp = seg_end

            if (left_end_bp  - seg_start) < min_length:
                continue
            if (right_end_bp - left_end_bp) < min_length:
                break   # subsequent k values only make left longer

            # Left: windows 0..k-1
            sl   = float(ps[k - 1])
            sql  = float(ps_sq[k - 1])
            varl = sql / k - (sl / k) ** 2

            # Right: windows k..n-1
            sr   = float(ps[n - 1]) - float(ps[k - 1])
            sqr  = float(ps_sq[n - 1]) - float(ps_sq[k - 1])
            nr   = n - k
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
            raw_std = float(np.std(cn_slice))
            cal_med = raw_med / gc_factor
            cal_std = raw_std / gc_factor
            sub = dict(parent)
            sub['start']     = s
            sub['end']       = e
            sub['cn_median'] = cal_med
            sub['cn_mean']   = float(np.mean(cn_slice)) / gc_factor
            sub['cn_std']    = cal_std
            sub['num_windows'] = len(cn_slice)
            sub['state']     = assign_state_from_cn(cal_med)
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
        cv   = cn_s / (cn_m + 1e-6)
        if cv <= cv_threshold:
            return [seg]

        gc_factor = seg.get('gc_bias_factor', 1.0)
        result = _find_best_split(seg['chrom'], seg['start'], seg['end'], gc_factor)
        if result is None:
            return [seg]

        left, right = result
        # Carry over keys that _make_sub left empty (chrom, avg_quality, etc.)
        for key, val in seg.items():
            left.setdefault(key, val)
            right.setdefault(key, val)
        # Correct start/end overridden by setdefault
        left['start'],  left['end']  = seg['start'], result[0]['end']
        right['start'], right['end'] = result[1]['start'], seg['end']

        return _recursive_split(left, depth + 1) + _recursive_split(right, depth + 1)

    result      = []
    n_candidate = 0
    n_extra     = 0

    for seg in segments:
        cn_m = seg.get('cn_median', 1.0)
        cn_s = seg.get('cn_std',    0.0)
        cv   = cn_s / (cn_m + 1e-6)

        if seg['state'] not in DUP_STATES or cv <= cv_threshold:
            result.append(seg)
            continue

        n_candidate += 1
        sub = _recursive_split(seg, 0)
        n_extra += len(sub) - 1
        result.extend(sub)

    elapsed = _t.time() - _t0
    print(f"[CV-SPLIT] {n_candidate:,} high-CV dup segments processed → "
          f"{n_extra:,} additional segments created ({elapsed:.1f}s)")
    print(f"[CV-SPLIT] Total: {len(segments):,} → {len(result):,} segments")
    return result


def filter_low_quality_segments(segments, threshold, entropy_threshold=1.5):
    """Hard quality + entropy filter: reclassify uncertain dup segments to Neutral.

    Two independent criteria (both applied):
      1. avg_quality < threshold  — k-mer quality too low (W1/W7 fix)
      2. max_entropy > entropy_threshold — HMM is uncertain about state boundary;
         high-entropy segments span multiple states → unreliable CN (W3 fix)
         entropy_threshold=1.5 bits (out of max 2.81 for 7 states).

    Unlike quality weighting (which biases CN values), this acts as a binary filter:
    segments that pass have intact CN values.

    Used in cn-accuracy mode instead of quality weighting.
    """
    if threshold <= 0 and entropy_threshold <= 0:
        return segments

    dup_states = {'LowDup', 'HighDup', 'Amp', 'HighAmp', 'MedAmp', 'ExtremeAmp'}
    reclassified_q = 0
    reclassified_e = 0
    result = []
    for seg in segments:
        if seg['state'] in dup_states:
            q = seg.get('avg_quality', 1.0)
            e = seg.get('max_entropy', 0.0)
            if threshold > 0 and q < threshold:
                seg = dict(seg)
                seg['state'] = 'Neutral'
                reclassified_q += 1
            elif entropy_threshold > 0 and e > entropy_threshold:
                seg = dict(seg)
                seg['state'] = 'Neutral'
                reclassified_e += 1
        result.append(seg)

    reclassified = reclassified_q + reclassified_e
    if reclassified_q > 0:
        print(f"[POST] Quality filter (threshold={threshold:.2f}): "
              f"reclassified {reclassified_q} low-quality dup segments to Neutral")
    if reclassified_e > 0:
        print(f"[POST] Entropy filter (threshold={entropy_threshold:.2f} bits): "
              f"reclassified {reclassified_e} high-entropy dup segments to Neutral")

    # Merge adjacent Neutral/NoisyNeutral segments (same-state only)
    if not result:
        return result
    merged = [result[0]]
    for seg in result[1:]:
        prev = merged[-1]
        if (seg['state'] == prev['state'] and
                seg['state'] in ('Neutral', 'NoisyNeutral') and
                seg['chrom'] == prev['chrom']):
            total_bp = (prev['end'] - prev['start']) + (seg['end'] - seg['start'])
            pw = (prev['end'] - prev['start']) / total_bp
            cw = (seg['end'] - seg['start']) / total_bp
            merged[-1] = {
                'chrom': prev['chrom'],
                'start': prev['start'],
                'end': seg['end'],
                'state': prev['state'],
                'cn_median': prev['cn_median'] * pw + seg['cn_median'] * cw,
                'cn_mean': prev['cn_mean'] * pw + seg['cn_mean'] * cw,
                'n_windows': prev['n_windows'] + seg['n_windows'],
                'avg_quality': prev.get('avg_quality', 1.0) * pw + seg.get('avg_quality', 1.0) * cw,
                'cn_std': _pooled_std(
                    prev.get('cn_std', 0.0), prev['n_windows'],
                    seg.get('cn_std', 0.0), seg['n_windows'],
                    prev['cn_median'], seg['cn_median'],
                ),
                'avg_entropy': prev.get('avg_entropy', 0.0) * pw + seg.get('avg_entropy', 0.0) * cw,
                'max_entropy': max(prev.get('max_entropy', 0.0), seg.get('max_entropy', 0.0)),
            }
        else:
            merged.append(seg)

    if len(merged) != len(result):
        print(f"[POST] Merged adjacent Neutral/NoisyNeutral after quality filter: "
              f"{len(result)} -> {len(merged)}")
    return merged


# =============================================================================
# MAIN
# =============================================================================

def apply_mode(mode):
    """Apply mode-specific configuration to global variables."""
    global STATES, SELF_TRANSITION, QUALITY_WEIGHT, MAX_MERGE_GAP, MIN_SEGMENT_LENGTHS
    global QUALITY_THRESHOLD, START_PROBS

    if mode not in MODE_CONFIGS:
        print(f"[ERROR] Unknown mode: {mode}. Use 'precision', 'sensitive', or 'cn-accuracy'.")
        sys.exit(1)

    cfg = MODE_CONFIGS[mode]
    STATES = cfg['states']
    SELF_TRANSITION = cfg['self_transition']
    QUALITY_WEIGHT = cfg['quality_weight']
    MAX_MERGE_GAP = cfg['max_merge_gap']
    MIN_SEGMENT_LENGTHS = cfg['min_segment_lengths']
    QUALITY_THRESHOLD = cfg.get('quality_threshold', 0.0)
    # cn-accuracy uses MedAmp/ExtremeAmp at positions 0,1 instead of HomDel/HetDel.
    # Start probabilities reflect that chromosomes almost always begin in Neutral.
    global START_PROBS
    if mode == 'cn-accuracy':
        START_PROBS = np.array([0.01, 0.01, 0.80, 0.03, 0.02, 0.01, 0.01, 0.11])
    else:
        START_PROBS = np.array([0.01, 0.05, 0.85, 0.05, 0.02, 0.01, 0.01])


def main():
    parser = argparse.ArgumentParser(
        description="7-State HMM CNV Segmentation (Log-Seg & Linear-Geno)"
    )
    parser.add_argument("--input", "-i", required=True,
                        help="Input BED file from Phase 1 preprocessing")
    parser.add_argument("--output", "-o", required=True,
                        help="Output segments BED file")
    parser.add_argument("--mode", "-m",
                        choices=['precision', 'sensitive', 'cn-accuracy'],
                        default='precision',
                        help="precision: best BISER F1 (QW=3.0, strict filters). "
                             "sensitive: high-recall candidates for ML (QW=1.0, liberal). "
                             "cn-accuracy: accurate CN values to CN=150 (QW=0, large segments, "
                             "hard quality filter, GC post-calibration).")
    parser.add_argument("--no-pomegranate", action="store_true",
                        help="Force use of simple Viterbi (skip pomegranate)")
    parser.add_argument("--extended", action="store_true",
                        help="Write extended columns (quality, cn_std, entropy)")
    parser.add_argument("--coarse-output",
                        help="Also run HMM on 5000bp coarse windows and write segments here "
                             "(for multi-resolution consensus features in the ML classifier)")
    parser.add_argument("--gc-content-bed",
                        help="GC content BED (chrom, start, end, gc_content) from "
                             "compute_gc_content.py. "
                             "In cn-accuracy mode: post-segmentation CN calibration (correct). "
                             "In other modes: pre-HMM log2_cn correction (use with caution).")
    parser.add_argument("--repeat-bed",
                        help="Repeat-annotated window BED from compute_repeat_annotation.py "
                             "(10-column: standard 8 cols + masked_fraction + repeat_class). "
                             "Adds masked_fraction and repeat_class columns to output segments. "
                             "Does NOT change HMM state assignment — metadata only.")
    parser.add_argument("--quality-threshold", type=float, default=None,
                        help="Override hard quality threshold for post-segmentation filter "
                             "(cn-accuracy default: 0.7; others: 0.0 = disabled). "
                             "Dup segments with avg_quality below threshold → reclassified Neutral.")
    parser.add_argument("--no-chrx-correction", action="store_true",
                        help="Disable chrX-specific HMM state correction. "
                             "Use for female samples (e.g., CHM13) where chrX is diploid "
                             "and should use the same model as autosomes.")
    parser.add_argument("--cn-reclassify-threshold", type=float, default=None,
                        help="Post-GC calibration CN threshold: LowDup segments with "
                             "cn_median below this value are reclassified to Neutral. "
                             "Default in cn-accuracy mode: 1.5 (midpoint between "
                             "Neutral CN=1.0 and LowDup CN=2.0).  Set to 0 to disable.")
    parser.add_argument("--cv-split-threshold", type=float, default=None,
                        help="CV threshold for post-merge segment splitting.  Dup/amp "
                             "segments with cn_std/cn_median above this value are split "
                             "at the variance-minimising breakpoint (recursive, max 3 "
                             "levels).  Default in cn-accuracy mode: 0.6.  "
                             "Set to 0 to disable.")
    parser.add_argument("--min-kmers", type=int, default=0,
                        help="Minimum num_kmers per 500bp window.  Windows below this "
                             "threshold have unreliable CN due to sparse k-mer coverage "
                             "(e.g. complex/GC-rich regions at k=72) and are set to "
                             "CN=1.0 (Neutral) before HMM to suppress false HetDel calls. "
                             "Recommended: 30 for CHM13 ONT k=72.  Default: 0 (disabled).")

    args = parser.parse_args()

    # Apply mode configuration
    apply_mode(args.mode)

    # Override quality threshold if explicitly provided
    if args.quality_threshold is not None:
        global QUALITY_THRESHOLD
        QUALITY_THRESHOLD = args.quality_threshold

    mode_labels = {
        'sensitive':   'HIGH-RECALL CANDIDATE GENERATION',
        'precision':   'PRECISION (Best BISER F1)',
        'cn-accuracy': 'CN ACCURACY (QW=0, large segments)',
    }
    mode_label = mode_labels.get(args.mode, args.mode.upper())
    print("=" * 60)
    print(f"KonuSeg {len(STATES)}-State HMM Segmentation — {mode_label}")
    print("=" * 60)
    print(f"Input:  {args.input}")
    print(f"Output: {args.output}")
    print(f"Mode:   {args.mode}")
    print()
    print("HMM Parameters:")
    print(f"  SELF_TRANSITION:    {SELF_TRANSITION}")
    print(f"  LowDup mean:        {STATES[3]['log2_mean']} (CN~{2**STATES[3]['log2_mean']:.2f})")
    print(f"  LowDup var:         {STATES[3]['log2_var']}")
    print(f"  Neutral var:        {STATES[2]['log2_var']}")
    if 7 in STATES:
        print(f"  NoisyNeutral mean:  {STATES[7]['log2_mean']} (CN~{2**STATES[7]['log2_mean']:.2f})")
        print(f"  NoisyNeutral var:   {STATES[7]['log2_var']}")
    print(f"  QUALITY_WEIGHT:     {QUALITY_WEIGHT}")
    print(f"  MIN_SEGMENT_LENGTHS:{MIN_SEGMENT_LENGTHS}")
    print(f"  MAX_MERGE_GAP:      {MAX_MERGE_GAP}")
    print(f"  QUALITY_THRESHOLD:  {QUALITY_THRESHOLD} "
          f"({'disabled' if QUALITY_THRESHOLD <= 0 else 'hard filter'})")
    print()

    # Load data
    df = load_data(args.input, min_kmers=args.min_kmers)

    # GC pre-HMM correction: only in non-cn-accuracy modes (known to harm SD detection)
    gc_df = None
    if args.gc_content_bed:
        gc_df = load_gc_content(args.gc_content_bed)
        if args.mode != 'cn-accuracy':
            df = apply_gc_correction(df, gc_df)
            gc_df = None  # already applied, don't do post-calibration
            print()
        else:
            print(f"[GC] cn-accuracy mode: GC content loaded for post-segmentation calibration")
            print()

    # Run segmentation
    use_pom = POMEGRANATE_AVAILABLE and not args.no_pomegranate
    apply_chrx = not args.no_chrx_correction
    segments = run_segmentation(df, use_pomegranate=use_pom, chrx_correction=apply_chrx)

    # Post-processing
    if MAX_MERGE_GAP > 0 or MIN_SEGMENT_LENGTHS:
        print(f"\n[POST] Post-processing: merge_gap={MAX_MERGE_GAP}, "
              f"min_lengths={MIN_SEGMENT_LENGTHS}")
        segments = filter_small_segments(segments, MIN_SEGMENT_LENGTHS)

    # Hard quality + entropy filter BEFORE merge (W7 fix: fragmented low-quality
    # segments must be removed first so they cannot be merged into composites that
    # then appear high-quality; W3 fix: high-entropy segments indicate HMM boundary
    # uncertainty and are reclassified to Neutral before merging propagates them)
    if QUALITY_THRESHOLD > 0:
        segments = filter_low_quality_segments(segments, QUALITY_THRESHOLD)

    if MAX_MERGE_GAP > 0 or MIN_SEGMENT_LENGTHS:
        segments = merge_nearby_dup_segments(segments, MAX_MERGE_GAP)

    # Post-segmentation GC CN calibration (cn-accuracy mode only)
    if gc_df is not None:
        print()
        segments = compute_segment_mean_gc(segments, df, gc_df)
        segments = apply_gc_cn_calibration(segments)

    # CN-based state reclassification (cn-accuracy mode):
    #   LowDup segments with cn_median < threshold → Neutral.
    #   Default threshold=1.5 (midpoint Neutral=1.0 / LowDup=2.0).
    #   Targets biologically expected ~45-55% Neutral fraction.
    if args.mode == 'cn-accuracy':
        cn_recl_threshold = args.cn_reclassify_threshold
        if cn_recl_threshold is None:
            cn_recl_threshold = 1.87   # default: NoisyNeutral/LowDup boundary
        if cn_recl_threshold > 0:
            print()
            segments = reclassify_by_cn_threshold(segments, cn_recl_threshold)

    # CV-based segment splitting (cn-accuracy mode):
    #   Dup/amp segments with high internal CV are split at the variance-minimising
    #   breakpoint.  Sub-segments get their state re-assigned from CN_STATE_BOUNDARIES.
    #   Reduces HighDup/Amp CV and improves per-state CN concordance.
    if args.mode == 'cn-accuracy':
        cv_thr = args.cv_split_threshold
        if cv_thr is None:
            cv_thr = 0.6   # default
        if cv_thr > 0:
            print()
            segments = split_high_cv_segments(segments, df, cv_threshold=cv_thr)

    # Repeat annotation (metadata only — does not change HMM state assignment)
    if args.repeat_bed:
        print()
        rm_df = load_repeat_annotation(args.repeat_bed)
        segments = compute_segment_repeat_annotation(segments, df, rm_df)

    # Print statistics
    print_statistics(segments)

    # Write output — sensitive and cn-accuracy always write extended columns
    do_extended = args.mode in ('sensitive', 'cn-accuracy') or args.extended
    write_output(segments, args.output, extended=do_extended)

    # Rec 1: coarse-scale (5000bp) segmentation for multi-resolution consensus
    if args.coarse_output:
        print(f"\n[COARSE] Generating 5000bp coarse-scale segmentation...")
        coarse_df = aggregate_coarse_windows(df, coarse_size=5000)
        coarse_segs = run_segmentation(coarse_df, use_pomegranate=use_pom, chrx_correction=apply_chrx)
        if MIN_SEGMENT_LENGTHS:
            coarse_segs = filter_small_segments(coarse_segs, MIN_SEGMENT_LENGTHS)
            coarse_segs = merge_nearby_dup_segments(coarse_segs, MAX_MERGE_GAP)
        write_output(coarse_segs, args.coarse_output, extended=False)
        print(f"[COARSE] Written {len(coarse_segs):,} coarse segments to {args.coarse_output}")

    print("\n" + "=" * 60)
    print("Segmentation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()

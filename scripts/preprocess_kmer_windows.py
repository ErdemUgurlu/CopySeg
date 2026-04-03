#!/usr/bin/env python3
"""
preprocess_kmer_windows.py — K-mer BED → KonuSeg window BED
============================================================

Converts k-mer depth data (sliding windows) into the 8-column window
format used by KonuSeg segmenters (HMM and PELT).

Works with any sample, any sequencing technology (ONT/PacBio), any k-mer
size, and any chromosome naming convention (chr prefix, bare numbers,
NCBI accessions).

Input format (4 columns, 1bp-step sliding windows):
  chrom  start  end  kmer_count

Output format (8 columns, preserving input chromosome names):
  chrom  start  end  cn  mean_count  log_ratio  num_kmers  num_filtered

Normalization:
  mean_count = mean raw k-mer count per window bin
  cn = mean_count / diploid_depth_median   (Neutral regions → CN ≈ 1.0)
  log_ratio  = log2(cn + epsilon)
  num_kmers  = number of contributing k-mer windows (coverage proxy)
  num_filtered = number of k-mers soft-capped by biological filter

Notes:
  - Chromosome names from input are preserved as-is (no remapping)
  - chrM/MT excluded (mitochondrial CN distorts normalization)
  - For XY samples: --sex XY excludes chrX/chrY from median computation
  - Memory-efficient: processes in 5M-line chunks, aggregates per bin

Usage:
  python3 scripts/preprocess_kmer_windows.py \\
    --input  /path/to/kmer_bloom_filter.bed \\
    --output output/my_run/cn_w500.bed \\
    --sex XX

Author: KonuSeg Team
Date: March 2026
"""

import argparse
import array as _array
import json
import math
import os
import sys
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd

from chrom_utils import (
    MITO_PATTERNS, CHRX_PATTERNS, CHRY_PATTERNS,
    is_mito, natural_chrom_key, resolve_chrom,
)

EPSILON = 1e-3

# ============================================================================
# PER-CLASS BIO-FILTER THRESHOLD FACTORS (× Gaussian peak)
# Applied when --rm-annotation-bed is provided.
# ============================================================================
# Rationale:
#   Satellite   30×  centromeric/pericentromeric arrays; extreme CN → aggressive
#   LINE/SINE  150×  SD bölgeleri LINE/SINE-rich olabilir; global default uygun
#   Simple_rep 150×  VNTR'lar Simple_repeat; global default uygun (LPA KIV-2)
#   LTR        100×  retrotransposons; mild inflation
#   DNA/other  150×  transposons; same as global default
THRESHOLD_FACTORS = {
    'Satellite':       30.0,
    'LINE':           150.0,   # 50→150: global default; LINE-rich SD k-merleri korunur
    'SINE':           150.0,   # 50→150: global default; SINE-rich SD k-merleri korunur
    'Simple_repeat':  150.0,   # 75→150: global default; VNTR (LPA KIV-2) k-merleri korunur
    'LTR':            100.0,
    'DNA':            150.0,
    'Low_complexity': 150.0,
    'Other':          150.0,
    'None':           150.0,
}


# ============================================================================
# MULTIPLICITY WEIGHT LOADER
# ============================================================================

class WeightLoader:
    """Lazy loader for per-chromosome k-mer multiplicity weight arrays.

    Weight arrays are produced by compute_multiplicity_weights.py (float16 .npy).
    Each array maps genomic position → weight = M_ref / M_target.

    Uses LRU cache to limit memory (max_cached chromosomes at a time).
    Handles chromosome name resolution (chr/bare/CM039) via chrom_utils.
    """

    def __init__(self, weight_dir, max_cached=3):
        self.weight_dir = weight_dir
        self.max_cached = max_cached
        self._cache = OrderedDict()
        self._chrom_map = {}  # input_name → npy_chrom_name

        manifest_path = os.path.join(weight_dir, 'manifest.json')
        if not os.path.exists(manifest_path):
            raise FileNotFoundError(f"manifest.json not found in {weight_dir}")
        with open(manifest_path) as f:
            self._manifest = json.load(f)

        self.k_target = self._manifest.get('k_target', '?')
        self.k_ref = self._manifest.get('k_ref', '?')
        self._weight_chroms = set(self._manifest.get('chromosomes', {}).keys())

    def _resolve(self, chrom):
        """Map input chromosome name to weight file chromosome name."""
        if chrom in self._chrom_map:
            return self._chrom_map[chrom]
        resolved = resolve_chrom(chrom, self._weight_chroms)
        self._chrom_map[chrom] = resolved
        return resolved

    def _ensure_loaded(self, chrom):
        """Load weight array for chromosome, return it (or None)."""
        resolved = self._resolve(chrom)
        if resolved is None:
            return None
        if resolved in self._cache:
            self._cache.move_to_end(resolved)
            return self._cache[resolved]
        npy_path = os.path.join(self.weight_dir, f"{resolved}.npy")
        if not os.path.exists(npy_path):
            return None
        arr = np.load(npy_path)
        self._cache[resolved] = arr
        self._cache.move_to_end(resolved)
        if len(self._cache) > self.max_cached:
            self._cache.popitem(last=False)
        return arr

    def get_weights_bulk(self, chrom, positions):
        """Return float32 weight array for the given numpy positions array."""
        arr = self._ensure_loaded(chrom)
        weights = np.ones(len(positions), dtype=np.float32)
        if arr is None:
            return weights
        valid = (positions >= 0) & (positions < len(arr))
        weights[valid] = arr[positions[valid]].astype(np.float32)
        return weights


# ============================================================================
# STEP 0a: RM annotation helpers
# ============================================================================

def load_rm_lookup(bed_path: str) -> dict:
    """
    Load a repeat_annotated_w500.bed file produced by compute_repeat_annotation.py.
    Returns {(chrom_cm039, bin_start): repeat_class} for all annotated windows.

    Expected columns (tab-separated, with optional '#' comment header):
      chrom  start  end  cn  mean_count  log_ratio  num_kmers  num_filtered
      masked_fraction  repeat_class
    repeat_class is column index 9.
    """
    print(f"[RM] Loading repeat annotation from {bed_path}...")
    df = pd.read_csv(
        bed_path, sep='\t', comment='#', header=None,
        usecols=[0, 1, 9],
        names=['chrom', 'start', 'repeat_class'],
        dtype={'chrom': str, 'start': np.int32, 'repeat_class': str},
    )
    lookup = {
        (row.chrom, int(row.start)): row.repeat_class
        for row in df.itertuples(index=False)
    }
    print(f"[RM] Loaded {len(lookup):,} window annotations")
    return lookup


def build_bin_threshold_map(rm_lookup: dict,
                             peak: float, global_factor: float) -> dict:
    """
    Build a {(chrom, bin_start): threshold} mapping.

    RM-annotated bins  → THRESHOLD_FACTORS[repeat_class] × peak
    Unannotated bins   → global_factor × peak  (fallback, not stored)

    Only RM-annotated bins with class-specific factors are stored explicitly;
    everything else uses the fallback (global_factor × peak) at query time.
    """
    threshold_map: dict = {}

    for (chrom, bin_start), repeat_class in rm_lookup.items():
        factor = THRESHOLD_FACTORS.get(repeat_class, global_factor)
        threshold_map[(chrom, bin_start)] = peak * factor

    n_reduced = sum(1 for v in threshold_map.values()
                    if v < peak * global_factor)
    print(f"[RM] Threshold map built: "
          f"{n_reduced:,} reduced-threshold bins, "
          f"rest → {global_factor:.0f}× peak ({peak * global_factor:.1f})")
    return threshold_map


# ============================================================================
# STEP 0: Histogram sampling → Gaussian peak → biological threshold
# ============================================================================

def build_histogram_sample(filepath, sample_chunks=10, chunk_size=5_000_000,
                           weight_loader=None):
    """
    Read the first sample_chunks chunks to build a k-mer count histogram.

    A 50M-line sample (default: 10 × 5M) is sufficient to identify the
    Gaussian peak (= single-copy sequencing depth) without reading the full
    2.5B-line file.

    When weight_loader is provided, multiplicity weights are applied to
    raw counts BEFORE histogramming, so the Gaussian peak reflects the
    corrected single-copy depth (not the multiplicity-inflated depth).

    NOTE (A4): For chromosome-sorted files this covers chr1 and part of chr2.
    For CHM13 and standard ONT datasets, chr1 is representative of genome-wide
    single-copy depth — low risk in practice.  For datasets where per-chromosome
    depth varies significantly (e.g. male haploid with chrX at 0.5×), consider
    increasing --hist-sample-chunks (e.g. 30) to span more chromosomes.
    True stratified middle-file sampling is not feasible here because pandas
    sequential CSV reading cannot efficiently seek to arbitrary byte offsets
    in multi-GB compressed/text files.

    Returns: numpy int64 array of length 65536 (index = raw k-mer count)
    """
    hist = np.zeros(65536, dtype=np.int64)

    # When weight_loader is provided, we need chrom + position columns
    if weight_loader is not None:
        usecols = [0, 1, 3]
        col_names = ['chrom', 'start', 'raw_cn']
        col_dtypes = {'chrom': str, 'start': np.int32, 'raw_cn': np.float32}
    else:
        usecols = [3]
        col_names = ['raw_cn']
        col_dtypes = {'raw_cn': np.float32}

    reader = pd.read_csv(
        filepath,
        sep='\t',
        header=None,
        comment='#',
        usecols=usecols,
        names=col_names,
        dtype=col_dtypes,
        chunksize=chunk_size,
    )

    n_chunks = 0
    for chunk in reader:
        counts = chunk['raw_cn'].values.copy()

        # Apply multiplicity weights before histogramming
        if weight_loader is not None:
            for chrom_name in chunk['chrom'].unique():
                mask = chunk['chrom'].values == chrom_name
                positions = chunk['start'].values[mask]
                weights = weight_loader.get_weights_bulk(chrom_name, positions)
                counts[mask] *= weights

        counts = np.clip(counts, 0, 65535).astype(np.int32)
        np.add.at(hist, counts, 1)
        n_chunks += 1
        if n_chunks >= sample_chunks:
            break

    print(f"[HIST] Histogram built from first {n_chunks} chunks "
          f"({n_chunks * chunk_size / 1e6:.0f}M lines sampled)")
    return hist


def find_gaussian_peak(histogram, search_min=2, search_max=500):
    """
    Find the biological single-copy peak in the k-mer count histogram.

    The histogram has a large noise spike at count=1 (sequencing errors /
    singleton k-mers) that decreases monotonically before the true
    single-copy biological peak (~genome-wide median depth).

    Strategy (K2 fix, 2026-03-30):
      1. Find all local maxima in the raw histogram.
      2. Identify the noise tail end (where histogram drops to <5% of its
         value at search_min — this separates the noise region from the
         biological signal region).
      3. Among local maxima after the noise tail, pick the tallest — this is
         the biological single-copy peak.

    Previous version used a "first-increase" valley detector that could
    misfire on a single noise spike. The local-maxima approach is immune
    to isolated spikes because it considers all peaks and picks the best one
    outside the noise region.

    Returns: int — estimated single-copy depth (e.g. ~28 for CHM13 ONT k=72)
    """
    h = histogram

    # Step 1: Find all local maxima in the histogram
    peaks = []
    for i in range(max(search_min, 2), min(search_max, len(h)) - 1):
        if h[i] > h[i - 1] and h[i] >= h[i + 1] and h[i] > 0:
            peaks.append(i)

    # Step 2: Identify noise tail end (where histogram drops to <5% of start)
    start_val = float(np.max(h[search_min:min(search_min + 10, search_max)]))
    noise_threshold = start_val * 0.05
    noise_tail_end = search_min
    for i in range(search_min, min(search_max, len(h))):
        if h[i] < noise_threshold:
            noise_tail_end = i
            break

    # Step 3: Among peaks after the noise tail, pick the tallest
    bio_peaks = [p for p in peaks if p >= noise_tail_end]

    if bio_peaks:
        peak = max(bio_peaks, key=lambda p: int(h[p]))
        # Find the valley (last decrease before the peak)
        valley = noise_tail_end
        for i in range(noise_tail_end, peak):
            if h[i] <= h[i + 1]:
                valley = i
                break
    elif peaks:
        # All peaks in noise zone — use the last peak as best guess
        peak = peaks[-1]
        valley = peak - 1 if peak > 0 else search_min
    else:
        # No peaks found — fallback: first increase method
        valley = search_min
        for i in range(search_min, min(search_max, len(h) - 2)):
            if h[i + 1] >= h[i]:
                valley = i
                break
        peak = int(np.argmax(h[valley:search_max + 1])) + valley

    if peak <= 0:
        peak = search_min

    print(f"[HIST] Noise tail valley at count={valley}  "
          f"(histogram[{valley}] = {histogram[valley]:,})")
    print(f"[HIST] Gaussian peak (single-copy depth): {peak}  "
          f"(histogram[{peak}] = {histogram[peak]:,})")
    print(f"[HIST] K2: {len(peaks)} local maxima found, "
          f"{len(bio_peaks)} after noise tail (end={noise_tail_end})")
    return peak


# ============================================================================
# STEP 1: Aggregate sliding windows into bins
# ============================================================================

def _winsorized_mean(values_buf, p_low=5.0, p_high=95.0):
    """
    Compute winsorized mean from an array.array('f') buffer.
    Clips values outside [p_low, p_high] percentile range, then averages.
    Falls back to plain mean when n < 3 (percentiles undefined).
    """
    arr = np.frombuffer(values_buf, dtype=np.float32)
    if len(arr) < 3:
        return float(arr.mean()) if len(arr) > 0 else 0.0
    lo = float(np.percentile(arr, p_low))
    hi = float(np.percentile(arr, p_high))
    clipped = arr[(arr >= lo) & (arr <= hi)]
    return float(clipped.mean()) if len(clipped) > 0 else float(arr.mean())


def aggregate_windows(filepath, window_size, bio_threshold=0,
                      bin_threshold_map=None, weight_loader=None):
    """
    Read 4-column sliding-window BED in 5M-line chunks.

    For each 500bp bin:
      - k-mers with raw_count > threshold are removed (biological filter).
      - Threshold is either global (bio_threshold) or per-bin (bin_threshold_map).
      - Valid k-mer counts are stored in a compact array.array('f') buffer.
      - Winsorized mean (5th–95th percentile) is computed at write-time.

    When weight_loader is provided, multiplicity weights are applied to
    raw counts BEFORE bio-filtering (corrected = raw * weight[pos]).

    Per-bin threshold mode (bin_threshold_map):
      - SINE/LINE bins   → class-specific × peak  (Alu/L1 inflation suppressed)
      - Satellite bins   → 30 × peak  (centromeric arrays suppressed)
      - Everything else  → global bio_threshold

    Memory: ~10 GB numpy-equivalent for 2.56B k-mers after bio-filtering.
    Feasible on the 64 GB cluster node.

    Returns:
      acc_values   : {(chrom, bin_start): array.array('f')}  — valid counts
      acc_filtered : {(chrom, bin_start): int}               — n_filtered count
    """
    CHUNK_SIZE = 5_000_000

    acc_values   = defaultdict(lambda: _array.array('f'))
    acc_filtered = defaultdict(int)
    skipped_chroms = set()
    n_total = 0
    n_skipped = 0
    n_bio_filtered_total = 0
    n_mult_corrected = 0

    use_per_bin = bin_threshold_map is not None

    print(f"[IO] Reading: {filepath}")
    print(f"[IO]   Chunk size: {CHUNK_SIZE:,} lines | bin size: {window_size}bp")
    print(f"[IO]   Aggregation: winsorized mean (p5–p95)")
    if weight_loader is not None:
        print(f"[IO]   Multiplicity correction: ENABLED "
              f"(k_target={weight_loader.k_target}, k_ref={weight_loader.k_ref})")
    if use_per_bin:
        sat_f  = THRESHOLD_FACTORS.get('Satellite',  30.0)
        line_f = THRESHOLD_FACTORS.get('LINE',       150.0)
        print(f"[IO]   Biological filter: per-bin RM-guided soft-cap "
              f"(LINE/SINE={line_f:.0f}×peak, Satellite={sat_f:.0f}×peak, "
              f"default={bio_threshold:.1f}; k-mers capped, not excluded)")
    elif bio_threshold > 0:
        print(f"[IO]   Biological filter: raw_count > {bio_threshold:.1f} → soft-capped "
              f"(k-mers attenuated to threshold, not excluded)")

    reader = pd.read_csv(
        filepath,
        sep='\t',
        header=None,
        comment='#',
        usecols=[0, 1, 3],
        names=['chrom', 'start', 'raw_cn'],
        dtype={'chrom': str, 'start': np.int32, 'raw_cn': np.float32},
        chunksize=CHUNK_SIZE,
    )

    chunk_num = 0
    for chunk in reader:
        chunk_num += 1
        if chunk_num % 20 == 0:
            print(f"[IO]   ...{n_total:,} lines processed", flush=True)

        # Filter out mitochondrial chromosomes — keep everything else
        mito_mask = chunk['chrom'].isin(MITO_PATTERNS)
        n_skipped += int(mito_mask.sum())
        skipped_chroms.update(chunk.loc[mito_mask, 'chrom'].unique())
        chunk = chunk[~mito_mask].copy()
        if chunk.empty:
            continue
        n_total += len(chunk)

        # === MULTIPLICITY WEIGHT CORRECTION ===
        # Applied BEFORE bio-filter so thresholds apply to corrected values.
        # weight[pos] = M_ref / M_target; corrected = raw * weight
        if weight_loader is not None:
            raw_cn_vals = chunk['raw_cn'].values.copy()
            for chrom_name in chunk['chrom'].unique():
                mask = chunk['chrom'].values == chrom_name
                positions = chunk['start'].values[mask]
                weights = weight_loader.get_weights_bulk(chrom_name, positions)
                raw_cn_vals[mask] *= weights
            n_corrected = int(np.sum(raw_cn_vals != chunk['raw_cn'].values))
            n_mult_corrected += n_corrected
            chunk['raw_cn'] = raw_cn_vals

        chunk['bin_start'] = (chunk['start'].values // window_size) * window_size

        # Soft-capping biological filter (A1 + A2 fix)
        # ---------------------------------------------------------------
        # OLD (binary exclusion): k-mers with raw_count > threshold are
        # discarded entirely.  Only below-threshold k-mers contribute to
        # the bin mean → systematic CN underestimation at high-CN loci
        # (A1).  Winsorized mean on top of binary exclusion compounds the
        # downward bias (A2).
        #
        # NEW (soft cap): k-mers above threshold are CAPPED to the
        # threshold value instead of discarded.  All k-mers contribute to
        # the bin mean; extreme values are attenuated to `threshold` rather
        # than removed entirely.  This eliminates the binary-exclusion bias
        # while still limiting bloom-filter artifact inflation.
        #
        # num_filtered continues to count above-threshold k-mers for the
        # downstream quality weight:
        #   quality = num_kmers / (num_kmers + num_filtered + 1)
        # ---------------------------------------------------------------
        if use_per_bin:
            unique_bins = chunk[['chrom', 'bin_start']].drop_duplicates().copy()
            unique_bins['_thr'] = [
                bin_threshold_map.get((r.chrom, int(r.bin_start)), bio_threshold)
                for r in unique_bins.itertuples(index=False)
            ]
            chunk = chunk.merge(unique_bins, on=['chrom', 'bin_start'], how='left')
            thr_vals = chunk['_thr'].fillna(bio_threshold).values
            over_mask = chunk['raw_cn'].values > thr_vals
            # Cap to per-bin threshold
            capped = np.where(over_mask,
                              thr_vals.astype(np.float32),
                              chunk['raw_cn'].values.astype(np.float32))
            chunk.drop(columns=['_thr'], inplace=True)
        elif bio_threshold > 0:
            over_mask = chunk['raw_cn'].values > bio_threshold
            capped = np.where(over_mask,
                              np.float32(bio_threshold),
                              chunk['raw_cn'].values.astype(np.float32))
        else:
            over_mask = np.zeros(len(chunk), dtype=bool)
            capped    = chunk['raw_cn'].values.astype(np.float32)

        n_bio_filtered_total += int(over_mask.sum())

        # Accumulate capped counts per bin (frombytes bulk-append for speed)
        chunk['_capped'] = capped
        for (chrom, bin_start), grp in chunk.groupby(
                ['chrom', 'bin_start'], sort=False)['_capped']:
            buf = acc_values[(chrom, int(bin_start))]
            buf.frombytes(grp.values.astype(np.float32).tobytes())
        chunk.drop(columns=['_capped'], inplace=True)

        # Track above-threshold k-mer count for quality metric (num_filtered)
        if over_mask.any():
            for (chrom, bin_start), n_f in (
                    chunk[over_mask].groupby(['chrom', 'bin_start'],
                                            sort=False)['raw_cn'].count().items()):
                acc_filtered[(chrom, int(bin_start))] += int(n_f)

    skipped_display = sorted(skipped_chroms - MITO_PATTERNS)
    print(f"[IO]   Total lines read:     {n_total:,}")
    print(f"[IO]   Skipped (chrM/etc):   {n_skipped:,}"
          f"{' — ' + str(skipped_display) if skipped_display else ''}")
    if weight_loader is not None:
        pct_mult = 100 * n_mult_corrected / max(n_total, 1)
        print(f"[IO]   Mult-corrected k-mers: {n_mult_corrected:,} ({pct_mult:.2f}%)")
    pct = 100 * n_bio_filtered_total / max(n_total, 1)
    print(f"[IO]   Bio-capped k-mers:    {n_bio_filtered_total:,} ({pct:.2f}%) "
          f"(soft-capped to threshold, not excluded)")
    total_bins = len(acc_values)
    print(f"[IO]   Total bins created:   {total_bins:,}")
    return acc_values, acc_filtered


# ============================================================================
# STEP 2: Compute diploid depth median for normalization
# ============================================================================


def compute_genome_median(acc, sex='XX', neutral_lo=0.4, neutral_hi=2.5, max_iter=5):
    """
    Compute genome-wide neutral-band depth median for CN normalization.

    K1 fix (Phase 2): Iterative neutral-band filtering.
    Problem: Including all bins (CN=1 to CN=150+) inflates the median
    by 3-7%, causing systematic CN underestimation genome-wide.
    Fix: Iteratively restrict to [neutral_lo×, neutral_hi×] × median,
    excluding high-CN duplicated regions from the baseline estimate.

    For XY samples: chrX and chrY are excluded from median computation
    (hemizygous coverage would bias the diploid median downward).

    Algorithm:
      1. Compute initial median from ALL diploid bins
      2. Keep only bins in [neutral_lo × median, neutral_hi × median]
      3. Recompute median from filtered bins
      4. Repeat until convergence (< 0.1% change) or max_iter

    This single value is used as the CN=1.0 baseline in write_output()
    (EK5 fix: genome-wide normalization for cross-chromosome consistency).
    """
    # For XY: exclude sex chromosomes from median (hemizygous depth)
    excluded_chroms = set()
    if sex == 'XY':
        excluded_chroms = CHRX_PATTERNS | CHRY_PATTERNS

    # Collect all winsorized means once
    all_means = []
    seen_chroms = set()
    for (chrom, bin_start), buf in acc.items():
        if chrom in excluded_chroms:
            continue
        if len(buf) < 3:
            continue
        all_means.append(_winsorized_mean(buf))
        seen_chroms.add(chrom)

    all_means = np.array(all_means, dtype=np.float64)
    n_total = len(all_means)

    if n_total == 0:
        print("[NORM] ERROR: No valid bins for genome median computation. "
              "Check input data and --sex setting.")
        sys.exit(1)

    # Initial unfiltered median
    genome_median = float(np.median(all_means))
    n_chroms = len(seen_chroms)
    if excluded_chroms:
        # Count how many actual chromosomes in the data were excluded
        all_data_chroms = set(chrom for (chrom, _) in acc.keys())
        n_actually_excluded = len(all_data_chroms - seen_chroms)
        sex_note = f" (sex={sex}, {n_actually_excluded} chroms excluded)"
    else:
        sex_note = ""
    print(f"[NORM] Initial genome-wide depth median: {genome_median:.4f}  "
          f"({n_total:,} bins across {n_chroms} chroms{sex_note})")

    # Iterative neutral-band filtering
    for i in range(max_iter):
        lo = neutral_lo * genome_median
        hi = neutral_hi * genome_median
        mask = (all_means >= lo) & (all_means <= hi)
        n_kept = int(mask.sum())
        new_median = float(np.median(all_means[mask]))
        change_pct = abs(new_median - genome_median) / genome_median * 100
        n_excluded = n_total - n_kept
        print(f"[NORM] Iter {i+1}: band [{lo:.2f}, {hi:.2f}] → "
              f"{n_kept:,} bins kept ({n_excluded:,} excluded, "
              f"{100*n_excluded/n_total:.1f}%), "
              f"median {genome_median:.4f} → {new_median:.4f} "
              f"(Δ={change_pct:.3f}%)")
        genome_median = new_median
        if change_pct < 0.1:
            print(f"[NORM] Converged after {i+1} iterations (Δ < 0.1%)")
            break

    print(f"[NORM] Final genome-wide depth median: {genome_median:.4f}")
    return genome_median


# ============================================================================
# STEP 3: Write 8-column output BED
# ============================================================================

def write_output(acc_values, acc_filtered, genome_median,
                 window_size, output_path):
    """
    Compute genome-wide normalized CN, log_ratio, and write 8-column BED.
    Output sorted by natural chromosome order then position.

    mean_count   = winsorized mean (p5–p95) of k-mer counts per bin (capped).
    cn           = mean_count / genome_median  (EK5 fix: genome-wide baseline)
    num_filtered = k-mers soft-capped by biological threshold.

    EK5 rationale: per-chromosome normalization causes the same element on
    different chromosomes to receive different CN values (SD-rich chr22 has
    inflated median → genuine CN=1 bins appear CN<1). Using the genome-wide
    neutral-band median eliminates this cross-chromosome inconsistency.
    Post-segmentation GC calibration (apply_gc_cn_calibration) handles
    residual per-chromosome GC bias without distorting cross-chrom CN.
    """
    rows = []

    for (chrom, bin_start), buf in acc_values.items():
        n_win = len(buf)
        if n_win < 1:
            continue
        mean_cnt = _winsorized_mean(buf)

        # EK5 fix: genome-wide median as baseline (cross-chrom consistent CN)
        cn = mean_cnt / genome_median if genome_median > 0 else 1.0
        cn = max(cn, EPSILON)

        log_ratio    = math.log2(cn)
        num_kmers    = n_win
        num_filtered = acc_filtered.get((chrom, bin_start), 0)

        rows.append((
            chrom,
            bin_start,
            bin_start + window_size,
            round(cn, 4),
            round(mean_cnt, 4),
            round(log_ratio, 6),
            num_kmers,
            num_filtered,
        ))

    rows.sort(key=lambda r: (natural_chrom_key(r[0]), r[1]))

    n_written = 0
    with open(output_path, 'w') as fh:
        fh.write(f"# KonuSeg preprocessed windows\n")
        fh.write(f"# window_size={window_size} | normalization=genome_wide\n")
        fh.write(f"# genome_depth_median={genome_median:.4f}\n")
        fh.write("# chrom\tstart\tend\tcn\tmean_count\tlog_ratio\tnum_kmers\tnum_filtered\n")
        for row in rows:
            fh.write('\t'.join(map(str, row)) + '\n')
            n_written += 1

    print(f"[IO] Written {n_written:,} windows → {output_path}")
    return n_written


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Preprocess k-mer BED → KonuSeg window BED"
    )
    parser.add_argument('--input', '-i', required=True,
                        help='K-mer BED file (4-col: chrom, start, end, count)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output 8-col window BED')
    parser.add_argument('--window-size', type=int, default=500,
                        help='Output bin size in bp (default: 500)')
    parser.add_argument('--sex', choices=['XX', 'XY'], default='XX',
                        help='Sample sex karyotype. XX: all chroms diploid '
                             '(default). XY: chrX/chrY excluded from genome-wide '
                             'median computation (hemizygous coverage).')
    parser.add_argument('--bio-threshold-factor', type=float, default=150.0,
                        help='Biological filter: k-mers with raw_count > factor × '
                             'Gaussian_peak are soft-capped. '
                             'Removes LINE/SINE/Satellite k-mer inflation. '
                             'Default: 150. Set to 0 to disable.')
    parser.add_argument('--rm-annotation-bed', default=None,
                        help='Repeat-annotated window BED produced by '
                             'compute_repeat_annotation.py (10 columns, repeat_class '
                             'in column 10). When provided, applies per-bin '
                             'RM-guided thresholds: Satellite=30×, LTR=100×, '
                             'default=bio-threshold-factor×. '
                             'Suppresses Alu/L1 inflation in segmental dup regions.')
    parser.add_argument('--hist-sample-chunks', type=int, default=10,
                        help='Number of 5M-line chunks to sample for histogram '
                             'peak detection (default: 10 = 50M lines). '
                             'Skipped when --bio-threshold-factor 0.')
    parser.add_argument('--weight-dir', default=None,
                        help='Per-chromosome multiplicity weight arrays (.npy) from '
                             'compute_multiplicity_weights.py. Applies per-position '
                             'correction: corrected = raw × weight[pos], where '
                             'weight = M_ref / M_target. Removes repeat element '
                             'k-mer inflation while preserving SD signal. '
                             'Backward compatible (no-op when not set).')
    args = parser.parse_args()

    use_rm = args.rm_annotation_bed is not None

    # Load multiplicity weights (before histogram, so peak reflects corrected values)
    weight_loader = None
    if args.weight_dir is not None:
        if not os.path.isdir(args.weight_dir):
            print(f"ERROR: --weight-dir not found: {args.weight_dir}", file=sys.stderr)
            sys.exit(1)
        weight_loader = WeightLoader(args.weight_dir)

    sex_desc = "all chroms diploid" if args.sex == 'XX' else "chrX/chrY excluded from median"
    print("=" * 60)
    print("KonuSeg — K-mer BED Preprocessing")
    print(f"  Sliding windows → {args.window_size}bp bins")
    print(f"  Sex: {args.sex} ({sex_desc})")
    print("  Normalization: genome-wide neutral-band median")
    if weight_loader is not None:
        print(f"  Multiplicity correction: {args.weight_dir}")
        print(f"    k_target={weight_loader.k_target}, k_ref={weight_loader.k_ref}")
    if args.bio_threshold_factor > 0:
        if use_rm:
            print(f"  Biological filter: RM-guided per-bin thresholds "
                  f"(Sat=30×, default={args.bio_threshold_factor}×)")
        else:
            print(f"  Biological filter: {args.bio_threshold_factor}× Gaussian peak "
                  "(global — use --rm-annotation-bed for per-bin thresholds)")
    print("=" * 60)
    print(f"Input:       {args.input}")
    print(f"Output:      {args.output}")
    print(f"Window size: {args.window_size}bp")
    if use_rm:
        print(f"RM annotation: {args.rm_annotation_bed}")
    if weight_loader is not None:
        print(f"Weight dir:  {args.weight_dir}")
    print()

    # Step 0: Build histogram sample → find Gaussian peak → set bio_threshold
    bio_threshold = 0.0
    gaussian_peak = 0
    if args.bio_threshold_factor > 0:
        print("[HIST] Sampling k-mer count histogram for Gaussian peak detection...")
        histogram = build_histogram_sample(args.input, args.hist_sample_chunks,
                                              weight_loader=weight_loader)
        gaussian_peak = find_gaussian_peak(histogram)
        bio_threshold = gaussian_peak * args.bio_threshold_factor
        print(f"[HIST] Biological threshold: {gaussian_peak} × "
              f"{args.bio_threshold_factor} = {bio_threshold:.1f}")
        print()

    # Step 0a (optional): Build per-bin threshold map from RM annotation
    bin_threshold_map = None
    if use_rm and gaussian_peak > 0:
        import os as _os
        if not _os.path.exists(args.rm_annotation_bed):
            print(f"[RM] WARNING: --rm-annotation-bed not found: {args.rm_annotation_bed}")
            print(f"[RM]          Falling back to global threshold {bio_threshold:.1f}")
        else:
            rm_lookup = load_rm_lookup(args.rm_annotation_bed)
            bin_threshold_map = build_bin_threshold_map(
                rm_lookup, gaussian_peak, args.bio_threshold_factor)
            print()

    # Step 1: Aggregate sliding-window k-mers into bins (with bio filter + winsorized mean)
    acc_values, acc_filtered = aggregate_windows(
        args.input, args.window_size,
        bio_threshold=bio_threshold,
        bin_threshold_map=bin_threshold_map,
        weight_loader=weight_loader)

    # Step 2: Compute genome-wide normalization factor
    print()
    genome_median = compute_genome_median(acc_values, sex=args.sex)

    # Step 3: Write normalized output
    print()
    n = write_output(acc_values, acc_filtered, genome_median,
                     args.window_size, args.output)

    print()
    print(f"Done. {n:,} windows written.")
    return 0


if __name__ == '__main__':
    sys.exit(main())

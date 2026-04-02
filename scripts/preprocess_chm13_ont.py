#!/usr/bin/env python3
"""
preprocess_chm13_ont.py — CHM13 ONT k-mer BED → KonuSeg 500bp window BED
=========================================================================

Converts CHM13 ONT k-mer depth data (chr prefix, 72bp sliding windows)
into the 8-column 500bp window format used by the KonuSeg HMM.

Input format (4 columns, chr prefix, 1bp-step sliding windows):
  chrN  start  end  kmer_count

Output format (8 columns, CM039 accessions):
  chrom  start  end  cn  mean_count  log_ratio  num_kmers  num_filtered

Normalization:
  mean_count = mean raw k-mer count in 500bp bin
  cn = mean_count / diploid_depth_median   (Neutral regions → CN ≈ 1.0)
  log_ratio  = log2(cn + epsilon)
  num_kmers  = number of contributing 72bp windows (coverage proxy)
  num_filtered = 0  (no bloom filter data for CHM13 ONT)

Notes:
  - CHM13 is female (XX) → chrX is diploid; included in normalization
  - chrM excluded (mitochondrial CN is thousands, distorts normalization)
  - Memory-efficient: processes in 5M-line chunks, aggregates per bin

Usage:
  python3 scripts/preprocess_chm13_ont.py \\
    --input  /home/klea/recomb/gra_bf/output/final_algo/chm13/\\
             chm13_ont_cpn_merged_algo_k72_minseg500.bed_noNoiseFiltering_bloom_filter_kmers.bed \\
    --output input/chm13_ont_cn_w500.bed

Author: KonuSeg Team
Date: March 2026
"""

import argparse
import array as _array
import math
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

# ============================================================================
# T2T-CHM13v2.0: chr prefix → CM039 GenBank accessions
# ============================================================================

CHR_TO_CM039 = {
    'chr1':  'CM039011.1',
    'chr2':  'CM039012.1',
    'chr3':  'CM039013.1',
    'chr4':  'CM039014.1',
    'chr5':  'CM039015.1',
    'chr6':  'CM039016.1',
    'chr7':  'CM039017.1',
    'chr8':  'CM039018.1',
    'chr9':  'CM039019.1',
    'chr10': 'CM039020.1',
    'chr11': 'CM039021.1',
    'chr12': 'CM039022.1',
    'chr13': 'CM039023.1',
    'chr14': 'CM039024.1',
    'chr15': 'CM039025.1',
    'chr16': 'CM039026.1',
    'chr17': 'CM039027.1',
    'chr18': 'CM039028.1',
    'chr19': 'CM039029.1',
    'chr20': 'CM039030.1',
    'chr21': 'CM039031.1',
    'chr22': 'CM039032.1',
    'chrX':  'CM039033.1',
}

# CM039 chromosome order for output sorting
CHROM_ORDER = list(CHR_TO_CM039.values())

# CHM13 is female (XX): all 23 chromosomes are diploid
DIPLOID_CHROMS = set(CHR_TO_CM039.values())

EPSILON = 1e-3

# ============================================================================
# NOR LOCUS COORDINATES (T2T-CHM13v2.0, chr prefix)
# ============================================================================
# 45S rDNA arrays sit on the five acrocentric NOR chromosomes.  All NOR loci
# share nearly-identical 45S units → their k-mers pool across chromosomes and
# exceed the global bio-filter threshold (>150× single-copy).  Applying the
# normal threshold here destroys rDNA signal.  Instead we set threshold = ∞
# so every rDNA k-mer is counted, regardless of its raw depth.
#
# Coordinates include ±500 kb flanking buffer around the annotated NOR array.
NOR_REGIONS = [
    ('chr13',  5_282_236,  9_808_776),   # chr13 NOR ± 500 kb
    ('chr14',  1_651_445,  3_317_142),   # chr14 NOR ± 500 kb
    ('chr15',  2_056_360,  5_206_717),   # chr15 NOR ± 500 kb
    ('chr21',  2_750_105,  6_072_749),   # chr21 NOR ± 500 kb
    ('chr22',  4_305_000,  6_220_000),   # chr22 NOR ± 500 kb
]

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
#   NOR regions  ∞   rDNA arrays; must not be filtered (see NOR_REGIONS above)
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
# STEP 0a: NOR bin set + RM annotation helpers
# ============================================================================

def build_nor_bins_set(window_size: int = 500) -> set:
    """
    Return the set of (cm039_chrom, bin_start) tuples covering all NOR loci
    (including ±500 kb flanking buffer defined in NOR_REGIONS).
    These bins will use threshold = ∞ (no bio-filter).
    """
    nor_bins: set = set()
    for (chrom_chr, start, end) in NOR_REGIONS:
        cm039 = CHR_TO_CM039.get(chrom_chr)
        if cm039 is None:
            continue
        first_bin = (start // window_size) * window_size
        last_bin  = (end   // window_size) * window_size
        for b in range(first_bin, last_bin + window_size, window_size):
            nor_bins.add((cm039, b))
    return nor_bins


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


def build_bin_threshold_map(rm_lookup: dict, nor_bins_set: set,
                             peak: float, global_factor: float) -> dict:
    """
    Build a {(chrom, bin_start): threshold} mapping.

    Priority:
      1. NOR bins           → threshold = ∞  (no filter — rDNA must not be clipped)
      2. RM-annotated bins  → THRESHOLD_FACTORS[repeat_class] × peak
      3. Unannotated bins   → global_factor × peak  (fallback)

    Only NOR and RM-annotated bins are stored explicitly; everything else uses
    the fallback (global_factor × peak) at query time.
    """
    threshold_map: dict = {}

    # 1. NOR bins → no filter
    for key in nor_bins_set:
        threshold_map[key] = float('inf')

    # 2. RM-annotated bins → class-specific factor (NOR takes priority)
    for (chrom, bin_start), repeat_class in rm_lookup.items():
        if (chrom, bin_start) in nor_bins_set:
            continue
        factor = THRESHOLD_FACTORS.get(repeat_class, global_factor)
        threshold_map[(chrom, bin_start)] = peak * factor

    n_nor  = sum(1 for v in threshold_map.values() if v == float('inf'))
    n_line = sum(1 for k, v in threshold_map.items()
                 if v != float('inf') and v < peak * global_factor)
    print(f"[RM] Threshold map built: "
          f"{n_nor:,} NOR bins (∞), "
          f"{n_line:,} reduced-threshold bins, "
          f"rest → {global_factor:.0f}× peak ({peak * global_factor:.1f})")
    return threshold_map


# ============================================================================
# STEP 0: Histogram sampling → Gaussian peak → biological threshold
# ============================================================================

def build_histogram_sample(filepath, sample_chunks=10, chunk_size=5_000_000):
    """
    Read the first sample_chunks chunks to build a k-mer count histogram.

    A 50M-line sample (default: 10 × 5M) is sufficient to identify the
    Gaussian peak (= single-copy sequencing depth) without reading the full
    2.5B-line file.

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

    reader = pd.read_csv(
        filepath,
        sep='\t',
        header=None,
        comment='#',
        usecols=[3],
        names=['raw_cn'],
        dtype={'raw_cn': np.float32},
        chunksize=chunk_size,
    )

    n_chunks = 0
    for chunk in reader:
        counts = np.clip(chunk['raw_cn'].values, 0, 65535).astype(np.int32)
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

def _winsorized_mean(values_buf, p_low=5.0, p_high=95.0, use_raw_mean=False):
    """
    Compute winsorized mean from an array.array('f') buffer.
    Clips values outside [p_low, p_high] percentile range, then averages.
    Falls back to plain mean when n < 3 (percentiles undefined).

    use_raw_mean: if True, return plain mean without percentile clipping.
    Used for NOR bins where right-skewed rDNA k-mer counts should not be
    truncated (p95 clipping causes 19-36% systematic underestimation).
    """
    arr = np.frombuffer(values_buf, dtype=np.float32)
    if len(arr) < 3:
        return float(arr.mean()) if len(arr) > 0 else 0.0
    if use_raw_mean:
        return float(arr.mean())
    lo = float(np.percentile(arr, p_low))
    hi = float(np.percentile(arr, p_high))
    clipped = arr[(arr >= lo) & (arr <= hi)]
    return float(clipped.mean()) if len(clipped) > 0 else float(arr.mean())


def aggregate_windows(filepath, window_size, bio_threshold=0,
                      bin_threshold_map=None):
    """
    Read 4-column sliding-window BED in 5M-line chunks.

    For each 500bp bin:
      - k-mers with raw_count > threshold are removed (biological filter).
      - Threshold is either global (bio_threshold) or per-bin (bin_threshold_map).
      - Valid k-mer counts are stored in a compact array.array('f') buffer.
      - Winsorized mean (5th–95th percentile) is computed at write-time.

    Per-bin threshold mode (bin_threshold_map):
      - NOR locus bins   → threshold = ∞  (rDNA k-mers kept regardless of count)
      - SINE/LINE bins   → threshold = 50 × peak  (Alu/L1 inflation suppressed)
      - Satellite bins   → threshold = 30 × peak  (centromeric arrays suppressed)
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

    use_per_bin = bin_threshold_map is not None

    print(f"[IO] Reading: {filepath}")
    print(f"[IO]   Chunk size: {CHUNK_SIZE:,} lines | bin size: {window_size}bp")
    print(f"[IO]   Aggregation: winsorized mean (p5–p95)")
    if use_per_bin:
        sat_f  = THRESHOLD_FACTORS.get('Satellite',  30.0)
        line_f = THRESHOLD_FACTORS.get('LINE',       150.0)
        print(f"[IO]   Biological filter: per-bin RM-guided soft-cap "
              f"(NOR=∞, LINE/SINE={line_f:.0f}×peak, Satellite={sat_f:.0f}×peak, "
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

        # Filter to known chromosomes
        known_mask = chunk['chrom'].isin(CHR_TO_CM039)
        n_skipped += int((~known_mask).sum())
        skipped_chroms.update(chunk.loc[~known_mask, 'chrom'].unique())
        chunk = chunk[known_mask].copy()
        if chunk.empty:
            continue
        n_total += len(chunk)

        chunk['chrom']     = chunk['chrom'].map(CHR_TO_CM039)
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
        # NOR bins (threshold=∞): over_mask is always False → no capping,
        # full rDNA signal preserved.
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
            # Cap to per-bin threshold (inf stays inf for NOR → never triggered)
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

    skipped_display = sorted(skipped_chroms - {'chrM'})
    print(f"[IO]   Total lines read:     {n_total:,}")
    print(f"[IO]   Skipped (chrM/etc):   {n_skipped:,}"
          f"{' — ' + str(skipped_display) if skipped_display else ''}")
    pct = 100 * n_bio_filtered_total / max(n_total, 1)
    print(f"[IO]   Bio-capped k-mers:    {n_bio_filtered_total:,} ({pct:.2f}%) "
          f"(soft-capped to threshold, not excluded)")
    total_bins = len(acc_values)
    print(f"[IO]   Total bins created:   {total_bins:,}")
    return acc_values, acc_filtered


# ============================================================================
# STEP 2: Compute diploid depth median for normalization
# ============================================================================


def compute_genome_median(acc, neutral_lo=0.4, neutral_hi=2.5, max_iter=5):
    """
    Compute genome-wide neutral-band depth median for CN normalization.

    K1 fix (Phase 2): Iterative neutral-band filtering.
    Problem: Including all bins (CN=1 to CN=150+) inflates the median
    by 3-7%, causing systematic CN underestimation genome-wide.
    Fix: Iteratively restrict to [neutral_lo×, neutral_hi×] × median,
    excluding high-CN duplicated regions from the baseline estimate.

    Algorithm:
      1. Compute initial median from ALL diploid bins
      2. Keep only bins in [neutral_lo × median, neutral_hi × median]
      3. Recompute median from filtered bins
      4. Repeat until convergence (< 0.1% change) or max_iter

    This single value is used as the CN=1.0 baseline in write_output()
    (EK5 fix: genome-wide normalization for cross-chromosome consistency).
    """
    # Collect all winsorized means once
    all_means = []
    for (chrom, bin_start), buf in acc.items():
        if chrom not in DIPLOID_CHROMS:
            continue
        if len(buf) < 3:
            continue
        all_means.append(_winsorized_mean(buf))

    all_means = np.array(all_means, dtype=np.float64)
    n_total = len(all_means)

    # Initial unfiltered median
    genome_median = float(np.median(all_means))
    print(f"[NORM] Initial genome-wide depth median: {genome_median:.4f}  "
          f"({n_total:,} bins across {len(CHROM_ORDER)} chroms)")

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
    Output sorted by chromosome order (chr1→chrX) then position.

    mean_count   = winsorized mean (p5–p95) of k-mer counts per bin (capped).
    cn           = mean_count / genome_median  (EK5 fix: genome-wide baseline)
    num_filtered = k-mers soft-capped by biological threshold.

    NOR bins use raw mean (no percentile clipping) to avoid systematic
    underestimation of rDNA CN caused by p95 truncation of right-skewed counts.

    EK5 rationale: per-chromosome normalization causes the same element on
    different chromosomes to receive different CN values (SD-rich chr22 has
    inflated median → genuine CN=1 bins appear CN<1). Using the genome-wide
    neutral-band median eliminates this cross-chromosome inconsistency.
    Post-segmentation GC calibration (apply_gc_cn_calibration) handles
    residual per-chromosome GC bias without distorting cross-chrom CN.
    """
    chrom_rank = {c: i for i, c in enumerate(CHROM_ORDER)}
    nor_bins = build_nor_bins_set(window_size)
    rows = []
    n_nor = 0

    for (chrom, bin_start), buf in acc_values.items():
        n_win = len(buf)
        if n_win < 1:
            continue
        is_nor = (chrom, bin_start) in nor_bins
        mean_cnt = _winsorized_mean(buf, use_raw_mean=is_nor)
        if is_nor:
            n_nor += 1

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

    rows.sort(key=lambda r: (chrom_rank.get(r[0], 99), r[1]))

    n_written = 0
    with open(output_path, 'w') as fh:
        fh.write(f"# KonuSeg CHM13 ONT preprocessed windows\n")
        fh.write(f"# window_size={window_size} | normalization=genome_wide\n")
        fh.write(f"# genome_depth_median={genome_median:.4f}\n")
        fh.write("# chrom\tstart\tend\tcn\tmean_count\tlog_ratio\tnum_kmers\tnum_filtered\n")
        for row in rows:
            fh.write('\t'.join(map(str, row)) + '\n')
            n_written += 1

    if n_nor:
        print(f"[NOR] {n_nor:,} NOR bins used raw mean (no winsorized clipping)")
    print(f"[IO] Written {n_written:,} windows → {output_path}")
    return n_written


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Preprocess CHM13 ONT k-mer BED → KonuSeg 500bp window BED"
    )
    parser.add_argument('--input', '-i', required=True,
                        help='CHM13 ONT k-mer BED (4-col, chr prefix, 72bp sliding windows)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output 8-col 500bp window BED (CM039 accessions)')
    parser.add_argument('--window-size', type=int, default=500,
                        help='Output bin size in bp (default: 500)')
    parser.add_argument('--bio-threshold-factor', type=float, default=150.0,
                        help='Biological filter: k-mers with raw_count > factor × '
                             'Gaussian_peak are excluded before window mean computation. '
                             'Removes LINE/SINE/Satellite k-mer inflation. '
                             'Default: 150. Set to 0 to disable.')
    parser.add_argument('--rm-annotation-bed', default=None,
                        help='Repeat-annotated 500bp window BED produced by '
                             'compute_repeat_annotation.py (10 columns, repeat_class '
                             'in column 10).  When provided, applies per-bin '
                             'RM-guided thresholds: NOR=∞, SINE/LINE=50×, '
                             'Satellite=30×, LTR=100×, default=bio-threshold-factor×. '
                             'Dramatically improves rDNA accuracy and suppresses '
                             'Alu/L1 inflation in segmental dup regions (e.g. SMN).')
    parser.add_argument('--hist-sample-chunks', type=int, default=10,
                        help='Number of 5M-line chunks to sample for histogram '
                             'peak detection (default: 10 = 50M lines). '
                             'Skipped when --bio-threshold-factor 0.')
    args = parser.parse_args()

    use_rm = args.rm_annotation_bed is not None

    print("=" * 60)
    print("KonuSeg — CHM13 ONT Preprocessing")
    print("  chr prefix → CM039 accessions")
    print("  72bp sliding windows → 500bp bins")
    print("  CHM13 is female: chrX included in diploid normalization")
    print("  Normalization: per-chromosome median")
    if args.bio_threshold_factor > 0:
        if use_rm:
            print(f"  Biological filter: RM-guided per-bin thresholds "
                  f"(NOR=∞, SINE/LINE=50×, Sat=30×, default={args.bio_threshold_factor}×)")
        else:
            print(f"  Biological filter: {args.bio_threshold_factor}× Gaussian peak "
                  "(global — use --rm-annotation-bed for per-bin thresholds)")
    print("=" * 60)
    print(f"Input:       {args.input}")
    print(f"Output:      {args.output}")
    print(f"Window size: {args.window_size}bp")
    if use_rm:
        print(f"RM annotation: {args.rm_annotation_bed}")
    print()

    # Step 0: Build histogram sample → find Gaussian peak → set bio_threshold
    bio_threshold = 0.0
    gaussian_peak = 0
    if args.bio_threshold_factor > 0:
        print("[HIST] Sampling k-mer count histogram for Gaussian peak detection...")
        histogram = build_histogram_sample(args.input, args.hist_sample_chunks)
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
            print("[RM] Building NOR bin set (rDNA threshold = ∞)...")
            nor_bins = build_nor_bins_set(args.window_size)
            print(f"[RM] NOR bins: {len(nor_bins):,} "
                  f"({len(nor_bins) * args.window_size / 1e6:.1f} Mb)")

            rm_lookup = load_rm_lookup(args.rm_annotation_bed)
            bin_threshold_map = build_bin_threshold_map(
                rm_lookup, nor_bins, gaussian_peak, args.bio_threshold_factor)
            print()

    # Step 1: Aggregate all 72bp windows into 500bp bins (with bio filter + winsorized mean)
    acc_values, acc_filtered = aggregate_windows(
        args.input, args.window_size,
        bio_threshold=bio_threshold,
        bin_threshold_map=bin_threshold_map)

    # Step 2: Compute genome-wide normalization factor
    print()
    genome_median = compute_genome_median(acc_values)

    # Step 3: Write normalized output
    print()
    n = write_output(acc_values, acc_filtered, genome_median,
                     args.window_size, args.output)

    print()
    print(f"Done. {n:,} windows written.")
    return 0


if __name__ == '__main__':
    sys.exit(main())

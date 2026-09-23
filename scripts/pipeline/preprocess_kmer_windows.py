#!/usr/bin/env python3

import argparse
import array as _array
import json
import math
import os
import sys
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd

import os, sys
_HERE = os.path.abspath(os.path.dirname(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from chrom_utils import (
    MITO_PATTERNS, CHRX_PATTERNS, CHRY_PATTERNS,
    is_mito, natural_chrom_key, resolve_chrom,
)

EPSILON = 1e-3

THRESHOLD_FACTORS = {
    'Satellite':       30.0,
    'LINE':           150.0,
    'SINE':           150.0,
    'Simple_repeat':  150.0,
    'LTR':            100.0,
    'DNA':            150.0,
    'Low_complexity': 150.0,
    'Other':          150.0,
    'None':           150.0,
}



class WeightLoader:

    def __init__(self, weight_dir, max_cached=3):
        self.weight_dir = weight_dir
        self.max_cached = max_cached
        self._cache = OrderedDict()
        self._chrom_map = {}

        manifest_path = os.path.join(weight_dir, 'manifest.json')
        if not os.path.exists(manifest_path):
            raise FileNotFoundError(f"manifest.json not found in {weight_dir}")
        with open(manifest_path) as f:
            self._manifest = json.load(f)

        self.k_target = self._manifest.get('k_target', '?')
        self.k_ref = self._manifest.get('k_ref', '?')
        self._weight_chroms = set(self._manifest.get('chromosomes', {}).keys())

    def _resolve(self, chrom):
        if chrom in self._chrom_map:
            return self._chrom_map[chrom]
        resolved = resolve_chrom(chrom, self._weight_chroms)
        self._chrom_map[chrom] = resolved
        return resolved

    def _ensure_loaded(self, chrom):
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
        arr = self._ensure_loaded(chrom)
        weights = np.ones(len(positions), dtype=np.float32)
        if arr is None:
            return weights
        valid = (positions >= 0) & (positions < len(arr))
        weights[valid] = arr[positions[valid]].astype(np.float32)
        return weights



REPEAT_KMER_WEIGHT = 0.01

RM_CLASS_NAME_TO_ID = {
    'Satellite': 1, 'Simple_repeat': 2, 'LINE': 3, 'SINE': 4,
    'LTR': 5, 'DNA': 6, 'Low_complexity': 7, 'Other': 8,
}

DEFAULT_CLASS_WEIGHTS_K32 = {
    1: 0.001,
    2: 0.30,
    3: 0.01,
    4: 0.001,
    5: 0.05,
    6: 0.10,
    7: 0.30,
    8: 0.05,
}


class RMMaskLoader:

    def __init__(self, mask_dir, repeat_weight=REPEAT_KMER_WEIGHT,
                 class_weights=None, max_cached=3):
        self.mask_dir = mask_dir
        self.repeat_weight = repeat_weight
        self.max_cached = max_cached
        self._cache = OrderedDict()
        self._chrom_map = {}

        if class_weights is not None:
            self._weight_lut = np.ones(9, dtype=np.float32)
            for cid, w in class_weights.items():
                if 1 <= cid <= 8:
                    self._weight_lut[cid] = float(w)
            self.repeat_weight = None
        else:
            self._weight_lut = None

        manifest_path = os.path.join(mask_dir, 'manifest.json')
        if not os.path.exists(manifest_path):
            raise FileNotFoundError(f"manifest.json not found in {mask_dir}")
        with open(manifest_path) as f:
            self._manifest = json.load(f)

        self._mask_chroms = set(self._manifest.get('chromosomes', {}).keys())

    def _resolve(self, chrom):
        if chrom in self._chrom_map:
            return self._chrom_map[chrom]
        resolved = resolve_chrom(chrom, self._mask_chroms)
        self._chrom_map[chrom] = resolved
        return resolved

    def _ensure_loaded(self, chrom):
        resolved = self._resolve(chrom)
        if resolved is None:
            return None
        if resolved in self._cache:
            self._cache.move_to_end(resolved)
            return self._cache[resolved]
        npy_path = os.path.join(self.mask_dir, f"{resolved}.npy")
        if not os.path.exists(npy_path):
            return None
        arr = np.load(npy_path)
        self._cache[resolved] = arr
        self._cache.move_to_end(resolved)
        if len(self._cache) > self.max_cached:
            self._cache.popitem(last=False)
        return arr

    def get_weights_bulk(self, chrom, positions):
        arr = self._ensure_loaded(chrom)
        weights = np.ones(len(positions), dtype=np.float32)
        if arr is None:
            return weights
        valid = (positions >= 0) & (positions < len(arr))
        mask_vals = arr[positions[valid]]

        if self._weight_lut is not None:
            weights[valid] = self._weight_lut[mask_vals]
        else:
            is_repeat = mask_vals > 0
            w = np.ones(int(valid.sum()), dtype=np.float32)
            w[is_repeat] = self.repeat_weight
            weights[valid] = w
        return weights

    def count_repeat_positions(self, chrom, positions):
        arr = self._ensure_loaded(chrom)
        if arr is None:
            return len(positions), 0
        valid = (positions >= 0) & (positions < len(arr))
        mask_vals = arr[positions[valid]]
        n_repeat = int((mask_vals > 0).sum())
        n_unique = int(valid.sum()) - n_repeat
        return n_unique, n_repeat

    def get_unique_mask_bulk(self, chrom, positions):
        arr = self._ensure_loaded(chrom)
        is_unique = np.zeros(len(positions), dtype=bool)
        if arr is None:
            return is_unique
        valid = (positions >= 0) & (positions < len(arr))
        is_unique[valid] = (arr[positions[valid]] == 0)
        return is_unique

    def get_weights_and_unique_bulk(self, chrom, positions):
        arr = self._ensure_loaded(chrom)
        weights   = np.ones(len(positions), dtype=np.float32)
        is_unique = np.zeros(len(positions), dtype=bool)
        if arr is None:
            return weights, is_unique
        valid = (positions >= 0) & (positions < len(arr))
        mask_vals = arr[positions[valid]]
        if self._weight_lut is not None:
            weights[valid] = self._weight_lut[mask_vals]
        else:
            is_repeat = mask_vals > 0
            w = np.ones(int(valid.sum()), dtype=np.float32)
            w[is_repeat] = self.repeat_weight
            weights[valid] = w
        is_unique[valid] = (mask_vals == 0)
        return weights, is_unique



def load_rm_lookup(bed_path: str) -> dict:
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



def build_histogram_sample(filepath, sample_chunks=10, chunk_size=5_000_000,
                           weight_loader=None, per_window_correct=False,
                           window_size=500, pw_percentile=25.0,
                           rm_mask_loader=None):
    hist = np.zeros(65536, dtype=np.int64)

    if weight_loader is not None or per_window_correct or rm_mask_loader is not None:
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

        if weight_loader is not None:
            for chrom_name in chunk['chrom'].unique():
                mask = chunk['chrom'].values == chrom_name
                positions = chunk['start'].values[mask]
                weights = weight_loader.get_weights_bulk(chrom_name, positions)
                counts[mask] *= weights

        if rm_mask_loader is not None:
            for chrom_name in chunk['chrom'].unique():
                mask = chunk['chrom'].values == chrom_name
                positions = chunk['start'].values[mask]
                rm_weights = rm_mask_loader.get_weights_bulk(chrom_name, positions)
                counts[mask] *= rm_weights

        if per_window_correct:
            bin_starts = (chunk['start'].values // window_size) * window_size
            df_tmp = pd.DataFrame({
                'chrom': chunk['chrom'].values, 'bin': bin_starts, 'c': counts,
            })
            bin_p5 = (df_tmp.groupby(['chrom', 'bin'], sort=False)['c']
                      .quantile(pw_percentile / 100.0).clip(lower=1.0)
                      .reset_index().rename(columns={'c': '_p5'}))
            df_tmp = df_tmp.merge(bin_p5, on=['chrom', 'bin'], how='left')
            p5 = df_tmp['_p5'].values.astype(np.float32)
            c_f = counts.astype(np.float32)
            counts = c_f * (p5 / np.maximum(c_f, p5))

        if rm_mask_loader is not None and rm_mask_loader.repeat_weight == 0.0:
            counts = counts[counts > 0]
        counts = np.clip(counts, 0, 65535).astype(np.int32)
        np.add.at(hist, counts, 1)
        n_chunks += 1
        if n_chunks >= sample_chunks:
            break

    print(f"[HIST] Histogram built from first {n_chunks} chunks "
          f"({n_chunks * chunk_size / 1e6:.0f}M lines sampled)")
    return hist


def find_gaussian_peak(histogram, search_min=2, search_max=500):
    h = histogram

    peaks = []
    for i in range(max(search_min, 2), min(search_max, len(h)) - 1):
        if h[i] > h[i - 1] and h[i] >= h[i + 1] and h[i] > 0:
            peaks.append(i)

    start_val = float(np.max(h[search_min:min(search_min + 10, search_max)]))
    noise_threshold = start_val * 0.05
    noise_tail_end = search_min
    for i in range(search_min, min(search_max, len(h))):
        if h[i] < noise_threshold:
            noise_tail_end = i
            break

    bio_peaks = [p for p in peaks if p >= noise_tail_end]

    if bio_peaks:
        peak = max(bio_peaks, key=lambda p: int(h[p]))
        valley = noise_tail_end
        for i in range(noise_tail_end, peak):
            if h[i] <= h[i + 1]:
                valley = i
                break
    elif peaks:
        peak = peaks[-1]
        valley = peak - 1 if peak > 0 else search_min
    else:
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


def fit_gaussian_to_histogram(histogram, p0_peak, search_max=500,
                              fit_window_sigmas=3.0,
                              min_r_squared=0.70):
    from scipy.optimize import curve_fit

    def gaussian(x, A, mu, sigma):
        return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

    h = np.asarray(histogram, dtype=np.float64)
    p0_peak = max(2, int(p0_peak))
    sigma_guess = max(1.5, 1.3 * np.sqrt(p0_peak))
    fit_lo = max(2, int(p0_peak - fit_window_sigmas * sigma_guess))
    fit_hi = min(search_max, int(p0_peak + fit_window_sigmas * sigma_guess))
    if fit_hi - fit_lo < 5:
        return float(p0_peak), float(sigma_guess), 0.0, False

    x = np.arange(fit_lo, fit_hi + 1)
    y = h[fit_lo:fit_hi + 1]
    if y.sum() == 0 or np.all(y == y[0]):
        return float(p0_peak), float(sigma_guess), 0.0, False

    p0 = [float(h[p0_peak]), float(p0_peak), float(sigma_guess)]
    lo_bounds = [0.0,           float(fit_lo), 0.5]
    hi_bounds = [float(h.max()) * 2.0, float(fit_hi), float(sigma_guess) * 3.0]
    try:
        popt, _ = curve_fit(gaussian, x, y, p0=p0,
                            bounds=(lo_bounds, hi_bounds),
                            maxfev=5000)
        A_fit, mu_fit, sigma_fit = popt
        y_pred = gaussian(x, *popt)
        ss_res = float(np.sum((y - y_pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        ok = r_squared >= min_r_squared
        return float(mu_fit), float(sigma_fit), float(r_squared), bool(ok)
    except (RuntimeError, ValueError) as e:
        print(f"[HIST] Gaussian fit failed: {e}. Falling back to p0_peak={p0_peak}")
        return float(p0_peak), float(sigma_guess), 0.0, False



def _winsorized_mean(values_buf, p_low=5.0, p_high=95.0, exclude_zeros=False):
    arr = np.frombuffer(values_buf, dtype=np.float32)
    n_total = len(arr)
    if exclude_zeros:
        arr = arr[arr > 0]
    n_used = len(arr)
    if n_used < 3:
        val = float(arr.mean()) if n_used > 0 else 0.0
        return (val, n_used) if exclude_zeros else val
    lo = float(np.percentile(arr, p_low))
    hi = float(np.percentile(arr, p_high))
    clipped = arr[(arr >= lo) & (arr <= hi)]
    val = float(clipped.mean()) if len(clipped) > 0 else float(arr.mean())
    return (val, n_used) if exclude_zeros else val


def aggregate_windows(filepath, window_size, bio_threshold=0,
                      bin_threshold_map=None, weight_loader=None,
                      per_window_correct=False, pw_percentile=25.0,
                      rm_mask_loader=None):
    CHUNK_SIZE = 5_000_000

    acc_values        = defaultdict(lambda: _array.array('f'))
    acc_filtered      = defaultdict(int)
    acc_unique_values = defaultdict(lambda: _array.array('f'))
    skipped_chroms = set()
    n_total = 0
    n_skipped = 0
    n_bio_filtered_total = 0
    n_mult_corrected = 0
    n_rm_downweighted = 0

    use_per_bin = bin_threshold_map is not None
    use_rm_mask = rm_mask_loader is not None

    print(f"[IO] Reading: {filepath}")
    print(f"[IO]   Chunk size: {CHUNK_SIZE:,} lines | bin size: {window_size}bp")
    print(f"[IO]   Aggregation: winsorized mean (p5–p95)")
    if weight_loader is not None:
        print(f"[IO]   Multiplicity correction: ENABLED "
              f"(k_target={weight_loader.k_target}, k_ref={weight_loader.k_ref})")
    if use_rm_mask:
        if rm_mask_loader.repeat_weight == 0.0:
            print(f"[IO]   RM k-mer EXCLUSION mode: ENABLED "
                  f"(repeat k-mers fully excluded from CN computation)")
        else:
            print(f"[IO]   RM k-mer down-weighting: ENABLED "
                  f"(repeat_weight={rm_mask_loader.repeat_weight})")
    if per_window_correct:
        print(f"[IO]   Per-window multiplicity correction: ENABLED (p{pw_percentile:g} normalization)")
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

        mito_mask = chunk['chrom'].isin(MITO_PATTERNS)
        n_skipped += int(mito_mask.sum())
        skipped_chroms.update(chunk.loc[mito_mask, 'chrom'].unique())
        chunk = chunk[~mito_mask].copy()
        if chunk.empty:
            continue
        n_total += len(chunk)

        if use_rm_mask:
            raw_orig_for_unique = chunk['raw_cn'].values.copy()

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

        if use_rm_mask:
            raw_cn_vals = chunk['raw_cn'].values.copy()
            is_unique_vec = np.zeros(len(chunk), dtype=bool)
            for chrom_name in chunk['chrom'].unique():
                chrom_sel = chunk['chrom'].values == chrom_name
                positions = chunk['start'].values[chrom_sel]
                rm_weights, is_unique_chrom = (
                    rm_mask_loader.get_weights_and_unique_bulk(chrom_name, positions))
                raw_cn_vals[chrom_sel] *= rm_weights
                n_rm_downweighted += int((rm_weights < 1.0).sum())
                is_unique_vec[chrom_sel] = is_unique_chrom
            chunk['raw_cn'] = raw_cn_vals

        if per_window_correct:
            raw_vals = chunk['raw_cn'].values.copy().astype(np.float32)
            bin_p5 = (chunk.groupby(['chrom', 'bin_start'], sort=False)['raw_cn']
                      .quantile(pw_percentile / 100.0).clip(lower=1.0)
                      .reset_index().rename(columns={'raw_cn': '_p5'}))
            chunk = chunk.merge(bin_p5, on=['chrom', 'bin_start'], how='left')
            p5_vals = chunk['_p5'].values.astype(np.float32)
            weights = p5_vals / np.maximum(raw_vals, p5_vals)
            chunk['raw_cn'] = raw_vals * weights
            chunk.drop(columns=['_p5'], inplace=True)
            n_mult_corrected += int(np.sum(weights < 0.999))

        if use_per_bin:
            unique_bins = chunk[['chrom', 'bin_start']].drop_duplicates().copy()
            unique_bins['_thr'] = [
                bin_threshold_map.get((r.chrom, int(r.bin_start)), bio_threshold)
                for r in unique_bins.itertuples(index=False)
            ]
            chunk = chunk.merge(unique_bins, on=['chrom', 'bin_start'], how='left')
            thr_vals = chunk['_thr'].fillna(bio_threshold).values
            over_mask = chunk['raw_cn'].values > thr_vals
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

        if use_rm_mask:
            if use_per_bin:
                over_mask_u = raw_orig_for_unique > thr_vals
                capped_unique = np.where(
                    over_mask_u,
                    thr_vals.astype(np.float32),
                    raw_orig_for_unique.astype(np.float32))
            elif bio_threshold > 0:
                over_mask_u = raw_orig_for_unique > bio_threshold
                capped_unique = np.where(
                    over_mask_u,
                    np.float32(bio_threshold),
                    raw_orig_for_unique.astype(np.float32))
            else:
                capped_unique = raw_orig_for_unique.astype(np.float32)

            if is_unique_vec.any():
                u_chrom = chunk['chrom'].values[is_unique_vec]
                u_bin   = chunk['bin_start'].values[is_unique_vec]
                u_vals  = capped_unique[is_unique_vec]
                u_df = pd.DataFrame({
                    'chrom':     u_chrom,
                    'bin_start': u_bin,
                    '_uc':       u_vals,
                })
                for (chrom, bin_start), grp in u_df.groupby(
                        ['chrom', 'bin_start'], sort=False)['_uc']:
                    buf = acc_unique_values[(chrom, int(bin_start))]
                    buf.frombytes(grp.values.astype(np.float32).tobytes())

        chunk['_capped'] = capped
        for (chrom, bin_start), grp in chunk.groupby(
                ['chrom', 'bin_start'], sort=False)['_capped']:
            buf = acc_values[(chrom, int(bin_start))]
            buf.frombytes(grp.values.astype(np.float32).tobytes())
        chunk.drop(columns=['_capped'], inplace=True)

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
    if use_rm_mask:
        pct_rm = 100 * n_rm_downweighted / max(n_total, 1)
        if rm_mask_loader.repeat_weight == 0.0:
            print(f"[IO]   RM-excluded k-mers:     {n_rm_downweighted:,} ({pct_rm:.2f}%) "
                  f"(repeat positions → zeroed, excluded from winsorized mean)")
        else:
            print(f"[IO]   RM-downweighted k-mers: {n_rm_downweighted:,} ({pct_rm:.2f}%) "
                  f"(repeat positions, weight={rm_mask_loader.repeat_weight})")
    pct = 100 * n_bio_filtered_total / max(n_total, 1)
    print(f"[IO]   Bio-capped k-mers:    {n_bio_filtered_total:,} ({pct:.2f}%) "
          f"(soft-capped to threshold, not excluded)")
    total_bins = len(acc_values)
    print(f"[IO]   Total bins created:   {total_bins:,}")
    if use_rm_mask:
        n_unique_total = sum(len(b) for b in acc_unique_values.values())
        n_unique_bins  = sum(1 for b in acc_unique_values.values() if len(b) > 0)
        pct_u = 100 * n_unique_total / max(n_total, 1)
        print(f"[IO]   Pass-2 unique-only buffer: "
              f"{n_unique_total:,} k-mers across {n_unique_bins:,} bins "
              f"({pct_u:.2f}% of input)")
    return acc_values, acc_filtered, acc_unique_values




def compute_genome_median(acc, sex='XX', neutral_lo=0.4, neutral_hi=2.5, max_iter=5,
                          exclude_zeros=False):
    excluded_chroms = set()
    if sex == 'XY':
        excluded_chroms = CHRX_PATTERNS | CHRY_PATTERNS

    all_means = []
    seen_chroms = set()
    n_empty_bins = 0
    for (chrom, bin_start), buf in acc.items():
        if chrom in excluded_chroms:
            continue
        if len(buf) < 3:
            continue
        if exclude_zeros:
            val, n_used = _winsorized_mean(buf, exclude_zeros=True)
            if n_used == 0:
                n_empty_bins += 1
                continue
            all_means.append(val)
        else:
            all_means.append(_winsorized_mean(buf))
        seen_chroms.add(chrom)

    all_means = np.array(all_means, dtype=np.float64)
    n_total = len(all_means)

    if n_total == 0:
        print("[NORM] ERROR: No valid bins for genome median computation. "
              "Check input data and --sex setting.")
        sys.exit(1)

    genome_median = float(np.median(all_means))
    n_chroms = len(seen_chroms)
    if excluded_chroms:
        all_data_chroms = set(chrom for (chrom, _) in acc.keys())
        n_actually_excluded = len(all_data_chroms - seen_chroms)
        sex_note = f" (sex={sex}, {n_actually_excluded} chroms excluded)"
    else:
        sex_note = ""
    if exclude_zeros and n_empty_bins > 0:
        print(f"[NORM] RM exclusion: {n_empty_bins:,} bins with no unique k-mers "
              f"(skipped from median computation)")
    print(f"[NORM] Initial genome-wide depth median: {genome_median:.4f}  "
          f"({n_total:,} bins across {n_chroms} chroms{sex_note})")

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



def write_output(acc_values, acc_filtered, genome_median,
                 window_size, output_path, rm_mask_active=False,
                 rm_exclude=False, acc_unique_values=None):
    if acc_unique_values is None:
        acc_unique_values = {}

    rows = []
    n_empty_bins = 0
    n_pass2_bins = 0

    for (chrom, bin_start), buf in acc_values.items():
        n_win = len(buf)
        if n_win < 1:
            continue

        if rm_exclude:
            mean_cnt, n_unique = _winsorized_mean(buf, exclude_zeros=True)
            if n_unique == 0:
                mean_cnt = genome_median
                n_empty_bins += 1
        else:
            mean_cnt = _winsorized_mean(buf)

        cn = mean_cnt / genome_median if genome_median > 0 else 1.0
        cn = max(cn, EPSILON)

        log_ratio    = math.log2(cn)
        num_kmers    = n_win
        num_filtered = acc_filtered.get((chrom, bin_start), 0)

        unique_buf = acc_unique_values.get((chrom, bin_start))
        if unique_buf is not None and len(unique_buf) > 0:
            unique_mean = _winsorized_mean(unique_buf)
            n_unique_k  = len(unique_buf)
            n_pass2_bins += 1
            unique_mean_cell = round(float(unique_mean), 4)
        else:
            unique_mean_cell = float('nan')
            n_unique_k = 0

        rows.append((
            chrom,
            bin_start,
            bin_start + window_size,
            round(cn, 4),
            round(mean_cnt, 4),
            round(log_ratio, 6),
            num_kmers,
            num_filtered,
            unique_mean_cell,
            n_unique_k,
        ))

    rows.sort(key=lambda r: (natural_chrom_key(r[0]), r[1]))

    n_written = 0
    with open(output_path, 'w') as fh:
        fh.write(f"# CopySeg preprocessed windows\n")
        fh.write(f"# window_size={window_size} | normalization=genome_wide\n")
        fh.write(f"# genome_depth_median={genome_median:.4f}\n")
        if rm_exclude:
            fh.write(f"# rm_mask_mode=exclude | repeat k-mers excluded from CN computation\n")
        elif rm_mask_active:
            fh.write(f"# rm_mask_mode=downweight | repeat k-mers down-weighted in raw counts\n")
        if rm_mask_active:
            fh.write(f"# pass2_unique_buffer=on | columns 9-10 carry unique-only stats\n")
        else:
            fh.write(f"# pass2_unique_buffer=off | columns 9-10 are nan/0 (no rm_mask_dir)\n")
        fh.write("# chrom\tstart\tend\tcn\tmean_count\tlog_ratio\tnum_kmers\tnum_filtered\t"
                 "unique_mean_count\tn_unique_kmers\n")
        for row in rows:
            fh.write('\t'.join(map(str, row)) + '\n')
            n_written += 1

    if rm_exclude and n_empty_bins > 0:
        print(f"[IO] RM exclusion: {n_empty_bins:,} bins with no unique k-mers "
              f"→ set to CN=1.0 (conservative default)")
    if rm_mask_active:
        print(f"[IO] Pass-2 unique-only stats: {n_pass2_bins:,}/{n_written:,} bins "
              f"({100*n_pass2_bins/max(n_written,1):.1f}%) have ≥1 unique k-mer")
    print(f"[IO] Written {n_written:,} windows → {output_path}")
    return n_written



def main():
    parser = argparse.ArgumentParser(
        description="Preprocess k-mer BED → CopySeg window BED"
    )
    parser.add_argument('--input', '-i', required=True,
                        help='K-mer BED file (4-col: chrom, start, end, count)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output 10-col window BED (cols 9-10 carry Pass-2 '
                             'unique-only stats; NaN/0 when --rm-mask-dir absent).')
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
    parser.add_argument('--per-window-correct', action='store_true', default=False,
                        help='Per-window minimum multiplicity correction (Solution B). '
                             'Within each window, normalizes k-mer counts by a low '
                             'percentile — repeat element spikes are suppressed while '
                             'uniform segmental duplication signal is preserved. '
                             'No external files or precomputation needed. '
                             'Mutually exclusive with --weight-dir.')
    parser.add_argument('--pw-percentile', type=float, default=25.0,
                        help='Percentile used by --per-window-correct (default: 25). '
                             'Lower = more aggressive correction. p5 over-corrects '
                             '(696 Mb false HetDel); p25 is gentler.')
    parser.add_argument('--rm-mask-dir', default=None,
                        help='Per-chromosome RM binary mask directory (.npy files) '
                             'from compute_rm_mask.py. When provided, k-mers at '
                             'RepeatMasker-annotated positions are down-weighted '
                             '(default: ×0.01) so unique-sequence k-mers dominate '
                             'the window CN. Primary correction for k<=50 repeat '
                             'k-mer sharing inflation. No Jellyfish dependency.')
    parser.add_argument('--repeat-kmer-weight', type=float, default=REPEAT_KMER_WEIGHT,
                        help=f'Weight for k-mers at repeat positions when using '
                             f'--rm-mask-dir (default: {REPEAT_KMER_WEIGHT}). '
                             f'0.0 = full exclusion, 0.01 = near-exclusion, '
                             f'1.0 = no effect.')
    parser.add_argument('--rm-exclude', action='store_true', default=False,
                        help='RM exclusion mode: fully exclude repeat-position '
                             'k-mers from CN computation (only unique-sequence '
                             'k-mers contribute to window CN). Requires '
                             '--rm-mask-dir. Sets --repeat-kmer-weight to 0.0 '
                             'and excludes zeroed values from winsorized mean. '
                             'Bins with no unique k-mers default to CN=1.0.')
    parser.add_argument('--rm-class-weights', default=None,
                        help='Per-repeat-class weights for RM down-weighting. '
                             'Format: "Satellite=0.001,Simple_repeat=0.3,LINE=0.01,'
                             'SINE=0.001,LTR=0.05,DNA=0.1,Low_complexity=0.3,Other=0.05". '
                             'Requires --rm-mask-dir. Mutually exclusive with '
                             '--rm-exclude and --repeat-kmer-weight. '
                             'Different repeat classes have different k-mer '
                             'multiplicities at small k — this allows preserving '
                             'VNTR signal (Simple_repeat) while aggressively '
                             'suppressing Alu inflation (SINE).')
    parser.add_argument('--peak-method', choices=['local_maxima', 'gaussian_fit'],
                        default='local_maxima',
                        help='Single-copy depth peak detection method. '
                             '"local_maxima" (default): the existing K2 algorithm '
                             '— tallest local maximum after the noise tail. '
                             '"gaussian_fit": scipy.optimize curve-fit a Gaussian '
                             'centered on the local_maxima output, refine peak as '
                             'the fitted mu. The fitted sigma is logged (useful '
                             'for noise modeling). Falls back to local_maxima '
                             'when the fit fails (r²<0.70).')
    args = parser.parse_args()

    use_rm = args.rm_annotation_bed is not None

    if args.rm_exclude:
        if args.rm_mask_dir is None:
            print("ERROR: --rm-exclude requires --rm-mask-dir.",
                  file=sys.stderr)
            sys.exit(1)
        args.repeat_kmer_weight = 0.0

    parsed_class_weights = None
    if args.rm_class_weights is not None:
        if args.rm_mask_dir is None:
            print("ERROR: --rm-class-weights requires --rm-mask-dir.",
                  file=sys.stderr)
            sys.exit(1)
        if args.rm_exclude:
            print("ERROR: --rm-class-weights and --rm-exclude are mutually exclusive.",
                  file=sys.stderr)
            sys.exit(1)
        parsed_class_weights = {}
        for pair in args.rm_class_weights.split(','):
            pair = pair.strip()
            if '=' not in pair:
                print(f"ERROR: Invalid --rm-class-weights format: '{pair}'. "
                      f"Expected 'ClassName=weight'.", file=sys.stderr)
                sys.exit(1)
            name, val = pair.split('=', 1)
            name = name.strip()
            if name not in RM_CLASS_NAME_TO_ID:
                print(f"ERROR: Unknown repeat class '{name}'. "
                      f"Valid: {list(RM_CLASS_NAME_TO_ID.keys())}", file=sys.stderr)
                sys.exit(1)
            cid = RM_CLASS_NAME_TO_ID[name]
            parsed_class_weights[cid] = float(val)

    if args.weight_dir is not None and args.per_window_correct:
        print("ERROR: --weight-dir and --per-window-correct are mutually exclusive.",
              file=sys.stderr)
        sys.exit(1)
    if args.weight_dir is not None and args.rm_mask_dir is not None:
        print("ERROR: --weight-dir and --rm-mask-dir are mutually exclusive. "
              "Both correct k-mer multiplicity — using both would double-correct.",
              file=sys.stderr)
        sys.exit(1)
    if args.rm_mask_dir is not None and args.per_window_correct:
        print("ERROR: --rm-mask-dir and --per-window-correct are mutually exclusive. "
              "Both address repeat k-mer inflation — using both would over-correct.",
              file=sys.stderr)
        sys.exit(1)

    weight_loader = None
    if args.weight_dir is not None:
        if not os.path.isdir(args.weight_dir):
            print(f"ERROR: --weight-dir not found: {args.weight_dir}", file=sys.stderr)
            sys.exit(1)
        weight_loader = WeightLoader(args.weight_dir)

    rm_mask_loader = None
    if args.rm_mask_dir is not None:
        if not os.path.isdir(args.rm_mask_dir):
            print(f"ERROR: --rm-mask-dir not found: {args.rm_mask_dir}", file=sys.stderr)
            sys.exit(1)
        rm_mask_loader = RMMaskLoader(args.rm_mask_dir,
                                      repeat_weight=args.repeat_kmer_weight,
                                      class_weights=parsed_class_weights)

    sex_desc = "all chroms diploid" if args.sex == 'XX' else "chrX/chrY excluded from median"
    print("=" * 60)
    print("CopySeg — K-mer BED Preprocessing")
    print(f"  Sliding windows → {args.window_size}bp bins")
    print(f"  Sex: {args.sex} ({sex_desc})")
    print("  Normalization: genome-wide neutral-band median")
    if weight_loader is not None:
        print(f"  Multiplicity correction: {args.weight_dir}")
        print(f"    k_target={weight_loader.k_target}, k_ref={weight_loader.k_ref}")
    if rm_mask_loader is not None:
        if args.rm_exclude:
            print(f"  RM k-mer EXCLUSION: {args.rm_mask_dir}")
            print(f"    repeat k-mers fully excluded — CN from unique k-mers only")
        elif parsed_class_weights is not None:
            id_to_name = {v: k for k, v in RM_CLASS_NAME_TO_ID.items()}
            print(f"  RM per-class weights: {args.rm_mask_dir}")
            for cid in sorted(parsed_class_weights.keys()):
                cname = id_to_name.get(cid, f"class_{cid}")
                print(f"    {cname}={parsed_class_weights[cid]}")
        else:
            print(f"  RM k-mer down-weighting: {args.rm_mask_dir}")
            print(f"    repeat_weight={args.repeat_kmer_weight}")
    if args.per_window_correct:
        print("  Per-window multiplicity correction: ENABLED")
        print(f"    Method: within-window p{args.pw_percentile:g} normalization (no external DB)")
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
    if rm_mask_loader is not None:
        print(f"RM mask dir: {args.rm_mask_dir}")
    print()

    bio_threshold = 0.0
    gaussian_peak = 0
    if args.bio_threshold_factor > 0:
        print("[HIST] Sampling k-mer count histogram for Gaussian peak detection...")
        histogram = build_histogram_sample(args.input, args.hist_sample_chunks,
                                              weight_loader=weight_loader,
                                              per_window_correct=args.per_window_correct,
                                              window_size=args.window_size,
                                              pw_percentile=args.pw_percentile,
                                              rm_mask_loader=rm_mask_loader)
        gaussian_peak = find_gaussian_peak(histogram)

        if args.peak_method == 'gaussian_fit':
            mu_fit, sigma_fit, r_sq, ok = fit_gaussian_to_histogram(
                histogram, p0_peak=gaussian_peak)
            print(f"[HIST] Gaussian curve fit: mu={mu_fit:.3f}, "
                  f"sigma={sigma_fit:.3f}, r²={r_sq:.4f} "
                  f"({'accepted' if ok else 'rejected — fallback to local_maxima'})")
            if ok:
                refined_peak = max(2, int(round(mu_fit)))
                if refined_peak != gaussian_peak:
                    print(f"[HIST] Peak refined: {gaussian_peak} → {refined_peak} "
                          f"(via Gaussian fit; sigma={sigma_fit:.2f} ≈ "
                          f"{100*sigma_fit/refined_peak:.1f}% CV)")
                gaussian_peak = refined_peak

        bio_threshold = gaussian_peak * args.bio_threshold_factor
        print(f"[HIST] Biological threshold: {gaussian_peak} × "
              f"{args.bio_threshold_factor} = {bio_threshold:.1f}")
        print()

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

    acc_values, acc_filtered, acc_unique_values = aggregate_windows(
        args.input, args.window_size,
        bio_threshold=bio_threshold,
        bin_threshold_map=bin_threshold_map,
        weight_loader=weight_loader,
        per_window_correct=args.per_window_correct,
        pw_percentile=args.pw_percentile,
        rm_mask_loader=rm_mask_loader)

    rm_exclude = (rm_mask_loader is not None and args.repeat_kmer_weight == 0.0)

    print()
    genome_median = compute_genome_median(acc_values, sex=args.sex,
                                          exclude_zeros=rm_exclude)

    print()
    n = write_output(acc_values, acc_filtered, genome_median,
                     args.window_size, args.output,
                     rm_mask_active=(rm_mask_loader is not None),
                     rm_exclude=rm_exclude,
                     acc_unique_values=acc_unique_values)

    print()
    print(f"Done. {n:,} windows written.")
    return 0


if __name__ == '__main__':
    sys.exit(main())

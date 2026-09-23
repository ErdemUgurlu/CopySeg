#!/usr/bin/env python3

import argparse
import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.abspath(os.path.dirname(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from segment_cnv_fused_lasso import assign_state_from_cn
from chrom_utils import is_chrx, is_chry


TANDEM_CLASSES = frozenset({'Simple_repeat', 'Satellite', 'Low_complexity', 'rRNA'})

DEFAULT_MIN_UNIQUE_KMERS_PER_BIN = 20
DEFAULT_MIN_TRUST_BINS           = 10
DEFAULT_MIN_TRUST_FRACTION       = 0.05

DEFAULT_MIN_BINS_FOR_GENOME_MEDIAN = 100



def compute_unique_genome_median(windows_df, sex='XX',
                                 neutral_lo=0.4, neutral_hi=2.5,
                                 max_iter=5,
                                 min_bins=DEFAULT_MIN_BINS_FOR_GENOME_MEDIAN):
    vals = windows_df['unique_mean_count'].values.astype(np.float64)
    n_unique = windows_df['n_unique_kmers'].values
    chroms = windows_df['chrom'].values

    valid = np.isfinite(vals) & (n_unique > 0) & (vals > 0)

    if sex == 'XY':
        for cname in np.unique(chroms):
            if is_chrx(str(cname)) or is_chry(str(cname)):
                valid &= (chroms != cname)

    arr = vals[valid]
    if arr.size == 0:
        return None
    if arr.size < min_bins:
        print(f"[PASS2-NORM] WARNING: only {arr.size} valid unique bins "
              f"(< {min_bins} required). Refusing to compute a normalizer "
              f"from this little data — all segments will route to fallback. "
              f"Check that --rm-mask-dir was set during Pass-1.")
        return None

    med = float(np.median(arr))
    print(f"[PASS2-NORM] Initial unique-only median: {med:.4f} "
          f"({arr.size:,} valid bins)")

    for i in range(max_iter):
        lo = neutral_lo * med
        hi = neutral_hi * med
        in_band = (arr >= lo) & (arr <= hi)
        n_kept = int(in_band.sum())
        if n_kept == 0:
            break
        new_med = float(np.median(arr[in_band]))
        change = abs(new_med - med) / med if med > 0 else 0.0
        print(f"[PASS2-NORM] Iter {i+1}: band [{lo:.2f}, {hi:.2f}] "
              f"→ {n_kept:,} bins kept, median {med:.4f} → {new_med:.4f} "
              f"(Δ={100*change:.3f}%)")
        med = new_med
        if change < 0.001:
            print(f"[PASS2-NORM] Converged after {i+1} iterations")
            break

    print(f"[PASS2-NORM] Final unique-only median: {med:.4f}")
    return med



def load_segments(seg_path):
    print(f"[IO] Loading segments: {seg_path}")
    header = None
    with open(seg_path) as fh:
        for line in fh:
            if line.startswith('#'):
                header = line.lstrip('#').rstrip('\n').split('\t')
                break
            else:
                break
    if header is None:
        raise ValueError(f"Segments BED missing #header line: {seg_path}")

    df = pd.read_csv(seg_path, sep='\t', comment='#', header=None,
                     names=header)
    print(f"[IO] Loaded {len(df):,} segments ({len(header)} columns: "
          f"{', '.join(header)})")
    required = {'chrom', 'start', 'end', 'cn_median'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Segments BED missing required columns: {missing}")
    return df, header


def load_windows(win_path):
    print(f"[IO] Loading windows: {win_path}")
    df = pd.read_csv(win_path, sep='\t', comment='#', header=None,
                     dtype={0: str})
    n_cols = len(df.columns)
    if n_cols == 8:
        df.columns = ['chrom', 'start', 'end', 'cn', 'mean_count',
                      'log_ratio', 'num_kmers', 'num_filtered']
        df['unique_mean_count'] = np.nan
        df['n_unique_kmers']    = 0
        print(f"[IO] Legacy 8-col file detected — Pass-2 columns synthesized "
              f"(unique_mean_count=NaN, n_unique_kmers=0). All segments will "
              f"route to fallback_too_few_unique.")
    elif n_cols >= 10:
        df.columns = (['chrom', 'start', 'end', 'cn', 'mean_count',
                       'log_ratio', 'num_kmers', 'num_filtered',
                       'unique_mean_count', 'n_unique_kmers']
                      + [f'col{i}' for i in range(10, n_cols)])
        df['unique_mean_count'] = pd.to_numeric(df['unique_mean_count'],
                                                errors='coerce')
        df['n_unique_kmers']    = pd.to_numeric(df['n_unique_kmers'],
                                                errors='coerce').fillna(0).astype(int)
    else:
        raise ValueError(f"Windows BED has unexpected column count: "
                         f"{n_cols} (need 8 or 10)")
    print(f"[IO] Loaded {len(df):,} windows  "
          f"(unique_mean_count: {df['unique_mean_count'].notna().sum():,} non-NaN; "
          f"n_unique_kmers>0: {(df['n_unique_kmers'] > 0).sum():,})")
    return df



def build_windows_index(windows_df):
    by_chrom = {}
    for chrom, grp in windows_df.groupby('chrom', sort=False):
        grp = grp.sort_values('start').reset_index(drop=True)
        by_chrom[chrom] = {
            'starts':            grp['start'].values,
            'unique_mean_count': grp['unique_mean_count'].values,
            'n_unique_kmers':    grp['n_unique_kmers'].values.astype(np.int64),
        }
    return by_chrom


def refine_segment(seg_row, windows_idx, unique_genome_median,
                   min_unique_per_bin=DEFAULT_MIN_UNIQUE_KMERS_PER_BIN,
                   min_trust_bins=DEFAULT_MIN_TRUST_BINS,
                   min_trust_fraction=DEFAULT_MIN_TRUST_FRACTION,
                   estimator_percentile=50.0):
    cn_pass1 = float(seg_row['cn_median'])
    repeat_class = str(seg_row.get('repeat_class', '') or '')

    chrom_data = windows_idx.get(seg_row['chrom'])
    if chrom_data is None:
        return cn_pass1, 'fallback_too_few_unique', 0, 0

    starts = chrom_data['starts']
    lo = int(np.searchsorted(starts, int(seg_row['start']), side='left'))
    hi = int(np.searchsorted(starts, int(seg_row['end']),   side='left'))
    n_total = hi - lo
    if n_total <= 0:
        return cn_pass1, 'fallback_too_few_unique', 0, n_total

    bin_unique     = chrom_data['unique_mean_count'][lo:hi]
    bin_n_unique   = chrom_data['n_unique_kmers'][lo:hi]
    trust_mask = (bin_n_unique >= min_unique_per_bin) & np.isfinite(bin_unique)
    n_trust = int(trust_mask.sum())

    if repeat_class in TANDEM_CLASSES:
        return cn_pass1, 'fallback_repeat_class', n_trust, n_total

    if (n_trust < min_trust_bins
            or (n_trust / max(n_total, 1)) < min_trust_fraction):
        return cn_pass1, 'fallback_too_few_unique', n_trust, n_total

    trust_vals = bin_unique[trust_mask]
    med = float(np.percentile(trust_vals, estimator_percentile))
    if not np.isfinite(med) or unique_genome_median is None or unique_genome_median <= 0:
        return cn_pass1, 'fallback_nan', n_trust, n_total

    cn_refined = float(med / unique_genome_median)
    cn_refined = max(cn_refined, 1e-3)
    return cn_refined, 'unique_median', n_trust, n_total


def refine_all(segments_df, windows_df, sex='XX',
               min_unique_per_bin=DEFAULT_MIN_UNIQUE_KMERS_PER_BIN,
               min_trust_bins=DEFAULT_MIN_TRUST_BINS,
               min_trust_fraction=DEFAULT_MIN_TRUST_FRACTION,
               lowdup_threshold=1.25,
               min_bins_for_genome_median=DEFAULT_MIN_BINS_FOR_GENOME_MEDIAN,
               estimator_percentile=50.0):
    print(f"\n[PASS2] Refining {len(segments_df):,} segments "
          f"(unique-only estimator = p{estimator_percentile:g})...")

    unique_genome_median = compute_unique_genome_median(
        windows_df, sex=sex, min_bins=min_bins_for_genome_median)
    if unique_genome_median is None:
        print("[PASS2-NORM] WARNING: No valid bins for unique-only normalizer. "
              "All segments will route to fallback.")

    windows_idx = build_windows_index(windows_df)

    cn_refined_arr   = np.empty(len(segments_df), dtype=np.float64)
    method_arr       = np.empty(len(segments_df), dtype=object)
    n_trust_arr      = np.empty(len(segments_df), dtype=np.int64)
    n_total_arr      = np.empty(len(segments_df), dtype=np.int64)

    for i, (_, seg) in enumerate(segments_df.iterrows()):
        cn_ref, method, n_trust, n_total = refine_segment(
            seg, windows_idx, unique_genome_median,
            min_unique_per_bin=min_unique_per_bin,
            min_trust_bins=min_trust_bins,
            min_trust_fraction=min_trust_fraction,
            estimator_percentile=estimator_percentile,
        )
        cn_refined_arr[i] = cn_ref
        method_arr[i]     = method
        n_trust_arr[i]    = n_trust
        n_total_arr[i]    = n_total

    out = segments_df.copy()
    out['cn_refined']    = cn_refined_arr
    out['refine_method'] = method_arr
    out['n_trust_bins']  = n_trust_arr
    out['n_total_bins']  = n_total_arr
    out['state_refined'] = [
        assign_state_from_cn(float(cn), lowdup_threshold=lowdup_threshold)
        for cn in cn_refined_arr
    ]

    method_counts = pd.Series(method_arr).value_counts()
    print(f"\n[PASS2] Method breakdown:")
    for m, c in method_counts.items():
        pct = 100 * c / len(segments_df)
        print(f"[PASS2]   {m:<28s} {c:>7,d}  ({pct:5.1f}%)")

    n_state_changed = int(
        (out['state_refined'] != out.get('state', out['state_refined'])).sum()
    ) if 'state' in out.columns else 0
    if 'state' in out.columns:
        print(f"[PASS2] State relabeled: {n_state_changed:,} segments "
              f"({100*n_state_changed/len(out):.1f}%)")
    return out



def write_refined(refined_df, output_path, original_columns):
    new_cols = ['cn_refined', 'state_refined', 'refine_method',
                'n_trust_bins', 'n_total_bins']
    col_order = original_columns + [c for c in new_cols
                                    if c not in original_columns]
    refined_df = refined_df[col_order]

    print(f"[IO] Writing {len(refined_df):,} refined segments → {output_path}")
    with open(output_path, 'w') as fh:
        fh.write('#' + '\t'.join(col_order) + '\n')
        for row in refined_df.itertuples(index=False, name=None):
            fh.write('\t'.join(_fmt(v) for v in row) + '\n')


def _fmt(v):
    if isinstance(v, float):
        if np.isnan(v):
            return 'nan'
        return f"{v:.4f}"
    return str(v)



def main():
    parser = argparse.ArgumentParser(
        description="Pass-2 segment CN refinement using unique-position k-mers")
    parser.add_argument('--segments', required=True,
                        help='Pass-1 segments BED (segs_cnacc_w500.bed)')
    parser.add_argument('--windows', required=True,
                        help='Pass-1 windows BED with unique columns (cn_w500.bed)')
    parser.add_argument('--output', required=True,
                        help='Output refined segments BED')
    parser.add_argument('--sex', choices=['XX', 'XY'], default='XX',
                        help='Sample sex karyotype (default XX). For XY, '
                             'chrX/chrY bins are excluded from unique_genome_median.')
    parser.add_argument('--min-unique-kmers-per-bin', type=int,
                        default=DEFAULT_MIN_UNIQUE_KMERS_PER_BIN,
                        help='Minimum n_unique_kmers per bin to count as "trust" '
                             f'(default {DEFAULT_MIN_UNIQUE_KMERS_PER_BIN})')
    parser.add_argument('--min-trust-bins', type=int,
                        default=DEFAULT_MIN_TRUST_BINS,
                        help='Absolute floor on number of trust bins per segment '
                             f'(default {DEFAULT_MIN_TRUST_BINS})')
    parser.add_argument('--min-trust-fraction', type=float,
                        default=DEFAULT_MIN_TRUST_FRACTION,
                        help='Fraction of segment bins that must be trust bins '
                             f'(default {DEFAULT_MIN_TRUST_FRACTION})')
    parser.add_argument('--lowdup-threshold', type=float, default=1.25,
                        help='Neutral/LowDup CN boundary for state assignment '
                             '(default 1.25, matches segmenter)')
    parser.add_argument('--min-bins-for-genome-median', type=int,
                        default=DEFAULT_MIN_BINS_FOR_GENOME_MEDIAN,
                        help='Minimum valid bin count for unique_genome_median. '
                             'Below this, the normalizer is too noisy and Pass-2 '
                             'refuses to compute it; all segments fall back to '
                             f'Pass-1 cn_median (default {DEFAULT_MIN_BINS_FOR_GENOME_MEDIAN}).')
    parser.add_argument('--percentile-estimator', type=float, default=50.0,
                        help='Percentile of unique-position counts used as the per-segment '
                             'central estimate (default 50 = median). Higher values (e.g. 67, 75) '
                             'shift toward high-count k-mers to recover high-CN/paralogous loci '
                             '(under-called by the median). Range 50-100.')
    args = parser.parse_args()

    segments_df, header = load_segments(args.segments)
    windows_df = load_windows(args.windows)

    refined = refine_all(
        segments_df, windows_df,
        sex=args.sex,
        min_unique_per_bin=args.min_unique_kmers_per_bin,
        min_trust_bins=args.min_trust_bins,
        min_trust_fraction=args.min_trust_fraction,
        lowdup_threshold=args.lowdup_threshold,
        min_bins_for_genome_median=args.min_bins_for_genome_median,
        estimator_percentile=args.percentile_estimator,
    )

    write_refined(refined, args.output, original_columns=header)
    print(f"\n[DONE] Wrote {len(refined):,} segments → {args.output}")
    return 0


if __name__ == '__main__':
    sys.exit(main())

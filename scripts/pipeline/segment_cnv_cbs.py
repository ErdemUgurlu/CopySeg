#!/usr/bin/env python3

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


EPSILON = 1e-3
WINDOW_BP = 500

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



def assign_state_from_cn(cn):
    for upper, state in CN_STATE_BOUNDARIES:
        if cn < upper:
            return state
    return 'HighAmp'



def load_data(bed_file, min_kmers=0):
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

        if min_kmers > 0 and 'num_kmers' in df.columns:
            low_cov = df['num_kmers'] < min_kmers
            n_low = int(low_cov.sum())
            if n_low > 0:
                df.loc[low_cov, 'cn'] = 1.0
                print(f"[IO] Low-coverage filter (num_kmers < {min_kmers}): "
                      f"{n_low:,} windows ({100 * n_low / len(df):.1f}%) set to CN=1.0")

        if 'num_kmers' in df.columns and 'num_filtered' in df.columns:
            df['quality'] = (df['num_kmers'] /
                             (df['num_kmers'] + df['num_filtered'] + 1))
        else:
            df['quality'] = 1.0

        df['log2_cn'] = np.log2(df['cn'].values + EPSILON)

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
    print(f"[GC] Loading GC content from {gc_bed_file}...")
    gc_df = pd.read_csv(gc_bed_file, sep='\t', comment='#', header=None)
    if len(gc_df.columns) < 4:
        raise ValueError(f"GC BED needs >= 4 columns, got {len(gc_df.columns)}")
    gc_df.columns = ['chrom', 'start', 'end', 'gc_content'][:len(gc_df.columns)]
    gc_df = gc_df[['chrom', 'start', 'gc_content']]
    print(f"[GC] Loaded {len(gc_df):,} windows (gc mean={gc_df['gc_content'].mean():.4f})")
    return gc_df


def load_repeat_annotation(repeat_bed_file):
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



def make_segment_from_windows(windows_df):
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
        'avg_entropy': 0.0,
        'max_entropy': 0.0,
    }


def run_pelt_chrom(df_chrom, penalty, min_size_windows):
    signal = df_chrom['log2_cn'].values
    n = len(signal)

    if n == 0:
        return []

    if n < max(min_size_windows * 2, 4):
        return [make_segment_from_windows(df_chrom)]

    algo = rpt.Pelt(model='l2', min_size=min_size_windows, jump=1)
    breakpoints = algo.fit_predict(signal, pen=penalty)

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
    import time as _t

    window_sizes = (df['end'] - df['start']).values
    window_bp = int(np.median(window_sizes)) if len(window_sizes) > 0 else WINDOW_BP
    if window_bp != WINDOW_BP:
        print(f"[CBS] Non-standard window size detected: {window_bp} bp (expected {WINDOW_BP})")

    min_size_windows = max(1, min_size // window_bp)
    print(f"[CBS] min_segment_length={min_size} bp → min_size_windows={min_size_windows}")
    print(f"[CBS] PELT penalty: {'BIC (auto)' if penalty is None else penalty}")

    chroms = list(df['chrom'].unique())
    all_segments = []
    total_t0 = _t.time()

    for chrom in chroms:
        t0 = _t.time()
        df_chrom = df[df['chrom'] == chrom].sort_values('start').reset_index(drop=True)
        n = len(df_chrom)

        pen = penalty if penalty is not None else (2.0 * np.log(max(n, 2)))

        segs = run_pelt_chrom(df_chrom, pen, min_size_windows)
        all_segments.extend(segs)

        elapsed = _t.time() - t0
        print(f"[CBS] {chrom}: {n:,} windows → {len(segs):,} segments ({elapsed:.1f}s)")

    total_elapsed = _t.time() - total_t0
    print(f"[CBS] Total: {len(df):,} windows → {len(all_segments):,} segments "
          f"({total_elapsed:.1f}s)")
    return all_segments



def compute_boundary_confidence(segments, df):
    if not segments:
        return segments

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

        lo_l = int(np.searchsorted(starts_arr, left['start'], side='left'))
        hi_l = int(np.searchsorted(starts_arr, left['end'],   side='left'))
        cn_left = cn_arr[lo_l:hi_l]

        lo_r = int(np.searchsorted(starts_arr, right['start'], side='left'))
        hi_r = int(np.searchsorted(starts_arr, right['end'],   side='left'))
        cn_right = cn_arr[lo_r:hi_r]

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



def _pooled_std(s1, n1, s2, n2, m1, m2):
    n1, n2 = max(int(n1), 1), max(int(n2), 1)
    within  = max(n1 - 1, 0) * s1 ** 2 + max(n2 - 1, 0) * s2 ** 2
    between = (n1 * n2) / (n1 + n2) * (m1 - m2) ** 2
    return float(np.sqrt(max(0.0, (within + between) / (n1 + n2 - 1))))



def filter_small_segments(segments, min_lengths):
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



def write_output(segments, output_file, extended=False):
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



def main():
    parser = argparse.ArgumentParser(
        description="CopySeg CBS/PELT CN Segmentation — drop-in alternative to HMM"
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

    if args.mode is not None:
        print(f"[INFO] --mode '{args.mode}' is HMM-only; CBS always runs in cn-accuracy equivalent mode")
    if args.no_pomegranate:
        print(f"[INFO] --no-pomegranate is a no-op for CBS (pomegranate not used)")
    if args.coarse_output:
        print(f"[INFO] --coarse-output is not implemented in CBS segmenter")

    print("=" * 60)
    print("CopySeg CBS/PELT Segmentation")
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

    df = load_data(args.input, min_kmers=args.min_kmers)

    gc_df = None
    if args.gc_content_bed:
        gc_df = load_gc_content(args.gc_content_bed)
        print()

    print()
    segments = run_segmentation(df, penalty=args.penalty,
                                min_size=args.min_segment_length)

    print()

    min_lengths = {
        'LowDup':  args.min_segment_length,
        'HighDup': args.min_segment_length,
        'Amp':     args.min_segment_length,
    }
    segments = filter_small_segments(segments, min_lengths)

    segments = filter_low_quality_segments(
        segments,
        threshold=args.quality_threshold,
        cv_filter_threshold=args.cv_filter_threshold,
    )

    segments = merge_nearby_dup_segments(segments, max_gap=10000)

    if gc_df is not None:
        print()
        segments = compute_segment_mean_gc(segments, df, gc_df)
        segments = apply_gc_cn_calibration(segments)

    if args.cn_reclassify_threshold > 0:
        print()
        segments = reclassify_by_cn_threshold(
            segments, lowdup_threshold=args.cn_reclassify_threshold)

    if args.cv_split_threshold > 0:
        print()
        segments = split_high_cv_segments(
            segments, df, cv_threshold=args.cv_split_threshold,
            min_length=args.min_segment_length)

    if args.repeat_bed:
        print()
        rm_df = load_repeat_annotation(args.repeat_bed)
        segments = compute_segment_repeat_annotation(segments, df, rm_df)

    print()
    segments = compute_boundary_confidence(segments, df)

    print_statistics(segments)
    do_extended = True
    write_output(segments, args.output, extended=do_extended)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
CN Accuracy Validation for KonuSeg HMM Output.
===============================================

Validates segment-level copy number accuracy using four metrics:

  A. State-CN concordance
       Fraction of dup segments where |cn_median - expected_cn| / expected_cn < 0.3.
       Expected CN: LowDup=2, HighDup=4, Amp=8, HighAmp=16.
       NOTE: For HighAmp (CN=16-150+), the median raw value IS the ground truth —
       the expected_cn=16 is just the lower bound.  A separate high-CN analysis
       shows the actual distribution.

  B. CN consistency (within-segment CV)
       cv = cn_std / cn_median per segment.
       Good: cv < 0.3.  Noisy / fragmented: cv > 0.5.
       Requires extended HMM output (--extended or cn-accuracy mode).

  C. AMY1 locus
       Region: CM039011.1:103,500,000-103,750,000 (T2T-CHM13v2.0).
       CHM13 carries 7 AMY1 copies. Assembly-validated (Jellyfish k=72): CN=7.0.
       Expected observed CN range: 5.0-9.0 (accounting for segmentation variation).

  D. ONT depth correlation (optional)
       Pearson r between per-window ONT depth and the CN assigned by HMM.
       Expected r > 0.85 for well-calibrated output.
       Requires --bam-features from extract_bam_features.py.

Usage:
  python3 scripts/validate_cn_accuracy.py \\
    --segments output/c5_my_run/segs_cnacc.bed \\
    [--bam-features output/c5_my_run/ont_bam_features_w500.tsv]

Author: KonuSeg Team
Date: March 2026
"""

import argparse
import sys
import numpy as np
import pandas as pd

# Expected CN for HMM dup states (haploid reference, CN=1 is normal).
# cn-accuracy mode uses MedAmp/ExtremeAmp instead of HomDel/HetDel.
EXPECTED_CN = {
    'LowDup':     2.0,
    'HighDup':    4.0,
    'Amp':        8.0,
    'MedAmp':    12.0,   # cn-accuracy: CN≈12 (log2=3.58); fills CN=9-22 gap
    'HighAmp':   22.0,   # cn-accuracy: CN≈22 (log2=4.5); was 16 floor
    'ExtremeAmp': 64.0,  # cn-accuracy: CN≈64 (log2=6.0); ultra-high amplification
}

# Concordance computed only on states where expected CN is well-defined.
# High-CN states (MedAmp/HighAmp/ExtremeAmp) have open-ended CN ranges
# and previously inflated concordance via auto-PASS (cn >= expected → error=0).
CONCORDANCE_STATES = {'LowDup': 2.0, 'HighDup': 4.0, 'Amp': 8.0}
HIGH_CN_STATES = {'MedAmp': 12.0, 'HighAmp': 22.0, 'ExtremeAmp': 64.0}

DUP_STATES = set(EXPECTED_CN.keys())

# AMY1A locus on chromosome 1 (same coordinates on CHM13v2.0 reference).
# Expected CN and range differ per sample.
_AMY1A_NAMES = {'chr1', 'CM039011.1'}
AMY1A_START = 103_500_000
AMY1A_END   = 103_750_000

# Per-sample AMY1 expected CN and acceptable range
_AMY1_SAMPLE_CONFIG = {
    'chm13': {'expected': 7.0, 'range': (5.0, 9.0), 'copies': 7,
              'label': 'CHM13, 7 copies, Jellyfish k=72 validated'},
    'hg002': {'expected': 3.0, 'range': (2.0, 4.5), 'copies': 3,
              'label': 'HG002, 3 copies, assembly-validated'},
}


# =============================================================================
# I/O
# =============================================================================

def load_segments(seg_file):
    """Load HMM segments BED file. Handles both minimal and extended formats."""
    print(f"[IO] Loading segments from {seg_file}...")
    df = pd.read_csv(seg_file, sep='\t', comment='#', header=None)

    base_cols  = ['chrom', 'start', 'end', 'state', 'cn_median', 'cn_mean', 'n_windows']
    ext_cols   = ['avg_quality', 'min_quality', 'cn_std', 'avg_repeats',
                  'avg_entropy', 'max_entropy', 'masked_fraction', 'repeat_class']
    # Extra columns from segment_cnv_fused_lasso.py extended output (18-col)
    extra_cols = ['gc_bias_factor', 'segment_iqr', 'boundary_conf']
    all_cols = base_cols + ext_cols + extra_cols  # 18 total
    ncols = len(df.columns)
    if ncols > len(all_cols):
        # Unknown extra columns: name them col_N to prevent crash
        all_cols = all_cols + [f'col_{i}' for i in range(len(all_cols), ncols)]
    df.columns = all_cols[:ncols]

    df['start']     = df['start'].astype(int)
    df['end']       = df['end'].astype(int)
    df['cn_median'] = df['cn_median'].astype(float)
    df['length']    = df['end'] - df['start']

    print(f"[IO] Loaded {len(df):,} segments")
    return df


def load_bam_features(bam_file):
    """Load ONT BAM feature TSV (chrom, start, end, ont_depth, ...)."""
    print(f"[IO] Loading BAM features from {bam_file}...")
    cols = ['chrom', 'start', 'end', 'ont_depth', 'ont_mapq_mean',
            'ont_unique_frac', 'ont_secondary_frac']
    df = pd.read_csv(bam_file, sep='\t', comment='#', header=None)
    df.columns = cols[:len(df.columns)]
    df['start']     = df['start'].astype(int)
    if 'end' in df.columns:
        df['end'] = df['end'].astype(int)
    df['ont_depth'] = df['ont_depth'].astype(float)
    print(f"[IO] Loaded {len(df):,} windows")
    return df


# =============================================================================
# METRIC A — State-CN concordance
# =============================================================================

def metric_a_concordance(df):
    """Report per-state CN distribution and overall concordance with expected CN.

    Concordance is computed only on CONCORDANCE_STATES (LowDup, HighDup, Amp)
    where expected CN is well-defined. High-CN states (MedAmp/HighAmp/ExtremeAmp)
    are reported separately without the cn >= expected → error=0 override that
    previously inflated overall concordance.
    """
    print("=" * 60)
    print("A. State-CN Concordance")
    print("=" * 60)
    print(f"  Concordance states: {CONCORDANCE_STATES}")
    print(f"  High-CN states (reported separately): {HIGH_CN_STATES}")
    print()

    # --- Concordance states (LowDup, HighDup, Amp) ---
    conc_errors = []

    for state, expected in CONCORDANCE_STATES.items():
        sdf = df[df['state'] == state]
        if len(sdf) == 0:
            print(f"  {state:10}: no segments")
            continue

        cn = sdf['cn_median'].values
        errors = np.abs(cn - expected) / expected
        frac_30 = (errors < 0.3).mean()
        frac_50 = (errors < 0.5).mean()
        conc_errors.extend(errors.tolist())

        p5, p50, p95 = np.percentile(cn, [5, 50, 95])
        print(f"  {state:10}: n={len(sdf):>6,}  "
              f"CN [p5={p5:.1f}, med={p50:.1f}, p95={p95:.1f}]  "
              f"within_30%={frac_30:.1%}  within_50%={frac_50:.1%}")

    if conc_errors:
        n_conc = len(df[df['state'].isin(CONCORDANCE_STATES)])
        frac_overall = np.mean(np.array(conc_errors) < 0.3)
        print(f"\n  CONCORDANCE: {frac_overall:.1%} of {n_conc:,} dup segments within 30% of expected CN")
        if frac_overall >= 0.60:
            print("  Status: GOOD")
        elif frac_overall >= 0.40:
            print("  Status: MODERATE")
        else:
            print("  Status: POOR — QW underestimation suspected (use --mode cn-accuracy)")

    # --- High-CN states (separate report, no auto-PASS) ---
    print()
    print("  High-CN Distribution (informational — not included in concordance):")
    for state, expected in HIGH_CN_STATES.items():
        sdf = df[df['state'] == state]
        if len(sdf) == 0:
            print(f"  {state:10}: no segments")
            continue
        cn = sdf['cn_median'].values
        p5, p50, p95 = np.percentile(cn, [5, 50, 95])
        errors = np.abs(cn - expected) / expected
        frac_30 = (errors < 0.3).mean()
        print(f"  {state:10}: n={len(sdf):>6,}  "
              f"CN [p5={p5:.1f}, med={p50:.1f}, p95={p95:.1f}]  "
              f"within_30%={frac_30:.1%}")
    print()


# =============================================================================
# METRIC B — CN consistency
# =============================================================================

def metric_b_consistency(df):
    """Report within-segment CN consistency (cv = cn_std / cn_median)."""
    print("=" * 60)
    print("B. CN Consistency (within-segment CV = cn_std / cn_median)")
    print("=" * 60)

    if 'cn_std' not in df.columns:
        print("  cn_std column not found — run HMM with --mode cn-accuracy or --extended")
        print()
        return

    dup_df = df[df['state'].isin(DUP_STATES) & (df['cn_median'] > 0.1)].copy()
    if len(dup_df) == 0:
        print("  No dup segments found.")
        print()
        return

    dup_df['cv'] = dup_df['cn_std'] / (dup_df['cn_median'] + 0.001)

    cv_med      = dup_df['cv'].median()
    frac_good   = (dup_df['cv'] < 0.3).mean()
    frac_noisy  = (dup_df['cv'] > 0.5).mean()

    print(f"  n dup segments:       {len(dup_df):,}")
    print(f"  CV median:            {cv_med:.3f}")
    print(f"  Fraction CV < 0.3 (good):   {frac_good:.1%}")
    print(f"  Fraction CV > 0.5 (noisy):  {frac_noisy:.1%}")

    # Per-state breakdown
    print()
    print("  Per-state breakdown:")
    for state in DUP_STATES:
        sdf = dup_df[dup_df['state'] == state]
        if len(sdf) == 0:
            continue
        print(f"    {state:10}: n={len(sdf):>5,}  cv_med={sdf['cv'].median():.3f}  "
              f"good={( sdf['cv'] < 0.3).mean():.1%}")

    if frac_good >= 0.70:
        print("\n  Status: GOOD (segments internally consistent)")
    elif frac_good >= 0.50:
        print("\n  Status: MODERATE")
    else:
        print("\n  Status: POOR (high within-segment CN variability — fragmentation suspected)")
    print()


# =============================================================================
# METRIC C — AMY1 locus
# =============================================================================

def metric_c_amy1(df, sample='chm13'):
    """Check AMY1A region for expected copy-number amplification.

    Expected CN and range are sample-specific (see _AMY1_SAMPLE_CONFIG).
    Large (>5kb) dup segments in the region represent the gene locus itself.
    Tiny (<5kb) segments are internal tandem repeats — excluded from CN estimate.
    """
    amy_cfg = _AMY1_SAMPLE_CONFIG.get(sample)
    if amy_cfg is None:
        print("=" * 60)
        print("C. AMY1 Locus Validation")
        print("=" * 60)
        print(f"  SKIP: no AMY1 GT available for sample '{sample}'")
        print()
        return

    exp_cn = amy_cfg['expected']
    cn_lo, cn_hi = amy_cfg['range']

    # Detect chromosome 1 naming convention from segment data
    seg_chroms = set(df['chrom'].unique())
    amy_chrom = next((c for c in _AMY1A_NAMES if c in seg_chroms), 'chr1')

    print("=" * 60)
    print("C. AMY1 Locus Validation")
    print("=" * 60)
    print(f"  Region: {amy_chrom}:{AMY1A_START:,}-{AMY1A_END:,}")
    print(f"  Expected CN ({amy_cfg['label']}): {cn_lo}-{cn_hi}")

    amy_df = df[
        (df['chrom'] == amy_chrom) &
        (df['start'] < AMY1A_END) &
        (df['end']   > AMY1A_START)
    ].copy()

    if len(amy_df) == 0:
        print("  No segments found in AMY1A region")
        print("  Status: FAIL (chromosome 1 missing or AMY1 region empty)")
        print()
        return

    dup_segs  = amy_df[amy_df['state'].isin(DUP_STATES)].copy()
    # Large dup segments (>5kb): represent actual gene copies, not internal tandem repeats
    large_dup = dup_segs[dup_segs['length'] > 5000]
    tiny_amp  = dup_segs[dup_segs['length'] <= 5000]

    print(f"  Total segments in region: {len(amy_df)}")
    print(f"  Dup/amp segments: {len(dup_segs)}  "
          f"(large >5kb: {len(large_dup)}, tiny ≤5kb: {len(tiny_amp)})")
    print()
    print(f"  {'State':<10}  {'CN_median':>10}  {'Length':>12}  {'Start':>12}")
    print(f"  {'-'*48}")
    for _, row in amy_df.sort_values('length', ascending=False).head(8).iterrows():
        print(f"  {row['state']:<10}  {row['cn_median']:>10.2f}  "
              f"{row['length']:>12,}  {row['start']:>12,}")

    print()
    total_dup_mb = dup_segs['length'].sum() / 1e6
    print(f"  Total dup coverage in region: {total_dup_mb:.2f} Mb")

    if len(large_dup) > 0:
        # Representative CN: weighted median of large dup segments
        rep_cn = (large_dup['cn_median'] * large_dup['length']).sum() / large_dup['length'].sum()
        max_tiny_cn = tiny_amp['cn_median'].max() if len(tiny_amp) > 0 else 0.0
        print(f"  Representative CN (large segs, length-weighted): {rep_cn:.2f}")
        if max_tiny_cn > rep_cn:
            print(f"  Note: tiny segment peak CN={max_tiny_cn:.1f} — internal tandem repeat, not gene CN")

        if cn_lo <= rep_cn <= cn_hi:
            print(f"  Status: PASS (CN={rep_cn:.1f} in expected range {cn_lo}-{cn_hi} "
                  f"for AMY1, assembly-validated CN={exp_cn})")
        elif rep_cn > cn_hi:
            print(f"  Status: MARGINAL (CN={rep_cn:.1f} > {cn_hi} — higher than expected; "
                  "possible k-mer inflation or segment boundary issue)")
        elif rep_cn >= cn_lo * 0.7:
            print(f"  Status: MARGINAL (CN={rep_cn:.1f} — below expected range {cn_lo}-{cn_hi} "
                  "but amplification detected)")
        else:
            print(f"  Status: FAIL (CN={rep_cn:.1f} < {cn_lo * 0.7:.1f} — AMY1 amplification "
                  "not detected at expected level)")
    else:
        max_cn = dup_segs['cn_median'].max() if len(dup_segs) > 0 else 0.0
        print(f"  No large (>5kb) dup segments found — only tiny repeats (max CN={max_cn:.1f})")
        print(f"  Status: FAIL (AMY1 gene locus not segmented as dup)")
    print()


# =============================================================================
# METRIC D — ONT depth correlation
# =============================================================================

def metric_d_ont_correlation(df, bam_file):
    """Pearson r between per-window ONT depth and HMM-assigned segment CN."""
    print("=" * 60)
    print("D. ONT Depth Correlation (independent depth validator)")
    print("=" * 60)

    bam_df = load_bam_features(bam_file)

    # For each window in bam_df, find its segment and get the CN
    # Build sorted segment table per chromosome for fast lookup
    seg_records = df[['chrom', 'start', 'end', 'cn_median']].copy()
    seg_records = seg_records.sort_values(['chrom', 'start']).reset_index(drop=True)

    matched_cn    = []
    matched_depth = []

    for chrom, wdf in bam_df.groupby('chrom'):
        chrom_segs = seg_records[seg_records['chrom'] == chrom]
        if len(chrom_segs) == 0:
            continue

        seg_starts = chrom_segs['start'].values
        seg_ends   = chrom_segs['end'].values
        seg_cns    = chrom_segs['cn_median'].values

        w_starts = wdf['start'].values
        w_ends   = wdf['end'].values if 'end' in wdf.columns else (wdf['start'] + 250).values
        depths   = wdf['ont_depth'].values

        # Vectorized interval assignment
        idx = np.searchsorted(seg_starts, w_starts, side='right') - 1
        idx = np.clip(idx, 0, len(chrom_segs) - 1)
        valid = (w_starts >= seg_starts[idx]) & (w_ends <= seg_ends[idx]) & (depths > 0)

        matched_cn.extend(seg_cns[idx[valid]].tolist())
        matched_depth.extend(depths[valid].tolist())

    n_matched = len(matched_cn)
    if n_matched < 100:
        print(f"  Too few matched windows ({n_matched}) — check BAM features alignment")
        print()
        return

    cns    = np.array(matched_cn)
    depths = np.array(matched_depth)

    try:
        from scipy.stats import pearsonr, spearmanr
        r_pearson,  p_pearson  = pearsonr(cns, depths)
        r_spearman, p_spearman = spearmanr(cns, depths)
    except ImportError:
        # fallback: manual Pearson
        r_pearson  = float(np.corrcoef(cns, depths)[0, 1])
        r_spearman = r_pearson
        p_pearson  = p_spearman = float('nan')

    print(f"  Matched windows:   {n_matched:,}")
    print(f"  Pearson  r:        {r_pearson:.4f}  (p={p_pearson:.2e})")
    print(f"  Spearman r:        {r_spearman:.4f}  (p={p_spearman:.2e})")

    # Per-CN-bin depth
    cn_bins = [0, 2, 4, 8, 16, 50, float('inf')]
    cn_labels = ['Neutral/Del', 'LowDup(~2)', 'HighDup(~4)', 'Amp(~8)',
                 'HighAmp(16-50)', 'HighAmp(50+)']
    print()
    print(f"  {'CN range':<18}  {'n':>7}  {'mean_cn':>9}  {'mean_ont_depth':>15}")
    print(f"  {'-'*52}")
    for lo, hi, label in zip(cn_bins[:-1], cn_bins[1:], cn_labels):
        mask = (cns >= lo) & (cns < hi)
        if mask.sum() == 0:
            continue
        print(f"  {label:<18}  {mask.sum():>7,}  "
              f"{cns[mask].mean():>9.2f}  {depths[mask].mean():>15.2f}")

    if r_pearson >= 0.85:
        print("\n  Status: GOOD (r≥0.85 — strong k-mer CN ↔ ONT depth agreement)")
    elif r_pearson >= 0.70:
        print(f"\n  Status: MODERATE (r={r_pearson:.3f}; target ≥0.85)")
    else:
        print(f"\n  Status: POOR (r={r_pearson:.3f}; systematic bias or coordinate mismatch)")
    print()


# =============================================================================
# SIZE DISTRIBUTION
# =============================================================================

def report_size_distribution(df):
    """Print segment length distribution for dup/amp states."""
    print("=" * 60)
    print("Segment Size Distribution (dup/amp states)")
    print("=" * 60)

    dup_df = df[df['state'].isin(DUP_STATES)].copy()
    if len(dup_df) == 0:
        print("  No dup segments found.")
        return

    bins   = [0, 1000, 5000, 10000, 50000, 100_000, 500_000, float('inf')]
    labels = ['<1kb', '1-5kb', '5-10kb', '10-50kb', '50-100kb', '100k-500kb', '>500kb']

    dup_df['size_bin'] = pd.cut(dup_df['length'], bins=bins, labels=labels)
    counts = dup_df['size_bin'].value_counts().reindex(labels, fill_value=0)

    total = len(dup_df)
    for label, count in counts.items():
        pct = 100.0 * count / max(1, total)
        bar = '█' * int(pct / 2)
        print(f"  {label:>10}: {count:>6,}  ({pct:5.1f}%)  {bar}")

    print()
    print(f"  Total dup segments:  {total:,}")
    print(f"  Median length:       {dup_df['length'].median():>10,.0f} bp")
    print(f"  Mean length:         {dup_df['length'].mean():>10,.0f} bp")
    print(f"  Max length:          {dup_df['length'].max():>10,.0f} bp")
    total_dup_bp = dup_df['length'].sum()
    print(f"  Total dup bases:     {total_dup_bp:>10,} bp ({total_dup_bp/1e6:.1f} Mb)")

    # Warn if too many tiny segments (fragmentation indicator)
    tiny_frac = counts.get('<1kb', 0) / max(1, total)
    if tiny_frac > 0.30:
        print(f"\n  WARNING: {tiny_frac:.0%} of segments are <1kb — "
              "consider increasing max_merge_gap or using cn-accuracy mode")
    print()


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Validate CN accuracy of KonuSeg HMM output"
    )
    parser.add_argument("--segments", required=True,
                        help="HMM segments BED file (preferably from cn-accuracy mode)")
    parser.add_argument("--bam-features",
                        help="ONT BAM features TSV from extract_bam_features.py "
                             "(enables metric D: depth correlation)")
    parser.add_argument("--output",
                        help="Write validation summary to this file (default: stdout only)")
    parser.add_argument("--sample", default="chm13",
                        help="Sample name for AMY1 locus validation "
                             "(chm13 or hg002; default: chm13)")
    args = parser.parse_args()

    print("=" * 60)
    print("KonuSeg CN Accuracy Validation")
    print("=" * 60)
    print(f"  Segments:     {args.segments}")
    if args.bam_features:
        print(f"  BAM features: {args.bam_features}")
    print()

    df = load_segments(args.segments)

    # Summary counts
    state_counts = df['state'].value_counts()
    print("Segment counts by state:")
    for state in ['Neutral', 'LowDup', 'HighDup', 'Amp', 'MedAmp', 'HighAmp', 'ExtremeAmp',
                  'HomDel', 'HetDel']:  # HomDel/HetDel kept for precision/sensitive mode outputs
        n = state_counts.get(state, 0)
        bases = df.loc[df['state'] == state, 'length'].sum()
        print(f"  {state:<10}: {n:>6,} segs  ({bases/1e6:6.2f} Mb)")
    print()

    # Metric A
    metric_a_concordance(df)

    # Metric B
    metric_b_consistency(df)

    # Metric C
    metric_c_amy1(df, sample=args.sample)

    # Metric D (optional)
    if args.bam_features:
        metric_d_ont_correlation(df, args.bam_features)

    # Size distribution
    report_size_distribution(df)

    print("=" * 60)
    print("Validation Complete")
    print("=" * 60)


if __name__ == "__main__":
    main()

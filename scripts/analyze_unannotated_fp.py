#!/usr/bin/env python3
"""
Analyze the 131 Mb of unannotated "FP" regions:
  - CN > 1.25 (CopySeg calls dup)
  - NOT in SEDEF SD
  - RepeatMasker class = nan (unannotated)

Goal: determine if these are pipeline artifacts or real signal.
"""

import sys
import numpy as np
from collections import defaultdict
import bisect

from chrom_utils import resolve_chrom, display_label

# ---------------------------------------------------------------------------
# Load SEDEF → merged intervals
# ---------------------------------------------------------------------------
def load_sedef_merged(path, target_chroms):
    raw = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue
            cols = line.rstrip('\n').split('\t')
            chrom = cols[0]
            resolved = resolve_chrom(chrom, target_chroms)
            if resolved is None:
                continue
            raw[resolved].append((int(cols[1]), int(cols[2])))
    merged = {}
    for chrom in sorted(raw):
        intervals = sorted(raw[chrom])
        result = []
        cs, ce = intervals[0]
        for s, e in intervals[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                result.append((cs, ce))
                cs, ce = s, e
        result.append((cs, ce))
        merged[chrom] = result
    return merged

def overlaps_sedef(chrom, start, end, sedef_merged):
    if chrom not in sedef_merged:
        return 0
    intervals = sedef_merged[chrom]
    starts = [iv[0] for iv in intervals]
    i = max(0, bisect.bisect_left(starts, start) - 1)
    total = 0
    while i < len(intervals):
        s, e = intervals[i]
        if s >= end:
            break
        ovlp = max(0, min(end, e) - max(start, s))
        total += ovlp
        i += 1
    return total

# ---------------------------------------------------------------------------
# Load segments
# ---------------------------------------------------------------------------
def load_segments(path):
    segs = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            seg = {
                'chrom': cols[0],
                'start': int(cols[1]),
                'end': int(cols[2]),
                'state': cols[3],
                'cn': float(cols[4]),
                'cn_mean': float(cols[5]) if len(cols) > 5 and cols[5] != '.' else None,
                'cn_std': float(cols[6]) if len(cols) > 6 and cols[6] != '.' else 0,
                'cv': float(cols[7]) if len(cols) > 7 and cols[7] != '.' else 0,
                'quality': float(cols[8]) if len(cols) > 8 and cols[8] != '.' else 1.0,
                'size': int(cols[2]) - int(cols[1]),
            }
            if len(cols) > 14:
                rc = cols[14]
                seg['repeat_class'] = None if rc in ('nan', '.', '') else rc
                seg['masked_fraction'] = float(cols[13]) if cols[13] not in ('.', '') else 0.0
            else:
                seg['repeat_class'] = None
                seg['masked_fraction'] = 0.0
            segs.append(seg)
    return segs

# ---------------------------------------------------------------------------
# Load raw windows for spot-checking
# ---------------------------------------------------------------------------
def load_windows(path, chrom_filter, start_filter, end_filter):
    """Load raw 500bp windows for a specific region."""
    windows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            chrom = cols[0]
            if chrom != chrom_filter:
                continue
            start = int(cols[1])
            end = int(cols[2])
            if start >= end_filter or end <= start_filter:
                continue
            cn = float(cols[3])
            windows.append({'start': start, 'end': end, 'cn': cn})
    return windows

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
def main():
    seg_path = "output/chm13_my_run_fl_v8_k1/segs_chm13_cnacc_w500.bed"
    sedef_path = "data/validation/chm13v2.0_SD.bed"
    windows_path = "output/chm13_my_run_fl_v8_k1/chm13_ont_cn_w500.bed"

    print("=" * 70)
    print("UNANNOTATED FP ANALYSIS")
    print("=" * 70)

    print("\n[1] Loading data...")
    segs = load_segments(seg_path)
    target_chroms = set(s['chrom'] for s in segs)
    sedef = load_sedef_merged(sedef_path, target_chroms)
    print(f"  Segments: {len(segs)}")
    print(f"  SEDEF chroms: {len(sedef)}")

    # ---------------------------------------------------------------------------
    # Classify segments
    # ---------------------------------------------------------------------------
    print("\n[2] Classifying segments...")

    categories = {
        'tp': [],           # CN>1.25, in SEDEF
        'fp_line': [],      # CN>1.25, not in SEDEF, repeat=LINE
        'fp_sine': [],      # CN>1.25, not in SEDEF, repeat=SINE
        'fp_sat': [],       # CN>1.25, not in SEDEF, repeat=Satellite
        'fp_ltr': [],       # CN>1.25, not in SEDEF, repeat=LTR
        'fp_simple': [],    # CN>1.25, not in SEDEF, repeat=Simple_repeat
        'fp_nan': [],       # CN>1.25, not in SEDEF, repeat=nan ← TARGET
        'fp_other': [],     # CN>1.25, not in SEDEF, other repeat
        'neutral': [],      # CN<=1.25
    }

    for seg in segs:
        if seg['cn'] <= 1.25:
            categories['neutral'].append(seg)
            continue

        ovlp = overlaps_sedef(seg['chrom'], seg['start'], seg['end'], sedef)
        ovlp_frac = ovlp / seg['size'] if seg['size'] > 0 else 0

        if ovlp_frac >= 0.5:
            categories['tp'].append(seg)
        else:
            rc = seg['repeat_class']
            if rc is None:
                categories['fp_nan'].append(seg)
            elif 'LINE' in rc:
                categories['fp_line'].append(seg)
            elif 'SINE' in rc:
                categories['fp_sine'].append(seg)
            elif 'Satellite' in rc or 'satellite' in rc:
                categories['fp_sat'].append(seg)
            elif 'LTR' in rc:
                categories['fp_ltr'].append(seg)
            elif 'Simple' in rc:
                categories['fp_simple'].append(seg)
            else:
                categories['fp_other'].append(seg)

    print("\n  Category breakdown:")
    for cat, segs_list in sorted(categories.items()):
        total_mb = sum(s['size'] for s in segs_list) / 1e6
        print(f"    {cat:12s}: {len(segs_list):>7,} segs  {total_mb:>8.1f} Mb")

    # ---------------------------------------------------------------------------
    # Deep analysis of fp_nan
    # ---------------------------------------------------------------------------
    fp_nan = categories['fp_nan']
    if not fp_nan:
        print("\nNo unannotated FP segments found.")
        return

    print("\n" + "=" * 70)
    print(f"[3] DEEP ANALYSIS: {len(fp_nan)} unannotated FP segments")
    print("=" * 70)

    cns = np.array([s['cn'] for s in fp_nan])
    sizes = np.array([s['size'] for s in fp_nan])
    quals = np.array([s['quality'] for s in fp_nan])
    cvs = np.array([s['cv'] for s in fp_nan])
    masked = np.array([s['masked_fraction'] for s in fp_nan])
    total_mb = sizes.sum() / 1e6

    # CN distribution
    print(f"\n  A. CN Distribution ({total_mb:.1f} Mb total)")
    print(f"     mean={cns.mean():.2f}  median={np.median(cns):.2f}  "
          f"p5={np.percentile(cns,5):.2f}  p95={np.percentile(cns,95):.2f}")

    cn_bins = [(1.25, 1.5, "borderline"), (1.5, 2.0, "low dup"),
               (2.0, 3.0, "clear dup"), (3.0, 6.0, "high dup"),
               (6.0, 50.0, "amp"), (50.0, 1e6, "extreme")]
    print(f"\n     CN Bin          Segments    Mb       %")
    print(f"     {'─'*50}")
    for lo, hi, label in cn_bins:
        mask = (cns >= lo) & (cns < hi)
        n = mask.sum()
        mb = sizes[mask].sum() / 1e6
        pct = mb / total_mb * 100 if total_mb > 0 else 0
        bar = "█" * int(pct / 2)
        print(f"     CN {lo:>5.1f}-{hi:>5.0f} ({label:>10s}): {n:>6,} segs  {mb:>7.1f} Mb  {pct:>5.1f}%  {bar}")

    # Quality distribution
    print(f"\n  B. Quality Distribution")
    print(f"     mean={quals.mean():.3f}  median={np.median(quals):.3f}  "
          f"p5={np.percentile(quals,5):.3f}")
    low_qual = (quals < 0.7).sum()
    print(f"     Low quality (<0.7): {low_qual} / {len(fp_nan)} ({low_qual/len(fp_nan)*100:.1f}%)")

    # Masked fraction
    print(f"\n  C. Masked Fraction (RepeatMasker coverage)")
    print(f"     mean={masked.mean():.3f}  median={np.median(masked):.3f}")
    print(f"     masked_frac = 0.0 (completely unannotated): "
          f"{(masked == 0).sum()} ({(masked == 0).sum()/len(fp_nan)*100:.1f}%)")
    print(f"     masked_frac > 0 but < 0.5 (partial repeat): "
          f"{((masked > 0) & (masked < 0.5)).sum()} "
          f"({((masked > 0) & (masked < 0.5)).sum()/len(fp_nan)*100:.1f}%)")
    print(f"     masked_frac >= 0.5 (majority repeat but class=nan): "
          f"{(masked >= 0.5).sum()} ({(masked >= 0.5).sum()/len(fp_nan)*100:.1f}%)")

    # Size distribution
    print(f"\n  D. Size Distribution")
    size_bins = [(0, 1000), (1000, 5000), (5000, 10000), (10000, 50000),
                 (50000, 100000), (100000, 1e9)]
    for lo, hi in size_bins:
        mask = (sizes >= lo) & (sizes < hi)
        n = mask.sum()
        mb = sizes[mask].sum() / 1e6
        label = f"{lo/1000:.0f}k-{hi/1000:.0f}k" if hi < 1e9 else f">{lo/1000:.0f}k"
        print(f"     {label:>12s}: {n:>6,} segs  {mb:>7.1f} Mb")

    # State distribution
    print(f"\n  E. State Distribution")
    state_counts = defaultdict(lambda: {'n': 0, 'bp': 0})
    for s in fp_nan:
        state_counts[s['state']]['n'] += 1
        state_counts[s['state']]['bp'] += s['size']
    for state in sorted(state_counts, key=lambda x: state_counts[x]['bp'], reverse=True):
        d = state_counts[state]
        print(f"     {state:>12s}: {d['n']:>6,} segs  {d['bp']/1e6:>7.1f} Mb")

    # Per-chromosome
    print(f"\n  F. Per-Chromosome Distribution (top 5)")
    chr_bp = defaultdict(int)
    for s in fp_nan:
        chr_bp[display_label(s['chrom'])] += s['size']
    for chrom, bp in sorted(chr_bp.items(), key=lambda x: -x[1])[:5]:
        print(f"     {chrom:>6s}: {bp/1e6:>7.1f} Mb")

    # ---------------------------------------------------------------------------
    # Borderline analysis (CN 1.25-1.5)
    # ---------------------------------------------------------------------------
    borderline = [s for s in fp_nan if 1.25 <= s['cn'] < 1.5]
    if borderline:
        print(f"\n  G. BORDERLINE Segments (CN 1.25-1.5): {len(borderline)} segs, "
              f"{sum(s['size'] for s in borderline)/1e6:.1f} Mb")
        bl_masked = np.array([s['masked_fraction'] for s in borderline])
        print(f"     masked_frac: mean={bl_masked.mean():.3f}  median={np.median(bl_masked):.3f}")
        print(f"     completely unannotated (mf=0): "
              f"{(bl_masked == 0).sum()} ({(bl_masked == 0).sum()/len(borderline)*100:.1f}%)")

    # ---------------------------------------------------------------------------
    # Clearly elevated analysis (CN >= 2.0)
    # ---------------------------------------------------------------------------
    elevated = [s for s in fp_nan if s['cn'] >= 2.0]
    if elevated:
        print(f"\n  H. CLEARLY ELEVATED Segments (CN >= 2.0): {len(elevated)} segs, "
              f"{sum(s['size'] for s in elevated)/1e6:.1f} Mb")
        el_masked = np.array([s['masked_fraction'] for s in elevated])
        el_cns = np.array([s['cn'] for s in elevated])
        print(f"     CN: mean={el_cns.mean():.2f}  median={np.median(el_cns):.2f}")
        print(f"     masked_frac: mean={el_masked.mean():.3f}  median={np.median(el_masked):.3f}")

        # Sample 10 largest clearly elevated segments
        elevated_sorted = sorted(elevated, key=lambda x: -x['size'])[:10]
        print(f"\n     Top 10 largest clearly elevated (CN>=2, nan repeat):")
        print(f"     {'Chrom':>14s}  {'Start':>12s}  {'End':>12s}  {'Size':>8s}  "
              f"{'CN':>6s}  {'Quality':>7s}  {'MaskedF':>7s}  State")
        print(f"     {'─'*85}")
        for s in elevated_sorted:
            chr_name = display_label(s['chrom'])
            print(f"     {chr_name:>14s}  {s['start']:>12,}  {s['end']:>12,}  "
                  f"{s['size']:>8,}  {s['cn']:>6.2f}  {s['quality']:>7.3f}  "
                  f"{s['masked_fraction']:>7.3f}  {s['state']}")

    # ---------------------------------------------------------------------------
    # Spot-check: look at raw windows for top elevated regions
    # ---------------------------------------------------------------------------
    print(f"\n" + "=" * 70)
    print(f"[4] RAW WINDOW SPOT-CHECK (top 3 elevated nan-FP regions)")
    print("=" * 70)

    if elevated:
        for seg in elevated_sorted[:3]:
            chr_name = display_label(seg['chrom'])
            print(f"\n  Region: {chr_name}:{seg['start']:,}-{seg['end']:,} "
                  f"({seg['size']:,} bp, CN={seg['cn']:.2f}, state={seg['state']})")

            windows = load_windows(windows_path, seg['chrom'],
                                   seg['start'], seg['end'])
            if windows:
                win_cns = [w['cn'] for w in windows]
                print(f"  Raw windows: n={len(windows)}, "
                      f"mean={np.mean(win_cns):.2f}, median={np.median(win_cns):.2f}, "
                      f"min={min(win_cns):.2f}, max={max(win_cns):.2f}, "
                      f"std={np.std(win_cns):.2f}")
                # Show distribution
                below_1 = sum(1 for c in win_cns if c < 1.0)
                near_1 = sum(1 for c in win_cns if 1.0 <= c < 1.25)
                above_125 = sum(1 for c in win_cns if 1.25 <= c < 2.0)
                above_2 = sum(1 for c in win_cns if c >= 2.0)
                print(f"  Window CN distribution: "
                      f"<1.0: {below_1} ({below_1/len(win_cns)*100:.0f}%), "
                      f"1.0-1.25: {near_1} ({near_1/len(win_cns)*100:.0f}%), "
                      f"1.25-2.0: {above_125} ({above_125/len(win_cns)*100:.0f}%), "
                      f">=2.0: {above_2} ({above_2/len(win_cns)*100:.0f}%)")
            else:
                print(f"  No raw windows found (file may need different path)")

    # ---------------------------------------------------------------------------
    # Summary verdict
    # ---------------------------------------------------------------------------
    print(f"\n" + "=" * 70)
    print(f"[5] SUMMARY")
    print("=" * 70)

    borderline_mb = sum(s['size'] for s in fp_nan if 1.25 <= s['cn'] < 1.5) / 1e6
    low_dup_mb = sum(s['size'] for s in fp_nan if 1.5 <= s['cn'] < 2.0) / 1e6
    clear_dup_mb = sum(s['size'] for s in fp_nan if s['cn'] >= 2.0) / 1e6

    print(f"\n  Total unannotated FP: {total_mb:.1f} Mb")
    print(f"  ├── CN 1.25-1.5 (borderline):    {borderline_mb:>7.1f} Mb  ← likely noise/normalization")
    print(f"  ├── CN 1.5-2.0 (low dup):        {low_dup_mb:>7.1f} Mb  ← ambiguous")
    print(f"  └── CN >= 2.0 (clearly elevated): {clear_dup_mb:>7.1f} Mb  ← likely real signal")
    print(f"\n  If borderline is noise → real FP = {low_dup_mb + clear_dup_mb:.1f} Mb")
    print(f"  Revised adjusted precision = {158.8 / (158.8 + low_dup_mb + clear_dup_mb + 28.3 + 4.8) * 100:.1f}%")
    print(f"  If borderline + low_dup are noise → real FP = {clear_dup_mb:.1f} Mb")
    print(f"  Best-case adjusted precision = {158.8 / (158.8 + clear_dup_mb + 28.3 + 4.8) * 100:.1f}%")

if __name__ == "__main__":
    main()

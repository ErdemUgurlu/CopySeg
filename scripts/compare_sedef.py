#!/usr/bin/env python3
"""
compare_sedef.py — Compare KonuSeg CN segments with SEDEF segmental duplication calls.

Metrics:
  1. Base-pair overlap (sensitivity / precision / F1) at CN>1.25 threshold
  2. CN distribution in SD vs non-SD regions
  3. State distribution in SD regions
  4. Per-chromosome breakdown
  5. Size-stratified SD detection rate

Usage:
  python3 scripts/compare_sedef.py \
      --segments output/chm13_my_run_fl_v8_k1/segs_chm13_cnacc_w500.bed \
      --sedef data/validation/chm13v2.0_SD.bed \
      --output output/chm13_my_run_fl_v8_k1/validation/sedef_comparison.md
"""

import argparse
import sys
from collections import defaultdict

import numpy as np

from chrom_utils import resolve_chrom, display_label, natural_chrom_key

# ---------------------------------------------------------------------------
# Load SEDEF BED → merged, non-redundant SD intervals per chromosome
# ---------------------------------------------------------------------------
def load_sedef_merged(path, target_chroms):
    """Load SEDEF BED9, merge overlapping entries → {chrom: [(start, end), ...]}

    Chromosome names are resolved to match `target_chroms` (the naming convention
    used in the KonuSeg segments).
    """
    raw = defaultdict(list)
    n_total = 0
    n_skipped = 0
    with open(path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue
            cols = line.rstrip('\n').split('\t')
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            resolved = resolve_chrom(chrom, target_chroms)
            n_total += 1
            if resolved is None:
                n_skipped += 1
                continue
            raw[resolved].append((start, end))
    print(f"[SEDEF] Loaded {n_total} entries, skipped {n_skipped} (chrY/chrM/unknown)")

    # Merge overlapping intervals
    merged = {}
    total_bases = 0
    total_regions = 0
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
        total_regions += len(result)
        total_bases += sum(e - s for s, e in result)
    print(f"[SEDEF] After merge: {total_regions} non-redundant regions, "
          f"{total_bases / 1e6:.1f} Mb total")
    return merged


# ---------------------------------------------------------------------------
# Load KonuSeg segments
# ---------------------------------------------------------------------------
def load_segments(path):
    """Load KonuSeg segments BED → list of dicts."""
    segments = []
    with open(path) as f:
        header = None
        for line in f:
            if line.startswith('#'):
                header = line.lstrip('#').strip().split('\t')
                continue
            cols = line.rstrip('\n').split('\t')
            seg = {
                'chrom': cols[0],
                'start': int(cols[1]),
                'end':   int(cols[2]),
                'state': cols[3],
                'cn':    float(cols[4]),   # cn_median
                'size':  int(cols[2]) - int(cols[1]),
            }
            if len(cols) > 14:
                seg['repeat_class'] = cols[14]
                seg['masked_fraction'] = float(cols[13]) if cols[13] != '.' else 0.0
            segments.append(seg)
    print(f"[SEGS] Loaded {len(segments)} segments")
    return segments


# ---------------------------------------------------------------------------
# Interval overlap helpers
# ---------------------------------------------------------------------------
def overlap_bp(s1, e1, s2, e2):
    """Base-pair overlap between two intervals."""
    return max(0, min(e1, e2) - max(s1, s2))


def query_overlap(intervals_sorted, qstart, qend):
    """Return total overlap (bp) of query with sorted interval list (binary search)."""
    import bisect
    # Find first interval that could overlap
    # intervals_sorted = [(start, end), ...] sorted by start
    idx = bisect.bisect_right([iv[0] for iv in intervals_sorted], qend) - 1
    total = 0
    # Scan backwards from idx
    i = max(0, bisect.bisect_left([iv[0] for iv in intervals_sorted], qstart) - 1)
    while i < len(intervals_sorted):
        s, e = intervals_sorted[i]
        if s >= qend:
            break
        total += overlap_bp(qstart, qend, s, e)
        i += 1
    return total


# ---------------------------------------------------------------------------
# Main comparison
# ---------------------------------------------------------------------------
def compare(sedef_merged, segments, output_path):
    lines = []
    def out(s=""):
        lines.append(s)
        print(s)

    out("# SEDEF vs KonuSeg Comparison Report")
    out()

    # -----------------------------------------------------------------------
    # 1. Genome-wide overlap statistics
    # -----------------------------------------------------------------------
    out("## 1. Base-pair Overlap Statistics")
    out()

    # Build SD lookup for fast query
    sd_bases_total = sum(sum(e - s for s, e in ivs) for ivs in sedef_merged.values())

    # Classify each segment
    dup_in_sd = 0        # our dup bases overlapping SD
    dup_not_in_sd = 0    # our dup bases NOT in SD
    neutral_in_sd = 0    # our neutral bases in SD (missed SDs)
    neutral_not_in_sd = 0
    cn_threshold = 1.25

    # Per-segment CN in SD overlap
    sd_seg_cns = []       # (overlap_bp, cn) for segments overlapping SD
    nonsd_seg_cns = []    # (size, cn) for segments NOT overlapping SD
    sd_states = defaultdict(int)   # state → bases in SD
    nonsd_states = defaultdict(int)

    # Per-chromosome stats
    chr_stats = defaultdict(lambda: {'sd_bases': 0, 'dup_in_sd': 0,
                                      'sd_detected': 0, 'total_dup': 0})

    for seg in segments:
        chrom = seg['chrom']
        sd_intervals = sedef_merged.get(chrom, [])
        ovlp = query_overlap(sd_intervals, seg['start'], seg['end'])
        non_ovlp = seg['size'] - ovlp
        is_dup = seg['cn'] > cn_threshold

        if is_dup:
            dup_in_sd += ovlp
            dup_not_in_sd += non_ovlp
            chr_stats[chrom]['total_dup'] += seg['size']
        else:
            neutral_in_sd += ovlp
            neutral_not_in_sd += non_ovlp

        if ovlp > 0:
            sd_seg_cns.append((ovlp, seg['cn']))
            sd_states[seg['state']] += ovlp
            if is_dup:
                chr_stats[chrom]['dup_in_sd'] += ovlp
                chr_stats[chrom]['sd_detected'] += ovlp
        if non_ovlp > 0:
            nonsd_seg_cns.append((non_ovlp, seg['cn']))
            nonsd_states[seg['state']] += non_ovlp

    # Compute SD bases per chromosome
    for chrom, ivs in sedef_merged.items():
        chr_stats[chrom]['sd_bases'] = sum(e - s for s, e in ivs)

    sd_detected = dup_in_sd
    sd_missed = neutral_in_sd
    total_sd_covered = sd_detected + sd_missed
    sensitivity = sd_detected / total_sd_covered * 100 if total_sd_covered > 0 else 0
    total_dup = dup_in_sd + dup_not_in_sd
    precision = dup_in_sd / total_dup * 100 if total_dup > 0 else 0
    f1 = 2 * sensitivity * precision / (sensitivity + precision) if (sensitivity + precision) > 0 else 0

    out(f"| Metric | Value |")
    out(f"|--------|-------|")
    out(f"| SEDEF SD footprint (merged) | {sd_bases_total / 1e6:.1f} Mb |")
    out(f"| KonuSeg dup bases (CN>{cn_threshold}) | {total_dup / 1e6:.1f} Mb |")
    out(f"| SD bases called as dup (TP) | {sd_detected / 1e6:.1f} Mb |")
    out(f"| SD bases called as neutral (FN) | {sd_missed / 1e6:.1f} Mb |")
    out(f"| Dup bases NOT in SD (FP*) | {dup_not_in_sd / 1e6:.1f} Mb |")
    out(f"| **Sensitivity** (SD detected) | **{sensitivity:.1f}%** |")
    out(f"| **Precision** (dup in SD) | **{precision:.1f}%** |")
    out(f"| **F1** | **{f1:.1f}%** |")
    out()
    out("> *FP note: \"FP\" bases may include real duplications not in SEDEF (e.g., tandem repeats,")
    out("> satellite arrays, rDNA — SEDEF has ≥1kb, ≥90% identity threshold). Also, k=72 detects")
    out("> only very high identity (>~98.6%) SDs, so low-identity SEDEF SDs are expected FN.*")
    out()

    # -----------------------------------------------------------------------
    # 2. CN distribution: SD vs non-SD
    # -----------------------------------------------------------------------
    out("## 2. CN Distribution in SD vs Non-SD Regions")
    out()

    if sd_seg_cns:
        sd_cns_weighted = np.array([cn for _, cn in sd_seg_cns])
        sd_weights = np.array([bp for bp, _ in sd_seg_cns], dtype=float)
        sd_weights /= sd_weights.sum()
        sd_mean = np.average(sd_cns_weighted, weights=sd_weights)
        # Weighted median
        sorted_idx = np.argsort(sd_cns_weighted)
        cumw = np.cumsum(sd_weights[sorted_idx])
        sd_median = sd_cns_weighted[sorted_idx[np.searchsorted(cumw, 0.5)]]
    else:
        sd_mean = sd_median = 0

    if nonsd_seg_cns:
        nonsd_cns_weighted = np.array([cn for _, cn in nonsd_seg_cns])
        nonsd_weights = np.array([bp for bp, _ in nonsd_seg_cns], dtype=float)
        nonsd_weights /= nonsd_weights.sum()
        nonsd_mean = np.average(nonsd_cns_weighted, weights=nonsd_weights)
        sorted_idx = np.argsort(nonsd_cns_weighted)
        cumw = np.cumsum(nonsd_weights[sorted_idx])
        nonsd_median = nonsd_cns_weighted[sorted_idx[np.searchsorted(cumw, 0.5)]]
    else:
        nonsd_mean = nonsd_median = 0

    out(f"| Region | Weighted Mean CN | Weighted Median CN |")
    out(f"|--------|-----------------|-------------------|")
    out(f"| In SEDEF SD | {sd_mean:.2f} | {sd_median:.2f} |")
    out(f"| Outside SD | {nonsd_mean:.2f} | {nonsd_median:.2f} |")
    out()

    # CN histogram bins
    cn_bins = [(0, 1.25, 'CN<1.25 (Neutral)'),
               (1.25, 3.0, 'CN 1.25-3 (LowDup)'),
               (3.0, 6.0, 'CN 3-6 (HighDup)'),
               (6.0, 12.0, 'CN 6-12 (Amp)'),
               (12.0, 22.0, 'CN 12-22 (MedAmp)'),
               (22.0, 50.0, 'CN 22-50 (HighAmp)'),
               (50.0, float('inf'), 'CN 50+ (ExtremeAmp)')]

    out("### CN Bin Distribution (base-pair weighted)")
    out()
    out(f"| CN Bin | In SD (Mb) | In SD (%) | Outside SD (Mb) | Outside SD (%) |")
    out(f"|--------|-----------|----------|----------------|---------------|")

    for lo, hi, label in cn_bins:
        sd_bp = sum(bp for bp, cn in sd_seg_cns if lo <= cn < hi)
        nonsd_bp = sum(bp for bp, cn in nonsd_seg_cns if lo <= cn < hi)
        total_sd_bp = sum(bp for bp, _ in sd_seg_cns)
        total_nonsd_bp = sum(bp for bp, _ in nonsd_seg_cns)
        sd_pct = sd_bp / total_sd_bp * 100 if total_sd_bp > 0 else 0
        nonsd_pct = nonsd_bp / total_nonsd_bp * 100 if total_nonsd_bp > 0 else 0
        out(f"| {label} | {sd_bp/1e6:.1f} | {sd_pct:.1f}% | {nonsd_bp/1e6:.1f} | {nonsd_pct:.1f}% |")
    out()

    # -----------------------------------------------------------------------
    # 3. State distribution in SD regions
    # -----------------------------------------------------------------------
    out("## 3. State Distribution in SEDEF SD Regions")
    out()
    out(f"| State | Bases in SD (Mb) | % of SD | Bases outside SD (Mb) | % outside |")
    out(f"|-------|-----------------|---------|----------------------|----------|")

    total_sd_state = sum(sd_states.values())
    total_nonsd_state = sum(nonsd_states.values())
    all_states = sorted(set(list(sd_states.keys()) + list(nonsd_states.keys())),
                       key=lambda s: -sd_states.get(s, 0))
    for state in all_states:
        sd_bp = sd_states.get(state, 0)
        nonsd_bp = nonsd_states.get(state, 0)
        sd_pct = sd_bp / total_sd_state * 100 if total_sd_state > 0 else 0
        nonsd_pct = nonsd_bp / total_nonsd_state * 100 if total_nonsd_state > 0 else 0
        out(f"| {state} | {sd_bp/1e6:.1f} | {sd_pct:.1f}% | {nonsd_bp/1e6:.1f} | {nonsd_pct:.1f}% |")
    out()

    # -----------------------------------------------------------------------
    # 4. Per-chromosome breakdown
    # -----------------------------------------------------------------------
    out("## 4. Per-Chromosome Breakdown")
    out()
    out(f"| Chr | SD (Mb) | Dup called (Mb) | SD detected (Mb) | Sensitivity |")
    out(f"|-----|---------|-----------------|-----------------|-------------|")

    for chrom in sorted(chr_stats, key=natural_chrom_key):
        cs = chr_stats[chrom]
        chr_label = display_label(chrom)
        sd_mb = cs['sd_bases'] / 1e6
        dup_mb = cs['total_dup'] / 1e6
        det_mb = cs['sd_detected'] / 1e6
        sens = cs['sd_detected'] / cs['sd_bases'] * 100 if cs['sd_bases'] > 0 else 0
        out(f"| {chr_label} | {sd_mb:.1f} | {dup_mb:.1f} | {det_mb:.1f} | {sens:.1f}% |")
    out()

    # -----------------------------------------------------------------------
    # 5. Size-stratified SD detection
    # -----------------------------------------------------------------------
    out("## 5. Size-Stratified SD Detection Rate")
    out()

    size_bins = [
        (1000, 5000, '1-5 kb'),
        (5000, 10000, '5-10 kb'),
        (10000, 50000, '10-50 kb'),
        (50000, 100000, '50-100 kb'),
        (100000, 500000, '100-500 kb'),
        (500000, float('inf'), '500+ kb'),
    ]

    out(f"| SD Size | Count | Detected (≥50% overlap dup) | Detection Rate |")
    out(f"|---------|-------|---------------------------|---------------|")

    for lo, hi, label in size_bins:
        count = 0
        detected = 0
        for chrom, intervals in sedef_merged.items():
            for s, e in intervals:
                sz = e - s
                if lo <= sz < hi:
                    count += 1
                    # Check if ≥50% of this SD is called as dup
                    sd_len = e - s
                    dup_ovlp = 0
                    for seg in segments:
                        if seg['chrom'] != chrom:
                            continue
                        if seg['cn'] > cn_threshold:
                            dup_ovlp += overlap_bp(s, e, seg['start'], seg['end'])
                    if dup_ovlp >= 0.5 * sd_len:
                        detected += 1
        rate = detected / count * 100 if count > 0 else 0
        out(f"| {label} | {count} | {detected} | {rate:.1f}% |")
    out()

    # -----------------------------------------------------------------------
    # 6. Caveats
    # -----------------------------------------------------------------------
    out("## 6. Important Caveats")
    out()
    out("1. **k-mer identity threshold**: k=72 requires ~98.6% identity over 72bp to share")
    out("   exact k-mers. SEDEF reports SDs at ≥90% identity. Low-identity SDs (90-98%)")
    out("   are expected to show CN≈1 in k-mer analysis — these are expected FN, not errors.")
    out()
    out("2. **SEDEF scope**: SEDEF detects pairwise segmental duplications (≥1kb, ≥90% identity).")
    out("   It does NOT cover: satellite arrays, rDNA repeats, or tandem repeats <1kb.")
    out("   Our pipeline detects elevated CN in all these regions → expected \"FP\" relative to SEDEF.")
    out()
    out("3. **CN magnitude**: SEDEF reports binary presence/absence of duplication.")
    out("   Our pipeline quantifies CN. A SEDEF SD with 2 copies should show CN≈2;")
    out("   one with 10 copies should show CN≈10. This comparison only checks binary overlap.")
    out()
    out("4. **CHM13 specifics**: No chrY (female). SEDEF chromosome names auto-resolved to match segments.")
    out()

    # Write output
    with open(output_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"\n[DONE] Report written to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Compare KonuSeg segments with SEDEF SD calls')
    parser.add_argument('--segments', required=True, help='KonuSeg segments BED')
    parser.add_argument('--sedef', required=True, help='SEDEF SD BED file')
    parser.add_argument('--output', required=True, help='Output report path')
    args = parser.parse_args()

    segments = load_segments(args.segments)
    target_chroms = set(seg['chrom'] for seg in segments)
    sedef_merged = load_sedef_merged(args.sedef, target_chroms)
    compare(sedef_merged, segments, args.output)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
compute_repeat_annotation.py — RepeatMasker .out → per-window repeat annotation BED
=====================================================================================

Parses a RepeatMasker .out file and annotates each 500bp genomic window with:
  - masked_fraction  : fraction of window bases covered by any RepeatMasker element
  - repeat_class     : dominant repeat class by covered bases
                       (Satellite | Simple_repeat | LINE | SINE | LTR | DNA |
                        Low_complexity | Other | None)

Input:
  RepeatMasker .out format (space-delimited, 3-line header):
    col 4  = query chromosome (chr prefix, e.g. chr1)
    col 5  = query begin (1-based, inclusive)
    col 6  = query end   (1-based, inclusive)
    col 10 = repeat class/family (e.g. "SINE/Alu", "LINE/L1", "Satellite/centr")

Output:
  10-column BED (same as KonuSeg window BED + 2 new columns):
    chrom, start, end, cn, mean_count, log_ratio, num_kmers, num_filtered,
    masked_fraction, repeat_class

Algorithm:
  Per-chromosome numpy bool mask (O(genome_size) memory).
  Processes one chromosome at a time to keep peak memory ≈ largest chromosome.
  For dominant class: tracks covered bases per class using a separate uint8 array
  with a class-ID encoding.

Usage:
  python3 scripts/compute_repeat_annotation.py \\
    --repeatmasker /home/klea/recomb/gra_bf/data/human/chm13/\\
                   chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out \\
    --windows      output/chm13_my_run/chm13_ont_cn_w500.bed \\
    --output       output/chm13_my_run/repeat_annotated_w500.bed

Author: KonuSeg Team
Date:   March 2026
"""

import argparse
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

# ============================================================================
# Constants
# ============================================================================

from chrom_utils import resolve_chrom, natural_chrom_key

# Repeat class IDs (used in uint8 mask; 0 = unmasked)
CLASS_IDS = {
    'Satellite':     1,
    'Simple_repeat': 2,
    'LINE':          3,
    'SINE':          4,
    'LTR':           5,
    'DNA':           6,
    'Low_complexity':7,
    'Other':         8,
}
ID_TO_CLASS = {v: k for k, v in CLASS_IDS.items()}
# Priority for ties: lower ID wins (Satellite > Simple_repeat > LINE > SINE …)


def normalize_class(class_family: str) -> str:
    """Map RepeatMasker class/family string to a broad category."""
    cls = class_family.split('/')[0].strip()
    if cls in ('Satellite', 'centr', 'rRNA', 'tRNA', 'snRNA', 'srpRNA', 'scRNA'):
        return 'Satellite'
    if cls in ('Simple_repeat', 'Tandem_repeat', 'TRF'):
        return 'Simple_repeat'
    if cls == 'LINE':
        return 'LINE'
    if cls == 'SINE':
        return 'SINE'
    if cls == 'LTR':
        return 'LTR'
    if cls == 'DNA':
        return 'DNA'
    if cls == 'Low_complexity':
        return 'Low_complexity'
    return 'Other'


# ============================================================================
# Step 1: Parse RepeatMasker .out
# ============================================================================

def parse_repeatmasker_out(filepath: str, target_chroms: set = None) -> dict:
    """
    Parse RepeatMasker .out file (space-delimited, 3-line header).

    Format:
      col 4: query_chr (chr prefix)   col 5: start (1-based)  col 6: end (1-based)
      col 10: class/family (e.g. "SINE/Alu", "LINE/L1", "Satellite/centr")
    Lines may have a trailing '*' (low-confidence call).

    target_chroms: set of chromosome names from the windows BED. Used to
    auto-detect naming convention (chr prefix vs CM039 accessions).

    Returns: dict {chrom: list of (start_0based, end_0based, class_id)}
    """
    records = defaultdict(list)
    n_parsed = n_skipped = 0

    print(f"[RM] Parsing .out: {filepath}")

    with open(filepath, 'r') as fh:
        for _ in range(3):          # skip 3-line header
            fh.readline()

        for line in fh:
            line = line.rstrip().rstrip('*').rstrip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 11:
                continue

            chrom_raw    = parts[4]
            class_family = parts[10]

            # Resolve chromosome name to match window BED naming
            if target_chroms is not None:
                resolved = resolve_chrom(chrom_raw, target_chroms)
                if resolved is None:
                    n_skipped += 1
                    continue
            else:
                resolved = chrom_raw

            try:
                start_0 = int(parts[5]) - 1    # 1-based → 0-based
                end_0   = int(parts[6])         # 1-based inclusive → exclusive
            except ValueError:
                n_skipped += 1
                continue

            if end_0 <= start_0:
                n_skipped += 1
                continue

            cls      = normalize_class(class_family)
            class_id = CLASS_IDS.get(cls, CLASS_IDS['Other'])
            records[resolved].append((start_0, end_0, class_id))
            n_parsed += 1

            if n_parsed % 1_000_000 == 0:
                print(f"[RM]   ...{n_parsed:,} records parsed", flush=True)

    print(f"[RM] Parsed {n_parsed:,} repeat records ({n_skipped:,} skipped/unrecognised)")
    return records


def parse_repeatmasker_bed(filepath: str, target_chroms: set = None) -> dict:
    """
    Parse RepeatMasker BED file (tab-delimited, no header, 0-based coordinates).

    Format (columns):
      0: chrom (chr prefix)  1: start (0-based)  2: end (0-based exclusive)
      3: repeat_name         4: SW_score          5: strand
      6: class               7: family            8: perc_div  9: ID

    target_chroms: set of chromosome names from the windows BED.

    Returns: dict {chrom: list of (start_0based, end_0based, class_id)}
    """
    records = defaultdict(list)
    n_parsed = n_skipped = 0

    print(f"[RM] Parsing BED: {filepath}")

    # Build resolution cache for chunk-level filtering
    _chrom_cache = {}
    def _resolve(name):
        if name not in _chrom_cache:
            if target_chroms is not None:
                _chrom_cache[name] = resolve_chrom(name, target_chroms)
            else:
                _chrom_cache[name] = name
        return _chrom_cache[name]

    CHUNK = 2_000_000
    reader = pd.read_csv(
        filepath,
        sep='\t',
        header=None,
        usecols=[0, 1, 2, 6],
        names=['chrom', 'start', 'end', 'cls'],
        dtype={'chrom': str, 'start': np.int32, 'end': np.int32, 'cls': str},
        chunksize=CHUNK,
    )

    for chunk in reader:
        chunk['_resolved'] = chunk['chrom'].map(_resolve)
        unknown = chunk['_resolved'].isna()
        n_skipped += int(unknown.sum())
        chunk = chunk[~unknown]
        if chunk.empty:
            continue

        for row in chunk.itertuples(index=False):
            if row.end <= row.start:
                n_skipped += 1
                continue
            cls      = normalize_class(row.cls)
            class_id = CLASS_IDS.get(cls, CLASS_IDS['Other'])
            records[row._resolved].append((int(row.start), int(row.end), class_id))
            n_parsed += 1

        print(f"[RM]   ...{n_parsed:,} records parsed", flush=True)

    print(f"[RM] Parsed {n_parsed:,} repeat records ({n_skipped:,} skipped/unrecognised)")
    return records


def parse_repeatmasker(filepath: str, target_chroms: set = None) -> dict:
    """Auto-detect format and parse RepeatMasker file."""
    if filepath.endswith('.bed') or filepath.endswith('.bed.gz'):
        return parse_repeatmasker_bed(filepath, target_chroms)
    return parse_repeatmasker_out(filepath, target_chroms)


# ============================================================================
# Step 2: Load windows BED
# ============================================================================

def load_windows(windows_bed: str) -> pd.DataFrame:
    """Load KonuSeg 8-column window BED into a DataFrame."""
    df = pd.read_csv(
        windows_bed,
        sep='\t',
        comment='#',
        header=None,
        names=['chrom', 'start', 'end', 'cn', 'mean_count',
               'log_ratio', 'num_kmers', 'num_filtered'],
        dtype={'chrom': str, 'start': np.int32, 'end': np.int32,
               'cn': np.float32, 'mean_count': np.float32,
               'log_ratio': np.float32, 'num_kmers': np.int32,
               'num_filtered': np.int32},
    )
    print(f"[WIN] Loaded {len(df):,} windows from {windows_bed}")
    return df


# ============================================================================
# Step 3: Annotate windows per chromosome
# ============================================================================

def annotate_windows(df: pd.DataFrame, records: dict, window_size: int) -> pd.DataFrame:
    """
    For each window, compute masked_fraction and dominant_repeat_class.

    Uses a per-chromosome uint8 numpy mask:
      0 = unmasked, 1-8 = repeat class IDs (lower ID = higher priority).
    For windows where multiple repeat classes overlap, the dominant class
    is the one covering the most bases (ties broken by class priority).
    """
    n_total = len(df)
    masked_frac  = np.zeros(n_total, dtype=np.float32)
    dom_class_id = np.zeros(n_total, dtype=np.uint8)

    for chrom in sorted(df['chrom'].unique(), key=natural_chrom_key):
        chrom_mask = df['chrom'] == chrom
        idx = np.where(chrom_mask)[0]
        if len(idx) == 0:
            continue

        chrom_recs = records.get(chrom, [])
        starts = df['start'].values[idx]
        ends   = df['end'].values[idx]

        # Chromosome length = last window end + some padding
        chrom_len = int(ends.max()) + 1

        if not chrom_recs:
            # No repeats on this chromosome — all zeros already
            continue

        # Build a per-class coverage count array for dominant-class assignment.
        # Shape: (chrom_len, n_classes) would be large; instead use a single
        # uint8 mask where we write class IDs with LOWER ID taking priority
        # (achieved by only overwriting 0 entries or entries with larger IDs).
        # A second pass gets the dominant class per window.

        # Simple approach: fill mask with class_id (highest priority = smallest id wins)
        # We write in reverse priority order so that higher-priority classes overwrite.
        mask = np.zeros(chrom_len, dtype=np.uint8)

        # Sort by class_id descending so high-priority (low ID) classes are written last
        chrom_recs_sorted = sorted(chrom_recs, key=lambda r: -r[2])
        for s, e, cid in chrom_recs_sorted:
            e = min(e, chrom_len)
            if s < e:
                mask[s:e] = cid

        # Also build a binary mask for masked_fraction
        any_mask = (mask > 0).view(np.uint8)

        print(f"[RM]   {chrom}: {len(chrom_recs):,} records, "
              f"{100 * any_mask.mean():.1f}% masked", flush=True)

        # Per-window annotation
        for local_i, global_i in enumerate(idx):
            s = int(starts[local_i])
            e = min(int(ends[local_i]), chrom_len)
            if e <= s:
                continue
            win_len = e - s
            win_mask   = any_mask[s:e]
            n_masked   = int(win_mask.sum())
            masked_frac[global_i] = n_masked / win_len

            if n_masked > 0:
                # Dominant class = most frequent non-zero value in mask[s:e]
                win_classes = mask[s:e]
                # numpy bincount for small arrays
                counts = np.bincount(win_classes, minlength=9)
                counts[0] = 0  # ignore unmasked
                dom_class_id[global_i] = int(counts.argmax())

    df = df.copy()
    df['masked_fraction'] = np.round(masked_frac, 4)
    df['repeat_class']    = [ID_TO_CLASS.get(int(c), 'None') if c > 0 else 'None'
                             for c in dom_class_id]
    return df


# ============================================================================
# Step 4: Write output & print summary
# ============================================================================

def write_output(df: pd.DataFrame, output_path: str):
    """Write annotated 10-column BED."""
    chrom_rank = {c: natural_chrom_key(c) for c in df['chrom'].unique()}
    df_sorted = df.sort_values(
        by=['chrom', 'start'],
        key=lambda col: col.map(chrom_rank) if col.name == 'chrom' else col
    )

    with open(output_path, 'w') as fh:
        fh.write("# KonuSeg repeat-annotated windows\n")
        fh.write("# chrom\tstart\tend\tcn\tmean_count\tlog_ratio\t"
                 "num_kmers\tnum_filtered\tmasked_fraction\trepeat_class\n")
        for row in df_sorted.itertuples(index=False):
            fh.write(f"{row.chrom}\t{row.start}\t{row.end}\t"
                     f"{row.cn}\t{row.mean_count}\t{row.log_ratio}\t"
                     f"{row.num_kmers}\t{row.num_filtered}\t"
                     f"{row.masked_fraction:.4f}\t{row.repeat_class}\n")

    print(f"[RM] Written {len(df_sorted):,} annotated windows → {output_path}")


def print_summary(df: pd.DataFrame):
    """Print repeat content summary."""
    n = len(df)
    print()
    print("=" * 60)
    print("REPEAT ANNOTATION SUMMARY")
    print("=" * 60)
    print(f"Total windows: {n:,}")
    print()

    # Masked fraction distribution
    for thr, label in [(0.0, 'any'), (0.5, '>50%'), (0.8, '>80%'), (0.95, '>95%')]:
        cnt = (df['masked_fraction'] > thr).sum()
        mb  = cnt * 500 / 1e6
        print(f"  masked_fraction > {thr:.2f}: {cnt:,} windows ({100*cnt/n:.1f}%, {mb:.0f} Mb)")

    print()
    print("By dominant repeat class:")
    by_cls = df.groupby('repeat_class').agg(
        n_windows=('chrom', 'count'),
        total_mb=('chrom', lambda x: len(x) * 500 / 1e6)
    ).sort_values('n_windows', ascending=False)

    for cls, row in by_cls.iterrows():
        pct = 100 * row['n_windows'] / n
        print(f"  {cls:<15s}: {row['n_windows']:>8,} windows  "
              f"({pct:5.1f}%,  {row['total_mb']:7.0f} Mb)")

    print()
    total_repeat_mb = (df['masked_fraction'] > 0.5).sum() * 500 / 1e6
    total_genome_mb = n * 500 / 1e6
    print(f"Total repetitive (>50% masked): {total_repeat_mb:.0f} Mb "
          f"/ {total_genome_mb:.0f} Mb = {100*total_repeat_mb/total_genome_mb:.1f}%")
    print("=" * 60)


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Annotate KonuSeg 500bp windows with RepeatMasker repeat classes"
    )
    parser.add_argument('--repeatmasker', '-r', required=True,
                        help='RepeatMasker .out file (space-delimited, 3-line header)')
    parser.add_argument('--windows', '-w', required=True,
                        help='KonuSeg 8-column 500bp window BED')
    parser.add_argument('--output', '-o', required=True,
                        help='Output 10-column annotated window BED')
    parser.add_argument('--window-size', type=int, default=500,
                        help='Window size in bp (default: 500)')
    args = parser.parse_args()

    print("=" * 60)
    print("KonuSeg — RepeatMasker Window Annotation")
    print("=" * 60)
    print(f"RepeatMasker: {args.repeatmasker}")
    print(f"Windows BED:  {args.windows}")
    print(f"Output:       {args.output}")
    print()

    # Load windows first to discover chromosome naming convention
    print()
    df = load_windows(args.windows)
    target_chroms = set(df['chrom'].unique())
    print(f"[WIN] Target chromosomes: {len(target_chroms)} ({', '.join(sorted(target_chroms, key=natural_chrom_key)[:5])}...)")

    # Parse RepeatMasker (.out or .bed — auto-detected from extension)
    print()
    records = parse_repeatmasker(args.repeatmasker, target_chroms)

    # Annotate
    print()
    print("[RM] Annotating windows (per-chromosome numpy mask)...")
    df = annotate_windows(df, records, args.window_size)

    # Summary
    print_summary(df)

    # Write
    write_output(df, args.output)

    return 0


if __name__ == '__main__':
    sys.exit(main())

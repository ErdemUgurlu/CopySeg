#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict

import numpy as np
import pandas as pd


import os, sys
_HERE = os.path.abspath(os.path.dirname(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from chrom_utils import resolve_chrom, natural_chrom_key

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


def normalize_class(class_family: str) -> str:
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



def parse_repeatmasker_out(filepath: str, target_chroms: set = None) -> dict:
    records = defaultdict(list)
    n_parsed = n_skipped = 0

    print(f"[RM] Parsing .out: {filepath}")

    with open(filepath, 'r') as fh:
        for _ in range(3):
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

            if target_chroms is not None:
                resolved = resolve_chrom(chrom_raw, target_chroms)
                if resolved is None:
                    n_skipped += 1
                    continue
            else:
                resolved = chrom_raw

            try:
                start_0 = int(parts[5]) - 1
                end_0   = int(parts[6])
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
    records = defaultdict(list)
    n_parsed = n_skipped = 0

    print(f"[RM] Parsing BED: {filepath}")

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
    if filepath.endswith('.bed') or filepath.endswith('.bed.gz'):
        return parse_repeatmasker_bed(filepath, target_chroms)
    return parse_repeatmasker_out(filepath, target_chroms)



def load_windows(windows_bed: str) -> pd.DataFrame:
    df = pd.read_csv(
        windows_bed,
        sep='\t',
        comment='#',
        header=None,
        usecols=range(8),
        names=['chrom', 'start', 'end', 'cn', 'mean_count',
               'log_ratio', 'num_kmers', 'num_filtered'],
        dtype={'chrom': str, 'start': np.int32, 'end': np.int32,
               'cn': np.float32, 'mean_count': np.float32,
               'log_ratio': np.float32, 'num_kmers': np.int32,
               'num_filtered': np.int32},
    )
    print(f"[WIN] Loaded {len(df):,} windows from {windows_bed}")
    return df



def annotate_windows(df: pd.DataFrame, records: dict, window_size: int) -> pd.DataFrame:
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

        chrom_len = int(ends.max()) + 1

        if not chrom_recs:
            continue


        mask = np.zeros(chrom_len, dtype=np.uint8)

        chrom_recs_sorted = sorted(chrom_recs, key=lambda r: -r[2])
        for s, e, cid in chrom_recs_sorted:
            e = min(e, chrom_len)
            if s < e:
                mask[s:e] = cid

        any_mask = (mask > 0).view(np.uint8)

        print(f"[RM]   {chrom}: {len(chrom_recs):,} records, "
              f"{100 * any_mask.mean():.1f}% masked", flush=True)

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
                win_classes = mask[s:e]
                counts = np.bincount(win_classes, minlength=9)
                counts[0] = 0
                dom_class_id[global_i] = int(counts.argmax())

    df = df.copy()
    df['masked_fraction'] = np.round(masked_frac, 4)
    df['repeat_class']    = [ID_TO_CLASS.get(int(c), 'None') if c > 0 else 'None'
                             for c in dom_class_id]
    return df



def write_output(df: pd.DataFrame, output_path: str):
    chrom_rank = {c: natural_chrom_key(c) for c in df['chrom'].unique()}
    df_sorted = df.sort_values(
        by=['chrom', 'start'],
        key=lambda col: col.map(chrom_rank) if col.name == 'chrom' else col
    )

    with open(output_path, 'w') as fh:
        fh.write("# CopySeg repeat-annotated windows\n")
        fh.write("# chrom\tstart\tend\tcn\tmean_count\tlog_ratio\t"
                 "num_kmers\tnum_filtered\tmasked_fraction\trepeat_class\n")
        for row in df_sorted.itertuples(index=False):
            fh.write(f"{row.chrom}\t{row.start}\t{row.end}\t"
                     f"{row.cn}\t{row.mean_count}\t{row.log_ratio}\t"
                     f"{row.num_kmers}\t{row.num_filtered}\t"
                     f"{row.masked_fraction:.4f}\t{row.repeat_class}\n")

    print(f"[RM] Written {len(df_sorted):,} annotated windows → {output_path}")


def print_summary(df: pd.DataFrame):
    n = len(df)
    print()
    print("=" * 60)
    print("REPEAT ANNOTATION SUMMARY")
    print("=" * 60)
    print(f"Total windows: {n:,}")
    print()

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



def main():
    parser = argparse.ArgumentParser(
        description="Annotate CopySeg 500bp windows with RepeatMasker repeat classes"
    )
    parser.add_argument('--repeatmasker', '-r', required=True,
                        help='RepeatMasker .out file (space-delimited, 3-line header)')
    parser.add_argument('--windows', '-w', required=True,
                        help='CopySeg 8-column 500bp window BED')
    parser.add_argument('--output', '-o', required=True,
                        help='Output 10-column annotated window BED')
    parser.add_argument('--window-size', type=int, default=500,
                        help='Window size in bp (default: 500)')
    args = parser.parse_args()

    print("=" * 60)
    print("CopySeg — RepeatMasker Window Annotation")
    print("=" * 60)
    print(f"RepeatMasker: {args.repeatmasker}")
    print(f"Windows BED:  {args.windows}")
    print(f"Output:       {args.output}")
    print()

    print()
    df = load_windows(args.windows)
    target_chroms = set(df['chrom'].unique())
    print(f"[WIN] Target chromosomes: {len(target_chroms)} ({', '.join(sorted(target_chroms, key=natural_chrom_key)[:5])}...)")

    print()
    records = parse_repeatmasker(args.repeatmasker, target_chroms)

    print()
    print("[RM] Annotating windows (per-chromosome numpy mask)...")
    df = annotate_windows(df, records, args.window_size)

    print_summary(df)

    write_output(df, args.output)

    return 0


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3

import argparse
import json
import os
import sys
from collections import defaultdict

import numpy as np

import os, sys
_HERE = os.path.abspath(os.path.dirname(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_HERE, _ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from compute_repeat_annotation import (
    parse_repeatmasker,
    CLASS_IDS,
    ID_TO_CLASS,
)
from chrom_utils import natural_chrom_key


def load_chrom_lengths_fai(fai_path: str) -> dict:
    lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip().split('\t')
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    print(f"[FAI] Loaded {len(lengths)} chromosome lengths from {fai_path}")
    return lengths


def load_chrom_lengths_fasta(fasta_path: str) -> dict:
    lengths = {}
    current_chrom = None
    current_len = 0
    print(f"[FASTA] Scanning chromosome lengths from {fasta_path} ...")

    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith('>'):
                if current_chrom is not None:
                    lengths[current_chrom] = current_len
                current_chrom = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.rstrip())
        if current_chrom is not None:
            lengths[current_chrom] = current_len

    print(f"[FASTA] Loaded {len(lengths)} chromosome lengths")
    return lengths


def load_chrom_lengths_from_records(records: dict) -> dict:
    lengths = {}
    for chrom, recs in records.items():
        max_end = max(e for _, e, _ in recs) if recs else 0
        lengths[chrom] = max_end + 10000
    print(f"[RM-infer] Inferred {len(lengths)} chromosome lengths from RM records")
    return lengths


def build_mask(records: list, chrom_len: int) -> np.ndarray:
    mask = np.zeros(chrom_len, dtype=np.uint8)

    sorted_recs = sorted(records, key=lambda r: -r[2])
    for s, e, cid in sorted_recs:
        e = min(e, chrom_len)
        if s < e:
            mask[s:e] = cid

    return mask


def erode_boundaries(mask: np.ndarray, kmer_size: int) -> np.ndarray:
    if kmer_size <= 1:
        return mask.copy()

    eroded = mask.copy()
    margin = kmer_size - 1

    is_repeat = mask > 0
    transitions = np.diff(is_repeat.astype(np.int8))
    change_idx = np.where(transitions != 0)[0]

    for idx in change_idx:
        lo = max(0, idx + 1 - margin)
        hi = min(len(mask), idx + 1 + margin)
        eroded[lo:hi] = 0

    if len(mask) > 0 and mask[0] > 0:
        eroded[:min(margin, len(mask))] = 0
    if len(mask) > 0 and mask[-1] > 0:
        eroded[max(0, len(mask) - margin):] = 0

    return eroded


def main():
    parser = argparse.ArgumentParser(
        description="Pre-compute per-chromosome RepeatMasker binary masks "
                    "for k-mer down-weighting in preprocessing"
    )
    parser.add_argument('--repeatmasker', '-r', required=True,
                        help='RepeatMasker .out or .bed file')
    parser.add_argument('--ref-fai', '-f', default=None,
                        help='Reference FASTA .fai index (for chromosome lengths)')
    parser.add_argument('--ref-fasta', default=None,
                        help='Reference FASTA file (alternative to .fai — scanned for lengths)')
    parser.add_argument('--output-dir', '-o', required=True,
                        help='Output directory for .npy mask files + manifest.json')
    parser.add_argument('--kmer-size', '-k', type=int, default=None,
                        help='K-mer size for boundary erosion. When provided, '
                             'positions within (k-1) bases of repeat/unique '
                             'boundaries are set to 0 (unique). This prevents '
                             'down-weighting of boundary k-mers that span both '
                             'repeat and unique sequence. Recommended for k<=50.')
    args = parser.parse_args()

    print("=" * 60)
    print("CopySeg — RepeatMasker Binary Mask Pre-computation")
    print("=" * 60)
    print(f"RepeatMasker: {args.repeatmasker}")
    print(f"Reference:    {args.ref_fai or args.ref_fasta or '(inferred from RM records)'}")
    print(f"Output dir:   {args.output_dir}")
    if args.kmer_size is not None:
        print(f"Boundary erosion: k={args.kmer_size} (margin={args.kmer_size - 1} bp)")
    print()

    chrom_lengths = None
    if args.ref_fai:
        chrom_lengths = load_chrom_lengths_fai(args.ref_fai)
    elif args.ref_fasta:
        chrom_lengths = load_chrom_lengths_fasta(args.ref_fasta)

    target_chroms = set(chrom_lengths.keys()) if chrom_lengths else None
    print()
    records = parse_repeatmasker(args.repeatmasker, target_chroms)

    if chrom_lengths is None:
        chrom_lengths = load_chrom_lengths_from_records(records)

    os.makedirs(args.output_dir, exist_ok=True)

    print()
    print("[MASK] Building per-chromosome binary masks...")
    ref_source = args.ref_fai or args.ref_fasta or 'inferred'
    manifest = {
        'source': os.path.abspath(args.repeatmasker),
        'ref': ref_source if ref_source == 'inferred' else os.path.abspath(ref_source),
        'class_ids': CLASS_IDS,
        'chromosomes': {},
    }

    total_bases = 0
    total_masked = 0

    for chrom in sorted(chrom_lengths.keys(), key=natural_chrom_key):
        chrom_len = chrom_lengths[chrom]
        chrom_recs = records.get(chrom, [])

        mask = build_mask(chrom_recs, chrom_len)

        if args.kmer_size is not None:
            n_before = int((mask > 0).sum())
            mask = erode_boundaries(mask, args.kmer_size)
            n_after = int((mask > 0).sum())
            n_eroded = n_before - n_after
        else:
            n_eroded = 0

        n_masked = int((mask > 0).sum())
        pct = 100 * n_masked / chrom_len if chrom_len > 0 else 0

        class_counts = {}
        for cid, cls_name in ID_TO_CLASS.items():
            cnt = int((mask == cid).sum())
            if cnt > 0:
                class_counts[cls_name] = cnt

        npy_path = os.path.join(args.output_dir, f"{chrom}.npy")
        np.save(npy_path, mask)

        manifest['chromosomes'][chrom] = {
            'length': chrom_len,
            'n_masked': n_masked,
            'pct_masked': round(pct, 2),
            'n_records': len(chrom_recs),
            'class_breakdown': class_counts,
        }

        total_bases += chrom_len
        total_masked += n_masked

        erode_str = f", {n_eroded:,} eroded" if n_eroded > 0 else ""
        print(f"[MASK]   {chrom}: {chrom_len:,} bp, "
              f"{n_masked:,} masked ({pct:.1f}%){erode_str}, "
              f"{len(chrom_recs):,} RM records → {npy_path}")

    manifest_path = os.path.join(args.output_dir, 'manifest.json')
    with open(manifest_path, 'w') as fh:
        json.dump(manifest, fh, indent=2)

    total_pct = 100 * total_masked / total_bases if total_bases > 0 else 0
    total_mb = total_bases / 1e6
    masked_mb = total_masked / 1e6
    n_chroms = len(manifest['chromosomes'])

    print()
    print("=" * 60)
    print("MASK PRE-COMPUTATION SUMMARY")
    print("=" * 60)
    print(f"Chromosomes:   {n_chroms}")
    print(f"Total bases:   {total_bases:,} ({total_mb:.0f} Mb)")
    print(f"Total masked:  {total_masked:,} ({masked_mb:.0f} Mb, {total_pct:.1f}%)")
    print(f"Manifest:      {manifest_path}")
    print(f"Output dir:    {args.output_dir}")
    print(f"Disk usage:    ~{total_bases / 1e9:.1f} GB (uint8 .npy files)")
    print("=" * 60)

    return 0


if __name__ == '__main__':
    sys.exit(main())

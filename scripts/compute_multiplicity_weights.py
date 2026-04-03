#!/usr/bin/env python3
"""
compute_multiplicity_weights.py — Per-position k-mer multiplicity weights
=========================================================================

Pre-computes correction weights for k-mer count normalization across
different k-mer sizes.  For each genomic position, the weight is:

    weight[pos] = M_ref[pos] / M_target[pos]

Where:
    M_target = reference multiplicity of target-size k-mer (e.g., k=32)
    M_ref    = reference multiplicity of reference-size k-mer (e.g., k=72)

This removes spurious k-mer sharing from repeat elements (Alu, LINE) while
preserving true segmental duplication signal.

Requires:
    - samtools (for FASTA extraction via faidx)
    - jellyfish (for reference k-mer counting)

Usage:
    # 1. Create Jellyfish DBs (one-time per reference + k)
    jellyfish count -m 32 -s 10G -C -t 16 reference.fa -o ref_k32.jf
    jellyfish count -m 72 -s 10G -C -t 16 reference.fa -o ref_k72.jf

    # 2. Compute per-position weights
    python3 scripts/compute_multiplicity_weights.py \\
        --ref-fasta reference.fa \\
        --k-target 32 --k-ref 72 \\
        --jf-target ref_k32.jf --jf-ref ref_k72.jf \\
        --output-dir output/mult_weights_k32_ref72/

Author: KonuSeg Team
Date: April 2026
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time

import numpy as np


# ============================================================================
# FASTA extraction (samtools faidx, consistent with query_kmer_cn.sh)
# ============================================================================

def get_chromosome_names(ref_fasta):
    """Get chromosome names and lengths from FASTA index."""
    fai_path = ref_fasta + ".fai"
    if not os.path.exists(fai_path):
        print(f"[FASTA] Creating index: {ref_fasta}")
        subprocess.run(["samtools", "faidx", ref_fasta], check=True)

    chroms = []
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            length = int(parts[1])
            chroms.append((name, length))
    return chroms


def extract_chromosome_sequence(ref_fasta, chrom, tmp_dir):
    """Extract full chromosome sequence using samtools faidx."""
    fa_path = os.path.join(tmp_dir, f"{chrom}.fa")
    subprocess.run(
        ["samtools", "faidx", ref_fasta, chrom],
        stdout=open(fa_path, "w"),
        check=True,
    )
    # Parse FASTA → uppercase string
    seq_parts = []
    with open(fa_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq_parts.append(line.strip().upper())
    os.remove(fa_path)
    return "".join(seq_parts)


# ============================================================================
# Jellyfish batch query
# ============================================================================

def batch_query_jellyfish(seq, k, jf_db, tmp_dir, batch_size=5_000_000):
    """Query Jellyfish DB for all k-mers at every position in seq.

    Returns: np.ndarray of uint32, length = len(seq) - k + 1.
    Positions with N-containing k-mers get count = 0.
    """
    n_positions = len(seq) - k + 1
    if n_positions <= 0:
        return np.zeros(0, dtype=np.uint32)

    counts = np.zeros(n_positions, dtype=np.uint32)

    for batch_start in range(0, n_positions, batch_size):
        batch_end = min(batch_start + batch_size, n_positions)

        # Write k-mers as FASTA
        fa_path = os.path.join(tmp_dir, "batch_kmers.fa")
        out_path = os.path.join(tmp_dir, "batch_counts.txt")
        valid_indices = []

        with open(fa_path, "w") as f:
            for i in range(batch_start, batch_end):
                kmer = seq[i : i + k]
                if "N" in kmer:
                    continue
                f.write(f">{i}\n{kmer}\n")
                valid_indices.append(i)

        if not valid_indices:
            continue

        # Query Jellyfish
        with open(out_path, "w") as out_fh:
            subprocess.run(
                ["jellyfish", "query", jf_db, "-s", fa_path],
                stdout=out_fh,
                check=True,
            )

        # Parse output: each line is "KMER_SEQ COUNT" or "ID COUNT"
        idx = 0
        with open(out_path) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    count = int(parts[-1])
                    if idx < len(valid_indices):
                        counts[valid_indices[idx]] = count
                    idx += 1

        # Cleanup temp files
        for p in (fa_path, out_path):
            if os.path.exists(p):
                os.remove(p)

    return counts


# ============================================================================
# Weight computation for one chromosome
# ============================================================================

def compute_weights_for_chrom(seq, k_target, k_ref, jf_target, jf_ref,
                              tmp_dir, batch_size):
    """Compute per-position weight = M_ref / M_target.

    Returns: np.float16 array of length len(seq) - k_target + 1.
    """
    n_target = len(seq) - k_target + 1
    n_ref = len(seq) - k_ref + 1

    print(f"  [1/3] Querying target DB (k={k_target}, {n_target:,} positions)...")
    m_target = batch_query_jellyfish(seq, k_target, jf_target, tmp_dir, batch_size)

    print(f"  [2/3] Querying reference DB (k={k_ref}, {n_ref:,} positions)...")
    m_ref_raw = batch_query_jellyfish(seq, k_ref, jf_ref, tmp_dir, batch_size)

    # Align arrays: m_ref is shorter (k_ref > k_target)
    # For positions where the k_ref-mer extends beyond chromosome end,
    # we cannot get M_ref → use 1/M_target (assume unique at k_ref).
    m_ref = np.ones(n_target, dtype=np.uint32)
    m_ref[:n_ref] = m_ref_raw

    print(f"  [3/3] Computing weights...")
    weights = np.ones(n_target, dtype=np.float32)

    # Compute weight = M_ref / M_target
    valid = m_target > 0
    weights[valid] = m_ref[valid].astype(np.float32) / m_target[valid].astype(np.float32)

    # Clamp weight to [0, 1] — M_ref should always be <= M_target
    # (a 72-mer match implies a 32-mer match, so M72 <= M32)
    weights = np.clip(weights, 0.0, 1.0)

    # Positions where M_target=0 (N-containing): weight=0 (will be skipped)
    weights[~valid] = 0.0

    return weights.astype(np.float16)


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Pre-compute per-position k-mer multiplicity correction weights",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create Jellyfish DBs (one-time):
  jellyfish count -m 32 -s 10G -C -t 16 reference.fa -o ref_k32.jf
  jellyfish count -m 72 -s 10G -C -t 16 reference.fa -o ref_k72.jf

  # Compute weights:
  python3 scripts/compute_multiplicity_weights.py \\
      --ref-fasta reference.fa \\
      --k-target 32 --k-ref 72 \\
      --jf-target ref_k32.jf --jf-ref ref_k72.jf \\
      --output-dir output/mult_weights_k32_ref72/
        """,
    )
    parser.add_argument("--ref-fasta", required=True,
                        help="Reference FASTA (must have .fai index or samtools available)")
    parser.add_argument("--k-target", type=int, required=True,
                        help="Target k-mer size (e.g., 32)")
    parser.add_argument("--k-ref", type=int, required=True,
                        help="Reference k-mer size (e.g., 72). Must be >= k-target.")
    parser.add_argument("--jf-target", required=True,
                        help="Jellyfish count DB for target k-mer size")
    parser.add_argument("--jf-ref", required=True,
                        help="Jellyfish count DB for reference k-mer size")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for weight .npy files and manifest.json")
    parser.add_argument("--batch-size", type=int, default=5_000_000,
                        help="Positions per Jellyfish query batch (default: 5M)")
    parser.add_argument("--chroms", default=None,
                        help="Comma-separated chromosome subset (default: all)")
    args = parser.parse_args()

    if args.k_ref < args.k_target:
        print("ERROR: --k-ref must be >= --k-target", file=sys.stderr)
        sys.exit(1)

    for path, label in [(args.ref_fasta, "Reference FASTA"),
                        (args.jf_target, "Target Jellyfish DB"),
                        (args.jf_ref, "Reference Jellyfish DB")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("KonuSeg — K-mer Multiplicity Weight Computation")
    print(f"  Reference FASTA: {args.ref_fasta}")
    print(f"  k_target:        {args.k_target}")
    print(f"  k_ref:           {args.k_ref}")
    print(f"  JF target:       {args.jf_target}")
    print(f"  JF ref:          {args.jf_ref}")
    print(f"  Output:          {args.output_dir}")
    print(f"  Batch size:      {args.batch_size:,}")
    print("=" * 60)
    print()

    # Get chromosome list from FASTA index
    all_chroms = get_chromosome_names(args.ref_fasta)
    print(f"[FASTA] Found {len(all_chroms)} sequences in reference")

    # Filter to requested subset
    if args.chroms:
        requested = set(args.chroms.split(","))
        all_chroms = [(n, l) for n, l in all_chroms if n in requested]
        print(f"[FASTA] Subset: {len(all_chroms)} chromosomes")

    # Skip mitochondrial
    from chrom_utils import is_mito
    all_chroms = [(n, l) for n, l in all_chroms if not is_mito(n)]

    manifest = {
        "k_target": args.k_target,
        "k_ref": args.k_ref,
        "reference": os.path.basename(args.ref_fasta),
        "ref_fasta": os.path.abspath(args.ref_fasta),
        "jf_target": os.path.abspath(args.jf_target),
        "jf_ref": os.path.abspath(args.jf_ref),
        "creation_date": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "chromosomes": {},
    }

    total_positions = 0
    total_corrected = 0

    with tempfile.TemporaryDirectory(prefix="mult_weights_") as tmp_dir:
        for i, (chrom, chrom_len) in enumerate(all_chroms, 1):
            print(f"\n[{i}/{len(all_chroms)}] Processing {chrom} "
                  f"({chrom_len:,} bp)...")
            t0 = time.time()

            seq = extract_chromosome_sequence(args.ref_fasta, chrom, tmp_dir)
            actual_len = len(seq)
            if actual_len < args.k_target:
                print(f"  SKIP: sequence too short ({actual_len} < {args.k_target})")
                continue

            weights = compute_weights_for_chrom(
                seq, args.k_target, args.k_ref,
                args.jf_target, args.jf_ref,
                tmp_dir, args.batch_size,
            )

            # Save
            npy_path = os.path.join(args.output_dir, f"{chrom}.npy")
            np.save(npy_path, weights)

            # Statistics
            n_pos = len(weights)
            w_f32 = weights.astype(np.float32)
            n_below_1 = int(np.sum(w_f32 < 0.999))
            n_zero = int(np.sum(w_f32 == 0.0))
            mean_w = float(np.mean(w_f32[w_f32 > 0])) if np.any(w_f32 > 0) else 0.0
            median_w = float(np.median(w_f32[w_f32 > 0])) if np.any(w_f32 > 0) else 0.0

            elapsed = time.time() - t0
            total_positions += n_pos
            total_corrected += n_below_1

            manifest["chromosomes"][chrom] = {
                "length": actual_len,
                "n_weights": n_pos,
                "file": f"{chrom}.npy",
                "mean_weight": round(mean_w, 4),
                "median_weight": round(median_w, 4),
                "frac_below_1": round(n_below_1 / max(n_pos, 1), 4),
                "n_zero": n_zero,
            }

            print(f"  Saved: {npy_path} ({n_pos:,} positions, "
                  f"{n_pos * 2 / 1e6:.1f} MB)")
            print(f"  Stats: mean_w={mean_w:.4f}, median_w={median_w:.4f}, "
                  f"frac_corrected={n_below_1/max(n_pos,1):.2%}, "
                  f"n_zero={n_zero:,}")
            print(f"  Time: {elapsed:.1f}s")

    # Write manifest
    manifest_path = os.path.join(args.output_dir, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\n[DONE] Manifest: {manifest_path}")

    # Summary
    frac = total_corrected / max(total_positions, 1)
    disk_mb = total_positions * 2 / 1e6
    print(f"\n{'=' * 60}")
    print(f"Summary: {len(manifest['chromosomes'])} chromosomes")
    print(f"  Total positions: {total_positions:,}")
    print(f"  Positions with weight < 1: {total_corrected:,} ({frac:.1%})")
    print(f"  Disk usage: {disk_mb:.1f} MB")
    print(f"{'=' * 60}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

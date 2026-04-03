#!/usr/bin/env python3
"""
plot_kmer_histogram.py — Raw k-mer count histogram visualization.

Samples the first N million lines from a k-mer BED file and produces
a histogram of k-mer counts, with the Gaussian peak (single-copy depth)
annotated.

Usage:
  python3 scripts/plot_kmer_histogram.py \
    --input /path/to/kmer_bloom_filter.bed \
    --output output/kmer_histogram.png

  # Zoom into biological peak region:
  python3 scripts/plot_kmer_histogram.py \
    --input /path/to/kmer.bed --output hist.png \
    --x-max 200 --sample-lines 100000000

  # Log-scale y-axis for full dynamic range:
  python3 scripts/plot_kmer_histogram.py \
    --input /path/to/kmer.bed --output hist_log.png --log-y
"""

import argparse
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # headless — no display needed
import matplotlib.pyplot as plt


# ── Histogram building (mirrors preprocess_kmer_windows.py logic) ────────

def build_histogram(filepath, max_lines=50_000_000, chunk_size=5_000_000):
    """Read the first max_lines lines and build a count histogram."""
    hist = np.zeros(65536, dtype=np.int64)

    reader = pd.read_csv(
        filepath, sep='\t', header=None, comment='#',
        usecols=[3], names=['count'],
        dtype={'count': np.float32},
        chunksize=chunk_size,
    )

    n_read = 0
    for chunk in reader:
        counts = np.clip(chunk['count'].values, 0, 65535).astype(np.int32)
        np.add.at(hist, counts, 1)
        n_read += len(chunk)
        if n_read >= max_lines:
            break

    print(f"[HIST] Read {n_read:,} lines from {filepath}")
    return hist, n_read


# ── Gaussian peak detection (same as preprocessing) ─────────────────────

def find_gaussian_peak(hist, search_min=2, search_max=500):
    """Find the biological single-copy peak."""
    h = hist
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
    elif peaks:
        peak = peaks[-1]
    else:
        peak = int(np.argmax(h[search_min:search_max + 1])) + search_min

    return max(peak, search_min), noise_tail_end


# ── Plotting ─────────────────────────────────────────────────────────────

def plot_histogram(hist, peak, noise_end, output_path,
                   x_max=500, log_y=False, n_lines=0, title=None):
    """Plot the k-mer count histogram with peak annotation."""
    x = np.arange(len(hist))
    y = hist.astype(np.float64)

    # Determine x range
    if x_max is None:
        # Auto: last bin with count > 0
        nonzero = np.nonzero(y)[0]
        x_max = int(nonzero[-1]) + 10 if len(nonzero) > 0 else 500
    x_max = min(x_max, len(hist))

    fig, ax = plt.subplots(figsize=(14, 6))

    # Bar plot
    ax.bar(x[1:x_max], y[1:x_max], width=1.0, color='#3b82f6',
           edgecolor='none', alpha=0.85, label='k-mer count frequency')

    # Noise region shading
    ax.axvspan(0, noise_end, alpha=0.10, color='red', label=f'Noise region (0–{noise_end})')

    # Peak annotation
    ax.axvline(peak, color='#ef4444', linewidth=2, linestyle='--',
               label=f'Gaussian peak = {peak}')
    ax.annotate(
        f'Peak = {peak}\n(single-copy depth)\nfreq = {hist[peak]:,}',
        xy=(peak, hist[peak]),
        xytext=(peak + max(x_max // 10, 15), hist[peak] * 0.85),
        fontsize=10, fontweight='bold', color='#ef4444',
        arrowprops=dict(arrowstyle='->', color='#ef4444', lw=1.5),
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='#ef4444', alpha=0.9),
    )

    # Multiples of peak for reference
    for mult in [2, 3, 4, 5, 7]:
        x_mult = peak * mult
        if x_mult < x_max:
            ax.axvline(x_mult, color='#9ca3af', linewidth=0.8,
                       linestyle=':', alpha=0.6)
            ax.text(x_mult, ax.get_ylim()[1] * 0.95, f'{mult}×',
                    ha='center', va='top', fontsize=8, color='#6b7280')

    if log_y:
        ax.set_yscale('log')
        ax.set_ylabel('Frequency (log scale)', fontsize=12)
    else:
        ax.set_ylabel('Frequency', fontsize=12)

    ax.set_xlabel('Raw k-mer count', fontsize=12)
    ax.set_xlim(0, x_max)

    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title(f'K-mer Count Histogram  (sampled {n_lines:,} lines)',
                     fontsize=14, fontweight='bold')

    ax.legend(loc='upper right', fontsize=9)
    ax.grid(axis='y', alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"[PLOT] Saved → {output_path}")


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot raw k-mer count histogram from BED file")
    parser.add_argument('--input', '-i', required=True,
                        help='K-mer BED file (4-col: chrom start end count)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output image path (.png or .pdf)')
    parser.add_argument('--sample-lines', type=int, default=50_000_000,
                        help='Max lines to sample (default: 50M)')
    parser.add_argument('--x-max', type=int, default=None,
                        help='Max x-axis value (default: auto)')
    parser.add_argument('--log-y', action='store_true',
                        help='Use log scale for y-axis')
    parser.add_argument('--title', default=None,
                        help='Custom plot title')
    args = parser.parse_args()

    hist, n_lines = build_histogram(args.input, max_lines=args.sample_lines)
    peak, noise_end = find_gaussian_peak(hist)
    print(f"[HIST] Gaussian peak (single-copy depth): {peak}")
    print(f"[HIST] Noise tail end: {noise_end}")
    print(f"[HIST] Histogram[{peak}] = {hist[peak]:,}")

    x_max = args.x_max
    if x_max is None:
        x_max = min(peak * 10, 500)  # default: show up to 10× peak

    plot_histogram(hist, peak, noise_end, args.output,
                   x_max=x_max, log_y=args.log_y, n_lines=n_lines,
                   title=args.title)
    return 0


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
"""CopySeg CN visualization — 4 diagnostic plots per run.

Usage:
    python3 scripts/visualize_cn_results.py \
        --windows  output/chm13_my_run_fl_v3/chm13_ont_cn_w500.bed \
        --segments output/chm13_my_run_fl_v3/segs_chm13_cnacc_w500.bed \
        --outdir   output/chm13_my_run_fl_v3/plots \
        --run-name fl_v3

Generates:
    Plot A: GT locus detail (raw windows + segment overlay) for each locus
    Plot B: Genome-wide CN histogram (preprocessing output)
    Plot C: Segment CN vs Length scatter (colored by state)
    Plot D: State distribution bar chart
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from chrom_utils import resolve_chrom, display_label

# ── GT loci (chr prefix — CHM13v2.0 coordinates) ──────────────────────────
GT_LOCI = [
    {'name': 'AMY1_Cluster', 'chrom': 'chr1', 'start': 103_504_385,
     'end': 103_513_421, 'expected_cn': 7.0, 'pad': 100_000},
    {'name': 'SMN1', 'chrom': 'chr5', 'start': 71_381_729,
     'end': 71_423_141, 'expected_cn': 2.0, 'pad': 50_000},
    {'name': 'SMN2', 'chrom': 'chr5', 'start': 70_791_126,
     'end': 70_837_821, 'expected_cn': 2.0, 'pad': 50_000},
    {'name': 'LPA_KIV2', 'chrom': 'chr6', 'start': 161_783_172,
     'end': 162_011_762, 'expected_cn': 9.0, 'pad': 50_000},
    {'name': 'TBC1D3', 'chrom': 'chr17', 'start': 39_044_723,
     'end': 39_055_625, 'expected_cn': 9.0, 'pad': 30_000},
    {'name': 'TP53_Control', 'chrom': 'chr17', 'start': 7_572_544,
     'end': 7_591_594, 'expected_cn': 1.0, 'pad': 30_000},
    {'name': 'NOTCH2NL', 'chrom': 'chr1', 'start': 145_265_708,
     'end': 145_345_897, 'expected_cn': 4.0, 'pad': 50_000},
]

# ── State colors ─────────────────────────────────────────────────────────────
STATE_COLORS = {
    'HomDel':     '#2c3e50',
    'HetDel':     '#8e44ad',
    'Neutral':    '#95a5a6',
    'LowDup':     '#3498db',
    'HighDup':    '#2ecc71',
    'Amp':        '#f39c12',
    'MedAmp':     '#e67e22',
    'HighAmp':    '#e74c3c',
    'ExtremeAmp': '#c0392b',
}

STATE_ORDER = ['HomDel', 'HetDel', 'Neutral', 'LowDup', 'HighDup',
               'Amp', 'MedAmp', 'HighAmp', 'ExtremeAmp']


def load_windows(path):
    """Load preprocessed windows BED (8-col or 10-col with repeat annotation)."""
    df = pd.read_csv(path, sep='\t', comment='#',
                     header=None, usecols=[0, 1, 2, 3],
                     names=['chrom', 'start', 'end', 'cn'],
                     dtype={'chrom': str, 'start': int, 'end': int, 'cn': float})
    return df


def load_segments(path):
    """Load segment BED (17-18 col extended PELT/HMM format)."""
    # Read with comment='#' to skip header lines starting with #
    df = pd.read_csv(path, sep='\t', comment='#', header=None,
                     dtype={0: str})
    # Assign standard column names based on actual column count
    std_cols = ['chrom', 'start', 'end', 'state', 'cn_median', 'cn_mean',
                'n_windows', 'avg_quality', 'min_quality', 'cn_std',
                'avg_repeats', 'avg_entropy', 'max_entropy',
                'masked_fraction', 'repeat_class', 'gc_bias_factor',
                'segment_iqr', 'boundary_conf']
    ncols = df.shape[1]
    df.columns = std_cols[:ncols] if ncols <= len(std_cols) else \
        std_cols + [f'col_{i}' for i in range(len(std_cols), ncols)]
    # Convert numeric columns
    for col in ['start', 'end', 'cn_median', 'cn_mean', 'cn_std',
                'n_windows', 'avg_quality']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna(subset=['start', 'end', 'cn_median'])
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['length'] = df['end'] - df['start']
    return df


# ═════════════════════════════════════════════════════════════════════════════
# PLOT A: GT Locus Detail
# ═════════════════════════════════════════════════════════════════════════════

def plot_locus_detail(windows, segments, locus, outdir, run_name,
                      data_chroms=None):
    """Plot raw window CN + segment overlay for a GT locus."""
    gt_chrom = locus['chrom']
    # Resolve GT chrom name to match the data's naming convention
    if data_chroms is not None:
        chrom = resolve_chrom(gt_chrom, data_chroms) or gt_chrom
    else:
        chrom = gt_chrom
    view_start = locus['start'] - locus['pad']
    view_end = locus['end'] + locus['pad']
    chr_label = display_label(chrom)

    # Filter windows
    wm = (windows['chrom'] == chrom) & \
         (windows['start'] >= view_start) & (windows['end'] <= view_end)
    w = windows[wm].copy()

    # Filter segments
    sm = (segments['chrom'] == chrom) & \
         (segments['end'] > view_start) & (segments['start'] < view_end)
    s = segments[sm].copy()

    if w.empty:
        print(f"  [WARN] No windows for {locus['name']} — skipping")
        return

    fig, ax = plt.subplots(figsize=(14, 5))

    # Raw windows as gray dots
    mid = (w['start'] + w['end']) / 2
    ax.scatter(mid / 1e6, w['cn'], s=3, c='#bdc3c7', alpha=0.5,
              edgecolors='none', label='Window CN', zorder=1)

    # Segments as colored bars
    for _, seg in s.iterrows():
        color = STATE_COLORS.get(seg['state'], '#95a5a6')
        seg_start = max(seg['start'], view_start)
        seg_end = min(seg['end'], view_end)
        ax.plot([seg_start / 1e6, seg_end / 1e6],
                [seg['cn_median'], seg['cn_median']],
                color=color, linewidth=3, solid_capstyle='butt', zorder=3)

    # GT region shading
    ax.axvspan(locus['start'] / 1e6, locus['end'] / 1e6,
               alpha=0.08, color='blue', zorder=0)

    # Expected CN line
    ax.axhline(y=locus['expected_cn'], color='red', linestyle='--',
               linewidth=1.5, alpha=0.7, label=f"Expected CN={locus['expected_cn']}")

    # CN=1 reference
    ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)

    # Y-axis: auto-scale with cap
    cn_max = min(w['cn'].quantile(0.98) * 1.3, w['cn'].max() * 1.1)
    cn_max = max(cn_max, locus['expected_cn'] * 1.5, 5)
    ax.set_ylim(-0.5, cn_max)

    ax.set_xlim(view_start / 1e6, view_end / 1e6)
    ax.set_xlabel(f'{chr_label} position (Mb)', fontsize=11)
    ax.set_ylabel('Copy Number', fontsize=11)
    ax.set_title(f"{locus['name']} — {run_name}\n"
                 f"Expected CN={locus['expected_cn']}  |  "
                 f"{len(s)} segments in view", fontsize=12)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.15)

    fig.tight_layout()
    fname = os.path.join(outdir, f"locus_{locus['name']}_{run_name}.png")
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  Saved: {fname}")


# ═════════════════════════════════════════════════════════════════════════════
# PLOT B: Genome-wide CN Histogram
# ═════════════════════════════════════════════════════════════════════════════

def plot_cn_histogram(windows, outdir, run_name):
    """Genome-wide CN distribution histogram."""
    cn = windows['cn'].values
    cn_clipped = cn[cn < 20]  # Focus on CN 0-20 range

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: CN 0-5 detail (Neutral/LowDup range)
    ax = axes[0]
    ax.hist(cn_clipped[cn_clipped < 5], bins=200, color='#3498db',
            alpha=0.7, edgecolor='none')
    ax.axvline(x=1.0, color='red', linestyle='--', linewidth=1.5,
               label='CN=1.0 (expected Neutral)')
    ax.axvline(x=1.25, color='orange', linestyle=':', linewidth=1,
               label='CN=1.25 (Neutral/LowDup boundary)')
    ax.set_xlabel('Copy Number', fontsize=11)
    ax.set_ylabel('Window count', fontsize=11)
    ax.set_title(f'CN Distribution (0-5) — {run_name}', fontsize=12)
    ax.legend(fontsize=9)
    ax.set_xlim(0, 5)
    ax.grid(True, alpha=0.15)

    # Right: Full range (log-y)
    ax = axes[1]
    bins = np.concatenate([
        np.arange(0, 10, 0.1),
        np.arange(10, 50, 1),
        np.arange(50, 200, 5),
        np.arange(200, 1001, 50),
    ])
    ax.hist(cn, bins=bins, color='#2ecc71', alpha=0.7, edgecolor='none')
    ax.set_xscale('symlog', linthresh=5)
    ax.set_yscale('log')
    ax.set_xlabel('Copy Number', fontsize=11)
    ax.set_ylabel('Window count (log)', fontsize=11)
    ax.set_title(f'CN Distribution (full range) — {run_name}', fontsize=12)
    ax.grid(True, alpha=0.15)

    # Summary stats
    neutral_frac = np.sum((cn > 0.7) & (cn < 1.25)) / len(cn) * 100
    dup_frac = np.sum(cn >= 1.25) / len(cn) * 100
    fig.text(0.5, -0.02,
             f"Total windows: {len(cn):,}  |  "
             f"Neutral (0.7-1.25): {neutral_frac:.1f}%  |  "
             f"Dup (>=1.25): {dup_frac:.1f}%  |  "
             f"Median CN: {np.median(cn):.3f}",
             ha='center', fontsize=10, style='italic')

    fig.tight_layout()
    fname = os.path.join(outdir, f"cn_histogram_{run_name}.png")
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {fname}")


# ═════════════════════════════════════════════════════════════════════════════
# PLOT C: Segment CN vs Length Scatter
# ═════════════════════════════════════════════════════════════════════════════

def plot_cn_vs_length(segments, outdir, run_name):
    """Segment CN vs segment length, colored by state."""
    fig, ax = plt.subplots(figsize=(12, 7))

    for state in STATE_ORDER:
        sdf = segments[segments['state'] == state]
        if sdf.empty:
            continue
        color = STATE_COLORS.get(state, '#95a5a6')
        alpha = 0.15 if state == 'Neutral' else 0.4
        size = 3 if state == 'Neutral' else 8
        ax.scatter(sdf['length'] / 1e3, sdf['cn_median'],
                   s=size, c=color, alpha=alpha, edgecolors='none',
                   label=f"{state} ({len(sdf):,})")

    ax.set_xscale('log')
    ax.set_yscale('symlog', linthresh=2)
    ax.set_xlabel('Segment length (kb)', fontsize=11)
    ax.set_ylabel('CN median', fontsize=11)
    ax.set_title(f'Segment CN vs Length — {run_name}', fontsize=12)
    ax.legend(loc='upper left', fontsize=8, ncol=2, markerscale=3)
    ax.grid(True, alpha=0.15)

    # Reference lines
    ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
    ax.axhline(y=1.25, color='orange', linestyle=':', linewidth=0.8, alpha=0.3)

    fig.tight_layout()
    fname = os.path.join(outdir, f"cn_vs_length_{run_name}.png")
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  Saved: {fname}")


# ═════════════════════════════════════════════════════════════════════════════
# PLOT D: State Distribution Bar Chart
# ═════════════════════════════════════════════════════════════════════════════

def plot_state_distribution(segments, outdir, run_name):
    """State distribution: segment count + total Mb side by side."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    counts = []
    mb = []
    labels = []
    colors = []

    for state in STATE_ORDER:
        sdf = segments[segments['state'] == state]
        if sdf.empty and state in ('HomDel', 'HetDel'):
            continue  # Skip empty deletion states for cleaner plot
        counts.append(len(sdf))
        mb.append(sdf['length'].sum() / 1e6 if not sdf.empty else 0)
        labels.append(state)
        colors.append(STATE_COLORS.get(state, '#95a5a6'))

    x = np.arange(len(labels))
    width = 0.7

    # Left: Segment counts
    ax = axes[0]
    bars = ax.bar(x, counts, width, color=colors, edgecolor='white', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Segment count', fontsize=11)
    ax.set_title(f'Segments by State — {run_name}', fontsize=12)
    for bar, c in zip(bars, counts):
        if c > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f'{c:,}', ha='center', va='bottom', fontsize=8)
    ax.grid(axis='y', alpha=0.15)

    # Right: Total Mb
    ax = axes[1]
    bars = ax.bar(x, mb, width, color=colors, edgecolor='white', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Total (Mb)', fontsize=11)
    ax.set_title(f'Coverage by State — {run_name}', fontsize=12)
    for bar, m in zip(bars, mb):
        if m > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f'{m:.0f}', ha='center', va='bottom', fontsize=8)
    ax.grid(axis='y', alpha=0.15)

    fig.tight_layout()
    fname = os.path.join(outdir, f"state_distribution_{run_name}.png")
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  Saved: {fname}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='CopySeg CN visualization — diagnostic plots per run')
    parser.add_argument('--windows', required=True,
                        help='Preprocessed windows BED (chm13_ont_cn_w500.bed '
                             'or repeat_annotated_w500.bed)')
    parser.add_argument('--segments', required=True,
                        help='Segment output BED (segs_chm13_cnacc_w500.bed)')
    parser.add_argument('--outdir', required=True,
                        help='Output directory for plots')
    parser.add_argument('--run-name', required=True,
                        help='Short run identifier (e.g. fl_v3)')
    parser.add_argument('--skip-windows', action='store_true',
                        help='Skip loading windows (faster — skips Plot A and B)')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    run = args.run_name

    # Load segments (always needed)
    print(f"[VIZ] Loading segments: {args.segments}")
    segments = load_segments(args.segments)
    print(f"[VIZ] Loaded {len(segments):,} segments")

    # Plot D: State distribution (no windows needed)
    print(f"\n[Plot D] State distribution")
    plot_state_distribution(segments, args.outdir, run)

    # Plot C: CN vs Length scatter (no windows needed)
    print(f"\n[Plot C] CN vs Length scatter")
    plot_cn_vs_length(segments, args.outdir, run)

    if not args.skip_windows:
        # Load windows
        print(f"\n[VIZ] Loading windows: {args.windows}")
        windows = load_windows(args.windows)
        print(f"[VIZ] Loaded {len(windows):,} windows")

        # Plot B: CN histogram
        print(f"\n[Plot B] CN histogram")
        plot_cn_histogram(windows, args.outdir, run)

        # Plot A: GT locus details
        data_chroms = set(windows['chrom'].unique()) | set(segments['chrom'].unique())
        for locus in GT_LOCI:
            print(f"\n[Plot A] {locus['name']}")
            plot_locus_detail(windows, segments, locus, args.outdir, run,
                              data_chroms=data_chroms)
    else:
        print("\n[VIZ] --skip-windows: skipping Plot A (locus detail) and Plot B (histogram)")

    print(f"\n[VIZ] Done — all plots saved to {args.outdir}/")


if __name__ == '__main__':
    main()

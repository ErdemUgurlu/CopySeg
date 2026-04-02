#!/usr/bin/env python3
"""Plot CV-split ablation comparison: fl_v8_k1 vs fl_v9a vs fl_v9b"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Data from runs
runs = ['fl_v8_k1\n(baseline)', 'fl_v9a\n(cv_split=0.6\nmin_seg=3k)', 'fl_v9b\n(cv_split=0.6\nmin_seg=10k)']
short_runs = ['fl_v8_k1', 'fl_v9a', 'fl_v9b']

# GT loci data
gt_loci = {
    'AMY1':    {'exp': 7.0, 'tol': 0.20, 'est': [7.06, 5.68, 6.07]},
    'SMN1':    {'exp': 2.0, 'tol': 0.25, 'est': [2.02, 2.06, 2.04]},
    'SMN2':    {'exp': 2.0, 'tol': 0.25, 'est': [2.10, 2.07, 2.50]},
    'LPA':     {'exp': 9.0, 'tol': 0.35, 'est': [6.52, 6.41, 6.44]},
    'TBC1D3':  {'exp': 9.0, 'tol': 0.25, 'est': [2.09, 2.04, 2.04]},
    'TP53':    {'exp': 1.0, 'tol': 0.15, 'est': [1.02, 1.02, 1.38]},
    'NOTCH2NL':{'exp': 4.0, 'tol': 0.30, 'est': [4.00, 3.36, 3.40]},
}

# Quality metrics
cv_median = [1.824, 1.537, 2.254]
cv_good_pct = [2.0, 2.4, 1.8]  # approx for fl_v8_k1
cv_noisy_pct = [88.0, 91.4, 94.5]
concordance = [82.0, 79.9, 79.6]

# State distribution (Mb)
states_data = {
    'Neutral':    [1830, 2427, 2233],
    'HetDel':     [40, 42, 40],
    'LowDup':     [665, 146, 183],
    'HighDup':    [107, 63, 92],
    'Amp':        [91, 38, 97],
    'MedAmp':     [54, 40, 88],
    'HighAmp':    [56, 76, 90],
    'ExtremeAmp': [212, 223, 232],
}

# Segment size distribution
seg_sizes = {
    '<1kb':      [0, 0, 0],
    '1-5kb':     [0, 60829, 0],
    '5-10kb':    [0, 4000, 31897],
    '10-50kb':   [0, 5661, 16424],
    '50-100kb':  [0, 832, 830],
    '>100kb':    [0, 353, 378],
}

# Total dup bases
total_dup_mb = [1185, 586, 782]
total_dup_segs = [4958, 71675, 49529]  # approx for fl_v8_k1

fig = plt.figure(figsize=(20, 16))
fig.suptitle('CV-Split Ablation: fl_v8_k1 vs fl_v9a vs fl_v9b', fontsize=16, fontweight='bold', y=0.98)

# --- Panel 1: GT Locus Accuracy (% error) ---
ax1 = fig.add_subplot(2, 3, 1)
loci_names = list(gt_loci.keys())
x = np.arange(len(loci_names))
width = 0.25
colors = ['#2196F3', '#FF9800', '#4CAF50']

for i, (run, color) in enumerate(zip(short_runs, colors)):
    pct_errs = []
    for locus in loci_names:
        d = gt_loci[locus]
        pct_err = (d['est'][i] - d['exp']) / d['exp'] * 100
        pct_errs.append(pct_err)
    bars = ax1.bar(x + i*width, pct_errs, width, label=run, color=color, alpha=0.85)

# Add tolerance bands (show for first locus as example)
ax1.axhline(y=0, color='black', linewidth=0.5)
ax1.set_xlabel('GT Locus')
ax1.set_ylabel('% Error')
ax1.set_title('GT Locus Accuracy')
ax1.set_xticks(x + width)
ax1.set_xticklabels(loci_names, rotation=45, ha='right', fontsize=9)
ax1.legend(fontsize=8, loc='lower left')
ax1.set_ylim(-85, 45)
ax1.grid(axis='y', alpha=0.3)

# --- Panel 2: CV Distribution ---
ax2 = fig.add_subplot(2, 3, 2)
x2 = np.arange(3)
bars_cv = ax2.bar(x2, cv_median, color=colors, alpha=0.85, edgecolor='black', linewidth=0.5)
ax2.set_xticks(x2)
ax2.set_xticklabels(short_runs)
ax2.set_ylabel('CV Median (dup segments)')
ax2.set_title('Within-Segment CV')
ax2.axhline(y=0.5, color='red', linestyle='--', alpha=0.6, label='Good threshold (0.5)')
ax2.axhline(y=0.3, color='green', linestyle='--', alpha=0.6, label='Excellent threshold (0.3)')
for bar, val in zip(bars_cv, cv_median):
    ax2.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.05,
             f'{val:.3f}', ha='center', va='bottom', fontweight='bold', fontsize=11)
ax2.legend(fontsize=8)
ax2.set_ylim(0, 3)
ax2.grid(axis='y', alpha=0.3)

# --- Panel 3: State Distribution (stacked bar) ---
ax3 = fig.add_subplot(2, 3, 3)
dup_states = ['LowDup', 'HighDup', 'Amp', 'MedAmp', 'HighAmp', 'ExtremeAmp']
dup_colors = ['#64B5F6', '#42A5F5', '#1E88E5', '#1565C0', '#0D47A1', '#0A2472']
bottom = np.zeros(3)
for state, color in zip(dup_states, dup_colors):
    vals = states_data[state]
    ax3.bar(x2, vals, bottom=bottom, label=state, color=color, alpha=0.85)
    bottom += np.array(vals)

ax3.set_xticks(x2)
ax3.set_xticklabels(short_runs)
ax3.set_ylabel('Dup Bases (Mb)')
ax3.set_title('Dup State Distribution')
ax3.legend(fontsize=7, loc='upper right')
for i, total in enumerate(total_dup_mb):
    ax3.text(i, total + 20, f'{total} Mb', ha='center', fontweight='bold', fontsize=10)
ax3.set_ylim(0, 1400)
ax3.grid(axis='y', alpha=0.3)

# --- Panel 4: Concordance & GT PASS ---
ax4 = fig.add_subplot(2, 3, 4)
gt_pass = [6, 6, 6]
x4 = np.arange(3)
width4 = 0.35

bars_conc = ax4.bar(x4 - width4/2, concordance, width4, label='Concordance %', color='#7E57C2', alpha=0.85)
ax4b = ax4.twinx()
bars_gt = ax4b.bar(x4 + width4/2, gt_pass, width4, label='GT PASS (/7)', color='#26A69A', alpha=0.85)

ax4.set_xticks(x4)
ax4.set_xticklabels(short_runs)
ax4.set_ylabel('Concordance %')
ax4b.set_ylabel('GT PASS count')
ax4.set_title('Concordance & GT Accuracy')
ax4.set_ylim(70, 90)
ax4b.set_ylim(0, 7.5)
for bar, val in zip(bars_conc, concordance):
    ax4.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.3,
             f'{val}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
for bar, val in zip(bars_gt, gt_pass):
    ax4b.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.1,
              f'{val}/7', ha='center', va='bottom', fontsize=10, fontweight='bold')
ax4.legend(loc='upper left', fontsize=8)
ax4b.legend(loc='upper right', fontsize=8)
ax4.grid(axis='y', alpha=0.3)

# --- Panel 5: Segment Size Comparison ---
ax5 = fig.add_subplot(2, 3, 5)
# Show fl_v9a and fl_v9b segment size distributions
sizes_v9a = [0, 60829, 4000, 5661, 832, 353]
sizes_v9b = [0, 0, 31897, 16424, 830, 378]
size_labels = ['<1kb', '1-5kb', '5-10kb', '10-50kb', '50-100kb', '>100kb']
x5 = np.arange(len(size_labels))
ax5.bar(x5 - 0.15, sizes_v9a, 0.3, label='v9a (min=3k)', color='#FF9800', alpha=0.85)
ax5.bar(x5 + 0.15, sizes_v9b, 0.3, label='v9b (min=10k)', color='#4CAF50', alpha=0.85)
ax5.set_xticks(x5)
ax5.set_xticklabels(size_labels, rotation=30, ha='right')
ax5.set_ylabel('Number of dup segments')
ax5.set_title('Dup Segment Size Distribution')
ax5.legend(fontsize=9)
ax5.grid(axis='y', alpha=0.3)

# --- Panel 6: Per-State CV Median ---
ax6 = fig.add_subplot(2, 3, 6)
cv_states = ['LowDup', 'HighDup', 'Amp', 'MedAmp', 'HighAmp', 'ExtremeAmp']
cv_v9a = [0.939, 1.407, 1.802, 2.152, 2.206, 1.094]
cv_v9b = [1.437, 2.476, 3.624, 3.394, 2.469, 1.378]
x6 = np.arange(len(cv_states))
ax6.bar(x6 - 0.15, cv_v9a, 0.3, label='v9a (min=3k)', color='#FF9800', alpha=0.85)
ax6.bar(x6 + 0.15, cv_v9b, 0.3, label='v9b (min=10k)', color='#4CAF50', alpha=0.85)
ax6.axhline(y=0.5, color='red', linestyle='--', alpha=0.6, label='Good threshold')
ax6.set_xticks(x6)
ax6.set_xticklabels(cv_states, rotation=30, ha='right')
ax6.set_ylabel('CV Median')
ax6.set_title('Per-State CV Median')
ax6.legend(fontsize=8)
ax6.grid(axis='y', alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.96])
out_path = '/Users/erdemugurlu/Desktop/KonuSeg_Phase2/output/cv_ablation_comparison.png'
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"Plot saved: {out_path}")

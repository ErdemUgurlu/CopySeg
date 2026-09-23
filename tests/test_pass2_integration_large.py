#!/usr/bin/env python3

import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, ".."))
_SCRIPTS = os.path.join(_REPO_ROOT, "scripts", "pipeline")
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

from refine_cn_unique import (
    refine_segment, refine_all, compute_unique_genome_median,
    DEFAULT_MIN_UNIQUE_KMERS_PER_BIN, DEFAULT_MIN_TRUST_BINS,
    DEFAULT_MIN_TRUST_FRACTION, TANDEM_CLASSES,
)



WINDOW_SIZE = 500
COVERAGE = 28
RM_CLASS_WEIGHTS_STR = (
    "Satellite=0.001,Simple_repeat=0.3,LINE=0.01,SINE=0.001,"
    "LTR=0.05,DNA=0.1,Low_complexity=0.3,Other=0.05")

M_BY_CLASS = {
    "Satellite":      10000,
    "Simple_repeat":  5,
    "LINE":           100,
    "SINE":           1000,
    "LTR":            50,
    "DNA":            20,
    "Low_complexity": 5,
    "Other":          50,
}

CLASS_TO_ID = {
    "Satellite":      1,
    "Simple_repeat":  2,
    "LINE":           3,
    "SINE":           4,
    "LTR":            5,
    "DNA":            6,
    "Low_complexity": 7,
    "Other":          8,
}


GENOME_LAYOUT = {
    "chr1": [
        ("ctrl_unique_left",       300_000,   1.0, 0.00,    None,             'unique_median',         1.0),
        ("TP53_archetype",         100_000,   1.0, 0.50,    "SINE",           'unique_median',         1.0),
        ("DenseAlu_singlecopy",     80_000,   1.0, 0.80,    "SINE",           'unique_median',         1.0),
        ("AMY1_archetype",          60_000,   7.0, 0.30,    "SINE",           'unique_median',         7.0),
        ("MidSD_50Alu",             80_000,   5.0, 0.50,    "SINE",           'unique_median',         5.0),
        ("LPA_VNTR_core",           50_000,   9.0, 1.00,    "Simple_repeat",  'fallback_repeat_class', None),
        ("Satellite_centromere",    40_000,   1.0, 1.00,    "Satellite",      'fallback_repeat_class', None),
        ("LowDup_CN3",              50_000,   3.0, 0.20,    "SINE",           'unique_median',         3.0),
        ("HighCN_AMP_CN15",         40_000,  15.0, 0.40,    "SINE",           'unique_median',        15.0),
        ("ctrl_unique_right",      200_000,   1.0, 0.00,    None,             'unique_median',         1.0),
    ],
    "chr5": [
        ("ctrl_unique",            200_000,   1.0, 0.00,    None,             'unique_median',         1.0),
        ("SMN1_archetype",          50_000,   2.0, 0.50,    "SINE",           'unique_median',         2.0),
        ("SMN2_archetype",          50_000,   2.0, 0.50,    "SINE",           'unique_median',         2.0),
        ("TBC1D3_archetype",        30_000,  11.0, 0.50,    "SINE",           'unique_median',        11.0),
        ("LowComplexity_seg",       40_000,   5.0, 0.95,    "Low_complexity", 'fallback_repeat_class', None),
        ("filler",                 500_000,   1.0, 0.00,    None,             'unique_median',         1.0),
    ],
    "chr17": [
        ("ctrl_allunique",       1_000_000,   1.0, 0.00,    None,             'unique_median',         1.0),
    ],
    "chr19": [
        ("ctrl",                   100_000,   1.0, 0.00,    None,             'unique_median',         1.0),
        ("Stress_CN11_90Alu",       50_000,  11.0, 0.90,    "SINE",           'unique_median',        11.0),
        ("Stress_CN11_95Alu",       50_000,  11.0, 0.95,    "SINE",           'unique_median',        11.0),
        ("Stress_CN11_98Alu",       50_000,  11.0, 0.98,    "SINE",           'fallback_too_few_unique', None),
        ("Stress_CN1_95Alu",        50_000,   1.0, 0.95,    "SINE",           'unique_median',         1.0),
        ("filler",                 300_000,   1.0, 0.00,    None,             'unique_median',         1.0),
    ],
}



def generate_synthetic_kmer_bed(genome_layout, output_path, seed=42):
    rng = np.random.default_rng(seed)
    total_positions = 0
    with open(output_path, 'w') as fh:
        for chrom, regions in genome_layout.items():
            offset = 0
            for (label, length, cn, rm_frac, rm_class,
                 _exp_method, _exp_cn) in regions:
                m_class = M_BY_CLASS.get(rm_class, 1) if rm_class else 1
                rm_threshold = int(rm_frac * 100)
                for i in range(length):
                    pos = offset + i
                    is_repeat = rm_class is not None and (pos % 100) < rm_threshold
                    if is_repeat:
                        count = COVERAGE * cn * m_class
                    else:
                        count = COVERAGE * cn
                    noise = 1.0 + rng.uniform(-0.05, 0.05)
                    count = max(1, int(count * noise))
                    fh.write(f"{chrom}\t{pos}\t{pos+1}\t{count}\n")
                offset += length
            total_positions += offset
    print(f"[SYNTH] k-mer BED: {output_path} ({total_positions:,} positions)")
    return total_positions


def generate_rm_mask(genome_layout, output_dir, skip_chroms=()):
    os.makedirs(output_dir, exist_ok=True)
    manifest_chroms = {}
    for chrom, regions in genome_layout.items():
        if chrom in skip_chroms:
            continue
        chrom_length = sum(r[1] for r in regions)
        mask = np.zeros(chrom_length, dtype=np.uint8)
        offset = 0
        for (label, length, cn, rm_frac, rm_class, _, _) in regions:
            if rm_class is not None:
                cls_id = CLASS_TO_ID[rm_class]
                rm_threshold = int(rm_frac * 100)
                for i in range(length):
                    pos = offset + i
                    if (pos % 100) < rm_threshold:
                        mask[pos] = cls_id
            offset += length
        np.save(os.path.join(output_dir, f"{chrom}.npy"), mask)
        n_repeat = int((mask > 0).sum())
        n_unique = chrom_length - n_repeat
        print(f"[SYNTH] RM mask {chrom}.npy: {n_unique:,} unique + "
              f"{n_repeat:,} repeat ({100*n_repeat/chrom_length:.1f}%)")
        manifest_chroms[chrom] = {"length": chrom_length}

    with open(os.path.join(output_dir, "manifest.json"), 'w') as fh:
        json.dump({"chromosomes": manifest_chroms}, fh, indent=2)
    print(f"[SYNTH] RM manifest: {os.path.join(output_dir, 'manifest.json')}")
    return output_dir


def generate_segments_bed(genome_layout, windows_bed, output_path,
                          cn_source='windows'):
    df = pd.read_csv(windows_bed, sep='\t', comment='#', header=None,
                     names=['chrom', 'start', 'end', 'cn', 'mean_count',
                            'log_ratio', 'num_kmers', 'num_filtered',
                            'unique_mean_count', 'n_unique_kmers'])
    df['start'] = df['start'].astype(int)
    df['end']   = df['end'].astype(int)
    df['cn']    = pd.to_numeric(df['cn'], errors='coerce')
    df['unique_mean_count'] = pd.to_numeric(df['unique_mean_count'],
                                            errors='coerce')

    rows = []
    for chrom, regions in genome_layout.items():
        offset = 0
        for (label, length, cn, rm_frac, rm_class,
             _exp_method, _exp_cn) in regions:
            seg_start = offset
            seg_end   = offset + length
            offset = seg_end

            mask = ((df['chrom'] == chrom) &
                    (df['start'] >= seg_start) &
                    (df['end']   <= seg_end))
            seg_windows = df[mask]
            if seg_windows.empty:
                cn_p1 = cn
            else:
                cn_p1 = float(seg_windows['cn'].median())

            if cn_p1 < 0.30:   state = 'HomDel'
            elif cn_p1 < 0.70: state = 'HetDel'
            elif cn_p1 < 1.25: state = 'Neutral'
            elif cn_p1 < 3.00: state = 'LowDup'
            elif cn_p1 < 6.00: state = 'HighDup'
            elif cn_p1 < 12.0: state = 'Amp'
            elif cn_p1 < 22.0: state = 'MedAmp'
            elif cn_p1 < 50.0: state = 'HighAmp'
            else:              state = 'ExtremeAmp'

            n_win = max(1, (seg_end - seg_start) // WINDOW_SIZE)
            rc = rm_class if rm_class is not None else 'None'
            rows.append({
                'chrom':           chrom,
                'start':           seg_start,
                'end':             seg_end,
                'state':           state,
                'cn_median':       cn_p1,
                'cn_mean':         cn_p1,
                'n_windows':       n_win,
                'avg_quality':     1.0,
                'min_quality':     1.0,
                'cn_std':          0.0,
                'avg_repeats':     rm_frac,
                'avg_entropy':     0.0,
                'max_entropy':     0.0,
                'masked_fraction': rm_frac,
                'repeat_class':    rc,
                '_label':          label,
                '_true_cn':        cn,
            })
    segs = pd.DataFrame(rows)

    hdr = ('#chrom\tstart\tend\tstate\tcn_median\tcn_mean\tn_windows\t'
           'avg_quality\tmin_quality\tcn_std\tavg_repeats\t'
           'avg_entropy\tmax_entropy\tmasked_fraction\trepeat_class\n')
    cols = ['chrom', 'start', 'end', 'state', 'cn_median', 'cn_mean',
            'n_windows', 'avg_quality', 'min_quality', 'cn_std',
            'avg_repeats', 'avg_entropy', 'max_entropy',
            'masked_fraction', 'repeat_class']
    with open(output_path, 'w') as fh:
        fh.write(hdr)
        for _, r in segs.iterrows():
            line = '\t'.join(
                f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c])
                for c in cols
            )
            fh.write(line + '\n')
    print(f"[SYNTH] Segments BED: {output_path} ({len(segs):,} segments)")
    return segs



def run_preprocess(kmer_bed, rm_mask_dir, output_bed, sex='XX',
                   use_class_weights=True):
    cmd = [sys.executable, os.path.join(_SCRIPTS, 'preprocess_kmer_windows.py'),
           '--input', kmer_bed,
           '--output', output_bed,
           '--window-size', str(WINDOW_SIZE),
           '--sex', sex,
           '--bio-threshold-factor', '150.0',
           '--rm-mask-dir', rm_mask_dir]
    if use_class_weights:
        cmd.extend(['--rm-class-weights', RM_CLASS_WEIGHTS_STR])
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:]); print(r.stderr[-1000:])
        raise RuntimeError(f"preprocess_kmer_windows failed: exit {r.returncode}")
    print(f"[RUN] preprocess OK → {output_bed}")
    return r.stdout


def run_refine(seg_bed, win_bed, out_bed, sex='XX'):
    cmd = [sys.executable, os.path.join(_SCRIPTS, 'refine_cn_unique.py'),
           '--segments', seg_bed,
           '--windows', win_bed,
           '--output', out_bed,
           '--sex', sex]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:]); print(r.stderr[-1000:])
        raise RuntimeError(f"refine_cn_unique failed: exit {r.returncode}")
    print(f"[RUN] refine OK → {out_bed}")
    return r.stdout


def load_refined(path):
    with open(path) as fh:
        hdr = fh.readline().lstrip('#').rstrip('\n').split('\t')
    df = pd.read_csv(path, sep='\t', comment='#', header=None, names=hdr)
    df['start']      = df['start'].astype(int)
    df['end']        = df['end'].astype(int)
    df['cn_median']  = pd.to_numeric(df['cn_median'], errors='coerce')
    df['cn_refined'] = pd.to_numeric(df['cn_refined'], errors='coerce')
    df['n_trust_bins'] = pd.to_numeric(df['n_trust_bins'], errors='coerce').astype(int)
    df['n_total_bins'] = pd.to_numeric(df['n_total_bins'], errors='coerce').astype(int)
    return df



def test_large_multichrom_integration(tolerance=0.30):
    td = tempfile.mkdtemp(prefix="copyseg_pass2_large_")
    print(f"\n[TEST] Large multi-chrom integration. Workdir: {td}\n")
    try:
        kmer_bed = os.path.join(td, "synthetic_k32_multichrom.bed")
        mask_dir = os.path.join(td, "rm_mask")
        win_bed  = os.path.join(td, "cn_w500.bed")
        seg_bed  = os.path.join(td, "segs.bed")
        ref_bed  = os.path.join(td, "segs_refined.bed")

        generate_synthetic_kmer_bed(GENOME_LAYOUT, kmer_bed)
        generate_rm_mask(GENOME_LAYOUT, mask_dir)
        run_preprocess(kmer_bed, mask_dir, win_bed, sex='XX')

        wdf = pd.read_csv(win_bed, sep='\t', comment='#', header=None)
        assert wdf.shape[1] == 10, f"Expected 10-col output, got {wdf.shape[1]}"
        unique_col_nan = pd.to_numeric(wdf[8], errors='coerce').isna().sum()
        n_unique_pos  = (pd.to_numeric(wdf[9], errors='coerce') > 0).sum()
        print(f"[VERIFY] {wdf.shape[1]}-col output, "
              f"{n_unique_pos:,} bins with unique k-mers, "
              f"{unique_col_nan:,} NaN unique_mean_count")

        segs_truth = generate_segments_bed(GENOME_LAYOUT, win_bed, seg_bed)
        run_refine(seg_bed, win_bed, ref_bed, sex='XX')

        refined = load_refined(ref_bed)
        truth_keyed = segs_truth.set_index(['chrom', 'start'])
        refined = refined.set_index(['chrom', 'start'])
        joined = refined.join(truth_keyed[['_label', '_true_cn']], how='left')

        print(f"\n{'Region':<28} {'TrueCN':>7} {'Pass1':>8} {'Pass2':>8} "
              f"{'Method':<26} {'NTrust':>6} {'NTotal':>6} {'Verdict'}")
        print("-" * 110)

        n_pass = 0
        n_fail = 0
        for (chrom, start), row in joined.iterrows():
            label    = row['_label']
            true_cn  = float(row['_true_cn'])
            cn_p1    = float(row['cn_median'])
            cn_p2    = float(row['cn_refined'])
            method   = row['refine_method']
            n_trust  = int(row['n_trust_bins'])
            n_total  = int(row['n_total_bins'])

            exp = None
            for r in GENOME_LAYOUT[chrom]:
                if r[0] == label:
                    exp = r; break

            exp_method = exp[5]
            exp_cn     = exp[6]

            verdict_parts = []
            ok_method = (method == exp_method)
            verdict_parts.append('m_ok' if ok_method else f'm_BAD(want {exp_method})')

            if method == 'unique_median':
                if exp_cn is None:
                    raise AssertionError(f"{label}: unique_median expected but no exp_cn given")
                rel_err = abs(cn_p2 - exp_cn) / exp_cn
                ok_cn = rel_err <= tolerance
                verdict_parts.append(f'cn_{"ok" if ok_cn else "BAD"}({100*rel_err:.1f}%)')
            elif method.startswith('fallback'):
                ok_cn = abs(cn_p2 - cn_p1) < 1e-3
                verdict_parts.append(f'cn_{"ok" if ok_cn else "BAD(!=p1)"}')
            else:
                ok_cn = True

            ok = ok_method and ok_cn
            verdict = ' / '.join(verdict_parts) + (' PASS' if ok else ' FAIL')
            if ok: n_pass += 1
            else:  n_fail += 1

            print(f"{label:<28} {true_cn:>7.2f} {cn_p1:>8.3f} {cn_p2:>8.3f} "
                  f"{method:<26} {n_trust:>6d} {n_total:>6d} {verdict}")

        print("-" * 110)
        print(f"INTEGRATION: {n_pass}/{n_pass+n_fail} regions passed routing+CN tests")
        if n_fail > 0:
            raise AssertionError(f"{n_fail} regions failed — investigate above")
        return True

    finally:
        if 'KEEP_ARTIFACTS' in os.environ:
            print(f"\n[TEST] Artifacts retained: {td}")
        else:
            shutil.rmtree(td)



def _windows_synth(bins, chrom="chr1"):
    rows = []
    for b in bins:
        s = b["start"]
        rows.append({
            "chrom": chrom,
            "start": s,
            "end":   s + WINDOW_SIZE,
            "cn":    b.get("cn", 1.0),
            "mean_count": b.get("mean_count", COVERAGE),
            "log_ratio":  np.log2(max(b.get("cn", 1.0), 1e-6)),
            "num_kmers": b.get("num_kmers", 469),
            "num_filtered": b.get("num_filtered", 0),
            "unique_mean_count": b.get("unique_mean_count", np.nan),
            "n_unique_kmers": b.get("n_unique_kmers", 0),
        })
    df = pd.DataFrame(rows)
    df['start'] = df['start'].astype(int)
    df['end']   = df['end'].astype(int)
    df['unique_mean_count'] = pd.to_numeric(df['unique_mean_count'], errors='coerce')
    df['n_unique_kmers']    = df['n_unique_kmers'].astype(int)
    return df


def _segs_synth(seg_dicts, chrom="chr1"):
    rows = []
    for s in seg_dicts:
        rows.append({
            "chrom": chrom, "start": s["start"], "end": s["end"],
            "state": s.get("state", "Neutral"),
            "cn_median": s["cn_median"], "cn_mean": s.get("cn_mean", s["cn_median"]),
            "n_windows": s.get("n_windows", (s["end"] - s["start"]) // WINDOW_SIZE),
            "avg_quality": 1.0, "min_quality": 1.0, "cn_std": 0.0,
            "avg_repeats": s.get("masked_fraction", 0.0),
            "avg_entropy": 0.0, "max_entropy": 0.0,
            "masked_fraction": s.get("masked_fraction", 0.0),
            "repeat_class": s.get("repeat_class", "None"),
        })
    return pd.DataFrame(rows)


def corner_t1_exact_threshold_trust_bins():
    bg = [{"start": i * WINDOW_SIZE, "cn": 1.0,
           "unique_mean_count": float(COVERAGE), "n_unique_kmers": 234}
          for i in range(500)]

    seg_bins_10 = []
    for i in range(20):
        seg_bins_10.append({
            "start": 1_000_000 + i * WINDOW_SIZE,
            "cn": 2.0,
            "unique_mean_count": 2.0 * COVERAGE if i < 10 else COVERAGE * 1.05,
            "n_unique_kmers": 234 if i < 10 else 5,
        })

    seg_bins_9 = []
    for i in range(20):
        seg_bins_9.append({
            "start": 2_000_000 + i * WINDOW_SIZE,
            "cn": 2.0,
            "unique_mean_count": 2.0 * COVERAGE if i < 9 else COVERAGE * 1.05,
            "n_unique_kmers": 234 if i < 9 else 5,
        })

    windows = pd.concat([
        _windows_synth(bg, chrom="chrBG"),
        _windows_synth(seg_bins_10, chrom="chr1"),
        _windows_synth(seg_bins_9, chrom="chr1"),
    ], ignore_index=True)

    segs = _segs_synth([
        {"start": 1_000_000, "end": 1_000_000 + 20 * WINDOW_SIZE,
         "cn_median": 2.0, "state": "LowDup", "repeat_class": "None",
         "masked_fraction": 0.50},
        {"start": 2_000_000, "end": 2_000_000 + 20 * WINDOW_SIZE,
         "cn_median": 2.0, "state": "LowDup", "repeat_class": "None",
         "masked_fraction": 0.55},
    ])
    refined = refine_all(segs, windows, sex="XX")
    assert refined.iloc[0]['refine_method'] == 'unique_median', (
        f"10 trust bins should pass: got {refined.iloc[0]['refine_method']}")
    assert refined.iloc[1]['refine_method'] == 'fallback_too_few_unique', (
        f"9 trust bins should fall back: got {refined.iloc[1]['refine_method']}")
    print(f"  CC1 PASS: 10 bins → {refined.iloc[0]['refine_method']}, "
          f"9 bins → {refined.iloc[1]['refine_method']}")
    return True


def corner_t2_empty_repeat_class():
    bg = [{"start": i * WINDOW_SIZE, "cn": 1.0,
           "unique_mean_count": float(COVERAGE), "n_unique_kmers": 234}
          for i in range(500)]
    seg_bins = [{"start": 1_000_000 + i * WINDOW_SIZE, "cn": 3.0,
                 "unique_mean_count": 3.0 * COVERAGE,
                 "n_unique_kmers": 234} for i in range(20)]
    windows = pd.concat([
        _windows_synth(bg, chrom="chrBG"),
        _windows_synth(seg_bins, chrom="chr1"),
    ], ignore_index=True)

    segs_rows = []
    for i, rc in enumerate([None, np.nan, '', 'None']):
        segs_rows.append({
            "start": 1_000_000 + i * 5_000_000,
            "end":   1_000_000 + i * 5_000_000 + 20 * WINDOW_SIZE,
            "cn_median": 3.0, "state": "LowDup",
            "repeat_class": rc, "masked_fraction": 0.05,
        })
    segs_rows[0]['start'] = 1_000_000
    segs_rows[0]['end']   = 1_000_000 + 20 * WINDOW_SIZE
    segs = _segs_synth(segs_rows)

    refined = refine_all(segs, windows, sex="XX")
    assert refined.iloc[0]['refine_method'] == 'unique_median', (
        f"empty repeat_class with trust bins should route unique_median: "
        f"got {refined.iloc[0]['refine_method']}")
    for i in range(1, 4):
        assert refined.iloc[i]['refine_method'] == 'fallback_too_few_unique', (
            f"row {i} empty class no bins should fallback_too_few_unique, "
            f"got {refined.iloc[i]['refine_method']}")
    print(f"  CC2 PASS: empty/None/NaN repeat_class never triggers fallback_repeat_class")
    return True


def corner_t3_degenerate_cn_median_zero():
    bg = [{"start": i * WINDOW_SIZE, "cn": 1.0,
           "unique_mean_count": float(COVERAGE), "n_unique_kmers": 234}
          for i in range(500)]
    seg_bins = [{"start": 1_000_000 + i * WINDOW_SIZE, "cn": 0.0,
                 "unique_mean_count": float(COVERAGE), "n_unique_kmers": 234}
                for i in range(20)]
    windows = pd.concat([
        _windows_synth(bg, chrom="chrBG"),
        _windows_synth(seg_bins, chrom="chr1"),
    ], ignore_index=True)
    segs = _segs_synth([
        {"start": 1_000_000, "end": 1_000_000 + 20 * WINDOW_SIZE,
         "cn_median": 0.0, "state": "HomDel", "repeat_class": "None",
         "masked_fraction": 0.0},
    ])
    refined = refine_all(segs, windows, sex="XX")
    row = refined.iloc[0]
    assert row['refine_method'] == 'unique_median', row.to_dict()
    assert abs(row['cn_refined'] - 1.0) < 0.30
    assert row['state_refined'] == 'Neutral', row['state_refined']
    print(f"  CC3 PASS: cn_median=0 → cn_refined={row['cn_refined']:.3f}, "
          f"state HomDel → {row['state_refined']}")
    return True


def corner_t4_boundary_kmer_outlier():
    bg = [{"start": i * WINDOW_SIZE, "cn": 1.0,
           "unique_mean_count": float(COVERAGE), "n_unique_kmers": 234}
          for i in range(500)]
    seg_bins = []
    for i in range(20):
        umc = COVERAGE * 10 if i == 10 else float(COVERAGE)
        seg_bins.append({"start": 1_000_000 + i * WINDOW_SIZE, "cn": 1.0,
                         "unique_mean_count": umc, "n_unique_kmers": 234})
    windows = pd.concat([
        _windows_synth(bg, chrom="chrBG"),
        _windows_synth(seg_bins, chrom="chr1"),
    ], ignore_index=True)
    segs = _segs_synth([
        {"start": 1_000_000, "end": 1_000_000 + 20 * WINDOW_SIZE,
         "cn_median": 1.0, "state": "Neutral", "repeat_class": "None",
         "masked_fraction": 0.20},
    ])
    refined = refine_all(segs, windows, sex="XX")
    row = refined.iloc[0]
    assert row['refine_method'] == 'unique_median'
    assert abs(row['cn_refined'] - 1.0) < 0.05, (
        f"Median should suppress 1-bin outlier, got cn_refined={row['cn_refined']}")
    print(f"  CC4 PASS: 1-bin 10x outlier suppressed by median → "
          f"cn_refined={row['cn_refined']:.3f}")
    return True


def corner_t5_missing_chrom_rm_mask():
    layout = {
        "chr1": [("ctrl", 500_000, 1.0, 0.00, None, 'unique_median', 1.0)],
        "chr2": [("orphan", 200_000, 2.0, 0.00, None, 'fallback_too_few_unique', None)],
    }
    td = tempfile.mkdtemp(prefix="copyseg_cc5_")
    try:
        kmer_bed = os.path.join(td, "kmer.bed")
        mask_dir = os.path.join(td, "rm_mask")
        win_bed  = os.path.join(td, "cn_w500.bed")
        seg_bed  = os.path.join(td, "segs.bed")
        ref_bed  = os.path.join(td, "ref.bed")

        generate_synthetic_kmer_bed(layout, kmer_bed)
        generate_rm_mask(layout, mask_dir, skip_chroms=("chr2",))
        run_preprocess(kmer_bed, mask_dir, win_bed, sex='XX')

        wdf = pd.read_csv(win_bed, sep='\t', comment='#', header=None,
                          names=['chrom', 'start', 'end', 'cn', 'mean_count',
                                 'log_ratio', 'num_kmers', 'num_filtered',
                                 'unique_mean_count', 'n_unique_kmers'])
        wdf['n_unique_kmers'] = pd.to_numeric(wdf['n_unique_kmers'], errors='coerce')
        chr2_n_unique = wdf.loc[wdf['chrom'] == 'chr2', 'n_unique_kmers']
        assert (chr2_n_unique == 0).all(), (
            f"chr2 should have 0 unique k-mers, got: {chr2_n_unique.describe()}")

        generate_segments_bed(layout, win_bed, seg_bed)
        run_refine(seg_bed, win_bed, ref_bed, sex='XX')
        refined = load_refined(ref_bed)

        chr2_row = refined[refined['chrom'] == 'chr2'].iloc[0]
        assert chr2_row['refine_method'] == 'fallback_too_few_unique', (
            f"chr2 (no mask) should be fallback, got {chr2_row['refine_method']}")
        print(f"  CC5 PASS: chr2 with no RM mask → "
              f"{chr2_row['refine_method']}, no crash, cn_refined preserved Pass-1")
        return True
    finally:
        shutil.rmtree(td)


def corner_t6_multi_segment_per_chrom_binary_search():
    bg_bins = [{"start": i * WINDOW_SIZE, "cn": 1.0,
                "unique_mean_count": float(COVERAGE), "n_unique_kmers": 234}
               for i in range(2000)]

    n_segs     = 30
    seg_len_bp = 10 * WINDOW_SIZE
    seg_bins   = []
    expected_cn = []
    for i in range(n_segs):
        cn_i = 1.0 if (i % 2 == 0) else 4.0
        expected_cn.append(cn_i)
        for j in range(10):
            seg_bins.append({
                "start": i * seg_len_bp + j * WINDOW_SIZE,
                "cn": cn_i,
                "unique_mean_count": cn_i * COVERAGE,
                "n_unique_kmers": 234,
            })

    windows = pd.concat([
        _windows_synth(bg_bins, chrom="chrBG"),
        _windows_synth(seg_bins, chrom="chr1"),
    ], ignore_index=True)

    segs = _segs_synth([{
        "start": i * seg_len_bp,
        "end":   (i + 1) * seg_len_bp,
        "cn_median": expected_cn[i],
        "state": "Neutral" if expected_cn[i] == 1.0 else "HighDup",
        "repeat_class": "None",
        "masked_fraction": 0.0,
    } for i in range(n_segs)])

    refined = refine_all(segs, windows, sex="XX")
    n_bad = 0
    for i in range(n_segs):
        cn_ref = refined.iloc[i]['cn_refined']
        if abs(cn_ref - expected_cn[i]) > 0.10:
            n_bad += 1
    assert n_bad == 0, f"{n_bad}/{n_segs} segments assigned wrong CN (binary-search off-by-one?)"
    print(f"  CC6 PASS: {n_segs} alternating segments, all correctly indexed by binary search")
    return True


def corner_t7_xy_chrx_exclusion_real_pipeline():
    layout = {
        "chr1": [
            ("ctrl",        400_000, 1.0, 0.00, None,     'unique_median', 1.0),
            ("SD_archetype", 50_000, 6.0, 0.30, "SINE",   'unique_median', 6.0),
            ("filler",      150_000, 1.0, 0.00, None,     'unique_median', 1.0),
        ],
        "chrX": [
            ("chrX_hemizygous", 600_000, 0.5, 0.00, None, 'unique_median', None),
        ],
    }
    td = tempfile.mkdtemp(prefix="copyseg_cc7_")
    try:
        kmer_bed = os.path.join(td, "kmer.bed")
        mask_dir = os.path.join(td, "rm_mask")
        win_bed  = os.path.join(td, "cn_w500.bed")
        seg_bed  = os.path.join(td, "segs.bed")
        ref_bed  = os.path.join(td, "ref.bed")

        generate_synthetic_kmer_bed(layout, kmer_bed)
        generate_rm_mask(layout, mask_dir)
        run_preprocess(kmer_bed, mask_dir, win_bed, sex='XY')
        generate_segments_bed(layout, win_bed, seg_bed)
        run_refine(seg_bed, win_bed, ref_bed, sex='XY')

        refined = load_refined(ref_bed)
        sd_row = refined[(refined['chrom'] == 'chr1') &
                         (refined['start'] == 400_000)].iloc[0]
        assert abs(sd_row['cn_refined'] - 6.0) < 1.0, (
            f"XY mode chr1 SD: expected ~6.0, got {sd_row['cn_refined']:.3f}")
        print(f"  CC7 PASS: XY mode chr1 SD_archetype cn_refined={sd_row['cn_refined']:.3f} "
              f"(chrX hemizygous bins excluded from median as designed)")
        return True
    finally:
        shutil.rmtree(td)


def corner_t8_legacy_8col_file():
    td = tempfile.mkdtemp(prefix="copyseg_cc8_")
    try:
        win_bed = os.path.join(td, "legacy_cn_w500.bed")
        with open(win_bed, 'w') as fh:
            fh.write("# legacy 8-col format\n")
            fh.write("# chrom\tstart\tend\tcn\tmean_count\tlog_ratio\tnum_kmers\tnum_filtered\n")
            for i in range(50):
                s = i * WINDOW_SIZE
                fh.write(f"chr1\t{s}\t{s+WINDOW_SIZE}\t1.0\t{COVERAGE:.4f}\t0.0\t450\t0\n")

        seg_bed = os.path.join(td, "segs.bed")
        with open(seg_bed, 'w') as fh:
            fh.write("#chrom\tstart\tend\tstate\tcn_median\tcn_mean\tn_windows\t"
                     "avg_quality\tmin_quality\tcn_std\tavg_repeats\t"
                     "avg_entropy\tmax_entropy\tmasked_fraction\trepeat_class\n")
            fh.write("chr1\t0\t10000\tNeutral\t1.0000\t1.0000\t20\t1.0\t1.0\t0.0\t0.0\t0.0\t0.0\t0.0\tNone\n")

        ref_bed = os.path.join(td, "ref.bed")
        run_refine(seg_bed, win_bed, ref_bed)
        refined = load_refined(ref_bed)
        row = refined.iloc[0]
        assert row['refine_method'] == 'fallback_too_few_unique', row.to_dict()
        assert abs(row['cn_refined'] - 1.0) < 1e-6
        print(f"  CC8 PASS: legacy 8-col → fallback_too_few_unique, no crash")
        return True
    finally:
        shutil.rmtree(td)


def corner_t9_segment_unaligned_to_bins():
    bg = [{"start": i * WINDOW_SIZE, "cn": 1.0,
           "unique_mean_count": float(COVERAGE), "n_unique_kmers": 234}
          for i in range(500)]
    seg_bins = [{"start": 1_000_000 + i * WINDOW_SIZE, "cn": 4.0,
                 "unique_mean_count": 4.0 * COVERAGE,
                 "n_unique_kmers": 234} for i in range(15)]
    windows = pd.concat([
        _windows_synth(bg, chrom="chrBG"),
        _windows_synth(seg_bins, chrom="chr1"),
    ], ignore_index=True)
    segs = _segs_synth([
        {"start": 1_000_137, "end": 1_007_321,
         "cn_median": 4.0, "state": "HighDup",
         "repeat_class": "None", "masked_fraction": 0.0},
    ])
    refined = refine_all(segs, windows, sex="XX")
    row = refined.iloc[0]
    assert row['n_total_bins'] > 5, f"n_total too small: {row['n_total_bins']}"
    assert row['refine_method'] == 'unique_median', row.to_dict()
    assert abs(row['cn_refined'] - 4.0) < 0.10
    print(f"  CC9 PASS: unaligned segment [1_000_137, 1_007_321) → "
          f"n_total={row['n_total_bins']}, cn_refined={row['cn_refined']:.3f}")
    return True



CORNER_TESTS = [
    ("CC1 — Exact-threshold trust bins (10 vs 9)",         corner_t1_exact_threshold_trust_bins),
    ("CC2 — Empty/None/NaN repeat_class never triggers tandem fallback",
                                                            corner_t2_empty_repeat_class),
    ("CC3 — Degenerate cn_median=0, recover via unique",   corner_t3_degenerate_cn_median_zero),
    ("CC4 — 1-bin extreme outlier suppressed by median",   corner_t4_boundary_kmer_outlier),
    ("CC5 — Missing chrom RM mask file → safe fallback",   corner_t5_missing_chrom_rm_mask),
    ("CC6 — 50 alternating segments, binary-search correct", corner_t6_multi_segment_per_chrom_binary_search),
    ("CC7 — XY mode chrX hemizygous excluded from median (end-to-end)",
                                                            corner_t7_xy_chrx_exclusion_real_pipeline),
    ("CC8 — Pass-2 on legacy 8-col windows file",          corner_t8_legacy_8col_file),
    ("CC9 — Segment boundaries unaligned to bin grid",     corner_t9_segment_unaligned_to_bins),
]


def main():
    print("=" * 78)
    print("Pass-2 — LARGE INTEGRATION + CORNER-CASE SUITE")
    print("=" * 78)

    n_pass = 0
    n_fail = 0
    failures = []

    print("\n--- Stage 1: Multi-chromosome integration ---")
    try:
        test_large_multichrom_integration()
        n_pass += 1
        print("  Multi-chrom integration: PASS")
    except Exception as e:
        n_fail += 1
        failures.append(("Multi-chrom integration", e))
        print(f"  Multi-chrom integration: FAIL — {e}")
        import traceback; traceback.print_exc()

    print("\n--- Stage 2: Corner cases ---")
    for name, fn in CORNER_TESTS:
        print(f"\n[{name}]")
        try:
            fn()
            n_pass += 1
        except Exception as e:
            n_fail += 1
            failures.append((name, e))
            print(f"  FAIL: {e}")
            import traceback; traceback.print_exc()

    print()
    print("=" * 78)
    n_total = n_pass + n_fail
    print(f"FINAL: {n_pass}/{n_total} tests PASSED, {n_fail} FAILED")
    if failures:
        for n, e in failures:
            print(f"  - {n}: {e}")
    print("=" * 78)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

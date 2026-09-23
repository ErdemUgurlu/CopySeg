#!/usr/bin/env python3

import json
import os
import sys
import tempfile
import shutil

import numpy as np

CHROM = "chr1"
GENOME_SIZE = 1_000_000
WINDOW_SIZE = 500
COVERAGE = 28
ALU_MULTIPLICITY = 1000

REGIONS = [
    (0,      200_000, 1.0, 0.50, "Single-copy_50%Alu"),
    (200_000, 400_000, 2.0, 0.50, "SD_CN2_50%Alu"),
    (400_000, 500_000, 9.0, 1.00, "SD_CN9_100%RM"),
    (500_000, 600_000, 5.0, 0.30, "SD_CN5_30%Alu"),
    (600_000, 800_000, 1.0, 0.80, "Single-copy_80%Alu"),
    (800_000, 900_000, 1.0, 1.00, "Pure_repeat_100%RM"),
    (900_000, 1_000_000, 1.0, 0.00, "Single-copy_unique"),
]


def create_synthetic_kmer_bed(output_path):
    with open(output_path, 'w') as f:
        for start, end, cn, rm_frac, label in REGIONS:
            for pos in range(start, end):
                is_repeat = (pos % 100) < int(rm_frac * 100)

                if is_repeat:
                    count = int(COVERAGE * cn * ALU_MULTIPLICITY)
                else:
                    count = int(COVERAGE * cn)

                np.random.seed(pos)
                noise = 1.0 + np.random.uniform(-0.10, 0.10)
                count = max(1, int(count * noise))

                f.write(f"{CHROM}\t{pos}\t{pos+1}\t{count}\n")

    n_lines = GENOME_SIZE
    print(f"[SYNTH] Created k-mer BED: {output_path} ({n_lines:,} lines)")
    return n_lines


def create_synthetic_rm_mask(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    mask = np.zeros(GENOME_SIZE, dtype=np.uint8)

    for start, end, cn, rm_frac, label in REGIONS:
        for pos in range(start, end):
            if (pos % 100) < int(rm_frac * 100):
                mask[pos] = 1

    npy_path = os.path.join(output_dir, f"{CHROM}.npy")
    np.save(npy_path, mask)

    n_repeat = int(mask.sum())
    n_unique = len(mask) - n_repeat
    pct_repeat = 100 * n_repeat / len(mask)
    print(f"[SYNTH] Created RM mask: {npy_path}")
    print(f"[SYNTH]   {n_unique:,} unique + {n_repeat:,} repeat "
          f"({pct_repeat:.1f}% repeat)")

    manifest = {
        "chromosomes": {CHROM: {"length": GENOME_SIZE}},
        "total_masked_bp": n_repeat,
        "total_genome_bp": GENOME_SIZE,
    }
    manifest_path = os.path.join(output_dir, "manifest.json")
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"[SYNTH] Created manifest: {manifest_path}")

    return output_dir


def create_synthetic_rm_mask_multiclass(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    mask = np.zeros(GENOME_SIZE, dtype=np.uint8)

    region_classes = {
        "Single-copy_50%Alu": 4,
        "SD_CN2_50%Alu": 4,
        "SD_CN9_100%RM": 2,
        "SD_CN5_30%Alu": 4,
        "Single-copy_80%Alu": 4,
        "Pure_repeat_100%RM": 1,
        "Single-copy_unique": 0,
    }

    for start, end, cn, rm_frac, label in REGIONS:
        cls_id = region_classes.get(label, 1)
        if cls_id == 0:
            continue
        for pos in range(start, end):
            if (pos % 100) < int(rm_frac * 100):
                mask[pos] = cls_id

    npy_path = os.path.join(output_dir, f"{CHROM}.npy")
    np.save(npy_path, mask)

    n_repeat = int((mask > 0).sum())
    print(f"[SYNTH] Created multiclass RM mask: {npy_path}")
    print(f"[SYNTH]   Unique: {GENOME_SIZE - n_repeat:,}, Repeat: {n_repeat:,}")
    for cid in sorted(set(mask) - {0}):
        cnt = int((mask == cid).sum())
        print(f"[SYNTH]   Class {cid}: {cnt:,} positions")

    manifest = {
        "chromosomes": {CHROM: {"length": GENOME_SIZE}},
        "total_masked_bp": n_repeat,
        "total_genome_bp": GENOME_SIZE,
    }
    with open(os.path.join(output_dir, "manifest.json"), 'w') as f:
        json.dump(manifest, f, indent=2)
    return output_dir


def run_preprocessing_classw(input_bed, rm_mask_dir, output_bed):
    script = os.path.join(os.path.dirname(__file__), '..', 'scripts', 'pipeline',
                          'preprocess_kmer_windows.py')
    script = os.path.abspath(script)

    class_weights = ("Satellite=0.001,Simple_repeat=0.3,LINE=0.01,"
                     "SINE=0.001,LTR=0.05,DNA=0.1,Low_complexity=0.3,Other=0.05")

    cmd = [
        sys.executable, script,
        '--input', input_bed,
        '--output', output_bed,
        '--window-size', str(WINDOW_SIZE),
        '--sex', 'XX',
        '--bio-threshold-factor', '150.0',
        '--rm-mask-dir', rm_mask_dir,
        '--rm-class-weights', class_weights,
    ]

    print(f"\n{'='*60}")
    print(f"Running: PER-CLASS WEIGHTS mode")
    print(f"{'='*60}")

    import subprocess
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"STDERR:\n{result.stderr}")
        print(f"STDOUT:\n{result.stdout}")
        raise RuntimeError(f"Preprocessing failed (exit {result.returncode})")
    print(result.stdout)
    return output_bed


def run_preprocessing(input_bed, rm_mask_dir, output_bed, exclude_mode=False):
    script = os.path.join(os.path.dirname(__file__), '..', 'scripts', 'pipeline',
                          'preprocess_kmer_windows.py')
    script = os.path.abspath(script)

    cmd = [
        sys.executable, script,
        '--input', input_bed,
        '--output', output_bed,
        '--window-size', str(WINDOW_SIZE),
        '--sex', 'XX',
        '--bio-threshold-factor', '150.0',
        '--rm-mask-dir', rm_mask_dir,
    ]
    if exclude_mode:
        cmd.append('--rm-exclude')
    else:
        cmd.extend(['--repeat-kmer-weight', '0.01'])

    print(f"\n{'='*60}")
    print(f"Running: {'EXCLUSION' if exclude_mode else 'DOWN-WEIGHT (Sol.C)'} mode")
    print(f"{'='*60}")

    import subprocess
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"STDERR:\n{result.stderr}")
        print(f"STDOUT:\n{result.stdout}")
        raise RuntimeError(f"Preprocessing failed (exit {result.returncode})")
    print(result.stdout)
    return output_bed


def analyze_results(output_bed, mode_name):
    print(f"\n{'='*60}")
    print(f"RESULTS: {mode_name}")
    print(f"{'='*60}")

    windows = []
    with open(output_bed) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            chrom, start, end, cn = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            windows.append((chrom, start, end, cn))

    if not windows:
        print("ERROR: No windows in output!")
        return {}

    results = {}
    print(f"\n{'Region':<30s} {'Expected':>8s} {'Median':>8s} {'Mean':>8s} "
          f"{'Err%':>8s} {'N_win':>6s} {'Status'}")
    print("-" * 85)

    for reg_start, reg_end, true_cn, rm_frac, label in REGIONS:
        region_cns = [cn for (c, s, e, cn) in windows
                      if s >= reg_start and e <= reg_end]
        if not region_cns:
            print(f"{label:<30s} {'N/A':>8s}")
            continue
        arr = np.array(region_cns)
        med = float(np.median(arr))
        mean = float(np.mean(arr))
        err_pct = (med - true_cn) / true_cn * 100 if true_cn > 0 else 0
        n_win = len(arr)

        if true_cn == 1.0:
            tolerance = 0.30
        elif true_cn <= 3.0:
            tolerance = 0.30
        else:
            tolerance = 0.35

        if abs(err_pct) <= tolerance * 100:
            status = "PASS"
        elif abs(err_pct) <= 50:
            status = "MARGINAL"
        else:
            status = "FAIL"

        results[label] = {
            'expected': true_cn, 'median': med, 'mean': mean,
            'err_pct': err_pct, 'n_win': n_win, 'status': status,
        }
        print(f"{label:<30s} {true_cn:>8.1f} {med:>8.2f} {mean:>8.2f} "
              f"{err_pct:>+7.1f}% {n_win:>6d} {status}")

    all_cns = np.array([cn for (_, _, _, cn) in windows])
    n_total = len(all_cns)
    n_neutral = int(np.sum((all_cns >= 0.85) & (all_cns <= 1.25)))
    n_lowdup = int(np.sum((all_cns > 1.25) & (all_cns <= 3.0)))
    n_highdup = int(np.sum((all_cns > 3.0) & (all_cns <= 6.0)))
    n_amp = int(np.sum(all_cns > 6.0))
    n_del = int(np.sum(all_cns < 0.85))

    pct_neutral = 100 * n_neutral / n_total
    pct_dup = 100 * (n_lowdup + n_highdup + n_amp) / n_total
    pct_del = 100 * n_del / n_total

    print(f"\nGenom Dağılımı:")
    print(f"  Neutral (0.85-1.25):  {n_neutral:>5d} ({pct_neutral:>5.1f}%)")
    print(f"  LowDup  (1.25-3.0):  {n_lowdup:>5d} ({100*n_lowdup/n_total:>5.1f}%)")
    print(f"  HighDup  (3.0-6.0):  {n_highdup:>5d} ({100*n_highdup/n_total:>5.1f}%)")
    print(f"  Amp        (>6.0):   {n_amp:>5d} ({100*n_amp/n_total:>5.1f}%)")
    print(f"  Deleted    (<0.85):  {n_del:>5d} ({pct_del:>5.1f}%)")
    print(f"  ---")
    print(f"  Neutral:  {pct_neutral:.1f}%")
    print(f"  Non-neutral: {pct_dup:.1f}%")

    results['_distribution'] = {
        'pct_neutral': pct_neutral, 'pct_dup': pct_dup, 'pct_del': pct_del,
    }

    return results


def main():
    print("=" * 60)
    print("RM Exclusion Mode — Synthetic Validation Test")
    print("=" * 60)
    print(f"\nParameters:")
    print(f"  Genome size: {GENOME_SIZE/1e6:.1f} Mb")
    print(f"  Window size: {WINDOW_SIZE} bp")
    print(f"  Coverage (peak): {COVERAGE}")
    print(f"  Alu multiplicity: {ALU_MULTIPLICITY}")
    print(f"\nRegions:")
    for s, e, cn, rm, label in REGIONS:
        print(f"  [{s/1000:.0f}-{e/1000:.0f}kb] CN={cn}, RM={rm*100:.0f}% — {label}")

    test_dir = tempfile.mkdtemp(prefix="copyseg_rm_excl_test_")
    print(f"\nTest directory: {test_dir}")

    try:
        kmer_bed = os.path.join(test_dir, "synthetic_k32.bed")
        rm_mask_dir = os.path.join(test_dir, "rm_mask")
        create_synthetic_kmer_bed(kmer_bed)
        create_synthetic_rm_mask(rm_mask_dir)

        out_downweight = os.path.join(test_dir, "cn_downweight.bed")
        run_preprocessing(kmer_bed, rm_mask_dir, out_downweight, exclude_mode=False)
        results_dw = analyze_results(out_downweight, "DOWN-WEIGHT (Sol.C, weight=0.01)")

        out_exclude = os.path.join(test_dir, "cn_exclude.bed")
        run_preprocessing(kmer_bed, rm_mask_dir, out_exclude, exclude_mode=True)
        results_ex = analyze_results(out_exclude, "EXCLUSION (weight=0.0)")

        print(f"\n{'='*60}")
        print("KARŞILAŞTIRMA: Down-weight vs Exclusion")
        print(f"{'='*60}")
        print(f"\n{'Region':<30s} {'Expected':>8s} {'DW_CN':>8s} {'EX_CN':>8s} "
              f"{'DW_Err':>8s} {'EX_Err':>8s} {'DW':>6s} {'EX':>6s}")
        print("-" * 95)

        for label in [r[4] for r in REGIONS]:
            dw = results_dw.get(label, {})
            ex = results_ex.get(label, {})
            exp = dw.get('expected', 0)
            print(f"{label:<30s} {exp:>8.1f} "
                  f"{dw.get('median', 0):>8.2f} {ex.get('median', 0):>8.2f} "
                  f"{dw.get('err_pct', 0):>+7.1f}% {ex.get('err_pct', 0):>+7.1f}% "
                  f"{dw.get('status', 'N/A'):>6s} {ex.get('status', 'N/A'):>6s}")

        dist_dw = results_dw.get('_distribution', {})
        dist_ex = results_ex.get('_distribution', {})
        print(f"\nGenom Dağılımı:")
        print(f"  {'Metric':<20s} {'Down-weight':>12s} {'Exclusion':>12s} {'Target':>12s}")
        print(f"  {'Neutral%':<20s} {dist_dw.get('pct_neutral',0):>11.1f}% "
              f"{dist_ex.get('pct_neutral',0):>11.1f}% {'55-65%':>12s}")
        print(f"  {'Non-neutral%':<20s} {dist_dw.get('pct_dup',0):>11.1f}% "
              f"{dist_ex.get('pct_dup',0):>11.1f}% {'30-40%':>12s}")

        n_pass_dw = sum(1 for v in results_dw.values()
                        if isinstance(v, dict) and v.get('status') == 'PASS')
        n_pass_ex = sum(1 for v in results_ex.values()
                        if isinstance(v, dict) and v.get('status') == 'PASS')
        n_regions = len(REGIONS)

        print(f"\n{'='*60}")
        print(f"SONUÇ: Down-weight {n_pass_dw}/{n_regions} PASS | "
              f"Exclusion {n_pass_ex}/{n_regions} PASS")

        neutral_ok = 50 <= dist_ex.get('pct_neutral', 0) <= 70
        print(f"Neutral dağılımı ({'OK' if neutral_ok else 'FAIL'}): "
              f"{dist_ex.get('pct_neutral', 0):.1f}% (hedef: 55-65%)")
        print(f"{'='*60}")

        if n_pass_ex >= 5 and neutral_ok:
            print("\n✓ EXCLUSION MODE TEST PASSED — gerçek veri ile test edilebilir")
            return 0
        else:
            print(f"\n✗ EXCLUSION MODE TEST FAILED — "
                  f"{n_pass_ex}/{n_regions} PASS, neutral={dist_ex.get('pct_neutral',0):.1f}%")
            return 1

    finally:
        print(f"\nTest artifacts: {test_dir}")


def test_per_class_weights():
    print(f"\n{'='*60}")
    print("TEST: Per-Class RM Weights")
    print(f"{'='*60}")

    scripts_dir = os.path.join(os.path.dirname(__file__), '..', 'scripts', 'pipeline')
    scripts_dir = os.path.abspath(scripts_dir)
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)

    from preprocess_kmer_windows import RMMaskLoader, RM_CLASS_NAME_TO_ID

    test_dir = tempfile.mkdtemp(prefix="copyseg_classw_test_")
    mask_dir = os.path.join(test_dir, "rm_mask")
    os.makedirs(mask_dir)

    mask = np.zeros(1000, dtype=np.uint8)
    mask[200:400] = 4
    mask[400:600] = 2
    mask[600:800] = 3
    mask[800:1000] = 1

    np.save(os.path.join(mask_dir, "chr1.npy"), mask)
    manifest = {"chromosomes": {"chr1": {"length": 1000}}}
    with open(os.path.join(mask_dir, "manifest.json"), 'w') as f:
        json.dump(manifest, f)

    class_weights = {1: 0.001, 2: 0.30, 3: 0.01, 4: 0.001}
    loader = RMMaskLoader(mask_dir, class_weights=class_weights)
    positions = np.arange(1000)
    weights = loader.get_weights_bulk("chr1", positions)

    n_pass = 0
    n_fail = 0

    unique_w = weights[0:200]
    assert np.all(unique_w == 1.0), f"Unique weights should be 1.0, got {unique_w[:5]}"
    print(f"  Unique (0-200):        weight={unique_w[0]:.3f} — PASS")
    n_pass += 1

    sine_w = weights[200:400]
    assert np.allclose(sine_w, 0.001), f"SINE weights should be 0.001, got {sine_w[0]}"
    print(f"  SINE (200-400):        weight={sine_w[0]:.3f} — PASS")
    n_pass += 1

    sr_w = weights[400:600]
    assert np.allclose(sr_w, 0.30), f"Simple_repeat weights should be 0.3, got {sr_w[0]}"
    print(f"  Simple_repeat (400-600): weight={sr_w[0]:.3f} — PASS")
    n_pass += 1

    line_w = weights[600:800]
    assert np.allclose(line_w, 0.01), f"LINE weights should be 0.01, got {line_w[0]}"
    print(f"  LINE (600-800):        weight={line_w[0]:.3f} — PASS")
    n_pass += 1

    sat_w = weights[800:1000]
    assert np.allclose(sat_w, 0.001), f"Satellite weights should be 0.001, got {sat_w[0]}"
    print(f"  Satellite (800-1000):  weight={sat_w[0]:.3f} — PASS")
    n_pass += 1

    shutil.rmtree(test_dir)
    print(f"\n  Per-class weights test: {n_pass}/{n_pass + n_fail} PASS")
    return n_pass == 5


def test_boundary_erosion():
    print(f"\n{'='*60}")
    print("TEST: Boundary Erosion (k=32)")
    print(f"{'='*60}")

    scripts_dir = os.path.join(os.path.dirname(__file__), '..', 'scripts', 'pipeline')
    scripts_dir = os.path.abspath(scripts_dir)
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)

    from compute_rm_mask import erode_boundaries

    n_pass = 0

    mask = np.zeros(500, dtype=np.uint8)
    mask[100:400] = 4
    eroded = erode_boundaries(mask, kmer_size=32)

    assert eroded[68] == 0, "Position well before repeat should stay unique"
    assert eroded[100] == 0, "Entry boundary position should be eroded"
    assert eroded[130] == 0, "Position within entry margin should be eroded"
    assert eroded[131] == 4, "First surviving interior position should remain repeat"
    assert eroded[368] == 4, "Last surviving interior position should remain repeat"
    assert eroded[369] == 0, "Position within exit margin should be eroded"
    assert eroded[399] == 0, "Last position of repeat within margin should be eroded"
    assert eroded[400] == 0, "Position after repeat should stay unique"
    n_surviving = int((eroded == 4).sum())
    assert n_surviving == 238, f"Expected 238 surviving, got {n_surviving}"
    print(f"  Large repeat block (300bp): margins eroded, {n_surviving} surviving — PASS")
    n_pass += 1

    mask2 = np.zeros(200, dtype=np.uint8)
    mask2[80:120] = 3
    eroded2 = erode_boundaries(mask2, kmer_size=32)
    assert np.all(eroded2 == 0), "Short repeat block should be entirely eliminated"
    print(f"  Short repeat block (40bp < 62): eliminated — PASS")
    n_pass += 1

    mask3 = np.zeros(100, dtype=np.uint8)
    mask3[20:80] = 2
    eroded3 = erode_boundaries(mask3, kmer_size=1)
    assert np.array_equal(mask3, eroded3), "k=1 should not erode anything"
    print(f"  k=1 (no erosion): unchanged — PASS")
    n_pass += 1

    mask4 = np.zeros(200, dtype=np.uint8)
    mask4[0:100] = 1
    eroded4 = erode_boundaries(mask4, kmer_size=32)
    assert eroded4[0] == 0, "Start of chromosome repeat should be eroded"
    assert eroded4[30] == 0, "Position within start margin should be eroded"
    assert eroded4[31] == 1, "Interior position after start margin should remain"
    assert eroded4[68] == 1, "Interior position before exit margin should remain"
    assert eroded4[69] == 0, "Position within exit margin should be eroded"
    print(f"  Repeat at chromosome start: boundary eroded — PASS")
    n_pass += 1

    print(f"\n  Boundary erosion test: {n_pass}/4 PASS")
    return n_pass == 4


def test_backward_compat():
    print(f"\n{'='*60}")
    print("TEST: Backward Compatibility (no class_weights)")
    print(f"{'='*60}")

    scripts_dir = os.path.join(os.path.dirname(__file__), '..', 'scripts', 'pipeline')
    scripts_dir = os.path.abspath(scripts_dir)
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)

    from preprocess_kmer_windows import RMMaskLoader

    test_dir = tempfile.mkdtemp(prefix="copyseg_compat_test_")
    mask_dir = os.path.join(test_dir, "rm_mask")
    os.makedirs(mask_dir)

    mask = np.zeros(500, dtype=np.uint8)
    mask[100:200] = 4
    mask[200:300] = 2
    mask[300:400] = 3
    np.save(os.path.join(mask_dir, "chr1.npy"), mask)
    manifest = {"chromosomes": {"chr1": {"length": 500}}}
    with open(os.path.join(mask_dir, "manifest.json"), 'w') as f:
        json.dump(manifest, f)

    loader_legacy = RMMaskLoader(mask_dir, repeat_weight=0.01, class_weights=None)
    positions = np.arange(500)
    w_legacy = loader_legacy.get_weights_bulk("chr1", positions)

    assert np.all(w_legacy[0:100] == 1.0)
    assert np.all(w_legacy[400:500] == 1.0)
    assert np.allclose(w_legacy[100:200], 0.01)
    assert np.allclose(w_legacy[200:300], 0.01)
    assert np.allclose(w_legacy[300:400], 0.01)

    print(f"  Legacy flat weight (w=0.01): all repeat classes → 0.01 — PASS")

    loader_excl = RMMaskLoader(mask_dir, repeat_weight=0.0, class_weights=None)
    w_excl = loader_excl.get_weights_bulk("chr1", positions)
    assert np.all(w_excl[100:400] == 0.0)
    print(f"  Exclusion mode (w=0.0): all repeat classes → 0.0 — PASS")

    shutil.rmtree(test_dir)
    print(f"\n  Backward compat test: 2/2 PASS")
    return True


if __name__ == '__main__':
    print("\n" + "=" * 60)
    print("UNIT TESTS")
    print("=" * 60)
    unit_pass = True
    unit_pass &= test_per_class_weights()
    unit_pass &= test_boundary_erosion()
    unit_pass &= test_backward_compat()

    if not unit_pass:
        print("\nUnit tests FAILED — skipping integration test")
        sys.exit(1)

    print("\nUnit tests PASSED — running integration test...")

    print("\n" + "=" * 60)
    print("INTEGRATION: Per-Class Weights Mode")
    print("=" * 60)
    test_dir = tempfile.mkdtemp(prefix="copyseg_classw_integ_")
    kmer_bed = os.path.join(test_dir, "synthetic_k32.bed")
    rm_mask_dir = os.path.join(test_dir, "rm_mask")

    create_synthetic_kmer_bed(kmer_bed)
    create_synthetic_rm_mask_multiclass(rm_mask_dir)

    out_classw = os.path.join(test_dir, "cn_classw.bed")
    run_preprocessing_classw(kmer_bed, rm_mask_dir, out_classw)
    results_cw = analyze_results(out_classw, "PER-CLASS WEIGHTS")

    n_pass_cw = sum(1 for v in results_cw.values()
                    if isinstance(v, dict) and v.get('status') == 'PASS')
    n_regions = len(REGIONS)
    print(f"\n{'='*60}")
    print(f"Per-class weights: {n_pass_cw}/{n_regions} PASS")
    print(f"{'='*60}")
    print(f"Test artifacts: {test_dir}")

    print("\n" + "=" * 60)
    print("INTEGRATION: Exclusion vs Down-Weight")
    print("=" * 60)
    sys.exit(main())

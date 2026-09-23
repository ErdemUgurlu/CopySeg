#!/usr/bin/env python3

import os
import sys

import numpy as np
import pandas as pd

_THIS_DIR = os.path.dirname(__file__)
_SCRIPTS = os.path.abspath(os.path.join(_THIS_DIR, "..", "scripts", "pipeline"))
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

from refine_cn_unique import (
    refine_segment, refine_all, compute_unique_genome_median,
    build_windows_index, load_segments, load_windows,
    DEFAULT_MIN_UNIQUE_KMERS_PER_BIN, DEFAULT_MIN_TRUST_BINS,
    DEFAULT_MIN_TRUST_FRACTION, TANDEM_CLASSES,
)



WINDOW_SIZE = 500
COVERAGE    = 28
N_PER_BIN   = 469


def _make_windows(bins, chrom="chr1"):
    rows = []
    for b in bins:
        start = b["start"]
        cn    = b.get("cn", 1.0)
        mean_count = b.get("mean_count", cn * COVERAGE)
        rows.append({
            "chrom":             chrom,
            "start":             start,
            "end":               start + WINDOW_SIZE,
            "cn":                cn,
            "mean_count":        mean_count,
            "log_ratio":         np.log2(max(cn, 1e-6)),
            "num_kmers":         b.get("num_kmers", N_PER_BIN),
            "num_filtered":      b.get("num_filtered", 0),
            "unique_mean_count": b.get("unique_mean_count", np.nan),
            "n_unique_kmers":    b.get("n_unique_kmers", 0),
        })
    df = pd.DataFrame(rows)
    df["start"]              = df["start"].astype(int)
    df["end"]                = df["end"].astype(int)
    df["unique_mean_count"]  = pd.to_numeric(df["unique_mean_count"], errors="coerce")
    df["n_unique_kmers"]     = df["n_unique_kmers"].astype(int)
    return df


def _make_segments(segs, chrom="chr1"):
    rows = []
    for s in segs:
        rows.append({
            "chrom":          chrom,
            "start":          s["start"],
            "end":            s["end"],
            "state":          s.get("state", "Neutral"),
            "cn_median":      s["cn_median"],
            "cn_mean":        s.get("cn_mean", s["cn_median"]),
            "n_windows":      s.get("n_windows", (s["end"] - s["start"]) // WINDOW_SIZE),
            "avg_quality":    s.get("avg_quality", 1.0),
            "min_quality":    s.get("min_quality", 1.0),
            "cn_std":         s.get("cn_std", 0.0),
            "avg_repeats":    s.get("avg_repeats", 0.0),
            "avg_entropy":    0.0,
            "max_entropy":    0.0,
            "masked_fraction": s.get("masked_fraction", 0.0),
            "repeat_class":   s.get("repeat_class", "None"),
        })
    df = pd.DataFrame(rows)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    return df



def test_t1_healthy_sd():
    sd_bins = [{"start": 100_000 + i * WINDOW_SIZE,
                "cn": 7.0,
                "mean_count": 7 * COVERAGE,
                "unique_mean_count": 7.0 * COVERAGE,
                "n_unique_kmers": 234}
               for i in range(20)]

    bg_bins = [{"start": i * WINDOW_SIZE,
                "cn": 1.0,
                "mean_count": COVERAGE,
                "unique_mean_count": float(COVERAGE),
                "n_unique_kmers": 234}
               for i in range(1000)]

    windows = pd.concat([
        _make_windows(bg_bins, chrom="chrBG"),
        _make_windows(sd_bins, chrom="chr1"),
    ], ignore_index=True)

    segments = _make_segments([
        {"start": 100_000, "end": 100_000 + 20 * WINDOW_SIZE,
         "cn_median": 6.4,
         "state": "Amp",
         "repeat_class": "None",
         "masked_fraction": 0.50},
    ])

    refined = refine_all(segments, windows, sex="XX")
    assert len(refined) == 1
    row = refined.iloc[0]
    assert row["refine_method"] == "unique_median", row.to_dict()
    assert abs(row["cn_refined"] - 7.0) < 0.5, (
        f"Expected cn_refined ≈ 7.0, got {row['cn_refined']:.3f}")
    assert row["state_refined"] in {"HighDup", "Amp"}, row["state_refined"]
    print(f"  T1 PASS: cn_refined={row['cn_refined']:.3f}, "
          f"method={row['refine_method']}, state={row['state_refined']}")
    return True



def test_t2_vntr_simple_repeat():
    bg_bins = [{"start": i * WINDOW_SIZE,
                "cn": 1.0,
                "mean_count": COVERAGE,
                "unique_mean_count": float(COVERAGE),
                "n_unique_kmers": 234}
               for i in range(1000)]

    vntr_bins = [{"start": 200_000 + i * WINDOW_SIZE,
                  "cn": 4.5,
                  "mean_count": 4.5 * COVERAGE,
                  "unique_mean_count": np.nan,
                  "n_unique_kmers": 0}
                 for i in range(10)]

    windows = pd.concat([
        _make_windows(bg_bins, chrom="chrBG"),
        _make_windows(vntr_bins, chrom="chr1"),
    ], ignore_index=True)

    segments = _make_segments([
        {"start": 200_000, "end": 200_000 + 10 * WINDOW_SIZE,
         "cn_median": 4.49,
         "state": "HighDup",
         "repeat_class": "Simple_repeat",
         "masked_fraction": 0.95},
    ])

    refined = refine_all(segments, windows, sex="XX")
    row = refined.iloc[0]
    assert row["refine_method"] == "fallback_repeat_class", row.to_dict()
    assert abs(row["cn_refined"] - 4.49) < 1e-6, row["cn_refined"]
    print(f"  T2 PASS: cn_refined={row['cn_refined']:.3f}, "
          f"method={row['refine_method']}")
    return True



def test_t3_no_trust_bins():
    bg_bins = [{"start": i * WINDOW_SIZE,
                "cn": 1.0,
                "mean_count": COVERAGE,
                "unique_mean_count": float(COVERAGE),
                "n_unique_kmers": 234}
               for i in range(1000)]

    repeat_bins = [{"start": 300_000 + i * WINDOW_SIZE,
                    "cn": 2.0,
                    "mean_count": 2 * COVERAGE,
                    "unique_mean_count": COVERAGE * 1.05,
                    "n_unique_kmers": 5}
                   for i in range(10)]

    windows = pd.concat([
        _make_windows(bg_bins, chrom="chrBG"),
        _make_windows(repeat_bins, chrom="chr1"),
    ], ignore_index=True)

    segments = _make_segments([
        {"start": 300_000, "end": 300_000 + 10 * WINDOW_SIZE,
         "cn_median": 2.0,
         "state": "LowDup",
         "repeat_class": "SINE",
         "masked_fraction": 0.92},
    ])

    refined = refine_all(segments, windows, sex="XX")
    row = refined.iloc[0]
    assert row["refine_method"] == "fallback_too_few_unique", row.to_dict()
    assert abs(row["cn_refined"] - 2.0) < 1e-6
    assert row["n_trust_bins"] == 0, row["n_trust_bins"]
    print(f"  T3 PASS: n_trust={row['n_trust_bins']}, "
          f"method={row['refine_method']}")
    return True



def test_t4_tp53_unique_recovery():
    bg_bins = [{"start": i * WINDOW_SIZE,
                "cn": 1.0,
                "mean_count": COVERAGE,
                "unique_mean_count": float(COVERAGE),
                "n_unique_kmers": 234}
               for i in range(1000)]

    tp53_bins = [{"start": 400_000 + i * WINDOW_SIZE,
                  "cn": 1.71,
                  "mean_count": 1.71 * COVERAGE,
                  "unique_mean_count": float(COVERAGE),
                  "n_unique_kmers": 234}
                 for i in range(12)]

    windows = pd.concat([
        _make_windows(bg_bins, chrom="chrBG"),
        _make_windows(tp53_bins, chrom="chr1"),
    ], ignore_index=True)

    segments = _make_segments([
        {"start": 400_000, "end": 400_000 + 12 * WINDOW_SIZE,
         "cn_median": 1.71,
         "state": "LowDup",
         "repeat_class": "None",
         "masked_fraction": 0.50},
    ])

    refined = refine_all(segments, windows, sex="XX")
    row = refined.iloc[0]
    assert row["refine_method"] == "unique_median", row.to_dict()
    assert abs(row["cn_refined"] - 1.0) < 0.20, (
        f"Expected cn_refined ≈ 1.0, got {row['cn_refined']:.3f}")
    assert row["state_refined"] == "Neutral", row["state_refined"]
    print(f"  T4 PASS: cn_refined={row['cn_refined']:.3f} (down from 1.71), "
          f"state Neutral (was LowDup), method={row['refine_method']}")
    return True



def test_t5_state_relabel():
    bg_bins = [{"start": i * WINDOW_SIZE,
                "cn": 1.0,
                "unique_mean_count": float(COVERAGE),
                "n_unique_kmers": 234}
               for i in range(1000)]

    seg_bins = [{"start": 500_000 + i * WINDOW_SIZE,
                 "cn": 1.10,
                 "unique_mean_count": 2.0 * COVERAGE,
                 "n_unique_kmers": 234}
                for i in range(15)]

    windows = pd.concat([
        _make_windows(bg_bins, chrom="chrBG"),
        _make_windows(seg_bins, chrom="chr1"),
    ], ignore_index=True)

    segments = _make_segments([
        {"start": 500_000, "end": 500_000 + 15 * WINDOW_SIZE,
         "cn_median": 1.10,
         "state": "Neutral",
         "repeat_class": "None",
         "masked_fraction": 0.40},
    ])

    refined = refine_all(segments, windows, sex="XX")
    row = refined.iloc[0]
    assert row["refine_method"] == "unique_median", row.to_dict()
    assert abs(row["cn_refined"] - 2.0) < 0.20, row["cn_refined"]
    assert row["state_refined"] == "LowDup", row["state_refined"]
    print(f"  T5 PASS: cn 1.10 → {row['cn_refined']:.3f}, "
          f"state Neutral → {row['state_refined']}")
    return True



def test_t6_backward_compat():
    bg_bins = [{"start": i * WINDOW_SIZE,
                "cn": 1.0,
                "unique_mean_count": np.nan,
                "n_unique_kmers": 0}
               for i in range(50)]

    seg_bins = [{"start": 600_000 + i * WINDOW_SIZE,
                 "cn": 3.5,
                 "unique_mean_count": np.nan,
                 "n_unique_kmers": 0}
                for i in range(10)]

    windows = pd.concat([
        _make_windows(bg_bins, chrom="chrBG"),
        _make_windows(seg_bins, chrom="chr1"),
    ], ignore_index=True)

    segments = _make_segments([
        {"start": 600_000, "end": 600_000 + 10 * WINDOW_SIZE,
         "cn_median": 3.5,
         "state": "HighDup",
         "repeat_class": "None",
         "masked_fraction": 0.20},
    ])

    median = compute_unique_genome_median(windows, sex="XX")
    assert median is None, f"Expected None (no valid bins), got {median}"

    refined = refine_all(segments, windows, sex="XX")
    row = refined.iloc[0]
    assert row["refine_method"] == "fallback_too_few_unique", row.to_dict()
    assert abs(row["cn_refined"] - 3.5) < 1e-6
    assert np.isfinite(row["cn_refined"])
    print(f"  T6 PASS: no-mask → fallback, cn_refined={row['cn_refined']:.3f}, "
          f"method={row['refine_method']}")
    return True



def test_t7_xy_chrx_exclusion():
    chrx_bins = [{"start": i * WINDOW_SIZE,
                  "cn": 0.5,
                  "unique_mean_count": COVERAGE * 0.5,
                  "n_unique_kmers": 234}
                 for i in range(500)]
    chr1_bins = [{"start": i * WINDOW_SIZE,
                  "cn": 1.0,
                  "unique_mean_count": float(COVERAGE),
                  "n_unique_kmers": 234}
                 for i in range(500)]

    windows = pd.concat([
        _make_windows(chrx_bins, chrom="chrX"),
        _make_windows(chr1_bins, chrom="chr1"),
    ], ignore_index=True)

    median_xy = compute_unique_genome_median(windows, sex="XY")
    assert abs(median_xy - COVERAGE) < 1.0, (
        f"XY median should ≈ {COVERAGE}, got {median_xy:.3f}")

    median_xx = compute_unique_genome_median(windows, sex="XX")
    print(f"  T7 PASS: median_xy={median_xy:.3f}, median_xx={median_xx:.3f}")
    return True



def test_t8_multi_segment_distribution():
    bg_bins = [{"start": i * WINDOW_SIZE,
                "cn": 1.0,
                "unique_mean_count": float(COVERAGE),
                "n_unique_kmers": 234}
               for i in range(1000)]

    sd_bins = [{"start": 1_000_000 + i * WINDOW_SIZE,
                "cn": 6.5,
                "unique_mean_count": 7.0 * COVERAGE,
                "n_unique_kmers": 234} for i in range(12)]
    vntr_bins = [{"start": 2_000_000 + i * WINDOW_SIZE,
                  "cn": 4.0,
                  "unique_mean_count": np.nan,
                  "n_unique_kmers": 0} for i in range(10)]
    repeat_bins = [{"start": 3_000_000 + i * WINDOW_SIZE,
                    "cn": 2.0,
                    "unique_mean_count": COVERAGE * 1.05,
                    "n_unique_kmers": 5} for i in range(10)]
    tp53_bins = [{"start": 4_000_000 + i * WINDOW_SIZE,
                  "cn": 1.5,
                  "unique_mean_count": float(COVERAGE),
                  "n_unique_kmers": 234} for i in range(12)]

    windows = pd.concat([
        _make_windows(bg_bins, chrom="chrBG"),
        _make_windows(sd_bins, chrom="chr1"),
        _make_windows(vntr_bins, chrom="chr1"),
        _make_windows(repeat_bins, chrom="chr1"),
        _make_windows(tp53_bins, chrom="chr1"),
    ], ignore_index=True)

    segments = _make_segments([
        {"start": 1_000_000, "end": 1_000_000 + 12 * WINDOW_SIZE,
         "cn_median": 6.5, "state": "Amp", "repeat_class": "None",
         "masked_fraction": 0.5},
        {"start": 2_000_000, "end": 2_000_000 + 10 * WINDOW_SIZE,
         "cn_median": 4.0, "state": "HighDup", "repeat_class": "Simple_repeat",
         "masked_fraction": 0.95},
        {"start": 3_000_000, "end": 3_000_000 + 10 * WINDOW_SIZE,
         "cn_median": 2.0, "state": "LowDup", "repeat_class": "SINE",
         "masked_fraction": 0.92},
        {"start": 4_000_000, "end": 4_000_000 + 12 * WINDOW_SIZE,
         "cn_median": 1.5, "state": "LowDup", "repeat_class": "None",
         "masked_fraction": 0.4},
    ])

    refined = refine_all(segments, windows, sex="XX")
    methods = list(refined["refine_method"].values)
    expected = ["unique_median", "fallback_repeat_class",
                "fallback_too_few_unique", "unique_median"]
    assert methods == expected, f"Methods: {methods}, expected: {expected}"
    print(f"  T8 PASS: methods distributed correctly: {methods}")
    return True



ALL_TESTS = [
    ("T1 — Healthy SD CN=7 → refined ≈ 7.0",                test_t1_healthy_sd),
    ("T2 — VNTR (Simple_repeat) → fallback_repeat_class",   test_t2_vntr_simple_repeat),
    ("T3 — All-repeat, no trust bins → fallback",           test_t3_no_trust_bins),
    ("T4 — TP53-like, unique recovery → CN ≈ 1.0",          test_t4_tp53_unique_recovery),
    ("T5 — State relabel Neutral → LowDup",                 test_t5_state_relabel),
    ("T6 — Backward-compat: no-mask file → no crash",       test_t6_backward_compat),
    ("T7 — XY sex-chrom exclusion from median",             test_t7_xy_chrx_exclusion),
    ("T8 — Multi-segment method distribution",              test_t8_multi_segment_distribution),
]


def main():
    print("=" * 70)
    print("Pass-2 Unique-Only Refinement — Unit Tests")
    print("=" * 70)
    n_pass = 0
    n_fail = 0
    for name, fn in ALL_TESTS:
        print(f"\n[{name}]")
        try:
            fn()
            n_pass += 1
        except Exception as e:
            print(f"  FAIL: {e}")
            import traceback
            traceback.print_exc()
            n_fail += 1
    print()
    print("=" * 70)
    print(f"RESULTS: {n_pass}/{len(ALL_TESTS)} PASSED, {n_fail} FAILED")
    print("=" * 70)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

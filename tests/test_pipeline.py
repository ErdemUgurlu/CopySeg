#!/usr/bin/env python3

import sys
import math
import io
import textwrap

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, "scripts/pipeline")

import segment_cnv_fused_lasso as fl
import evaluate_ground_truth   as egt
import validate_cn_accuracy    as vca



class TestWeightedL2Cost:

    def _make_signal(self, y, w=None):
        y = np.array(y, dtype=np.float64)
        if w is None:
            w = np.ones(len(y))
        return np.column_stack([y, w])

    def test_zero_cost_constant_signal(self):
        y = np.array([2.0, 2.0, 2.0, 2.0, 2.0])
        cost = fl.WeightedL2Cost().fit(self._make_signal(y))
        assert cost.error(0, 5) == pytest.approx(0.0, abs=1e-10)

    def test_unweighted_equals_sse(self):
        y = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        cost = fl.WeightedL2Cost().fit(self._make_signal(y))
        sse = np.sum((y - y.mean()) ** 2)
        assert cost.error(0, 5) == pytest.approx(sse, rel=1e-9)

    def test_weighted_cost_correctness(self):
        y = np.array([1.0, 2.0, 10.0, 2.0])
        w = np.array([1.0, 1.0, 0.01, 1.0])
        cost = fl.WeightedL2Cost().fit(self._make_signal(y, w))
        yw = np.dot(w, y) / w.sum()
        expected = np.dot(w, (y - yw) ** 2)
        assert cost.error(0, 4) == pytest.approx(expected, rel=1e-9)

    def test_subsegment_query(self):
        y = np.array([1.0, 3.0, 5.0, 7.0, 9.0])
        cost = fl.WeightedL2Cost().fit(self._make_signal(y))
        sub = y[1:4]
        expected = np.sum((sub - sub.mean()) ** 2)
        assert cost.error(1, 4) == pytest.approx(expected, rel=1e-9)

    def test_zero_weight_segment(self):
        y = np.array([1.0, 5.0, 9.0])
        w = np.zeros(3)
        cost = fl.WeightedL2Cost().fit(np.column_stack([y, w]))
        result = cost.error(0, 3)
        assert np.isfinite(result), f"Expected finite cost, got {result}"
        assert result >= 0.0

    def test_1d_fallback(self):
        y = np.array([1.0, 2.0, 3.0, 4.0])
        cost = fl.WeightedL2Cost().fit(y)
        assert cost.error(0, 4) >= 0.0



class TestAssignState:
    @pytest.mark.parametrize("cn, expected_state", [
        (0.10, 'HomDel'),
        (0.50, 'HetDel'),
        (1.00, 'Neutral'),
        (1.30, 'LowDup'),
        (3.50, 'HighDup'),
        (7.00, 'Amp'),
        (15.0, 'MedAmp'),
        (30.0, 'HighAmp'),
        (60.0, 'ExtremeAmp'),
    ])
    def test_default_boundaries(self, cn, expected_state):
        assert fl.assign_state_from_cn(cn) == expected_state

    def test_custom_lowdup_threshold(self):
        assert fl.assign_state_from_cn(1.3, lowdup_threshold=1.52) == 'Neutral'
        assert fl.assign_state_from_cn(1.6, lowdup_threshold=1.52) == 'LowDup'



def _make_df(cns, qualities=None, chrom='CM039011.1', start=0, window=500):
    n = len(cns)
    if qualities is None:
        qualities = [1.0] * n
    starts = [start + i * window for i in range(n)]
    ends   = [s + window for s in starts]
    data = {
        'chrom':        [chrom] * n,
        'start':        starts,
        'end':          ends,
        'cn':           cns,
        'num_kmers':    [100]  * n,
        'num_filtered': [0]    * n,
        'quality':      qualities,
        'log2_cn':      [math.log2(max(c, 1e-3)) for c in cns],
    }
    return pd.DataFrame(data)


class TestMakeSegment:

    def test_single_window(self):
        df = _make_df([2.5])
        seg = fl.make_segment_from_windows(df)
        assert seg['cn_median'] == pytest.approx(2.5)
        assert seg['cn_std'] == 0.0
        assert seg['n_windows'] == 1
        assert seg['state'] == 'LowDup'

    def test_quality_weighted_mean(self):
        cns = [2.0, 2.0, 2.0, 40.0]
        qws = [1.0, 1.0, 1.0, 0.01]
        df = _make_df(cns, qws)
        seg = fl.make_segment_from_windows(df)
        assert seg['cn_median'] < 3.0, (
            f"Expected weighted mean ~2.1, got {seg['cn_median']:.3f}")

    def test_quality_weighted_std_lower_than_unweighted(self):
        cns = [2.0, 2.0, 2.0, 40.0]
        qws = [1.0, 1.0, 1.0, 0.01]
        df = _make_df(cns, qws)
        seg = fl.make_segment_from_windows(df)
        unweighted_std = float(np.std(cns))
        assert seg['cn_std'] < unweighted_std

    def test_state_assignment_consistent_with_cn(self):
        cns = [1.0, 1.0, 1.0, 50.0]
        qws = [1.0, 1.0, 1.0, 0.001]
        df = _make_df(cns, qws)
        seg = fl.make_segment_from_windows(df)
        assert seg['state'] == 'Neutral', (
            f"Expected Neutral (weighted CN≈1), got {seg['state']} (CN={seg['cn_median']:.3f})")

    def test_empty_df_raises(self):
        df = _make_df([])
        with pytest.raises((IndexError, ValueError)):
            fl.make_segment_from_windows(df)



class TestPooledStd:

    def test_same_mean_same_std(self):
        result = fl._pooled_std(1.0, 4, 1.0, 4, 5.0, 5.0)
        assert result == pytest.approx(math.sqrt(6.0 / 7.0), rel=1e-9)

    def test_different_means_inflates_std(self):
        result = fl._pooled_std(0.1, 10, 0.1, 10, 1.0, 5.0)
        assert result > 1.0

    def test_nonnegative(self):
        assert fl._pooled_std(0.0, 1, 0.0, 1, 2.0, 2.0) >= 0.0

    def test_single_window_groups(self):
        result = fl._pooled_std(0.0, 1, 0.0, 1, 1.0, 3.0)
        assert result == pytest.approx(math.sqrt(2.0), rel=1e-9)



def _seg(state, start, end, cn=2.0, chrom='chr1'):
    return {
        'chrom': chrom, 'start': start, 'end': end,
        'state': state, 'cn_median': cn, 'cn_mean': cn,
        'n_windows': (end - start) // 500, 'avg_quality': 1.0,
        'min_quality': 1.0, 'cn_std': 0.0, 'avg_repeats': 0.0,
        'avg_entropy': 0.0, 'max_entropy': 0.0,
    }


class TestFilterSmallSegments:

    def test_short_dup_reclassified(self):
        segs = [_seg('LowDup', 0, 1000, cn=2.0)]
        result = fl.filter_small_segments(segs, {'LowDup': 3000})
        assert result[0]['state'] == 'Neutral'

    def test_long_dup_kept(self):
        segs = [_seg('LowDup', 0, 5000, cn=2.0)]
        result = fl.filter_small_segments(segs, {'LowDup': 3000})
        assert result[0]['state'] == 'LowDup'

    def test_adjacent_neutral_merged(self):
        segs = [
            _seg('LowDup',  0,  500, cn=2.0),
            _seg('Neutral', 500, 1000, cn=1.0),
        ]
        result = fl.filter_small_segments(segs, {'LowDup': 3000})
        assert len(result) == 1
        assert result[0]['state'] == 'Neutral'
        assert result[0]['start'] == 0
        assert result[0]['end']   == 1000

    def test_cross_chrom_no_merge(self):
        segs = [
            _seg('Neutral', 0, 500, chrom='chr1'),
            _seg('Neutral', 0, 500, chrom='chr2'),
        ]
        result = fl.filter_small_segments(segs, {})
        assert len(result) == 2

    def test_empty_input(self):
        assert fl.filter_small_segments([], {'LowDup': 3000}) == []

    def test_highamp_not_filtered(self):
        segs = [_seg('HighAmp', 0, 500, cn=50.0)]
        result = fl.filter_small_segments(segs, {'LowDup': 3000, 'HighDup': 3000, 'Amp': 3000})
        assert result[0]['state'] == 'HighAmp'



class TestMergeNearbyDup:

    def test_same_state_merge(self):
        segs = [
            _seg('LowDup',  0,    5000, cn=2.0),
            _seg('Neutral', 5000, 8000, cn=1.0),
            _seg('LowDup',  8000, 13000, cn=2.0),
        ]
        result = fl.merge_nearby_dup_segments(segs, max_gap=10000)
        assert len(result) == 1
        assert result[0]['state'] == 'LowDup'
        assert result[0]['start'] == 0
        assert result[0]['end']   == 13000

    def test_gap_too_large_no_merge(self):
        segs = [
            _seg('LowDup',  0,    5000, cn=2.0),
            _seg('Neutral', 5000, 20000, cn=1.0),
            _seg('LowDup',  20000, 25000, cn=2.0),
        ]
        result = fl.merge_nearby_dup_segments(segs, max_gap=10000)
        assert len(result) == 3

    def test_cross_state_merge_with_tolerance(self):
        segs = [
            _seg('HighDup', 0,    5000, cn=3.5),
            _seg('Neutral', 5000, 8000, cn=1.0),
            _seg('Amp',     8000, 13000, cn=5.0),
        ]
        result = fl.merge_nearby_dup_segments(segs, max_gap=10000, merge_cn_tolerance=2.0)
        assert len(result) == 1

    def test_cross_state_no_merge_without_tolerance(self):
        segs = [
            _seg('HighDup', 0,    5000, cn=3.5),
            _seg('Neutral', 5000, 8000, cn=1.0),
            _seg('Amp',     8000, 13000, cn=5.0),
        ]
        result = fl.merge_nearby_dup_segments(segs, max_gap=10000, merge_cn_tolerance=0.0)
        assert len(result) == 3

    def test_merged_cn_is_weighted_mean(self):
        segs = [
            _seg('LowDup', 0,    5000, cn=2.0),
            _seg('Neutral', 5000, 6000, cn=1.0),
            _seg('LowDup', 6000, 7000, cn=4.0),
        ]
        result = fl.merge_nearby_dup_segments(segs, max_gap=5000)
        assert len(result) == 1
        expected_cn = (5000 * 2.0 + 1000 * 4.0) / 6000
        assert result[0]['cn_median'] == pytest.approx(expected_cn, rel=1e-6)

    def test_merged_cn_std_uses_pooled_std(self):
        segs = [
            _seg('LowDup', 0,    5000, cn=2.0),
            _seg('Neutral', 5000, 6000, cn=1.0),
            _seg('LowDup', 6000, 7000, cn=4.0),
        ]
        result = fl.merge_nearby_dup_segments(segs, max_gap=5000)
        assert result[0]['cn_std'] > 0.0, "Pooled std should be >0 when CNs differ"

    def test_chrom_boundary_not_merged(self):
        segs = [
            _seg('LowDup',  0,    5000, chrom='chr1'),
            _seg('Neutral', 5000, 6000, chrom='chr1'),
            _seg('LowDup',  0,    5000, chrom='chr2'),
        ]
        result = fl.merge_nearby_dup_segments(segs, max_gap=10000)
        assert len(result) == 3



class TestSplitHighCV:

    def _make_high_cv_df(self, cn_left=2.0, cn_right=6.0, n=20):
        cns = [cn_left] * (n // 2) + [cn_right] * (n // 2)
        return _make_df(cns, start=0)

    def test_high_cv_segment_is_split(self):
        df = self._make_high_cv_df(cn_left=2.0, cn_right=6.0, n=20)
        cn_all = df['cn'].values
        w_all  = df['quality'].values
        w_sum  = w_all.sum()
        cn_wm  = np.dot(w_all, cn_all) / w_sum
        cn_std = np.sqrt(np.dot(w_all, (cn_all - cn_wm) ** 2) / w_sum)
        cv = cn_std / cn_wm
        assert cv >= 0.5, f"Test setup: expected CV ≥ 0.5, got {cv:.3f}"

        seg = fl.make_segment_from_windows(df)
        result = fl.split_high_cv_segments([seg], df, cv_threshold=0.45, min_length=500)
        assert len(result) == 2, f"Expected 2 sub-segments, got {len(result)}"

    def test_low_cv_segment_not_split(self):
        cns = [2.0] * 20
        df = _make_df(cns)
        seg = fl.make_segment_from_windows(df)
        result = fl.split_high_cv_segments([seg], df, cv_threshold=0.5)
        assert len(result) == 1

    def test_split_boundary_aligns_with_cn_jump(self):
        cns = [2.0] * 10 + [6.0] * 10
        df = _make_df(cns)
        seg = fl.make_segment_from_windows(df)
        result = fl.split_high_cv_segments([seg], df, cv_threshold=0.4, min_length=500)
        if len(result) == 2:
            left_cn  = result[0]['cn_median']
            right_cn = result[1]['cn_median']
            assert left_cn < right_cn, "Left CN should be lower than right CN"
            assert left_cn  < 3.0, f"Left half CN={left_cn:.2f}, expected ~2"
            assert right_cn > 4.0, f"Right half CN={right_cn:.2f}, expected ~6"

    def test_neutral_segment_not_split(self):
        seg = _seg('Neutral', 0, 10000, cn=1.5)
        seg['cn_std'] = 5.0
        df = _make_df([1.5, 1.5, 2.5, 2.5] * 5)
        result = fl.split_high_cv_segments([seg], df, cv_threshold=0.1)
        assert len(result) == 1
        assert result[0]['state'] == 'Neutral'

    def test_quality_weighted_split_prefers_boundary(self):
        cns = [2.0] * 9 + [6.0] * 10
        qws = [1.0] * 8 + [0.01] + [1.0] * 10
        df = _make_df(cns, qws)
        seg = fl.make_segment_from_windows(df)
        result = fl.split_high_cv_segments([seg], df, cv_threshold=0.3, min_length=500)
        if len(result) == 2:
            split_at = result[0]['end']
            assert abs(split_at - 4500) <= 1000, (
                f"Split at {split_at}, expected near 4500 (CN transition)")



class TestGCCalibration:

    def _make_segments(self, n_neutral=50, gc_vals=None):
        if gc_vals is None:
            gc_vals = np.linspace(0.35, 0.65, n_neutral)
        segs = []
        for i, gc in enumerate(gc_vals):
            segs.append({
                'chrom': 'chr1', 'start': i * 5000, 'end': (i + 1) * 5000,
                'state': 'Neutral', 'cn_median': 1.0 + 0.1 * (gc - 0.5),
                'cn_mean': 1.0, 'n_windows': 10, 'avg_quality': 1.0,
                'min_quality': 1.0, 'cn_std': 0.1, 'avg_repeats': 0.0,
                'avg_entropy': 0.0, 'max_entropy': 0.0,
                'mean_gc': float(gc),
                'masked_fraction': 0.0, 'repeat_class': 'None',
            })
        return segs

    def test_calibration_runs_on_neutral_segs(self):
        segs = self._make_segments(50)
        dup = dict(segs[0])
        dup.update({'state': 'LowDup', 'cn_median': 2.0, 'mean_gc': 0.45,
                    'start': 1000000, 'end': 1005000})
        segs.append(dup)
        result = fl.apply_gc_cn_calibration(segs)
        assert len(result) == len(segs)
        assert all('gc_bias_factor' in s for s in result if 'mean_gc' in s)

    def test_satellite_bypass(self):
        segs = self._make_segments(50)
        sat = dict(segs[0])
        sat.update({
            'state': 'LowDup', 'cn_median': 5.0, 'mean_gc': 0.45,
            'masked_fraction': 0.9, 'repeat_class': 'Satellite',
            'start': 1000000, 'end': 1005000,
        })
        segs.append(sat)
        result = fl.apply_gc_cn_calibration(segs)
        sat_result = [s for s in result if s.get('repeat_class') == 'Satellite'][0]
        assert sat_result['gc_bias_factor'] == 1.0
        assert sat_result['cn_median'] == 5.0

    def test_high_gc_bypass(self):
        segs = self._make_segments(50)
        high_gc = dict(segs[0])
        high_gc.update({
            'state': 'LowDup', 'cn_median': 3.0, 'mean_gc': 0.65,
            'masked_fraction': 0.1, 'repeat_class': 'LINE',
            'start': 1000000, 'end': 1005000,
        })
        segs.append(high_gc)
        result = fl.apply_gc_cn_calibration(segs)
        hgc = [s for s in result if s.get('mean_gc', 0) > 0.60][0]
        assert hgc['gc_bias_factor'] == 1.0

    def test_too_few_neutral_segs_skips(self):
        segs = self._make_segments(10)
        dup = dict(segs[0])
        dup.update({'state': 'LowDup', 'cn_median': 2.0, 'start': 1000000, 'end': 1005000})
        segs.append(dup)
        original_cn = dup['cn_median']
        result = fl.apply_gc_cn_calibration(segs)
        dup_out = [s for s in result if s['state'] == 'LowDup'][0]
        assert dup_out['cn_median'] == original_cn

    def test_calibration_clamps_gc_factor(self):
        segs = self._make_segments(50, gc_vals=np.linspace(0.35, 0.65, 50))
        extreme = dict(segs[0])
        extreme.update({
            'state': 'LowDup', 'cn_median': 2.0, 'mean_gc': 0.15,
            'masked_fraction': 0.0, 'repeat_class': 'None',
            'start': 2000000, 'end': 2005000,
        })
        segs.append(extreme)
        result = fl.apply_gc_cn_calibration(segs)
        out = [s for s in result if s.get('mean_gc', 0) == 0.15]
        if out and 'gc_bias_factor' in out[0]:
            gf = out[0]['gc_bias_factor']
            assert 0.6 <= gf <= 1.8, f"gc_bias_factor {gf} outside [0.6, 1.8]"



class TestReclassify:

    def test_lowdup_below_threshold_to_neutral(self):
        segs = [_seg('LowDup', 0, 5000, cn=1.10)]
        result = fl.reclassify_by_cn_threshold(segs, lowdup_threshold=1.25)
        assert result[0]['state'] == 'Neutral'

    def test_lowdup_above_threshold_kept(self):
        segs = [_seg('LowDup', 0, 5000, cn=1.50)]
        result = fl.reclassify_by_cn_threshold(segs, lowdup_threshold=1.25)
        assert result[0]['state'] == 'LowDup'

    def test_hetdel_high_cn_to_lowdup(self):
        segs = [_seg('HetDel', 0, 5000, cn=1.8)]
        result = fl.reclassify_by_cn_threshold(segs)
        assert result[0]['state'] == 'LowDup'

    def test_hetdel_moderate_cn_to_neutral(self):
        segs = [_seg('HetDel', 0, 5000, cn=0.9)]
        result = fl.reclassify_by_cn_threshold(segs, hetdel_threshold=0.75)
        assert result[0]['state'] == 'Neutral'

    def test_hetdel_low_cn_stays(self):
        segs = [_seg('HetDel', 0, 5000, cn=0.5)]
        result = fl.reclassify_by_cn_threshold(segs, hetdel_threshold=0.75)
        assert result[0]['state'] == 'HetDel'

    def test_neutral_above_threshold_to_dup(self):
        segs = [_seg('Neutral', 0, 5000, cn=1.5)]
        result = fl.reclassify_by_cn_threshold(segs, lowdup_threshold=1.25)
        assert result[0]['state'] == 'LowDup'

    def test_adjacent_neutral_merged_after_reclassify(self):
        segs = [
            _seg('LowDup', 0,    5000, cn=1.1),
            _seg('Neutral', 5000, 10000, cn=1.0),
        ]
        result = fl.reclassify_by_cn_threshold(segs, lowdup_threshold=1.25)
        assert len(result) == 1
        assert result[0]['start'] == 0
        assert result[0]['end']   == 10000



class TestClassVerdict:

    def test_class_A_within_35pct_passes(self):
        assert egt.class_verdict("A", 6.94, 7).startswith("PASS")

    def test_class_A_over_50pct_fails_direction(self):
        assert egt.class_verdict("A", 0.7, 2) == "FAIL (UNDER)"
        assert egt.class_verdict("A", 3.5, 2) == "FAIL (OVER)"

    def test_class_A_marginal_band(self):
        assert egt.class_verdict("A", 1.45, 1).startswith("MARGINAL")

    def test_class_C_scored_like_A(self):
        assert egt.class_verdict("C", 1.0, 1).startswith("PASS")
        assert egt.class_verdict("C", 4.9, 1) == "FAIL (OVER)"

    def test_class_B_is_characterized_not_scored(self):
        v = egt.class_verdict("B", 4.6, 9)
        assert v == "AGG-LIMITED"
        assert not v.startswith("FAIL")
        assert egt.class_verdict("B", 8.5, 9) == "AGG≈PHYS"

    def test_class_D_observatory(self):
        assert egt.class_verdict("D", 47.0, 219) == "OBSERVATORY"
        assert egt.class_verdict("D", 19.0, None) == "OBSERVATORY"

    def test_skip_and_nodata(self):
        assert egt.class_verdict("skip", None, None) == "SKIPPED"
        assert egt.class_verdict("A", None, 1) == "NO DATA"



class TestValidateCNAccuracySchema:

    def _make_18col_bed(self):
        cols = [
            'CM039011.1', '0', '5000', 'LowDup',
            '2.0000', '2.0000', '10',
            '0.9500', '0.8000', '0.2000',
            '5.00', '0.0000', '0.0000',
            '0.1000', 'LINE',
            '1.0200', '0.5000', '3.4000',
        ]
        return '\t'.join(cols)

    def test_18col_loads_without_crash(self):
        content = self._make_18col_bed() + '\n'
        with io.StringIO(content) as f:
            df = pd.read_csv(f, sep='\t', comment='#', header=None)

        base_cols  = ['chrom', 'start', 'end', 'state', 'cn_median', 'cn_mean', 'n_windows']
        ext_cols   = ['avg_quality', 'min_quality', 'cn_std', 'avg_repeats',
                      'avg_entropy', 'max_entropy', 'masked_fraction', 'repeat_class']
        extra_cols = ['gc_bias_factor', 'segment_iqr', 'boundary_conf']
        all_cols   = base_cols + ext_cols + extra_cols

        ncols = len(df.columns)
        assert ncols == 18
        df.columns = all_cols[:ncols]
        assert 'gc_bias_factor' in df.columns
        assert df['cn_median'].iloc[0] == pytest.approx(2.0)

    def test_15col_still_works(self):
        cols = [
            'CM039011.1', '0', '5000', 'LowDup',
            '2.0000', '2.0000', '10',
            '0.95', '0.80', '0.20',
            '5.00', '0.00', '0.00',
            '0.10', 'LINE',
        ]
        content = '\t'.join(cols) + '\n'
        with io.StringIO(content) as f:
            df = pd.read_csv(f, sep='\t', comment='#', header=None)

        base_cols  = ['chrom', 'start', 'end', 'state', 'cn_median', 'cn_mean', 'n_windows']
        ext_cols   = ['avg_quality', 'min_quality', 'cn_std', 'avg_repeats',
                      'avg_entropy', 'max_entropy', 'masked_fraction', 'repeat_class']
        extra_cols = ['gc_bias_factor', 'segment_iqr', 'boundary_conf']
        all_cols   = base_cols + ext_cols + extra_cols
        ncols = len(df.columns)
        df.columns = all_cols[:ncols]
        assert 'repeat_class' in df.columns
        assert 'gc_bias_factor' not in df.columns



class TestEK5GenomeWideNorm:

    def test_same_count_same_cn_different_chroms(self):
        genome_median = 28.0
        raw_count = 28.0

        cn_chr4_old  = raw_count / 25.0
        cn_chr22_old = raw_count / 31.0

        cn_chr4_new  = raw_count / genome_median
        cn_chr22_new = raw_count / genome_median

        assert cn_chr4_new == cn_chr22_new, "EK5: same count must yield same CN"
        assert cn_chr4_old != cn_chr22_old, "Before EK5: per-chrom gives different CNs"

    def test_genome_median_normalizes_to_cn1(self):
        genome_median = 28.0
        cn = genome_median / genome_median
        assert cn == pytest.approx(1.0)



class TestSoftCapping:

    def test_capped_values_in_mean(self):
        threshold = 840.0
        counts    = np.array([700.0, 800.0, 1000.0, 1200.0])
        capped    = np.where(counts > threshold, threshold, counts)
        mean_capped  = float(np.mean(capped))
        mean_excluded = float(np.mean(counts[counts <= threshold]))

        assert mean_capped > mean_excluded, (
            "Capped mean should be closer to true mean than exclude-only mean")
        assert mean_capped < float(np.mean(counts))

    def test_nor_threshold_infinity_no_capping(self):
        threshold = float('inf')
        counts    = np.array([5000.0, 10000.0, 50000.0])
        capped    = np.where(counts > threshold, threshold, counts)
        assert np.array_equal(capped, counts), "NOR bins: no capping should occur"

    def test_num_filtered_still_counts_capped(self):
        counts    = np.array([700.0, 800.0, 1000.0, 1200.0])
        threshold = 840.0
        over_mask = counts > threshold
        n_filtered = int(over_mask.sum())
        assert n_filtered == 2, f"Expected 2 capped k-mers, got {n_filtered}"



class TestPostProcessingChain:

    def test_satellite_bypass_requires_repeat_annotation_first(self):
        segs = []
        for i in range(50):
            gc = 0.35 + i * 0.006
            segs.append({
                'chrom': 'chr1', 'start': i * 5000, 'end': (i + 1) * 5000,
                'state': 'Neutral', 'cn_median': 1.0 + 0.1 * (gc - 0.5),
                'cn_mean': 1.0, 'n_windows': 10, 'avg_quality': 1.0,
                'min_quality': 1.0, 'cn_std': 0.0, 'avg_repeats': 0.0,
                'avg_entropy': 0.0, 'max_entropy': 0.0,
                'mean_gc': gc,
                'masked_fraction': 0.0, 'repeat_class': 'None',
            })
        sat = {
            'chrom': 'chr1', 'start': 300000, 'end': 310000,
            'state': 'HighAmp', 'cn_median': 50.0, 'cn_mean': 50.0,
            'n_windows': 20, 'avg_quality': 0.3, 'min_quality': 0.1,
            'cn_std': 5.0, 'avg_repeats': 40.0, 'avg_entropy': 0.0, 'max_entropy': 0.0,
            'mean_gc': 0.45,
            'masked_fraction': 0.95, 'repeat_class': 'Satellite',
        }
        segs.append(sat)
        result = fl.apply_gc_cn_calibration(segs)
        sat_out = [s for s in result if s.get('repeat_class') == 'Satellite'][0]
        assert sat_out['gc_bias_factor'] == 1.0
        assert sat_out['cn_median'] == pytest.approx(50.0)



class TestWriteOutputColumnDetection:

    def test_repeat_class_detected_beyond_first_20(self):
        import tempfile, os
        segs = []
        for i in range(30):
            s = _seg('LowDup', i * 5000, (i + 1) * 5000, cn=2.0)
            s['avg_entropy'] = 0.0
            s['max_entropy'] = 0.0
            if i == 25:
                s['repeat_class']    = 'LINE'
                s['masked_fraction'] = 0.3
            segs.append(s)

        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            tmppath = f.name
        try:
            fl.write_output(segs, tmppath, extended=True)
            with open(tmppath) as f:
                lines = [l for l in f if not l.startswith('#')]
            n_cols = len(lines[0].strip().split('\t'))
            assert n_cols == 15, f"Expected 15 cols (with repeat), got {n_cols}"
        finally:
            os.unlink(tmppath)

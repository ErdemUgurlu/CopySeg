#!/usr/bin/env python3

import argparse
import sys
import os
import pandas as pd
import numpy as np
from datetime import datetime


_CHM13_VALIDATED_K = 72
_CHM13_RDNA_TOTAL  = 219

GROUND_TRUTH_CHM13 = [
    ("chr17",   7572544,   7591594, "TP53_Control",     1,   1.0,   "A",
     "Single-Copy Control",
     "TP53; single-copy, no paralogs; every 72-mer unique. The only fully "
     "k-robust CN=1 anchor. Source: NCBI Gene 7157"),
    ("chr1",  103504385, 103513421, "AMY1_Cluster",     7,   0.999, "A",
     "Segmental Dup",
     "AMY1; 7 copies >99.9% identical → aggregate≈physical=7. CHM13=H7.3 "
     "haplotype (NOT GRCh38's 3). Source: Yilmaz et al. 2024 Science; "
     "Megablast: 7 copies @>=99.9%"),
    ("chr16",    170722,    171564, "HBA1",             2,   1.0,   "A",
     "Segmental Dup",
     "HBA1+HBA2 alpha-globin, coding identical → aggregate≈2. NB: 842bp≈2 "
     "windows, subtelomeric → resolution/GC floor may under-call. "
     "Source: NCBI Gene 3039/3040"),

    ("chr12",   6544873,   6548727, "GAPDH",            1,   1.0,   "C",
     "Single-Copy Control",
     "GAPDH; 1 functional copy + 60+ processed pseudogenes (>80%). Exon k-mers "
     "inflate at k32; ≈1 at k72 (72-mers cross splice junctions). NOT a clean "
     "control. Source: Liu 2009; PLoS ONE 2012"),
    ("chr7",    5644671,   5648124, "ACTB",             1,   1.0,   "C",
     "Single-Copy Control",
     "ACTB; 1 functional copy + pseudogenes + ACTG1 paralog (~89-90% CDS). "
     "k32-inflated, ≈1 at k72. Source: PLoS ONE 2012"),
    ("chr7",   75860000,  75962983, "GTF2I",            1,   0.995, "C",
     "Single-Copy Control",
     "GTF2I; 1 functional copy in 7q11.23 Williams-Beuren LCRs; block-B "
     "duplicons ~99.5% → aggregate ~1.5-2.5 in LCR-overlap. Source: Cuscó 2008"),
    ("chr4",  193305857, 193328177, "FRG1",             1,   0.97,  "C",
     "Single-Copy Control",
     "FRG1; 1 functional copy but 23 dispersed paralogs in CHM13 (vs 9 GRCh38). "
     "Largest, most k-dependent inflation (~2.5-3.4 k72, ~6-9 k32). NOT a clean "
     "control. Source: Nurk 2022; Aganezov 2022 Science"),
    ("chr7",    6088690,   6126814, "PMS2",             1,   0.99,  "C",
     "Single-Copy Control",
     "PMS2; 1 functional copy; PMS2CL ~98-100% over exons 12-15 → local ~2, "
     "k-irreducible there. Mostly unique gene body → reads ~1. Source: De Vos 2004"),
    ("chr22",  42605990,  42610301, "CYP2D6",           1,   0.94,  "C",
     "Single-Copy Control",
     "CYP2D6; 1 copy (+CYP2D7/2D8P pseudogenes ~94%); CYP2D7 merges at k32. "
     "Source: NCBI Gene 1565 CHM13v2.0"),

    ("chr6",  161783172, 162011762, "LPA_KIV2_Array",  22,   0.96,  "B",
     "VNTR Array",
     "LPA KIV-2; 22 units in T2T-CHM13 (NOT GRCh38's 6). 95-98% identity → "
     "aggregate<physical at k72 (k-mer floor). Source: Behera 2024 (PMC11515395)"),
    ("chr17",  39044723,  39055625, "TBC1D3_Cluster",   9,   0.992, "B",
     "Segmental Dup",
     "TBC1D3; 9 copies >99% (gene conversion). aggregate ~identity^72 attenuated "
     "(~5 at k72) < physical 9; shorter k raises it. Pass-2 unique-only WORSENS "
     "(no unique k-mers). Source: Yoo 2024 Genome Res"),
    ("chr16",  22785024,  22814310, "NPIP",            26,   0.97,  "B",
     "Segmental Dup",
     "NPIP; 26 copies on chr16 (27 total), identity %97-99.6 (wide). Divergent "
     "old members carry unique k-mers → aggregate<<26. Source: Del Rosario/Eichler "
     "2025 Cell Genomics"),
    ("chr1",  145265708, 145345897, "NOTCH2NL_Family",  3,   0.997, "B",
     "Segmental Dup",
     "NOTCH2NL; scored q-arm window = 3 tandem copies (A/B/C). Full family 4-5 "
     "(NOTCH2/NOTCH2NLR on p-arm, geographically split). aggregate mixes p-arm. "
     "Source: Fiddes 2018 Cell; Vollger 2022"),
    ("chr1",  121194339, 121402237, "SRGAP2C",          4,   0.995, "B",
     "Segmental Dup",
     "SRGAP2C; 4 family paralogs (A/B/C/D) but INCOMPLETE dup: promoter+first ~9 "
     "exons ×4, gene body ×1 → position-dependent step, not flat 4. "
     "Source: Dennis 2012 Cell; Vollger 2022"),
    ("chrX",  152456258, 152470523, "OPN1MW",           3,   0.98,  "B",
     "Segmental Dup",
     "OPN1 array: 1 OPN1LW + 2 OPN1MW = 3 (CHM13 XX, haploid). >98% identity "
     "tandem; introns divergent. Source: medRxiv 2026 opsin assembly"),
    ("chr5",   71381729,  71423141, "SMN1",             1,   0.999, "B",
     "Segmental Dup",
     "SMN1; physical 1. SMN1↔SMN2 >99.9% (15 PSVs) → aggregate ~2, NOT k-mer "
     "splittable (no unique k-mers). measured ~2 = correct k-mer behavior, not "
     "error. Source: Blackburn 2024 (PMC11927269)"),
    ("chr5",   70791126,  70837821, "SMN2",             1,   0.999, "B",
     "Segmental Dup",
     "SMN2; physical 1. Same near-identical pool as SMN1 → aggregate ~2, not "
     "recoverable. Source: Blackburn 2024"),
    ("chr1",  103466373, 103474647, "AMY2A",            1,   0.99,  "B",
     "Segmental Dup",
     "AMY2A; physical 1 but shares exonic k-mers with 7 AMY1 + AMY2B → aggregate "
     ">physical (~4). NOT over-call. 8.3kb engulfed by AMY1 segment. "
     "Source: Yilmaz 2024; Bolognini 2024 Nature"),

    ("chr13",   5817416,   9348041, "rDNA_chr13",     219,   0.995, "D",
     "rDNA Array",
     "45S rDNA NOR chr13 (largest, ~80 units). physical_cn=219 is the 5-NOR "
     "TOTAL (Nurk 2022), NOT per-NOR — inter-NOR %99.4-99.7 → k-mers can't split. "
     "ddPCR 409±9 diploid (~205 haploid) corroborates total. Observatory"),
    ("chr14",   2099537,   2817811, "rDNA_chr14",     219,   0.995, "D",
     "rDNA Array",
     "45S rDNA NOR chr14 (smallest, ~15 units). 219 = 5-NOR total. Observatory"),
    ("chr15",   2506442,   4707485, "rDNA_chr15",     219,   0.995, "D",
     "rDNA Array",
     "45S rDNA NOR chr15. 219 = 5-NOR total. Observatory"),
    ("chr21",   3108298,   5612715, "rDNA_chr21",     219,   0.995, "D",
     "rDNA Array",
     "45S rDNA NOR chr21. 219 = 5-NOR total. Observatory"),
    ("chr22",   4793794,   5720650, "rDNA_chr22",     219,   0.995, "D",
     "rDNA Array",
     "45S rDNA NOR chr22. 219 = 5-NOR total. Observatory"),
    ("chr1",  227746662, 228024151, "5S_rDNA_chr1",   128,   0.99,  "D",
     "Tandem Array",
     "5S rDNA; 128 composite units (5S+AluY+2 subunits) %98-100. "
     "Source: Hoyt 2022 Science (PMC9301658). Observatory (compression)"),
    ("chr4",  193541579, 193543650, "DUX4",            65,   0.98,  "D",
     "Tandem Array",
     "D4Z4/DUX4; chr4q35≈33 + chr10q26≈32-33 = ~65 combined. 4q/10q ~98% → "
     "k-mer pools both arrays (paralog-collapse like TBC1D3). "
     "Source: Huang 2024 (PMC11092085). Observatory"),
    ("chr1",   16004103,  16075615, "NBPF1",         None,   0.96,  "D",
     "Segmental Dup",
     "NBPF1; physical GENE = 1, but signal is DUF1220/Olduvai family dosage "
     "(~16 paralogs %95-98). NO published T2T aggregate DUF1220 total; target is "
     "inherently k-dependent (k72~3, k32~16). Observatory. Source: O'Bleness 2012"),
    ("chr15",  32184804,  32243499, "GOLGA8A",       None,   0.98,  "D",
     "Segmental Dup",
     "GOLGA8A; ~15 GOLGA8 subfamily copies (GRCh38-era; no clean T2T count). "
     "core-duplicon → uncertain aggregate direction. Observatory (low confidence). "
     "Source: PMC6920530"),

    ("chrY",   6100000,   7031000, "TSPY_Array",     None,   None,  "skip",
     "Array (chrY)",
     "TSPY; chrY absent in CHM13 (XX) → SKIPPED"),
]

_HG002_VALIDATED_K = 72

GROUND_TRUTH_HG002 = [
    ("chr17",  7572544,   7591594, "TP53_Control",     1,   1.0,   "A",
     "Single-Copy Control",
     "TP53; single-copy, haplotype-invariant. Liftover +4887bp (original coords ok)"),
    ("chr1", 103828019, 103837055, "AMY1_Cluster",     3,   0.999, "A",
     "Segmental Dup",
     "AMY1; PROVISIONAL physical=3 (HG002 haplotype; CHM13=7). Needs independent "
     "HG002 physical GT. LIFTOVER from chr1:103504385 (offset +323634bp)"),
    ("chr1", 103790008, 103798282, "AMY2A",            1,   0.99,  "B",
     "Segmental Dup",
     "AMY2A; physical 1, aggregate>physical (AMY1 sharing). LIFTOVER offset +323635bp"),
    ("chr1", 149121462, 149204681, "NOTCH2NL_Family",  3,   0.997, "B",
     "Segmental Dup",
     "NOTCH2NL; PROVISIONAL physical=3 (HG002). LIFTOVER offset +3855754bp; "
     "old CHM13 coords gave CN=93 artefact"),
    ("chr6", 161783172, 162011762, "LPA_KIV2_Array",   1,   0.96,  "B",
     "VNTR Array",
     "LPA KIV-2; PROVISIONAL — HG002 minimal KIV-2 expansion (flanking unique "
     "dominates GT window). Original CHM13 coords (liftover hits non-LPA dups)"),
    ("chr17", 37556735,  37567643, "TBC1D3_Cluster",  11,   0.992, "B",
     "Segmental Dup",
     "TBC1D3; PROVISIONAL physical=11 (HG002 has more copies than CHM13's 9). "
     "aggregate<physical (id^k). LIFTOVER offset -1487988bp"),
    ("chr1",  16920193,  16985605, "NBPF1",         None,   0.96,  "D",
     "Segmental Dup",
     "NBPF1; physical gene=1, DUF1220 family dosage k-dependent → observatory. "
     "LIFTOVER offset +916090bp"),
    ("chr5",  71381729,  71423141, "SMN1",             1,   0.999, "D",
     "Segmental Dup",
     "SMN1; OBSERVATORY in HG002: liftover shows SMN1+SMN2 overlap at "
     "CM039015.1:72.42M — cannot evaluate separately (combined ≈6.4)"),
    ("chr5",  70791126,  70837821, "SMN2",             1,   0.999, "D",
     "Segmental Dup",
     "SMN2; OBSERVATORY: same SMN1+SMN2 overlap in HG002 assembly"),
    ("chr1", 231302374, 231545229, "5S_rDNA_chr1",    91,   0.99,  "D",
     "Tandem Array",
     "5S rDNA; PROVISIONAL ~91 units in HG002 (CHM13=128). Observatory. "
     "LIFTOVER offset +3555712bp"),
    ("chr13", 5411627,   5578556, "rDNA_chr13",     None,   None,  "skip",
     "rDNA Array", "rDNA; HG002 assembly incomplete → SKIPPED"),
    ("chr14", 1879084,   2069240, "rDNA_chr14",     None,   None,  "skip",
     "rDNA Array", "rDNA; HG002 assembly incomplete → SKIPPED"),
    ("chr15", 2325498,   2537498, "rDNA_chr15",     None,   None,  "skip",
     "rDNA Array", "rDNA; HG002 assembly incomplete → SKIPPED"),
    ("chr21", 3105606,   3341782, "rDNA_chr21",     None,   None,  "skip",
     "rDNA Array", "rDNA; HG002 assembly incomplete → SKIPPED"),
    ("chr22", 4860006,   5063528, "rDNA_chr22",     None,   None,  "skip",
     "rDNA Array", "rDNA; HG002 assembly incomplete → SKIPPED"),
]

SAMPLE_CONFIG = {
    'chm13': {
        'gt': GROUND_TRUTH_CHM13,
        'validated_k': _CHM13_VALIDATED_K,
        'rdna_total': _CHM13_RDNA_TOTAL,
        'label': 'CHM13 (female, XX)',
    },
    'hg002': {
        'gt': GROUND_TRUTH_HG002,
        'validated_k': _HG002_VALIDATED_K,
        'rdna_total': None,
        'label': 'HG002 (male, XY) — physical GT PROVISIONAL',
    },
}

_CHR_TO_CM039 = {
    "chr1":  "CM039011.1", "chr2":  "CM039012.1", "chr3":  "CM039013.1",
    "chr4":  "CM039014.1", "chr5":  "CM039015.1", "chr6":  "CM039016.1",
    "chr7":  "CM039017.1", "chr8":  "CM039018.1", "chr9":  "CM039019.1",
    "chr10": "CM039020.1", "chr11": "CM039021.1", "chr12": "CM039022.1",
    "chr13": "CM039023.1", "chr14": "CM039024.1", "chr15": "CM039025.1",
    "chr16": "CM039026.1", "chr17": "CM039027.1", "chr18": "CM039028.1",
    "chr19": "CM039029.1", "chr20": "CM039030.1", "chr21": "CM039031.1",
    "chr22": "CM039032.1", "chrX":  "CM039033.1",
}

_FULL_SCHEMA = [
    "chrom", "start", "end", "state", "cn_median", "cn_mean",
    "n_windows", "avg_quality", "min_quality", "cn_std",
    "avg_repeats", "avg_entropy", "max_entropy",
    "masked_fraction", "repeat_class",
    "gc_bias_factor", "segment_iqr", "boundary_conf",
]

def load_windows(bed_path: str) -> pd.DataFrame:
    names = ["chrom", "start", "end", "cn", "mean_count", "log_ratio",
             "num_kmers", "num_filtered"]
    df = pd.read_csv(bed_path, sep="\t", comment="#", header=None,
                     usecols=range(8), names=names)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    df["cn"]    = pd.to_numeric(df["cn"], errors="coerce")
    return df


def load_segments(bed_path: str, cn_column: str = "cn_median") -> pd.DataFrame:
    header_cols = None
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                header_cols = line.lstrip("#").rstrip("\n").split("\t")
                break
            else:
                break

    if header_cols is not None and 'cn_median' in header_cols:
        df = pd.read_csv(bed_path, sep="\t", comment="#", header=None,
                         names=header_cols)
    else:
        ncols = len(_FULL_SCHEMA)
        with open(bed_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                ncols = len(line.split("\t"))
                break
        names = _FULL_SCHEMA[:ncols]
        df = pd.read_csv(bed_path, sep="\t", comment="#", header=None,
                         names=names)

    df["start"]     = df["start"].astype(int)
    df["end"]       = df["end"].astype(int)
    df["cn_median"] = pd.to_numeric(df["cn_median"], errors="coerce")

    if cn_column != "cn_median":
        if cn_column not in df.columns:
            raise ValueError(
                f"--cn-column '{cn_column}' not present in segments BED "
                f"(have: {list(df.columns)}). Run refine_cn_unique.py first?")
        df["cn_median"] = pd.to_numeric(df[cn_column], errors="coerce")
        if cn_column == "cn_refined" and "state_refined" in df.columns:
            df["state"] = df["state_refined"]
    return df

def intersect_region(segs_df: pd.DataFrame, chrom_cm039: str,
                     gt_start: int, gt_end: int) -> pd.DataFrame:
    mask = (
        (segs_df["chrom"] == chrom_cm039) &
        (segs_df["end"]   > gt_start) &
        (segs_df["start"] < gt_end)
    )
    hits = segs_df[mask].copy()
    if hits.empty:
        return hits
    hits["ovlp_start"] = hits["start"].clip(lower=gt_start)
    hits["ovlp_end"]   = hits["end"].clip(upper=gt_end)
    hits["ovlp_len"]   = hits["ovlp_end"] - hits["ovlp_start"]
    return hits

def weighted_cn(hits: pd.DataFrame, outlier_percentile: float = 95.0) -> tuple:
    if hits.empty:
        return (None, None, 0, "NO_COVERAGE")

    total_bp  = hits["ovlp_len"].sum()
    peak_cn   = hits["cn_median"].max()
    dom_state = hits.loc[hits["ovlp_len"].idxmax(), "state"]

    cn_threshold = np.percentile(hits["cn_median"], outlier_percentile)
    robust_hits  = hits[hits["cn_median"] <= cn_threshold]

    if robust_hits.empty or robust_hits["ovlp_len"].sum() == 0:
        w_cn = (hits["cn_median"] * hits["ovlp_len"]).sum() / total_bp
    else:
        robust_bp = robust_hits["ovlp_len"].sum()
        w_cn = (robust_hits["cn_median"] * robust_hits["ovlp_len"]).sum() / robust_bp

    return (round(w_cn, 4), round(peak_cn, 4), int(total_bp), dom_state)

def window_cn(win_df: pd.DataFrame, chrom: str,
              gt_start: int, gt_end: int) -> float:
    mask = (
        (win_df["chrom"] == chrom) &
        (win_df["end"]   > gt_start) &
        (win_df["start"] < gt_end)
    )
    hits = win_df[mask]
    if hits.empty:
        return None
    return round(float(hits["cn"].median()), 4)

PASS_ERR = 0.35
MARGINAL_ERR = 0.50

def class_verdict(cn_class, measured, physical):
    if cn_class == "skip":
        return "SKIPPED"
    if measured is None:
        return "NO DATA"
    if cn_class == "D":
        return "OBSERVATORY"
    if physical is None or physical <= 0:
        return "OBSERVATORY"
    err = abs(measured - physical) / physical
    if cn_class in ("A", "C"):
        if err <= PASS_ERR:
            return "PASS (≤35%)"
        if err <= MARGINAL_ERR:
            return "MARGINAL (<50%)"
        return "FAIL (OVER)" if measured > physical else "FAIL (UNDER)"
    if cn_class == "B":
        if err <= PASS_ERR:
            return "AGG≈PHYS"
        return "AGG-LIMITED"
    return "OBSERVATORY"

def evaluate(segments_path: str, kmer_size: int = 72,
             sample: str = 'chm13', windows_path: str = None,
             cn_column: str = "cn_median") -> pd.DataFrame:
    segs = load_segments(segments_path, cn_column=cn_column)

    if cn_column == "cn_refined":
        if windows_path is not None:
            print("[EVAL] --cn-column cn_refined: suppressing --windows path "
                  "(window CN would read Pass-1 `cn`, shadowing Pass-2).",
                  file=sys.stderr)
            windows_path = None
        outlier_pct = 100.0
    else:
        outlier_pct = 95.0

    win_df = load_windows(windows_path) if windows_path else None
    cfg = SAMPLE_CONFIG[sample]
    gt_table = cfg['gt']

    seg_chroms = set(segs['chrom'].unique())
    _use_cm039 = any(c.startswith('CM0') for c in seg_chroms)
    if win_df is not None:
        win_chroms = set(win_df['chrom'].unique())
        _win_use_cm039 = any(c.startswith('CM0') for c in win_chroms)
    else:
        _win_use_cm039 = False

    rows = []
    for (chrom, start, end, gene, physical_cn, paralog_id, cn_class,
         category, notes) in gt_table:
        chrom_seg = _CHR_TO_CM039.get(chrom, chrom) if _use_cm039 else chrom
        chrom_win = _CHR_TO_CM039.get(chrom, chrom) if _win_use_cm039 else chrom

        gt_len = end - start
        hits = intersect_region(segs, chrom_seg, start, end)
        est_cn_seg, peak_cn, cov_bp, dom_state = weighted_cn(hits, outlier_percentile=outlier_pct)

        win_cn_val = window_cn(win_df, chrom_win, start, end) if win_df is not None else None
        est_cn = win_cn_val if win_cn_val is not None else est_cn_seg
        cn_source = "window" if win_cn_val is not None else "segment"

        verdict = class_verdict(cn_class, est_cn, physical_cn)

        if est_cn is None or physical_cn is None or (isinstance(physical_cn, (int, float)) and physical_cn <= 0):
            abs_err = None
            pct_err_str = "—"
        else:
            abs_err = round(abs(est_cn - physical_cn), 4)
            rel_err = round((est_cn - physical_cn) / physical_cn * 100, 1)
            pct_err_str = f"{rel_err:+.1f}%"

        n_segs  = len(hits) if not hits.empty else 0
        cov_pct = round(cov_bp / gt_len * 100, 1) if gt_len > 0 else 0
        id_str  = f"{paralog_id*100:.1f}%" if paralog_id is not None else "—"

        rows.append({
            "Region":       gene,
            "Category":     category,
            "Class":        cn_class,
            "Chr":          chrom,
            "GT_Start":     start,
            "GT_End":       end,
            "GT_Len_kb":    round(gt_len / 1e3, 1),
            "Physical_CN":  physical_cn,
            "Paralog_Id":   id_str,
            "Estimated_CN": est_cn,
            "Seg_CN":       est_cn_seg,
            "Peak_CN":      peak_cn,
            "Abs_Error":    abs_err,
            "Pct_Error":    pct_err_str,
            "N_Segs":       n_segs,
            "Cov_bp":       cov_bp,
            "Cov_Pct":      cov_pct,
            "Dom_State":    dom_state,
            "Verdict":      verdict,
            "CN_Source":    cn_source,
            "Notes":        notes,
        })
    return pd.DataFrame(rows)

_CLASS_DESC = {
    "A": "aggregate ≈ physical (scored directly)",
    "C": "aggregate > physical, Pass-2 recovers (scored)",
    "B": "aggregate ≠ physical, k-mer floor (characterized)",
    "D": "observatory (high-CN array / no clean target)",
}

def _fmt(v, nd=2):
    return f"{v:.{nd}f}" if isinstance(v, (int, float)) and pd.notna(v) else "—"

def build_markdown(df: pd.DataFrame, segments_path: str, kmer_size: int = 72,
                   sample: str = 'chm13', cn_column: str = "cn_median") -> str:
    run_name = os.path.basename(os.path.dirname(segments_path))
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    cfg = SAMPLE_CONFIG[sample]
    validated_k = cfg['validated_k']

    scored = df[df["Class"].isin(["A", "C"]) & df["Estimated_CN"].notna()]
    agglim = df[df["Class"] == "B"]
    obs    = df[df["Class"] == "D"]
    skipped = df[df["Class"] == "skip"]

    n_pass     = scored["Verdict"].str.startswith("PASS").sum()
    n_marginal = scored["Verdict"].str.startswith("MARGINAL").sum()
    n_fail     = scored["Verdict"].str.startswith("FAIL").sum()

    lines = []
    def ln(s=""): lines.append(s)

    ln("# CopySeg Ground Truth Evaluation — Independent Physical CN")
    ln()
    ln(f"**Run:** `{run_name}`  ")
    ln(f"**Sample:** {cfg['label']}  ")
    ln(f"**Segments:** `{segments_path}`  ")
    ln(f"**Measured column:** `{cn_column}` "
       f"({'aggregate' if cn_column == 'cn_median' else 'Pass-2 de-convolved (physical estimate)'})  ")
    ln(f"**k-mer length:** k={kmer_size} | **Reference:** T2T-CHM13v2.0 | **Generated:** {now}  ")
    ln()
    ln("> **GT = INDEPENDENT physical CN** (literature/assembly/alignment, NOT k-mer/Jellyfish "
       "— avoids circularity). The pipeline measures *aggregate* k-mer CN; whether it equals "
       "*physical* CN is a per-locus property of paralog identity, encoded as the achievability "
       "**Class** (A/B/C/D). Only A & C carry PASS/FAIL. See `docs/ground_truth_independent.md`.")
    ln()

    ln("## Summary")
    ln()
    ln("| Metric | Value |")
    ln("|--------|-------|")
    ln(f"| **PASS** (class A+C) | **{n_pass}** / {len(scored)} |")
    ln(f"| MARGINAL | {n_marginal} |")
    ln(f"| **FAIL** | **{n_fail}** |")
    ln(f"| Aggregate-limited (class B) | {len(agglim)} (characterized, not scored) |")
    ln(f"| Observatory (class D) | {len(obs)} |")
    ln(f"| Skipped | {len(skipped)} |")
    if not scored.empty and scored["Abs_Error"].notna().any():
        ln(f"| Median abs. error (scored) | {scored['Abs_Error'].median():.2f} CN |")
    ln()

    ln("## Class A & C — scored against physical CN")
    ln()
    ln("_A: aggregate≈physical (high-identity or single-copy). "
       "C: aggregate>physical from pseudogene sharing; Pass-2 unique-only should recover physical._")
    ln()
    ln("| Region | Cls | Phys.CN | Measured | Peak | %Err | Cov% | Dom.State | Verdict |")
    ln("|--------|-----|---------|----------|------|------|------|-----------|---------|")
    for _, r in df[df["Class"].isin(["A", "C"])].iterrows():
        ln(f"| {r['Region']} | {r['Class']} | {r['Physical_CN']} | "
           f"{_fmt(r['Estimated_CN'])} | {_fmt(r['Peak_CN'])} | {r['Pct_Error']} | "
           f"{r['Cov_Pct']}% | `{r['Dom_State']}` | {r['Verdict']} |")
    ln()

    ln("## Class B — aggregate-limited (characterized, not pass/fail)")
    ln()
    ln("_aggregate ≠ physical is EXPECTED here (near-identical sisters share all k-mers → "
       "aggregate>physical; divergent families share none → aggregate<physical). k-mers cannot "
       "recover physical. Reported for honesty; the measured value is the correct aggregate._")
    ln()
    ln("| Region | Phys.CN | Identity | Measured (agg) | Gap | Note |")
    ln("|--------|---------|----------|----------------|-----|------|")
    for _, r in agglim.iterrows():
        gap = "—"
        if pd.notna(r["Estimated_CN"]) and isinstance(r["Physical_CN"], (int, float)):
            gap = "under" if r["Estimated_CN"] < r["Physical_CN"] else "over"
            gap = f"{gap} ({r['Pct_Error']})"
        short = r["Notes"].split(".")[0]
        ln(f"| {r['Region']} | {r['Physical_CN']} | {r['Paralog_Id']} | "
           f"{_fmt(r['Estimated_CN'])} | {gap} | {short} |")
    ln()

    if not obs.empty:
        ln("## Class D — observatory (monitor only)")
        ln()
        rdna_total = cfg.get('rdna_total')
        if rdna_total:
            ln(f"_45S rDNA physical = {rdna_total} TOTAL across 5 NORs (Nurk 2022), not per-NOR "
               f"— inter-NOR ~99.5% identity aggregates them._")
            ln()
        ln("| Region | Phys.CN | Measured | Peak | Note |")
        ln("|--------|---------|----------|------|------|")
        for _, r in obs.iterrows():
            pcn = "k-dep" if pd.isna(r["Physical_CN"]) else r["Physical_CN"]
            short = r["Notes"].split(".")[0]
            ln(f"| {r['Region']} | {pcn} | {_fmt(r['Estimated_CN'])} | "
               f"{_fmt(r['Peak_CN'])} | {short} |")
        ln()

    ln("## Per-class rollup")
    ln()
    ln("| Class | Meaning | N | Median Phys | Median Meas |")
    ln("|-------|---------|---|-------------|-------------|")
    for cls in ["A", "C", "B", "D"]:
        grp = df[(df["Class"] == cls) & df["Estimated_CN"].notna()]
        if grp.empty:
            continue
        phys = grp["Physical_CN"].apply(lambda x: x if isinstance(x, (int, float)) else np.nan)
        mp = f"{phys.median():.1f}" if phys.notna().any() else "—"
        mm = f"{grp['Estimated_CN'].median():.1f}"
        ln(f"| {cls} | {_CLASS_DESC[cls]} | {len(grp)} | {mp} | {mm} |")
    ln()

    ln("---")
    ln(f"*Auto-generated by evaluate_ground_truth.py — {now}. GT source: "
       f"docs/ground_truth_independent.md (independent physical CN, non-circular).*")
    return "\n".join(lines)

def main():
    parser = argparse.ArgumentParser(description="CopySeg Ground Truth Evaluator (independent physical CN)")
    parser.add_argument("--segments", required=True, help="CopySeg BED file")
    parser.add_argument("--output", default=None, help="Output Markdown path")
    parser.add_argument("--kmer-size", type=int, default=72,
                        help="K-mer size of the run (for report labeling; default 72)")
    parser.add_argument("--sample", choices=list(SAMPLE_CONFIG.keys()), default="chm13",
                        help=f"Sample → GT table (choices: {', '.join(SAMPLE_CONFIG.keys())}; default chm13)")
    parser.add_argument("--windows", default=None,
                        help="Window-level CN BED (cn_w500.bed); uses window median within GT region.")
    parser.add_argument("--cn-column", default="cn_median",
                        help="Segment column carrying the CN to score. 'cn_median' = aggregate "
                             "(default); 'cn_refined' = Pass-2 de-convolved physical estimate (k<=50). "
                             "cn_refined auto-disables --windows.")
    args = parser.parse_args()

    if not os.path.exists(args.segments):
        print(f"ERROR: segments file not found: {args.segments}", file=sys.stderr)
        sys.exit(1)
    if args.windows and not os.path.exists(args.windows):
        print(f"ERROR: windows file not found: {args.windows}", file=sys.stderr)
        sys.exit(1)

    cfg = SAMPLE_CONFIG[args.sample]
    _use_windows = bool(args.windows) and args.cn_column != "cn_refined"
    cn_mode = "window-level" if _use_windows else f"segment-level ({args.cn_column})"
    print(f"[EVAL] Loading segments: {args.segments} "
          f"(k={args.kmer_size}, sample={args.sample}, CN mode: {cn_mode})", file=sys.stderr)
    if args.windows:
        print(f"[EVAL] Loading windows: {args.windows}", file=sys.stderr)
    df = evaluate(args.segments, kmer_size=args.kmer_size, sample=args.sample,
                  windows_path=args.windows, cn_column=args.cn_column)

    print("\n" + "=" * 80)
    print(f"CopySeg GT Evaluation — SUMMARY (CN: {cn_mode})")
    print("=" * 80)
    summary_cols = ["Region", "Class", "Physical_CN", "Estimated_CN",
                    "Peak_CN", "Pct_Error", "Verdict"]
    print(df[summary_cols].to_string(index=False))
    print("=" * 80)

    md = build_markdown(df, args.segments, kmer_size=args.kmer_size,
                        sample=args.sample, cn_column=args.cn_column)
    if args.output:
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        with open(args.output, "w") as f:
            f.write(md)
        print(f"\n[EVAL] Markdown report → {args.output}", file=sys.stderr)
    else:
        print("\n" + md)

    sys.exit(0)

if __name__ == "__main__":
    main()

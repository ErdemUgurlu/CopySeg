#!/bin/bash
#SBATCH --job-name=kmer_query
#SBATCH --partition=hi_end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/kmer_query_%j.out
#SBATCH --error=logs/kmer_query_%j.err

set -euo pipefail

# =============================================================================
# CopySeg — Single-Locus k-mer CN Query
#
# Mevcut Jellyfish DB'yi kullanarak tek bir lokusun k-mer CN'ini ölçer.
# Jellyfish DB yoksa hata verir (önce run_kmer_gt_validation.sh çalıştır).
#
# Kullanım:
#   sbatch scripts/query_kmer_cn.sh NPEPPS chr17:48392837-48485055
#   sbatch scripts/query_kmer_cn.sh SMN1 chr5:71381729-71423141
#   sbatch scripts/query_kmer_cn.sh MYC chr8:127735434-127742951
#
#   # k=32 ile:
#   KMER_SIZE=32 sbatch scripts/query_kmer_cn.sh NPEPPS chr17:48392837-48485055
#
#   # HG002 DB ile:
#   JF_DB=output/k32_gt_validation_hg002/hg002_diploid_k32.jf DIPLOID=2 \
#     sbatch scripts/query_kmer_cn.sh NPEPPS chr17:48392837-48485055
#
# Gereksinimler:
#   - Jellyfish DB (run_kmer_gt_validation.sh ile önceden oluşturulmuş)
#   - CHM13v2.0 FASTA (lokus sekansı çıkarmak için)
# =============================================================================

if [ $# -lt 2 ]; then
    echo "Usage: $0 <locus_name> <chr:start-end>"
    echo ""
    echo "Examples:"
    echo "  sbatch scripts/query_kmer_cn.sh NPEPPS chr17:48392837-48485055"
    echo "  sbatch scripts/query_kmer_cn.sh TP53 chr17:7572544-7591594"
    echo "  KMER_SIZE=32 sbatch scripts/query_kmer_cn.sh MYC chr8:127735434-127742951"
    exit 1
fi

LOCUS_NAME="$1"
REGION="$2"

KMER_SIZE="${KMER_SIZE:-72}"
DIPLOID="${DIPLOID:-1}"

echo "=================================================="
echo "CopySeg — k-mer CN Query: ${LOCUS_NAME}"
echo "=================================================="
echo "Region:    ${REGION}"
echo "K-mer:     ${KMER_SIZE}"
echo "Diploid:   ${DIPLOID} (1=haploid DB, 2=diploid DB → divide by 2)"
echo "Job ID:    ${SLURM_JOB_ID:-local}"
echo "Started:   $(date)"
echo "=================================================="
echo ""

# =============================================================================
# ENVIRONMENT
# =============================================================================

if command -v conda &>/dev/null; then
    eval "$(conda shell.bash hook)"
elif [ -f "/home/erdem.ugurlu/miniconda3/etc/profile.d/conda.sh" ]; then
    source "/home/erdem.ugurlu/miniconda3/etc/profile.d/conda.sh"
fi
conda activate konuseg

# =============================================================================
# PATHS
# =============================================================================

CHM13_FASTA="${CHM13_FASTA:-}"
if [ -z "${CHM13_FASTA}" ]; then
    for candidate in \
        "/home/erdem.ugurlu/KonuSeg_Phase2/data/chm13v2.0.fa" \
        "/home/klea/recomb/gra_bf/data/human/chm13/chm13v2.0.fa"; do
        if [ -f "${candidate}" ]; then
            CHM13_FASTA="${candidate}"
            break
        fi
    done
fi

if [ -z "${CHM13_FASTA}" ] || [ ! -f "${CHM13_FASTA}" ]; then
    echo "ERROR: CHM13 FASTA not found"
    exit 1
fi

# Jellyfish DB — auto-detect
JF_DB="${JF_DB:-}"
if [ -z "${JF_DB}" ]; then
    for candidate in \
        "output/k${KMER_SIZE}_gt_validation/chm13v2_k${KMER_SIZE}.jf" \
        "output/k${KMER_SIZE}_gt_validation_hg002/hg002_diploid_k${KMER_SIZE}.jf" \
        "output/k72_gt_validation/chm13v2_k72.jf"; do
        if [ -f "${candidate}" ]; then
            JF_DB="${candidate}"
            break
        fi
    done
fi

if [ -z "${JF_DB}" ] || [ ! -f "${JF_DB}" ]; then
    echo "ERROR: Jellyfish DB not found."
    echo "       Önce run_kmer_gt_validation.sh çalıştır."
    exit 1
fi

echo "CHM13 FASTA: ${CHM13_FASTA}"
echo "Jellyfish DB: ${JF_DB}"
echo ""

# =============================================================================
# QUERY
# =============================================================================

TMPDIR=$(mktemp -d)
REGION_FA="${TMPDIR}/${LOCUS_NAME}_region.fa"
KMERS_FA="${TMPDIR}/${LOCUS_NAME}_kmers.fa"
COUNTS_TXT="${TMPDIR}/${LOCUS_NAME}_counts.txt"

# 1. Extract region
echo "[1/4] Extracting region sequence..."
samtools faidx "${CHM13_FASTA}" "${REGION}" > "${REGION_FA}"
REGION_LEN=$(awk 'NR>1' "${REGION_FA}" | tr -d '\n' | wc -c)
echo "      Region length: ${REGION_LEN} bp"

# 2. Extract k-mers
echo "[2/4] Extracting ${KMER_SIZE}-mers..."
python3 -c "
seq = ''
with open('${REGION_FA}') as f:
    for line in f:
        if not line.startswith('>'): seq += line.strip().upper()
k = ${KMER_SIZE}
n = 0
with open('${KMERS_FA}', 'w') as out:
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if 'N' not in kmer:
            out.write(f'>{n}\n{kmer}\n')
            n += 1
print(f'      Extracted {n:,} valid {k}-mers')
"

# 3. Query
echo "[3/4] Querying Jellyfish DB..."
jellyfish query "${JF_DB}" -s "${KMERS_FA}" > "${COUNTS_TXT}"

# 4. Statistics
echo "[4/4] Computing statistics..."
echo ""
python3 -c "
import numpy as np

DIPLOID = ${DIPLOID}
counts = []
with open('${COUNTS_TXT}') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            counts.append(int(parts[1]))

arr = np.array(counts)
median = float(np.median(arr))
mean = float(np.mean(arr))
norm_median = median / DIPLOID
norm_mean = mean / DIPLOID

print('=' * 50)
print(f'  RESULT: ${LOCUS_NAME} (${REGION})')
print('=' * 50)
print(f'  K-mer size:      ${KMER_SIZE}')
print(f'  K-mers queried:  {len(arr):,}')
print(f'  Region length:   ${REGION_LEN} bp')
print()
print(f'  Raw median:      {median:.0f}')
print(f'  Raw mean:        {mean:.1f}')
print(f'  Raw IQR:         [{np.percentile(arr,25):.0f}, {np.percentile(arr,75):.0f}]')
print(f'  Raw range:       [{arr.min()}, {arr.max()}]')
print(f'  P5-P95:          [{np.percentile(arr,5):.0f}, {np.percentile(arr,95):.0f}]')
print()
if DIPLOID > 1:
    print(f'  Normalized CN:   {norm_median:.1f}  (raw / {DIPLOID})')
    print(f'  Normalized mean: {norm_mean:.1f}')
else:
    print(f'  ▶ k-mer CN = {median:.0f}')
print('=' * 50)
"

# Cleanup
rm -rf "${TMPDIR}"

echo ""
echo "Completed: $(date)"

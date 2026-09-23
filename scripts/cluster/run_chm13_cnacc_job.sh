#!/bin/bash
#SBATCH --job-name=copyseg_chm13
#SBATCH --partition=hi_end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/copyseg_chm13_%j.out
#SBATCH --error=logs/copyseg_chm13_%j.err


export SEX="${SEX:-XX}"
export SAMPLE="${SAMPLE:-chm13}"

exec scripts/cluster/run_copyseg_pipeline.sh "chm13_${1:-run_$(date +%Y%m%d_%H%M)}"

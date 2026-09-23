#!/bin/bash
#SBATCH --job-name=copyseg_hg002
#SBATCH --partition=hi_end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=logs/copyseg_hg002_%j.out
#SBATCH --error=logs/copyseg_hg002_%j.err


export SEX="${SEX:-XY}"
export SAMPLE="${SAMPLE:-hg002}"

exec scripts/cluster/run_copyseg_pipeline.sh "hg002_${1:-run_$(date +%Y%m%d_%H%M)}"

#!/bin/bash
#SBATCH --job-name=lc025_frac
#SBATCH --output=log/%x_%A_%a.log
#SBATCH --error=log/%x_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --partition=icpu-stocker
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --array=1992-2020
#SBATCH --mail-type=END,FAIL

set -euo pipefail

mkdir -p log

module purge
module load R/4.4.2-gfbf-2024a UDUNITS/2.2.28-GCCcore-13.3.0 \
PROJ/9.4.1-GCCcore-13.3.0 GDAL/3.10.0-foss-2024a CDO/2.4.4-gompi-2024a

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export GDAL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

export LC_YEAR="${SLURM_ARRAY_TASK_ID}"
export REMAKE_ALL=FALSE

Rscript 12_make_lc025_fractions.R

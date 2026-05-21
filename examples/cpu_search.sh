#!/usr/bin/env bash
#SBATCH -J PAIPAI_search_cpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=31
#SBATCH -t 72:00:00
#SBATCH --partition=cpu
#SBATCH --output=paipai_search_cpu.log
#SBATCH --mem=80G

set -euo pipefail

module purge
module load WebProxy
module load Anaconda3/2024.02-1

source /sw/eb/sw/Anaconda3/2024.02-1/etc/profile.d/conda.sh
conda activate materialsframework-main
export PAIPAI_PYTHON="$CONDA_PREFIX/bin/python"
export PATH="$CONDA_PREFIX/bin:$PATH"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV"
echo "CONDA_PREFIX=$CONDA_PREFIX"
echo "PAIPAI_PYTHON=$PAIPAI_PYTHON"

paipai \
  --input struc.in \
  --mode search \
  --device cpu \
  --fast 15 \
  --slow 15 \
  --pool-cap 128 \
  --steps 2000 \
  --temp 700 \
  --p-swap-metal 70 \
  --p-swap-inter 30 \
  --p-cluster-inter 0 \
  --p-exch-metal 0 \
  --p-exch-inter 0 \
  --intsite-neighbor-cutoff 3.5

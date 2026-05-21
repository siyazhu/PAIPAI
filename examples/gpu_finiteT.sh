#!/usr/bin/env bash
#SBATCH -J PAIPAI_finiteT_gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 72:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --output=paipai_finiteT.log
#SBATCH --mem=80G

set -euo pipefail

module purge
module load WebProxy
module load Anaconda3/2024.02-1
module load cuDNN/8.9.2.26-CUDA-12.1.1

source /sw/eb/sw/Anaconda3/2024.02-1/etc/profile.d/conda.sh
conda activate materialsframework-main
export PAIPAI_PYTHON="$CONDA_PREFIX/bin/python"
export PATH="$CONDA_PREFIX/bin:$PATH"

export XLA_FLAGS=--xla_gpu_cuda_data_dir=/sw/eb/sw/CUDA/12.1.1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
echo "CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV"
echo "PAIPAI_PYTHON=$PAIPAI_PYTHON"
echo "which python=$(which python)"

nvidia-smi

paipai \
  --input struc.in \
  --mode finiteT \
  --device cuda \
  --ngpu 1 \
  --fast 0 \
  --slow 1 \
  --steps 200000 \
  --temp 700 \
  --p-swap-metal 70 \
  --p-swap-inter 30 \
  --p-cluster-inter 0 \
  --p-exch-metal 0 \
  --p-exch-inter 0 \
  --intsite-neighbor-cutoff 3.5

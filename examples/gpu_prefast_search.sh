#!/usr/bin/env bash
#SBATCH -J PAIPAI_prefast_gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 72:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:4
#SBATCH --output=paipai_prefast_search.log
#SBATCH --mem=200G

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

nvidia-smi

paipai \
  --input struc.in \
  --mode search \
  --root search_prefast \
  --device cuda \
  --ngpu 4 \
  --fast 2 \
  --slow 2 \
  --pool-cap 128 \
  --steps 200000 \
  --temp 50 \
  --fmax-screen 0.1 \
  --max-steps-screen 100 \
  --fmax-refine 0.01 \
  --max-steps-refine 500 \
  --p-swap-metal 60 \
  --p-swap-inter 10 \
  --p-hop-inter 0 \
  --p-cluster-inter 30 \
  --p-exch-metal 0 \
  --p-exch-inter 0 \
  --intsite-neighbor-cutoff 2.5 \
  --intsite-hop-cutoff 2.0 \
  --prefast on \
  --prefast-warmup-steps 1000 \
  --prefast-candidates-per-slot 6 \
  --prefast-diagnostics summary

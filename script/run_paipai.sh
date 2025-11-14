#!/usr/bin/env bash
set -euo pipefail

########################################
# User-configurable section
########################################

# Initial structure file (readable by your Structure class)
INPUT_STRUC="struc.in"

# C++ Monte Carlo master executable
MASTER="./mc_paipai"

# Python interpreter for workers
PYTHON="python"        # or /scratch/.../orbenv/bin/python

# ML potential model and device
MODEL="GRACE-2L-OMAT"
DEVICE="cpu"           # use "cuda" if running on a GPU node

# Number of fast and slow workers
NUM_FAST=15
NUM_SLOW=15

# Monte Carlo parameters
MC_STEPS=2000
MC_TEMP=0.001

# Maximum capacity of the waiting pool (shared candidate buffer)
POOL_CAP=128

# Root working directory (all subfolders will be created here)
ROOT_DIR="$(pwd)"

########################################
# Initialization and directory setup
########################################
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[INFO] ROOT_DIR = $ROOT_DIR"
echo "[INFO] INPUT_STRUC = $INPUT_STRUC"
echo "[INFO] NUM_FAST = $NUM_FAST, NUM_SLOW = $NUM_SLOW"
echo "[INFO] MODEL = $MODEL, DEVICE = $DEVICE"
echo "[INFO] MC_STEPS = $MC_STEPS, MC_TEMP = $MC_TEMP"
echo "[INFO] POOL_CAP = $POOL_CAP"

# Create required subdirectories
mkdir -p "$ROOT_DIR/fast" \
         "$ROOT_DIR/waiting_pool" \
         "$ROOT_DIR/waiting_work" \
         "$ROOT_DIR/refine_outbox" \
         "$ROOT_DIR/reports" \
         "$ROOT_DIR/counters" \
         "$ROOT_DIR/mcprocess"

# Cleanup handler – kill all background jobs on exit or interruption
cleanup() {
    echo "[INFO] Cleaning up: killing all background workers..."
    kill 0 2>/dev/null || true
}
trap cleanup EXIT INT TERM

########################################
# Launch fast workers
########################################

echo "[INFO] Launching $NUM_FAST fast workers..."
for (( i=1; i<=NUM_FAST; ++i )); do
    LOG_FAST="fast_${i}.log"
    echo "  -> fast_worker slot $i (log: $LOG_FAST)"
    $PYTHON fast_worker.py \
        --slot "$i" \
        --root "$ROOT_DIR" \
        --model "$MODEL" \
        --device "$DEVICE" \
        --pool_cap "$POOL_CAP" \
        >"$LOG_FAST" 2>&1 &
done

########################################
# Launch slow workers
########################################

echo "[INFO] Launching $NUM_SLOW slow workers..."
for (( j=1; j<=NUM_SLOW; ++j )); do
    ID=$(printf "slow-%02d" "$j")
    LOG_SLOW="slow_${ID}.log"
    echo "  -> slow_worker $ID (log: $LOG_SLOW)"
    $PYTHON slow_worker.py \
        --worker-id "$ID" \
        --root "$ROOT_DIR" \
        --model "$MODEL" \
        --device "$DEVICE" \
        >"$LOG_SLOW" 2>&1 &
done

########################################
# Run the C++ Monte Carlo master
########################################

echo "[INFO] Running C++ master ..."
echo "       $MASTER $INPUT_STRUC --workers $NUM_FAST --steps $MC_STEPS --temp $MC_TEMP"
echo "       (MC progress will be written to mc.log)"

# Note: --workers must match the number of fast workers (slots)
$MASTER "$INPUT_STRUC" \
    --workers "$NUM_FAST" \
    --steps "$MC_STEPS" \
    --temp "$MC_TEMP"

echo "[INFO] C++ master finished. Terminating workers..."

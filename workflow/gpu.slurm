#!/bin/bash
#SBATCH --job-name=sc-gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpuk[001-018]
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.log

# GPU job template for this repo. Usage: sbatch workflow/gpu.slurm src/<script>.py [args]
# Amarel notes (verified 2026-09-22, see docs/01):
#   --exclude gpuk*   : those nodes read home so slowly that `import torch` times out
#   LD_LIBRARY_PATH   : system libstdc++ is too old for conda's scipy once torch is loaded
#   CUDA_VISIBLE_DEVICES: some nodes do not isolate cards; the echo confirms we got exactly one

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

source /cache/home/xw347/miniforge3/etc/profile.d/conda.sh
conda activate sc
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

echo "node: $(hostname)  CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}  start: $(date)"
python "$@"
echo "end: $(date)"

#!/bin/bash
# Launch the 256^3 3D flagship on the HPC (32 FFT threads), detached.
# Part of the user-approved 48-core budget: 16 matrix workers + 32 here.
# Run on the remote from the project root.
set -eu
cd "$(dirname "$0")/.."
PY=/home/rc/anaconda3/bin/python3
[ -x "$PY" ] || PY=$(command -v python3)
mkdir -p results_hpc
{ PFC_FFT_THREADS=32 OMP_NUM_THREADS=1 nohup "$PY" hpc/hpc_run_one.py \
    --kind poly3d --n 256 --r -0.25 --relax 200 --seed 7 \
    --strain-to 0.06 --out results_hpc/p3d_flagship_256 \
    > flagship_256.log 2>&1 < /dev/null & }
sleep 5
echo "flagship launched:"
head -2 flagship_256.log

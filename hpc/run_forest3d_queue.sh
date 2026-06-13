#!/bin/bash
# HPC dilute-forest 3D queue: 256^3, ~20 cells, dilute NF so dislocation lines
# are spatially resolved (the clean junction-hardening verdict the 96^3 local
# box could not give). Run detached from the project root on the HPC.
#   nohup bash hpc/run_forest3d_queue.sh > forest3d.log 2>&1 < /dev/null &
set -u
cd "$(dirname "$0")/.."
PY=/home/rc/anaconda3/bin/python3
[ -x "$PY" ] || PY=$(command -v python3)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# 256^3 is heavy: give each worker several FFT threads, few concurrent.
THREADS=${THREADS:-8}
mkdir -p results_hpc logs
run_cell() {
  local nf=$1 seed=$2
  local tag="f3d_n256c20_nf${nf}_s${seed}"
  [ -f "results_hpc/${tag}/summary.json" ] && { echo "skip $tag"; return; }
  echo "start $tag"
  FOREST3D_N=256 FOREST3D_CELLS=20 FOREST3D_NF=$nf FOREST3D_SEED=$seed \
    FOREST3D_TAG=$tag PFC_FFT_THREADS=$THREADS \
    OUT_BASE=results_hpc "$PY" src/run_forest_3d.py > "logs/${tag}.log" 2>&1
  # relocate output (driver writes results/forest_3d/<tag>; move to results_hpc)
  [ -d "results/forest_3d/${tag}" ] && mv "results/forest_3d/${tag}" "results_hpc/${tag}"
  echo "end $tag rc=$?"
}
# 6 cells; run 6 concurrent at THREADS each (<=48 cores)
for nf in 0 4 8; do for s in 7 11; do run_cell $nf $s & done; done
wait
echo "forest3d queue done: $(date)"

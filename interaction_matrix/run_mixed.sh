#!/bin/bash
# Phase 2b: direct measurement of the STEM-inventory forest alpha (validate the
# 0.69 prediction). 4 densities x 8 seeds = 32 jobs, big-first, saturating.
set -u
cd ~/BO/interaction_matrix
export PYTHONPATH=~/BO/exadis_src/python
R=mixed; rm -rf "$R"; mkdir -p "$R"; TH=3; MAXJ=32
SEEDS="1234 5678 2222 3333 4444 5555 6666 7777"
COMMON="LBOX=5000 NPROBE=2 ERATE=1e4 MAX_STRAIN=0.0012 FLOW_LO=0.0007 MAXSEG=300"
: > "$R/jobs.txt"
for nf in 32 24 16 8; do
  for s in $SEEDS; do
    echo "OMP_NUM_THREADS=$TH SEED=$s NFOREST=$nf $COMMON OUT=$R/n${nf}_s${s} python3 -u build_mixed.py > $R/n${nf}_s${s}.log 2>&1" >> "$R/jobs.txt"
  done
done
echo "MIXED JOBS: $(wc -l < $R/jobs.txt)  (concurrency=$MAXJ x $TH = $((MAXJ*TH)) cores)"
while IFS= read -r line; do
  eval "$line" &
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do wait -n; done
done < "$R/jobs.txt"
wait
echo "MIXED ALL_DONE"

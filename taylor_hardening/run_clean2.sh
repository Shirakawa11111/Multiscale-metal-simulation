#!/bin/bash
# High-utilization 8-seed clean-protocol campaign.
# Lesson: small/medium DDD jobs scale to only ~few cores, and heterogeneous job
# durations cause stragglers that idle the node. Fix: a big saturated pool,
# big-jobs-FIRST (so long NL=300 runs overlap instead of straggling alone),
# moderate threads, high concurrency.
set -u
cd ~/BO/taylor_hardening
export PYTHONPATH=~/BO/exadis_src/python
R=clean2
rm -rf "$R"; mkdir -p "$R"
TH=6; MAXJ=18
SEEDS="1234 5678 2222 3333 4444 5555 6666 7777"
COMMON="LBOX=16000 N_PROBE=8 ERATE=1e4 MAX_STRAIN=0.0008 FLOW_LO=0.0005 NGRID=48 MAXSEG=600 MINSEG=150"
: > "$R/jobs.txt"
for nl in 300 200 100 50; do
  for s in $SEEDS; do
    echo "OMP_NUM_THREADS=$TH SEED=$s NUM_LINES=$nl $COMMON OUT=$R/F_n${nl}_s${s} python3 -u forest_probe.py > $R/F_n${nl}_s${s}.log 2>&1" >> "$R/jobs.txt"
  done
done
for s in $SEEDS; do
  echo "OMP_NUM_THREADS=$TH SEED=$s NUM_LINES=8 $COMMON OUT=$R/B_s${s} python3 -u forest_probe.py > $R/B_s${s}.log 2>&1" >> "$R/jobs.txt"
done
echo "CLEAN2 JOBS: $(wc -l < "$R/jobs.txt")  (concurrency=$MAXJ x $TH threads = $((MAXJ*TH)) cores)"
while IFS= read -r line; do
  eval "$line" &
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do wait -n; done
done < "$R/jobs.txt"
wait
echo "CLEAN2 ALL_DONE"

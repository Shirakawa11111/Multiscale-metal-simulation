#!/bin/bash
# Rate-extrapolation of the clean (K=8) forest-hardening coefficient toward
# quasi-static, to pin down alpha and remove the residual rate drift
# (0.74@1e4 vs 0.85@5e3). High-utilization: many jobs x few threads, big-first
# (slow low-rate/high-density jobs launched first so they overlap, no straggler).
set -u
cd ~/BO/taylor_hardening
export PYTHONPATH=~/BO/exadis_src/python
R=rateK8
rm -rf "$R"; mkdir -p "$R"
TH=4; MAXJ=28
SEEDS="1234 5678 2222"
RATES="3e3 5e3 1e4"
COMMON="LBOX=16000 N_PROBE=8 MAX_STRAIN=0.0008 FLOW_LO=0.0005 NGRID=48 MAXSEG=600 MINSEG=150"
: > "$R/jobs.txt"
# big-first: low rate + high density first (slowest), descending
for e in $RATES; do
  for nl in 200 100 50; do
    for s in $SEEDS; do
      echo "OMP_NUM_THREADS=$TH SEED=$s NUM_LINES=$nl ERATE=$e $COMMON OUT=$R/F_e${e}_n${nl}_s${s} python3 -u forest_probe.py > $R/F_e${e}_n${nl}_s${s}.log 2>&1" >> "$R/jobs.txt"
    done
  done
done
for e in $RATES; do for s in $SEEDS; do
  echo "OMP_NUM_THREADS=$TH SEED=$s NUM_LINES=8 ERATE=$e $COMMON OUT=$R/B_e${e}_s${s} python3 -u forest_probe.py > $R/B_e${e}_s${s}.log 2>&1" >> "$R/jobs.txt"
done; done
echo "RATEK8 JOBS: $(wc -l < "$R/jobs.txt")  (concurrency=$MAXJ x $TH = $((MAXJ*TH)) cores)"
while IFS= read -r line; do
  eval "$line" &
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do wait -n; done
done < "$R/jobs.txt"
wait
echo "RATEK8 ALL_DONE"

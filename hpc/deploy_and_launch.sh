#!/bin/bash
# One-shot deployment: verify >=64 free CPUs on the HPC (user policy
# 2026-06-12: start only when 64 CPUs are available, use 64), rsync code,
# launch the matrix queue detached. Local usage (macOS bash 3.2 safe):
#   bash hpc/deploy_and_launch.sh
#
# Exit codes: 0 launched · 1 probe/ssh failure · 2 fewer than 64 CPUs free
#             3 matrix already running
set -eu
HOST=rc@192.168.196.40
ROOT=/home/rc/BO/pfc_matrix
SSH="ssh -o ConnectTimeout=8 -o BatchMode=yes $HOST"

# --- free-CPU gate -----------------------------------------------------
# Instantaneous occupancy, not ps lifetime averages: sample /proc/stat
# twice 3s apart for true idle fraction AND take 1-min loadavg; gate on
# the conservative minimum, rounding busy UP.
PROBE='n=$(nproc)
read -r _ u1 n1 s1 i1 w1 _ < /proc/stat
t1=$((u1+n1+s1+i1+w1))
sleep 3
read -r _ u2 n2 s2 i2 w2 _ < /proc/stat
t2=$((u2+n2+s2+i2+w2))
dt=$((t2-t1)); di=$(( (i2+w2)-(i1+w1) ))
idle_cores=$(( dt>0 ? n*di/dt : 0 ))
load=$(awk "{print int(\$1+0.999)}" /proc/loadavg)
free_by_load=$(( n-load ))
free=$(( idle_cores<free_by_load ? idle_cores : free_by_load ))
echo "$n $free"'
OUT=$($SSH "$PROBE") || { echo "ABORT: ssh probe to $HOST failed (rc=$?)"; exit 1; }
read -r NPROC FREE <<< "$OUT"
case "$NPROC" in ''|*[!0-9]*) echo "ABORT: bad probe output: '$OUT'"; exit 1;; esac
case "$FREE" in ''|*[!0-9-]*) echo "ABORT: bad probe output: '$OUT'"; exit 1;; esac
echo "HPC: $NPROC cores, ~$FREE free (instantaneous, conservative)"
if [ "$FREE" -lt 64 ]; then
    echo "ABORT: fewer than 64 free CPUs (policy: start only at >=64)"
    exit 2
fi

if $SSH "pgrep -f hpc_run_one.py > /dev/null"; then
    echo "ABORT: matrix workers already running on HPC"
    exit 3
fi

# --- deploy + launch ---------------------------------------------------
$SSH "mkdir -p $ROOT"
rsync -az -e "ssh -o ConnectTimeout=8 -o BatchMode=yes" \
    --exclude results --exclude '*.npz' src hpc "$HOST:$ROOT/"
# Background ONLY the nohup command ({ ... & } grouping); cd/chmod failures
# must stay foreground and fatal, and tail must run in $ROOT.
$SSH "cd $ROOT && chmod +x hpc/run_queue.sh \
  && { nohup bash hpc/run_queue.sh > queue.log 2>&1 < /dev/null & } \
  && sleep 2 && echo launched && tail -3 queue.log"
echo "monitor:  ssh $HOST 'tail -5 $ROOT/queue.log; ls $ROOT/results_hpc | wc -l'"
echo "pull:     rsync -az --include='*/' --include='summary.json' --include='*.log' --exclude='*' $HOST:$ROOT/results_hpc/ results_hpc_pulled/"

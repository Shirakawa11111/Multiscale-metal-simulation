#!/bin/bash
# One-shot deployment: verify >=64 free CPUs on the HPC (user policy
# 2026-06-12: start only when 64 CPUs are available, use 64), rsync code,
# launch the matrix queue detached. Local usage:
#   bash hpc/deploy_and_launch.sh
set -eu
HOST=rc@192.168.196.40
ROOT=/home/rc/BO/pfc_matrix
SSH="ssh -o ConnectTimeout=8 -o BatchMode=yes $HOST"

read -r NPROC BUSY <<< "$($SSH "echo \$(nproc) \$(ps -eo pcpu --no-headers | awk '{s+=\$1} END {printf \"%d\", s/100}')")"
FREE=$((NPROC - BUSY))
echo "HPC: $NPROC cores, ~$BUSY busy, ~$FREE free"
if [ "$FREE" -lt 64 ]; then
    echo "ABORT: fewer than 64 free CPUs (policy: start only at >=64)"
    exit 2
fi

if $SSH "pgrep -f hpc_run_one.py > /dev/null"; then
    echo "ABORT: matrix workers already running on HPC"
    exit 3
fi

$SSH "mkdir -p $ROOT"
rsync -az --exclude results --exclude '*.npz' src hpc "$HOST:$ROOT/"
$SSH "cd $ROOT && chmod +x hpc/run_queue.sh && nohup bash hpc/run_queue.sh > queue.log 2>&1 < /dev/null & sleep 1; echo launched; tail -3 queue.log"
echo "monitor:  ssh $HOST 'tail -5 $ROOT/queue.log; ls $ROOT/results_hpc | wc -l'"
echo "pull:     rsync -az --include='*/' --include='summary.json' --include='*.log' --exclude='*' $HOST:$ROOT/results_hpc/ results_hpc_pulled/"

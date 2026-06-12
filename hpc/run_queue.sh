#!/bin/bash
# Remote queue runner: executes manifest lines with at most $MAXP concurrent
# single-threaded workers. Run from the project root on the HPC:
#   nohup bash hpc/run_queue.sh > queue.log 2>&1 < /dev/null &
set -u
cd "$(dirname "$0")/.."

MAXP=${MAXP:-64}
export PFC_FFT_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

PY=/home/rc/anaconda3/bin/python3
[ -x "$PY" ] || PY=$(command -v python3)
"$PY" -c "import numpy, scipy" || { echo "FATAL: no numpy/scipy in $PY"; exit 1; }
export PY

mkdir -p results_hpc logs
echo "queue start: $(wc -l < hpc/manifest.txt) runs, maxp=$MAXP, py=$PY"

run_one() {
    local line="$1"
    local name="${line%%|*}"
    local args="${line#*|}"
    if [ -f "results_hpc/${name}/summary.json" ]; then
        echo "skip $name (done)"
        return 0
    fi
    echo "start $name"
    # shellcheck disable=SC2086
    "$PY" hpc/hpc_run_one.py $args > "logs/${name}.log" 2>&1
    echo "end $name rc=$?"
}
export -f run_one

xargs -d '\n' -P "$MAXP" -I{} bash -c 'run_one "$@"' _ {} < "${MANIFEST:-hpc/manifest.txt}"
echo "queue done: $(date)"

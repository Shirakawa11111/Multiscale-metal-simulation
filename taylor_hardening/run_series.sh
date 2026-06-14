#!/bin/bash
# Drive the Taylor density series: ONE fresh python process per density
# (creating multiple ExaDiS networks in one process triggers a pyexadis
# double-free). Then aggregate into the Taylor fit.
set -u
cd ~/BO/taylor_hardening
export PYTHONPATH=~/BO/exadis_src/python
NUM_LINES_LIST="${NUM_LINES_LIST:-10 20 40 80}"
ERATE="${ERATE:-3e4}"
MAX_STRAIN="${MAX_STRAIN:-0.003}"
FLOW_LO="${FLOW_LO:-0.002}"
OUT="${OUT:-taylor_out}"
OMP="${OMP:-48}"
rm -rf "$OUT"; mkdir -p "$OUT"
echo "DRIVER START: densities=[$NUM_LINES_LIST] erate=$ERATE max_strain=$MAX_STRAIN"
for NL in $NUM_LINES_LIST; do
  echo "===== DENSITY num_lines=$NL ====="
  OMP_NUM_THREADS=$OMP NUM_LINES=$NL ERATE=$ERATE MAX_STRAIN=$MAX_STRAIN \
    FLOW_LO=$FLOW_LO NGRID=32 MAXSEG=600 MINSEG=150 OUT=$OUT \
    python3 -u taylor_series.py 2>&1 | grep -E "DENSITY|rho0=|FAILED|TAYLOR|free\(\)|Error"
done
echo "===== AGGREGATE ====="
python3 -u aggregate_taylor.py "$OUT"
echo "DRIVER ALL_DONE"

#!/bin/bash
# Forest-probe Taylor series: vary num_lines -> vary pinned forest density
# rho_f. One fresh process per density (avoids the pyexadis double-free). The
# forest is PINNED so rho stays ~ rho_f (no annihilation) -> the density lever
# across runs survives, unlike the free-line initial-density series.
set -u
cd ~/BO/taylor_hardening
export PYTHONPATH=~/BO/exadis_src/python
NUM_LINES_LIST="${NUM_LINES_LIST:-20 40 80 160}"
ERATE="${ERATE:-3e4}"
MAX_STRAIN="${MAX_STRAIN:-0.003}"
FLOW_LO="${FLOW_LO:-0.0015}"
PROBE_FRAC="${PROBE_FRAC:-0.25}"
OUT="${OUT:-forest_out}"
OMP="${OMP:-48}"
rm -rf "$OUT"; mkdir -p "$OUT"
echo "FOREST DRIVER START: num_lines=[$NUM_LINES_LIST] erate=$ERATE probe_frac=$PROBE_FRAC"
for NL in $NUM_LINES_LIST; do
  echo "===== FOREST num_lines=$NL ====="
  OMP_NUM_THREADS=$OMP NUM_LINES=$NL ERATE=$ERATE MAX_STRAIN=$MAX_STRAIN \
    FLOW_LO=$FLOW_LO PROBE_FRAC=$PROBE_FRAC NGRID=32 MAXSEG=600 MINSEG=150 \
    OUT=$OUT python3 -u forest_probe.py 2>&1 | \
    grep -E "num_lines=|rho_forest|Error|Traceback|free\(\)|FAILED"
done
echo "===== AGGREGATE ====="
python3 -u aggregate_forest.py "$OUT"
echo "FOREST ALL_DONE"

#!/bin/bash
# Faithful-matrix attempt with PINNED-END bowing probe (Frank-Read-like) + low
# rate + peak/critical stress. One representative pair per junction type x
# densities x seeds. Tests whether collinear becomes dominant (canonical).
set -u; cd ~/BO/interaction_matrix; export PYTHONPATH=~/BO/exadis_src/python
R=matrix_src; TH=3; MAXJ=36
rm -rf "$R"; mkdir -p "$R"
# representative (MSYS,FSYS,type) pairs in ExaDiS ordering
PAIRS="0,3,collinear 0,8,Lomer 0,6,Hirth 0,4,glissile 0,1,coplanar 0,0,self"
COMMON="LBOX=5000 N_PROBE=1 ERATE=3e3 MAX_STRAIN=0.0005 FLOW_LO=0.0003 MAXSEG=300"
: > "$R/jobs.txt"
for nf in 32 16 8; do
  for p in $PAIRS; do
    m=${p%%,*}; rest=${p#*,}; f=${rest%%,*}; t=${rest#*,}
    for s in 1234 5678 2222; do
      echo "OMP_NUM_THREADS=$TH MSYS=$m FSYS=$f NFOREST=$nf SEED=$s $COMMON OUT=$R/${t}_n${nf}_s${s} python3 -u build_src.py > $R/${t}_n${nf}_s${s}.log 2>&1" >> "$R/jobs.txt"
    done
  done
done
echo "SRC JOBS: $(wc -l < $R/jobs.txt) (conc=$MAXJ x $TH=$((MAXJ*TH)) cores)"
while IFS= read -r line; do
  eval "$line" &
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do wait -n; done
done < "$R/jobs.txt"
wait
echo "SRC ALL_DONE"

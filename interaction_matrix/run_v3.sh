#!/bin/bash
set -u; cd ~/BO/interaction_matrix; export PYTHONPATH=~/BO/exadis_src/python
R=matrix_v3; TH=3; MAXJ=36
rm -rf "$R"; mkdir -p "$R"
python3 gen_jobs_v3.py "$R" "$TH" > "$R/jobs.txt" 2>"$R/jobs.count"
echo "V3 $(cat $R/jobs.count) (concurrency=$MAXJ x $TH = $((MAXJ*TH)) cores)"
while IFS= read -r line; do
  eval "$line" &
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do wait -n; done
done < "$R/jobs.txt"
wait
echo "V3 ALL_DONE"

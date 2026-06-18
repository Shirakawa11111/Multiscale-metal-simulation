#!/bin/bash
# Large-scale FCC interaction-matrix campaign, saturating the 128-core node.
# Many small jobs x few threads, big-first ordering (from gen_jobs.py).
set -u
cd ~/BO/interaction_matrix
export PYTHONPATH=~/BO/exadis_src/python
R=matrix; TH=3; MAXJ=40
rm -rf "$R"; mkdir -p "$R"
python3 gen_jobs.py "$R" "$TH" > "$R/jobs.txt" 2> "$R/jobs.count"
echo "MATRIX $(cat $R/jobs.count)  (concurrency=$MAXJ x $TH = $((MAXJ*TH)) cores)"
while IFS= read -r line; do
  eval "$line" &
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do wait -n; done
done < "$R/jobs.txt"
wait
echo "MATRIX ALL_DONE"

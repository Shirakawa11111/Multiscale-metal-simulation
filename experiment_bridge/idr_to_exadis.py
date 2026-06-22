"""IDR -> ExaDiS network (CLI + thin wrapper).

Lowers a defect_idr_v1 document to an ExaDiS manual-network JSON under a chosen assignment/cell policy.
The implementation lives in defect_ir.adapters.to_exadis; this module gives the roadmap-named entry point
and a CLI for sweeping policies (the BO/UQ hook).

  python3 experiment_bridge/idr_to_exadis.py results_exadis/cu_stem_idr.json \
      --assignment top1|sample --cell as_is|thickened_periodic --zbox 5 --seed 0 -o out.json
"""

import os, sys, json, argparse

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from defect_ir.adapters.to_exadis import idr_to_exadis_network  # re-exported


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("idr_json")
    ap.add_argument(
        "--assignment",
        default="top1",
        choices=["top1", "sample_linewise", "sample_edgewise"],
        help="top1=deterministic; sample_linewise=physical Monte-Carlo (default for UQ); "
        "sample_edgewise=artifact stress-test only",
    )
    ap.add_argument("--cell", default="as_is", choices=["as_is", "thickened_periodic"])
    ap.add_argument("--zbox", type=float, default=5.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()
    doc = json.load(open(a.idr_json))
    net = idr_to_exadis_network(
        doc,
        assignment_policy=a.assignment,
        cell_policy=a.cell,
        zbox=a.zbox,
        seed=a.seed,
    )
    out = a.out or a.idr_json.replace(".json", f"_exadis_{a.assignment}_{a.cell}.json")
    json.dump(net, open(out, "w"), indent=1)
    print(
        f"{net['network_counts']} periodic={net['cell']['is_periodic']} policy={a.assignment}/{a.cell} -> {out}"
    )


if __name__ == "__main__":
    main()

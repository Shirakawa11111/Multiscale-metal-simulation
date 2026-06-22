"""Synthetic g.b experiment [g.b-ready interface, LOCAL, no DDD].

Demonstrates how the IDR collapses assignment ambiguity when diffraction-contrast (g.b) data arrives.
For each reconstructed line we take its candidate set, pick ONE synthetic 'true' Burgers per parent line
(line-coherent, NOT per segment), generate g.b invisibility observations (visible <=> |g.b|>tol) under a
chosen set of reflections, share that truth across all the line's segments, and apply
`defect_ir.uncertainty.apply_gb_constraints`. Reports how assignment entropy collapses:
  no g.b  -> 3 candidates, ~1.58 bits
  1 g     -> partial, entropy drops
  2 g     -> typically resolved, entropy ~0
This is the path from a geometry-only IDR to a gb_validated one; plug in real reflections to use it.

  python3 experiment_bridge/synthetic_gb.py  -> results_exadis/synthetic_gb.{json,md}
"""

import os, sys, json, random
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
from defect_ir.uncertainty import apply_gb_constraints, assignment_entropy

OUT = os.path.join(os.path.dirname(__file__), "results_exadis")
IDR = os.path.join(OUT, "cu_stem_idr.json")
TOL = 0.1
G_SETS = {
    "no_gb": [],
    "1_reflection_g200": [[2, 0, 0]],
    "2_reflections_g200_g020": [[2, 0, 0], [0, 2, 0]],
    "3_reflections_g200_g020_g002": [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
}


def visible(g, b, tol=TOL):
    g = np.array(g, float)
    b = np.array(b, float)
    return abs(float(g @ b) / (np.linalg.norm(g) * np.linalg.norm(b))) > tol


def main():
    doc = json.load(open(IDR))
    labels = doc["topology"]["edge_labels"]
    rng = random.Random(0)
    # LINE-COHERENT synthetic ground truth: ONE random true Burgers per parent_line_id, shared by all its
    # edges (consistent with the IDR v1.1 line-coherent principle -- avoids re-introducing the edgewise artifact).
    eid_to_plid = {e["id"]: e.get("parent_line_id") for e in doc["geometry"]["edges"]}
    line_true = {}
    for lab in labels:
        pid = eid_to_plid.get(lab["edge_id"])
        if pid not in line_true:
            cands = lab["slip_system_candidates"]
            line_true[pid] = cands[rng.randrange(len(cands))]["b"]
    true_bs = [line_true[eid_to_plid.get(lab["edge_id"])] for lab in labels]
    rep = {
        "source": os.path.basename(IDR),
        "tol": TOL,
        "n_edges": len(labels),
        "scenarios": {},
    }
    for name, gs in G_SETS.items():
        ents, ncand, nresolved = [], [], 0
        statuses = {}
        for lab, true_b in zip(labels, true_bs):
            cands = lab["slip_system_candidates"]
            obs = [{"g": g, "visible": visible(g, true_b)} for g in gs]
            kept, status = apply_gb_constraints(cands, obs, tol=TOL)
            ents.append(assignment_entropy(kept))
            ncand.append(len(kept))
            nresolved += len(kept) == 1
            statuses[status] = statuses.get(status, 0) + 1
        rep["scenarios"][name] = {
            "n_reflections": len(gs),
            "mean_entropy_bits": round(float(np.mean(ents)), 4),
            "mean_candidates": round(float(np.mean(ncand)), 2),
            "frac_resolved": round(nresolved / len(labels), 3),
            "status_counts": statuses,
        }
    json.dump(rep, open(os.path.join(OUT, "synthetic_gb.json"), "w"), indent=1)
    rows = "\n".join(
        f"| {k} | {v['n_reflections']} | {v['mean_entropy_bits']} | {v['mean_candidates']} | {v['frac_resolved']} |"
        for k, v in rep["scenarios"].items()
    )
    md = f"""# Synthetic g·b ambiguity collapse (g·b-ready interface)

`{rep['source']}`, {rep['n_edges']} edges, invisibility tol={TOL}. Synthetic ground truth: ONE random true
Burgers per **parent line** (line-coherent, consistent with IDR v1.1), shared by all that line's segments.

| scenario | # reflections | mean entropy (bits) | mean candidates | frac resolved |
|--|--|--|--|--|
{rows}

**Reading.** With no g·b the assignment is ~3-way ambiguous (~1.58 bits). Each added reflection that obeys
the invisibility criterion (|g·b|≈0) removes incompatible candidates; ~2 well-chosen reflections collapse
most lines to a single slip system (`gb_validated`). This is the interface to plug in **real** diffraction
contrast: populate each line's `uncertainty.gb_constraints` with observed (g, visible) pairs and call
`apply_gb_constraints` to upgrade the IDR from `geometry_only_pending_gb` to `gb_validated`.
"""
    open(os.path.join(OUT, "synthetic_gb.md"), "w").write(md)
    for k, v in rep["scenarios"].items():
        print(
            f"{k:>32}: entropy={v['mean_entropy_bits']} bits  mean_cand={v['mean_candidates']}  resolved={v['frac_resolved']}"
        )
    print("-> results_exadis/synthetic_gb.{json,md}")


if __name__ == "__main__":
    main()

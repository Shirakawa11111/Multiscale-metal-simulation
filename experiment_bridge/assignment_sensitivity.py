"""Monte-Carlo assignment / cell-policy sensitivity audit (LOCAL, no DDD run).  [M3/M4 seed]

Quantifies how much the *lowered* ExaDiS network depends on the slip-system assignment policy and the
cell policy -- i.e. the structural uncertainty the IDR exposes, BEFORE spending any HPC. The DDD-outcome
part (relaxation/loading survival, topology events) is M3-on-HPC and is deferred.

  python3 experiment_bridge/assignment_sensitivity.py [cu_stem_idr.json] [N]
  -> results_exadis/assignment_sensitivity.{json,md}
"""

import os, sys, json, math
from collections import Counter

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
from defect_ir.adapters.to_exadis import idr_to_exadis_network

OUT = os.path.join(os.path.dirname(__file__), "results_exadis")
IDR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(OUT, "cu_stem_idr.json")
N = int(sys.argv[2]) if len(sys.argv) > 2 else 50


def seg_system_key(seg):
    # identify the chosen slip system of a segment by its (b, n) signature
    return tuple(round(x, 3) for x in seg[2:8])


def main():
    doc = json.load(open(IDR))
    n_edges = doc["topology"]["counts"]["n_edges"]

    # --- Monte-Carlo over assignment policy=sample ---
    per_seg = [
        Counter() for _ in range(n_edges)
    ]  # system signature distribution per segment
    inventories = []  # per-sample slip-system inventory
    for s in range(N):
        net = idr_to_exadis_network(
            doc, assignment_policy="sample", cell_policy="as_is", seed=s
        )
        inv = Counter()
        for i, seg in enumerate(net["segs"]):
            key = seg_system_key(seg)
            per_seg[i][key] += 1
            inv[key] += 1
        inventories.append(inv)

    # per-segment assignment entropy across MC samples (bits)
    def entropy(counter, total):
        return -sum(
            (c / total) * math.log2(c / total) for c in counter.values() if c > 0
        )

    seg_entropy = [entropy(c, N) for c in per_seg]
    mean_seg_entropy = sum(seg_entropy) / len(seg_entropy) if seg_entropy else 0.0
    frac_multi = (
        sum(1 for c in per_seg if len(c) > 1) / len(per_seg) if per_seg else 0.0
    )

    # slip-system inventory variability across samples (how unstable is the population?)
    all_keys = set().union(*[set(inv) for inv in inventories]) if inventories else set()
    inv_stats = {}
    for k in all_keys:
        vals = [inv.get(k, 0) for inv in inventories]
        mu = sum(vals) / len(vals)
        sd = (sum((v - mu) ** 2 for v in vals) / len(vals)) ** 0.5
        inv_stats[str(k)] = {
            "mean": round(mu, 1),
            "sd": round(sd, 1),
            "cv": round(sd / mu, 2) if mu > 0 else None,
        }

    # --- cell-policy contrast (deterministic) ---
    net_foil = idr_to_exadis_network(doc, cell_policy="as_is")
    cells = {}
    for zb in (3.0, 5.0, 10.0):
        netp = idr_to_exadis_network(doc, cell_policy="thickened_periodic", zbox=zb)
        cells[f"thickened_zbox{zb:g}"] = {
            "is_periodic": netp["cell"]["is_periodic"],
            "z_h": netp["cell"]["h_b"][2][2],
        }
    cells["as_is_foil"] = {
        "is_periodic": net_foil["cell"]["is_periodic"],
        "z_h": net_foil["cell"]["h_b"][2][2],
    }

    rep = {
        "n_samples": N,
        "n_edges": n_edges,
        "assignment_policy_sample": {
            "mean_segment_entropy_bits": round(mean_seg_entropy, 4),
            "frac_segments_multivalued": round(frac_multi, 4),
            "n_distinct_systems_used": len(all_keys),
            "inventory_mean_cv": round(
                sum(v["cv"] for v in inv_stats.values() if v["cv"] is not None)
                / max(1, sum(1 for v in inv_stats.values() if v["cv"] is not None)),
                3,
            ),
        },
        "cell_policy_contrast": cells,
        "interpretation": (
            "Because every segment's geometric assignment is ~3-way degenerate "
            "(confidence ~0.33), sampled lowerings reshuffle the slip-system inventory "
            "substantially -- so any downstream DDD metric must be reported as a "
            "distribution over assignment samples, not a single number. Cell policy "
            "(foil non-periodic vs thickened periodic) changes z-periodicity and box "
            "height, a separate, deterministic knob."
        ),
    }
    json.dump(
        rep, open(os.path.join(OUT, "assignment_sensitivity.json"), "w"), indent=1
    )
    a = rep["assignment_policy_sample"]
    md = (
        f"""# Assignment / cell-policy sensitivity audit (local, no DDD)

`{os.path.basename(IDR)}`, {N} Monte-Carlo assignment samples, {n_edges} segments.

## Assignment policy = sample
| metric | value |
|--|--|
| mean per-segment assignment entropy (bits) | {a['mean_segment_entropy_bits']} |
| fraction of segments that take >1 system across samples | {a['frac_segments_multivalued']} |
| distinct slip systems used | {a['n_distinct_systems_used']} |
| inventory coefficient-of-variation (mean) | {a['inventory_mean_cv']} |

## Cell policy contrast
| policy | is_periodic | z box (b) |
|--|--|--|
"""
        + "\n".join(
            f"| {k} | {v['is_periodic']} | {v['z_h']:.0f} |" for k, v in cells.items()
        )
        + f"""

**Interpretation.** {rep['interpretation']}
"""
    )
    open(os.path.join(OUT, "assignment_sensitivity.md"), "w").write(md)
    print(json.dumps(rep["assignment_policy_sample"]))
    print("cell policies:", json.dumps(cells))
    print("-> results_exadis/assignment_sensitivity.{json,md}")


if __name__ == "__main__":
    main()

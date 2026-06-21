"""Lower a Defect-IDR document to an ExaDiS manual-network dict (the stem_network.json shape).

This is the `idr_to_exadis` stage: it takes the uncertainty-aware IDR and commits to ONE concrete
network given an assignment policy + cell policy. Because the policy is an argument, BO/UQ can sweep
it (top-1 vs sampled slip system, foil vs thickened-periodic cell, endpoint policy) without touching
the upstream reconstruction.
"""
import random


def _pick_system(label, policy, rng):
    cands = label.get("slip_system_candidates", [])
    if not cands:
        return None
    if policy == "sample":
        r = rng.random(); acc = 0.0
        for c in cands:
            acc += c.get("prior", 0.0)
            if r <= acc:
                return c
        return cands[-1]
    # default top-1 (the chosen_system, else highest prior)
    cs = label.get("chosen_system")
    return next((c for c in cands if c["system_id"] == cs), max(cands, key=lambda c: c.get("prior", 0.0)))


def idr_to_exadis_network(doc, assignment_policy="top1", cell_policy="as_is",
                          zbox=1.0, seed=0):
    """Return an ExaDiS manual-network dict.

    assignment_policy: 'top1' (chosen_system) | 'sample' (Monte-Carlo from priors).
    cell_policy: 'as_is' (keep IDR cell + periodicity) |
                 'thickened_periodic' (scale h[2,2] by zbox, force is_periodic all True --
                  the hardening-pilot policy; see experiment_bridge/CELL_POLICY.md).
    """
    rng = random.Random(seed)
    g = doc["geometry"]; topo = doc["topology"]
    vid_to_idx = {v["id"]: i for i, v in enumerate(g["vertices"])}
    labels = {l["edge_id"]: l for l in topo.get("edge_labels", [])}

    nodes = []
    for v in g["vertices"]:
        pos = list(v["pos"])
        if len(pos) == 2:
            pos = pos + [0.0]
        con = "PINNED_NODE" if v.get("constraint") in ("pinned", "surface") else "UNCONSTRAINED"
        nodes.append([float(pos[0]), float(pos[1]), float(pos[2]), con])

    segs = []
    used_policy_counts = {"with_system": 0, "no_system": 0}
    for e in g["edges"]:
        i1, i2 = vid_to_idx[e["v1"]], vid_to_idx[e["v2"]]
        lab = labels.get(e["id"])
        sysc = _pick_system(lab, assignment_policy, rng) if lab else None
        if sysc is None:
            used_policy_counts["no_system"] += 1
            continue
        b, n = sysc["b"], sysc["n"]
        segs.append([i1, i2, float(b[0]), float(b[1]), float(b[2]),
                     float(n[0]), float(n[1]), float(n[2])])
        used_policy_counts["with_system"] += 1

    cell = dict(g["cell"])
    h = cell.get("h")
    is_per = list(cell.get("is_periodic", [True, True, False]))
    if cell_policy == "thickened_periodic" and h is not None:
        h = [list(row) for row in h]
        h[2][2] = h[2][2] * float(zbox)
        is_per = [True, True, True]

    return {
        "format": "exadis_python_manual_network_v0",
        "template_id": doc["provenance"].get("input_ref") or "idr_lowered",
        "template_type": "defect_idr_lowered",
        "units": doc["units"],
        "cell": {"h_b": h, "box_size_m": cell.get("box_size_m"), "is_periodic": is_per},
        "nodes": nodes,
        "segs": segs,
        "node_columns": ["x_b", "y_b", "z_b", "constraint"],
        "seg_columns": ["node1", "node2", "burg_x", "burg_y", "burg_z", "plane_x", "plane_y", "plane_z"],
        "network_counts": {"nodes": len(nodes), "segments": len(segs)},
        "provenance": {"lowered_from": "defect_idr_v1",
                       "assignment_policy": assignment_policy,
                       "cell_policy": cell_policy, "zbox": zbox, "seed": seed,
                       **used_policy_counts},
    }

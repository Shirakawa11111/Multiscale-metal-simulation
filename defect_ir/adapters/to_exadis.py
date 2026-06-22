"""Lower a Defect-IDR document to an ExaDiS manual-network dict (the stem_network.json shape).

This is the `idr_to_exadis` stage: it takes the uncertainty-aware IDR and commits to ONE concrete
network given an assignment policy + cell policy. Because the policy is an argument, BO/UQ can sweep
it (top-1 vs sampled slip system, foil vs thickened-periodic cell, endpoint policy) without touching
the upstream reconstruction.
"""

import random


def _sample(cands, rng):
    r = rng.random()
    acc = 0.0
    for c in cands:
        acc += c.get("prior", 0.0)
        if r <= acc:
            return c
    return cands[-1]


def _pick_system(label, policy, rng):
    """policy in {top1, sample/sample_edgewise, sample_linewise}.
    NB: sample_linewise is resolved upstream (one draw per parent line); if it reaches here it means the
    edge has no parent_line_id, so we fall back to an edgewise draw."""
    cands = label.get("slip_system_candidates", [])
    if not cands:
        return None
    if policy in ("sample", "sample_edgewise", "sample_linewise"):
        return _sample(cands, rng)
    cs = label.get("chosen_system")  # top1
    return next(
        (c for c in cands if c["system_id"] == cs),
        max(cands, key=lambda c: c.get("prior", 0.0)),
    )


def idr_to_exadis_network(
    doc,
    assignment_policy="top1",
    cell_policy="as_is",
    zbox=1.0,
    seed=0,
    endpoint_policy="pinned",
):
    """Return an ExaDiS manual-network dict.

    assignment_policy:
      'top1'            -> the chosen_system (deterministic baseline).
      'sample_linewise' -> PHYSICAL Monte-Carlo default: ONE draw per parent reconstructed line,
                           applied to all its segments (no within-line Burgers discontinuities).
      'sample_edgewise' -> per-edge draw; creates within-line discontinuities -> use ONLY as an
                           artifact stress-test / upper bound (see REAL_NETWORK_AUDIT.md v1.1).
      'sample'          -> DEPRECATED ambiguous alias; normalized to 'sample_linewise' with a warning.
    cell_policy: 'as_is' (keep IDR cell + periodicity) |
                 'thickened_periodic' (scale h[2,2] by zbox, force is_periodic all True --
                  the hardening-pilot policy; see experiment_bridge/CELL_POLICY.md).
    endpoint_policy: 'pinned' (honor the IDR vertex constraints) |
                 'free' (release pinned endpoints -> all UNCONSTRAINED, tests anchor sensitivity).
    """
    if assignment_policy == "sample":
        import warnings

        warnings.warn(
            "assignment_policy='sample' is ambiguous; using the physical default 'sample_linewise'. "
            "Pass 'sample_edgewise' explicitly only for the within-line-discontinuity stress-test.",
            stacklevel=2,
        )
        assignment_policy = "sample_linewise"
    rng = random.Random(seed)
    g = doc["geometry"]
    topo = doc["topology"]
    vid_to_idx = {v["id"]: i for i, v in enumerate(g["vertices"])}
    labels = {l["edge_id"]: l for l in topo.get("edge_labels", [])}

    nodes = []
    for v in g["vertices"]:
        pos = list(v["pos"])
        if len(pos) == 2:
            pos = pos + [0.0]
        pinned = endpoint_policy != "free" and v.get("constraint") in (
            "pinned",
            "surface",
        )
        con = "PINNED_NODE" if pinned else "UNCONSTRAINED"
        nodes.append([float(pos[0]), float(pos[1]), float(pos[2]), con])

    # LINE-COHERENT sampling: one draw per parent reconstructed line, applied to all its segments
    # (prevents adjacent segments of the SAME physical line getting different Burgers -> artificial junctions).
    line_choice = {}
    if assignment_policy == "sample_linewise":
        for e in g["edges"]:
            pid = e.get("parent_line_id")
            if pid is None:
                continue
            if pid not in line_choice:
                lab = labels.get(e["id"])
                cands = lab.get("slip_system_candidates", []) if lab else []
                line_choice[pid] = _sample(cands, rng) if cands else None

    segs = []
    used_policy_counts = {"with_system": 0, "no_system": 0}
    for e in g["edges"]:
        i1, i2 = vid_to_idx[e["v1"]], vid_to_idx[e["v2"]]
        lab = labels.get(e["id"])
        pid = e.get("parent_line_id")
        if assignment_policy == "sample_linewise" and pid in line_choice:
            sysc = line_choice[pid]  # coherent across the whole line
        else:
            sysc = _pick_system(lab, assignment_policy, rng) if lab else None
        if sysc is None:
            used_policy_counts["no_system"] += 1
            continue
        b, n = sysc["b"], sysc["n"]
        segs.append(
            [
                i1,
                i2,
                float(b[0]),
                float(b[1]),
                float(b[2]),
                float(n[0]),
                float(n[1]),
                float(n[2]),
            ]
        )
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
        "seg_columns": [
            "node1",
            "node2",
            "burg_x",
            "burg_y",
            "burg_z",
            "plane_x",
            "plane_y",
            "plane_z",
        ],
        "network_counts": {"nodes": len(nodes), "segments": len(segs)},
        "provenance": {
            "lowered_from": "defect_idr_v1",
            "assignment_policy": assignment_policy,
            "cell_policy": cell_policy,
            "zbox": zbox,
            "seed": seed,
            **used_policy_counts,
        },
    }

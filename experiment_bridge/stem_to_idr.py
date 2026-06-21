"""STEM reconstruction -> Defect-IR (uncertainty-aware).  [M2]

Replaces the single-assignment stage of `stem_to_exadis.py`: instead of forcing one slip system per
reconstructed line via argmin|n.t|, this emits a full `defect_idr_v1` document where every line carries
the **top-k slip-system candidates with confidence**, plus explicit z-depth / Burgers-source /
endpoint-policy uncertainty. Downstream, `defect_ir.adapters.to_exadis.idr_to_exadis_network` lowers it
to an ExaDiS network under a chosen assignment_policy + cell_policy (the BO/UQ hook).

  python3 experiment_bridge/stem_to_idr.py          # -> results_exadis/cu_stem_idr.json (+ report)
"""

import os, sys, glob, json
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
from defect_ir import schema as S
from defect_ir.validators import validate_idr
from defect_ir.uncertainty import slip_system_candidates, network_assignment_summary
from defect_ir.adapters.to_exadis import idr_to_exadis_network

RECON = os.path.join(os.path.dirname(__file__), "recon_data")
OUT = os.path.join(os.path.dirname(__file__), "results_exadis")
B_CU = 2.556e-10
LXY_B, LZ_B = 8000.0, 600.0


def fcc_catalog():
    B = [
        [0, 1, -1],
        [1, 0, -1],
        [1, -1, 0],
        [0, 1, -1],
        [1, 0, 1],
        [1, 1, 0],
        [0, 1, 1],
        [1, 0, -1],
        [1, 1, 0],
        [0, 1, 1],
        [1, 0, 1],
        [1, -1, 0],
    ]
    N = [[1, 1, 1]] * 3 + [[-1, 1, 1]] * 3 + [[1, -1, 1]] * 3 + [[1, 1, -1]] * 3
    return [{"system_id": i, "b": B[i], "n": N[i]} for i in range(12)]


def load_lines():
    lines = []
    for fp in sorted(glob.glob(os.path.join(RECON, "points_3d*.txt"))):
        pts = [
            [float(x) for x in ln.split()] for ln in open(fp) if len(ln.split()) == 3
        ]
        if len(pts) >= 2:
            lines.append(np.array(pts))
    return lines


def build_idr():
    cat = fcc_catalog()
    lines = load_lines()
    allp = np.vstack(lines)
    lo, hi = allp.min(0), allp.max(0)
    sp_xy = max(hi[0] - lo[0], hi[1] - lo[1], 1e-9)
    sp_z = max(hi[2] - lo[2], 1e-9)

    def to_b(p):  # identical mapping to stem_to_exadis (fidelity)
        return [
            0.05 * LXY_B + 0.9 * LXY_B * (p[0] - lo[0]) / sp_xy,
            0.05 * LXY_B + 0.9 * LXY_B * (p[1] - lo[1]) / sp_xy,
            0.1 * LZ_B + 0.8 * LZ_B * (p[2] - lo[2]) / sp_z,
        ]

    doc = S.empty_idr("fcc_cu", 3, "stem_3d_reconstruction", "b", B_CU)
    doc["provenance"].update(
        tool="stem_to_idr",
        input_ref="experiment_bridge/recon_data/points_3d*.txt",
        notes="Burgers geometry-only (plane containment); g.b not done; z weakly constrained (stereo).",
    )
    doc["geometry"]["cell"] = {
        "h": [[LXY_B, 0, 0], [0, LXY_B, 0], [0, 0, LZ_B]],
        "box_size_m": [LXY_B * B_CU, LXY_B * B_CU, LZ_B * B_CU],
        "is_periodic": [True, True, False],
    }
    doc["topology"]["slip_system_catalog"] = cat

    vid = 0
    edges = []
    edge_labels = []
    for ln in lines:
        bln = [to_b(p) for p in ln]
        first = vid
        for j, p in enumerate(bln):
            endp = j in (0, len(bln) - 1)
            doc["geometry"]["vertices"].append(
                {
                    "id": vid,
                    "pos": [float(p[0]), float(p[1]), float(p[2])],
                    "role": "endpoint" if endp else "interior",
                    "constraint": "pinned" if endp else "free",
                }
            )
            vid += 1
        t = [bln[-1][k] - bln[0][k] for k in range(3)]
        cands = slip_system_candidates(t, cat, T=0.15, topk=3)
        for j in range(len(bln) - 1):
            eid = len(edges)
            edges.append(
                {
                    "id": eid,
                    "v1": first + j,
                    "v2": first + j + 1,
                    "kind": "dislocation_segment",
                }
            )
            edge_labels.append(
                S.edge_label(eid, cands, assignment_status="geometry_only_pending_gb")
            )
    doc["geometry"]["edges"] = edges
    doc["topology"]["edge_labels"] = edge_labels
    doc["topology"]["counts"] = {
        "n_vertices": vid,
        "n_edges": len(edges),
        "n_lines": len(lines),
    }
    summ = network_assignment_summary(edge_labels)
    doc["uncertainty"] = {
        "z_depth": {
            "weakly_constrained": True,
            "sigma_nm": 30,
            "note": "STEM stereo weakly constrains z; mapped to foil thickness ~150 nm",
        },
        "burgers_assignment": {
            "method": "geometric_plane_containment",
            "validated_by_gb": False,
        },
        "endpoint_policy": {"policy": "pinned_due_to_truncated_reconstruction"},
        "system_size_caveat": f"{len(lines)} reconstructed lines -> limited hardening statistics",
        "assignment_summary": summ,
    }
    doc["simulation_targets"] = {
        "engine": "exadis",
        "cell_policy": "foil_nonperiodic_z",
        "loading_mode": "strain_rate",
        "edir": [0, 0, 1],
        "erate": 1e4,
        "force_model": "SUBCYCLING_MODEL",
        "mobility": "FCC_0",
        "note": "hardening pilot may switch to thickened_periodic; see CELL_POLICY.md",
    }
    return doc, summ


def write_report(doc, summ, ok, errs, warns):
    """Human + machine readable audit report (the auditability deliverable)."""
    c = doc["topology"]["counts"]
    # top slip-system histogram (chosen)
    hist = {}
    for lab in doc["topology"]["edge_labels"]:
        hist[lab.get("chosen_system")] = hist.get(lab.get("chosen_system"), 0) + 1
    rep = {
        "source": doc["provenance"]["input_ref"],
        "valid": ok,
        "n_errors": len(errs),
        "n_warnings": len(warns),
        "n_lines": c.get("n_lines"),
        "n_vertices": c.get("n_vertices"),
        "n_edges": c.get("n_edges"),
        "assignment": summ,
        "burgers_source": "geometry_only_pending_gb",
        "cell_policy": doc["simulation_targets"].get("cell_policy"),
        "endpoint_policy": doc["uncertainty"]["endpoint_policy"]["policy"],
        "z_uncertainty_nm": doc["uncertainty"]["z_depth"]["sigma_nm"],
        "chosen_system_histogram": {
            str(k): v for k, v in sorted(hist.items(), key=lambda x: -x[1])
        },
    }
    json.dump(rep, open(os.path.join(OUT, "cu_stem_idr_report.json"), "w"), indent=1)
    md = f"""# STEM -> IDR audit report (Cu)

Source: `{rep['source']}`  |  valid={ok}  errors={len(errs)} warnings={len(warns)}

| quantity | value |
|--|--|
| reconstructed lines | {rep['n_lines']} |
| vertices / edges | {rep['n_vertices']} / {rep['n_edges']} |
| Burgers source | {rep['burgers_source']} |
| endpoint policy | {rep['endpoint_policy']} |
| z uncertainty | ~{rep['z_uncertainty_nm']} nm (stereo-weak) |
| cell policy | {rep['cell_policy']} |

## Assignment uncertainty (the key audit result)
| metric | value |
|--|--|
| mean confidence | {summ['mean_confidence']} |
| min confidence | {summ['min_confidence']} |
| frac low-confidence (<0.5) | {summ['frac_low_confidence(<0.5)']} |
| mean assignment entropy (bits) | {summ['mean_assignment_entropy_bits']} |
| ambiguous edges (entropy>0.8) | {summ['n_ambiguous(entropy>0.8)']} / {rep['n_edges']} |

**Highlight.** The IDR did not merely reformat the network — it **exposed** that the STEM->DDD
slip-system assignment is intrinsically ambiguous without g·b: geometry fixes the {{111}} *plane*, but
the three ⟨110⟩ Burgers on it are near-degenerate (mean confidence ~{summ['mean_confidence']},
~{summ['mean_assignment_entropy_bits']} bits), so essentially every segment is assignment-ambiguous.
The legacy `stem_to_exadis.py` forced one of these and hid it. A real g·b analysis would set
`assignment_status = gb_validated` and collapse the priors.
"""
    open(os.path.join(OUT, "cu_stem_idr_report.md"), "w").write(md)


def main():
    os.makedirs(OUT, exist_ok=True)
    doc, summ = build_idr()
    ok, errs, warns = validate_idr(doc)
    json.dump(doc, open(os.path.join(OUT, "cu_stem_idr.json"), "w"), indent=1)
    # also lower to an ExaDiS network (top-1 policy) so the old pipeline still gets its input
    net = idr_to_exadis_network(doc, assignment_policy="top1", cell_policy="as_is")
    json.dump(net, open(os.path.join(OUT, "stem_network_from_idr.json"), "w"), indent=1)
    write_report(doc, summ, ok, errs, warns)
    print(
        f"cu_stem_idr.json: valid={ok} errs={len(errs)} warns={len(warns)} "
        f"verts={doc['topology']['counts']['n_vertices']} edges={doc['topology']['counts']['n_edges']} "
        f"lines={doc['topology']['counts']['n_lines']}"
    )
    print("assignment uncertainty:", json.dumps(summ))
    print(f"lowered -> stem_network_from_idr.json: {net['network_counts']}")
    print("report -> cu_stem_idr_report.{json,md}")
    for e in errs:
        print("  ERROR:", e)


if __name__ == "__main__":
    main()

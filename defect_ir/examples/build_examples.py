"""Build + validate the two worked IDR examples (Cu STEM dislocation graph, graphene defect graph).
Doubles as the smoke test for the whole defect_ir package. Run:  python3 -m defect_ir.examples.build_examples
"""
import json, os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from defect_ir import schema as S
from defect_ir.validators import validate_idr
from defect_ir.uncertainty import slip_system_candidates, network_assignment_summary, assignment_entropy
from defect_ir.adapters.to_exadis import idr_to_exadis_network

HERE = os.path.dirname(__file__)


def fcc_catalog():
    # 12 FCC systems: 4 {111} planes x 3 <110> in-plane Burgers (matches experiment_bridge convention)
    B = [[0, 1, -1], [1, 0, -1], [1, -1, 0], [0, 1, -1], [1, 0, 1], [1, 1, 0],
         [0, 1, 1], [1, 0, -1], [1, 1, 0], [0, 1, 1], [1, 0, 1], [1, -1, 0]]
    N = [[1, 1, 1]] * 3 + [[-1, 1, 1]] * 3 + [[1, -1, 1]] * 3 + [[1, 1, -1]] * 3
    return [{"system_id": i, "b": B[i], "n": N[i]} for i in range(12)]


# ---------------- Example 1: Cu STEM dislocation graph (3D) ----------------
def build_cu_stem():
    cat = fcc_catalog()
    doc = S.empty_idr("fcc_cu", 3, "stem_3d_reconstruction", "b", 2.556e-10)
    doc["provenance"].update(tool="stem_to_idr", input_ref="recon_data/points_3d*.txt",
                             notes="Burgers geometry-only (plane containment); g.b not done. z weakly constrained.")
    LXY, LZ = 8000.0, 600.0
    doc["geometry"]["cell"] = {"h": [[LXY, 0, 0], [0, LXY, 0], [0, 0, LZ]],
                               "box_size_m": [LXY * 2.556e-10, LXY * 2.556e-10, LZ * 2.556e-10],
                               "is_periodic": [True, True, False]}   # foil: z non-periodic
    doc["topology"]["slip_system_catalog"] = cat
    # two short reconstructed-like polylines (3 nodes each)
    polylines = [
        [[1000, 1000, 200], [2200, 1500, 230], [3400, 2000, 260]],
        [[5000, 4000, 300], [5300, 5200, 320], [5600, 6400, 350]],
    ]
    vid = 0; edges = []; edge_labels = []
    for li, pl in enumerate(polylines):
        first = vid
        for j, p in enumerate(pl):
            role = "endpoint" if j in (0, len(pl) - 1) else "interior"
            con = "pinned" if j in (0, len(pl) - 1) else "free"
            doc["geometry"]["vertices"].append({"id": vid, "pos": [float(p[0]), float(p[1]), float(p[2])],
                                                 "role": role, "constraint": con})
            vid += 1
        # tangent of the line -> slip-system candidates (keep top-3 with confidence)
        t = [pl[-1][k] - pl[0][k] for k in range(3)]
        cands = slip_system_candidates(t, cat, T=0.15, topk=3)
        for j in range(len(pl) - 1):
            eid = len(edges)
            edges.append({"id": eid, "v1": first + j, "v2": first + j + 1, "kind": "dislocation_segment"})
            edge_labels.append(S.edge_label(eid, cands, assignment_status="geometry_only_pending_gb"))
    doc["geometry"]["edges"] = edges
    doc["topology"]["edge_labels"] = edge_labels
    doc["topology"]["counts"] = {"n_vertices": vid, "n_edges": len(edges)}
    doc["uncertainty"] = {
        "z_depth": {"weakly_constrained": True, "sigma_nm": 30,
                    "note": "STEM stereo weakly constrains z; mapped heuristically to foil thickness"},
        "burgers_assignment": {"method": "geometric_plane_containment", "validated_by_gb": False},
        "endpoint_policy": {"policy": "pinned_due_to_truncated_reconstruction"},
        "system_size_caveat": "few reconstructed lines -> limited hardening statistics",
        "assignment_summary": network_assignment_summary(edge_labels),
    }
    doc["simulation_targets"] = {"engine": "exadis", "cell_policy": "foil_nonperiodic_z",
                                 "loading_mode": "strain_rate", "edir": [0, 0, 1], "erate": 1e4,
                                 "force_model": "SUBCYCLING_MODEL", "mobility": "FCC_0",
                                 "note": "hardening pilot may use thickened_periodic cell; see CELL_POLICY.md"}
    return doc


# ---------------- Example 2: graphene defect graph (2D, 5|7 cores) ----------------
def build_graphene():
    doc = S.empty_idr("graphene_2d", 2, "pfc_defect_detection", "a0", 2.46e-10)
    doc["provenance"].update(tool="defect_analysis", input_ref="results/*.npz",
                             notes="5|7 coordination cores from periodic Delaunay; PFC density field.")
    L = 100.0
    doc["geometry"]["cell"] = {"h": [[L, 0, 0], [0, L, 0], [0, 0, 0]],
                               "box_size_m": [L * 2.46e-10, L * 2.46e-10, 0.0],
                               "is_periodic": [True, True, False]}
    # one 5-7 dislocation core = a 5-coordinated atom paired with a 7-coordinated atom
    atoms = [{"id": 0, "pos": [48.0, 50.0], "role": "atom", "coordination": 5},
             {"id": 1, "pos": [52.0, 50.0], "role": "atom", "coordination": 7},
             {"id": 2, "pos": [50.0, 53.0], "role": "atom", "coordination": 6},
             {"id": 3, "pos": [50.0, 47.0], "role": "atom", "coordination": 6}]
    for a in atoms:
        doc["geometry"]["vertices"].append({"id": a["id"], "pos": a["pos"], "role": "atom",
                                            "constraint": "free"})
    doc["geometry"]["edges"] = [{"id": 0, "v1": 0, "v2": 1, "kind": "core_pair"}]
    doc["topology"]["vertex_labels"] = [{"vertex_id": a["id"], "coordination": a["coordination"],
                                         "defect_type": ("5" if a["coordination"] == 5 else
                                                         "7" if a["coordination"] == 7 else "none")}
                                        for a in atoms]
    doc["topology"]["edge_labels"] = [{"edge_id": 0, "kind": "5-7_core",
                                       "slip_system_candidates": [], "chosen_system": None,
                                       "assignment_confidence": None,
                                       "assignment_status": "topology_inferred"}]
    doc["topology"]["counts"] = {"n_vertices": len(atoms), "n_edges": 1,
                                 "n_cores": 1, "density_m2": 1.0 / ((L * 2.46e-10) ** 2)}
    doc["uncertainty"] = {"burgers_assignment": {"method": "topology_5_7_core", "validated_by_gb": False},
                          "note": "2D Burgers implied by 5-7 separation; not a 3D vector"}
    doc["simulation_targets"] = {"engine": "pfc", "loading_mode": "affine_tension",
                                 "note": "PFC density-field seeding"}
    return doc


def main():
    out = {}
    for name, doc in [("cu_stem_idr", build_cu_stem()), ("graphene_defect_idr", build_graphene())]:
        ok, errs, warns = validate_idr(doc)
        path = os.path.join(HERE, name + ".json")
        json.dump(doc, open(path, "w"), indent=1)
        print(f"[{name}] valid={ok}  errors={len(errs)} warnings={len(warns)} -> {os.path.basename(path)}")
        for e in errs:
            print("   ERROR:", e)
        for w in warns:
            print("   warn:", w)
        out[name] = doc
    # demo: lower the Cu IDR to ExaDiS under two assignment policies (UQ hook)
    cu = out["cu_stem_idr"]
    net_top1 = idr_to_exadis_network(cu, assignment_policy="top1", cell_policy="as_is")
    net_smpl = idr_to_exadis_network(cu, assignment_policy="sample", cell_policy="thickened_periodic", zbox=5.0, seed=1)
    print(f"\nIDR->ExaDiS lowering (Cu): top1 segs={net_top1['network_counts']['segments']} "
          f"(periodic {net_top1['cell']['is_periodic']}); "
          f"sampled+thickened segs={net_smpl['network_counts']['segments']} "
          f"(periodic {net_smpl['cell']['is_periodic']})")
    print("assignment uncertainty (Cu):", json.dumps(cu["uncertainty"]["assignment_summary"]))


if __name__ == "__main__":
    main()

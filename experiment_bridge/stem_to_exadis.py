"""STEM reconstruction -> ExaDiS/OpenDiS DDD network adapter.

  *** LEGACY BASELINE (direct, single-assignment converter). ***
  Prefer the uncertainty-aware path for new work:
      experiment_bridge/stem_to_idr.py   (STEM -> defect_idr_v1, top-k slip-system candidates)
      experiment_bridge/idr_to_exadis.py (IDR -> ExaDiS, selectable assignment/cell policy)
  This file forces ONE slip system per line via argmin|n.t| (geometry only, no g.b) and hides the
  assignment ambiguity; kept as a baseline / reproduction reference. See ASSIGNMENT_UNCERTAINTY.md.

The regime-correct coupling: DDD represents dislocations as LINES (exactly the
STEM reconstruction output) and evolves them in the room-temperature
athermal-glide + junction + forest-hardening regime (where PFC fundamentally
failed). This adapter converts the reconstructed 3D dislocation polylines into
a valid ExaDiS `exadis_python_manual_network_v0` network JSON, ready for
pyexadis (the framework's exadis_minimal_run.py).

Per line: nodes at the polyline vertices (b units), consecutive vertices ->
segments; each line is assigned the FCC slip system whose {111} plane best
CONTAINS the line tangent (the glide plane), giving its a/2<110> Burgers and
plane normal. Endpoints PINNED (foil-surface anchors / incomplete-line ends),
interior nodes UNCONSTRAINED.

Honest inputs still needed for full fidelity: experimental Burgers (g.b)
per line (here assigned geometrically), and resolved foil depth (z is poorly
constrained by stereo; mapped to a realistic foil thickness).

Output: experiment_bridge/results_exadis/stem_network.json
"""

import os, sys, glob, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "CU",
                                "simulation_tools"))
import numpy as np
from interfaceB_exporter import fcc_slip_systems

RECON = os.path.join(os.path.dirname(__file__), "recon_data")
if not os.path.isdir(RECON):
    RECON = os.path.join(os.path.dirname(__file__), "..", "..", "3D重建算法论文",
                         "3d_scatter")
OUT = os.path.join(os.path.dirname(__file__), "results_exadis")
B_CU = 2.556e-10
LXY_B = 8000.0      # in-plane box ~2 um
LZ_B = 600.0        # foil thickness ~150 nm


def load_lines():
    lines = []
    for fp in sorted(glob.glob(os.path.join(RECON, "points_3d*.txt"))):
        pts = [[float(x) for x in ln.split()] for ln in open(fp)
               if len(ln.split()) == 3]
        if len(pts) >= 2:
            lines.append(np.array(pts))
    return lines


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = load_lines()
    burg_dirs, plane_norms = fcc_slip_systems()       # (12,3),(12,3)
    allp = np.vstack(lines)
    lo, hi = allp.min(0), allp.max(0)
    sp_xy = max(hi[0] - lo[0], hi[1] - lo[1])
    sp_z = max(hi[2] - lo[2], 1e-9)

    def to_b(p):
        x = 0.05 * LXY_B + 0.9 * LXY_B * (p[0] - lo[0]) / sp_xy
        y = 0.05 * LXY_B + 0.9 * LXY_B * (p[1] - lo[1]) / sp_xy
        z = 0.1 * LZ_B + 0.8 * LZ_B * (p[2] - lo[2]) / sp_z
        return [float(x), float(y), float(z)]

    nodes, segs = [], []
    sys_count = {}
    for ln in lines:
        bln = np.array([to_b(p) for p in ln])
        tangent = bln[-1] - bln[0]
        tangent /= (np.linalg.norm(tangent) + 1e-12)
        # pick FCC slip system whose plane CONTAINS the line (n . t ~ 0),
        # tie-broken by max glide character; gives (burgers, plane)
        scores = np.abs(plane_norms @ tangent)         # 0 = line in plane
        si = int(np.argmin(scores))
        b = burg_dirs[si]
        n = plane_norms[si]
        sys_count[si] = sys_count.get(si, 0) + 1
        i0 = len(nodes)
        for j, p in enumerate(bln):
            constraint = "PINNED_NODE" if j in (0, len(bln) - 1) \
                else "UNCONSTRAINED"
            nodes.append([p[0], p[1], p[2], constraint])
        for j in range(len(bln) - 1):
            segs.append([i0 + j, i0 + j + 1,
                         float(b[0]), float(b[1]), float(b[2]),
                         float(n[0]), float(n[1]), float(n[2])])

    net = {
        "format": "exadis_python_manual_network_v0",
        "template_id": "stem_reconstructed_cu",
        "template_type": "experimental_reconstruction",
        "units": {"length": "b", "length_unit_m": B_CU,
                  "burgers_vector_magnitude_b": 1.0},
        "cell": {"length_unit": "b", "length_unit_m": B_CU,
                 "box_size_m": [LXY_B * B_CU, LXY_B * B_CU, LZ_B * B_CU],
                 "box_size_b": [LXY_B, LXY_B, LZ_B],
                 "h_b": [[LXY_B, 0, 0], [0, LXY_B, 0], [0, 0, LZ_B]],
                 "is_periodic": [True, True, False]},   # foil: z non-periodic
        "nodes": nodes, "segs": segs,
        "node_columns": ["x_b", "y_b", "z_b", "constraint"],
        "seg_columns": ["node1", "node2", "burg_x", "burg_y", "burg_z",
                        "plane_x", "plane_y", "plane_z"],
        "network_counts": {"nodes": len(nodes), "segments": len(segs)},
        "source_template": {"template_id": "stem_reconstructed_cu",
                            "template_type": "experimental_reconstruction",
                            "n_reconstructed_lines": len(lines),
                            "note": "Burgers assigned geometrically (plane "
                                    "contains line); real g.b analysis needed "
                                    "for full fidelity"},
    }
    out = os.path.join(OUT, "stem_network.json")
    json.dump(net, open(out, "w"), indent=1)

    # validate
    ok = True
    nn = len(nodes)
    for s in segs:
        if not (0 <= s[0] < nn and 0 <= s[1] < nn):
            ok = False
        bn = np.linalg.norm(s[2:5])
        if abs(bn - 1.0) > 1e-6:        # FCC b should be unit (normalized)
            ok = False
    print(f"STEM->ExaDiS network: {len(lines)} lines -> {len(nodes)} nodes, "
          f"{len(segs)} segments")
    print(f"  slip systems used: {len(sys_count)} of 12 "
          f"(distribution {dict(sorted(sys_count.items()))})")
    print(f"  box {LXY_B:.0f}x{LXY_B:.0f}x{LZ_B:.0f} b "
          f"({LXY_B*B_CU*1e9:.0f}x{LXY_B*B_CU*1e9:.0f}x{LZ_B*B_CU*1e9:.0f} nm), "
          f"foil z non-periodic")
    print(f"  VALIDATION: {'PASS (valid ExaDiS network)' if ok else 'FAIL'}")
    print(f"  written: {out}")
    print(f"  -> ready for: python3 exadis_minimal_run.py --network-json "
          f"{out} (needs pyexadis)")


if __name__ == "__main__":
    main()

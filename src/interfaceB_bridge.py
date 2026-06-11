"""Interface B (CU framework) -> PFC bridge.

Reads the DAMASK ROI exports in CU/outputs/interfaceB_main/*.h5 and computes
the scale mapping between DAMASK dislocation-density state variables and a
PFC simulation box:

  - PFC length unit: triangular lattice constant a0 = 4*pi/sqrt(3) maps to
    the Cu Burgers vector b = 2.556e-10 m (frozen in 单晶铜拉伸模拟
    docs/units_mapping.md), so 1 PFC unit = b/a0 = 3.523e-11 m.
  - 2D dislocation areal density rho_2D [m^-2] is compared directly with the
    DAMASK line density rho [m^-2] (line length per volume == intersection
    count per area for straight parallel lines).

Outputs per ROI: density at t*, sigma history range, the physical size of
our PFC boxes, how many dislocations the ROI density implies in such a box,
and the minimum box needed to host >= 2 dislocations at ROI density.
Writes results/interfaceB_bridge/report.json.
"""

import os, glob, json
import numpy as np
import h5py

B_CU = 2.556e-10                      # m, Cu Burgers vector
A0 = 4 * np.pi / np.sqrt(3.0)         # PFC lattice constant (dimensionless)
SCALE = B_CU / A0                     # m per PFC length unit
IB_GLOB = "/Users/bojingkai/Desktop/CU/outputs/interfaceB_main/*.h5"
OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   "interfaceB_bridge")


def analyze(path):
    with h5py.File(path) as f:
        roi_id = f["roi/roi_id"][()].decode()
        roi_type = f["roi/roi_type"][()].decode()
        tstar_inc = int(f["roi/observed_tstar_increment"][()])
        incs = f["driving/increments"][:]
        i_t = int(np.argmin(np.abs(incs - tstar_inc)))
        rho_avg = float(f["driving/rho_total_avg"][i_t])
        rho_p95 = float(f["driving/rho_total_p95"][i_t])
        sig = f["driving/sigma_avg"][:]          # (T,3,3) Pa
        svm = np.sqrt(1.5 * np.einsum("tij,tij->t",
                                      sig - sig.mean(axis=(1, 2),
                                                     keepdims=True) *
                                      np.eye(3), sig))
        rho_cell_max = float(f["cp/rho_total_cp"][i_t].max())
        ddd_range = f["ddd/target_density_range_m2"][:].tolist()
        rho_seed_rule = f["mapping/rho_seed_rule"][()].decode()

    boxes = {}
    for npix, dx in (("512^2", np.pi / 4), ("2048^2", np.pi / 4)):
        n = int(npix.split("^")[0])
        L_phys = n * dx * SCALE
        boxes[npix] = dict(
            L_nm=L_phys * 1e9,
            n_dislocations_at_rho_avg=rho_avg * L_phys ** 2,
            n_dislocations_at_rho_cell_max=rho_cell_max * L_phys ** 2,
            rho_of_2_cores_m2=2.0 / L_phys ** 2,
        )
    L_min_2 = np.sqrt(2.0 / rho_avg)
    return dict(
        file=os.path.basename(path), roi_id=roi_id, roi_type=roi_type,
        tstar_increment=tstar_inc,
        rho_avg_m2=rho_avg, rho_p95_m2=rho_p95, rho_cell_max_m2=rho_cell_max,
        sigma_vm_max_MPa=float(svm.max() / 1e6),
        ddd_target_density_range_m2=ddd_range,
        rho_seed_rule=rho_seed_rule,
        pfc_boxes=boxes,
        min_box_for_2_cores_at_rho_avg_um=L_min_2 * 1e6,
        min_grid_at_dx_pi4=int(np.ceil(L_min_2 / SCALE / (np.pi / 4))),
    )


def main():
    os.makedirs(OUT, exist_ok=True)
    files = sorted(glob.glob(IB_GLOB))
    report = dict(scale_m_per_pfc_unit=SCALE, a0_pfc=A0, b_cu_m=B_CU,
                  rois=[analyze(p) for p in files])
    with open(os.path.join(OUT, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    for r in report["rois"]:
        b = r["pfc_boxes"]["512^2"]
        print(f"ROI {r['roi_id']} ({r['roi_type']}): rho_avg={r['rho_avg_m2']:.2e} "
              f"rho_cell_max={r['rho_cell_max_m2']:.2e} m^-2, "
              f"sigma_vm_max={r['sigma_vm_max_MPa']:.1f} MPa")
        print(f"   512^2 PFC box = {b['L_nm']:.1f} nm -> "
              f"{b['n_dislocations_at_rho_avg']:.3f} cores at rho_avg; "
              f"2 cores in this box ~ rho={b['rho_of_2_cores_m2']:.2e} m^-2")
        print(f"   box for >=2 cores at rho_avg: "
              f"{r['min_box_for_2_cores_at_rho_avg_um']:.2f} um "
              f"(grid {r['min_grid_at_dx_pi4']}^2)")


if __name__ == "__main__":
    main()

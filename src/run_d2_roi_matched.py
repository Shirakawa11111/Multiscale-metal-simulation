"""D2 flagship: ROI density-matched PFC tension run.

Physical correspondence (see interfaceB_bridge.py):
  1536^2 grid @ dx=pi/4  ->  L = 42.5 nm (a0 <-> b_Cu = 2.556 A)
  4 seeded edge cores (quadrupole) -> rho_2D = 2.2e15 m^-2
  ~ DAMASK hotspot ROI cell-max density (2.05e15 m^-2,
    state_hotspot_80000_component_2, increment 80000)

Loading: area-conserving tension, 0.25%/step to 10%, 400 relax/step.
Outputs DDD-event-proxy-style metrics: core count, density trace,
annihilation/nucleation events with field snapshots.
Output dir: results/d2_roi_matched_1536/
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   "d2_roi_matched_1536")

N = 1536
DEPS = 0.0025
N_STRAIN = 40
RELAX = 400
DT = 0.5
QUAD = [(0.33, 0.25, +1), (0.33, 0.75, -1),
        (0.67, 0.25, -1), (0.67, 0.75, +1)]
ROI_REF = dict(
    roi_id="state_hotspot_80000_component_2",
    rho_target_m2=2.05e15,
    scale_m_per_unit=2.556e-10 / (4 * np.pi / np.sqrt(3.0)),
)


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    m = PFC2D(N, N, r=-0.25, psi_bar=-0.25)
    m.init_dislocations(QUAD)
    m.step(DT, n=600)

    L_m = m.lx * ROI_REF["scale_m_per_unit"]
    rho_sim = 4.0 / L_m ** 2
    print(f"box {L_m*1e9:.1f} nm, rho_sim={rho_sim:.2e} m^-2 "
          f"(target {ROI_REF['rho_target_m2']:.2e})", flush=True)

    pts = find_peaks(m.psi, m.dx, m.dy)
    d = find_dislocations(pts, m.lx, m.ly)
    m.save(os.path.join(OUT, "initial.npz"))
    rows = [dict(exx=0.0, F=m.free_energy(), sigma=m.stress(),
                 cores=len(d["cores"]), rho=d["rho"],
                 core_xy=d["cores"].tolist())]
    print(f"eps=0 cores={len(d['cores'])}", flush=True)

    n_events = 0
    for i in range(N_STRAIN):
        m.apply_strain(DEPS, area_conserving=True)
        m.step(DT, n=RELAX)
        pts = find_peaks(m.psi, m.dx, m.dy)
        d = find_dislocations(pts, m.lx, m.ly)
        nc = len(d["cores"])
        rows.append(dict(exx=m.exx, F=m.free_energy(), sigma=m.stress(),
                         cores=nc, rho=d["rho"],
                         core_xy=d["cores"].tolist()))
        print(f"exx={m.exx*100:5.2f}% sigma={rows[-1]['sigma']:+.5f} "
              f"cores={nc}", flush=True)
        if nc != rows[-2]["cores"] and n_events < 8:
            m.save(os.path.join(OUT, f"event_{m.exx*100:.2f}pct.npz"))
            n_events += 1
        if i % 10 == 9:
            m.save(os.path.join(OUT, f"snap_{m.exx*100:.2f}pct.npz"))

    m.save(os.path.join(OUT, "final.npz"))
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, roi_ref=ROI_REF, rho_sim_m2=rho_sim,
                       L_nm=L_m * 1e9, n=N, deps=DEPS, relax=RELAX, dt=DT,
                       wall_s=time.time() - t0), f, indent=1)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

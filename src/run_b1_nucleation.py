"""B1: dislocation nucleation in a perfect crystal under uniaxial tension.

512x512 box, area-conserving tension to 10% in 0.25% steps, 400 relax steps
per increment. Records F(eps), core count, dislocation density; saves field
snapshots at key strains. Output: results/b1_nucleation_512/
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   "b1_nucleation_512")
os.makedirs(OUT, exist_ok=True)

DEPS = 0.0025
N_STRAIN = 40           # -> 10%
RELAX = 400
DT = 0.5


def main():
    t0 = time.time()
    m = PFC2D(512, 512, r=-0.25, psi_bar=-0.25)
    m.init_crystal()
    m.step(DT, n=300)

    rows = []
    nucleated_at = None
    for i in range(N_STRAIN):
        m.apply_strain(DEPS, area_conserving=True)
        m.step(DT, n=RELAX)
        F = m.free_energy()
        pts = find_peaks(m.psi, m.dx, m.dy)
        d = find_dislocations(pts, m.lx, m.ly)
        nc = len(d["cores"])
        rows.append(dict(exx=m.exx, F=F, cores=nc, rho=d["rho"],
                         n_peaks=len(pts)))
        print(f"exx={m.exx*100:5.2f}%  F={F:.6f}  cores={nc:3d}  "
              f"rho={d['rho']:.2e}", flush=True)
        if nc > 0 and nucleated_at is None:
            nucleated_at = m.exx
            m.save(os.path.join(OUT, f"snap_nucleation_{m.exx*100:.2f}pct.npz"))
        if i % 8 == 7:
            m.save(os.path.join(OUT, f"snap_{m.exx*100:.2f}pct.npz"))

    m.save(os.path.join(OUT, "final.npz"))
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, nucleated_at=nucleated_at,
                       wall_s=time.time() - t0, n=512, deps=DEPS,
                       relax=RELAX, dt=DT), f, indent=1)
    print(f"done, nucleation at exx={nucleated_at}, "
          f"wall {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

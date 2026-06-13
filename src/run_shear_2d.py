"""2D simple-shear mechanism study (complements the tensile M-series).

Shear is the natural resolved driver for dislocation glide and avoids the
amorphization path of uniaxial tension. Cases:
  shear_quad  pre-seeded quadrupole -> glide / annihilation under shear
  shear_poly  melt-quenched polycrystal -> shear yield, GB sliding, flow

Output: results/shear_2d/<case>/  (per-step tau, cores, core positions)
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D
from defect_analysis import find_peaks, find_dislocations

BASE = os.path.join(os.path.dirname(__file__), "..", "results", "shear_2d")
DGAMMA = 0.0025
N_STEPS = 64          # -> gamma = 16%
RELAX = 300
DT = 0.5
QUAD = [(0.33, 0.25, +1), (0.33, 0.75, -1), (0.67, 0.25, -1), (0.67, 0.75, +1)]


def make(case):
    m = PFC2D(512, 512, r=-0.25, psi_bar=-0.25)
    if case == "shear_quad":
        m.init_dislocations(QUAD)
        m.step(DT, n=600)
    else:
        m.init_random(noise=0.05, seed=7)
        m.step(DT, n=3000)
    return m


def run(case):
    out = os.path.join(BASE, case)
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    m = make(case)
    pts = find_peaks(m.psi, m.dx, m.dy)
    d = find_dislocations(pts, m.lx, m.ly)
    m.save(os.path.join(out, "initial.npz"))
    rows = [dict(gamma=0.0, tau=m.shear_stress(), F=m.free_energy(),
                 cores=len(d["cores"]), core_xy=d["cores"].tolist())]
    print(f"[{case}] gamma=0 cores={len(d['cores'])}", flush=True)
    for i in range(N_STEPS):
        m.apply_shear(DGAMMA)
        m.step(DT, n=RELAX)
        pts = find_peaks(m.psi, m.dx, m.dy)
        d = find_dislocations(pts, m.lx, m.ly)
        rows.append(dict(gamma=m.gamma, tau=m.shear_stress(),
                         F=m.free_energy(), cores=len(d["cores"]),
                         core_xy=d["cores"].tolist()))
        if i % 4 == 3:
            print(f"[{case}] gamma={m.gamma*100:5.2f}% tau={rows[-1]['tau']:+.5f} "
                  f"cores={len(d['cores'])}", flush=True)
        if i % 16 == 15:
            m.save(os.path.join(out, f"snap_{m.gamma*100:.1f}pct.npz"))
    m.save(os.path.join(out, "final.npz"))
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, dgamma=DGAMMA, relax=RELAX,
                       wall_s=time.time() - t0), f, indent=1)
    print(f"[{case}] done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    for c in (sys.argv[1:] or ["shear_quad", "shear_poly"]):
        run(c)

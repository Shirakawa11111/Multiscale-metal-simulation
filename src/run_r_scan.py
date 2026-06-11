"""r-scan: effective-temperature dependence of the multiplication threshold.

512^2 quadrupole (4 cores) strained to 16% at r in {-0.35, -0.25(ref), -0.15}.
(-0.25 reference behavior is interpolated from D2 at 1536^2; this scan keeps
box size fixed at 512^2 so only r varies within the scan.)
Output: results/r_scan_512/<tag>/
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D
from defect_analysis import find_peaks, find_dislocations

BASE = os.path.join(os.path.dirname(__file__), "..", "results", "r_scan_512")
QUAD = [(0.33, 0.25, +1), (0.33, 0.75, -1), (0.67, 0.25, -1), (0.67, 0.75, +1)]
DEPS, N_STRAIN, RELAX, DT = 0.0025, 64, 400, 0.5  # -> 16%


def run(r):
    tag = f"r{r:+.2f}".replace("+", "p").replace("-", "m").replace(".", "")
    out = os.path.join(BASE, tag)
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    m = PFC2D(512, 512, r=r, psi_bar=-0.25)
    m.init_dislocations(QUAD)
    m.step(DT, n=600)
    rows = []
    for i in range(N_STRAIN):
        m.apply_strain(DEPS, area_conserving=True)
        m.step(DT, n=RELAX)
        pts = find_peaks(m.psi, m.dx, m.dy)
        d = find_dislocations(pts, m.lx, m.ly)
        rows.append(dict(exx=m.exx, sigma=m.stress(), cores=len(d["cores"])))
        print(f"[{tag}] exx={m.exx*100:5.2f}% sigma={rows[-1]['sigma']:+.5f} "
              f"cores={len(d['cores'])}", flush=True)
        if i % 16 == 15:
            m.save(os.path.join(out, f"snap_{m.exx*100:.1f}pct.npz"))
    m.save(os.path.join(out, "final.npz"))
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, r=r, wall_s=time.time() - t0), f, indent=1)
    print(f"[{tag}] done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    for r in ([float(x) for x in sys.argv[1:]] or [-0.35, -0.15]):
        run(r)

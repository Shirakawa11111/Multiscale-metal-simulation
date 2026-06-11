"""D2 continuation: resume the ROI-matched flagship from 10.5% and strain
to 20% to map the multiplication cascade / density saturation rho(eps).
Output: results/d2_roi_matched_1536/ (rows appended in summary_part2.json)
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   "d2_roi_matched_1536")
DEPS = 0.0025
N_STRAIN = 38   # 10.5% -> ~20%
RELAX = 400
DT = 0.5


def main():
    t0 = time.time()
    m = PFC2D.load(os.path.join(OUT, "final.npz"))
    print(f"resumed at exx={m.exx*100:.2f}%", flush=True)
    rows = []
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
        prev = rows[-2]["cores"] if len(rows) > 1 else 8
        if nc != prev and n_events < 10:
            m.save(os.path.join(OUT, f"event2_{m.exx*100:.2f}pct.npz"))
            n_events += 1
        if i % 10 == 9:
            m.save(os.path.join(OUT, f"snap2_{m.exx*100:.2f}pct.npz"))

    m.save(os.path.join(OUT, "final_20pct.npz"))
    with open(os.path.join(OUT, "summary_part2.json"), "w") as f:
        json.dump(dict(rows=rows, wall_s=time.time() - t0,
                       relax=RELAX, deps=DEPS), f, indent=1)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

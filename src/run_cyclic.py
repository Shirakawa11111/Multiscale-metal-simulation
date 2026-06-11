"""Cyclic loading (fatigue protocol, cf. 单晶铜拉伸模拟): symmetric strain
cycles eps_xx in [-amp, +amp], triangle wave, quadrupole-seeded 512^2.
Tracks per-step dislocation count and stress -> hysteresis loops, and
per-cycle dislocation retention (shakedown vs accumulation).
Output: results/cyclic_512/
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "cyclic_512")
QUAD = [(0.33, 0.25, +1), (0.33, 0.75, -1), (0.67, 0.25, -1), (0.67, 0.75, +1)]
AMP = 0.03          # +/-3% strain amplitude
DEPS = 0.0025
N_CYCLES = 8
RELAX = 300
DT = 0.5


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    m = PFC2D(512, 512, r=-0.25, psi_bar=-0.25)
    m.init_dislocations(QUAD)
    m.step(DT, n=600)

    rows = []
    n_steps_quarter = int(round(AMP / DEPS))
    # triangle wave: 0 -> +amp -> -amp -> 0 per cycle
    pattern = ([+1] * n_steps_quarter + [-1] * 2 * n_steps_quarter
               + [+1] * n_steps_quarter)
    for cyc in range(N_CYCLES):
        for sgn in pattern:
            m.apply_strain(sgn * DEPS, area_conserving=True)
            m.step(DT, n=RELAX)
            pts = find_peaks(m.psi, m.dx, m.dy)
            d = find_dislocations(pts, m.lx, m.ly)
            rows.append(dict(cycle=cyc, exx=m.exx, sigma=m.stress(),
                             cores=len(d["cores"])))
        print(f"cycle {cyc}: end cores={rows[-1]['cores']} "
              f"exx={m.exx*100:.2f}% ({time.time()-t0:.0f}s)", flush=True)
        m.save(os.path.join(OUT, f"cycle_{cyc}.npz"))

    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, amp=AMP, n_cycles=N_CYCLES, relax=RELAX,
                       wall_s=time.time() - t0), f, indent=1)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

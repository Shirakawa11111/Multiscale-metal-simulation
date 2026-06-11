"""C3: 3D BCC polycrystal tension at 128^3 (production version of the
64^3 smoke). Quench -> anneal -> volume-conserving tension to 6%.
Output: results/c3_poly3d_128/
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc3d import PFC3D, A_BCC, find_peaks_3d

OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   "c3_poly3d_128")
DX = 12 * A_BCC / 128  # same resolution as the 64^3 smoke (commensurate)


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    m = PFC3D(128, 128, 128, dx=DX, r=-0.25, psi_bar=-0.25)
    m.init_random(noise=0.05, seed=11)
    m.step(0.5, n=3000)
    n0 = len(find_peaks_3d(m.psi, m.dx, m.dy, m.dz))
    print(f"quenched: peaks={n0} F={m.free_energy():.6f} "
          f"({time.time()-t0:.0f}s)", flush=True)
    m.save(os.path.join(OUT, "initial.npz"))

    rows = []
    for i in range(24):
        m.apply_strain(0.0025, volume_conserving=True)
        m.step(0.5, n=200)
        s = m.stress()
        n = len(find_peaks_3d(m.psi, m.dx, m.dy, m.dz)) if i % 4 == 3 else None
        rows.append(dict(exx=m.exx, sigma=s, F=m.free_energy(), peaks=n))
        print(f"exx={m.exx*100:5.2f}% sigma={s:+.6f}"
              + (f" peaks={n}" if n else ""), flush=True)
        if i % 8 == 7:
            m.save(os.path.join(OUT, f"snap_{m.exx*100:.2f}pct.npz"))

    m.save(os.path.join(OUT, "final.npz"))
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, n0_peaks=n0, wall_s=time.time() - t0),
                  f, indent=1)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

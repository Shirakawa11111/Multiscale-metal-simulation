"""Grain-coarsening study (tests the M15 finding directly): anneal a quenched
128^3 BCC polycrystal for increasing times and measure whether grains coarsen
(disordered_frac drops, GB network breaks into resolvable pieces).

If frac is flat with anneal time -> grain size is frozen by the quench
(M15 confirmed). If frac drops -> grains coarsen and individual dislocation
lines should eventually resolve.

Output: results/grain_coarsening/ (per-checkpoint metrics + figure)
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc3d import PFC3D, A_BCC, find_peaks_3d
from defect_analysis_3d import find_dislocation_lines

OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   "grain_coarsening")
DX = 12 * A_BCC / 128
CHECKPOINTS = [3000, 7000, 15000, 30000]   # cumulative anneal steps
SEED = 11


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    m = PFC3D(128, 128, 128, dx=DX, r=-0.25, psi_bar=-0.25)
    m.init_random(noise=0.05, seed=SEED)
    box = np.array([m.lx, m.ly, m.lz])

    rows = []
    done = 0
    for target in CHECKPOINTS:
        m.step(0.5, n=target - done)
        done = target
        pts = find_peaks_3d(m.psi, m.dx, m.dy, m.dz)
        r = find_dislocation_lines(pts, box)
        rows.append(dict(anneal_steps=done, atoms=len(pts),
                         disordered_frac=r["disordered_frac"],
                         n_lines=r["n_lines"], F=m.free_energy(),
                         top_sizes=sorted(r["line_sizes"], reverse=True)[:5]))
        print(f"anneal={done}: atoms={len(pts)} frac={r['disordered_frac']:.3f} "
              f"lines={r['n_lines']} F={m.free_energy():.6f} "
              f"({time.time()-t0:.0f}s)", flush=True)
        m.save(os.path.join(OUT, f"anneal_{done}.npz"))

    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, seed=SEED, n=128, wall_s=time.time() - t0),
                  f, indent=1)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        steps = [r["anneal_steps"] for r in rows]
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
        ax[0].plot(steps, [r["disordered_frac"] for r in rows], "o-")
        ax[0].set_xlabel("anneal steps")
        ax[0].set_ylabel("disordered fraction")
        ax[0].set_title("grain coarsening: GB content vs anneal")
        ax[0].grid(alpha=0.3)
        ax[1].plot(steps, [r["F"] for r in rows], "s-", color="tab:red")
        ax[1].set_xlabel("anneal steps")
        ax[1].set_ylabel("free energy")
        ax[1].set_title("energy relaxation")
        ax[1].grid(alpha=0.3)
        fig.savefig(os.path.join(OUT, "coarsening.png"), dpi=140,
                    bbox_inches="tight")
    except Exception as ex:
        print("plot skipped:", ex)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

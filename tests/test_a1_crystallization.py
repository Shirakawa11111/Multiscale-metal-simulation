"""Gate A1: crystallization smoke tests.

A1a (melt quench, polycrystal expected):
  - free energy decreases monotonically
  - density-based lattice spacing within 3% of 4pi/sqrt(3)
  - >=85% atoms 6-coordinated (grain boundaries account for the rest)
A1b (seeded perfect crystal):
  - >=99% atoms 6-coordinated after relaxation
  - spacing within 2% of theory
  - final energy <= polycrystal final energy
"""

import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from pfc2d import PFC2D, A_LATTICE, FFT_BACKEND
from defect_analysis import find_peaks, coordination, density_spacing

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "a1_crystallization")
os.makedirs(OUT, exist_ok=True)


def analyze(m):
    pts = find_peaks(m.psi, m.dx, m.dy)
    a_d = density_spacing(pts, m.lx, m.ly)
    coord, _ = coordination(pts, m.lx, m.ly)
    frac6 = (coord == 6).mean() if len(coord) else 0.0
    return pts, a_d, coord, frac6


def main():
    t0 = time.time()

    # --- A1a: melt quench ---
    m = PFC2D(256, 256, r=-0.25, psi_bar=-0.25)
    m.init_random(noise=0.05, seed=42)
    energies = [m.free_energy()]
    for _ in range(40):
        m.step(0.5, n=50)
        energies.append(m.free_energy())
    e = np.array(energies)
    rises = int((np.diff(e) > 1e-10).sum())
    pts, a_d, coord, frac6 = analyze(m)
    dev = abs(a_d - A_LATTICE) / A_LATTICE
    print(f"[A1a melt] F {e[0]:.6f}->{e[-1]:.6f} rises={rises} "
          f"peaks={len(pts)} a_density={a_d:.3f} (dev {100*dev:.1f}%) frac6={frac6:.3f}")
    a1a = rises == 0 and len(pts) > 500 and dev < 0.03 and frac6 >= 0.85
    m.save(os.path.join(OUT, "a1a_melt_final.npz"))

    # --- A1b: seeded perfect crystal ---
    mc = PFC2D(256, 256, r=-0.25, psi_bar=-0.25)
    mc.init_crystal()
    e0c = mc.free_energy()
    mc.step(0.5, n=500)
    e1c = mc.free_energy()
    ptsc, a_dc, coordc, frac6c = analyze(mc)
    devc = abs(a_dc - A_LATTICE) / A_LATTICE
    print(f"[A1b crystal] F {e0c:.6f}->{e1c:.6f} peaks={len(ptsc)} "
          f"a_density={a_dc:.3f} (dev {100*devc:.1f}%) frac6={frac6c:.3f}")
    a1b = frac6c >= 0.99 and devc < 0.02 and e1c <= e[-1] + 1e-9
    mc.save(os.path.join(OUT, "a1b_crystal_final.npz"))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
        ax[0].imshow(m.psi, origin="lower", cmap="viridis")
        ax[0].set_title(f"A1a melt-grown, t={m.time:.0f}")
        sc = ax[1].scatter(pts[:, 0], pts[:, 1], c=coord, s=6, cmap="coolwarm",
                           vmin=4, vmax=8)
        plt.colorbar(sc, ax=ax[1], label="coordination")
        ax[1].set_aspect("equal")
        ax[1].set_title("A1a coordination")
        ax[2].imshow(mc.psi, origin="lower", cmap="viridis")
        ax[2].set_title("A1b seeded crystal")
        fig.savefig(os.path.join(OUT, "summary.png"), dpi=130, bbox_inches="tight")
    except Exception as ex:
        print("plot skipped:", ex)

    print(f"wall {time.time()-t0:.1f}s ({FFT_BACKEND})")
    print("GATE A1a:", "PASS" if a1a else "FAIL")
    print("GATE A1b:", "PASS" if a1b else "FAIL")
    return 0 if (a1a and a1b) else 1


if __name__ == "__main__":
    sys.exit(main())

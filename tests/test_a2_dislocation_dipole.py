"""Gate A2: edge-dislocation dipole seeding and relaxation.

Seeds a +b/-b edge dipole via isotropic elastic displacement fields, relaxes,
and checks that the defect detector finds a stable, small, even number of
5|7 dislocation cores (>=2) while the rest of the lattice stays 6-coordinated.
Also records the dipole's slow attraction (climb/glide) as a sanity signal.
"""

import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "a2_dipole")
os.makedirs(OUT, exist_ok=True)


def detect(m):
    pts = find_peaks(m.psi, m.dx, m.dy)
    return pts, find_dislocations(pts, m.lx, m.ly)


def main():
    t0 = time.time()
    m = PFC2D(256, 256, r=-0.25, psi_bar=-0.25)
    m.init_dislocation_dipole()

    history = []
    m.step(0.2, n=200)  # initial hard relax of the seeded field
    for ck in range(10):
        m.step(0.5, n=100)
        pts, d = detect(m)
        history.append((m.time, d["n5"], d["n7"], len(d["cores"])))
        print(f"t={m.time:7.1f}  n5={d['n5']:3d} n7={d['n7']:3d} "
              f"cores={len(d['cores'])} rho={d['rho']:.2e}")

    pts, d = detect(m)
    n_cores = len(d["cores"])
    frac_defective = (d["n5"] + d["n7"]) / max(len(pts), 1)
    cores_stable = all(h[3] == 2 for h in history[-4:])
    ok = (n_cores == 2 and frac_defective < 0.02 and cores_stable)

    m.save(os.path.join(OUT, "final_state.npz"))
    np.save(os.path.join(OUT, "history.npy"), np.array(history))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        ax[0].imshow(m.psi, origin="lower", cmap="viridis",
                     extent=[0, m.lx, 0, m.ly])
        if len(d["cores"]):
            ax[0].plot(d["cores"][:, 0], d["cores"][:, 1], "rx", ms=12, mew=3)
        ax[0].set_title(f"psi + cores, t={m.time:.0f}")
        ax[1].plot(pts[:, 0], pts[:, 1], ".", ms=2, color="gray")
        if len(d["fives"]):
            ax[1].plot(d["fives"][:, 0], d["fives"][:, 1], "b^", ms=7, label="5")
        if len(d["sevens"]):
            ax[1].plot(d["sevens"][:, 0], d["sevens"][:, 1], "rv", ms=7, label="7")
        ax[1].legend()
        ax[1].set_aspect("equal")
        ax[1].set_title("5|7 coordination map")
        fig.savefig(os.path.join(OUT, "summary.png"), dpi=130, bbox_inches="tight")
    except Exception as ex:
        print("plot skipped:", ex)

    print(f"final: cores={n_cores}, defective-frac={frac_defective:.4f}, "
          f"wall {time.time()-t0:.1f}s")
    print("GATE A2:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

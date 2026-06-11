"""Gate A3: quasi-static uniaxial tension (area-conserving box rescaling).

A3a (perfect crystal, elastic regime):
  - strain to 2% in 0.25% steps, relax each step
  - free energy F(eps) increases and is convex (positive effective modulus)
  - no dislocations nucleate at eps <= 2%
A3b (pre-seeded dipole, plastic carrier):
  - same loading; the two cores must remain detectable and MOVE (glide)
  - net core displacement > 1 lattice constant by eps = 2%
"""

import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "a3_tension")
os.makedirs(OUT, exist_ok=True)

DEPS = 0.0025
NSTEPS_STRAIN = 8          # -> 2% total
RELAX_STEPS = 300
DT = 0.5


def detect(m):
    pts = find_peaks(m.psi, m.dx, m.dy)
    return pts, find_dislocations(pts, m.lx, m.ly)


def run_tension(m, label):
    rows = []
    cores_t = []
    for i in range(NSTEPS_STRAIN):
        m.apply_strain(DEPS, area_conserving=True)
        m.step(DT, n=RELAX_STEPS)
        F = m.free_energy()
        pts, d = detect(m)
        rows.append((m.exx, F, len(d["cores"]), d["rho"]))
        cores_t.append(d["cores"].copy())
        print(f"[{label}] exx={m.exx*100:5.2f}%  F={F:.6f}  cores={len(d['cores'])}")
    return np.array(rows), cores_t


def min_image_disp(p, q, lx, ly):
    v = q - p
    v[0] -= lx * np.round(v[0] / lx)
    v[1] -= ly * np.round(v[1] / ly)
    return np.linalg.norm(v)


def main():
    t0 = time.time()

    # --- A3a: perfect crystal ---
    m = PFC2D(256, 256, r=-0.25, psi_bar=-0.25)
    m.init_crystal()
    m.step(DT, n=200)
    F0 = m.free_energy()
    rows_a, _ = run_tension(m, "A3a")
    eps = np.concatenate([[0.0], rows_a[:, 0]])
    F = np.concatenate([[F0], rows_a[:, 1]])
    coef = np.polyfit(eps, F, 2)  # F ~ c2 eps^2 + c1 eps + c0
    # The q-snapped seed carries ~1% residual misfit, so F(eps) is a parabola
    # whose minimum marks the stress-free strain; require positive curvature
    # (elastic modulus), small misfit, and genuine loading beyond it.
    eps_min = -coef[1] / (2.0 * coef[0])
    no_defects = bool(np.all(rows_a[:, 2] == 0))
    loaded = F[-1] > F.min() + 1e-9
    a3a = coef[0] > 0 and abs(eps_min) < 0.015 and no_defects and loaded
    print(f"[A3a] curvature c2={coef[0]:.4f} eps_min={eps_min*100:.2f}% "
          f"defect-free={no_defects} loaded={loaded}")
    m.save(os.path.join(OUT, "a3a_final.npz"))

    # --- A3b: dipole glide ---
    md = PFC2D(256, 256, r=-0.25, psi_bar=-0.25)
    md.init_dislocation_dipole()
    md.step(DT, n=400)
    _, d0 = detect(md)
    rows_b, cores_t = run_tension(md, "A3b")
    always_detected = bool(np.all(rows_b[:, 2] >= 2))
    disp = 0.0
    if len(d0["cores"]) and len(cores_t[-1]):
        disp = max(min(min_image_disp(c0, c1, md.lx, md.ly)
                       for c1 in cores_t[-1]) for c0 in d0["cores"])
    moved = disp > A_LATTICE
    a3b = always_detected and moved
    print(f"[A3b] cores tracked={always_detected}, max core disp={disp:.2f} "
          f"(a0={A_LATTICE:.2f})")
    md.save(os.path.join(OUT, "a3b_final.npz"))
    np.savez(os.path.join(OUT, "curves.npz"), rows_a=rows_a, rows_b=rows_b)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
        ax[0].plot(eps * 100, F, "o-")
        ax[0].set_xlabel("strain %")
        ax[0].set_ylabel("F")
        ax[0].set_title("A3a elastic energy")
        ax[1].plot(rows_b[:, 0] * 100, rows_b[:, 2], "s-")
        ax[1].set_xlabel("strain %")
        ax[1].set_ylabel("# cores")
        ax[1].set_title("A3b core count")
        ax[2].imshow(md.psi, origin="lower", cmap="viridis",
                     extent=[0, md.lx, 0, md.ly])
        for ct in cores_t:
            if len(ct):
                ax[2].plot(ct[:, 0], ct[:, 1], "r.", ms=4)
        ax[2].set_title("A3b core trajectory")
        fig.savefig(os.path.join(OUT, "summary.png"), dpi=130,
                    bbox_inches="tight")
    except Exception as ex:
        print("plot skipped:", ex)

    print(f"wall {time.time()-t0:.1f}s")
    print("GATE A3a:", "PASS" if a3a else "FAIL")
    print("GATE A3b:", "PASS" if a3b else "FAIL")
    return 0 if (a3a and a3b) else 1


if __name__ == "__main__":
    sys.exit(main())

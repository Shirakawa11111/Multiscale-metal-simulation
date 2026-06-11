"""Gate C1: 3D PFC BCC crystallization smoke test (64^3).

C1a (seeded one-mode BCC): relax, energy below liquid, peak count matches
  BCC site density (2 atoms / a^3 cell) within 10%, NN distance ~ sqrt(3)/2 a.
C1b (melt quench): crystallizes from noise (peaks form, energy drops below
  the uniform-liquid value).
"""

import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from scipy.spatial import cKDTree
from pfc3d import PFC3D, A_BCC, find_peaks_3d, FFT_BACKEND

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "c1_bcc_3d")
os.makedirs(OUT, exist_ok=True)


def nn_distance(pts, lx, ly, lz):
    if len(pts) < 9:
        return np.nan
    tree = cKDTree(pts, boxsize=[lx, ly, lz])
    d, _ = tree.query(pts, k=2)
    return float(np.median(d[:, 1]))


DX_C = 6 * A_BCC / 64  # 6 BCC cells across 64 grid points (commensurate)


def main():
    t0 = time.time()

    # liquid reference energy (uniform psi_bar)
    m0 = PFC3D(64, 64, 64, dx=DX_C, r=-0.25, psi_bar=-0.25)
    f_liquid = m0.free_energy()

    # --- C1a seeded BCC ---
    m = PFC3D(64, 64, 64, dx=DX_C, r=-0.25, psi_bar=-0.25)
    m.init_crystal()
    e0 = m.free_energy()
    m.step(0.5, n=400)
    e1 = m.free_energy()
    pts = find_peaks_3d(m.psi, m.dx, m.dy, m.dz)
    v = m.lx * m.ly * m.lz
    n_expect = 2.0 * v / A_BCC ** 3
    nn = nn_distance(pts, m.lx, m.ly, m.lz)
    nn_theory = np.sqrt(3.0) / 2.0 * A_BCC
    dev_n = abs(len(pts) - n_expect) / n_expect
    dev_nn = abs(nn - nn_theory) / nn_theory
    print(f"[C1a] F {e0:.6f}->{e1:.6f} (liquid {f_liquid:.6f}) "
          f"peaks={len(pts)} expect={n_expect:.0f} (dev {100*dev_n:.1f}%) "
          f"NN={nn:.3f} theory={nn_theory:.3f} (dev {100*dev_nn:.1f}%)")
    c1a = e1 < f_liquid and dev_n < 0.10 and dev_nn < 0.05
    m.save(os.path.join(OUT, "c1a_bcc.npz"))

    # --- C1b melt quench ---
    mq = PFC3D(64, 64, 64, dx=DX_C, r=-0.25, psi_bar=-0.25)
    mq.init_random(noise=0.05, seed=3)
    mq.step(0.5, n=1500)
    eq = mq.free_energy()
    ptsq = find_peaks_3d(mq.psi, mq.dx, mq.dy, mq.dz)
    print(f"[C1b] F={eq:.6f} (liquid {f_liquid:.6f}) peaks={len(ptsq)}")
    c1b = eq < f_liquid and len(ptsq) > 0.5 * n_expect
    mq.save(os.path.join(OUT, "c1b_quench.npz"))

    print(f"wall {time.time()-t0:.1f}s ({FFT_BACKEND})")
    print("GATE C1a:", "PASS" if c1a else "FAIL")
    print("GATE C1b:", "PASS" if c1b else "FAIL")
    return 0 if (c1a and c1b) else 1


if __name__ == "__main__":
    sys.exit(main())

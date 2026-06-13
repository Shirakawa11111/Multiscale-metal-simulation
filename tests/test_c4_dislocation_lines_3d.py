"""Gate C4: 3D dislocation-line detector validation (CSP metric).

C4a (perfect BCC): disordered fraction < 2%, zero lines.
C4b (BCC with a melted cylindrical void): the void rim is flagged as a
  connected disordered region (>=1 line cluster, disordered_frac jumps).
C4c (strained perfect crystal, 2%): detector reports ~0 lines.
C4d (synthetic affine invariance): an ideal BCC lattice under a known 15%
  volume-conserving affine F must give CSP ~ 0 (the property that makes the
  metric valid under large strain, unlike a fixed-distance metric).
"""

import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from pfc3d import PFC3D, A_BCC, find_peaks_3d
from defect_analysis_3d import find_dislocation_lines, disorder_metric

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "c4_lines_3d")
os.makedirs(OUT, exist_ok=True)
DX = 6 * A_BCC / 64


def main():
    t0 = time.time()

    # --- C4a perfect ---
    m = PFC3D(64, 64, 64, dx=DX, r=-0.25, psi_bar=-0.25)
    m.init_crystal()
    m.step(0.5, n=300)
    pts = find_peaks_3d(m.psi, m.dx, m.dy, m.dz)
    box = np.array([m.lx, m.ly, m.lz])
    r = find_dislocation_lines(pts, box)
    print(f"[C4a perfect] frac={r['disordered_frac']:.4f} lines={r['n_lines']}")
    c4a = r["disordered_frac"] < 0.02 and r["n_lines"] == 0

    # --- C4b void ---
    mv = PFC3D(64, 64, 64, dx=DX, r=-0.25, psi_bar=-0.25)
    mv.init_crystal()
    # carve a cylindrical liquid void along z
    x = np.arange(64) * mv.dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rr = np.sqrt((X - mv.lx / 2) ** 2 + (Y - mv.ly / 2) ** 2)
    void = rr < 1.3 * A_BCC   # thin cylinder along z
    # deep mass depletion keeps the pore open (shallow voids dissolve by
    # Gibbs-Thomson, cf. M7); short relax leaves a persistent rim defect line
    mv.psi[void] = mv.psi_bar - 0.6
    mv.step(0.5, n=150)
    ptsv = find_peaks_3d(mv.psi, mv.dx, mv.dy, mv.dz)
    rv = find_dislocation_lines(ptsv, box)
    print(f"[C4b void]    frac={rv['disordered_frac']:.4f} lines={rv['n_lines']} "
          f"sizes={rv['line_sizes'][:5]}")
    c4b = rv["n_lines"] >= 1 and rv["disordered_frac"] > r["disordered_frac"]
    np.savez(os.path.join(OUT, "void_labels.npz"), pts=ptsv,
             labels=rv["labels"], disorder=rv["disorder"])

    # --- C4c strained perfect (strain invariance) ---
    ms = PFC3D(64, 64, 64, dx=DX, r=-0.25, psi_bar=-0.25)
    ms.init_crystal()
    ms.step(0.5, n=200)
    for _ in range(8):
        ms.apply_strain(0.0025, volume_conserving=True)
        ms.step(0.5, n=150)
    pts_s = find_peaks_3d(ms.psi, ms.dx, ms.dy, ms.dz)
    box_s = np.array([ms.lx, ms.ly, ms.lz])
    rs = find_dislocation_lines(pts_s, box_s)
    print(f"[C4c strained] exx={ms.exx*100:.1f}% frac={rs['disordered_frac']:.4f} "
          f"lines={rs['n_lines']}")
    c4c = rs["n_lines"] == 0 and rs["disordered_frac"] < 0.03

    # --- C4d synthetic affine invariance ---
    from defect_analysis_3d import disorder_metric
    nc = 5
    base = []
    for i in range(nc):
        for j in range(nc):
            for k in range(nc):
                base.append([i, j, k])
                base.append([i + 0.5, j + 0.5, k + 0.5])
    bp = np.array(base, dtype=float)
    bbox = np.array([nc, nc, nc], dtype=float)
    ex = 0.15
    fy = 1.0 / np.sqrt(1 + ex)
    F = np.diag([1 + ex, fy, fy])
    csp_strained = disorder_metric(bp @ F.T, bbox * np.array([1 + ex, fy, fy]),
                                   a0=1.0)
    print(f"[C4d affine] CSP max at 15% affine strain = {csp_strained.max():.2e}")
    c4d = csp_strained.max() < 1e-3

    print(f"wall {time.time()-t0:.1f}s")
    for tag, ok in [("C4a", c4a), ("C4b", c4b), ("C4c", c4c), ("C4d", c4d)]:
        print(f"GATE {tag}:", "PASS" if ok else "FAIL")
    return 0 if (c4a and c4b and c4c and c4d) else 1


if __name__ == "__main__":
    sys.exit(main())

"""Gate C2: 3D BCC elastic tension smoke (64^3, volume-conserving).

Seeded BCC strained to 2% in 0.25% steps: stress-strain slope positive,
field stays finite and crystalline (peak count within 5% of initial),
stress at 2% positive (loaded past any snap-misfit minimum).
"""

import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from pfc3d import PFC3D, A_BCC, find_peaks_3d

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "c2_bcc_tension")
os.makedirs(OUT, exist_ok=True)
DX_C = 6 * A_BCC / 64


def main():
    t0 = time.time()
    m = PFC3D(64, 64, 64, dx=DX_C, r=-0.25, psi_bar=-0.25)
    m.init_crystal()
    m.step(0.5, n=300)
    n0 = len(find_peaks_3d(m.psi, m.dx, m.dy, m.dz))

    rows = []
    for i in range(8):
        m.apply_strain(0.0025, volume_conserving=True)
        m.step(0.5, n=200)
        s = m.stress()
        rows.append((m.exx, m.free_energy(), s))
        print(f"exx={m.exx*100:5.2f}%  F={rows[-1][1]:.6f}  sigma={s:+.6f}")

    rows = np.array(rows)
    n1 = len(find_peaks_3d(m.psi, m.dx, m.dy, m.dz))
    slope = np.polyfit(rows[:, 0], rows[:, 2], 1)[0]
    finite = np.isfinite(rows).all()
    peaks_ok = abs(n1 - n0) / max(n0, 1) < 0.05
    print(f"peaks {n0}->{n1}, modulus={slope:.4f}, final sigma={rows[-1,2]:+.5f}")
    m.save(os.path.join(OUT, "final.npz"))
    np.save(os.path.join(OUT, "curve.npy"), rows)

    ok = finite and peaks_ok and slope > 0 and rows[-1, 2] > 0
    print(f"wall {time.time()-t0:.1f}s")
    print("GATE C2:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

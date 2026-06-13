"""Proper mobility re-measurement (addresses the adversarial-review fixes for
the flagged-weak mobility result):
  - LARGE box (768^2) so a glide of MANY Burgers vectors is captured (not the
    ~1b quantization staircase of the earlier 256^2 attempt)
  - sub-grid core tracking via the refined peak finder
  - >=6 stress levels
  - fit v against the TOTAL resolved Peach-Koehler stress
       tau_eff = tau_applied - tau_mutual,  tau_mutual = mu*b/(2 pi (1-nu) L)
    (the dipole's own attraction, subtracted analytically) — so the slope is
    the intrinsic mobility, not contaminated by the pair interaction.
Reports v(tau_eff): linear (overdamped) form, mobility M, and an honest R^2.
Output: results/mobility_proper/
"""

import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "mobility_proper")
N = 768
DT = 0.5
MU = 0.0545          # measured PFC shear modulus
NU = 1.0 / 3.0
B = A_LATTICE


def core_x_near(m, y_target):
    d = find_dislocations(find_peaks(m.psi, m.dx, m.dy), m.lx, m.ly)
    c = d["cores"]
    if len(c) == 0:
        return None, 0
    dy = np.abs((c[:, 1] - y_target + m.ly / 2) % m.ly - m.ly / 2)
    return float(c[np.argmin(dy)][0]), len(c)


def glide_at(gamma_hold, mpfc=True):
    m = PFC2D(N, N, r=-0.25, psi_bar=-0.25)
    # wide Y separation (0.25 ly) so mutual attraction is weak & glide is long
    m.init_dislocations([(0.5, 0.32, +1), (0.5, 0.68, -1)])
    (m.step_mpfc(DT, n=700, beta=10.0) if mpfc else m.step(DT, n=700))
    y_lo, y_hi = 0.32 * m.ly, 0.68 * m.ly
    for _ in range(int(round(gamma_hold / 0.0025))):
        m.apply_shear(0.0025)
        (m.step_mpfc(DT, n=40, beta=10.0) if mpfc else m.step(DT, n=40))
    tau_app = m.shear_stress()
    xl0, _ = core_x_near(m, y_lo)
    xh0, _ = core_x_near(m, y_hi)
    if xl0 is None or xh0 is None:
        return None
    track = []
    t0 = m.time
    for i in range(30):
        (m.step_mpfc(DT, n=60, beta=10.0) if mpfc else m.step(DT, n=60))
        xl, nc = core_x_near(m, y_lo)
        xh, _ = core_x_near(m, y_hi)
        if xl is None or xh is None or nc < 2:
            break
        dlo = (xl - xl0 + m.lx / 2) % m.lx - m.lx / 2
        dhi = (xh - xh0 + m.lx / 2) % m.lx - m.lx / 2
        track.append((m.time - t0, 0.5 * (abs(dlo) + abs(dhi))))
    track = np.array(track)
    if len(track) < 8 or track[-1, 1] < 2 * B:    # require >2b of glide
        return None
    v = float(np.polyfit(track[:, 0], track[:, 1], 1)[0])
    # mutual attraction at the (fixed) Y separation L = 0.36 ly
    L = 0.36 * m.ly
    tau_mutual = MU * B / (2 * np.pi * (1 - NU) * L)
    return dict(tau_app=float(tau_app), tau_mutual=float(tau_mutual),
                tau_eff=float(tau_app - tau_mutual), v=v,
                glide_total=float(track[-1, 1]), n_pts=len(track))


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    pts = []
    for gh in (0.02, 0.03, 0.04, 0.05, 0.06, 0.08):
        r = glide_at(gh, mpfc=True)
        if r:
            pts.append(r)
            print(f"  gamma={gh:.3f}: tau_app={r['tau_app']:.5f} "
                  f"tau_eff={r['tau_eff']:.5f} v={r['v']:.5f} "
                  f"glide={r['glide_total']/B:.1f}b", flush=True)
    res = {"points": pts}
    if len(pts) >= 4:
        te = np.array([p["tau_eff"] for p in pts])
        v = np.array([p["v"] for p in pts])
        M = float(np.sum(v * te) / np.sum(te ** 2) / B)   # through-origin in tau_eff
        r2 = 1 - np.var(v - M * B * te) / (np.var(v) + 1e-12)
        # also a free-intercept fit to expose any residual threshold
        slope, intc = np.polyfit(te, v, 1)
        res.update(mobility_M=M, linear_r2=float(r2),
                   free_fit_slope_over_b=float(slope / B),
                   free_fit_intercept=float(intc))
        print(f"\nPROPER MOBILITY (vs tau_eff): M={M:.3f}, R^2={r2:.3f}; "
              f"free-fit intercept={intc:.5f} (≈0 confirms no spurious threshold)",
              flush=True)
    with open(os.path.join(OUT, "mobility_proper.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

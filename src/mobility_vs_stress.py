"""Sub-problem B step 1: measure the PFC glide mobility law v(tau) across
several resolved-shear levels, test for linear (Newtonian/overdamped) drag
v = M*tau*b — the form DDD/ExaDiS assumes — and convert M to physical units
via the Interface-B scale map, yielding a transferable dislocation drag
coefficient for the downstream discrete-dislocation model.

Holds a sequence of fixed shear levels; at each, glides a Y-separated edge
dipole and fits the steady velocity. Output: results/mobility_law/
"""

import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D, A_LATTICE

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "mobility_law")
# Interface-B scale map (from interfaceB_bridge): 1 PFC length = b_Cu/a0
B_CU = 2.556e-10
SCALE_L = B_CU / A_LATTICE          # m per PFC length unit


def core_x(m, y_target):
    from defect_analysis import find_peaks, find_dislocations
    d = find_dislocations(find_peaks(m.psi, m.dx, m.dy), m.lx, m.ly)
    c = d["cores"]
    if len(c) == 0:
        return None, len(c)
    dy = np.abs((c[:, 1] - y_target + m.ly / 2) % m.ly - m.ly / 2)
    return float(c[np.argmin(dy)][0]), len(c)


def glide_velocity_at(gamma_hold, mpfc=False, n=256, r=-0.25):
    m = PFC2D(n, n, r=r, psi_bar=-0.25)
    m.init_dislocations([(0.5, 0.3, +1), (0.5, 0.7, -1)])
    (m.step_mpfc(0.5, n=600, beta=10.0) if mpfc else m.step(0.5, n=600))
    y_lo, y_hi = 0.3 * m.ly, 0.7 * m.ly
    for _ in range(int(round(gamma_hold / 0.0025))):
        m.apply_shear(0.0025)
        (m.step_mpfc(0.5, n=40, beta=10.0) if mpfc else m.step(0.5, n=40))
    tau = m.shear_stress()
    xl0, _ = core_x(m, y_lo)
    xh0, _ = core_x(m, y_hi)
    if xl0 is None or xh0 is None:
        return None
    track = []
    t0 = m.time
    for i in range(20):
        (m.step_mpfc(0.5, n=60, beta=10.0) if mpfc else m.step(0.5, n=60))
        xl, nc = core_x(m, y_lo)
        xh, _ = core_x(m, y_hi)
        if xl is None or xh is None or nc < 2:
            break
        dlo = (xl - xl0 + m.lx / 2) % m.lx - m.lx / 2
        dhi = (xh - xh0 + m.lx / 2) % m.lx - m.lx / 2
        track.append((m.time - t0, 0.5 * (abs(dlo) + abs(dhi))))
    track = np.array(track)
    if len(track) < 6:
        return None
    v = float(np.polyfit(track[:, 0], track[:, 1], 1)[0])
    return dict(tau=float(tau), v=v, gamma_hold=gamma_hold)


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    res = {}
    for mpfc, label in [(False, "diffusive"), (True, "mpfc")]:
        pts = []
        for gh in (0.015, 0.025, 0.035, 0.05, 0.07):
            r = glide_velocity_at(gh, mpfc=mpfc)
            if r:
                pts.append(r)
                print(f"[{label}] gamma={gh:.3f} tau={r['tau']:.5f} v={r['v']:.5f}",
                      flush=True)
        if len(pts) >= 3:
            tau = np.array([p["tau"] for p in pts])
            v = np.array([p["v"] for p in pts])
            # linear-drag fit v = M*b*tau (through origin) + check curvature
            b = A_LATTICE
            M = float(np.sum(v * tau) / np.sum(tau ** 2) / b)   # LS through 0
            r2 = 1 - np.var(v - M * b * tau) / (np.var(v) + 1e-12)
            res[label] = dict(points=pts, mobility_M_pfc=M, linear_r2=r2)
            print(f"[{label}] linear-drag M = {M:.3f} PFC units, R^2={r2:.3f}",
                  flush=True)
    # physical conversion (PFC mobility -> note: needs a PFC time<->s map to be
    # fully physical; here we report the dimensionless M and the length scale,
    # and the drag in PFC units, flagging the time-calibration as the open item)
    res["scale"] = dict(length_m_per_pfc=SCALE_L, b_cu_m=B_CU,
                        note="M is dimensionless PFC mobility; physical drag "
                             "B=1/M needs a PFC-time<->s calibration (open).")
    with open(os.path.join(OUT, "mobility_law.json"), "w") as f:
        json.dump(res, f, indent=1)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(7, 5.5))
        for label in ("diffusive", "mpfc"):
            if label in res:
                pts = res[label]["points"]
                tau = np.array([p["tau"] for p in pts])
                v = np.array([p["v"] for p in pts])
                M = res[label]["mobility_M_pfc"]
                plt.plot(tau, v, "o", label=f"{label} data")
                tt = np.linspace(0, tau.max(), 20)
                plt.plot(tt, M * A_LATTICE * tt, "--",
                         label=f"{label} fit M={M:.2f} (R2={res[label]['linear_r2']:.2f})")
        plt.xlabel("resolved shear stress tau (PFC units)")
        plt.ylabel("glide velocity v (PFC units)")
        plt.title("PFC dislocation mobility law v(tau) — linear-drag test")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.savefig(os.path.join(OUT, "mobility_law.png"), dpi=140,
                    bbox_inches="tight")
    except Exception as ex:
        print("plot skipped:", ex)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

"""Decisive test for the 'hardening forest PFC' program: can 2D triangular
PFC realize Taylor forest hardening (positive alpha) with (a) a proper
intersecting-slip-system forest and (b) elastically-relaxed (MPFC) dynamics?

Earlier Taylor attempt FAILED (alpha=-0.287, softening) because it seeded
annihilating dipoles on the SAME glide system. Here:
  - mobile population: edge dipoles with Burgers along x (0deg) — driven by
    the applied epsilon_xy shear (max Schmid factor)
  - FOREST: dipoles with Burgers at 60deg and 120deg (the other two triangular
    glide systems) — low resolved shear under epsilon_xy, so they act as a
    relatively-immobile forest the mobile dislocations must cut through
Vary forest density; measure flow stress under MPFC. If tau_flow rises with
sqrt(rho_forest) (alpha>0), forest hardening is realized in 2D and the heavy
XPFC-FCC upgrade may be unnecessary; if it still softens, FCC/two-mode is
required.

Output: results/forest_hardening/
"""

import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   "forest_hardening")
N = 384
DT = 0.5
DGAMMA = 0.0025
N_SHEAR = 32          # -> 8% shear
RELAX = 150
BETA = 10.0


def measure_mu():
    m = PFC2D(256, 256, r=-0.25, psi_bar=-0.25)
    m.init_crystal()
    m.step_mpfc(DT, n=300, beta=BETA)
    g, t = [], []
    for _ in range(5):
        m.apply_shear(0.0025)
        m.step_mpfc(DT, n=250, beta=BETA)
        g.append(m.gamma)
        t.append(m.shear_stress())
    return float(np.polyfit(g, t, 1)[0])


def build(n_forest, seed):
    """2 mobile x-dipoles + n_forest dipoles split over 60/120-deg systems."""
    rng = np.random.default_rng(seed)
    cores = [(0.3, 0.35, +1, 0.0), (0.3, 0.65, -1, 0.0),
             (0.7, 0.35, +1, 0.0), (0.7, 0.65, -1, 0.0)]
    for i in range(n_forest):
        ang = 60.0 if i % 2 == 0 else 120.0
        x = rng.uniform(0.1, 0.9)
        y = rng.uniform(0.15, 0.45)
        cores.append((x, y, +1, ang))
        cores.append((x, min(0.88, y + 0.4), -1, ang))
    m = PFC2D(N, N, r=-0.25, psi_bar=-0.25); m.glide_kc=float(os.environ.get("GLIDE_KC","0"))
    m.init_dislocations(cores)
    m.step_mpfc(DT, n=700, beta=BETA)
    return m


def flow_stress(m):
    taus = []
    for i in range(N_SHEAR):
        m.apply_shear(DGAMMA)
        m.step_mpfc(DT, n=RELAX, beta=BETA)
        taus.append(m.shear_stress())
    taus = np.array(taus)
    return float(np.mean(taus[-N_SHEAR // 3:])), taus.tolist()


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    mu = measure_mu()
    b = A_LATTICE
    print(f"mu (MPFC) = {mu:.4f}", flush=True)
    rows = []
    for n_forest in (0, 4, 8, 16, 24):
        for seed in (7, 11):
            m = build(n_forest, seed)
            d = find_dislocations(find_peaks(m.psi, m.dx, m.dy), m.lx, m.ly)
            n0 = len(d["cores"])
            rho = n0 / (m.lx * m.ly)
            tf, hist = flow_stress(m)
            rows.append(dict(n_forest=n_forest, seed=seed, n_cores=n0,
                             rho=rho, flow_stress=tf))
            print(f"  n_forest={n_forest} seed={seed}: cores={n0} "
                  f"rho={rho:.2e} tau_flow={tf:.5f}", flush=True)
    rho = np.array([r["rho"] for r in rows])
    tf = np.array([r["flow_stress"] for r in rows])
    x = np.sqrt(rho)
    slope, tau0 = np.polyfit(x, tf, 1)
    alpha = slope / (mu * b)
    r2 = 1 - np.var(tf - (slope * x + tau0)) / (np.var(tf) + 1e-12)
    result = dict(mu=mu, alpha=float(alpha), tau0=float(tau0),
                  taylor_r2=float(r2), rows=rows,
                  verdict=("HARDENS (alpha>0): forest hardening realized in 2D"
                           if alpha > 0 else
                           "still SOFTENS (alpha<0): needs XPFC-FCC / two-mode"))
    with open(os.path.join(OUT, "forest_hardening.json"), "w") as f:
        json.dump(result, f, indent=1)
    print(f"\nFOREST HARDENING: alpha = {alpha:.3f}, R^2 = {r2:.3f}\n"
          f"  -> {result['verdict']}", flush=True)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

"""DUAL-FIX culminating test: two-mode crystallography (square lattice, supplies
intersecting slip systems + junction geometry) + climb-suppressed dynamics
(glide_kc) — does the flow stress finally RISE with dislocation density
(Taylor alpha > 0)?

Density = SEEDED forest count (lattice-agnostic; sidesteps the triangular-only
5|7 detector which is ill-defined on a square lattice). Flow stress is the
lattice-agnostic shear_stress plateau. Sweep forest density x glide_kc.

Square slip systems: <10> Burgers at 0/90 deg, <11> at 45/135 deg.
Output: hardening_program/results_twomode/
"""

import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from pfc_twomode import TwoModePFC

OUT = os.path.join(os.path.dirname(__file__), "results_twomode")
N = 384
DT = 0.5
DG = 0.0025
N_SHEAR = 32
RELAX = 150
BETA = 10.0
A_SQ = 2 * np.pi   # square lattice constant (q=1)


def build(n_forest, kc, seed):
    rng = np.random.default_rng(seed)
    # mobile <10> dipoles (Burgers 0deg) + forest on <11> (45/135) and <10>(90)
    cores = [(0.3, 0.30, +1, 0.0), (0.3, 0.70, -1, 0.0),
             (0.7, 0.30, +1, 0.0), (0.7, 0.70, -1, 0.0)]
    fang = [45.0, 135.0, 90.0]
    for i in range(n_forest):
        a = fang[i % 3]
        x = rng.uniform(0.12, 0.88)
        y = rng.uniform(0.12, 0.42)
        cores.append((x, y, +1, a))
        cores.append((x, min(0.88, y + 0.45), -1, a))
    m = TwoModePFC(N, N)
    m.glide_kc = kc
    m.init_dislocations_square(cores)
    m.step_mpfc(DT, n=700, beta=BETA)
    return m


def flow_stress(m):
    taus = []
    for i in range(N_SHEAR):
        m.apply_shear(DG)
        m.step_mpfc(DT, n=RELAX, beta=BETA)
        taus.append(m.shear_stress())
    taus = np.array(taus)
    return float(np.mean(taus[-N_SHEAR // 3:]))


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    mu_m = TwoModePFC(192, 192)
    mu_m.init_square(amp=0.35, amp2=0.18)
    mu_m.step(DT, n=300)
    g, t = [], []
    for _ in range(5):
        mu_m.apply_shear(0.0025)
        mu_m.step(DT, n=250)
        g.append(mu_m.gamma)
        t.append(mu_m.shear_stress())
    mu = float(np.polyfit(g, t, 1)[0])
    print(f"two-mode square shear modulus mu={mu:.4f}", flush=True)

    res = {}
    for kc in (0.0, 1.0):
        rows = []
        for nf in (0, 6, 12, 20):
            tf = np.mean([flow_stress(build(nf, kc, s)) for s in (7, 11)])
            rho = (4 + 2 * nf) / (N * mu_m.dx0 * (1 + 0) * N * mu_m.dx0)  # rough areal density (seeded)
            # use seeded core count directly as density proxy (box fixed)
            rho = (4 + 2 * nf)
            rows.append(dict(n_forest=nf, seeded_cores=4 + 2 * nf,
                             flow_stress=float(tf)))
            print(f"  kc={kc} nf={nf}: seeded={4+2*nf} tau_flow={tf:.5f}",
                  flush=True)
        sc = np.array([r["seeded_cores"] for r in rows], float)
        tf = np.array([r["flow_stress"] for r in rows])
        # Taylor: tau ~ sqrt(rho); rho proportional to seeded count
        slope, intc = np.polyfit(np.sqrt(sc), tf, 1)
        res[f"kc{kc}"] = dict(rows=rows, taylor_slope=float(slope),
                              hardens=bool(slope > 0))
        print(f"  kc={kc}: d(tau)/d(sqrt(N)) = {slope:+.6f} "
              f"=> {'HARDENS' if slope > 0 else 'softens'}", flush=True)
    res["mu"] = mu
    with open(os.path.join(OUT, "twomode_forest.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(f"\nDUAL-FIX VERDICT:")
    for kc in (0.0, 1.0):
        r = res[f"kc{kc}"]
        print(f"  two-mode + climb_kc={kc}: "
              f"{'HARDENS (alpha>0)!' if r['hardens'] else 'still softens'} "
              f"(slope={r['taylor_slope']:+.6f})")
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

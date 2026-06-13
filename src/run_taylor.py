"""Sub-problem B loop-closure (dimensionless, no time-calibration needed):
extract the Taylor forest-hardening coefficient alpha from PFC.

Taylor relation: tau_flow = tau_0 + alpha * mu * b * sqrt(rho), where alpha is
the dimensionless coefficient DAMASK's dislotwin law normally FITS/ASSUMES
(~0.3-0.5 for metals). We seed 2D crystals with a controlled forest density
rho = N_dislocations / area, shear each, read the flow stress, and fit
tau_flow vs sqrt(rho). The slope / (mu*b) = alpha is a transferable,
dimensionless PFC->crystal-plasticity constitutive parameter — closing a real
upscaling loop without needing the PFC-time<->seconds map.

mu (PFC shear modulus) is measured independently (perfect-crystal shear).
Output: results/taylor/  (flow stress vs rho, alpha fit)
"""

import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "taylor")
N = 512
DT = 0.5
DGAMMA = 0.0025
N_SHEAR = 40       # -> 10% shear
RELAX = 200


def measure_mu():
    """PFC shear modulus from a perfect crystal (dtau/dgamma)."""
    m = PFC2D(256, 256, r=-0.25, psi_bar=-0.25)
    m.init_crystal()
    m.step(DT, n=300)
    g, t = [], []
    for _ in range(6):
        m.apply_shear(0.0025)
        m.step(DT, n=300)
        g.append(m.gamma)
        t.append(m.shear_stress())
    return float(np.polyfit(g, t, 1)[0])


def seed_forest(n_dipoles, seed):
    """Place n_dipoles random ±b edge dipoles (forest density 2*n_dipoles/area)."""
    rng = np.random.default_rng(seed)
    cores = []
    for _ in range(n_dipoles):
        x = rng.uniform(0.1, 0.9)
        y = rng.uniform(0.15, 0.45)
        cores.append((x, y, +1))
        cores.append((x, min(0.85, y + 0.4), -1))   # partner on a parallel plane
    # enforce net-zero already balanced (+/- pairs)
    m = PFC2D(N, N, r=-0.25, psi_bar=-0.25)
    m.init_dislocations(cores)
    m.step(DT, n=600)
    return m


def flow_stress(m):
    """Shear to N_SHEAR; flow stress = mean tau over the last third (plateau)."""
    taus = []
    for i in range(N_SHEAR):
        m.apply_shear(DGAMMA)
        m.step(DT, n=RELAX)
        taus.append(m.shear_stress())
    taus = np.array(taus)
    return float(np.mean(taus[-N_SHEAR // 3:])), taus.tolist()


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    mu = measure_mu()
    print(f"PFC shear modulus mu = {mu:.4f}", flush=True)
    b = A_LATTICE
    rows = []
    for n_dip in (2, 4, 8, 16, 24):
        for seed in (7, 11):
            m = seed_forest(n_dip, seed)
            pts = find_peaks(m.psi, m.dx, m.dy)
            d = find_dislocations(pts, m.lx, m.ly)
            n0 = len(d["cores"])
            area = m.lx * m.ly
            rho = n0 / area
            tf, hist = flow_stress(m)
            rows.append(dict(n_dip=n_dip, seed=seed, n_cores=n0,
                             rho=rho, flow_stress=tf))
            print(f"  n_dip={n_dip} seed={seed}: cores={n0} rho={rho:.2e} "
                  f"tau_flow={tf:.5f}", flush=True)

    # fit tau_flow = tau0 + alpha*mu*b*sqrt(rho)
    rho = np.array([r["rho"] for r in rows])
    tf = np.array([r["flow_stress"] for r in rows])
    x = np.sqrt(rho)
    slope, tau0 = np.polyfit(x, tf, 1)
    alpha = slope / (mu * b)
    r2 = 1 - np.var(tf - (slope * x + tau0)) / (np.var(tf) + 1e-12)
    result = dict(mu=mu, b=b, alpha=float(alpha), tau0=float(tau0),
                  taylor_r2=float(r2), rows=rows,
                  note="alpha = Taylor forest-hardening coefficient "
                       "(dimensionless); literature metals ~0.3-0.5")
    with open(os.path.join(OUT, "taylor.json"), "w") as f:
        json.dump(result, f, indent=1)
    print(f"\nTAYLOR FIT: tau_flow = {tau0:.5f} + alpha*mu*b*sqrt(rho), "
          f"alpha = {alpha:.3f}, R^2 = {r2:.3f}  (metals: alpha~0.3-0.5)",
          flush=True)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(7, 5.5))
        plt.plot(np.sqrt(rho), tf, "o", label="PFC flow stress")
        xx = np.linspace(0, np.sqrt(rho).max(), 20)
        plt.plot(xx, slope * xx + tau0, "--",
                 label=f"Taylor fit: alpha={alpha:.2f}, R2={r2:.2f}")
        plt.xlabel(r"$\sqrt{\rho}$  (PFC units)")
        plt.ylabel(r"flow stress $\tau$ (PFC units)")
        plt.title("Taylor forest hardening from PFC: tau ~ alpha*mu*b*sqrt(rho)")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.savefig(os.path.join(OUT, "taylor.png"), dpi=140,
                    bbox_inches="tight")
    except Exception as ex:
        print("plot skipped:", ex)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

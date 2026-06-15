"""Scale dependence of the Taylor coefficient: small (Lbox=10000 b) vs large
(Lbox=16000 b) multi-slip forest-probe series, fixed 2 probes.

Finding: the DIFFERENTIAL forest-hardening slope d(tau)/d(mu b sqrt(rho)) drops
from ~1.2 (small, strong-obstacle/Orowan limit) to ~0.5 (large, BULK-Cu range)
as the forest grows and all 12 slip systems are better represented. The large
box also carries a density-independent offset (longer probes / fewer carriers
under strain-rate control), so a naive through-origin fit looks worse even
though the slope is more bulk-like.

  python3 plot_scale_comparison.py
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
MU, B = 54.6e9, 2.55e-10


def load(fp):
    P = json.load(open(fp))["points"]
    P.sort(key=lambda p: p["rho_forest"])
    rho = np.array([p["rho_forest"] for p in P])
    tau = np.array([p["flow_stress"] for p in P])
    return rho, tau


def fits(rho, tau):
    x = MU * B * np.sqrt(rho)
    a0 = np.sum(x * tau) / np.sum(x * x)            # through origin
    A = np.vstack([x, np.ones_like(x)]).T
    (a1, c1), _, _, _ = np.linalg.lstsq(A, tau, rcond=None)   # with intercept
    r2 = lambda p: 1 - np.sum((tau - p) ** 2) / np.sum((tau - tau.mean()) ** 2)
    return a0, r2(a0 * x), a1, c1, r2(a1 * x + c1)


small = load(os.path.join(HERE, "forest_taylor_fixed.json"))
big = load(os.path.join(HERE, "bulk_taylor.json"))

fig, ax = plt.subplots(figsize=(8, 6.2))
for (rho, tau), name, col in [(small, "small box (L=10000 b)", "tab:red"),
                              (big, "large multi-slip box (L=16000 b)", "tab:purple")]:
    a0, r2_0, a1, c1, r2_1 = fits(rho, tau)
    x = MU * B * np.sqrt(rho) / 1e6
    ax.scatter(x, tau / 1e6, s=95, color=col, zorder=6,
               label=f"{name}")
    xx = np.linspace(0, x.max() * 1.1, 50)
    ax.plot(xx, (a1 * (xx * 1e6) + c1) / 1e6, "--", color=col, alpha=0.9,
            label=fr"  fit $\tau$={c1/1e6:.0f}+{a1:.2f}$\mu b\sqrt{{\rho}}$ (R$^2$={r2_1:.2f})")
    print(f"{name}: through-origin alpha={a0:.2f} (R2={r2_0:.2f}); "
          f"slope alpha={a1:.2f}, offset={c1/1e6:.1f} MPa (R2={r2_1:.2f})")

xx = np.linspace(0, (MU * B * np.sqrt(big[0].max()) / 1e6) * 1.1, 50)
ax.fill_between(xx, 0.3 * xx, 0.5 * xx, color="tab:green", alpha=0.15,
                label=r"Cu bulk $\alpha$=0.3-0.5 (Taylor slope)")
ax.set_xlabel(r"$\mu\, b\, \sqrt{\rho_f}$ (MPa)")
ax.set_ylabel(r"forest-probe flow stress $\tau_c$ (MPa)")
ax.set_title("Taylor slope moves toward bulk Cu as the multi-slip forest grows")
ax.legend(fontsize=8.2, loc="upper left")
ax.grid(alpha=0.3)
ax.set_xlim(0, xx.max())
ax.set_ylim(0, 115)
fig.tight_layout()
fp = os.path.join(HERE, "taylor_scale_comparison.png")
fig.savefig(fp, dpi=130, bbox_inches="tight")
print("saved", fp)

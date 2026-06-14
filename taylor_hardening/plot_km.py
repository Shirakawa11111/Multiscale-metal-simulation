"""Kocks-Mecking dynamic steady-state convergence: regardless of INITIAL
density, all DDD configs converge to a common (rho_ss, sigma_ss) set by the
strain rate (multiplication-annihilation balance). This is the hallmark of
dislocation plasticity that PFC categorically lacks (PFC softens monotonically).

  python3 plot_km.py /tmp/taylor_curves.json
"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/taylor_curves.json"
OUT = os.path.dirname(os.path.abspath(__file__))
c = json.load(open(SRC))
MU, B = 54.6e9, 2.55e-10

fig, ax = plt.subplots(1, 3, figsize=(16, 5))
colors = plt.cm.viridis(np.linspace(0, 0.85, len(c)))
rss, sss = [], []
for col, nl in zip(colors, sorted(c, key=int)):
    r = np.array(c[nl])
    strain, stress, rho = r[:, 1] * 100, r[:, 2] / 1e6, r[:, 3]
    ax[0].plot(strain, rho, "-", color=col, lw=1.6,
               label=fr"$\rho_0$={rho[0]:.1e}")
    ax[1].plot(strain, stress, "-", color=col, lw=1.6,
               label=fr"$\rho_0$={rho[0]:.1e}")
    m = r[:, 1] >= 0.7 * r[:, 1].max()
    rss.append(np.mean(rho[m]))
    sss.append(np.mean(stress[m]))

ax[0].axhline(np.mean(rss), color="k", ls=":", lw=1,
              label=fr"$\rho_{{ss}}\approx${np.mean(rss):.1e}")
ax[0].set_xlabel("strain (%)")
ax[0].set_ylabel(r"dislocation density $\rho$ (m$^{-2}$)")
ax[0].set_yscale("log")
ax[0].set_title("Density converges to a common steady state\n(dense annihilates DOWN to $\\rho_{ss}$)")
ax[0].legend(fontsize=8)
ax[0].grid(alpha=0.3, which="both")

ax[1].set_xlabel("strain (%)")
ax[1].set_ylabel(r"flow stress $\sigma$ (MPa)")
ax[1].set_title("Flow stress converges to a common steady state\n(Kocks-Mecking dynamic recovery)")
ax[1].legend(fontsize=8)
ax[1].grid(alpha=0.3)

# steady-state point(s) vs Taylor line
rss, sss = np.array(rss), np.array(sss)
x = MU * B * np.sqrt(rss) / 1e6
ax[2].scatter(x, sss, s=90, color="tab:red", zorder=5,
              label=r"DDD steady states ($\rho_{ss},\sigma_{ss}$)")
alpha = float(np.sum(x * sss) / np.sum(x * x))
xx = np.linspace(0, x.max() * 1.3, 50)
ax[2].plot(xx, alpha * xx, "k--", label=fr"$\alpha$={alpha:.2f} (this regime)")
ax[2].fill_between(xx, 0.3 * xx, 0.5 * xx, color="tab:green", alpha=0.15,
                   label=r"Cu bulk $\alpha$=0.3-0.5")
ax[2].set_xlabel(r"$\mu\, b\, \sqrt{\rho_{ss}}$ (MPa)")
ax[2].set_ylabel(r"$\sigma_{ss}$ (MPa)")
ax[2].set_title(r"Steady states cluster at one $\rho_{ss}$" "\n(no $\\sqrt{\\rho}$ lever from initial density)")
ax[2].legend(fontsize=8)
ax[2].grid(alpha=0.3)
ax[2].set_xlim(0, x.max() * 1.3)
ax[2].set_ylim(0, max(sss.max(), 0.5 * x.max()) * 1.2)

fig.tight_layout()
fp = os.path.join(OUT, "km_steady_state.png")
fig.savefig(fp, dpi=130, bbox_inches="tight")
print(f"rho_ss mean = {np.mean(rss):.2e} (range {rss.min():.2e}-{rss.max():.2e})")
print(f"sigma_ss mean = {np.mean(sss):.1f} MPa")
print(f"alpha (steady-state cluster) = {alpha:.2f}")
print("saved", fp)

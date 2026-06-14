"""Thematic closure: does the REAL STEM-reconstructed Cu network obey the same
DDD Taylor forest-hardening law as the synthetic forest-probe series?

Overlays the reconstructed config's free-evolution (sigma, rho) trajectory (from
experiment_bridge/results_exadis/stress_strain_dens.dat) on the forest-probe
Taylor line tau_c = alpha mu b sqrt(rho) (alpha~1.4, R2=0.87).

Honest caveat: the reconstructed run is a FREE evolution (pinned foil anchors,
multiplication 3.2e12->7.7e12), not a controlled fixed-forest depinning test, so
this is a consistency check (does it sit in the same Taylor regime?), not an
exact alpha measurement.

  python3 overlay_stem_taylor.py
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
MU, B = 54.6e9, 2.55e-10

fj = json.load(open(os.path.join(HERE, "forest_taylor_fixed.json")))
rho_f = np.array([p["rho_forest"] for p in fj["points"]])
tau_f = np.array([p["flow_stress"] for p in fj["points"]])
alpha = fj["taylor_alpha"]

stem = np.loadtxt(os.path.join(HERE, "..", "experiment_bridge", "results_exadis",
                                "stress_strain_dens.dat"), ndmin=2)
s_strain, s_stress, s_rho = stem[:, 1], stem[:, 2], stem[:, 3]
# plastic-flow portion (past the elastic loading / first yield)
pf = s_strain >= 3e-4
xs = MU * B * np.sqrt(s_rho) / 1e6
# representative reconstructed point: mean over plastic flow
rho_rec = float(np.mean(s_rho[pf]))
tau_rec = float(np.mean(np.abs(s_stress[pf])))
alpha_rec = tau_rec / (MU * B * np.sqrt(rho_rec))

fig, ax = plt.subplots(figsize=(7.5, 6))
xf = MU * B * np.sqrt(rho_f) / 1e6
ax.scatter(xf, tau_f / 1e6, s=110, color="tab:red", zorder=6,
           label=fr"synthetic forest-probe ($\alpha$={alpha:.2f}, R$^2$={fj['r2']:.2f})")
xx = np.linspace(0, max(xf.max(), xs[pf].max()) * 1.1, 50)
ax.plot(xx, alpha * xx, "k--", zorder=4, label=r"Taylor fit $\tau=\alpha\mu b\sqrt{\rho}$")
ax.fill_between(xx, 0.3 * xx, 0.5 * xx, color="tab:green", alpha=0.15,
                label=r"Cu bulk $\alpha$=0.3-0.5")
# reconstructed-config trajectory (plastic flow)
ax.scatter(xs[pf], np.abs(s_stress[pf]) / 1e6, s=28, color="tab:blue", alpha=0.55,
           zorder=5, label="real STEM-reconstructed config (flow trajectory)")
ax.scatter([MU * B * np.sqrt(rho_rec) / 1e6], [tau_rec / 1e6], marker="*",
           s=420, color="navy", edgecolor="white", zorder=7,
           label=fr"reconstructed mean ($\alpha_{{eff}}$={alpha_rec:.2f})")
ax.set_xlabel(r"$\mu\, b\, \sqrt{\rho}$ (MPa)")
ax.set_ylabel(r"flow stress $\tau$ (MPa)")
ax.set_title("Real reconstructed Cu network obeys the same DDD Taylor regime")
ax.legend(fontsize=8.5, loc="upper left")
ax.grid(alpha=0.3)
ax.set_xlim(0, xx.max())
ax.set_ylim(0, max(tau_f.max() / 1e6, (np.abs(s_stress[pf]) / 1e6).max()) * 1.25)
fig.tight_layout()
fp = os.path.join(HERE, "stem_on_taylor_line.png")
fig.savefig(fp, dpi=130, bbox_inches="tight")
print(f"reconstructed config: rho_mean={rho_rec:.2e}, tau_mean={tau_rec/1e6:.1f} MPa, "
      f"alpha_eff={alpha_rec:.2f}  (forest-probe alpha={alpha:.2f})")
print("saved", fp)

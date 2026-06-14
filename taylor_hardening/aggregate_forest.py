"""Aggregate forest-probe results into the Taylor fit + figure.

Reads OUT/n*/result.json (one per pinned-forest density), fits
tau_c = alpha * mu * b * sqrt(rho_forest) through the origin, and plots the
per-density stress-strain curves + the Taylor line vs the Cu literature band.

  python3 aggregate_forest.py forest_out
"""
import os, sys, glob, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else "forest_out"
MU, B = 54.6e9, 2.55e-10

R = []
for fp in sorted(glob.glob(os.path.join(OUT, "n*", "result.json"))):
    R.append(json.load(open(fp)))
R.sort(key=lambda r: r["rho_forest"])
if not R:
    print("no results found in", OUT)
    sys.exit(1)

rho = np.array([r["rho_forest"] for r in R])
sig = np.array([r["flow_stress"] for r in R])
err = np.array([r.get("flow_std", 0.0) for r in R])
x = MU * B * np.sqrt(rho)
alpha = float(np.sum(x * sig) / np.sum(x * x))
if len(sig) > 1:
    pred = alpha * x
    ss_res = float(np.sum((sig - pred) ** 2))
    ss_tot = float(np.sum((sig - sig.mean()) ** 2))
    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else 0.0
else:
    r2 = 0.0

summary = dict(taylor_alpha=alpha, r2=r2, n_points=len(R), mu=MU, b=B,
               points=[dict(rho_forest=r["rho_forest"],
                            rho_flow=r["rho_flow"],
                            flow_stress=r["flow_stress"],
                            num_lines=r["num_lines"], k_probe=r["k_probe"])
                       for r in R],
               note="Forest-probe Taylor: pinned forest at controlled rho_f + "
                    "mobile probes; flow stress = tau_c = alpha mu b sqrt(rho_f).")
with open(os.path.join(OUT, "forest_taylor.json"), "w") as f:
    json.dump(summary, f, indent=1)

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
colors = plt.cm.plasma(np.linspace(0, 0.85, len(R)))
for col, r in zip(colors, R):
    c = np.array(r["curve"])
    ax[0].plot(c[:, 1] * 100, c[:, 2] / 1e6, "-", color=col, lw=1.5,
               label=fr"$\rho_f$={r['rho_forest']:.1e} (N={r['num_lines']})")
ax[0].set_xlabel("strain (%)")
ax[0].set_ylabel("stress (MPa)")
ax[0].set_title("Forest-probe stress-strain (pinned forest + mobile probes)")
ax[0].legend(fontsize=8)
ax[0].grid(alpha=0.3)

xp = x / 1e6
ax[1].errorbar(xp, sig / 1e6, yerr=err / 1e6, fmt="o", ms=9, capsize=4,
               color="tab:red", zorder=5, label=r"$\tau_c$ (forest-probe)")
xx = np.linspace(0, xp.max() * 1.15, 50)
ax[1].plot(xx, alpha * xx, "k--",
           label=fr"Taylor fit $\alpha$={alpha:.2f} (R$^2$={r2:.2f})")
ax[1].fill_between(xx, 0.3 * xx, 0.5 * xx, color="tab:green", alpha=0.15,
                   label=r"Cu bulk $\alpha$=0.3-0.5")
ax[1].set_xlabel(r"$\mu\, b\, \sqrt{\rho_f}$ (MPa)")
ax[1].set_ylabel(r"flow stress $\tau_c$ (MPa)")
ax[1].set_title(r"Taylor law $\tau_c=\alpha\,\mu\,b\,\sqrt{\rho_f}$ (pinned forest)")
ax[1].legend(fontsize=8)
ax[1].grid(alpha=0.3)
ax[1].set_xlim(0, xp.max() * 1.15)
ax[1].set_ylim(0, max((sig / 1e6).max(), alpha * xx.max()) * 1.2)

fig.tight_layout()
fp = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                  "forest_taylor.png")
fig.savefig(fp, dpi=130, bbox_inches="tight")
print(f"n_points={len(R)}")
for r in R:
    print(f"  rho_f={r['rho_forest']:.3e}  tau_c={r['flow_stress']/1e6:.1f} MPa"
          f"  (N={r['num_lines']}, k_probe={r['k_probe']})")
print(f"=== FOREST-PROBE TAYLOR alpha={alpha:.3f}, R^2={r2:.3f} ===")
print("saved", fp)

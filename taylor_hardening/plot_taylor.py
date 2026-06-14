"""Plot the Taylor hardening density series: flow stress vs sqrt(rho), with the
per-run stress-strain curves, and the fitted Taylor coefficient alpha.

  python3 plot_taylor.py taylor_out
"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else "taylor_out"
S = json.load(open(os.path.join(OUT, "taylor_series.json")))
curves = json.load(open(os.path.join(OUT, "taylor_curves.json")))
mu, b = S["mu"], S["b"]
series = S["series"]

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))

# (1) stress-strain curves per density
for s in series:
    res = np.array(curves[str(s["num_lines"])])
    ax[0].plot(res[:, 1] * 100, res[:, 2] / 1e6, "-",
               label=fr"$\rho_0$={s['rho0']:.1e} m$^{{-2}}$")
ax[0].set_xlabel("strain (%)")
ax[0].set_ylabel("flow stress (MPa)")
ax[0].set_title("FCC Cu DDD (junctions on): stress-strain vs initial density")
ax[0].legend(fontsize=8)
ax[0].grid(alpha=0.3)

# (2) Taylor: sigma_flow vs mu b sqrt(rho)
rr = np.array([s["rho_flow"] for s in series])
ss = np.array([s["flow_stress"] for s in series]) / 1e6
err = np.array([s.get("flow_std", 0.0) for s in series]) / 1e6
x = mu * b * np.sqrt(rr) / 1e6
ax[1].errorbar(x, ss, yerr=err, fmt="o", ms=9, capsize=4, color="tab:red",
               label="DDD flow stress")
alpha = S["taylor_alpha"]
xx = np.linspace(0, x.max() * 1.1, 50)
ax[1].plot(xx, alpha * xx, "k--",
           label=fr"Taylor fit $\alpha$={alpha:.3f} (R$^2$={S['r2']:.2f})")
# literature band alpha 0.3-0.5
ax[1].fill_between(xx, 0.3 * xx, 0.5 * xx, color="tab:green", alpha=0.15,
                   label=r"Cu literature $\alpha$=0.3-0.5")
ax[1].set_xlabel(r"$\mu\, b\, \sqrt{\rho}$ (MPa)")
ax[1].set_ylabel(r"flow stress $\sigma$ (MPa)")
ax[1].set_title(r"Taylor law $\sigma=\alpha\,\mu\,b\,\sqrt{\rho}$ from DDD")
ax[1].legend(fontsize=8)
ax[1].grid(alpha=0.3)

fig.tight_layout()
fp = os.path.join(OUT, "taylor_hardening.png")
fig.savefig(fp, dpi=130, bbox_inches="tight")
print("alpha =", alpha, "R2 =", S["r2"])
print("saved", fp)

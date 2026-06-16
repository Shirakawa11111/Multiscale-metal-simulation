"""Multi-seed robustness of the forest-probe Taylor coefficient.

Three independent random seeds x four forest densities (large box, line-length
rho_f), baseline-subtracted with the measured no-forest tau0 = 72.1 MPa. The
result is an honest negative on robustness: the coefficient has large
seed-to-seed scatter, so alpha is bracketed, not determined, at this system size.

  python3 plot_multiseed.py
"""
import os, json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MU, B, TAU0 = 54.6e9, 2.55e-10, 72.1e6
rows = [json.load(open(os.path.join(d, "result.json")))
        for d in sorted(glob.glob("/tmp/ms_out/n*"))
        if os.path.exists(os.path.join(d, "result.json"))]
seeds = sorted(set(r["seed"] for r in rows))

fig, ax = plt.subplots(figsize=(8, 6.2))
cols = {1234: "tab:blue", 5678: "tab:green", 9012: "tab:red"}
alphas = []
for s in seeds:
    forest = sorted([r for r in rows if r["seed"] == s and r["rho_forest"] > 0],
                    key=lambda r: r["rho_forest"])
    rho = np.array([r["rho_forest"] for r in forest])
    tau = np.array([r["flow_stress"] for r in forest])
    std = np.array([r["flow_std"] for r in forest])
    x = MU * B * np.sqrt(rho) / 1e6
    y = (tau - TAU0) / 1e6
    a = float(np.sum((x * 1e6) * (y * 1e6)) / np.sum((x * 1e6) ** 2))
    alphas.append(a)
    noisy = std / tau > 0.1
    ax.errorbar(x, y, yerr=std / 1e6, fmt="o", ms=8, color=cols[s], capsize=3,
                label=fr"seed {s}: $\alpha$={a:.2f}" + (" (noisy)" if noisy.sum() >= 2 else ""))
    xx = np.linspace(0, x.max() * 1.1, 50)
    ax.plot(xx, a * xx, "--", color=cols[s], alpha=0.6, lw=1)

alphas = np.array(alphas)
xx = np.linspace(0, 75, 50)
ax.fill_between(xx, 0.3 * xx, 0.5 * xx, color="0.6", alpha=0.25,
                label=r"Cu bulk $\alpha$=0.3-0.5")
ax.plot(xx, alphas.mean() * xx, "k-", lw=2,
        label=fr"mean $\alpha$={alphas.mean():.2f}$\pm${alphas.std():.2f} (n=3 seeds)")
ax.set_xlabel(r"$\mu\, b\, \sqrt{\rho_f}$ (MPa)")
ax.set_ylabel(r"$\tau_c-\tau_0$ (measured baseline removed, MPa)")
ax.set_title("Forest-probe Taylor coefficient: large seed-to-seed scatter\n"
             fr"$\alpha$ = {alphas.mean():.2f} $\pm$ {alphas.std():.2f} "
             fr"(range {alphas.min():.2f}-{alphas.max():.2f}); not robust at n=4/seed")
ax.legend(fontsize=8.5, loc="upper left")
ax.grid(alpha=0.3)
ax.set_xlim(0, 75)
fig.tight_layout()
fp = os.path.dirname(os.path.abspath(__file__)) + "/multiseed_taylor.png"
fig.savefig(fp, dpi=130, bbox_inches="tight")
print(f"alpha per seed: {[round(a,3) for a in alphas]}")
print(f"mean={alphas.mean():.3f} SD={alphas.std():.3f}")
print("saved", fp)

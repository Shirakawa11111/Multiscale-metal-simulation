"""Directly-measured carrier baseline -> baseline-subtracted forest hardening.

Peer review (correctly) flagged that the large-box through-origin alpha=1.70 is a
rejected fit and that the offset interpretation was asserted, not measured. So we
MEASURED it: a no-forest run (NUM_LINES=2, N_PROBE=2, rho_forest=0) gives the pure
2-probe carrier/drag baseline tau0:
    small box (L=10000): tau0 = 32.8 MPa
    large box (L=16000): tau0 = 72.1 MPa
The baseline rises strongly with box size (same 2 probes, larger cell -> higher
probe velocity to carry the imposed rate -> higher drag), confirming the
carrier-baseline interpretation directly rather than by assertion.

Subtracting the MEASURED tau0 from each forest point isolates the forest
contribution; fitting (tau - tau0) = alpha * mu b sqrt(rho_f) through the origin
then gives the forest-hardening coefficient WITHOUT the carrier offset.

  python3 baseline_subtracted.py
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
MU, B = 54.6e9, 2.55e-10
TAU0 = {"small": 32.8e6, "large": 72.1e6}   # measured no-forest 2-probe baselines


def load(fp):
    P = json.load(open(fp))["points"]
    P.sort(key=lambda p: p["rho_forest"])
    return (np.array([p["rho_forest"] for p in P]),
            np.array([p["flow_stress"] for p in P]))


def fit0(x, y):
    a = float(np.sum(x * y) / np.sum(x * x))
    r2 = 1 - np.sum((y - a * x) ** 2) / np.sum((y - y.mean()) ** 2)
    return a, r2


small = load(os.path.join(HERE, "forest_taylor_fixed.json"))
big = load(os.path.join(HERE, "bulk_taylor.json"))

fig, ax = plt.subplots(figsize=(8, 6.2))
for (rho, tau), key, col in [(small, "small", "tab:red"), (big, "large", "tab:purple")]:
    x = MU * B * np.sqrt(rho)
    yr = tau - TAU0[key]                      # baseline-subtracted forest stress
    a, r2 = fit0(x, yr)
    perpt = yr / x
    lbl = "small box (L=10000 b)" if key == "small" else "large multi-slip box (L=16000 b)"
    ax.scatter(x / 1e6, yr / 1e6, s=95, color=col, zorder=6,
               label=fr"{lbl}: $\alpha$={a:.2f} (R$^2$={r2:.2f})")
    xx = np.linspace(0, x.max() / 1e6 * 1.1, 50)
    ax.plot(xx, a * xx, "--", color=col, alpha=0.85)
    print(f"{key} box: measured tau0={TAU0[key]/1e6:.1f} MPa; "
          f"baseline-subtracted through-origin alpha={a:.3f} (R2={r2:.2f}); "
          f"per-point alpha={[round(v,2) for v in perpt]}")

xx = np.linspace(0, (MU * B * np.sqrt(big[0].max()) / 1e6) * 1.1, 50)
ax.fill_between(xx, 0.3 * xx, 0.5 * xx, color="tab:green", alpha=0.18,
                label=r"Cu bulk $\alpha$=0.3-0.5")
ax.set_xlabel(r"$\mu\, b\, \sqrt{\rho_f}$ (MPa)")
ax.set_ylabel(r"forest stress $\tau_c-\tau_0$ (measured baseline removed, MPa)")
ax.set_title("Forest hardening with the MEASURED carrier baseline removed\n"
             "(large box lands in the bulk-Cu band)")
ax.legend(fontsize=8.5, loc="upper left")
ax.grid(alpha=0.3)
ax.set_xlim(0, xx.max())
ax.set_ylim(-5, None)
fig.tight_layout()
fp = os.path.join(HERE, "baseline_subtracted.png")
fig.savefig(fp, dpi=130, bbox_inches="tight")
print("saved", fp)

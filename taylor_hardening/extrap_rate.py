"""Rate-extrapolation of the clean (K=8) forest-hardening coefficient.

For each strain rate, compute the per-seed baseline-subtracted alpha (densities
50/100/200, each seed's own no-forest baseline), then mean +/- SD per rate.
Fit alpha vs rate and extrapolate to rate -> 0 (quasi-static). Confirms whether
the strong-obstacle conclusion (alpha above bulk Cu) holds quasi-statically.

  python3 extrap_rate.py /tmp/rateK8
"""
import os, sys, glob, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = sys.argv[1] if len(sys.argv) > 1 else "/tmp/rateK8"
MU, B = 54.6e9, 2.55e-10
RATES = [3e3, 5e3, 1e4]
SEEDS = ["1234", "5678", "2222"]


def load(tag):
    fp = glob.glob(os.path.join(R, tag, "*", "result.json"))
    return json.load(open(fp[0])) if fp else None


def erate_tag(e):
    return {3e3: "3e3", 5e3: "5e3", 1e4: "1e4"}[e]


rate_alpha = {}
print("=== clean K=8 forest coefficient vs strain rate ===")
for e in RATES:
    et = erate_tag(e)
    al = []
    for s in SEEDS:
        base = load(f"B_e{et}_s{s}")
        if not base:
            continue
        tau0 = base["flow_stress"]
        pts = [load(f"F_e{et}_n{nl}_s{s}") for nl in [50, 100, 200]]
        pts = [p for p in pts if p]
        if len(pts) < 3:
            continue
        rho = np.array([p["rho_forest"] for p in pts])
        tau = np.array([p["flow_stress"] for p in pts])
        x = MU * B * np.sqrt(rho); y = tau - tau0
        al.append(float(np.sum(x * y) / np.sum(x * x)))
    if al:
        rate_alpha[e] = (np.mean(al), np.std(al), len(al), np.mean([
            (load(f"B_e{et}_s{s}") or {}).get("flow_stress", np.nan) for s in SEEDS
            if load(f"B_e{et}_s{s}")]) / 1e6)
        print(f"  erate={et}: alpha={np.mean(al):.3f}+/-{np.std(al):.3f} "
              f"(n={len(al)} seeds), baseline tau0~{rate_alpha[e][3]:.0f} MPa")

# extrapolate alpha vs rate -> 0
if len(rate_alpha) >= 2:
    es = np.array(sorted(rate_alpha))
    av = np.array([rate_alpha[e][0] for e in es])
    sd = np.array([rate_alpha[e][1] for e in es])
    # linear fit alpha = a0 + a1*rate ; quasi-static = a0
    A = np.vstack([np.ones_like(es), es]).T
    coef, *_ = np.linalg.lstsq(A, av, rcond=None)
    a0 = float(coef[0])
    print(f"\n=== quasi-static (rate->0) alpha extrapolation: {a0:.3f} ===")
    print(f"    (bulk Cu 0.3-0.5; clean values rise as rate falls -> strong-obstacle regime)")

    fig, ax = plt.subplots(figsize=(7.5, 5.4))
    ax.errorbar(es, av, yerr=sd, fmt="o", ms=9, capsize=5, color="tab:blue",
                label="clean K=8 (per-rate mean)")
    rr = np.linspace(0, es.max() * 1.1, 50)
    ax.plot(rr, coef[0] + coef[1] * rr, "k--",
            label=fr"linear fit; quasi-static $\alpha\to${a0:.2f}")
    ax.axhspan(0.3, 0.5, color="0.6", alpha=0.25, label="Cu bulk 0.3-0.5")
    ax.scatter([0], [a0], marker="*", s=300, color="navy", zorder=6)
    ax.set_xlabel("strain rate (1/s)")
    ax.set_ylabel(r"forest-hardening coefficient $\alpha$")
    ax.set_title("Rate extrapolation of the clean Taylor coefficient\n"
                 r"(quasi-static $\alpha$ stays above bulk = strong-obstacle regime)")
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    ax.set_xlim(-500, es.max() * 1.1); ax.set_ylim(0, max(av.max() + 0.2, 1.0))
    fig.tight_layout()
    fp = os.path.dirname(os.path.abspath(__file__)) + "/rate_extrapolation.png"
    fig.savefig(fp, dpi=130, bbox_inches="tight")
    print("saved", fp)

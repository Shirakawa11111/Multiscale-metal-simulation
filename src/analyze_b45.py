"""B4/B5 combined analysis: strain-rate scan (3-point, with b3 as the
RELAX=400 midpoint) + nucleation-case overview.
Output: results/b45_series_512/analysis_b45.png + metrics_b45.json
"""

import os, json
import numpy as np

BASE = os.path.join(os.path.dirname(__file__), "..", "results")
RATE_CASES = [  # (label, summary path, relax steps)
    ("fast (RELAX=100)", "b45_series_512/b5_rate_fast/summary.json", 100),
    ("base (RELAX=400)", "b_series_512/b3_polycrystal/summary.json", 400),
    ("slow (RELAX=1600)", "b45_series_512/b5_rate_slow/summary.json", 1600),
]
NUC_CASES = [
    ("pore v1 (dissolved)", "b45_series_512/b4_pore/summary.json"),
    ("noise 0.02", "b45_series_512/b4_noise/summary.json"),
]
EPS_CMP = 0.085  # common strain for rate comparison


def load(rel):
    with open(os.path.join(BASE, rel)) as f:
        rows = json.load(f)["rows"]
    eps = np.array([r["exx"] for r in rows])
    sig = np.array([r["sigma"] if r["sigma"] is not None else np.nan
                    for r in rows])
    cores = np.array([r["cores"] for r in rows])
    return eps, sig, cores


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    metrics = {}
    for label, path, relax in RATE_CASES:
        eps, sig, cores = load(path)
        axes[0].plot(eps * 100, sig, "o-", ms=3, label=label)
        axes[1].plot(eps * 100, cores, "s-", ms=3, label=label)
        i = int(np.nanargmin(np.abs(eps - EPS_CMP)))
        metrics[label] = dict(relax=relax,
                              flow_stress_at_8p5=float(sig[i]),
                              cores_at_8p5=int(cores[i]),
                              cores_initial=int(cores[0]))
    axes[0].set_xlabel("strain (%)")
    axes[0].set_ylabel("stress (PFC units)")
    axes[0].set_title("polycrystal stress-strain vs strain rate")
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    axes[1].set_xlabel("strain (%)")
    axes[1].set_ylabel("# dislocation cores")
    axes[1].set_title("dislocation retention vs strain rate")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    # rate-sensitivity summary: flow stress & cores at common strain vs rate
    rr = np.array([[1.0 / m["relax"], m["flow_stress_at_8p5"],
                    m["cores_at_8p5"]] for m in metrics.values()])
    ax2 = axes[2]
    ax2.semilogx(rr[:, 0], rr[:, 1], "o-", color="tab:red",
                 label="flow stress @8.5%")
    ax2.set_xlabel("effective strain rate ~ 1/RELAX")
    ax2.set_ylabel("flow stress", color="tab:red")
    ax2b = ax2.twinx()
    ax2b.semilogx(rr[:, 0], rr[:, 2], "s--", color="tab:blue",
                  label="cores @8.5%")
    ax2b.set_ylabel("# cores", color="tab:blue")
    ax2.set_title("rate sensitivity at 8.5% strain")
    ax2.grid(alpha=0.3)

    # log-log slope: sigma ~ rate^m (viscoplastic rate exponent)
    msens = np.polyfit(np.log(rr[:, 0]), np.log(rr[:, 1]), 1)[0]
    metrics["rate_exponent_m"] = float(msens)
    ax2.text(0.05, 0.92, f"m = dln(sigma)/dln(rate) = {msens:.3f}",
             transform=ax2.transAxes)

    for label, path in NUC_CASES:
        eps, sig, cores = load(path)
        metrics[label] = dict(final_cores=int(cores[-1]),
                              max_cores=int(cores.max()))

    out = os.path.join(BASE, "b45_series_512")
    fig.savefig(os.path.join(out, "analysis_b45.png"), dpi=140,
                bbox_inches="tight")
    with open(os.path.join(out, "metrics_b45.json"), "w") as f:
        json.dump(metrics, f, indent=1)
    print(json.dumps(metrics, indent=1))


if __name__ == "__main__":
    main()

"""Analysis of B-series v2: stress-strain curves, dislocation activity,
mechanism extraction. Produces results/b_series_512/analysis.png + metrics
JSON (modulus, peak stress, toughness — same metric set as
Topology research/analyze_tensile_scan.py).
"""

import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np

BASE = os.path.join(os.path.dirname(__file__), "..", "results", "b_series_512")
CASES = ["b1b_void", "b2_dipole", "b3_polycrystal"]
LABELS = {"b1b_void": "void notch", "b2_dipole": "dipole",
          "b3_polycrystal": "polycrystal"}


def metrics(eps, sig, elastic_max=0.02):
    el = (eps > 0) & (eps <= elastic_max)
    E = float(np.polyfit(eps[el], sig[el], 1)[0]) if el.sum() >= 2 else np.nan
    ipk = int(np.argmax(sig))
    toughness = float(np.trapz(sig[:ipk + 1], eps[:ipk + 1]))
    return dict(modulus=E, peak_stress=float(sig[ipk]),
                peak_strain=float(eps[ipk]), toughness=toughness)


def main():
    data = {}
    out_metrics = {}
    for c in CASES:
        with open(os.path.join(BASE, c, "summary.json")) as f:
            rows = json.load(f)["rows"]
        eps = np.array([r["exx"] for r in rows])
        sig = np.array([r["sigma"] for r in rows])
        cores = np.array([r["cores"] for r in rows])
        rho = np.array([r["rho"] for r in rows])
        data[c] = (eps, sig, cores, rho)
        out_metrics[c] = metrics(eps, sig)
        # event extraction
        events = []
        for i in range(1, len(cores)):
            if cores[i] != cores[i - 1]:
                events.append(dict(exx=float(eps[i]),
                                   d_cores=int(cores[i] - cores[i - 1]),
                                   cores=int(cores[i])))
        out_metrics[c]["events"] = events[:40]
        out_metrics[c]["initial_cores"] = int(cores[0])
        out_metrics[c]["final_cores"] = int(cores[-1])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    ax = axes[0, 0]
    for c in CASES:
        eps, sig, *_ = data[c]
        ax.plot(eps * 100, sig, "o-", ms=3, label=LABELS[c])
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("stress (PFC units)")
    ax.set_title("stress-strain")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    for c in CASES:
        eps, _, cores, _ = data[c]
        ax.plot(eps * 100, cores, "s-", ms=3, label=LABELS[c])
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("# dislocation cores")
    ax.set_yscale("symlog")
    ax.set_title("dislocation count")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    for c in CASES:
        eps, _, _, rho = data[c]
        ax.plot(eps * 100, rho, "^-", ms=3, label=LABELS[c])
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("rho (cores / area)")
    ax.set_title("dislocation density evolution")
    ax.legend()
    ax.grid(alpha=0.3)

    # zoom: stress around events for b2 (annihilation) and b3 (plastic flow)
    ax = axes[1, 1]
    for c in ("b2_dipole", "b3_polycrystal"):
        eps, sig, cores, _ = data[c]
        dsig = np.gradient(sig, eps)
        ax.plot(eps * 100, dsig, "-", label=f"d(sigma)/d(eps) {LABELS[c]}")
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("tangent modulus")
    ax.set_title("tangent modulus (softening events)")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.savefig(os.path.join(BASE, "analysis.png"), dpi=140,
                bbox_inches="tight")
    with open(os.path.join(BASE, "metrics.json"), "w") as f:
        json.dump(out_metrics, f, indent=1)
    print(json.dumps(out_metrics, indent=1)[:3000])


if __name__ == "__main__":
    main()

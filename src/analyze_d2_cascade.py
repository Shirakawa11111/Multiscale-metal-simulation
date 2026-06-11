"""D2 full-range analysis (part 1 + part 2): yield -> avalanche -> flow.
Outputs: results/d2_roi_matched_1536/analysis_cascade.png + cascade.json
"""

import os, json
import numpy as np

OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   "d2_roi_matched_1536")
SCALE = 2.556e-10 / (4 * np.pi / np.sqrt(3.0))  # m per PFC unit


def main():
    rows = []
    for fn in ("summary.json", "summary_part2.json"):
        with open(os.path.join(OUT, fn)) as f:
            rows += json.load(f)["rows"]
    eps = np.array([r["exx"] for r in rows])
    sig = np.array([r["sigma"] for r in rows])
    cores = np.array([r["cores"] for r in rows])
    L_m = 1536 * np.pi / 4 * SCALE
    rho = cores / L_m ** 2

    i_peak = int(np.argmax(sig))
    flow = float(np.mean(sig[-5:]))
    cascade = dict(
        sigma_peak=float(sig[i_peak]), eps_peak=float(eps[i_peak]),
        sigma_flow_late=flow, softening_ratio=flow / float(sig[i_peak]),
        cores_max=int(cores.max()), rho_max_m2=float(rho.max()),
        rho_initial_m2=float(rho[0]),
        multiplication_onset_eps=float(eps[np.argmax(cores > 4)]))
    with open(os.path.join(OUT, "cascade.json"), "w") as f:
        json.dump(cascade, f, indent=1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    ax = axes[0]
    ax.plot(eps * 100, sig, "b-", lw=1.5)
    ax.axvline(eps[i_peak] * 100, color="gray", ls=":")
    ax.annotate(f"yield {sig[i_peak]:.4f}\n@{eps[i_peak]*100:.1f}%",
                (eps[i_peak] * 100, sig[i_peak]),
                textcoords="offset points", xytext=(10, -30))
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("stress (PFC units)")
    ax.set_title("yield -> softening -> flow")
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.semilogy(eps * 100, np.maximum(cores, 0.5), "r-", lw=1.5)
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("# dislocation cores (log)")
    ax.set_title("multiplication cascade 4 -> ~190")
    ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(eps * 100, rho / 1e15, "m-", lw=1.5)
    ax.axhline(2.05, color="gray", ls="--", label="DAMASK ROI cell max")
    ax.set_xlabel("strain (%)")
    ax.set_ylabel(r"rho ($10^{15}$ m$^{-2}$)")
    ax.set_title("density: ROI-matched start -> saturation")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(OUT, "analysis_cascade.png"), dpi=140,
                bbox_inches="tight")
    print(json.dumps(cascade, indent=1))


if __name__ == "__main__":
    main()

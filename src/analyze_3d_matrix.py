"""3D BCC polycrystal matrix analysis (16 cells: r x rate x seed, 128^3).

Compares 3D mechanical response to the 2D maps: flow stress & modulus over
(r, strain-rate), seed-averaged; 3D rate exponent m_3D(r).
Outputs: results_hpc_3d/ANALYSIS/matrix3d.png, metrics3d.json
"""

import os, json
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..", "results_hpc_3d")
OUT = os.path.join(ROOT, "ANALYSIS")
RS = (-0.35, -0.30, -0.25, -0.20)
RELAXES = (100, 400)
SEEDS = (7, 11)
EPS_FLOW = 0.05


def load(name):
    with open(os.path.join(ROOT, name, "summary.json")) as f:
        s = json.load(f)
    rows = s["rows"]
    return (np.array([x["exx"] for x in rows]),
            np.array([x["sigma"] for x in rows]))


def main():
    os.makedirs(OUT, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    metrics = {}

    # modulus (early slope) and flow stress maps
    mod = np.full((len(RS), len(RELAXES)), np.nan)
    flow = np.full_like(mod, np.nan)
    for i, r in enumerate(RS):
        for j, rx in enumerate(RELAXES):
            ms, fs = [], []
            for s in SEEDS:
                try:
                    eps, sig = load(f"p3d_r{r}_x{rx}_s{s}")
                except FileNotFoundError:
                    continue
                el = (eps > 0) & (eps <= 0.015)
                if el.sum() >= 2:
                    ms.append(np.polyfit(eps[el], sig[el], 1)[0])
                fs.append(sig[int(np.argmin(np.abs(eps - EPS_FLOW)))])
            if ms:
                mod[i, j] = np.mean(ms)
            if fs:
                flow[i, j] = np.mean(fs)
    metrics["modulus_3d"] = mod.tolist()
    metrics["flow_stress_3d_at_5pct"] = flow.tolist()

    ax = axes[0]
    for j, rx in enumerate(RELAXES):
        ax.plot(RS, flow[:, j], "o-", label=f"RELAX={rx}")
    ax.set_xlabel("r (quench depth)")
    ax.set_ylabel(f"flow stress @ {EPS_FLOW*100:.0f}%")
    ax.set_title("3D polycrystal flow stress vs r")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1]
    rates = 1.0 / np.array(RELAXES)
    ms_r = []
    for i, r in enumerate(RS):
        if np.isfinite(flow[i]).all() and (flow[i] > 0).all():
            m = np.polyfit(np.log(rates), np.log(flow[i]), 1)[0]
        else:
            m = np.nan
        ms_r.append(m)
        ax.plot(rates, flow[i], "s-", label=f"r={r} (m={m:.2f})")
    metrics["rate_exponent_3d_vs_r"] = dict(zip(map(str, RS), ms_r))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("strain rate ~ 1/RELAX")
    ax.set_ylabel("flow stress")
    ax.set_title("3D rate sensitivity (2-point m)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # representative stress-strain curves
    ax = axes[2]
    for r in RS:
        try:
            eps, sig = load(f"p3d_r{r}_x400_s7")
            ax.plot(eps * 100, sig, "-", label=f"r={r}")
        except FileNotFoundError:
            pass
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("stress")
    ax.set_title("3D stress-strain (RELAX=400, seed 7)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.savefig(os.path.join(OUT, "matrix3d.png"), dpi=140,
                bbox_inches="tight")
    with open(os.path.join(OUT, "metrics3d.json"), "w") as f:
        json.dump(metrics, f, indent=1)
    print(json.dumps(metrics, indent=1))


if __name__ == "__main__":
    main()

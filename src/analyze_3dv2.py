"""3D v2 analysis: (1) GB-distance decomposition of the Σ5 bicrystal size
series — proves disorder is GB-localized (no intragranular emission);
(2) in-window (<=8%) poly3d mechanism maps with CSP.
Outputs: results_hpc_3dv2/ANALYSIS/{bicrystal_sizeeffect,poly3d_window}.png + json
"""

import sys, os, json, glob
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from pfc3d import PFC3D, find_peaks_3d
from defect_analysis_3d import find_dislocation_lines

ROOT = os.path.join(os.path.dirname(__file__), "..", "results_hpc_3dv2")
OUT = os.path.join(ROOT, "ANALYSIS")
A = 8.886


def gb_distance_profile(npz):
    m = PFC3D.load(npz)
    pts = find_peaks_3d(m.psi, m.dx, m.dy, m.dz)
    box = np.array([m.lx, m.ly, m.lz])
    r = find_dislocation_lines(pts, box)
    zc = pts[:, 2]
    dgb = np.minimum(np.abs(zc - m.lz / 2), np.minimum(zc, m.lz - zc))
    dis = r["disorder"] > 0.3
    shells = [(0, 2), (2, 4), (4, 6), (6, 99)]
    prof = []
    for lo, hi in shells:
        sel = (dgb >= lo * A) & (dgb < hi * A)
        prof.append(dis[sel].mean() if sel.sum() else 0.0)
    return m.exx, prof, r["n_lines"]


def main():
    os.makedirs(OUT, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    metrics = {}

    # (1) GB-distance decomposition at final state, all sizes
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    sizes = sorted(glob.glob(os.path.join(ROOT, "bic3d_*_r-0.3")))
    labels = ["0-2", "2-4", "4-6", ">6"]
    x = np.arange(4)
    w = 0.8 / max(len(sizes), 1)
    bic = {}
    for k, d in enumerate(sizes):
        name = os.path.basename(d)
        exx, prof, nl = gb_distance_profile(os.path.join(d, "final.npz"))
        bic[name] = dict(exx=float(exx), profile=prof, n_lines=nl)
        ax[0].bar(x + k * w, np.array(prof) * 100, w,
                  label=name.replace("bic3d_", "").replace("_r-0.3", ""))
    metrics["bicrystal_gb_distance"] = bic
    ax[0].set_xticks(x + w * len(sizes) / 2)
    ax[0].set_xticklabels([f"{l} a0" for l in labels])
    ax[0].set_xlabel("distance from GB plane")
    ax[0].set_ylabel("% atoms disordered (CSP>0.3)")
    ax[0].set_title("Σ5 bicrystal @9.4%: disorder is GB-localized\n"
                    "(grain interior >4a0 stays clean -> no emission)")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3, axis="y")

    # onset trace for the 160 case
    d = os.path.join(ROOT, "bic3d_n160c10_r-0.3")
    trace = []
    for fp in sorted(glob.glob(os.path.join(d, "snap_*.npz"))
                     + [os.path.join(d, "initial.npz"),
                        os.path.join(d, "final.npz")]):
        exx, prof, nl = gb_distance_profile(fp)
        trace.append((exx * 100, prof[3]))   # interior (>6 a0) disorder
    trace.sort()
    tr = np.array(trace)
    ax[1].plot(tr[:, 0], tr[:, 1] * 100, "ro-")
    ax[1].set_xlabel("strain (%)")
    ax[1].set_ylabel("grain-interior (>6a0) disordered %")
    ax[1].set_title("160^3: interior stays ~0 -> GB thickens, no emission")
    ax[1].grid(alpha=0.3)
    fig.savefig(os.path.join(OUT, "bicrystal_sizeeffect.png"), dpi=140,
                bbox_inches="tight")

    # (2) in-window poly3d maps
    RS = (-0.20, -0.25, -0.30, -0.35)
    RELAXES = (100, 400)
    SEEDS = (7, 11)
    mod = np.full((4, 2), np.nan)
    flow = np.full((4, 2), np.nan)
    for i, r in enumerate(RS):
        for j, rx in enumerate(RELAXES):
            ms, fs = [], []
            for s in SEEDS:
                try:
                    dd = json.load(open(os.path.join(
                        ROOT, f"p3dw_r{r}_x{rx}_s{s}", "summary.json")))
                except FileNotFoundError:
                    continue
                rows = dd["rows"]
                eps = np.array([x["exx"] for x in rows])
                sig = np.array([x["sigma"] for x in rows])
                el = (eps > 0) & (eps <= 0.02)
                if el.sum() >= 2:
                    ms.append(np.polyfit(eps[el], sig[el], 1)[0])
                fs.append(sig[int(np.argmin(np.abs(eps - 0.07)))])
            if ms:
                mod[i, j] = np.mean(ms)
            if fs:
                flow[i, j] = np.mean(fs)
    metrics["poly3d_window"] = dict(r=list(RS), relax=list(RELAXES),
                                    modulus=mod.tolist(), flow7=flow.tolist())
    fig2, ax2 = plt.subplots(1, 2, figsize=(12, 5))
    for j, rx in enumerate(RELAXES):
        ax2[0].plot(RS, mod[:, j], "o-", label=f"RELAX={rx}")
        ax2[1].plot(RS, flow[:, j] * 1e3, "s-", label=f"RELAX={rx}")
    ax2[0].set_xlabel("r"); ax2[0].set_ylabel("modulus")
    ax2[0].set_title("3D poly modulus vs r (<=2%, valid)")
    ax2[0].legend(); ax2[0].grid(alpha=0.3)
    ax2[1].set_xlabel("r"); ax2[1].set_ylabel("flow stress @7% (x1e-3)")
    ax2[1].set_title("3D poly flow stress (within validity window)")
    ax2[1].legend(); ax2[1].grid(alpha=0.3)
    fig2.savefig(os.path.join(OUT, "poly3d_window.png"), dpi=140,
                 bbox_inches="tight")

    with open(os.path.join(OUT, "metrics_3dv2.json"), "w") as f:
        json.dump(metrics, f, indent=1)
    print("interior-disorder by size (final, >6a0 shell):")
    for n, b in metrics["bicrystal_gb_distance"].items():
        print(f"  {n}: interior={b['profile'][3]*100:.1f}% "
              f"(0-2a0={b['profile'][0]*100:.0f}%) exx={b['exx']*100:.1f}%")


if __name__ == "__main__":
    main()

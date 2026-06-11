"""D2 flagship analysis: core trajectories, stress-strain, multiplication
event. Output: results/d2_roi_matched_1536/analysis_d2.png + events JSON.
"""

import os, json
import numpy as np

OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   "d2_roi_matched_1536")


def main():
    with open(os.path.join(OUT, "summary.json")) as f:
        s = json.load(f)
    rows = s["rows"]
    eps = np.array([r["exx"] for r in rows])
    sig = np.array([r["sigma"] for r in rows])
    cores = np.array([r["cores"] for r in rows])
    trajs = [np.array(r["core_xy"]) for r in rows]

    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from pfc2d import PFC2D
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6))

    ax = axes[0]
    ax.plot(eps * 100, sig, "o-", ms=3, color="tab:blue")
    ax2 = ax.twinx()
    ax2.plot(eps * 100, cores, "s--", ms=4, color="tab:red")
    ax2.set_ylabel("# cores", color="tab:red")
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("stress", color="tab:blue")
    ax.set_title(f"D2 ROI-matched: {s['L_nm']:.1f} nm box, "
                 f"rho={s['rho_sim_m2']:.1e} m$^{{-2}}$")
    ax.grid(alpha=0.3)

    ax = axes[1]
    cmap = plt.cm.viridis
    for i, t in enumerate(trajs):
        if len(t):
            ax.scatter(t[:, 0], t[:, 1], s=14, color=cmap(i / len(trajs)))
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(0, eps[-1] * 100))
    plt.colorbar(sm, ax=ax, label="strain (%)")
    ax.set_aspect("equal")
    ax.set_title("core trajectories (glide under tension)")

    # multiplication event field
    ax = axes[2]
    ev_files = sorted(f for f in os.listdir(OUT) if f.startswith("event_"))
    if ev_files:
        m = PFC2D.load(os.path.join(OUT, ev_files[-1]))
        from defect_analysis import find_peaks, find_dislocations
        pts = find_peaks(m.psi, m.dx, m.dy)
        d = find_dislocations(pts, m.lx, m.ly)
        ax.imshow(m.psi, origin="lower", cmap="viridis",
                  extent=[0, m.lx, 0, m.ly])
        if len(d["cores"]):
            ax.plot(d["cores"][:, 0], d["cores"][:, 1], "rx", ms=14, mew=3)
        ax.set_title(f"{ev_files[-1]}: {len(d['cores'])} cores")
        events = dict(file=ev_files[-1], n_cores=len(d["cores"]),
                      cores_xy=d["cores"].tolist())
        with open(os.path.join(OUT, "events.json"), "w") as f:
            json.dump(events, f, indent=1)

    fig.savefig(os.path.join(OUT, "analysis_d2.png"), dpi=140,
                bbox_inches="tight")

    # glide distances
    if len(trajs[0]) and len(trajs[-2]):
        d0 = trajs[0]
        print("initial cores:", np.round(d0, 1).tolist())
    print(f"stress at 10.5%: {sig[-1]:.5f}, cores: {cores[0]}->{cores[-1]}")
    i_mult = int(np.argmax(cores > cores[0]))
    if cores.max() > cores[0]:
        print(f"multiplication at exx={eps[i_mult]*100:.2f}% "
              f"({cores[i_mult-1]}->{cores[i_mult]} cores)")


if __name__ == "__main__":
    main()

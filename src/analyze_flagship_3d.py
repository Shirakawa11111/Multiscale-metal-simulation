"""256^3 flagship analysis: stress-strain + 3D dislocation-line evolution
under tension (the scale where individual lines resolve, unlike 128^3 where
the GB network percolates).

Outputs: results_hpc_3d/ANALYSIS/flagship3d.png, flagship3d.json
"""

import os, sys, json, glob
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc3d import PFC3D, find_peaks_3d
from defect_analysis_3d import find_dislocation_lines

DIR = os.path.join(os.path.dirname(__file__), "..", "results_hpc_3d",
                   "p3d_flagship_256")
OUT = os.path.join(os.path.dirname(__file__), "..", "results_hpc_3d",
                   "ANALYSIS")


def main():
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(DIR, "summary.json")) as f:
        s = json.load(f)
    eps = np.array([r["exx"] for r in s["rows"]])
    sig = np.array([r["sigma"] for r in s["rows"]])

    # line evolution across available snapshots
    snaps = [("initial.npz", 0.0)]
    for fp in sorted(glob.glob(os.path.join(DIR, "snap_*.npz"))):
        pct = float(os.path.basename(fp).split("_")[1].replace("pct.npz", ""))
        snaps.append((os.path.basename(fp), pct / 100))
    snaps.append(("final.npz", eps[-1]))

    line_rows = []
    for fn, e in snaps:
        m = PFC3D.load(os.path.join(DIR, fn))
        pts = find_peaks_3d(m.psi, m.dx, m.dy, m.dz)
        box = np.array([m.lx, m.ly, m.lz])
        r = find_dislocation_lines(pts, box)
        line_rows.append(dict(exx=e, atoms=len(pts), n_lines=r["n_lines"],
                              disordered_frac=r["disordered_frac"],
                              total_line_length=r.get("total_line_length", 0),
                              top_sizes=sorted(r["line_sizes"],
                                               reverse=True)[:5]))
        print(f"{fn}: exx={e*100:.1f}% atoms={len(pts)} lines={r['n_lines']} "
              f"frac={r['disordered_frac']:.3f} "
              f"L_tot={r.get('total_line_length',0):.0f}", flush=True)

    metrics = dict(stress_strain=dict(exx=eps.tolist(), sigma=sig.tolist()),
                   line_evolution=line_rows,
                   box_cells=round(float(np.load(
                       os.path.join(DIR, "initial.npz"))["dx0"]), 4))
    with open(os.path.join(OUT, "flagship3d.json"), "w") as f:
        json.dump(metrics, f, indent=1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    ax[0].plot(eps * 100, sig, "b-o", ms=3)
    ax[0].set_xlabel("strain (%)")
    ax[0].set_ylabel("stress")
    ax[0].set_title("256³ polycrystal stress-strain")
    ax[0].grid(alpha=0.3)
    le = line_rows
    ax[1].plot([r["exx"] * 100 for r in le],
               [r["disordered_frac"] for r in le], "r-s")
    ax[1].set_xlabel("strain (%)")
    ax[1].set_ylabel("disordered fraction")
    ax[1].set_title("defected-atom fraction vs strain")
    ax[1].grid(alpha=0.3)
    ax[2].plot([r["exx"] * 100 for r in le],
               [r["total_line_length"] for r in le], "m-^")
    ax[2].set_xlabel("strain (%)")
    ax[2].set_ylabel("total line length (proxy)")
    ax[2].set_title("dislocation content vs strain")
    ax[2].grid(alpha=0.3)
    fig.savefig(os.path.join(OUT, "flagship3d.png"), dpi=140,
                bbox_inches="tight")
    print("saved figure")


if __name__ == "__main__":
    main()

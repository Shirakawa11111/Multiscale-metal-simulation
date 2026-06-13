"""Kinematic carrier decomposition for the polycrystal rate-sensitivity story.

Splits the defect population into:
  - GB network  : 5|7 atoms in the single largest connected cluster (the
                  percolating grain-boundary network) -> GB-sliding carrier
  - intragranular: 5|7 atoms in small isolated clusters (mobile glide
                  dislocations inside grains) -> dislocation-glide carrier

Tracks the split vs strain for cold (r=-0.35) and near-spinodal (r=-0.25)
to test whether the carrier shifts (glide->GB) with r, or whether plasticity
is GB-mediated at ALL r (which would re-label the m(r) trend as GB-sliding
flow softening toward the spinodal, not a glide->GBS crossover).

Output: results_hpc/ANALYSIS/kinematic_decomp.{png,json}
"""

import sys, os, glob, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from scipy.spatial import cKDTree
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, coordination

ROOT = os.path.join(os.path.dirname(__file__), "..", "results_hpc")
OUT = os.path.join(ROOT, "ANALYSIS")


def carrier_split(m):
    pts = find_peaks(m.psi, m.dx, m.dy)
    coord, _ = coordination(pts, m.lx, m.ly)
    defect = pts[coord != 6]
    n_def = len(defect)
    if n_def < 2:
        return dict(n_defect=n_def, gb_frac=0.0, intra_frac=0.0,
                    n_clusters=0, largest=0)
    # connected components of defect atoms (periodic), link 1.6*a0
    tiled = np.vstack([defect + np.array([sx * m.lx, sy * m.ly])
                       for sx in (-1, 0, 1) for sy in (-1, 0, 1)])
    src = np.concatenate([np.arange(n_def)] * 9)
    tree = cKDTree(tiled)
    parent = list(range(n_def))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n_def):
        for j in tree.query_ball_point(defect[i], r=1.6 * A_LATTICE):
            j0 = int(src[j])
            if j0 != i:
                ri, rj = find(i), find(j0)
                if ri != rj:
                    parent[ri] = rj
    from collections import Counter
    sizes = Counter(find(i) for i in range(n_def))
    sz = sorted(sizes.values(), reverse=True)
    largest = sz[0]
    # GB network = largest percolating cluster; intragranular = the rest
    gb = largest
    intra = n_def - largest
    return dict(n_defect=int(n_def), gb_frac=float(gb / len(pts)),
                intra_frac=float(intra / len(pts)),
                n_clusters=len(sz), largest=int(largest),
                intra_clusters=int(sum(1 for s in sz if s < 6)))


def main():
    os.makedirs(OUT, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    res = {}
    for r, color in [(-0.35, "tab:blue"), (-0.25, "tab:red")]:
        d = os.path.join(ROOT, f"poly_r{r}_x400_s7")
        files = {}
        for fp in (glob.glob(os.path.join(d, "snap_*.npz"))
                   + [os.path.join(d, "initial.npz")]):
            mm = PFC2D.load(fp)
            files[round(mm.exx * 100, 1)] = fp
        traj = []
        for pct in sorted(files):
            m = PFC2D.load(files[pct])
            s = carrier_split(m)
            traj.append((pct, s["gb_frac"], s["intra_frac"],
                         s["intra_clusters"]))
        traj = np.array(traj)
        res[str(r)] = dict(strain_pct=traj[:, 0].tolist(),
                           gb_frac=traj[:, 1].tolist(),
                           intra_frac=traj[:, 2].tolist(),
                           intra_clusters=traj[:, 3].tolist())
        ax[0].plot(traj[:, 0], traj[:, 1], "o-", color=color,
                   label=f"r={r} GB-network")
        ax[0].plot(traj[:, 0], traj[:, 2], "s--", color=color,
                   label=f"r={r} intragranular")
        ax[1].plot(traj[:, 0], traj[:, 3], "^-", color=color, label=f"r={r}")
    ax[0].set_xlabel("strain (%)")
    ax[0].set_ylabel("defect fraction")
    ax[0].set_title("carrier split: GB-network vs intragranular")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)
    ax[1].set_xlabel("strain (%)")
    ax[1].set_ylabel("# isolated intragranular clusters")
    ax[1].set_title("mobile glide-dislocation count (small clusters)")
    ax[1].legend()
    ax[1].grid(alpha=0.3)
    fig.savefig(os.path.join(OUT, "kinematic_decomp.png"), dpi=140,
                bbox_inches="tight")
    with open(os.path.join(OUT, "kinematic_decomp.json"), "w") as f:
        json.dump(res, f, indent=1)
    for r, v in res.items():
        print(f"r={r}: GB-frac {v['gb_frac'][0]:.3f}->{v['gb_frac'][-1]:.3f}, "
              f"intra-frac {v['intra_frac'][0]:.3f}->{v['intra_frac'][-1]:.3f}, "
              f"isolated clusters {v['intra_clusters'][0]}->{v['intra_clusters'][-1]}")


if __name__ == "__main__":
    main()

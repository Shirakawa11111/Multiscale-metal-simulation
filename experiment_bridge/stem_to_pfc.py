"""Experiment->simulation bridge (capstone): seed a PFC simulation from REAL
STEM-reconstructed 3D dislocation lines, then evolve under load.

This is the end-to-end "image -> simulation" data-assimilation demo that closes
the gap the pure-simulation work lacked (no experimental anchor): the
phase-winding defect-introduction algorithm (which migrated cleanly 2D->3D)
ingests the experimentally reconstructed dislocation geometry directly.

Input: 3D重建算法论文/3d_scatter/points_3d*.txt — each file is one
reconstructed dislocation line (polyline of normalized 3D points; z is the thin
foil-normal direction, ~planar). We map (x,y) into a PFC box and seed an edge
dislocation per line.

Honest simplifications (flagged for the real pipeline):
  - Burgers vectors are assigned uniform (along x) with alternating sign to net
    zero (PBC requirement); REAL signs/Burgers need experimental Burgers
    analysis (g.b) per line.
  - 2D PFC (the reconstruction is ~planar); a full 3D line seeding uses
    pfc3d.init_dislocation_lines.

Output: experiment_bridge/results/
"""

import os, sys, glob, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

RECON = os.path.join(os.path.dirname(__file__), "..", "..", "3D重建算法论文",
                     "3d_scatter")
OUT = os.path.join(os.path.dirname(__file__), "results")
N = 512
DT = 0.5


def load_recon_lines():
    lines = []
    for fp in sorted(glob.glob(os.path.join(RECON, "points_3d*.txt"))):
        pts = []
        for ln in open(fp):
            v = ln.split()
            if len(v) == 3:
                pts.append([float(v[0]), float(v[1]), float(v[2])])
        if len(pts) >= 2:
            lines.append(np.array(pts))
    return lines


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    lines = load_recon_lines()
    print(f"loaded {len(lines)} reconstructed dislocation lines", flush=True)
    # global (x,y) bbox -> map into the box interior [0.1,0.9]
    allpts = np.vstack([ln[:, :2] for ln in lines])
    lo, hi = allpts.min(0), allpts.max(0)
    span = (hi - lo).max()

    def to_frac(xy):
        return 0.1 + 0.8 * (xy - lo) / span

    # one edge dislocation per line at its centroid; uniform Burgers (x),
    # alternating sign to net-zero (use an even count).
    n_use = len(lines) - (len(lines) % 2)
    cores = []
    for i in range(n_use):
        c = to_frac(lines[i][:, :2].mean(0))
        cores.append((float(c[0]), float(c[1]), +1 if i % 2 == 0 else -1, 0.0))
    print(f"seeding {len(cores)} edge dislocations from reconstruction "
          f"(net-zero, uniform Burgers x)", flush=True)

    m = PFC2D(N, N, r=-0.25, psi_bar=-0.25)
    m.init_dislocations(cores)
    m.save(os.path.join(OUT, "seeded.npz"))
    d0 = find_dislocations(find_peaks(m.psi, m.dx, m.dy), m.lx, m.ly)
    print(f"  t=0 seeded: detected cores={len(d0['cores'])} (target {len(cores)})",
          flush=True)
    # relax (assimilation: let the lattice accommodate the real configuration)
    m.step(DT, n=400)
    d1 = find_dislocations(find_peaks(m.psi, m.dx, m.dy), m.lx, m.ly)
    print(f"  after relax: cores={len(d1['cores'])} (config stable if ~retained)",
          flush=True)
    m.save(os.path.join(OUT, "relaxed.npz"))

    # evolve under area-conserving tension (assimilation -> prediction)
    rows = []
    for i in range(16):
        m.apply_strain(0.0025, area_conserving=True)
        m.step(DT, n=200)
        d = find_dislocations(find_peaks(m.psi, m.dx, m.dy), m.lx, m.ly)
        rows.append(dict(exx=m.exx, sigma=m.stress(), cores=len(d["cores"])))
        if i % 4 == 3:
            print(f"  exx={m.exx*100:.1f}% sigma={rows[-1]['sigma']:+.5f} "
                  f"cores={len(d['cores'])}", flush=True)
    m.save(os.path.join(OUT, "evolved.npz"))

    result = dict(n_recon_lines=len(lines), n_seeded=len(cores),
                  cores_t0=len(d0["cores"]), cores_relaxed=len(d1["cores"]),
                  tension=rows, wall_s=time.time() - t0,
                  note="REAL STEM-reconstructed dislocations seeded into PFC and "
                       "evolved under tension — end-to-end image->simulation "
                       "data-assimilation pipeline.")
    with open(os.path.join(OUT, "bridge.json"), "w") as f:
        json.dump(result, f, indent=1)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        rel = PFC2D.load(os.path.join(OUT, "relaxed.npz"))
        pts = find_peaks(rel.psi, rel.dx, rel.dy)
        d = find_dislocations(pts, rel.lx, rel.ly)
        fig, ax = plt.subplots(1, 2, figsize=(14, 6))
        # reconstructed lines (mapped) overlaid on PFC core positions
        for ln in lines:
            fr = to_frac(ln[:, :2])
            ax[0].plot(fr[:, 0] * rel.lx, fr[:, 1] * rel.ly, "-", lw=1.5,
                       color="tab:blue", alpha=0.7)
        if len(d["cores"]):
            ax[0].plot(d["cores"][:, 0], d["cores"][:, 1], "rx", ms=8, mew=2,
                       label="PFC cores (relaxed)")
        ax[0].set_title("STEM-reconstructed lines (blue) -> seeded PFC cores (red)")
        ax[0].legend()
        ax[0].set_aspect("equal")
        ax[1].imshow(rel.psi, origin="lower", cmap="viridis",
                     extent=[0, rel.lx, 0, rel.ly])
        ax[1].set_title("PFC density (real-config seeded, relaxed)")
        fig.savefig(os.path.join(OUT, "bridge.png"), dpi=130,
                    bbox_inches="tight")
        print("figure saved", flush=True)
    except Exception as ex:
        print("plot skipped:", ex)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

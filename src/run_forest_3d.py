"""3D forest-hardening test — the redirect after 2D decisively could not
harden (junctions are 3D). Mobile dislocation lines glide under shear through
a forest of lines on intersecting axes; if 3D PFC forms stabilizing junctions,
flow stress should RISE with forest density (alpha>0) — unlike 2D.

Loading: simple xz shear (3D BCC is shear-stable to >=15%, M18, dodging the
tension amorphization ceiling). Mobile lines: along z, Burgers x (glide in x
under eps_xz). Forest: lines along x and y (thread the x-z glide plane) ->
the mobile line must cut/junction with them.

Env: FOREST3D_N (grid), FOREST3D_CELLS, FOREST3D_NF (forest line pairs).
Output: results/forest_3d/<tag>/
"""

import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc3d import PFC3D, A_BCC, find_peaks_3d
from defect_analysis_3d import find_dislocation_lines

N = int(os.environ.get("FOREST3D_N", "96"))
CELLS = int(os.environ.get("FOREST3D_CELLS", "9"))
NF = int(os.environ.get("FOREST3D_NF", "0"))
SEED = int(os.environ.get("FOREST3D_SEED", "7"))
TAG = os.environ.get("FOREST3D_TAG", f"n{N}c{CELLS}_nf{NF}_s{SEED}")
OUT = os.path.join(os.path.dirname(__file__), "..", "results", "forest_3d", TAG)
DX = CELLS * A_BCC / N
DT = 0.5
DG = 0.0025
N_SHEAR = 32          # -> 8% shear
RELAX = 150
BETA = 10.0


def build():
    rng = np.random.default_rng(SEED)
    lines = [
        # mobile dipole: z-lines, Burgers +-x, on two glide planes (y offset)
        dict(axis="z", pos=(0.5, 0.35), burgers="x", sign=+1),
        dict(axis="z", pos=(0.5, 0.65), burgers="x", sign=-1),
    ]
    # forest: pairs of lines along x and y, Burgers along z (thread x-z plane),
    # net Burgers cancels within each pair
    for i in range(NF):
        ax = "x" if i % 2 == 0 else "y"
        a = rng.uniform(0.15, 0.85)
        b = rng.uniform(0.15, 0.45)
        lines.append(dict(axis=ax, pos=(a, b), burgers="z", sign=+1))
        lines.append(dict(axis=ax, pos=(a, min(0.85, b + 0.4)),
                          burgers="z", sign=-1))
    m = PFC3D(N, N, N, dx=DX, r=-0.25, psi_bar=-0.25)
    m.init_dislocation_lines(lines)
    m.step_mpfc(DT, n=500, beta=BETA)
    return m


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    m = build()
    box = np.array([m.lx, m.ly, m.lz])
    pts = find_peaks_3d(m.psi, m.dx, m.dy, m.dz)
    r0 = find_dislocation_lines(pts, box)
    rows = [dict(gxz=0.0, tau=m.shear_stress(), frac=r0["disordered_frac"],
                 n_lines=r0["n_lines"])]
    print(f"[{TAG}] seeded: frac={r0['disordered_frac']:.3f} "
          f"lines={r0['n_lines']} atoms={len(pts)}", flush=True)
    for i in range(N_SHEAR):
        m.apply_shear(DG)
        m.step_mpfc(DT, n=RELAX, beta=BETA)
        if i % 4 == 3 or i == N_SHEAR - 1:
            pts = find_peaks_3d(m.psi, m.dx, m.dy, m.dz)
            rr = find_dislocation_lines(pts, box)
            rows.append(dict(gxz=m.gxz, tau=m.shear_stress(),
                             frac=rr["disordered_frac"], n_lines=rr["n_lines"]))
            print(f"[{TAG}] gxz={m.gxz*100:.1f}% tau={rows[-1]['tau']:+.6f} "
                  f"frac={rr['disordered_frac']:.3f} lines={rr['n_lines']}",
                  flush=True)
    tau_flow = float(np.mean([r["tau"] for r in rows[-4:]]))
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, NF=NF, N=N, CELLS=CELLS, seed=SEED,
                       tau_flow=tau_flow, wall_s=time.time() - t0), f, indent=1)
    m.save(os.path.join(OUT, "final.npz"))
    print(f"[{TAG}] tau_flow={tau_flow:.5f} done {time.time()-t0:.0f}s",
          flush=True)


if __name__ == "__main__":
    main()

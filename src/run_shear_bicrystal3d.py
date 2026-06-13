"""3D Σ5 bicrystal under SIMPLE SHEAR (xz) — the protocol that avoids the
~10-12% uniaxial amorphization ceiling (3D BCC shear is stable to >=15%,
verified). Tests whether a clean twist GB glides / emits dislocations under
resolved shear, which uniaxial tension could not resolve (M16).

Output: results/shear_bic3d_<tag>/
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc3d import PFC3D, A_BCC, find_peaks_3d
from defect_analysis_3d import find_dislocation_lines

N = int(os.environ.get("SBIC_N", "128"))
CELLS = int(os.environ.get("SBIC_CELLS", "10"))
R = float(os.environ.get("SBIC_R", "-0.25"))
TAG = os.environ.get("SBIC_TAG", f"n{N}c{CELLS}_r{R}")
OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   f"shear_bic3d_{TAG}")
NZ = N * 2
DX = CELLS * A_BCC / N
DGAMMA = 0.0025
N_STEPS = 60          # -> gamma_xz = 15%
RELAX = 200
DT = 0.5


def analyze(m, poly=False):
    pts = find_peaks_3d(m.psi, m.dx, m.dy, m.dz)
    box = np.array([m.lx, m.ly, m.lz])
    r = find_dislocation_lines(pts, box)
    if poly:
        # no GB plane in a polycrystal: report total CSP-disordered fraction
        return r, len(pts), float(r["disordered_frac"])
    zc = pts[:, 2]
    dgb = np.minimum(np.abs(zc - m.lz / 2), np.minimum(zc, m.lz - zc))
    interior = dgb > 6 * A_BCC
    bulk_dis = ((r["disorder"][interior] > 0.3).mean()
                if interior.sum() else 0.0)
    return r, len(pts), float(bulk_dis)


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    if os.environ.get("SBIC_KIND") == "poly":
        m = PFC3D(N, N, N, dx=DX, r=R, psi_bar=-0.25)
        m.init_random(noise=0.05, seed=int(os.environ.get("SBIC_SEED", "7")))
        m.step(DT, n=2000)
    else:
        m = PFC3D(N, N, NZ, dx=DX, r=R, psi_bar=-0.25)
        m.init_bicrystal_csl()
        m.step(DT, n=800)
    m.save(os.path.join(OUT, "initial.npz"))
    r, npts, bulk = analyze(m, poly=(os.environ.get("SBIC_KIND")=="poly"))
    rows = [dict(gxz=0.0, tau=m.shear_stress(), F=m.free_energy(),
                 atoms=npts, n_lines=r["n_lines"], interior_dis=bulk)]
    print(f"[{TAG}] g=0 atoms={npts} lines={r['n_lines']} interior_dis={bulk:.3f}",
          flush=True)
    for i in range(N_STEPS):
        m.apply_shear(DGAMMA)
        m.step(DT, n=RELAX)
        r, npts, bulk = analyze(m, poly=(os.environ.get("SBIC_KIND")=="poly"))
        rows.append(dict(gxz=m.gxz, tau=m.shear_stress(), F=m.free_energy(),
                         atoms=npts, n_lines=r["n_lines"], interior_dis=bulk))
        if i % 4 == 3:
            print(f"[{TAG}] gxz={m.gxz*100:5.2f}% tau={rows[-1]['tau']:+.6f} "
                  f"lines={r['n_lines']} interior_dis={bulk:.3f} atoms={npts}",
                  flush=True)
        if i % 12 == 11:
            m.save(os.path.join(OUT, f"snap_{m.gxz*100:.1f}pct.npz"))
    m.save(os.path.join(OUT, "final.npz"))
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, N=N, CELLS=CELLS, R=R, relax=RELAX,
                       wall_s=time.time() - t0), f, indent=1)
    print(f"[{TAG}] done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

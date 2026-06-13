"""M16 done properly: tension of a clean Σ5 twist bicrystal (3D).

The Σ5 CSL twist GB has defects ONLY at the boundary (zero bulk mismatch,
verified). Under volume-conserving tension the GB acts as a dislocation
source: emitted lines appear as NEW line clusters in the previously clean
bulk, and disordered_frac rises above its GB baseline. This isolates
GB-mediated dislocation nucleation — the 3D analogue of the 2D pore
"dislocation factory" (M6) but at a crystallographically clean interface.

Output: results/bicrystal_csl/  (per-step metrics + snapshots + figure)
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc3d import PFC3D, A_BCC, find_peaks_3d
from defect_analysis_3d import find_dislocation_lines

OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   "bicrystal_csl")
NX = 80          # 5 in-plane cells (Σ5 commensurate)
NZ = 160         # 10 cells along z (two 5-cell grains)
DX = 5 * A_BCC / NX
DEPS = 0.0025
N_STRAIN = 32    # -> 8%
RELAX = 200
DT = 0.5


def analyze(m):
    pts = find_peaks_3d(m.psi, m.dx, m.dy, m.dz)
    box = np.array([m.lx, m.ly, m.lz])
    r = find_dislocation_lines(pts, box)
    zc = pts[:, 2]
    near_gb = ((np.abs(zc - m.lz / 2) < 3 * A_BCC)
               | (np.minimum(zc, m.lz - zc) < 3 * A_BCC))
    bulk_dis = int((r["disorder"][~near_gb] > 0.06).sum())
    return r, len(pts), bulk_dis


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    m = PFC3D(NX, NX, NZ, dx=DX, r=-0.25, psi_bar=-0.25)
    m.init_bicrystal_csl()
    m.step(DT, n=800)            # relax the GB
    m.save(os.path.join(OUT, "initial.npz"))

    rows = []
    r, npts, bulk = analyze(m)
    rows.append(dict(exx=0.0, sigma=m.stress(), F=m.free_energy(),
                     atoms=npts, n_lines=r["n_lines"],
                     disordered_frac=r["disordered_frac"],
                     bulk_disordered=bulk))
    print(f"eps=0: atoms={npts} lines={r['n_lines']} "
          f"frac={r['disordered_frac']:.3f} bulk_dis={bulk}", flush=True)

    base_bulk = bulk
    for i in range(N_STRAIN):
        m.apply_strain(DEPS, volume_conserving=True)
        m.step(DT, n=RELAX)
        r, npts, bulk = analyze(m)
        rows.append(dict(exx=m.exx, sigma=m.stress(), F=m.free_energy(),
                         atoms=npts, n_lines=r["n_lines"],
                         disordered_frac=r["disordered_frac"],
                         bulk_disordered=bulk))
        print(f"exx={m.exx*100:5.2f}% sigma={rows[-1]['sigma']:+.6f} "
              f"lines={r['n_lines']} bulk_dis={bulk}", flush=True)
        if bulk > base_bulk + 5 and not os.path.exists(
                os.path.join(OUT, "emission.npz")):
            m.save(os.path.join(OUT, "emission.npz"))
            print(f"  >>> GB emission at exx={m.exx*100:.2f}% "
                  f"(bulk {base_bulk}->{bulk})", flush=True)
        if i % 8 == 7:
            m.save(os.path.join(OUT, f"snap_{m.exx*100:.2f}pct.npz"))

    m.save(os.path.join(OUT, "final.npz"))
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, NX=NX, NZ=NZ, relax=RELAX,
                       wall_s=time.time() - t0), f, indent=1)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        eps = [r["exx"] * 100 for r in rows]
        fig, ax = plt.subplots(1, 3, figsize=(17, 5))
        ax[0].plot(eps, [r["sigma"] for r in rows], "b-o", ms=3)
        ax[0].set_xlabel("strain (%)"); ax[0].set_ylabel("stress")
        ax[0].set_title("Σ5 bicrystal stress-strain"); ax[0].grid(alpha=0.3)
        ax[1].plot(eps, [r["bulk_disordered"] for r in rows], "r-s", ms=3)
        ax[1].set_xlabel("strain (%)")
        ax[1].set_ylabel("bulk disordered atoms")
        ax[1].set_title("GB dislocation emission into bulk")
        ax[1].grid(alpha=0.3)
        ax[2].plot(eps, [r["n_lines"] for r in rows], "m-^", ms=3)
        ax[2].set_xlabel("strain (%)"); ax[2].set_ylabel("# line clusters")
        ax[2].set_title("dislocation line count"); ax[2].grid(alpha=0.3)
        fig.savefig(os.path.join(OUT, "bicrystal_csl.png"), dpi=140,
                    bbox_inches="tight")
    except Exception as ex:
        print("plot skipped:", ex)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

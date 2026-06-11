"""B4/B5 mechanism runs (512^2 tension, same protocol as run_b_series).

B4 heterogeneous/homogeneous nucleation:
  b4_pore    perfect crystal + mass-depleted pore (persistent stress
             concentrator) -> dislocation emission from the pore rim
  b4_noise   perfect crystal + thermal noise during straining
             (noise_amp=0.02) -> noise-assisted nucleation
B5 strain-rate effect on the polycrystal (b3 @ RELAX=400 is the midpoint):
  b5_rate_fast   RELAX=100  (4x faster effective strain rate)
  b5_rate_slow   RELAX=1600 (4x slower)

Output: results/b45_series_512/<case>/
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D
from defect_analysis import find_peaks, find_dislocations

BASE = os.path.join(os.path.dirname(__file__), "..", "results",
                    "b45_series_512")

DEPS = 0.0025
N_STRAIN = 40
DT = 0.5
N = 512

CASES = {
    "b4_pore":      dict(kind="pore", relax=400, noise=0.0),
    "b4_noise":     dict(kind="perfect", relax=400, noise=0.02),
    "b5_rate_fast": dict(kind="poly", relax=100, noise=0.0),
    "b5_rate_slow": dict(kind="poly", relax=1600, noise=0.0),
    # v1 pore (depth 0.2, 4a0) dissolved via Gibbs-Thomson by ~6% strain;
    # depth 0.5 / 8a0 verified persistent to t=3000
    "b4_pore_v2":   dict(kind="pore2", relax=400, noise=0.0),
    "b4_noise_v2":  dict(kind="perfect", relax=400, noise=0.06),
}


def make_case(kind):
    m = PFC2D(N, N, r=-0.25, psi_bar=-0.25)
    if kind == "pore":
        m.init_crystal()
        m.step(DT, n=200)
        m.add_void(radius=4 * 7.255, depth=0.2)
        m.step(DT, n=600)
    elif kind == "pore2":
        m.init_crystal()
        m.step(DT, n=200)
        m.add_void(radius=8 * 7.255, depth=0.5)
        m.step(DT, n=600)
    elif kind == "perfect":
        m.init_crystal()
        m.step(DT, n=300)
    elif kind == "poly":
        m.init_random(noise=0.05, seed=7)
        m.step(DT, n=3000)
    return m


def run_case(name, cfg):
    out = os.path.join(BASE, name)
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    m = make_case(cfg["kind"])
    pts = find_peaks(m.psi, m.dx, m.dy)
    d = find_dislocations(pts, m.lx, m.ly)
    base_cores = len(d["cores"])
    m.save(os.path.join(out, "initial.npz"))
    rows = [dict(exx=m.exx, F=m.free_energy(), sigma=m.stress(),
                 cores=base_cores, rho=d["rho"])]
    print(f"[{name}] eps=0 cores={base_cores}", flush=True)

    n_events = 0
    for i in range(N_STRAIN):
        m.apply_strain(DEPS, area_conserving=True)
        m.step(DT, n=cfg["relax"], noise_amp=cfg["noise"], seed=1000 + i)
        pts = find_peaks(m.psi, m.dx, m.dy)
        d = find_dislocations(pts, m.lx, m.ly)
        nc = len(d["cores"])
        rows.append(dict(exx=m.exx, F=m.free_energy(), sigma=m.stress(),
                         cores=nc, rho=d["rho"]))
        print(f"[{name}] exx={m.exx*100:5.2f}% sigma={rows[-1]['sigma']:+.5f} "
              f"cores={nc}", flush=True)
        if nc != rows[-2]["cores"] and n_events < 6:
            m.save(os.path.join(out, f"event_{m.exx*100:.2f}pct.npz"))
            n_events += 1
        if i % 8 == 7:
            m.save(os.path.join(out, f"snap_{m.exx*100:.2f}pct.npz"))

    m.save(os.path.join(out, "final.npz"))
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, wall_s=time.time() - t0, n=N, deps=DEPS,
                       cfg=cfg, base_cores=base_cores), f, indent=1)
    print(f"[{name}] done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    names = sys.argv[1:] or list(CASES)
    for nm in names:
        run_case(nm, CASES[nm])

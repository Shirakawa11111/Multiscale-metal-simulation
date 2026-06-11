"""B-series mechanism runs (512^2, area-conserving tension, sigma(eps) curves).

B1b void:        perfect crystal + melted disc notch -> heterogeneous nucleation
B2  dipole:      pre-seeded +b/-b dipole strained to 10% -> multiplication?
B3  polycrystal: melt-quenched polycrystal -> GB-mediated plasticity

Each case logs strain, free energy, stress, core count, dislocation density,
and saves snapshots around nucleation/first-multiplication events.
Output: results/b_series_512/<case>/
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D
from defect_analysis import find_peaks, find_dislocations

BASE = os.path.join(os.path.dirname(__file__), "..", "results", "b_series_512")

DEPS = 0.0025
N_STRAIN = 40            # -> 10%
RELAX = 400
DT = 0.5
N = 512


def make_case(name):
    m = PFC2D(N, N, r=-0.25, psi_bar=-0.25)
    if name == "b1b_void":
        m.init_crystal()
        m.step(DT, n=200)
        m.add_void(radius=4 * 7.255)
        m.step(DT, n=400)
    elif name == "b2_dipole":
        m.init_dislocation_dipole()
        m.step(DT, n=400)
    elif name == "b3_polycrystal":
        m.init_random(noise=0.05, seed=7)
        m.step(DT, n=3000)  # grow + anneal grains
    return m


def run_case(name):
    out = os.path.join(BASE, name)
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    m = make_case(name)
    pts = find_peaks(m.psi, m.dx, m.dy)
    d = find_dislocations(pts, m.lx, m.ly)
    base_cores = len(d["cores"])
    m.save(os.path.join(out, "initial.npz"))
    rows = [dict(exx=m.exx, F=m.free_energy(), sigma=m.stress(),
                 cores=base_cores, rho=d["rho"])]
    print(f"[{name}] eps=0 cores={base_cores}", flush=True)

    event_logged = False
    for i in range(N_STRAIN):
        m.apply_strain(DEPS, area_conserving=True)
        m.step(DT, n=RELAX)
        pts = find_peaks(m.psi, m.dx, m.dy)
        d = find_dislocations(pts, m.lx, m.ly)
        nc = len(d["cores"])
        rows.append(dict(exx=m.exx, F=m.free_energy(), sigma=m.stress(),
                         cores=nc, rho=d["rho"]))
        print(f"[{name}] exx={m.exx*100:5.2f}% sigma={rows[-1]['sigma']:+.5f} "
              f"cores={nc}", flush=True)
        if not event_logged and nc > base_cores:
            m.save(os.path.join(out, f"event_{m.exx*100:.2f}pct.npz"))
            event_logged = True
        if i % 8 == 7:
            m.save(os.path.join(out, f"snap_{m.exx*100:.2f}pct.npz"))

    m.save(os.path.join(out, "final.npz"))
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, wall_s=time.time() - t0, n=N,
                       deps=DEPS, relax=RELAX, dt=DT,
                       base_cores=base_cores), f, indent=1)
    print(f"[{name}] done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    cases = sys.argv[1:] or ["b1b_void", "b2_dipole", "b3_polycrystal"]
    for c in cases:
        run_case(c)

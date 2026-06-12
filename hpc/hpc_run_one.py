"""Unified single-run driver for the HPC parameter matrix.

One process = one matrix cell, single-threaded FFT (PFC_FFT_THREADS=1 set by
the queue runner); 64 concurrent processes fill the 64-CPU allowance.

Kinds:
  poly   polycrystal tension      (quench seed, RELAX = inverse strain rate)
  quad   quadrupole tension       (multiplication threshold)
  cyc    polycrystal cyclic       (amp, 8 cycles)

Usage:
  python3 hpc_run_one.py --kind poly --r -0.25 --relax 400 --seed 7 \
      --n 512 --strain-to 0.16 --out results_hpc/poly_r-0.25_x400_s7
"""

import argparse, json, os, sys, time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "src"))

import numpy as np
from pfc2d import PFC2D, FFT_BACKEND
from defect_analysis import find_peaks, find_dislocations

DT = 0.5
DEPS = 0.0025


def detect(m):
    pts = find_peaks(m.psi, m.dx, m.dy)
    return find_dislocations(pts, m.lx, m.ly)


def make_model(a):
    m = PFC2D(a.n, a.n, r=a.r, psi_bar=-0.25)
    if a.kind in ("poly", "cyc"):
        m.init_random(noise=0.05, seed=a.seed)
        m.step(DT, n=3000)
    else:  # quad
        m.init_dislocations([(0.33, 0.25, +1), (0.33, 0.75, -1),
                             (0.67, 0.25, -1), (0.67, 0.75, +1)])
        m.step(DT, n=600)
    return m


def run_tension(m, a, rows):
    n_steps = int(round(a.strain_to / DEPS))
    for i in range(n_steps):
        m.apply_strain(DEPS, area_conserving=True)
        m.step(DT, n=a.relax)
        d = detect(m)
        rows.append(dict(exx=m.exx, sigma=m.stress(), F=m.free_energy(),
                         cores=len(d["cores"]), rho=d["rho"]))
        if i % 8 == 7:
            print(f"exx={m.exx*100:.2f}% cores={rows[-1]['cores']}",
                  flush=True)
        if i % 16 == 15:
            m.save(os.path.join(a.out, f"snap_{m.exx*100:.2f}pct.npz"))


def run_cyclic(m, a, rows):
    nq = int(round(a.amp / DEPS))
    pattern = [+1] * nq + [-1] * 2 * nq + [+1] * nq
    for cyc in range(a.cycles):
        for sgn in pattern:
            m.apply_strain(sgn * DEPS, area_conserving=True)
            m.step(DT, n=a.relax)
            d = detect(m)
            rows.append(dict(cycle=cyc, exx=m.exx, sigma=m.stress(),
                             cores=len(d["cores"]), rho=d["rho"]))
        print(f"cycle {cyc}: cores={rows[-1]['cores']}", flush=True)
        m.save(os.path.join(a.out, f"cycle_{cyc}.npz"))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--kind", choices=["poly", "quad", "cyc"], required=True)
    p.add_argument("--r", type=float, default=-0.25)
    p.add_argument("--relax", type=int, default=400)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--n", type=int, default=512)
    p.add_argument("--strain-to", type=float, default=0.16)
    p.add_argument("--amp", type=float, default=0.01)
    p.add_argument("--cycles", type=int, default=8)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    print(f"backend={FFT_BACKEND} cfg={vars(a)}", flush=True)
    m = make_model(a)
    d0 = detect(m)
    m.save(os.path.join(a.out, "initial.npz"))
    rows = [dict(exx=0.0, sigma=m.stress(), F=m.free_energy(),
                 cores=len(d0["cores"]), rho=d0["rho"])]

    if a.kind == "cyc":
        run_cyclic(m, a, rows)
    else:
        run_tension(m, a, rows)

    m.save(os.path.join(a.out, "final.npz"))
    # atomic write: summary.json existence is the queue's resume-skip marker,
    # so a partial file from a killed worker must never be left behind
    tmp = os.path.join(a.out, "summary.json.tmp")
    with open(tmp, "w") as f:
        json.dump(dict(rows=rows, cfg=vars(a), backend=FFT_BACKEND,
                       wall_s=time.time() - t0), f, indent=1)
    os.replace(tmp, os.path.join(a.out, "summary.json"))
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

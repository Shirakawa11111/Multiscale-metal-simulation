"""Decisive diagnostic: does a PERFECTLY SESSILE (pinned) forest harden?

This isolates the dynamics question from everything else. We seed a forest +
mobile dislocations, then PIN the forest core regions (blend psi back to its
captured reference value each sub-step), making the forest truly immobile
(no climb, no glide, no annihilation) — the ideal limit a true glide-only
dynamics would approximate. The mobile dislocations glide freely under shear
and must interact with the held forest. Because PFC applies shear via the
k-metric (psi stays on the reference grid), a reference-frame pin co-shears
with the lattice.

DECISION:
  - pinned forest HARDENS (tau_flow rises with pinned density) -> keeping the
    forest sessile is the lever -> a true glide-only dynamics is worth building
  - pinned forest still SOFTENS / mobile bypasses (Orowan) -> hardening is
    blocked by crystallography/fundamentals, not just forest mobility (matches
    Berry 2014 pre-built LC lock -> Orowan bypass) -> dynamics fix won't suffice

Compares PINNED vs UNPINNED at matched seeded density.
Output: hardening_program/results_pinned/
"""

import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "results_pinned")
N = 384
DT = 0.5
DG = 0.0025
N_SHEAR = 28
RELAX = 150
MU = 0.0545


def build(n_forest, seed):
    rng = np.random.default_rng(seed)
    cores = [(0.3, 0.30, +1, 0.0), (0.3, 0.70, -1, 0.0),
             (0.7, 0.30, +1, 0.0), (0.7, 0.70, -1, 0.0)]
    forest_xy = []
    for i in range(n_forest):
        ang = 60.0 if i % 2 == 0 else 120.0
        x = rng.uniform(0.12, 0.88)
        y = rng.uniform(0.12, 0.42)
        cores.append((x, y, +1, ang))
        cores.append((x, min(0.88, y + 0.45), -1, ang))
        forest_xy += [(x, y), (x, min(0.88, y + 0.45))]
    m = PFC2D(N, N, r=-0.25, psi_bar=-0.25)
    m.init_dislocations(cores)
    m.step(DT, n=600)
    return m, forest_xy


def forest_mask(m, forest_xy, radius=2.0 * A_LATTICE):
    x = (np.arange(m.nx) * m.dx)[None, :]
    y = (np.arange(m.ny) * m.dy)[:, None]
    mask = np.zeros((m.ny, m.nx), bool)
    for fx, fy in forest_xy:
        cx, cy = fx * m.lx, fy * m.ly
        rr = np.sqrt(((x - cx + m.lx / 2) % m.lx - m.lx / 2) ** 2
                     + ((y - cy + m.ly / 2) % m.ly - m.ly / 2) ** 2)
        mask |= rr < radius
    return mask


def flow_stress(m, pin_mask=None, pin_target=None):
    taus = []
    for i in range(N_SHEAR):
        m.apply_shear(DG)
        for _ in range(5):                      # relax in sub-chunks
            m.step(DT, n=RELAX // 5)
            if pin_mask is not None:            # re-impose the pin
                m.psi = np.where(pin_mask, pin_target, m.psi)
        taus.append(m.shear_stress())
    return float(np.mean(taus[-N_SHEAR // 3:]))


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    res = {}
    for pinned in (False, True):
        rows = []
        for nf in (4, 12, 20):
            tfs = []
            for seed in (7, 11):
                m, fxy = build(nf, seed)
                if pinned:
                    mask = forest_mask(m, fxy)
                    target = m.psi.copy()
                    tf = flow_stress(m, mask, target)
                else:
                    tf = flow_stress(m)
                tfs.append(tf)
            rows.append(dict(n_forest=nf, seeded=4 + 2 * nf,
                             flow_stress=float(np.mean(tfs))))
            print(f"  pinned={pinned} nf={nf}: tau_flow={np.mean(tfs):.5f}",
                  flush=True)
        sc = np.array([r["seeded"] for r in rows], float)
        tf = np.array([r["flow_stress"] for r in rows])
        slope = float(np.polyfit(np.sqrt(sc), tf, 1)[0])
        res["pinned" if pinned else "free"] = dict(
            rows=rows, taylor_slope=slope, hardens=bool(slope > 0))
        print(f"  pinned={pinned}: d(tau)/d(sqrt(N))={slope:+.6f} "
              f"=> {'HARDENS' if slope > 0 else 'softens'}", flush=True)
    with open(os.path.join(OUT, "pinned_forest.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(f"\nVERDICT: pinned forest "
          f"{'HARDENS -> dynamics is the lever, build true glide-only' if res['pinned']['hardens'] else 'still SOFTENS/bypassed -> crystallography/fundamental block, dynamics fix insufficient'}")
    print(f"  free slope={res['free']['taylor_slope']:+.6f}, "
          f"pinned slope={res['pinned']['taylor_slope']:+.6f}")
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

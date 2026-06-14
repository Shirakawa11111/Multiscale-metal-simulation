"""DIAGNOSTIC v3: cleanest annihilation/mobility signal, tracking-free.

The count-based signal (number of 5|7 dislocation cores) is immune to the
greedy-matching ambiguity that inflates the v1/v2 'path length'. Here we:
  - seed the same 2 mobile (0deg) + 2 forest (60/120deg) dipoles (8 cores)
  - relax, then report n5, n7 and core count
  - apply 8% x-shear, reporting at each step the TOTAL 5-coord and 7-coord
    atom counts and the paired-core count. A monotone/strong drop = the
    seeded dislocations are annihilating, i.e. mobile soft modes, NOT sessile
    obstacles.
  - CONTROL: a no-drive baseline of equal duration to show whether the count
    is stable absent shear (sessile) or already decays (mobile/diffusive).
"""

import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results",
                   "diag_forest_mobility")
N = 256
DT = 0.5
RELAX0 = 700
DGAMMA = 0.0025
N_SHEAR = 32
RELAX = 150


def diag(m):
    pts = find_peaks(m.psi, m.dx, m.dy)
    d = find_dislocations(pts, m.lx, m.ly)
    return d["n5"], d["n7"], len(d["cores"])


def build():
    cores = [(0.30, 0.35, +1, 0.0), (0.30, 0.65, -1, 0.0),
             (0.70, 0.35, +1, 0.0), (0.70, 0.65, -1, 0.0),
             (0.50, 0.20, +1, 60.0), (0.50, 0.80, -1, 60.0),
             (0.20, 0.50, +1, 120.0), (0.80, 0.50, -1, 120.0)]
    m = PFC2D(N, N, r=-0.25, psi_bar=-0.25)
    m.init_dislocations(cores)
    return m


def main():
    os.makedirs(OUT, exist_ok=True)

    # immediately after seeding (pre-relax) -> expect ~8 cores
    m = build()
    pre = diag(m)
    print(f"seeded (pre-relax): n5={pre[0]} n7={pre[1]} cores={pre[2]} "
          f"(seeded 8)", flush=True)

    m.step_mpfc(DT, n=RELAX0, beta=10.0)
    post = diag(m)
    print(f"after relax (no shear): n5={post[0]} n7={post[1]} cores={post[2]}",
          flush=True)

    # ---- CONTROL: no-drive baseline, same number of relax blocks as driven ----
    mc = build()
    mc.step_mpfc(DT, n=RELAX0, beta=10.0)
    ctrl = [diag(mc)[2]]
    for _ in range(N_SHEAR):
        mc.step_mpfc(DT, n=RELAX, beta=10.0)  # NO shear
        ctrl.append(diag(mc)[2])
    print(f"CONTROL core count (no drive): {ctrl}", flush=True)

    # ---- DRIVEN: identical, but with x-shear each block ----
    driven = [post[2]]
    n5s = [post[0]]
    n7s = [post[1]]
    for _ in range(N_SHEAR):
        m.apply_shear(DGAMMA)
        m.step_mpfc(DT, n=RELAX, beta=10.0)
        a, b, c = diag(m)
        n5s.append(a)
        n7s.append(b)
        driven.append(c)
    print(f"DRIVEN core count (8% shear): {driven}", flush=True)

    summary = dict(
        seeded_pre_relax_cores=pre[2],
        after_relax_cores=post[2],
        control_core_series=ctrl,
        driven_core_series=driven,
        driven_n5_series=n5s, driven_n7_series=n7s,
        control_min=int(min(ctrl)), control_final=ctrl[-1],
        driven_min=int(min(driven)), driven_final=driven[-1],
        shear_pct=DGAMMA*N_SHEAR*100,
    )
    with open(os.path.join(OUT, "diag_forest_count.json"), "w") as f:
        json.dump(summary, f, indent=1)

    print("\n==== COUNT DIAGNOSTIC SUMMARY ====", flush=True)
    print(f"seeded 8 -> after relax {post[2]} cores "
          f"(loss during relaxation alone: {pre[2]-post[2]} from pre-relax "
          f"detect {pre[2]})", flush=True)
    print(f"CONTROL (no drive): start={ctrl[0]} min={min(ctrl)} end={ctrl[-1]}",
          flush=True)
    print(f"DRIVEN  (8% shear): start={driven[0]} min={min(driven)} "
          f"end={driven[-1]}", flush=True)


if __name__ == "__main__":
    main()

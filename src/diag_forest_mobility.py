"""DIAGNOSTIC: are the 'forest' dislocations in 2D triangular PFC immobile
obstacles (sessile) or mobile/annihilating soft modes?

This directly tests diagnosis (a) of the forest-softening result: forest
hardening requires the mobile dislocation to cut through SESSILE forest
junctions. If the forest cores themselves glide/climb or annihilate under the
applied shear, they cannot act as obstacles and the model softens.

Setup (256^2, r=-0.25, triangular):
  - 2 MOBILE dipoles: Burgers along x (0 deg) -> max Schmid factor under x-shear
  - FOREST: 2 dipoles on the 60deg and 120deg glide systems (the canonical
    'intersecting slip system' forest, low resolved shear under x-shear)
Relax, then apply incremental x-shear. At each strain:
  - locate all dislocation cores (5|7 pairs)
  - MATCH each forest core to its nearest core from the previous frame and
    accumulate per-core path length (how far it moved) -> mobility
  - record total core count -> annihilation

Reports measured numbers: forest-core displacement vs strain, and core-count
trajectory. Verdict: sessile obstacle vs mobile/annihilating soft mode.
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
RELAX0 = 700        # initial relaxation of the seeded structure
DGAMMA = 0.0025
N_SHEAR = 32        # -> 8% engineering shear
RELAX = 150         # relaxation steps per shear increment
A0 = A_LATTICE      # lattice spacing in this nondimensional unit


def min_image(d, L):
    """Wrap a coordinate difference into [-L/2, L/2)."""
    return d - L * np.round(d / L)


def core_distance(a, b, lx, ly):
    dx = min_image(a[0] - b[0], lx)
    dy = min_image(a[1] - b[1], ly)
    return float(np.hypot(dx, dy))


def detect_cores(m):
    pts = find_peaks(m.psi, m.dx, m.dy)
    d = find_dislocations(pts, m.lx, m.ly)
    return d["cores"], d


def match_tracks(prev, cur, lx, ly, max_jump):
    """Greedy nearest-neighbour matching of cores between two frames.
    Returns list of (i_prev, j_cur, dist) for matched pairs and the set of
    unmatched indices on each side. max_jump caps a plausible per-step move."""
    matches = []
    used_cur = set()
    if len(prev) == 0 or len(cur) == 0:
        return matches, set(range(len(prev))), set(range(len(cur)))
    # build all pairwise periodic distances
    for i, p in enumerate(prev):
        best_j, best_d = -1, 1e30
        for j, c in enumerate(cur):
            if j in used_cur:
                continue
            dd = core_distance(p, c, lx, ly)
            if dd < best_d:
                best_d, best_j = dd, j
        if best_j >= 0 and best_d <= max_jump:
            used_cur.add(best_j)
            matches.append((i, best_j, best_d))
    matched_prev = {i for i, _, _ in matches}
    unmatched_prev = set(range(len(prev))) - matched_prev
    unmatched_cur = set(range(len(cur))) - used_cur
    return matches, unmatched_prev, unmatched_cur


def main():
    os.makedirs(OUT, exist_ok=True)

    # ---- seed a small fixed set: 2 mobile (0deg) + forest (60/120deg) ----
    # mobile dipoles, Burgers along x, separated in y (glide plane horizontal)
    cores = [
        (0.30, 0.35, +1, 0.0), (0.30, 0.65, -1, 0.0),   # mobile dipole A
        (0.70, 0.35, +1, 0.0), (0.70, 0.65, -1, 0.0),   # mobile dipole B
        # FOREST: 60deg system dipole
        (0.50, 0.20, +1, 60.0), (0.50, 0.80, -1, 60.0),
        # FOREST: 120deg system dipole
        (0.20, 0.50, +1, 120.0), (0.80, 0.50, -1, 120.0),
    ]
    # tag which seeded cores are forest (60/120) for reporting the seeded
    # forest fractional positions
    forest_seed_fracs = [(0.50, 0.20), (0.50, 0.80), (0.20, 0.50), (0.80, 0.50)]

    m = PFC2D(N, N, r=-0.25, psi_bar=-0.25)
    m.init_dislocations(cores)
    # relax the seeded structure (no shear yet) using MPFC dynamics, same as
    # the forest-hardening production runs
    m.step_mpfc(DT, n=RELAX0, beta=10.0)

    lx, ly = m.lx, m.ly
    a0_phys = A0  # nondim spacing reference; max plausible per-increment jump
    max_jump = 3.0 * a0_phys  # cores shouldn't teleport > ~3 lattice spacings

    cur_cores, d0 = detect_cores(m)
    print(f"after relax: {len(cur_cores)} cores (n5={d0['n5']}, n7={d0['n7']})",
          flush=True)

    # identify which detected cores are the forest ones: nearest detected core
    # to each seeded forest fractional position
    forest_idx = []
    for fx, fy in forest_seed_fracs:
        target = np.array([fx * lx, fy * ly])
        if len(cur_cores) == 0:
            break
        dists = [core_distance(target, c, lx, ly) for c in cur_cores]
        j = int(np.argmin(dists))
        if dists[j] < 6.0 * a0_phys and j not in forest_idx:
            forest_idx.append(j)
    forest_anchor = [cur_cores[j].copy() for j in forest_idx]
    print(f"matched {len(forest_idx)} forest cores at relax", flush=True)

    # per-forest-core accumulated path length and net displacement from anchor
    n_forest = len(forest_anchor)
    path_len = np.zeros(n_forest)         # integrated |displacement| (mobility)
    track_pos = [a.copy() for a in forest_anchor]
    track_alive = [True] * n_forest

    history = []
    history.append(dict(gamma=float(m.gamma), n_cores=int(len(cur_cores)),
                        forest_path=path_len.tolist(),
                        forest_net=[0.0] * n_forest,
                        forest_alive=list(track_alive)))

    prev_cores = cur_cores
    for it in range(N_SHEAR):
        m.apply_shear(DGAMMA)
        m.step_mpfc(DT, n=RELAX, beta=10.0)
        cur_cores, dd = detect_cores(m)

        # advance each tracked forest core by nearest-match in current frame
        for k in range(n_forest):
            if not track_alive[k]:
                continue
            if len(cur_cores) == 0:
                track_alive[k] = False
                continue
            dists = [core_distance(track_pos[k], c, lx, ly) for c in cur_cores]
            j = int(np.argmin(dists))
            if dists[j] <= max_jump:
                step = dists[j]
                path_len[k] += step
                track_pos[k] = cur_cores[j].copy()
            else:
                # no nearby core -> this forest core annihilated / vanished
                track_alive[k] = False

        forest_net = [core_distance(track_pos[k], forest_anchor[k], lx, ly)
                      if track_alive[k] else np.nan for k in range(n_forest)]
        history.append(dict(gamma=float(m.gamma),
                            n_cores=int(len(cur_cores)),
                            forest_path=path_len.tolist(),
                            forest_net=forest_net,
                            forest_alive=list(track_alive)))
        if it % 4 == 0 or it == N_SHEAR - 1:
            alive = sum(track_alive)
            mean_path = np.mean(path_len[:n_forest]) if n_forest else 0.0
            print(f"  gamma={m.gamma:.3f}: total_cores={len(cur_cores)} "
                  f"forest_alive={alive}/{n_forest} "
                  f"mean_forest_path={mean_path/a0_phys:.2f} a0", flush=True)
        prev_cores = cur_cores

    # ---- summary metrics ----
    n_cores_series = [h["n_cores"] for h in history]
    final = history[-1]
    final_paths_a0 = [p / a0_phys for p in final["forest_path"]]
    final_net_a0 = [(v / a0_phys if np.isfinite(v) else None)
                    for v in final["forest_net"]]
    alive_final = sum(final["forest_alive"])

    # mobile-vs-forest contrast: total core count is the cleanest annihilation
    # signal; forest path length is the cleanest mobility signal
    summary = dict(
        N=N, r=-0.25, dgamma=DGAMMA, n_shear=N_SHEAR, relax=RELAX,
        n_cores_initial=n_cores_series[0],
        n_cores_final=n_cores_series[-1],
        n_cores_min=int(min(n_cores_series)),
        n_cores_series=n_cores_series,
        n_forest_tracked=n_forest,
        forest_alive_final=alive_final,
        forest_final_path_len_a0=final_paths_a0,
        forest_final_net_disp_a0=final_net_a0,
        mean_forest_path_a0=float(np.mean(final_paths_a0)) if n_forest else 0.0,
        max_forest_path_a0=float(np.max(final_paths_a0)) if n_forest else 0.0,
    )
    with open(os.path.join(OUT, "diag_forest_mobility.json"), "w") as f:
        json.dump(dict(summary=summary, history=history), f, indent=1)

    print("\n==== DIAGNOSTIC SUMMARY ====", flush=True)
    print(f"cores: initial={summary['n_cores_initial']} "
          f"min={summary['n_cores_min']} final={summary['n_cores_final']}",
          flush=True)
    print(f"forest cores tracked={n_forest} alive_at_end={alive_final}",
          flush=True)
    print(f"forest per-core path length (a0): "
          f"{[round(x,2) for x in final_paths_a0]}", flush=True)
    print(f"forest per-core NET displacement (a0): "
          f"{[None if x is None else round(x,2) for x in final_net_a0]}",
          flush=True)
    print(f"mean forest path = {summary['mean_forest_path_a0']:.2f} a0, "
          f"max = {summary['max_forest_path_a0']:.2f} a0 "
          f"over {DGAMMA*N_SHEAR*100:.0f}% shear", flush=True)


if __name__ == "__main__":
    main()

"""DIAGNOSTIC v2: refined mobility/annihilation test with controls.

Adds to v1:
  (1) a NOISE-FLOOR control: how much do detected core positions jitter from
      detection alone, on a relaxed structure held at FIXED shear (no driving)?
      Subtract this to know real motion.
  (2) cleaner separation: track BOTH the mobile (0deg) and forest (60/120deg)
      cores so the contrast (mobile should move a lot; forest is the question)
      is explicit.
  (3) robust matching with a per-step cap and reporting of net vs path so a
      core that wanders back is not counted as 'pinned'.

256^2, r=-0.25, MPFC dynamics, 8% x-shear.
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
A0 = A_LATTICE


def mi(d, L):
    return d - L * np.round(d / L)


def cdist(a, b, lx, ly):
    return float(np.hypot(mi(a[0]-b[0], lx), mi(a[1]-b[1], ly)))


def cores_of(m):
    pts = find_peaks(m.psi, m.dx, m.dy)
    d = find_dislocations(pts, m.lx, m.ly)
    return d["cores"]


def track_one(start, frames, lx, ly, max_jump):
    """Follow a single core from `start` through a list of core-arrays.
    Returns (path_len, net_disp, alive, final_pos)."""
    pos = np.array(start, float)
    path = 0.0
    alive = True
    for cur in frames:
        if not alive or len(cur) == 0:
            alive = False
            continue
        dists = [cdist(pos, c, lx, ly) for c in cur]
        j = int(np.argmin(dists))
        if dists[j] <= max_jump:
            path += dists[j]
            pos = np.array(cur[j], float)
        else:
            alive = False
    net = cdist(pos, np.array(start, float), lx, ly) if alive else np.nan
    return path, net, alive, pos


def main():
    os.makedirs(OUT, exist_ok=True)

    seeded = dict(
        mobile=[(0.30, 0.35, +1, 0.0), (0.30, 0.65, -1, 0.0),
                (0.70, 0.35, +1, 0.0), (0.70, 0.65, -1, 0.0)],
        forest=[(0.50, 0.20, +1, 60.0), (0.50, 0.80, -1, 60.0),
                (0.20, 0.50, +1, 120.0), (0.80, 0.50, -1, 120.0)],
    )
    cores = seeded["mobile"] + seeded["forest"]
    m = PFC2D(N, N, r=-0.25, psi_bar=-0.25)
    m.init_dislocations(cores)
    m.step_mpfc(DT, n=RELAX0, beta=10.0)
    lx, ly = m.lx, m.ly
    max_jump = 3.0 * A0

    c0 = cores_of(m)
    print(f"after relax: {len(c0)} cores detected", flush=True)

    # ---- (1) NOISE FLOOR: hold shear fixed, relax more, measure jitter ----
    floor_frames = []
    mtmp_psi = m.psi.copy()
    for _ in range(8):
        m.step_mpfc(DT, n=RELAX, beta=10.0)   # same relax cadence, NO shear
        floor_frames.append(cores_of(m))
    # track every detected core through the no-drive frames
    floor_paths = []
    for c in c0:
        p, _, alive, _ = track_one(c, floor_frames, lx, ly, max_jump)
        if alive:
            floor_paths.append(p / A0)
    noise_floor = float(np.mean(floor_paths)) if floor_paths else 0.0
    noise_floor_max = float(np.max(floor_paths)) if floor_paths else 0.0
    print(f"NOISE FLOOR (no drive, 8 frames): mean path={noise_floor:.2f} a0, "
          f"max={noise_floor_max:.2f} a0 over comparable relaxation", flush=True)

    # restore the pre-floor field so the driven run starts clean
    m.psi = mtmp_psi
    m._psi_prev = mtmp_psi.copy()
    m._update_k()

    # anchors for mobile + forest from current detection
    def nearest_detected(fx, fy, used):
        target = np.array([fx*lx, fy*ly])
        ds = [cdist(target, c, lx, ly) for c in c0]
        order = np.argsort(ds)
        for j in order:
            if j not in used and ds[j] < 6.0*A0:
                used.add(int(j))
                return c0[j].copy()
        return None

    used = set()
    mob_anchors, for_anchors = [], []
    for (fx, fy, *_ ) in seeded["mobile"]:
        a = nearest_detected(fx, fy, used)
        if a is not None:
            mob_anchors.append(a)
    for (fx, fy, *_ ) in seeded["forest"]:
        a = nearest_detected(fx, fy, used)
        if a is not None:
            for_anchors.append(a)
    print(f"anchored {len(mob_anchors)} mobile + {len(for_anchors)} forest "
          f"cores", flush=True)

    # ---- driven shear, collect frames ----
    frames = []
    n_series = [len(c0)]
    for it in range(N_SHEAR):
        m.apply_shear(DGAMMA)
        m.step_mpfc(DT, n=RELAX, beta=10.0)
        cc = cores_of(m)
        frames.append(cc)
        n_series.append(len(cc))

    # track mobile + forest anchors through driven frames
    def summarize(anchors):
        res = []
        for a in anchors:
            p, net, alive, pos = track_one(a, frames, lx, ly, max_jump)
            res.append(dict(path_a0=p/A0,
                            net_a0=(net/A0 if np.isfinite(net) else None),
                            alive=alive))
        return res

    mob = summarize(mob_anchors)
    forr = summarize(for_anchors)

    def stat(rs, key):
        v = [r[key] for r in rs if r[key] is not None]
        return (float(np.mean(v)), float(np.max(v))) if v else (0.0, 0.0)

    mob_path = stat(mob, "path_a0")
    for_path = stat(forr, "path_a0")
    mob_net = stat(mob, "net_a0")
    for_net = stat(forr, "net_a0")

    summary = dict(
        noise_floor_mean_a0=noise_floor, noise_floor_max_a0=noise_floor_max,
        n_cores_initial=n_series[0], n_cores_final=n_series[-1],
        n_cores_min=int(min(n_series)), n_cores_series=n_series,
        mobile=dict(n=len(mob), path_mean_a0=mob_path[0], path_max_a0=mob_path[1],
                    net_mean_a0=mob_net[0], net_max_a0=mob_net[1], detail=mob),
        forest=dict(n=len(forr), path_mean_a0=for_path[0], path_max_a0=for_path[1],
                    net_mean_a0=for_net[0], net_max_a0=for_net[1], detail=forr),
        shear_pct=DGAMMA*N_SHEAR*100,
    )
    with open(os.path.join(OUT, "diag_forest_mobility2.json"), "w") as f:
        json.dump(summary, f, indent=1)

    print("\n==== DIAGNOSTIC v2 SUMMARY ====", flush=True)
    print(f"shear applied: {summary['shear_pct']:.0f}%", flush=True)
    print(f"NOISE FLOOR (detection jitter, no drive): "
          f"mean={noise_floor:.2f} a0  max={noise_floor_max:.2f} a0", flush=True)
    print(f"core count: init={n_series[0]} min={min(n_series)} "
          f"final={n_series[-1]}  series={n_series}", flush=True)
    print(f"MOBILE (0deg): path mean={mob_path[0]:.2f} max={mob_path[1]:.2f} a0 | "
          f"net mean={mob_net[0]:.2f} max={mob_net[1]:.2f} a0", flush=True)
    print(f"FOREST (60/120deg): path mean={for_path[0]:.2f} max={for_path[1]:.2f} a0 | "
          f"net mean={for_net[0]:.2f} max={for_net[1]:.2f} a0", flush=True)
    print(f"forest net displacement vs noise floor: "
          f"{for_net[0]:.2f} a0 vs {noise_floor:.2f} a0 jitter", flush=True)


if __name__ == "__main__":
    main()

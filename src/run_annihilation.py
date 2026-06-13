"""Loop-closure (redirected, dimensionless): measure the dislocation
ANNIHILATION DISTANCE y_ann from PFC — the critical glide-plane separation
below which two opposite-sign edge dislocations spontaneously annihilate.

y_ann (in units of b) is an explicit recovery parameter in Kocks-Mecking /
DAMASK-dislotwin density-evolution laws (drho/dgamma has a -2 y_ann rho v / b
annihilation term). It is geometric and dimensionless, so PFC supplies it to
crystal plasticity WITHOUT any PFC-time<->seconds calibration — a clean
upscaling transfer.

Method: seed a dipole with cores on two glide planes separated by d (along y),
at the SAME x (so the glide-plane offset is the annihilation-relevant
distance). Relax with NO applied stress; if the pair annihilates (cores -> 0)
the separation was within y_ann. Bisect over d to find the critical distance.
Output: results/annihilation/
"""

import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "annihilation")
N = 256
DT = 0.5


def annihilates(d_over_b, r=-0.25, relax=1500, seed=0):
    """Seed an edge dipole with glide-plane (y) separation d = d_over_b * b at
    the same x; relax unstressed; return True if it annihilates."""
    b = A_LATTICE
    m = PFC2D(N, N, r=r, psi_bar=-0.25)
    dy = d_over_b * b / m.ly
    yc = 0.5
    # same x, opposite sign, separated by dy along the glide-plane normal
    m.init_dislocations([(0.5, yc - dy / 2, +1), (0.5, yc + dy / 2, -1)])
    m.step(DT, n=relax)
    d = find_dislocations(find_peaks(m.psi, m.dx, m.dy), m.lx, m.ly)
    return len(d["cores"]) == 0, len(d["cores"])


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    b = A_LATTICE
    # scan separation in units of b
    scan = []
    for dob in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
        ann, nc = annihilates(dob)
        scan.append(dict(d_over_b=dob, annihilated=ann, n_cores=nc))
        print(f"  d/b={dob:.1f}: annihilated={ann} (cores={nc})", flush=True)
    # critical y_ann = midpoint between largest annihilating and smallest surviving
    ann_d = [s["d_over_b"] for s in scan if s["annihilated"]]
    surv_d = [s["d_over_b"] for s in scan if not s["annihilated"]]
    if ann_d and surv_d:
        y_ann = 0.5 * (max(ann_d) + min(d for d in surv_d if d > max(ann_d)))
    elif ann_d:
        y_ann = max(ann_d)
    else:
        y_ann = float("nan")
    result = dict(y_ann_over_b=float(y_ann), scan=scan, b=b,
                  note="annihilation distance in units of Burgers vector; "
                       "Kocks-Mecking/dislotwin recovery parameter "
                       "(typical metals y_ann ~ 1-6 b). Dimensionless -> "
                       "transferable to CP with no time calibration.")
    with open(os.path.join(OUT, "annihilation.json"), "w") as f:
        json.dump(result, f, indent=1)
    print(f"\nANNIHILATION DISTANCE y_ann = {y_ann:.1f} b "
          f"(metals: ~1-6 b) -- transferable CP recovery parameter", flush=True)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

"""Stage-1 cheapest decisive test: is the forest-hardening blocker DYNAMIC
(climb) or CRYSTALLOGRAPHIC?  Binary two-line junction, climb-suppressed vs
climb-active, in the EXISTING 3D BCC code — no XPFC build.

Two intersecting dislocation lines (a mobile z-line, Burgers x, gliding under
xz shear; a forest y-line threading its glide plane, Burgers z). Under shear
they meet and may react. We probe junction survival / glide resistance in two
kinetic limits, both reaching the SAME accumulated shear:
  - CLIMB-SUPPRESSED proxy: high beta (fast elastic) + FEW MPFC sub-steps per
    shear increment -> little diffusive (climb) advance between increments
  - CLIMB-ACTIVE: many sub-steps per increment -> climb fully operative
DECISION: if the junction stays put / raises tau and the forest line stays
sessile in the climb-suppressed limit but unzips/annihilates and softens in
the climb-active limit -> dynamics is the lever (-> Stage 2). If bypassed in
both -> crystallography also blocks (-> XPFC/two-mode co-required).

Output: results/junction_test/
"""

import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc3d import PFC3D, A_BCC, find_peaks_3d
from defect_analysis_3d import find_dislocation_lines

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "junction_test")
N = 96
CELLS = 9
DX = CELLS * A_BCC / N
DT = 0.5
DG = 0.0025
TOTAL_GAMMA = 0.06
BETA = 10.0


def forest_displacement(m, y_axis_x0, y_axis_z0):
    """Track the threading y-line's core position drift (proxy for whether it
    stays sessile). Returns mean |disorder-centroid drift| is hard in 3D, so we
    report the disordered-atom count near the original forest line plane."""
    pts = find_peaks_3d(m.psi, m.dx, m.dy, m.dz)
    box = np.array([m.lx, m.ly, m.lz])
    r = find_dislocation_lines(pts, box)
    return r["disordered_frac"], r["n_lines"], len(pts)


def run(label, sub_steps):
    """sub_steps = MPFC sub-steps per shear increment (few=climb-suppressed)."""
    m = PFC3D(N, N, N, dx=DX, r=-0.25, psi_bar=-0.25)
    # mobile z-line dipole (Burgers x) + forest y-line dipole (Burgers z)
    # threading the glide plane; net Burgers per component cancels within pairs.
    m.init_dislocation_lines([
        dict(axis="z", pos=(0.5, 0.30), burgers="x", sign=+1),
        dict(axis="z", pos=(0.5, 0.70), burgers="x", sign=-1),
        dict(axis="y", pos=(0.30, 0.5), burgers="z", sign=+1),
        dict(axis="y", pos=(0.70, 0.5), burgers="z", sign=-1),
    ])
    m.step_mpfc(DT, n=500, beta=BETA)
    f0, nl0, np0 = forest_displacement(m, 0, 0)
    rows = [dict(gxz=0.0, tau=m.shear_stress(), frac=f0, n_lines=nl0)]
    n_inc = int(round(TOTAL_GAMMA / DG))
    for i in range(n_inc):
        m.apply_shear(DG)
        m.step_mpfc(DT, n=sub_steps, beta=BETA)
        if i % 3 == 2 or i == n_inc - 1:
            f, nl, npp = forest_displacement(m, 0, 0)
            rows.append(dict(gxz=m.gxz, tau=m.shear_stress(), frac=f,
                             n_lines=nl, atoms=npp))
    tau_flow = float(np.mean([r["tau"] for r in rows[-3:]]))
    return dict(label=label, sub_steps=sub_steps, rows=rows,
                tau_flow=tau_flow, frac0=f0, frac_end=rows[-1]["frac"],
                nlines0=nl0, nlines_end=rows[-1]["n_lines"])


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    res = {}
    for label, sub in [("climb_suppressed", 20), ("climb_active", 300)]:
        r = run(label, sub)
        res[label] = r
        print(f"[{label}] sub_steps={sub}: tau_flow={r['tau_flow']:.5f} "
              f"frac {r['frac0']:.3f}->{r['frac_end']:.3f} "
              f"lines {r['nlines0']}->{r['nlines_end']}", flush=True)
    # decision proxy: does climb-suppressed retain more defect content (junction
    # survives) and higher tau than climb-active?
    cs, ca = res["climb_suppressed"], res["climb_active"]
    retains = cs["frac_end"] > ca["frac_end"] + 0.02
    harder = cs["tau_flow"] > ca["tau_flow"] * 1.05
    res["verdict"] = (
        "DYNAMICS IS THE LEVER: climb-suppressed retains junction/defects "
        "and/or is harder -> proceed to Stage 2 (climb-suppressed + two-mode FCC)"
        if (retains or harder) else
        "climb suppression alone insufficient -> crystallography also blocks "
        "(XPFC/two-mode co-required); OR junction never formed (check seeding)")
    with open(os.path.join(OUT, "junction_test.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(f"\nVERDICT: {res['verdict']}")
    print(f"  climb-suppressed: tau_flow={cs['tau_flow']:.5f} "
          f"frac_end={cs['frac_end']:.3f}")
    print(f"  climb-active:     tau_flow={ca['tau_flow']:.5f} "
          f"frac_end={ca['frac_end']:.3f}")
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

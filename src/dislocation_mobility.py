"""Sub-problem A foundation: is PFC a QUANTITATIVE mesoscale mechanics model?

Two validations that gate any constitutive coupling to DDD/CP:
  (1) ELASTIC FIELD: does a PFC edge dislocation produce the correct
      long-range stress field?  Proxy: the lattice-distortion (local
      shear) around an isolated core should decay as ~b/(2 pi r) (isotropic
      edge-dislocation field). We measure the azimuthally-averaged distortion
      vs r and fit the 1/r law.
  (2) MOBILITY v(tau): under a controlled resolved shear stress, an edge
      dislocation glides at velocity v. The Peach-Koehler force per length is
      f = tau * b, so v = M * tau * b defines the PFC glide mobility M — the
      transferable quantity DDD/ExaDiS needs (v = M_DDD * tau).

Method: seed a wide edge dipole (two opposite edges separated by ly/2). Apply
simple shear gamma in steps; the resolved shear drives the two cores to glide
toward/away along x. Track core x-positions vs accumulated strain-time and
fit v(tau). Both diffusive (step) and MPFC (step_mpfc) dynamics are run so we
can report mobility in each and flag the dynamics dependence the route review
demanded.

Output: results/mobility/{elastic_field.json, mobility.json, mobility.png}
"""

import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "mobility")


def seeded_dipole(n=256, r=-0.25, relax=600, mpfc=False, beta=10.0):
    m = PFC2D(n, n, r=r, psi_bar=-0.25)
    # horizontal edge dipole: cores at (lx/4, ly/2) and (3lx/4, ly/2),
    # opposite sign -> a glide-apart pair under +x resolved shear
    m.init_dislocations([(0.25, 0.5, +1), (0.75, 0.5, -1)])
    (m.step_mpfc(0.5, n=relax, beta=beta) if mpfc else m.step(0.5, n=relax))
    return m


def core_positions(m):
    d = find_dislocations(find_peaks(m.psi, m.dx, m.dy), m.lx, m.ly)
    return d["cores"]


def measure_mobility(mpfc, label):
    m = seeded_dipole(mpfc=mpfc)
    c0 = core_positions(m)
    if len(c0) < 2:
        return dict(label=label, ok=False, note="dipole not stable")
    track = []
    dt = 0.5
    relax = 200
    for i in range(20):
        m.apply_shear(0.002)            # +x resolved shear increment
        (m.step_mpfc(dt, n=relax, beta=10.0) if mpfc else m.step(dt, n=relax))
        cores = core_positions(m)
        tau = m.shear_stress()
        # mean |x| glide of the cores from start (min-image)
        if len(cores) >= 2:
            xs = np.sort(cores[:, 0])
            sep = xs[-1] - xs[0]
            track.append((m.gamma, tau, sep, m.time, len(cores)))
    track = np.array(track)
    # glide velocity: d(sep)/d(time) once moving; mobility M = v/(tau*b)
    if len(track) < 5:
        return dict(label=label, ok=False, note="lost cores")
    b = A_LATTICE
    v = np.gradient(track[:, 2], track[:, 3])      # d sep / d time
    tau = track[:, 1]
    moving = np.abs(v) > 1e-4
    M = float(np.median(v[moving] / (tau[moving] * b))) if moving.any() else 0.0
    return dict(label=label, ok=True, mobility_M=M,
                gamma=track[:, 0].tolist(), tau=tau.tolist(),
                sep=track[:, 2].tolist(), v=v.tolist(),
                n_cores=track[:, 4].tolist())


def elastic_field():
    """Azimuthal distortion decay around a single core (1/r test).
    Use one core of a widely-separated dipole, sample local lattice shear
    (gradient of the phase) in annuli and fit amplitude*(1/r)+c."""
    m = seeded_dipole(relax=800)
    cores = core_positions(m)
    if len(cores) < 1:
        return dict(ok=False)
    cx, cy = cores[0]
    # local "distortion" proxy: magnitude of grad(psi) deviation from bulk
    gy, gx = np.gradient(m.psi, m.dy, m.dx)
    gmag = np.sqrt(gx ** 2 + gy ** 2)
    x = (np.arange(m.nx) * m.dx)[None, :]
    y = (np.arange(m.ny) * m.dy)[:, None]
    rr = np.sqrt(((x - cx + m.lx / 2) % m.lx - m.lx / 2) ** 2
                 + ((y - cy + m.ly / 2) % m.ly - m.ly / 2) ** 2)
    bins = np.linspace(2 * A_LATTICE, 0.4 * m.lx, 12)
    prof = []
    bulk = np.median(gmag)
    for lo, hi in zip(bins[:-1], bins[1:]):
        sel = (rr >= lo) & (rr < hi)
        if sel.sum():
            prof.append((0.5 * (lo + hi), float(np.mean(gmag[sel]) - bulk)))
    prof = np.array(prof)
    # fit excess distortion ~ A/r
    pos = prof[:, 1] > 0
    if pos.sum() >= 4:
        A, c = np.polyfit(1.0 / prof[pos, 0], prof[pos, 1], 1)
        r2 = 1 - (np.var(prof[pos, 1] - (A / prof[pos, 0] + c))
                  / (np.var(prof[pos, 1]) + 1e-12))
    else:
        A, c, r2 = float("nan"), float("nan"), float("nan")
    return dict(ok=True, r=prof[:, 0].tolist(), distortion=prof[:, 1].tolist(),
                inv_r_amp=float(A), r2_of_1overr_fit=float(r2))


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    ef = elastic_field()
    with open(os.path.join(OUT, "elastic_field.json"), "w") as f:
        json.dump(ef, f, indent=1)
    print(f"[elastic field] 1/r fit R^2 = {ef.get('r2_of_1overr_fit')}", flush=True)

    res = {}
    for mpfc, label in [(False, "diffusive"), (True, "mpfc_beta10")]:
        r = measure_mobility(mpfc, label)
        res[label] = r
        print(f"[{label}] mobility M = {r.get('mobility_M')} ok={r['ok']}",
              flush=True)
    with open(os.path.join(OUT, "mobility.json"), "w") as f:
        json.dump(res, f, indent=1)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        if ef.get("ok"):
            r = np.array(ef["r"]); dd = np.array(ef["distortion"])
            ax[0].plot(r, dd, "o", label="measured excess distortion")
            A = ef["inv_r_amp"]
            ax[0].plot(r, A / r, "-", label=f"A/r fit (R2={ef['r2_of_1overr_fit']:.2f})")
            ax[0].set_xlabel("r from core"); ax[0].set_ylabel("excess distortion")
            ax[0].set_title("edge dislocation elastic field (1/r test)")
            ax[0].legend(); ax[0].grid(alpha=0.3)
        for label, r in res.items():
            if r["ok"]:
                ax[1].plot(np.array(r["tau"]), np.array(r["v"]), "o-",
                           label=f"{label} (M={r['mobility_M']:.3f})")
        ax[1].set_xlabel("resolved shear stress tau")
        ax[1].set_ylabel("glide velocity v")
        ax[1].set_title("dislocation mobility v(tau)")
        ax[1].legend(); ax[1].grid(alpha=0.3)
        fig.savefig(os.path.join(OUT, "mobility.png"), dpi=140,
                    bbox_inches="tight")
    except Exception as ex:
        print("plot skipped:", ex)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

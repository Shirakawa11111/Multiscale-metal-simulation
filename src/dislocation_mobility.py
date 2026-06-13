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
    # Edge dipole separated along Y (perpendicular to the x-glide direction),
    # so the two opposite-sign edges sit on DIFFERENT horizontal glide planes
    # (y=0.3 ly, y=0.7 ly). Under applied epsilon_xy shear they glide in
    # OPPOSITE x-directions on their own planes; their mutual interaction at
    # large Delta-y is weak and mostly climb (slow), so the measured x-glide
    # velocity reflects the APPLIED resolved stress, not pair annihilation.
    m.init_dislocations([(0.5, 0.3, +1), (0.5, 0.7, -1)])
    (m.step_mpfc(0.5, n=relax, beta=beta) if mpfc else m.step(0.5, n=relax))
    return m


def core_positions(m):
    d = find_dislocations(find_peaks(m.psi, m.dx, m.dy), m.lx, m.ly)
    return d["cores"]


def _match_core_x(cores, y_target, ly):
    """x-position of the core nearest a given glide-plane y (periodic)."""
    if len(cores) == 0:
        return None
    dy = np.abs((cores[:, 1] - y_target + ly / 2) % ly - ly / 2)
    return float(cores[np.argmin(dy)][0])


def measure_mobility(mpfc, label, gamma_hold=0.03):
    """Apply a FIXED shear gamma_hold, then hold it and let the two edges glide
    at constant resolved stress; track each core's absolute x vs time and fit
    v. Mobility M = v / (tau * b)."""
    m = seeded_dipole(mpfc=mpfc)
    c0 = core_positions(m)
    if len(c0) < 2:
        return dict(label=label, ok=False, note="dipole not stable")
    y_lo, y_hi = 0.3 * m.ly, 0.7 * m.ly
    # ramp to a fixed shear quickly, then hold
    for _ in range(int(gamma_hold / 0.0025)):
        m.apply_shear(0.0025)
        (m.step_mpfc(0.5, n=40, beta=10.0) if mpfc else m.step(0.5, n=40))
    tau = m.shear_stress()
    # now hold gamma fixed and watch glide
    x_lo0 = _match_core_x(core_positions(m), y_lo, m.ly)
    x_hi0 = _match_core_x(core_positions(m), y_hi, m.ly)
    track = []
    t0 = m.time
    for i in range(24):
        (m.step_mpfc(0.5, n=60, beta=10.0) if mpfc else m.step(0.5, n=60))
        cores = core_positions(m)
        xl = _match_core_x(cores, y_lo, m.ly)
        xh = _match_core_x(cores, y_hi, m.ly)
        if xl is None or xh is None:
            break
        # unwrapped glide distance of each core from its start (min-image)
        dlo = ((xl - x_lo0 + m.lx / 2) % m.lx - m.lx / 2)
        dhi = ((xh - x_hi0 + m.lx / 2) % m.lx - m.lx / 2)
        glide = 0.5 * (abs(dlo) + abs(dhi))   # mean |glide| of the two
        track.append((m.time - t0, glide, len(cores)))
    track = np.array(track)
    if len(track) < 6:
        return dict(label=label, ok=False, note="lost cores", tau=float(tau))
    b = A_LATTICE
    # steady glide velocity = slope of glide-distance vs time (robust fit)
    v = float(np.polyfit(track[:, 0], track[:, 1], 1)[0])
    M = v / (tau * b) if tau != 0 else float("nan")
    return dict(label=label, ok=True, mobility_M=float(M), v=v,
                tau=float(tau), gamma_hold=gamma_hold,
                t=track[:, 0].tolist(), glide=track[:, 1].tolist(),
                n_cores=track[:, 2].tolist())


def _free_energy_density(m):
    from pfc2d import _rfft2, _irfft2
    psi_h = _rfft2(m.psi)
    lin_term = _irfft2(m.lin * psi_h, m.psi.shape)
    return 0.5 * m.psi * lin_term + 0.25 * m.psi ** 4


def elastic_field():
    """Elastic-energy decay around a core (1/r^2 test). A dislocation's elastic
    strain ~ b/(2 pi r), so its excess free-energy DENSITY ~ strain^2 ~ 1/r^2.
    We coarse-grain the local f over a unit cell, take the azimuthal mean of the
    excess-over-bulk vs r, and fit A/r^2. (Cleaner & more physical than a
    grad-psi proxy.)"""
    from scipy.ndimage import uniform_filter
    m = seeded_dipole(relax=800)
    cores = core_positions(m)
    if len(cores) < 1:
        return dict(ok=False)
    cx, cy = cores[0]
    f = _free_energy_density(m)
    w = max(3, int(round(A_LATTICE / m.dx)))      # coarse-grain over a cell
    fcg = uniform_filter(f, size=w, mode="wrap")
    bulk = np.median(fcg)
    x = (np.arange(m.nx) * m.dx)[None, :]
    y = (np.arange(m.ny) * m.dy)[:, None]
    rr = np.sqrt(((x - cx + m.lx / 2) % m.lx - m.lx / 2) ** 2
                 + ((y - cy + m.ly / 2) % m.ly - m.ly / 2) ** 2)
    bins = np.linspace(1.5 * A_LATTICE, 0.35 * m.lx, 12)
    prof = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        sel = (rr >= lo) & (rr < hi)
        if sel.sum() > 5:
            prof.append((0.5 * (lo + hi), float(np.mean(fcg[sel]) - bulk)))
    prof = np.array(prof)
    pos = prof[:, 1] > 0
    if pos.sum() >= 4:
        rp, ep = prof[pos, 0], prof[pos, 1]
        # fit log(excess) = log A - p log r ; report exponent p (elastic -> 2)
        p_fit, logA = np.polyfit(np.log(rp), np.log(ep), 1)
        pexp = -p_fit
        pred = np.exp(logA) * rp ** p_fit
        r2 = 1 - np.var(ep - pred) / (np.var(ep) + 1e-12)
    else:
        pexp, r2 = float("nan"), float("nan")
    return dict(ok=True, r=prof[:, 0].tolist(),
                excess_energy=prof[:, 1].tolist(),
                decay_exponent_p=float(pexp), r2_of_powerlaw=float(r2),
                note="elastic dislocation -> p~2 (energy density ~ 1/r^2)")


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    ef = elastic_field()
    with open(os.path.join(OUT, "elastic_field.json"), "w") as f:
        json.dump(ef, f, indent=1)
    print(f"[elastic field] energy-density decay exponent p = "
          f"{ef.get('decay_exponent_p')} (elastic->2), R^2="
          f"{ef.get('r2_of_powerlaw')}", flush=True)

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
            r = np.array(ef["r"]); dd = np.array(ef["excess_energy"])
            ax[0].loglog(r, np.abs(dd), "o", label="excess energy density")
            p = ef["decay_exponent_p"]
            ax[0].set_xlabel("r from core"); ax[0].set_ylabel("excess f")
            ax[0].set_title(f"elastic field: energy~1/r^{p:.1f} "
                            f"(R2={ef['r2_of_powerlaw']:.2f}, elastic->2)")
            ax[0].legend(); ax[0].grid(alpha=0.3, which="both")
        for label, r in res.items():
            if r.get("ok"):
                ax[1].plot(np.array(r["t"]), np.array(r["glide"]), "o-",
                           label=f"{label}: M={r['mobility_M']:.3f}, tau={r['tau']:.4f}")
        ax[1].set_xlabel("time (held shear)")
        ax[1].set_ylabel("mean glide distance")
        ax[1].set_title("dislocation glide under constant resolved stress")
        ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
        fig.savefig(os.path.join(OUT, "mobility.png"), dpi=140,
                    bbox_inches="tight")
    except Exception as ex:
        print("plot skipped:", ex)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

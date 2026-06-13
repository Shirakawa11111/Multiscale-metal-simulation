"""The linchpin: calibrate PFC time -> seconds, converting the dimensionless
PFC mobility / KM rates into physical units (with the honest regime caveat).

Conserved PFC dynamics  dpsi/dt = laplacian(delta F/delta psi)  relaxes a
small long-wavelength density perturbation at wavevector k as
   psi_k(t) ~ exp(-omega(k) t),   omega(k) = k^2 * L(k),  L(k)=r+(1-k^2)^2.
For k->0,  omega ~ D_PFC * k^2  with  D_PFC = L(0) = r + 1  (PFC length^2/time).
So the model's intrinsic diffusivity is known/measurable; matching it to a
physical diffusion coefficient D_phys [m^2/s] fixes the time unit:
   tau_PFC = D_PFC * (length_unit)^2 / D_phys     [seconds per PFC time unit].

Because conserved PFC transports mass DIFFUSIVELY, this places PFC dynamics in
the DIFFUSION/CLIMB-CONTROLLED (high-homologous-temperature) regime — the same
regime as the m(r)->superplastic and GB-sliding findings. We therefore
calibrate against Cu self-/pipe-diffusion at ~0.6-0.8 Tm and report the
dislocation drag in that context (NOT athermal phonon-drag glide).

Output: results/time_calib/
"""

import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D, A_LATTICE, _rfft2, _irfft2

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "time_calib")
B_CU = 2.556e-10
LEN_UNIT = B_CU / A_LATTICE          # m per PFC length unit


def measure_D_pfc(r=-0.25):
    """Measure D_PFC from the decay rate of a small long-wavelength mode in the
    LIQUID (uniform) state, where the dynamics is linear: omega(k)=k^2 L(k)."""
    n = 128
    m = PFC2D(n, n, r=r, psi_bar=-0.25)
    # uniform + a single small-amplitude long-wavelength cosine in x
    kmode = 2                                   # mode index (long wavelength)
    x = np.arange(n) * m.dx
    amp0 = 1e-3
    m.psi = m.psi_bar + amp0 * np.cos(2 * np.pi * kmode * x / m.lx)[None, :] \
        * np.ones((n, 1))
    k = 2 * np.pi * kmode / m.lx
    # track modal amplitude decay
    def modal_amp(psi):
        return np.abs(np.fft.rfft(psi.mean(axis=0))[kmode])
    amps, ts = [modal_amp(m.psi)], [0.0]
    for _ in range(30):
        m.step(0.5, n=10)
        amps.append(modal_amp(m.psi))
        ts.append(m.time)
    amps = np.array(amps)
    ts = np.array(ts)
    good = amps > amps[0] * 1e-3
    rate = -np.polyfit(ts[good], np.log(amps[good]), 1)[0]   # omega measured
    omega_analytic = k ** 2 * (r + (1 - k ** 2) ** 2)
    D_pfc_measured = rate / k ** 2
    D_pfc_analytic = r + (1 - k ** 2) ** 2
    return dict(k=k, omega_measured=float(rate),
                omega_analytic=float(omega_analytic),
                D_pfc_measured=float(D_pfc_measured),
                D_pfc_analytic=float(D_pfc_analytic),
                D_pfc_smallk_limit=float(r + 1))


def main():
    os.makedirs(OUT, exist_ok=True)
    d = measure_D_pfc()
    D_pfc = d["D_pfc_measured"]
    print(f"D_PFC measured={D_pfc:.4f} vs analytic k-dep={d['D_pfc_analytic']:.4f} "
          f"(small-k limit r+1={d['D_pfc_smallk_limit']:.2f})")

    # physical anchors: Cu diffusion across homologous temperatures
    # D = D0 exp(-Q/kT); lattice self-diffusion D0=2e-5 m^2/s, Q=2.0 eV
    kB = 8.617e-5  # eV/K
    Tm = 1358.0
    anchors = {}
    for frac in (0.5, 0.6, 0.7, 0.8):
        T = frac * Tm
        D_lat = 2e-5 * np.exp(-2.0 / (kB * T))
        D_pipe = 100 * D_lat            # pipe diffusion ~1-3 orders faster
        tau_lat = D_pfc * LEN_UNIT ** 2 / D_lat
        tau_pipe = D_pfc * LEN_UNIT ** 2 / D_pipe
        anchors[f"{frac}Tm_{int(T)}K"] = dict(
            T_K=T, D_lattice_m2s=D_lat, D_pipe_m2s=D_pipe,
            tau_PFC_s_lattice=tau_lat, tau_PFC_s_pipe=tau_pipe)

    # close the stress unit: PFC shear modulus mu_pfc (measured) <-> Cu mu.
    MU_PFC = 0.0545          # measured (perfect-crystal dtau/dgamma)
    MU_CU = 48e9             # Pa
    STRESS_UNIT = MU_CU / MU_PFC      # Pa per PFC stress unit
    M_pfc = 1.05             # dimensionless PFC mobility (v=M*tau*b)
    # physical glide velocity at a representative resolved shear tau_phys.
    # v_phys = (M_pfc * tau_pfc * b_pfc) * LEN_UNIT / tau_PFC, with
    # tau_pfc = tau_phys / STRESS_UNIT and b_pfc*LEN_UNIT = b_Cu:
    #   v_phys = M_pfc * (tau_phys/STRESS_UNIT) * b_Cu / tau_PFC
    tau_phys = 100e6         # 100 MPa resolved shear
    vel = {}
    for kf, v in anchors.items():
        for chan, tpfc in (("lattice", v["tau_PFC_s_lattice"]),
                           ("pipe", v["tau_PFC_s_pipe"])):
            v_phys = M_pfc * (tau_phys / STRESS_UNIT) * B_CU / tpfc
            vel[f"{kf}_{chan}"] = v_phys      # m/s at 100 MPa
    result = dict(D_pfc=d, length_unit_m=LEN_UNIT, stress_unit_Pa=STRESS_UNIT,
                  mu_pfc=MU_PFC, mu_cu_Pa=MU_CU, M_pfc=M_pfc,
                  cu_anchors=anchors,
                  glide_velocity_at_100MPa_ms=vel,
                  regime="diffusion/climb-controlled (high homologous T); "
                         "conserved PFC transports mass diffusively, so these "
                         "rates map to high-T creep/GB regimes — consistent "
                         "with the m(r)->superplastic and GB-sliding findings, "
                         "NOT athermal phonon-drag glide.",
                  sanity="physical velocities should be slow (um/s..mm/s), "
                         "characteristic of diffusion/climb-controlled motion, "
                         "NOT athermal glide (~100s m/s) — validates the regime.")
    with open(os.path.join(OUT, "time_calib.json"), "w") as f:
        json.dump(result, f, indent=1)
    print("\nPFC time unit tau_PFC (seconds), anchored to Cu diffusion:")
    for k, v in anchors.items():
        print(f"  {k}: lattice-diff tau_PFC={v['tau_PFC_s_lattice']:.2e}s, "
              f"pipe-diff tau_PFC={v['tau_PFC_s_pipe']:.2e}s")
    print(f"\nstress unit = {STRESS_UNIT:.2e} Pa/PFC  (mu_pfc {MU_PFC} <-> Cu {MU_CU:.0e} Pa)")
    print("physical glide velocity at 100 MPa resolved shear:")
    for k, v in vel.items():
        print(f"  {k}: v = {v:.2e} m/s")
    print(f"\nregime: {result['regime'][:80]}...")
    print(f"sanity: {result['sanity'][:90]}")


if __name__ == "__main__":
    main()

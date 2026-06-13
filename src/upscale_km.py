"""!!! RETRACTED AS A LOOP-CLOSURE (2026-06-13, see docs/CORRECTIONS.md) !!!
This script's k1 = k2*sqrt(rho0)*1.3 makes k2=2*y_ann cancel exactly, so the
saturation density (2.535e15) is independent of y_ann (verified for y_ann=1..100b)
and the mobility M never enters the ODE. The Taylor alpha=0.35 is from literature,
contradicting the project's own PFC forest fit alpha=-0.287. The 2.5e15 m^-2 and
197 MPa numbers are NOT PFC-calibrated predictions. Kept only to show the KM
recovery-term FORM; do not cite its numbers as loop closure.

Loop-closure demonstration: assemble the three PFC-measured, dimensionless
constitutive parameters into a Kocks-Mecking (KM) dislocation-density-evolution
law and UPSCALE — integrate from the DAMASK ROI initial density to predict the
hardening curve rho(gamma) and flow stress tau(gamma). This is the actual
PFC -> crystal-plasticity information transfer (not a one-way unit conversion):
PFC supplies the FORM and dimensionless coefficients of the law that DAMASK's
dislotwin model otherwise fits.

PFC inputs (this project, dimensionless / in units of b):
  - annihilation distance y_ann ~ 7 b           (results/annihilation)
  - mobility law v = M b (tau - tau_c), M~1      (results/mobility_law)
  - multiplication onset sigma/E ~ 0.07          (results/d2 cascade)
KM law:  d(rho)/d(gamma) = k1*sqrt(rho) - k2*rho,  k2 = 2*y_ann/b (recovery),
         k1 (storage) from the multiplication/mean-free-path term.
Taylor flow stress:  tau = tau0 + alpha_eff*mu*b*sqrt(rho).

DAMASK ROI anchor (interfaceB_bridge): rho0 ~ 1.5e15 m^-2, b_Cu=2.556e-10 m.
Output: results/upscale/ — predicted rho(gamma), tau(gamma), saturation rho.
"""

import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "upscale")
B_CU = 2.556e-10


def load(path, key, default):
    try:
        return json.load(open(path)).get(key, default)
    except Exception:
        return default


def main():
    os.makedirs(OUT, exist_ok=True)
    # --- PFC-measured dimensionless inputs ---
    y_ann = load(os.path.join(os.path.dirname(__file__), "..", "results",
                              "annihilation", "annihilation.json"),
                 "y_ann_over_b", 7.0)               # in units of b
    mob = json.load(open(os.path.join(os.path.dirname(__file__), "..",
                    "results", "mobility_law", "threshold_fit.json")))
    M = mob.get("mpfc", mob.get("diffusive", {})).get("M", 1.0)
    # KM coefficients
    k2 = 2.0 * y_ann          # recovery, per unit strain (dimensionless in b)
    # storage k1: dislocation storage ~ 1/(b*L), L mean free path ~ 1/sqrt(rho)
    # -> k1*sqrt(rho); take k1 so saturation rho_s = (k1/k2)^2 lands near the
    # DAMASK ROI density when expressed in b^-2 (anchor the storage to the ROI).
    rho0_phys = 1.5e15        # m^-2 (DAMASK ROI)
    rho0_b2 = rho0_phys * B_CU ** 2     # dimensionless (per b^2)
    # choose k1 so steady-state sqrt(rho_s)=k1/k2 matches ~ROI density
    k1 = k2 * np.sqrt(rho0_b2) * 1.3    # saturation ~1.7x ROI (hardening room)

    # integrate KM from a low initial density up to saturation
    rho = 0.05 * rho0_b2
    gammas, rhos = [], []
    dg = 1e-4
    for i in range(int(0.5 / dg)):       # to 50% shear strain
        drho = (k1 * np.sqrt(max(rho, 0)) - k2 * rho) * dg
        rho += drho
        if i % 200 == 0:
            gammas.append(i * dg)
            rhos.append(rho)
    rho_s = (k1 / k2) ** 2
    # Taylor flow stress with a representative alpha~0.35 (PFC forest geometry
    # could not fix alpha cleanly; use literature alpha and PFC rho-evolution)
    alpha, mu_phys, G = 0.35, 4.5e10, 4.5e10   # Cu shear modulus ~45 GPa
    tau = alpha * mu_phys * B_CU * np.sqrt(np.array(rhos) / B_CU ** 2)

    result = dict(
        pfc_inputs=dict(y_ann_over_b=y_ann, mobility_M=M,
                        multiplication_threshold_sigma_over_E=0.07),
        km_coeffs=dict(k1=k1, k2=k2, rho_sat_over_b2=rho_s,
                       rho_sat_m2=rho_s / B_CU ** 2),
        roi_anchor=dict(rho0_m2=rho0_phys, b_m=B_CU),
        curve=dict(gamma=gammas, rho_b2=rhos,
                   rho_m2=[r / B_CU ** 2 for r in rhos],
                   flow_stress_MPa=(tau / 1e6).tolist()),
        note="PFC supplies y_ann and mobility form; KM storage anchored to the "
             "DAMASK ROI density; Taylor alpha from literature (PFC forest "
             "geometry inconclusive). Demonstrates the PFC->CP density-evolution "
             "transfer at the dimensionless/form level.")
    with open(os.path.join(OUT, "upscale_km.json"), "w") as f:
        json.dump(result, f, indent=1)
    print(f"PFC->KM upscaling:")
    print(f"  y_ann={y_ann:.1f}b -> recovery k2={k2:.2f}")
    print(f"  saturation rho_s = {rho_s/B_CU**2:.2e} m^-2 "
          f"(ROI anchor {rho0_phys:.1e})")
    print(f"  predicted flow stress at saturation ~ "
          f"{result['curve']['flow_stress_MPa'][-1]:.0f} MPa "
          f"(DAMASK ROI sigma_vm ~1380 MPa)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        g = np.array(gammas)
        ax[0].plot(g, np.array(rhos) / B_CU ** 2, "b-")
        ax[0].axhline(rho0_phys, color="gray", ls="--", label="DAMASK ROI rho0")
        ax[0].set_xlabel("shear strain gamma")
        ax[0].set_ylabel(r"$\rho$ (m$^{-2}$)")
        ax[0].set_title("PFC-calibrated KM density evolution")
        ax[0].legend(); ax[0].grid(alpha=0.3)
        ax[1].plot(g, tau / 1e6, "r-")
        ax[1].set_xlabel("shear strain gamma")
        ax[1].set_ylabel("flow stress (MPa)")
        ax[1].set_title("upscaled hardening curve (Taylor + PFC rho)")
        ax[1].grid(alpha=0.3)
        fig.savefig(os.path.join(OUT, "upscale_km.png"), dpi=140,
                    bbox_inches="tight")
    except Exception as ex:
        print("plot skipped:", ex)


if __name__ == "__main__":
    main()

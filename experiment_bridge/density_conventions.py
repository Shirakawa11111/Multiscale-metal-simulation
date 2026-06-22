"""Density reporting conventions for a STEM-foil network [LOCAL, no DDD].

The cell-policy audit showed apparent volume density is a convention (rho ~ 1/zbox). This makes the
distinction explicit and foil-aware:
  - projected-area line density  Lambda_A = total_line_length / projected_area   [m^-1]  <- foil-NATIVE
  - volume density               rho      = total_line_length / volume = Lambda_A / z   [m^-2]  <- convention
A STEM foil directly constrains Lambda_A (length per projected area); the bulk-equivalent rho requires a
declared effective thickness / zbox. Computed from the as-built lowered network (deterministic).

  python3 experiment_bridge/density_conventions.py  -> results_exadis/density_conventions.{json,md}
"""
import os, sys, json
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
from defect_ir.adapters.to_exadis import idr_to_exadis_network

OUT = os.path.join(os.path.dirname(__file__), "results_exadis")
IDR = os.path.join(OUT, "cu_stem_idr.json")


def main():
    doc = json.load(open(IDR))
    b = doc["units"]["length_unit_m"]
    net = idr_to_exadis_network(doc, assignment_policy="top1", cell_policy="as_is")
    h = np.array(net["cell"]["h_b"], float)
    nodes = np.array([[n[0], n[1], n[2]] for n in net["nodes"]], float)
    # total built line length (b units -> m)
    Lb = sum(np.linalg.norm(nodes[int(s[0])] - nodes[int(s[1])]) for s in net["segs"])
    L_m = Lb * b
    LXY_b = h[0, 0]
    proj_area_m2 = (LXY_b * b) ** 2                       # foil projected (xy) area
    lambda_area_m1 = L_m / proj_area_m2                   # foil-native line density [m^-1]
    z_conv = {"foil_z600": 600.0, "thickened_z3_1800": 1800.0,
              "thickened_z5_3000": 3000.0, "thickened_z10_6000": 6000.0}
    vol = {k: lambda_area_m1 / (z * b) for k, z in z_conv.items()}    # rho = Lambda_A / z_eff (AS-BUILT)

    # ---- relaxed layer: fold in the measured DDD line-length relaxation (cell-policy audit M5) ----
    # The as-built geometry is the foil-native reconstruction state; the DDD-legalized/relaxed state is the
    # actual simulation initial condition. relaxed rho per cell config comes from cell_policy_audit_summary;
    # relaxed Lambda_A = rho_relaxed * z * b, and the relaxation fraction = relaxed/as-built.
    relaxed = None
    audit_path = os.path.join(OUT, "cell_policy_audit_summary.json")
    if os.path.exists(audit_path):
        rho_relaxed = json.load(open(audit_path)).get("rho_app_m2", {})
        cfg_map = {"foil_z600": ("foil", 600.0), "thickened_z3_1800": ("z3LT", 1800.0),
                   "thickened_z5_3000": ("z5LT", 3000.0), "thickened_z10_6000": ("z10LT", 6000.0)}
        relaxed = {}
        for k, (akey, z) in cfg_map.items():
            if akey in rho_relaxed:
                rr = rho_relaxed[akey]
                la_r = rr * z * b
                relaxed[k] = {
                    "volume_density_m2": float(rr),
                    "lambda_area_m1": float(la_r),
                    "line_length_relaxation_fraction": round(la_r / lambda_area_m1, 4),
                }
    rep = {
        "source": os.path.basename(IDR),
        "as_built": {
            "total_line_length_m": L_m,
            "projected_area_m2": proj_area_m2,
            "lambda_area_m1_foil_native": lambda_area_m1,
            "volume_density_m2_by_convention": {k: float(v) for k, v in vol.items()},
        },
        "relaxed_from_cell_policy_audit": relaxed,
        "note": ("Lambda_A (projected line density) is the foil-native observable (independent of z). "
                 "Volume density rho = Lambda_A / z_eff is convention-dependent (rho ~ 1/zbox). Report TWO "
                 "states: AS-BUILT (the STEM reconstruction geometry, foil-native) and RELAXED (the "
                 "DDD-legalized simulation initial condition; line length relaxes to ~0.59-0.69). Always "
                 "quote rho with its z_eff/zbox convention AND its state (as-built vs relaxed)."),
    }
    json.dump(rep, open(os.path.join(OUT, "density_conventions.json"), "w"), indent=1)
    rows = "\n".join(f"| {k} | {v:.3e} |" for k, v in vol.items())
    if relaxed:
        rlx_rows = "\n".join(
            f"| {k} | {v['volume_density_m2']:.3e} | {v['lambda_area_m1']:.3e} | {v['line_length_relaxation_fraction']} |"
            for k, v in relaxed.items())
        fr = [v["line_length_relaxation_fraction"] for v in relaxed.values()]
        relaxed_md = f"""
## Relaxed state (DDD-legalized initial condition)
After zero-stress DDD relaxation the line length contracts to {min(fr)}–{max(fr)} of as-built (cell-dependent:
foil relaxes most). The relaxed network is the actual *simulation* initial condition; the as-built is the
*microscope reconstruction* state. Report both — they answer different questions.

| convention | relaxed ρ (m⁻²) | relaxed Λ_A (m⁻¹) | relaxation fraction |
|--|--|--|--|
{rlx_rows}
"""
    else:
        relaxed_md = ""
    md = f"""# Density reporting conventions (foil-aware, as-built vs relaxed)

`{rep['source']}`, top1 lowering.

## As-built state (STEM reconstruction geometry, foil-native)
| quantity | value |
|--|--|
| total line length | {L_m:.3e} m |
| projected (xy) area | {proj_area_m2:.3e} m² |
| **Λ_A (projected line density, foil-native)** | **{lambda_area_m1:.3e} m⁻¹** |

Bulk-equivalent volume density ρ = Λ_A / z_eff (convention-dependent, ∝ 1/z):

| convention | ρ (m⁻²) |
|--|--|
{rows}
{relaxed_md}
**Reading.** A STEM foil directly constrains **Λ_A** (line length per projected area) — report this. The
bulk-equivalent **ρ** is a derived quantity that requires a declared effective thickness / zbox and scales
as 1/z (consistent with `CELL_POLICY_AUDIT.md`). And report the **state**: *as-built* (microscope geometry)
vs *relaxed* (DDD-legalized simulation initial condition, line length ×0.59–0.69). So: projected line density
is foil-native; volume density is convention-dependent; both come in an as-built and a relaxed flavor.
"""
    open(os.path.join(OUT, "density_conventions.md"), "w").write(md)
    rtxt = ("; relaxed foil Lambda_A=%.3e (frac %.3f)"
            % (relaxed["foil_z600"]["lambda_area_m1"], relaxed["foil_z600"]["line_length_relaxation_fraction"])
            ) if relaxed else ""
    print(f"Lambda_A={lambda_area_m1:.3e} m^-1 (foil-native, as-built); "
          f"rho: foil={vol['foil_z600']:.2e}, z5={vol['thickened_z5_3000']:.2e} m^-2{rtxt}")
    print("-> results_exadis/density_conventions.{json,md}")


if __name__ == "__main__":
    main()

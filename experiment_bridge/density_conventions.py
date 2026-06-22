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
    vol = {k: lambda_area_m1 / (z * b) for k, z in z_conv.items()}    # rho = Lambda_A / z_eff
    rep = {
        "source": os.path.basename(IDR),
        "total_line_length_m": L_m,
        "projected_area_m2": proj_area_m2,
        "lambda_area_m1_foil_native": lambda_area_m1,
        "volume_density_m2_by_convention": {k: float(v) for k, v in vol.items()},
        "note": ("Lambda_A (projected line density) is the foil-native observable (independent of z). "
                 "Volume density rho = Lambda_A / z_eff is convention-dependent (rho ~ 1/zbox), matching "
                 "the cell-policy audit. Always quote rho with its z_eff/zbox convention."),
    }
    json.dump(rep, open(os.path.join(OUT, "density_conventions.json"), "w"), indent=1)
    rows = "\n".join(f"| {k} | {v:.3e} |" for k, v in vol.items())
    md = f"""# Density reporting conventions (foil-aware)

`{rep['source']}` (as-built, top1 lowering).

| quantity | value |
|--|--|
| total line length | {L_m:.3e} m |
| projected (xy) area | {proj_area_m2:.3e} m² |
| **Λ_A (projected line density, foil-native)** | **{lambda_area_m1:.3e} m⁻¹** |

Bulk-equivalent volume density ρ = Λ_A / z_eff (convention-dependent):

| convention | ρ (m⁻²) |
|--|--|
{rows}

**Reading.** A STEM foil directly constrains **Λ_A** (line length per projected area) — report this. The
bulk-equivalent **ρ** is a derived quantity that requires a declared effective thickness / zbox and scales
as 1/z (consistent with `CELL_POLICY_AUDIT.md`). So: *projected line density is foil-native; volume density
is a convention-dependent derived quantity.*
"""
    open(os.path.join(OUT, "density_conventions.md"), "w").write(md)
    print(f"Lambda_A={lambda_area_m1:.3e} m^-1 (foil-native); "
          f"rho: foil={vol['foil_z600']:.2e}, z5={vol['thickened_z5_3000']:.2e} m^-2")
    print("-> results_exadis/density_conventions.{json,md}")


if __name__ == "__main__":
    main()

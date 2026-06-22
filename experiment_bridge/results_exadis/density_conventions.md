# Density reporting conventions (foil-aware)

`cu_stem_idr.json` (as-built, top1 lowering).

| quantity | value |
|--|--|
| total line length | 1.364e-05 m |
| projected (xy) area | 4.181e-12 m² |
| **Λ_A (projected line density, foil-native)** | **3.262e+06 m⁻¹** |

Bulk-equivalent volume density ρ = Λ_A / z_eff (convention-dependent):

| convention | ρ (m⁻²) |
|--|--|
| foil_z600 | 2.127e+13 |
| thickened_z3_1800 | 7.091e+12 |
| thickened_z5_3000 | 4.254e+12 |
| thickened_z10_6000 | 2.127e+12 |

**Reading.** A STEM foil directly constrains **Λ_A** (line length per projected area) — report this. The
bulk-equivalent **ρ** is a derived quantity that requires a declared effective thickness / zbox and scales
as 1/z (consistent with `CELL_POLICY_AUDIT.md`). So: *projected line density is foil-native; volume density
is a convention-dependent derived quantity.*

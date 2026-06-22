# Density reporting conventions (foil-aware, as-built vs relaxed)

`cu_stem_idr.json`, top1 lowering.

## As-built state (STEM reconstruction geometry, foil-native)
| quantity | value |
|--|--|
| total line length | 1.364e-05 m |
| projected (xy) area | 4.181e-12 m² |
| **Λ_A (projected line density, foil-native)** | **3.262e+06 m⁻¹** |

Bulk-equivalent volume density ρ = Λ_A / z_eff (convention-dependent, ∝ 1/z):

| convention | ρ (m⁻²) |
|--|--|
| foil_z600 | 2.127e+13 |
| thickened_z3_1800 | 7.091e+12 |
| thickened_z5_3000 | 4.254e+12 |
| thickened_z10_6000 | 2.127e+12 |

## Relaxed state (DDD-legalized initial condition)
After zero-stress DDD relaxation the line length contracts to 0.59–0.6856 of as-built (cell-dependent:
foil relaxes most). The relaxed network is the actual *simulation* initial condition; the as-built is the
*microscope reconstruction* state. Report both — they answer different questions.

| convention | relaxed ρ (m⁻²) | relaxed Λ_A (m⁻¹) | relaxation fraction |
|--|--|--|--|
| foil_z600 | 1.255e+13 | 1.925e+06 | 0.59 |
| thickened_z3_1800 | 4.363e+12 | 2.007e+06 | 0.6153 |
| thickened_z5_3000 | 2.718e+12 | 2.084e+06 | 0.6389 |
| thickened_z10_6000 | 1.458e+12 | 2.237e+06 | 0.6856 |

**Reading.** A STEM foil directly constrains **Λ_A** (line length per projected area) — report this. The
bulk-equivalent **ρ** is a derived quantity that requires a declared effective thickness / zbox and scales
as 1/z (consistent with `CELL_POLICY_AUDIT.md`). And report the **state**: *as-built* (microscope geometry)
vs *relaxed* (DDD-legalized simulation initial condition, line length ×0.59–0.69). So: projected line density
is foil-native; volume density is convention-dependent; both come in an as-built and a relaxed flavor.

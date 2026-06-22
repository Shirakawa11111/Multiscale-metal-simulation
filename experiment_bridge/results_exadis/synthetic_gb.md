# Synthetic g·b ambiguity collapse (g·b-ready interface)

`cu_stem_idr.json`, 243 edges, invisibility tol=0.1. Synthetic ground truth: ONE random true
Burgers per **parent line** (line-coherent, consistent with IDR v1.1), shared by all that line's segments.

| scenario | # reflections | mean entropy (bits) | mean candidates | frac resolved |
|--|--|--|--|--|
| no_gb | 0 | 1.5849 | 3.0 | 0.0 |
| 1_reflection_g200 | 1 | 0.7037 | 1.7 | 0.296 |
| 2_reflections_g200_g020 | 2 | 0.0 | 1.0 | 1.0 |
| 3_reflections_g200_g020_g002 | 3 | 0.0 | 1.0 | 1.0 |

**Reading.** With no g·b the assignment is ~3-way ambiguous (~1.58 bits). Each added reflection that obeys
the invisibility criterion (|g·b|≈0) removes incompatible candidates; ~2 well-chosen reflections collapse
most lines to a single slip system (`gb_validated`). This is the interface to plug in **real** diffraction
contrast: populate each line's `uncertainty.gb_constraints` with observed (g, visible) pairs and call
`apply_gb_constraints` to upgrade the IDR from `geometry_only_pending_gb` to `gb_validated`.

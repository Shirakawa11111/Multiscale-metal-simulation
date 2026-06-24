# STEM-to-DDD v2 — release notes

**Version:** `STEM-to-DDD v2 / IDR v1.1`
**Manifest:** `results_exadis/audit_package_manifest.json` · **Index:** `AUDIT_PACKAGE_INDEX.md`

## Scope
An **auditable microscope-to-DDD initial-condition package** — a real STEM dislocation reconstruction carried
into ExaDiS/DDD with every modeling assumption (slip-system assignment, cell policy, density convention,
endpoint policy) exposed, quantified, and self-checked. **It is not a hardening prediction.**

## Core claims (each backed by a result file)
1. **Geometry-only assignment exposes a 3-way Burgers ambiguity.** {111} plane fixed; the 3 ⟨110⟩ Burgers are
   near-degenerate → mean confidence 0.333, entropy 1.585 bits, 243/243 edges ambiguous. (`cu_stem_idr_report`)
2. **Edgewise sampling is self-falsified as unphysical.** Per-edge draws inject 142.8/216 within-line Burgers
   discontinuities (artificial junctions); the old "assignment → topology 5.4×" is **retracted**.
   (`assignment_sensitivity`, `bo_uq_pilot_summary_v0_superseded`)
3. **`sample_linewise` is the physical default.** One draw per parent line → 0 within-line discontinuities;
   assignment becomes a *minor* topology knob (junctions ≈ top-1). (`v11_linewise_summary`)
4. **Cell policy controls the density reporting convention, not physics.** ρ_app ∝ 1/zbox; force model ~8%;
   survival / line-length relaxation / junction count are cell-robust. Foil-native observable is **Λ_A**
   (3.26e6 m⁻¹); ρ_vol needs a declared z_eff and state (as-built vs relaxed). (`CELL_POLICY_AUDIT`,
   `density_conventions`)
5. **The IDR is g·b-ready.** Synthetic line-coherent reflections collapse entropy 1.585 → 0.704 → 0 over
   0/1/2 reflections; real (g, visible) data plugs into `gb_constraints`. (`synthetic_gb`, `GB_DATA_REQUIREMENTS`)

## Reproduce
```
python3 experiment_bridge/run_local_audit_package.py     # -> LOCAL AUDIT PACKAGE: PASS
```
Evidence committed: `results_exadis/local_audit_pass.txt` (run log) and
`results_exadis/py_compile_pass.txt` (fresh-clone `py_compile` — all package `.py` are valid LF multi-line
Python; a GitHub-raw viewer may collapse line breaks, which is a rendering artifact, not the file content).

## Key invariants (regression guards)
- within-line discontinuities: linewise **0**, edgewise artifact **142.8**/216
- assignment entropy: no-g·b **1.585 bits** → 2 reflections **0**
- Λ_A foil-native (as-built) **3.26e6 m⁻¹**; ρ_vol ∝ 1/z_eff

## Limitations (stated honestly)
- **No real g·b yet** — assignment resolution is shown synthetically; real diffraction contrast is the next input.
- **Small network** — 27 lines / 270 vertices / 243 edges.
- **Short DDD loading** — sensitivity ratios and stability, **not** converged hardening numbers.

## Sealed / deferred
- **Sealed:** collinear / DDD interaction-kernel validation (drift-limited bounded negative,
  `interaction_matrix/multislip_flow/CONCLUSION.md`).
- **Deferred:** full Bayesian optimization (UQ objectives frozen first, `BO_UQ_PILOT.md`).

## Next stage (not started)
Real / partial g·b integration: compare **geometry-only → partial-g·b → gb-validated** IDRs on the same
network across assignment entropy, slip-system inventory, DDD stability, junction count, and density
convention. Experimental input spec: `GB_DATA_REQUIREMENTS.md`.

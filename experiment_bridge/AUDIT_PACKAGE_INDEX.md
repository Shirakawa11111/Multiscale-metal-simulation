# STEM-to-DDD v2 — audit package index

Navigation for the `experiment_bridge/` audit package (IDR v1.1). The package is an *auditable* path from a
real STEM reconstruction to a DDD initial condition — not a hardening curve. One-line reproduce of every
deterministic step: `python3 experiment_bridge/run_local_audit_package.py`.

## Core report
- **[`STEM_TO_DDD_V2_AUDIT.md`](STEM_TO_DDD_V2_AUDIT.md)** — the integrated argument (input → IDR uncertainty
  → lowering correction → DDD stability → cell policy → g·b path) + paper-ready conclusion.
- **`results_exadis/audit_summary_figure.png`** — 4-panel main result (`make_audit_figure.py`).

## Input → IDR
- `stem_to_idr.py` → `results_exadis/cu_stem_idr.json` + `cu_stem_idr_report.{md,json}`
  (27 lines / 270 vertices / 243 edges; mean confidence 0.333; entropy 1.585 bits; 243/243 ambiguous).

## Assignment (lowering policy)
- `assignment_sensitivity.py` → `results_exadis/assignment_sensitivity.{md,json}`
- **`sample_linewise`** = physical default (0 within-line discontinuities).
- `sample_edgewise` = **deprecated** artifact stress-test (142.8/216 within-line Burgers discontinuities).

## Density & cell policy
- **[`CELL_POLICY_AUDIT.md`](CELL_POLICY_AUDIT.md)** — ρ_app ∝ 1/zbox (convention); survival/topology robust.
- `density_conventions.py` → `results_exadis/density_conventions.{md,json}` — foil-native **Λ_A** (z-invariant,
  as-built 3.26e6 m⁻¹ / relaxed 1.92e6 m⁻¹) vs convention-dependent **ρ_vol** (as-built + relaxed states).
- JSON schema: `{ source, as_built{...}, relaxed_from_cell_policy_audit{...}, note }`.

## g·b-ready path
- `synthetic_gb.py` → `results_exadis/synthetic_gb.{md,json}` — **line-coherent** (one true Burgers per
  parent line); entropy collapse 1.585 → 0.704 → 0 bits over 0/1/2 reflections.
- Interface: `defect_ir.uncertainty.apply_gb_constraints` (real diffraction data plugs in here).

## Reproducibility
- **`run_local_audit_package.py`** — one-click: runs all local steps + figure, checks outputs, asserts
  invariants (edgewise 142.8 / linewise 0; g·b 1.585→0).
- `tests/test_defect_ir.py` — gate test (examples valid, lowering round-trips, linewise coherent, g·b collapse).
- [`AUDIT_MANIFEST.md`](AUDIT_MANIFEST.md) — authoritative record of the DDD/HPC runs (not rerun locally).

## Superseded (do not read as current)
- `results_exadis/bo_uq_pilot_summary_v0_superseded.json` — v0 edgewise ranking (5.4× topology). **RETRACTED**;
  canonical is `results_exadis/v11_linewise_summary.json` + [`BO_UQ_PILOT.md`](BO_UQ_PILOT.md).
- `stem_to_exadis.py` — LEGACY direct single-assignment converter (kept as baseline).

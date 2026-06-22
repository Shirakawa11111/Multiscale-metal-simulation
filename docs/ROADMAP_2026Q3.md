# Roadmap — 2026 Q3

> Supersedes the PFC-centric `RESEARCH_PLAN.md` (kept as historical record, frozen).
> Honest top-level status: [`../PROJECT_STATUS.md`](../PROJECT_STATUS.md).

The project has grown from a PFC starting point into a **defect-graph → multiscale-simulation framework**.
The core contribution is not a single solver but a **unified defect intermediate representation (IDR)** that
carries 2D (graphene/h-BN) and 3D (Cu dislocation / STEM) defects into simulation *with quantified
uncertainty* — turning each route from "it runs" into "it is auditable."

## 1. Sealed branches (done — not mainline)
- **PFC regime boundary**: PFC (conservative diffusive variant) does not forest-harden under room-temperature
  monotonic tension — a regime result, not a bug. Reframed as a 2D method prototype + high-T/creep branch.
- **DDD interaction-kernel validation** (`../interaction_matrix/`): controlled binary reactions confirm the
  collinear *mechanism* + remobilization length scaling (R²=0.89), but local strength, pairwise MFP, and the
  co-driven multi-slip flow-stress density-lever **do not reproduce** canonical collinear dominance.
  **Drift-limited bounded negative**; sealed as a protocol-dependence benchmark
  (`../interaction_matrix/multislip_flow/CONCLUSION.md`).

## 2. Current mainline
- **IDR** — unified defect intermediate representation (the structural backbone / "central nervous system").
- **STEM-to-DDD v2** — uncertainty-aware pipeline: `stem_to_idr` → `idr_to_exadis`, with slip-system
  *candidates + confidence* (not a single forced assignment), explicit cell-policy and endpoint-policy.
- **BO/UQ calibration layer** — start with sensitivity/UQ pilot (not full Bayesian optimization).

## 3. Near-term milestones
- **M1 — IDR schema + validators** (`defect_ir/`) — **DONE ✓**: 6-section `defect_idr_v1`
  (`schema.py`), structural + physics validators (`validators.py`: unique ids, b·n=0, box consistency,
  candidate priors), `uncertainty.py` (candidate priors + entropy), `adapters/to_exadis.py`,
  `IDR_SPEC.md`, 2 validated examples, gate test.
- **M2 — STEM-to-DDD uncertainty envelope** — **IN PROGRESS**: `stem_to_idr.py` (real recon → IDR,
  top-k candidates + audit report `cu_stem_idr_report.{json,md}`) + `idr_to_exadis.py` (CLI, selectable
  assignment/cell policy); `stem_to_exadis.py` marked LEGACY; `CELL_POLICY.md`, `ASSIGNMENT_UNCERTAINTY.md`.
  *Remaining:* formal sensitivity sweep harness (local), then M3.
- **M3 — real-network DDD audit report** — **PILOT DONE, v1.1-corrected** (experiment_bridge/REAL_NETWORK_AUDIT.md):
  import is stable & auditable; the v0 "assignment→topology 5.4×" was a per-edge sampling artifact
  (within-line Burgers discontinuities), corrected with line-coherent `sample_linewise`. Robust result:
  **cell policy dominates apparent density (~5.2×, deconfounded from force)**; assignment ambiguity is a
  minor topology knob once propagated per line.
- **M4 — BO/UQ pilot** — **PILOT DONE, v1.1-corrected** (experiment_bridge/BO_UQ_PILOT.md): corrected knob
  ranking — cell policy → density (dominant); force → density (1.17×); line-coherent assignment → minor
  topology; endpoint minor; survival robust. (v0 edgewise ranking retracted to an appendix.)
  **Default lowering policy: `sample_linewise`.**

Default pipeline: `stem_to_idr.py → idr_to_exadis.py`. Legacy: `stem_to_exadis.py`.

## 4. No longer investing
- chasing canonical collinear coefficient / random probe-matrix recovery
- PFC room-temperature forest hardening
- any new DDD interaction-kernel protocol beyond reproduction + documentation

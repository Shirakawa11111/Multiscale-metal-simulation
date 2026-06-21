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
- **M1 — IDR schema + validators** (`defect_ir/`): geometry / topology / physics_labels / uncertainty /
  provenance / simulation_targets; 2 worked examples (Cu STEM dislocation graph, graphene defect graph).
- **M2 — STEM-to-DDD uncertainty envelope**: split `stem_to_exadis.py` into `stem_to_idr.py` +
  `idr_to_exadis.py`; per-line slip-system candidate set + assignment confidence; `CELL_POLICY.md`,
  `ASSIGNMENT_UNCERTAINTY.md`.
- **M3 — real-network DDD audit report**: zero-stress relaxation stability, loading response, density
  evolution, topology events, slip-system inventory, assignment/z-scaling/endpoint sensitivity → audit JSON + md.
- **M4 — BO/UQ pilot**: sensitivity over image/recon, crystallography-assignment, DDD-legalization, loading
  knobs; objectives = stability + interpretability (network survival, density-growth plausibility,
  topology-event rate, assignment sensitivity, agreement with observed density increment).

## 4. No longer investing
- chasing canonical collinear coefficient / random probe-matrix recovery
- PFC room-temperature forest hardening
- any new DDD interaction-kernel protocol beyond reproduction + documentation

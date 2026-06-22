# Defect-centric Multiscale Simulation Framework

**From 2D graphene/h-BN defect graphs to 3D STEM-to-DDD dislocation networks** — a unified,
uncertainty-aware pipeline that carries experimental/simulated defects through multiscale simulation
and **quantifies** the assumptions at each step (not just "it runs", but "it is auditable").

Developed and tested by Claude Code in a self-driven iteration loop. Honest top-level status:
[PROJECT_STATUS.md](PROJECT_STATUS.md) · roadmap: [docs/ROADMAP_2026Q3.md](docs/ROADMAP_2026Q3.md).

## Current mainline
- **`defect_ir/`** — the unified Defect Intermediate Representation (IDR): one schema for 2D atomic
  defect graphs and 3D dislocation networks, with per-assignment **candidates + confidence** and
  explicit uncertainty. Spec: [defect_ir/IDR_SPEC.md](defect_ir/IDR_SPEC.md). **[M1 ✓]**
- **`experiment_bridge/`** — STEM-to-DDD v2: `stem_to_idr.py` (reconstruction → IDR, top-k slip-system
  candidates) → `idr_to_exadis.py` (IDR → ExaDiS, selectable assignment/cell policy). Consolidated audit:
  `STEM_TO_DDD_V2_AUDIT.md` (one-click reproduce via `run_local_audit_package.py`, main figure
  `make_audit_figure.py`). Legacy direct converter `stem_to_exadis.py` kept as baseline. Transparency:
  `CELL_POLICY.md`, `ASSIGNMENT_UNCERTAINTY.md`. **[M2 ✓]**
- **BO/UQ** — sensitivity/UQ over assignment / cell / endpoint policy (next).

## Structure
```
defect_ir/          unified IDR schema + validators + uncertainty + adapters (+ examples, spec)
experiment_bridge/  STEM->IDR / IDR->ExaDiS / real-network DDD; cell & assignment policy docs
interaction_matrix/ SEALED: DDD interaction-kernel validation (collinear bounded negative)
taylor_hardening/   DDD Taylor/Kocks-Mecking regime-correct hardening; 5-reviewer peer review
md_rung/            MD elastic-constant anchoring of the Taylor prefactor (mu b)
src/                PFC branch: 2D/3D phase-field-crystal mechanism + regime-boundary prototype
validation_4d/      Interface-B / PFC validation boundary (FINDINGS.md)
docs/               ROADMAP_2026Q3.md (current) · RESEARCH_PLAN.md (frozen, PFC-era)
```

## Quick start
```bash
python3 -m defect_ir.examples.build_examples     # build + validate the 2D & 3D IDR examples
python3 experiment_bridge/stem_to_idr.py         # real STEM recon -> IDR + audit report
python3 tests/test_defect_ir.py                  # IDR gate test
```

## Sealed branches (results, no longer mainline)
- **DDD interaction-kernel validation** — collinear mechanism + remobilization length scaling are real
  (R²=0.89), but canonical collinear *dominance* does not reproduce in the multi-slip flow-stress
  density-lever (drift-limited bounded negative). [interaction_matrix/multislip_flow/CONCLUSION.md](interaction_matrix/multislip_flow/CONCLUSION.md)
- **PFC regime boundary** — PFC (conservative diffusive variant) does not forest-harden under
  room-temperature monotonic tension (a regime result). Retained as a 2D method prototype + high-T branch.

---

## PFC branch — gate status & mechanism results (historical, all gates pass)
<details><summary>expand</summary>

| gate | content | key number |
|--|--|--|
| A1a/b | 2D crystallization (melt/seed) | lattice-const dev 0.3%/0.5%, perfect-crystal frac6=1.000 |
| A2 | edge-dipole insertion + relax | exactly 2 stable cores, defect frac 0.44% |
| A3a/b | elastic tension / slip | modulus 0.211 (energy) vs 0.209 (stress); slip ~4a₀ @2% |
| C1/C2 | 3D BCC crystallize / tension | 2 peaks/cell, NN dev 4.1%; linear stress, peaks conserved |
| D1/D2 | Interface-B density match + cascade | ROI ρ≈2.2e15 → yield σ/E≈7% @11.6% → avalanche (316 cores) → flow |

`tests/`: `test_a1_crystallization.py`, `test_a2_dislocation_dipole.py`, `test_a3_tension.py`,
`test_c1_bcc_3d.py`. Deps: numpy, scipy, matplotlib, pyfftw (optional). Figures under `results/`.
</details>

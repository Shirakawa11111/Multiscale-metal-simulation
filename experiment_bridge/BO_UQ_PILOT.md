# BO/UQ sensitivity pilot (M4, v1.1)

A *sensitivity* pilot (not full Bayesian optimization, not stress-curve matching): sweep the IDR's exposed
knobs and rank which one drives which **objective** (stability + interpretability), using the physical
**line-coherent** assignment policy. Data: `results_exadis/v11_linewise_summary.json`,
`results_exadis/assignment_sensitivity.json`; figure `results_exadis/v11_linewise.png`. Builds on
`REAL_NETWORK_AUDIT.md` (v1.1).

## Corrected knob → objective ranking
| rank | knob | objective | magnitude |
|--|--|--|--|
| 1 | **cell policy** (foil ↔ thickened periodic) | **apparent density** | **~5.2×** (deconfounded: same-force foil→thickened) |
| 2 | **force model** at fixed thickened cell | density | ~1.17× (DDD_FFT vs LineTension) |
| 3 | **assignment** (`sample_linewise`) | topology | **minor / residual** (~top-1 ± seed scatter) |
| 4 | endpoint policy (pinned/free) | topology | minor (~7%) |
| — | network survival | **robust** to all knobs (import does not collapse) |

## Reading it
- **Cell policy is the dominant uncertainty** — and it is the cell *volume / periodicity*, not the force
  model (force only 1.17× at a fixed thickened cell). Any density number must be reported with its cell policy.
- **Assignment ambiguity is a minor topology knob** once propagated *per line* (`sample_linewise`). The
  earlier "assignment dominates topology ~5×" was an edgewise artifact (see Appendix + `REAL_NETWORK_AUDIT.md`).
- **Survival is robust** — the STEM→DDD import is a stable initial condition under every knob.

## Decision value
- The highest-leverage modeling choice is the **cell policy** (foil vs thickened-periodic) for density
  normalization → the natural next experiment is a **cell-policy density audit** (sweep zbox + force at fixed
  line-coherent assignment).
- Resolving the slip-system assignment via experimental **g·b** is still worthwhile — it upgrades a
  geometry-only IDR to a physics-validated one (collapses the assignment entropy) — but, per v1.1, it is
  **not** expected to produce a large topology swing; its value is correctness/auditability, not a 5× effect.

## UQ objective freeze (pre-BO)
Before any *full* Bayesian optimization the objective is frozen into three classes, so BO optimizes a defined
target rather than a vague stress curve. The pilot above tells us which class each pipeline knob actually
moves (see `STEM_TO_DDD_V2_AUDIT.md` §"three classes"):

| class | observables | what the audit found |
|--|--|--|
| **A — stability** | network survival, segment-count growth, line-length relaxation (0.59–0.69) | **robust** to all knobs |
| **B — topology** | junction count, topology-event proxy, **within-line discontinuity ≡ 0** | minor; line-coherent-sensitive |
| **C — reporting** | **Λ_A projected line density (foil-native)**, ρ_vol *under declared cell convention* | a **convention**, not physics (ρ ∝ 1/zbox) |

Within-line discontinuity is **not** a tunable objective — it is a hard legality invariant (must be 0;
`sample_edgewise` violates it). The first BO/UQ question is therefore *which knobs move a reportable
observable vs which only change a reporting convention*, not "minimize flow stress." Full BO stays **deferred**
until (a) a real g·b assignment collapses class B's residual ambiguity and (b) converged loading makes the
class-A density-growth objective discriminating.

*Caveats:* small system (270 nodes), short loading (density-growth objective not yet discriminating); these
are sensitivity ratios, not converged hardening numbers.

---

## Appendix — SUPERSEDED v0 (edgewise) ranking
The v0 pilot (35 runs, `m4/`, `bo_uq_pilot.png`, `bo_uq_pilot_summary.json`) used the edgewise `sample`
policy and reported: assignment → topology **5.4×** (top1 ~5 vs sample ~27); cell/zbox → density 4.4×;
endpoint minor; survival robust. The 5.4× was an edgewise within-line-discontinuity artifact (143/216
discontinuities); it is retracted. Density/survival conclusions from v0 stand (and are refined above).

# BO/UQ sensitivity pilot (M4)

Goal (per the roadmap): a *sensitivity* pilot — not full Bayesian optimization, and not stress-curve
matching. Sweep the IDR's exposed knobs and rank which one drives which **objective** (stability +
interpretability). 35 ExaDiS runs on the real 270-node Cu network: 5 configs (cell/zbox/endpoint) ×
{top-1 + 6 assignment samples}, conservative ≤30 cores. Figure: `results_exadis/bo_uq_pilot.png`,
data: `results_exadis/bo_uq_pilot_summary.json` (builds on the M3 audit, `REAL_NETWORK_AUDIT.md`).

## Knobs × objectives — the ranking
| knob | objective it drives | magnitude |
|--|--|--|
| **slip-system assignment** (top-1 vs sampling the 3-way ambiguity) | **topology** (junction formation) | **5.4×** (top-1 ~5 vs sample ~27 junctions) + within-sample CV ~0.15 |
| **cell policy / zbox** (foil vs thickened periodic) | **density** | **4.4×** (foil 1.4e13 vs thickened 3.2e12); topology ~flat (z3/z5/z10 → 29/30/30) |
| **endpoint policy** (pinned vs free) | topology | minor, ~1.07× (30 vs 28) |
| — | **network survival** | **robust** (~0.41 across *all* knobs — import does not collapse) |

## Reading it
- **Assignment policy is the dominant uncertainty for topology**, by far: the legacy single-assignment
  hides a 5.4× swing, and even among valid samples topology varies ~12–17%. Density, by contrast, is
  *insensitive* to assignment (CV ≈ 0).
- **Cell policy is the dominant knob for density** (4.4×), but barely touches topology or survival.
- **Endpoint policy** is a minor knob here (short loading); **survival is robust** to every knob —
  the STEM→DDD import is a stable initial condition regardless.

## Decision value
The pilot tells you where to spend effort: **resolving the slip-system assignment ambiguity (experimental
g·b) is the single highest-value next step**, because it controls the topology objective by ~5×; cell
policy must always be *reported* (it sets density 4.4×) but is a modeling choice, not an unknown; endpoint
and zbox-within-periodic are second-order. This is the framework doing its job — quantifying and ranking
the experimental→simulation uncertainty so downstream DDD numbers come with an audit trail.

*Caveats:* small system (270 nodes), short relaxation/loading; density-growth was ~1.0 over 300 load steps
(no multiplication captured at this scale) so the density-growth-plausibility objective is not yet
discriminating — a longer-loading variant is the natural extension. These are sensitivity ratios, not
converged hardening numbers.

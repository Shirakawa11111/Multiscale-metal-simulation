# Real-network DDD audit (M3, v1.1)

Question: is the STEM → IDR → ExaDiS import **stable and auditable**, and how do downstream DDD outcomes
depend on the policies the IDR exposes? Runs on the real 27-line Cu network (`cu_stem_idr.json`, 270 nodes):
zero-stress relaxation + short stress-controlled loading. Harness: `real_network_audit.py`; reproducibility
in `AUDIT_MANIFEST.md`. This document gives the **corrected (v1.1)** conclusions; the v0 edgewise result is
recorded as a falsified artifact at the bottom.

## 1. Corrected results (line-coherent sampling)
- **The import is stable.** Across all policies the network survives zero-stress relaxation (segment count
  grows 1.8–2.6× via remesh; it does not delete away); line length relaxes to 0.53–0.84 of built value as
  the pinned reconstructed lines straighten under line tension. A viable DDD initial condition.
- **Assignment ambiguity is a MINOR knob for topology.** Junction nodes after loading, top1 vs
  `sample_edgewise` vs **`sample_linewise`** (the physical default), per cell/force:

  | cell/force | top1 | sample_edgewise | **sample_linewise** | linewise/top1 |
  |--|--|--|--|--|
  | foil + LineTension | 10 | 25.5 | **6.3** | 0.63 |
  | thickened + DDD_FFT | 2 | 29.5 | **7.2** | 3.58 |
  | thickened + LineTension | 6 | 26.3 | **5.5** | 0.92 |

  Line-coherent junction counts sit at ~top-1 level (ratio scattered 0.6–3.6 around 1, not 5×). So the
  geometric assignment ambiguity does **not** strongly amplify topology in this network.
- **Cell policy DOMINATES apparent density (~5.2×), deconfounded from force.** ρ after relax:
  foil+LineTension 1.16e13, thickened+LineTension 2.25e12, thickened+DDD_FFT 2.64e12 → the **cell** drives
  density 5.2× (foil→thickened at the *same* LineTension force); the **force** model is minor (1.17×).

## 2. The falsification (why v0 was wrong) — found & fixed
A review flagged that the v0 `sample` policy drew a slip system **independently per edge**, so adjacent
segments of the *same* reconstructed line could get different Burgers → artificial within-line
discontinuities → artificial junctions. Verified (`results_exadis/assignment_sensitivity.{json,md}`):
**`sample_edgewise` injects ~143/216 (66%) within-line Burgers discontinuities; `sample_linewise` injects 0**
(and `top1` 0). The v0 "5.4× topology swing" was therefore **predominantly an edgewise-sampling artifact**.
Fix: edges carry `parent_line_id`; `sample_linewise` makes one draw per line; `sample` is deprecated;
`sample_edgewise` is retained only as a discontinuity stress-test / upper bound.
Figure: `results_exadis/v11_linewise.png`; data: `results_exadis/v11_linewise_summary.json`; 39-run rerun in `v11/`.

## 3. What stands, what is bounded
- **Stands:** the import is auditable & stable; the assignment is genuinely ambiguous *as a label*
  (geometry fixes the {111} plane, the 3 ⟨110⟩ Burgers are near-degenerate: mean confidence ~0.33,
  ~1.58 bits — `cu_stem_idr_report.md`), but that ambiguity propagates **per line**, not per segment, and
  is a minor topology driver; **cell policy is the dominant density knob.**
- **Bounded / caveats:** small system (270 nodes), short relaxation/loading (density-growth ~1.0, no
  multiplication captured); these are sensitivity ratios, not converged hardening numbers. `as_is`+DDD_FFT
  is unsupported (FFT needs full PBC).

## v0 (SUPERSEDED — historical)
The original M3 pilot reported top1 ≈ 5 vs `sample` ≈ 27 junctions (~5.4×) and concluded "assignment policy
dominates topology." That used edgewise sampling and is retracted per §2 — kept here only as the record of
the self-correction.

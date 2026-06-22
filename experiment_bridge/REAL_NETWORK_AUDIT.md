# Real-network DDD audit (M3)

Question: is the STEM → IDR → ExaDiS import **stable and auditable**, and how do downstream DDD outcomes
depend on the assignment / cell policy that the IDR exposes? 10 runs on the real 27-line Cu network
(`cu_stem_idr.json`, 270 nodes), each: zero-stress relaxation + short stress-controlled loading.
Harness: `real_network_audit.py`. Figure: `results_exadis/real_network_audit.png`,
data: `results_exadis/real_network_audit_summary.json` + per-run `audit/*/audit.json`.

## 1. The import is stable (no collapse)
Across all policies the network **survives** zero-stress relaxation — segment count grows (1.8–2.6×, from
remesh splitting), it does not delete away. Line length relaxes to 0.53–0.84 of the built value: the
pinned-end reconstructed lines **straighten** under line tension between their anchors (physical, expected).
So the ingested experimental network is a viable DDD initial condition, not a degenerate one.

> **⚠️ SUPERSEDED by the v1.1 correction (see bottom).** Section 2's "5.4×" used *edgewise* sampling,
> which is now known to be a sampling artifact. Line-coherent sampling shows topology is **not** strongly
> assignment-sensitive. Read §"v1.1 correction" for the corrected conclusion.

## 2. [edgewise, superseded] Downstream topology under edgewise sampling
| assignment policy | junction nodes after loading |
|--|--|
| **top1** (legacy single argmin\|n·t\|) | **~5** |
| **sample** (the real 3-way candidate ambiguity) | **~27 (5.4×)** |

This is the central audit result: because every line's Burgers is geometrically 3-way degenerate, the
legacy top-1 converter places many lines on the *same* (first-listed) system → few cross-system reactions
→ it **under-counts junctions ~5×**. Sampling the true ambiguity produces a realistically reactive forest.
**Any junction / hardening number from a single top-1 assignment is not trustworthy as a point estimate.**

## 3. Apparent density is dominated by the CELL policy
| cell policy | ρ after relax | force |
|--|--|--|
| `as_is` (foil, z non-periodic) | 1.4e13 m⁻² | LineTension |
| `thickened_periodic` (zbox=5) | 3.2e12 m⁻² | DDD_FFT | 

A ~4.4× difference — the thickened periodic cell dilutes the density (5× z-volume) and uses full N-body
elastic force. So density must always be reported *with* its cell policy (see `CELL_POLICY.md`).

## 4. Sensitivity at fixed policy (the UQ envelope)
Across 4 assignment samples at a fixed cell policy: relaxed line length is stable (CV ≈ 0.01–0.06),
density CV ≈ 0.02–0.08, but junction count varies ±3–4 (~15%). So density is robust to the assignment
ambiguity; **topology is not** — it must be reported as a distribution.

## Verdict
The STEM-to-DDD import is **auditable and stable**, but two policy knobs the IDR makes explicit —
slip-system assignment (→ topology, 5×) and cell policy (→ density, 4.4×) — dominate the downstream DDD
outcome. This is exactly what the IDR was built to expose: the framework now **quantifies** where the
experimental-to-simulation uncertainty lives instead of hiding it behind a single forced choice.
*Caveat:* small system (270 nodes, short relaxation/loading); these are sensitivity ratios, not converged
hardening numbers.

---

## v1.1 correction — line-coherent assignment (the artifact, found & fixed)
A review flagged that `sample` drew a slip system **independently per edge**, so adjacent segments of the
*same* reconstructed line could get different Burgers → artificial within-line discontinuities → artificial
junctions. Verified: edgewise sampling produces **142/216 (66%)** within-line Burgers discontinuities;
`top1` and the new **`sample_linewise`** (one draw per parent line) produce **0**. Rerun (39 ExaDiS runs,
`v11/`, figure `results_exadis/v11_linewise.png`, data `results_exadis/v11_linewise_summary.json`):

| cell/force | top1 | sample_edgewise | sample_linewise | linewise/top1 |
|--|--|--|--|--|
| foil + LineTension | 10 | 25.5 | **6.3** | 0.63 |
| thickened + DDD_FFT | 2 | 29.5 | **7.2** | 3.58 |
| thickened + LineTension | 6 | 26.3 | **5.5** | 0.92 |

**Corrected conclusion.** The ~5.4× junction increase was **predominantly an edgewise-sampling artifact**.
With physically-coherent line-wise sampling, junction counts return to ~top-1 level (scattered around it,
ratio 0.6–3.6, not 5×). So **slip-system assignment ambiguity does NOT strongly amplify topology** once
sampled coherently — there is only a modest residual seed variation. `sample_linewise` is the correct UQ
protocol; `sample_edgewise` is retained only as a discontinuity stress-test / upper bound.

**Cell vs force deconfounded** (separating the two knobs that were coupled in §3): density after relax is
foil+LineTension 1.16e13, thickened+LineTension 2.25e12, thickened+DDD_FFT 2.64e12 → the **cell** drives
density **5.2×** (foil→thickened at the *same* LineTension force); the **force** model is minor (~1.17×).
So "cell policy drives density" stands, and it is the cell volume/periodicity, not the force model.

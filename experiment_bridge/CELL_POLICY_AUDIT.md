# Cell-policy density audit (M2/M3 follow-up, line-coherent)

Maps ρ_app and the robust observables vs the cell/force policy, at fixed physical `sample_linewise`
assignment. 35 ExaDiS runs (5 cell/force configs × top1 + 6 line-coherent seeds), conservative ≤30 cores.
Figure `results_exadis/cell_policy_audit.png`; data `results_exadis/cell_policy_audit_summary.json`.

## ρ_app after relaxation
| config | force | ρ (m⁻²) |
|--|--|--|
| foil (z=600) | LineTension | 1.25e13 |
| thickened z3 | LineTension | 4.36e12 |
| thickened z5 | LineTension | 2.72e12 |
| thickened z10 | LineTension | 1.46e12 |
| thickened z5 | DDD_FFT | 2.50e12 |

## Findings
- **Apparent density is a cell-NORMALIZATION artifact: ρ_app ∝ 1/zbox.** The line length is fixed; thickening
  the cell only inflates the volume, so ρ = length/volume scales as 1/zbox (z3→z5→z10: 4.36→2.72→1.46e12).
  foil vs thickened-z5 = **4.6×** at the *same* force. So a foil-network density is meaningless without its
  cell convention.
- **Force model is minor for density** (~8%: LT 2.72e12 vs FFT 2.50e12 at z5).
- **Survival and topology are ~cell-independent (robust):** segment survival ~2.0 and junction count ~5–8
  across all cells (DDD_FFT adds a few junctions vs LineTension); line-length relaxation 0.59–0.69.

## Verdict
The cell policy controls only the **density reporting convention** (ρ ∝ 1/volume), not physics. The
**cell-robust, reportable observables** are network survival, line-length relaxation, and junction count.
Any apparent-density number from a STEM foil network must be quoted with its cell convention (foil vs
thickened-periodic + zbox). *Caveat:* small system, short loading — sensitivity ratios, not converged numbers.

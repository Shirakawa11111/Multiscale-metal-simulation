# Multi-slip COLL on/off cell — conclusion (bounded negative)

> **One sentence:** In a co-driven, cross-slip-enabled, multiplying-FR-source multi-slip DDD cell
> (ExaDiS / FCC_0 / DDD_FFT), the **canonical collinear interaction dominance does NOT reproduce** as a
> flow-stress or storage/MFP effect — the collinear mechanism is real, but the collective coefficient is
> not recovered. This is a bounded negative, not a failure.

## What was tested
Per the adversarial scrutiny ([../evolving_forest_mfp/SCRUTINY_VERDICT.md](../evolving_forest_mfp/SCRUTINY_VERDICT.md)),
the decisive test: a primary slip system with **multiplying Frank–Read sources** gliding through an
**evolvable forest** of a partner system, under **strain-rate control** (flow stress = output), with
**CrossSlip** enabled, comparing junction-type *pairs* at matched density. Decision instrument
(`multislip_flow.py`) gates: **co-drive** (`EDIR_MODE=opt_pair`, both systems Schmid ≈ 0.41), per-system
density ledger, partner-forest drift, plateau, ambiguous-fraction, carrier-starvation. The verdict uses
**resolved shear** τ_RSS = σ_flow·|S_primary| (NOT raw scalar — `opt_pair` gives each pair a different
Schmid), and the **measured** forest density, not the target.

## Result — density lever (measured ρ_f 1.65e13 → 1.5e14, ~9×)

`density_lever.png`, `density_lever_results.csv`, `density_lever_summary_XSLIP{on,off}.json`.

| R_RSS = τ_RSS(coll)/τ_RSS(glissile) | ρ_f≈1.6e13 | ρ_f≈1.5e14 | density gain G_τ |
|--|--|--|--|
| XSLIP on  | 0.963 | 0.989 | 1.03 |
| XSLIP off | 0.954 | 1.008 | 1.06 |

- **coll/glissile ≈ 1.0 at both densities, both cross-slip settings, FLAT** (no density-threshold growth).
  Canonical target is ≈ √(0.62/0.12) ≈ **2.3 and growing** — not observed.
- **coll_opp / coll_same = 1.00** everywhere (the annihilation single-bit toggle is null even co-driven).
- **coll / hirth ≈ 1.1–1.19** (Hirth weakest — canonical *direction* correct, but small, not dominant).
- **XSLIP on vs off: identical** — the "cross-slip replenishment" mechanism does not manifest.
- L_mf / storage agrees: collinear MFP is *longer* (weaker storage), not shorter.

## Five independent measurements, all consistent
1. Local remobilization strength (`../binary_reaction_matrix`, SEALED): collinear not dominant (mid-pack).
2. Pairwise MFP (`../evolving_forest_mfp`, RETRACTED): sampling-starved; geometry can't host the coefficient.
3. Multi-slip flow, RSS, density lever, XSLIP on: coll/glissile ≈ 1, no gain.
4. Same, XSLIP off: identical.
5. Annihilation toggle (coll_opp vs coll_same): null at every condition.

## Honest limitations (why "bounded", not absolute)
- **Drag-dominance:** at fixed erate the *absolute* flow stress is kinematic (∝1/ρ_mobile — it *dropped*
  51→9 MPa as density rose 9×, backwards for Taylor). Only the type-*ratio* (cancels common drag) is usable;
  a true Taylor hardening measurement needs rate-extrapolation / quasi-static (cf. `../../taylor_hardening`).
- **Collinear drift gate:** the co-driven collinear partner structurally deforms 8–22% (it cannot be held as
  a fixed-density forest — collinear systems co-deform strongly), so the *formal* gate reads AMBIGUOUS at low
  density. The conclusion is robust to this: the result is identical at the readable high density and across
  XSLIP, the low-density drift is mixed-sign and small under XSLIP off, and where positive it *inflates*
  collinear yet collinear stays ≤ glissile.
- Single mobility law (FCC_0), single FR-source model, finite cell, erate=1e6 (rate-dominated).

## Interpretation
The canonical Madec/Devincre collinear coefficient appears to depend on specifics not captured here —
quasi-static loading, the particular nodal/mobility/cross-slip rules, or fully-developed multi-slip
population statistics — rather than being a generic, readily-reproducible flow-stress observable in this
ExaDiS protocol. The collinear *mechanism* (partial annihilation + length scaling τ_c=K(μb/l)ln(l/b),
R²=0.89) is solid; the *collective dominance* is not reproduced by any of the local, pairwise, or
multi-slip measures attempted.

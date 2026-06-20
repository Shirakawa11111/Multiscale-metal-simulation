# evolving_forest_mfp — progress & honest retraction

## Goal
Test whether canonical collinear DOMINANCE is a collective mean-free-path / storage effect
(since the collinear LOCAL strength was shown NOT dominant in `../binary_reaction_matrix`).

## Pairwise assay (attempts 1–3) — SUPERSEDED, do not extend
Built a pair-resolved evolving-forest assay (mobile system m + evolvable forest system f,
DDD_FFT+Collision+Topology, full PBC). Forest-alone density-stable (<1% drift) ✓. Mechanism
(2-segment collinear annihilation, opp vs same) validated ✓.

- **Attempt 1** (FR sources, τ=120): mobile multiplied ~12×; opp consumed only ~8% more mobile,
  late and post-wrap. Multiplication swamps annihilation.
- **Attempt 2** (FR, τ=40 density scan): mobile barely flowed (γ~5e-6); high ρ_f forest self-annihilates.
- **Attempt 3** (non-multiplying mobile LINES, τ=60, ρ_f=3e12): opp==same (rel d_mob≈0%).

### ⚠️ RETRACTION (4-expert adversarial scrutiny, `SCRUTINY_VERDICT.md`)
The attempt-3 "null" is **NOT a valid negative** and is retracted: the mobile glided across only
~0.31 forest cells (−1.8% decay = essentially no glide) → **sampling-starved**, not a real null.
The pairwise geometry is **structurally incapable** of hosting the multi-slip collinear coefficient
(collinear b3=0 → no storage channel, no carrier replenishment, no co-driven coefficient). Every
pairwise re-tune trades one confound for another. **Do not run more pairwise scans.**

## Pivot (current mainline): multi-slip COLL on/off cell — `../multislip_flow/`
The single test that cannot return ambiguous (per scrutiny):
- Small FCC cell, full PBC, DDD_FFT + Collision + Topology + **CrossSlip** (required for replenishment),
  **strain-rate control** (flow stress is an OUTPUT), **multiplying FR sources** (insert_frank_read_src).
- Compare steady-state τ_flow and L_mf = 1/(b·dρ_stored/dγ) across forest types at matched density:
  **collinear vs glissile vs Hirth duplex**, ρ_f ∈ {3e12, 3e13}.
- **CONFIRM** collinear dominance: collinear τ_flow highest, gap grows with ρ, τ_flow ratio →
  canonical √(0.62/0.12) ≈ 2.3× over glissile.
- **REFUTE**: at a flow plateau (γ>1e-3), all gates passed (density-drift<5%, carrier-symmetry,
  plateau reached), collinear == glissile/others at both densities → canonical a_coll≈0.62 does not
  reproduce in this ExaDiS/FCC_0 setup; publish the bounded negative.

Gates before reading any number: density-drift <5%; carrier-symmetry; γ>1e-3 steady plateau.

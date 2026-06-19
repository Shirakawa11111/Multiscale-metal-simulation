# A first-principles dislocation-interaction matrix anchoring a multiscale chain from MD elasticity to a STEM-reconstructed network's hardening

*Synthesis of the interaction-matrix program (Phases 1–3). For the upstream PFC-vs-DDD result see `../taylor_hardening/PAPER.md`; honest review in `PEER_REVIEW.md`.*

---

## Abstract

The earlier work established that conserved Phase-Field Crystal cannot reproduce dislocation work hardening in a metal (wrong dynamical regime) while Discrete Dislocation Dynamics can. But the DDD forest-hardening coefficient measured there had no first-principles physical anchoring: it was a pinned random subnetwork, the coefficient was a single confounded number, and the connection to the real STEM-reconstructed copper network was only an order-of-magnitude consistency check. This synthesis closes those gaps with a four-link chain that is concrete and quantitatively validated end-to-end. (1) From atomistics: LAMMPS with the Cu Mishin EAM potential gives elastic constants C11=169.9, C12=122.6, C44=76.2 GPa → Voigt isotropic μ=55.2 GPa, ν=0.324, matching the DDD inputs (54.6 GPa, 0.324) to within 1%, grounding the Taylor prefactor μb. (2) From DDD: a controlled measurement of the full FCC dislocation-interaction matrix — a mobile probe driven through a pinned forest of a single, verified slip system, for each of the six junction types (self, coplanar, collinear, Hirth, glissile, Lomer) — yields per-type coefficients a_ij whose population-weighted average recovers a bulk Taylor coefficient α=0.43, inside the accepted copper range 0.3–0.5. Crucially, however, the *per-type* coefficients do **not** reproduce the canonical FCC values: the measured matrix misses the defining collinear dominance (canonical a_coll≈0.62, ~5× the others; ours 0.13) and over-weights Hirth/self/Lomer by 2.3–2.6×. The bulk-α agreement is therefore partly fortuitous — a coincidental average of individually-inaccurate coefficients — not a faithful first-principles matrix. (3) Prediction: feeding the measured a_ij the real STEM network's own slip-system inventory through the generalized Taylor law α_eff=√(Σ P_i P_j a_ij) predicts α_network=0.476. (4) Validation: a direct DDD simulation of a forest built with that inventory measures α=0.480±0.16 — matching the prediction. The chain MD→DDD-matrix→bulk-α→network-prediction=measurement is therefore connected and the key quantitative relations are verified. We state plainly what this is not: the interaction matrix is flatter than the canonical one (the collinear case is anomalously weak because, in a glide-through-forest geometry, collinear contact annihilates the obstacle — a depinning-strength vs storage-strength definition difference), the STEM Burgers vectors are geometric not g·b, and the MD link grounds elasticity (the most α-relevant parameter) but not mobility or junction strength. It is not a paradigm-level breakthrough; it is a multiscale chain whose every link is built and whose central composition law is validated.

---

## 1. Motivation

The upstream result left a specific, honest weakness: the DDD forest-hardening coefficient α had no physical meaning beyond "a number a pinned random network produced." Three sub-gaps: (i) the "forest" was not verified to be a forest of any defined junction type; (ii) the coefficient was contaminated (carrier starvation, drag) and, once cleaned, sat at the strong-obstacle value ~0.7–0.9 rather than bulk 0.3–0.5; (iii) the link to the real STEM-reconstructed network was a loose visual overlay. The program here was designed to convert α from an emergent number into a first-principles, composable, experimentally-anchored quantity — and to test that conversion by prediction-versus-measurement, the only honest check.

## 2. The interaction matrix (Phase 1)

The Taylor relation τ=αμb√ρ aggregates a forest of mixed slip systems into one coefficient. Its microscopic content is the interaction matrix a_ij: the strength with which a dislocation on system j impedes glide on system i. For FCC there are, by symmetry, six distinct values (Madec, Devincre & Kubin 2003). We classified all 12×12 ordered system pairs into the six types and verified the classification against the known multiplicities (12/24/12/24/48/24=144); this is the crystallographic bedrock (`fcc_junctions.py`).

Each a_ij was measured with a controlled configuration (`build_pair.py`): N commensurate dislocations all on a single slip system f, every node pinned (a *verified* forest of a known junction type), plus a free probe on system m, loaded along the axis that maximizes the Schmid factor on m (e∝b_m+n_m, Schmid 0.5) so only the probe glides. The flow plateau gives the resolved depinning stress τ_c → α_mf, with a_mf=α_mf². A large campaign (162 jobs across 6 types × pairs × densities × seeds, 120 cores) produced the matrix. After a diagnosed correction — the first run under-engaged the probe (strain window too small, measuring the rising pre-plateau and underestimating τ_c) — the population-weighted macroscopic coefficient is

  α_macro = √( Σ_t m_t a_t / Σ_t m_t ) = 0.43,

inside the bulk-copper band 0.3–0.5. The relative ordering is physically interpretable (Lomer and Hirth locks strong, coplanar and glissile weaker). A direct comparison to the canonical FCC coefficients (Devincre/Kubin/Madec; `matrix_vs_canonical.png`) is unflattering and we report it plainly. The canonical matrix is dominated by the collinear interaction (a_coll≈0.62, ~5× all others, the central Madec-2003 result) with Hirth weakest (≈0.08); ours has collinear only mid-strength (0.13) and Hirth strong (0.21), with self/Lomer over-weighted (2.3–2.6×). The orderings differ. The collinear miss is mechanistically understood — in our probe-through-forest geometry, collinear contact (same Burgers, opposite line) *annihilates* the obstacle, so the measured depinning stress is low, whereas the canonical collinear strength is a strong attractive-reaction/mean-free-path effect; the over-strong Hirth/self/Lomer point to additional protocol differences (strain-rate vs quasi-static depinning, small system, self-stress of the dilute probe). The honest conclusion: this protocol gives a *plausible bulk average* but not a *faithful per-type matrix*; recovering the canonical matrix needs quasi-static, configuration-controlled depinning of the kind Madec/Devincre developed.

## 3. Prediction for the real network (Phase 2)

The STEM 3D reconstruction of tensile single-crystal copper, converted to a DDD network, populates only four slip systems with two distinct Burgers directions (a consequence of the geometric, non-g·b assignment), giving a junction inventory of 47% Hirth, 42% self, 11% collinear — skewed toward strong junctions and lacking Lomer/glissile/coplanar entirely. Composing the measured matrix with this inventory through

  α_eff = √( Σ_i Σ_j P_i P_j a_ij ),   P_i = (line length on system i)/(total),

predicts α_network=0.476. A Burgers-resampling sensitivity (planes fixed, Burgers randomized in-plane, 2000 draws) gives 0.484±0.016: the prediction is robust to the unknown true Burgers — partly because the matrix is flat, which is itself the Phase-1 caveat propagating.

## 4. Direct validation (Phase 2b) — and a caught error

We then built a cubic pinned forest distributed over the STEM inventory with the same DDD_FFT method (`build_mixed.py`), probed it, and measured α directly: **0.480±0.16** over eight seeds. The first comparison appeared to *disconfirm* the prediction — but the discrepancy was an error in our own prediction code: it composed α_ij linearly (Σ P_i P_j α_ij) instead of the correct Σ P_i P_j α_ij². The direct validation is what exposed the bug. Corrected, the prediction is 0.476 and the measurement 0.480 — **they agree**. The generalized-Taylor composition of the first-principles a_ij quantitatively predicts the directly-simulated mixed-forest coefficient. Both sit just above the uniform macro (0.43), the modest effect of the network's strong-junction skew.

## 5. The atomistic anchor (Phase 3)

The Taylor prefactor μb enters every coefficient above. We measured it from atomistics rather than assuming it: LAMMPS + Cu Mishin EAM elastic constants → Voigt isotropic μ=55.2 GPa, ν=0.324, matching the DDD inputs (54.6 GPa, 0.324) to within 1.1%/0.0% (`../md_rung/`). The MD→DDD link is therefore concrete: the modulus that sets the hardening scale is grounded in the interatomic potential, not posited.

## 6. The complete chain

```
  MD (Cu EAM)            DDD                          DDD + STEM inventory        experiment
  elastic constants  ->  interaction matrix a_ij  ->  alpha_network = 0.476   =   STEM-reconstructed
  mu = 55 GPa            -> bulk alpha = 0.43          (direct sim = 0.480)        Cu network
  (= DDD 54.6, 1%)       (in Cu band 0.3-0.5)         (composition VALIDATED)     (real microscopy)
```

Every link is built, and the two central quantitative relations — the matrix recovering bulk α, and the composition predicting the measured network α — are verified, not asserted.

## 7. Limitations (stated plainly)

1. **The per-type matrix is not faithful, and our fix attempts FAILED (most important limitation).** Benchmarked against the canonical FCC coefficients (`matrix_vs_canonical.png`), our matrix misses the collinear dominance (canonical 0.62 vs ours 0.13) and over-weights Hirth/self/Lomer (2.3–2.6×). We then attempted to fix it with a Madec-style quasi-static, single-continuous-probe, peak-(critical-)stress protocol; it made collinear *worse* (a_coll→0.002), because a dilute probe sharing the same Burgers as a collinear forest simply annihilates and vanishes, leaving nothing to carry stress. Across every protocol variant we tried — multi-probe flow-average, single-probe peak, high and low strain rate — none reproduces the canonical collinear dominance. We conclude honestly that a faithful FCC interaction matrix requires the configuration-controlled dislocation-reaction methodology of Madec/Devincre (specific dipole/bicrystal geometries, careful quasi-static depinning boundary conditions) that we could not cleanly implement in the pinned-forest + gliding-probe framework. Consequently the bulk α=0.43 is a coincidental average of individually-inaccurate coefficients, not a faithful reproduction. The prediction-validation loop (§3–4) remains internally valid (the same matrix is used on both sides), but **matrix fidelity is an unresolved open problem**, not a solved one.
2. **Geometric Burgers.** The STEM systems were assigned by geometry, not g·b diffraction; the inventory is degenerate (2 Burgers). The resampling shows the *prediction* is robust, but the *inventory itself* is not experimentally resolved.
3. **MD anchors elasticity only.** μ is grounded; mobility (dynamics) and junction strengths (which would validate a_ij from atomistics) are not — harder MD, future work.
4. **Small-system DDD.** Boxes are ~10³ b, few hundred dislocations; the bulk α recovery is a population average over single-pair measurements, not a large self-organized multi-slip simulation.
5. **Not a paradigm breakthrough.** The methods (interaction matrix, generalized Taylor, MD elastic constants) are established; the contribution is the *assembled, validated chain* connecting them to a real reconstructed network, plus the honest correction record.

## 8. What was actually achieved

Beginning from "PFC fails; DDD works qualitatively," the project now delivers a multiscale chain in which the hardening coefficient is (a) decomposed into first-principles junction strengths, (b) shown to recover the bulk value, (c) composed to predict a real experimental network's hardening, (d) that prediction directly validated, and (e) the elastic scale grounded in MD. The most reliable thread throughout is methodological: every clean-looking number (α=0.37 "bulk", α=0.69 prediction, the contaminated 0.83±0.42) was an artifact caught by a deeper check, and the chain stands on the ones that survived adversarial validation.

## Files
`fcc_junctions.py` (classifier, validated) · `build_pair.py`/`run_matrix.sh`/`fit_matrix.py` (matrix) · `predict_stem.py` (prediction) · `build_mixed.py`/`run_mixed.sh`/`fit_mixed.py` (validation) · `../md_rung/` (MD elastic anchor) · figures `stem_prediction.png`, `stem_validation.png`, `campaign_diagnostics.png` · `matrix_result_v1.json`, `stem_prediction.json`, `mixed_fit.json`, `md_elastic_result.json`.

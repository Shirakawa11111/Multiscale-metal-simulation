# Local binary-reaction benchmark (SEALED)

> **One sentence**: Controlled binary reactions confirm the correct mechanisms and the
> collinear remobilization length scaling, **but the local remobilization matrix does NOT
> reproduce collinear dominance** — so collinear dominance is a *collective forest / mean-free-path*
> effect, pursued in [`../evolving_forest_mfp/`](../evolving_forest_mfp/). This benchmark is sealed;
> it is a mechanistic benchmark / negative discriminator, **not** the mainline.

## What is solid (control-validated)

| Result | Value | Evidence |
|--|--|--|
| collinear mechanism | partial annihilation, line→0.71, no stable junction | `binary_collinear.py`, observer |
| collinear remobilization scaling | τ_c=K(μb/l)ln(l/b), K≈0.75, R²=0.89 | `collinear_scaling.png` |
| local ranking (6 types, mechanisms confirmed) | self 225 > collinear≈glissile≈coplanar 175 > Lomer≈Hirth 125 MPa | `valid_ranking.png` |
| **collinear is NOT locally dominant** | mid-pack, not ~5× | `valid_ranking.png` |

## Four controls (all pass) — `local_controls.png`

1. **LineTension null** — without pairwise elastic interaction the **junctions do not form**
   (glissile/Hirth anneal→1.00 = no reaction); collinear still annihilates (collision-driven) but
   τ_c falls 175→100. ⇒ the junction mechanisms genuinely require DDD_FFT elastic interaction; not artifacts.
2. **rann / maxseg convergence** — τ_c stable for MAXSEG≥80 (collinear 180, Hirth 100) and robust to
   rann=5/10/20; only the coarse MAXSEG=40 edge differs. ⇒ converged at our standard MAXSEG=80.
3. **box convergence** — collinear τ_c=180 at L_box/L_seg = 4, 6, 8 (flat). ⇒ τ_c is **not** an
   artifact of the `spanB>0.45·LBOX` remobilization criterion.
4. **orientation family** — 6 crystallographically distinct collinear families give τ_c=160±20 (12% spread),
   all annihilate. ⇒ "collinear not locally dominant" is geometry-robust, not a single-pair fluke.

## Why sealed (not mainline)
The local strength matrix is now a clean, converged, robust quantity — and it **demonstrably cannot**
explain canonical collinear hardening. Continuing to refine `a_ij_local` would only make a settled
conclusion prettier. The missing physics is **storage / mean-free-path**, measured in the new mainline.

Artifacts: `binary_local_summary.json`, `valid_ranking.png`, `collinear_scaling.png`, `local_controls.png`.

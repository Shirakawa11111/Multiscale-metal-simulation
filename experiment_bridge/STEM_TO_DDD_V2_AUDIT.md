# STEM-to-DDD v2 — audit package (IDR v1.1)

A single, reportable, reproducible audit of the STEM→IDR→DDD pipeline on the real Cu reconstruction. The
contribution is **not** a hardening curve — it is an *auditable* path from a microscope reconstruction to a
DDD initial condition, with every assumption (slip-system assignment, cell policy, density convention,
endpoint policy) exposed, quantified, and — where it was wrong — self-corrected.

> One-line thesis: *geometric slip-system assignment is highly uncertain, but propagated in a physically
> coherent (line-wise) way it is a minor topology knob on this network; the cell policy sets only the
> reported density convention, while survival / line-length relaxation / junction count are the robust,
> cell-independent observables.*

---

## 1. Input — real STEM reconstruction
- 27 reconstructed dislocation lines → **270 vertices, 243 edges** (`cu_stem_idr.json`).
- Each edge carries `parent_line_id` (which physical line it belongs to) — the hook for line-coherence.
- Cell: foil, z non-periodic (`foil_nonperiodic_z`); z-depth uncertainty ±30 nm; endpoints pinned
  (truncated reconstruction). Source: `experiment_bridge/recon_data/points_3d*.txt`.

## 2. IDR uncertainty (geometry-only, pre-g·b)
Geometry fixes the {111} glide plane (`argmin|n·t|`), but the three ⟨110⟩ Burgers in that plane are
near-degenerate without diffraction contrast:

| metric | value |
|--|--|
| mean assignment confidence | **0.333** (≈ 1/3 — three-way tie) |
| mean assignment entropy | **1.585 bits** (= log₂3) |
| edges ambiguous (entropy > 0.8) | **243 / 243** |
| Burgers source | `geometry_only_pending_gb` |

→ The IDR does not hide this: every edge ships a 3-candidate set with priors, not one forced choice.

## 3. Lowering correction — the self-falsified artifact
The v0 pilots sampled assignments **per edge**. Quantified on this network (`assignment_sensitivity.json`):

| policy | within-line Burgers discontinuities | of intra-line adjacencies |
|--|--|--|
| top1 | 0 | 216 |
| **sample_edgewise (deprecated)** | **142.8** | 216 |
| **sample_linewise (default)** | **0** | 216 |

Edgewise sampling manufactured ~143 Burgers/plane discontinuities *inside single physical lines* → artificial
junctions → the v0 "assignment → topology **5.4×**" result. **Line-coherent** sampling (`sample_linewise`,
one draw per `parent_line_id`) injects **zero** within-line discontinuities and is now the default. The
ambiguity is still propagated (inventory CV 0.62 across seeds) — but per line, not per segment.

## 4. DDD stability (real-network audit, v1.1-corrected)
With line-coherent lowering, the imported network is a **stable** ExaDiS initial condition:
- **Network survival** ≈ 2.0 (segment count grows under remesh, does not collapse).
- **Line-length relaxation** to **0.59–0.69** of the as-built length under zero-stress relax — relaxes, no
  collapse.
- **Junction count** ≈ **5–8** (line-coherent), i.e. **≈ top-1**, not the 27 of the edgewise artifact.

→ Corrected conclusion (`REAL_NETWORK_AUDIT.md`): assignment ambiguity is a **minor topology knob** once
propagated per line; it does **not** drive a large topology swing.

## 5. Cell policy & density convention
The dominant knob is the **cell policy**, and it controls only the **density reporting convention**
(`CELL_POLICY_AUDIT.md`, `density_conventions.md`):

| config | ρ_app (relaxed) | | convention | foil-aware density |
|--|--|--|--|--|
| foil z600 / LineTension | 1.25e13 | | **Λ_A (projected, foil-native)** | **3.26e6 m⁻¹** |
| thick z3 / LT | 4.36e12 | | ρ_vol foil (z600) | 2.13e13 m⁻² |
| thick z5 / LT | 2.72e12 | | ρ_vol thick-z5 | 4.25e12 m⁻² |
| thick z10 / LT | 1.46e12 | | (ρ_vol = Λ_A / z_eff) | ∝ 1/z |
| thick z5 / DDD_FFT | 2.50e12 | | | |

- **ρ_app ∝ 1/zbox** — a pure cell-volume normalization (line length fixed, volume thickened): foil vs
  thick-z5 = **4.6×** at the *same* force. **Force model is minor** (LT vs FFT = 0.92×, ~8%).
- **Foil-native observable is Λ_A** = line length / projected area = **3.26e6 m⁻¹** (z-independent). The
  bulk-equivalent ρ = Λ_A / z_eff is a derived, convention-dependent number — always quote it with z_eff/zbox.
- Cross-check: as-built ρ_foil 2.13e13 × relax_len_frac(~0.59) ≈ **1.26e13** ≈ the relaxed foil 1.25e13. ✓
- **Cell-robust observables:** survival (~2.0), line-length relaxation (0.59–0.69), junction count (5–8).

## 6. g·b-ready path (line-coherent)
The IDR upgrades from geometry-only to physics-validated when diffraction-contrast data arrives, via
`uncertainty.apply_gb_constraints` (invisibility |g·b|≈0). Synthetic demo (`synthetic_gb.md`) with **one
line-coherent** ground-truth Burgers per parent line:

| reflections | mean entropy | frac resolved |
|--|--|--|
| none | 1.585 bits | 0% |
| 1 (g=200) | 0.704 bits | 30% |
| 2 (g=200,020) | 0.0 bits | 100% |

→ ~2 well-chosen reflections collapse the assignment to a single slip system (`gb_validated`). Real data
plugs into each line's `uncertainty.gb_constraints` as observed (g, visible) pairs — same code path.

---

## Reportable observables (the three classes)
This audit fixes what BO/UQ should be defined over (see `BO_UQ_PILOT.md` → UQ objective freeze):

- **A — stability:** network survival, segment-count growth, line-length relaxation. *(robust to all knobs)*
- **B — topology:** junction count, topology-event proxy, **within-line discontinuity ≡ 0** (a hard
  legality invariant, not a free knob). *(minor, line-coherent-sensitive)*
- **C — reporting:** **Λ_A projected line density (foil-native)**; ρ_vol only with a declared cell convention.
  *(a convention, not physics)*

## Reproduce
```
python3 experiment_bridge/stem_to_idr.py            # recon -> cu_stem_idr.json + report
python3 experiment_bridge/assignment_sensitivity.py # edgewise(142.8) vs linewise(0) within-line discont
python3 experiment_bridge/density_conventions.py    # Lambda_A (foil-native) + rho_vol by convention
python3 experiment_bridge/synthetic_gb.py           # line-coherent g.b entropy collapse 1.58->0.70->0
# DDD (HPC, <=30 cores): real_network_audit.py (M3), cell_policy_audit (M5) — see those docs
```
All local steps are deterministic; DDD steps record seed + config in their summary JSON (`AUDIT_MANIFEST.md`).

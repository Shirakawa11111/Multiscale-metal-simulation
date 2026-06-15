# Editorial Decision & Multi-Reviewer Report

**Manuscript:** *Regime-correct mesoscale plasticity: why Phase-Field Crystal cannot work-harden a 3D metal, and how Discrete Dislocation Dynamics seeded from STEM tomography quantitatively can*

**Review mode:** full (5 independent reviewers: EIC + Methodology + Domain + Perspective + Devil's Advocate)

**Editorial Decision: MAJOR REVISION** (unanimous across all five reviewers)

---

## 1. Score summary

| Dimension | EIC | Methodology | Domain | Perspective | DA |
|---|---|---|---|---|---|
| Originality | 7 | 7 | 6.5 | 5–7 | — |
| Methodological rigor | 6 | 4 | 4.5 | 5–6 | — |
| Evidence sufficiency | 5 | 4 | — | — | — |
| Argument coherence | 8 | 7 | 5.5 | 7 | — |
| Writing quality | 9 | 8 | — | — | — |

All five recommend **Major Revision**. None recommends reject; all judge the conceptual core (regime-mismatch thesis + carrier-confound finding + image→DDD route) genuinely worth publishing. The fixes are evidentiary and claim-calibration, not fatal.

**Data integrity note (raised independently by EIC, Methodology, Domain, DA):** every headline number reproduces *exactly* from the committed JSON files. The author does not fudge — even the embarrassing R²=−4.5 large-box through-origin fit is reported faithfully. The problems are over-interpretation and experimental design, not honesty of reporting.

---

## 2. CRITICAL issues (consensus; Devil's Advocate CRITICAL → Accept is blocked)

**[CR-1] The fit framing is backwards — the intercept-fit slopes are the bulk-band coefficients; the headline through-origin α=1.39 is the mis-specified fit.** *(Domain MAJOR-2, Methodology C2, DA C1/C2, EIC W1)*
The community defines the Taylor coefficient as the √ρ slope with the athermal/friction offset τ₀ separated (Devincre–Hoc–Kubin 2008; Kubin–Devincre–Hoc 2008). The intercept fits already give: small box α′=1.20 (τ₀=11 MPa), large box α′=0.52 (τ₀=64 MPa), and — decisively — the dismissed constant-fraction "control," re-fit with an intercept, gives **α′=0.31, τ₀=46 MPa, R²=0.97** (Domain reviewer recompute), i.e. a clean bulk-Cu slope. The paper leads with the wrong (through-origin) construction and explains its inflated value by the Orowan limit, when its own data already sit at 0.3–0.5 once τ₀ is separated.

**[CR-2] The "convergence into the bulk band" claim rests on a rejected fit + a dominant offset.** *(DA C1, Methodology C2, Perspective M3)*
The large-box through-origin fit gives R²=−4.5 (worse than the mean); the celebrated α′=0.52 exists only under a 2-parameter fit to 4 points whose 64 MPa intercept is 63–80% of the signal, while τ varies only 1.28× as ρ varies 6.2×. Per-point through-origin α falls monotonically (2.81→1.44), opposite to constant √ρ scaling. As written, the single most important inference ("α magnitude is scale, not principle") is not supported by that dataset alone.

**[CR-3] The constant-fraction "control" is mislabeled, and its number is misquoted.** *(Domain C-1)*
The abstract/§3.4 report R²=−0.94; the stored value is **−9.40**. More importantly, the series is *not* flat — with an intercept it yields a bulk-band slope (CR-1). The "constant fraction = carrier confound = meaningless flat response" narrative draws the opposite of the correct conclusion from the data.

**[CR-4] The STEM "validation" is partly circular and oversold.** *(Perspective C1/C2, DA M4)*
The hardening run loads along `edir=[0,0,1]` — the foil-normal, *synthetic* axis whose z-coordinate carries ~1% of the in-plane signal and is stretched ~700× by the adapter; the project's own route notes record that x-tension gives Schmid≈0 on the geometric slip systems, so the load axis was effectively chosen to make them resolve. "Closes the loop from experiment to a validated law" must be demoted to a feasibility/ingestion demonstration; α_eff≈1.0 is an order-of-magnitude sanity check, not validation.

---

## 3. MAJOR issues (consensus)

- **[M-A] No replication; n=4, single hard-coded seed (1234).** No error bars on α; one point (α=1.73) carries the small-box R². "Validated law" overstates this evidence base. *(All reviewers)*
- **[M-B] ρ_forest computed by node fraction (`rho_total·pinned/n_nodes`), not pinned line length; the pinned set is not verified to thread the probe glide planes** — so "ρ_f" is total pinned density, not forest density, and α is not the forest-interaction coefficient. *(Domain C-2, Methodology M3)*
- **[M-C] "PFC fundamentally cannot harden" overgeneralizes from enumeration to a universal.** Must engage Skaugen–Angheluta–Viñals (PRL 121, 255504, 2018; glide/climb timescale separation) and APFC (Salvalaglio–Elder); soften to a scoped, mechanism-specific claim about the absence of dislocation *storage/junction-locking* under the conserved variants tested. *(Domain MAJOR-1, DA M1)*
- **[M-D] KM ρ_ss(ε̇) is never measured (single strain rate); two-sided convergence mixes incompatible runs** (annihilation arm at erate=3e4→ρ_ss≈2.4e12 vs the STEM multiplication arm at erate=1e4→7.7e12, ~3× apart). Soften to "convergence at fixed rate, consistent with KM," or measure ≥2 rates. *(Domain MAJOR-3/4, Perspective M2)*
- **[M-E] Strain-rate drag may contaminate the flow-stress plateau itself**, not just the offset (2 probes carry the whole rate at v≈ε̇/(ρ_mobile·b); τ_drag=Bv/b is forest-independent and rate-dependent). Run a rate sweep at fixed ρ_f. *(Methodology M2, DA alt-explanation 4)*
- **[M-F] PFC negative result shown with no in-text figures**; PFC numerics under-specified (no grid, Δt, ψ̄, r, no resolution/timestep convergence showing softening is not a discretization artifact). *(EIC W3, Methodology M4)*
- **[M-G] Missing references** defining the very quantities measured. *(Domain)*

---

## 4. New data resolving part of the review (added this revision cycle)

The author ran the decisive control the reviewers requested (Methodology M1, "the one control that would settle it is named but not run"): a **no-forest, 2-probe run (ρ_f=0)** directly measuring the carrier baseline τ₀:
- small box τ₀ = 32.8 MPa; large box τ₀ = 72.1 MPa (rises with box size → confirms the carrier/drag interpretation **by measurement**, not assertion).
- Subtracting the *measured* τ₀ from the forest series: large box → α = 0.37 (R²=0.89), per-point [0.25, 0.29, 0.37, 0.41] — **in the bulk-Cu band**; small box → α = 0.79 (noisier; its baseline ≈ its lowest forest point).

This converts CR-1/CR-2/M-A(partly) from "asserted" to "measured," and supports the *corrected* headline the reviewers point to: **DDD recovers a bulk-Cu forest-hardening coefficient (≈0.3–0.4) once the measured carrier/friction baseline is separated.** It does not fix n=4/single-seed (M-A) or the ρ_f definition (M-B).

---

## 5. Disagreements among reviewers

- **Significance of the DDD half.** Domain reviewer scores it "incremental, below state of the art" (reproduces Arsenlis 2007 / Devincre 2008 on a small scale); EIC/Perspective see the *synthesis + STEM route* as the real, publishable contribution. **Editorial arbitration:** both are right — reframe the contribution as (i) the scoped PFC negative result and (ii) the image→DDD feasibility route, not as "DDD quantitatively reproduces Taylor" (which is partial and confounded).
- **Whether α≈1 is tautological.** DA argues pinning the forest makes τ∝√ρ near-geometric; Domain argues the coefficient still encodes junction statistics if ρ_f is the true forest density. **Arbitration:** the baseline-subtracted bulk-band coefficient (§4) is harder to dismiss as tautology, but the ρ_f-definition fix (M-B) is required before the coefficient is claimed to be the forest-interaction α.

---

## 6. Revision Roadmap (prioritized)

**Must-fix before resubmission (gates acceptance):**
1. **Re-frame the Taylor analysis (CR-1, CR-3, M-B).** Lead with intercept fits + the measured baseline (§4); report α=slope and τ₀=offset separately; correct −9.4; retract the "constant-fraction = flat confound" claim and report that series' bulk-band slope; recompute ρ_f as pinned line length and report the pinned set's slip-system composition.
2. **Rescope the convergence claim (CR-2).** State the large-box through-origin fit is rejected (R²=−4.5); present the baseline-subtracted α=0.37 as the bulk-band evidence, with the n=4/single-seed caveat front and centre.
3. **Demote the STEM claim (CR-4).** "Feasibility/ingestion route," not "validation/closing the loop"; disclose the loading-axis/z-fabrication; present α_eff≈1.0 as order-of-magnitude only.
4. **Soften the PFC claim (M-C).** Engage Skaugen 2018 PRL + APFC; change "fundamentally cannot" (title included) to a scoped statement about absent storage/junction-locking under the tested conserved variants.
5. **Add ≥1 PFC softening figure + PFC numerical parameters / a resolution check (M-F).**

**Strongly recommended:**
6. Replicate seeds (≥5/density) + error bars on α (M-A); rate sweep at fixed ρ_f (M-E); measure ρ_ss at ≥2 rates or soften KM language (M-D); add the missing references (M-G; Devincre 2006/2008, Kubin 2008, Madec 2002/2017, Franciosi–Zaoui 1982, Skaugen 2018, Berry 2006, Salvalaglio–Elder 2022, a STEM-tomography ref, the ExaDiS paper).

**Minor:** resolve all "[DOI to verify]"; reconcile 27/33 line count and "6×" density-span numbers; vendor the `interfaceB_exporter` dependency or document it.

---

## 7. Editorial note

This is a well-written, refreshingly honest manuscript whose claims currently run ahead of its shown evidence — and, in the Taylor analysis, are framed in a way that *understates* the genuinely good result hiding in the data (bulk-band forest-hardening slopes once τ₀ is separated). The fix is not spin; it is leading with the right fit, measuring (now done) rather than asserting the baseline, demoting the STEM claim to feasibility, and scoping the PFC negative result. Do that and this becomes a clean, publishable methods/scoping contribution.

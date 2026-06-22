# Slip-system assignment uncertainty (STEM → DDD)

`stem_to_exadis.py` assigns each reconstructed line a single slip system by `argmin|n·t|` (the {111}
plane that best contains the line tangent). This is **geometry-only**: it uses no experimental g·b
diffraction-contrast Burgers determination, and the z-tangent is weakly constrained by STEM stereo.
Two or three systems are often nearly equally compatible.

**Upgrade (IDR v2 path).** Instead of one forced choice, `defect_ir` keeps, per line, the **top-k
slip-system candidates with normalized priors** (`defect_ir.uncertainty.slip_system_candidates`):
- `plane_containment_score = |n̂·t̂|` (0 = line lies in the glide plane = best),
- `line_character_score    = |b̂·t̂|` (1 = screw, 0 = edge),
- `prior = softmax(quality / T)` over the kept candidates,
- `assignment_confidence` = prior of the chosen system; `assignment_status = geometry_only_pending_gb`.

`network_assignment_summary` then reports mean confidence, per-edge entropy (bits), and the number of
ambiguous edges. On the worked Cu example the lines come out **genuinely ambiguous** (mean confidence
~0.33, ~1.6 bits) — which the old single-assignment silently hid. Downstream, `idr_to_exadis_network`
can take the top-1 OR **sample** an assignment, enabling a Monte-Carlo assignment-sensitivity audit
(planned M3). A real g·b analysis would set `assignment_status = gb_validated` and collapse the priors.

## v1.1 clarification (ambiguity exists; topology amplification is minor)
Two claims must be kept separate:
- **Assignment ambiguity EXISTS** (unchanged, true): without g·b, geometry fixes the {111} plane but the 3
  ⟨110⟩ Burgers on it are near-degenerate — mean confidence ~0.33, ~1.58 bits, 243/243 edges ambiguous
  (`cu_stem_idr_report.md`).
- **Assignment ambiguity does NOT strongly amplify topology** (corrected): the v0 "~5×" junction increase
  was an *edgewise-sampling artifact* — per-edge draws break Burgers continuity within a single physical
  line (~143/216 within-line discontinuities). The physical **`sample_linewise`** policy (one draw per
  `parent_line_id`) gives junction counts at ~top-1 level. So the ambiguity must be propagated **per line**;
  its real impact on topology is minor in this network. See `REAL_NETWORK_AUDIT.md` v1.1.
Resolving it via experimental g·b is still valuable for *correctness* (collapsing the entropy), not for a
large topology swing.

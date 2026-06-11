# Paper Outline (academic-paper skill, outline-only mode)

> Generated 2026-06-11 by the autonomous iteration loop.
> ⚠️ Paper Configuration Record awaits user confirmation (IRON RULE checkpoint) —
> adjust any field below and the outline will be revised accordingly.

## Paper Configuration Record

| Field | Value |
|-------|-------|
| Working title | Phase-field crystal simulation of tensile dislocation mechanisms: from glide to multiplication avalanche at crystal-plasticity-informed densities |
| Paper type | Computational methods + mechanism study (IMRaD, adapted) |
| Discipline | Computational materials science / metal plasticity |
| Target venue tier | Computational Materials Science; Modelling Simul. Mater. Sci. Eng. (alt: Materialia) |
| Citation format | Vancouver (numbered) — CMS house style |
| Output format | LaTeX (.tex + .bib) when drafting begins |
| Language | English main text; bilingual abstract (EN + 中文) |
| Target length | 7,000–8,500 words, 8–10 figures |
| Source corpus | docs/STAGE_REPORT.md + results/ figure set (self-generated data) |
| AI disclosure | Required (simulations, analysis, and draft produced with Claude Code assistance; user-directed) |

## Detailed Outline with Evidence Map

### 1. Introduction (~1,100 words)
- 1.1 Multiscale modeling of metal plasticity: CP-FEM (DAMASK) ↔ DDD ↔ MD chains;
  the mesoscale gap at diffusive timescales. [cite: Roters 2019 DAMASK; LeSar DDD review — verify DOIs]
- 1.2 PFC as the atomic-resolution/diffusive-timescale bridge. [cite: Elder 2002 PRL;
  Elder & Grant 2004 PRE; Berry 2006 PRE; avalanches: Chan 2010 PRL — verify DOIs]
- 1.3 Gap: PFC branches in multiscale frameworks are designed but rarely
  implemented with density-matched coupling to CP state variables.
- 1.4 This work: implements the PFC branch (Layer 4b) of a DAMASK→DDD→MD
  framework; Interface B data contract; mechanism zoo M1–M7 + cascade.
- 1.5 Contributions list (numerics lessons, mechanism quantification,
  density-matched coupling feasibility).

### 2. Model and Methods (~1,800 words)
- 2.1 PFC free energy & conserved dynamics; parameters (r=-0.25, ψ̄=-0.25).
- 2.2 Numerics: rfft spectral semi-implicit + stabilized splitting (C=3maxψ²);
  per-step real-field round trip — include the phantom-mode stripe instability
  as a documented pitfall (Fig. S1: debug_blowup_fields).
  [Evidence: PROGRESS.md iteration 2; killer-probe validation]
- 2.3 Loading: affine box rescaling (k-rescaling), area/volume-conserving
  tension; effective strain rate = 1/RELAX. [cite: Hirouchi 2009 — verify]
- 2.4 Stress measurement: virtual affine deformation dF/dε; cross-validation
  vs energy curvature (0.209 vs 0.211). [Evidence: A3a + stress validation]
- 2.5 Defect detection: sub-grid peak refinement → periodic Delaunay → 5|7
  pairing; phase-winding dislocation seeding (dipole/quadrupole).
  [cite: Skaugen 2018 — verify; Evidence: gates A1/A2]
- 2.6 3D extension: one-mode BCC; gates C1/C2. [Evidence: c1/c2 results]
- 2.7 Scale mapping to DAMASK ROI (Interface B): a₀↔b_Cu; density-matched box
  sizing. [Evidence: interfaceB_bridge report.json]

### 3. Results (~2,600 words)
- 3.1 Verification gate suite (Table 1: A1–A3, C1–C3 with pass metrics).
- 3.2 M1/M2 glide & annihilation (dipole, 9.68% annihilation). [analysis.png]
- 3.3 M4 polycrystal plasticity & rate sensitivity m=0.28 (Fig: analysis_b45).
- 3.4 M5 noise-assisted nucleation phase diagram (threshold band 0.03–0.045).
  [nucleation_phase_diagram.png]
- 3.5 M6/M7 pore mechanics: Gibbs-Thomson dissolution vs persistent-pore
  rim emission ("dislocation factory", 18-dislocation cluster).
  [b4_pore_evolution.png, b4_pore_v2_overview.png]
- 3.6 Flagship: ROI-density-matched cascade — yield σ/E≈7% @11.6%,
  avalanche to 316 cores (ρ=1.75e17 m⁻²), flow at 49% of peak.
  [analysis_d2.png, analysis_cascade.png, cascade.json]
- 3.7 3D polycrystal tension (modulus 75% of single crystal). [analysis_iter7]

### 4. Discussion (~1,400 words)
- 4.1 Mechanism map vs DDD/MD expectations; what PFC adds at diffusive timescales.
- 4.2 Rate exponent m=0.28 ↔ diffusion/GB-dominated creep & superplasticity.
- 4.3 Density-matched coupling: 40 nm boxes reach DAMASK hotspot densities —
  implications for Interface B round-tripping (event metrics → DDD priors).
- 4.4 Limitations: PFC-unit stress (no absolute MPa without elastic-constant
  calibration); diffusive vs inertial timescales; 2D-dominant evidence;
  detector pairing at high density; single-parameter-point (r, ψ̄).
- 4.5 Outlook: r×rate×density matrix (HPC), 3D dislocation-line tracking,
  cyclic loading vs 单晶铜 fatigue protocol, XPFC for FCC Cu.

### 5. Conclusion (~350 words)

### Mandatory back matter
Data availability (GitHub repo + npz archives), AI disclosure, CRediT,
COI, funding, acknowledgments.

## Figure Plan (8 main + 2 SI)

| # | Figure | Source |
|---|--------|--------|
| 1 | Framework schematic: DAMASK→Interface B→PFC (Layer 4b) | to draw |
| 2 | Gate suite collage (crystallization, dipole, elastic validation) | a1/a2/a3 |
| 3 | Mechanism stress-strain + core counts (single/dipole/poly) | analysis.png |
| 4 | Rate sensitivity 3-panel + m fit | analysis_b45.png |
| 5 | Nucleation phase diagram | nucleation_phase_diagram.png |
| 6 | Pore evolution: dissolution vs factory | pore figs |
| 7 | Cascade 3-panel (yield/avalanche/density) | analysis_cascade.png |
| 8 | Core trajectories at ROI-matched density | analysis_d2.png |
| S1 | Phantom-mode stripe instability | debug_blowup_fields.png |
| S2 | 3D BCC gates + polycrystal curve | c1/c2/c3 |

## Key References (all to be DOI-verified in Phase 5a — IRON RULE: no unverified citations)

1. Elder, Katakowski, Haataja, Grant, PRL 88, 245701 (2002) — PFC introduction
2. Elder & Grant, PRE 70, 051605 (2004) — PFC elasticity/plasticity
3. Berry, Grant, Elder, PRE 73, 031609 (2006) — PFC dislocation dynamics
4. Stefanovic, Haataja, Provatas, PRL 96, 225504 (2006) — modified PFC
5. Chan, Tsekenis, Dantzig, Dahmen, Goldenfeld, PRL 105, 015502 (2010) — PFC plasticity avalanches
6. Skaugen, Angheluta, Viñals, PRB 97, 054113 (2018) — dislocation kinematics in PFC
7. Hirouchi, Takaki, Tomita, Comput. Mater. Sci. 44 (2009) — PFC deformation technique
8. Greenwood, Provatas, Rottler, PRL 105, 045702 (2010) — XPFC
9. Roters et al., Comput. Mater. Sci. 158 (2019) — DAMASK
10. Bertin, Aubry, Arsenlis, Cai, MSMSE 27 (2019) / ExaDiS — DDD (verify exact ref)

## Next checkpoints (per skill IRON RULES)
- [ ] User confirms/edits Paper Configuration Record
- [ ] Literature pass with DOI verification (Phase 1 + 5a)
- [ ] Outline approval → Phase 3 argumentation → Phase 4 drafting

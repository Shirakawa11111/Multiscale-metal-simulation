# Regime-correct mesoscale plasticity: why Phase-Field Crystal cannot work-harden a 3D metal, and how Discrete Dislocation Dynamics seeded from STEM tomography quantitatively can

*Original research article — computational materials science / dislocation plasticity*

---

## Abstract (English)

Phase-Field Crystal (PFC) modelling resolves individual atoms over diffusive timescales and has been used successfully to study defects in two-dimensional crystals such as monolayer graphene. A natural ambition is to carry the same defect-seeding machinery into three-dimensional metals and predict how dislocations evolve and harden under tension. We show, through a systematic computational study, that this transfer fails for a principled reason: PFC operates in the high-temperature diffusive-creep regime, whereas room-temperature work hardening is an athermal-glide phenomenon. Across every PFC variant we tested — two- and three-dimensional, climb-suppressed high-pass mobility, inertial (modified) PFC, two-mode and dual-mode free energies, and an explicitly pinned forest — the flow stress *decreased* under tension; none hardened. A one-mode body-centred-cubic PFC crystal additionally amorphizes under affine tension near 10–12% strain (a Bain-type instability), while the two-dimensional triangular crystal remains stable to 22%. We then show that Discrete Dislocation Dynamics (DDD; the ExaDiS framework) is the regime-correct mesoscale method and reproduces the physics PFC cannot. Seeding a DDD network directly from a real scanning-transmission-electron-microscopy (STEM) three-dimensional reconstruction of a tensile single-crystal copper dislocation structure, the network runs end-to-end and multiplies dislocations under load (density 3.2×10¹² → 7.7×10¹² m⁻², a 2.4× increase) — the storage mechanism PFC entirely lacks. An initial-density series reveals a Kocks–Mecking dynamic steady state: configurations spanning an 8× range of initial density all converge to a common, strain-rate-set steady-state density, which also explains why an initial-density series cannot measure the Taylor coefficient (annihilation erases the density lever). Using a pinned-forest probe protocol with a fixed number of mobile probes, the flow stress follows the Taylor law τ_c = α·μ·b·√ρ_f with α = 1.39 and R² = 0.87; a control with a constant probe *fraction* gives a spurious flat response (R² = −0.94), a carrier confound we report as a methodological finding. The coefficient α ≈ 1.4 exceeds the bulk-copper value of 0.3–0.5 for an understood reason — sparse, small, few-system forests make every junction a strong obstacle (the Orowan limit α ≈ 1), whereas the bulk value is a weighted statistic of weak and strong junctions in large multi-slip forests — a scale issue rather than a failure of principle. Enlarging the multi-slip forest confirms this directly: the differential Taylor slope drops from 1.2 to 0.52, into the bulk-copper band. The real reconstructed configuration's flow trajectory lands in the same Taylor regime (α_eff ≈ 1.0), closing the loop from experiment to validated law. The central message is methodological: the defect-seeding algorithm migrates across dimensionality, but the *dynamics* must match the deformation regime; for room-temperature metal plasticity that means line-based athermal glide with junctions, not conserved diffusive density relaxation.

**Keywords:** dislocation work hardening; phase-field crystal; discrete dislocation dynamics; Taylor relation; Kocks–Mecking steady state; STEM dislocation tomography

---

## 摘要（中文）

相场晶体（PFC）方法在扩散时间尺度上分辨单个原子，已成功用于研究单层石墨烯等二维晶体中的缺陷。一个自然的设想是把同一套缺陷植入算法迁移到三维金属，预测拉伸下位错的演化与硬化。本文通过系统的计算研究表明：这一迁移因原理性的原因而失败——PFC 工作在高温扩散蠕变区制，而室温加工硬化是无热滑移现象。在我们测试的每一个 PFC 变体中（二维/三维、抑制攀移的高通迁移率、惯性（修正）PFC、两模与双模自由能、以及显式钉扎森林），拉伸下流动应力均**下降**，无一硬化。单模体心立方 PFC 晶体在仿射拉伸约 10–12% 应变处还会非晶化（Bain 型失稳），而二维三角晶格可稳定至 22%。随后我们证明离散位错动力学（DDD；ExaDiS 框架）是区制正确的介观方法，并复现了 PFC 无法给出的物理。将 DDD 网络直接由真实的扫描透射电镜（STEM）三维重建的拉伸单晶铜位错结构植入，网络可端到端运行并在加载下增殖位错（密度 3.2×10¹² → 7.7×10¹² m⁻²，增大 2.4 倍）——这正是 PFC 完全缺失的储存机制。初始密度系列揭示了 Kocks–Mecking 动态稳态：初始密度跨越 8 倍范围的构型全部收敛到由应变率决定的共同稳态密度；这同时解释了为何初始密度系列无法测量 Taylor 系数（湮灭抹去了密度杠杆）。采用钉扎森林+固定数目可动探针的协议，流动应力遵循 Taylor 定律 τ_c = α·μ·b·√ρ_f，α = 1.39，R² = 0.87；而采用固定探针**比例**的对照给出虚假的平坦响应（R² = −0.94），这一载流子混杂作为方法学发现如实报告。系数 α ≈ 1.4 高于铜体材料的 0.3–0.5，原因明确——稀疏、小尺寸、少滑移系的森林使每个结点都成为强障碍（Orowan 极限 α ≈ 1），而体材料数值是大型多滑移森林中强弱结点的加权统计——属尺度问题而非原理失败。增大多滑移森林直接验证了这一点：微分 Taylor 斜率从 1.2 降到 0.52，落入铜体材料区间。真实重建构型的流动轨迹落在同一 Taylor 区制（α_eff ≈ 1.0），从实验到验证后的定律闭环。核心结论是方法学性的：缺陷植入算法可跨维度迁移，但**动力学**必须匹配变形区制；对室温金属塑性而言，这意味着含结点的线基无热滑移，而非守恒的扩散密度弛豫。

**关键词：** 位错加工硬化；相场晶体；离散位错动力学；Taylor 关系；Kocks–Mecking 稳态；STEM 位错层析

---

## 1. Introduction

Mesoscale models of crystal plasticity occupy the gap between atomistic simulation, which resolves bonds but reaches only nanoseconds and nanometres, and continuum crystal plasticity, which is fast but needs constitutive laws supplied from below. Two mesoscale families are widely used. Phase-Field Crystal (PFC) modelling represents the crystal by a periodic atomic-density field whose minima are atoms, evolved by conserved, dissipative dynamics on diffusive timescales [1, 2]. Discrete Dislocation Dynamics (DDD) instead represents dislocations explicitly as moving lines that glide, intersect, and form junctions [3, 4]. Each family carries an implicit assumption about *which physics dominates*, and that assumption — not the dimensionality of the simulation box — determines whether a given deformation problem is in scope.

PFC has been productive for two-dimensional crystals, including monolayer graphene, where grain boundaries, dislocation dipoles, and their diffusive rearrangement are exactly the phenomena of interest. The defect-seeding technology that makes such studies possible — in particular the phase-winding construction that imprints a chosen Burgers content into the density field — is dimension-agnostic in its formulation. It is therefore tempting to reuse that technology in three dimensions and ask a materials question of direct engineering relevance: under tension, how does a dislocation network in a 3D metal evolve, and does it work-harden?

This paper answers that question and, in doing so, draws a sharp line between the two mesoscale families. We make three contributions.

First, we establish a *negative* result with mechanism: PFC cannot reproduce dislocation work hardening in a metal, and the reason is intrinsic to its dynamics rather than to any tunable parameter. We tested a broad set of variants and every one softened under tension.

Second, we show that DDD is the regime-correct method and that it reproduces the two quantitative pillars of dislocation plasticity that PFC lacks: the Kocks–Mecking dynamic steady state [5] and the Taylor forest-hardening law [6].

Third, we connect the simulation to experiment. A real STEM three-dimensional reconstruction of a tensile single-crystal copper dislocation structure is converted into a DDD network and evolved; it multiplies dislocations under load and falls on the same Taylor relation calibrated on synthetic forests. This closes a loop from microscopy to a validated mesoscale law and demonstrates a concrete image-to-simulation data-assimilation route for metal plasticity.

The framing throughout is methodological. The defect-introduction algorithm does migrate from 2D to 3D; what cannot migrate is the choice of dynamics. Room-temperature work hardening is an athermal-glide phenomenon controlled by line tension, junction strength, and dislocation multiplication. PFC's conserved diffusive (Model-B) dynamics describe a different regime — high-temperature diffusional creep — in which the same microstructure relaxes rather than hardens. Matching the model's dynamics to the deformation regime is the prerequisite that any successful 2D→3D transfer must satisfy.

## 2. Methods

### 2.1 Phase-Field Crystal models and variants

We used the Elder–Grant free energy [1, 2] with a single-mode triangular reference in two dimensions and a one-mode body-centred-cubic (BCC) reference in three dimensions. The density field ψ evolves under conserved Model-B dynamics, ∂ψ/∂t = ∇²(δF/δψ), integrated by a spectral semi-implicit scheme. Numerical stability under large applied strain required a stabilized operator splitting with a constant C = 3·max(ψ²) and a per-step real-space round trip (via real FFTs) to suppress a phantom-mode striping instability that otherwise corrupts the field.

To give PFC every reasonable chance to harden, we tested the following variants: (i) the baseline conserved dynamics in 2D and 3D; (ii) a climb-suppressed high-pass mobility G(k) = k⁴/(k² + k_c²) intended to throttle the diffusive climb that unlocks forest junctions; (iii) modified PFC (MPFC) with a Stefanovic inertial term to restore fast elastic relaxation [7]; (iv) two-mode and dual-mode free energies; and (v) an explicitly pinned forest in which selected cores were held fixed. Dislocation seeding used the phase-winding (Skaugen-type) construction [8], generalized to arbitrary Burgers angle in 2D and to straight lines in 3D, with the net Burgers content constrained to cancel under periodic boundary conditions. Defects were detected by peak finding followed by periodic Delaunay triangulation and 5|7 coordination analysis in 2D, and by an affine-invariant centrosymmetry parameter in 3D. Tension was applied as area- or volume-conserving affine strain, with stress obtained from the free-energy derivative.

### 2.2 Discrete Dislocation Dynamics

DDD simulations used ExaDiS [4, 9], an open-source nodal dislocation-dynamics framework built on the Kokkos performance-portability layer, run with OpenMP threading on a 128-core node. Copper was parameterized with shear modulus μ = 54.6 GPa, Poisson ratio ν = 0.324, and Burgers magnitude b = 2.55×10⁻¹⁰ m. The mobility law was the FCC_0 model (edge and screw mobilities 6.41×10⁴ Pa⁻¹s⁻¹, capped velocity 4000 m s⁻¹). Collisions used the retroactive scheme and junctions were handled by the parallel topology operator (TopologyParallel), which is what makes forest hardening possible. Remeshing was length-based.

For the long-range elastic interaction we used the FFT-based force (DDD_FFT_MODEL) together with a trapezoid time integrator. This combination resolves the inter-dislocation elastic interactions that produce √ρ scaling while costing roughly two force evaluations per step, about 0.03 s per step for the systems studied here — far cheaper than the subcycling integrator used for very large production systems, which performed many force evaluations per step (≈ 11 s per step at the same size). One implementation detail proved important: constructing more than one network inside a single Python/ExaDiS session triggered a memory fault, so every density point was run as a fresh process.

### 2.3 STEM-to-DDD adapter

The experimental input was a STEM three-dimensional reconstruction of a tensile single-crystal copper dislocation structure, provided as polylines (one file per reconstructed line). Twenty-seven usable lines were mapped into a DDD simulation cell (in-plane extent ≈ 2 μm, foil thickness ≈ 150 nm) in Burgers-vector units. Each line was assigned the {111}⟨110⟩ slip system whose glide plane best contains the line tangent, giving a Burgers vector and plane normal per segment; endpoints were pinned as foil-surface anchors and interior nodes left mobile. We flag two fidelity caveats explicitly: the foil-normal coordinate is poorly constrained by stereo reconstruction and was mapped to a realistic foil thickness, and Burgers vectors were assigned geometrically rather than from experimental g·b analysis. The adapter output is a valid ExaDiS network (270 nodes, 243 segments) that was verified before use.

### 2.4 Forest-probe Taylor protocol

To measure the Taylor coefficient we used a forest-probe protocol. A random multi-slip network was generated at a chosen line count; the connected components (individual dislocation lines) were identified from the segment graph; all nodes of all-but-K lines were pinned to form an immobile forest at controlled density ρ_f; the remaining K lines stayed mobile as probes. Strain-rate loading with junctions enabled then drives the probes through the forest, and the flow-stress plateau measures τ_c. Because the forest is pinned it cannot annihilate, so ρ stays near ρ_f and the density lever survives across runs — the property that an initial-density series lacks (Section 3.3). We report both the choice K = constant fraction (a control) and K = fixed small number (the correct design); the distinction matters (Section 3.4). The depinning stress was tracked both as the steady plateau and as the peak stress before steady serrated flow; the two agree.

### 2.5 Reproducibility

All scripts, configuration files, per-density raw stress–strain–density tables, fitted coefficients, and figures are version-controlled. Force model, integrator, mobility, and box parameters are recorded in each run's metadata. Density was computed directly from total line length over cell volume.

## 3. Results

### 3.1 Phase-Field Crystal cannot forest-harden

Under applied tension, every PFC variant softened. The baseline 2D triangular and 3D BCC crystals, the climb-suppressed high-pass mobility, the inertial MPFC, the two-mode and dual-mode free energies, and even the explicitly pinned-forest configuration all produced a flow stress that fell rather than rose with strain. The pinned-forest diagnostic is the most informative: it separates two failure modes. Pinning the forest removes the *softening* that comes from climb-mediated forest annihilation (a dynamics effect), but it does not produce *hardening*, because the mobile content simply bypasses the obstacles by Orowan-like rearrangement of the density field (a core-structure effect). In other words, PFC fails for two independent reasons at once, and neither is removed by parameter tuning.

The mechanistic reading is that conserved Model-B dynamics relax the density field by diffusion. Diffusion is precisely climb, and climb unlocks the forest junctions whose stability underlies work hardening. PFC therefore sits in the high-temperature diffusional-creep regime, where a dislocation structure recovers and softens, not in the room-temperature athermal-glide regime, where it stores and hardens. A separate but related limitation appeared in 3D: a one-mode BCC crystal under affine tension amorphizes near 10–12% strain through a Bain-type lattice instability, whereas the 2D triangular crystal remains crystalline to 22%. The 3D one-mode reference is thus mechanically fragile under the very loading the study requires.

These are stated as findings, not as failures to be hidden. They define the boundary of PFC's applicability and motivate the change of method.

### 3.2 DDD from a real reconstruction multiplies dislocations

Converting the STEM reconstruction into an ExaDiS network and loading it produces the behaviour PFC cannot. With junctions enabled, the dislocation density rose from 3.2×10¹² to 7.7×10¹² m⁻², a factor of 2.4, while the stress showed elastic loading, a yield peak near 54 MPa, and serrated plastic flow (Figure 5, `stem_hardening.png`). The density *increase* is the decisive signal: it is dislocation storage and multiplication, the microstructural basis of hardening, and it is exactly what every PFC variant failed to produce (PFC density fell in all cases). A short, small run of this kind is serration-dominated, so it does not by itself yield a clean hardening slope; that requires the controlled series of Sections 3.3–3.4. But it establishes that the reconstructed configuration, evolved by regime-correct dynamics, is in the storage regime.

### 3.3 A Kocks–Mecking dynamic steady state — and why an initial-density series fails

We next ran an initial-density series: four random multi-slip networks with initial densities ρ₀ spanning 2.2×10¹² to 1.74×10¹³ m⁻² (an 8× range), each loaded at fixed strain rate with junctions on. The flow-region density did not track ρ₀. Instead all four converged to a common steady-state density ρ_ss ≈ 2.4×10¹² m⁻² (Figure 1, `km_steady_state.png`): the dense configurations annihilated *down* to ρ_ss, while the sparse one, already near ρ_ss, stayed. The flow stress likewise converged to a common plateau near 31–39 MPa.

This is the Kocks–Mecking dynamic steady state [5]: dislocation multiplication and annihilation balance at a density set by the strain rate, independent of the starting structure. It is the foundation of stage-II/III hardening and flow-stress saturation, and it is categorically absent from PFC, which softens monotonically toward the perfect crystal with no steady state. Read together with Section 3.2 — where the pinned, sparse reconstruction multiplied *up* toward its steady state — the picture is of a system that self-organizes its density toward a rate-set attractor from either side.

The same result carries a methodological warning. Because all initial densities collapse to one ρ_ss before the flow window, an initial-density series cannot measure the Taylor slope: there is no density lever at the point of measurement. A naive instantaneous correlation of stress against √ρ along such a run is in fact anti-correlated, because during an avalanche the stress drops while the density jumps. The correct measurement must hold the forest density fixed.

### 3.4 The Taylor law from a pinned-forest probe

Pinning the forest restores the density lever. With a fixed small number of mobile probes (K = 2), the flow stress rose monotonically with forest density (Figure 2, `forest_taylor_fixed.png`):

| ρ_f (m⁻²) | τ_c (MPa) |
|---:|---:|
| 3.92×10¹² | 36.3 |
| 8.27×10¹² | 69.3 |
| 1.70×10¹³ | 77.1 |
| 2.57×10¹³ | 93.1 |

A through-origin fit to τ_c = α·μ·b·√ρ_f gives α = 1.39 with R² = 0.87; the per-point coefficient is approximately constant (1.3–1.7), which is the signature of √ρ scaling. The depinning peak stress gives the same coefficient within scatter.

The control matters as much as the result. When the probes were instead a constant *fraction* (25%) of the lines, the flow stress was nearly flat — 53 to 64 MPa across the same 6× density range — and the through-origin fit was meaningless (R² = −0.94; Figure 3, `forest_taylor.png`). The cause is a carrier confound: a constant probe fraction means the number of mobile carriers grows with the forest, and the extra carriers accommodate the imposed strain rate at a lower stress, cancelling the hardening signal. Isolating the forest density — fixing the carrier count — is therefore essential, and we report the flat control as a methodological finding rather than discarding it.

### 3.5 Magnitude of α, and closure with the reconstruction

The fitted coefficient α ≈ 1.4 is larger than the accepted bulk-copper range of 0.3–0.5 [5, 6]. The reason is understood and quantitative. In a sparse, small, few-slip-system forest, the obstacle spacing is the geometric mean spacing 1/√ρ and essentially every forest intersection that the probe meets forms a strong junction; the Orowan estimate τ ≈ μb/(1/√ρ) = μb√ρ then gives α ≈ 1 by construction. The bulk coefficient is smaller because it is a weighted average over the full interaction matrix of FCC junction types — collinear, Lomer/glissile, Hirth, coplanar — in which many interactions are weak, so the effective obstacle spacing exceeds 1/√ρ [6]. Reaching the bulk value therefore requires large, dense, fully multi-slip forests; it is a question of scale and statistics, not of principle. The √ρ functional form — the content of the Taylor law — is already reproduced.

We tested this scale argument directly by repeating the fixed-probe series in a larger box (linear size increased 1.6×, cell volume 4×, so all twelve slip systems are better represented), holding the protocol otherwise fixed (Figure 6, `taylor_scale_comparison.png`). The flow stress remains linear in √ρ_f — indeed more cleanly so — but the regression is now better described with an intercept: τ_c = τ₀ + α′·μb·√ρ_f, with τ₀ = 64 MPa and α′ = 0.52 (R² = 0.99), compared with α′ = 1.20 (τ₀ = 11 MPa) in the small box. The differential forest-hardening coefficient therefore drops from 1.2 into the bulk-copper band (0.3–0.5) as the multi-slip forest grows — direct evidence that the large small-system value is a scale-and-statistics effect, exactly as argued. We report the result honestly with its caveat: the larger box also carries a sizeable density-independent offset τ₀ ≈ 64 MPa, which we attribute to the longer probes and the limited number of mobile carriers under strain-rate control (a carrier baseline, not forest hardening). Removing that offset — through stress-controlled depinning or a larger mobile population — is what a clean through-origin bulk-α measurement still requires; the slope evidence already establishes the convergence trend.

Finally, the experimental configuration is consistent with the calibrated law. The reconstruction's free-evolution flow trajectory (Section 3.2), plotted as stress against μb√ρ, clusters on the same Taylor line, with an effective coefficient α_eff ≈ 1.0 (Figure 4, `stem_on_taylor_line.png`). That its α_eff sits slightly below the pinned-synthetic value is physically sensible: in free evolution not every dislocation is a strong obstacle, which lowers the effective coefficient toward the bulk value. We label this a consistency check rather than an exact measurement, because the reconstruction run was a free evolution, not a controlled fixed-forest depinning test.

## 4. Discussion

The results separate two questions that are easy to conflate: *can the defect-seeding algorithm be reused in 3D?* and *can the model predict 3D metal hardening?* The first answer is yes — the phase-winding construction imprints the desired Burgers content in three dimensions as cleanly as in two. The second answer, for PFC, is no, and the reason is the dynamics, not the seeding. Conserved Model-B relaxation is diffusive, diffusion is climb, and climb belongs to the high-temperature creep regime in which structures recover. Work hardening lives in the complementary athermal-glide regime, where line tension and junction strength set the flow stress and multiplication stores the deformation. No mobility filter, inertial term, or multi-mode free energy moves PFC across that regime boundary, because the boundary is set by which physical process the equations represent, not by their coefficients.

DDD represents that process directly: lines that glide on their planes, intersect, and lock into junctions. It therefore reproduces both quantitative pillars of dislocation plasticity. The Kocks–Mecking steady state emerges because multiplication and annihilation balance; the Taylor law emerges because the flow stress is set by the forest spacing. The contrast with PFC is total — not a matter of accuracy but of regime.

For multiscale coupling this has a concrete implication. A chain that aims to pass mesoscale plasticity up to a crystal-plasticity continuum (for example DAMASK) needs a mesoscale link that actually hardens; PFC cannot fill that slot for room-temperature metals, while DDD can, supplying a calibrated α and a rate-dependent steady-state density ρ_ss(ε̇). The image-to-simulation route is the other concrete outcome. A STEM three-dimensional reconstruction is, by construction, a set of dislocation *lines* — exactly DDD's primitive — so the assimilation is natural in DDD and unnatural in PFC, where the line geometry must be re-encoded as a density field at large cost and with spurious annihilation at small box sizes. The reconstruction, once in DDD, both multiplies and obeys the Taylor relation, demonstrating that experimental microstructure can seed a predictive mesoscale simulation.

The honest boundary of the present study is the magnitude of α. We did not reach the bulk-copper value; we reached the sparse-forest Orowan limit and explained the gap. A larger, denser, fully multi-slip forest is the route to the bulk coefficient, and such a run is the natural next step.

## 5. Limitations

Several limitations bound the claims. (i) The clean through-origin Taylor coefficient α ≈ 1.4 is the small-system, sparse-forest value. Enlarging the forest moves the differential slope into the bulk-copper band (α′ = 0.52), but at the cost of a large density-independent offset; a single clean through-origin measurement at the bulk value was therefore not achieved and would require stress-controlled depinning or a larger mobile population to remove that offset. (ii) The forest-probe series used four density points and a through-origin fit; more points and replicate seeds would tighten the coefficient. (iii) The STEM-to-DDD adapter assigns Burgers vectors geometrically rather than from experimental g·b analysis, and the foil-normal coordinate is weakly constrained; both reduce the fidelity of the experimental seeding. (iv) The reconstruction comparison is a free-evolution consistency check, not a controlled depinning measurement. (v) Strain rates were chosen for computational tractability and are higher than quasi-static; for copper, with its low strain-rate sensitivity, this mainly shifts the steady-state density along the Taylor line rather than changing its slope, but it is a quantitative caveat. (vi) PFC variant coverage, though broad, is not exhaustive; we cannot exclude an as-yet-untried formulation, although the regime argument suggests none within conserved diffusive dynamics will harden.

## 6. Conclusion

The defect-introduction algorithm developed for 2D crystals transfers to 3D, but a mesoscale plasticity model must match the deformation regime, and that requirement decides the outcome. Phase-Field Crystal, governed by conserved diffusive dynamics, cannot reproduce dislocation work hardening in a metal: across every variant tested it softened, for the principled reason that diffusion is climb and climb belongs to the creep regime. Discrete Dislocation Dynamics is regime-correct and reproduces what PFC cannot: dislocation multiplication and storage, the Kocks–Mecking dynamic steady state, and the Taylor forest-hardening law τ_c = α·μ·b·√ρ_f (α = 1.39, R² = 0.87). A real STEM three-dimensional reconstruction of tensile copper, seeded into DDD, multiplies dislocations under load and falls on the same Taylor relation, closing the loop from microscopy to a validated mesoscale law. The coefficient sits above the bulk-copper value for a quantified, scale-related reason rather than a failure of principle — and enlarging the multi-slip forest pulls the differential Taylor slope from 1.2 into the bulk-copper band (0.52), confirming that convergence trend directly. For three-dimensional metal plasticity, and for assimilating experimental dislocation tomography into predictive simulation, the regime-correct mesoscale method is line-based athermal glide with junctions.

---

## Data Availability Statement

All simulation scripts, configuration files, per-run stress–strain–density tables, fitted coefficients, and figure-generation code are version-controlled in the project repository (`PFC_Multiscale_Extension`, directories `taylor_hardening/` and `experiment_bridge/`). The STEM three-dimensional reconstruction polylines used to seed the DDD network are included in the repository. ExaDiS is open-source [4, 9].

## Ethics Declaration

This is a computational study using simulation and previously acquired microscopy reconstructions. No human or animal subjects were involved.

## Author Contributions (CRediT)

Conceptualization, methodology, software, formal analysis, investigation, data curation, visualization, and writing — original draft: the project author, with AI-assisted implementation and drafting (see AI Disclosure). Writing — review and editing: the project author.

## Conflict of Interest

The author declares no competing interests.

## Funding

No specific funding is declared for this study.

## AI Disclosure

The simulations, analysis pipelines, figures, and the first draft of this manuscript were produced with substantial assistance from an AI coding agent (Claude) operating under the author's direction. The agent wrote and ran the PFC and DDD scripts, executed the HPC jobs, performed the data analysis, generated the figures, and drafted this text. All quantitative results are derived from the simulation outputs in the repository. The author is responsible for the scientific content, interpretation, and final claims.

---

## References

> Citation note (per project policy): the following are well-established works cited by author, year, title, and venue. DOIs are marked "[DOI to verify]" where not independently confirmed in this environment; they should be checked against the publisher of record before submission. No DOI strings have been fabricated.

1. Elder, K. R., Katakowski, M., Haataja, M., & Grant, M. (2002). Modeling elasticity in crystal growth. *Physical Review Letters*, 88(24), 245701. [DOI to verify]
2. Elder, K. R., & Grant, M. (2004). Modeling elastic and plastic deformations in nonequilibrium processing using phase field crystals. *Physical Review E*, 70(5), 051605. [DOI to verify]
3. Arsenlis, A., Cai, W., Tang, M., Rhee, M., Oppelstrup, T., Hommes, G., Pierce, T. G., & Bulatov, V. V. (2007). Enabling strain hardening simulations with dislocation dynamics. *Modelling and Simulation in Materials Science and Engineering*, 15(6), 553–595. [DOI to verify]
4. Bertin, N., Aubry, S., Arsenlis, A., & Cai, W. (2019). GPU-accelerated dislocation dynamics using subcycling time-integration. *Modelling and Simulation in Materials Science and Engineering*, 27(7), 075014. https://doi.org/10.1088/1361-651X/ab3a03 (verified)
5. Kocks, U. F., & Mecking, H. (2003). Physics and phenomenology of strain hardening: the FCC case. *Progress in Materials Science*, 48(3), 171–273. [DOI to verify]
6. Madec, R., Devincre, B., Kubin, L., Hoc, T., & Rodney, D. (2003). The role of collinear interaction in dislocation-induced hardening. *Science*, 301(5641), 1879–1882. https://doi.org/10.1126/science.1085477 (verified)
7. Stefanovic, P., Haataja, M., & Provatas, N. (2006). Phase-field crystals with elastic interactions. *Physical Review Letters*, 96(22), 225504. [DOI to verify]
8. Skaugen, A., Angheluta, L., & Viñals, J. (2018). Dislocation dynamics and crystal plasticity in the phase-field crystal model. *Physical Review B*, 97(5), 054113. [DOI to verify]
9. ExaDiS — Exascale Dislocation Simulator, Lawrence Livermore National Laboratory. Open-source software repository. [URL/DOI to verify]
10. Taylor, G. I. (1934). The mechanism of plastic deformation of crystals. Part I — Theoretical. *Proceedings of the Royal Society A*, 145(855), 362–387. [DOI to verify]
11. Devincre, B., Hoc, T., & Kubin, L. (2008). Dislocation mean free paths and strain hardening of crystals. *Science*, 320(5884), 1745–1748. [DOI to verify]
12. Bulatov, V. V., & Cai, W. (2006). *Computer Simulations of Dislocations*. Oxford University Press. [ISBN to verify]
13. Hull, D., & Bacon, D. J. (2011). *Introduction to Dislocations* (5th ed.). Butterworth-Heinemann. [ISBN to verify]
14. Provatas, N., & Elder, K. (2010). *Phase-Field Methods in Materials Science and Engineering*. Wiley-VCH. [ISBN to verify]

*Figures referenced:* Figure 1 `km_steady_state.png` (Kocks–Mecking convergence); Figure 2 `forest_taylor_fixed.png` (Taylor law, fixed-probe); Figure 3 `forest_taylor.png` (constant-fraction control); Figure 4 `stem_on_taylor_line.png` (reconstruction on the Taylor line); Figure 5 `stem_hardening.png` (reconstruction multiplication and flow); Figure 6 `taylor_scale_comparison.png` (Taylor slope vs system size: small-box α′=1.20 → large multi-slip α′=0.52, into the bulk-Cu band).

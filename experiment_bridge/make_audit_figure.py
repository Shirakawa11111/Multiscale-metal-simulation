"""Four-panel main-result figure for the STEM-to-DDD v2 audit package.

Every QUANTITATIVE label is read from the result JSONs at runtime -> numbers cannot drift from the data.
(Panel A counts/confidence/entropy from cu_stem_idr_report.json + assignment_sensitivity.json; B-D from their
own JSONs. The box LAYOUT in A is a fixed schematic, and the downstream-stability descriptor is qualitative.)
  A  pipeline schematic: STEM reconstruction -> IDR (candidates+uncertainty) -> line-coherent ExaDiS lowering
  B  assignment entropy collapse under g.b reflections           (synthetic_gb.json)
  C  edgewise vs line-coherent within-line Burgers discontinuities (assignment_sensitivity.json)
  D  density convention: Lambda_A is z-invariant, rho_vol ~ 1/z_eff (density_conventions.json)

  python3 experiment_bridge/make_audit_figure.py  -> results_exadis/audit_summary_figure.png
"""
import os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = os.path.join(os.path.dirname(__file__), "results_exadis")
gb = json.load(open(os.path.join(OUT, "synthetic_gb.json")))
asg = json.load(open(os.path.join(OUT, "assignment_sensitivity.json")))
den = json.load(open(os.path.join(OUT, "density_conventions.json")))
idr = json.load(open(os.path.join(OUT, "cu_stem_idr_report.json")))

fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
fig.suptitle("STEM-to-DDD v2 — uncertainty-aware, self-auditing pipeline (IDR v1.1)",
             fontsize=13, fontweight="bold")

# ---- Panel A: pipeline schematic ----
a = ax[0, 0]
a.set_xlim(0, 10); a.set_ylim(0, 10); a.axis("off")
a.set_title("A  Auditable path: microscope → DDD initial condition", fontsize=11, loc="left")
asgm = idr["assignment"]
disc_lw = int(round(asg["policies"]["sample_linewise"]["within_line_discontinuities_mean"]))
boxes = [
    (1.0, 6.3, f"STEM\nreconstruction\n{idr['n_lines']} lines / {idr['n_vertices']} v / {idr['n_edges']} e", "#cfe8ff"),
    (1.0, 2.0, f"IDR\ncandidates + uncertainty\n{asgm['mean_confidence']:.3g} conf · {asgm['mean_assignment_entropy_bits']:.3g} bits", "#d8f0d8"),
    (6.2, 2.0, f"line-coherent\nlowering (sample_linewise)\nwithin-line discont = {disc_lw}", "#ffe6cc"),
    (6.2, 6.3, "ExaDiS / DDD\nstable IC\nrobust topology + density", "#f0d8f0"),
]
for x, y, t, c in boxes:
    a.add_patch(FancyBboxPatch((x, y), 2.8, 2.6, boxstyle="round,pad=0.1",
                               fc=c, ec="#333", lw=1.2))
    a.text(x + 1.4, y + 1.3, t, ha="center", va="center", fontsize=8.3)
arr = dict(arrowstyle="-|>", color="#333", lw=1.6, mutation_scale=16)
a.add_patch(FancyArrowPatch((2.4, 6.3), (2.4, 4.6), **arr))
a.add_patch(FancyArrowPatch((3.8, 3.3), (6.2, 3.3), **arr))
a.add_patch(FancyArrowPatch((7.6, 4.6), (7.6, 6.3), **arr))
a.text(5.0, 3.55, "g·b-ready", ha="center", fontsize=7.5, style="italic", color="#555")

# ---- Panel B: entropy collapse ----
b = ax[0, 1]
sc = gb["scenarios"]
order = ["no_gb", "1_reflection_g200", "2_reflections_g200_g020", "3_reflections_g200_g020_g002"]
labels = ["no g·b", "1 reflection", "2 reflections", "3 reflections"]
ents = [sc[k]["mean_entropy_bits"] for k in order]
res = [sc[k]["frac_resolved"] * 100 for k in order]
bars = b.bar(labels, ents, color="#4a90d9", width=0.6)
b.axhline(1.585, ls="--", lw=1, color="#999"); b.text(2.6, 1.62, "log₂3 = 1.585", fontsize=8, color="#666")
for bar, e, r in zip(bars, ents, res):
    b.text(bar.get_x() + bar.get_width() / 2, e + 0.03, f"{e:.2f} b\n{r:.0f}% resolved",
           ha="center", va="bottom", fontsize=8)
b.set_ylabel("mean assignment entropy (bits)"); b.set_ylim(0, 1.9)
b.set_title("B  g·b collapses assignment ambiguity (line-coherent)", fontsize=11, loc="left")
b.tick_params(axis="x", labelsize=9)

# ---- Panel C: edgewise vs linewise ----
c = ax[1, 0]
pol = asg["policies"]
names = ["top1", "sample_edgewise\n(deprecated)", "sample_linewise\n(default)"]
disc = [pol["top1"]["within_line_discontinuities_mean"],
        pol["sample_edgewise"]["within_line_discontinuities_mean"],
        pol["sample_linewise"]["within_line_discontinuities_mean"]]
tot = pol["sample_edgewise"]["intra_line_adjacencies"]
bars = c.bar(names, disc, color=["#888", "#d9534f", "#5cb85c"], width=0.6)
for bar, d in zip(bars, disc):
    c.text(bar.get_x() + bar.get_width() / 2, d + 2, f"{d}", ha="center", va="bottom", fontsize=9)
c.axhline(tot, ls=":", lw=1, color="#999"); c.text(0.0, tot - 12, f"{tot} intra-line adjacencies", fontsize=8, color="#666")
c.set_ylabel("within-line Burgers discontinuities")
c.set_title("C  Edgewise sampling is an artifact; line-coherent is clean", fontsize=11, loc="left")
c.set_ylim(0, tot + 20); c.tick_params(axis="x", labelsize=8.5)

# ---- Panel D: density convention ----
d = ax[1, 1]
zmap = {"foil_z600": 600, "thickened_z3_1800": 1800, "thickened_z5_3000": 3000, "thickened_z10_6000": 6000}
ab = den["as_built"]["volume_density_m2_by_convention"]
zs = [zmap[k] for k in zmap]; rho_ab = [ab[k] for k in zmap]
d.loglog(zs, rho_ab, "o-", color="#d9534f", label="ρ_vol as-built (∝ 1/z)")
rlx = den.get("relaxed_from_cell_policy_audit") or {}
if rlx:
    rho_r = [rlx[k]["volume_density_m2"] for k in zmap if k in rlx]
    zr = [zmap[k] for k in zmap if k in rlx]
    d.loglog(zr, rho_r, "s--", color="#f0ad4e", label="ρ_vol relaxed (∝ 1/z)")
la = den["as_built"]["lambda_area_m1_foil_native"]
d.set_xlabel("declared effective thickness z_eff (b)")
d.set_ylabel("ρ_vol (m⁻²)")
d.set_title("D  ρ_vol is a 1/z convention; Λ_A is the invariant", fontsize=11, loc="left")
d.text(0.5, 0.06, f"Λ_A (foil-native, z-invariant) = {la:.2e} m⁻¹",
       transform=d.transAxes, fontsize=8.5, color="#2a6", fontweight="bold")
d.legend(fontsize=8, loc="upper right"); d.grid(True, which="both", ls=":", alpha=0.4)

fig.tight_layout(rect=[0, 0, 1, 0.96])
path = os.path.join(OUT, "audit_summary_figure.png")
fig.savefig(path, dpi=150)
print(f"-> {os.path.relpath(path, os.path.dirname(__file__))}")

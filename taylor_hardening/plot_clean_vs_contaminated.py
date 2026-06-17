"""Clean vs contaminated forest-probe Taylor coefficient.

Contaminated (K=2 probes, drag-inflated baseline): alpha = 0.83 +/- 0.42 (3 seeds,
one config never converged). Clean (K=8 carriers, large box, measured baseline):
alpha = 0.735 +/- 0.063 (7 seeds). The diagnostic-driven fix (more carriers,
removing carrier starvation) collapses the seed scatter ~7x and pins the
coefficient firmly in the strong-obstacle regime, ABOVE bulk Cu (0.3-0.5).

  python3 plot_clean_vs_contaminated.py
"""
import os, glob, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MU, B = 54.6e9, 2.55e-10


def alpha_of(root, seed, baseTag, fTag, dens):
    def load(t):
        fp = glob.glob(os.path.join(root, t, "*", "result.json"))
        return json.load(open(fp[0])) if fp else None
    base = load(baseTag(seed))
    if not base:
        return None
    tau0 = base["flow_stress"]
    pts = [load(fTag(seed, nl)) for nl in dens]
    pts = [p for p in pts if p]
    if len(pts) < 3:
        return None
    rho = np.array([p["rho_forest"] for p in pts])
    tau = np.array([p["flow_stress"] for p in pts])
    x = MU * B * np.sqrt(rho)
    y = tau - tau0
    return float(np.sum(x * y) / np.sum(x * x))


# clean (K=8): /tmp/clean2, tags B_s{seed}, F_n{nl}_s{seed}
clean = []
for s in ["1234", "5678", "2222", "3333", "4444", "5555", "6666", "7777"]:
    a = alpha_of("/tmp/clean2", s, lambda s: f"B_s{s}",
                 lambda s, nl: f"F_n{nl}_s{s}", [50, 100, 200])
    if a is not None:
        clean.append(a)

# contaminated (K=2, drag-inflated 72.1 MPa baseline): per-seed baseline-
# subtracted alphas from the committed multi-seed run (RESULTS.md step 4).
cont = [0.366, 0.742, 1.373]   # seeds 1234, 5678, 9012

clean = np.array(clean); cont = np.array(cont)
fig, ax = plt.subplots(figsize=(8, 5.6))
ax.axhspan(0.3, 0.5, color="0.6", alpha=0.25, label="Cu bulk alpha 0.3-0.5")
for i, a in enumerate(cont):
    ax.plot(0 + 0.04 * (i - 1), a, "o", ms=9, color="tab:red")
for i, a in enumerate(clean):
    ax.plot(1 + 0.03 * (i - 3.5), a, "o", ms=9, color="tab:blue")
if len(cont):
    ax.errorbar(0, cont.mean(), yerr=cont.std(), fmt="_", ms=40, color="tab:red",
                capsize=8, lw=2, label=f"contaminated K=2: {cont.mean():.2f}+/-{cont.std():.2f}")
if len(clean):
    ax.errorbar(1, clean.mean(), yerr=clean.std(), fmt="_", ms=40, color="tab:blue",
                capsize=8, lw=2, label=f"clean K=8: {clean.mean():.2f}+/-{clean.std():.2f}")
ax.set_xticks([0, 1]); ax.set_xticklabels(["contaminated\n(K=2, drag baseline)", "clean\n(K=8, measured baseline)"])
ax.set_ylabel(r"forest-hardening coefficient $\alpha$")
ax.set_title("Removing carrier starvation collapses the scatter ~7x\n"
             r"and pins $\alpha\approx0.74$ — above bulk Cu (strong-obstacle regime)")
ax.set_xlim(-0.5, 1.5); ax.set_ylim(0, 1.6)
ax.legend(fontsize=9, loc="upper right"); ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fp = os.path.dirname(os.path.abspath(__file__)) + "/clean_vs_contaminated.png"
fig.savefig(fp, dpi=130, bbox_inches="tight")
print(f"clean alpha = {clean.mean():.3f} +/- {clean.std():.3f} (n={len(clean)})")
print(f"contaminated alpha = {cont.mean():.3f} +/- {cont.std():.3f} (n={len(cont)})")
print("saved", fp)

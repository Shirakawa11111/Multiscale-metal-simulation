"""4-panel density-lever figure for the collinear-dominance bounded negative.
Usage: python3 make_dlev_figure.py <xslip_on_dir> <xslip_off_dir> [out.png]
"""
import json, glob, re, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ON = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dlev/dlev"
OFF = sys.argv[2] if len(sys.argv) > 2 else "/tmp/dlev2/dlev2"
OUT = sys.argv[3] if len(sys.argv) > 3 else "density_lever.png"


def load(root):
    D = {}
    for f in glob.glob(root + "/*/flow.json"):
        D[re.search(r"/([^/]+)/flow.json", f).group(1)] = json.load(open(f))
    return D


def rss(D, t):
    return D[t]["tau_flow_MPa"] * abs(D[t]["schmid_primary"])


def seeds(D, lv, P):
    return [k for k in D if re.match(rf"{lv}_{P}_s\d+$", k)]


def mean_over_seeds(D, lv, P, fn):
    ks = seeds(D, lv, P)
    return float(np.mean([fn(D, k) for k in ks])) if ks else np.nan


on, off = load(ON), load(OFF)
fig, ax = plt.subplots(2, 2, figsize=(13, 10))

# Panel A: R_RSS coll/glissile vs measured rho_f
for D, lab, c, mk in [(on, "XSLIP on", "tab:blue", "o"), (off, "XSLIP off", "tab:red", "s")]:
    xs, ys = [], []
    for lv in ("lo", "hi"):
        rho = mean_over_seeds(D, lv, "glissile", lambda d, k: d[k]["rho_forest_settled"])
        co = mean_over_seeds(D, lv, "coll_opp", rss); gl = mean_over_seeds(D, lv, "glissile", rss)
        xs.append(rho); ys.append(co / gl)
    ax[0, 0].plot(xs, ys, mk + "-", color=c, ms=11, label=lab)
ax[0, 0].axhline(2.3, ls="--", color="green", lw=2, label="canonical ~2.3 (Madec)")
ax[0, 0].axhline(1.0, ls=":", color="gray")
ax[0, 0].set_xscale("log"); ax[0, 0].set_ylim(0, 2.6)
ax[0, 0].set_xlabel("measured forest density $\\rho_f$ (m$^{-2}$)")
ax[0, 0].set_ylabel("$R_{RSS}$ = $\\tau_{RSS}^{coll}/\\tau_{RSS}^{gliss}$"); ax[0, 0].legend(loc="upper left")
ax[0, 0].set_title("(A) Collinear vs glissile (RSS): ~1 and FLAT\n(dominance would be ~2.3 & growing)")

# Panel B: R_RSS opp/same
for D, lab, c, mk in [(on, "XSLIP on", "tab:blue", "o"), (off, "XSLIP off", "tab:red", "s")]:
    xs, ys = [], []
    for lv in ("lo", "hi"):
        rho = mean_over_seeds(D, lv, "glissile", lambda d, k: d[k]["rho_forest_settled"])
        co = mean_over_seeds(D, lv, "coll_opp", rss); cs = mean_over_seeds(D, lv, "coll_same", rss)
        xs.append(rho); ys.append(co / cs)
    ax[0, 1].plot(xs, ys, mk + "-", color=c, ms=11, label=lab)
ax[0, 1].axhline(1.0, ls=":", color="gray")
ax[0, 1].set_xscale("log"); ax[0, 1].set_ylim(0.8, 1.4)
ax[0, 1].set_xlabel("measured forest density $\\rho_f$ (m$^{-2}$)")
ax[0, 1].set_ylabel("$\\tau_{RSS}^{opp}/\\tau_{RSS}^{same}$")
ax[0, 1].set_title("(B) Annihilation single-bit toggle: ~1\n(opposite vs same sense -> no effect)")
ax[0, 1].legend()

# Panel C: collinear forest drift (the gate limitation) vs density
for D, lab, c, mk in [(on, "XSLIP on", "tab:blue", "o"), (off, "XSLIP off", "tab:red", "s")]:
    xs, ys = [], []
    for lv in ("lo", "hi"):
        rho = mean_over_seeds(D, lv, "glissile", lambda d, k: d[k]["rho_forest_settled"])
        dr = float(np.max([abs(D[k]["forest_drift"]) for k in seeds(D, lv, "coll_opp")]))
        xs.append(rho); ys.append(dr * 100)
    ax[1, 0].plot(xs, ys, mk + "-", color=c, ms=11, label=lab)
ax[1, 0].axhline(5, ls="--", color="orange", lw=2, label="5% gate")
ax[1, 0].set_xscale("log")
ax[1, 0].set_xlabel("measured forest density $\\rho_f$ (m$^{-2}$)")
ax[1, 0].set_ylabel("collinear partner forest drift (%)")
ax[1, 0].set_title("(C) Co-driven collinear partner is NOT a stable forest\n(drift-limited -> strict gate AMBIGUOUS)")
ax[1, 0].legend()

# Panel D: protocol ladder (qualitative)
ax[1, 1].axis("off")
steps = [
    ("free-probe (glide-through)", "carrier extinction -> underestimate", "fail"),
    ("FR-source probe", "source bow-out baseline dominates", "fail"),
    ("local binary reaction", "mechanism + length scaling REAL;\nlocal strength NOT dominant (mid-pack)", "ok-neg"),
    ("pairwise evolving-forest MFP", "sampling-starved; geometry can't host", "retract"),
    ("multi-slip density-lever (this)", "R_RSS coll/gliss ~1, flat, XSLIP-indep", "ok-neg"),
]
y = 0.92
for name, note, kind in steps:
    col = {"fail": "tab:red", "retract": "tab:gray", "ok-neg": "tab:green"}[kind]
    ax[1, 1].text(0.02, y, "▶ " + name, fontsize=11, fontweight="bold", color=col, transform=ax[1, 1].transAxes)
    ax[1, 1].text(0.07, y - 0.058, note, fontsize=9, color="black", transform=ax[1, 1].transAxes)
    y -= 0.185
ax[1, 1].set_title("(D) Observable ladder: every probe -> NO collinear dominance", loc="center")

fig.suptitle("Bounded negative for canonical collinear dominance in the present ExaDiS/FCC_0 protocol",
             fontweight="bold", fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT, dpi=130, bbox_inches="tight")
print("wrote", OUT)

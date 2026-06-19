"""Capstone: does the first-principles matrix PREDICT the STEM network's
hardening consistently with direct measurement? Compare:
  - alpha predicted from (measured a_ij) x (STEM junction inventory)
  - alpha measured directly (clean forest-probe, rate-extrapolated)
  - uniform-population bulk macro alpha
  - STEM free-evolution effective alpha
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
p = json.load(open(os.path.join(HERE, "stem_prediction.json")))

labels = ["uniform-population\nbulk macro",
          "predicted for STEM\n(matrix x inventory)",
          "direct clean\nforest-probe",
          "STEM free-evolution\n(effective)"]
vals = [p["macro_uniform"], p["alpha_resampled_mean"], 0.85, 1.0]
errs = [0.0, p["alpha_resampled_std"], 0.11, 0.0]   # 0.85 = midpoint of 0.74-0.96
colors = ["0.6", "tab:red", "tab:blue", "tab:purple"]

fig, ax = plt.subplots(figsize=(8.5, 5.2))
x = np.arange(len(labels))
ax.bar(x, vals, 0.6, yerr=errs, capsize=5, color=colors)
ax.axhspan(0.3, 0.5, color="tab:green", alpha=0.15, label="Cu bulk 0.3-0.5")
for i, v in enumerate(vals):
    ax.text(i, v + 0.03, f"{v:.2f}", ha="center", fontsize=11)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel(r"forest-hardening coefficient $\alpha$")
ax.set_title("Experiment ↔ first-principles matrix: the prediction loop closes\n"
             "(STEM inventory is strong-junction-skewed → predicts α above bulk, "
             "consistent with direct measurement)")
ax.set_ylim(0, 1.2); ax.legend(fontsize=9, loc="upper left"); ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fp = os.path.join(HERE, "stem_prediction.png")
fig.savefig(fp, dpi=130, bbox_inches="tight")
print(f"predicted alpha_network = {p['alpha_resampled_mean']:.3f} +/- {p['alpha_resampled_std']:.3f}")
print("saved", fp)

"""Plot the first-principles FCC interaction matrix and the recovered macroscopic
Taylor coefficient against the bulk-Cu band.

  python3 plot_matrix.py /tmp/matrix
"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = sys.argv[1] if len(sys.argv) > 1 else "/tmp/matrix"
d = json.load(open(os.path.join(R, "matrix_result.json")))
amat = d["a_matrix"]; macro = d["alpha_macro"]
order = ["coplanar", "collinear", "glissile", "Hirth", "self", "Lomer"]
order = [t for t in order if t in amat]
alpha = [np.sqrt(amat[t]) for t in order]
# Madec/Devincre-type literature interaction coefficients a_ij (approx, for ref)
lit_a = {"self": 0.09, "coplanar": 0.05, "collinear": 0.57, "Hirth": 0.06,
         "glissile": 0.10, "Lomer": 0.11}
lit_alpha = [np.sqrt(lit_a[t]) for t in order]

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
x = np.arange(len(order))
ax[0].bar(x - 0.2, alpha, 0.4, color="tab:blue", label="DDD (this work)")
ax[0].bar(x + 0.2, lit_alpha, 0.4, color="0.6", label="literature ~a_ij (ref)")
ax[0].set_xticks(x); ax[0].set_xticklabels(order, rotation=30, ha="right")
ax[0].set_ylabel(r"$\alpha_{ij}=\sqrt{a_{ij}}$")
ax[0].set_title("First-principles FCC interaction matrix per junction type")
ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3, axis="y")

ax[1].axhspan(0.3, 0.5, color="tab:green", alpha=0.2, label="Cu bulk 0.3-0.5")
ax[1].bar([0], [macro], 0.5, color="tab:red")
ax[1].text(0, macro + 0.02, f"{macro:.3f}", ha="center", fontsize=12)
ax[1].set_xticks([0]); ax[1].set_xticklabels(["population-weighted\nmacroscopic α"])
ax[1].set_ylabel(r"$\alpha_{macro}=\sqrt{\langle a_{ij}\rangle}$")
ax[1].set_ylim(0, 0.7)
ax[1].set_title("Recovered bulk Taylor coefficient\n(from the measured matrix)")
ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3, axis="y")
fig.tight_layout()
fp = os.path.dirname(os.path.abspath(__file__)) + "/interaction_matrix.png"
fig.savefig(fp, dpi=130, bbox_inches="tight")
print(f"macro alpha = {macro:.3f}")
print("saved", fp)

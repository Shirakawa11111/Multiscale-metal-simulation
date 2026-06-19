"""Fit the directly-measured STEM-inventory forest alpha (Phase 2b validation).
Per seed, fit resolved tau vs sqrt(rho_f) across densities (slope = alpha,
baseline removed); average over seeds. Compare to the matrix prediction (0.69)
and the uniform-population macro (0.43).

  python3 fit_mixed.py /tmp/mixed
"""
import os, sys, glob, json
import numpy as np

R = sys.argv[1] if len(sys.argv) > 1 else "/tmp/mixed"
MU, B = 54.6e9, 2.55e-10
rows = [json.load(open(f)) for f in glob.glob(os.path.join(R, "*", "mixed_result.json"))]
print(f"loaded {len(rows)} runs")
byseed = {}
for r in rows:
    byseed.setdefault(r["seed"], []).append(r)

alphas = []
for s, rs in byseed.items():
    if len(rs) < 2:
        continue
    x = MU * B * np.sqrt(np.array([r["rho_f"] for r in rs]))
    y = np.array([r["tau_resolved"] for r in rs])
    A = np.vstack([x, np.ones_like(x)]).T
    (slope, icpt), *_ = np.linalg.lstsq(A, y, rcond=None)
    alphas.append(float(slope))
alphas = np.array(alphas)
print(f"per-seed alpha (slope, baseline-removed): {[round(a,3) for a in alphas]}")
if len(alphas):
    m, sd = float(alphas.mean()), float(alphas.std())
    print(f"\n=== directly-measured STEM-inventory forest alpha = {m:.3f} +/- {sd:.3f} "
          f"(n={len(alphas)} seeds) ===")
    print(f"    matrix prediction (matrix x inventory) = 0.69")
    print(f"    uniform-population macro              = 0.43")
    verdict = "CONSISTENT with prediction" if abs(m - 0.69) < 2 * (sd + 0.05) else "DEVIATES from prediction"
    print(f"    -> {verdict}")
    json.dump(dict(alpha_direct=m, alpha_direct_std=sd, n_seeds=len(alphas),
                   predicted=0.69, uniform=0.43),
              open(os.path.join(R, "mixed_fit.json"), "w"), indent=1)

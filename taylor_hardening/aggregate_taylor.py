"""Aggregate per-density results (one process per density) into the Taylor fit.

Reads OUT/d*/result.json (each written by a fresh taylor_series.py process),
fits sigma_flow = alpha * mu * b * sqrt(rho) through the origin, and writes the
combined taylor_series.json + taylor_curves.json that plot_taylor.py consumes.

  python3 aggregate_taylor.py taylor_out
"""
import os, sys, glob, json
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else "taylor_out"
MU, B_CU = 54.6e9, 2.55e-10

results = []
for fp in sorted(glob.glob(os.path.join(OUT, "d*", "result.json"))):
    results.append(json.load(open(fp)))
results.sort(key=lambda r: r["rho_flow"])

series, allres = [], {}
for r in results:
    curve = r.pop("curve", None)
    series.append(r)
    if curve is not None:
        allres[str(r["num_lines"])] = curve

rr = np.array([s["rho_flow"] for s in series])
ss = np.array([s["flow_stress"] for s in series])
x = MU * B_CU * np.sqrt(rr)
alpha = float(np.sum(x * ss) / np.sum(x * x)) if len(ss) else 0.0
if len(ss) > 1:
    pred = alpha * x
    ss_res = float(np.sum((ss - pred) ** 2))
    ss_tot = float(np.sum((ss - ss.mean()) ** 2))
    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else 0.0
else:
    r2 = 0.0

out = dict(series=series, taylor_alpha=alpha, r2=r2, mu=MU, b=B_CU,
           n_densities=len(series),
           note="Taylor flow-stress density series (one process per density), "
                "FCC Cu, ExaDiS DDD with junctions. sigma_flow vs mu*b*sqrt(rho).")
with open(os.path.join(OUT, "taylor_series.json"), "w") as f:
    json.dump(out, f, indent=1)
with open(os.path.join(OUT, "taylor_curves.json"), "w") as f:
    json.dump(allres, f)

print(f"n_densities={len(series)}")
for s in series:
    print(f"  rho_flow={s['rho_flow']:.3e}  flow={s['flow_stress']/1e6:.1f} MPa")
print(f"=== TAYLOR alpha={alpha:.3f}, R^2={r2:.3f} (Cu bulk lit ~0.3-0.5; "
      f"sparse small-box DDD ~1) ===")

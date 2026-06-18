"""Fit the FCC interaction matrix from the large-scale campaign and recover the
macroscopic Taylor coefficient.

For each (mobile, forest) pair, fit resolved tau_c vs sqrt(rho_f) across forest
densities (slope = alpha_mf, intercept = probe self-baseline). Average alpha_mf
(and a_mf = alpha_mf^2) per junction type. Then the macroscopic forest-hardening
coefficient for an equal-population multi-slip forest is
        alpha_macro = sqrt(<a_ij>)  (population-weighted over junction types),
which for FCC Cu should land near the literature 0.3-0.5 if the strong junctions
(collinear/Lomer) are correctly down-weighted by the many weak ones.

  python3 fit_matrix.py /tmp/matrix
"""
import os, sys, glob, json
import numpy as np

R = sys.argv[1] if len(sys.argv) > 1 else "/tmp/matrix"
MU, B = 54.6e9, 2.55e-10
# FCC ordered-pair multiplicities (population weights for the macro average)
MULT = {"self": 12, "coplanar": 24, "collinear": 12, "Hirth": 24,
        "glissile": 48, "Lomer": 24}

rows = []
for fp in glob.glob(os.path.join(R, "*", "pair_result.json")):
    try:
        rows.append(json.load(open(fp)))
    except Exception:
        pass
print(f"loaded {len(rows)} pair runs")

# group by (m,f); fit resolved tau vs sqrt(rho_f) across densities+seeds
bypair = {}
for r in rows:
    bypair.setdefault((r["msys"], r["fsys"]), []).append(r)

pair_alpha = {}
for (m, f), rs in bypair.items():
    x = MU * B * np.sqrt(np.array([r["rho_f"] for r in rs]))
    y = np.array([r["tau_resolved"] for r in rs])
    if len(rs) >= 2 and np.ptp(x) > 0:
        A = np.vstack([x, np.ones_like(x)]).T
        (slope, icpt), *_ = np.linalg.lstsq(A, y, rcond=None)
        a = float(slope)
    else:
        a = float(np.mean(y / x))           # through-origin fallback
    pair_alpha[(m, f)] = (a, rs[0]["jtype"], len(rs))

# average per junction type
bytype = {}
for (m, f), (a, t, n) in pair_alpha.items():
    bytype.setdefault(t, []).append(a)

print(f"\n{'junction':>10} {'alpha_ij':>9} {'a_ij':>8} {'npairs':>7} {'mult':>5}")
amat = {}
for t in ["self", "coplanar", "Hirth", "glissile", "Lomer", "collinear"]:
    if t in bytype:
        al = np.array(bytype[t]); a = float(np.mean(al))
        amat[t] = a
        print(f"{t:>10} {a:>9.3f} {a**2:>8.3f} {len(al):>7} {MULT[t]:>5}")

# macroscopic alpha = sqrt( sum(mult_t * a_t^2) / sum(mult_t) ) over measured types
if amat:
    num = sum(MULT[t] * amat[t] ** 2 for t in amat)
    den = sum(MULT[t] for t in amat)
    alpha_macro = float(np.sqrt(num / den))
    print(f"\n=== population-weighted macroscopic alpha = {alpha_macro:.3f} "
          f"(bulk Cu literature 0.3-0.5) ===")
    out = dict(a_matrix=amat, alpha_macro=alpha_macro,
               n_pairs=len(pair_alpha), n_runs=len(rows),
               note="FCC interaction matrix from controlled single-system forest "
                    "depinning; alpha_macro = sqrt(pop-weighted <a_ij>).")
    json.dump(out, open(os.path.join(R, "matrix_result.json"), "w"), indent=1)
    print("saved", os.path.join(R, "matrix_result.json"))

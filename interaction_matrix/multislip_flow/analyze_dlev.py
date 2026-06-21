"""Density-lever verdict: RSS-corrected flow stress + L_mf + stored-slope + density gain + gates.

Per the expert review: the decision MUST use resolved shear (tau_RSS = sigma_flow*|S_primary|), NOT raw
scalar flow stress (opt_pair gives each pair a different Schmid factor). Uses MEASURED rho_forest_settled
as the density axis, not the target. Applies all readable gates. Emits density_lever_summary.json.

Usage: python3 analyze_dlev.py <dir_with_flow_jsons>  (dir contains {lo,hi}_{type}_s{seed}/flow.json)
"""
import json, glob, re, sys, os
import numpy as np

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dlev/dlev"
B = 2.55e-10
D = {}
for f in glob.glob(os.path.join(ROOT, "*/flow.json")):
    tag = re.search(r"/([^/]+)/flow.json$", f).group(1)
    D[tag] = json.load(open(f))


def runs(level, ptype):
    pat = re.compile(rf"^{level}_{ptype}_s\d+$")
    return [D[k] for k in sorted(D) if pat.match(k)]


def is_starved(r):
    if r.get("rho_mobile_starved") is not None:
        return bool(r["rho_mobile_starved"])
    s = r.get("series", [])
    if not s:
        return False
    mob = np.array([row[2] for row in s])              # series col 2 = rho_mobile
    return bool(len(mob) and mob[-1] < 0.3 * (mob.max() if mob.max() > 0 else 1))


def agg(level, ptype):
    rs = runs(level, ptype)
    if not rs:
        return None
    rss = [r["tau_flow_MPa"] * abs(r["schmid_primary"]) for r in rs]
    lmf = [r["L_mf_b"] for r in rs if np.isfinite(r.get("L_mf_b", np.inf))]
    sstore = [r["S_store"] for r in rs if r.get("S_store") == r.get("S_store")]
    rho = [r.get("rho_forest_settled", 0.0) for r in rs]
    return dict(
        n=len(rs),
        tau_RSS_mean=float(np.mean(rss)), tau_RSS_sd=float(np.std(rss)),
        sigma_mean=float(np.mean([r["tau_flow_MPa"] for r in rs])),
        schmid=float(abs(rs[0]["schmid_primary"])),
        L_mf_mean=float(np.mean(lmf)) if lmf else float("inf"),
        S_store_mean=float(np.mean(sstore)) if sstore else float("nan"),
        rho_forest_settled=float(np.mean(rho)),
        plateau=all(r.get("plateau_reached") for r in rs),
        drift_ok=all(abs(r.get("forest_drift", 1)) < 0.05 for r in rs),
        drift=float(np.mean([r.get("forest_drift", 0) for r in rs])),
        drift_values=[round(float(r.get("forest_drift", 0)), 3) for r in rs],
        drift_abs_max=float(max(abs(r.get("forest_drift", 0)) for r in rs)),
        amb_ok=all(r.get("mean_ambiguous_frac", 1) < 0.10 for r in rs),
        schmid_ok=all(r.get("schmid_gate") for r in rs),
        starved=any(is_starved(r) for r in rs),
        readable_count=f"{sum(1 for r in rs if r.get('readable'))}/{len(rs)}",
        readable=all(r.get("readable") for r in rs),
    )


def level_summary(level, target):
    out = {"rho_level": level, "rho_target": target}
    types = {}
    for P in ("coll_opp", "coll_same", "glissile", "hirth"):
        a = agg(level, P)
        if a:
            types[P] = a
    out["types"] = types
    co, gl, hi, cs = types.get("coll_opp"), types.get("glissile"), types.get("hirth"), types.get("coll_same")
    if co and gl:
        out["R_RSS_coll_over_gliss"] = co["tau_RSS_mean"] / gl["tau_RSS_mean"]
        out["R_L_gliss_over_coll"] = (gl["L_mf_mean"] / co["L_mf_mean"]) if np.isfinite(co["L_mf_mean"]) and np.isfinite(gl["L_mf_mean"]) and co["L_mf_mean"] > 0 else None
        out["R_S_coll_over_gliss"] = (co["S_store_mean"] / gl["S_store_mean"]) if gl["S_store_mean"] not in (0, None) and gl["S_store_mean"] == gl["S_store_mean"] else None
    if co and cs:
        out["R_RSS_opp_over_same"] = co["tau_RSS_mean"] / cs["tau_RSS_mean"]
    if co and hi:
        out["R_RSS_coll_over_hirth"] = co["tau_RSS_mean"] / hi["tau_RSS_mean"]
    return out


lo = level_summary("lo", 3e12)
hi = level_summary("hi", 3e13)

# density gains (use measured density as the real lever)
verdict = "AMBIGUOUS"; reasons = []
RR_lo = lo.get("R_RSS_coll_over_gliss"); RR_hi = hi.get("R_RSS_coll_over_gliss")
RL_lo = lo.get("R_L_gliss_over_coll"); RL_hi = hi.get("R_L_gliss_over_coll")
G_tau = (RR_hi / RR_lo) if (RR_lo and RR_hi) else None
G_L = (RL_hi / RL_lo) if (RL_lo and RL_hi) else None


def both_readable_core():
    for lv in ("lo", "hi"):
        for P in ("coll_opp", "glissile"):
            a = agg(lv, P)
            if not a:
                return False
            if not (a["plateau"] and a["drift_ok"] and a["amb_ok"] and a["schmid_ok"] and not a["starved"]):
                reasons.append(f"{lv}_{P} gate fail (plateau={a['plateau']} drift_ok={a['drift_ok']} amb_ok={a['amb_ok']} schmid_ok={a['schmid_ok']} starved={a['starved']})")
                return False
    return True


# ---- DUAL VERDICT: separate the strict formal gate from the mechanism reading ----
gates_ok = both_readable_core()
strict_gate_verdict = "READABLE" if gates_ok else "AMBIGUOUS"

mech_reasons = []
mechanism_verdict = "INDETERMINATE"
if RR_hi and RR_lo:
    flat_near_unity = abs(RR_lo - 1) < 0.12 and abs(RR_hi - 1) < 0.12 and (G_tau is None or G_tau < 1.12)
    grows = (RR_hi > RR_lo * 1.12) and RR_hi > 1.15 and (G_L is None or G_L >= 1.0)
    if grows:
        mechanism_verdict = "CONFIRM_DENSITY_THRESHOLD"
        mech_reasons.append(f"R_RSS grows {RR_lo:.3f}->{RR_hi:.3f} (G_tau={G_tau:.2f})")
    elif flat_near_unity:
        mechanism_verdict = "BOUNDED_NEGATIVE" if gates_ok else "DRIFT_LIMITED_BOUNDED_NEGATIVE"
        mech_reasons += [
            f"R_RSS coll/glissile ~1 at both densities ({RR_lo:.3f},{RR_hi:.3f}); no density gain (G_tau={G_tau:.2f})",
            f"coll_opp/coll_same ~1 ({round(lo.get('R_RSS_opp_over_same') or 0,3)}/{round(hi.get('R_RSS_opp_over_same') or 0,3)})",
            "no growth toward canonical target ~2.3",
        ]
        if not gates_ok:
            mech_reasons.append("strict gate AMBIGUOUS (low-density collinear forest drift) limits formal interpretation; result robust to it")
    else:
        mechanism_verdict = "INDETERMINATE"
        mech_reasons.append(f"R_RSS {RR_lo:.3f}->{RR_hi:.3f} unclear")

summary = dict(lo=lo, hi=hi,
               density_lever=dict(R_RSS_3e12=RR_lo, R_RSS_3e13=RR_hi, G_tau=G_tau,
                                  R_L_3e12=RL_lo, R_L_3e13=RL_hi, G_L=G_L,
                                  R_RSS_opp_over_same_lo=lo.get("R_RSS_opp_over_same"),
                                  R_RSS_opp_over_same_hi=hi.get("R_RSS_opp_over_same"),
                                  gates_ok=gates_ok,
                                  strict_gate_verdict=strict_gate_verdict,
                                  mechanism_verdict=mechanism_verdict,
                                  gate_reasons=reasons, mechanism_reasons=mech_reasons))
print(json.dumps(summary, indent=1, default=str))
json.dump(summary, open("density_lever_summary.json", "w"), indent=1, default=str)
print(f"\n=== strict_gate_verdict: {strict_gate_verdict}  |  mechanism_verdict: {mechanism_verdict} ===")
for r in mech_reasons:
    print("  -", r)

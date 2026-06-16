"""Analyze the DDD diagnostic campaign:
  Batch A (rate sweep): is the forest flow stress (flow - baseline) rate-stable
                        (physical forest hardening) or rising with rate (drag)?
  Batch B1 (carrier sweep, no forest): does tau0 fall as mobile probes K rise?
  Batch B2 (fixed total, vary K): does flow stress drop with more carriers?

  python3 analyze_campaign.py /tmp/camp
"""
import os, sys, json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = sys.argv[1] if len(sys.argv) > 1 else "/tmp/camp"


def load(tag):
    fp = glob.glob(os.path.join(R, tag, "*", "result.json"))
    return json.load(open(fp[0])) if fp else None


def erate_val(s):
    return float(s)


# ---- Batch A: rate sweep ----
rates = ["1e4", "2e4", "4e4", "7e4", "1.2e5"]
seeds = ["1234", "5678"]
print("=== Batch A: rate sweep (fixed forest NL=100) ===")
print(f"{'erate':>8} {'seed':>5} {'flowF':>7} {'base':>6} {'excess':>7} (MPa)")
A = {}
for e in rates:
    exc = []
    for s in seeds:
        f = load(f"rF_e{e}_s{s}"); b = load(f"rB_e{e}_s{s}")
        if f and b:
            fe = (f["flow_stress"] - b["flow_stress"]) / 1e6
            exc.append(fe)
            print(f"{e:>8} {s:>5} {f['flow_stress']/1e6:7.1f} {b['flow_stress']/1e6:6.1f} {fe:7.1f}")
    if exc:
        A[erate_val(e)] = (np.mean(exc), np.std(exc))
print()

# ---- Batch B1: carrier sweep tau0 vs K ----
print("=== Batch B1: carrier sweep, no forest, tau0 vs K ===")
B1 = {}
for k in [2, 4, 8, 16, 32]:
    r = load(f"cB_k{k}")
    if r:
        B1[k] = r["flow_stress"] / 1e6
        print(f"  K={k:>2}: tau0 = {r['flow_stress']/1e6:.1f} MPa")
print()

# ---- Batch B2: fixed total, vary K ----
print("=== Batch B2: fixed total NL=100, vary mobile K ===")
B2 = {}
for k in [2, 5, 10, 20]:
    r = load(f"cF_k{k}")
    if r:
        B2[k] = r["flow_stress"] / 1e6
        print(f"  K={k:>2}: flow = {r['flow_stress']/1e6:.1f} MPa (k_probe={r.get('k_probe')})")

# ---- figure ----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
if A:
    es = sorted(A); m = [A[e][0] for e in es]; sd = [A[e][1] for e in es]
    ax[0].errorbar(es, m, yerr=sd, fmt="o-", capsize=4, color="tab:blue")
    ax[0].set_xscale("log"); ax[0].set_xlabel("strain rate (1/s)")
    ax[0].set_ylabel("forest excess  flow - baseline (MPa)")
    ax[0].set_title("Rate sweep: is forest stress rate-stable?\n(flat = physical, rising = drag)")
    ax[0].grid(alpha=0.3); ax[0].set_ylim(0, None)
if B1:
    ks = sorted(B1)
    ax[1].plot(ks, [B1[k] for k in ks], "s-", color="tab:red")
    ax[1].set_xscale("log", base=2); ax[1].set_xlabel("mobile probes K (no forest)")
    ax[1].set_ylabel(r"baseline $\tau_0$ (MPa)")
    ax[1].set_title("Carrier sweep: does drag baseline fall\nas carriers increase?")
    ax[1].grid(alpha=0.3); ax[1].set_ylim(0, None)
if B2:
    ks = sorted(B2)
    ax[2].plot(ks, [B2[k] for k in ks], "^-", color="tab:green")
    ax[2].set_xlabel("mobile probes K (fixed total=100)")
    ax[2].set_ylabel("flow stress (MPa)")
    ax[2].set_title("Fixed total, more carriers:\nflow stress drops?")
    ax[2].grid(alpha=0.3); ax[2].set_ylim(0, None)
fig.tight_layout()
fp = os.path.dirname(os.path.abspath(__file__)) + "/campaign_diagnostics.png"
fig.savefig(fp, dpi=130, bbox_inches="tight")
# verdict
if A:
    es = sorted(A); vals = [A[e][0] for e in es]
    drift = (vals[-1] - vals[0]) / vals[0] * 100 if vals[0] else 0
    print(f"\nRATE VERDICT: forest excess from {vals[0]:.1f} (e={es[0]:.0e}) to "
          f"{vals[-1]:.1f} MPa (e={es[-1]:.0e}) = {drift:+.0f}% over the rate range")
    print("  -> small drift => forest hardening is rate-stable (physical); "
          "large rise => drag-contaminated")
print("saved", fp)

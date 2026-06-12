"""Merged v1+v2 matrix analysis: m(r) over 7 solid r values and the
multiplication-threshold curve with ensemble statistics (9 r x 3 seeds).
Outputs: results_hpc/ANALYSIS/matrix_v2.png, matrix_v2_metrics.json
"""

import os, json
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..", "results_hpc")
OUT = os.path.join(ROOT, "ANALYSIS")
RS = (-0.35, -0.30, -0.25, -0.23, -0.21, -0.20)   # solid rows only
RELAXES = (100, 400, 1600)
SEEDS = (7, 11, 13)
EPS_CMP = 0.14
QUAD_RS = (-0.34, -0.32, -0.30, -0.28, -0.26, -0.24, -0.22, -0.20, -0.18)


def load(name):
    with open(os.path.join(ROOT, name, "summary.json")) as f:
        s = json.load(f)
    rows = s["rows"]
    return (np.array([x["exx"] for x in rows]),
            np.array([x["sigma"] for x in rows]),
            np.array([x["cores"] for x in rows]))


def main():
    os.makedirs(OUT, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2))
    metrics = {}

    # ---- m(r) over 7 r values ----
    rates = 1.0 / np.array(RELAXES)
    ms, ms_lo, ms_hi = [], [], []
    for r in RS:
        flows = []
        for rx in RELAXES:
            vals = []
            for s in SEEDS:
                try:
                    eps, sig, cores = load(f"poly_r{r}_x{rx}_s{s}")
                    vals.append(sig[int(np.argmin(np.abs(eps - EPS_CMP)))])
                except FileNotFoundError:
                    pass
            flows.append(vals)
        mean_f = np.array([np.mean(v) for v in flows])
        if (mean_f > 0).all():
            m = np.polyfit(np.log(rates), np.log(mean_f), 1)[0]
            # seed-bootstrap spread of m
            boots = []
            for _ in range(200):
                rng = np.random.default_rng(_)
                f_b = [np.mean(rng.choice(v, len(v))) for v in flows]
                if all(x > 0 for x in f_b):
                    boots.append(np.polyfit(np.log(rates), np.log(f_b), 1)[0])
            ms.append(m)
            ms_lo.append(np.percentile(boots, 16))
            ms_hi.append(np.percentile(boots, 84))
        else:
            ms.append(np.nan)
            ms_lo.append(np.nan)
            ms_hi.append(np.nan)
    metrics["m_vs_r"] = dict(r=list(RS), m=ms, lo=ms_lo, hi=ms_hi)
    ax = axes[0]
    ax.errorbar(RS, ms, yerr=[np.array(ms) - ms_lo, np.array(ms_hi) - ms],
                fmt="ko-", capsize=4)
    ax.axhline(0.5, color="r", ls=":", label="classic superplasticity m=0.5")
    ax.axvline(-0.1875, color="b", ls=":", label="liquid spinodal r=-3ψ̄²")
    ax.set_xlabel("r (quench depth)")
    ax.set_ylabel("rate exponent m")
    ax.set_title("rate sensitivity vs quench depth (7 r, bootstrap 1σ)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # ---- multiplication threshold with ensemble stats ----
    on_mean, on_sd, sc_mean, sc_sd = [], [], [], []
    for r in QUAD_RS:
        ons, scs = [], []
        for s in SEEDS:
            try:
                eps, sig, cores = load(f"quadf_r{r}_s{s}")
            except FileNotFoundError:
                continue
            if cores.max() > cores[0]:
                i_on = int(np.argmax(cores > cores[0]))
                ons.append(eps[i_on] * 100)
                scs.append(sig[i_on])   # critical multiplication STRESS
        on_mean.append(np.mean(ons) if ons else np.nan)
        on_sd.append(np.std(ons) if ons else np.nan)
        sc_mean.append(np.mean(scs) if scs else np.nan)
        sc_sd.append(np.std(scs) if scs else np.nan)
    # linear fit of the critical stress and its zero crossing
    cf = np.polyfit(QUAD_RS, sc_mean, 1)
    r_zero = -cf[1] / cf[0]
    metrics["quad_threshold"] = dict(r=list(QUAD_RS), onset_mean=on_mean,
                                     onset_sd=on_sd,
                                     sigma_crit_mean=sc_mean,
                                     sigma_crit_sd=sc_sd,
                                     sigma_crit_linear_fit=list(cf),
                                     sigma_crit_zero_at_r=float(r_zero))
    ax = axes[1]
    ax.errorbar(QUAD_RS, on_mean, yerr=on_sd, fmt="ko-", capsize=4)
    ax.set_xlabel("r")
    ax.set_ylabel("multiplication onset strain (%)")
    ax.set_title("onset(r), 3-realization ensemble")
    ax.grid(alpha=0.3)
    ax = axes[2]
    ax.errorbar(QUAD_RS, sc_mean, yerr=sc_sd, fmt="rs-", capsize=4,
                label="sigma_crit (at onset)")
    rr = np.linspace(min(QUAD_RS), r_zero, 50)
    ax.plot(rr, np.polyval(cf, rr), "k--", alpha=0.6,
            label=f"linear fit, ->0 at r={r_zero:.3f}")
    ax.set_xlabel("r")
    ax.set_ylabel("critical multiplication stress")
    ax.set_title("sigma_crit(r): monotone; strain-step is a modulus artifact")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.savefig(os.path.join(OUT, "matrix_v2.png"), dpi=140,
                bbox_inches="tight")
    with open(os.path.join(OUT, "matrix_v2_metrics.json"), "w") as f:
        json.dump(metrics, f, indent=1)
    print(json.dumps(metrics, indent=1))


if __name__ == "__main__":
    main()

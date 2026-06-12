"""Full matrix analysis (57 cells, results_hpc/):

  poly (45): flow stress & dislocation retention maps over (r, rate),
             seed-averaged; rate exponent m(r)
  quad  (4): multiplication onset strain / peak stress vs r
  cyc   (8): per-cycle recovery rate & dissipation vs amplitude

Outputs: results_hpc/ANALYSIS/matrix_maps.png, matrix_metrics.json
"""

import os, json, glob
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..", "results_hpc")
OUT = os.path.join(ROOT, "ANALYSIS")
RS = (-0.35, -0.30, -0.25, -0.20, -0.15)
RELAXES = (100, 400, 1600)
SEEDS = (7, 11, 13)
EPS_CMP = 0.14


def load(name):
    with open(os.path.join(ROOT, name, "summary.json")) as f:
        s = json.load(f)
    rows = s["rows"]
    return (np.array([x["exx"] for x in rows]),
            np.array([x["sigma"] for x in rows]),
            np.array([x["cores"] for x in rows]), s)


def main():
    os.makedirs(OUT, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(19, 11))
    metrics = {}

    # ---- poly: flow stress / retention maps ----
    flow = np.full((len(RS), len(RELAXES)), np.nan)
    flow_sd = np.full_like(flow, np.nan)
    keep = np.full_like(flow, np.nan)
    for i, r in enumerate(RS):
        for j, rx in enumerate(RELAXES):
            vals, kept = [], []
            for s in SEEDS:
                try:
                    eps, sig, cores, _ = load(f"poly_r{r}_x{rx}_s{s}")
                    k = int(np.argmin(np.abs(eps - EPS_CMP)))
                    vals.append(sig[k])
                    kept.append(cores[k] / max(cores[0], 1))
                except FileNotFoundError:
                    pass
            if vals:
                flow[i, j] = np.mean(vals)
                flow_sd[i, j] = np.std(vals)
                keep[i, j] = np.mean(kept)
    metrics["poly_flow_stress_at_14pct"] = flow.tolist()
    metrics["poly_flow_stress_sd"] = flow_sd.tolist()
    metrics["poly_core_retention"] = keep.tolist()

    ax = axes[0, 0]
    im = ax.imshow(flow, origin="lower", aspect="auto", cmap="viridis")
    ax.set_xticks(range(3), [f"1/{x}" for x in RELAXES])
    ax.set_yticks(range(5), [str(r) for r in RS])
    ax.set_xlabel("strain rate ~ 1/RELAX")
    ax.set_ylabel("r (quench depth)")
    ax.set_title(f"flow stress @ {EPS_CMP*100:.0f}% (seed-avg)")
    plt.colorbar(im, ax=ax)
    for i in range(5):
        for j in range(3):
            if np.isfinite(flow[i, j]):
                ax.text(j, i, f"{flow[i,j]*1e3:.1f}", ha="center",
                        va="center", color="w", fontsize=8)

    # rate exponent m(r)
    ax = axes[0, 1]
    ms = []
    rates = 1.0 / np.array(RELAXES)
    for i, r in enumerate(RS):
        if np.isfinite(flow[i]).all() and (flow[i] > 0).all():
            m = np.polyfit(np.log(rates), np.log(flow[i]), 1)[0]
        else:
            m = np.nan
        ms.append(m)
        ax.plot(rates, flow[i], "o-", label=f"r={r} (m={m:.2f})")
    metrics["rate_exponent_m_vs_r"] = dict(zip(map(str, RS), ms))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("strain rate")
    ax.set_ylabel("flow stress")
    ax.set_title("rate sensitivity by quench depth")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[0, 2]
    im = ax.imshow(keep, origin="lower", aspect="auto", cmap="magma")
    ax.set_xticks(range(3), [f"1/{x}" for x in RELAXES])
    ax.set_yticks(range(5), [str(r) for r in RS])
    ax.set_title("dislocation retention cores(14%)/cores(0)")
    plt.colorbar(im, ax=ax)

    # ---- quad: multiplication threshold vs r ----
    ax = axes[1, 0]
    qr, onset, speak = [], [], []
    for r in (-0.35, -0.30, -0.25, -0.20):
        try:
            eps, sig, cores, _ = load(f"quad_r{r}")
        except FileNotFoundError:
            continue
        qr.append(r)
        i_on = int(np.argmax(cores > cores[0]))
        onset.append(eps[i_on] * 100 if cores.max() > cores[0] else np.nan)
        speak.append(sig.max())
    metrics["quad_multiplication"] = dict(r=qr, onset_pct=onset,
                                          sigma_peak=speak)
    ax.plot(qr, onset, "ko-", label="multiplication onset %")
    ax2 = ax.twinx()
    ax2.plot(qr, speak, "rs--", label="peak stress")
    ax2.set_ylabel("peak stress", color="r")
    ax.set_xlabel("r")
    ax.set_ylabel("onset strain (%)")
    ax.set_title("multiplication threshold vs quench depth")
    ax.grid(alpha=0.3)

    # ---- cyc: amplitude scan ----
    ax = axes[1, 1]
    axd = axes[1, 2]
    amps = (0.005, 0.01, 0.02, 0.03)
    rec_rate, diss_late = [], []
    for amp in amps:
        rr, dd = [], []
        for s in (7, 11):
            try:
                with open(os.path.join(ROOT, f"cyc_a{amp}_s{s}",
                                       "summary.json")) as f:
                    sm = json.load(f)
            except FileNotFoundError:
                continue
            rows = sm["rows"]
            cyc = np.array([x.get("cycle", 0) for x in rows])
            cores = np.array([x["cores"] for x in rows])
            eps = np.array([x["exx"] for x in rows])
            sig = np.array([x["sigma"] for x in rows])
            ncyc = int(cyc.max()) + 1
            end = [cores[cyc == c][-1] for c in range(ncyc)]
            if end[0] > 0:
                rr.append(100 * (end[0] - end[-1]) / end[0] / ncyc)
            dd.append(abs(np.trapz(sig[cyc == ncyc - 1],
                                   eps[cyc == ncyc - 1])))
        rec_rate.append(np.mean(rr) if rr else np.nan)
        diss_late.append(np.mean(dd) if dd else np.nan)
    metrics["cyclic_amplitude_scan"] = dict(
        amps=list(amps), recovery_pct_per_cycle=rec_rate,
        last_cycle_dissipation=diss_late)
    ax.plot(np.array(amps) * 100, rec_rate, "go-")
    ax.set_xlabel("strain amplitude (%)")
    ax.set_ylabel("recovery (%/cycle)")
    ax.set_title("cyclic recovery vs amplitude")
    ax.grid(alpha=0.3)
    axd.semilogy(np.array(amps) * 100, diss_late, "mo-")
    axd.set_xlabel("strain amplitude (%)")
    axd.set_ylabel("last-cycle dissipation")
    axd.set_title("stabilized dissipation vs amplitude")
    axd.grid(alpha=0.3)

    fig.savefig(os.path.join(OUT, "matrix_maps.png"), dpi=140,
                bbox_inches="tight")
    with open(os.path.join(OUT, "matrix_metrics.json"), "w") as f:
        json.dump(metrics, f, indent=1)
    print(json.dumps(metrics, indent=1)[:2400])


if __name__ == "__main__":
    main()

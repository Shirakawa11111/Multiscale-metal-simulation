"""Iteration-7 analysis: pore emission (rim-filtered), noise-assisted
nucleation, 3D polycrystal stress curve.
Outputs: results/analysis_iter7.png + results/metrics_iter7.json
"""

import os, sys, json, glob
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D, A_LATTICE
from defect_analysis import find_peaks, find_dislocations

BASE = os.path.join(os.path.dirname(__file__), "..", "results")


def pore_geometry(m, thresh_amp=0.3):
    """Locate the pore: connected low-amplitude region around box center.
    Returns (cx, cy, r_eff). Amplitude proxy: |psi - psi_bar| smoothed."""
    from scipy.ndimage import uniform_filter, label
    amp = uniform_filter(np.abs(m.psi - m.psi_bar), size=9, mode="wrap")
    solid_level = np.percentile(amp, 80)
    mask = amp < thresh_amp * solid_level
    lbl, n = label(mask)
    if n == 0:
        return None
    # component containing the box center (pore was seeded there)
    cmp_id = lbl[m.ny // 2, m.nx // 2]
    if cmp_id == 0:
        sizes = np.bincount(lbl.ravel())[1:]
        cmp_id = int(np.argmax(sizes)) + 1
    ys, xs = np.where(lbl == cmp_id)
    cx, cy = xs.mean() * m.dx, ys.mean() * m.dy
    r_eff = np.sqrt(len(xs) * m.dx * m.dy / np.pi)
    return cx, cy, r_eff


def pore_emission_analysis():
    """Per-snapshot: pore size + emitted cores (outside rim + 2 a0)."""
    out = []
    for fp in sorted(glob.glob(os.path.join(
            BASE, "b45_series_512", "b4_pore_v2", "*.npz"))):
        m = PFC2D.load(fp)
        g = pore_geometry(m)
        pts = find_peaks(m.psi, m.dx, m.dy)
        d = find_dislocations(pts, m.lx, m.ly)
        if g is None:
            emitted = len(d["cores"])
            r_eff = 0.0
        else:
            cx, cy, r_eff = g
            if len(d["cores"]):
                rr = np.sqrt((d["cores"][:, 0] - cx) ** 2
                             + (d["cores"][:, 1] - cy) ** 2)
                emitted = int((rr > r_eff + 2 * A_LATTICE).sum())
            else:
                emitted = 0
        out.append(dict(file=os.path.basename(fp), exx=m.exx,
                        pore_r_over_a0=r_eff / A_LATTICE,
                        cores_total=len(d["cores"]), cores_emitted=emitted))
    return sorted(out, key=lambda r: r["exx"])


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    metrics = {}

    # --- pore emission ---
    pe = pore_emission_analysis()
    metrics["pore_v2"] = pe
    ax = axes[0]
    exx = [r["exx"] * 100 for r in pe]
    ax.plot(exx, [r["cores_emitted"] for r in pe], "ro-", label="emitted")
    ax.plot(exx, [r["cores_total"] for r in pe], "k.--", alpha=0.5,
            label="total detected")
    ax2 = ax.twinx()
    ax2.plot(exx, [r["pore_r_over_a0"] for r in pe], "b^--", alpha=0.6)
    ax2.set_ylabel("pore radius / a0", color="b")
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("# cores")
    ax.set_title("pore_v2: rim-filtered emission")
    ax.legend()
    ax.grid(alpha=0.3)

    # --- noise nucleation ---
    with open(os.path.join(BASE, "b45_series_512", "b4_noise_v2",
                           "summary.json")) as f:
        rows = json.load(f)["rows"]
    eps = np.array([r["exx"] for r in rows])
    cores = np.array([r["cores"] for r in rows])
    sig = np.array([r["sigma"] for r in rows])
    first = float(eps[np.argmax(cores > 0)]) if cores.max() > 0 else None
    metrics["noise_v2"] = dict(noise_amp=0.06, first_nucleation_exx=first,
                               max_cores=int(cores.max()),
                               final_cores=int(cores[-1]))
    ax = axes[1]
    ax.plot(eps * 100, sig, "b-o", ms=3)
    ax2 = ax.twinx()
    ax2.plot(eps * 100, cores, "r-s", ms=4)
    ax2.set_ylabel("# cores", color="r")
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("stress", color="b")
    ax.set_title(f"noise 0.06: first nucleation at "
                 f"{first*100:.1f}%" if first else "noise 0.06: none")
    ax.grid(alpha=0.3)

    # --- 3D polycrystal ---
    with open(os.path.join(BASE, "c3_poly3d_128", "summary.json")) as f:
        c3 = json.load(f)
    eps3 = np.array([r["exx"] for r in c3["rows"]])
    sig3 = np.array([r["sigma"] for r in c3["rows"]])
    el = eps3 <= 0.015
    mod3 = float(np.polyfit(eps3[el], sig3[el], 1)[0])
    metrics["c3_128"] = dict(modulus_early=mod3, n_atoms=c3["n0_peaks"],
                             single_crystal_modulus=0.0353,
                             ratio=mod3 / 0.0353)
    ax = axes[2]
    ax.plot(eps3 * 100, sig3, "g-o", ms=3, label="3D polycrystal 128^3")
    ax.plot(eps3 * 100, 0.0353 * eps3, "k--", alpha=0.5,
            label="single-crystal slope (C2)")
    ax.set_xlabel("strain (%)")
    ax.set_ylabel("stress")
    ax.set_title(f"3D polycrystal: modulus {mod3:.4f} "
                 f"({100*mod3/0.0353:.0f}% of single crystal)")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.savefig(os.path.join(BASE, "analysis_iter7.png"), dpi=140,
                bbox_inches="tight")
    with open(os.path.join(BASE, "metrics_iter7.json"), "w") as f:
        json.dump(metrics, f, indent=1)
    print(json.dumps(metrics, indent=1)[:2200])


if __name__ == "__main__":
    main()

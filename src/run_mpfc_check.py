"""Decisive cross-check: does the rate-sensitivity ladder m(r) survive when
elastic strain is relaxed FAST (MPFC inertial term) instead of only by slow
diffusion?

In MPFC with beta=10, dt=0.5 the acoustic time ~sqrt(beta)~3 << every RELAX
budget (100/400/1600), so the elastic field is always equilibrated and any
remaining rate dependence reflects plastic kinetics, not diffusion. If m still
rises from cold r to near-spinodal r, the glide->GB-creep trend is physical;
if it collapses, the original diffusive m(r) was an artifact.

Endpoints r=-0.35 (cold) and r=-0.21 (near spinodal) x 3 rates x 2 seeds.
Output: results/mpfc_check/<r>_<relax>_<seed>/summary.json
"""

import sys, os, time, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from pfc2d import PFC2D
from defect_analysis import find_peaks, find_dislocations

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "mpfc_check")
DEPS = 0.0025
N_STRAIN = 48          # -> 12% (within 2D validity)
DT = 0.5
BETA = 10.0
N = 512


def run(r, relax, seed):
    tag = f"r{r}_x{relax}_s{seed}"
    out = os.path.join(OUT, tag)
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    m = PFC2D(N, N, r=r, psi_bar=-0.25)
    m.init_random(noise=0.05, seed=seed)
    m.step_mpfc(DT, n=3000, beta=BETA)        # MPFC quench/anneal
    rows = []
    for i in range(N_STRAIN):
        m.apply_strain(DEPS, area_conserving=True)
        m.step_mpfc(DT, n=relax, beta=BETA)
        d = find_dislocations(find_peaks(m.psi, m.dx, m.dy), m.lx, m.ly)
        rows.append(dict(exx=m.exx, sigma=m.stress(), cores=len(d["cores"])))
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(dict(rows=rows, r=r, relax=relax, seed=seed, beta=BETA,
                       wall_s=time.time() - t0), f, indent=1)
    print(f"[{tag}] done {time.time()-t0:.0f}s, sigma@8%={rows[31]['sigma']:.5f}",
          flush=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    for r in (-0.35, -0.21):
        for relax in (100, 400, 1600):
            for seed in (7, 11):
                run(r, relax, seed)
    # analyze m(r)
    rates = 1.0 / np.array([100, 400, 1600])
    res = {}
    for r in (-0.35, -0.21):
        flow = []
        for relax in (100, 400, 1600):
            v = []
            for seed in (7, 11):
                d = json.load(open(os.path.join(OUT, f"r{r}_x{relax}_s{seed}",
                                                "summary.json")))
                rows = d["rows"]
                eps = np.array([x["exx"] for x in rows])
                sig = np.array([x["sigma"] for x in rows])
                v.append(sig[int(np.argmin(np.abs(eps - 0.08)))])
            flow.append(np.mean(v))
        flow = np.array(flow)
        m_exp = (np.polyfit(np.log(rates), np.log(flow), 1)[0]
                 if (flow > 0).all() else float("nan"))
        res[str(r)] = dict(flow_at_8pct=flow.tolist(), m=m_exp)
    res["diffusive_reference"] = {"-0.35": 0.186, "-0.21": 0.402}
    with open(os.path.join(OUT, "m_compare.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1), flush=True)


if __name__ == "__main__":
    main()

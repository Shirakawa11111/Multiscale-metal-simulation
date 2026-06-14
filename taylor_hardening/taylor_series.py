"""Taylor hardening law via a DDD density series — the QUANTITATIVE capstone.

Why this experiment exists
--------------------------
PFC could not forest-harden (every variant softened). DDD (ExaDiS) then showed
the REAL STEM-reconstructed Cu config MULTIPLIES dislocations (2.4x density).
But a single small run is serration-dominated, so an *instantaneous*
sigma-vs-sqrt(rho) correlation is meaningless: during an avalanche the stress
DROPS while the density JUMPS, so the two are anti-correlated point-by-point
(this is exactly why the first attempt gave a spurious negative alpha).

The textbook measurement of the Taylor law is the FLOW STRESS across a DENSITY
SERIES: independent FCC networks at increasing rho0, each run to steady plastic
flow with junctions on; then

        sigma_flow = alpha * mu * b * sqrt(rho)

gives the Taylor coefficient alpha. For Cu the literature value is ~0.3-0.5
(e.g. Bertin et al., MSMSE 2019; Madec/Devincre). A positive alpha in this
range is the quantitative work-hardening signature PFC fundamentally lacks.

Run on HPC (pyexadis built):
  OMP_NUM_THREADS=16 PYTHONPATH=~/BO/exadis_src/python python3 taylor_series.py
"""
import os, sys, json
import numpy as np

EX = os.path.expanduser("~/BO/exadis_src/python")
if EX not in sys.path:
    sys.path.append(EX)
import pyexadis
from pyexadis_base import (ExaDisNet, DisNetManager, SimulateNetworkPerf,
                           CalForce, MobilityLaw, TimeIntegration, Collision,
                           Topology, Remesh)
from pyexadis_utils import dislocation_density

LBOX = float(os.environ.get("LBOX", "10000"))         # b units (~2.55 um)
NUM_LINES = [int(x) for x in os.environ.get("NUM_LINES", "10,20,40,80").split(",")]
ERATE = float(os.environ.get("ERATE", "1e4"))
MAX_STRAIN = float(os.environ.get("MAX_STRAIN", "0.008"))
FLOW_LO = float(os.environ.get("FLOW_LO", "0.004"))   # plastic window for flow stress
OUT = os.environ.get("OUT", "taylor_out")
MU = 54.6e9
B_CU = 2.55e-10


NGRID = int(os.environ.get("NGRID", "32"))
MAXSEG = float(os.environ.get("MAXSEG", "600"))
MINSEG = float(os.environ.get("MINSEG", "150"))


def run_one(num_lines, seed, wdir):
    # DDD_FFT_MODEL (long-range elastic forest interactions, the source of the
    # sqrt(rho) Taylor scaling) + Trapezoid (2 force evals/step, far faster than
    # Subcycling's many) -> tractable for a multi-density demonstration series.
    state = {"crystal": "fcc", "burgmag": B_CU, "mu": MU, "nu": 0.324, "a": 6.0,
             "maxseg": MAXSEG, "minseg": MINSEG, "rtol": 10.0, "rann": 10.0,
             "nextdt": 1e-10, "maxdt": 1e-9}
    G = ExaDisNet()
    G.generate_line_config("fcc", LBOX, num_lines, maxseg=MAXSEG, seed=seed)
    net = DisNetManager(G)
    rho0 = dislocation_density(net, state["burgmag"])

    calforce = CalForce(force_mode="DDD_FFT_MODEL", state=state, Ngrid=NGRID,
                        cell=net.cell)
    mobility = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0,
                           Mscrew=64103.0, vmax=4000.0)
    timeint = TimeIntegration(integrator="Trapezoid", state=state,
                              force=calforce, mobility=mobility)
    collision = Collision(collision_mode="Retroactive", state=state)
    topology = Topology(topology_mode="TopologyParallel", state=state,
                        force=calforce, mobility=mobility)   # JUNCTIONS ON
    remesh = Remesh(remesh_rule="LengthBased", state=state)
    os.makedirs(wdir, exist_ok=True)
    sim = SimulateNetworkPerf(calforce=calforce, mobility=mobility,
                              timeint=timeint, collision=collision,
                              topology=topology, remesh=remesh, cross_slip=None,
                              vis=None, loading_mode="strain_rate", erate=ERATE,
                              edir=np.array([0., 0., 1.]), max_strain=MAX_STRAIN,
                              burgmag=state["burgmag"], state=state,
                              print_freq=10, plot_freq=10**9, write_freq=10**9,
                              write_dir=wdir)
    sim.run(net, state)
    rhof = dislocation_density(net, state["burgmag"])
    # parse the reliably-written stress-strain table:
    # columns [istep, strain, stress, density, elapsed]
    dat = os.path.join(wdir, "stress_strain_dens.dat")
    res = np.loadtxt(dat, ndmin=2)
    return res, float(rho0), float(rhof)


def main():
    pyexadis.initialize()
    os.makedirs(OUT, exist_ok=True)
    series, allres = [], {}

    def save(series, allres):
        # Taylor fit through origin: sigma = alpha * mu * b * sqrt(rho)
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
                   erate=ERATE, max_strain=MAX_STRAIN, flow_window_lo=FLOW_LO,
                   lbox_b=LBOX,
                   note="Taylor flow-stress density series, FCC Cu, ExaDiS DDD "
                        "with junctions (TopologyParallel). sigma_flow vs "
                        "mu*b*sqrt(rho).")
        with open(os.path.join(OUT, "taylor_series.json"), "w") as f:
            json.dump(out, f, indent=1)
        with open(os.path.join(OUT, "taylor_curves.json"), "w") as f:
            json.dump(allres, f)
        return alpha, r2

    for k, nl in enumerate(NUM_LINES):
        print(f"\n==== density run {k+1}/{len(NUM_LINES)}: num_lines={nl} ====",
              flush=True)
        try:
            res, rho0, rhof = run_one(nl, seed=1234 + k,
                                      wdir=os.path.join(OUT, f"d{nl}"))
        except Exception as e:
            print(f"  density {nl} FAILED: {e}", flush=True)
            continue
        strain, stress, dens = res[:, 1], res[:, 2], res[:, 3]
        mask = strain >= FLOW_LO
        flow = float(np.mean(np.abs(stress[mask]))) if mask.any() \
            else float(np.abs(stress[-1]))
        flow_std = float(np.std(np.abs(stress[mask]))) if mask.any() else 0.0
        rho_flow = float(np.mean(dens[mask])) if mask.any() else rhof
        series.append(dict(num_lines=nl, rho0=rho0, rho_final=rhof,
                           rho_flow=rho_flow, flow_stress=flow,
                           flow_std=flow_std, n_pts=int(len(strain))))
        allres[str(nl)] = res.tolist()
        # per-density result (one process per density avoids a pyexadis
        # double-free when multiple networks are built in one process)
        with open(os.path.join(OUT, f"d{nl}", "result.json"), "w") as f:
            json.dump(dict(series[-1], curve=res.tolist()), f)
        a, r2 = save(series, allres)   # incremental save after each density
        print(f"  rho0={rho0:.3e} rho_flow={rho_flow:.3e} "
              f"flow_stress={flow/1e6:.1f}+/-{flow_std/1e6:.1f} MPa "
              f"[running alpha={a:.3f}]", flush=True)

    alpha, r2 = save(series, allres)
    print(f"\n=== TAYLOR: alpha={alpha:.3f}, R^2={r2:.3f} "
          f"(Cu literature ~0.3-0.5) ===")
    pyexadis.finalize()


if __name__ == "__main__":
    main()

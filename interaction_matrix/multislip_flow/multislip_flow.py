"""Multi-slip COLL on/off flow-stress cell -- the decisive collinear-dominance test.

Per the 4-expert adversarial scrutiny: the pairwise MFP assay structurally cannot host the
multi-slip collinear coefficient. The single test that cannot return ambiguous:

  Primary mobile system m (MULTIPLYING Frank-Read sources) gliding through an EVOLVABLE forest
  of a partner system f, under STRAIN-RATE control (flow stress is an OUTPUT), with CrossSlip ON
  (required for carrier replenishment). Measure steady-state flow stress tau_flow and the stored
  density. Compare forest TYPES at matched density:
     coll_opp  : collinear partner, opposite sense (bf=-bm)  -> ANNIHILATION on
     coll_same : collinear partner, same sense (bf=+bm)      -> ANNIHILATION off (single-bit toggle)
     glissile / hirth : other junction partners (cross-type, canonical ratio reference)

  CONFIRM collinear dominance: coll_opp flow stress highest, gap grows with rho_f, and
     tau_flow(coll_opp)/tau_flow(glissile) -> canonical sqrt(0.62/0.12) ~ 2.3x.
  REFUTE: at a flow plateau (strain>1e-3), all gates passed, coll_opp == coll_same == others.

  JTYPE=coll_opp RHO_F=3e12 ERATE=1e3 LBOX=16000 NSTEPS=4000 SEED=1 OUT=ms/coll_opp_s1 \
    PYTHONPATH=~/BO/exadis_src/python python3 multislip_flow.py
"""
import os, sys, json
import numpy as np

EX = os.path.expanduser("~/BO/exadis_src/python")
if EX not in sys.path:
    sys.path.append(EX)
import pyexadis
from pyexadis_base import (ExaDisNet, DisNetManager, SimulateNetwork, CalForce,
                           MobilityLaw, TimeIntegration, Collision, Topology, Remesh,
                           CrossSlip, NodeConstraints)
from pyexadis_utils import insert_infinite_line, insert_frank_read_src, dislocation_density

# verified ExaDiS FCC 12-system table
BI = np.array([[0, 1, -1], [1, 0, -1], [1, -1, 0], [0, 1, -1], [1, 0, 1], [1, 1, 0],
               [0, 1, 1], [1, 0, -1], [1, 1, 0], [0, 1, 1], [1, 0, 1], [1, -1, 0]], float)
NI = np.array([[1, 1, 1]] * 3 + [[-1, 1, 1]] * 3 + [[1, -1, 1]] * 3 + [[1, 1, -1]] * 3, float)
B_CU, MU, NU = 2.55e-10, 54.6e9, 0.324


def envf(k, d): return float(os.environ.get(k, d))
def envi(k, d): return int(os.environ.get(k, d))


JTYPE  = os.environ.get("JTYPE", "coll_opp")
LBOX   = envf("LBOX", "16000")
RHO_F  = envf("RHO_F", "3e12")
NFOREST = envi("NFOREST", "0")
ERATE  = envf("ERATE", "1e3")
KM     = envi("KM", "8")             # FR source count on primary
LSRC   = envf("LSRC", "0")           # FR source length (0 -> 0.18*LBOX)
NSTEPS = envi("NSTEPS", "4000")
NREL   = envi("NREL", "200")
REC    = envi("REC", "25")
SEED   = envi("SEED", "1")
MAXSEG = envf("MAXSEG", "100")
NGRID  = envi("NGRID", "32")
RANN   = envf("RANN", "10")
XSLIP  = os.environ.get("XSLIP", "1") == "1"
OUT    = os.environ.get("OUT", "ms_out")
rng = np.random.default_rng(SEED)


def hat(v):
    v = np.asarray(v, float); n = np.linalg.norm(v)
    return v / n if n > 0 else v
def cos_par(u, v):
    return abs(float(np.dot(hat(u), hat(v))))


def pick_pair():
    bm, nm = BI[0].copy(), NI[0].copy()      # primary m = [0,1,-1](111)
    j = JTYPE.lower()
    if j in ("coll_opp", "coll_same"):
        nf = NI[3].copy()                    # 2nd {111} sharing bm
        bf = (-1.0 if j == "coll_opp" else 1.0) * bm.copy()
    elif j == "glissile":
        bf, nf = BI[4].copy(), NI[4].copy()
    elif j == "hirth":
        bf, nf = BI[6].copy(), NI[6].copy()
    elif j == "lomer":
        bf, nf = BI[5].copy(), NI[5].copy()
    else:
        raise ValueError("JTYPE %s" % JTYPE)
    return bm, nm, bf, nf


bm, nm, bf, nf = pick_pair()
b3a, b3b = bm + bf, bm - bf
FAM_J = [v for v in (b3a, b3b) if np.linalg.norm(v) > 1e-6 and cos_par(v, bm) < 0.94 and cos_par(v, bf) < 0.94]
# loading axis: uniaxial along (bhat+nhat)/sqrt2 -> Schmid 0.5 on the primary system
EDIR = hat(hat(bm) + hat(nm))


def nforest_for_rho():
    if NFOREST > 0:
        return NFOREST
    perline = 1.0 / (LBOX * B_CU) ** 2
    return max(2, int(round(RHO_F / perline)))


def build():
    cell = pyexadis.Cell(h=LBOX * np.eye(3), is_periodic=[True, True, True])
    C = np.array(cell.center())
    nodes, segs = [], []
    # forest: evolvable infinite lines of partner system f (NOT pinned)
    nF = nforest_for_rho()
    for _ in range(nF):
        o = C + (rng.random(3) - 0.5) * 0.85 * LBOX
        th = float(rng.choice([0, 30, 60, 90]))
        try:
            ok = insert_infinite_line(cell, nodes, segs, bf, nf, o, theta=th, maxseg=MAXSEG, trial=True)
        except TypeError:
            ok = 1.0
        if ok and ok > 0:
            insert_infinite_line(cell, nodes, segs, bf, nf, o, theta=th, maxseg=MAXSEG)
    # primary mobile: multiplying Frank-Read sources
    Lsrc = LSRC if LSRC > 0 else 0.18 * LBOX
    placed = 0
    for k in range(KM):
        c = C + (rng.random(3) - 0.5) * 0.7 * LBOX
        th = float(rng.choice([30, 60]))
        try:
            insert_frank_read_src(cell, nodes, segs, bm, nm, Lsrc, c, theta=th)
            placed += 1
        except Exception as e:
            if k == 0:
                print("FR insert err:", e, flush=True)
    if placed == 0:
        raise RuntimeError("0 FR sources placed")
    nodes = np.array(nodes); segs = np.array(segs)
    return cell, nodes, segs, nF, placed


def modules(state, cell):
    cf = CalForce(force_mode="DDD_FFT_MODEL", state=state, Ngrid=NGRID, cell=cell)
    mob = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    ti = TimeIntegration(integrator="Trapezoid", state=state, force=cf, mobility=mob)
    col = Collision(collision_mode="Retroactive", state=state)
    topo = Topology(topology_mode="TopologyParallel", state=state, force=cf, mobility=mob)
    rm = Remesh(remesh_rule="LengthBased", state=state)
    xs = CrossSlip(state=state, cross_slip_mode="ForceBasedParallel", force=cf, mobility=mob) if XSLIP else None
    return cf, mob, ti, col, topo, rm, xs


def main():
    pyexadis.initialize()
    os.makedirs(OUT, exist_ok=True)
    cell, nodes, segs, nF, nsrc = build()
    net = DisNetManager(ExaDisNet(cell, nodes, segs))
    state = {"crystal": "fcc", "burgmag": B_CU, "mu": MU, "nu": NU, "a": 6.0,
             "maxseg": MAXSEG, "minseg": MAXSEG / 4, "rtol": 10.0, "rann": RANN,
             "nextdt": 1e-12, "maxdt": 1e-10}
    rho_built = dislocation_density(net, B_CU)

    # PHASE A: brief zero-stress settle
    cf, mob, ti, col, topo, rm, xs = modules(state, net.cell)
    SimulateNetwork(calforce=cf, mobility=mob, timeint=ti, collision=col, topology=topo,
                    remesh=rm, cross_slip=xs, vis=None, state=state, burgmag=B_CU,
                    loading_mode="stress", applied_stress=np.zeros(6), max_step=NREL,
                    print_freq=10**9, plot_freq=10**9, write_freq=10**9, write_dir=OUT).run(net, state)
    rho_settled = dislocation_density(net, B_CU)

    # PHASE B: strain-rate loading along EDIR, flow stress is the output
    cf2, mob2, ti2, col2, topo2, rm2, xs2 = modules(state, net.cell)
    state["edir"] = EDIR.copy()
    sim = SimulateNetwork(calforce=cf2, mobility=mob2, timeint=ti2, collision=col2, topology=topo2,
                          remesh=rm2, cross_slip=xs2, vis=None, state=state, burgmag=B_CU,
                          loading_mode="strain_rate", erate=ERATE, edir=EDIR.copy(),
                          max_step=NSTEPS, print_freq=REC, plot_freq=10**9, write_freq=10**9,
                          write_dir=OUT)
    sim.run(net, state)

    # results: [istep, strain, stress, density, elapsed]
    R = np.array(sim.results) if len(sim.results) else np.zeros((0, 5))
    rho_end = dislocation_density(net, B_CU)
    out = dict(jtype=JTYPE, rho_f_target=RHO_F, nforest=nF, nsrc=nsrc, erate=ERATE,
               lbox=LBOX, ngrid=NGRID, xslip=XSLIP, seed=SEED, edir=EDIR.tolist(),
               rho_built=rho_built, rho_settled=rho_settled, rho_end=rho_end,
               forest_drift=float((rho_end - rho_settled) / rho_settled) if rho_settled else 0.0)
    if len(R):
        strain = R[:, 1]; stress = R[:, 2]; dens = R[:, 3]
        # flow stress = mean stress over the plateau (last 40% of strain, strain>1e-3 reached?)
        plateau = strain > max(0.4 * strain.max(), 1e-3)
        tau_flow = float(np.mean(stress[plateau])) if plateau.any() else float(stress[-1])
        out.update(strain_max=float(strain.max()), tau_flow_Pa=tau_flow,
                   tau_flow_MPa=tau_flow / 1e6, stress_end_MPa=float(stress[-1] / 1e6),
                   density_end=float(dens[-1]), plateau_reached=bool(strain.max() > 1e-3),
                   series=[[float(s), float(t / 1e6), float(d)] for s, t, d in zip(strain[::2], stress[::2], dens[::2])])
    json.dump(out, open(os.path.join(OUT, "flow.json"), "w"), indent=1)
    sm = out.get("strain_max", 0); tf = out.get("tau_flow_MPa", float("nan"))
    print(f"FLOW {JTYPE} s{SEED} rho_f={RHO_F}: tau_flow={tf:.1f}MPa strain_max={sm:.2e} "
          f"plateau={out.get('plateau_reached')} drift={out['forest_drift']:+.2f} nF={nF} nsrc={nsrc}", flush=True)
    pyexadis.finalize()


if __name__ == "__main__":
    main()

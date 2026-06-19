"""Phase 1 step 4-5: measure the collinear REMOBILIZATION critical stress.

Build the controlled collinear annihilation pair (confirmed mechanism), relax at
ZERO stress to form the partial-annihilation residual, measure the residual
length l_res, then apply a FIXED resolved shear tau on the mobile system
(sigma = tau (b^n + n^b)) and run quasi-statically. Detect whether the residual
REMOBILIZES (the unpinned nodes sweep large distances / plastic strain takes off)
or stays trapped. Sweeping tau brackets the critical remobilization stress tau_c;
sweeping segment length L gives the scaling tau_c ~ (mu b / l_res) ln(l_res/b).

  LSEG=<b> TAU_MPA=<MPa> PYTHONPATH=~/BO/exadis_src/python python3 collinear_strength.py
"""
import os, sys, json
import numpy as np

EX = os.path.expanduser("~/BO/exadis_src/python")
if EX not in sys.path:
    sys.path.append(EX)
import pyexadis
from pyexadis_base import (ExaDisNet, DisNetManager, SimulateNetwork, CalForce,
                           MobilityLaw, TimeIntegration, Collision, Topology, Remesh,
                           NodeConstraints)

B_CU, MU = 2.55e-10, 54.6e9
LBOX = float(os.environ.get("LBOX", "6000"))
LSEG = float(os.environ.get("LSEG", "1500"))
PHI = float(os.environ.get("PHI", "45"))
GAP = float(os.environ.get("GAP", "6"))
MAXSEG = float(os.environ.get("MAXSEG", "80"))
NREL = int(os.environ.get("NREL", "500"))         # zero-stress relaxation steps
NLOAD = int(os.environ.get("NLOAD", "600"))       # loading steps
TAU_MPA = float(os.environ.get("TAU_MPA", "100"))
OUT = os.environ.get("OUT", "cs_out")

b = np.array([0., 1., -1.])
n1 = np.array([1., 1., 1.])
n2 = np.array([-1., 1., 1.])
bh = b / np.linalg.norm(b)
n1h = n1 / np.linalg.norm(n1)


def add_segment(center, n, nodes, segs, burg):
    e = np.cross(n, b); e = e / np.linalg.norm(e)
    xi = np.cos(np.radians(PHI)) * bh + np.sin(np.radians(PHI)) * e
    xi = xi / np.linalg.norm(xi)
    nn = max(4, int(round(LSEG / MAXSEG)))
    i0 = len(nodes)
    nhat = n / np.linalg.norm(n)
    for j in range(nn + 1):
        p = center + (j / nn - 0.5) * LSEG * xi
        con = int(NodeConstraints.PINNED_NODE) if j in (0, nn) else int(NodeConstraints.UNCONSTRAINED)
        nodes.append(np.concatenate((p, [con])))
    for j in range(nn):
        segs.append(np.concatenate(([i0 + j, i0 + j + 1], burg, nhat)))
    return i0, len(nodes)


def line_length(pos, nid):
    return float(sum(np.linalg.norm(pos[int(a)] - pos[int(c)]) for a, c in nid))


def make_modules(state, cell):
    cf = CalForce(force_mode="DDD_FFT_MODEL", state=state, Ngrid=32, cell=cell)
    mob = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    ti = TimeIntegration(integrator="Trapezoid", state=state, force=cf, mobility=mob)
    col = Collision(collision_mode="Retroactive", state=state)
    topo = Topology(topology_mode="TopologyParallel", state=state, force=cf, mobility=mob)
    rm = Remesh(remesh_rule="LengthBased", state=state)
    return cf, mob, ti, col, topo, rm


def main():
    pyexadis.initialize()
    os.makedirs(OUT, exist_ok=True)
    cell = pyexadis.Cell(h=LBOX * np.eye(3), is_periodic=[True, True, True])
    C = np.array(cell.center())
    neutral = (n1 / np.linalg.norm(n1) + n2 / np.linalg.norm(n2))
    neutral = neutral / np.linalg.norm(neutral)
    nodes, segs = [], []
    add_segment(C + 0.5 * GAP * neutral, n1, nodes, segs, b)
    add_segment(C - 0.5 * GAP * neutral, n2, nodes, segs, -b)
    nodes = np.array(nodes); segs = np.array(segs)
    L0 = line_length(nodes[:, :3], segs[:, :2])

    net = DisNetManager(ExaDisNet(cell, nodes, segs))
    state = {"crystal": "fcc", "burgmag": B_CU, "mu": MU, "nu": 0.324, "a": 6.0,
             "maxseg": MAXSEG, "minseg": MAXSEG / 4, "rtol": 10.0, "rann": 10.0,
             "nextdt": 1e-11, "maxdt": 1e-9}
    cf, mob, ti, col, topo, rm = make_modules(state, net.cell)

    # ---- Phase A: zero-stress relaxation (form residual) ----
    simA = SimulateNetwork(calforce=cf, mobility=mob, timeint=ti, collision=col,
                           topology=topo, remesh=rm, vis=None, state=state, max_step=NREL,
                           loading_mode="stress", applied_stress=np.zeros(6),
                           print_freq=10**9, plot_freq=10**9, write_freq=10**9, write_dir=OUT)
    simA.run(net, state)
    dA = net.get_disnet(ExaDisNet).export_data()
    posA = dA["nodes"]["positions"]; nidA = dA["segs"]["nodeids"]; conA = dA["nodes"]["constraints"]
    L_res = line_length(posA, nidA)
    unpin = (conA[:, 0] != int(NodeConstraints.PINNED_NODE))
    posA_unpin = posA[unpin].copy()

    # ---- Phase B: fixed resolved shear tau on (b, n1) ----
    tau = TAU_MPA * 1e6
    sig = tau * (np.outer(bh, n1h) + np.outer(n1h, bh))   # resolved shear tensor
    voigt = np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[1, 2], sig[0, 2], sig[0, 1]])
    cf2, mob2, ti2, col2, topo2, rm2 = make_modules(state, net.cell)
    simB = SimulateNetwork(calforce=cf2, mobility=mob2, timeint=ti2, collision=col2,
                           topology=topo2, remesh=rm2, vis=None, state=state, max_step=NLOAD,
                           loading_mode="stress", applied_stress=voigt,
                           print_freq=10**9, plot_freq=10**9, write_freq=10**9, write_dir=OUT)
    simB.run(net, state)
    dB = net.get_disnet(ExaDisNet).export_data()
    posB = dB["nodes"]["positions"]; nidB = dB["segs"]["nodeids"]; conB = dB["nodes"]["constraints"]
    L_fin = line_length(posB, nidB)
    unpinB = (conB[:, 0] != int(NodeConstraints.PINNED_NODE))

    # remobilization signal: how far did the unpinned dislocation sweep?
    # (mean displacement of unpinned nodes, in b; and line-length growth)
    nB = posB[unpinB]
    sweep = float(np.mean(np.linalg.norm(nB - nB.mean(0), axis=1))) if len(nB) else 0.0
    # robust proxy: max coordinate span growth of the unpinned cloud B vs A
    span = lambda P: float(np.linalg.norm(P.max(0) - P.min(0))) if len(P) else 0.0
    spanA, spanB = span(posA_unpin), span(nB)
    length_growth = (L_fin - L_res) / L_res if L_res else 0.0
    remobilized = bool((spanB - spanA) / max(spanA, 1.0) > 0.5 or length_growth > 0.5)

    out = dict(protocol="collinear_strength", LSEG=LSEG, tau_MPa=TAU_MPA,
               L0=float(L0), L_res=float(L_res), L_fin=float(L_fin),
               anneal_fraction=float(L_res / L0),
               span_unpin_before=spanA, span_unpin_after=spanB,
               length_growth=float(length_growth), remobilized=remobilized,
               n_unpin_after=int(unpinB.sum()))
    json.dump(out, open(os.path.join(OUT, "strength.json"), "w"), indent=1)
    print(f"COLLINEAR L={LSEG} tau={TAU_MPA}MPa: L_res/L0={L_res/L0:.2f} "
          f"span {spanA:.0f}->{spanB:.0f} growth={length_growth:+.2f} "
          f"remobilized={remobilized}", flush=True)
    pyexadis.finalize()


if __name__ == "__main__":
    main()

"""Generalized controlled binary-reaction strength for ANY junction type.

Two finite pinned-end segments on two chosen slip systems (B1,N1)/(B2,N2),
crossing within rann at the box centre, relaxed at zero stress (forms the
reaction product: junction / annihilation residual / dipole), then a fixed
resolved shear on the mobile system (B1,N1) tests remobilization. Used to rank
the six FCC junction types with one consistent protocol.

  B1=0,1,-1 N1=1,1,1 B2=... N2=... LSEG=<b> TAU_MPA=<MPa> python3 binary_strength.py
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


def vec(name, default):
    return np.array([float(x) for x in os.environ.get(name, default).split(",")])


B1, N1, S1 = vec("B1", "0,1,-1"), vec("N1", "1,1,1"), float(os.environ.get("S1", "1"))
B2, N2, S2 = vec("B2", "0,1,-1"), vec("N2", "-1,1,1"), float(os.environ.get("S2", "-1"))
LBOX = float(os.environ.get("LBOX", "6000"))
LSEG = float(os.environ.get("LSEG", "1500"))
PHI = float(os.environ.get("PHI", "45"))
GAP = float(os.environ.get("GAP", "6"))
MAXSEG = float(os.environ.get("MAXSEG", "80"))
NREL = int(os.environ.get("NREL", "500"))
NLOAD = int(os.environ.get("NLOAD", "1500"))
TAU_MPA = float(os.environ.get("TAU_MPA", "100"))
FORCE = os.environ.get("FORCE", "DDD_FFT_MODEL")
JTYPE = os.environ.get("JTYPE", "?")
OUT = os.environ.get("OUT", "bs_out")


def add_segment(center, burg, n, sgn, nodes, segs):
    bh = burg / np.linalg.norm(burg)
    e = np.cross(n, burg); e = e / np.linalg.norm(e)
    xi = sgn * (np.cos(np.radians(PHI)) * bh + np.sin(np.radians(PHI)) * e)
    xi = xi / np.linalg.norm(xi)
    nn = max(4, int(round(LSEG / MAXSEG)))
    i0 = len(nodes); nhat = n / np.linalg.norm(n)
    for j in range(nn + 1):
        p = center + (j / nn - 0.5) * LSEG * xi
        con = int(NodeConstraints.PINNED_NODE) if j in (0, nn) else int(NodeConstraints.UNCONSTRAINED)
        nodes.append(np.concatenate((p, [con])))
    for j in range(nn):
        segs.append(np.concatenate(([i0 + j, i0 + j + 1], burg, nhat)))


def llen(pos, nid):
    return float(sum(np.linalg.norm(pos[int(a)] - pos[int(c)]) for a, c in nid))


def modules(state, cell):
    cf = CalForce(force_mode=FORCE, state=state, Ngrid=32, cell=cell)
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
    neutral = N1 / np.linalg.norm(N1) + N2 / np.linalg.norm(N2)
    neutral = neutral / np.linalg.norm(neutral)
    nodes, segs = [], []
    add_segment(C + 0.5 * GAP * neutral, B1, N1, S1, nodes, segs)
    add_segment(C - 0.5 * GAP * neutral, B2, N2, S2, nodes, segs)
    nodes = np.array(nodes); segs = np.array(segs)
    L0 = llen(nodes[:, :3], segs[:, :2])

    net = DisNetManager(ExaDisNet(cell, nodes, segs))
    state = {"crystal": "fcc", "burgmag": B_CU, "mu": MU, "nu": 0.324, "a": 6.0,
             "maxseg": MAXSEG, "minseg": MAXSEG / 4, "rtol": 10.0, "rann": 10.0,
             "nextdt": 1e-11, "maxdt": 1e-9}
    cf, mob, ti, col, topo, rm = modules(state, net.cell)
    simA = SimulateNetwork(calforce=cf, mobility=mob, timeint=ti, collision=col, topology=topo,
                           remesh=rm, vis=None, state=state, max_step=NREL, loading_mode="stress",
                           applied_stress=np.zeros(6), print_freq=10**9, plot_freq=10**9,
                           write_freq=10**9, write_dir=OUT)
    simA.run(net, state)
    dA = net.get_disnet(ExaDisNet).export_data()
    posA, nidA, conA = dA["nodes"]["positions"], dA["segs"]["nodeids"], dA["nodes"]["constraints"]
    L_res = llen(posA, nidA)
    unpinA = posA[conA[:, 0] != int(NodeConstraints.PINNED_NODE)]
    spanA = float(np.linalg.norm(unpinA.max(0) - unpinA.min(0))) if len(unpinA) else 0.0

    bh1, nh1 = B1 / np.linalg.norm(B1), N1 / np.linalg.norm(N1)
    sig = TAU_MPA * 1e6 * (np.outer(bh1, nh1) + np.outer(nh1, bh1))
    voigt = np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[1, 2], sig[0, 2], sig[0, 1]])
    cf2, mob2, ti2, col2, topo2, rm2 = modules(state, net.cell)
    simB = SimulateNetwork(calforce=cf2, mobility=mob2, timeint=ti2, collision=col2, topology=topo2,
                           remesh=rm2, vis=None, state=state, max_step=NLOAD, loading_mode="stress",
                           applied_stress=voigt, print_freq=10**9, plot_freq=10**9,
                           write_freq=10**9, write_dir=OUT)
    simB.run(net, state)
    dB = net.get_disnet(ExaDisNet).export_data()
    posB, nidB, conB = dB["nodes"]["positions"], dB["segs"]["nodeids"], dB["nodes"]["constraints"]
    L_fin = llen(posB, nidB)
    unpinB = posB[conB[:, 0] != int(NodeConstraints.PINNED_NODE)]
    spanB = float(np.linalg.norm(unpinB.max(0) - unpinB.min(0))) if len(unpinB) else 0.0
    length_growth = (L_fin - L_res) / L_res if L_res else 0.0
    remobilized = bool(spanB > 0.45 * LBOX or length_growth > 1.5)

    out = dict(jtype=JTYPE, B1=B1.tolist(), N1=N1.tolist(), B2=B2.tolist(), N2=N2.tolist(),
               LSEG=LSEG, tau_MPa=TAU_MPA, L0=float(L0), L_res=float(L_res),
               anneal_fraction=float(L_res / L0), span_before=spanA, span_after=spanB,
               length_growth=float(length_growth), remobilized=remobilized)
    json.dump(out, open(os.path.join(OUT, "strength.json"), "w"), indent=1)
    print(f"{JTYPE:>10} tau={TAU_MPA:.0f}MPa: L_res/L0={L_res/L0:.2f} span {spanA:.0f}->{spanB:.0f} "
          f"growth={length_growth:+.2f} remob={remobilized}", flush=True)
    pyexadis.finalize()


if __name__ == "__main__":
    main()

"""Per-type MECHANISM verification (zero-stress, no loading).

Build two finite pinned-end segments on (B1,N1)/(B2,N2), each line at PHI to its
OWN Burgers in its own plane, with an explicit per-segment line SENSE (S1,S2=+/-1).
Relax at zero stress with DDD_FFT + Collision + Topology and classify the reaction
from a topology observer:
  partial_annihilation  (length drops)        -> collinear / self annihilation
  junction              (degree>=3 node, len ~kept) -> Hirth/glissile/Lomer/coplanar
  pass_through / fully_deleted / none
Scanning S1,S2 finds the geometry that produces the canonical reaction for each
junction type, BEFORE any strength measurement.

  B1=.. N1=.. S1=1 B2=.. N2=.. S2=-1 JTYPE=collinear python3 binary_observe.py
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


def vec(name, d):
    return np.array([float(x) for x in os.environ.get(name, d).split(",")])


B1, N1, S1 = vec("B1", "0,1,-1"), vec("N1", "1,1,1"), float(os.environ.get("S1", "1"))
B2, N2, S2 = vec("B2", "0,1,-1"), vec("N2", "-1,1,1"), float(os.environ.get("S2", "-1"))
LBOX = float(os.environ.get("LBOX", "6000"))
LSEG = float(os.environ.get("LSEG", "1500"))
PHI = float(os.environ.get("PHI", "45"))
GAP = float(os.environ.get("GAP", "6"))
MAXSEG = float(os.environ.get("MAXSEG", "80"))
NREL = int(os.environ.get("NREL", "500"))
JTYPE = os.environ.get("JTYPE", "?")
OUT = os.environ.get("OUT", "bo_out")


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
    cf = CalForce(force_mode="DDD_FFT_MODEL", state=state, Ngrid=32, cell=net.cell)
    mob = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    ti = TimeIntegration(integrator="Trapezoid", state=state, force=cf, mobility=mob)
    col = Collision(collision_mode="Retroactive", state=state)
    topo = Topology(topology_mode="TopologyParallel", state=state, force=cf, mobility=mob)
    rm = Remesh(remesh_rule="LengthBased", state=state)
    sim = SimulateNetwork(calforce=cf, mobility=mob, timeint=ti, collision=col, topology=topo,
                          remesh=rm, vis=None, state=state, max_step=NREL, loading_mode="stress",
                          applied_stress=np.zeros(6), print_freq=10**9, plot_freq=10**9,
                          write_freq=10**9, write_dir=OUT)
    sim.run(net, state)
    d = net.get_disnet(ExaDisNet).export_data()
    pos, nid = d["nodes"]["positions"], d["segs"]["nodeids"]
    Lf, sf = llen(pos, nid), len(nid)
    deg = {}
    for a, c in nid:
        deg[int(a)] = deg.get(int(a), 0) + 1
        deg[int(c)] = deg.get(int(c), 0) + 1
    njunc = sum(1 for v in deg.values() if v >= 3)
    frac = Lf / L0 if L0 else 0
    if sf == 0:
        mech = "fully_deleted"
    elif frac < 0.85:
        mech = "partial_annihilation"
    elif njunc > 0:
        mech = "junction"
    elif abs(frac - 1) < 0.06:
        mech = "pass_through_none"
    else:
        mech = "other"
    out = dict(jtype=JTYPE, S1=S1, S2=S2, length_fraction=float(frac),
               n_junction_nodes=int(njunc), n_segs_after=int(sf), mechanism=mech)
    json.dump(out, open(os.path.join(OUT, "obs.json"), "w"), indent=1)
    print(f"{JTYPE:>10} S=({S1:+.0f},{S2:+.0f}): len_frac={frac:.2f} junc={njunc} -> {mech}", flush=True)
    pyexadis.finalize()


if __name__ == "__main__":
    main()

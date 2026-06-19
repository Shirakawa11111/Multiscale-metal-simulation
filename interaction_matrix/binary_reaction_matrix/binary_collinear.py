"""Phase 1 (mechanism): controlled binary COLLINEAR reaction in ExaDiS.

Two FINITE dislocation segments with the SAME Burgers vector on two intersecting
{111} planes, both pinned at the ends, free in the middle, crossing near the box
centre at ~PHI to the common Burgers. At ZERO applied stress, with DDD_FFT +
Collision + Topology, we let the system relax and OBSERVE the mechanism (do they
attract, partially annihilate, leave a residual segment, form a junction, pass
through, or get fully deleted?). We do NOT measure strength yet — first confirm
the physics is the canonical collinear reaction (Madec/Devincre).

A topology observer JSON records line-length and node/segment changes so the
mechanism is diagnosed, not eyeballed.

  LSEG=<b> PHI=<deg> GAP=<b> PYTHONPATH=~/BO/exadis_src/python python3 binary_collinear.py
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
LBOX = float(os.environ.get("LBOX", "4000"))
LSEG = float(os.environ.get("LSEG", "1500"))      # segment length (b)
PHI = float(os.environ.get("PHI", "45"))          # angle of line to common Burgers (deg)
GAP = float(os.environ.get("GAP", "40"))          # initial separation of the two segments (b)
MAXSEG = float(os.environ.get("MAXSEG", "80"))
NSTEPS = int(os.environ.get("NSTEPS", "400"))
FORCE = os.environ.get("FORCE", "DDD_FFT_MODEL")
TOPO = os.environ.get("TOPO", "1")                # 1 = TopologyParallel on
COLL = os.environ.get("COLL", "1")                # 1 = collision on
OUT = os.environ.get("OUT", "bc_out")

# collinear pair: SAME Burgers, two intersecting {111} planes (b in both)
b = np.array([0., 1., -1.])
n1 = np.array([1., 1., 1.])
n2 = np.array([-1., 1., 1.])
bh = b / np.linalg.norm(b)


def add_segment(center, n, nodes, segs, burg):
    e = np.cross(n, b); e = e / np.linalg.norm(e)          # edge dir in plane
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


def line_length(nodes, segs):
    tot = 0.0
    for s in segs:
        a, c = int(s[0]), int(s[1])
        tot += np.linalg.norm(nodes[a, :3] - nodes[c, :3])
    return tot


def main():
    pyexadis.initialize()
    os.makedirs(OUT, exist_ok=True)
    # DDD_FFT needs full PBC; box is large vs the segment so periodic images are
    # far, and the pinned ends anchor the segments (effectively finite reaction).
    cell = pyexadis.Cell(h=LBOX * np.eye(3), is_periodic=[True, True, True])
    C = np.array(cell.center())
    # offset the two segments along a neutral direction so they start close but distinct
    neutral = (n1 / np.linalg.norm(n1) + n2 / np.linalg.norm(n2))
    neutral = neutral / np.linalg.norm(neutral)
    nodes, segs = [], []
    # collinear ANNIHILATION pair: same Burgers line but OPPOSITE sense (+b, -b)
    # so the reaction product where they meet is b+(-b)=0 (annihilation+residual).
    add_segment(C + 0.5 * GAP * neutral, n1, nodes, segs, b)
    add_segment(C - 0.5 * GAP * neutral, n2, nodes, segs, -b)
    nodes = np.array(nodes); segs = np.array(segs)
    n0, s0 = len(nodes), len(segs)
    L0 = line_length(nodes, segs)

    net = DisNetManager(ExaDisNet(cell, nodes, segs))
    state = {"crystal": "fcc", "burgmag": B_CU, "mu": MU, "nu": 0.324, "a": 6.0,
             "maxseg": MAXSEG, "minseg": MAXSEG / 4, "rtol": 10.0, "rann": 10.0,
             "nextdt": 1e-11, "maxdt": 1e-9}
    calforce = CalForce(force_mode=FORCE, state=state, Ngrid=32, cell=net.cell)
    mobility = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    timeint = TimeIntegration(integrator="Trapezoid", state=state, force=calforce, mobility=mobility)
    collision = Collision(collision_mode="Retroactive", state=state) if COLL == "1" else None
    topology = Topology(topology_mode="TopologyParallel", state=state, force=calforce,
                        mobility=mobility) if TOPO == "1" else None
    remesh = Remesh(remesh_rule="LengthBased", state=state)
    sim = SimulateNetwork(calforce=calforce, mobility=mobility, timeint=timeint,
                          collision=collision, topology=topology, remesh=remesh, vis=None,
                          state=state, max_step=NSTEPS, loading_mode="stress",
                          applied_stress=np.zeros(6), print_freq=50, plot_freq=10**9,
                          write_freq=10**9, write_dir=OUT)
    sim.run(net, state)

    data = net.get_disnet(ExaDisNet).export_data()
    pos = data["nodes"]["positions"]
    nid = data["segs"]["nodeids"]
    nf, sf = pos.shape[0], len(nid)
    Lf = sum(np.linalg.norm(pos[int(a)] - pos[int(c)]) for a, c in nid)
    # connectivity -> detect junction nodes (degree >= 3)
    deg = {}
    for a, c in nid:
        deg[int(a)] = deg.get(int(a), 0) + 1
        deg[int(c)] = deg.get(int(c), 0) + 1
    n_junction = sum(1 for d in deg.values() if d >= 3)

    frac = Lf / L0 if L0 else 0
    if sf == 0:
        mech = "fully_deleted"
    elif n_junction > 0 and frac < 0.97:
        mech = "junction_or_partial_annihilation"
    elif frac < 0.85:
        mech = "partial_annihilation_residual"
    elif abs(frac - 1) < 0.05 and n_junction == 0:
        mech = "pass_through_or_no_reaction"
    else:
        mech = "reacted_other"

    obs = dict(protocol="binary_collinear_relax", force_mode=FORCE,
               topology=bool(TOPO == "1"), collision=bool(COLL == "1"),
               L_over_b=LSEG, phi_deg=PHI, gap=GAP,
               n_nodes_before=n0, n_nodes_after=int(nf),
               n_segs_before=s0, n_segs_after=int(sf),
               line_length_before=float(L0), line_length_after=float(Lf),
               length_fraction=float(frac), n_junction_nodes=int(n_junction),
               segment_deleted=bool(sf < s0 * 0.5),
               mechanism=mech)
    json.dump(obs, open(os.path.join(OUT, "observer.json"), "w"), indent=1)
    print("OBSERVER:", json.dumps(obs), flush=True)
    pyexadis.finalize()


if __name__ == "__main__":
    main()

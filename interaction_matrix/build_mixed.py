"""Phase 2b: direct validation of the STEM prediction. Build a cubic pinned
forest distributed over the STEM network's actual slip-system inventory (the
4 systems at their measured populations), probe on the network's dominant
system, and measure the effective forest-hardening alpha with the SAME
DDD_FFT method used for the interaction matrix. Compare to the matrix
prediction (alpha_network ~= 0.69) and the uniform-forest macro (0.43).

  MSYS via dominant STEM system; NFOREST/SEED/density from env.
"""
import os, sys, json
import numpy as np

EX = os.path.expanduser("~/BO/exadis_src/python")
if EX not in sys.path:
    sys.path.append(EX)
import pyexadis
from pyexadis_base import (ExaDisNet, DisNetManager, SimulateNetworkPerf, CalForce,
                           MobilityLaw, TimeIntegration, Collision, Topology, Remesh,
                           NodeConstraints)
from pyexadis_utils import insert_infinite_line, dislocation_density

B_CU, MU = 2.55e-10, 54.6e9
# STEM network inventory (geometric assignment): (burg, plane, population)
STEM = [([0, 1, -1], [-1, 1, 1], 0.58),
        ([0, 1, 1], [1, -1, 1], 0.25),
        ([0, 1, 1], [1, 1, -1], 0.12),
        ([0, 1, -1], [1, 1, 1], 0.04)]
PROBE_B, PROBE_N = [0, 1, -1], [-1, 1, 1]      # dominant network system

LBOX = float(os.environ.get("LBOX", "5000"))
NFOREST = int(os.environ.get("NFOREST", "24"))
NPROBE = int(os.environ.get("NPROBE", "2"))
ERATE = float(os.environ.get("ERATE", "1e4"))
MAX_STRAIN = float(os.environ.get("MAX_STRAIN", "0.0012"))
FLOW_LO = float(os.environ.get("FLOW_LO", "0.0007"))
MAXSEG = float(os.environ.get("MAXSEG", "300"))
SEED = int(os.environ.get("SEED", "1234"))
OUT = os.environ.get("OUT", "mixed_out")


def main():
    pyexadis.initialize()
    os.makedirs(OUT, exist_ok=True)
    cell = pyexadis.Cell(LBOX)
    h = np.array(cell.h); origin = np.array(cell.origin)
    np.random.seed(SEED)
    nodes, segs = [], []
    franges = []
    # distribute NFOREST forest lines over the STEM systems by population
    counts = np.random.multinomial(NFOREST, [s[2] for s in STEM])
    for (b, n, _), c in zip(STEM, counts):
        for k in range(c):
            pos = origin + np.random.rand(3) @ h.T
            i0 = len(nodes)
            insert_infinite_line(cell, nodes, segs, np.array(b, float), np.array(n, float),
                                 pos, theta=float(np.random.choice([0, 30, 60, 90])), maxseg=MAXSEG)
            franges.append((i0, len(nodes)))
    for k in range(NPROBE):
        pos = origin + np.random.rand(3) @ h.T
        insert_infinite_line(cell, nodes, segs, np.array(PROBE_B, float),
                             np.array(PROBE_N, float), pos, theta=0.0, maxseg=MAXSEG)

    nodes = np.array(nodes); segs = np.array(segs)
    for a, b in franges:
        nodes[a:b, 3] = int(NodeConstraints.PINNED_NODE)
    net = DisNetManager(ExaDisNet(cell, nodes, segs))
    rho_total = dislocation_density(net, B_CU)
    pin = nodes[:, 3] == int(NodeConstraints.PINNED_NODE)
    tot = pinl = 0.0
    for s in segs:
        a, b = int(s[0]), int(s[1])
        L = np.linalg.norm(nodes[a, :3] - nodes[b, :3])
        tot += L
        if pin[a] and pin[b]:
            pinl += L
    rho_f = rho_total * (pinl / tot if tot else 0)

    bm = np.array(PROBE_B, float); nm = np.array(PROBE_N, float)
    em = bm / np.linalg.norm(bm) + nm / np.linalg.norm(nm); em = em / np.linalg.norm(em)
    schmid = abs(np.dot(em, bm / np.linalg.norm(bm)) * np.dot(em, nm / np.linalg.norm(nm)))

    state = {"crystal": "fcc", "burgmag": B_CU, "mu": MU, "nu": 0.324, "a": 6.0,
             "maxseg": MAXSEG, "minseg": MAXSEG / 4, "rtol": 10.0, "rann": 10.0,
             "nextdt": 1e-10, "maxdt": 1e-9}
    calforce = CalForce(force_mode="DDD_FFT_MODEL", state=state, Ngrid=32, cell=net.cell)
    mobility = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    timeint = TimeIntegration(integrator="Trapezoid", state=state, force=calforce, mobility=mobility)
    collision = Collision(collision_mode="Retroactive", state=state)
    topology = Topology(topology_mode="TopologyParallel", state=state, force=calforce, mobility=mobility)
    remesh = Remesh(remesh_rule="LengthBased", state=state)
    sim = SimulateNetworkPerf(calforce=calforce, mobility=mobility, timeint=timeint,
                              collision=collision, topology=topology, remesh=remesh,
                              cross_slip=None, vis=None, loading_mode="strain_rate",
                              erate=ERATE, edir=em, max_strain=MAX_STRAIN,
                              burgmag=state["burgmag"], state=state, print_freq=50,
                              plot_freq=10**9, write_freq=10**9, write_dir=OUT)
    sim.run(net, state)
    dat = np.loadtxt(os.path.join(OUT, "stress_strain_dens.dat"), ndmin=2)
    strain, stress = dat[:, 1], dat[:, 2]
    m = strain >= FLOW_LO
    axial = float(np.mean(np.abs(stress[m]))) if m.any() else float(np.abs(stress[-1]))
    tau = schmid * axial
    alpha = tau / (MU * B_CU * np.sqrt(rho_f)) if rho_f > 0 else 0.0
    res = dict(kind="stem_inventory_forest", rho_f=float(rho_f), schmid=float(schmid),
               axial_flow=axial, tau_resolved=float(tau), alpha=float(alpha),
               n_forest=NFOREST, seed=SEED, counts=[int(c) for c in counts])
    json.dump(res, open(os.path.join(OUT, "mixed_result.json"), "w"), indent=1)
    print(f"  STEM-inventory forest: rho_f={rho_f:.2e} tau={tau/1e6:.1f} MPa alpha={alpha:.3f} "
          f"(predicted 0.69; uniform 0.43)", flush=True)
    pyexadis.finalize()


if __name__ == "__main__":
    main()

"""DEFINITIVE TEST: does the REAL STEM-reconstructed Cu dislocation config
WORK-HARDEN under full FCC DDD (with junctions) — the thing PFC fundamentally
could not do?

Adapts ExaDiS's canonical fcc_Cu hardening setup (SUBCYCLING N-body elastic
force + TopologyParallel JUNCTIONS + FCC_0 mobility + strain_rate loading ->
stress-strain curve) to seed from the STEM-reconstructed network instead of a
synthetic config. If the flow stress RISES with strain (work hardening), DDD
from the reconstructed config hardens (PFC could not).

Run on HPC (pyexadis built):
  PYTHONPATH=<exadis>/python python3 run_stem_hardening.py
"""

import os, sys, json
import numpy as np

EX = os.path.expanduser("~/BO/exadis_src/python")
if EX not in sys.path:
    sys.path.append(EX)
import pyexadis
from pyexadis_base import (ExaDisNet, DisNetManager, SimulateNetworkPerf,
                           CalForce, MobilityLaw, TimeIntegration, Collision,
                           Topology, Remesh, NodeConstraints)

NET = os.environ.get("STEM_NET", "stem_network.json")
MAX_STRAIN = float(os.environ.get("MAX_STRAIN", "0.0015"))
ERATE = float(os.environ.get("ERATE", "1e4"))
ZBOX = float(os.environ.get("ZBOX", "5.0"))  # thicken thin foil box for FFT
OUT = os.environ.get("OUT", "stem_hardening_out")


def main():
    net_j = json.loads(open(NET).read())
    pyexadis.initialize()
    state = {"crystal": "fcc", "burgmag": 2.556e-10, "mu": 54.6e9,
             "nu": 0.324, "a": 6.0,
             "maxseg": 300.0, "minseg": 50.0, "rtol": 10.0, "rann": 10.0,
             "nextdt": 1e-10, "maxdt": 1e-9}
    cmap = {"PINNED_NODE": NodeConstraints.PINNED_NODE,
            "UNCONSTRAINED": NodeConstraints.UNCONSTRAINED}
    nodes = np.asarray([[float(r[0]), float(r[1]), float(r[2]), cmap[r[3]]]
                        for r in net_j["nodes"]])
    segs = np.asarray(net_j["segs"], dtype=float)
    # fully periodic box for the N-body FFT force (standard hardening setup)
    h = np.asarray(net_j["cell"]["h_b"], dtype=float)
    h[2, 2] *= ZBOX   # thicken z so FFT voxels are not high-aspect
    cell = pyexadis.Cell(h=h, is_periodic=[True, True, True])
    net = DisNetManager(ExaDisNet(cell, nodes, segs))

    calforce = CalForce(force_mode="SUBCYCLING_MODEL", state=state,
                        Ngrid=32, cell=net.cell)
    mobility = MobilityLaw(mobility_law="FCC_0", state=state,
                           Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    timeint = TimeIntegration(integrator="Subcycling",
                              rgroups=[0.0, 50.0, 200.0, 600.0], state=state,
                              force=calforce, mobility=mobility)
    collision = Collision(collision_mode="Retroactive", state=state)
    topology = Topology(topology_mode="TopologyParallel", state=state,
                        force=calforce, mobility=mobility)   # JUNCTIONS ON
    remesh = Remesh(remesh_rule="LengthBased", state=state)
    os.makedirs(OUT, exist_ok=True)
    sim = SimulateNetworkPerf(calforce=calforce, mobility=mobility,
                              timeint=timeint, collision=collision,
                              topology=topology, remesh=remesh, cross_slip=None,
                              vis=None, loading_mode="strain_rate", erate=ERATE,
                              edir=np.array([0., 0., 1.]),
                              max_strain=MAX_STRAIN, burgmag=state["burgmag"],
                              state=state, print_freq=10, plot_freq=10**9,
                              write_freq=10**9, write_dir=OUT)
    sim.run(net, state)
    pyexadis.finalize()
    print("STEM HARDENING RUN DONE -> stress-strain in", OUT)


if __name__ == "__main__":
    main()

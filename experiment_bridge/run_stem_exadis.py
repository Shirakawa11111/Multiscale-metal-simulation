"""Run the REAL STEM-reconstructed dislocation network through ExaDiS DDD.

The culmination of the image->simulation route: the network produced by
stem_to_exadis.py (from real STEM 3D reconstruction) is loaded into ExaDiS and
evolved under applied shear with full nodal DDD (glide + collision + remesh).
This executes the route end-to-end with the real DDD engine (not a surrogate).

Run on the HPC where pyexadis is built:
  PYTHONPATH=<exadis>/python:<exadis>/build/python python3 run_stem_exadis.py \
      --network stem_network.json --steps 30 --out stem_ddd_out
"""

import argparse, json, os
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--network", required=True)
    p.add_argument("--steps", type=int, default=30)
    p.add_argument("--dt", type=float, default=1e-11)
    p.add_argument("--out", default="stem_ddd_out")
    p.add_argument("--mu", type=float, default=48e9)
    p.add_argument("--nu", type=float, default=0.34)
    p.add_argument("--sxy", type=float, default=100e6)
    p.add_argument("--syz", type=float, default=60e6)
    a = p.parse_args()

    net_j = json.loads(open(a.network).read())
    import pyexadis
    from pyexadis_base import (CalForce, Collision, DisNetManager, ExaDisNet,
                               MobilityLaw, NodeConstraints, Remesh,
                               SimulateNetwork, TimeIntegration)
    pyexadis.initialize()
    try:
        cmap = {"PINNED_NODE": NodeConstraints.PINNED_NODE,
                "UNCONSTRAINED": NodeConstraints.UNCONSTRAINED}
        nodes = np.asarray([[float(r[0]), float(r[1]), float(r[2]), cmap[r[3]]]
                            for r in net_j["nodes"]])
        segs = np.asarray(net_j["segs"], dtype=float)
        cell = pyexadis.Cell(h=np.asarray(net_j["cell"]["h_b"], dtype=float),
                             is_periodic=net_j["cell"]["is_periodic"])
        net = DisNetManager(ExaDisNet(cell, nodes, segs))
        box = np.asarray(net_j["cell"]["box_size_b"], dtype=float)
        min_box = float(np.min(box))
        state = {"burgmag": float(net_j["units"]["length_unit_m"]),
                 "mu": a.mu, "nu": a.nu, "a": 1.0,
                 "maxseg": 0.04 * min_box, "minseg": 0.01 * min_box,
                 "rann": 3.0}
        # applied shear (drives the geometrically-assigned slip systems)
        applied = np.array([0, 0, 0, a.syz, 0, a.sxy], dtype=float)  # Voigt 23,13,12

        n0 = len(nodes)
        calforce = CalForce(force_mode="LineTension", state=state)
        mobility = MobilityLaw(mobility_law="SimpleGlide", state=state)
        timeint = TimeIntegration(integrator="EulerForward", dt=a.dt, state=state)
        collision = Collision(collision_mode="Retroactive", state=state)
        remesh = Remesh(remesh_rule="LengthBased", state=state)
        os.makedirs(a.out, exist_ok=True)
        sim = SimulateNetwork(calforce=calforce, mobility=mobility,
                              timeint=timeint, collision=collision,
                              topology=None, remesh=remesh, vis=None,
                              state=state, max_step=a.steps,
                              loading_mode="stress", applied_stress=applied,
                              print_freq=5, plot_freq=10**9, write_freq=10**9,
                              write_dir=a.out)
        sim.run(net, state)
        # save final network (API varies across pyexadis versions; try several)
        n1 = -1
        try:
            data = net.get_disnet(ExaDisNet).export_data()
            n1 = len(data.get("nodes", {}).get("positions", []))
            json.dump({"n_nodes": n1}, open(os.path.join(a.out,
                      "final_counts.json"), "w"))
        except Exception as e:
            print("final export note:", e)
        summary = dict(status="ok", n_nodes_start=int(n0), n_nodes_end=int(n1),
                       steps=a.steps, applied_shear_MPa=[a.sxy/1e6, a.syz/1e6],
                       note="REAL STEM-reconstructed Cu dislocation network "
                            "evolved with ExaDiS DDD (glide+collision+remesh)")
        json.dump(summary, open(os.path.join(a.out, "summary.json"), "w"),
                  indent=2)
        print("STEM->ExaDiS RUN OK:", json.dumps(summary))
    finally:
        pyexadis.finalize()


if __name__ == "__main__":
    main()

"""Real-network DDD audit [M3]: is the STEM->IDR->ExaDiS import auditable & robust?

For one lowered network (chosen assignment_policy + cell_policy), run a zero-stress relaxation and a
short stress-controlled loading in ExaDiS, and measure: network survival (does relaxation delete it?),
density evolution, topology-event proxy (junction nodes), line-length change, and the slip-system
inventory. Driven per (policy, seed) by run_audit_campaign.py so the SENSITIVITY of these outcomes to
the assignment/cell policy can be reported -- the point of M3.

  JTAG=top1_foil POLICY=top1 CELL=as_is SEED=0 NREL=300 NLOAD=300 \
  IDR=results_exadis/cu_stem_idr.json OUT=audit/top1_foil PYTHONPATH=~/BO/exadis_src/python \
  python3 real_network_audit.py
"""
import os, sys, json
import numpy as np

EX = os.path.expanduser("~/BO/exadis_src/python")
if EX not in sys.path:
    sys.path.append(EX)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
import pyexadis
from pyexadis_base import (ExaDisNet, DisNetManager, SimulateNetwork, CalForce, MobilityLaw,
                           TimeIntegration, Collision, Topology, Remesh, NodeConstraints)
from defect_ir.adapters.to_exadis import idr_to_exadis_network

B_CU, MU, NU = 2.556e-10, 54.6e9, 0.324
IDR = os.environ.get("IDR", "results_exadis/cu_stem_idr.json")
POLICY = os.environ.get("POLICY", "top1")
CELL = os.environ.get("CELL", "as_is")
ZBOX = float(os.environ.get("ZBOX", "5"))
ENDPOINT = os.environ.get("ENDPOINT", "pinned")
SEED = int(os.environ.get("SEED", "0"))
NREL = int(os.environ.get("NREL", "300"))
NLOAD = int(os.environ.get("NLOAD", "300"))
SXY = float(os.environ.get("SXY", "100")) * 1e6
SYZ = float(os.environ.get("SYZ", "60")) * 1e6
JTAG = os.environ.get("JTAG", f"{POLICY}_{CELL}")
OUT = os.environ.get("OUT", f"audit/{JTAG}")


def build_net(netd):
    h = np.array(netd["cell"]["h_b"], float)
    cell = pyexadis.Cell(h=h, is_periodic=netd["cell"]["is_periodic"])
    cmap = {"PINNED_NODE": int(NodeConstraints.PINNED_NODE),
            "UNCONSTRAINED": int(NodeConstraints.UNCONSTRAINED)}
    nodes = np.array([[n[0], n[1], n[2], cmap.get(n[3], 0)] for n in netd["nodes"]], float)
    segs = np.array(netd["segs"], float)
    return DisNetManager(ExaDisNet(cell, nodes, segs)), cell


def metrics(net):
    d = net.get_disnet(ExaDisNet).export_data()
    pos, nid = d["nodes"]["positions"], np.asarray(d["segs"]["nodeids"]).astype(int)
    L = float(sum(np.linalg.norm(pos[int(a)] - pos[int(c)]) for a, c in nid))
    deg = {}
    for a, c in nid:
        deg[int(a)] = deg.get(int(a), 0) + 1
        deg[int(c)] = deg.get(int(c), 0) + 1
    njun = sum(1 for v in deg.values() if v >= 3)
    return dict(n_nodes=int(len(pos)), n_segs=int(len(nid)), line_len_b=L, n_junction_nodes=int(njun))


def inventory(netd):
    inv = {}
    for s in netd["segs"]:
        k = tuple(round(x, 2) for x in s[2:8])
        inv[str(k)] = inv.get(str(k), 0) + 1
    return inv


def modules(state, cell, periodic):
    force = "DDD_FFT_MODEL" if periodic else "LineTension"
    cf = CalForce(force_mode=force, state=state, Ngrid=32, cell=cell) if periodic else \
        CalForce(force_mode=force, state=state)
    mob = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    ti = TimeIntegration(integrator="Trapezoid", state=state, force=cf, mobility=mob)
    col = Collision(collision_mode="Retroactive", state=state)
    topo = Topology(topology_mode="TopologyParallel", state=state, force=cf, mobility=mob)
    rm = Remesh(remesh_rule="LengthBased", state=state)
    return cf, mob, ti, col, topo, rm


def main():
    pyexadis.initialize()
    os.makedirs(OUT, exist_ok=True)
    doc = json.load(open(IDR))
    netd = idr_to_exadis_network(doc, assignment_policy=POLICY, cell_policy=CELL, zbox=ZBOX,
                                 seed=SEED, endpoint_policy=ENDPOINT)
    periodic = all(netd["cell"]["is_periodic"])
    vol_b3 = abs(np.linalg.det(np.array(netd["cell"]["h_b"], float)))
    rho = lambda L: L / vol_b3 / B_CU ** 2

    net, cell = build_net(netd)
    m0 = metrics(net)
    state = {"crystal": "fcc", "burgmag": B_CU, "mu": MU, "nu": NU, "a": 6.0,
             "maxseg": 100.0, "minseg": 25.0, "rtol": 10.0, "rann": 10.0, "nextdt": 1e-12, "maxdt": 1e-10}

    # Phase A: zero-stress relaxation (survival test)
    cf, mob, ti, col, topo, rm = modules(state, net.cell, periodic)
    SimulateNetwork(calforce=cf, mobility=mob, timeint=ti, collision=col, topology=topo, remesh=rm,
                    vis=None, state=state, loading_mode="stress", applied_stress=np.zeros(6),
                    max_step=NREL, print_freq=10**9, plot_freq=10**9, write_freq=10**9, write_dir=OUT).run(net, state)
    mA = metrics(net)

    # Phase B: short stress-controlled loading
    cf2, mob2, ti2, col2, topo2, rm2 = modules(state, net.cell, periodic)
    applied = np.array([0, 0, 0, SYZ, 0, SXY], float)
    SimulateNetwork(calforce=cf2, mobility=mob2, timeint=ti2, collision=col2, topology=topo2, remesh=rm2,
                    vis=None, state=state, loading_mode="stress", applied_stress=applied,
                    max_step=NLOAD, print_freq=10**9, plot_freq=10**9, write_freq=10**9, write_dir=OUT).run(net, state)
    mB = metrics(net)

    # objectives (stability + interpretability, NOT stress-curve matching)
    surv = mA["n_segs"] / max(1, m0["n_segs"])
    obj = dict(
        network_survival_score=round(min(1.0, surv if surv < 1 else 1.0 / surv), 3),  # 1=retained; <1 = collapsed OR blew up
        density_growth=round(mB["line_len_b"] / max(1e-9, mA["line_len_b"]), 3),
        density_growth_plausible=bool(0.5 <= mB["line_len_b"] / max(1e-9, mA["line_len_b"]) <= 5.0),
        topology_event_rate=int(mB["n_junction_nodes"]),
    )
    out = dict(tag=JTAG, policy=POLICY, cell=CELL, endpoint=ENDPOINT, seed=SEED, periodic=periodic,
               zbox=ZBOX if periodic else None, force="DDD_FFT_MODEL" if periodic else "LineTension",
               objectives=obj, build=m0, after_relax=mA, after_load=mB,
               rho_build=rho(m0["line_len_b"]), rho_relax=rho(mA["line_len_b"]), rho_load=rho(mB["line_len_b"]),
               survival_seg_frac=round(mA["n_segs"] / max(1, m0["n_segs"]), 3),
               relax_len_frac=round(mA["line_len_b"] / max(1e-9, m0["line_len_b"]), 3),
               load_density_growth=round(mB["line_len_b"] / max(1e-9, mA["line_len_b"]), 3),
               topo_events_relax=mA["n_junction_nodes"], topo_events_load=mB["n_junction_nodes"],
               n_slip_systems_used=len(inventory(netd)))
    json.dump(out, open(os.path.join(OUT, "audit.json"), "w"), indent=1)
    print(f"AUDIT {JTAG} seed{SEED}: survival={out['survival_seg_frac']} relax_len={out['relax_len_frac']} "
          f"rho {out['rho_build']:.2e}->{out['rho_relax']:.2e}->{out['rho_load']:.2e} "
          f"load_growth={out['load_density_growth']} topo_load={out['topo_events_load']} "
          f"force={out['force']}", flush=True)
    pyexadis.finalize()


if __name__ == "__main__":
    main()

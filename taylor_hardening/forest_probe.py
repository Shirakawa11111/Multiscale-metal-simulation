"""Forest-probe Taylor measurement (Madec/Devincre style) — the RIGOROUS test.

The initial-density series failed to give a sqrt(rho) lever because free
dislocations ANNIHILATE down to a common Kocks-Mecking steady state before the
flow window (all runs end at the same rho_ss). The fix: PIN the forest.

Construction (no fragile manual crystallography — pin a subset of a
generate_line_config network):
  1. generate_line_config(num_lines) -> commensurate FCC lines on random {111}.
  2. find connected components (= individual dislocation lines) from the segs.
  3. PIN every node of all-but-K lines  (the immobile forest at density rho_f).
     The remaining K lines stay mobile (the probes that carry plastic flow).
  4. strain-rate load with junctions on; the probes glide and get pinned by
     forest junctions; the flow-stress plateau = tau_c = alpha*mu*b*sqrt(rho_f).

Because the forest is PINNED it cannot annihilate, so rho stays ~ rho_f and the
density lever across runs (different num_lines) SURVIVES -> a real Taylor series.

  OMP_NUM_THREADS=48 NUM_LINES=<N> PYTHONPATH=~/BO/exadis_src/python \
      python3 forest_probe.py
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
from pyexadis_utils import dislocation_density

LBOX = float(os.environ.get("LBOX", "10000"))
NUM_LINES = int(os.environ.get("NUM_LINES", "40"))
N_PROBE = int(os.environ.get("N_PROBE", "0"))          # 0 -> use PROBE_FRAC
PROBE_FRAC = float(os.environ.get("PROBE_FRAC", "0.25"))
ERATE = float(os.environ.get("ERATE", "3e4"))
MAX_STRAIN = float(os.environ.get("MAX_STRAIN", "0.003"))
FLOW_LO = float(os.environ.get("FLOW_LO", "0.0015"))
NGRID = int(os.environ.get("NGRID", "32"))
MAXSEG = float(os.environ.get("MAXSEG", "600"))
MINSEG = float(os.environ.get("MINSEG", "150"))
OUT = os.environ.get("OUT", "forest_out")
MU, B_CU = 54.6e9, 2.55e-10


def connected_components(n_nodes, seg_nodeids):
    """Union-find over segment endpoints -> list of node-index sets (lines)."""
    parent = list(range(n_nodes))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for s in seg_nodeids:
        a, b = int(s[0]), int(s[1])
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comps = {}
    for i in range(n_nodes):
        comps.setdefault(find(i), []).append(i)
    return list(comps.values())


def main():
    pyexadis.initialize()
    os.makedirs(OUT, exist_ok=True)
    wdir = os.path.join(OUT, f"n{NUM_LINES}_s{os.environ.get('SEED','1234')}")
    os.makedirs(wdir, exist_ok=True)
    state = {"crystal": "fcc", "burgmag": B_CU, "mu": MU, "nu": 0.324, "a": 6.0,
             "maxseg": MAXSEG, "minseg": MINSEG, "rtol": 10.0, "rann": 10.0,
             "nextdt": 1e-10, "maxdt": 1e-9}

    SEED = int(os.environ.get("SEED", "1234"))
    G = ExaDisNet()
    G.generate_line_config("fcc", LBOX, NUM_LINES, maxseg=MAXSEG, seed=SEED)
    data = G.export_data()
    nodeids = data["segs"]["nodeids"]
    n_nodes = data["nodes"]["positions"].shape[0]
    comps = connected_components(n_nodes, nodeids)
    n_lines = len(comps)
    k_probe = N_PROBE if N_PROBE > 0 else max(1, int(round(PROBE_FRAC * n_lines)))
    # smallest components are likely artifacts; keep the LARGEST lines as probes
    comps_sorted = sorted(comps, key=len, reverse=True)
    probe_nodes = set()
    for comp in comps_sorted[:k_probe]:
        probe_nodes.update(comp)

    cons = data["nodes"]["constraints"].copy()
    pinned = 0
    for i in range(n_nodes):
        if i not in probe_nodes:
            cons[i, 0] = int(NodeConstraints.PINNED_NODE)
            pinned += 1
    data["nodes"]["constraints"] = cons

    G2 = ExaDisNet()
    G2.import_data(data)
    net = DisNetManager(G2)
    rho_total = dislocation_density(net, state["burgmag"])
    # forest density by pinned LINE LENGTH (not node fraction): sum lengths of
    # segments whose both endpoints are pinned, as a fraction of total length
    pos = data["nodes"]["positions"]
    pin_mask = (data["nodes"]["constraints"][:, 0] == int(NodeConstraints.PINNED_NODE))
    tot_len = pinned_len = 0.0
    for s in nodeids:
        a, b = int(s[0]), int(s[1])
        L = float(np.linalg.norm(pos[a] - pos[b]))
        tot_len += L
        if pin_mask[a] and pin_mask[b]:
            pinned_len += L
    len_frac = pinned_len / tot_len if tot_len > 0 else pinned / n_nodes
    rho_forest = rho_total * len_frac          # line-length forest density
    rho_forest_nodefrac = rho_total * pinned / n_nodes   # old proxy, for comparison

    calforce = CalForce(force_mode="DDD_FFT_MODEL", state=state, Ngrid=NGRID,
                        cell=net.cell)
    mobility = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0,
                           Mscrew=64103.0, vmax=4000.0)
    timeint = TimeIntegration(integrator="Trapezoid", state=state,
                              force=calforce, mobility=mobility)
    collision = Collision(collision_mode="Retroactive", state=state)
    topology = Topology(topology_mode="TopologyParallel", state=state,
                        force=calforce, mobility=mobility)
    remesh = Remesh(remesh_rule="LengthBased", state=state)
    sim = SimulateNetworkPerf(calforce=calforce, mobility=mobility,
                              timeint=timeint, collision=collision,
                              topology=topology, remesh=remesh, cross_slip=None,
                              vis=None, loading_mode="strain_rate", erate=ERATE,
                              edir=np.array([0., 0., 1.]), max_strain=MAX_STRAIN,
                              burgmag=state["burgmag"], state=state,
                              print_freq=20, plot_freq=10**9, write_freq=10**9,
                              write_dir=wdir)
    sim.run(net, state)

    dat = os.path.join(wdir, "stress_strain_dens.dat")
    res = np.loadtxt(dat, ndmin=2)
    strain, stress, dens = res[:, 1], res[:, 2], res[:, 3]
    mask = strain >= FLOW_LO
    flow = float(np.mean(np.abs(stress[mask]))) if mask.any() \
        else float(np.abs(stress[-1]))
    flow_std = float(np.std(np.abs(stress[mask]))) if mask.any() else 0.0
    rho_flow = float(np.mean(dens[mask])) if mask.any() else rho_total
    # depinning proxy: peak stress (the stress at which the probes first break
    # through the forest, before steady serrated flow)
    tau_peak = float(np.max(np.abs(stress)))
    alpha_pt = flow / (MU * B_CU * np.sqrt(rho_forest)) if rho_forest > 0 else 0.0

    out = dict(num_lines=NUM_LINES, n_lines=int(n_lines), k_probe=int(k_probe),
               seed=SEED, n_nodes=int(n_nodes), n_pinned=int(pinned),
               rho_total=float(rho_total), rho_forest=float(rho_forest),
               rho_forest_nodefrac=float(rho_forest_nodefrac),
               rho_flow=rho_flow, flow_stress=flow, flow_std=flow_std,
               tau_peak=tau_peak, alpha_point=float(alpha_pt), erate=ERATE,
               max_strain=MAX_STRAIN, curve=res.tolist())
    with open(os.path.join(wdir, "result.json"), "w") as f:
        json.dump(out, f)
    print(f"  num_lines={NUM_LINES} n_lines={n_lines} k_probe={k_probe} "
          f"rho_forest={rho_forest:.3e} rho_flow={rho_flow:.3e} "
          f"flow={flow/1e6:.1f}+/-{flow_std/1e6:.1f} MPa alpha_pt={alpha_pt:.3f}",
          flush=True)
    pyexadis.finalize()


if __name__ == "__main__":
    main()

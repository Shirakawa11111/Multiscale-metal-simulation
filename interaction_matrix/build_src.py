"""Phase 1b: controlled single-(mobile, forest)-system DDD configuration to
measure an FCC interaction-matrix coefficient a_ij per junction type.

Builds a forest of pinned, commensurate dislocations all on ONE slip system f
(so the forest is a VERIFIED forest of a known junction type, not a random
pinned subset), plus a free probe on slip system m. Tension is applied along the
axis that maximizes the Schmid factor on the mobile system m (e ∝ b_m + n_m ->
Schmid 0.5), so only the probe glides (forest is pinned) and the flow plateau
measures the depinning resistance of the (m,f) junction type:
        tau_c = alpha_mf * mu * b * sqrt(rho_f),   a_mf = alpha_mf^2.

Run on HPC: MSYS=<m> FSYS=<f> PYTHONPATH=~/BO/exadis_src/python python3 build_pair.py
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

# ExaDiS FCC 12 slip systems (integer b, n), matching generate_line_config order
BI = np.array([[0, 1, -1], [1, 0, -1], [1, -1, 0], [0, 1, -1], [1, 0, 1], [1, 1, 0],
               [0, 1, 1], [1, 0, -1], [1, 1, 0], [0, 1, 1], [1, 0, 1], [1, -1, 0]], float)
NI = np.array([[1, 1, 1]] * 3 + [[-1, 1, 1]] * 3 + [[1, -1, 1]] * 3 + [[1, 1, -1]] * 3, float)

B_CU, MU = 2.55e-10, 54.6e9
MSYS = int(os.environ.get("MSYS", "0"))
FSYS = int(os.environ.get("FSYS", "8"))
LBOX = float(os.environ.get("LBOX", "8000"))
NFOREST = int(os.environ.get("NFOREST", "16"))
NPROBE = int(os.environ.get("NPROBE", "2"))
ERATE = float(os.environ.get("ERATE", "2e3"))
MAX_STRAIN = float(os.environ.get("MAX_STRAIN", "0.001"))
FLOW_LO = float(os.environ.get("FLOW_LO", "0.0006"))
MAXSEG = float(os.environ.get("MAXSEG", "400"))
SEED = int(os.environ.get("SEED", "1234"))
OUT = os.environ.get("OUT", "pair_out")


def jtype(i, j):
    b1, n1, b2, n2 = BI[i], NI[i], BI[j], NI[j]
    par = lambda u, v: np.all(np.cross(u, v) == 0)
    if par(b1, b2) and par(n1, n2):
        return "self"
    if par(n1, n2):
        return "coplanar"
    if par(b1, b2):
        return "collinear"
    if int(round(np.dot(b1, b2))) == 0:
        return "Hirth"
    for s in (1, -1):
        b3 = b1 + s * b2
        if np.all(b3 == 0):
            continue
        if sorted(np.abs(b3)) == [0, 1, 1] and (np.dot(b3, n1) == 0 or np.dot(b3, n2) == 0):
            return "glissile"
    return "Lomer"


def main():
    pyexadis.initialize()
    os.makedirs(OUT, exist_ok=True)
    cell = pyexadis.Cell(LBOX)
    h = np.array(cell.h); origin = np.array(cell.origin)
    np.random.seed(SEED)
    nodes, segs = [], []
    forest_ranges = []
    bf, nf = BI[FSYS], NI[FSYS]
    for k in range(NFOREST):
        pos = origin + np.random.rand(3) @ h.T
        i0 = len(nodes)
        insert_infinite_line(cell, nodes, segs, bf, nf, pos,
                             theta=float(np.random.choice([0, 30, 60, 90])), maxseg=MAXSEG)
        forest_ranges.append((i0, len(nodes)))
    n_forest_nodes = len(nodes)
    bm, nm = BI[MSYS], NI[MSYS]
    # PROBE = pinned-end finite bowing segment (Frank-Read-like) so it cannot
    # simply annihilate and vanish; collinear reaction traps the remaining arms.
    xi = np.cross(nm, bm); xi = xi/np.linalg.norm(xi)        # edge line dir in plane
    Lseg = 0.5*LBOX
    nseg = max(4, int(Lseg/MAXSEG))
    for k in range(NPROBE):
        c = origin + (0.2+0.6*np.random.rand(3)) @ h.T
        i0 = len(nodes)
        for j in range(nseg+1):
            p = c + (j/nseg - 0.5)*Lseg*xi
            con = int(NodeConstraints.PINNED_NODE) if j in (0, nseg) else int(NodeConstraints.UNCONSTRAINED)
            nodes.append(np.concatenate((p, [con])))
        for j in range(nseg):
            segs.append(np.concatenate(([i0+j, i0+j+1], bm, nm)))

    nodes = np.array(nodes); segs = np.array(segs)
    for a, b in forest_ranges:
        nodes[a:b, 3] = int(NodeConstraints.PINNED_NODE)   # pin the whole forest
    net = DisNetManager(ExaDisNet(cell, nodes, segs))
    rho_total = dislocation_density(net, B_CU)
    # forest density by pinned line length
    pin = nodes[:, 3] == int(NodeConstraints.PINNED_NODE)
    tot = pinl = 0.0
    for s in segs:
        a, b = int(s[0]), int(s[1])
        L = np.linalg.norm(nodes[a, :3] - nodes[b, :3])
        tot += L
        if pin[a] and pin[b]:
            pinl += L
    rho_f = rho_total * (pinl / tot if tot else 0)

    # tension axis maximizing Schmid on mobile system m (Schmid -> 0.5)
    em = bm / np.linalg.norm(bm) + nm / np.linalg.norm(nm)
    em = em / np.linalg.norm(em)
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
    axial_flow = float(np.mean(np.abs(stress[m]))) if m.any() else float(np.abs(stress[-1]))
    tau_resolved = schmid * axial_flow           # resolved shear on system m
    alpha = tau_resolved / (MU * B_CU * np.sqrt(rho_f)) if rho_f > 0 else 0.0
    axial_peak = float(np.max(np.abs(stress)))   # critical/depinning threshold
    tau_peak = schmid * axial_peak
    alpha_peak = tau_peak / (MU * B_CU * np.sqrt(rho_f)) if rho_f > 0 else 0.0
    res = dict(msys=MSYS, fsys=FSYS, jtype=jtype(MSYS, FSYS), schmid=float(schmid),
               rho_f=float(rho_f), rho_total=float(rho_total), n_forest=NFOREST,
               axial_flow=axial_flow, tau_resolved=float(tau_resolved),
               alpha=float(alpha), a_ij=float(alpha ** 2), erate=ERATE, seed=SEED,
               n_probe=NPROBE, alpha_peak=float(alpha_peak), tau_peak=float(tau_peak))
    json.dump(res, open(os.path.join(OUT, "pair_result.json"), "w"), indent=1)
    print(f"  ({MSYS},{FSYS}) {jtype(MSYS,FSYS):>9}: schmid={schmid:.3f} rho_f={rho_f:.2e} "
          f"axial={axial_flow/1e6:.1f} tau={tau_resolved/1e6:.1f} MPa alpha={alpha:.3f} a_ij={alpha**2:.3f}",
          flush=True)
    pyexadis.finalize()


if __name__ == "__main__":
    main()

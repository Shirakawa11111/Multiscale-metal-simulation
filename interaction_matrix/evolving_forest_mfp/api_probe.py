"""Resolve the evolving_forest_mfp API assumptions on THIS ExaDiS build (cheap, no physics)."""
import os, sys, inspect
import numpy as np
EX = os.path.expanduser("~/BO/exadis_src/python")
if EX not in sys.path:
    sys.path.append(EX)
import pyexadis
from pyexadis_base import (ExaDisNet, DisNetManager, SimulateNetwork, CalForce,
                           MobilityLaw, TimeIntegration, Collision, Topology, Remesh,
                           NodeConstraints)

pyexadis.initialize()
B_CU, MU = 2.55e-10, 54.6e9
LBOX = 6000.0
cell = pyexadis.Cell(h=LBOX*np.eye(3), is_periodic=[True, True, True])
C = np.array(cell.center())
b = np.array([0., 1., -1.]); n = np.array([1., 1., 1.])
bh = b/np.linalg.norm(b); e = np.cross(n, b); e = e/np.linalg.norm(e)
xi = (np.cos(np.radians(45))*bh + np.sin(np.radians(45))*e); xi /= np.linalg.norm(xi)
nodes, segs = [], []
for j in range(6):
    p = C + (j/5-0.5)*1500*xi
    con = int(NodeConstraints.PINNED_NODE) if j in (0, 5) else int(NodeConstraints.UNCONSTRAINED)
    nodes.append(np.concatenate((p, [con])))
for j in range(5):
    segs.append(np.concatenate(([j, j+1], b, n/np.linalg.norm(n))))
nodes = np.array(nodes); segs = np.array(segs)
net = DisNetManager(ExaDisNet(cell, nodes, segs))

print("=== export_data structure ===")
d = net.get_disnet(ExaDisNet).export_data()
print("top keys:", list(d.keys()))
print("nodes keys:", list(d["nodes"].keys()))
print("segs keys:", list(d["segs"].keys()))
for k, v in d["segs"].items():
    try:
        print(f"  segs[{k!r}] shape={np.asarray(v).shape}")
    except Exception as ex:
        print(f"  segs[{k!r}] (unsizable: {ex})")
print("nodes has 'tags'?", "tags" in d["nodes"])
print("nodes has 'constraints'?", "constraints" in d["nodes"])

print("=== SimulateNetwork hooks ===")
for h in ("step_begin", "step_end", "step", "step_update_response", "run"):
    print(f"  SimulateNetwork.{h}: {'yes' if hasattr(SimulateNetwork, h) else 'NO'}")

print("=== ExaDisNet velocity accessor ===")
G = net.get_disnet(ExaDisNet)
print("  get_velocities:", "yes" if hasattr(G, "get_velocities") else "NO")
print("  ExaDisNet methods sample:", [m for m in dir(G) if not m.startswith("_")][:25])

print("=== cell.closest_image ===")
print("  cell.closest_image:", "yes" if hasattr(cell, "closest_image") else "NO")

print("=== insert_infinite_line / dislocation_density ===")
try:
    from pyexadis_utils import insert_infinite_line, dislocation_density
    print("  imported OK")
    print("  insert_infinite_line sig:", str(inspect.signature(insert_infinite_line)))
    print("  dislocation_density sig:", str(inspect.signature(dislocation_density)))
except Exception as ex:
    print("  IMPORT FAIL:", ex)

pyexadis.finalize()
print("PROBE DONE")

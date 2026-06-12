"""Export PFC density-field snapshots to .xyz for OVITO (3D: enables DXA
dislocation-line analysis; 2D fields export with z=0).

Usage: python3 export_xyz.py <snapshot.npz> [out.xyz]
Auto-detects 2D vs 3D from the psi array shape.
"""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np


def export(npz_path, out_path=None):
    d = np.load(npz_path)
    psi = d["psi"]
    out_path = out_path or npz_path.replace(".npz", ".xyz")
    if psi.ndim == 3:
        from pfc3d import PFC3D, find_peaks_3d
        m = PFC3D.load(npz_path)
        pts = find_peaks_3d(m.psi, m.dx, m.dy, m.dz)
        box = (m.lx, m.ly, m.lz)
        pts3 = pts
    else:
        from pfc2d import PFC2D
        from defect_analysis import find_peaks
        m = PFC2D.load(npz_path)
        pts = find_peaks(m.psi, m.dx, m.dy)
        box = (m.lx, m.ly, 1.0)
        pts3 = np.column_stack([pts, np.zeros(len(pts))])
    with open(out_path, "w") as f:
        f.write(f"{len(pts3)}\n")
        f.write(f'Lattice="{box[0]} 0 0 0 {box[1]} 0 0 0 {box[2]}" '
                f'Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for p in pts3:
            f.write(f"Cu {p[0]:.4f} {p[1]:.4f} {p[2]:.4f}\n")
    print(f"{out_path}: {len(pts3)} atoms, box {box}")
    return out_path


if __name__ == "__main__":
    export(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)

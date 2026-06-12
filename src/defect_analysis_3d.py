"""3D dislocation-line detection for PFC BCC density fields.

Strategy (no OVITO dependency): per-atom common-neighbor-style disorder via
coordination deviation from the perfect-BCC value, then cluster the
disordered atoms into connected components = dislocation cores / lines.

For BCC, the first+second coordination shell has 14 neighbors (8 at
sqrt(3)/2 a, 6 at a). A perfect-lattice atom sees 14 neighbors within
1.15 a; dislocation cores deviate. We flag atoms whose 14-neighbor shell is
geometrically distorted (RMS deviation of neighbor distances from the ideal
two-shell pattern) and cluster them by proximity into lines, reporting line
count, total disordered fraction, and per-cluster extent (a proxy for line
length).
"""

import numpy as np
from scipy.spatial import cKDTree


A_BCC = 2.0 * np.pi * np.sqrt(2.0)


def _tiled(pts, box):
    """27-image tiling for periodic neighbour search."""
    imgs, idx = [], []
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for sz in (-1, 0, 1):
                imgs.append(pts + np.array([sx, sy, sz]) * box)
                idx.append(np.arange(len(pts)))
    return np.vstack(imgs), np.concatenate(idx)


def disorder_metric(pts, box, a0=A_BCC):
    """Per-atom structural disorder: RMS gap between the sorted 14 nearest
    neighbour distances and the ideal BCC two-shell distances
    (8 x sqrt(3)/2 a, 6 x a), normalised by a0. Robust to global strain via
    per-atom rescaling to the atom's own mean first-shell distance."""
    if len(pts) < 15:
        return np.zeros(len(pts))
    tiled, src = _tiled(pts, box)
    tree = cKDTree(tiled)
    d, _ = tree.query(pts, k=15)        # self + 14
    d = d[:, 1:]                         # drop self
    ideal = np.array([np.sqrt(3) / 2] * 8 + [1.0] * 6) * a0
    # normalise each atom by its own scale (first 8 neighbours' mean / ideal)
    scale = d[:, :8].mean(axis=1) / (np.sqrt(3) / 2 * a0)
    scale[scale == 0] = 1.0
    dn = d / scale[:, None]
    return np.sqrt(np.mean((dn - ideal) ** 2, axis=1)) / a0


# Absolute disorder cutoff: a perfect (sub-grid-refined) BCC crystal sits at
# ~0.031 (max ~0.042, from residual grid quantization); a melted void rim is
# ~0.13. 0.06 cleanly separates them and, being absolute, does not flag the
# tail of a defect-free crystal the way an adaptive mean+2std threshold does.
DISORDER_CUTOFF = 0.06


def find_dislocation_lines(pts, box, thresh=DISORDER_CUTOFF, link=0.9,
                           a0=A_BCC):
    """Cluster disordered atoms into dislocation lines.

    thresh: absolute disorder cutoff (default DISORDER_CUTOFF = 0.06).
    link:   atoms within link*a0 join the same line cluster.
    Returns dict(n_lines, disordered_frac, line_sizes, line_extents, labels,
                 disorder).
    """
    dis = disorder_metric(pts, box, a0)
    mask = dis > thresh
    core = pts[mask]
    if len(core) == 0:
        return dict(n_lines=0, disordered_frac=0.0, line_sizes=[],
                    line_extents=[], labels=np.full(len(pts), -1),
                    disorder=dis, thresh=float(thresh))

    # connected components on the periodic proximity graph of core atoms
    tiled, src = _tiled(core, box)
    tree = cKDTree(tiled)
    pairs = tree.query_ball_point(core, r=link * a0)
    parent = list(range(len(core)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, neigh in enumerate(pairs):
        for j in neigh:
            j0 = int(src[j])
            if j0 != i:
                ri, rj = find(i), find(j0)
                if ri != rj:
                    parent[ri] = rj

    roots = {}
    for i in range(len(core)):
        roots.setdefault(find(i), []).append(i)
    clusters = [v for v in roots.values() if len(v) >= 3]  # drop point noise

    sizes, extents = [], []
    for c in clusters:
        cp = core[c]
        # periodic extent: max pairwise min-image distance proxy via PCA span
        cc = cp - cp.mean(axis=0)
        cc -= box * np.round(cc / box)
        span = np.sqrt(np.linalg.eigvalsh(np.cov(cc.T) + 1e-9 * np.eye(3)))
        sizes.append(len(c))
        extents.append(float(2 * span.max()))  # ~ line length proxy

    labels = np.full(len(pts), -1)
    core_idx = np.where(mask)[0]
    for cid, c in enumerate(clusters):
        labels[core_idx[c]] = cid
    return dict(n_lines=len(clusters),
                disordered_frac=float(mask.mean()),
                line_sizes=sizes, line_extents=extents,
                total_line_length=float(sum(extents)),
                labels=labels, disorder=dis, thresh=float(thresh))

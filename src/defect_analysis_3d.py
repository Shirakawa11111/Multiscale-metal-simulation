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


def disorder_metric(pts, box, a0=A_BCC, n_neigh=8):
    """Per-atom centrosymmetry parameter (Kelchner et al. 1998).

    CSP = sum over the n_neigh/2 best antipodal neighbour pairs of
    |r_a + r_b|^2, normalised by the mean squared bond length. For ANY
    affine deformation F, an antipodal pair (r, -r) maps to (F r, -F r) whose
    sum is still zero, so CSP is exactly 0 for a perfectly crystalline atom at
    arbitrary (even anisotropic, volume-conserving) strain — unlike a
    fixed-ideal distance metric, which spuriously flags strained-but-perfect
    lattices. CSP is nonzero only at genuine defects (cores, surfaces, GBs).

    n_neigh=8 uses the BCC first shell (8 atoms = 4 antipodal pairs)."""
    n = len(pts)
    if n < n_neigh + 1:
        return np.zeros(n)
    tiled, src = _tiled(pts, box)
    tree = cKDTree(tiled)
    d, idx = tree.query(pts, k=n_neigh + 1)   # self + n_neigh
    csp = np.zeros(n)
    for i in range(n):
        nbr = tiled[idx[i, 1:]] - pts[i]      # n_neigh relative vectors
        mean_sq = np.mean(np.sum(nbr ** 2, axis=1)) + 1e-12
        used = np.zeros(len(nbr), dtype=bool)
        s = 0.0
        order = np.argsort(-np.sum(nbr ** 2, axis=1))  # longest first
        for a in order:
            if used[a]:
                continue
            # best antipode: minimise |r_a + r_b|
            cand = np.where(~used)[0]
            cand = cand[cand != a]
            if len(cand) == 0:
                break
            sums = np.sum((nbr[a] + nbr[cand]) ** 2, axis=1)
            b = cand[np.argmin(sums)]
            used[a] = used[b] = True
            s += sums.min()
        csp[i] = s / mean_sq
    return csp


# Absolute CSP cutoff: a perfect BCC crystal sits at CSP ~0.04 (max ~0.06,
# from sub-grid peak-refinement noise) AT ANY AFFINE STRAIN (CSP is
# affine-invariant — verified synthetically to 1e-5 at 15% strain); genuine
# defects (void rims, GB cores) read ~1-5. 0.3 separates them with wide
# margin and, being affine-invariant, does NOT flag a strained-but-perfect
# lattice (the failure mode of the earlier fixed-distance metric).
DISORDER_CUTOFF = 0.3


def find_dislocation_lines(pts, box, thresh=DISORDER_CUTOFF, link=1.2,
                           a0=A_BCC):
    """Cluster disordered atoms into dislocation lines.

    thresh: absolute CSP cutoff (default DISORDER_CUTOFF = 0.3).
    link:   atoms within link*a0 join the same line cluster (1.2 connects
            defect atoms across the slightly dilated cores/void rims, where
            spacing exceeds the bulk nearest-neighbour distance).
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

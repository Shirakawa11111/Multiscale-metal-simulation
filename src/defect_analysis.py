"""Defect detection for 2D PFC density fields.

Pipeline: density peaks -> Delaunay triangulation -> coordination number ->
5|7 coordination pairs = edge dislocations in the triangular lattice.
Mirrors the 5-7 topology identification used in Topology research/ and the
peak-reconstruction idea from BJK_hBN_Multiscale.
"""

import numpy as np
from scipy.ndimage import maximum_filter, label, center_of_mass
from scipy.spatial import Delaunay, cKDTree


def find_peaks(psi, dx, dy, min_distance_frac=0.5, a0=4 * np.pi / np.sqrt(3)):
    """Return (N,2) array of atom positions [x, y] from density maxima."""
    size = max(3, int(min_distance_frac * a0 / max(dx, dy)))
    mx = maximum_filter(psi, size=size, mode="wrap")
    thresh = psi.mean() + 0.3 * psi.std()
    mask = (psi >= mx) & (psi > thresh)
    lbl, n = label(mask)
    if n == 0:
        return np.zeros((0, 2))
    coms = center_of_mass(psi - psi.min(), lbl, range(1, n + 1))
    pts = np.array([[c[1] * dx, c[0] * dy] for c in coms])
    return pts


def _periodic_delaunay(pts, lx, ly, pad_frac=0.15):
    """Delaunay on a periodically padded copy; returns triangulation plus the
    index map back to original points (or -1 for ghost)."""
    images = [pts]
    idx = [np.arange(len(pts))]
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            if sx == 0 and sy == 0:
                continue
            shifted = pts + np.array([sx * lx, sy * ly])
            keep = ((shifted[:, 0] > -pad_frac * lx) & (shifted[:, 0] < (1 + pad_frac) * lx)
                    & (shifted[:, 1] > -pad_frac * ly) & (shifted[:, 1] < (1 + pad_frac) * ly))
            images.append(shifted[keep])
            idx.append(np.where(keep)[0])
    allpts = np.vstack(images)
    src = np.concatenate(idx)
    tri = Delaunay(allpts)
    return tri, allpts, src


def coordination(pts, lx, ly):
    """Coordination number per point via periodic Delaunay neighbors."""
    tri, allpts, src = _periodic_delaunay(pts, lx, ly)
    n = len(pts)
    neigh = [set() for _ in range(n)]
    indptr, indices = tri.vertex_neighbor_vertices
    for i in range(n):  # first n entries of allpts are the originals
        for j in indices[indptr[i]:indptr[i + 1]]:
            neigh[i].add(int(src[j]))
        neigh[i].discard(i)
    return np.array([len(s) for s in neigh]), neigh


def find_dislocations(pts, lx, ly, a0=4 * np.pi / np.sqrt(3)):
    """Pair each 5-coordinated atom with the nearest 7-coordinated atom.
    Returns dict with positions of 5s, 7s, paired dipole midpoints (= dislocation
    core positions), and areal dislocation density."""
    if len(pts) < 10:
        return dict(n5=0, n7=0, cores=np.zeros((0, 2)), rho=0.0,
                    fives=np.zeros((0, 2)), sevens=np.zeros((0, 2)))
    coord, _ = coordination(pts, lx, ly)
    fives = pts[coord == 5]
    sevens = pts[coord == 7]
    cores = []
    if len(fives) and len(sevens):
        tree = cKDTree(sevens, boxsize=None)
        used = set()
        for f in fives:
            d, j = tree.query(f, k=min(len(sevens), 4))
            d = np.atleast_1d(d)
            j = np.atleast_1d(j)
            for dd, jj in zip(d, j):
                if jj not in used and dd < 2.0 * a0:
                    used.add(int(jj))
                    cores.append(0.5 * (f + sevens[jj]))
                    break
    cores = np.array(cores) if cores else np.zeros((0, 2))
    rho = len(cores) / (lx * ly)
    return dict(n5=int((coord == 5).sum()), n7=int((coord == 7).sum()),
                cores=cores, rho=rho, fives=fives, sevens=sevens)


def lattice_spacing(pts, lx, ly):
    """Median nearest-neighbor distance (periodic, approximate via KDTree
    on tiled points). Biased low in polycrystals; prefer density_spacing."""
    if len(pts) < 4:
        return np.nan
    tiled = []
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            tiled.append(pts + np.array([sx * lx, sy * ly]))
    tree = cKDTree(np.vstack(tiled))
    d, _ = tree.query(pts, k=2)
    return float(np.median(d[:, 1]))


def density_spacing(pts, lx, ly):
    """Lattice constant inferred from areal peak density assuming a perfect
    triangular lattice: area/atom = (√3/2) a²  =>  a = sqrt(2A/(√3 N)).
    Robust against grain boundaries (unbiased by defective neighborhoods)."""
    if len(pts) < 4:
        return np.nan
    return float(np.sqrt(2.0 * lx * ly / (np.sqrt(3.0) * len(pts))))

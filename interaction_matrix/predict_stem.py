"""Phase 2: predict the real STEM-reconstructed network's forest-hardening
coefficient from the first-principles interaction matrix + the network's own
slip-system / junction inventory. Connects the measured a_ij (Phase 1) to the
experimental anchor.

Generalized Taylor law for a population of slip systems:
    alpha_network = sqrt( sum_i sum_j P_i P_j a_{type(i,j)} ),
where P_i = (line length on system i)/(total length) and a_{type} is the
measured interaction coefficient for the junction type of the (i,j) pair.

Because the experimental Burgers vectors were assigned GEOMETRICALLY (a known
weakness), we also run a Burgers-resampling sensitivity: for each line keep its
(geometric) glide plane but resample its Burgers among the <110> directions in
that plane, recompute alpha_network over many samples -> a prediction
DISTRIBUTION that quantifies the uncertainty from the unknown true Burgers.

  python3 predict_stem.py
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NET = os.path.join(HERE, "..", "experiment_bridge", "results_exadis", "stem_network.json")
MATRIX = os.path.join(HERE, "matrix_result_v2.json")

# FCC <110> Burgers in each {111} plane (integer)
PLANES = {(1, 1, 1): [(0, 1, -1), (1, 0, -1), (1, -1, 0)],
          (-1, 1, 1): [(0, 1, -1), (1, 0, 1), (1, 1, 0)],
          (1, -1, 1): [(0, 1, 1), (1, 0, -1), (1, 1, 0)],
          (1, 1, -1): [(0, 1, 1), (1, 0, 1), (1, -1, 0)]}


def plane_key(n):
    n = np.round(np.array(n) / (np.abs(n[np.argmax(np.abs(n))]) or 1)).astype(int)
    for k in PLANES:
        if np.all(np.cross(k, n) == 0):
            return k
    return None


def jtype(b1, n1, b2, n2):
    b1, n1, b2, n2 = map(lambda v: np.round(np.array(v)).astype(int), (b1, n1, b2, n2))
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


def load_lines():
    net = json.load(open(NET))
    nodes = np.array([n[:3] for n in net["nodes"]], float)
    lines = {}    # (burg_tuple, plane_tuple) -> total length
    perline = []  # (burg int, plane int, length)
    # group segs by their (burg,plane); each contributes its length
    segdata = {}
    for s in net["segs"]:
        a, b = int(s[0]), int(s[1])
        L = float(np.linalg.norm(nodes[a] - nodes[b]))
        burg = tuple(np.round(np.array(s[2:5]) * np.sqrt(2)).astype(int))   # -> <110> int
        pk = plane_key(s[5:8])
        segdata.setdefault((burg, pk), 0.0)
        segdata[(burg, pk)] += L
    return segdata


def alpha_from_inventory(systems, amat):
    """systems: list of (P_i, burg_i, plane_i). Returns alpha_network."""
    tot = 0.0
    for Pi, bi, ni in systems:
        for Pj, bj, nj in systems:
            t = jtype(bi, ni, bj, nj)
            tot += Pi * Pj * amat.get(t, 0.0)
    return float(np.sqrt(max(tot, 0.0)))


def main():
    amat = json.load(open(MATRIX))["a_matrix"]
    seg = load_lines()
    total = sum(seg.values())
    systems = [(L / total, list(b), list(p)) for (b, p), L in seg.items()]
    print("STEM network slip-system inventory (geometric assignment):")
    for P, b, p in sorted(systems, key=lambda x: -x[0]):
        print(f"  P={P:.3f}  b={b}  n={p}")
    # junction inventory among the populated systems
    inv = {}
    for Pi, bi, ni in systems:
        for Pj, bj, nj in systems:
            t = jtype(bi, ni, bj, nj)
            inv[t] = inv.get(t, 0.0) + Pi * Pj
    print("\njunction-type population weights P_i*P_j:")
    for t, w in sorted(inv.items(), key=lambda x: -x[1]):
        print(f"  {t:>10}: {w:.3f}   (a_ij={amat.get(t,0):.3f})")

    a_geom = alpha_from_inventory(systems, amat)
    print(f"\nPredicted alpha_network (geometric Burgers) = {a_geom:.3f}")

    # Burgers-resampling sensitivity: keep planes, resample Burgers in-plane
    rng = np.random.default_rng(0)
    planes_lengths = {}
    for (b, p), L in seg.items():
        planes_lengths.setdefault(p, 0.0)
        planes_lengths[p] += L
    samples = []
    for it in range(2000):
        syslist = []
        for p, L in planes_lengths.items():
            burg = PLANES[p][rng.integers(0, 3)]
            syslist.append((L / total, list(burg), list(p)))
        samples.append(alpha_from_inventory(syslist, amat))
    samples = np.array(samples)
    print(f"\nBurgers-resampling sensitivity (planes fixed, Burgers random in-plane, n=2000):")
    print(f"  alpha_network = {samples.mean():.3f} +/- {samples.std():.3f}  "
          f"(5-95%: {np.percentile(samples,5):.3f}-{np.percentile(samples,95):.3f})")
    print(f"  uniform-population macro alpha (reference) = "
          f"{np.sqrt(sum(amat.values())*0+ json.load(open(MATRIX))['alpha_macro']**2):.3f}")

    out = dict(alpha_geometric=a_geom, junction_inventory=inv,
               alpha_resampled_mean=float(samples.mean()),
               alpha_resampled_std=float(samples.std()),
               alpha_resampled_p5=float(np.percentile(samples, 5)),
               alpha_resampled_p95=float(np.percentile(samples, 95)),
               n_systems=len(systems), macro_uniform=json.load(open(MATRIX))["alpha_macro"],
               note="alpha_network from measured a_ij + STEM junction inventory; "
                    "resampling quantifies the geometric-Burgers uncertainty.")
    json.dump(out, open(os.path.join(HERE, "stem_prediction.json"), "w"), indent=1)
    print("\nsaved", os.path.join(HERE, "stem_prediction.json"))


if __name__ == "__main__":
    main()

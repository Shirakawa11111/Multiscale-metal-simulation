"""FCC slip systems and dislocation-junction classification (the crystallographic
bedrock for a first-principles interaction-matrix measurement).

12 FCC slip systems = 4 {111} planes x 3 <110> Burgers in each plane. For an
ordered pair of systems (mobile m, forest f) the junction type is one of six
(Madec, Devincre & Kubin, Science 2003): self, coplanar, collinear, Hirth,
glissile, Lomer(-Cottrell). We classify geometrically and VALIDATE against the
known multiplicities of the 12x12=144 ordered pairs:
    self 12, coplanar 24, collinear 12, Hirth 24, glissile 48, Lomer 24.
Per system (11 partners): 2 coplanar, 1 collinear, 2 Hirth, 4 glissile, 2 Lomer.
"""
import numpy as np
import itertools


def fcc_systems():
    """Return list of (b, n) integer vectors: b = <110> Burgers (in-plane),
    n = {111} plane normal, with b . n = 0."""
    planes = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]
    burgers = [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]
    systems = []
    for n in planes:
        n = np.array(n)
        for b in burgers:
            b = np.array(b)
            if np.dot(b, n) == 0:               # Burgers lies in the glide plane
                systems.append((b, n))
    return systems


def _parallel(u, v):
    return np.all(np.cross(u, v) == 0)


def junction_type(s1, s2):
    b1, n1 = s1
    b2, n2 = s2
    if _parallel(b1, b2) and _parallel(n1, n2):
        return "self"
    if _parallel(n1, n2):
        return "coplanar"          # same glide plane, different Burgers
    if _parallel(b1, b2):
        return "collinear"         # same Burgers, different plane (annihilation)
    # different plane AND different Burgers -> junction by b3 = b1 +/- b2
    d = int(np.dot(b1, b2))        # for <110>: 0 (perp) or +/-1 (60 deg)
    if d == 0:
        return "Hirth"             # perpendicular Burgers -> Hirth lock
    # d = +/-1: glissile if the junction Burgers lies in a glide plane of either
    # parent (it can glide); else Lomer-Cottrell sessile (junction on {100}).
    for sign in (1, -1):
        b3 = b1 + sign * b2
        if np.all(b3 == 0):
            continue
        # is b3 a <110> lattice glide direction lying in n1 or n2?
        is110 = sorted(np.abs(b3)) == [0, 1, 1]
        if is110 and (np.dot(b3, n1) == 0 or np.dot(b3, n2) == 0):
            return "glissile"
    return "Lomer"                 # junction Burgers is <100>-type / sessile


def classify_all():
    S = fcc_systems()
    counts = {}
    pairs = {}
    for i, j in itertools.product(range(len(S)), range(len(S))):
        t = junction_type(S[i], S[j])
        counts[t] = counts.get(t, 0) + 1
        pairs.setdefault(t, []).append((i, j))
    return S, counts, pairs


if __name__ == "__main__":
    S, counts, pairs = classify_all()
    print(f"FCC slip systems: {len(S)}  (expect 12)")
    expected = {"self": 12, "coplanar": 24, "collinear": 12,
                "Hirth": 24, "glissile": 48, "Lomer": 24}
    print(f"\n{'type':>10} {'measured':>9} {'expected':>9}  {'ok':>3}")
    ok_all = True
    for t in ["self", "coplanar", "collinear", "Hirth", "glissile", "Lomer"]:
        m = counts.get(t, 0); e = expected[t]; ok = (m == e)
        ok_all &= ok
        print(f"{t:>10} {m:>9} {e:>9}  {'OK' if ok else 'XX':>3}")
    print(f"\ntotal pairs: {sum(counts.values())} (expect 144)")
    print(f"CLASSIFICATION VALIDATION: {'PASS' if ok_all else 'FAIL'}")
    # one representative pair per type (for the controlled measurement)
    print("\nrepresentative (mobile, forest) pairs per junction type:")
    for t in ["self", "coplanar", "collinear", "Hirth", "glissile", "Lomer"]:
        if pairs.get(t):
            i, j = next((a, b) for a, b in pairs[t] if a != b) if t != "self" else pairs[t][0]
            print(f"  {t:>10}: m=sys{i} b={S[i][0]} n={S[i][1]}  |  f=sys{j} b={S[j][0]} n={S[j][1]}")

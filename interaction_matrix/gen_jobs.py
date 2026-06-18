"""Emit the job list for the large-scale FCC interaction-matrix campaign.
Enumerates (mobile, forest) system pairs per junction type (ExaDiS ordering),
selects K diverse pairs per type, and crosses with forest densities and seeds.
Big-first (highest NFOREST first) so long jobs overlap. Prints one env-prefixed
command per line for the semaphore driver.

  python3 gen_jobs.py <ROOT> <TH>
"""
import sys
import numpy as np
from collections import defaultdict

BI = np.array([[0, 1, -1], [1, 0, -1], [1, -1, 0], [0, 1, -1], [1, 0, 1], [1, 1, 0],
               [0, 1, 1], [1, 0, -1], [1, 1, 0], [0, 1, 1], [1, 0, 1], [1, -1, 0]], float)
NI = np.array([[1, 1, 1]] * 3 + [[-1, 1, 1]] * 3 + [[1, -1, 1]] * 3 + [[1, 1, -1]] * 3, float)


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


ROOT = sys.argv[1] if len(sys.argv) > 1 else "matrix"
TH = sys.argv[2] if len(sys.argv) > 2 else "3"
K_PAIRS = 4
SEEDS = [1234, 5678, 2222, 3333]
NFORESTS = [24, 12, 6]            # big-first (3 densities for the slope alpha_mf)
# fast jobs (the matrix is RELATIVE across junction types; the shared probe
# self-baseline cancels in the slope, so a higher rate is acceptable here)
COMMON = "LBOX=7000 NPROBE=2 ERATE=1e4 MAX_STRAIN=0.0003 FLOW_LO=0.0002 MAXSEG=400"

bytype = defaultdict(list)
for m in range(12):
    for f in range(12):
        bytype[jtype(m, f)].append((m, f))

# select K diverse pairs per type (evenly spaced)
sel = {}
for t, prs in bytype.items():
    idx = np.linspace(0, len(prs) - 1, min(K_PAIRS, len(prs))).round().astype(int)
    sel[t] = [prs[i] for i in sorted(set(idx))]

lines = []
for nf in NFORESTS:                 # big-first: outer loop over density desc
    for t, prs in sel.items():
        for (m, f) in prs:
            for s in SEEDS:
                tag = f"{t}_m{m}_f{f}_n{nf}_s{s}"
                lines.append(
                    f"OMP_NUM_THREADS={TH} MSYS={m} FSYS={f} NFOREST={nf} SEED={s} "
                    f"{COMMON} OUT={ROOT}/{tag} python3 -u build_pair.py "
                    f"> {ROOT}/{tag}.log 2>&1")
for ln in lines:
    print(ln)
import sys as _s
print(f"# TOTAL {len(lines)} jobs", file=_s.stderr)

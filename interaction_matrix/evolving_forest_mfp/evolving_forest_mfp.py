"""
evolving_forest_mfp.py  --  Pair-resolved EVOLVING-forest mean-free-path / storage assay.

GOAL: measure the COLLECTIVE dislocation MFP / storage rate of a mobile slip system m
flowing through an EVOLVABLE forest slip system f, to test whether collinear OPPOSITE-sense
annihilation shortens the MFP far more than other junction types (collective effect, NOT a
single-junction local strength).

Modules (all verified against binary_strength.py / binary_collinear.py / build_pair.py that RUN):
  DDD_FFT_MODEL (full pairwise elastic) + full PBC + Collision(Retroactive)
  + Topology(TopologyParallel) + Remesh(LengthBased) + FCC_0 mobility, Trapezoid integrator.

KEY DESIGN POINTS (grafted from the 3 designs, with every blocker/major fix applied):
  * Forest built with insert_infinite_line and NOT re-pinned -> EVOLVABLE (departure from
    build_pair.py which pins the forest). Optional single light anchor node per line, gated
    by a forest-alone stability control (run separately with NO mobile sources).
  * Mobile = a FIXED number of pinned-end Frank-Read-like sources (verified add_segment idiom
    from binary_strength.py / build_src.py). Fixed probe count is mandatory (PROGRESS lesson:
    fractional/varying carriers give flat noisy signal).
  * Sense control by Burgers SIGN on the 2nd {111}: opposite (b_f=-b_m) annihilates (validated
    binary_collinear.py -> length_fraction 0.71); same (b_f=+b_m) does not. The canonical
    reaction MUST be pre-validated with binary_observe.py-style S1/S2 scan before trusting the
    population run (env PAIR_SMOKE=1 runs that 2-segment check).
  * gamma_m = MOBILE-ONLY swept area (Burgers-filtered to FAM_m), accumulated EVERY step against
    positions cached at the START of the step (override step hook), PBC-folded. Global dEp is a
    cross-check only (it sums ALL systems and cannot isolate m, fatal for collinear).
  * Ledger classifies by PLANE (primary; collinear shares Burgers so Burgers cannot separate
    mobile from forest) + node degree + glide velocity. ALL lengths PBC-folded (closest_image).
    Screw segments (null plane) handled explicitly. Total length cross-checked vs
    dislocation_density() (not the vacuous self-closure check).
  * REMOVED/annihilated length computed two independent ways and BOTH logged:
      (1) per-step negative jumps in rho_total that coincide with a collision/topology event
          (segment-count drop), integrated -> rho_removed_events;
      (2) length balance against an explicit source term:
          rho_removed_bal = (L_built + L_created_by_sources) - L_present - L_into_junction,
          where L_created_by_sources is tracked as the cumulative POSITIVE rho_total jumps that
          coincide with NO collision (pure bow-out / multiplication / remesh growth).
    Storage slope and annihilation slope are reported SEPARATELY (never summed into one rho_stop),
    because rho_stored ~ gamma/(bL) assumes RETAINED length, while removed is destroyed length.
  * Guards: spanB>0.45*LBOX hard-flag (verified binary_strength.py guard); carrier-survival
    window 0.3-0.95 for the fit; fixed tau (no auto-ladder); one OS process == one network ==
    one (relax + load) -> avoids the documented pyexadis double-free.

ENV (parameterized like our other scripts):
  JTYPE in {coll_opp, coll_same, glissile, hirth, lomer, coplanar, self}
  SENSE  -1 (opposite, annihilate) | +1 (same)   [only used to override coll sign]
  B_M,N_M,B_F,N_F  comma vectors (optional explicit override of the pair)
  RHO_F (target m^-2, sets NFOREST if NFOREST unset) | NFOREST (direct line count)
  TAU_MPA   resolved shear on m only
  LBOX (b)  NSTEPS  SEED  OUT
  Controls: FORCE=DDD_FFT_MODEL|LINE_TENSION  TOPO=1|0  COLL=1|0  KM (mobile sources, default 6)
            FOREST_ALONE=1 (no mobile sources -> forest self-evolution control)
            PAIR_SMOKE=1   (2-segment binary_observe-style sense check, no forest)
            ANCHOR=0|1 (light 1-node-per-forest-line face anchor)

  Example:
    JTYPE=coll_opp RHO_F=3e12 TAU_MPA=45 LBOX=20000 NSTEPS=6000 SEED=1 \
    OUT=mfp/coll_opp_s1 PYTHONPATH=~/BO/exadis_src/python python3 evolving_forest_mfp.py
"""
import os, sys, json
import numpy as np

EX = os.path.expanduser("~/BO/exadis_src/python")
if EX not in sys.path:
    sys.path.append(EX)
import pyexadis
from pyexadis_base import (ExaDisNet, DisNetManager, SimulateNetwork, CalForce,
                           MobilityLaw, TimeIntegration, Collision, Topology, Remesh,
                           NodeConstraints)
from pyexadis_utils import insert_infinite_line, dislocation_density

# ---- verified ExaDiS FCC 12-system table (build_pair.py) ----
BI = np.array([[0, 1, -1], [1, 0, -1], [1, -1, 0], [0, 1, -1], [1, 0, 1], [1, 1, 0],
               [0, 1, 1], [1, 0, -1], [1, 1, 0], [0, 1, 1], [1, 0, 1], [1, -1, 0]], float)
NI = np.array([[1, 1, 1]] * 3 + [[-1, 1, 1]] * 3 + [[1, -1, 1]] * 3 + [[1, 1, -1]] * 3, float)
B_CU, MU = 2.55e-10, 54.6e9

# ----------------------------------------------------------------------------------------
# ENV
# ----------------------------------------------------------------------------------------
def envf(k, d): return float(os.environ.get(k, d))
def envi(k, d): return int(os.environ.get(k, d))
def envvec(k):
    v = os.environ.get(k, "")
    return np.array([float(x) for x in v.split(",")]) if v else None

JTYPE   = os.environ.get("JTYPE", "coll_opp")
SENSE   = envi("SENSE", -1)                    # -1 opposite (annihilate), +1 same
LBOX    = envf("LBOX", "20000")
TAU     = envf("TAU_MPA", "45") * 1e6
RHO_F   = envf("RHO_F", "3e12")
NFOREST = envi("NFOREST", "0")                 # 0 -> derive from RHO_F
KM      = envi("KM", "6")                      # FIXED mobile source/line count
MOBILE_MODE = os.environ.get("MOBILE_MODE", "fr")  # 'fr' = FR pinned sources (multiply);
                                                   # 'lines' = infinite gliding lines (NO multiplication)
NSTEPS  = envi("NSTEPS", "6000")
NREL    = envi("NREL", "300")                  # zero-stress settle steps
REC     = envi("REC", "20")                    # ledger record cadence
SEED    = envi("SEED", "1")
MAXSEG  = envf("MAXSEG", "80")
NGRID   = envi("NGRID", "64")                  # convergence knob (verify 32/64/128)
RANN    = envf("RANN", "10")
FORCE   = os.environ.get("FORCE", "DDD_FFT_MODEL")
TOPO    = os.environ.get("TOPO", "1") == "1"
COLL    = os.environ.get("COLL", "1") == "1"
ANCHOR  = os.environ.get("ANCHOR", "0") == "1"
FOREST_ALONE = os.environ.get("FOREST_ALONE", "0") == "1"
PAIR_SMOKE   = os.environ.get("PAIR_SMOKE", "0") == "1"
OUT     = os.environ.get("OUT", "mfp_out")
rng = np.random.default_rng(SEED)

def hat(v):
    v = np.asarray(v, float); n = np.linalg.norm(v)
    return v / n if n > 0 else v
def cos_par(u, v):
    u, v = hat(u), hat(v)
    return abs(float(np.dot(u, v)))          # |cos|, sense-agnostic

# ----------------------------------------------------------------------------------------
# PAIR DEFINITION  (m fixed = system 0; forest f per JTYPE; sense by Burgers sign)
# Validated reactions: binary_collinear.py uses b_f = -b_m on a 2nd {111} for OPPOSITE.
# ----------------------------------------------------------------------------------------
def pick_pair():
    bm, nm = BI[0].copy(), NI[0].copy()       # m = [0,1,-1] on (111)
    eb, en = envvec("B_M"), envvec("N_M")
    if eb is not None: bm = eb
    if en is not None: nm = en
    bf = envvec("B_F"); nf = envvec("N_F")
    if bf is not None and nf is not None:
        return bm, nm, bf, nf
    j = JTYPE.lower()
    if j in ("coll_opp", "coll_same", "collinear"):
        nf = NI[3].copy()                     # 2nd {111}: (-1,1,1); same |b| family
        s = -1.0 if (j == "coll_opp" or (j == "collinear" and SENSE < 0)) else +1.0
        bf = s * bm.copy()
    elif j == "glissile":
        bf, nf = BI[4].copy(), NI[4].copy()   # shares a plane family; b3 = <110> glissile
    elif j == "hirth":
        bf, nf = BI[6].copy(), NI[6].copy()   # dot(b_m,b_f)=0 -> Hirth lock
    elif j == "lomer":
        bf, nf = BI[5].copy(), NI[5].copy()
    elif j == "coplanar":
        bf, nf = BI[1].copy(), NI[1].copy()   # same plane (111), different <110>
    elif j == "self":
        bf, nf = bm.copy(), nm.copy()
    else:
        raise ValueError("unknown JTYPE %s" % JTYPE)
    return bm, nm, bf, nf

bm, nm, bf, nf = pick_pair()
b3a, b3b = bm + bf, bm - bf
# A GENUINE junction Burgers must be non-null AND not parallel to either parent.
# (collinear bf=-bm: b3a=0, b3b=2bm||bm -> BOTH excluded -> FAM_J empty, as physics requires.)
FAM_J = [v for v in (b3a, b3b)
         if np.linalg.norm(v) > 1e-6 and cos_par(v, bm) < 0.94 and cos_par(v, bf) < 0.94]
# m Schmid tensor (for stress + gamma resolution); ExaDiS Voigt order [xx,yy,zz,yz,xz,xy]
A_m = 0.5 * (np.outer(hat(bm), hat(nm)) + np.outer(hat(nm), hat(bm)))
A2  = np.array([A_m[0, 0], A_m[1, 1], A_m[2, 2], 2 * A_m[1, 2], 2 * A_m[0, 2], 2 * A_m[0, 1]])

# ----------------------------------------------------------------------------------------
# segs burgers/planes accessor  (robust: prefer named keys, fall back to raw matrix)
# segs were built as rows [nA, nB, bx,by,bz, nx,ny,nz].
# ----------------------------------------------------------------------------------------
def seg_arrays(data):
    S = data["segs"]
    nid = np.asarray(S["nodeids"]).astype(int)
    if "burgers" in S and "planes" in S:
        return nid, np.asarray(S["burgers"], float), np.asarray(S["planes"], float)
    raw = None
    for k in ("data", "array", "segs"):
        if k in S:
            raw = np.asarray(S[k], float); break
    if raw is not None and raw.shape[1] >= 8:
        return nid, raw[:, 2:5], raw[:, 5:8]
    # last resort: no per-seg burgers/planes -> plane-less; classifier degrades to degree only
    return nid, None, None

def closest_vec(cell, ra, rb):
    """PBC minimum-image vector rb-ra using cell.closest_image if available."""
    try:
        rb_img = np.array(cell.closest_image(Rref=ra, R=rb))
        return rb_img - ra
    except Exception:
        d = rb - ra
        H = np.array(cell.h)
        L = np.array([H[0, 0], H[1, 1], H[2, 2]])
        return d - L * np.round(d / L)        # orthorhombic fallback

# ----------------------------------------------------------------------------------------
# COARSE density ledger.  Classify each segment by PLANE (primary) + Burgers-family +
# node degree.  PLANE is primary because collinear shares the Burgers vector.
# ----------------------------------------------------------------------------------------
class Ledger:
    def __init__(self, cell):
        self.cell = cell
        self.vol = abs(np.linalg.det(np.array(cell.h)))   # b^3
        self.rows = []
        self.L_built = None
        self.mob_peak = 1e-30
        self.rho_tot_prev = None
        self.nseg_prev = None
        self.cum_removed_events = 0.0     # rho-units, from collision-coincident drops
        self.cum_created = 0.0            # rho-units, from non-collision growth (source/bow/remesh)

    def _to_rho(self, L_b):
        return L_b / self.vol / B_CU ** 2  # 1/m^2  (L in b, vol in b^3, /b^2)

    def classify_one(self, b_i, p_i, deg_a, deg_c, vmag):
        """Return class key. PLANE-first for collinear; Burgers-family otherwise."""
        anchored = (deg_a >= 3 or deg_c >= 3)
        on_jun_b = any(cos_par(b_i, r) > 0.94 for r in FAM_J) if FAM_J else False
        if on_jun_b:
            return "junction"
        # plane test (guard null/screw plane: norm ~ 0)
        pl_ok = (p_i is not None) and (np.linalg.norm(p_i) > 1e-6)
        on_nm = pl_ok and cos_par(p_i, nm) > 0.97
        on_nf = pl_ok and cos_par(p_i, nf) > 0.97
        is_m_b = cos_par(b_i, bm) > 0.94
        is_f_b = cos_par(b_i, bf) > 0.94
        # COLLINEAR: bm and bf are parallel -> Burgers cannot separate; use PLANE.
        collinear = cos_par(bm, bf) > 0.94
        if collinear:
            if anchored:
                return "residual"
            if on_nm:
                return "mobile"
            if on_nf:
                return "forest"
            # screw / off-plane / ambiguous: use glide-velocity as tiebreak
            return "mobile" if vmag > VMOB_THR else "residual"
        # NON-collinear: Burgers family is decisive, plane corroborates
        if is_m_b:
            if anchored: return "residual"
            return "mobile" if (on_nm or not pl_ok) else "residual"
        if is_f_b:
            if anchored: return "residual"
            return "forest" if (on_nf or not pl_ok) else "residual"
        return "offsys"

    def snapshot(self, net, gamma_m, istep, vmag_node, span_frac=0.0):
        d = net.get_disnet(ExaDisNet).export_data()
        pos = d["nodes"]["positions"]
        nid, burg, plane = seg_arrays(d)
        deg = {}
        for a, c in nid:
            deg[int(a)] = deg.get(int(a), 0) + 1
            deg[int(c)] = deg.get(int(c), 0) + 1
        pool = {"mobile": 0., "forest": 0., "junction": 0., "residual": 0., "offsys": 0.}
        ambiguous = 0.0
        for i, (a, c) in enumerate(nid):
            a, c = int(a), int(c)
            L = float(np.linalg.norm(closest_vec(self.cell, pos[a], pos[c])))
            b_i = burg[i] if burg is not None else bm
            p_i = plane[i] if plane is not None else None
            vm = 0.5 * (vmag_node.get(a, 0.0) + vmag_node.get(c, 0.0))
            cls = self.classify_one(b_i, p_i, deg.get(a, 0), deg.get(c, 0), vm)
            pool[cls] += L
            if cls == "offsys":
                ambiguous += L
        Ltot = sum(pool.values())
        rho = {k: self._to_rho(v) for k, v in pool.items()}
        rho_total = self._to_rho(Ltot)
        if self.L_built is None:
            self.L_built = Ltot
        self.mob_peak = max(self.mob_peak, rho["mobile"])
        nseg = len(nid)
        # ---- direct removed/created accounting via per-record rho_total deltas ----
        rho_removed_events = self.cum_removed_events
        if self.rho_tot_prev is not None:
            d_rho = rho_total - self.rho_tot_prev
            seg_dropped = (self.nseg_prev is not None and nseg < self.nseg_prev)
            if d_rho < 0 and seg_dropped:
                self.cum_removed_events += (-d_rho)     # length lost at a topo/collision event
            elif d_rho > 0:
                self.cum_created += d_rho                # source bow-out / multiplication / remesh
            rho_removed_events = self.cum_removed_events
        self.rho_tot_prev = rho_total
        self.nseg_prev = nseg
        # ---- length-balance removed (independent cross-check, with explicit source term) ----
        rho_built = self._to_rho(self.L_built)
        rho_removed_bal = max(0.0, rho_built + self.cum_created - rho_total - rho["junction"])
        # ---- stored vs removed kept SEPARATE ----
        rho_stored = rho["junction"] + rho["residual"]     # RETAINED immobile product
        row = dict(
            istep=int(istep), gamma_m=float(gamma_m),
            rho_mobile=rho["mobile"], rho_forest=rho["forest"],
            rho_junction=rho["junction"], rho_residual=rho["residual"],
            rho_offsys=rho["offsys"], rho_total=rho_total,
            rho_stored=rho_stored,
            rho_removed_events=float(rho_removed_events),
            rho_removed_balance=float(rho_removed_bal),
            mobile_survival=rho["mobile"] / self.mob_peak,
            ambiguous_frac=float(ambiguous / Ltot) if Ltot > 0 else 0.0,
            n_seg=int(nseg),
            mobile_span_frac=float(span_frac),
        )
        self.rows.append(row)
        return row

# velocity threshold for "moving" (b/s); FCC_0 vmax=4000 m/s scaled to b -> set generously low
VMOB_THR = envf("VMOB_THR", "1.0")   # m/s; below this a mobile-Burgers arm counts as stored

# ----------------------------------------------------------------------------------------
# SimulateNetwork subclass:
#   - cache node positions + velocities at the START of every step (step hook fires before
#     topology/remesh in this build's loop, mirroring the swept-area xold convention);
#   - accumulate MOBILE-ONLY swept-area gamma_m every step;
#   - record the ledger every REC steps.
# ----------------------------------------------------------------------------------------
class MFPSim(SimulateNetwork):
    def attach(self, ledger):
        self.led = ledger
        self.cell = ledger.cell
        self.gamma_m = 0.0
        self.gamma_m_global = 0.0
        self._old = None          # dict tag-> position at step start
        self.flagged_wrap = False

    def _node_state(self, N):
        d = N.get_disnet(ExaDisNet).export_data()
        pos = d["nodes"]["positions"]
        tags = d["nodes"].get("tags", None)
        if tags is not None:
            keys = [tuple(int(x) for x in t) for t in np.asarray(tags)]
        else:
            keys = list(range(len(pos)))   # fallback: row index (only safe intra-step)
        # velocities for the stored/mobile velocity tiebreak (best-effort)
        vmag = {}
        try:
            G = N.get_disnet(ExaDisNet)
            vel = G.get_velocities() if hasattr(G, "get_velocities") else None
            if vel is not None:
                vel = np.asarray(vel)
                for i in range(len(pos)):
                    vmag[i] = float(np.linalg.norm(vel[i]))
        except Exception:
            pass
        return d, pos, keys, vmag

    def step_begin(self, N, state):
        # cache OLD positions keyed by stable node tag (verified: nodes carry 'tags')
        d, pos, keys, _ = self._node_state(N)
        self._old = {keys[i]: pos[i].copy() for i in range(len(pos))}
        # also keep current burgers/nodeids snapshot for mobile-seg filtering
        nid, burg, plane = seg_arrays(d)
        self._old_segs = (nid, burg, plane)
        self._old_keys = keys

    def _mobile_swept_gamma(self, N):
        """gamma increment from MOBILE (FAM_m, plane nm) segments only, swept-area kernel,
        evaluated against positions cached at step_begin. PBC-folded. Burgers-filtered so
        forest/junction motion (and collinear annihilation of the FOREST) does not leak in."""
        if self._old is None:
            return 0.0
        d, pos, keys, _ = self._node_state(N)
        nid, burg, plane = seg_arrays(d)
        keypos_new = {keys[i]: pos[i] for i in range(len(pos))}
        dE = np.zeros((3, 3))
        for i, (a, c) in enumerate(nid):
            a, c = int(a), int(c)
            b_i = burg[i] if burg is not None else bm
            if cos_par(b_i, bm) <= 0.94:
                continue
            p_i = plane[i] if plane is not None else None
            if (p_i is not None) and np.linalg.norm(p_i) > 1e-6 and cos_par(p_i, nm) <= 0.90:
                continue                      # only m-plane arms (forest collinear excluded)
            ka, kc = keys[a], keys[c]
            r1, r2 = pos[a], pos[c]
            r3 = self._old.get(ka, r1)        # old A (step start); new node -> zero swept
            r4 = self._old.get(kc, r2)        # old B
            # swept-area kernel (pyexadis plastic_strain): n = 0.5*cross(r2-r3, r1-r4), PBC-fold
            v1 = closest_vec(self.cell, r3, r2)   # r2 - r3
            v2 = closest_vec(self.cell, r4, r1)   # r1 - r4
            nA = 0.5 * np.cross(v1, v2)
            P = np.outer(nA, bm)
            dE += 0.5 * (P + P.T)
        dE /= self.led.vol                     # /V (b^3); gamma is dimensionless
        return float(np.dot(np.array([dE[0,0],dE[1,1],dE[2,2],2*dE[1,2],2*dE[0,2],2*dE[0,1]]), A2))

    def step_end(self, N, state):
        # PRIMARY gamma_m: mobile-only swept area (every step)
        self.gamma_m += self._mobile_swept_gamma(N)
        # CROSS-CHECK gamma_m: global dEp resolved on m (sums ALL systems; upper bound only)
        dEp = np.asarray(state.get("dEp", np.zeros(6)), float)
        if dEp.shape == (6,):
            self.gamma_m_global += float(np.dot(dEp, A2))
        # span-wrap guard: MOBILE arms only (plane||nm, Burgers||bm, unpinned). The forest
        # infinite lines legitimately span the whole box, so they must NOT trip this flag.
        if (state.get("istep", 0) % REC) == 0:
            d, pos, keys, vmag = self._node_state(N)
            con = d["nodes"]["constraints"][:, 0]
            nid, burg, plane = seg_arrays(d)
            mob_nodes = set()
            for i, (a, c) in enumerate(nid):
                a, c = int(a), int(c)
                b_i = burg[i] if burg is not None else bm
                p_i = plane[i] if plane is not None else None
                pl_ok = (p_i is not None) and np.linalg.norm(p_i) > 1e-6
                if cos_par(b_i, bm) > 0.94 and (cos_par(p_i, nm) > 0.90 if pl_ok else True):
                    if con[a] != int(NodeConstraints.PINNED_NODE): mob_nodes.add(a)
                    if con[c] != int(NodeConstraints.PINNED_NODE): mob_nodes.add(c)
            span_frac = 0.0
            if mob_nodes:
                mp = pos[list(mob_nodes)]
                # per-axis span; KM sources are spread across ~0.6*LBOX at build, so only a
                # near-full-box span (genuine PBC wrap of a single glide front) should flag.
                span_frac = float(np.max(mp.max(0) - mp.min(0)) / LBOX)
                if span_frac > 0.85:
                    self.flagged_wrap = True
            self.led.snapshot(N, self.gamma_m, state.get("istep", 0), vmag, span_frac)

# ----------------------------------------------------------------------------------------
# Builders
# ----------------------------------------------------------------------------------------
def add_pinned_source(center, burg, n, nodes, segs, Lseg, line_ref=None):
    """Frank-Read-like pinned-end finite source (verified add_segment / build_src idiom).
    line_ref: Burgers used to ORIENT the line (default = burg). For a collinear pair pass a
    COMMON reference (e.g. bm) for both arms so that burg=-bm gives an annihilating (+b,xi)/(-b,xi)
    crossing instead of a flipped parallel line (binary_observe.py line-sense lesson)."""
    lr = burg if line_ref is None else line_ref
    bh = hat(lr)
    e = np.cross(n, lr); e = hat(e)
    xi = hat(np.cos(np.radians(45.0)) * bh + np.sin(np.radians(45.0)) * e)
    nn = max(4, int(round(Lseg / MAXSEG)))
    i0 = len(nodes); nhat = hat(n)
    for j in range(nn + 1):
        p = center + (j / nn - 0.5) * Lseg * xi
        con = int(NodeConstraints.PINNED_NODE) if j in (0, nn) else int(NodeConstraints.UNCONSTRAINED)
        nodes.append(np.concatenate((p, [con])))
    for j in range(nn):
        segs.append(np.concatenate(([i0 + j, i0 + j + 1], burg, nhat)))

def nforest_for_rho(cell):
    """Choose NFOREST so the as-built forest density ~ RHO_F. Each infinite line ~ one box
    period (~ LBOX*b). rho ~ NFOREST*LBOX*b / (LBOX*b)^3 = NFOREST / (LBOX*b)^2."""
    if NFOREST > 0:
        return NFOREST
    perline = 1.0 / (LBOX * B_CU) ** 2            # rough density per single wrapping line
    return max(2, int(round(RHO_F / perline)))

def build():
    cell = pyexadis.Cell(h=LBOX * np.eye(3), is_periodic=[True, True, True])
    C = np.array(cell.center())
    nodes, segs = [], []
    forest_ranges = []
    nF = 0 if FOREST_ALONE is False and False else nforest_for_rho(cell)
    nF = nforest_for_rho(cell)
    for _ in range(nF):
        o = C + (rng.random(3) - 0.5) * 0.85 * LBOX
        th = float(rng.choice([0, 30, 60, 90]))
        i0 = len(nodes)
        try:
            ok = insert_infinite_line(cell, nodes, segs, bf, nf, o, theta=th, maxseg=MAXSEG, trial=True)
        except TypeError:
            ok = 1.0                              # older signature: no trial kwarg
        if ok and ok > 0:
            insert_infinite_line(cell, nodes, segs, bf, nf, o, theta=th, maxseg=MAXSEG)
            forest_ranges.append((i0, len(nodes)))
    # EVOLVABLE forest: do NOT re-pin.  Optional ONE light anchor per line (suppresses rigid
    # drift only; interior stays free to react).  Gated by ANCHOR env + forest-alone control.
    nodes = np.array(nodes) if len(nodes) else np.zeros((0, 4))
    if ANCHOR and len(forest_ranges):
        for a, b in forest_ranges:
            # anchor the single node nearest a box face
            seg = nodes[a:b, :3]
            face_dist = np.min(np.minimum(seg, LBOX - seg), axis=1)
            nodes[a + int(np.argmin(face_dist)), 3] = int(NodeConstraints.PINNED_NODE)
    nodes = list(nodes)
    # mobile carriers (skip if forest-alone control)
    if not FOREST_ALONE:
        if MOBILE_MODE == "lines":
            # NON-multiplying mobile: infinite gliding lines on system m. Their density changes
            # ONLY by reaction with the forest (annihilation/junction), so opp-vs-same isolates
            # collinear consumption with NO FR-multiplication confound. line sense via common ref bm.
            placed = 0
            for k in range(KM):
                o = C + (rng.random(3) - 0.5) * 0.7 * LBOX
                th = float(rng.choice([30, 60]))   # commensurate (like forest) AND off-screw so they glide
                try:
                    ok = insert_infinite_line(cell, nodes, segs, bm, nm, o, theta=th, maxseg=MAXSEG, trial=True)
                except TypeError:
                    ok = 1.0
                if ok and ok > 0:
                    insert_infinite_line(cell, nodes, segs, bm, nm, o, theta=th, maxseg=MAXSEG)
                    placed += 1
            if placed == 0:
                raise RuntimeError("MOBILE_MODE=lines placed 0 mobile lines (theta non-commensurate?)")
        else:
            # FR pinned-end sources (multiply under stress)
            Lseg = 0.25 * LBOX                        # SHORTER than 0.5 LBOX -> avoid span-wrap
            for k in range(KM):
                c = C + (rng.random(3) - 0.5) * 0.6 * LBOX
                add_pinned_source(c, bm, nm, nodes, segs, Lseg)
    nodes = np.array(nodes); segs = np.array(segs)
    return cell, nodes, segs

def modules(state, cell):
    cf = CalForce(force_mode=FORCE, state=state, Ngrid=NGRID, cell=cell)
    mob = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    ti = TimeIntegration(integrator="Trapezoid", state=state, force=cf, mobility=mob)
    col = Collision(collision_mode="Retroactive", state=state) if COLL else None
    topo = Topology(topology_mode="TopologyParallel", state=state, force=cf, mobility=mob) if TOPO else None
    rm = Remesh(remesh_rule="LengthBased", state=state)
    return cf, mob, ti, col, topo, rm

# ----------------------------------------------------------------------------------------
# PAIR_SMOKE: 2-segment binary_observe-style sense check (NO forest population).
#   confirms coll_opp annihilates (length_fraction < 0.85) and coll_same does NOT (~1.0)
#   at the chosen plane geometry/rann BEFORE any population run.
# ----------------------------------------------------------------------------------------
def run_pair_smoke():
    pyexadis.initialize(); os.makedirs(OUT, exist_ok=True)
    cell = pyexadis.Cell(h=LBOX * np.eye(3), is_periodic=[True, True, True])
    C = np.array(cell.center())
    neutral = hat(hat(nm) + hat(nf))
    nodes, segs = [], []
    GAP, LSEG = 6.0, 1500.0   # GAP<=rann so the two segments actually contact (binary_collinear.py lesson)
    # collinear pair shares the Burgers line family -> orient BOTH arms by a common ref (bm) so
    # bf=-bm yields an annihilating crossing, not a flipped parallel line.
    lref = bm if cos_par(bm, bf) > 0.94 else None
    add_pinned_source(C + 0.5 * GAP * neutral, bm, nm, nodes, segs, LSEG, line_ref=lref)
    add_pinned_source(C - 0.5 * GAP * neutral, bf, nf, nodes, segs, LSEG, line_ref=lref)
    nodes = np.array(nodes); segs = np.array(segs)
    L0 = float(sum(np.linalg.norm(nodes[int(a), :3] - nodes[int(c), :3]) for a, c in segs[:, :2]))
    net = DisNetManager(ExaDisNet(cell, nodes, segs))
    state = dict(crystal="fcc", burgmag=B_CU, mu=MU, nu=0.324, a=6.0, maxseg=MAXSEG,
                 minseg=MAXSEG / 4, rtol=10.0, rann=RANN, nextdt=1e-11, maxdt=1e-9)
    cf, mob, ti, col, topo, rm = modules(state, net.cell)
    sim = SimulateNetwork(calforce=cf, mobility=mob, timeint=ti, collision=col, topology=topo,
                          remesh=rm, vis=None, state=state, max_step=NREL, loading_mode="stress",
                          applied_stress=np.zeros(6), print_freq=10**9, plot_freq=10**9,
                          write_freq=10**9, write_dir=OUT)
    sim.run(net, state)
    d = net.get_disnet(ExaDisNet).export_data()
    pos, nid = d["nodes"]["positions"], np.asarray(d["segs"]["nodeids"]).astype(int)
    Lf = float(sum(np.linalg.norm(closest_vec(net.cell, pos[a], pos[c])) for a, c in nid))
    deg = {}
    for a, c in nid:
        deg[int(a)] = deg.get(int(a), 0) + 1; deg[int(c)] = deg.get(int(c), 0) + 1
    njun = sum(1 for v in deg.values() if v >= 3)
    frac = Lf / L0 if L0 else 0.0
    mech = ("partial_annihilation" if frac < 0.85 else
            "junction" if njun > 0 else
            "pass_through" if abs(frac - 1) < 0.06 else "other")
    out = dict(mode="pair_smoke", jtype=JTYPE, bm=bm.tolist(), nm=nm.tolist(),
               bf=bf.tolist(), nf=nf.tolist(), length_fraction=frac,
               n_junction=njun, mechanism=mech)
    json.dump(out, open(os.path.join(OUT, "pair_smoke.json"), "w"), indent=1)
    print("PAIR_SMOKE:", json.dumps(out), flush=True)
    pyexadis.finalize()

# ----------------------------------------------------------------------------------------
# MAIN production run: settle (zero stress) -> load (stress on m only) -> ledger -> fit.
# ----------------------------------------------------------------------------------------
def main():
    if PAIR_SMOKE:
        return run_pair_smoke()
    pyexadis.initialize(); os.makedirs(OUT, exist_ok=True)
    cell, nodes, segs = build()
    net = DisNetManager(ExaDisNet(cell, nodes, segs))
    state = dict(crystal="fcc", burgmag=B_CU, mu=MU, nu=0.324, a=6.0, maxseg=MAXSEG,
                 minseg=MAXSEG / 4, rtol=10.0, rann=RANN, nextdt=1e-11, maxdt=1e-9)
    rho_f_built = dislocation_density(net, B_CU)   # whole-net density right after build

    # PHASE A: zero-stress settle (let forest + sources relax; DO NOT record gamma)
    cfA, mobA, tiA, colA, topoA, rmA = modules(state, net.cell)
    SimulateNetwork(calforce=cfA, mobility=mobA, timeint=tiA, collision=colA, topology=topoA,
                    remesh=rmA, vis=None, state=state, max_step=NREL, loading_mode="stress",
                    applied_stress=np.zeros(6), print_freq=10**9, plot_freq=10**9,
                    write_freq=10**9, write_dir=OUT).run(net, state)
    rho_f_settled = dislocation_density(net, B_CU)

    if FOREST_ALONE:
        # forest self-evolution control: report intrinsic density drift, no loading.
        cfB, mobB, tiB, colB, topoB, rmB = modules(state, net.cell)
        led = Ledger(net.cell)
        sim = MFPSim(calforce=cfB, mobility=mobB, timeint=tiB, collision=colB, topology=topoB,
                     remesh=rmB, vis=None, state=state, max_step=NSTEPS, loading_mode="stress",
                     applied_stress=np.zeros(6), print_freq=REC, plot_freq=10**9,
                     write_freq=10**9, write_dir=OUT)
        sim.attach(led); sim.run(net, state)
        rho0 = led.rows[0]["rho_forest"] if led.rows else 0.0
        rhoT = led.rows[-1]["rho_forest"] if led.rows else 0.0
        out = dict(mode="forest_alone", jtype=JTYPE, seed=SEED, rho_f_built=rho_f_built,
                   rho_f_settled=rho_f_settled, rho_forest_0=rho0, rho_forest_end=rhoT,
                   forest_drift_frac=(rhoT - rho0) / rho0 if rho0 else 0.0,
                   rows=led.rows)
        json.dump(out, open(os.path.join(OUT, "ledger.json"), "w"), indent=1)
        print("FOREST_ALONE drift=%.3f" % out["forest_drift_frac"], flush=True)
        pyexadis.finalize(); return

    # PHASE B: stress on m only (resolved shear tau); accumulate gamma_m + ledger
    sig = TAU * (np.outer(hat(bm), hat(nm)) + np.outer(hat(nm), hat(bm)))
    voigt = np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[1, 2], sig[0, 2], sig[0, 1]])
    cfB, mobB, tiB, colB, topoB, rmB = modules(state, net.cell)
    led = Ledger(net.cell)
    sim = MFPSim(calforce=cfB, mobility=mobB, timeint=tiB, collision=colB, topology=topoB,
                 remesh=rmB, vis=None, state=state, max_step=NSTEPS, loading_mode="stress",
                 applied_stress=voigt, print_freq=REC, plot_freq=10**9,
                 write_freq=10**9, write_dir=OUT)
    sim.attach(led); sim.run(net, state)

    # -------- ANALYSIS (in-script): fit storage and removal slopes SEPARATELY --------
    R = led.rows
    g = np.array([r["gamma_m"] for r in R]) if R else np.array([])
    res = analyze(R, g, rho_f_built, rho_f_settled, sim)
    res.update(dict(jtype=JTYPE, sense=SENSE, tau_MPa=TAU / 1e6, LBOX=LBOX,
                    nforest=int(nforest_for_rho(cell)), KM=KM, seed=SEED,
                    force=FORCE, topo=TOPO, coll=COLL, anchor=ANCHOR,
                    ngrid=NGRID, rann=RANN, span_wrap_flagged=sim.flagged_wrap,
                    rho_f_built=rho_f_built, rho_f_settled=rho_f_settled,
                    bm=bm.tolist(), nm=nm.tolist(), bf=bf.tolist(), nf=nf.tolist(),
                    gamma_m_final=float(g[-1]) if len(g) else 0.0,
                    gamma_m_global_final=float(sim.gamma_m_global)))
    res["rows"] = R
    json.dump(res, open(os.path.join(OUT, "ledger.json"), "w"), indent=1)
    print("RESULT %s s%d tau%.0f: S_arrest=%.2e(R2=%.2f) S_store=%.2e S_remove=%.2e "
          "L_mf=%.0fb Lambda=%.3f fit_n=%d win_span=%.2f wrap=%s" %
          (JTYPE, SEED, TAU / 1e6, res["S_arrest"], res["R2_arrest"], res["S_store"],
           res["S_remove"], res["L_mf_b"], res["Lambda"], res["fit_n"],
           res["win_span_max"], sim.flagged_wrap), flush=True)
    pyexadis.finalize()

def linfit(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 4 or np.ptp(x) <= 0:
        return float("nan"), float("nan"), 0.0
    try:
        s, b = np.polyfit(x, y, 1)
    except Exception:
        return float("nan"), float("nan"), 0.0
    yh = b + s * x
    ss = 1 - np.sum((y - yh) ** 2) / max(1e-30, np.sum((y - y.mean()) ** 2))
    return float(s), float(b), float(ss)

def analyze(R, g, rho_f_built, rho_f_settled, sim):
    if not R:
        return dict(S_store=float("nan"), S_remove=float("nan"), L_mf_b=float("inf"),
                    Lambda=float("inf"), R2_store=0.0, R2_remove=0.0, fit_n=0,
                    gamma_global_consistency=float("nan"))
    surv = np.array([r["mobile_survival"] for r in R])
    span = np.array([r.get("mobile_span_frac", 0.0) for r in R])
    stored = np.array([r["rho_stored"] for r in R])           # junction + residual (RETAINED arrest)
    removed = np.array([r["rho_removed_balance"] for r in R])  # junction-SUBTRACTED removal (annihilation)
    arrest = stored + removed                                  # total mobile arrest (the MFP driver)
    rhoF = np.array([r["rho_forest"] for r in R])
    # PRE-WRAP window only: once the mobile glide front wraps the PBC box the MFP is corrupted.
    gmin = g[max(1, len(g) // 8)] if len(g) else 0.0
    win = (span < 0.70) & (g > gmin) & np.r_[True, np.diff(g) >= 0]
    if win.sum() < 4:
        win = (span < 0.85) & np.r_[True, np.diff(g) >= 0]
    if win.sum() < 4:
        win = np.r_[True, np.diff(g) >= 0]
    gw = g[win]
    S_store, c_store, R2_store = linfit(gw, stored[win])
    S_remove, c_remove, R2_remove = linfit(gw, removed[win])
    S_arrest, c_arrest, R2_arrest = linfit(gw, arrest[win])
    S_tot = S_arrest if S_arrest == S_arrest else 0.0
    L_mf_b = (1.0 / (B_CU * (S_tot / B_CU))) if S_tot > 0 else float("inf")
    # NOTE: rho is in 1/m^2, gamma dimensionless -> S has units 1/m^2. L = 1/(b*S) in meters.
    L_mf_m = 1.0 / (B_CU * S_tot) if S_tot > 0 else float("inf")
    L_mf_b = L_mf_m / B_CU if np.isfinite(L_mf_m) else float("inf")
    rho_f_fit = float(np.mean(rhoF[win])) if win.sum() else (rho_f_settled or rho_f_built)
    Lambda = L_mf_m * np.sqrt(rho_f_fit) if np.isfinite(L_mf_m) and rho_f_fit > 0 else float("inf")
    # gamma cross-check (global should be >= mobile-only; flag if mobile >> global)
    gg = float(sim.gamma_m_global)
    cons = (g[-1] / gg) if gg != 0 else float("nan")
    win_span_max = float(np.max(span[win])) if win.sum() else 1.0
    return dict(S_store=float(S_store), S_remove=float(S_remove), S_arrest=float(S_arrest),
                S_tot=float(S_tot), c_store=float(c_store), c_remove=float(c_remove),
                R2_store=float(R2_store), R2_remove=float(R2_remove), R2_arrest=float(R2_arrest),
                L_mf_b=float(L_mf_b), L_mf_m=float(L_mf_m), Lambda=float(Lambda),
                rho_f_fit=float(rho_f_fit), fit_n=int(win.sum()), win_span_max=win_span_max,
                gamma_global_consistency=float(cons),
                mean_ambiguous_frac=float(np.mean([r["ambiguous_frac"] for r in R])))

if __name__ == "__main__":
    main()

"""Multi-slip COLL on/off DECISION cell -- collinear-dominance test (upgraded per expert review).

Upgrades over the v1 scaffold (so it can DECIDE, not just run):
  (1) CO-DRIVEN loading: EDIR_MODE=opt_pair finds the uniaxial axis maximizing min(|S_m|,|S_f|)
      so BOTH the primary and the partner system are near-equally driven (collinear carriers
      replenish bidirectionally). v1's EDIR drove the partner at only ~1/3 of primary -> a
      coll_opp==coll_same there was inconclusive. Schmid factors are output + gated.
  (2) PER-SYSTEM density ledger: rho_mobile_primary / rho_forest_partner / rho_junction /
      rho_residual / rho_removed_proxy / ambiguous_frac, classified by plane (primary for
      collinear, since it shares the Burgers family) + Burgers + node degree. forest_drift is the
      PARTNER-forest drift (not total density, which FR multiplication contaminates).
  (3) GATES in JSON: plateau (plastic gamma>1e-3 & d(tau)/d(gamma)->0), forest_drift<5%,
      ambiguous_frac<0.1 -> readable. Non-readable runs are 'diagnostic', never a conclusion.
  Baselines: RHO_F=0 (source-only), FOREST_ALONE=1, XSLIP=0 ablation.

  JTYPE=coll_opp EDIR_MODE=opt_pair RHO_F=3e12 ERATE=1e5 LBOX=10000 NSTEPS=12000 SEED=1 \
    OUT=ms/coll_opp_s1 PYTHONPATH=~/BO/exadis_src/python python3 multislip_flow.py
"""
import os, sys, json
import numpy as np

EX = os.path.expanduser("~/BO/exadis_src/python")
if EX not in sys.path:
    sys.path.append(EX)
import pyexadis
from pyexadis_base import (ExaDisNet, DisNetManager, SimulateNetwork, CalForce,
                           MobilityLaw, TimeIntegration, Collision, Topology, Remesh,
                           CrossSlip, NodeConstraints)
from pyexadis_utils import insert_infinite_line, insert_frank_read_src, dislocation_density

BI = np.array([[0, 1, -1], [1, 0, -1], [1, -1, 0], [0, 1, -1], [1, 0, 1], [1, 1, 0],
               [0, 1, 1], [1, 0, -1], [1, 1, 0], [0, 1, 1], [1, 0, 1], [1, -1, 0]], float)
NI = np.array([[1, 1, 1]] * 3 + [[-1, 1, 1]] * 3 + [[1, -1, 1]] * 3 + [[1, 1, -1]] * 3, float)
B_CU, MU, NU = 2.55e-10, 54.6e9, 0.324


def envf(k, d): return float(os.environ.get(k, d))
def envi(k, d): return int(os.environ.get(k, d))


JTYPE  = os.environ.get("JTYPE", "coll_opp")
EDIR_MODE = os.environ.get("EDIR_MODE", "opt_pair")   # 'opt_pair' (co-drive) | 'primary'
LBOX   = envf("LBOX", "10000")
RHO_F  = envf("RHO_F", "3e12")             # 0 -> no forest (source-only baseline)
NFOREST = envi("NFOREST", "-1")
ERATE  = envf("ERATE", "1e5")
KM     = envi("KM", "14")
LSRC   = envf("LSRC", "0")                  # 0 -> 0.35*LBOX (low activation)
NSTEPS = envi("NSTEPS", "12000")
NREL   = envi("NREL", "150")
REC    = envi("REC", "40")
SEED   = envi("SEED", "1")
MAXSEG = envf("MAXSEG", "100")
NGRID  = envi("NGRID", "32")
RANN   = envf("RANN", "10")
XSLIP  = os.environ.get("XSLIP", "1") == "1"
FOREST_ALONE = os.environ.get("FOREST_ALONE", "0") == "1"
OUT    = os.environ.get("OUT", "ms_out")
rng = np.random.default_rng(SEED)


def hat(v):
    v = np.asarray(v, float); n = np.linalg.norm(v)
    return v / n if n > 0 else v
def cos_par(u, v):
    return abs(float(np.dot(hat(u), hat(v))))
def schmid(b, n, e):
    return float(np.dot(hat(b), e) * np.dot(hat(n), e))


def pick_pair():
    bm, nm = BI[0].copy(), NI[0].copy()
    j = JTYPE.lower()
    if j in ("coll_opp", "coll_same"):
        nf = NI[3].copy()
        bf = (-1.0 if j == "coll_opp" else 1.0) * bm.copy()
    elif j == "glissile":
        bf, nf = BI[4].copy(), NI[4].copy()
    elif j == "hirth":
        bf, nf = BI[6].copy(), NI[6].copy()
    elif j == "lomer":
        bf, nf = BI[5].copy(), NI[5].copy()
    else:
        raise ValueError(JTYPE)
    return bm, nm, bf, nf


bm, nm, bf, nf = pick_pair()
b3a, b3b = bm + bf, bm - bf
FAM_J = [v for v in (b3a, b3b) if np.linalg.norm(v) > 1e-6 and cos_par(v, bm) < 0.94 and cos_par(v, bf) < 0.94]
COLLINEAR = cos_par(bm, bf) > 0.94


def opt_edir():
    """Fibonacci-sphere search for the uniaxial axis maximizing min(|S_m|,|S_f|)."""
    N = 4000; idx = np.arange(N)
    phi = np.pi * (3 - np.sqrt(5)) * idx; z = 1 - 2 * idx / (N - 1); r = np.sqrt(np.clip(1 - z * z, 0, 1))
    pts = np.c_[r * np.cos(phi), r * np.sin(phi), z]
    sm = (pts @ hat(bm)) * (pts @ hat(nm))
    sf = (pts @ hat(bf)) * (pts @ hat(nf))
    return pts[np.argmax(np.minimum(np.abs(sm), np.abs(sf)))]


EDIR = hat(hat(bm) + hat(nm)) if EDIR_MODE == "primary" else hat(opt_edir())
S_PRIM = schmid(bm, nm, EDIR); S_PART = schmid(bf, nf, EDIR)
S_RATIO = abs(S_PART / S_PRIM) if S_PRIM != 0 else 0.0
A_m = 0.5 * (np.outer(hat(bm), hat(nm)) + np.outer(hat(nm), hat(bm)))
A2_m = np.array([A_m[0, 0], A_m[1, 1], A_m[2, 2], 2 * A_m[1, 2], 2 * A_m[0, 2], 2 * A_m[0, 1]])


def nforest():
    if RHO_F <= 0:
        return 0
    if NFOREST >= 0:
        return NFOREST
    perline = 1.0 / (LBOX * B_CU) ** 2
    return max(2, int(round(RHO_F / perline)))


# ---------------- per-system density ledger ----------------
class Ledger:
    def __init__(self, cell):
        self.cell = cell
        self.vol = abs(np.linalg.det(np.array(cell.h)))
        self.rows = []
    def _rho(self, Lb): return Lb / self.vol / B_CU ** 2
    def classify(self, b_i, p_i, anchored):
        if any(cos_par(b_i, r) > 0.94 for r in FAM_J):
            return "junction"
        pl = (p_i is not None) and np.linalg.norm(p_i) > 1e-6
        if COLLINEAR:
            if anchored: return "residual"
            if pl and cos_par(p_i, nm) > 0.97: return "mobile"
            if pl and cos_par(p_i, nf) > 0.97: return "forest"
            return "residual"
        if cos_par(b_i, bm) > 0.94:
            return "residual" if anchored else "mobile"
        if cos_par(b_i, bf) > 0.94:
            return "residual" if anchored else "forest"
        return "offsys"
    def snapshot(self, net, gamma_p, istep, tau, strain=0.0):
        d = net.get_disnet(ExaDisNet).export_data()
        pos = d["nodes"]["positions"]; S = d["segs"]
        nid = np.asarray(S["nodeids"]).astype(int)
        burg = np.asarray(S["burgers"], float) if "burgers" in S else None
        plane = np.asarray(S["planes"], float) if "planes" in S else None
        deg = {}
        for a, c in nid:
            deg[int(a)] = deg.get(int(a), 0) + 1; deg[int(c)] = deg.get(int(c), 0) + 1
        pool = {"mobile": 0., "forest": 0., "junction": 0., "residual": 0., "offsys": 0.}
        for i, (a, c) in enumerate(nid):
            a, c = int(a), int(c)
            L = float(np.linalg.norm(pos[a] - pos[c]))
            b_i = burg[i] if burg is not None else bm
            p_i = plane[i] if plane is not None else None
            cls = self.classify(b_i, p_i, deg.get(a, 0) >= 3 or deg.get(c, 0) >= 3)
            pool[cls] += L
        Ltot = sum(pool.values()) or 1e-30
        row = dict(istep=int(istep), gamma_p=float(gamma_p), strain=float(strain), tau_MPa=float(tau / 1e6),
                   rho_mobile=self._rho(pool["mobile"]), rho_forest=self._rho(pool["forest"]),
                   rho_junction=self._rho(pool["junction"]), rho_residual=self._rho(pool["residual"]),
                   rho_stored=self._rho(pool["junction"] + pool["residual"]),
                   rho_total=self._rho(Ltot), ambiguous_frac=float(pool["offsys"] / Ltot))
        self.rows.append(row)
        return row


class MSSim(SimulateNetwork):
    def attach(self, led):
        self.led = led; self.gamma_p = 0.0
    def step_end(self, N, state):
        dEp = np.asarray(state.get("dEp", np.zeros(6)), float)
        if dEp.shape == (6,):
            self.gamma_p += float(np.dot(dEp, A2_m))
        if (state.get("istep", 0) % REC) == 0:
            self.led.snapshot(N, self.gamma_p, state.get("istep", 0),
                              float(state.get("stress", 0.0)), float(state.get("strain", 0.0)))


def build():
    cell = pyexadis.Cell(h=LBOX * np.eye(3), is_periodic=[True, True, True])
    C = np.array(cell.center()); nodes, segs = [], []
    nF = nforest()
    for _ in range(nF):
        o = C + (rng.random(3) - 0.5) * 0.85 * LBOX
        th = float(rng.choice([0, 30, 60, 90]))
        try:
            ok = insert_infinite_line(cell, nodes, segs, bf, nf, o, theta=th, maxseg=MAXSEG, trial=True)
        except TypeError:
            ok = 1.0
        if ok and ok > 0:
            insert_infinite_line(cell, nodes, segs, bf, nf, o, theta=th, maxseg=MAXSEG)
    nsrc = 0
    if not FOREST_ALONE:
        Lsrc = LSRC if LSRC > 0 else 0.35 * LBOX
        for k in range(KM):
            c = C + (rng.random(3) - 0.5) * 0.7 * LBOX
            try:
                insert_frank_read_src(cell, nodes, segs, bm, nm, Lsrc, c, theta=float(rng.choice([30, 60])))
                nsrc += 1
            except Exception as e:
                if k == 0: print("FR err", e, flush=True)
        if nsrc == 0:
            raise RuntimeError("0 FR sources")
    return cell, np.array(nodes), np.array(segs), nF, nsrc


def modules(state, cell):
    cf = CalForce(force_mode="DDD_FFT_MODEL", state=state, Ngrid=NGRID, cell=cell)
    mob = MobilityLaw(mobility_law="FCC_0", state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
    ti = TimeIntegration(integrator="Trapezoid", state=state, force=cf, mobility=mob)
    col = Collision(collision_mode="Retroactive", state=state)
    topo = Topology(topology_mode="TopologyParallel", state=state, force=cf, mobility=mob)
    rm = Remesh(remesh_rule="LengthBased", state=state)
    xs = CrossSlip(state=state, cross_slip_mode="ForceBasedParallel", force=cf, mobility=mob) if XSLIP else None
    return cf, mob, ti, col, topo, rm, xs


def fwd(s, e): return float(s) > 0


def main():
    pyexadis.initialize(); os.makedirs(OUT, exist_ok=True)
    cell, nodes, segs, nF, nsrc = build()
    net = DisNetManager(ExaDisNet(cell, nodes, segs))
    state = {"crystal": "fcc", "burgmag": B_CU, "mu": MU, "nu": NU, "a": 6.0,
             "maxseg": MAXSEG, "minseg": MAXSEG / 4, "rtol": 10.0, "rann": RANN,
             "nextdt": 1e-12, "maxdt": 1e-10}
    led = Ledger(net.cell)
    cf, mob, ti, col, topo, rm, xs = modules(state, net.cell)
    SimulateNetwork(calforce=cf, mobility=mob, timeint=ti, collision=col, topology=topo, remesh=rm,
                    cross_slip=xs, vis=None, state=state, burgmag=B_CU, loading_mode="stress",
                    applied_stress=np.zeros(6), max_step=NREL, print_freq=10**9, plot_freq=10**9,
                    write_freq=10**9, write_dir=OUT).run(net, state)
    settled = led.snapshot(net, 0.0, 0, 0.0)

    cf2, mob2, ti2, col2, topo2, rm2, xs2 = modules(state, net.cell)
    state["edir"] = EDIR.copy()
    sim = MSSim(calforce=cf2, mobility=mob2, timeint=ti2, collision=col2, topology=topo2, remesh=rm2,
                cross_slip=xs2, vis=None, state=state, burgmag=B_CU, loading_mode="strain_rate",
                erate=ERATE, edir=EDIR.copy(), max_step=NSTEPS, print_freq=REC, plot_freq=10**9,
                write_freq=10**9, write_dir=OUT)
    sim.attach(led); sim.run(net, state)

    R = led.rows
    eps = np.array([r.get("strain", 0.0) for r in R]); tau = np.array([r["tau_MPa"] for r in R])
    rf = np.array([r["rho_forest"] for r in R]); amb = np.array([r["ambiguous_frac"] for r in R])
    rstore = np.array([r["rho_stored"] for r in R]); rmob = np.array([r["rho_mobile"] for r in R])
    Eyoung = 2 * MU * (1 + NU) / 1e6                  # MPa
    eps_p = eps - tau / Eyoung                         # plastic strain = total - elastic
    eps_p = np.maximum.accumulate(np.clip(eps_p, 0, None)) if len(eps_p) else eps_p
    # plateau = yielding: d(tau)/d(eps) falls well below E (plastic flow caps the stress)
    plateau = False; dtau_de_over_E = float("nan")
    if len(eps) >= 8 and np.ptp(eps) > 0:
        late = slice(max(1, 2 * len(eps) // 3), len(eps))
        if np.ptp(eps[late]) > 0:
            dtau_de = float(np.polyfit(eps[late], tau[late], 1)[0])
            dtau_de_over_E = dtau_de / Eyoung
            plateau = bool(eps_p.max() > 1e-4 and abs(dtau_de) < 0.20 * Eyoung)
    # flow stress = mean tau over the post-yield window
    yld = eps_p > max(0.4 * eps_p.max(), 1e-5) if len(eps_p) else np.zeros(0, bool)
    tau_flow = float(np.mean(tau[yld])) if yld.any() else float(tau[-1] if len(tau) else 0)
    # forest-partner drift (settled -> end), the RIGHT density-drift gate
    rf_settled = settled["rho_forest"]; rf_end = R[-1]["rho_forest"] if R else 0.0
    forest_drift = float((rf_end - rf_settled) / rf_settled) if rf_settled > 1e6 else 0.0
    mean_amb = float(np.mean(amb)) if len(amb) else 1.0
    # storage MFP from stored-density rise per PLASTIC strain
    L_mf_b = float("inf"); S_store = float("nan")
    if yld.sum() >= 4 and np.ptp(eps_p[yld]) > 0:
        S_store = float(np.polyfit(eps_p[yld], rstore[yld], 1)[0])
        if S_store > 0:
            L_mf_b = (1.0 / (B_CU * S_store)) / B_CU
    g = eps_p   # for downstream reporting (plastic strain)

    drift_gate = abs(forest_drift) < 0.05
    amb_gate = mean_amb < 0.10
    schmid_gate = (S_RATIO > 0.8) if COLLINEAR else (S_RATIO > 0.5)
    readable = bool(plateau and drift_gate and amb_gate)

    out = dict(jtype=JTYPE, edir_mode=EDIR_MODE, edir=EDIR.tolist(),
               schmid_primary=S_PRIM, schmid_partner=S_PART, schmid_ratio=S_RATIO, schmid_gate=schmid_gate,
               rho_f_target=RHO_F, nforest=nF, nsrc=nsrc, erate=ERATE, lbox=LBOX, xslip=XSLIP,
               forest_alone=FOREST_ALONE, seed=SEED, collinear=COLLINEAR,
               eps_p_max=float(g.max()) if len(g) else 0.0, strain_max=float(eps.max()) if len(eps) else 0.0,
               strain_reached=bool(len(g) and g.max() > 1e-4),
               tau_flow_MPa=tau_flow, tau_RSS_MPa=float(tau_flow * abs(S_PRIM)),
               dtau_deps_over_E=dtau_de_over_E if dtau_de_over_E == dtau_de_over_E else None,
               rho_mobile_starved=bool(len(rmob) and rmob[-1] < 0.3 * (max(rmob) if len(rmob) else 1)),
               plateau_reached=plateau, forest_drift=forest_drift, mean_ambiguous_frac=mean_amb,
               rho_mobile_end=float(rmob[-1]) if len(rmob) else 0.0,
               rho_forest_settled=rf_settled, rho_forest_end=rf_end,
               rho_stored_end=float(rstore[-1]) if len(rstore) else 0.0,
               S_store=S_store, L_mf_b=L_mf_b,
               drift_gate=drift_gate, ambiguous_gate=amb_gate, readable=readable,
               series=[[r.get("strain", 0.0), r["tau_MPa"], r["rho_mobile"], r["rho_forest"], r["rho_stored"]] for r in R[::2]])
    json.dump(out, open(os.path.join(OUT, "flow.json"), "w"), indent=1)
    print(f"FLOW {JTYPE} s{SEED} edir={EDIR_MODE}(Sm={S_PRIM:.2f},Sf={S_PART:.2f},r={S_RATIO:.2f}) "
          f"rho_f={RHO_F}: tau_flow={tau_flow:.1f}MPa eps_p={out['eps_p_max']:.2e} eps={out['strain_max']:.2e} "
          f"plateau={plateau} drift={forest_drift:+.2f} L_mf={L_mf_b:.0f}b READABLE={readable}", flush=True)
    pyexadis.finalize()


if __name__ == "__main__":
    main()

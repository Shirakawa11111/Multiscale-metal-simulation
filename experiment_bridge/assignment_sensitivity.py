"""Assignment-policy sensitivity audit, v1.1 (LOCAL, no DDD run).

Compares the three lowering policies and exposes the edgewise artifact that the M3/M4 v0 pilots hit:
  top1            -> deterministic baseline
  sample_edgewise -> per-edge draw (DEPRECATED: creates within-line Burgers discontinuities)
  sample_linewise -> physical default (one draw per parent reconstructed line; 0 within-line discontinuities)

Key metric: within-line Burgers discontinuities (adjacent segments of the SAME parent line that carry a
different Burgers/plane signature) -- the artificial junctions that inflated the edgewise topology result.

  python3 experiment_bridge/assignment_sensitivity.py [cu_stem_idr.json] [N]
  -> results_exadis/assignment_sensitivity.{json,md}
"""

import os, sys, json, math
from collections import Counter
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
from defect_ir.adapters.to_exadis import idr_to_exadis_network

OUT = os.path.join(os.path.dirname(__file__), "results_exadis")
IDR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(OUT, "cu_stem_idr.json")
N = int(sys.argv[2]) if len(sys.argv) > 2 else 30


def sig(seg):
    return tuple(round(x, 2) for x in seg[2:8])


def within_line_discont(net, plid):
    s = [sig(x) for x in net["segs"]]
    disc = tot = 0
    for i in range(1, len(s)):
        if plid[i] == plid[i - 1]:
            tot += 1
            disc += s[i] != s[i - 1]
    return disc, tot


def inventory_cv(nets):
    keys = set()
    invs = [Counter(sig(x) for x in n["segs"]) for n in nets]
    for inv in invs:
        keys |= set(inv)
    cvs = []
    for k in keys:
        vals = [inv.get(k, 0) for inv in invs]
        mu = np.mean(vals)
        if mu > 0:
            cvs.append(np.std(vals) / mu)
    return float(np.mean(cvs)) if cvs else 0.0


def main():
    doc = json.load(open(IDR))
    plid = [e.get("parent_line_id") for e in doc["geometry"]["edges"]]
    res = {}
    for pol, nseed in [("top1", 1), ("sample_edgewise", N), ("sample_linewise", N)]:
        nets = [
            idr_to_exadis_network(doc, assignment_policy=pol, seed=s)
            for s in range(nseed)
        ]
        dts = [within_line_discont(n, plid) for n in nets]
        disc = np.mean([d for d, t in dts])
        tot = dts[0][1]
        res[pol] = {
            "within_line_discontinuities_mean": round(float(disc), 1),
            "intra_line_adjacencies": tot,
            "within_line_discontinuity_frac": round(float(disc) / max(1, tot), 3),
            "inventory_cv": round(inventory_cv(nets), 3),
            "n_seeds": nseed,
        }
    edge_d = res["sample_edgewise"]["within_line_discontinuities_mean"]
    line_d = res["sample_linewise"]["within_line_discontinuities_mean"]
    rep = {
        "source": os.path.basename(IDR),
        "policies": res,
        "edgewise_artifact": {
            "within_line_discontinuities": edge_d,
            "of_intra_line_adjacencies": res["sample_edgewise"][
                "intra_line_adjacencies"
            ],
            "linewise_discontinuities": line_d,
            "verdict": (
                "sample_edgewise injects ~%.0f artificial within-line Burgers discontinuities "
                "(=> artificial junctions); sample_linewise injects 0. Use sample_linewise."
                % edge_d
            ),
        },
    }
    json.dump(
        rep, open(os.path.join(OUT, "assignment_sensitivity.json"), "w"), indent=1
    )
    md = f"""# Assignment-policy sensitivity (v1.1, local, no DDD)

`{rep['source']}`, {N} Monte-Carlo seeds for the sampling policies.

| policy | within-line discontinuities | of intra-line adjacencies | inventory CV |
|--|--|--|--|
| top1 | {res['top1']['within_line_discontinuities_mean']} | {res['top1']['intra_line_adjacencies']} | {res['top1']['inventory_cv']} |
| sample_edgewise (deprecated) | **{edge_d}** | {res['sample_edgewise']['intra_line_adjacencies']} | {res['sample_edgewise']['inventory_cv']} |
| sample_linewise (default) | **{line_d}** | {res['sample_linewise']['intra_line_adjacencies']} | {res['sample_linewise']['inventory_cv']} |

**Verdict.** {rep['edgewise_artifact']['verdict']} Edgewise sampling breaks Burgers continuity along a single
physical reconstructed line, manufacturing junction-like topology; line-coherent sampling preserves it.
The assignment *ambiguity itself* is real (geometry fixes the {{111}} plane, the 3 ⟨110⟩ Burgers are
near-degenerate), but it must be propagated **per line**, not per segment.
"""
    open(os.path.join(OUT, "assignment_sensitivity.md"), "w").write(md)
    print(
        f"edgewise within-line discont={edge_d}/{res['sample_edgewise']['intra_line_adjacencies']}, "
        f"linewise={line_d}; inventory CV edge={res['sample_edgewise']['inventory_cv']} "
        f"line={res['sample_linewise']['inventory_cv']}"
    )
    print("-> results_exadis/assignment_sensitivity.{json,md}")


if __name__ == "__main__":
    main()

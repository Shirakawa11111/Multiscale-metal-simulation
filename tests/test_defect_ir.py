"""Gate test for the defect_ir package: examples build + validate, and lowering round-trips."""

import os, sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from defect_ir.validators import validate_idr
from defect_ir.adapters.to_exadis import idr_to_exadis_network
from defect_ir.examples.build_examples import build_cu_stem, build_graphene


def test_examples_valid():
    for doc in (build_cu_stem(), build_graphene()):
        ok, errs, warns = validate_idr(doc)
        assert ok, f"invalid IDR: {errs}"


def test_lowering_roundtrip():
    cu = build_cu_stem()
    net = idr_to_exadis_network(cu, assignment_policy="top1")
    assert net["network_counts"]["segments"] == len(cu["geometry"]["edges"])
    assert net["cell"]["is_periodic"] == [True, True, False]
    netp = idr_to_exadis_network(cu, cell_policy="thickened_periodic", zbox=5.0)
    assert netp["cell"]["is_periodic"] == [True, True, True]


def test_invalid_rejected():
    bad = build_cu_stem()
    bad["geometry"]["edges"][0]["v1"] = 99999
    ok, errs, _ = validate_idr(bad)
    assert not ok and any("not a known vertex" in e for e in errs)


def test_linewise_coherent():
    from defect_ir.adapters.to_exadis import idr_to_exadis_network

    cu = build_cu_stem()
    plid = [e["parent_line_id"] for e in cu["geometry"]["edges"]]
    net = idr_to_exadis_network(cu, assignment_policy="sample_linewise", seed=1)
    sig = [tuple(round(x, 2) for x in s[2:8]) for s in net["segs"]]
    for i in range(1, len(sig)):
        if plid[i] == plid[i - 1]:
            assert (
                sig[i] == sig[i - 1]
            ), "linewise must keep one Burgers per parent line"


if __name__ == "__main__":
    test_examples_valid()
    test_lowering_roundtrip()
    test_invalid_rejected()
    test_linewise_coherent()
    test_gb_collapse()
    print("defect_ir gate test: PASS")


def test_gb_collapse():
    from defect_ir.uncertainty import apply_gb_constraints, assignment_entropy

    cu = build_cu_stem()
    cands = cu["topology"]["edge_labels"][0]["slip_system_candidates"]
    true_b = cands[0]["b"]
    import numpy as np

    obs = [
        {"g": g, "visible": abs(np.dot(g, true_b)) > 1e-6}
        for g in ([2, 0, 0], [0, 2, 0], [0, 0, 2])
    ]
    kept, status = apply_gb_constraints(cands, obs)
    assert len(kept) == 1 and status == "gb_validated"
    assert assignment_entropy(kept) == 0.0
    # no observations -> unchanged
    k2, st2 = apply_gb_constraints(cands, [])
    assert len(k2) == len(cands) and st2 == "geometry_only_pending_gb"

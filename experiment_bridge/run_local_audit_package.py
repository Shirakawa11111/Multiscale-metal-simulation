"""One-click LOCAL reproducer for the STEM-to-DDD v2 audit package.

Runs every DETERMINISTIC (non-DDD, non-HPC) step in order and checks its outputs exist:
  stem_to_idr.py          -> cu_stem_idr.json + report           (recon -> IDR)
  assignment_sensitivity  -> assignment_sensitivity.{json,md}    (edgewise 142.8 vs linewise 0)
  density_conventions     -> density_conventions.{json,md}       (Lambda_A foil-native + as-built/relaxed)
  synthetic_gb            -> synthetic_gb.{json,md}               (line-coherent g.b entropy collapse)
  tests/test_defect_ir    -> gate test PASS
  make_audit_figure       -> audit_summary_figure.png            (4-panel main result, if matplotlib)
The DDD/HPC steps (real_network_audit, cell_policy_audit) are NOT run here; their authoritative records are
the summary JSONs + AUDIT_MANIFEST.md. This script verifies those summaries are present, not that they rerun.

  python3 experiment_bridge/run_local_audit_package.py   -> exit 0 iff every local step + output is OK
"""
import os, sys, subprocess, json

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(HERE, "results_exadis")
PY = sys.executable


def run(label, argv, outputs):
    print(f"\n=== {label} ===")
    r = subprocess.run([PY] + argv, cwd=ROOT, capture_output=True, text=True)
    tail = (r.stdout.strip().splitlines() or [""])[-1]
    print(f"  $ {' '.join(argv)}\n  -> {tail}")
    if r.returncode != 0:
        print(f"  !! FAILED (rc={r.returncode})\n{r.stderr[-800:]}")
        return False
    missing = [o for o in outputs if not os.path.exists(os.path.join(OUT, o))]
    if missing:
        print(f"  !! MISSING OUTPUTS: {missing}")
        return False
    print(f"  OK ({len(outputs)} outputs present)")
    return True


def check_present(label, files):
    print(f"\n=== {label} (presence-only; authoritative DDD/HPC records) ===")
    missing = [f for f in files if not os.path.exists(os.path.join(OUT, f))]
    if missing:
        print(f"  !! MISSING: {missing}")
        return False
    print(f"  OK ({len(files)} summary files present)")
    return True


def main():
    steps = [
        ("STEM -> IDR", ["experiment_bridge/stem_to_idr.py"],
         ["cu_stem_idr.json", "cu_stem_idr_report.json", "cu_stem_idr_report.md"]),
        ("assignment sensitivity (edgewise vs linewise)", ["experiment_bridge/assignment_sensitivity.py"],
         ["assignment_sensitivity.json", "assignment_sensitivity.md"]),
        ("density conventions (as-built + relaxed)", ["experiment_bridge/density_conventions.py"],
         ["density_conventions.json", "density_conventions.md"]),
        ("synthetic g.b (line-coherent)", ["experiment_bridge/synthetic_gb.py"],
         ["synthetic_gb.json", "synthetic_gb.md"]),
        ("gate test", ["tests/test_defect_ir.py"], []),
    ]
    ok = all(run(*s) for s in steps)
    # optional figure (needs matplotlib)
    fig = os.path.join(HERE, "make_audit_figure.py")
    if os.path.exists(fig):
        ok = run("main-result figure", ["experiment_bridge/make_audit_figure.py"],
                 ["audit_summary_figure.png"]) and ok
    # DDD/HPC summaries: presence-only
    ok = check_present("DDD/HPC audit records",
                       ["cell_policy_audit_summary.json"]) and ok

    print("\n" + "=" * 60)
    if ok:
        # tiny invariant cross-check: edgewise artifact present, linewise clean
        a = json.load(open(os.path.join(OUT, "assignment_sensitivity.json")))
        ew = a["policies"]["sample_edgewise"]["within_line_discontinuities_mean"]
        lw = a["policies"]["sample_linewise"]["within_line_discontinuities_mean"]
        g = json.load(open(os.path.join(OUT, "synthetic_gb.json")))
        e0 = g["scenarios"]["no_gb"]["mean_entropy_bits"]
        e2 = g["scenarios"]["2_reflections_g200_g020"]["mean_entropy_bits"]
        print(f"LOCAL AUDIT PACKAGE: PASS")
        print(f"  invariants: edgewise within-line discont={ew} (artifact), linewise={lw} (clean); "
              f"g.b entropy {e0}->{e2} bits")
        sys.exit(0)
    print("LOCAL AUDIT PACKAGE: FAIL (see above)")
    sys.exit(1)


if __name__ == "__main__":
    main()

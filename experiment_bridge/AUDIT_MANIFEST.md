# Audit reproducibility manifest

The DDD audits are driven by shell job-lists generated on the HPC and run with a small semaphore
(`real_network_audit.py` per (policy, cell, force, endpoint, seed)); there is no `run_audit_campaign.py`
— this manifest is the authoritative record.

## IDR v1.1 decisive rerun (`v11/`) — top1 vs edgewise vs linewise
- IDR: `results_exadis/cu_stem_idr.json` (27 lines, 270 verts, 243 edges, `parent_line_id` per edge)
- assignment: `top1` (seed 0); `sample_edgewise`, `sample_linewise` (seeds 1–6 each)
- cell/force: `as_is`+LineTension (foil); `thickened_periodic`(zbox=5)+DDD_FFT_MODEL;
  `thickened_periodic`(zbox=5)+LineTension (deconfounding control). `as_is`+DDD_FFT is unsupported (FFT needs PBC).
- endpoint: `pinned`; NREL=300, NLOAD=300; stress σxy=100, σyz=60 MPa; mobility FCC_0; OMP=3, ≤10 concurrent.
- engine: pyexadis at `~/BO/exadis_src/python`. Outputs: `v11/<tag>/audit.json`.
- analysis: `results_exadis/v11_linewise_summary.json`, figure `results_exadis/v11_linewise.png`.

## Earlier pilots (superseded by v1.1 for the topology claim)
- M3 `audit/` (10 runs) and M4 `m4/` (35 runs): used the edgewise `sample` policy; their topology
  numbers are an upper bound (edgewise artifact). Density/survival conclusions stand.

Commit hash and exact job lists are in git history; the generators are the inline `emit()` loops in the
campaign commits (`run_v11.sh`, `run_m4.sh`, `run_audit.sh`).

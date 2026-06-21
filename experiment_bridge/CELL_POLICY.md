# Cell policy — STEM foil → DDD (transparency note)

The reconstructed STEM specimen is a **foil** (z = thin direction, weakly constrained). Different
stages use different cells *on purpose*; this file makes the transformation an explicit, auditable
assumption rather than a hidden one.

| stage | cell | periodicity | why |
|--|--|--|--|
| `reconstruction_cell` | raw STEM/tomography coords, z = foil thickness | z non-periodic | physical foil geometry |
| `exadis_ingestion_cell` (`stem_network.json`, `run_stem_exadis.py`) | `h_b` from ingestion | **[T, T, F]** | preserve foil boundary; LineTension/SimpleGlide short evolution |
| `hardening_simulation_cell` (`run_stem_hardening.py`) | `h_b` with `h[2,2] *= ZBOX` (default 5) | **[T, T, T]** | DDD_FFT / SUBCYCLING N-body elastic force needs a periodic box with non-degenerate z voxels |

**Honest caveat.** The thickened-periodic cell is used ONLY for the bulk-like *hardening pilot*; it
does **not** represent the true foil boundary. Any flow-stress / density-evolution number from that
cell is a bulk-like estimate, not a foil measurement. The `defect_ir` lowering exposes this as a
`cell_policy` argument (`as_is` vs `thickened_periodic`) so its effect can be swept by BO/UQ rather
than baked in.

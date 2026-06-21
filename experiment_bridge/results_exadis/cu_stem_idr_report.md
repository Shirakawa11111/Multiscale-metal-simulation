# STEM -> IDR audit report (Cu)

Source: `experiment_bridge/recon_data/points_3d*.txt`  |  valid=True  errors=0 warnings=0

| quantity | value |
|--|--|
| reconstructed lines | 27 |
| vertices / edges | 270 / 243 |
| Burgers source | geometry_only_pending_gb |
| endpoint policy | pinned_due_to_truncated_reconstruction |
| z uncertainty | ~30 nm (stereo-weak) |
| cell policy | foil_nonperiodic_z |

## Assignment uncertainty (the key audit result)
| metric | value |
|--|--|
| mean confidence | 0.3333 |
| min confidence | 0.3333 |
| frac low-confidence (<0.5) | 1.0 |
| mean assignment entropy (bits) | 1.5849 |
| ambiguous edges (entropy>0.8) | 243 / 243 |

**Highlight.** The IDR did not merely reformat the network — it **exposed** that the STEM->DDD
slip-system assignment is intrinsically ambiguous without g·b: geometry fixes the {111} *plane*, but
the three ⟨110⟩ Burgers on it are near-degenerate (mean confidence ~0.3333,
~1.5849 bits), so essentially every segment is assignment-ambiguous.
The legacy `stem_to_exadis.py` forced one of these and hid it. A real g·b analysis would set
`assignment_status = gb_validated` and collapse the priors.

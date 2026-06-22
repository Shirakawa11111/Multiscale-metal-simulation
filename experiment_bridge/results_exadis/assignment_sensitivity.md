# Assignment-policy sensitivity (v1.1, local, no DDD)

`cu_stem_idr.json`, 30 Monte-Carlo seeds for the sampling policies.

| policy | within-line discontinuities | of intra-line adjacencies | inventory CV |
|--|--|--|--|
| top1 | 0.0 | 216 | 0.0 |
| sample_edgewise (deprecated) | **142.8** | 216 | 0.206 |
| sample_linewise (default) | **0.0** | 216 | 0.622 |

**Verdict.** sample_edgewise injects ~143 artificial within-line Burgers discontinuities (=> artificial junctions); sample_linewise injects 0. Use sample_linewise. Edgewise sampling breaks Burgers continuity along a single
physical reconstructed line, manufacturing junction-like topology; line-coherent sampling preserves it.
The assignment *ambiguity itself* is real (geometry fixes the {111} plane, the 3 ⟨110⟩ Burgers are
near-degenerate), but it must be propagated **per line**, not per segment.

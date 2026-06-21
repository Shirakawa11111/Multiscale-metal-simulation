# Assignment / cell-policy sensitivity audit (local, no DDD)

`cu_stem_idr.json`, 50 Monte-Carlo assignment samples, 243 segments.

## Assignment policy = sample
| metric | value |
|--|--|
| mean per-segment assignment entropy (bits) | 1.5548 |
| fraction of segments that take >1 system across samples | 1.0 |
| distinct slip systems used | 12 |
| inventory coefficient-of-variation (mean) | 0.212 |

## Cell policy contrast
| policy | is_periodic | z box (b) |
|--|--|--|
| thickened_zbox3 | [True, True, True] | 1800 |
| thickened_zbox5 | [True, True, True] | 3000 |
| thickened_zbox10 | [True, True, True] | 6000 |
| as_is_foil | [True, True, False] | 600 |

**Interpretation.** Because every segment's geometric assignment is ~3-way degenerate (confidence ~0.33), sampled lowerings reshuffle the slip-system inventory substantially -- so any downstream DDD metric must be reported as a distribution over assignment samples, not a single number. Cell policy (foil non-periodic vs thickened periodic) changes z-periodicity and box height, a separate, deterministic knob.

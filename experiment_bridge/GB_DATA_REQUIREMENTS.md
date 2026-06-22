# g·b data requirements — turning a geometry-only IDR into a gb-validated one

What diffraction-contrast data the STEM-to-DDD v2 pipeline needs to *resolve* the slip-system assignment,
and exactly how it enters the IDR. This is the experimental-input spec for the next stage; the interface and
the entropy-collapse mechanism are already validated synthetically (`synthetic_gb.py`).

## Why g·b is needed
Geometry alone fixes the **glide plane** (the {111} with `argmin|n·t|`) but leaves the **Burgers vector**
3-way degenerate: the three ⟨110⟩ directions lying in that {111} are near-equally compatible with the line
tangent. Result on the real Cu network: every edge is ~3-way ambiguous — mean confidence 0.333, entropy
**1.585 bits** (= log₂3), 243/243 edges ambiguous. No amount of geometry removes this; it requires contrast.

## The criterion (what the code applies)
A dislocation is (to first order) **invisible** under reflection **g** when **g·b = 0**. Each observation is
a pair `(g, visible)`. `defect_ir.uncertainty.apply_gb_constraints` keeps a candidate only if its predicted
visibility matches every observation, using the normalized criterion |ĝ·b̂| > tol (tol = angle tolerance,
scale-invariant proxy for g·b = 0). It returns `gb_validated` (collapsed to one), `gb_partial` (narrowed),
or `gb_inconsistent` (no candidate fits — flags a contradiction).

## How many reflections per line?
From the line-coherent synthetic study (`synthetic_gb.md`), per reconstructed line:

| reflections | mean entropy | frac fully resolved | meaning |
|--|--|--|--|
| 0 (geometry only) | 1.585 bits | 0% | 3-way degenerate |
| 1 well-chosen | 0.704 bits | ~30% | one candidate ruled in/out; usually a 2-way residue |
| **2 well-chosen** | **0.0 bits** | **100%** | **single slip system (`gb_validated`)** |

**Recommendation: acquire 2 well-chosen reflections per line.** One reflection roughly halves the entropy
but rarely uniquely resolves the 3 ⟨110⟩; two non-collinear reflections (e.g. a `g=200` type and a `g=020`
type) partition the three candidates uniquely. A third reflection is useful only as a **consistency check**
(should return `gb_validated` with no new information) or where the first two are nearly collinear.

**"Well-chosen"** = reflections whose invisibility pattern *differs* across the three in-plane ⟨110⟩
candidates. A reflection with g·b ≠ 0 for all three carries no discriminating power for that line. In
practice pick g's that make *different* candidates invisible (g·b = 0 for one but not the others).

## How it enters the IDR
Per line (or per edge, line-coherent), populate `uncertainty.gb_constraints` with the observed reflections:
```json
{
  "uncertainty": {
    "gb_constraints": {
      "<parent_line_id>": {
        "observations": [
          {"g": [2, 0, 0], "visible": true},
          {"g": [0, 2, 0], "visible": false}
        ]
      }
    }
  }
}
```
The lowering pass calls `apply_gb_constraints(candidates, observations)` → narrows each line's candidate set
and updates `assignment_status` (`geometry_only_pending_gb` → `gb_partial` / `gb_validated`). Downstream
DDD lowering and every audit observable are **unchanged** — only the assignment confidence improves.

## Partial coverage (the realistic case)
Real campaigns rarely image every line under two reflections. The IDR handles a **mixed** network natively:
- lines with ≥2 discriminating reflections → `gb_validated` (entropy 0),
- lines with 1 reflection → `gb_partial` (entropy reduced, candidate set narrowed),
- lines with none → `geometry_only_pending_gb` (unchanged 3-candidate set).

Report the network as a **coverage histogram** (fraction validated / partial / pending) plus the
network-mean entropy before vs after. The natural validation experiment is then:
**geometry-only IDR → partial-g·b IDR → (mostly) gb-validated IDR**, comparing assignment entropy,
slip-system inventory, DDD network stability, junction count, and the density convention across the three.

## Acquisition checklist (for the microscopist)
1. For each reconstructed line, note its geometry-assigned {111} plane (from `cu_stem_idr_report`).
2. Choose 2 non-collinear reflections that discriminate the three in-plane ⟨110⟩ (different invisibility).
3. Record `visible` / `invisible` for each (g, line).
4. Enter as `(g, visible)` pairs in `uncertainty.gb_constraints` keyed by `parent_line_id`.
5. Re-run lowering; check `assignment_status` and the coverage histogram.

> Bottom line: **2 well-chosen reflections per line collapse the assignment to a single slip system.** One is
> partially informative; partial coverage is supported and yields a mixed validated/pending network. This is
> the highest-value *experimental* input for upgrading STEM-to-DDD v2 from auditable to physics-validated.

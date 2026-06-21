# defect_ir — unified Defect Intermediate Representation (IDR)

The structural backbone of the framework: **one** uncertainty-aware representation for both
- 3D crystalline dislocation networks (Cu, STEM-reconstructed → DDD/ExaDiS), and
- 2D atomic defect graphs (graphene / h-BN, 5|7 coordination cores → PFC).

Before IDR, the 2D and 3D defect work were two parallel pipelines. IDR makes them one framework and,
crucially, turns each route from "it runs" into "it is **auditable**" — every Burgers/slip-system
assignment carries **candidates + confidence**, not a single forced choice, so BO/UQ can act on the
representation itself.

## Layout
```
defect_ir/
  schema.py        schema v1: 6 sections (provenance, units, geometry, topology, uncertainty,
                   simulation_targets); enums; empty_idr(); edge_label()
  validators.py    validate_idr(doc, strict) -> (ok, errors, warnings); assert_valid()
  uncertainty.py   slip_system_candidates(tangent, catalog) -> top-k w/ softmax priors;
                   assignment_entropy(); network_assignment_summary()
  adapters/
    to_exadis.py   idr_to_exadis_network(doc, assignment_policy, cell_policy, zbox, seed)
                   -> ExaDiS manual-network dict (the stem_network.json shape)
  examples/
    build_examples.py        builds + validates both examples (also the smoke test)
    cu_stem_idr.json         3D Cu STEM dislocation graph (foil, z non-periodic)
    graphene_defect_idr.json 2D graphene 5|7 defect graph
```

## The six sections
- **provenance** — material, dimension, source_method, tool, input_ref, notes.
- **units** — length unit + `length_unit_m` (so 2D and 3D share a metric).
- **geometry** — `cell` (h / box_size_m / is_periodic) + `vertices` (`id`, `pos`, role, constraint) +
  `edges` (`id`, `v1`, `v2`, kind). Works for polylines (3D lines) and atom graphs (2D).
- **topology** — `slip_system_catalog`; per-edge `edge_labels` with **`slip_system_candidates`**
  (each with plane-containment / line-character scores + `prior`), `chosen_system`,
  `assignment_confidence`, `assignment_status`; per-vertex `vertex_labels` (coordination, defect_type).
- **uncertainty** — named, explicit: z-depth, Burgers source (geometry vs g·b), endpoint policy,
  system-size caveat, and an `assignment_summary` (mean confidence, entropy, # ambiguous edges).
- **simulation_targets** — engine + cell/loading/force policy for downstream legalization.

## Why this matters
The old `experiment_bridge/stem_to_exadis.py` forced one slip system per line via `argmin|n·t|`.
The IDR keeps the top-k candidates with priors; `idr_to_exadis_network(..., assignment_policy=...)`
can then take the top-1 or **sample** an assignment, and `network_assignment_summary` reports how
ambiguous the network is. That is the hook the BO/UQ layer needs.

## Run
```bash
python3 -m defect_ir.examples.build_examples   # build + validate examples, demo lowering
python3 tests/test_defect_ir.py                # gate test
```

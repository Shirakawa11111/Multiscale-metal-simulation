# Defect-IR spec — `defect_idr_v1` (normative)

A defect-IDR document is a JSON object with `format_version = "defect_idr_v1"` and six required
sections. Validated by `defect_ir.validators.validate_idr(doc, strict=False) -> (ok, errors, warnings)`.

## Top level (required)
`format_version, provenance, units, geometry, topology, uncertainty, simulation_targets`

## provenance  (required: material, dimension, source_method)
| field | type | notes |
|--|--|--|
| material | enum | `fcc_cu, graphene_2d, hbn_2d, bcc_fe, other` |
| dimension | 2 \| 3 | controls vertex `pos` length |
| source_method | enum | `stem_3d_reconstruction, pfc_defect_detection, md_defect_detection, synthetic, other` |
| tool, input_ref, created, notes | str | optional |

## units  (required: length, length_unit_m)
`length` (e.g. `"b"`, `"a0"`), `length_unit_m` (>0), `burgers_magnitude_b` (optional).

## geometry  (required: cell, vertices, edges)
- **cell**: `is_periodic` (3-list of bool, required); `h` (3×3) and/or `box_size_m` (3-list).
  Invariant (warn): `box_size_m[k] ≈ h[k][k] · length_unit_m`.
- **vertices**: each `{id (unique), pos ([x,y] if dim=2 else [x,y,z]), role?, constraint?}`.
  `role ∈ {endpoint, interior, core, atom, junction}`, `constraint ∈ {pinned, free, surface, unknown}`.
- **edges**: each `{id (unique), v1, v2, kind?}`; `v1,v2` must be known vertex ids.
  `kind ∈ {dislocation_segment, core_pair, bond, other}`. Warn on zero-length edges.

## topology
- **slip_system_catalog**: `[{system_id, b:[3], n:[3]}]` (crystallographic reference).
- **edge_labels**: per edge `{edge_id, slip_system_candidates, chosen_system, assignment_confidence,
  assignment_status}`. Each candidate: `{system_id, b:[3], n:[3], plane_containment_score,
  line_character_score, gb_visibility_score|null, prior}`. Invariants:
  - candidate `prior` values sum to ≈ 1 (warn otherwise);
  - `chosen_system` ∈ candidates (error otherwise); should be the max-prior unless a `policy` is declared (warn);
  - **physics (error):** for `dimension=3`, every candidate must satisfy `b·n = 0` (Burgers lies in glide plane).
  - `assignment_status ∈ {geometry_only_pending_gb, gb_validated, topology_inferred, assumed, unknown}`.
- **vertex_labels**: per vertex `{vertex_id, coordination, defect_type}` (2D: `5|7|none`, etc.).
- **counts**: free-form summary (`n_vertices, n_edges, n_lines, density_m2, ...`).

## uncertainty  (named, explicit)
Recommended keys: `z_depth {weakly_constrained, sigma_nm}`, `burgers_assignment {method, validated_by_gb}`,
`endpoint_policy {policy}`, `system_size_caveat`, `assignment_summary` (from
`uncertainty.network_assignment_summary`). Warn: a `stem_3d_reconstruction` source with empty uncertainty.

## simulation_targets  (required: engine)
`engine ∈ {exadis, pfc, lammps, none}`, plus engine-specific policy: `cell_policy, loading_mode, edir,
erate, force_model, mobility, ...` (consumed by the adapters; see `experiment_bridge/CELL_POLICY.md`).

## Lowering
`defect_ir.adapters.to_exadis.idr_to_exadis_network(doc, assignment_policy, cell_policy, zbox, seed)`
commits the uncertainty-aware IDR to one concrete ExaDiS network. `assignment_policy ∈ {top1, sample}`,
`cell_policy ∈ {as_is, thickened_periodic}` — both are the BO/UQ sweep knobs.

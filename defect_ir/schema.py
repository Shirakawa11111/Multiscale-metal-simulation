"""Defect Intermediate Representation (IDR) — schema v1.

A single, uncertainty-aware representation that carries BOTH
  - 3D crystalline dislocation networks (Cu, STEM-reconstructed -> DDD/ExaDiS), and
  - 2D atomic defect graphs (graphene / h-BN, 5|7 coordination cores -> PFC),
so that 2D and 3D defect work become one framework instead of two parallel pipelines.

Design rule: **expose** domain-specific assumptions and uncertainty, never hide them.
The novelty is not a solver — it is this intermediate layer, on which BO/UQ can act
(slip-system assignment policy, z-scaling, endpoint policy, cell policy, ... ).

A document is a plain dict (JSON-serializable) with six sections:

  provenance          where the defects came from, by what method, with what caveats
  units               length unit + physical scale (so 2D and 3D share a metric)
  geometry            cell + vertices + edges (the graph; works for lines AND atom graphs)
  topology            crystallographic / topological labels per edge & vertex,
                      with slip-system CANDIDATES + confidence (not a single forced choice)
  uncertainty         explicit, named uncertainty fields (z depth, Burgers source, ...)
  simulation_targets  engine + cell/loading/force policy for downstream legalization

See `validators.validate_idr` for the enforced contract and `examples/` for one worked
instance per domain. Adapters in `adapters/` lower an IDR doc to a concrete engine input.
"""

FORMAT_VERSION = "defect_idr_v1"

# ---- enumerations (kept small + explicit) ----
MATERIALS = {"fcc_cu", "graphene_2d", "hbn_2d", "bcc_fe", "other"}
SOURCE_METHODS = {"stem_3d_reconstruction", "pfc_defect_detection",
                  "md_defect_detection", "synthetic", "other"}
VERTEX_ROLES = {"endpoint", "interior", "core", "atom", "junction"}
CONSTRAINTS = {"pinned", "free", "surface", "unknown"}
EDGE_KINDS = {"dislocation_segment", "core_pair", "bond", "other"}
BURGERS_SOURCES = {"geometry", "gb_contrast", "topology", "assumed", "unknown"}
ASSIGNMENT_STATUS = {"geometry_only_pending_gb", "gb_validated",
                     "topology_inferred", "assumed", "unknown"}
ENGINES = {"exadis", "pfc", "lammps", "none"}

# ---- required keys per section (the enforced contract) ----
REQUIRED_TOP = ["format_version", "provenance", "units", "geometry",
                "topology", "uncertainty", "simulation_targets"]
REQUIRED_PROVENANCE = ["material", "dimension", "source_method"]
REQUIRED_UNITS = ["length", "length_unit_m"]
REQUIRED_GEOMETRY = ["cell", "vertices", "edges"]
REQUIRED_CELL = ["is_periodic"]                       # h or box_size also expected (checked softly)
REQUIRED_VERTEX = ["id", "pos"]                       # pos = [x, y] (2D) or [x, y, z] (3D)
REQUIRED_EDGE = ["id", "v1", "v2"]
REQUIRED_SIM = ["engine"]


def empty_idr(material, dimension, source_method, length_unit="b", length_unit_m=2.556e-10):
    """Return a minimal, valid-shaped IDR skeleton to fill in."""
    return {
        "format_version": FORMAT_VERSION,
        "provenance": {"material": material, "dimension": int(dimension),
                       "source_method": source_method, "input_ref": None,
                       "tool": None, "created": None, "notes": ""},
        "units": {"length": length_unit, "length_unit_m": float(length_unit_m),
                  "burgers_magnitude_b": 1.0},
        "geometry": {
            "cell": {"h": None, "box_size_m": None, "is_periodic": [True, True, dimension == 3]},
            "vertices": [],   # {id, pos:[..], role, constraint}
            "edges": [],      # {id, v1, v2, kind}
        },
        "topology": {
            "slip_system_catalog": [],   # [{system_id, b:[..], n:[..]}]
            "edge_labels": [],           # see edge_label() below
            "vertex_labels": [],         # {vertex_id, coordination, defect_type}
            "counts": {},
        },
        "uncertainty": {},               # named fields; see examples
        "simulation_targets": {"engine": "none"},
    }


def edge_label(edge_id, candidates, chosen_system=None,
               assignment_status="geometry_only_pending_gb"):
    """Build an uncertainty-aware edge label.

    candidates: list of dicts, each
        {system_id, b:[3], n:[3], plane_containment_score, line_character_score,
         gb_visibility_score (or None), prior}
    `prior` is the normalized confidence over candidates (see uncertainty.candidate_priors).
    chosen_system: the system_id selected by the current policy (top-1 by default).
    """
    if chosen_system is None and candidates:
        chosen_system = max(candidates, key=lambda c: c.get("prior", 0.0))["system_id"]
    conf = next((c.get("prior") for c in candidates if c["system_id"] == chosen_system), None)
    return {"edge_id": int(edge_id),
            "slip_system_candidates": candidates,
            "chosen_system": chosen_system,
            "assignment_confidence": conf,
            "assignment_status": assignment_status}

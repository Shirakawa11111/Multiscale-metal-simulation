"""Validate a Defect-IDR document against schema v1. No third-party deps."""
from . import schema as S


def _is_num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def validate_idr(doc, strict=False):
    """Return (ok: bool, errors: list[str], warnings: list[str]).

    strict=True promotes warnings (soft/recommended fields) to errors.
    """
    errs, warns = [], []

    def need(cond, msg):
        if not cond:
            errs.append(msg)

    if not isinstance(doc, dict):
        return False, ["IDR must be a dict"], []
    for k in S.REQUIRED_TOP:
        need(k in doc, f"missing top-level section: {k}")
    if errs:
        return False, errs, warns
    if doc.get("format_version") != S.FORMAT_VERSION:
        warns.append(f"format_version != {S.FORMAT_VERSION} (got {doc.get('format_version')})")

    # provenance
    p = doc["provenance"]
    for k in S.REQUIRED_PROVENANCE:
        need(k in p, f"provenance.{k} missing")
    if p.get("material") not in S.MATERIALS:
        warns.append(f"provenance.material '{p.get('material')}' not in known set")
    if p.get("dimension") not in (2, 3):
        errs.append(f"provenance.dimension must be 2 or 3 (got {p.get('dimension')})")
    if p.get("source_method") not in S.SOURCE_METHODS:
        warns.append(f"provenance.source_method '{p.get('source_method')}' not in known set")
    dim = p.get("dimension")

    # units
    u = doc["units"]
    for k in S.REQUIRED_UNITS:
        need(k in u, f"units.{k} missing")
    need(_is_num(u.get("length_unit_m")) and u.get("length_unit_m", 0) > 0,
         "units.length_unit_m must be a positive number")

    # geometry
    g = doc["geometry"]
    for k in S.REQUIRED_GEOMETRY:
        need(k in g, f"geometry.{k} missing")
    cell = g.get("cell", {})
    need(isinstance(cell.get("is_periodic"), list) and len(cell["is_periodic"]) == 3,
         "geometry.cell.is_periodic must be a 3-list of bools")
    if cell.get("h") is None and cell.get("box_size_m") is None:
        warns.append("geometry.cell has neither 'h' nor 'box_size_m' (size unknown)")

    verts = g.get("vertices", [])
    edges = g.get("edges", [])
    need(isinstance(verts, list), "geometry.vertices must be a list")
    need(isinstance(edges, list), "geometry.edges must be a list")
    ids = set()
    for i, v in enumerate(verts):
        for k in S.REQUIRED_VERTEX:
            need(k in v, f"vertex[{i}].{k} missing")
        if "id" in v:
            ids.add(v["id"])
        pos = v.get("pos")
        if not (isinstance(pos, list) and len(pos) == dim and all(_is_num(c) for c in pos)):
            errs.append(f"vertex[{i}].pos must be a {dim}-list of numbers")
        if v.get("role") and v["role"] not in S.VERTEX_ROLES:
            warns.append(f"vertex[{i}].role '{v['role']}' unknown")
        if v.get("constraint") and v["constraint"] not in S.CONSTRAINTS:
            warns.append(f"vertex[{i}].constraint '{v['constraint']}' unknown")
    for i, e in enumerate(edges):
        for k in S.REQUIRED_EDGE:
            need(k in e, f"edge[{i}].{k} missing")
        for end in ("v1", "v2"):
            if e.get(end) not in ids:
                errs.append(f"edge[{i}].{end}={e.get(end)} not a known vertex id")
        if e.get("kind") and e["kind"] not in S.EDGE_KINDS:
            warns.append(f"edge[{i}].kind '{e['kind']}' unknown")

    # topology: edge labels & candidate confidence
    topo = doc["topology"]
    edge_ids = {e.get("id") for e in edges}
    for i, lab in enumerate(topo.get("edge_labels", [])):
        if lab.get("edge_id") not in edge_ids:
            errs.append(f"topology.edge_labels[{i}].edge_id not a known edge")
        cands = lab.get("slip_system_candidates", [])
        if cands:
            psum = sum(c.get("prior", 0.0) for c in cands)
            if abs(psum - 1.0) > 0.02:
                warns.append(f"edge_labels[{i}] candidate priors sum to {psum:.3f} (expect ~1.0)")
            if lab.get("chosen_system") not in {c.get("system_id") for c in cands}:
                errs.append(f"edge_labels[{i}].chosen_system not among candidates")
        if lab.get("assignment_status") and lab["assignment_status"] not in S.ASSIGNMENT_STATUS:
            warns.append(f"edge_labels[{i}].assignment_status '{lab['assignment_status']}' unknown")

    # simulation_targets
    st = doc["simulation_targets"]
    for k in S.REQUIRED_SIM:
        need(k in st, f"simulation_targets.{k} missing")
    if st.get("engine") not in S.ENGINES:
        warns.append(f"simulation_targets.engine '{st.get('engine')}' unknown")

    # uncertainty must be present (may be sparse) but should be non-trivial for experimental sources
    if p.get("source_method") == "stem_3d_reconstruction" and not doc.get("uncertainty"):
        warns.append("STEM source but empty uncertainty section (z-depth / Burgers source should be stated)")

    if strict:
        errs += [f"[strict] {w}" for w in warns]
        warns = []
    return (len(errs) == 0), errs, warns


def assert_valid(doc, strict=False):
    ok, errs, warns = validate_idr(doc, strict=strict)
    if not ok:
        raise ValueError("Invalid IDR:\n  - " + "\n  - ".join(errs))
    return warns

"""Uncertainty helpers: turn geometric slip-system scores into a CANDIDATE SET + confidence.

The old `stem_to_exadis.py` forced a single slip system per line via argmin|n.t|. That throws
away the fact that two or three systems are often nearly equally compatible with a reconstructed
line, and that the Burgers vector is geometry-only (no g.b contrast). Here we keep the top-k
candidates with normalized priors so downstream BO/UQ can (a) sample assignments, (b) report an
assignment-sensitivity, instead of trusting one hard choice.
"""

import math


def _norm(v):
    n = math.sqrt(sum(c * c for c in v))
    return [c / n for c in v] if n > 0 else list(v)


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def slip_system_candidates(
    tangent, slip_systems, T=0.15, topk=3, w_plane=1.0, w_char=0.0
):
    """Score every slip system against a line tangent and return top-k candidates w/ priors.

    tangent: 3-vector (line direction, need not be unit).
    slip_systems: list of {system_id, b:[3], n:[3]}.
    Scoring:
      plane_containment_score = |n . t_hat|        (0 = line lies in the glide plane = best)
      line_character_score    = |b . t_hat|        (1 = screw, 0 = edge; diagnostic, not used by default)
      quality = w_plane*(1-containment) + w_char*line_character
      prior   = softmax(quality / T)   over the kept candidates
    Returns candidates sorted by prior (desc), each:
      {system_id, b, n, plane_containment_score, line_character_score, gb_visibility_score=None, prior}
    """
    t = _norm(tangent)
    scored = []
    for s in slip_systems:
        n = _norm(s["n"])
        b = _norm(s["b"])
        containment = abs(_dot(n, t))  # lower is better
        character = abs(_dot(b, t))  # 1=screw, 0=edge
        quality = w_plane * (1.0 - containment) + w_char * character
        scored.append(
            {
                "system_id": s["system_id"],
                "b": list(s["b"]),
                "n": list(s["n"]),
                "plane_containment_score": round(containment, 4),
                "line_character_score": round(character, 4),
                "gb_visibility_score": None,
                "_q": quality,
            }
        )
    scored.sort(key=lambda c: c["_q"], reverse=True)
    keep = scored[:topk]
    mx = max(c["_q"] for c in keep)
    exps = [math.exp((c["_q"] - mx) / max(T, 1e-6)) for c in keep]
    Z = sum(exps)
    for c, e in zip(keep, exps):
        c["prior"] = round(e / Z, 4)
        del c["_q"]
    return keep


def assignment_entropy(candidates):
    """Shannon entropy (bits) of the candidate priors — a per-edge ambiguity score.
    0 = unambiguous (one system), higher = more ambiguous (assignment-sensitive)."""
    ps = [c.get("prior", 0.0) for c in candidates if c.get("prior", 0.0) > 0]
    return round(-sum(p * math.log2(p) for p in ps), 4) if ps else 0.0


def network_assignment_summary(edge_labels):
    """Aggregate confidence/entropy across all edges -> an auditable assignment-uncertainty report."""
    confs = [
        l.get("assignment_confidence")
        for l in edge_labels
        if l.get("assignment_confidence") is not None
    ]
    ents = [
        assignment_entropy(l.get("slip_system_candidates", [])) for l in edge_labels
    ]
    n = len(edge_labels) or 1
    return {
        "n_edges": len(edge_labels),
        "mean_confidence": round(sum(confs) / len(confs), 4) if confs else None,
        "min_confidence": round(min(confs), 4) if confs else None,
        "frac_low_confidence(<0.5)": round(sum(c < 0.5 for c in confs) / n, 4),
        "mean_assignment_entropy_bits": round(sum(ents) / n, 4),
        "n_ambiguous(entropy>0.8)": sum(e > 0.8 for e in ents),
    }


# ---- g.b (diffraction-contrast) constraint: collapse assignment ambiguity with real/synthetic data ----
def apply_gb_constraints(candidates, observations, tol=0.1):
    """Filter + renormalize slip-system candidates using g.b invisibility observations.

    A dislocation with Burgers b is (to first order) INVISIBLE under reflection g when g.b = 0. Each
    observation is {g:[3], visible:bool}; a candidate is kept only if its predicted visibility
    (|g.b_hat| > tol) matches every observation. Surviving candidates are renormalized; each gets a
    `gb_compatible=True` and `gb_visibility_score` = min over observations of |g.b_hat|.
    Returns (kept_candidates, assignment_status) -- status 'gb_validated' if it collapses to one, else
    'gb_partial'. With no observations, returns the input unchanged ('geometry_only_pending_gb').
    """
    if not observations:
        return list(candidates), "geometry_only_pending_gb"
    kept = []
    for c in candidates:
        bh = _norm(c["b"])
        scores = []
        ok = True
        for o in observations:
            gb = abs(_dot(_norm(o["g"]), bh))
            predicted_visible = gb > tol
            scores.append(gb)
            if predicted_visible != bool(o["visible"]):
                ok = False
                break
        if ok:
            cc = dict(c)
            cc["gb_compatible"] = True
            cc["gb_visibility_score"] = round(min(scores), 4)
            kept.append(cc)
    if not kept:  # observations inconsistent with all -> keep all, flag
        return list(candidates), "gb_inconsistent"
    Z = sum(c.get("prior", 0.0) for c in kept) or 1.0
    for c in kept:
        c["prior"] = round(c.get("prior", 0.0) / Z, 4)
    status = "gb_validated" if len(kept) == 1 else "gb_partial"
    return kept, status

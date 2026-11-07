"""
Volume selection shared by bki-plot, bki-distrib and notebook_plot.

Spec grammar:
  pick / surface   NAME_RE#SEL
  highlight        [LABEL=]NAME_RE[#SEL]
where SEL is an index list (``0,2,5-7``, radio ``meta.index``) or ``@REGEX``
matched against the ``_PICK_KEYS`` metadata fields.
"""
from __future__ import annotations

import re

from .ir import IRCollection, VolumeEntry

# metadata fields an @REGEX selector is matched against (those present)
_PICK_KEYS = ("surface", "node", "nuclide", "path", "part", "pv", "pvpath")

# LABEL= prefix of a highlight spec: no regex characters before the '='
_LABEL_RE = re.compile(r"^([^=#\\()\[\]{}^$|*+?.]+)=(.+)$")


# ---------------------------------------------------------------------------
# Spec parsing
# ---------------------------------------------------------------------------

def _parse_index_list(sel: str) -> set[int]:
    """'0,2,5-7' → {0, 2, 5, 6, 7}."""
    out: set[int] = set()
    for tok in sel.split(","):
        tok = tok.strip()
        if not tok:
            continue
        a, sep, b = tok.partition("-")
        if sep:
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(a))
    return out


def parse_pick_spec(spec: str) -> tuple:
    """
    Parse a ``NAME_RE#SEL`` selector (``--pick`` / ``--surface``).

    ``SEL`` is either an index list (``0,2,5-7``, radio ``meta.index``) or
    ``@REGEX`` matched against ``meta.surface``, ``meta.node``,
    ``meta.nuclide`` (radio) and ``meta.path``, ``meta.part``, ``meta.pv``,
    ``meta.pvpath`` (geometry).

    Returns (name_re, indices | None, match_re | None).
    """
    name, sep, sel = spec.rpartition("#")
    if not sep or not name or not sel:
        raise ValueError(f"bad spec '{spec}' (expected NAME_RE#SEL)")
    if sel.startswith("@"):
        return re.compile(name), None, re.compile(sel[1:])
    try:
        return re.compile(name), _parse_index_list(sel), None
    except ValueError:
        raise ValueError(f"bad index list '{sel}' in spec '{spec}'") from None


def parse_highlight_spec(spec: str) -> tuple:
    """
    Parse a ``--highlight`` spec ``[LABEL=]NAME_RE[#SEL]``.

    ``SEL`` is as in :func:`parse_pick_spec`; without it every placement of
    the matching volumes is highlighted.  ``LABEL`` (no regex characters)
    names the group in the legend; it defaults to the spec itself.

    Returns (label, name_re, indices | None, match_re | None).
    """
    m = _LABEL_RE.match(spec)
    label, body = ((m.group(1).strip(), m.group(2).lstrip()) if m
                   else (spec, spec))
    if "#" in body:
        return (label,) + parse_pick_spec(body)
    if not body:
        raise ValueError(f"bad highlight spec '{spec}' (expected [LABEL=]NAME_RE[#SEL])")
    return label, re.compile(body), None, None


def validate_picks(surface: list[str] | None = None,
                   pick: list[str] | None = None,
                   highlight: list[str] | None = None) -> None:
    """Parse every ``--surface``, ``--pick`` and ``--highlight`` spec, in
    that order.  Raises ValueError("--OPT: reason") for the first bad one
    (a bad regex included), for the CLI to report before any selection."""
    for opt, specs, parse in (("--surface", surface, parse_pick_spec),
                              ("--pick", pick, parse_pick_spec),
                              ("--highlight", highlight, parse_highlight_spec)):
        for spec in specs or ():
            try:
                parse(spec)
            except (ValueError, re.error) as exc:
                raise ValueError(f"{opt}: {exc}") from None


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def _spec_hits(v: VolumeEntry, indices: set[int] | None, match,
               verb: str = "picked", need_path: bool = True) -> bool:
    """True if the SEL part of one spec (`indices` or `match`) selects `v`;
    the spec's NAME_RE is assumed to match already.

    An index list on a geometry volume raises ValueError ("geometry volumes
    are <verb> by path").  With `need_path`, so does an @REGEX on a geometry
    volume without ``meta.path`` (an IR from an older bki-extract)."""
    m = v.meta
    if indices is not None:
        if v.kind != "radio":
            raise ValueError(
                f"index selector matches {v.kind} volume '{v.name}'; "
                f"geometry volumes are {verb} by path (NAME_RE#@REGEX)")
        if m.get("index") in indices:
            return True
    if match is not None:
        if need_path and v.kind != "radio" and "path" not in m:
            raise ValueError(
                f"{v.kind} volume '{v.name}' has no placement path; "
                f"re-run bki-extract to pick geometry volumes")
        if any(match.search(str(m[k])) for k in _PICK_KEYS
               if m.get(k) is not None):
            return True
    return False


def _picked(v: VolumeEntry, specs: list[tuple]) -> bool:
    """Pick/surface filter: True if some spec naming `v` selects it.
    Volumes not named by any spec are kept."""
    applicable = [s for s in specs if s[0].search(v.name)]
    if not applicable:
        return True
    # evaluate every spec (no short-circuit): an index selector on a geometry
    # volume must raise whatever the spec order
    hits = [_spec_hits(v, indices, match) for _, indices, match in applicable]
    return any(hits)


def match_highlight(v: VolumeEntry, specs: list[tuple]) -> int | None:
    """Index of the first highlight spec matching `v`, or None.

    Unlike the pick filter, a volume not named by any spec is *not* matched."""
    for k, (_, name_re, indices, match) in enumerate(specs):
        if not name_re.search(v.name):
            continue
        if indices is None and match is None:
            return k
        if _spec_hits(v, indices, match, verb="highlighted", need_path=False):
            return k
    return None


def highlight_groups(vols: list[VolumeEntry], highlight: list[str]) -> tuple:
    """Split `vols` by highlight specs.

    Returns (labels, groups, context): one volume list per spec, in spec
    order, and the volumes no spec matches."""
    specs  = [parse_highlight_spec(h) for h in highlight]
    groups: list[list[VolumeEntry]] = [[] for _ in specs]
    context: list[VolumeEntry] = []
    for v in vols:
        k = match_highlight(v, specs)
        (context if k is None else groups[k]).append(v)
    return [sp[0] for sp in specs], groups, context


def _compile(pat: str) -> re.Pattern:
    try:
        return re.compile(pat)
    except re.error as exc:
        raise ValueError(f"bad regex '{pat}': {exc}") from None


def select_volumes(
    ir: IRCollection,
    kinds: list[str] | None = None,
    name_patterns: list[str] | None = None,
    kind_patterns: dict[str, list[str]] | None = None,
    surfaces: list[str] | None = None,
    picks: list[str] | None = None,
) -> list[VolumeEntry]:
    """
    Apply the bki-plot selectors to ``ir.volumes`` (see
    :func:`duneradgeo.plot.plot_ir` for their meaning).
    Raises ValueError when a pick spec cannot be applied to a volume or a
    name/kind pattern is not a valid regex.
    """
    vols = ir.volumes
    if kinds is not None:
        kind_set = set(kinds)
        vols = [v for v in vols if v.kind in kind_set]
    if name_patterns:
        pats = [_compile(p) for p in name_patterns]
        vols = [v for v in vols if any(p.search(v.name) for p in pats)]
    if kind_patterns:
        compiled = {k: [_compile(p) for p in ps]
                    for k, ps in kind_patterns.items()}
        vols = [v for v in vols
                if v.kind not in compiled
                or any(p.search(v.name) for p in compiled[v.kind])]
    if surfaces:
        specs = [parse_pick_spec(s) for s in surfaces]
        vols = [v for v in vols if v.kind != "radio" or _picked(v, specs)]
    if picks:
        specs = [parse_pick_spec(s) for s in picks]
        vols = [v for v in vols if _picked(v, specs)]
    return vols

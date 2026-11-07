"""
bki-plot  —  render a JSON IR as a matplotlib PDF

Produces the same 3D + 3×2D layout as the original plot_gdml_detectors.py:
                 --orientation hd (default)   --orientation vd (x vertical)
  top-left       3D ortho view                3D ortho view
  top-right      Z–Y  beam face               Z–X  beam face
  bot-left       Z–X  top view                Z–Y  top view
  bot-right      X–Y  side view               Y–X  side view

Usage:
  bki-plot detector.json --out detector.pdf
  bki-plot detector.json --out detector.pdf --view 20 -60
  bki-plot detector.json --highlight 'bottom planes=^volTPCPlane#@volTPC_PV_Bot'
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from ._draw import BOX_FACES, box_corners_world, project_polygon
from .ir import IRCollection, VolumeEntry


# ---------------------------------------------------------------------------
# Colour/style constants
# ---------------------------------------------------------------------------

_ANODE_COL   = "royalblue"
_OPDET_COL   = "darkorange"
_CATHODE_COL = "teal"
_FSHAPER_COL = "slategray"
_RADIO_COL   = "crimson"
_FACE_ALPHA  = 0.25
_EDGE_ALPHA  = 0.6

_STRUCT_PALETTE: dict[str, tuple] = {
    "volTPCActive":      ("#4CAF50", 0.08),   # green
    "volTPCActiveInner": ("#4CAF50", 0.08),
    "volGaseousArgon":   ("#90CAF9", 0.18),   # light blue
    "volTPC":            ("#FFF176", 0.06),   # pale yellow
    "volCryostat":       ("#CE93D8", 0.05),   # pale purple
}
_STRUCT_DEFAULT = ("#BDBDBD", 0.08)

# Per-surface colours used when individual radio surfaces are selected;
# above this many groups the colouring collapses to one colour per producer.
_SURFACE_PALETTE = ("tab:red", "tab:blue", "tab:green", "tab:purple",
                    "tab:brown", "tab:pink", "tab:olive", "tab:cyan",
                    "tab:orange", "tab:gray")
_MAX_SURFACE_GROUPS = 20

# --highlight: highlighted groups take _SURFACE_PALETTE colours, everything
# else is drawn as grey context with face opacity `dim`.
_CONTEXT_COL   = "#9E9E9E"
# low: faces of overlapping context boxes stack up along the view direction
_DIM_DEFAULT   = 0.02
_HL_FACE_ALPHA = 0.45
_HL_LINEWIDTH  = 1.2
# a highlighted box whose 2D projection is thinner than this fraction of the
# axis span is drawn as a line, so e.g. 0.02 cm wire planes stay visible
_THIN_FRAC     = 0.004

# Axis layouts, keyed by orientation.  World axes: 0 = x (drift),
# 1 = y, 2 = z (beam).  "order3d" maps world axes onto the 3D plot's
# (X, Y, Z=up) axes; "views" lists the (horizontal, vertical, title) world
# axes of the top-right, bottom-left and bottom-right 2D panels.
#   hd: drift horizontal, y vertical.
#   vd: drift vertical (x up); each panel keeps its role, with x and y
#       swapped.  The 3D Y axis (world y) is inverted to stay right-handed.
_LAYOUTS: dict[str, dict] = {
    "hd": dict(order3d=(2, 0, 1), flip3d_y=False,
               labels=("x (drift)", "y (vert)", "z (beam)"),
               views=((2, 1, "Z–Y  (beam face)"),
                      (2, 0, "Z–X  (top view)"),
                      (0, 1, "X–Y  (side view)"))),
    "vd": dict(order3d=(2, 1, 0), flip3d_y=True,
               labels=("x (drift, vert)", "y", "z (beam)"),
               views=((2, 0, "Z–X  (beam face)"),
                      (2, 1, "Z–Y  (top view)"),
                      (1, 0, "Y–X  (side view)"))),
}


# ---------------------------------------------------------------------------
# Volume selection
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
    ``meta.nuclide`` (radio) and ``meta.path``, ``meta.part`` (geometry).

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


_PICK_KEYS = ("surface", "node", "nuclide", "path", "part", "pv", "pvpath")


def _picked(v: VolumeEntry, specs: list[tuple]) -> bool:
    applicable = [s for s in specs if s[0].search(v.name)]
    if not applicable:
        return True          # volumes not named by any spec are untouched
    m = v.meta
    for _, indices, match in applicable:
        if indices is not None:
            if v.kind != "radio":
                raise ValueError(
                    f"index selector matches {v.kind} volume '{v.name}'; "
                    f"geometry volumes are picked by path (NAME_RE#@REGEX)")
            if m.get("index") in indices:
                return True
        if match is not None:
            if v.kind != "radio" and "path" not in m:
                raise ValueError(
                    f"{v.kind} volume '{v.name}' has no placement path; "
                    f"re-run bki-extract to pick geometry volumes")
            if any(match.search(str(m[k])) for k in _PICK_KEYS
                   if m.get(k) is not None):
                return True
    return False


_LABEL_RE = re.compile(r"^([^=#\\()\[\]{}^$|*+?.]+)=(.+)$")


def parse_highlight_spec(spec: str) -> tuple:
    """
    Parse a ``--highlight`` spec ``[LABEL=]NAME_RE[#SEL]``.

    ``SEL`` is as in :func:`parse_pick_spec`; without it every placement of
    the matching volumes is highlighted.  ``LABEL`` (no regex characters)
    names the group in the legend; it defaults to the spec itself.

    Returns (label, name_re, indices | None, match_re | None).
    """
    m = _LABEL_RE.match(spec)
    label, body = (m.group(1).strip(), m.group(2)) if m else (spec, spec)
    if "#" in body:
        return (label,) + parse_pick_spec(body)
    if not body:
        raise ValueError(f"bad highlight spec '{spec}' (expected [LABEL=]NAME_RE[#SEL])")
    return label, re.compile(body), None, None


def match_highlight(v: VolumeEntry, specs: list[tuple]) -> int | None:
    """Index of the first highlight spec matching `v`, or None.

    Unlike :func:`_picked`, a volume not named by any spec is *not* matched."""
    m = v.meta
    for k, (_, name_re, indices, match) in enumerate(specs):
        if not name_re.search(v.name):
            continue
        if indices is None and match is None:
            return k
        if indices is not None:
            if v.kind != "radio":
                raise ValueError(
                    f"index selector matches {v.kind} volume '{v.name}'; "
                    f"geometry volumes are highlighted by path (NAME_RE#@REGEX)")
            if m.get("index") in indices:
                return k
        if match is not None and any(match.search(str(m[key])) for key in _PICK_KEYS
                                     if m.get(key) is not None):
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


def select_volumes(
    ir: IRCollection,
    kinds: list[str] | None = None,
    name_patterns: list[str] | None = None,
    kind_patterns: dict[str, list[str]] | None = None,
    surfaces: list[str] | None = None,
    picks: list[str] | None = None,
) -> list[VolumeEntry]:
    """
    Apply the bki-plot selectors to ``ir.volumes`` (see :func:`plot_ir`).
    Raises ValueError when a pick spec cannot be applied to a volume.
    """

    vols = ir.volumes
    if kinds is not None:
        kind_set = set(kinds)
        vols = [v for v in vols if v.kind in kind_set]
    if name_patterns:
        pats = [re.compile(p) for p in name_patterns]
        vols = [v for v in vols if any(p.search(v.name) for p in pats)]
    if kind_patterns:
        compiled = {k: [re.compile(p) for p in ps]
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


def _world_bounds(vols: list[VolumeEntry], pad_frac: float = 0.03) -> tuple:
    """
    Padded world-frame (lo, hi) [x, y, z] enclosing all corners of `vols`,
    or (None, None) when empty.  The 2D views draw projections of these same
    corners, so the box also encloses every drawn 2D polygon.
    """
    if not vols:
        return None, None
    arr = np.vstack([box_corners_world(v.pos, v.ext, v.rot) for v in vols])
    lo, hi = arr.min(axis=0), arr.max(axis=0)
    pad = pad_frac * np.maximum(hi - lo, 1.0)
    return lo - pad, hi + pad


def _surface_tag(v: VolumeEntry) -> str:
    """'producer#idx' for radio entries, 'volume#copy' for parts of
    multi-part solids, plain name otherwise."""
    idx = v.meta.get("index")
    if v.kind == "radio" and idx is not None:
        return f"{v.name}#{idx}"
    copy = v.meta.get("copy")
    return f"{v.name}#{copy}" if copy is not None else v.name


def _n_placements(vols: list[VolumeEntry]) -> int:
    """Placements in `vols`, counting the parts of a multi-part solid once."""
    return len({(v.kind, v.name, v.meta.get("copy", id(v))) for v in vols})


# ---------------------------------------------------------------------------
# Core plotting function (importable from notebooks etc.)
# ---------------------------------------------------------------------------

def plot_ir(
    ir: IRCollection,
    out_path: str,
    title: str = "",
    view: tuple = (20.0, -60.0),
    kinds: list[str] | None = None,
    name_patterns: list[str] | None = None,
    kind_patterns: dict[str, list[str]] | None = None,
    surfaces: list[str] | None = None,
    picks: list[str] | None = None,
    proportional: bool = True,
    orientation: str = "hd",
    highlight: list[str] | None = None,
    dim: float = _DIM_DEFAULT,
    zoom_highlight: bool = False,
) -> None:
    """
    Render an IRCollection to a PDF/PNG file.

    Parameters
    ----------
    kinds : list of str, optional
        Only plot volumes whose ``kind`` is in this list.  Default: all.
    name_patterns : list of str, optional
        Only plot volumes whose ``name`` matches at least one regex
        (re.search).  Applied across all kinds after the ``kinds`` filter.
    kind_patterns : dict[str, list[str]], optional
        Per-kind name filters, e.g. ``{"radio": ["GenInCathode"],
        "struct": ["volGaseousArgon"]}``.  For each kind present in this
        dict only volumes whose name matches at least one pattern are kept;
        kinds absent from the dict are left unfiltered.
    surfaces : list of str, optional
        Radio generation-surface selectors ``PRODUCER_RE#SEL``, where SEL is
        an index list (``0,2,5-7``, matched against ``meta.index``) or
        ``@REGEX`` (matched against ``meta.surface``, ``meta.node`` and
        ``meta.nuclide``).  Only producers matched by some spec are
        narrowed.  When given, each selected surface gets its own colour.
    picks : list of str, optional
        Like ``surfaces`` but for every kind: ``NAME_RE#@REGEX`` keeps the
        placements of matching volumes whose ``meta.path`` (e.g.
        ``.../volTPC:37/volTPCPlaneV:0``) or ``meta.part`` matches; index
        lists apply to radio entries only.  Volumes not matched by any
        ``NAME_RE`` are untouched.
    proportional : bool, optional
        Draw the 2D views to scale (1 cm == 1 cm on both axes; default).
        When False each 2D panel is stretched to fill its grid cell.
    orientation : {"hd", "vd", "auto"}, optional
        ``"hd"`` (default) draws the drift axis x horizontally with y
        vertical; ``"vd"`` draws x vertically, in every 2D view that
        contains it and in the 3D view.  ``"auto"`` picks ``"vd"`` when
        ``ir.detector == "vd"``, else ``"hd"``.
    highlight : list of str, optional
        ``[LABEL=]NAME_RE[#SEL]`` specs (SEL as for ``picks``, optional).
        Applied after the selectors: volumes matched by a spec are drawn on
        top in that spec's colour, with one legend entry per spec; all other
        selected volumes are drawn as grey context.  Highlighted boxes too
        thin to see in a 2D view are drawn as lines.
    dim : float, optional
        Face opacity of the context volumes when ``highlight`` is given.
    zoom_highlight : bool, optional
        Fit the axes to the highlighted volumes instead of all drawn ones.

    Axes convention: x = drift, y = vertical (HD) / horizontal (VD),
    z = beam.  In HD mode the 3D axes are re-ordered as (z, x, y) to match
    the original script; in VD mode as (z, y, x).
    """
    if orientation == "auto":
        orientation = "vd" if ir.detector == "vd" else "hd"
    if orientation not in _LAYOUTS:
        raise ValueError(f"unknown orientation '{orientation}' "
                         f"(expected one of: hd, vd, auto)")
    layout = _LAYOUTS[orientation]
    o3     = layout["order3d"]

    try:
        import matplotlib
    except ImportError as exc:
        raise ImportError("plotting needs matplotlib: "
                          "pip install matplotlib") from exc
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    vols = select_volumes(ir, kinds=kinds, name_patterns=name_patterns,
                          kind_patterns=kind_patterns, surfaces=surfaces,
                          picks=picks)
    all_vols = vols

    hl_labels: list[str] = []
    hl_groups: list[list[VolumeEntry]] = []
    if highlight:
        hl_labels, hl_groups, vols = highlight_groups(vols, highlight)
        for lbl, g in zip(hl_labels, hl_groups):
            if not g:
                print(f"[warn] --highlight '{lbl}' matches no selected volume",
                      file=sys.stderr)

    filtered = IRCollection(source=ir.source, detector=ir.detector, volumes=vols,
                            frame=ir.frame)

    anodes  = filtered.by_kind("anode")
    opdets  = filtered.by_kind("opdet")
    cathodes = filtered.by_kind("cathode")
    fshapers = filtered.by_kind("fieldshaper")
    structs = filtered.by_kind("struct")
    radios  = filtered.by_kind("radio")

    fig = plt.figure(figsize=(16, 10))
    if title:
        fig.suptitle(title, fontsize=11, fontweight="bold")

    ax3 = fig.add_subplot(2, 2, 1, projection="3d")
    # (axes, horizontal world axis, vertical world axis, title) per 2D panel
    views2d = [(fig.add_subplot(2, 2, k), ah, av, ttl)
               for k, (ah, av, ttl) in zip((2, 3, 4), layout["views"])]

    # ---- axis limits (enclose every drawn volume, or the highlights) ----
    hl_vols = [v for g in hl_groups for v in g]
    lo, hi = _world_bounds(hl_vols if zoom_highlight and hl_vols else all_vols)

    # ---- batched box drawing (context and highlights) ----
    from matplotlib.collections import LineCollection, PolyCollection
    from matplotlib.colors import to_rgba

    def _draw_batch(bvols: list[VolumeEntry], col: str, label: str,
                    face_alpha: float, edge_alpha: float, lw: float,
                    zorder: int, thin_as_lines: bool = False):
        if not bvols:
            return
        face, edge = to_rgba(col, face_alpha), to_rgba(col, edge_alpha)
        corners = [box_corners_world(v.pos, v.ext, v.rot) for v in bvols]
        faces3d = [[[c[j][o3[0]], c[j][o3[1]], c[j][o3[2]]] for j in fi]
                   for c in corners for fi in BOX_FACES]
        pc3 = Poly3DCollection(faces3d, facecolor=face, edgecolor=edge, linewidth=lw)
        pc3.set_zorder(zorder)
        ax3.add_collection3d(pc3)
        ax3.plot([], [], [], color=col, label=label, linewidth=3)
        for ax2d, ah, av, _ in views2d:
            polys, lines = [], []
            for c in corners:
                poly = project_polygon(c, ah, av)
                if thin_as_lines and lo is not None:
                    pmin, pmax = poly.min(axis=0), poly.max(axis=0)
                    w, h = pmax - pmin
                    mid = (pmin + pmax) / 2
                    if w < _THIN_FRAC * (hi[ah] - lo[ah]):
                        lines.append([(mid[0], pmin[1]), (mid[0], pmax[1])])
                        continue
                    if h < _THIN_FRAC * (hi[av] - lo[av]):
                        lines.append([(pmin[0], mid[1]), (pmax[0], mid[1])])
                        continue
                polys.append(poly)
            if polys:
                ax2d.add_collection(PolyCollection(polys, facecolor=face, edgecolor=edge,
                                                   linewidth=lw, zorder=zorder))
            if lines:
                ax2d.add_collection(LineCollection(lines, colors=[to_rgba(col, 1.0)],
                                                   linewidth=1.5, zorder=zorder))
            ax2d.plot([], [], color=col, label=label, linewidth=3)

    if highlight:
        _draw_batch(vols, _CONTEXT_COL, f"other ({_n_placements(vols)})",
                    dim, min(1.0, 5 * dim), 0.2, zorder=1)
        for k, (lbl, g) in enumerate(zip(hl_labels, hl_groups)):
            _draw_batch(g, _SURFACE_PALETTE[k % len(_SURFACE_PALETTE)],
                        f"{lbl} ({_n_placements(g)})",
                        _HL_FACE_ALPHA, 1.0, _HL_LINEWIDTH, zorder=10 + k,
                        thin_as_lines=True)
        anodes = opdets = cathodes = fshapers = structs = radios = []

    # ---- structural volumes (batched per name) ----
    s_face_groups: dict[str, list]       = defaultdict(list)
    s_poly_groups: dict[str, dict]       = {}

    for v in structs:
        corners = box_corners_world(v.pos, v.ext, v.rot)
        for fi in BOX_FACES:
            s_face_groups[v.name].append(
                [[corners[j][o3[0]], corners[j][o3[1]], corners[j][o3[2]]] for j in fi])
        if v.name not in s_poly_groups:
            s_poly_groups[v.name] = {ax2d: [] for ax2d, *_ in views2d}
        for ax2d, ah, av, _ in views2d:
            s_poly_groups[v.name][ax2d].append(project_polygon(corners, ah, av))

    for vname, faces in s_face_groups.items():
        col, alpha = _STRUCT_PALETTE.get(vname, _STRUCT_DEFAULT)
        ax3.add_collection3d(Poly3DCollection(
            faces, alpha=alpha, facecolor=col, edgecolor=col, linewidth=0.3))
        ax3.plot([], [], [], color=col, label=vname, linewidth=3)

    for vname, ax_polys in s_poly_groups.items():
        col, alpha = _STRUCT_PALETTE.get(vname, _STRUCT_DEFAULT)
        for ax2d, polys in ax_polys.items():
            for k, poly2d in enumerate(polys):
                ax2d.add_patch(MplPolygon(poly2d, closed=True,
                                          facecolor=col, edgecolor=col,
                                          alpha=alpha, linewidth=0.3,
                                          label=vname if k == 0 else "_", zorder=1))

    # ---- generic box-drawing helper ----
    def _draw_boxes(vols: list[VolumeEntry], col: str, label_pfx: str, zorder3d: int,
                    show_count: bool = True):
        for i, v in enumerate(vols):
            corners   = box_corners_world(v.pos, v.ext, v.rot)
            lbl       = (f"{label_pfx} ({len(vols)})" if show_count else label_pfx) \
                        if i == 0 else "_"
            faces_3d  = [[[corners[j][o3[0]], corners[j][o3[1]], corners[j][o3[2]]]
                           for j in fi] for fi in BOX_FACES]
            ax3.add_collection3d(Poly3DCollection(
                faces_3d, alpha=_FACE_ALPHA,
                facecolor=col, edgecolor=col, linewidth=0.4))
            if i == 0:
                ax3.plot([], [], [], color=col, label=lbl, linewidth=3)
            for ax2d, ah, av, _ in views2d:
                poly2d = project_polygon(corners, ah, av)
                ax2d.add_patch(MplPolygon(poly2d, closed=True,
                                          facecolor=col, edgecolor=col,
                                          alpha=_FACE_ALPHA, linewidth=0.4,
                                          label=lbl, zorder=zorder3d))
                lbl = "_"

    _draw_boxes(anodes, _ANODE_COL, "Anodes",    zorder3d=3)
    _draw_boxes(cathodes, _CATHODE_COL, "Cathodes", zorder3d=3)
    # one ring = several boxes: count placements, not boxes
    _draw_boxes(fshapers, _FSHAPER_COL,
                f"Field shapers ({_n_placements(fshapers)})", zorder3d=2,
                show_count=False)
    _draw_boxes(opdets, _OPDET_COL, "Opt. det.", zorder3d=4)
    picked_radio = picks and any(parse_pick_spec(s)[0].search(v.name)
                                 for s in picks for v in radios)
    if (surfaces or picked_radio) and radios:
        # one colour per selected surface (producer#index), or per producer
        # when there are too many surfaces for a readable legend
        groups: dict[str, list[VolumeEntry]] = defaultdict(list)
        for v in radios:
            groups[_surface_tag(v)].append(v)
        if len(groups) > _MAX_SURFACE_GROUPS:
            groups = defaultdict(list)
            for v in radios:
                groups[v.name].append(v)
        per_surface = all(len(g) == 1 for g in groups.values())
        for k, (lbl, gvols) in enumerate(groups.items()):
            _draw_boxes(gvols, _SURFACE_PALETTE[k % len(_SURFACE_PALETTE)],
                        lbl, zorder3d=2, show_count=not per_surface)
    else:
        _draw_boxes(radios, _RADIO_COL, "Radio vol.", zorder3d=2)

    if lo is not None:
        span = hi - lo
        ax3.set_xlim(lo[o3[0]], hi[o3[0]])
        if layout["flip3d_y"]:
            ax3.set_ylim(hi[o3[1]], lo[o3[1]])
        else:
            ax3.set_ylim(lo[o3[1]], hi[o3[1]])
        ax3.set_zlim(lo[o3[2]], hi[o3[2]])
        ax3.set_box_aspect((span[o3[0]], span[o3[1]], span[o3[2]]))

    names = "xyz"
    ax3.set_xlabel(f"{names[o3[0]]} [cm]")
    ax3.set_ylabel(f"{names[o3[1]]} [cm]")
    ax3.set_zlabel(f"{names[o3[2]]} [cm]")
    ax3.view_init(elev=view[0], azim=view[1])
    ax3.set_proj_type("ortho")
    ax3.set_title("3D (ortho, proportional)")
    ax3.legend(fontsize=7, loc="upper right", frameon=False)

    # All 2D views share the same world bounding box, so matching axes line
    # up across panels; when proportional they are also drawn to scale
    # (1 cm == 1 cm on both axes), otherwise they fill their grid cell.
    labels = layout["labels"]
    for ax2d, ah, av, ttl in views2d:
        ax2d.set_xlabel(f"{labels[ah]} [cm]")
        ax2d.set_ylabel(f"{labels[av]} [cm]")
        ax2d.set_title(ttl)
        if lo is not None:
            ax2d.set_xlim(lo[ah], hi[ah])
            ax2d.set_ylim(lo[av], hi[av])
        else:
            ax2d.autoscale()
        ax2d.set_aspect("equal" if proportional else "auto", adjustable="box")
        # above the highlight collections (zorder 10+)
        ax2d.legend(fontsize=7, frameon=False).set_zorder(100)
        ax2d.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"[info] plot written to {out_path}", file=sys.stderr)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Listing helpers
# ---------------------------------------------------------------------------

def _hl_mark(v, marks: dict | None) -> str:
    """'*k' for a volume highlighted by spec k (1-based), else ''."""
    k = marks.get(id(v)) if marks else None
    return f"*{k + 1}" if k is not None else ""


def _print_list(volumes: list, marks: dict | None = None) -> None:
    """Compact summary: kinds → unique names → placement count.

    `marks` maps id(volume) → highlight spec index; highlighted placements
    are counted per name and radio surfaces are flagged *k."""
    names_by_kind: dict[str, list[str]] = {}
    for v in volumes:
        names_by_kind.setdefault(v.kind, [])
        if v.name not in names_by_kind[v.kind]:
            names_by_kind[v.kind].append(v.name)
    for kind in sorted(names_by_kind):
        kvols = [v for v in volumes if v.kind == kind]
        print(f"{kind}  ({_n_placements(kvols)} placements)")
        for name in names_by_kind[kind]:
            vs = [v for v in kvols if v.name == name]
            n  = _n_placements(vs)
            hl = [v for v in vs if _hl_mark(v, marks)]
            tag = ""
            if hl:
                specs = sorted({_hl_mark(v, marks) for v in hl})
                tag = f"  {_n_placements(hl)} highlighted {' '.join(specs)}"
            print(f"  {name}  ({n}" + (f", {len(vs)} boxes)" if len(vs) != n else ")")
                  + tag)
            if kind == "radio":
                # one line per generation surface, for use with --surface
                for v in vs:
                    m = v.meta
                    print(f"  {_hl_mark(v, marks):>2}#{m.get('index', '?'):<4}"
                          f" {m.get('surface', '-'):<22}"
                          f" {m.get('nuclide') or '-':<16}"
                          f" {m.get('node') or m.get('source', '-')}")


def _print_list_detail(volumes: list, marks: dict | None = None) -> None:
    """One line per placement: kind, name, and x/y/z ranges in the IR frame
    (column hl: *k when highlighted by spec k).  The path column is the
    physvol-name path (meta.pvpath) when the GDML names its physvols, else
    the logical-volume path (meta.path)."""
    # header
    print(f"{'#':>5} {'hl':>3}  {'kind':<11}  {'name':<44}  {'surface':<22}"
          f"  {'x_min':>10}  {'x_max':>10}"
          f"  {'y_min':>10}  {'y_max':>10}"
          f"  {'z_min':>10}  {'z_max':>10}  cm  path")
    print("-" * 167)
    for idx, v in enumerate(volumes):
        corners = box_corners_world(v.pos, v.ext, v.rot)
        lo, hi  = corners.min(axis=0), corners.max(axis=0)
        print(f"{idx:5d} {_hl_mark(v, marks):>3}  {v.kind:<11}  {_surface_tag(v):<44}"
              f"  {v.meta.get('surface', v.meta.get('part', '')):<22}"
              f"  {lo[0]:10.3f}  {hi[0]:10.3f}"
              f"  {lo[1]:10.3f}  {hi[1]:10.3f}"
              f"  {lo[2]:10.3f}  {hi[2]:10.3f}      {v.meta.get('pvpath') or v.meta.get('path', '')}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("ir_file", help="JSON IR produced by bki-extract")
    ap.add_argument("--list", "-l", action="store_true",
                    help="print volume kinds and names grouped by kind (radio "
                         "producers with one line per surface), then exit "
                         "(respects selectors)")
    ap.add_argument("--list-detail", "-L", action="store_true",
                    help="print every placement with its x/y/z ranges in the IR "
                         "frame and its path (physvol names when the GDML has "
                         "them), then exit (respects selectors)")
    ap.add_argument("--out", "-o", default=None,
                    help="output plot file (default: <stem>.pdf)")
    ap.add_argument("--view", nargs=2, type=float, default=(20.0, -60.0),
                    metavar=("ELEV", "AZIM"),
                    help="3D view elevation and azimuth in degrees (default: 20 -60)")
    ap.add_argument("--title", default=None,
                    help="plot title (default: derived from IR source)")
    ap.add_argument("--proportional", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="draw the 2D views to scale (1 cm == 1 cm on both "
                         "axes); --no-proportional stretches each view to "
                         "fill its panel (default: proportional)")
    ap.add_argument("--orientation", choices=("hd", "vd", "auto"), default="hd",
                    help="hd: drift axis x horizontal, y vertical; vd: x "
                         "vertical in every view; auto: vd for VD IRs, hd "
                         "otherwise (default: hd)")
    # ---- volume selectors ----
    grp = ap.add_argument_group("volume selectors")
    grp.add_argument("--no-detector", action="store_true",
                     help="exclude detector volumes (anode, opdet, cathode, "
                          "fieldshaper, struct)")
    grp.add_argument("--no-radio", action="store_true",
                     help="exclude radiological background volumes (radio)")
    grp.add_argument("--kind", nargs="+", default=None,
                     metavar="KIND",
                     help="show only these kinds: anode opdet cathode fieldshaper "
                          "struct radio "
                          "(applied after --no-detector / --no-radio)")
    grp.add_argument("--name", nargs="+", default=None,
                     metavar="PATTERN",
                     help="show only volumes whose name matches at least one "
                          "regex across all kinds (re.search)")
    grp.add_argument("--anode", nargs="+", default=None, metavar="PATTERN",
                     help="regex filter on anode volume names")
    grp.add_argument("--opdet", nargs="+", default=None, metavar="PATTERN",
                     help="regex filter on opdet volume names")
    grp.add_argument("--cathode", nargs="+", default=None, metavar="PATTERN",
                     help="regex filter on cathode volume names")
    grp.add_argument("--fieldshaper", nargs="+", default=None, metavar="PATTERN",
                     help="regex filter on field-shaper volume names")
    grp.add_argument("--struct", nargs="+", default=None, metavar="PATTERN",
                     help="regex filter on struct volume names")
    grp.add_argument("--radio", nargs="+", default=None, metavar="PATTERN",
                     help="regex filter on radio producer names")
    grp.add_argument("--surface", nargs="+", default=None, metavar="SPEC",
                     help="radio generation-surface selector PRODUCER_RE#SEL; "
                          "SEL is an index list (0,2,5-7) or @REGEX matched "
                          "against surface label, node and nuclide.  Producers "
                          "not matched by any spec are left unfiltered; "
                          "selected surfaces are coloured individually")
    grp.add_argument("--pick", nargs="+", default=None, metavar="SPEC",
                     help="placement selector NAME_RE#@REGEX for any kind: "
                          "keep placements of matching volumes whose path "
                          "(shown by -L, e.g. .../volTPC:37/volTPCPlaneV:0) "
                          "or part matches REGEX; for radio entries it works "
                          "as --surface.  Volumes not matched by any NAME_RE "
                          "are left unfiltered")
    hgrp = ap.add_argument_group("highlighting (applied after the selectors)")
    hgrp.add_argument("--highlight", nargs="+", default=None, metavar="SPEC",
                      help="[LABEL=]NAME_RE[#SEL]: draw the matching volumes on top, "
                           "one colour and legend entry per spec, and the rest as "
                           "grey context.  SEL as for --pick (optional; @REGEX "
                           "also matches the physvol name pv and pvpath)")
    hgrp.add_argument("--dim", type=float, default=_DIM_DEFAULT, metavar="F",
                      help=f"face opacity of the context volumes (default: {_DIM_DEFAULT})")
    hgrp.add_argument("--zoom-highlight", action="store_true",
                      help="fit the axes to the highlighted volumes")
    args = ap.parse_args()

    ir_path = Path(args.ir_file)
    if not ir_path.is_file():
        sys.exit(f"[error] '{args.ir_file}' not found")

    if ir_path.suffix.lower() != ".json":
        sys.exit(f"[error] '{ir_path.name}' is not a JSON IR file — "
                 f"run bki-extract first to produce one")

    ir = IRCollection.load(ir_path)
    print(f"[info] loaded {len(ir.volumes)} volumes from {ir_path.name}", file=sys.stderr)

    # build effective kinds list from coarse + fine selectors
    _DET_KINDS   = {"anode", "opdet", "cathode", "fieldshaper", "struct"}
    _RADIO_KINDS = {"radio"}
    active: set[str] = _DET_KINDS | _RADIO_KINDS
    if args.no_detector:
        active -= _DET_KINDS
    if args.no_radio:
        active -= _RADIO_KINDS
    # --kind refines further (intersection with active, or exact list if --no-* not used)
    if args.kind:
        requested = set(args.kind)
        unknown   = requested - (_DET_KINDS | _RADIO_KINDS)
        if unknown:
            print(f"[warn] unknown kind(s): {', '.join(sorted(unknown))}; "
                  f"valid kinds are: anode opdet cathode fieldshaper struct radio", file=sys.stderr)
        active &= requested
    kinds = sorted(active) if active != (_DET_KINDS | _RADIO_KINDS) else None

    out_path = args.out or str(ir_path.with_suffix(".pdf"))
    title    = args.title or (
        f"{Path(ir.source).name}  —  {ir.detector.upper()} detector layout"
        if ir.source else ir_path.stem
    )

    kind_patterns: dict[str, list[str]] = {}
    for kind, flag in (("anode", args.anode), ("opdet", args.opdet),
                       ("cathode", args.cathode),
                       ("fieldshaper", args.fieldshaper), ("struct", args.struct),
                       ("radio", args.radio)):
        if flag:
            kind_patterns[kind] = flag

    for opt, specs, parse in (("--surface", args.surface, parse_pick_spec),
                              ("--pick", args.pick, parse_pick_spec),
                              ("--highlight", args.highlight, parse_highlight_spec)):
        for spec in specs or ():
            try:
                parse(spec)
            except (ValueError, re.error) as exc:
                sys.exit(f"[error] {opt}: {exc}")

    # same filters for the listings and the plot
    try:
        vols = select_volumes(ir, kinds=kinds, name_patterns=args.name,
                              kind_patterns=kind_patterns or None,
                              surfaces=args.surface, picks=args.pick)
    except ValueError as exc:
        sys.exit(f"[error] {exc}")

    if args.list or args.list_detail:
        marks = None
        if args.highlight:
            try:
                _, groups, _ = highlight_groups(vols, args.highlight)
            except ValueError as exc:
                sys.exit(f"[error] {exc}")
            marks = {id(v): k for k, g in enumerate(groups) for v in g}
        (_print_list_detail if args.list_detail else _print_list)(vols, marks)
        sys.exit(0)

    try:
        plot_ir(ir, out_path, title=title, view=tuple(args.view),
                kinds=kinds, name_patterns=args.name,
                kind_patterns=kind_patterns or None,
                surfaces=args.surface, picks=args.pick,
                proportional=args.proportional,
                orientation=args.orientation,
                highlight=args.highlight, dim=args.dim,
                zoom_highlight=args.zoom_highlight)
    except ValueError as exc:
        sys.exit(f"[error] {exc}")


if __name__ == "__main__":
    main()

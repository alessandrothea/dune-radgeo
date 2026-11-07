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
import sys
from collections import defaultdict
from typing import Any, NamedTuple
from pathlib import Path

import numpy as np

from ._draw import BOX_FACES, aabb, box_corners_world, project_polygon
from ._select import highlight_groups, parse_pick_spec, select_volumes, validate_picks
from ._style import (ANODE_COL, CATHODE_COL, CONTEXT_COL, DIM_DEFAULT, FACE_ALPHA,
                     FSHAPER_COL, HL_FACE_ALPHA, HL_LINEWIDTH, MAX_SURFACE_GROUPS,
                     OPDET_COL, RADIO_COL, STRUCT_DEFAULT, STRUCT_PALETTE,
                     SURFACE_PALETTE, THIN_FRAC, import_mpl)
from .ir import IRCollection, VolumeEntry


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
# Helpers
# ---------------------------------------------------------------------------

def _world_bounds(vols: list[VolumeEntry], pad_frac: float = 0.03) -> tuple:
    """
    Padded (lo, hi) [x, y, z] in the IR frame enclosing all corners of
    `vols`, or (None, None) when empty.  The 2D views draw projections of
    these same corners, so the box also encloses every drawn 2D polygon.
    """
    if not vols:
        return None, None
    boxes = [aabb(v.pos, v.ext, v.rot) for v in vols]
    lo = np.min([b[0] for b in boxes], axis=0)
    hi = np.max([b[1] for b in boxes], axis=0)
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


def _group_by(vols: list[VolumeEntry], key) -> dict[str, list[VolumeEntry]]:
    """`vols` grouped by key(v), groups in order of first appearance."""
    groups: dict[str, list[VolumeEntry]] = defaultdict(list)
    for v in vols:
        groups[key(v)].append(v)
    return groups


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

class _Panels(NamedTuple):
    """The four panels of a plot_ir figure and their shared bounds."""
    ax3: Any         # 3D view
    views2d: list    # (axes, horizontal world axis, vertical world axis, title)
    o3: tuple        # world axes shown as the 3D view's X, Y, Z
    lo: Any          # _world_bounds() of the figure; None when nothing is drawn
    hi: Any


def _faces3d(corners: list, o3: tuple) -> list:
    """The six faces of every box in `corners` ((8,3) arrays), as vertex
    lists in the 3D view's axis order."""
    return [[[c[j][o3[0]], c[j][o3[1]], c[j][o3[2]]] for j in fi]
            for c in corners for fi in BOX_FACES]


def _draw_batch(pn: _Panels, vols: list[VolumeEntry], col: str, label: str,
                face_alpha: float, edge_alpha: float, lw: float, zorder: int, *,
                patches: bool = False, per_box: bool = False,
                thin_as_lines: bool = False) -> None:
    """
    Draw the boxes of `vols` in colour `col` on all four panels, with a
    legend entry `label` in each.

    By default each panel gets one collection and the 2D legend entries are
    lines.  ``patches``: one Polygon per box in the 2D views, the legend
    entry then shows a patch.  ``per_box``: one 3D collection per box, and
    the 2D legend entry in the first 2D view only.  ``thin_as_lines``: boxes
    thinner than THIN_FRAC of a 2D view are drawn as lines (collections only).
    """
    if not vols:
        return
    from matplotlib.collections import LineCollection, PolyCollection
    from matplotlib.colors import to_rgba
    from matplotlib.patches import Polygon
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    face, edge = to_rgba(col, face_alpha), to_rgba(col, edge_alpha)
    corners = [box_corners_world(v.pos, v.ext, v.rot) for v in vols]
    for i, cs in enumerate([[c] for c in corners] if per_box else [corners]):
        pn.ax3.add_collection3d(Poly3DCollection(
            _faces3d(cs, pn.o3), facecolor=face, edgecolor=edge, linewidth=lw,
            zorder=zorder))
        if i == 0:
            pn.ax3.plot([], [], [], color=col, label=label, linewidth=3)

    lo, hi = pn.lo, pn.hi
    for n, (ax2d, ah, av, _) in enumerate(pn.views2d):
        polys = [project_polygon(c, ah, av) for c in corners]
        if patches:
            lbl = label if n == 0 or not per_box else "_"
            for poly in polys:
                ax2d.add_patch(Polygon(poly, closed=True, facecolor=face,
                                       edgecolor=edge, linewidth=lw,
                                       label=lbl, zorder=zorder))
                lbl = "_"
            continue
        lines = []
        if thin_as_lines and lo is not None:
            boxes, polys = polys, []
            for poly in boxes:
                pmin, pmax = poly.min(axis=0), poly.max(axis=0)
                w, h = pmax - pmin
                mid = (pmin + pmax) / 2
                if w < THIN_FRAC * (hi[ah] - lo[ah]):
                    lines.append([(mid[0], pmin[1]), (mid[0], pmax[1])])
                elif h < THIN_FRAC * (hi[av] - lo[av]):
                    lines.append([(pmin[0], mid[1]), (pmax[0], mid[1])])
                else:
                    polys.append(poly)
        if polys:
            ax2d.add_collection(PolyCollection(polys, facecolor=face, edgecolor=edge,
                                               linewidth=lw, zorder=zorder))
        if lines:
            ax2d.add_collection(LineCollection(lines, colors=[to_rgba(col, 1.0)],
                                               linewidth=1.5, zorder=zorder))
        ax2d.plot([], [], color=col, label=label, linewidth=3)


def _draw_highlight(pn: _Panels, context: list[VolumeEntry], labels: list[str],
                    groups: list[list[VolumeEntry]], dim: float) -> None:
    """Grey context volumes, then each highlight group on top in its colour."""
    _draw_batch(pn, context, CONTEXT_COL, f"other ({_n_placements(context)})",
                dim, min(1.0, 5 * dim), 0.2, 1)
    for k, (lbl, g) in enumerate(zip(labels, groups)):
        _draw_batch(pn, g, SURFACE_PALETTE[k % len(SURFACE_PALETTE)],
                    f"{lbl} ({_n_placements(g)})", HL_FACE_ALPHA, 1.0,
                    HL_LINEWIDTH, 10 + k, thin_as_lines=True)


def _draw_kinds(pn: _Panels, vols: list[VolumeEntry], colour_surfaces: bool) -> None:
    """Volumes coloured by kind: struct per name, then anodes, cathodes,
    field shapers, optical detectors and radio volumes (one colour per
    surface, or per producer, when `colour_surfaces`)."""
    by_kind = _group_by(vols, lambda v: v.kind)

    for name, g in _group_by(by_kind["struct"], lambda v: v.name).items():
        col, alpha = STRUCT_PALETTE.get(name, STRUCT_DEFAULT)
        _draw_batch(pn, g, col, name, alpha, alpha, 0.3, 1, patches=True)

    def boxes(g, col, label, zorder):
        _draw_batch(pn, g, col, label, FACE_ALPHA, FACE_ALPHA, 0.4, zorder,
                    patches=True, per_box=True)

    boxes(by_kind["anode"], ANODE_COL, f"Anodes ({len(by_kind['anode'])})", 3)
    boxes(by_kind["cathode"], CATHODE_COL,
          f"Cathodes ({len(by_kind['cathode'])})", 3)
    # one ring = several boxes: count placements, not boxes
    boxes(by_kind["fieldshaper"], FSHAPER_COL,
          f"Field shapers ({_n_placements(by_kind['fieldshaper'])})", 2)
    boxes(by_kind["opdet"], OPDET_COL, f"Opt. det. ({len(by_kind['opdet'])})", 4)

    radios = by_kind["radio"]
    if not colour_surfaces:
        boxes(radios, RADIO_COL, f"Radio vol. ({len(radios)})", 2)
        return
    # one colour per selected surface (producer#index), or per producer when
    # there are too many surfaces for a readable legend
    groups = _group_by(radios, _surface_tag)
    if len(groups) > MAX_SURFACE_GROUPS:
        groups = _group_by(radios, lambda v: v.name)
    per_surface = all(len(g) == 1 for g in groups.values())
    for k, (lbl, g) in enumerate(groups.items()):
        boxes(g, SURFACE_PALETTE[k % len(SURFACE_PALETTE)],
              lbl if per_surface else f"{lbl} ({len(g)})", 2)


def _setup_axes(pn: _Panels, layout: dict, view: tuple, proportional: bool) -> None:
    """Limits, labels, titles and legends of the four panels."""
    ax3, o3, lo, hi = pn.ax3, pn.o3, pn.lo, pn.hi
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
    for ax2d, ah, av, ttl in pn.views2d:
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
    dim: float = DIM_DEFAULT,
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

    plt = import_mpl()

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

    fig = plt.figure(figsize=(16, 10))
    if title:
        fig.suptitle(title, fontsize=11, fontweight="bold")
    # axis limits enclose every drawn volume, or the highlights
    hl_vols = [v for g in hl_groups for v in g]
    lo, hi = _world_bounds(hl_vols if zoom_highlight and hl_vols else all_vols)
    pn = _Panels(ax3=fig.add_subplot(2, 2, 1, projection="3d"),
                 views2d=[(fig.add_subplot(2, 2, k), ah, av, ttl)
                          for k, (ah, av, ttl) in zip((2, 3, 4), layout["views"])],
                 o3=layout["order3d"], lo=lo, hi=hi)

    if highlight:
        _draw_highlight(pn, vols, hl_labels, hl_groups, dim)
    else:
        picked_radio = picks and any(parse_pick_spec(s)[0].search(v.name)
                                     for s in picks for v in vols
                                     if v.kind == "radio")
        _draw_kinds(pn, vols, colour_surfaces=bool(surfaces or picked_radio))
    _setup_axes(pn, layout, view, proportional)

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
        lo, hi = aabb(v.pos, v.ext, v.rot)
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
    hgrp.add_argument("--dim", type=float, default=DIM_DEFAULT, metavar="F",
                      help=f"face opacity of the context volumes (default: {DIM_DEFAULT})")
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

    # same filters for the listings and the plot
    try:
        validate_picks(surface=args.surface, pick=args.pick, highlight=args.highlight)
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

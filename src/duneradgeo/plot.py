"""
bki-plot  —  render a JSON IR as a matplotlib PDF

Produces the same 3D + 3×2D layout as the original plot_gdml_detectors.py:
  top-left   3D ortho view (proportional axes)
  top-right  Z–Y  beam face
  bot-left   Z–X  top view
  bot-right  X–Y  side view

Usage:
  bki-plot detector.json --out detector.pdf
  bki-plot detector.json --out detector.pdf --view 20 -60
"""
from __future__ import annotations

import argparse
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

    Axes convention: x = drift, y = vertical, z = beam.
    The 3D axes are re-ordered as (z, x, y) to match the original script.
    """
    import re as _re

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    # ---- apply filters ----
    vols = ir.volumes
    if kinds is not None:
        kind_set = set(kinds)
        vols = [v for v in vols if v.kind in kind_set]
    if name_patterns:
        pats = [_re.compile(p) for p in name_patterns]
        vols = [v for v in vols if any(p.search(v.name) for p in pats)]
    if kind_patterns:
        compiled = {k: [_re.compile(p) for p in ps]
                    for k, ps in kind_patterns.items()}
        vols = [v for v in vols
                if v.kind not in compiled
                or any(p.search(v.name) for p in compiled[v.kind])]

    filtered = IRCollection(source=ir.source, detector=ir.detector, volumes=vols)

    anodes  = filtered.by_kind("anode")
    opdets  = filtered.by_kind("opdet")
    structs = filtered.by_kind("struct")
    radios  = filtered.by_kind("radio")

    fig = plt.figure(figsize=(16, 10))
    if title:
        fig.suptitle(title, fontsize=11, fontweight="bold")

    ax3   = fig.add_subplot(2, 2, 1, projection="3d")
    ax_zy = fig.add_subplot(2, 2, 2)
    ax_zx = fig.add_subplot(2, 2, 3)
    ax_xy = fig.add_subplot(2, 2, 4)

    # ---- structural volumes (batched per name) ----
    s_face_groups: dict[str, list]       = defaultdict(list)
    s_poly_groups: dict[str, dict]       = {}

    for v in structs:
        corners = box_corners_world(v.pos, v.ext, v.rot)
        for fi in BOX_FACES:
            s_face_groups[v.name].append(
                [[corners[j][2], corners[j][0], corners[j][1]] for j in fi])
        if v.name not in s_poly_groups:
            s_poly_groups[v.name] = {ax_zy: [], ax_zx: [], ax_xy: []}
        for ax2d, ah, av in ((ax_zy, 2, 1), (ax_zx, 2, 0), (ax_xy, 0, 1)):
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
    def _draw_boxes(vols: list[VolumeEntry], col: str, label_pfx: str, zorder3d: int):
        for i, v in enumerate(vols):
            corners   = box_corners_world(v.pos, v.ext, v.rot)
            lbl       = f"{label_pfx} ({len(vols)})" if i == 0 else "_"
            faces_3d  = [[[corners[j][2], corners[j][0], corners[j][1]]
                           for j in fi] for fi in BOX_FACES]
            ax3.add_collection3d(Poly3DCollection(
                faces_3d, alpha=_FACE_ALPHA,
                facecolor=col, edgecolor=col, linewidth=0.4))
            if i == 0:
                ax3.plot([], [], [], color=col, label=lbl, linewidth=3)
            for ax2d, ah, av in ((ax_zy, 2, 1), (ax_zx, 2, 0), (ax_xy, 0, 1)):
                poly2d = project_polygon(corners, ah, av)
                ax2d.add_patch(MplPolygon(poly2d, closed=True,
                                          facecolor=col, edgecolor=col,
                                          alpha=_FACE_ALPHA, linewidth=0.4,
                                          label=lbl, zorder=zorder3d))
                lbl = "_"

    _draw_boxes(anodes, _ANODE_COL, "Anodes",    zorder3d=3)
    _draw_boxes(opdets, _OPDET_COL, "Opt. det.", zorder3d=4)
    _draw_boxes(radios, _RADIO_COL, "Radio vol.", zorder3d=2)

    # ---- axis limits (based on anode+opdet corners) ----
    ref_vols = anodes + opdets or structs or radios
    all_corners = [box_corners_world(v.pos, v.ext, v.rot) for v in ref_vols]
    if all_corners:
        arr = np.vstack(all_corners)
        lo, hi = arr.min(axis=0), arr.max(axis=0)
        pad = 0.03 * np.maximum(hi - lo, 1.0)
        lo, hi = lo - pad, hi + pad
        span = hi - lo
        ax3.set_xlim(lo[2], hi[2])
        ax3.set_ylim(lo[0], hi[0])
        ax3.set_zlim(lo[1], hi[1])
        ax3.set_box_aspect((span[2], span[0], span[1]))

    ax3.set_xlabel("z [cm]"); ax3.set_ylabel("x [cm]"); ax3.set_zlabel("y [cm]")
    ax3.view_init(elev=view[0], azim=view[1])
    ax3.set_proj_type("ortho")
    ax3.set_title("3D (ortho, proportional)")
    ax3.legend(fontsize=7, loc="upper right", frameon=False)

    for ax2d, xl, yl, ttl in (
        (ax_zy, "z (beam) [cm]",  "y (vert) [cm]",  "Z–Y  (beam face)"),
        (ax_zx, "z (beam) [cm]",  "x (drift) [cm]", "Z–X  (top view)"),
        (ax_xy, "x (drift) [cm]", "y (vert) [cm]",  "X–Y  (side view)"),
    ):
        ax2d.set_xlabel(xl); ax2d.set_ylabel(yl); ax2d.set_title(ttl)
        ax2d.autoscale()
        ax2d.legend(fontsize=7, frameon=False)
        ax2d.grid(True, alpha=0.25)
    ax_zy.set_aspect("equal")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"[info] plot written to {out_path}", file=sys.stderr)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Listing helpers
# ---------------------------------------------------------------------------

def _print_list(volumes: list) -> None:
    """Compact summary: kinds → unique names → placement count."""
    from collections import Counter
    counts      = Counter(v.kind for v in volumes)
    names_by_kind: dict[str, list[str]] = {}
    for v in volumes:
        names_by_kind.setdefault(v.kind, [])
        if v.name not in names_by_kind[v.kind]:
            names_by_kind[v.kind].append(v.name)
    for kind in sorted(names_by_kind):
        print(f"{kind}  ({counts[kind]} placements)")
        for name in names_by_kind[kind]:
            n = sum(1 for v in volumes if v.kind == kind and v.name == name)
            print(f"  {name}  ({n})")


def _print_list_detail(volumes: list) -> None:
    """One line per placement: kind, name, and world-frame x/y/z ranges."""
    from ._draw import box_corners_world
    # header
    print(f"{'#':>5}  {'kind':<8}  {'name':<44}"
          f"  {'x_min':>10}  {'x_max':>10}"
          f"  {'y_min':>10}  {'y_max':>10}"
          f"  {'z_min':>10}  {'z_max':>10}  cm")
    print("-" * 130)
    for idx, v in enumerate(volumes):
        corners = box_corners_world(v.pos, v.ext, v.rot)
        lo, hi  = corners.min(axis=0), corners.max(axis=0)
        print(f"{idx:5d}  {v.kind:<8}  {v.name:<44}"
              f"  {lo[0]:10.3f}  {hi[0]:10.3f}"
              f"  {lo[1]:10.3f}  {hi[1]:10.3f}"
              f"  {lo[2]:10.3f}  {hi[2]:10.3f}")


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
                    help="print volume kinds and names grouped by kind, then exit")
    ap.add_argument("--list-detail", "-L", action="store_true",
                    help="print every placement with pos, ext and 8 world-frame "
                         "corner coordinates, then exit (respects selectors)")
    ap.add_argument("--out", "-o", default=None,
                    help="output plot file (default: <stem>.pdf)")
    ap.add_argument("--view", nargs=2, type=float, default=(20.0, -60.0),
                    metavar=("ELEV", "AZIM"),
                    help="3D view elevation and azimuth in degrees (default: 20 -60)")
    ap.add_argument("--title", default=None,
                    help="plot title (default: derived from IR source)")
    # ---- volume selectors ----
    grp = ap.add_argument_group("volume selectors")
    grp.add_argument("--no-detector", action="store_true",
                     help="exclude detector volumes (anode, opdet, struct)")
    grp.add_argument("--no-radio", action="store_true",
                     help="exclude radiological background volumes (radio)")
    grp.add_argument("--kind", nargs="+", default=None,
                     metavar="KIND",
                     help="show only these kinds: anode opdet struct radio "
                          "(applied after --no-detector / --no-radio)")
    grp.add_argument("--name", nargs="+", default=None,
                     metavar="PATTERN",
                     help="show only volumes whose name matches at least one "
                          "regex across all kinds (re.search, applied last)")
    grp.add_argument("--anode", nargs="+", default=None, metavar="PATTERN",
                     help="regex filter on anode volume names")
    grp.add_argument("--opdet", nargs="+", default=None, metavar="PATTERN",
                     help="regex filter on opdet volume names")
    grp.add_argument("--struct", nargs="+", default=None, metavar="PATTERN",
                     help="regex filter on struct volume names")
    grp.add_argument("--radio", nargs="+", default=None, metavar="PATTERN",
                     help="regex filter on radio producer names")
    args = ap.parse_args()

    ir_path = Path(args.ir_file)
    if not ir_path.is_file():
        sys.exit(f"[error] '{args.ir_file}' not found")

    ir = IRCollection.load(ir_path)
    print(f"[info] loaded {len(ir.volumes)} volumes from {ir_path.name}", file=sys.stderr)

    if args.list:
        _print_list(ir.volumes)
        sys.exit(0)

    # build effective kinds list from coarse + fine selectors
    _DET_KINDS   = {"anode", "opdet", "struct"}
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
                  f"valid kinds are: anode opdet struct radio", file=sys.stderr)
        active &= requested
    kinds = sorted(active) if active != (_DET_KINDS | _RADIO_KINDS) else None

    out_path = args.out or str(ir_path.with_suffix(".pdf"))
    title    = args.title or (
        f"{Path(ir.source).name}  —  {ir.detector.upper()} detector layout"
        if ir.source else ir_path.stem
    )

    kind_patterns: dict[str, list[str]] = {}
    for kind, flag in (("anode", args.anode), ("opdet", args.opdet),
                       ("struct", args.struct), ("radio", args.radio)):
        if flag:
            kind_patterns[kind] = flag

    if args.list_detail:
        # apply the same filters as the plot would, then print detail
        import re as _re
        vols = ir.volumes
        if kinds is not None:
            vols = [v for v in vols if v.kind in set(kinds)]
        if args.name:
            pats = [_re.compile(p) for p in args.name]
            vols = [v for v in vols if any(p.search(v.name) for p in pats)]
        if kind_patterns:
            compiled = {k: [_re.compile(p) for p in ps]
                        for k, ps in kind_patterns.items()}
            vols = [v for v in vols
                    if v.kind not in compiled
                    or any(p.search(v.name) for p in compiled[v.kind])]
        _print_list_detail(vols)
        sys.exit(0)

    plot_ir(ir, out_path, title=title, view=tuple(args.view),
            kinds=kinds, name_patterns=args.name,
            kind_patterns=kind_patterns or None)


if __name__ == "__main__":
    main()

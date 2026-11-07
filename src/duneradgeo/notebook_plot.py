"""
Plotly-based interactive plotting for Jupyter notebooks.

Public API:
  make_figure(ir)          → plotly Figure with 3D + 3×2D sub-plots
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from ._draw import BOX_EDGES, BOX_FACES, box_corners_world, project_polygon
from ._select import highlight_groups
from ._style import CONTEXT_COL, DIM_DEFAULT, HL_FACE_ALPHA, KIND_COLOUR, SURFACE_PALETTE
from .ir import IRCollection

# ---------------------------------------------------------------------------
# Colour maps
# ---------------------------------------------------------------------------

# bki-plot's kind colours; struct volumes, coloured per name there, get one colour
_KIND_COLOUR: dict[str, str] = dict(KIND_COLOUR, struct="mediumseagreen")
_KIND_OPACITY: dict[str, float] = {
    "anode":  0.25,
    "opdet":  0.20,
    "cathode": 0.25,
    "fieldshaper": 0.30,
    "struct": 0.10,
    "radio":  0.20,
}


def _hex_to_rgba(col: str, alpha: float) -> str:
    """Convert named/hex colour to 'rgba(r,g,b,a)' accepted by Plotly."""
    import matplotlib.colors as mc
    r, g, b, _ = mc.to_rgba(col)
    return f"rgba({int(r*255)},{int(g*255)},{int(b*255)},{alpha:.2f})"


# ---------------------------------------------------------------------------
# Single-box Plotly traces
# ---------------------------------------------------------------------------

def _box_mesh_trace(center: np.ndarray, extents: np.ndarray, rot: np.ndarray,
                    colour: str, opacity: float, name: str,
                    legendgroup: str, showlegend: bool) -> Any:
    """Return a go.Mesh3d trace for one box."""
    import plotly.graph_objects as go

    c   = box_corners_world(center, extents, rot)
    # Mesh3d triangles from quads
    i_idx, j_idx, k_idx = [], [], []
    for face in BOX_FACES:
        a, b, cc, d = face
        i_idx += [a, a]
        j_idx += [b, cc]
        k_idx += [cc, d]

    return go.Mesh3d(
        x=c[:, 2], y=c[:, 0], z=c[:, 1],   # (z, x, y) → beam, drift, vert
        i=i_idx, j=j_idx, k=k_idx,
        color=colour, opacity=opacity,
        flatshading=True,
        name=name, legendgroup=legendgroup,
        showlegend=showlegend,
        hovertemplate=(
            f"<b>{name}</b><br>"
            f"pos  (x,y,z) = ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}) cm<br>"
            f"ext  (dx,dy,dz) = ({extents[0]:.1f}, {extents[1]:.1f}, {extents[2]:.1f}) cm"
            "<extra></extra>"
        ),
    )


def _box_line_trace(center: np.ndarray, extents: np.ndarray, rot: np.ndarray,
                    colour: str, legendgroup: str) -> Any:
    """Return a go.Scatter3d edge-wireframe trace for one box."""
    import plotly.graph_objects as go

    c    = box_corners_world(center, extents, rot)
    xs, ys, zs = [], [], []
    for a, b in BOX_EDGES:
        xs += [c[a][2], c[b][2], None]
        ys += [c[a][0], c[b][0], None]
        zs += [c[a][1], c[b][1], None]

    return go.Scatter3d(
        x=xs, y=ys, z=zs,
        mode="lines",
        line=dict(color=colour, width=1),
        legendgroup=legendgroup,
        showlegend=False,
        hoverinfo="skip",
    )


# ---------------------------------------------------------------------------
# 2D projection traces
# ---------------------------------------------------------------------------

def _box_proj_trace(center: np.ndarray, extents: np.ndarray, rot: np.ndarray,
                    colour: str, opacity: float,
                    ax_h: int, ax_v: int,
                    name: str, legendgroup: str, showlegend: bool) -> Any:
    """Return a go.Scatter (filled polygon) for one box projected onto a 2D plane."""
    import plotly.graph_objects as go

    pts   = project_polygon(box_corners_world(center, extents, rot), ax_h, ax_v)
    xs    = list(pts[:, 0]) + [pts[0, 0]]
    ys    = list(pts[:, 1]) + [pts[0, 1]]
    rgba  = _hex_to_rgba(colour, opacity)
    return go.Scatter(
        x=xs, y=ys,
        fill="toself",
        fillcolor=rgba,
        line=dict(color=colour, width=0.5),
        name=name, legendgroup=legendgroup,
        showlegend=showlegend,
        hoverinfo="skip",
        mode="lines",
    )


# ---------------------------------------------------------------------------
# Main figure builder
# ---------------------------------------------------------------------------

def make_figure(ir: IRCollection, height: int = 900,
                highlight: list[str] | None = None, dim: float = DIM_DEFAULT) -> Any:
    """
    Build an interactive Plotly figure with:
      - sub-plot (1,1)  →  3D view  (z=beam, x=drift, y=vert)
      - sub-plot (1,2)  →  Z–Y  (beam face)
      - sub-plot (2,1)  →  Z–X  (top view)
      - sub-plot (2,2)  →  X–Y  (side view)

    All sub-plots share one legend (volumes grouped by kind + name).

    With ``highlight`` (``[LABEL=]NAME_RE[#SEL]`` specs, as for bki-plot
    --highlight) the matching volumes get one colour and legend group per
    spec and all others are drawn in grey with opacity ``dim``.
    """
    from plotly.subplots import make_subplots

    fig = make_subplots(
        rows=2, cols=2,
        column_widths=[0.55, 0.45],
        row_heights=[0.5, 0.5],
        specs=[
            [{"type": "scene"}, {"type": "xy"}],
            [{"type": "xy"},    {"type": "xy"}],
        ],
        subplot_titles=["3D (z=beam, x=drift, y=vert)",
                        "Z–Y  (beam face)",
                        "Z–X  (top view)",
                        "X–Y  (side view)"],
    )

    seen_groups: set[str] = set()

    # (colour, opacity, legend name/group) per volume
    styles: dict[int, tuple] = {}
    if highlight:
        from matplotlib.colors import to_hex   # plotly has no "tab:" colours
        labels, groups, context = highlight_groups(ir.volumes, highlight)
        for v in context:
            styles[id(v)] = (CONTEXT_COL, dim, "other")
        for k, (lbl, g) in enumerate(zip(labels, groups)):
            col = to_hex(SURFACE_PALETTE[k % len(SURFACE_PALETTE)])
            for v in g:
                styles[id(v)] = (col, HL_FACE_ALPHA, lbl)
        # highlights last, so they are drawn on top
        order = context + [v for g in groups for v in g]
    else:
        order = ir.volumes

    for v in order:
        if highlight:
            col, alpha, group = styles[id(v)]
            name = group
        else:
            col     = _KIND_COLOUR.get(v.kind, "grey")
            alpha   = _KIND_OPACITY.get(v.kind, 0.20)
            group   = f"{v.kind}:{v.name}"
            name    = v.name
        is_first = group not in seen_groups
        seen_groups.add(group)

        # 3D
        fig.add_trace(_box_mesh_trace(v.pos, v.ext, v.rot, col, alpha,
                                      name, group, is_first),
                      row=1, col=1)
        fig.add_trace(_box_line_trace(v.pos, v.ext, v.rot, col, group),
                      row=1, col=1)

        # 2D projections: (ax_h, ax_v, row, col)
        proj_specs = [
            (2, 1, 1, 2),  # Z–Y
            (2, 0, 2, 1),  # Z–X
            (0, 1, 2, 2),  # X–Y
        ]
        for ah, av, r, c in proj_specs:
            fig.add_trace(_box_proj_trace(v.pos, v.ext, v.rot, col, alpha,
                                          ah, av, name, group, False),
                          row=r, col=c)

    # ---- 2D axis labels ----
    # scaleanchor uses short-form axis id ('y2', 'y3', …), not 'yaxis2'
    _AXIS_LABEL = {0: "x [cm]", 1: "y [cm]", 2: "z [cm]"}
    proj_axes = [
        ("xaxis2", "yaxis2", "y2", 2, 1),
        ("xaxis3", "yaxis3", "y3", 2, 0),
        ("xaxis4", "yaxis4", "y4", 0, 1),
    ]
    for xk, yk, y_anchor, ah, av in proj_axes:
        fig.update_layout(**{
            xk: dict(title_text=_AXIS_LABEL[ah],
                     scaleanchor=y_anchor, scaleratio=1, showgrid=True),
            yk: dict(title_text=_AXIS_LABEL[av], showgrid=True),
        })

    # ---- 3D scene axes ----
    fig.update_scenes(
        xaxis_title="z [cm]",
        yaxis_title="x [cm]",
        zaxis_title="y [cm]",
        aspectmode="data",
    )

    fig.update_layout(
        height=height,
        title_text=f"{Path(ir.source).name if ir.source else 'bkginspector'}"
                   f"  —  {ir.detector.upper() if ir.detector else 'DUNE'} layout",
        legend=dict(itemsizing="constant", tracegroupgap=2),
        margin=dict(t=60, b=20, l=20, r=20),
    )
    return fig

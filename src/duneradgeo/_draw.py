"""Shared box geometry helpers used by both the matplotlib and Plotly backends."""
from __future__ import annotations

import numpy as np

# 8 unit-cube corners: bottom face (z-) then top face (z+), CCW from -x,-y
BOX_CORNERS = np.array([
    [-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1],
    [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1],
], dtype=float)

# 6 faces, each a list of 4 corner indices
BOX_FACES = [
    [0, 1, 2, 3],   # z-
    [4, 5, 6, 7],   # z+
    [0, 1, 5, 4],   # y-
    [2, 3, 7, 6],   # y+
    [0, 3, 7, 4],   # x-
    [1, 2, 6, 5],   # x+
]

# 12 edges
BOX_EDGES = [
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7),
]


def box_corners_world(center: np.ndarray, extents: np.ndarray,
                      rot: np.ndarray) -> np.ndarray:
    """(8,3) world-frame corners of a box."""
    local = BOX_CORNERS * (extents / 2)
    return (rot @ local.T).T + center


def project_polygon(corners_world: np.ndarray, ax_h: int, ax_v: int) -> np.ndarray:
    """Project 8 corners onto a 2D plane; return CCW convex-hull vertices."""
    pts = corners_world[:, [ax_h, ax_v]]
    c   = pts.mean(axis=0)
    order = np.argsort(np.arctan2(pts[:, 1] - c[1], pts[:, 0] - c[0]))
    return pts[order]

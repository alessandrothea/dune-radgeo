"""Result type, measurement registry and box helpers shared by measurements."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .._draw import box_corners_world
from ..ir import IRCollection, VolumeEntry

X  = 0          # drift axis
YZ = (1, 2)     # plane transverse to the drift (anodes, cathode)


# ---------------------------------------------------------------------------
# Result + registry
# ---------------------------------------------------------------------------

@dataclass
class Result:
    """Outcome of one measurement."""
    name:    str
    title:   str
    summary: list[str]  = field(default_factory=list)   # human-readable lines
    values:  dict       = field(default_factory=dict)   # headline numbers, None if n/a
    columns: list[str]  = field(default_factory=list)   # per-row table (--detail)
    rows:    list[dict] = field(default_factory=list)


@dataclass
class Measurement:
    name: str
    help: str
    func: Callable[[IRCollection], Result]


MEASUREMENTS: dict[str, Measurement] = {}


def measurement(name: str, help: str):
    """Register ``func(ir) -> Result`` as measurement `name`."""
    def deco(func):
        if name in MEASUREMENTS:
            raise ValueError(f"duplicate measurement '{name}'")
        MEASUREMENTS[name] = Measurement(name, help, func)
        return func
    return deco


def run(ir: IRCollection, names: list[str] | None = None) -> list[Result]:
    """Run the named measurements (default: all, in registration order)."""
    return [MEASUREMENTS[n].func(ir) for n in (names or list(MEASUREMENTS))]


# ---------------------------------------------------------------------------
# Helpers shared by measurements
# ---------------------------------------------------------------------------

def aabb(v: VolumeEntry) -> tuple[np.ndarray, np.ndarray]:
    """World-frame (lo, hi) corners of a volume's box."""
    corners = box_corners_world(v.pos, v.ext, v.rot)
    return corners.min(axis=0), corners.max(axis=0)


def tags(vols: list[VolumeEntry]) -> list[str]:
    """'name#k' labels: meta.copy for multi-part solids, else per-name order."""
    seen: dict[str, int] = {}
    out = []
    for v in vols:
        k = v.meta.get("copy")
        if k is None:
            k = seen[v.name] = seen.get(v.name, -1) + 1
        out.append(f"{v.name}#{k}")
    return out


def overlap_area(lo1, hi1, lo2, hi2, axes: tuple[int, int]) -> float:
    """Overlap area of two boxes projected on `axes`; 0 if disjoint."""
    d = [min(hi1[a], hi2[a]) - max(lo1[a], lo2[a]) for a in axes]
    return d[0] * d[1] if d[0] > 0 and d[1] > 0 else 0.0


def struct_box(ir: IRCollection, pattern: str) -> tuple | None:
    """World-frame (lo, hi) enclosing every struct volume whose name matches
    `pattern` (re.search), or None if there is none (e.g. r"^volCryostat$")."""
    boxes = [aabb(v) for v in ir.by_kind("struct") if re.search(pattern, v.name)]
    if not boxes:
        return None
    return (np.min([b[0] for b in boxes], axis=0),
            np.max([b[1] for b in boxes], axis=0))


def box_distances(lo, hi, los: np.ndarray, his: np.ndarray) -> np.ndarray:
    """Shortest distance from box (lo, hi) to each box in (los, his), shape
    (N, 3); 0 where they touch or overlap."""
    d = np.maximum(0.0, np.maximum(los - hi, lo - his))
    return np.sqrt((d ** 2).sum(axis=1))

"""pyg4ometry GDML backend: read the GDML with ``pyg4ometry.gdml.Reader``.

pyg4ometry evaluates the GDML itself (constants, variables, expressions,
units, loops, scales, assemblies of ``<physvol>``), and the result is copied
into the plain tables of :class:`GeometryBase`, so the search, placement
paths and box decomposition are exactly those of the other backends; only the
parsing differs.  Nothing refers to pyg4ometry once the constructor returns.

Optional: ``pip install "duneradgeo[pyg4ometry]"``.  pyg4ometry imports
``vtk`` (about 560 MB) when it is imported, although neither the reader nor
the mesher needs it; when ``vtk`` is not installed a stand-in module is
supplied for the duration of that import (see :func:`_import_pyg4ometry`).

Solids: box, tube and torus take their bounding box from their parameters
(the same code as the built-in parser), booleans keep their operands.  Any
other solid type gets the bounding box of the mesh pyg4ometry builds for it,
so ``volume_gen`` can sample it; such solids yield no box decomposition.
"""
from __future__ import annotations

import re
import sys
import types
import xml.etree.ElementTree as ET

import numpy as np

from ._gdml import GeometryBase, torus_bbox, tube_bbox

_BOOLEANS = {"Union": "union", "Subtraction": "subtraction",
             "Intersection": "intersection"}

_pg = None


class Pyg4ometryUnavailable(RuntimeError):
    """pyg4ometry cannot be imported."""


def _import_pyg4ometry(feature: str = "the 'pyg4ometry' GDML backend"):
    """The pyg4ometry module, imported once.

    pyg4ometry's ``visualisation`` package imports ``vtk`` at module level.
    Without VTK a placeholder module is put in ``sys.modules`` while
    pyg4ometry is imported and removed again afterwards: the reader, the
    writer and the meshes never call it, and a later real use of VTK fails
    with its own ImportError."""
    global _pg
    if _pg is not None:
        return _pg
    stub = None
    try:
        import vtk                                           # noqa: F401
    except ImportError:
        stub = types.ModuleType("vtk")
        stub.__getattr__ = lambda name: (
            (_ for _ in ()).throw(AttributeError(name)) if name.startswith("__")
            else type(name, (), {}))
        sys.modules["vtk"] = stub
    try:
        import pyg4ometry
    except ImportError as exc:
        raise Pyg4ometryUnavailable(
            f"{feature} needs pyg4ometry, which cannot be imported ({exc}): "
            "install it with 'pip install \"duneradgeo[pyg4ometry]\"'") from exc
    finally:
        if stub is not None and sys.modules.get("vtk") is stub:
            del sys.modules["vtk"]
    _pg = pyg4ometry
    return _pg


def _matrix(pg, rotation: list[float]) -> np.ndarray:
    """3×3 matrix of a pyg4ometry rotation (tait-bryan x-y-z angles [rad]):
    the GDML rotation matrix, i.e. the one the other backends call ``R``."""
    return np.array(pg.transformation.tbxyz2matrix(rotation), dtype=float)


class Pyg4GDMLGeometry(GeometryBase):
    """GDML geometry read by pyg4ometry, same API as
    :class:`~duneradgeo._gdml.GDMLGeometry`.  Lengths are converted to cm."""

    backend = "pyg4ometry"

    def __init__(self, path: str, verbose: bool = False):
        super().__init__()
        pg = _import_pyg4ometry()
        self._pg = pg
        self.skipped: dict[str, int] = {}       # unsupported element → count
        try:
            reg = pg.gdml.Reader(path, skipMaterials=True).getRegistry()
            world = reg.getWorldVolume()
        except Exception as exc:               # pyg4ometry raises many types
            raise RuntimeError(
                f"pyg4ometry cannot read {path}: {type(exc).__name__}: {exc}") from exc
        if world is None:
            raise RuntimeError(f"{path} has no <world> in its <setup>")

        self._pvnames = self._scan_physvol_names(path)
        for name, solid in reg.solidDict.items():
            self._add_solid(name, solid, verbose)
        for name, lv in reg.logicalVolumeDict.items():
            self._add_volume(name, lv)
        self._world = world.name

        if self.skipped:
            what = ", ".join(f"{n} {k}" for k, n in sorted(self.skipped.items()))
            print(f"[pyg4ometry] skipped (not supported): {what}", file=sys.stderr)

    # -- table filling -------------------------------------------------------

    def _unit(self, name: str) -> float:
        """Factor of a GDML unit to cm (pyg4ometry's are relative to mm)."""
        return float(self._pg.gdml.Units.unit(name)) / 10.0

    def _add_solid(self, name: str, s, verbose: bool) -> None:
        kind = type(s).__name__
        ev = s.evaluateParameter
        if kind == "Box":
            u = self._unit(s.lunit)
            ext = np.array([ev(s.pX), ev(s.pY), ev(s.pZ)], dtype=float) * u
            self._boxes[name] = ext
            self._prims[name] = (np.zeros(3), ext)
        elif kind == "Tubs":
            u, au = self._unit(s.lunit), float(self._pg.gdml.Units.unit(s.aunit))
            self._prims[name] = tube_bbox(
                ev(s.pRMin) * u, ev(s.pRMax) * u, ev(s.pDz) * u / 2,
                ev(s.pSPhi) * au, ev(s.pDPhi) * au)
        elif kind == "Torus":
            u, au = self._unit(s.lunit), float(self._pg.gdml.Units.unit(s.aunit))
            self._prims[name] = torus_bbox(
                ev(s.pRtor) * u, ev(s.pRmax) * u,
                ev(s.pSPhi) * au, ev(s.pDPhi) * au)
        elif kind in _BOOLEANS:
            # <firstposition>/<firstrotation> are not kept by pyg4ometry's
            # booleans either: the first operand sits at the origin
            pos = np.array(s.tra2[1].eval(), dtype=float) / 10.0
            rot = _matrix(self._pg, s.tra2[0].eval())
            self._bools[name] = (_BOOLEANS[kind], s.obj1.name, s.obj2.name,
                                 pos, rot, np.zeros(3), np.eye(3))
        else:
            self._add_meshed(name, s, kind, verbose)

    def _add_meshed(self, name: str, s, kind: str, verbose: bool) -> None:
        """A solid of another type: its bounding box is that of its mesh."""
        try:
            verts = np.array(s.mesh().toVerticesAndPolygons()[0], dtype=float) / 10.0
            lo, hi = verts.min(axis=0), verts.max(axis=0)
        except Exception as exc:
            self.skipped[kind] = self.skipped.get(kind, 0) + 1
            if verbose:
                print(f"[pyg4ometry] {kind} {name}: no bounding box ({exc})",
                      file=sys.stderr)
            return
        self._shape_bbox[name] = ((lo + hi) / 2, hi - lo)

    def _add_volume(self, name: str, lv) -> None:
        daughters = []
        names = self._pvnames.get(name)
        phys = [pv for pv in lv.daughterVolumes if type(pv).__name__ == "PhysicalVolume"]
        if names is not None and len(names) != len(phys):
            names = None                       # not the same list: use the fallback
        k = -1
        for pv in lv.daughterVolumes:
            sub = getattr(pv, "logicalVolume", None)
            if sub is None or type(pv).__name__ != "PhysicalVolume":
                kind = type(pv).__name__
                self.skipped[kind] = self.skipped.get(kind, 0) + 1
                continue
            k += 1
            pos = np.array(pv.position.eval(), dtype=float) / 10.0
            rot = _matrix(self._pg, pv.rotation.eval())
            # the GDML rotation is the inverse of the placement (as in the
            # built-in parser); pyg4ometry's scale is not applied
            daughters.append((sub.name, pos, rot.T,
                              names[k] if names is not None
                              else self._guess_pv_name(pv, sub.name)))
        solid = lv.solid
        self._vols[name] = {"solid": solid.name if solid is not None else None,
                            "daughters": daughters}

    @staticmethod
    def _scan_physvol_names(path: str) -> dict:
        """logical volume → the ``name`` of each of its ``<physvol>`` in order
        (None when the GDML gives none), read from the XML: pyg4ometry names
        an unnamed physvol itself.  {} when the file cannot be scanned (e.g.
        ``<!ENTITY>`` includes)."""
        try:
            root = ET.parse(path).getroot()
        except ET.ParseError:
            return {}
        out = {}
        for structure in root.iter("structure"):
            for vol in structure.findall("volume"):
                out[vol.get("name")] = [pv.get("name") for pv in vol.findall("physvol")]
        return out

    @staticmethod
    def _guess_pv_name(pv, vol_name: str) -> str | None:
        """Fallback when the XML was not scanned: pyg4ometry names an unnamed
        physvol ``<volume>_PV`` (``<volume>_<n>_PV`` for the n-th repetition),
        so a physvol named like that in the GDML is taken as unnamed."""
        if re.fullmatch(re.escape(vol_name) + r"(_\d+)?_PV", pv.name):
            return None
        return pv.name

"""ROOT GDML backend: read the GDML with ``TGeoManager::Import``.

ROOT parses the file (units, expressions, constants, rotations) and the
result is copied into the plain tables of :class:`GeometryBase`, so the
search, placement paths and box decomposition are exactly those of the
pure-Python backend; only the parsing differs.  Nothing refers to ROOT once
the constructor returns.

Needs a ROOT built with the ``geom`` and ``gdml`` features: the dunesw/UPS
ROOT in the container has them; the ``root`` pip wheels (0.1a12) do not yet.
"""
from __future__ import annotations

import ctypes
import os
import re
import sys
import tempfile
from contextlib import contextmanager

import numpy as np

from ._gdml import GeometryBase
from ._root import ROOTUnavailable, import_root

# TGeo classes taken as primitives, i.e. what the pure-Python backend handles
# (GDML box, tube, torus).  Their local bounding box is ROOT's own.
_PRIM_CLASSES = {"TGeoBBox", "TGeoTube", "TGeoTubeSeg", "TGeoTorus"}

# TGeoBoolNode::EGeoBoolType
_BOOL_OPS = {0: "union", 1: "intersection", 2: "subtraction"}

# Import chatter that carries no information about the geometry
_QUIET = re.compile(r"Unsupported GDML Tag Used|^Info in <")


def _import_root():
    ROOT = import_root("the 'root' GDML backend")
    if not hasattr(ROOT, "TGeoManager") or ROOT.gSystem.Load("libGdml") < 0:
        raise ROOTUnavailable(
            f"ROOT {ROOT.gROOT.GetVersion()} has no TGeo/GDML support "
            "(libGeom/libGdml missing); use --gdml-backend python")
    return ROOT


@contextmanager
def _captured_output(ROOT, verbose: bool):
    """Collect C/C++ stdout+stderr of the block; replay it, filtered, to stderr.

    TGDMLParse writes "Unsupported GDML Tag" for every <color>/<extension>
    straight to std::cout, so the file descriptors are redirected."""
    libc = ctypes.CDLL(None)
    sys.stdout.flush()
    sys.stderr.flush()
    saved = [os.dup(1), os.dup(2)]
    with tempfile.TemporaryFile(mode="w+b") as tmp:
        os.dup2(tmp.fileno(), 1)
        os.dup2(tmp.fileno(), 2)
        try:
            yield
        finally:
            ROOT.std.cout.flush()
            ROOT.std.cerr.flush()
            libc.fflush(None)
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            for fd in saved:
                os.close(fd)
            tmp.seek(0)
            text = tmp.read().decode(errors="replace")
        n_quiet = 0
        for line in text.splitlines():
            if not line.strip():
                continue
            if not verbose and _QUIET.search(line):
                n_quiet += 1
                continue
            print(f"[root] {line}", file=sys.stderr)
        if n_quiet:
            print(f"[root] ({n_quiet} info / unsupported-tag lines hidden; "
                  "--verbose shows them)", file=sys.stderr)


def _matrix(m) -> tuple[np.ndarray, np.ndarray]:
    """(translation, 3×3 rotation) of a TGeoMatrix: x_mother = R @ x + t."""
    t = m.GetTranslation()
    r = m.GetRotationMatrix()
    return (np.array([t[0], t[1], t[2]], dtype=float),
            np.array([r[i] for i in range(9)], dtype=float).reshape(3, 3))


class RootGDMLGeometry(GeometryBase):
    """GDML geometry read by ROOT (``TGeoManager::Import``), same API as
    :class:`~duneradgeo._gdml.GDMLGeometry`."""

    backend = "root"

    def __init__(self, path: str, verbose: bool = False):
        super().__init__()
        ROOT = _import_root()
        ROOT.TGeoManager.SetVerboseLevel(1 if verbose else 0)
        with _captured_output(ROOT, verbose):
            geom = ROOT.TGeoManager.Import(path)
        if not geom:
            raise RuntimeError(f"TGeoManager::Import failed for {path}")

        for shape in geom.GetListOfShapes():
            self._add_shape(shape)
        for vol in geom.GetListOfVolumes():
            self._add_volume(vol)
        self._world = geom.GetTopVolume().GetName()

    # -- table filling -------------------------------------------------------

    def _add_shape(self, shape) -> None:
        name = shape.GetName()
        cls  = shape.ClassName()
        if hasattr(shape, "GetDX"):
            # ROOT's own bounding box, for every shape class (composites
            # included): the box BaseRadioGen samples volume_gen nodes in
            o = shape.GetOrigin()
            self._shape_bbox[name] = (
                np.array([o[0], o[1], o[2]], dtype=float),
                2 * np.array([shape.GetDX(), shape.GetDY(), shape.GetDZ()]))
        if cls in _PRIM_CLASSES:
            o   = shape.GetOrigin()
            ext = 2 * np.array([shape.GetDX(), shape.GetDY(), shape.GetDZ()])
            self._prims[name] = (np.array([o[0], o[1], o[2]], dtype=float), ext)
            if cls == "TGeoBBox":
                self._boxes[name] = ext
        elif shape.IsComposite():
            node   = shape.GetBoolNode()
            first  = node.GetLeftShape().GetName()
            second = node.GetRightShape().GetName()
            fpos, frot = _matrix(node.GetLeftMatrix())
            pos,  rot  = _matrix(node.GetRightMatrix())
            op = _BOOL_OPS.get(int(node.GetBooleanOperator()), "unknown")
            self._bools[name] = (op, first, second, pos, rot, fpos, frot)

    def _add_volume(self, vol) -> None:
        shape = vol.GetShape()
        daughters = []
        for i in range(vol.GetNdaughters()):
            node  = vol.GetNode(i)
            dname = node.GetVolume().GetName()
            pos, rot = _matrix(node.GetMatrix())
            daughters.append((dname, pos, rot, self._pv_name(node, dname)))
        mat = vol.GetMaterial()
        if mat:
            self._vol_material[vol.GetName()] = mat.GetName()
            # TGeoMaterial density is in g/cm3
            self._densities[mat.GetName()] = float(mat.GetDensity())
        self._vols[vol.GetName()] = {
            "solid": shape.GetName() if shape else None,
            "daughters": daughters,
        }

    @staticmethod
    def _pv_name(node, vol_name: str) -> str | None:
        """The GDML physvol name, or None when the physvol was unnamed:
        TGDMLParse then names the node ``<volume>_<copy number>``."""
        name = node.GetName()
        return None if name == f"{vol_name}_{node.GetNumber()}" else name

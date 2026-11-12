"""GDML geometry: world-frame positions and extents of named volumes.

Two backends share the queries of :class:`GeometryBase`:

* ``python`` — :class:`GDMLGeometry`, a pure-Python parser (numpy only);
* ``root``   — :class:`~duneradgeo._gdml_root.RootGDMLGeometry`, which reads
  the GDML with ROOT's ``TGeoManager::Import`` (needs ROOT with TGeo/GDML).

:func:`load_geometry` picks one by name.
"""
from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np

from ._draw import aabb


def _rot_matrix_rad(rx: float, ry: float, rz: float) -> np.ndarray:
    """3×3 rotation matrix from GDML Euler angles [rad]: R = Rz @ Ry @ Rx.

    R maps a boolean operand's frame into the boolean solid's frame.  A
    physvol is placed with the inverse, R.T (see :meth:`_parse_physvol`)."""
    cx, sx = np.cos(rx), np.sin(rx)
    cy, sy = np.cos(ry), np.sin(ry)
    cz, sz = np.cos(rz), np.sin(rz)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx,  cx]])
    Ry = np.array([[cy, 0, sy], [0,  1,  0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz,  0], [0,   0,  1]])
    return Rz @ Ry @ Rx


def _arc_points(r: float, phi0: float, dphi: float) -> np.ndarray:
    """(N, 2) points on a circle arc of radius r that include its x/y extremes."""
    if dphi >= 2 * np.pi - 1e-9:
        phis = np.arange(4) * np.pi / 2
    else:
        quads = np.arange(np.ceil(phi0 / (np.pi / 2)),
                          np.floor((phi0 + dphi) / (np.pi / 2)) + 1) * np.pi / 2
        phis = np.concatenate([[phi0, phi0 + dphi], quads])
    return r * np.column_stack([np.cos(phis), np.sin(phis)])

_MASS_UNIT = {"g": 1.0, "mg": 1e-3, "kg": 1e3, "ug": 1e-6}               # → g
_VOLUME_UNIT = {"cm3": 1.0, "mm3": 1e-3, "m3": 1e6, "dm3": 1e3, "l": 1e3}  # → cm3


def density_factor(unit: str) -> float | None:
    """Factor of a GDML density unit (``g/cm3``, ``kg/m3``, …) to g/cm3, or
    None when the unit is not understood.  GDML's default is g/cm3."""
    unit = unit.strip().lower().replace(" ", "")
    if not unit:
        return 1.0
    mass, sep, vol = unit.partition("/")
    if not sep or mass not in _MASS_UNIT or vol not in _VOLUME_UNIT:
        return None
    return _MASS_UNIT[mass] / _VOLUME_UNIT[vol]


def _placed_aabb(bbox: tuple, pos: np.ndarray, rot: np.ndarray) -> tuple:
    """(lo, hi) of box `bbox` = (centre, ext) mapped by x -> rot @ x + pos:
    the envelope of its 8 corners."""
    centre, ext = bbox
    return aabb(pos + rot @ centre, ext, rot)


class GeometryBase:
    """
    Backend-independent queries on a parsed geometry.

    A backend fills these tables (lengths in cm) and inherits the search:

    ``_vols``
        logical volume → ``{"solid": name | None, "daughters": [(vol, pos,
        rot, pvname), ...]}`` in placement order, where ``rot``/``pos`` map
        the daughter frame into the mother frame (``x_m = rot @ x_d + pos``)
        and ``pvname`` is the physvol name, or None when the GDML has none.
    ``_prims``
        primitive solid (box, tube, torus) → ``(centre, ext)``, its local
        axis-aligned bounding box.
    ``_bools``
        boolean solid → ``(op, first, second, pos, rot, fpos, frot)``, ``op``
        one of union / subtraction / intersection; ``pos, rot`` place the
        second operand and ``fpos, frot`` the first in the boolean's frame.
    ``_boxes``
        box solid → full extents ``(dx, dy, dz)``.
    ``_shape_bbox``
        solid → ``(centre, ext)``, its local bounding box as the reader
        computed it (the ROOT backend: ROOT's own, for every shape).  Solids
        not listed get theirs from :meth:`solid_bbox`.
    ``_vol_material``
        logical volume → name of its material (volumes without one are absent).
    ``_densities``
        material → density [g/cm3] (materials without a density, e.g. NIST
        ones that the GDML does not define, are absent).
    ``_world``
        name of the world volume.
    """

    backend = ""

    def __init__(self) -> None:
        self._boxes: dict[str, np.ndarray] = {}
        self._prims: dict[str, tuple]      = {}
        self._bools: dict[str, tuple]      = {}
        self._vols:  dict[str, dict]       = {}
        self._world: str | None            = None
        self._vol_material: dict[str, str]  = {}
        self._densities: dict[str, float]   = {}
        self._shape_bbox: dict[str, tuple] = {}

    # -- tree search -------------------------------------------------------

    def find_positions(
        self,
        vol_patterns: list[str],
        world_vol: str | None = None,
        verbose: bool = False,
        with_path: bool = False,
    ) -> list[tuple]:
        """
        Return (name, world_pos, world_rot) for every placement of volumes
        whose name matches any pattern in vol_patterns.

        A backward BFS from the matching volumes collects every volume that
        can contain one; a DFS from the world volume then composes the
        placements, pruning subtrees outside that set.

        With ``with_path`` each tuple also carries the placement path below
        the world volume, ``volA:0/volB:3/name:1``, where ``:k`` counts the
        placements of that volume within its parent in GDML order.  It does
        not depend on the patterns, so it identifies a placement across
        extractions of the same (or a locally edited) GDML.  Two more items
        follow: the physvol name of the placement and the physvol names
        along the path, ``pvA/pvB/pv`` (None and "" when the GDML leaves
        physvols unnamed, as the refactored geometries do; unnamed steps of
        a partly named path are empty).
        """
        patterns = [re.compile(p) for p in vol_patterns]
        world    = world_vol or self._world
        if world is None:
            raise ValueError("No world volume; pass world_vol explicitly")

        target_vols = {
            name for name in self._vols
            if any(p.search(name) for p in patterns)
        }
        if verbose:
            print(f"  [dbg] {len(target_vols)} target vols: "
                  f"{sorted(target_vols)[:6]}", file=sys.stderr)
        if not target_vols:
            return []

        # BFS backward to find all ancestors of targets
        parent_of: dict[str, set[str]] = defaultdict(set)
        for pname, pdata in self._vols.items():
            for (dname, *_) in pdata["daughters"]:
                parent_of[dname].add(pname)

        can_contain: set[str] = set(target_vols)
        frontier: set[str]   = set(target_vols)
        while frontier:
            nxt: set[str] = set()
            for v in frontier:
                for p in parent_of[v]:
                    if p not in can_contain:
                        can_contain.add(p)
                        nxt.add(p)
            frontier = nxt

        if verbose:
            print(f"  [dbg] {len(can_contain)} vols in ancestor set", file=sys.stderr)

        results: list[tuple] = []

        def walk(vol_name, w_pos, w_rot, ancestry, path, pvname, pvpath):
            if vol_name in ancestry:
                return
            if vol_name in target_vols:
                hit = (vol_name, w_pos.copy(), w_rot.copy())
                pvp = pvpath if pvpath.strip("/") else ""
                results.append(hit + (path, pvname, pvp) if with_path else hit)
                # keep recursing: a target may contain other targets
            if vol_name not in self._vols or vol_name not in can_contain:
                return
            new_anc = ancestry | {vol_name}
            seen: dict[str, int] = defaultdict(int)
            for dname, local_pos, local_rot, dpv in self._vols[vol_name]["daughters"]:
                # count every sibling, pruned or not, so k is pattern-independent
                k = seen[dname]
                seen[dname] += 1
                if dname not in can_contain:
                    continue
                step = f"{dname}:{k}"
                walk(dname, w_pos + w_rot @ local_pos, w_rot @ local_rot, new_anc,
                     f"{path}/{step}" if path else step, dpv,
                     f"{pvpath}/{dpv or ''}" if path else (dpv or ""))

        walk(world, np.zeros(3), np.eye(3), frozenset(), "", None, "")
        return results

    # -- utilities ---------------------------------------------------------

    @property
    def volume_names(self) -> set[str]:
        return set(self._vols.keys())

    def solid_parts(self, solid: str, _seen: frozenset = frozenset()) -> list[tuple]:
        """
        Primitive pieces of `solid` as (name, centre, rot, ext) in the solid's
        own frame, each piece being the local bounding box of a box, tube or
        torus.  Both operands of a union are kept; subtractions and
        intersections keep only their first operand (the outer envelope).
        Unsupported primitives yield no pieces.
        """
        if solid in self._prims:
            centre, ext = self._prims[solid]
            return [(solid, centre, np.eye(3), ext)]
        node = self._bools.get(solid)
        if node is None or solid in _seen:
            return []
        op, first, second, pos, rot, fpos, frot = node
        seen  = _seen | {solid}
        parts = [(n, fpos + frot @ c, frot @ r, e)
                 for n, c, r, e in self.solid_parts(first, seen)]
        if op == "union":
            parts += [(n, pos + rot @ c, rot @ r, e)
                      for n, c, r, e in self.solid_parts(second, seen)]
        return parts

    def volume_material(self, vol_name: str) -> str | None:
        """Name of the material of logical volume `vol_name`, or None."""
        return self._vol_material.get(vol_name)

    def volume_density(self, vol_name: str) -> float | None:
        """Density [g/cm3] of the material of logical volume `vol_name`, or
        None when the volume has no material or the material no density."""
        return self._densities.get(self._vol_material.get(vol_name))

    def volume_parts(self, vol_name: str) -> list[tuple]:
        """:meth:`solid_parts` of the solid of logical volume `vol_name`."""
        vol = self._vols.get(vol_name)
        if vol is None or not vol.get("solid"):
            return []
        return self.solid_parts(vol["solid"])

    def solid_bbox(self, solid: str, _seen: frozenset = frozenset()) -> tuple | None:
        """
        ``(centre, ext)``: the local axis-aligned bounding box of `solid`,
        the box ROOT gives the shape (TGeoBBox origin and DX/DY/DZ), which
        BaseRadioGen samples a volume_gen node in.  A primitive's own box; a
        boolean's as TGeoUnion / TGeoSubtraction / TGeoIntersection::
        ComputeBBox build it from the operands' boxes, each placed by its 8
        corners: the envelope of both operands (union), of the first one
        (subtraction), or the overlap of both (intersection).  None when a
        solid involved is not supported (cycle-safe).
        """
        bb = self._shape_bbox.get(solid)
        if bb is None:
            bb = self._prims.get(solid)
        if bb is not None:
            return bb
        node = self._bools.get(solid)
        if node is None or solid in _seen:
            return None
        op, first, second, pos, rot, fpos, frot = node
        seen = _seen | {solid}
        a = self.solid_bbox(first, seen)
        b = None if op == "subtraction" else self.solid_bbox(second, seen)
        if a is None or (b is None and op != "subtraction"):
            return None
        lo, hi = _placed_aabb(a, fpos, frot)
        if op == "union":
            blo, bhi = _placed_aabb(b, pos, rot)
            lo, hi = np.minimum(lo, blo), np.maximum(hi, bhi)
        elif op == "intersection":
            blo, bhi = _placed_aabb(b, pos, rot)
            lo, hi = np.maximum(lo, blo), np.minimum(hi, bhi)
            hi = np.maximum(hi, lo)              # disjoint operands: empty box
        elif op != "subtraction":
            return None
        return (lo + hi) / 2, hi - lo

    def volume_bbox(self, vol_name: str) -> tuple | None:
        """:meth:`solid_bbox` of the solid of logical volume `vol_name`."""
        vol = self._vols.get(vol_name)
        if vol is None or not vol.get("solid"):
            return None
        return self.solid_bbox(vol["solid"])


class GDMLGeometry(GeometryBase):
    """
    GDML geometry read in pure Python (``xml.etree``): fills the
    :class:`GeometryBase` tables from the define, solids, structure and
    setup sections; the search is :meth:`GeometryBase.find_positions`.
    The full FD 1×8×14 geometry loads and searches in about 1 s.
    """

    backend = "python"

    def __init__(self, path: str):
        super().__init__()
        tree = ET.parse(path)
        root = tree.getroot()
        for elem in root.iter():           # strip namespace prefix if present
            if "}" in elem.tag:
                elem.tag = elem.tag.split("}", 1)[1]

        self._defs:  dict[str, np.ndarray] = {}   # named positions [cm]
        self._rots:  dict[str, np.ndarray] = {}   # named rotations (3×3)

        self._parse_defines(root)
        self._parse_materials(root)
        self._parse_solids(root)
        self._parse_structure(root)
        self._parse_setup(root)

    # -- unit / expression helpers ----------------------------------------

    @staticmethod
    def _ufact(unit: str) -> float:
        return {"cm": 1.0, "mm": 0.1, "m": 100.0, "um": 1e-4,
                "centimeter": 1.0, "millimeter": 0.1, "meter": 100.0,
                "micrometer": 1e-4,
                "in": 2.54, "ft": 30.48}.get(unit.lower().strip(), 1.0)

    @staticmethod
    def _aufact(unit: str) -> float:
        """Angle unit → radians.  GDML writers spell these differently
        (``deg`` in the refactored geometries, ``radian``/``degree`` in ggd)."""
        return {"rad": 1.0, "radian": 1.0, "mrad": 1e-3, "milliradian": 1e-3,
                "deg": np.pi / 180.0, "degree": np.pi / 180.0,
                }.get(unit.lower().strip(), 1.0)

    @staticmethod
    def _fval(s: str) -> float:
        """Parse a GDML attribute that may be a simple arithmetic expression."""
        s = s.strip()
        if not s:
            return 0.0
        try:
            return float(s)
        except ValueError:
            clean = re.sub(r"[^0-9.+\-*/() eE]", "", s)
            if not clean:
                return 0.0
            return float(eval(clean, {"__builtins__": {}}))  # noqa: S307

    # -- section parsers --------------------------------------------------

    def _parse_defines(self, root):
        defines = root.find("define")
        if defines is None:
            return
        for child in defines:
            name = child.get("name", "")
            if child.tag == "position":
                u = self._ufact(child.get("unit", "mm"))
                self._defs[name] = np.array([
                    self._fval(child.get("x", "0")) * u,
                    self._fval(child.get("y", "0")) * u,
                    self._fval(child.get("z", "0")) * u,
                ])
            elif child.tag == "rotation":
                u = self._aufact(child.get("unit", "rad"))
                self._rots[name] = _rot_matrix_rad(
                    self._fval(child.get("x", "0")) * u,
                    self._fval(child.get("y", "0")) * u,
                    self._fval(child.get("z", "0")) * u,
                )

    def _parse_materials(self, root):
        """Material densities, in g/cm3; a material without a readable <D>
        has none."""
        materials = root.find("materials")
        if materials is None:
            return
        for mat in materials.findall("material"):
            d = mat.find("D")
            if d is None or not mat.get("name"):
                continue
            factor = density_factor(d.get("unit", ""))
            if factor is None:
                print(f"  [warn] material {mat.get('name')}: density unit "
                      f"{d.get('unit')!r} not understood, density skipped",
                      file=sys.stderr)
                continue
            try:
                self._densities[mat.get("name")] = self._fval(d.get("value", "")) * factor
            except (SyntaxError, ValueError, TypeError):
                print(f"  [warn] material {mat.get('name')}: density "
                      f"{d.get('value')!r} not understood, density skipped",
                      file=sys.stderr)

    def _parse_solids(self, root):
        solids = root.find("solids")
        if solids is None:
            return
        for child in solids:
            name = child.get("name", "")
            u    = self._ufact(child.get("lunit", "mm"))
            au   = self._aufact(child.get("aunit", "rad"))
            if child.tag == "box":
                self._boxes[name] = np.array([
                    self._fval(child.get("x", "0")) * u,
                    self._fval(child.get("y", "0")) * u,
                    self._fval(child.get("z", "0")) * u,
                ])
                self._prims[name] = (np.zeros(3), self._boxes[name])
            elif child.tag == "tube":
                phi0 = self._fval(child.get("startphi", "0")) * au
                dphi = self._fval(child.get("deltaphi", str(2 * np.pi / au))) * au
                pts  = np.vstack([
                    _arc_points(self._fval(child.get(k, "0")) * u, phi0, dphi)
                    for k in ("rmin", "rmax")])
                hz   = self._fval(child.get("z", "0")) * u / 2
                self._prims[name] = self._aabb(pts, -hz, hz)
            elif child.tag == "torus":
                phi0 = self._fval(child.get("startphi", "0")) * au
                dphi = self._fval(child.get("deltaphi", str(2 * np.pi / au))) * au
                rmax = self._fval(child.get("rmax", "0")) * u
                rtor = self._fval(child.get("rtor", "0")) * u
                if rtor >= rmax:
                    # the cross-section at angle phi spans radii rtor ± rmax
                    # along (cos phi, sin phi): in x–y the torus covers the
                    # annular sector of a tube with those radii
                    pts = np.vstack([_arc_points(rtor - rmax, phi0, dphi),
                                     _arc_points(rtor + rmax, phi0, dphi)])
                else:
                    # self-intersecting torus: pad the centre arc by rmax
                    arc = _arc_points(rtor, phi0, dphi)
                    pts = np.vstack([arc - rmax, arc + rmax])
                self._prims[name] = self._aabb(pts, -rmax, rmax)
        for child in solids:
            if child.tag in ("subtraction", "union", "intersection"):
                first  = child.find("first")
                second = child.find("second")
                if first is not None and second is not None:
                    pos, rot = self._parse_transform(child)
                    # the first operand is taken as placed at the origin
                    # (<firstposition>/<firstrotation> are not read)
                    self._bools[child.get("name", "")] = (
                        child.tag, first.get("ref", ""), second.get("ref", ""), pos, rot,
                        np.zeros(3), np.eye(3))

    @staticmethod
    def _aabb(xy: np.ndarray, z0: float, z1: float) -> tuple:
        lo = np.array([xy[:, 0].min(), xy[:, 1].min(), z0])
        hi = np.array([xy[:, 0].max(), xy[:, 1].max(), z1])
        return (lo + hi) / 2, hi - lo

    def _parse_transform(self, elem) -> tuple:
        """(pos, rot) from the position/rotation children of a physvol or boolean."""
        pos = np.zeros(3)
        rot = np.eye(3)
        for child in elem:
            tag = child.tag
            if tag == "position":
                u = self._ufact(child.get("unit", "mm"))
                pos = np.array([
                    self._fval(child.get("x", "0")) * u,
                    self._fval(child.get("y", "0")) * u,
                    self._fval(child.get("z", "0")) * u,
                ])
            elif tag == "positionref":
                pos = self._defs.get(child.get("ref", ""), np.zeros(3))
            elif tag == "rotation":
                u = self._aufact(child.get("unit", "rad"))
                rot = _rot_matrix_rad(
                    self._fval(child.get("x", "0")) * u,
                    self._fval(child.get("y", "0")) * u,
                    self._fval(child.get("z", "0")) * u,
                )
            elif tag == "rotationref":
                rot = self._rots.get(child.get("ref", ""), np.eye(3))
        return pos, rot

    def _parse_physvol(self, pv) -> tuple:
        volref = None
        for child in pv:
            if child.tag == "volumeref":
                volref = child.get("ref")
        pos, rot = self._parse_transform(pv)
        # Geant4 places daughters with the inverse of the GDML rotation (and
        # boolean operands without it).  Checked on the VD geometries: only
        # this combination centres the field-shaper rings on the TPC, joins
        # their corner tori to the tube ends and lays the cathode X-ARAPUCAs
        # flat in the cathode plane.
        return volref, pos, rot.T, pv.get("name")

    def _parse_structure(self, root):
        structure = root.find("structure")
        if structure is None:
            return
        for vol in structure:
            if vol.tag != "volume":
                continue
            name      = vol.get("name", "")
            solid     = None
            daughters = []
            for child in vol:
                if child.tag == "solidref":
                    solid = child.get("ref")
                elif child.tag == "materialref":
                    self._vol_material[name] = child.get("ref")
                elif child.tag == "physvol":
                    vref, pos, rot, pvname = self._parse_physvol(child)
                    if vref:
                        daughters.append((vref, pos, rot, pvname))
            self._vols[name] = {"solid": solid, "daughters": daughters}

    def _parse_setup(self, root):
        setup = root.find("setup")
        if setup is not None:
            world = setup.find("world")
            if world is not None:
                self._world = world.get("ref")


BACKENDS = ("python", "root")


def load_geometry(path: str, backend: str = "python",
                  verbose: bool = False) -> GeometryBase:
    """Parse the GDML at `path` with the named backend (see :data:`BACKENDS`)."""
    if backend == "python":
        return GDMLGeometry(path)
    if backend == "root":
        from ._gdml_root import RootGDMLGeometry
        return RootGDMLGeometry(path, verbose=verbose)
    raise ValueError(f"unknown GDML backend {backend!r}; choose from {BACKENDS}")

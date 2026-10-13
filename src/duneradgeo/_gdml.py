"""Pure-Python GDML parser: world-frame positions and extents of named volumes."""
from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np


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

class GDMLGeometry:
    """
    Parse a DUNE GDML file and locate world-frame positions of named volumes.

    Algorithm:
      1. BFS backward from target volumes to build `can_contain` ancestor set.
      2. DFS from the world volume, composing translations and rotations, pruning
         any subtree whose root is not in `can_contain`.

    Full FD 1×8×14 geometry runs in ~1 s with this approach.
    """

    def __init__(self, path: str):
        tree = ET.parse(path)
        root = tree.getroot()
        for elem in root.iter():           # strip namespace prefix if present
            if "}" in elem.tag:
                elem.tag = elem.tag.split("}", 1)[1]

        self._defs:  dict[str, np.ndarray] = {}   # named positions [cm]
        self._rots:  dict[str, np.ndarray] = {}   # named rotations (3×3)
        self._boxes: dict[str, np.ndarray] = {}   # solid → (dx,dy,dz) [cm]
        self._prims: dict[str, tuple]      = {}   # box/tube/torus → local AABB (centre, ext)
        self._bools: dict[str, tuple]      = {}   # boolean → (op, first, second, pos, rot)
        self._vols:  dict[str, dict]       = {}   # vol → {solid, daughters}
        self._world: str | None            = None
        self._solid_first: dict[str, str]  = {}   # boolean-solid → outer ref

        self._parse_defines(root)
        self._parse_solids(root)
        self._parse_structure(root)
        self._parse_setup(root)

    # -- unit / expression helpers ----------------------------------------

    @staticmethod
    def _ufact(unit: str) -> float:
        return {"cm": 1.0, "mm": 0.1, "m": 100.0, "um": 1e-4,
                "in": 2.54, "ft": 30.48}.get(unit.lower().strip(), 1.0)

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
                unit = child.get("unit", "rad")
                u = 1.0 if unit == "rad" else np.pi / 180.0
                self._rots[name] = _rot_matrix_rad(
                    self._fval(child.get("x", "0")) * u,
                    self._fval(child.get("y", "0")) * u,
                    self._fval(child.get("z", "0")) * u,
                )

    def _parse_solids(self, root):
        solids = root.find("solids")
        if solids is None:
            return
        for child in solids:
            name = child.get("name", "")
            u    = self._ufact(child.get("lunit", "mm"))
            au   = 1.0 if child.get("aunit", "rad") == "rad" else np.pi / 180.0
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
                arc  = _arc_points(self._fval(child.get("rtor", "0")) * u, phi0, dphi)
                # swept disk of radius rmax around the arc
                pts  = np.vstack([arc - rmax, arc + rmax])
                self._prims[name] = self._aabb(pts, -rmax, rmax)
        for child in solids:
            if child.tag in ("subtraction", "union", "intersection"):
                first  = child.find("first")
                second = child.find("second")
                if first is not None:
                    self._solid_first[child.get("name", "")] = first.get("ref", "")
                if first is not None and second is not None:
                    pos, rot = self._parse_transform(child)
                    self._bools[child.get("name", "")] = (
                        child.tag, first.get("ref", ""), second.get("ref", ""), pos, rot)

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
                unit = child.get("unit", "rad")
                u = 1.0 if unit == "rad" else np.pi / 180.0
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
        return volref, pos, rot.T

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
                elif child.tag == "physvol":
                    vref, pos, rot = self._parse_physvol(child)
                    if vref:
                        daughters.append((vref, pos, rot))
            self._vols[name] = {"solid": solid, "daughters": daughters}

    def _parse_setup(self, root):
        setup = root.find("setup")
        if setup is not None:
            world = setup.find("world")
            if world is not None:
                self._world = world.get("ref")

    # -- tree search -------------------------------------------------------

    def find_positions(
        self,
        vol_patterns: list[str],
        world_vol: str | None = None,
        verbose: bool = False,
    ) -> list[tuple[str, np.ndarray, np.ndarray]]:
        """
        Return (name, world_pos, world_rot) for every placement of volumes
        whose name matches any pattern in vol_patterns.
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
            for (dname, _, _) in pdata["daughters"]:
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

        results: list[tuple[str, np.ndarray, np.ndarray]] = []

        def walk(vol_name, w_pos, w_rot, ancestry):
            if vol_name in ancestry:
                return
            if vol_name in target_vols:
                results.append((vol_name, w_pos.copy(), w_rot.copy()))
                # keep recursing: a target may contain other targets
            if vol_name not in self._vols or vol_name not in can_contain:
                return
            new_anc = ancestry | {vol_name}
            for dname, local_pos, local_rot in self._vols[vol_name]["daughters"]:
                if dname not in can_contain:
                    continue
                walk(dname, w_pos + w_rot @ local_pos, w_rot @ local_rot, new_anc)

        walk(world, np.zeros(3), np.eye(3), frozenset())
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
        op, first, second, pos, rot = node
        seen  = _seen | {solid}
        parts = self.solid_parts(first, seen)
        if op == "union":
            parts += [(n, pos + rot @ c, rot @ r, e)
                      for n, c, r, e in self.solid_parts(second, seen)]
        return parts

    def volume_parts(self, vol_name: str) -> list[tuple]:
        """:meth:`solid_parts` of the solid of logical volume `vol_name`."""
        vol = self._vols.get(vol_name)
        if vol is None or not vol.get("solid"):
            return []
        return self.solid_parts(vol["solid"])

    def box_extents(self, vol_name: str) -> np.ndarray | None:
        """Full (dx,dy,dz) extents [cm] for vol_name; None if not a box.

        For boolean solids (subtraction/union/intersection), follows the
        first-reference chain until a box is found (cycle-safe).
        """
        vol = self._vols.get(vol_name)
        if vol is None:
            return None
        solid = vol.get("solid", "")
        seen: set[str] = set()
        while solid and solid not in seen:
            ext = self._boxes.get(solid)
            if ext is not None:
                return ext
            seen.add(solid)
            solid = self._solid_first.get(solid, "")
        return None

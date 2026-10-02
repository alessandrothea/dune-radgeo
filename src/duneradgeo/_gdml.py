"""Pure-Python GDML parser: world-frame positions and extents of named volumes."""
from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np


def _rot_matrix_rad(rx: float, ry: float, rz: float) -> np.ndarray:
    """3×3 rotation matrix from GDML Euler angles [rad].
    Convention: intrinsic ZYX = extrinsic XYZ → R = Rz @ Ry @ Rx."""
    cx, sx = np.cos(rx), np.sin(rx)
    cy, sy = np.cos(ry), np.sin(ry)
    cz, sz = np.cos(rz), np.sin(rz)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx,  cx]])
    Ry = np.array([[cy, 0, sy], [0,  1,  0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz,  0], [0,   0,  1]])
    return Rz @ Ry @ Rx


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
            if child.tag == "box":
                name = child.get("name", "")
                u = self._ufact(child.get("lunit", "mm"))
                self._boxes[name] = np.array([
                    self._fval(child.get("x", "0")) * u,
                    self._fval(child.get("y", "0")) * u,
                    self._fval(child.get("z", "0")) * u,
                ])
        for child in solids:
            if child.tag in ("subtraction", "union", "intersection"):
                first = child.find("first")
                if first is not None:
                    self._solid_first[child.get("name", "")] = first.get("ref", "")

    def _parse_physvol(self, pv) -> tuple:
        volref = None
        pos    = np.zeros(3)
        rot    = np.eye(3)
        for child in pv:
            tag = child.tag
            if tag == "volumeref":
                volref = child.get("ref")
            elif tag == "position":
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
        return volref, pos, rot

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

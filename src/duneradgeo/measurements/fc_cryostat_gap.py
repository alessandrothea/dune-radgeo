"""fc-cryostat-gap: field-cage side walls ↔ cryostat walls (y−, y+, z−, z+)."""
from __future__ import annotations

import numpy as np

from ..ir import IRCollection
from ._core import YZ, Result, aabb, measurement, struct_box

# (wall label, axis, side): side 0 = low face, 1 = high face
_WALLS = tuple((f"{'xyz'[a]}{'-+'[s]}", a, s) for a in YZ for s in (0, 1))


def _names(solids: np.ndarray, gaps: np.ndarray, value: float) -> str:
    """Piece types whose gap equals `value` (ties listed together)."""
    return ", ".join(sorted(set(solids[np.abs(gaps - value) < 1e-6])))


@measurement("fc-cryostat-gap",
             "field-cage side walls ↔ cryostat walls (y−, y+, z−, z+)")
def fc_cryostat_gap(ir: IRCollection) -> Result:
    """
    For each of the four field-cage side walls (the top and bottom, along
    the drift, are not walls), the distance from the cryostat wall to the
    field-cage surface facing it.  The cryostat wall is the face of the
    volCryostat box, i.e. the inner surface in contact with the argon.
      gap_min  cryostat wall ↔ the closest field-shaper piece on that wall
      gap_max  cryostat wall ↔ the most recessed piece on that wall
               (e.g. slim tubes of the photon-detector window)
    """
    res = Result("fc-cryostat-gap", "field-cage side wall ↔ cryostat wall",
                 columns=["wall", "fc_solid", "n_parts", "gap_min", "gap_max",
                          "cryostat_wall", "fc_surface"])
    fshapers = ir.by_kind("fieldshaper")
    cryo = struct_box(ir, r"^volCryostat$")
    missing = [k for k, ok in (("fieldshaper", fshapers), ("struct volCryostat", cryo))
               if not ok]
    if missing:
        res.summary.append(f"no {' or '.join(missing)} volumes in the IR")
        return res

    # 1. Boxes.  Field-shaper pieces (one box per tube / corner of each
    #    ring) as (N, 3) arrays of world-frame corners, and the cryostat:
    #    the volCryostat box, whose faces are the inner cryostat walls.
    fboxes   = [aabb(v) for v in fshapers]
    flo, fhi = np.array([b[0] for b in fboxes]), np.array([b[1] for b in fboxes])
    clo, chi = cryo

    # 2. Field-cage envelope: the bounding box of all pieces.  Its y / z
    #    faces are the outermost field-cage surfaces, one per side wall.
    elo, ehi = flo.min(axis=0), fhi.max(axis=0)

    # 3. Each piece's outward face per wall and how far it sits inside the
    #    envelope face of that wall.  A piece belongs to the wall it is
    #    closest to (a ring's y-running tubes to the z walls, its z-running
    #    tubes to the y walls; corners to whichever is nearer, ties to y).
    face  = np.column_stack([fhi[:, a] if s else flo[:, a] for _, a, s in _WALLS])
    inset = np.column_stack([ehi[a] - fhi[:, a] if s else flo[:, a] - elo[a]
                             for _, a, s in _WALLS])
    owner = np.argmin(inset, axis=1)

    solids = np.array([v.meta.get("part", v.name) for v in fshapers])
    for k, (wall, a, s) in enumerate(_WALLS):
        mine = owner == k
        if not mine.any():
            continue
        # 4. Gap to the cryostat wall, measured along the wall normal:
        #      low face  (y−, z−):  piece face − cryostat face
        #      high face (y+, z+):  cryostat face − piece face
        cwall = chi[a] if s else clo[a]
        gaps  = (cwall - face[mine, k]) if s else (face[mine, k] - cwall)
        sol   = solids[mine]
        gmin, gmax = float(gaps.min()), float(gaps.max())
        res.values.update({f"{wall}_gap_min": gmin, f"{wall}_gap_max": gmax})
        res.summary.append(
            f"{wall:<3} gap {gmin:9.3f}  ({_names(sol, gaps, gmin)})   "
            f"up to {gmax:.3f}  ({_names(sol, gaps, gmax)})   "
            f"cryostat wall {'xyz'[a]}={cwall:.3f}")
        # 5. Detail rows: per wall and piece type, so e.g. standard and slim
        #    tubes show up separately.
        for name in sorted(set(sol)):
            g = gaps[sol == name]
            fs = face[mine, k][sol == name]
            res.rows.append(dict(
                wall=wall, fc_solid=str(name), n_parts=int(g.size),
                gap_min=float(g.min()), gap_max=float(g.max()),
                cryostat_wall=float(cwall),
                fc_surface=float(fs.max() if s else fs.min())))
    return res

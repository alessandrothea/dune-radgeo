"""opdet-fc-gap: wall optical detectors ↔ closest field-cage surface."""
from __future__ import annotations

import statistics

import numpy as np

from ..ir import IRCollection
from ._core import YZ, Result, aabb, box_distances, measurement, tags


@measurement("opdet-fc-gap",
             "wall optical detectors ↔ closest field-cage surface")
def opdet_fc_gap(ir: IRCollection) -> Result:
    """
    Wall optical detectors are the opdet boxes whose centre lies outside the
    field cage's y–z envelope (the cathode-mounted ones lie inside it); the
    wall is the envelope face they are beyond (y−, y+, z−, z+).  For each,
    the closest field-shaper part and the shortest distance between the two
    boxes.  Field-shaper tubes are stored as their bounding boxes, so for a
    detector facing a straight tube the gap to the tube surface is exact;
    for a corner (torus) part it is a lower bound.
    """
    res = Result("opdet-fc-gap", "wall optical detector ↔ closest field-cage surface",
                 columns=["wall", "opdet", "gap", "fc_part", "fc_solid",
                          "opdet_x", "opdet_y", "opdet_z"])
    opdets, fshapers = ir.by_kind("opdet"), ir.by_kind("fieldshaper")
    missing = [k for k, vs in (("opdet", opdets), ("fieldshaper", fshapers)) if not vs]
    if missing:
        res.summary.append(f"no {' or '.join(missing)} volumes in the IR")
        return res

    # 1. Field-cage boxes.  Each field-shaper ring is stored as one box per
    #    tube / corner piece (8 or 12 per ring); collect their world-frame
    #    (lo, hi) corners into arrays so distances can be vectorised.
    fboxes   = [aabb(v) for v in fshapers]
    flo, fhi = np.array([b[0] for b in fboxes]), np.array([b[1] for b in fboxes])
    # 2. Field-cage envelope: the bounding box of all pieces.  Its y and z
    #    faces are the outer surfaces of the cage walls.
    elo, ehi = flo.min(axis=0), fhi.max(axis=0)
    ftags    = tags(fshapers)

    n_inner = 0
    for otag, v in zip(tags(opdets), opdets):
        lo, hi = aabb(v)
        c = (lo + hi) / 2
        # 3. Wall detector or not.  How far the detector centre lies beyond
        #    each y / z face of the envelope (positive = outside that face).
        #    The largest value names the wall it is mounted on; if no value is
        #    positive the detector is inside the cage (cathode-mounted) and is
        #    skipped.  x is ignored: the cage is open towards anode and cathode.
        beyond = {}
        for a in YZ:
            beyond[f"{'xyz'[a]}-"] = elo[a] - c[a]
            beyond[f"{'xyz'[a]}+"] = c[a] - ehi[a]
        wall = max(beyond, key=beyond.get)
        if beyond[wall] <= 0:
            n_inner += 1
            continue
        # 4. Closest field-cage surface.  Shortest distance from the detector
        #    box to every field-shaper piece box; the minimum is the gap and
        #    its piece the closest one.  Per axis the separation is the space
        #    between the facing faces (0 if the boxes overlap on that axis);
        #    a detector facing a tube overlaps it in the other two axes, so
        #    the gap reduces to face-to-tube-surface along the wall normal.
        d = box_distances(lo, hi, flo, fhi)
        j = int(np.argmin(d))
        res.rows.append(dict(
            wall=wall, opdet=otag, gap=float(d[j]), fc_part=ftags[j],
            fc_solid=fshapers[j].meta.get("part", fshapers[j].name),
            opdet_x=float(c[0]), opdet_y=float(c[1]), opdet_z=float(c[2])))

    if not res.rows:
        res.summary.append("no optical detectors outside the field cage")
    # 5. Per-wall summary: mean / min / max gap over the wall's detectors, and
    #    which kind of piece was closest (e.g. slim vs. standard tubes).
    for wall in ("y-", "y+", "z-", "z+"):
        rows = [r for r in res.rows if r["wall"] == wall]
        if not rows:
            continue
        gaps = [r["gap"] for r in rows]
        mean = statistics.mean(gaps)
        res.values.update({f"{wall}_gap_mean": mean, f"{wall}_gap_min": min(gaps)})
        res.summary.append(
            f"{wall:<3} gap {mean:9.3f}  (min {min(gaps):.3f}, max {max(gaps):.3f}, "
            f"{len(rows)} detectors)")
        solids: dict[str, int] = {}
        for r in rows:
            solids[r["fc_solid"]] = solids.get(r["fc_solid"], 0) + 1
        res.summary.append("      closest: " + ", ".join(
            f"{n} ({k})" for n, k in sorted(solids.items(), key=lambda t: -t[1])))
    if n_inner:
        res.summary.append(f"{n_inner} optical detector(s) inside the field cage "
                           f"(cathode-mounted) skipped")
    return res

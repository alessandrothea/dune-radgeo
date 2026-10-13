"""anode-cryostat-gap: top / bottom anodes ↔ cryostat roof / floor along x."""
from __future__ import annotations

import statistics

from ..ir import IRCollection
from ._core import X, Result, aabb, measurement, struct_box, tags


@measurement("anode-cryostat-gap",
             "top / bottom anodes ↔ cryostat top / bottom face along the drift axis")
def anode_cryostat_gap(ir: IRCollection) -> Result:
    """
    Distance along x from each anode to the cryostat face behind it: top
    anodes to the cryostat top face (roof), bottom anodes to the bottom face
    (floor).  The cryostat faces are those of the volCryostat box (inner
    cryostat surface; on top it encloses the gas ullage).
      gap_inner  from the anode's inner face (the one facing the drift volume)
      gap_outer  from the anode's outer face (the one facing the cryostat)
    On the top side the liquid surface (bottom of volGaseousArgon) is shown
    for reference.
    """
    res = Result("anode-cryostat-gap", "anode ↔ cryostat top / bottom face along x",
                 columns=["side", "anode", "gap_inner", "gap_outer",
                          "inner_face", "outer_face", "cryostat_face"])
    anodes = ir.by_kind("anode")
    cryo = struct_box(ir, r"^volCryostat$")
    missing = [k for k, ok in (("anode", anodes), ("struct volCryostat", cryo)) if not ok]
    if missing:
        res.summary.append(f"no {' or '.join(missing)} volumes in the IR")
        return res
    clo, chi = cryo

    # 1. Top / bottom split.  The plane dividing the two drift regions is the
    #    cathode (mean centre x of its modules); without cathode entries the
    #    cryostat's mid-height is used instead.
    cathodes = ir.by_kind("cathode")
    if cathodes:
        split = statistics.mean(float(sum(aabb(c)[i][X] for i in (0, 1)) / 2)
                                for c in cathodes)
    else:
        split = float(clo[X] + chi[X]) / 2

    for tag, a in zip(tags(anodes), anodes):
        lo, hi = aabb(a)
        top = (lo[X] + hi[X]) / 2 > split
        # 2. Faces.  A top anode faces the drift volume with its lower face
        #    and the cryostat roof with its upper face; a bottom anode the
        #    other way round.
        inner, outer = (lo[X], hi[X]) if top else (hi[X], lo[X])
        # 3. Gaps along x to the cryostat face behind the anode:
        #      top:     cryostat top face − anode face
        #      bottom:  anode face − cryostat bottom face
        cface = chi[X] if top else clo[X]
        sign = 1 if top else -1
        res.rows.append(dict(
            side="top" if top else "bottom", anode=tag,
            gap_inner=float(sign * (cface - inner)),
            gap_outer=float(sign * (cface - outer)),
            inner_face=float(inner), outer_face=float(outer),
            cryostat_face=float(cface)))

    # 4. Liquid surface for reference: the bottom of the gas ullage.
    gas = struct_box(ir, r"^volGaseousArgon$")

    # 5. Per-side summary: mean / min / max over the side's anodes (equal for
    #    a flat anode plane parallel to the cryostat face).
    for side in ("top", "bottom"):
        rows = [r for r in res.rows if r["side"] == side]
        if not rows:
            res.summary.append(f"{side:<6}  no anodes")
            res.values.update({f"{side}_gap_inner": None, f"{side}_gap_outer": None})
            continue
        gi = [r["gap_inner"] for r in rows]
        go = [r["gap_outer"] for r in rows]
        res.values.update({f"{side}_gap_inner": statistics.mean(gi),
                           f"{side}_gap_outer": statistics.mean(go)})
        # each line names both surfaces of the gap: the anode face (inner =
        # towards the drift volume, outer = towards the cryostat) and the
        # cryostat face behind it
        cname = f"cryostat {side}"
        res.summary.append(
            f"{side:<6}  anode inner face ↔ {cname:<15} {statistics.mean(gi):9.3f}  "
            f"(min {min(gi):.3f}, max {max(gi):.3f}, {len(rows)} anodes)")
        res.summary.append(
            f"        anode outer face ↔ {cname:<15} {statistics.mean(go):9.3f}")
        res.summary.append(f"        {cname} face at x={rows[0]['cryostat_face']:.3f}")
        if side == "top" and gas is not None:
            surf = float(gas[0][X])
            below = statistics.mean(r["inner_face"] for r in rows)
            res.values["liquid_surface_x"] = surf
            res.summary.append(
                f"        liquid surface at x={surf:.3f}  (gas ullage "
                f"{float(gas[1][X]) - surf:.3f} thick); anode inner face "
                f"{surf - below:.3f} below it")
    return res

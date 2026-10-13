"""drift-gaps: cathode ↔ top / bottom anode distance along the drift axis."""
from __future__ import annotations

import re
import statistics

import numpy as np

from ..ir import IRCollection
from ._core import X, YZ, Result, aabb, measurement, overlap_area, tags

GAP_TOLERANCE = 0.01   # cm: pairs further than this from their side's median are flagged


def _facing(lo, hi, boxes: list[tuple]) -> int | None:
    """Index of the box in `boxes` with the largest y–z overlap, or None."""
    areas = [overlap_area(lo, hi, blo, bhi, YZ) for blo, bhi in boxes]
    best = int(np.argmax(areas))
    return best if areas[best] > 0 else None


@measurement("drift-gaps",
             "cathode ↔ top / bottom anode distance along the drift axis")
def drift_gaps(ir: IRCollection) -> Result:
    """
    Pair every anode with the cathode module it faces (largest y–z overlap)
    and measure along x:
      gap          facing surfaces: anode bottom − cathode top (top side),
                   cathode bottom − anode top (bottom side)
      centre_dist  |anode centre − cathode centre|
    Anodes above the cathode are "top", below "bottom".  The volTPCActive
    x-extent on each side is reported as a cross-check.
    """
    res = Result("drift-gaps", "cathode ↔ anode distance along x (drift)",
                 columns=["side", "anode", "cathode", "gap", "centre_dist",
                          "anode_x", "cathode_x"])
    anodes, cathodes = ir.by_kind("anode"), ir.by_kind("cathode")
    missing = [k for k, vs in (("anode", anodes), ("cathode", cathodes)) if not vs]
    if missing:
        res.summary.append(f"no {' or '.join(missing)} volumes in the IR")
        return res

    # 1. World-frame boxes.  Every IR entry is a (possibly rotated) box;
    #    aabb() turns it into axis-aligned (lo, hi) corners so faces along
    #    x can be compared directly.  The cathode is the outer box of the
    #    perforated CathodeGrid solid, i.e. its full thickness (~4 cm).
    cboxes, ctags = [aabb(c) for c in cathodes], tags(cathodes)
    unpaired = []
    for atag, a in zip(tags(anodes), anodes):
        alo, ahi = aabb(a)
        # 2. Pairing.  The anode's partner is the cathode module whose y–z
        #    footprint overlaps it most, i.e. the one directly across the
        #    drift.  Several anodes may share one cathode module (v7: 80
        #    cathodes for 80 top + 160 bottom anodes).  No overlap → unpaired.
        j = _facing(alo, ahi, cboxes)
        if j is None:
            unpaired.append(atag)
            continue
        clo, chi = cboxes[j]
        # 3. Side.  Decided by position, not by volume name: an anode whose
        #    centre is above its cathode's centre (larger x) is "top".  This
        #    also works for v4/v5, which have top anodes only.
        ax, cx = (alo[X] + ahi[X]) / 2, (clo[X] + chi[X]) / 2
        top = ax > cx
        # 4. Distances along x.
        #    gap:          between the facing surfaces —
        #                    top:    anode lower face − cathode upper face
        #                    bottom: cathode lower face − anode upper face
        #    centre_dist:  between centres; exceeds gap by half the cathode
        #                  plus half the anode thickness (~2 cm in VD).
        res.rows.append(dict(
            side="top" if top else "bottom", anode=atag, cathode=ctags[j],
            gap=float(alo[X] - chi[X] if top else clo[X] - ahi[X]),
            centre_dist=float(abs(ax - cx)),
            anode_x=float(ax), cathode_x=float(cx)))

    # 5. Cross-check inputs: the mean cathode plane, used to split the
    #    active volumes into top / bottom drift regions, and the
    #    volTPCActive(Inner) boxes themselves.
    cathode_x = statistics.mean((lo[X] + hi[X]) / 2 for lo, hi in cboxes)
    active = [aabb(v) for v in ir.by_kind("struct")
              if re.match(r"^volTPCActive", v.name)]

    # 6. Per-side summary.
    for side in ("top", "bottom"):
        rows = [r for r in res.rows if r["side"] == side]
        if not rows:
            res.summary.append(f"{side:<6}  no anodes")
            res.values.update({f"{side}_gap": None, f"{side}_centre_dist": None})
            continue
        # mean / min / max over all anode–cathode pairs on this side; for
        # flat, parallel planes they coincide
        gaps = [r["gap"] for r in rows]
        mean = statistics.mean(gaps)
        cmean = statistics.mean(r["centre_dist"] for r in rows)
        res.values.update({f"{side}_gap": mean, f"{side}_centre_dist": cmean})
        res.summary.append(
            f"{side:<6}  gap {mean:9.3f}  (min {min(gaps):.3f}, max {max(gaps):.3f}, "
            f"{len(rows)} anodes)   centre-to-centre {cmean:9.3f}")
        # pairs that deviate from the side's typical (median) gap point to a
        # tilted or misplaced module
        med = statistics.median(gaps)
        off = [r for r in rows if abs(r["gap"] - med) > GAP_TOLERANCE]
        if off:
            res.summary.append(f"        {len(off)} pair(s) differ from the median gap "
                               f"by > {GAP_TOLERANCE} cm (see --detail)")
        # active-volume x-extent on the same side of the cathode: the
        # simulated drift length, expected to be close to the gap
        above = side == "top"
        ext = [hi[X] - lo[X] for lo, hi in active
               if ((lo[X] + hi[X]) / 2 > cathode_x) == above]
        if ext:
            res.values[f"{side}_active_x"] = statistics.mean(ext)
            res.summary.append(f"        volTPCActive x-extent {statistics.mean(ext):9.3f}  "
                               f"({len(ext)} volumes)")
    if unpaired:
        res.summary.append(f"{len(unpaired)} anode(s) face no cathode: "
                           f"{', '.join(unpaired[:5])}{' …' if len(unpaired) > 5 else ''}")
    return res

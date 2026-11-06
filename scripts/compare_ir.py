#!/usr/bin/env python3
"""
Compare two JSON IR files volume by volume (regression check).

Usage:
    compare_ir.py ref.json new.json
    compare_ir.py ref.json new.json --tol 0.1 --ignore-rot '^volTPC'
    compare_ir.py ref.json new.json --only-common -v

Entries are grouped by (kind, name, part) and each reference entry is paired
with the nearest entry of the same group in the other IR.  A group passes
when both IRs have the same number of entries, the pairing is one-to-one,
and every pair agrees in centre (--tol), extents (--tol) and rotation
(--rot-tol).  Extract both IRs with the same --frame (e.g. volEnclosureTPC)
when the GDMLs nest the detector differently.

Exits with status 1 if any group fails.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict

import numpy as np

from duneradgeo.ir import IRCollection


def _groups(ir: IRCollection, kinds: set[str] | None) -> dict[tuple, list]:
    g: dict[tuple, list] = defaultdict(list)
    for v in ir.volumes:
        if kinds and v.kind not in kinds:
            continue
        g[(v.kind, v.name, v.meta.get("part", ""))].append(v)
    return g


def _compare(a: list, b: list, tol: float, rot_tol: float, check_rot: bool) -> tuple[bool, str]:
    pa = np.array([v.pos for v in a])
    pb = np.array([v.pos for v in b])
    d2 = ((pa[:, None, :] - pb[None, :, :]) ** 2).sum(axis=2)
    j  = d2.argmin(axis=1)
    dpos = np.sqrt(d2[np.arange(len(a)), j])
    dext = np.array([np.abs(a[i].ext - b[j[i]].ext).max() for i in range(len(a))])
    drot = np.array([np.abs(a[i].rot - b[j[i]].rot).max() for i in range(len(a))])
    n_pos = int((dpos > tol).sum())
    n_ext = int((dext > tol).sum())
    n_rot = int((drot > rot_tol).sum()) if check_rot else 0
    unique = len(set(j.tolist())) == len(a)
    ok = len(a) == len(b) and unique and not (n_pos or n_ext or n_rot)
    msg = (f"n={len(a)}/{len(b)}  max|dpos|={dpos.max():.3f}  max|dext|={dext.max():.3f}"
           f"  max|drot|={drot.max():.2g}")
    bad = []
    if len(a) != len(b):
        bad.append("count")
    if not unique:
        bad.append("pairing not one-to-one")
    if n_pos:
        shifts = Counter(tuple(np.round(pb[j[i]] - pa[i], 2)) for i in np.where(dpos > tol)[0])
        top = ", ".join(f"{n}×({s[0]:g},{s[1]:g},{s[2]:g})" for s, n in shifts.most_common(3))
        bad.append(f"{n_pos} pos [{top}]")
    if n_ext:
        bad.append(f"{n_ext} ext")
    if n_rot:
        bad.append(f"{n_rot} rot")
    if bad:
        msg += "  ✗ " + "; ".join(bad)
    elif not check_rot and (drot > rot_tol).any():
        msg += f"  ({int((drot > rot_tol).sum())} rot, ignored)"
    return ok, msg


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("ref", help="reference JSON IR")
    ap.add_argument("new", help="JSON IR to check against the reference")
    ap.add_argument("--tol", type=float, default=0.01,
                    help="position / extent tolerance in cm (default: 0.01)")
    ap.add_argument("--rot-tol", type=float, default=1e-6,
                    help="max abs difference of rotation-matrix elements (default: 1e-6)")
    ap.add_argument("--kinds", nargs="+", default=None,
                    help="only compare these kinds (anode opdet cathode fieldshaper struct radio)")
    ap.add_argument("--ignore-rot", nargs="+", default=[], metavar="NAME_RE",
                    help="do not fail on rotation differences of these volume names")
    ap.add_argument("--only-common", action="store_true",
                    help="do not fail on groups present in only one IR")
    ap.add_argument("--verbose", "-v", action="store_true",
                    help="also list the groups that pass")
    args = ap.parse_args()

    ra, rb = IRCollection.load(args.ref), IRCollection.load(args.new)
    fa = (ra.frame or {}).get("volume", "world")
    fb = (rb.frame or {}).get("volume", "world")
    print(f"[info] ref: {len(ra.volumes)} entries, frame {fa}   new: {len(rb.volumes)} entries, frame {fb}")
    if fa != fb:
        print(f"[warn] frames differ ({fa} vs {fb}); positions are not comparable",
              file=sys.stderr)

    kinds = set(args.kinds) if args.kinds else None
    ga, gb = _groups(ra, kinds), _groups(rb, kinds)
    ignore_rot = [re.compile(p) for p in args.ignore_rot]

    failed = 0
    for key in sorted(set(ga) | set(gb)):
        kind, name, part = key
        label = f"{kind:<11} {name}{'/' + part if part else ''}"
        if key not in ga or key not in gb:
            side = "ref" if key in ga else "new"
            print(f"{'skip' if args.only_common else 'FAIL'}  {label:<50} only in {side}")
            failed += not args.only_common
            continue
        check_rot = not any(p.search(name) for p in ignore_rot)
        ok, msg = _compare(ga[key], gb[key], args.tol, args.rot_tol, check_rot)
        failed += not ok
        if not ok or args.verbose:
            print(f"{'ok  ' if ok else 'FAIL'}  {label:<50} {msg}")

    n = len(set(ga) | set(gb))
    print(f"[info] {n - failed}/{n} groups agree")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()

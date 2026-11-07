#!/usr/bin/env python3
"""
Compare the python and root GDML backends on the same GDML files.

Usage (the root backend needs ROOT with TGeo/GDML, e.g. in the dunesw
container: dsw py scripts/compare_backends.py ...):

    compare_backends.py geometry.gdml [more.gdml ...]
    compare_backends.py geometry.gdml --patterns '^volCryostat$' '^volArapuca' -v

For each file both backends are loaded and compared on

  * the logical volume names and the world volume;
  * the solid tables: primitive bounding boxes (centre, extents), box
    extents and boolean operands with their placements;
  * every placement matched by the detector profile's patterns (anode,
    opdet, cathode, fieldshaper, struct and passive, or --patterns):
    paired by placement path, compared in world position, rotation and
    physvol names, and in the box parts of the volume.

Exits with status 1 if anything differs beyond --tol / --rot-tol.
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np

from duneradgeo._gdml import load_geometry
from duneradgeo._profiles import detect_profile, get_profile

_KINDS = ("anode", "opdet", "cathode", "fieldshaper", "struct", "passive")


def _load(path: str, backend: str, verbose: bool):
    t0 = time.time()
    geo = load_geometry(path, backend=backend, verbose=verbose)
    print(f"[info] {backend:<6} {len(geo.volume_names)} volumes  "
          f"({time.time() - t0:.1f}s)")
    return geo


class _Report:
    def __init__(self, verbose: bool):
        self.verbose = verbose
        self.failed  = 0

    def line(self, ok: bool, label: str, msg: str, details: list[str] = ()) -> None:
        self.failed += not ok
        if ok and not self.verbose:
            return
        print(f"{'ok  ' if ok else 'FAIL'}  {label:<34} {msg}")
        for d in list(details)[:8]:
            print(f"        {d}")
        if len(details) > 8:
            print(f"        … {len(details) - 8} more")


def _maxdiff(a, b) -> float:
    return float(np.abs(np.asarray(a, float) - np.asarray(b, float)).max())


def _cmp_solids(py, rt, rep: _Report, tol: float, rot_tol: float) -> None:
    for attr, what in (("_prims", "primitive bboxes"), ("_boxes", "box extents")):
        a, b = getattr(py, attr), getattr(rt, attr)
        only = sorted(set(a) ^ set(b))
        bad = []
        for n in sorted(set(a) & set(b)):
            if attr == "_prims":
                (ca, ea), (cb, eb) = a[n], b[n]
                d = max(_maxdiff(ca, cb), _maxdiff(ea, eb))
                if d > tol:
                    bad.append(f"{n}: python c={np.round(ca, 3).tolist()} ext={np.round(ea, 3).tolist()}"
                               f"  root c={np.round(cb, 3).tolist()} ext={np.round(eb, 3).tolist()}")
            elif _maxdiff(a[n], b[n]) > tol:
                bad.append(f"{n}: python {np.round(a[n], 3).tolist()}  root {np.round(b[n], 3).tolist()}")
        details = [f"only in {'python' if n in a else 'root'}: {n}" for n in only] + bad
        rep.line(not details, what, f"n={len(a)}/{len(b)}  {len(bad)} differ, {len(only)} unpaired",
                 details)

    a, b = py._bools, rt._bools
    only = sorted(set(a) ^ set(b))
    bad = []
    for n in sorted(set(a) & set(b)):
        opa, fa, sa, pa, ra, fpa, fra = a[n]
        opb, fb, sb, pb, rb, fpb, frb = b[n]
        why = []
        if (opa, fa, sa) != (opb, fb, sb):
            why.append(f"operands python {(opa, fa, sa)} root {(opb, fb, sb)}")
        if _maxdiff(pa, pb) > tol or _maxdiff(fpa, fpb) > tol:
            why.append(f"pos python {np.round(pa, 3).tolist()} root {np.round(pb, 3).tolist()}"
                       + (f" first python {np.round(fpa, 3).tolist()} root {np.round(fpb, 3).tolist()}"
                          if _maxdiff(fpa, fpb) > tol else ""))
        if _maxdiff(ra, rb) > rot_tol or _maxdiff(fra, frb) > rot_tol:
            why.append(f"rot |d|={max(_maxdiff(ra, rb), _maxdiff(fra, frb)):.2g}")
        if why:
            bad.append(f"{n}: " + "; ".join(why))
    details = [f"only in {'python' if n in a else 'root'}: {n}" for n in only] + bad
    rep.line(not details, "boolean solids", f"n={len(a)}/{len(b)}  {len(bad)} differ, {len(only)} unpaired",
             details)


def _cmp_parts(pa: list, pb: list, tol: float, rot_tol: float) -> str | None:
    if [p[0] for p in pa] != [p[0] for p in pb]:
        return f"parts python {[p[0] for p in pa]} root {[p[0] for p in pb]}"
    d = max((max(_maxdiff(a[1], b[1]), _maxdiff(a[3], b[3])) for a, b in zip(pa, pb)), default=0.0)
    r = max((_maxdiff(a[2], b[2]) for a, b in zip(pa, pb)), default=0.0)
    if d > tol or r > rot_tol:
        worst = max(zip(pa, pb), key=lambda ab: max(_maxdiff(ab[0][1], ab[1][1]),
                                                   _maxdiff(ab[0][3], ab[1][3])))
        return (f"part {worst[0][0]}: |dc,dext|={d:.3f} |drot|={r:.2g}  python ext="
                f"{np.round(worst[0][3], 3).tolist()} root ext={np.round(worst[1][3], 3).tolist()}")
    return None


def _cmp_placements(py, rt, kind: str, patterns: list[str], rep: _Report,
                    tol: float, rot_tol: float) -> None:
    ha = {h[3]: h for h in py.find_positions(patterns, with_path=True)}
    hb = {h[3]: h for h in rt.find_positions(patterns, with_path=True)}
    only = sorted(set(ha) ^ set(hb))
    common = [p for p in ha if p in hb]
    dpos = drot = 0.0
    bad_pv: list[str] = []
    bad_xf: list[str] = []
    for p in common:
        a, b = ha[p], hb[p]
        dp, dr = _maxdiff(a[1], b[1]), _maxdiff(a[2], b[2])
        dpos, drot = max(dpos, dp), max(drot, dr)
        if dp > tol or dr > rot_tol:
            bad_xf.append(f"{p}: |dpos|={dp:.3f} |drot|={dr:.2g}")
        if (a[4], a[5]) != (b[4], b[5]):
            bad_pv.append(f"{p}: python pv={a[4]!r} pvpath={a[5]!r}  root pv={b[4]!r} pvpath={b[5]!r}")
    bad_parts: list[str] = []
    for name in sorted({ha[p][0] for p in common}):
        why = _cmp_parts(py.volume_parts(name), rt.volume_parts(name), tol, rot_tol)
        ea, eb = py.box_extents(name), rt.box_extents(name)
        if (ea is None) != (eb is None) or (ea is not None and _maxdiff(ea, eb) > tol):
            why = (why + "; " if why else "") + f"box_extents python {ea} root {eb}"
        if why:
            bad_parts.append(f"{name}: {why}")
    details = ([f"only in {'python' if p in ha else 'root'}: {p}" for p in only]
               + bad_xf + bad_pv + bad_parts)
    msg = (f"n={len(ha)}/{len(hb)}  max|dpos|={dpos:.3g}  max|drot|={drot:.2g}"
           f"  pv diff={len(bad_pv)}  parts diff={len(bad_parts)}")
    rep.line(not details, f"placements {kind}", msg, details)


def compare(path: str, patterns: list[str] | None, tol: float, rot_tol: float,
            verbose: bool) -> int:
    print(f"== {path}")
    py = _load(path, "python", False)
    rt = _load(path, "root", verbose)
    rep = _Report(verbose)

    va, vb = py.volume_names, rt.volume_names
    rep.line(va == vb and py._world == rt._world, "volumes / world",
             f"n={len(va)}/{len(vb)}  world {py._world}/{rt._world}",
             [f"only in python: {n}" for n in sorted(va - vb)]
             + [f"only in root: {n}" for n in sorted(vb - va)])
    _cmp_solids(py, rt, rep, tol, rot_tol)

    if patterns:
        groups = [("custom", patterns)]
    else:
        key = detect_profile(py.volume_names)
        _, prof = get_profile(key)
        print(f"[info] profile {key or 'generic'}")
        groups = [(k, prof.get(f"{k}_patterns", [])) for k in _KINDS]
    for kind, pats in groups:
        if pats:
            _cmp_placements(py, rt, kind, pats, rep, tol, rot_tol)

    print(f"[info] {'agree' if not rep.failed else f'{rep.failed} check(s) differ'}")
    return rep.failed


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gdml", nargs="+", help="GDML file(s)")
    ap.add_argument("--patterns", nargs="+", default=None, metavar="RE",
                    help="volume-name patterns to compare placements of "
                         "(default: the detected profile's)")
    ap.add_argument("--tol", type=float, default=1e-3,
                    help="position / extent tolerance in cm (default: 1e-3)")
    ap.add_argument("--rot-tol", type=float, default=1e-6,
                    help="max abs difference of rotation-matrix elements (default: 1e-6)")
    ap.add_argument("--verbose", "-v", action="store_true",
                    help="also list the checks that pass, and show ROOT's import output")
    args = ap.parse_args()

    failed = sum(compare(p, args.patterns, args.tol, args.rot_tol, args.verbose) > 0
                 for p in args.gdml)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()

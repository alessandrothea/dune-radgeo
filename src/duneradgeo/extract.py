"""
bki-extract  —  extract detector & radiological volumes → JSON IR

Usage:
  bki-extract geometry.gdml --out ir.json
  bki-extract job.fcl       --out ir.json
  bki-extract geometry.gdml --out ir.json \\
      --anode-pattern 'volAnodePlate' --opdet-pattern '^volArapuca'
  bki-extract geometry.gdml --out ir.json --frame volEnclosureTPC

If a FCL file is given:
  * #include directives are resolved automatically via $FHICL_FILE_PATH
  * radiological producers (RadioGen / Decay0Gen) → kind="radio" volumes,
    and each producer's full FHiCL table → the IR's top-level "producers"
  * the GDML filename is extracted from AuxDetGeometry.GDML and used to
    populate anode / opdet / cathode / fieldshaper / struct volumes as well.
  * producers with distrib_x/y/z: the formulas are evaluated with ROOT (TF1),
    so PyROOT is required; the box on that axis becomes the effective range
    of the distribution TF1::GetRandom generates (--npx, --containment).

The refactored and the ggd-generated VD GDMLs use the same volume names but
nest the cryostat differently (ggd adds volSupportEnc/volEnclosureCryostat),
so their world coordinates differ by a global offset.  --frame expresses
every volume in the frame of one placed volume (e.g. volEnclosureTPC) and
records its world placement in the IR, so IRs of both kinds overlay.
--with-passive adds the profile's passive material (steel shell, insulation,
support structure, shielding, cavern) as struct volumes.

--gdml-backend root reads the GDML with ROOT (TGeoManager::Import) instead
of the built-in parser; both feed the same search, so IRs from the two
backends can be compared with scripts/compare_ir.py.

Non-box solids are approximated by boxes: tubes and tori by their bounding
box, unions by one box per operand (e.g. each straight tube and corner of a
field-shaper ring), subtractions by their outer (first) operand.  A radio
volume_rand / volume_gen box is the bounding box of the volume's solid (a
union: the whole composite) in the world, as BaseRadioGen samples it.

Output is a JSON file readable by bki-plot and the nb_B notebook.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
import time
from pathlib import Path

import numpy as np

from ._draw import box_corners_world
from ._fcl import (DEFAULT_CONTAINMENT, DISTRIB_AXES, DISTRIB_COMMON,
                   FhiclParser, default_fcl_cache_dir, distrib_keys, expand_fcl,
                   producer_volumes, select_producers)
from ._gdml import BACKENDS, GeometryBase, load_geometry
from ._root import ROOTUnavailable
from ._tf1 import GENERATOR_NPX
from ._profiles import detect_profile, get_profile, profile_detector
from ._search import gdml_from_fcl, resolve_search_path
from .ir import IRCollection, VolumeEntry


class ExtractError(Exception):
    """An extraction step cannot go on; bki-extract exits with its message."""


# radio row keys copied into the IR meta (when set)
_RADIO_META = (("source", "node", "nuclide", "module_type", "BqPercc",
                "rate", "surface", "T0", "T1")
               + tuple(k for axis in DISTRIB_AXES for k in distrib_keys(axis).values())
               + DISTRIB_COMMON)

# ---------------------------------------------------------------------------
# GDML → VolumeEntry list
# ---------------------------------------------------------------------------

def _gdml_volumes(gdml_path: Path,
                  anode_pats: list[str] | None,
                  opdet_pats: list[str] | None,
                  struct_pats: list[str] | None,
                  cathode_pats: list[str] | None = None,
                  fieldshaper_pats: list[str] | None = None,
                  passive: bool = False,
                  verbose: bool = False,
                  backend: str = "python") -> tuple[list[VolumeEntry], str, GeometryBase]:
    """Parse a GDML and return (entries, detector_key, geo).

    Raises ExtractError when the backend cannot read it."""
    print(f"[info] parsing {gdml_path.name} ({backend} backend) …", file=sys.stderr)
    t0  = time.time()
    try:
        geo = load_geometry(str(gdml_path), backend=backend, verbose=verbose)
    except RuntimeError as exc:
        raise ExtractError(str(exc)) from exc
    print(f"[info] {len(geo.volume_names)} logical volumes  ({time.time()-t0:.1f}s)",
          file=sys.stderr)

    profile_key = detect_profile(geo.volume_names)
    _, prof      = get_profile(profile_key)
    if profile_key:
        print(f"[info] detected {prof['label']}  (profile {profile_key})", file=sys.stderr)
    else:
        print("[warn] unknown detector type; using generic patterns", file=sys.stderr)

    ap = anode_pats  or prof["anode_patterns"]
    op = opdet_pats  or prof["opdet_patterns"]
    cp = cathode_pats or prof.get("cathode_patterns", [])
    fp = fieldshaper_pats or prof.get("fieldshaper_patterns", [])
    sp = struct_pats or prof.get("struct_patterns", [])

    entries: list[VolumeEntry] = []

    def _collect(patterns, kind):
        hits = geo.find_positions(patterns, verbose=verbose, with_path=True)
        ok = 0
        copies: dict[str, int] = {}
        for name, pos, rot, path, pvname, pvpath in hits:
            parts = geo.volume_parts(name)
            if not parts:
                print(f"[warn] no box/tube/torus solid for {kind} {name}, skipping",
                      file=sys.stderr)
                continue
            copy = copies[name] = copies.get(name, -1) + 1
            corners = []
            for part, centre, prot, ext in parts:
                wpos, wrot = pos + rot @ centre, rot @ prot
                # multi-part volumes: tag each box with its solid and placement
                meta = {"part": part, "copy": copy} if len(parts) > 1 else {}
                meta["path"] = path
                if pvname:
                    meta["pv"] = pvname
                if pvpath:
                    meta["pvpath"] = pvpath
                entries.append(VolumeEntry(kind=kind, name=name,
                                           pos=wpos, ext=ext, rot=wrot, meta=meta))
                corners.append(box_corners_world(wpos, ext, wrot))
            corners = np.vstack(corners)
            lo, hi  = corners.min(axis=0), corners.max(axis=0)
            c, d    = (lo + hi) / 2, hi - lo
            nparts  = f"  ({len(parts)} parts)" if len(parts) > 1 else ""
            print(f"  {kind:<8} {name:<40s}"
                  f"  x={c[0]:8.2f}  y={c[1]:8.2f}  z={c[2]:8.2f}"
                  f"  dx={d[0]:8.2f}  dy={d[1]:8.2f}  dz={d[2]:8.2f}{nparts}")
            ok += 1
        print(f"[info] {ok}/{len(hits)} {kind} placements", file=sys.stderr)

    _collect(ap, "anode")
    _collect(op, "opdet")
    if cp:
        _collect(cp, "cathode")
    if fp:
        _collect(fp, "fieldshaper")
    if sp:
        _collect(sp, "struct")
    if passive:
        pp = prof.get("passive_patterns", [])
        if pp:
            _collect(pp, "struct")
        else:
            print("[warn] --with-passive: profile has no passive patterns", file=sys.stderr)

    return entries, profile_detector(profile_key), geo


def _to_frame(entries: list[VolumeEntry], geo: GeometryBase, volume: str) -> dict:
    """Re-express `entries` (world frame) in the frame of `volume`, in place.

    `volume` must be placed exactly once (else ExtractError).  Returns the IR
    frame record."""
    hits = geo.find_positions([f"^{re.escape(volume)}$"], with_path=True)
    if len(hits) != 1:
        raise ExtractError(f"--frame {volume}: expected exactly one placement, "
                           f"found {len(hits)}")
    _, t, R, path, *_ = hits[0]
    for e in entries:
        e.pos = R.T @ (e.pos - t)
        e.rot = R.T @ e.rot
    print(f"[info] frame {volume} at x={t[0]:.2f} y={t[1]:.2f} z={t[2]:.2f}"
          f"  ({path})", file=sys.stderr)
    return {"volume": volume, "path": path, "pos": t.tolist(), "rot": R.tolist()}


# ---------------------------------------------------------------------------
# FCL → radio VolumeEntry list
# ---------------------------------------------------------------------------

def _parse_fcl(fcl_path: Path, cache_dir: str | None = None) -> tuple[dict, list]:
    """Expand (via fhicl-dump or inline includes) and parse a FCL once.

    Returns (doc, top_order) so callers can reuse the parsed document instead
    of re-running the expensive expansion step.  When `cache_dir` is given the
    expanded text is persisted there and reused by later invocations.
    """
    print(f"[info] parsing FCL {fcl_path.name} …", file=sys.stderr)
    parser = FhiclParser(expand_fcl(str(fcl_path), cache_dir=cache_dir))
    doc    = parser.parse()
    return doc, parser.order


def _fcl_radio_volumes(doc: dict,
                       top_order: list,
                       geo: GeometryBase | None = None,
                       prefix: str = "",
                       all_configured: bool = False,
                       select_re: str | None = None,
                       containment: float = DEFAULT_CONTAINMENT,
                       npx: int = GENERATOR_NPX) -> tuple[list[VolumeEntry], dict]:
    """Build radio VolumeEntry rows from an already-parsed FCL document.

    Returns (entries, producers): the rows, and the full FHiCL table of every
    selected producer keyed by label (kept even when no row resolves).
    Raises ExtractError when a producer's distrib_* formula needs PyROOT and
    it is missing; other producer errors are warnings."""
    table, labels = select_producers(doc, top_order,
                                     prefix=prefix,
                                     all_configured=all_configured,
                                     select_re=select_re)
    producers = {name: _json_safe(table[name]) for name in labels}
    entries: list[VolumeEntry] = []
    for name in labels:
        try:
            rows = producer_volumes(name, table[name],
                                    geo=geo, containment=containment, npx=npx)
        except ROOTUnavailable as exc:
            raise ExtractError(f"{name}: {exc}") from exc
        except Exception as exc:
            print(f"[warn] {name}: {exc}", file=sys.stderr)
            continue

        for r in rows:
            lo, hi = r.get("lo"), r.get("hi")
            if lo is None or hi is None or None in lo or None in hi:
                print(f"[warn] {name}#{r['index']}: unresolved box (needs --gdml), skipping",
                      file=sys.stderr)
                continue
            lo_a, hi_a = np.array(lo, dtype=float), np.array(hi, dtype=float)
            pos = (lo_a + hi_a) / 2
            ext = hi_a - lo_a
            rot = np.eye(3)
            meta = {k: r.get(k) for k in _RADIO_META if r.get(k) is not None}
            meta["index"] = r["index"]
            entries.append(VolumeEntry(kind="radio", name=name,
                                       pos=pos, ext=ext, rot=rot, meta=meta))
    print(f"[info] {len(entries)} radio volume entries", file=sys.stderr)
    return entries, producers


def _json_safe(v):
    """Deep copy of a parsed FHiCL value with ±inf/nan as strings (strict JSON)."""
    if isinstance(v, dict):
        return {k: _json_safe(x) for k, x in v.items()}
    if isinstance(v, list):
        return [_json_safe(x) for x in v]
    if isinstance(v, float) and not math.isfinite(v):
        return str(v)
    return v


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    try:
        _main()
    except ExtractError as exc:
        sys.exit(f"[error] {exc}")


def _main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("input",
                    help="GDML geometry file or expanded FHiCL job file")
    ap.add_argument("--out", "-o", required=True,
                    help="output JSON file (IR)")
    # geometry patterns
    ap.add_argument("--anode-pattern",  nargs="+", default=None)
    ap.add_argument("--opdet-pattern",  nargs="+", default=None)
    ap.add_argument("--cathode-pattern", nargs="+", default=None)
    ap.add_argument("--fieldshaper-pattern", nargs="+", default=None)
    ap.add_argument("--struct-pattern", nargs="+", default=None)
    ap.add_argument("--struct", nargs="+", default=None, metavar="VOLNAME",
                    help="exact structural volume names (anchored, overrides --struct-pattern)")
    ap.add_argument("--with-passive", action="store_true",
                    help="also extract the profile's passive material "
                         "(steel shell, insulation, support, shielding, cavern) as struct")
    ap.add_argument("--frame", default=None, metavar="VOLNAME",
                    help="express positions in the frame of this (singly placed) "
                         "volume, e.g. volCryostat or volEnclosureTPC, instead of "
                         "the world")
    # FCL options
    ap.add_argument("--prefix",         default="",
                    help="only producers whose label starts with this")
    ap.add_argument("--select",         default=None,
                    help="regex filter on producer labels")
    ap.add_argument("--all-configured", action="store_true",
                    help="include producers not on any trigger path")
    ap.add_argument("--containment",    type=float, default=DEFAULT_CONTAINMENT,
                    help="distrib_x/y/z containment fraction for the effective "
                         "range (default %(default)s)")
    ap.add_argument("--npx",            type=int, default=GENERATOR_NPX,
                    help="TF1 Npx of the generator: the effective range is that of "
                         "the distribution TF1::GetRandom draws from an Npx-point "
                         "table; larsimrad BaseRadioGen keeps ROOT's default "
                         "(default: %(default)s)")
    ap.add_argument("--no-geometry",    action="store_true",
                    help="skip GDML parsing even when a FCL is given")
    ap.add_argument("--gdml",           default=None,
                    help="override GDML path (when input is FCL)")
    ap.add_argument("--gdml-backend",   choices=BACKENDS, default="python",
                    help="GDML reader: the pure-Python parser, or ROOT's "
                         "TGeoManager::Import (needs ROOT with TGeo/GDML, "
                         "e.g. in the dunesw container) (default: %(default)s)")
    ap.add_argument("--fcl-cache",      default=default_fcl_cache_dir(), metavar="DIR",
                    help="directory where expanded FCL files are cached and "
                         "reused across runs (default: %(default)s)")
    ap.add_argument("--no-fcl-cache",   action="store_true",
                    help="always re-run the FCL expansion, ignoring the cache")
    ap.add_argument("--verbose", "-v",  action="store_true")
    args = ap.parse_args()

    inp = Path(args.input)
    is_fcl = inp.suffix.lower() == ".fcl"
    if not inp.is_file():
        search_var = "FHICL_FILE_PATH" if is_fcl else "FW_SEARCH_PATH"
        resolved = resolve_search_path(args.input, var=search_var)
        if resolved is None:
            raise ExtractError(f"'{args.input}' not found on disk or in ${search_var}")
        inp = resolved
        is_fcl = inp.suffix.lower() == ".fcl"  # re-check after resolution

    entries:     list[VolumeEntry] = []
    detector_key = ""
    source       = str(inp)

    # ---- Expand + parse the FCL exactly once (fhicl-dump is expensive) ----
    fcl_doc:   dict | None = None
    fcl_order: list        = []
    if is_fcl:
        cache_dir = None if args.no_fcl_cache else args.fcl_cache
        fcl_doc, fcl_order = _parse_fcl(inp, cache_dir=cache_dir)

    # ---- Resolve GDML path (needed before radio extraction for volume_rand/gen) ----
    gdml_path: Path | None = None
    if not args.no_geometry:
        if args.gdml:
            gdml_path = Path(args.gdml)
            if not gdml_path.is_file():
                gdml_path = resolve_search_path(args.gdml)
                if gdml_path is None:
                    raise ExtractError(f"--gdml '{args.gdml}' not found")
        elif is_fcl:
            name = gdml_from_fcl(fcl_doc)
            if name:
                gdml_path = resolve_search_path(name)
                if gdml_path is None:
                    print(f"[warn] '{name}' not found in $FW_SEARCH_PATH; "
                          "radio volume_rand/gen boxes will be unresolved",
                          file=sys.stderr)
        else:
            gdml_path = inp  # input is already GDML

    # ---- Parse GDML (shared by detector volumes and FCL radio resolution) ----
    geo: GeometryBase | None = None
    if gdml_path is not None:
        struct_pats = (
            [f"^{n}$" for n in args.struct] if args.struct
            else args.struct_pattern
        )
        geo_entries, detector_key, geo = _gdml_volumes(
            gdml_path,
            anode_pats=args.anode_pattern,
            opdet_pats=args.opdet_pattern,
            struct_pats=struct_pats,
            cathode_pats=args.cathode_pattern,
            fieldshaper_pats=args.fieldshaper_pattern,
            passive=args.with_passive,
            verbose=args.verbose,
            backend=args.gdml_backend,
        )
        entries += geo_entries
        if not is_fcl:
            source = str(gdml_path)

    # ---- FCL input: radio volumes (geo available for volume_rand/gen) ----
    producers: dict = {}
    if is_fcl:
        radio, producers = _fcl_radio_volumes(
            fcl_doc, fcl_order,
            geo=geo,
            prefix=args.prefix,
            all_configured=args.all_configured,
            select_re=args.select,
            containment=args.containment,
            npx=args.npx,
        )
        entries += radio

    if not entries:
        raise ExtractError("nothing extracted — check patterns or use --verbose")

    frame = None
    if args.frame and args.frame != "world":
        if geo is None:
            raise ExtractError("--frame needs the GDML (drop --no-geometry)")
        frame = _to_frame(entries, geo, args.frame)

    ir = IRCollection(source=source, detector=detector_key, volumes=entries,
                      frame=frame, producers=producers)
    ir.save(args.out)
    print(f"[info] wrote {len(entries)} volume entries"
          + (f" and {len(producers)} producer configs" if producers else "")
          + f" to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()

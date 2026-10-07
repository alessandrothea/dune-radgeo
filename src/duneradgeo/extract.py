"""
bki-extract  —  extract detector & radiological volumes → JSON IR

Usage:
  bki-extract geometry.gdml --out ir.json
  bki-extract job.fcl       --out ir.json
  bki-extract geometry.gdml --out ir.json \\
      --anode-pattern 'volAnodePlate' --opdet-pattern '^volArapuca'

If a FCL file is given:
  * #include directives are resolved automatically via $FHICL_FILE_PATH
  * radiological producers (RadioGen / Decay0Gen) → kind="radio" volumes
  * the GDML filename is extracted from AuxDetGeometry.GDML and used to
    populate anode / opdet / struct volumes as well.

Output is a JSON file readable by bki-plot and the nb_B notebook.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from ._draw import box_corners_world
from ._fcl import (FhiclParser, default_fcl_cache_dir, expand_fcl,
                   producer_volumes, select_producers)
from ._gdml import GDMLGeometry
from ._profiles import detect_profile, get_profile
from ._search import gdml_from_fcl, resolve_search_path
from .ir import IRCollection, VolumeEntry


# ---------------------------------------------------------------------------
# GDML → VolumeEntry list
# ---------------------------------------------------------------------------

def _gdml_volumes(gdml_path: Path,
                  anode_pats: list[str] | None,
                  opdet_pats: list[str] | None,
                  struct_pats: list[str] | None,
                  verbose: bool = False) -> tuple[list[VolumeEntry], str, GDMLGeometry]:
    """Parse a GDML and return (entries, detector_key, geo)."""
    print(f"[info] parsing {gdml_path.name} …", file=sys.stderr)
    t0  = time.time()
    geo = GDMLGeometry(str(gdml_path))
    print(f"[info] {len(geo.volume_names)} logical volumes  ({time.time()-t0:.1f}s)",
          file=sys.stderr)

    profile_key = detect_profile(geo.volume_names)
    _, prof      = get_profile(profile_key)
    if profile_key:
        print(f"[info] detected {prof['label']}", file=sys.stderr)
    else:
        print("[warn] unknown detector type; using generic patterns", file=sys.stderr)

    ap = anode_pats  or prof["anode_patterns"]
    op = opdet_pats  or prof["opdet_patterns"]
    sp = struct_pats or prof.get("struct_patterns", [])

    entries: list[VolumeEntry] = []

    def _collect(patterns, kind):
        hits = geo.find_positions(patterns, verbose=verbose)
        ok = 0
        for name, pos, rot in hits:
            ext = geo.box_extents(name)
            if ext is None:
                print(f"[warn] no box extents for {kind} {name}, skipping", file=sys.stderr)
                continue
            print(f"  {kind:<8} {name:<40s}"
                  f"  x={pos[0]:8.2f}  y={pos[1]:8.2f}  z={pos[2]:8.2f}"
                  f"  dx={ext[0]:8.2f}  dy={ext[1]:8.2f}  dz={ext[2]:8.2f}")
            entries.append(VolumeEntry(kind=kind, name=name,
                                       pos=pos, ext=ext, rot=rot, meta={}))
            ok += 1
        print(f"[info] {ok}/{len(hits)} {kind} placements", file=sys.stderr)

    _collect(ap, "anode")
    _collect(op, "opdet")
    if sp:
        _collect(sp, "struct")

    return entries, (profile_key or ""), geo


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
                       geo: GDMLGeometry | None = None,
                       prefix: str = "",
                       all_configured: bool = False,
                       select_re: str | None = None,
                       containment: float = 0.99) -> list[VolumeEntry]:
    """Build radio VolumeEntry rows from an already-parsed FCL document."""
    table, labels = select_producers(doc, top_order,
                                     prefix=prefix,
                                     all_configured=all_configured,
                                     select_re=select_re)
    entries: list[VolumeEntry] = []
    for name in labels:
        try:
            rows = producer_volumes(name, table[name],
                                    geo=geo, containment=containment)
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
            meta = {k: r.get(k) for k in
                    ("source", "node", "nuclide", "module_type", "BqPercc",
                     "rate", "surface", "T0", "T1", "distrib_x", "x_eff")
                    if r.get(k) is not None}
            meta["index"] = r["index"]
            entries.append(VolumeEntry(kind="radio", name=name,
                                       pos=pos, ext=ext, rot=rot, meta=meta))
    print(f"[info] {len(entries)} radio volume entries", file=sys.stderr)
    return entries


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
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
    ap.add_argument("--struct-pattern", nargs="+", default=None)
    ap.add_argument("--struct", nargs="+", default=None, metavar="VOLNAME",
                    help="exact structural volume names (anchored, overrides --struct-pattern)")
    # FCL options
    ap.add_argument("--prefix",         default="",
                    help="only producers whose label starts with this")
    ap.add_argument("--select",         default=None,
                    help="regex filter on producer labels")
    ap.add_argument("--all-configured", action="store_true",
                    help="include producers not on any trigger path")
    ap.add_argument("--containment",    type=float, default=0.99,
                    help="distrib_x containment fraction (default 0.99)")
    ap.add_argument("--no-geometry",    action="store_true",
                    help="skip GDML parsing even when a FCL is given")
    ap.add_argument("--gdml",           default=None,
                    help="override GDML path (when input is FCL)")
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
            sys.exit(f"[error] '{args.input}' not found on disk or in ${search_var}")
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
                    sys.exit(f"[error] --gdml '{args.gdml}' not found")
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
    geo: GDMLGeometry | None = None
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
            verbose=args.verbose,
        )
        entries += geo_entries
        if not is_fcl:
            source = str(gdml_path)

    # ---- FCL input: radio volumes (geo available for volume_rand/gen) ----
    if is_fcl:
        entries += _fcl_radio_volumes(
            fcl_doc, fcl_order,
            geo=geo,
            prefix=args.prefix,
            all_configured=args.all_configured,
            select_re=args.select,
            containment=args.containment,
        )

    if not entries:
        sys.exit("[error] nothing extracted — check patterns or use --verbose")

    ir = IRCollection(source=source, detector=detector_key, volumes=entries)
    ir.save(args.out)
    print(f"[info] wrote {len(entries)} volume entries to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()

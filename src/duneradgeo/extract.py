"""
bki-extract  —  extract detector & radiological volumes → JSON IR

\b
Usage:
  bki-extract geometry.gdml --out ir.json
  bki-extract job.fcl       --out ir.json
  bki-extract geometry.gdml --out ir.json \\
      --anode-pattern 'volAnodePlate' --opdet-pattern '^volArapuca'
  bki-extract geometry.gdml --out ir.json --frame volEnclosureTPC

\b
From a FCL file (#includes expanded by fhicl-dump, else via $FHICL_FILE_PATH):
  * RadioGen / Decay0Gen producers → kind="radio" volumes, and each
    producer's full FHiCL table → the IR's top-level "producers";
  * the GDML named in AuxDetGeometry.GDML → the detector volumes
    (anode / opdet / cathode / fieldshaper / struct) and the radio boxes;
  * a distrib_x/y/z formula (needs PyROOT) narrows the box on that axis to
    the effective range of what TF1::GetRandom generates (--npx,
    --containment; see duneradgeo._tf1).

--frame expresses every volume in the frame of one singly placed volume and
records its world placement in the IR: the perl and the ggd
VD GDMLs nest the cryostat differently (a global offset), so their IRs in
the frame of volEnclosureTPC overlay.  --with-passive adds the profile's
passive material (steel shell, insulation, support, shielding, cavern) as
struct volumes.  --gdml-backend root reads the GDML with ROOT
(TGeoManager::Import); scripts/compare_ir.py compares the two backends' IRs.

Non-box solids become boxes: tubes and tori their bounding box, unions one
box per primitive operand (e.g. each straight tube and corner of a
field-shaper ring), subtractions and intersections their first operand.  A
radio volume_rand / volume_gen box is the ROOT bounding box of the volume's
solid (a union: the whole composite), as BaseRadioGen samples it.

Output: a JSON IR, read by bki-plot, bki-distrib and the nb_B notebook.
"""
from __future__ import annotations

import math
import re
import time
from pathlib import Path

import click
import numpy as np

from . import _cli
from ._cli import CONTEXT_SETTINGS, CliError, VariadicOption
from ._draw import aabb
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
    """An extraction step cannot go on; the command reports its message."""


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
    _cli.info(f"parsing {gdml_path.name} ({backend} backend) …")
    t0  = time.time()
    try:
        geo = load_geometry(str(gdml_path), backend=backend, verbose=verbose)
    except RuntimeError as exc:
        raise ExtractError(str(exc)) from exc
    _cli.info(f"{len(geo.volume_names)} logical volumes  ({time.time()-t0:.1f}s)")

    profile_key = detect_profile(geo.volume_names)
    _, prof      = get_profile(profile_key)
    if profile_key:
        _cli.info(f"detected {prof['label']}  (profile {profile_key})")
    else:
        _cli.warn("unknown detector type; using generic patterns")

    ap = anode_pats  or prof["anode_patterns"]
    op = opdet_pats  or prof["opdet_patterns"]
    cp = cathode_pats or prof.get("cathode_patterns", [])
    fp = fieldshaper_pats or prof.get("fieldshaper_patterns", [])
    sp = struct_pats or prof.get("struct_patterns", [])

    entries: list[VolumeEntry] = []
    rows: list[tuple] = []            # (kind, name, centre, size, parts) per summary line

    def _collect(patterns, kind):
        hits = geo.find_positions(patterns, verbose=verbose, with_path=True)
        ok = 0
        copies: dict[str, int] = {}
        for name, pos, rot, path, pvname, pvpath in hits:
            parts = geo.volume_parts(name)
            if not parts:
                _cli.warn(f"no box/tube/torus solid for {kind} {name}, skipping")
                continue
            copy = copies[name] = copies.get(name, -1) + 1
            boxes = []
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
                boxes.append(aabb(wpos, ext, wrot))
            lo = np.min([b[0] for b in boxes], axis=0)
            hi = np.max([b[1] for b in boxes], axis=0)
            c, d    = (lo + hi) / 2, hi - lo
            nparts  = f"  ({len(parts)} parts)" if len(parts) > 1 else ""
            rows.append((kind, name, c, d, nparts.strip()))
            ok += 1
        _cli.info(f"{ok}/{len(hits)} {kind} placements")

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
            _cli.warn("--with-passive: profile has no passive patterns")

    _print_summary(rows)
    return entries, profile_detector(profile_key), geo


def _print_summary(rows: list[tuple]) -> None:
    """The per-volume summary (bounding box centre and size, cm) as a table."""
    if not rows:
        return
    from rich.table import Table
    table = Table(box=None, header_style="bold", pad_edge=False)
    table.add_column("kind")
    table.add_column("name", overflow="fold")
    for h in ("x", "y", "z", "dx", "dy", "dz"):
        table.add_column(h, justify="right")
    table.add_column("parts", style="dim")
    for kind, name, c, d, nparts in rows:
        table.add_row(kind, name, *(f"{v:.2f}" for v in (*c, *d)), nparts)
    _cli.console().print(table)


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
    _cli.info(f"frame {volume} at x={t[0]:.2f} y={t[1]:.2f} z={t[2]:.2f}"
              f"  ({path})")
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
    _cli.info(f"parsing FCL {fcl_path.name} …")
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
            _cli.warn(f"{name}: {exc}")
            continue

        for r in rows:
            lo, hi = r.get("lo"), r.get("hi")
            if lo is None or hi is None or None in lo or None in hi:
                _cli.warn(f"{name}#{r['index']}: unresolved box (needs --gdml), skipping")
                continue
            lo_a, hi_a = np.array(lo, dtype=float), np.array(hi, dtype=float)
            pos = (lo_a + hi_a) / 2
            ext = hi_a - lo_a
            rot = np.eye(3)
            meta = {k: r.get(k) for k in _RADIO_META if r.get(k) is not None}
            meta["index"] = r["index"]
            entries.append(VolumeEntry(kind="radio", name=name,
                                       pos=pos, ext=ext, rot=rot, meta=meta))
    _cli.info(f"{len(entries)} radio volume entries")
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

def main(argv: list | None = None) -> None:
    cli(args=argv, prog_name="bki-extract")


def _patterns(values) -> list[str] | None:
    """The values of a variadic option as a list, None when it was not given."""
    return list(values) or None


class _Args:
    """The options as attributes (args.out, args.frame, ...)."""

    def __init__(self, kw: dict):
        self.__dict__.update(kw)


@click.command(context_settings=CONTEXT_SETTINGS,
               help=__doc__.split("\n", 2)[2].strip())
@click.argument("input", metavar="INPUT")
@click.option("--out", "-o", required=True, metavar="OUT",
              help="output JSON file (IR)")
@click.option("--anode-pattern", cls=VariadicOption, metavar="ANODE_PATTERN...",
              help="regexes of the anode volume names (default: the profile's)")
@click.option("--opdet-pattern", cls=VariadicOption, metavar="OPDET_PATTERN...",
              help="regexes of the optical detector volume names (default: the profile's)")
@click.option("--cathode-pattern", cls=VariadicOption, metavar="CATHODE_PATTERN...",
              help="regexes of the cathode volume names (default: the profile's)")
@click.option("--fieldshaper-pattern", cls=VariadicOption,
              metavar="FIELDSHAPER_PATTERN...",
              help="regexes of the field shaper volume names (default: the profile's)")
@click.option("--struct-pattern", cls=VariadicOption, metavar="STRUCT_PATTERN...",
              help="regexes of the structural volume names (default: the profile's)")
@click.option("--struct", cls=VariadicOption, metavar="VOLNAME...",
              help="exact structural volume names (anchored, overrides --struct-pattern)")
@click.option("--with-passive", is_flag=True,
              help="also extract the profile's passive material "
                   "(steel shell, insulation, support, shielding, cavern) as struct")
@click.option("--frame", default=None, metavar="VOLNAME",
              help="express positions in the frame of this (singly placed) "
                   "volume, e.g. volCryostat or volEnclosureTPC, instead of "
                   "the world")
@click.option("--prefix", default="",
              help="only producers whose label starts with this")
@click.option("--select", "select_re", default=None,
              help="regex filter on producer labels")
@click.option("--all-configured", is_flag=True,
              help="include producers not on any trigger path")
@click.option("--containment", type=float, default=DEFAULT_CONTAINMENT,
              show_default=True,
              help="distrib_x/y/z containment fraction for the effective range")
@click.option("--npx", type=int, default=GENERATOR_NPX, show_default=True,
              help="TF1 Npx of the generator: the effective range is that of "
                   "the distribution TF1::GetRandom draws from an Npx-point "
                   "table; larsimrad BaseRadioGen keeps ROOT's default")
@click.option("--no-geometry", is_flag=True,
              help="skip GDML parsing even when a FCL is given")
@click.option("--gdml", default=None, metavar="GDML",
              help="override GDML path (when input is FCL)")
@click.option("--gdml-backend", type=click.Choice(list(BACKENDS)), default="python",
              show_default=True,
              help="GDML reader: the pure-Python parser, or ROOT's "
                   "TGeoManager::Import (needs ROOT with TGeo/GDML, "
                   "e.g. in the dunesw container)")
@click.option("--fcl-cache", default=default_fcl_cache_dir(), metavar="DIR",
              show_default=True,
              help="directory where expanded FCL files are cached and "
                   "reused across runs")
@click.option("--no-fcl-cache", is_flag=True,
              help="always re-run the FCL expansion, ignoring the cache")
@click.option("--verbose", "-v", is_flag=True)
def cli(**kw) -> None:
    try:
        _run(_Args(kw))
    except (ExtractError, ROOTUnavailable, ValueError) as exc:
        raise CliError(str(exc)) from None


def _run(args) -> None:
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
                    _cli.warn(f"'{name}' not found in $FW_SEARCH_PATH; "
                              "radio volume_rand/gen boxes will be unresolved")
        else:
            gdml_path = inp  # input is already GDML

    # ---- Parse GDML (shared by detector volumes and FCL radio resolution) ----
    geo: GeometryBase | None = None
    if gdml_path is not None:
        struct_pats = (
            [f"^{n}$" for n in args.struct] if args.struct
            else _patterns(args.struct_pattern)
        )
        geo_entries, detector_key, geo = _gdml_volumes(
            gdml_path,
            anode_pats=_patterns(args.anode_pattern),
            opdet_pats=_patterns(args.opdet_pattern),
            struct_pats=struct_pats,
            cathode_pats=_patterns(args.cathode_pattern),
            fieldshaper_pats=_patterns(args.fieldshaper_pattern),
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
            select_re=args.select_re,
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
    _cli.info(f"wrote {len(entries)} volume entries"
              + (f" and {len(producers)} producer configs" if producers else "")
              + f" to {args.out}")


if __name__ == "__main__":
    main()

"""
bki-measure  —  geometric measurements on a JSON IR

Usage:
  bki-measure detector.json                    # every measurement, summaries
  bki-measure detector.json drift-gaps -d      # one measurement, with per-row table
  bki-measure detector.json --json out.json    # also write the results as JSON
  bki-measure --list                           # available measurements

Lengths are in cm, in the world frame (x = drift, z = beam).

The measurements live in the duneradgeo.measurements package, one module
each; adding a module there adds a measurement (see its docstring).
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

from .ir import IRCollection
from .measurements import MEASUREMENTS, Result, run


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _fmt(x) -> str:
    return f"{x:.3f}" if isinstance(x, float) else str(x)


def print_result(res: Result, detail: bool = False) -> None:
    print(f"== {res.name}: {res.title}  [cm]")
    for line in res.summary:
        print(f"  {line}")
    if detail and res.rows:
        cells  = [[_fmt(r.get(c, "")) for c in res.columns] for r in res.rows]
        widths = [max(len(c), *(len(row[i]) for row in cells))
                  for i, c in enumerate(res.columns)]
        print()
        print("  " + "  ".join(c.rjust(w) for c, w in zip(res.columns, widths)))
        for row in cells:
            print("  " + "  ".join(v.rjust(w) for v, w in zip(row, widths)))
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("ir_file", nargs="?", help="JSON IR produced by bki-extract")
    ap.add_argument("names", nargs="*", metavar="MEASUREMENT",
                    help="measurements to run (default: all; see --list)")
    ap.add_argument("--list", action="store_true",
                    help="list the available measurements and exit")
    ap.add_argument("--detail", "-d", action="store_true",
                    help="also print each measurement's per-row table")
    ap.add_argument("--json", metavar="FILE", default=None,
                    help="write the results (summary, values, rows) as JSON")
    args = ap.parse_args()

    if args.list:
        width = max(len(n) for n in MEASUREMENTS)
        for m in MEASUREMENTS.values():
            print(f"{m.name:<{width}}  {m.help}")
        return
    if not args.ir_file:
        ap.error("ir_file is required (unless --list)")

    unknown = [n for n in args.names if n not in MEASUREMENTS]
    if unknown:
        sys.exit(f"[error] unknown measurement(s): {', '.join(unknown)}; "
                 f"available: {', '.join(MEASUREMENTS)}")

    ir_path = Path(args.ir_file)
    if not ir_path.is_file():
        sys.exit(f"[error] '{args.ir_file}' not found")
    ir = IRCollection.load(ir_path)
    print(f"[info] loaded {len(ir.volumes)} volumes from {ir_path.name}", file=sys.stderr)

    results = run(ir, args.names or None)
    for res in results:
        print_result(res, detail=args.detail)

    if args.json:
        Path(args.json).write_text(json.dumps(
            {"source": ir.source, "detector": ir.detector,
             "results": [asdict(r) for r in results]}, indent=2))
        print(f"[info] results written to {args.json}", file=sys.stderr)


if __name__ == "__main__":
    main()

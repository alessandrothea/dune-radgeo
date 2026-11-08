"""
bki-measure  —  geometric measurements on a JSON IR

Usage:
  bki-measure detector.json                    # every measurement, summaries
  bki-measure detector.json drift-gaps -d      # one measurement, with per-row table
  bki-measure detector.json --json out.json    # also write the results as JSON
  bki-measure --list                           # available measurements

Lengths are in cm, in the IR frame (x = drift, z = beam).

The measurements live in the duneradgeo.measurements package, one module
each; adding a module there adds a measurement (see its docstring).
"""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import click

from ._cli import CONTEXT_SETTINGS, CliError, console, info
from .ir import IRCollection
from .measurements import MEASUREMENTS, Result, run


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _fmt(x) -> str:
    return f"{x:.3f}" if isinstance(x, float) else str(x)


def render_list(con) -> None:
    """The available measurements: name and one-line help."""
    from rich.table import Table
    t = Table(box=None, pad_edge=False, header_style="bold")
    t.add_column("measurement", no_wrap=True)
    t.add_column("help")
    for m in MEASUREMENTS.values():
        t.add_row(m.name, m.help)
    con.print(t)


def render_result(con, res: Result, detail: bool = False) -> None:
    """A panel with the summary lines (spacing kept) and, with `detail`, the
    per-row table below it."""
    from rich.box import SIMPLE_HEAD
    from rich.console import Group
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text

    parts = [Text("\n".join(res.summary) if res.summary else "(no summary)",
                  no_wrap=True, overflow="ignore")]
    if detail and res.rows:
        t = Table(box=SIMPLE_HEAD, header_style="bold", pad_edge=False)
        for c in res.columns:
            t.add_column(c, justify="right", no_wrap=True)
        for r in res.rows:
            t.add_row(*[Text(_fmt(r.get(c, ""))) for c in res.columns])
        parts.append(t)
    con.print(Panel(Group(*parts), title=Text(f"{res.name}: {res.title}  [cm]", style="bold"),
                    title_align="left", expand=False))


# ---------------------------------------------------------------------------
# command line (click)
# ---------------------------------------------------------------------------

_HELP = """Geometric measurements on a JSON IR.

\b
Usage:
  bki-measure detector.json                    # every measurement, summaries
  bki-measure detector.json drift-gaps -d      # one measurement, with per-row table
  bki-measure detector.json --json out.json    # also write the results as JSON
  bki-measure --list                           # available measurements

Lengths are in cm, in the IR frame (x = drift, z = beam).

The measurements live in the duneradgeo.measurements package, one module
each; adding a module there adds a measurement (see its docstring)."""


@click.command(context_settings=CONTEXT_SETTINGS, help=_HELP)
@click.argument("ir_file", required=False, metavar="[IR_FILE]",
                type=click.Path(dir_okay=False))
@click.argument("names", nargs=-1, metavar="[MEASUREMENT]...")
@click.option("--list", "list_", is_flag=True,
              help="List the available measurements and exit.")
@click.option("--detail", "-d", is_flag=True,
              help="Also print each measurement's per-row table.")
@click.option("--json", "json_file", type=click.Path(dir_okay=False), default=None,
              metavar="FILE", help="Write the results (summary, values, rows) as JSON.")
def cli(ir_file, names: tuple, list_: bool, detail: bool, json_file) -> None:
    # IR_FILE is a JSON IR produced by bki-extract; MEASUREMENT: the
    # measurements to run (default: all; see --list)
    if list_:
        render_list(console())
        return
    if not ir_file:
        raise click.UsageError("IR_FILE is required (unless --list)")

    unknown = [n for n in names if n not in MEASUREMENTS]
    if unknown:
        raise CliError(f"unknown measurement(s): {', '.join(unknown)}; "
                       f"available: {', '.join(MEASUREMENTS)}")

    ir_path = Path(ir_file)
    if not ir_path.is_file():
        raise CliError(f"'{ir_file}' not found")
    try:
        ir = IRCollection.load(ir_path)
        info(f"loaded {len(ir.volumes)} volumes from {ir_path.name}")
        results = run(ir, list(names) or None)
    except ValueError as exc:
        raise CliError(str(exc)) from None

    con = console()
    for res in results:
        render_result(con, res, detail=detail)

    if json_file:
        Path(json_file).write_text(json.dumps(
            {"source": ir.source, "detector": ir.detector,
             "results": [asdict(r) for r in results]}, indent=2))
        info(f"results written to {json_file}")


def main(argv: list | None = None) -> None:
    cli(args=argv, prog_name="bki-measure")


if __name__ == "__main__":
    main()

"""
radgeo  —  one command line for the duneradgeo tools.

  radgeo extract      GDML / FHiCL  ->  JSON IR
  radgeo plot         render a JSON IR as a PDF / PNG
  radgeo measure      geometric measurements on a JSON IR
  radgeo profiles     draw the distrib_x/y/z position profiles of radio producers
  radgeo backgrounds  inspect the radiological backgrounds of a JSON IR (alias: bkg)

``radgeo COMMAND --help`` shows the options of each command.
"""
from __future__ import annotations

import click

from . import bkg, distrib, extract, measure, plot
from ._cli import CONTEXT_SETTINGS

ALIASES = {"bkg": "backgrounds"}


class _RadgeoGroup(click.Group):
    """A group that also accepts the aliases of ALIASES (not listed in --help)."""

    def get_command(self, ctx, cmd_name):
        return super().get_command(ctx, ALIASES.get(cmd_name, cmd_name))

    def resolve_command(self, ctx, args):
        name, cmd, rest = super().resolve_command(ctx, args)
        return (cmd.name if cmd is not None else name), cmd, rest


@click.group(cls=_RadgeoGroup, context_settings=CONTEXT_SETTINGS,
             help="Inspect DUNE background sources: detector geometry (GDML) and "
                  "radiological volumes (FHiCL).\n\n"
                  "\b\nExamples:\n"
                  "  radgeo extract --help\n"
                  "  radgeo plot detector.json --out detector.pdf\n"
                  "  radgeo backgrounds ir.json list")
def cli() -> None:
    pass


_SUBCOMMANDS = (
    ("extract", extract.cli, "Extract a GDML geometry and FHiCL parameters into a JSON IR."),
    ("plot", plot.cli, "Render a JSON IR as a matplotlib PDF or PNG."),
    ("measure", measure.cli, "Geometric measurements on a JSON IR."),
    ("profiles", distrib.cli, "Draw the distrib_x/y/z position profiles of radio producers."),
    ("backgrounds", bkg.cli, "Inspect the radiological backgrounds of a JSON IR."),
)
for _name, _cmd, _short in _SUBCOMMANDS:
    _cmd.name = _name            # the name shown in 'Usage: radgeo <name> ...'
    _cmd.short_help = _short
    cli.add_command(_cmd, _name)


def main(argv: list | None = None) -> None:
    cli.main(args=argv, prog_name="radgeo")

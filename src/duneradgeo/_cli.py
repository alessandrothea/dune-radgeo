"""
Shared pieces of the radgeo command lines, built on click and rich.

- ``CONTEXT_SETTINGS``: -h as well as --help.
- ``CliError``: a message for the user, printed as ``[error] ...`` with exit
  code 1 (click's own usage errors exit with 2).
- ``info`` / ``warn`` / ``error``: the ``[info] ...`` / ``[warn] ...`` /
  ``[error] ...`` lines on stderr, coloured when stderr is a terminal.
  The text is exactly what was printed before, so logs can still be grepped.
- ``console``: a rich Console for tables on stdout, 160 columns wide when the
  output is not a terminal (unless COLUMNS is set) so piped tables stay readable.
- ``VariadicOption``: an option that takes one or more values up to the next
  option, as argparse's ``nargs="+"`` did: ``--name A B C``.
"""
from __future__ import annotations

import os

import click

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}
PIPE_WIDTH = 160            # table width when stdout is not a terminal

_STYLES = {"info": "dim", "warn": "yellow", "error": "bold red"}


class CliError(click.ClickException):
    """A message for the user: printed as '[error] ...', exit code 1."""

    def show(self, file=None) -> None:
        error(self.message)


def _rich_console():
    try:
        from rich.console import Console
    except ImportError:
        raise CliError("the 'rich' package is needed (pip install rich)") from None
    return Console


def _err_console():
    return _rich_console()(stderr=True, highlight=False, markup=False, soft_wrap=True)


def _say(level: str, msg: str) -> None:
    try:
        console = _err_console()
    except CliError:                       # no rich: plain text, so errors still show
        click.echo(f"[{level}] {msg}", err=True)
        return
    console.print(f"[{level}] {msg}", style=_STYLES[level])


def info(msg: str) -> None:
    _say("info", msg)


def warn(msg: str) -> None:
    _say("warn", msg)


def error(msg: str) -> None:
    _say("error", msg)


def console():
    """A rich Console on stdout for tables (wide when piped)."""
    Console = _rich_console()
    c = Console(highlight=False)
    if not c.is_terminal and "COLUMNS" not in os.environ:
        c = Console(width=PIPE_WIDTH, highlight=False)
    return c


class VariadicOption(click.Option):
    """An option taking one or more values, up to the next option.

    ``--name A B C --kind radio`` gives ``name == ("A", "B", "C")``.  A token
    that starts with '-' ends the values unless it is a number (so ``--x -5``
    style values stay values); write ``--opt=VALUE`` for any other value that
    starts with '-' (bare words after it are still taken as further values,
    so put another option or the positional argument before them).  Declare it as
    ``@click.option("--name", cls=VariadicOption, metavar="PATTERN...")``.
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("multiple", True)
        super().__init__(*args, **kwargs)

    def add_to_parser(self, parser, ctx):
        super().add_to_parser(parser, ctx)
        names = [n for n in self.opts]                       # '--name', '-n'
        found = None
        for n in names:
            found = (parser._long_opt.get(n) if n.startswith("--")
                     else parser._short_opt.get(n)) or found
        if found is None:                                    # pragma: no cover
            return
        previous = found.process

        def process(value, state):
            # `value` is the first value; take the following ones while they
            # are not options
            while state.rargs and not _looks_like_option(state.rargs[0]):
                previous(value, state)
                value = state.rargs.pop(0)
            previous(value, state)

        found.process = process


def _looks_like_option(token: str) -> bool:
    if not token.startswith("-") or token == "-":
        return False
    try:
        float(token)
    except ValueError:
        return True
    return False

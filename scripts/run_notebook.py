#!/usr/bin/env python3
"""
Run nb_B.ipynb non-interactively against a given JSON IR file (IR_FILE: the
JSON IR produced by bki-extract).

Examples:

\b
    run_notebook.py detector.json
    run_notebook.py detector.json --out executed.ipynb
    run_notebook.py detector.json --timeout 300

The executed notebook is written to <ir_stem>_nb_B.ipynb in the current
directory unless --out is given.  Raises SystemExit(1) on cell errors.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import click

from duneradgeo._cli import CONTEXT_SETTINGS, CliError, info


@click.command(context_settings=CONTEXT_SETTINGS, help=__doc__)
@click.argument("ir_file")
@click.option("--out", "-o", default=None,
              help="path for the executed output notebook (default: <ir_stem>_nb_B.ipynb in CWD)")
@click.option("--timeout", type=int, default=120, show_default=True,
              help="per-cell execution timeout in seconds")
def cli(ir_file, out, timeout):
    ir_path = Path(ir_file).resolve()
    if not ir_path.is_file():
        raise CliError(f"IR file not found: {ir_path}")

    nb_path = Path(__file__).parent.parent / "notebooks" / "nb_B.ipynb"
    if not nb_path.is_file():
        raise CliError(f"notebook not found: {nb_path}")

    out_path = Path(out) if out else Path(f"{ir_path.stem}_nb_B.ipynb")

    env = os.environ.copy()
    env["BKI_IR_FILE"] = str(ir_path)

    cmd = [
        sys.executable, "-m", "jupyter", "nbconvert",
        "--to", "notebook",
        "--execute",
        f"--ExecutePreprocessor.timeout={timeout}",
        "--output", str(out_path.resolve()),
        str(nb_path),
    ]

    info(f"executing {nb_path.name} with IR={ir_path.name} …")
    result = subprocess.run(cmd, env=env)

    if result.returncode != 0:
        raise CliError(f"notebook execution failed (exit {result.returncode})")

    info(f"executed notebook written to {out_path}")


def main(argv: list | None = None) -> None:
    cli.main(args=argv, prog_name="run_notebook.py")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Run nb_B.ipynb non-interactively against a given JSON IR file.

Usage:
    run_notebook.py detector.json
    run_notebook.py detector.json --out executed.ipynb
    run_notebook.py detector.json --timeout 300

The executed notebook is written to <ir_stem>_nb_B.ipynb in the current
directory unless --out is given.  Raises SystemExit(1) on cell errors.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("ir_file", help="JSON IR produced by bki-extract")
    ap.add_argument("--out", "-o", default=None,
                    help="path for the executed output notebook "
                         "(default: <ir_stem>_nb_B.ipynb in CWD)")
    ap.add_argument("--timeout", type=int, default=120,
                    help="per-cell execution timeout in seconds (default: 120)")
    args = ap.parse_args()

    ir_path = Path(args.ir_file).resolve()
    if not ir_path.is_file():
        sys.exit(f"[error] IR file not found: {ir_path}")

    nb_path = Path(__file__).parent.parent / "notebooks" / "nb_B.ipynb"
    if not nb_path.is_file():
        sys.exit(f"[error] notebook not found: {nb_path}")

    out_path = Path(args.out) if args.out else Path(f"{ir_path.stem}_nb_B.ipynb")

    env = os.environ.copy()
    env["BKI_IR_FILE"] = str(ir_path)

    cmd = [
        sys.executable, "-m", "jupyter", "nbconvert",
        "--to", "notebook",
        "--execute",
        f"--ExecutePreprocessor.timeout={args.timeout}",
        "--output", str(out_path.resolve()),
        str(nb_path),
    ]

    print(f"[info] executing {nb_path.name} with IR={ir_path.name} …", flush=True)
    result = subprocess.run(cmd, env=env)

    if result.returncode != 0:
        sys.exit(f"[error] notebook execution failed (exit {result.returncode})")

    print(f"[info] executed notebook written to {out_path}")


if __name__ == "__main__":
    main()

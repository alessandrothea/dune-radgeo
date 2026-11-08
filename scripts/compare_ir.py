#!/usr/bin/env python3
"""
Compare two JSON IR files volume by volume (regression check): REF is the
reference JSON IR, NEW the one to check against it.

Examples:

\b
    compare_ir.py ref.json new.json
    compare_ir.py ref.json new.json --tol 0.1 --ignore-rot '^volTPC'
    compare_ir.py ref.json new.json --only-common -v

Entries are grouped by (kind, name, part) and each reference entry is paired
with the nearest entry of the same group in the other IR.  A group passes
when both IRs have the same number of entries, the pairing is one-to-one,
and every pair agrees in centre (--tol), extents (--tol) and rotation
(--rot-tol).  Extract both IRs with the same --frame (e.g. volEnclosureTPC)
when the GDMLs nest the detector differently.

Exits with status 1 if any group fails.
"""
from __future__ import annotations

import re
from collections import Counter, defaultdict
from pathlib import Path

import click
import numpy as np
from rich.text import Text

from duneradgeo._cli import CONTEXT_SETTINGS, CliError, VariadicOption, console, warn
from duneradgeo.ir import IRCollection


def _groups(ir: IRCollection, kinds: set[str] | None) -> dict[tuple, list]:
    g: dict[tuple, list] = defaultdict(list)
    for v in ir.volumes:
        if kinds and v.kind not in kinds:
            continue
        g[(v.kind, v.name, v.meta.get("part", ""))].append(v)
    return g


def _compare(a: list, b: list, tol: float, rot_tol: float, check_rot: bool) -> tuple[bool, str]:
    pa = np.array([v.pos for v in a])
    pb = np.array([v.pos for v in b])
    d2 = ((pa[:, None, :] - pb[None, :, :]) ** 2).sum(axis=2)
    j  = d2.argmin(axis=1)
    dpos = np.sqrt(d2[np.arange(len(a)), j])
    dext = np.array([np.abs(a[i].ext - b[j[i]].ext).max() for i in range(len(a))])
    drot = np.array([np.abs(a[i].rot - b[j[i]].rot).max() for i in range(len(a))])
    n_pos = int((dpos > tol).sum())
    n_ext = int((dext > tol).sum())
    n_rot = int((drot > rot_tol).sum()) if check_rot else 0
    unique = len(set(j.tolist())) == len(a)
    ok = len(a) == len(b) and unique and not (n_pos or n_ext or n_rot)
    msg = (f"n={len(a)}/{len(b)}  max|dpos|={dpos.max():.3f}  max|dext|={dext.max():.3f}"
           f"  max|drot|={drot.max():.2g}")
    bad = []
    if len(a) != len(b):
        bad.append("count")
    if not unique:
        bad.append("pairing not one-to-one")
    if n_pos:
        shifts = Counter(tuple(np.round(pb[j[i]] - pa[i], 2)) for i in np.where(dpos > tol)[0])
        top = ", ".join(f"{n}×({s[0]:g},{s[1]:g},{s[2]:g})" for s, n in shifts.most_common(3))
        bad.append(f"{n_pos} pos [{top}]")
    if n_ext:
        bad.append(f"{n_ext} ext")
    if n_rot:
        bad.append(f"{n_rot} rot")
    if bad:
        msg += "  ✗ " + "; ".join(bad)
    elif not check_rot and (drot > rot_tol).any():
        msg += f"  ({int((drot > rot_tol).sum())} rot, ignored)"
    return ok, msg


def _status(ok: bool | None) -> Text:
    """'ok  ' green, 'FAIL' bold red, 'skip' dim (words unchanged when piped)."""
    if ok is None:
        return Text("skip", style="dim")
    return Text("ok  ", style="green") if ok else Text("FAIL", style="bold red")


def _row(cons, status: Text, label: str, msg: str, width: int = 50) -> None:
    line = Text.assemble(status, "  ", f"{label:<{width}} ", msg)
    cons.print(line, soft_wrap=True)


def _load(path: str) -> IRCollection:
    if not Path(path).is_file():
        raise CliError(f"'{path}' not found")
    try:
        return IRCollection.load(path)
    except (OSError, ValueError, KeyError) as exc:
        raise CliError(f"cannot read '{path}': {exc}") from None


@click.command(context_settings=CONTEXT_SETTINGS, help=__doc__)
@click.argument("ref")
@click.argument("new")
@click.option("--tol", type=float, default=0.01, show_default=True,
              help="position / extent tolerance in cm")
@click.option("--rot-tol", type=float, default=1e-6, show_default=True,
              help="max abs difference of rotation-matrix elements")
@click.option("--kinds", cls=VariadicOption, metavar="KIND...", default=None,
              help="only compare these kinds (anode opdet cathode fieldshaper struct radio)")
@click.option("--ignore-rot", cls=VariadicOption, metavar="NAME_RE...", default=(),
              help="do not fail on rotation differences of these volume names")
@click.option("--only-common", is_flag=True,
              help="do not fail on groups present in only one IR")
@click.option("--verbose", "-v", is_flag=True, help="also list the groups that pass")
def cli(ref, new, tol, rot_tol, kinds, ignore_rot, only_common, verbose):
    ra, rb = _load(ref), _load(new)
    cons = console()
    fa = (ra.frame or {}).get("volume", "world")
    fb = (rb.frame or {}).get("volume", "world")
    cons.print(Text(f"[info] ref: {len(ra.volumes)} entries, frame {fa}   "
                    f"new: {len(rb.volumes)} entries, frame {fb}", style="dim"), soft_wrap=True)
    if fa != fb:
        warn(f"frames differ ({fa} vs {fb}); positions are not comparable")

    kinds = set(kinds) if kinds else None
    ga, gb = _groups(ra, kinds), _groups(rb, kinds)
    try:
        ignore = [re.compile(p) for p in ignore_rot]
    except re.error as exc:
        raise CliError(f"bad --ignore-rot pattern: {exc}") from None

    failed = 0
    for key in sorted(set(ga) | set(gb)):
        kind, name, part = key
        label = f"{kind:<11} {name}{'/' + part if part else ''}"
        if key not in ga or key not in gb:
            side = "ref" if key in ga else "new"
            _row(cons, _status(None if only_common else False), label, f"only in {side}")
            failed += not only_common
            continue
        check_rot = not any(p.search(name) for p in ignore)
        ok, msg = _compare(ga[key], gb[key], tol, rot_tol, check_rot)
        failed += not ok
        if not ok or verbose:
            _row(cons, _status(ok), label, msg)

    n = len(set(ga) | set(gb))
    cons.print(Text(f"[info] {n - failed}/{n} groups agree", style="dim"), soft_wrap=True)
    if failed:
        raise SystemExit(1)


def main(argv: list | None = None) -> None:
    cli.main(args=argv, prog_name="compare_ir.py")


if __name__ == "__main__":
    main()

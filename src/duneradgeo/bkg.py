"""
radgeo backgrounds  —  inspect the radiological backgrounds of a JSON IR

A background is one radiological producer: the radio volumes of the IR that
share its name (its sampling boxes) and its FHiCL parameters (IR "producers").

Usage:
  radgeo backgrounds ir.json list                         # every background
  radgeo backgrounds ir.json show Po218 Cathode           # parameters and volumes (regex on names)
  radgeo backgrounds ir.json select nuclide=Po218 'BqPercc<0.01'
  radgeo backgrounds ir.json region --slice x=-328        # active on a plane
  radgeo backgrounds ir.json region --x -330 -300 --z 0 500
  radgeo backgrounds ir.json region --box -330 -300 -700 700 0 2000 --volumes
  radgeo backgrounds ir.json suggest                      # other selections worth trying

select EXPR:  FIELD=REGEX (re.search on any value of the field), FIELD!=REGEX
(no value matches), FIELD<N, <=, >, >= (any numeric value), or a bare FIELD
(the background has it); several EXPRs must all hold.  FIELDs are the FHiCL
parameters of the producer plus the metadata of its volumes (nuclide, surface,
node, source, module_type, BqPercc, T0, ...), `name` and `nvolumes`;
`suggest` lists them.  Every command takes --json for machine-readable output.

Coordinates are in the frame stored in the IR (world unless it was extracted
with --frame).  A background is active in a region when the axis-aligned box
of at least one of its volumes overlaps it (touching counts).

Needs `rich` for the tables (pip install rich), and nothing else but numpy.
"""
from __future__ import annotations

import difflib
import json
import re
import shlex
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import click
import numpy as np

from ._cli import CONTEXT_SETTINGS, CliError, console
from ._draw import aabb
from .ir import IRCollection

AXES = "xyz"
_MAX_CHOICES = 12            # a field with more distinct values is not suggested
_NOT_SUGGESTED = {"name", "nvolumes", "index"}   # plus every distrib_* field
_BOX_FIELD = re.compile(r"^[XYZ][01]$")        # sampling boxes: see `region`


# ---------------------------------------------------------------------------
# model
# ---------------------------------------------------------------------------

def _flat(v):
    """The scalars inside a value: a list's items, a dict's values."""
    if isinstance(v, dict):
        for x in v.values():
            yield from _flat(x)
    elif isinstance(v, (list, tuple)):
        for x in v:
            yield from _flat(x)
    else:
        yield v


def _s(v) -> str:
    """A value as text, the form that regexes and tables see."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def _num(v):
    """float(v) for a number (not a bool); None otherwise."""
    if isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v)
    return None


@dataclass
class Background:
    """One radiological producer: its FHiCL parameters and its volumes."""
    name:    str
    params:  dict = field(default_factory=dict)
    volumes: list = field(default_factory=list)       # list[VolumeEntry]

    @cached_property
    def fields(self) -> dict:
        """field -> distinct scalar values (producer parameters, then the
        metadata of the volumes), plus `name` and `nvolumes`."""
        out: dict = {"name": [self.name], "nvolumes": [len(self.volumes)]}

        def add(key, value):
            lst = out.setdefault(key, [])
            for x in _flat(value):
                if x is not None and x not in lst:
                    lst.append(x)

        for k, v in self.params.items():
            add(k, v)
        for vol in self.volumes:
            for k, v in vol.meta.items():
                add(k, v)
        return out

    def texts(self, key: str) -> list:
        return [_s(x) for x in self.fields.get(key, [])]

    @cached_property
    def boxes(self) -> list:
        """(lo, hi) of every volume, in the IR frame."""
        return [aabb(v.pos, v.ext, v.rot) for v in self.volumes]


def backgrounds(ir: IRCollection) -> list:
    """The backgrounds of an IR, sorted by name: one per distinct radio volume
    name.  Producers without volumes are not listed."""
    by_name: dict = defaultdict(list)
    for v in ir.volumes:
        if v.kind == "radio":
            by_name[v.name].append(v)
    return [Background(n, ir.producers.get(n, {}), vols)
            for n, vols in sorted(by_name.items())]


# ---------------------------------------------------------------------------
# select
# ---------------------------------------------------------------------------

_EXPR = re.compile(r"^\s*([A-Za-z_][\w.]*)\s*(?:(!=|>=|<=|=|>|<)\s*(.*))?$", re.S)


@dataclass
class Expr:
    field: str
    op:    str | None        # None: the field exists
    value: str = ""
    regex: object = None     # compiled value, for = and !=
    number: float = 0.0      # parsed value, for the numeric operators

    def __str__(self) -> str:
        return self.field if self.op is None else f"{self.field}{self.op}{self.value}"


def parse_expr(text: str) -> Expr:
    """'FIELD', 'FIELD=REGEX', 'FIELD!=REGEX' or 'FIELD<N' (<=, >, >=).
    Raises ValueError, with a message for the user, if it is malformed."""
    m = _EXPR.match(text)
    if not m:
        raise ValueError(f"bad expression {text!r} (expected FIELD, FIELD=REGEX, "
                         f"FIELD!=REGEX or FIELD<N, <=, >, >=)")
    key, op, value = m.group(1), m.group(2), m.group(3) or ""
    e = Expr(key, op, value)
    if op in ("=", "!="):
        try:
            e.regex = re.compile(value)
        except re.error as exc:
            raise ValueError(f"bad regex {value!r} in {text!r}: {exc}") from None
    elif op is not None:
        try:
            e.number = float(value)
        except ValueError:
            raise ValueError(f"{text!r}: {op} needs a number, got {value!r}") from None
    return e


def known_fields(bgs: list) -> list:
    return sorted({k for b in bgs for k in b.fields})


def check_fields(exprs: list, bgs: list) -> None:
    """ValueError naming the unknown field of an expression, with the closest
    known fields."""
    known = known_fields(bgs)
    for e in exprs:
        if e.field not in known:
            close = difflib.get_close_matches(e.field, known, n=3)
            hint = (f"; did you mean {', '.join(close)}?" if close
                    else f"; known fields: {', '.join(known)}")
            raise ValueError(f"unknown field '{e.field}'{hint}")


def matches(bg: Background, e: Expr) -> bool:
    values = bg.fields.get(e.field)
    if e.op is None:
        return values is not None
    if e.op == "!=":                                 # also true when absent
        return not any(e.regex.search(_s(x)) for x in values or [])
    if values is None:
        return False
    if e.op == "=":
        return any(e.regex.search(_s(x)) for x in values)
    cmp = {"<": lambda a, b: a < b, "<=": lambda a, b: a <= b,
           ">": lambda a, b: a > b, ">=": lambda a, b: a >= b}[e.op]
    return any(cmp(n, e.number) for n in map(_num, values) if n is not None)


def select(bgs: list, exprs: list) -> list:
    return [b for b in bgs if all(matches(b, e) for e in exprs)]


# ---------------------------------------------------------------------------
# region
# ---------------------------------------------------------------------------

@dataclass
class Region:
    lo: np.ndarray
    hi: np.ndarray

    def text(self) -> str:
        parts = []
        for a, lo, hi in zip(AXES, self.lo, self.hi):
            if np.isinf(lo) and np.isinf(hi):
                continue
            parts.append(f"{a} = {lo:g}" if lo == hi else
                         f"{a} in [{'-inf' if np.isinf(lo) else f'{lo:g}'}, "
                         f"{'inf' if np.isinf(hi) else f'{hi:g}'}]")
        return ", ".join(parts) or "everywhere"

    def overlaps(self, box: tuple) -> bool:
        lo, hi = box
        return bool(np.all(hi >= self.lo) and np.all(lo <= self.hi))


def make_region(box=None, ranges=None, slices=None, thickness: float = 0.0) -> Region:
    """The intersection of: `box` = (x0, x1, y0, y1, z0, z1); `ranges` =
    {axis: (lo, hi)}; `slices` = {axis: value}, a plane of the given
    `thickness` (0: a plane).  ValueError if empty or malformed."""
    lo, hi = np.full(3, -np.inf), np.full(3, np.inf)

    def limit(i, a, b):
        if a > b:
            raise ValueError(f"empty range for {AXES[i]}: {a:g} > {b:g}")
        lo[i], hi[i] = max(lo[i], a), min(hi[i], b)
        if lo[i] > hi[i]:
            raise ValueError(f"the constraints on {AXES[i]} do not intersect")

    if box:
        for i in range(3):
            limit(i, box[2 * i], box[2 * i + 1])
    for a, (x0, x1) in (ranges or {}).items():
        limit(AXES.index(a), x0, x1)
    for a, v in (slices or {}).items():
        limit(AXES.index(a), v - thickness / 2, v + thickness / 2)
    return Region(lo, hi)


def in_region(bgs: list, region: Region) -> list:
    """(background, [indices of its volumes overlapping the region]) for every
    background with at least one."""
    out = []
    for b in bgs:
        hits = [i for i, box in enumerate(b.boxes) if region.overlaps(box)]
        if hits:
            out.append((b, hits))
    return out


# ---------------------------------------------------------------------------
# summaries (data for both the tables and --json)
# ---------------------------------------------------------------------------

def _range_text(box: tuple, axis: int) -> str:
    return f"[{box[0][axis]:.2f}, {box[1][axis]:.2f}]"


def surface_counts(b: Background) -> list:
    """[(surface, n volumes)] in order of appearance."""
    c = Counter(v.meta.get("surface", "?") for v in b.volumes)
    return list(c.items())


def summary(b: Background) -> dict:
    """The row of `list`: what a background is, in a few columns."""
    activity = b.texts("BqPercc") or b.texts("rate")
    return {
        "name": b.name,
        "module_type": ", ".join(b.texts("module_type")),
        "nuclide": ", ".join(b.texts("nuclide")),
        "activity": ", ".join(activity),
        "activity_field": "BqPercc" if b.texts("BqPercc") else ("rate" if activity else ""),
        "material": ", ".join(b.texts("material") or b.texts("Material")),
        "sampling": ", ".join(b.texts("source")),
        "volumes": len(b.volumes),
        "surfaces": [{"surface": s, "volumes": n} for s, n in surface_counts(b)],
        "distrib": [a for a in AXES if f"distrib_{a}" in b.fields],
    }


def detail(b: Background, limit: int = 25) -> dict:
    """What `show` prints: the parameters, the surfaces and the first `limit`
    volumes (0: all)."""
    vols = []
    for v, box in zip(b.volumes, b.boxes):
        vols.append({"index": v.meta.get("index"), "node": v.meta.get("node"),
                     "surface": v.meta.get("surface"),
                     "lo": [round(float(x), 3) for x in box[0]],
                     "hi": [round(float(x), 3) for x in box[1]]})
    shown = vols if not limit else vols[:limit]
    return {"name": b.name, "params": b.params, "surfaces": summary(b)["surfaces"],
            "volumes": shown, "n_volumes": len(vols)}


def region_rows(hits: list) -> list:
    """The rows of `region`: per background, how many of its volumes overlap
    and the box enclosing those."""
    rows = []
    for b, idx in hits:
        lo = np.min([b.boxes[i][0] for i in idx], axis=0)
        hi = np.max([b.boxes[i][1] for i in idx], axis=0)
        rows.append({"name": b.name, "overlapping": len(idx), "volumes": len(b.volumes),
                     "lo": [round(float(x), 3) for x in lo],
                     "hi": [round(float(x), 3) for x in hi],
                     "indices": idx})
    return rows


# ---------------------------------------------------------------------------
# suggestions
# ---------------------------------------------------------------------------

_NAME_TOKEN = re.compile(r"[A-Z][a-z]+|[A-Z]+(?![a-z])|\d+(?:x\d+)+")
_NAME_GLUE = {"Gen", "In", "At", "From"}      # words that join name parts


def _plain_or_anchored(value: str, others: list) -> str:
    """A regex for `value` that does not also match the other values."""
    plain = re.escape(value)
    if any(o != value and re.search(plain, o) for o in others):
        return f"^{plain}$"
    return plain


def _command(*words: str) -> str:
    """The arguments of a radgeo backgrounds command, shell-quoted (the IR file goes first)."""
    return " ".join(shlex.quote(w) for w in words)


def _round_threshold(values: list) -> float:
    """A round number near the median of `values`."""
    med = float(np.median(values))
    return float(f"{med:.2g}")


def suggestions(bgs: list) -> dict:
    """Selections and regions that split the backgrounds: the values of the
    fields that two or more backgrounds share; thresholds for the numeric
    fields; the CamelCase words of the names; the planes the volumes lie on."""
    fields_out, numeric_out = [], []
    for key in known_fields(bgs):
        if key in _NOT_SUGGESTED or key.startswith("distrib_") or _BOX_FIELD.match(key):
            continue
        raw = [x for b in bgs for x in b.fields.get(key, [])]
        if raw and all(_num(x) is not None for x in raw):
            per_bg = [max(_num(x) for x in b.fields[key]) for b in bgs if key in b.fields]
            if len(set(per_bg)) >= 2:
                thr = _round_threshold(per_bg)
                expr = f"{key}>{thr:g}"
                numeric_out.append({
                    "field": key, "backgrounds": len(per_bg),
                    "min": min(per_bg), "max": max(per_bg),
                    "matching": len(select(bgs, [parse_expr(expr)])),
                    "command": _command("select", expr)})
            continue
        count: Counter = Counter()
        for b in bgs:
            for t in b.texts(key):
                count[t] += 1
        if len(count) < 2 or len(count) > 4 * _MAX_CHOICES:
            continue
        values = []
        for t, c in sorted(count.items(), key=lambda kv: (-kv[1], kv[0])):
            if c < 2:
                continue
            expr = f"{key}={_plain_or_anchored(t, list(count))}"
            n = len(select(bgs, [parse_expr(expr)]))
            values.append({"value": t, "backgrounds": n,
                           "command": _command("select", expr)})
        if values and len(values) <= _MAX_CHOICES and any(
                0 < v["backgrounds"] < len(bgs) for v in values):
            fields_out.append({"field": key, "values": values})

    tokens: Counter = Counter()
    for b in bgs:
        tokens.update(set(_NAME_TOKEN.findall(b.name)))
    names = []
    for tok, _ in sorted(tokens.items(), key=lambda kv: (-kv[1], kv[0])):
        if len(tok) < 3 or tok in _NAME_GLUE:
            continue
        n = len(select(bgs, [parse_expr(f"name={re.escape(tok)}")]))
        if 1 < n < len(bgs):
            names.append({"token": tok, "backgrounds": n,
                          "command": _command("select", f"name={re.escape(tok)}")})
    names.sort(key=lambda d: (-d["backgrounds"], d["token"]))

    planes: dict = {}
    for b in bgs:
        for s in b.texts("surface"):
            m = re.fullmatch(r"plane_([xyz])@(.+)", s)
            if m:
                planes[(m.group(1), float(m.group(2)))] = s
    regions = []
    for (a, c), _ in sorted(planes.items()):
        hits = in_region(bgs, make_region(slices={a: c}))
        label = planes[(a, c)]
        on = [(b, dict(surface_counts(b))[label]) for b in bgs
              if label in dict(surface_counts(b))]
        regions.append({"slice": f"{a}={c:g}", "on_plane": len(on),
                        "volumes_on_plane": sum(n for _, n in on),
                        "crossing": len(hits),
                        "command": _command("region", "--slice", f"{a}={c:g}")})
    allvols = [box for b in bgs for box in b.boxes]
    extent = {}
    if allvols:
        lo = np.min([b[0] for b in allvols], axis=0)
        hi = np.max([b[1] for b in allvols], axis=0)
        extent = {a: [round(float(lo[i]), 3), round(float(hi[i]), 3)]
                  for i, a in enumerate(AXES)}
    return {"fields": fields_out, "numeric": numeric_out, "names": names,
            "regions": regions, "extent": extent}


# ---------------------------------------------------------------------------
# rendering (rich)
# ---------------------------------------------------------------------------

def _surfaces_text(items: list) -> str:
    return ", ".join(f"{d['surface']}" + (f" ×{d['volumes']}" if d["volumes"] > 1 else "")
                     for d in items)


def render_list(console, rows: list, title: str = "", extra: dict | None = None,
                extra_title: str = "matched") -> None:
    from rich.table import Table
    t = Table(title=title or None, title_justify="left", header_style="bold")
    for col, kw in (("name", {"style": "bold", "no_wrap": True}),
                    ("module", {}), ("nuclide", {}), ("activity", {"justify": "right"}),
                    ("material", {}), ("sampling", {}), ("vols", {"justify": "right"}),
                    ("surfaces", {}), ("distrib", {})):
        t.add_column(col, **kw)
    if extra:
        t.add_column(extra_title)
    for r in rows:
        act = r["activity"] + (f" ({r['activity_field']})" if r["activity"] and
                               r["activity_field"] != "BqPercc" else "")
        cells = [r["name"], r["module_type"], r["nuclide"], act, r["material"],
                 r["sampling"], str(r["volumes"]), _surfaces_text(r["surfaces"]),
                 ",".join(r["distrib"])]
        if extra:
            cells.append(extra.get(r["name"], ""))
        t.add_row(*cells)
    console.print(t)
    console.print(f"[dim]{len(rows)} background(s)[/dim]")


def render_show(console, d: dict) -> None:
    from rich.panel import Panel
    from rich.table import Table
    p = Table.grid(padding=(0, 2))
    for k, v in d["params"].items():
        text = (", ".join(_s(x) for x in v) if isinstance(v, list)
                else json.dumps(v) if isinstance(v, dict) else _s(v))
        p.add_row(f"[bold]{k}[/bold]", text)
    if not d["params"]:
        p.add_row("[dim]no producer parameters in the IR[/dim]", "")
    console.print(Panel(p, title=f"[bold]{d['name']}[/bold]", title_align="left"))
    s = Table(title="surfaces", title_justify="left", header_style="bold")
    s.add_column("surface")
    s.add_column("volumes", justify="right")
    for x in d["surfaces"]:
        s.add_row(x["surface"], str(x["volumes"]))
    console.print(s)
    v = Table(title="volumes", title_justify="left", header_style="bold")
    for col, kw in (("#", {"justify": "right"}), ("node", {}), ("surface", {}),
                    ("x", {"justify": "right"}), ("y", {"justify": "right"}),
                    ("z", {"justify": "right"})):
        v.add_column(col, **kw)
    for r in d["volumes"]:
        v.add_row(_s(r["index"]), r["node"] or "-", r["surface"] or "-",
                  *(f"[{r['lo'][i]:.2f}, {r['hi'][i]:.2f}]" for i in range(3)))
    console.print(v)
    if d["n_volumes"] > len(d["volumes"]):
        console.print(f"[dim]… {d['n_volumes'] - len(d['volumes'])} more "
                      f"(--limit 0 shows all)[/dim]")


def render_region(console, region: Region, rows: list, total: int, frame: str,
                  volumes: dict | None = None) -> None:
    from rich.table import Table
    t = Table(title=f"active in {region.text()}  ({frame} frame)", title_justify="left",
              header_style="bold")
    t.add_column("name", style="bold", no_wrap=True)
    t.add_column("overlapping", justify="right")
    for a in AXES:
        t.add_column(f"{a} extent", justify="right")
    for r in rows:
        t.add_row(r["name"], f"{r['overlapping']} / {r['volumes']}",
                  *(f"[{r['lo'][i]:.2f}, {r['hi'][i]:.2f}]" for i in range(3)))
    console.print(t)
    console.print(f"[dim]{len(rows)} of {total} background(s)[/dim]")
    if volumes:
        v = Table(title="overlapping volumes", title_justify="left", header_style="bold")
        for col in ("background", "#", "node", "surface", "x", "y", "z"):
            v.add_column(col)
        for name, items in volumes.items():
            for r in items:
                v.add_row(name, _s(r["index"]), r["node"] or "-", r["surface"] or "-",
                          *(f"[{r['lo'][i]:.2f}, {r['hi'][i]:.2f}]" for i in range(3)))
        console.print(v)


def render_suggestions(console, d: dict, total: int) -> None:
    from rich.table import Table
    console.print(f"[bold]{total}[/bold] background(s).  Selections that split them; "
                  f"run as [bold]radgeo backgrounds IR_FILE <command>[/bold]:")
    t = Table(header_style="bold", title="by field", title_justify="left")
    t.add_column("field", style="bold")
    t.add_column("value")
    t.add_column("n", justify="right")
    t.add_column("command", overflow="fold")
    for f in d["fields"]:
        for k, v in enumerate(f["values"]):
            t.add_row(f["field"] if k == 0 else "", v["value"], str(v["backgrounds"]),
                      v["command"], end_section=(k == len(f["values"]) - 1))
    console.print(t)
    if d["numeric"]:
        u = Table(header_style="bold", title="numeric fields", title_justify="left")
        u.add_column("field", style="bold")
        u.add_column("n", justify="right")
        u.add_column("min", justify="right")
        u.add_column("max", justify="right")
        u.add_column("above the median", justify="right")
        u.add_column("command", overflow="fold")
        for x in d["numeric"]:
            u.add_row(x["field"], str(x["backgrounds"]), f"{x['min']:g}", f"{x['max']:g}",
                      str(x["matching"]), x["command"])
        console.print(u)
    if d["names"]:
        n = Table(header_style="bold", title="by words in the name", title_justify="left")
        n.add_column("word", style="bold")
        n.add_column("n", justify="right")
        n.add_column("command", overflow="fold")
        for x in d["names"]:
            n.add_row(x["token"], str(x["backgrounds"]), x["command"])
        console.print(n)
    if d["regions"]:
        r = Table(header_style="bold", title="planes the volumes lie on",
                  title_justify="left")
        r.add_column("plane", style="bold")
        r.add_column("on the plane", justify="right")
        r.add_column("volumes on it", justify="right")
        r.add_column("crossing it", justify="right")
        r.add_column("command", overflow="fold")
        for x in d["regions"]:
            r.add_row(x["slice"], str(x["on_plane"]), str(x["volumes_on_plane"]),
                      str(x["crossing"]), x["command"])
        console.print(r)
    if d["extent"]:
        console.print("[dim]the volumes span " + ",  ".join(
            f"{a} [{lo:g}, {hi:g}]" for a, (lo, hi) in d["extent"].items())
            + "  (try region --x / --y / --z / --box)[/dim]")


# ---------------------------------------------------------------------------
# command line (click)
# ---------------------------------------------------------------------------

BkgError = CliError


class _Group(click.Group):
    """A group that reports the ValueErrors of the commands (bad expression,
    unknown field, empty region) as errors."""

    def invoke(self, ctx):
        try:
            return super().invoke(ctx)
        except ValueError as exc:
            raise BkgError(str(exc)) from None


@dataclass
class _Session:
    """What the commands share: the IR file, loaded when a command needs it."""
    path:    str
    as_json: bool

    def load(self) -> tuple:
        """(ir, backgrounds); an error if the file or its backgrounds are missing."""
        p = Path(self.path)
        if not p.is_file():
            raise BkgError(f"'{self.path}' not found")
        ir = IRCollection.load(p)
        bgs = backgrounds(ir)
        if not bgs:
            raise BkgError(f"no radiological background in '{self.path}'")
        return ir, bgs

    def console(self):
        return console()

    def emit(self, data, render) -> None:
        """JSON of `data`, or render(console) for the tables."""
        if self.as_json:
            click.echo(json.dumps(data, indent=2, default=str))
        else:
            render(self.console())


def _parse_slices(items) -> dict:
    out = {}
    for s in items:
        m = re.fullmatch(r"\s*([xyzXYZ])\s*=\s*(\S+)", s)
        if not m:
            raise ValueError(f"bad --slice {s!r} (expected AXIS=VALUE, e.g. x=-328)")
        try:
            out[m.group(1).lower()] = float(m.group(2))
        except ValueError:
            raise ValueError(f"bad --slice {s!r}: {m.group(2)!r} is not a number") from None
    return out


def _vol_rows(b: Background, idx: list) -> list:
    d = detail(b, limit=0)["volumes"]
    return [d[i] for i in idx]


def _json_option(fn):
    """--json on a command: the same switch as the group's, so it can follow
    the command name (radgeo backgrounds IR_FILE select ... --json)."""
    def callback(ctx, param, value):
        if value:
            ctx.obj.as_json = True
    return click.option("--json", "json_flag", is_flag=True, expose_value=False,
                        callback=callback, help="Print JSON instead of tables.")(fn)


@click.group(cls=_Group, context_settings=CONTEXT_SETTINGS,
             help="Inspect the radiological backgrounds of a JSON IR (radgeo extract).\n\n"
                  "A background is one radiological producer: its radio volumes "
                  "(sampling boxes) and its FHiCL parameters.\n\n"
                  "\b\nExamples:\n"
                  "  radgeo backgrounds ir.json list\n"
                  "  radgeo backgrounds ir.json show Po218 Cathode\n"
                  "  radgeo backgrounds ir.json select nuclide=Po218 'BqPercc<0.01'\n"
                  "  radgeo backgrounds ir.json region --slice x=-328\n"
                  "  radgeo backgrounds ir.json region --x -330 -300 --z 0 500 --volumes\n"
                  "  radgeo backgrounds ir.json select material=LAr --json")
@click.argument("ir_file", metavar="IR_FILE", type=click.Path(dir_okay=False))
@click.option("--json", "as_json", is_flag=True,
              help="Print JSON instead of tables (also accepted after the command; "
                   "needs no rich).")
@click.pass_context
def cli(ctx, ir_file: str, as_json: bool) -> None:
    ctx.obj = _Session(ir_file, as_json)


@cli.command("list")
@click.option("--sort", type=click.Choice(["name", "volumes"]), default="name",
              show_default=True, help="Order of the rows (volumes: most first).")
@_json_option
@click.pass_obj
def list_cmd(sess: _Session, sort: str) -> None:
    """Every background: module, nuclide, activity, material, sampling, volumes."""
    _, bgs = sess.load()
    rows = [summary(b) for b in bgs]
    if sort == "volumes":
        rows.sort(key=lambda r: -r["volumes"])
    sess.emit(rows, lambda c: render_list(c, rows))


@cli.command("show")
@click.argument("names", nargs=-1, required=True, metavar="NAME_RE...")
@click.option("--limit", type=int, default=25, show_default=True,
              help="Volumes shown per background (0: all).")
@_json_option
@click.pass_obj
def show_cmd(sess: _Session, names: tuple, limit: int) -> None:
    """The FHiCL parameters, surfaces and volumes of the backgrounds whose name
    matches a regex (re.search)."""
    _, bgs = sess.load()
    try:
        pats = [re.compile(n) for n in names]
    except re.error as exc:
        raise BkgError(f"bad regex: {exc}") from None
    chosen = [b for b in bgs if any(p.search(b.name) for p in pats)]
    if not chosen:
        raise BkgError(f"no background matches {', '.join(names)}; "
                       f"see 'radgeo backgrounds {sess.path} list'")
    data = [detail(b, limit) for b in chosen]

    def render(console):
        for d in data:
            render_show(console, d)
    sess.emit(data, render)


@cli.command("select", help="""The backgrounds for which every EXPR holds.

\b
EXPR is
  FIELD=REGEX    re.search on any value of the field
  FIELD!=REGEX   no value matches (also true if the field is absent)
  FIELD<N        also <=, >, >=: any numeric value
  FIELD          the background has the field

A FIELD is a FHiCL parameter of the producer or a metadata key of its volumes
(nuclide, surface, node, source, module_type, BqPercc, T0, ...), plus name and
nvolumes; `suggest` lists them.""")
@click.argument("exprs", nargs=-1, required=True, metavar="EXPR...")
@_json_option
@click.pass_obj
def select_cmd(sess: _Session, exprs: tuple) -> None:
    _, bgs = sess.load()
    parsed = [parse_expr(e) for e in exprs]
    check_fields(parsed, bgs)
    chosen = select(bgs, parsed)
    rows = [summary(b) for b in chosen]
    shown = {b.name: "; ".join(
        f"{e.field}={','.join(b.texts(e.field)[:3])}" for e in parsed
        if e.op is not None and e.field in b.fields) for b in chosen}
    sess.emit({"expressions": [str(e) for e in parsed], "backgrounds": rows},
              lambda c: render_list(c, rows, title="  AND  ".join(str(e) for e in parsed),
                                    extra=shown if any(shown.values()) else None))


@cli.command("region", help="""The backgrounds active in a region.

The region is the intersection of the constraints, in the frame stored in the
IR; a background is active when the box of at least one of its volumes
overlaps it (touching counts).""")
@click.option("--x", "x", type=(float, float), default=None, metavar="LO HI")
@click.option("--y", "y", type=(float, float), default=None, metavar="LO HI")
@click.option("--z", "z", type=(float, float), default=None, metavar="LO HI")
@click.option("--box", type=click.Tuple([float] * 6), default=None,
              metavar="X0 X1 Y0 Y1 Z0 Z1", help="A 3D box.")
@click.option("--slice", "slices", multiple=True, metavar="AXIS=VALUE",
              help="The plane AXIS = VALUE, e.g. x=-328; repeat for a line or a point.")
@click.option("--thickness", type=float, default=0.0, show_default=True, metavar="T",
              help="Thickness of the --slice planes.")
@click.option("--volumes", is_flag=True, help="Also list the overlapping volumes.")
@_json_option
@click.pass_obj
def region_cmd(sess: _Session, x, y, z, box, slices, thickness, volumes) -> None:
    ir, bgs = sess.load()
    ranges = {a: v for a, v in (("x", x), ("y", y), ("z", z)) if v}
    region = make_region(box=box, ranges=ranges, slices=_parse_slices(slices),
                         thickness=thickness)
    hits = in_region(bgs, region)
    rows = region_rows(hits)
    frame = (ir.frame or {}).get("volume", "world")
    vols = {b.name: _vol_rows(b, idx) for b, idx in hits} if volumes else None
    data = {"region": {"lo": [None if np.isinf(v) else float(v) for v in region.lo],
                       "hi": [None if np.isinf(v) else float(v) for v in region.hi]},
            "frame": frame, "backgrounds": rows, **({"volumes": vols} if vols else {})}
    sess.emit(data, lambda c: render_region(c, region, rows, len(bgs), frame, vols))


@cli.command("suggest")
@_json_option
@click.pass_obj
def suggest_cmd(sess: _Session) -> None:
    """Selections and regions that split the backgrounds, with the commands to
    run them."""
    _, bgs = sess.load()
    data = suggestions(bgs)
    sess.emit(data, lambda c: render_suggestions(c, data, len(bgs)))


def main(argv: list | None = None) -> None:
    cli(args=argv, prog_name="radgeo backgrounds")


if __name__ == "__main__":
    main()

"""
bki-distrib  —  draw the distrib_x/y/z position profiles stored in a JSON IR

RadioGen/Decay0Gen producers with a distrib_<axis> parameter draw that
coordinate from a TF1 formula in x over the box range (TF1::GetRandom,
larsimrad BaseRadioGen).  Each page shows the formula (TF1::Eval) and the
distribution as generated, with the effective range of each marked (the
generated one, shaded, is the box bki-extract writes).
GetRandom samples from an Npx-point table (default 100), so features narrower
than range/Npx are lost; see the duneradgeo._tf1 docstring.
Producers sharing a formula and range share one plot: a page of the PDF, or
one file per plot for other formats.  Needs PyROOT.  Not modelled: the redraw
of points outside the volume_gen / material volumes.  Coordinates are world
coordinates, as the generator samples them, also for IRs written with
bki-extract --frame.

Usage:
  bki-distrib detector.json                     # -> detector_distrib.pdf
  bki-distrib detector.json --out distrib.png --log   # distrib_0.png, …
  bki-distrib detector.json --list              # text only, no matplotlib
  bki-distrib detector.json --radio 'Mesh|Cathode'

Reads the distrib_* ranges, containment and Npx that bki-extract writes.
"""
from __future__ import annotations

import argparse
import dataclasses
import re
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import _tf1
from ._fcl import DISTRIB_AXES, DISTRIB_COMMON, distrib_keys
from ._root import ROOTUnavailable
from .ir import IRCollection, VolumeEntry
from ._select import select_volumes, validate_picks
from ._style import import_mpl

_NPTS       = 100001
_EFF_COL    = "tab:orange"
_CURVE_COL  = "tab:blue"        # as generated (TF1::GetRandom)
_FORMULA_COL = "0.1"            # the formula


@dataclass
class DistribGroup:
    """One formula on one sampling range, and the producers using it."""
    axis:      str
    formula:   str
    range_:    tuple                       # (lo, hi) sampling range [cm]
    eff:       tuple                       # effective range, as generated [cm]
    formula_eff: tuple                     # effective range of the formula [cm]
    containment: float                     # fraction in either effective range
    npx:       int                         # TF1 Npx of the generated distribution
    eff_computed: bool = False             # eff recomputed for another --npx
    surfaces:  list[str]  = field(default_factory=list)
    sources:   list[str]  = field(default_factory=list)   # volume_rand:…, explicit, …
    # producer name -> [nuclide, number of placements]
    producers: dict       = field(default_factory=dict)

    def label(self) -> str:
        return f"distrib_{self.axis}"

    @property
    def dist(self) -> _tf1.Distribution:
        """The formula on the sampling range with the group's Npx (cached;
        raises ValueError when ROOT cannot parse the formula)."""
        return _tf1.distribution(self.formula, *self.range_, self.npx)


def _key_range(r) -> tuple:
    return tuple(round(float(v), 3) for v in r)


_RANGES = ("range", "eff", "formula_eff")      # distrib_keys of the three ranges


def collect(vols: list[VolumeEntry]) -> list[DistribGroup]:
    """Group the distrib_* entries of the radio volumes in `vols`.

    Raises ValueError when an entry lacks the keys bki-extract writes."""
    groups: dict[tuple, DistribGroup] = {}
    for v in vols:
        if v.kind != "radio":
            continue
        m = v.meta
        for axis in DISTRIB_AXES:
            keys = distrib_keys(axis)
            formula = m.get(keys["formula"])
            if not formula:
                continue
            need = [keys[k] for k in _RANGES] + list(DISTRIB_COMMON)
            missing = [k for k in need if k not in m]
            if missing:
                raise ValueError(f"{v.name}: IR entry lacks {', '.join(missing)}; "
                                 f"re-run bki-extract")
            rng, eff, feff = (tuple(map(float, m[keys[k]])) for k in _RANGES)
            cont, npx = float(m[DISTRIB_COMMON[0]]), int(m[DISTRIB_COMMON[1]])
            key = (axis, formula, _key_range(rng), _key_range(eff), _key_range(feff),
                   cont, npx)
            g = groups.get(key)
            if g is None:
                g = groups[key] = DistribGroup(axis=axis, formula=formula, range_=rng,
                                               eff=eff, formula_eff=feff,
                                               containment=cont, npx=npx)
            surf = m.get("surface")
            if surf and surf not in g.surfaces:
                g.surfaces.append(surf)
            src = m.get("source")
            if src and src not in g.sources:
                g.sources.append(src)
            p = g.producers.setdefault(v.name, [m.get("nuclide", ""), 0])
            p[1] += 1
    return list(groups.values())


def with_npx(g: DistribGroup, npx: int) -> DistribGroup:
    """A copy of `g` for another TF1 Npx, its generated effective range
    recomputed (needs PyROOT); `g` itself when `npx` is its own."""
    if npx == g.npx:
        return g
    return dataclasses.replace(g, eff=g.dist.with_npx(npx).interval(g.containment),
                               npx=npx, eff_computed=True)


@dataclass
class Profile:
    """A group's distributions on a grid: the formula and what is generated."""
    x:           np.ndarray
    formula:     np.ndarray      # formula pdf, normalised over the full range
    generated:   np.ndarray      # GetRandom pdf, normalised over the full range


def profile(g: DistribGroup, lo: float | None = None, hi: float | None = None,
            npts: int = _NPTS) -> Profile:
    """Both pdfs of a resolved group on [lo, hi] (default: the full range)."""
    lo = g.range_[0] if lo is None else lo
    hi = g.range_[1] if hi is None else hi
    d = g.dist
    x = np.linspace(lo, hi, npts)
    return Profile(x, d.formula_pdf(x), d.pdf(x))


def fractions(g: DistribGroup) -> tuple:
    """(generated, formula): the measured fraction of the generated
    distribution inside g.eff, and of the formula pdf inside g.formula_eff."""
    d = g.dist
    c = d.cdf(list(g.eff))
    return float(c[1] - c[0]), d.formula_fraction(*g.formula_eff)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _fmt_range(r, w: int = 0) -> str:
    """'[lo, hi]' with 3 decimals, each number right-aligned in `w` characters."""
    return f"[{r[0]:{w}.3f}, {r[1]:{w}.3f}]"


_LABEL_LINK = re.compile(r"^(?:Gen)?(?:In|At|From)(?=[A-Z0-9])")


def _common_suffix(names: list[str]) -> str:
    s = names[0]
    for name in names[1:]:
        n = 0
        while n < min(len(s), len(name)) and s[-1 - n] == name[-1 - n]:
            n += 1
        s = s[len(s) - n:]
    return s


def group_label(g: DistribGroup) -> str:
    """
    Short name of the group, from its producer names: the name itself for a
    single producer, else the part the names end with, from a CamelCase word
    on and without a leading Gen/In/At/From (``UpperMesh1x8x14`` for
    ``K42From42ArGenInUpperMesh1x8x14``, ``Rn222ChainFromPo218GenInUpperMesh1x8x14``,
    ...).  "N producers" when the names share no such part.
    """
    names = list(g.producers)
    if len(names) == 1:
        return names[0]
    m = re.search(r"[A-Z]", _common_suffix(names)) if names else None
    label = _LABEL_LINK.sub("", m.string[m.start():]) if m else ""
    return label if len(label) >= 3 else f"{len(names)} producers"


def nuclides(g: DistribGroup) -> list[str]:
    """The producers' nuclides, each once, in producer order."""
    out = []
    for nucl, _ in g.producers.values():
        if nucl and nucl not in out:
            out.append(nucl)
    return out


def page_title(g: DistribGroup, nmax: int | None = None) -> str:
    """'<label> · x profile, N producers (nuclides)'; at most `nmax`
    nuclides, then '…'."""
    label, nucl, n = group_label(g), nuclides(g), len(g.producers)
    if nmax is not None and len(nucl) > nmax:
        nucl = nucl[:nmax] + ["…"]
    what = f"{g.axis} profile"
    if n > 1:
        what += f", {n} producers"
    if nucl:
        what += f" ({', '.join(nucl)})"
    return what if label == f"{n} producers" else f"{label} · {what}"


def page_subtitle(g: DistribGroup, source: str = "") -> str:
    """Where the distribution sits, its sampling volume(s) and `source`,
    joined by ' · '."""
    pct = f"{100 * g.containment:g}%"
    lo, hi = g.eff
    plane = [s for s in g.surfaces if s.startswith(f"plane_{g.axis}@")]
    if plane:
        where = (f"layer at {g.axis} ≈ {float(plane[0].split('@')[1]):.1f} cm, "
                 f"{pct} in {hi - lo:.1f} cm")
    else:
        where = f"{pct} in [{lo:.1f}, {hi:.1f}] cm"
        other = [_surface_text(s) for s in g.surfaces]
        if other:
            where = f"{', '.join(other)}, {where}"
    parts = [where] + list(g.sources) + ([source] if source else [])
    return "  ·  ".join(parts)


def print_groups(groups: list[DistribGroup]) -> None:
    for i, g in enumerate(groups):
        print(f"[{i}] {group_label(g)} · {g.label()}  range {_fmt_range(g.range_)}"
              f"  containment {g.containment:g}  npx {g.npx}")
        try:
            frac_gen, frac_formula = fractions(g)
            print(f"    eff_range (generated) {_fmt_range(g.eff)}  (measured {frac_gen:.4f})"
                  + ("  [recomputed for --npx]" if g.eff_computed else ""))
            print(f"    formula eff range     {_fmt_range(g.formula_eff)}"
                  f"  (measured {frac_formula:.4f})")
        except ValueError as exc:
            print(f"    [cannot evaluate: {exc}]")
        if g.surfaces:
            print(f"    surface: {', '.join(g.surfaces)}")
        print(f"    formula: {g.formula}")
        for name, (nucl, n) in g.producers.items():
            print(f"    {name}  ({nucl})" + (f"  x{n} placements" if n > 1 else ""))
        print()


# ---- formula pretty-printing ---------------------------------------------

_EXP_SIGN = re.compile(r"[0-9.][eE]$")       # '1.2e' before a sign: exponent


def _strip_outer(s: str) -> str:
    """'(a*b)' -> 'a*b' when the parentheses enclose the whole term."""
    while s.startswith("(") and s.endswith(")"):
        depth = 0
        for i, ch in enumerate(s):
            depth += (ch == "(") - (ch == ")")
            if depth == 0 and i < len(s) - 1:
                return s
        s = s[1:-1]
    return s


def _tidy(s: str) -> str:
    s = s.replace("TMath::", "")
    s = re.sub(r"\bAbs\(([^()]*)\)", r"|\1|", s)
    s = re.sub(r"\bPower\((\|[^|]*\||[^,()]+),\s*([^()]+)\)", r"\1^\2", s)
    s = re.sub(r"\b(Exp|Sqrt|Log|Gaus)\(", lambda m: m.group(1).lower() + "(", s)
    s = re.sub(r"(\d)\.(?!\d)", r"\1", s)          # '325.' -> '325'
    s = re.sub(r"[eE]([+-]?)0+(\d)", r"e\1\2", s)   # 'e-007' -> 'e-7'
    s = s.replace("e+", "e")
    return s.replace("*", "·")


def _terms(formula: str) -> list[str]:
    """The formula split before each top-level '+' / '-' (not an exponent sign);
    the pieces concatenate back to the formula exactly (spaces included)."""
    terms, cur, depth = [], "", 0
    for ch in formula:
        if ch in "+-" and depth == 0 and cur.strip() and not _EXP_SIGN.search(cur.rstrip()):
            terms.append(cur)
            cur = ch
            continue
        depth += (ch == "(") - (ch == ")")
        cur += ch
    terms.append(cur)
    return terms


def raw_formula_lines(formula: str, axis: str = "x", width: int = 86) -> list[str]:
    """The formula exactly as in the FHiCL, ``distrib_x: "…"``, wrapped at
    top-level terms (no characters added or dropped)."""
    lead = f'distrib_{axis}: "'
    pad = " " * len(lead)
    avail = max(20, width - len(lead))
    lines, cur, last = [], "", ""
    for term in _terms(formula):
        while len(term) > avail:                       # a single over-long term
            if cur:
                lines.append(cur)
                cur = ""
            lines.append(term[:avail])
            term = term[avail:]
        if cur and len(cur) + len(term) > avail:
            lines.append(cur)
            cur = ""
        cur += term
        last = term
    if len(cur) + 1 > avail:                           # room for the closing quote
        head = cur[:len(cur) - len(last)]
        if head:
            lines.append(head)
            cur = last
        if len(cur) + 1 > avail:
            lines.append(cur[:avail - 1])
            cur = cur[avail - 1:]
    lines.append(cur + '"')
    return [(lead if k == 0 else pad) + s for k, s in enumerate(lines)]


def formula_lines(formula: str, width: int = 70) -> list[str]:
    """
    The TF1 formula split into its top-level terms, one per line:
    ``f(x) = a`` / ``     + b`` …, with TMath:: dropped, Abs(a) as |a|,
    Power(a,n) as a^n and '*' as '·'.  Monospace layout, for display only.
    """
    terms = _terms(formula.replace(" ", ""))

    lead, pad = "f(x) = ", " " * 7
    lines = []
    for k, term in enumerate(terms):
        sign, body = ("", term) if term[0] not in "+-" else (term[0], term[1:])
        body = _strip_outer(body)
        keep = sign == "-" and len(_terms(body)) > 1   # -(a+b) is not -a+b
        body = _tidy(body)
        if keep:
            body = f"({body})"
        if k == 0:
            first, prefix = lead + ("-" if sign == "-" else ""), pad
        else:
            first, prefix = pad + f"{sign} ", pad + "  "
        wrapped = textwrap.wrap(body, width=max(20, width - len(first)),
                                break_on_hyphens=False) or [""]
        lines.append(first + wrapped[0])
        lines += [prefix + "  " + w for w in wrapped[1:]]
    return lines


def _surface_text(s: str) -> str:
    m = re.fullmatch(r"plane_([xyz])@(.+)", s)
    return f"plane {m.group(1)} = {m.group(2)} cm" if m else s


# ---- page layout ----------------------------------------------------------
# Positions are figure fractions; sizes are points (1/72 inch).

_MONO   = "DejaVu Sans Mono"
_GREY   = "0.35"
_MAX_PRODUCERS = 16
_LEFT   = 0.065          # left edge of the title, the plots and the formula text
_RIGHT  = 0.985          # right margin: page number, side panel text
_TOP    = 0.875          # top of the plot and of the side panel
_AX_W   = 0.565          # plot width
_XLABEL_H = 0.085        # room for the x label (and tick labels) under a plot
_ZOOM_SPLIT = 0.5        # share of the plot height for the full-range plot
_TITLE_Y, _SOURCE_Y, _FOOT_Y = 0.950, 0.915, 0.015
_PANEL  = 0.665          # side panel left edge
_HEADING_H   = 0.035     # side panel heading, down to the text under it
_SECTION_GAP = 0.035     # side panel: notes to the next heading
_NOTE_GAP    = 0.012     # side panel: sampling table to the notes
_BOX_X, _BOX_Y = _LEFT, 0.04   # formula box: text left edge, box bottom
_BOX_PAD, _BOX_GAP = 0.012, 0.012   # inside the box: top/bottom, between blocks
_BOX_INSET   = 0.006     # box edge to text, left and right; corner radius
_BOX_SIZE    = 10.0      # largest formula font size
_CAP_SIZE    = 8.5       # caption font size
_HEADING_SIZE = 11.5
_NOTE_SIZE   = 8
_CHAR_W      = 0.61      # mono character width / font size
_CAP_CHAR_W  = 0.5       # mean caption (sans) character width / font size
_LS_MONO     = 1.45      # line spacings (× font size)
_LS_TABLE    = 1.5
_LS_NOTE     = 1.4
_LS_CAP      = 1.6       # caption line plus the gap under it
_ZOOM_FRAC   = 0.5       # zoom when the eff range is narrower than this × range
_ZOOM_MARGIN = 0.25      # zoom window: eff range ± this × its width


def _text_h(nlines: int, size: float, ls: float, fig) -> float:
    """Height of `nlines` lines of `size`-point text with line spacing `ls`,
    in figure units."""
    return nlines * size * ls / 72 / fig.get_figheight()


def _mono_size(lines: list[str], fig, x0: float, base: float,
               right: float = _RIGHT) -> float:
    """Largest size <= `base` at which the longest mono line fits in [x0, right]."""
    chars = max((len(s) for s in lines), default=1)
    avail = (right - x0) * fig.get_figwidth() * 72          # points
    return min(base, avail / (_CHAR_W * chars))


def _plot_pdf(ax, pr: Profile, g: DistribGroup, log: bool, legend: bool) -> None:
    c = f"{100 * g.containment:g}%"
    ax.plot(pr.x, pr.formula, color=_FORMULA_COL, lw=1.4, label="formula (TF1::Eval)")
    ax.plot(pr.x, pr.generated, color=_CURVE_COL, lw=1.8, alpha=0.9,
            label=f"as generated (TF1::GetRandom, Npx={g.npx})")
    ax.axvspan(*g.eff, color=_EFF_COL, alpha=0.22, lw=0,
               label=f"effective range, as generated ({c})")
    for e in g.eff:
        ax.axvline(e, color=_EFF_COL, lw=1.0, ls="--")
    for k, e in enumerate(g.formula_eff):
        ax.axvline(e, color=_FORMULA_COL, lw=1.0, ls=":",
                   label=f"effective range of the formula ({c})" if k == 0 else None)
    if log:
        ax.set_yscale("log")
    if legend:
        ax.legend(loc="best", fontsize=10, framealpha=0.9)
    ax.set_xlim(pr.x[0], pr.x[-1])
    ax.set_xlabel(f"{g.axis}  [cm, world]", fontsize=11)
    ax.set_ylabel("probability density  [1/cm]", fontsize=11)
    ax.tick_params(labelsize=10)
    ax.grid(alpha=0.3)


def zoom_window(g: DistribGroup) -> tuple | None:
    """The zoom window around both effective ranges, or None when no zoom is
    useful (they cover more than _ZOOM_FRAC of the range)."""
    a = min(g.eff[0], g.formula_eff[0])
    b = max(g.eff[1], g.formula_eff[1])
    w, W = b - a, g.range_[1] - g.range_[0]
    if not (W > 0 and w < _ZOOM_FRAC * W):
        return None
    return (max(g.range_[0], a - _ZOOM_MARGIN * w), min(g.range_[1], b + _ZOOM_MARGIN * w))


def _draw_formula_box(fig, g: DistribGroup) -> float:
    """The formula box at the bottom left: the FHiCL parameter, then the
    compact form, each under a caption (sans, grey).  Returns its top."""
    from matplotlib.patches import FancyBboxPatch
    raw, comp = raw_formula_lines(g.formula, g.axis), formula_lines(g.formula)
    blocks = [("FHiCL · evaluated by ROOT", raw), ("compact form · display only", comp)]
    lines = raw + comp
    fsize = min(_BOX_SIZE, _mono_size(lines, fig, _BOX_X, _BOX_SIZE, right=_BOX_X + _AX_W))
    mono_h = _text_h(1, fsize, _LS_MONO, fig)                # one mono line
    cap_h = _text_h(1, _CAP_SIZE, _LS_CAP, fig)
    box_h = (2 * _BOX_PAD + len(blocks) * cap_h + len(lines) * mono_h
             + (len(blocks) - 1) * _BOX_GAP)
    box_top = _BOX_Y + box_h

    width = max(max(len(s) for s in lines) * _CHAR_W * fsize,
                max(len(c) for c, _ in blocks) * _CAP_CHAR_W * _CAP_SIZE)
    width = min(width / 72 / fig.get_figwidth(), _AX_W) + 2 * _BOX_INSET
    fig.add_artist(FancyBboxPatch(
        (_BOX_X - _BOX_INSET, _BOX_Y), width, box_h, transform=fig.transFigure,
        boxstyle=f"round,pad=0,rounding_size={_BOX_INSET}", fc="0.965", ec="0.8",
        lw=0.8, zorder=0))
    y = box_top - _BOX_PAD
    for k, (caption, block) in enumerate(blocks):
        if k:
            y -= _BOX_GAP
        fig.text(_BOX_X, y, caption, fontsize=_CAP_SIZE, color=_GREY, style="italic",
                 va="top", ha="left")
        y -= cap_h
        fig.text(_BOX_X, y, "\n".join(block), family=_MONO, fontsize=fsize,
                 va="top", ha="left", linespacing=_LS_MONO)
        y -= len(block) * mono_h
    return box_top


def _draw_plots(fig, g: DistribGroup, log: bool, zoom: bool,
                bottom: float) -> tuple | None:
    """Both pdfs over the full range, between `bottom` and _TOP, with a zoom
    on the effective ranges under it when `zoom` and zoom_window(g) allow.
    Returns fractions(g), or None when the formula cannot be evaluated."""
    try:
        pr, frac = profile(g), fractions(g)
    except ValueError as exc:
        ax = fig.add_axes((_LEFT, bottom, _AX_W, _TOP - bottom))
        ax.text(0.5, 0.5, f"cannot evaluate the formula:\n{exc}", fontsize=11,
                ha="center", va="center", transform=ax.transAxes, color="tab:red")
        ax.set_xlim(*g.range_)
        return None
    win = zoom_window(g) if zoom else None
    if win is None:
        ax = fig.add_axes((_LEFT, bottom, _AX_W, _TOP - bottom))
        _plot_pdf(ax, pr, g, log, legend=True)
        return frac
    h = _TOP - bottom - _XLABEL_H                 # both plots, without the top x label
    h_main = _ZOOM_SPLIT * h
    ax = fig.add_axes((_LEFT, _TOP - h_main, _AX_W, h_main))
    ax_z = fig.add_axes((_LEFT, bottom, _AX_W, h - h_main))
    _plot_pdf(ax, pr, g, log, legend=True)
    _plot_pdf(ax_z, profile(g, *win), g, log, legend=False)
    ax_z.set_title(f"zoom: effective ranges ± {100 * _ZOOM_MARGIN:g}% of their span",
                   fontsize=10, loc="left", color=_GREY)
    return frac


def _draw_side_panel(fig, g: DistribGroup, frac: tuple | None) -> None:
    """The right-hand panel: the sampling table (with the measured fractions
    `frac`, if given), notes, and the producers."""
    cont = f"{100 * g.containment:g}%"
    if frac is not None:
        cont += (f"  (measured: generated {100 * frac[0]:.2f}%,"
                 f" formula {100 * frac[1]:.2f}%)")
    rows = [("range", _fmt_range(g.range_, 9) + " cm"),
            ("eff range", _fmt_range(g.eff, 9) + " cm  generated"
             + ("*" if g.eff_computed else "")),
            ("", _fmt_range(g.formula_eff, 9) + " cm  formula"),
            ("containment", cont),
            ("Npx", f"{g.npx}  (TF1::GetRandom table)"),
            ("volume", ", ".join(g.sources) or "—")]
    kw = max(len(k) for k, _ in rows)
    sampling = [f"{k:<{kw}}  {v}" for k, v in rows]

    names = list(g.producers.items())
    shown = names[:_MAX_PRODUCERS]
    nw = max(len(n) for n, _ in shown)
    prods = [f"{n:<{nw}}  {nucl}" + (f"  ×{k}" if k > 1 else "") for n, (nucl, k) in shown]
    if len(names) > len(shown):
        prods.append(f"… and {len(names) - len(shown)} more (see --list)")

    notes = ["not modelled: redraw outside the volume_gen / material volumes"]
    if g.eff_computed:
        notes.insert(0, f"* recomputed for --npx {g.npx} (the IR has another Npx)")

    x0, y = _PANEL, _TOP
    fig.text(x0, y, "Sampling", fontsize=_HEADING_SIZE, fontweight="bold", va="top")
    y -= _HEADING_H
    fs = _mono_size(sampling, fig, x0, 9.0)
    fig.text(x0, y, "\n".join(sampling), family=_MONO, fontsize=fs, va="top",
             linespacing=_LS_TABLE)
    y -= _text_h(len(sampling), fs, _LS_TABLE, fig) + _NOTE_GAP
    fig.text(x0, y, "\n".join(notes), fontsize=_NOTE_SIZE, color=_GREY, va="top",
             linespacing=_LS_NOTE)
    y -= _text_h(len(notes), _NOTE_SIZE, _LS_NOTE, fig) + _SECTION_GAP
    fig.text(x0, y, f"Producers ({len(names)})", fontsize=_HEADING_SIZE,
             fontweight="bold", va="top")
    y -= _HEADING_H
    fig.text(x0, y, "\n".join(prods), family=_MONO, fontsize=_mono_size(prods, fig, x0, 8.5),
             va="top", linespacing=_LS_MONO)


def _fit_text(fig, y: float, make, n: int, **kw) -> None:
    """Draw make(k) at (_LEFT, y) for the largest k <= n whose text ends
    left of _RIGHT (make(None) is the full text, tried first)."""
    renderer = fig.canvas.get_renderer()
    for k in [None] + list(range(n - 1, -1, -1)):
        t = fig.text(_LEFT, y, make(k), va="baseline", **kw)
        if (k == 0 or t.get_window_extent(renderer).x1
                <= _RIGHT * fig.bbox.width):
            return
        t.remove()


def _draw_page(fig, g: DistribGroup, log: bool, source: str,
               page: int, npages: int, zoom: bool = True) -> None:
    """One distribution on `fig`: title, plot (+ zoom), formula box, side panel."""
    _fit_text(fig, _TITLE_Y, lambda k: page_title(g, k), len(nuclides(g)),
              fontsize=15, fontweight="bold")
    fig.text(_LEFT, _SOURCE_Y, page_subtitle(g, source), fontsize=10,
             color=_GREY, va="baseline")
    fig.text(_RIGHT, _FOOT_Y, f"{page} / {npages}", fontsize=9, color=_GREY, ha="right")
    box_top = _draw_formula_box(fig, g)
    frac = _draw_plots(fig, g, log, zoom, bottom=box_top + _XLABEL_H)
    _draw_side_panel(fig, g, frac)


def plot_groups(groups: list[DistribGroup], out_path: str, title: str = "",
                log: bool = False, zoom: bool = True,
                doc_title: str = "") -> list[str]:
    """
    One page per distribution: a multi-page PDF when `out_path` ends in .pdf,
    otherwise one file per distribution (<stem>_NN<suffix>).
    `title` (e.g. the source FHiCL) ends the subtitle of each page;
    `doc_title` is the PDF's Title metadata.  Returns the files written.
    """
    plt = import_mpl("--list works without it")
    from matplotlib.backends.backend_pdf import PdfPages

    out = Path(out_path)
    n   = len(groups)

    def _page(i: int, g: DistribGroup):
        fig = plt.figure(figsize=(14.0, 8.0))
        _draw_page(fig, g, log=log, source=title, page=i + 1, npages=n, zoom=zoom)
        return fig

    if out.suffix.lower() == ".pdf":
        with PdfPages(out) as pdf:
            for i, g in enumerate(groups):
                fig = _page(i, g)
                pdf.savefig(fig)
                plt.close(fig)
            pdf.infodict()["Title"] = doc_title or "distrib_* position profiles"
        written = [str(out)]
        print(f"[info] wrote {n} page(s) to {out}", file=sys.stderr)
    else:
        width = len(str(n - 1))
        written = []
        for i, g in enumerate(groups):
            path = out.with_name(f"{out.stem}_{i:0{width}d}{out.suffix}")
            fig = _page(i, g)
            fig.savefig(path, dpi=110)
            plt.close(fig)
            written.append(str(path))
        print(f"[info] wrote {n} file(s): {written[0]} … {written[-1]}"
              if n > 1 else f"[info] wrote {written[0]}", file=sys.stderr)
    return written


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("ir_file", help="JSON IR produced by bki-extract")
    ap.add_argument("--out", "-o", default=None,
                    help="output file (default: <stem>_distrib.pdf): a PDF gets "
                         "one page per distribution, other formats one file "
                         "per distribution, <stem>_NN<suffix>")
    ap.add_argument("--list", "-l", action="store_true",
                    help="print the distributions, effective ranges and producers, "
                         "then exit (needs ROOT, not matplotlib)")
    ap.add_argument("--radio", nargs="+", default=None, metavar="PATTERN",
                    help="regex filter on radio producer names")
    ap.add_argument("--pick", nargs="+", default=None, metavar="SPEC",
                    help="PRODUCER_RE#SEL as in bki-plot --surface (index list "
                         "or @REGEX on surface, node, nuclide)")
    ap.add_argument("--axis", nargs="+", choices=DISTRIB_AXES, default=None,
                    help="only these axes (default: all present)")
    ap.add_argument("--npx", type=int, default=None,
                    help="TF1 Npx for the generated distribution (default: the "
                         "distrib_eff_npx stored in the IR); another value recomputes the "
                         "generated effective range")
    ap.add_argument("--log", action="store_true", help="logarithmic y axis")
    ap.add_argument("--no-zoom", action="store_true",
                    help="no zoomed view on the effective range (by default it is "
                         "added below the plot when the effective range is "
                         "narrower than half the full range)")
    ap.add_argument("--title", default=None,
                    help="text ending each page subtitle (default: the source FHiCL name)")
    args = ap.parse_args()

    ir_path = Path(args.ir_file)
    if not ir_path.is_file():
        sys.exit(f"[error] '{args.ir_file}' not found")
    ir = IRCollection.load(ir_path)

    try:
        validate_picks(pick=args.pick)
        vols = select_volumes(ir, kinds=["radio"],
                              kind_patterns={"radio": args.radio} if args.radio else None,
                              picks=args.pick)
    except ValueError as exc:
        sys.exit(f"[error] {exc}")

    try:
        groups = collect(vols)
    except ValueError as exc:
        sys.exit(f"[error] {exc}")
    if args.axis:
        groups = [g for g in groups if g.axis in args.axis]
    if not groups:
        sys.exit("[error] no radio volumes with a distrib_* formula "
                 "(after the selectors)")
    print(f"[info] {len(groups)} distribution(s) from "
          f"{sum(len(g.producers) for g in groups)} producer entr(ies)", file=sys.stderr)

    # every group is evaluated by ROOT; drop those whose formula it cannot parse
    # or that integrates to 0 over its range (checked on the dense grid, so a
    # narrow peak is kept)
    kept = []
    for g in groups:
        try:
            if args.npx:
                g = with_npx(g, args.npx)
            g.dist.norm                       # raises ValueError
            kept.append(g)
        except ROOTUnavailable as exc:
            sys.exit(f"[error] {exc}")
        except ValueError as exc:
            print(f"[warn] {g.label()} {g.formula!r}: {exc}; skipped", file=sys.stderr)
    groups = kept
    if not groups:
        sys.exit("[error] no distribution could be evaluated")

    if args.list:
        print_groups(groups)
        return

    out_path = args.out or str(ir_path.with_name(ir_path.stem + "_distrib.pdf"))
    title = args.title or (Path(ir.source).name if ir.source else ir_path.stem)
    try:
        plot_groups(groups, out_path, title=title, log=args.log, zoom=not args.no_zoom,
                    doc_title=f"{ir_path.stem}: distrib_* profiles")
    except ImportError as exc:
        sys.exit(f"[error] {exc}")


if __name__ == "__main__":
    main()

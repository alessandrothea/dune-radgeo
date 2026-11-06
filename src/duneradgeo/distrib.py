"""
bki-distrib  —  draw the distrib_x/y/z position profiles stored in a JSON IR

RadioGen/Decay0Gen producers with a distrib_<axis> parameter do not sample
that coordinate uniformly.  The parameter is a TF1 formula in the variable x,
drawn with TF1::GetRandom over the box range on that axis (larsimrad
BaseRadioGen).  This command draws each formula as a pdf over that range,
with the effective range bki-extract used for the box (<axis>_eff_range, the
shortest interval holding eff_containment of the pdf) shaded.  Producers sharing a formula and range share one plot; each
plot is a page of the PDF (other formats: one file per plot).
Coordinates are world coordinates, as the generator samples them, also for
IRs written with bki-extract --frame.

Usage:
  bki-distrib detector.json                     # -> detector_distrib.pdf
  bki-distrib detector.json --out distrib.png --log   # distrib_0.png, …
  bki-distrib detector.json --list              # text only, no matplotlib
  bki-distrib detector.json --radio 'Mesh|Cathode'

IRs written before bki-extract stored <axis>_range only have the containment
interval: the profile is then drawn over that interval only (or over
--range LO HI).  Re-run bki-extract to get the full sampling range.
"""
from __future__ import annotations

import argparse
import re
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ._fcl import DISTRIB_AXES, eval_tf1
from .ir import IRCollection, VolumeEntry
from .plot import parse_pick_spec, select_volumes

_NPTS       = 100001
_EFF_COL    = "tab:orange"
_CURVE_COL  = "tab:blue"


@dataclass
class DistribGroup:
    """One formula on one sampling range, and the producers using it."""
    axis:      str
    formula:   str
    rng:       tuple                       # (lo, hi) sampling range [cm]
    full:      bool                        # rng is the stored <axis>_range
    eff:       tuple | None = None         # effective range (containment interval) [cm]
    containment: float | None = None       # eff_containment stored by bki-extract
    surfaces:  list[str]  = field(default_factory=list)
    sources:   list[str]  = field(default_factory=list)   # volume_rand:…, explicit, …
    # producer name -> [nuclide, number of placements]
    producers: dict       = field(default_factory=dict)

    def label(self) -> str:
        return f"distrib_{self.axis}"


def _key_range(r) -> tuple:
    return tuple(round(float(v), 3) for v in r)


def collect(vols: list[VolumeEntry], fallback: tuple | None = None) -> list[DistribGroup]:
    """Group the distrib_* entries of the radio volumes in `vols`."""
    groups: dict[tuple, DistribGroup] = {}
    for v in vols:
        if v.kind != "radio":
            continue
        m = v.meta
        for ax, axis in enumerate(DISTRIB_AXES):
            formula = m.get(f"distrib_{axis}")
            if not formula:
                continue
            # <axis>_eff_range; IRs before the rename call it <axis>_eff
            eff  = m.get(f"{axis}_eff_range", m.get(f"{axis}_eff"))
            rng  = m.get(f"{axis}_range")
            full = rng is not None
            if rng is None:
                if fallback is not None:
                    rng = fallback
                else:
                    # the box extent on this axis is the containment interval
                    rng = eff if eff else (v.pos[ax] - v.ext[ax] / 2,
                                           v.pos[ax] + v.ext[ax] / 2)
            cont = m.get("eff_containment")     # absent in IRs before it was stored
            key = (axis, formula, _key_range(rng), _key_range(eff) if eff else None, cont)
            g = groups.get(key)
            if g is None:
                g = groups[key] = DistribGroup(axis=axis, formula=formula,
                                               rng=tuple(map(float, rng)), full=full,
                                               eff=tuple(map(float, eff)) if eff else None,
                                               containment=cont)
            surf = m.get("surface")
            if surf and surf not in g.surfaces:
                g.surfaces.append(surf)
            src = m.get("source")
            if src and src not in g.sources:
                g.sources.append(src)
            p = g.producers.setdefault(v.name, [m.get("nuclide", ""), 0])
            p[1] += 1
    return list(groups.values())


def _integrate(f: np.ndarray, x: np.ndarray) -> float:
    """Trapezoid rule (np.trapz is gone in numpy 2, np.trapezoid absent in 1.x)."""
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(x)))


def evaluate(g: DistribGroup, npts: int = _NPTS) -> tuple:
    """
    (x, pdf, fraction in eff, norm): the formula normalised to 1 over g.rng.
    The fraction is integrated on its own grid over the effective range, so
    it does not depend on where the range edges fall on the full-range grid.
    """
    x = np.linspace(g.rng[0], g.rng[1], npts)
    f = eval_tf1(g.formula, x)
    norm = _integrate(f, x)
    if not norm > 0:
        raise ValueError("formula integrates to 0 over the range")
    frac = None
    if g.eff is not None:
        xe, pe = pdf_on(g, g.eff[0], g.eff[1], norm, npts)
        frac = _integrate(pe, xe)
    return x, f / norm, frac, norm


def pdf_on(g: DistribGroup, lo: float, hi: float, norm: float,
           npts: int = _NPTS) -> tuple:
    """(x, pdf) on [lo, hi], with the normalisation `norm` of the full range."""
    x = np.linspace(lo, hi, npts)
    return x, eval_tf1(g.formula, x) / norm


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _fmt_rng(r) -> str:
    return f"[{r[0]:.3f}, {r[1]:.3f}]"


def print_groups(groups: list[DistribGroup]) -> None:
    for i, g in enumerate(groups):
        try:
            _, _, frac, _ = evaluate(g)
            ftxt = ""
            if g.containment is not None:
                ftxt += f"  containment {g.containment:g}"
            if frac is not None:
                ftxt += f"  (measured {frac:.4f})"
        except Exception as exc:                      # noqa: BLE001
            ftxt = f"  [cannot evaluate: {exc}]"
        print(f"[{i}] {g.label()}  range {_fmt_rng(g.rng)}"
              f"{'' if g.full else ' (containment interval only)'}"
              f"  eff_range {_fmt_rng(g.eff) if g.eff else '-'}{ftxt}")
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


def formula_lines(formula: str, width: int = 70) -> list[str]:
    """
    The TF1 formula split into its top-level terms, one per line:
    ``f(x) = a`` / ``     + b`` …, with TMath:: dropped, Abs(a) as |a|,
    Power(a,n) as a^n and '*' as '·'.  Monospace layout.
    """
    s = formula.replace(" ", "")
    terms, cur, depth = [], "", 0
    for ch in s:
        if ch in "+-" and depth == 0 and cur and not _EXP_SIGN.search(cur):
            terms.append(cur)
            cur = ch
            continue
        depth += (ch == "(") - (ch == ")")
        cur += ch
    terms.append(cur)

    lead, pad = "f(x) = ", " " * 7
    lines = []
    for k, term in enumerate(terms):
        sign, body = ("", term) if term[0] not in "+-" else (term[0], term[1:])
        body = _tidy(_strip_outer(body))
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

_MONO   = "DejaVu Sans Mono"
_GREY   = "0.35"
_MAX_PRODUCERS = 16
_TOP    = 0.875          # top of the plot and of the side panel (figure units)
_AX_W   = 0.565          # plot width
_PANEL  = 0.665          # side panel left edge
_ZOOM_FRAC   = 0.5       # zoom when the eff range is narrower than this × range
_ZOOM_MARGIN = 0.25      # zoom window: eff range ± this × its width


def _mono_size(lines: list[str], fig, x0: float, base: float) -> float:
    """Largest size <= `base` at which the longest mono line fits in [x0, 0.985]."""
    chars = max((len(s) for s in lines), default=1)
    avail = (0.985 - x0) * fig.get_figwidth() * 72          # points
    return min(base, avail / (0.61 * chars))


def _eff_label(g: DistribGroup, frac: float | None) -> str:
    if g.containment is not None:
        return f"effective range ({100 * g.containment:g}% containment)"
    if frac is not None:
        return f"effective range ({100 * frac:.1f}% of the pdf)"
    return "effective range"


def _plot_pdf(ax_, x, pdf, g: DistribGroup, frac, log: bool, legend: bool) -> None:
    ax_.plot(x, pdf, color=_CURVE_COL, lw=1.8, label="pdf")
    if g.eff is not None:
        ax_.axvspan(*g.eff, color=_EFF_COL, alpha=0.22, lw=0, label=_eff_label(g, frac))
        for e in g.eff:
            ax_.axvline(e, color=_EFF_COL, lw=1.0, ls="--")
    if log:
        ax_.set_yscale("log")
    if legend:
        ax_.legend(loc="best", fontsize=10, framealpha=0.9)
    ax_.set_xlim(x[0], x[-1])
    ax_.set_xlabel(f"{g.axis}  [cm, world]", fontsize=11)
    ax_.set_ylabel("probability density  [1/cm]", fontsize=11)
    ax_.tick_params(labelsize=10)
    ax_.grid(alpha=0.3)


def zoom_window(g: DistribGroup) -> tuple | None:
    """The zoom window around the eff range, or None when no zoom is useful."""
    if g.eff is None:
        return None
    w, W = g.eff[1] - g.eff[0], g.rng[1] - g.rng[0]
    if not (W > 0 and w < _ZOOM_FRAC * W):
        return None
    return (max(g.rng[0], g.eff[0] - _ZOOM_MARGIN * w),
            min(g.rng[1], g.eff[1] + _ZOOM_MARGIN * w))


def _draw_page(fig, g: DistribGroup, log: bool, source: str,
               page: int, npages: int, zoom: bool = True) -> None:
    """One distribution on `fig`: title, plot (+ zoom), formula box, side panel."""
    flines = formula_lines(g.formula)
    line_h = 10.5 * 1.45 / 72 / fig.get_figheight()          # one mono line, fig units
    box_top = 0.045 + len(flines) * line_h + 0.025
    ax_bottom = box_top + 0.085

    # ---- title ----
    where = ", ".join(_surface_text(s) for s in g.surfaces) or "—"
    fig.text(0.065, 0.950, f"distrib_{g.axis}   ·   {where}",
             fontsize=15, fontweight="bold", va="baseline")
    if source:
        fig.text(0.065, 0.915, source, fontsize=10, color=_GREY, va="baseline")
    fig.text(0.985, 0.015, f"{page} / {npages}", fontsize=9, color=_GREY, ha="right")

    # ---- plot (+ zoom on the effective range) ----
    try:
        x, pdf, frac, norm = evaluate(g)
        err = None
    except Exception as exc:                      # noqa: BLE001
        x = pdf = frac = norm = None
        err = exc
    win = zoom_window(g) if (zoom and err is None) else None
    if win is None:
        ax_ = fig.add_axes((0.065, ax_bottom, _AX_W, _TOP - ax_bottom))
    else:
        h = _TOP - ax_bottom - 0.085                       # 0.085: x label of the top plot
        h_main = 0.5 * h
        ax_ = fig.add_axes((0.065, _TOP - h_main, _AX_W, h_main))
        ax_z = fig.add_axes((0.065, ax_bottom, _AX_W, h - h_main))
    if err is not None:
        ax_.text(0.5, 0.5, f"cannot evaluate the formula:\n{err}", fontsize=11,
                 ha="center", va="center", transform=ax_.transAxes, color="tab:red")
        ax_.set_xlim(*g.rng)
    else:
        _plot_pdf(ax_, x, pdf, g, frac, log, legend=True)
    if win is not None:
        xz, pz = pdf_on(g, win[0], win[1], norm)
        _plot_pdf(ax_z, xz, pz, g, frac, log, legend=False)
        ax_z.set_title(f"zoom: effective range ± {100 * _ZOOM_MARGIN:g}% of its width",
                       fontsize=10, loc="left", color=_GREY)

    # ---- formula (monospace) ----
    fig.text(0.065, box_top - 0.012, "\n".join(flines), family=_MONO, fontsize=10.5,
             va="top", ha="left", linespacing=1.45,
             bbox=dict(boxstyle="round,pad=0.5", fc="0.96", ec="0.82"))

    # ---- side panel ----
    def rng_txt(r):
        return f"[{r[0]:9.3f}, {r[1]:9.3f}]"

    rows = [("range", rng_txt(g.rng) + " cm")]
    if g.eff is not None:
        rows.append(("eff range", rng_txt(g.eff) + " cm"))
        cont = [f"{100 * g.containment:g}%"] if g.containment is not None else []
        if frac is not None:
            cont.append(f"(measured {100 * frac:.2f}%)" if cont
                        else f"{100 * frac:.2f}% of the pdf (measured)")
        if cont:
            rows.append(("containment", " ".join(cont)))
    rows.append(("volume", ", ".join(g.sources) or "—"))
    kw = max(len(k) for k, _ in rows)
    sampling = "\n".join(f"{k:<{kw}}  {v}" for k, v in rows)

    names = list(g.producers.items())
    shown = names[:_MAX_PRODUCERS]
    nw = max(len(n) for n, _ in shown)
    prods = "\n".join(f"{n:<{nw}}  {nucl}" + (f"  ×{k}" if k > 1 else "")
                      for n, (nucl, k) in shown)
    if len(names) > len(shown):
        prods += f"\n… and {len(names) - len(shown)} more (see --list)"

    x0, y = _PANEL, _TOP
    fig.text(x0, y, "Sampling", fontsize=11.5, fontweight="bold", va="top")
    y -= 0.035
    fs = _mono_size(sampling.split("\n"), fig, x0, 9.0)
    fig.text(x0, y, sampling, family=_MONO, fontsize=fs, va="top", linespacing=1.5)
    y -= (sampling.count("\n") + 1) * fs * 1.5 / 72 / fig.get_figheight() + 0.04
    if not g.full:
        fig.text(x0, y, "Full range not stored in this IR: drawn over the\n"
                        "containment interval only. Re-run bki-extract.",
                 fontsize=9, color="tab:red", va="top")
        y -= 0.07
    fig.text(x0, y, f"Producers ({len(names)})", fontsize=11.5, fontweight="bold", va="top")
    y -= 0.035
    fig.text(x0, y, prods, family=_MONO, fontsize=_mono_size(prods.split("\n"), fig, x0, 8.5),
             va="top", linespacing=1.45)


def plot_groups(groups: list[DistribGroup], out_path: str, title: str = "",
                log: bool = False, zoom: bool = True) -> list[str]:
    """
    One page per distribution: a multi-page PDF when `out_path` ends in .pdf,
    otherwise one file per distribution (<stem>_NN<suffix>).
    `title` is printed under each page heading.  Returns the files written.
    """
    try:
        import matplotlib
    except ImportError as exc:
        raise ImportError("plotting needs matplotlib (pip install matplotlib); "
                          "--list works without it") from exc
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
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
            pdf.infodict()["Title"] = title or "distrib_* position profiles"
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
                    help="print the distributions and their producers, then exit")
    ap.add_argument("--radio", nargs="+", default=None, metavar="PATTERN",
                    help="regex filter on radio producer names")
    ap.add_argument("--pick", nargs="+", default=None, metavar="SPEC",
                    help="PRODUCER_RE#SEL as in bki-plot --surface (index list "
                         "or @REGEX on surface, node, nuclide)")
    ap.add_argument("--axis", nargs="+", choices=DISTRIB_AXES, default=None,
                    help="only these axes (default: all present)")
    ap.add_argument("--range", nargs=2, type=float, default=None,
                    metavar=("LO", "HI"),
                    help="sampling range [cm] for IR entries without a stored "
                         "<axis>_range (older IRs)")
    ap.add_argument("--log", action="store_true", help="logarithmic y axis")
    ap.add_argument("--no-zoom", action="store_true",
                    help="no zoomed view on the effective range (by default it is "
                         "added below the plot when the effective range is "
                         "narrower than half the full range)")
    ap.add_argument("--title", default=None, help="figure title")
    args = ap.parse_args()

    ir_path = Path(args.ir_file)
    if not ir_path.is_file():
        sys.exit(f"[error] '{args.ir_file}' not found")
    ir = IRCollection.load(ir_path)

    for spec in args.pick or ():
        try:
            parse_pick_spec(spec)
        except (ValueError, re.error) as exc:
            sys.exit(f"[error] --pick: {exc}")
    try:
        vols = select_volumes(ir, kinds=["radio"],
                              kind_patterns={"radio": args.radio} if args.radio else None,
                              picks=args.pick)
    except ValueError as exc:
        sys.exit(f"[error] {exc}")

    groups = collect(vols, fallback=tuple(args.range) if args.range else None)
    if args.axis:
        groups = [g for g in groups if g.axis in args.axis]
    if not groups:
        sys.exit("[error] no radio volumes with a distrib_* formula "
                 "(after the selectors)")
    print(f"[info] {len(groups)} distribution(s) from "
          f"{sum(len(g.producers) for g in groups)} producer entr(ies)", file=sys.stderr)
    if any(not g.full for g in groups) and args.range is None:
        print("[warn] IR has no <axis>_range: drawing over the containment "
              "interval only; re-run bki-extract or pass --range",
              file=sys.stderr)

    if args.list:
        print_groups(groups)
        return

    out_path = args.out or str(ir_path.with_name(ir_path.stem + "_distrib.pdf"))
    title = args.title or (f"source: {Path(ir.source).name}"
                           if ir.source else ir_path.stem)
    try:
        plot_groups(groups, out_path, title=title, log=args.log, zoom=not args.no_zoom)
    except ImportError as exc:
        sys.exit(f"[error] {exc}")


if __name__ == "__main__":
    main()

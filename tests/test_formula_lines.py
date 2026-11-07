"""ROOT-free checks of the distrib_x/y/z formula text layout of bki-distrib
(duneradgeo.distrib.raw_formula_lines and formula_lines).

Run with pytest, or directly: python tests/test_formula_lines.py.  Needs no
PyROOT: the layout is string processing only.
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _runner import run_tests  # noqa: E402

from duneradgeo.distrib import (DistribGroup, formula_lines, group_label,  # noqa: E402
                                page_subtitle, page_title, raw_formula_lines)

# the upper-mesh and bulk K42 profiles of tests/test_tf1.py
MESH = "(7.2257e+4*TMath::Exp(-(x+320.)*(x+320.)/(2.*0.5*0.5)))+1."
BULK = ("10255.2+(41.3788*TMath::Abs(325.-x))-(0.135462*(325.-x)*(325.-x))"
        "+(0.000203379*TMath::Power(TMath::Abs(325.-x),3))"
        "-(1.14356e-007*TMath::Power(TMath::Abs(325.-x),4))")
# synthetic: 24 top-level terms with exponent signs; one long top-level term;
# signs in exponents and arguments, and a leading minus
POLY = "+".join(f"{1.5 * k:g}e-0{k % 9}*TMath::Power(x,{k})" for k in range(1, 25))
NEST = "(" + POLY + ")*TMath::Exp(-x/1e+3)"
SIGNS = "-3.2E-5*x*x-(x-1.e-2)+TMath::Gaus(x,-12.5,2.e+1)-TMath::Sqrt(TMath::Abs(x))"
FORMULAS = {"MESH": MESH, "BULK": BULK, "POLY": POLY, "NEST": NEST, "SIGNS": SIGNS}
N_TERMS = {"MESH": 2, "BULK": 5, "POLY": 24, "NEST": 1, "SIGNS": 4}
WIDTHS = (86, 60, 40, 32, 20)


def unwrap(lines: list, axis: str = "x") -> list:
    """The formula text of each line of raw_formula_lines, without the
    ``distrib_x: "`` lead, the continuation indent and the closing quote."""
    lead = f'distrib_{axis}: "'
    pad = " " * len(lead)
    assert lines[0].startswith(lead), lines[0]
    assert all(s.startswith(pad) for s in lines[1:]), lines
    assert lines[-1].endswith('"'), lines[-1]
    body = [lines[0][len(lead):]] + [s[len(pad):] for s in lines[1:]]
    body[-1] = body[-1][:-1]
    return body


def balanced(s: str) -> bool:
    depth = 0
    for ch in s:
        depth += (ch == "(") - (ch == ")")
        if depth < 0:
            return False
    return depth == 0


# ---------------------------------------------------------------------------
# raw_formula_lines
# ---------------------------------------------------------------------------

def test_raw_round_trip():
    """Joining the lines gives back the formula, at any width."""
    for name, f in FORMULAS.items():
        for width in WIDTHS + (120, 200):
            for axis in ("x", "z"):
                text = "".join(unwrap(raw_formula_lines(f, axis=axis, width=width), axis))
                assert text == f, (name, width, axis, text)


def test_raw_single_line():
    assert raw_formula_lines(MESH) == [f'distrib_x: "{MESH}"']
    assert raw_formula_lines("x", axis="y") == ['distrib_y: "x"']


def test_raw_wraps_at_terms():
    """When every term fits on a line, the lines break before a top-level
    '+'/'-' only, never inside parentheses or an exponent."""
    for name, width in (("BULK", 86), ("POLY", 86), ("POLY", 60), ("SIGNS", 60),
                        ("SIGNS", 40)):
        body = unwrap(raw_formula_lines(FORMULAS[name], width=width))
        assert len(body) > 1, (name, width)
        assert all(balanced(s) for s in body), (name, width, body)
        assert all(s[0] in "+-" for s in body[1:]), (name, width, body)
        assert not any(re.search(r"[0-9.][eE]$", s) for s in body[:-1]), (name, width, body)
    # BULK at the default width: its 5 terms on 3 lines
    body = unwrap(raw_formula_lines(BULK))
    assert len(body) == 3 and body[1].startswith("+(0.000203379*") \
        and body[2].startswith("-(1.14356e-007*"), body


def test_raw_width_bound():
    """No line is wider than `width`, which is raised to the lead plus 20
    characters of formula."""
    for name, f in FORMULAS.items():
        for width in WIDTHS:
            lines = raw_formula_lines(f, width=width)
            bound = max(width, len('distrib_x: "') + 20)
            assert max(map(len, lines)) <= bound, (name, width, lines)
    # a long term is cut to fit
    assert len(raw_formula_lines(NEST, width=40)) > len(NEST) // 40


def test_raw_width_bound_closing_quote():
    for f, width in ((POLY, 120), ("1" * 74, 86), ("2" * 148, 86)):
        lines = raw_formula_lines(f, width=width)
        assert max(map(len, lines)) <= width, (width, [len(s) for s in lines])


def test_raw_round_trip_spaces():
    f = "TMath::Gaus(x, 100, 1) + 2 * x"
    assert "".join(unwrap(raw_formula_lines(f))) == f


# ---------------------------------------------------------------------------
# formula_lines
# ---------------------------------------------------------------------------

_TERM_LINE = re.compile(r"^ {7}[+-] \S")     # the first line of a term after the first


def test_formula_lines_mesh():
    lines = formula_lines(MESH)
    assert len(lines) == 2, lines
    assert lines[0].startswith("f(x) = 7.2257e4·exp(")      # outer parentheses dropped
    assert "(x+320)" in lines[0] and "320." not in lines[0]
    assert lines[1] == "       + 1"


def test_formula_lines_bulk():
    lines = formula_lines(BULK)
    assert len(lines) == 5, lines
    assert lines[0] == "f(x) = 10255.2"
    assert [s[7] for s in lines[1:]] == ["+", "-", "+", "-"]
    assert "41.3788·|325-x|" in lines[1]
    assert "(325-x)·(325-x)" in lines[2]
    assert lines[3].endswith("|325-x|^3") and lines[4].endswith("|325-x|^4")
    assert "1.14356e-7" in lines[4]                          # e-007 -> e-7
    text = "".join(lines)
    assert not any(tok in text for tok in ("TMath", "Abs", "Power", "*"))


def test_formula_lines_layout():
    """One term per line (plus wrapped continuations, indented deeper), and
    wrapping only adds whitespace."""
    for name, f in FORMULAS.items():
        wide = formula_lines(f, width=1000)
        assert len(wide) == N_TERMS[name], (name, wide)
        for width in (70, 40, 20):
            lines = formula_lines(f, width=width)
            assert lines[0].startswith("f(x) = "), (name, lines)
            terms = [s for s in lines[1:] if _TERM_LINE.match(s)]
            assert len(terms) == N_TERMS[name] - 1, (name, width, lines)
            others = [s for s in lines[1:] if not _TERM_LINE.match(s)]
            assert all(s.startswith(" " * 9) for s in others), (name, width, others)
            assert "".join("".join(lines).split()) == "".join("".join(wide).split())
            # the first line of each term fits in `width` (at least 20 + indent)
            for s in [lines[0]] + terms:
                assert len(s) <= max(width, 29), (name, width, s)


def test_formula_lines_signs():
    assert formula_lines("-x+1") == ["f(x) = -x", "       + 1"]
    # exponent signs do not start a term
    assert formula_lines("2.5e-3*x-1e+5") == ["f(x) = 2.5e-3·x", "       - 1e5"]
    assert len(formula_lines(SIGNS)) == 4


def test_formula_lines_subtracted_sum():
    assert formula_lines("1-(x-2)") == ["f(x) = 1", "       - (x-2)"]
    assert formula_lines("-(x+1)") == ["f(x) = -(x+1)"]


def _group(producers, surfaces=("bulk",), eff=(-425.12, 415.96)):
    return DistribGroup(axis="x", formula="1", range_=(-425.12, 425.2), eff=eff,
                        formula_eff=eff, containment=0.99, npx=100,
                        surfaces=list(surfaces), sources=["volume_rand:volCryostat"],
                        producers={n: [nu, 1] for n, nu in producers})


MESH_PRODUCERS = [("K42From42ArGenInUpperMesh1x8x14", "K42"),
                  ("Rn222ChainFromPo218GenInUpperMesh1x8x14", "Po218"),
                  ("Rn222ChainFromBi210GenInUpperMesh1x8x14", "Bi210+Po210")]


def test_group_label():
    assert group_label(_group(MESH_PRODUCERS)) == "UpperMesh1x8x14"
    assert group_label(_group([("Rn222ChainPo218GenInLAr", "Po218")])) == \
        "Rn222ChainPo218GenInLAr"
    # a shared part that is not a word: fall back to the count
    assert group_label(_group([("FooA1", "K40"), ("BarB1", "K40")])) == "2 producers"
    # the common suffix starts mid-word: cut to the next CamelCase word
    assert group_label(_group([("AlphaCathodeGrid", "K40"),
                               ("BetaCathodeGrid", "U238")])) == "CathodeGrid"


def test_page_title():
    g = _group(MESH_PRODUCERS)
    assert page_title(g) == ("UpperMesh1x8x14 · x profile, 3 producers "
                             "(K42, Po218, Bi210+Po210)")
    assert page_title(g, 1) == "UpperMesh1x8x14 · x profile, 3 producers (K42, …)"
    assert page_title(_group([("Rn222ChainPo218GenInLAr", "Po218")])) == \
        "Rn222ChainPo218GenInLAr · x profile (Po218)"
    assert page_title(_group([("FooA1", "K40"), ("BarB1", "K40")])) == \
        "x profile, 2 producers (K40)"


def test_page_subtitle():
    layer = _group(MESH_PRODUCERS, surfaces=["plane_x@-320.917"],
                   eff=(-323.082, -318.753))
    assert page_subtitle(layer, "f.fcl") == ("layer at x ≈ -320.9 cm, 99% in 4.3 cm"
                                             "  ·  volume_rand:volCryostat  ·  f.fcl")
    assert page_subtitle(_group(MESH_PRODUCERS)) == \
        "bulk, 99% in [-425.1, 416.0] cm  ·  volume_rand:volCryostat"


if __name__ == "__main__":
    run_tests(globals())

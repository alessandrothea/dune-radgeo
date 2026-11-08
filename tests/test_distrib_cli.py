"""Checks of the command line of bki-distrib (duneradgeo.distrib, click + rich):
option syntax, error exits, --help and the --list blocks.

Run with pytest, or directly: python tests/test_distrib_cli.py.  Needs click
and rich (skipped without them).  The option syntax and the error exits need
neither ROOT nor matplotlib; the --list output is checked only when PyROOT is
importable.  IRs are small temporary files.
"""
import importlib.util
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
if importlib.util.find_spec("click") is None or importlib.util.find_spec("rich") is None:
    try:
        import pytest
        pytest.skip("click and rich are needed", allow_module_level=True)
    except ImportError:
        print("skipped: click and rich are needed")
        sys.exit(0)

from click.testing import CliRunner  # noqa: E402
from _runner import run_tests  # noqa: E402

from duneradgeo import distrib  # noqa: E402
from duneradgeo._fcl import distrib_keys  # noqa: E402
from duneradgeo.ir import IRCollection, VolumeEntry  # noqa: E402

HAVE_ROOT = importlib.util.find_spec("ROOT") is not None
LONG = "+".join(f"{k}.5*TMath::Gaus(x,{-300 + k},{k}.)" for k in range(1, 25))


def _radio(name, nuclide, axis, formula, surface="bulk"):
    k = distrib_keys(axis)
    meta = {"nuclide": nuclide, "surface": surface, "source": "volume_rand:volCryostat",
            k["formula"]: formula, k["range"]: [-100.0, 100.0],
            k["eff"]: [-50.0, 50.0], k["formula_eff"]: [-60.0, 60.0],
            "distrib_eff_containment": 0.99, "distrib_eff_npx": 100}
    return VolumeEntry("radio", name, np.zeros(3), np.ones(3), np.eye(3), meta)


def make_ir(formula="TMath::Gaus(x,0,20)", axis="x") -> IRCollection:
    vols = [_radio("Ar39GenInLAr", "Ar39", axis, formula),
            _radio("K42GenInLAr", "K42", axis, formula),
            _radio("Po218GenInMesh", "Po218", "z", "1.", "plane_z@5")]
    return IRCollection(source="a.fcl", detector="vd", volumes=vols)


def run(*args):
    res = CliRunner().invoke(distrib.cli, [str(a) for a in args])
    return res.exit_code, res.output


def with_ir(fn, **kw):
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "ir.json")
        make_ir(**kw).save(path)
        fn(path, d)


def captured(*args):
    """Run with distrib._run replaced by a recorder of its arguments."""
    seen = []
    orig = distrib._run
    distrib._run = lambda *a: seen.append(a)
    try:
        code, out = run(*args)
    finally:
        distrib._run = orig
    return code, out, (seen[0] if seen else None)


# --- option syntax -----------------------------------------------------------

def test_defaults():
    code, out, a = captured("ir.json")
    assert code == 0, out
    # ir_file, out, list_, radio, pick, axis, npx, log, no_zoom, title
    assert a == ("ir.json", None, False, (), (), (), None, False, False, None)


def test_every_option():
    code, out, a = captured("ir.json", "-o", "o.png", "-l", "--radio", "Mesh", "Cath",
                            "--pick", "^A#0", "^B#@x", "--axis", "x", "z",
                            "--npx", "1000", "--log", "--no-zoom", "--title", "my.fcl")
    assert code == 0, out
    assert a == ("ir.json", "o.png", True, ("Mesh", "Cath"), ("^A#0", "^B#@x"),
                 ("x", "z"), 1000, True, True, "my.fcl")
    code, out, a = captured("--out", "o.pdf", "--list", "ir.json")
    assert code == 0 and a[:3] == ("ir.json", "o.pdf", True), out


def test_variadic_stops_at_the_next_option():
    code, out, a = captured("ir.json", "--radio", "A|B", "C", "--log")
    assert code == 0 and a[3] == ("A|B", "C") and a[7] is True, out
    code, out, a = captured("ir.json", "--axis", "y", "--pick", "p#1")
    assert code == 0 and a[5] == ("y",) and a[4] == ("p#1",), out


def test_negative_npx():
    code, out, a = captured("ir.json", "--npx", "-5")
    assert code == 0 and a[6] == -5, out


def test_usage_errors_exit_2():
    for args in (["ir.json", "--axis", "w"], ["ir.json", "--axis", "x", "q"],
                 ["ir.json", "--npx", "many"], ["ir.json", "--nope"], []):
        code, out = run(*args)
        assert code == 2, (args, code, out)
    assert "Invalid value" in run("ir.json", "--axis", "w")[1]


# --- help --------------------------------------------------------------------

def test_help():
    for flag in ("--help", "-h"):
        code, out = run(flag)
        assert code == 0
        for text in ("bki-distrib", "IR_FILE", "--out", "-o", "--list", "-l",
                     "--radio PATTERN...", "--pick SPEC...", "--axis {x,y,z}...",
                     "--npx N", "--log", "--no-zoom", "--title",
                     "<stem>_NN<suffix>", "needs ROOT, not matplotlib",
                     "PRODUCER_RE#SEL", "distrib_eff_npx", "narrower than half",
                     "bki-distrib detector.json --list", "Needs PyROOT"):
            assert text in " ".join(out.split()) or text in out, text


# --- error exits (before any ROOT use) ----------------------------------------

def test_missing_file():
    code, out = run("/nonexistent/ir.json")
    assert code == 1 and "[error] '/nonexistent/ir.json' not found" in out


def test_bad_pick():
    def go(path, d):
        code, out = run(path, "--pick", "nohash")
        assert code == 1 and out.startswith("[error] ") and "nohash" in out, out
    with_ir(go)


def test_no_distribution_after_selectors():
    def go(path, d):
        code, out = run(path, "--radio", "NoSuchProducer")
        assert code == 1 and "[error] " in out, out
        code, out = run(path, "--axis", "y")
        assert code == 1 and "no radio volumes with a distrib_* formula" in out, out
    with_ir(go)


def test_missing_distrib_keys():
    with tempfile.TemporaryDirectory() as d:
        ir = make_ir()
        del ir.volumes[0].meta["distrib_eff_npx"]
        path = os.path.join(d, "ir.json")
        ir.save(path)
        code, out = run(path)
        assert code == 1 and "[error] Ar39GenInLAr: IR entry lacks" in out, out


# --- --list (needs ROOT) -------------------------------------------------------

def test_list_blocks():
    if not HAVE_ROOT:
        return

    def go(path, d):
        code, out = run(path, "--list")
        assert code == 0, out
        assert "[0] " in out and "[1] " in out
        assert "distrib_x" in out and "distrib_z" in out
        assert "range [-100.000, 100.000]" in out
        assert "containment 0.99" in out and "npx 100" in out
        assert "eff_range (generated) [-50.000, 50.000]" in out
        assert "formula eff range" in out and "[-60.000, 60.000]" in out
        assert "(measured " in out
        assert "surface: bulk" in out and "surface: plane_z@5" in out
        assert "formula: TMath::Gaus(x,0,20)" in out
        for name in ("Ar39GenInLAr", "K42GenInLAr", "Po218GenInMesh"):
            assert name in out
        assert "producer" in out and "nuclide" in out and "placements" in out
        assert "Ar39" in out and "K42" in out
        code, out = run(path, "--list", "--axis", "z")
        assert "distrib_x" not in out and "Po218GenInMesh" in out
        code, out = run(path, "--list", "--radio", "K42")
        assert "K42GenInLAr" in out and "Ar39GenInLAr" not in out
        code, out = run(path, "--list", "--npx", "1000")
        assert "npx 1000" in out and "[recomputed for --npx]" in out
    with_ir(go)


def test_list_formula_is_never_truncated():
    if not HAVE_ROOT:
        return

    def go(path, d):
        code, out = run(path, "--list", "--axis", "x")
        assert code == 0, out
        assert "…" not in out.split("formula:")[1]
        # folded over lines: stripping the panel borders and newlines gives it back
        body = out.split("formula:", 1)[1].split("producer", 1)[0]
        flat = "".join(ch for ch in body if ch not in "│\n ")
        assert LONG.replace(" ", "") in flat
    with_ir(go, formula=LONG)


def test_unparsable_formula_is_skipped():
    if not HAVE_ROOT:
        return

    def go(path, d):
        code, out = run(path, "--list", "--axis", "x")
        assert code == 1 and "[warn] distrib_x" in out and "skipped" in out
        assert "[error] no distribution could be evaluated" in out, out
    with_ir(go, formula="TMath::Gaus(x,0,")


def test_print_groups_without_root():
    """print_groups still works on a group ROOT cannot evaluate."""
    g = distrib.collect(make_ir().volumes)[0]
    import io
    from unittest import mock
    from rich.console import Console
    buf = io.StringIO()
    with mock.patch.object(distrib._cli, "console",
                           lambda: Console(file=buf, width=120, highlight=False)), \
            mock.patch.object(distrib, "fractions", side_effect=ValueError("no ROOT")):
        distrib.print_groups([g])
    out = buf.getvalue()
    assert "[0] " in out and "distrib_x" in out and "[cannot evaluate: no ROOT]" in out
    assert "formula: TMath::Gaus(x,0,20)" in out and "Ar39GenInLAr" in out


if __name__ == "__main__":
    run_tests(globals())

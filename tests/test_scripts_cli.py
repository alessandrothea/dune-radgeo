"""Checks of the click command lines of scripts/compare_ir.py,
scripts/compare_backends.py and scripts/run_notebook.py.

compare_ir is run on two small in-memory IRs; compare_backends (needs TGeo)
and run_notebook (needs jupyter) are only checked for --help and for a clean
error on a missing input.

Run with pytest, or directly: python tests/test_scripts_cli.py.  Needs click
and rich (skipped without them).
"""
import importlib.util
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
if importlib.util.find_spec("click") is None or importlib.util.find_spec("rich") is None:
    try:
        import pytest
        pytest.skip("click and rich are needed", allow_module_level=True)
    except ImportError:
        print("skipped: click and rich are needed")
        sys.exit(0)

import numpy as np  # noqa: E402
from click.testing import CliRunner  # noqa: E402
from _runner import run_tests  # noqa: E402

from duneradgeo.ir import IRCollection, VolumeEntry  # noqa: E402

SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")


def load_script(name):
    spec = importlib.util.spec_from_file_location(f"_script_{name}",
                                                  os.path.join(SCRIPTS, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def entry(kind, name, x, part=None):
    meta = {"part": part} if part else {}
    return VolumeEntry(kind=kind, name=name, pos=np.array([x, 0.0, 0.0]),
                       ext=np.array([4.0, 10.0, 10.0]), rot=np.eye(3), meta=meta)


def make_ir(shift=0.0):
    vols = [entry("anode", "volAnodePlate", 325.0 + 20 * i + (shift if i == 1 else 0))
            for i in range(3)]
    vols.append(entry("cathode", "volCathodeGrid", -328.0))
    return IRCollection(source="x.gdml", detector="vd", volumes=vols)


def write_irs(tmp, shift):
    a, b = os.path.join(tmp, "ref.json"), os.path.join(tmp, "new.json")
    make_ir().save(a)
    make_ir(shift).save(b)
    return a, b


def run(script, *args):
    mod = load_script(script)
    return CliRunner().invoke(mod.cli, list(args))


# --- compare_ir -------------------------------------------------------------

def test_compare_ir_identical():
    with tempfile.TemporaryDirectory() as tmp:
        a, b = write_irs(tmp, 0.0)
        res = run("compare_ir", a, b, "-v")
    assert res.exit_code == 0, res.output
    assert "ok    anode" in res.output and "FAIL" not in res.output
    assert "2/2 groups agree" in res.output
    quiet = None
    with tempfile.TemporaryDirectory() as tmp:
        a, b = write_irs(tmp, 0.0)
        quiet = run("compare_ir", a, b)
    assert quiet.exit_code == 0 and "ok    anode" not in quiet.output
    assert "2/2 groups agree" in quiet.output


def test_compare_ir_changed_entry():
    with tempfile.TemporaryDirectory() as tmp:
        a, b = write_irs(tmp, 0.5)
        res = run("compare_ir", a, b)
        assert res.exit_code != 0
        assert "FAIL  anode" in res.output and "1 pos" in res.output
        assert "1/2 groups agree" in res.output
        # within a large enough tolerance the same pair agrees
        res = run("compare_ir", a, b, "--tol", "1")
        assert res.exit_code == 0 and "2/2 groups agree" in res.output


def test_compare_ir_only_common_and_kinds():
    with tempfile.TemporaryDirectory() as tmp:
        a, b = write_irs(tmp, 0.0)
        ir = make_ir()
        ir.volumes.append(entry("opdet", "volOpDetSensitive", 1.0))
        ir.save(b)
        res = run("compare_ir", a, b)
        assert res.exit_code != 0 and "FAIL  opdet" in res.output and "only in new" in res.output
        res = run("compare_ir", a, b, "--only-common")
        assert res.exit_code == 0 and "skip  opdet" in res.output
        res = run("compare_ir", a, b, "--kinds", "anode", "cathode")
        assert res.exit_code == 0 and "2/2 groups agree" in res.output


def test_compare_ir_ignore_rot():
    with tempfile.TemporaryDirectory() as tmp:
        a, b = write_irs(tmp, 0.0)
        ir = make_ir()
        ir.volumes[0].rot = np.array([[0.0, -1, 0], [1, 0, 0], [0, 0, 1]])
        ir.save(b)
        res = run("compare_ir", a, b)
        assert res.exit_code != 0 and "FAIL  anode" in res.output
        res = run("compare_ir", a, b, "--ignore-rot", "^volAnode", "^volOther", "-v")
        assert res.exit_code == 0 and "ignored" in res.output
        res = run("compare_ir", a, b, "--rot-tol", "2")
        assert res.exit_code == 0


def test_compare_ir_missing_file():
    with tempfile.TemporaryDirectory() as tmp:
        a, _ = write_irs(tmp, 0.0)
        res = run("compare_ir", a, os.path.join(tmp, "nope.json"))
    assert res.exit_code == 1 and "[error]" in res.output and "not found" in res.output
    assert run("compare_ir", "only_one.json").exit_code == 2


# --- --help of the three scripts ----------------------------------------------

def _flat(text):
    return " ".join(text.split())


def test_help_options():
    wanted = {
        "compare_ir": ["REF", "NEW", "--tol", "--rot-tol", "--kinds KIND...",
                       "--ignore-rot NAME_RE...", "--only-common", "-v, --verbose",
                       "compare_ir.py ref.json new.json --only-common -v", "-h, --help"],
        "compare_backends": ["GDML...", "--patterns RE...", "--tol", "--rot-tol",
                             "-v, --verbose", "[default: 0.001]", "-h, --help",
                             "compare_backends.py geometry.gdml [more.gdml ...]"],
        "run_notebook": ["IR_FILE", "-o, --out", "--timeout", "[default: 120]",
                         "-h, --help", "run_notebook.py detector.json --timeout 300"],
    }
    for script, items in wanted.items():
        for flag in ("--help", "-h"):
            res = run(script, flag)
            assert res.exit_code == 0, (script, res.output)
            text = _flat(res.output)
            for item in items:
                assert item in text, (script, item, res.output)


def test_missing_input_is_a_clean_error():
    res = run("compare_backends", "/nonexistent/x.gdml")
    assert res.exit_code == 1 and res.output.startswith("[error]") and "not found" in res.output
    res = run("run_notebook", "/nonexistent/x.json")
    assert res.exit_code == 1 and res.output.startswith("[error] IR file not found")
    assert run("compare_backends").exit_code == 2
    assert run("run_notebook").exit_code == 2


def test_option_values_are_checked():
    assert run("compare_backends", "x.gdml", "--tol", "abc").exit_code == 2
    assert run("run_notebook", "x.json", "--timeout", "1.5").exit_code == 2
    assert run("compare_ir", "a", "b", "--rot-tol", "z").exit_code == 2


if __name__ == "__main__":
    run_tests(globals())

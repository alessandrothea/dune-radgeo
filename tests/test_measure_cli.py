"""Checks of the bki-measure command line (duneradgeo.measure) on a small
in-memory IR written to a temporary file; needs numpy, click and rich
(skipped without them).  Run with pytest, or: python tests/test_measure_cli.py.
"""
import importlib.util
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
if any(importlib.util.find_spec(m) is None for m in ("click", "rich", "numpy")):
    try:
        import pytest
        pytest.skip("click, rich and numpy are needed", allow_module_level=True)
    except ImportError:
        print("skipped: click, rich and numpy are needed")
        sys.exit(0)

import numpy as np  # noqa: E402
from click.testing import CliRunner  # noqa: E402
from _runner import run_tests  # noqa: E402

from duneradgeo import measure  # noqa: E402
from duneradgeo.ir import IRCollection, VolumeEntry  # noqa: E402


def box(kind, name, lo, hi, **meta) -> VolumeEntry:
    lo, hi = np.array(lo, float), np.array(hi, float)
    return VolumeEntry(kind=kind, name=name, pos=(lo + hi) / 2, ext=hi - lo,
                       rot=np.eye(3), meta=meta)


def write_ir(d: str) -> str:
    vols = [box("cathode", "volCathode", (-2, -50, 0), (2, 50, 100)),
            box("anode", "volAnodePlate", (300, -50, 0), (304, 50, 100)),
            box("anode", "volAnodePlate", (-304, -50, 0), (-300, 50, 100)),
            box("struct", "volCryostat", (-400, -300, -50), (400, 300, 200))]
    for i in range(3):
        vols.append(box("fieldshaper", "FieldShaperLongtube",
                        (-100, 200 + 5 * i, 0), (100, 202 + 5 * i, 100)))
    path = os.path.join(d, "ir.json")
    IRCollection(source="x.gdml", detector="vd", volumes=vols).save(path)
    return path


def run(*args):
    res = CliRunner().invoke(measure.cli, list(args))
    return res.exit_code, res.output


def test_list_needs_no_ir():
    code, out = run("--list")
    assert code == 0
    for name in ("drift-gaps", "fc-cryostat-gap", "anode-cryostat-gap", "opdet-fc-gap"):
        assert name in out
    assert "cathode" in out and "help" in out


def test_named_measurement_and_detail():
    with tempfile.TemporaryDirectory() as d:
        ir = write_ir(d)
        code, out = run(ir, "drift-gaps")
        assert code == 0 and "drift-gaps" in out and "top" in out and "bottom" in out
        assert "[cm]" in out and "gap" in out
        assert "centre_dist" not in out            # the per-row table needs --detail
        code, out = run(ir, "drift-gaps", "-d")
        assert code == 0 and "centre_dist" in out and "volAnodePlate#0" in out
        code, out2 = run(ir, "drift-gaps", "--detail")
        assert out2 == out


def test_default_runs_all():
    with tempfile.TemporaryDirectory() as d:
        code, out = run(write_ir(d))
        assert code == 0
        for name in ("drift-gaps", "fc-cryostat-gap", "anode-cryostat-gap", "opdet-fc-gap"):
            assert name in out


def test_json_file():
    with tempfile.TemporaryDirectory() as d:
        ir, out = write_ir(d), os.path.join(d, "out.json")
        code, _ = run(ir, "drift-gaps", "--json", out)
        assert code == 0
        data = json.load(open(out))
        assert data["detector"] == "vd" and data["source"] == "x.gdml"
        assert [r["name"] for r in data["results"]] == ["drift-gaps"]
        assert len(data["results"][0]["rows"]) == 2
        assert data["results"][0]["values"]["top_gap"] is not None


def test_errors():
    with tempfile.TemporaryDirectory() as d:
        ir = write_ir(d)
        code, out = run(ir, "nope")
        assert code == 1 and "[error] unknown measurement(s): nope; available: " in out
        assert "drift-gaps" in out
        code, out = run(os.path.join(d, "missing.json"))
        assert code == 1 and "[error]" in out and "not found" in out
        code, out = run()
        assert code == 2 and "IR_FILE is required" in out
        code, out = run(ir, "--bogus")
        assert code == 2


def test_help():
    code, out = run("--help")
    assert code == 0
    for s in ("--list", "--detail", "-d", "--json FILE", "-h", "IR_FILE", "MEASUREMENT",
              "Lengths are in cm", "bki-measure detector.json drift-gaps -d"):
        assert s in out
    assert run("-h")[1] == out


if __name__ == "__main__":
    run_tests(globals())

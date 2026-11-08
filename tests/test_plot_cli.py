"""Checks of the radgeo plot command line (duneradgeo.plot): option syntax,
--list / --list-detail tables and the error exits.

Run with pytest, or directly: python tests/test_plot_cli.py.  The IR is a
small in-memory IRCollection saved to a temporary file; plots are written to a
temporary directory with the Agg backend and skipped when matplotlib is missing.
"""
import importlib.util
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
for _mod in ("click", "rich"):
    if importlib.util.find_spec(_mod) is None:
        try:
            import pytest
            pytest.skip(f"{_mod} is not installed", allow_module_level=True)
        except ImportError:
            print(f"skipped: {_mod} is not installed")
            sys.exit(0)

from click.testing import CliRunner  # noqa: E402
from _runner import run_tests  # noqa: E402

from duneradgeo import plot  # noqa: E402
from duneradgeo.ir import IRCollection, VolumeEntry  # noqa: E402

HAVE_MPL = importlib.util.find_spec("matplotlib") is not None
_CRYO = "volDetEnclosure:0/volCryostat:0"


def vol(kind, name, x=0.0, **meta) -> VolumeEntry:
    return VolumeEntry(kind=kind, name=name, pos=np.array([x, 0.0, 0.0]),
                       ext=np.ones(3), rot=np.eye(3), meta=meta)


def make_ir() -> IRCollection:
    vols = [
        vol("anode", "volAnodePlate", -5.0, path=f"{_CRYO}/volAnodePlate:0"),
        vol("anode", "volAnodePlate", 5.0, path=f"{_CRYO}/volAnodePlate:1"),
        vol("cathode", "volCathodeGrid", path=f"{_CRYO}/volCathodeGrid:0"),
        vol("radio", "K42Gen", 1.0, nuclide="K42", surface="cathode", index=0),
        vol("radio", "K42Gen", 2.0, nuclide="K42", surface="anode", index=1),
        vol("radio", "Ar39GenInLAr", 3.0, nuclide="Ar39", surface="bulk",
            source="volume_rand:volCryostat", index=0),
    ]
    return IRCollection(source="test.gdml", detector="vd", volumes=vols)


def run(tmp, *args):
    path = os.path.join(tmp, "ir.json")
    if not os.path.exists(path):
        make_ir().save(path)
    res = CliRunner().invoke(plot.cli, [path, *args], env={"COLUMNS": "200"})
    return res.exit_code, res.output


def test_list():
    with tempfile.TemporaryDirectory() as tmp:
        code, out = run(tmp, "--list")
        assert code == 0, out
        assert "anode  (2 placements)" in out and "volAnodePlate" in out
        assert "radio  (3 placements)" in out
        assert "volCathodeGrid" in out and "Ar39GenInLAr" in out
        # one row per radio surface, with producer name, nuclide, node/source
        assert "cathode" in out and "bulk" in out and "volume_rand:volCryostat" in out
        code, out = run(tmp, "--list", "--kind", "anode", "cathode", "--name", "Anode")
        assert code == 0 and "volAnodePlate" in out
        assert "volCathodeGrid" not in out and "K42Gen" not in out


def test_list_detail():
    with tempfile.TemporaryDirectory() as tmp:
        code, out = run(tmp, "-L", "--anode", "Plate", "--no-radio")
        assert code == 0, out
        assert "x_min" in out and "z_max" in out and "path" in out
        assert f"{_CRYO}/volAnodePlate:1" in out
        assert "-5.500" in out and "-4.500" in out        # x range of the first anode
        assert "volCathodeGrid" in out
        code, out = run(tmp, "-L", "--highlight", "a=Anode", "--no-radio")
        assert code == 0 and "*1" in out


def test_errors():
    with tempfile.TemporaryDirectory() as tmp:
        code, out = run(tmp, "--list", "--pick", "volAnode")
        assert code == 1 and "[error]" in out and "Traceback" not in out, out
        code, out = run(tmp, "--list", "--name", "(")
        assert code == 1 and "[error]" in out and "Traceback" not in out, out
        res = CliRunner().invoke(plot.cli, [os.path.join(tmp, "missing.json")])
        assert res.exit_code == 1 and "[error]" in res.output and "not found" in res.output
        res = CliRunner().invoke(plot.cli, [__file__])
        assert res.exit_code == 1 and "not a JSON IR file" in res.output
        code, out = run(tmp, "--orientation", "sideways")
        assert code == 2


def test_help():
    res = CliRunner().invoke(plot.cli, ["-h"])
    assert res.exit_code == 0
    for s in ("--list-detail", "--view", "--no-proportional", "--zoom-highlight",
              "--highlight", "--pick", "hd|vd|auto", "radgeo plot detector.json"):
        assert s in res.output, s


def test_plots():
    if not HAVE_MPL:
        return
    with tempfile.TemporaryDirectory() as tmp:
        png = os.path.join(tmp, "a.png")
        code, out = run(tmp, "--out", png, "--view", "20", "-60")
        assert code == 0, out
        assert os.path.getsize(png) > 1000 and "plot written to" in out
        png = os.path.join(tmp, "b.png")
        code, out = run(tmp, "--highlight", "a=^volAnode", "--kind", "anode", "cathode",
                        "--out", png, "--dim", "0.1", "--no-proportional",
                        "--orientation", "vd")
        assert code == 0 and os.path.exists(png), out
        png = os.path.join(tmp, "c.png")
        code, out = run(tmp, "--highlight", "z=^nothing", "--out", png)
        assert code == 0 and "matches no selected volume" in out, out
        png = os.path.join(tmp, "d.png")
        code, out = run(tmp, "--pick", "volAnode#@:1", "--radio", "K42", "--out", png)
        assert code == 0 and os.path.exists(png), out


if __name__ == "__main__":
    run_tests(globals())

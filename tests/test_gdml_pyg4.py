"""Checks of the pyg4ometry GDML backend (optional extra ``pyg4ometry``).

Its tables must equal those of the built-in parser on the example
geometries (examples/gdml), compared by scripts/compare_backends.py; on a
small GDML with constants and expressions, which the built-in parser cannot
read, it must evaluate them; and without pyg4ometry it fails with a clean
error.

Run with pytest, or directly: python tests/test_gdml_pyg4.py.  The tests that
need pyg4ometry are skipped when it is not installed.
"""
import glob
import importlib.util
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _runner import run_tests  # noqa: E402

from duneradgeo import _gdml_pyg4  # noqa: E402
from duneradgeo._gdml import load_geometry  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = sorted(glob.glob(os.path.join(ROOT, "examples", "gdml", "*.gdml")))

# constants and expressions in lengths and positions; a cone, which the
# built-in parser does not know; a union with a rotated second operand
_GDML = """<?xml version="1.0"?>
<gdml>
 <define>
  <constant name="Lx" value="200"/>
  <position name="p0" x="Lx/2" y="0" z="0" unit="mm"/>
  <rotation name="r90" z="90" unit="deg"/>
 </define>
 <materials><element name="H" formula="H" Z="1"><atom value="1"/></element>
  <material name="m" Z="1"><D value="1"/><atom value="1"/></material></materials>
 <solids>
  <box name="wb" x="10*Lx" y="10*Lx" z="10*Lx" lunit="mm"/>
  <box name="b1" x="Lx" y="Lx" z="Lx" lunit="mm"/>
  <box name="b2" x="Lx/4" y="Lx/2" z="Lx/4" lunit="mm"/>
  <cone name="c" rmax1="10" rmax2="20" z="40" deltaphi="360" aunit="deg" lunit="mm"/>
  <union name="u"><first ref="b1"/><second ref="b2"/>
   <position x="Lx" y="0" z="0" unit="mm"/><rotationref ref="r90"/></union>
 </solids>
 <structure>
  <volume name="vu"><materialref ref="m"/><solidref ref="u"/></volume>
  <volume name="vc"><materialref ref="m"/><solidref ref="c"/></volume>
  <volume name="vw"><materialref ref="m"/><solidref ref="wb"/>
   <physvol name="pv"><volumeref ref="vu"/><positionref ref="p0"/></physvol>
   <physvol><volumeref ref="vc"/></physvol></volume>
 </structure>
 <setup name="Default" version="1.0"><world ref="vw"/></setup>
</gdml>
"""


def have_pyg4ometry() -> bool:
    try:
        _gdml_pyg4._import_pyg4ometry()
    except _gdml_pyg4.Pyg4ometryUnavailable:
        return False
    return True


def skip_without_pyg4ometry():
    if not have_pyg4ometry():
        try:
            import pytest
            pytest.skip("pyg4ometry is not installed")
        except ImportError:
            print("  skipped: pyg4ometry is not installed")
            return True
    return False


def same(a, b, tol=1e-6) -> bool:
    return np.allclose(a, b, atol=tol)


def load_compare_backends():
    spec = importlib.util.spec_from_file_location(
        "_script_compare_backends", os.path.join(ROOT, "scripts", "compare_backends.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_agrees_with_the_builtin_parser():
    """scripts/compare_backends.py: volumes, solid tables and every profile
    placement (paths, positions, rotations, physvol names, box parts)."""
    if skip_without_pyg4ometry():
        return
    assert EXAMPLES
    compare = load_compare_backends().compare
    for path in EXAMPLES:
        assert compare(path, None, 1e-3, 1e-6, False, ("python", "pyg4ometry")) == 0, path


def test_constants_expressions_and_other_solids():
    if skip_without_pyg4ometry():
        return
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "g.gdml")
        with open(path, "w") as f:
            f.write(_GDML)
        g = load_geometry(path, backend="pyg4ometry")
    assert same(g._boxes["wb"], [200, 200, 200])               # 10*Lx mm, in cm
    (name, pos, rot, path_, pv, pvpath), = g.find_positions(["^vu$"], with_path=True)
    assert same(pos, [10, 0, 0]) and pv == "pv"                # Lx/2 mm
    # the unnamed physvol of vc has no name; the cone has no primitive
    # decomposition but a bounding box (mesh: rmax 20 mm, z 40 mm)
    (_, _, _, _, pv_c, _), = g.find_positions(["^vc$"], with_path=True)
    assert pv_c is None
    assert g.volume_parts("vc") == []
    c, e = g.volume_bbox("vc")
    assert same(c, [0, 0, 0], 1e-3) and same(e, [4, 4, 4], 1e-3)
    op, first, second, bpos, brot, _, _ = g._bools["u"]
    assert (op, first, second) == ("union", "b1", "b2")
    assert same(bpos, [20, 0, 0]) and same(abs(brot), [[0, 1, 0], [1, 0, 0], [0, 0, 1]])


def test_missing_pyg4ometry_is_a_clean_error():
    from click.testing import CliRunner
    from duneradgeo import extract
    saved = sys.modules.get("pyg4ometry", ...), _gdml_pyg4._pg
    sys.modules["pyg4ometry"] = None
    _gdml_pyg4._pg = None
    try:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "g.gdml")
            with open(path, "w") as f:
                f.write(_GDML)
            res = CliRunner().invoke(extract.cli, [
                path, "-o", os.path.join(tmp, "ir.json"),
                "--gdml-backend", "pyg4ometry"])
    finally:
        if saved[0] is ...:
            sys.modules.pop("pyg4ometry", None)
        else:
            sys.modules["pyg4ometry"] = saved[0]
        _gdml_pyg4._pg = saved[1]
    assert res.exit_code == 1, res.output
    assert "[error] " in res.output and "pyg4ometry" in res.output
    assert "Traceback" not in res.output
    assert isinstance(res.exception, SystemExit)


if __name__ == "__main__":
    run_tests(globals())

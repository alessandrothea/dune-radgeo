"""Checks of the radgeo extract command line (click + rich) on a tiny GDML of boxes
and small FHiCL files; no ROOT needed (the ROOT-missing error is forced).

Run with pytest, or directly: python tests/test_extract_cli.py.  Needs click
and rich (skipped without them).
"""
import contextlib
import importlib.util
import json
import os
import sys
import tempfile

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

from duneradgeo import _root, extract  # noqa: E402

# volEnclosureTPC (placed once, at z = +1000 mm in the world) holds an anode
# plate, two Arapuca, a cathode grid and a volExtra; volSteelShell is a passive
# volume of the world.  volAnodePlate is the hint of the "vd-perl" profile.
_GDML = """<?xml version="1.0"?>
<gdml>
 <materials>
  <element name="Al" formula="Al" Z="13"><atom value="26.98"/></element>
  <material name="Alu" Z="13"><D value="2.699" unit="g/cm3"/><atom value="26.98"/></material>
  <material name="Steel" Z="26"><D value="7850" unit="kg/m3"/><atom value="55.8"/></material>
  <material name="LAr" Z="18"><D value="1.396"/><atom value="39.95"/></material>
  <material name="Odd" Z="1"><D value="3" unit="furlong"/><atom value="1"/></material>
 </materials>
 <solids>
  <box name="World" x="100000" y="100000" z="100000" lunit="mm"/>
  <box name="Enc" x="4000" y="4000" z="4000" lunit="mm"/>
  <box name="Anode" x="20" y="1000" z="1000" lunit="mm"/>
  <box name="Arapuca" x="10" y="100" z="200" lunit="mm"/>
  <box name="Cathode" x="30" y="2000" z="2000" lunit="mm"/>
  <box name="Extra" x="50" y="60" z="70" lunit="mm"/>
  <box name="Shell" x="8000" y="8000" z="8000" lunit="mm"/>
 </solids>
 <structure>
  <volume name="volAnodePlate"><materialref ref="Alu"/><solidref ref="Anode"/></volume>
  <volume name="volArapuca_0"><solidref ref="Arapuca"/></volume>
  <volume name="volArapuca_1"><solidref ref="Arapuca"/></volume>
  <volume name="volCathodeGrid"><materialref ref="G4_Cu"/><solidref ref="Cathode"/></volume>
  <volume name="volExtra"><materialref ref="Odd"/><solidref ref="Extra"/></volume>
  <volume name="volSteelShell"><materialref ref="Steel"/><solidref ref="Shell"/></volume>
  <volume name="volEnclosureTPC"><materialref ref="LAr"/><solidref ref="Enc"/>
   <physvol><volumeref ref="volAnodePlate"/>
    <position name="pa" x="-1000" y="0" z="0" unit="mm"/></physvol>
   <physvol><volumeref ref="volArapuca_0"/>
    <position name="po0" x="-900" y="0" z="100" unit="mm"/></physvol>
   <physvol><volumeref ref="volArapuca_1"/>
    <position name="po1" x="-900" y="0" z="-100" unit="mm"/></physvol>
   <physvol><volumeref ref="volCathodeGrid"/>
    <position name="pc" x="0" y="0" z="0" unit="mm"/></physvol>
   <physvol><volumeref ref="volExtra"/>
    <position name="pe" x="500" y="0" z="0" unit="mm"/></physvol>
  </volume>
  <volume name="volWorld"><solidref ref="World"/>
   <physvol><volumeref ref="volEnclosureTPC"/>
    <position name="pt" x="0" y="0" z="1000" unit="mm"/></physvol>
   <physvol><volumeref ref="volSteelShell"/>
    <position name="ps" x="0" y="0" z="0" unit="mm"/></physvol>
  </volume>
 </structure>
 <setup name="Default" version="1.0"><world ref="volWorld"/></setup>
</gdml>
"""

# cm; the first producer is an explicit box, no distrib_*
_FCL = """
bulk: {
  module_type: "RadioGen"
  Nuclide: "Ar39"
  BqPercc: 1.0
  X0: -10  X1: 10  Y0: -20  Y1: 20  Z0: 0  Z1: 50
  T0: -1e6  T1: 1e6
}
"""
_FCL_DISTRIB = _FCL.replace("Z1: 50", 'Z1: 50\n  distrib_x: "TMath::Gaus(x,0,1)"')


@contextlib.contextmanager
def workdir():
    with tempfile.TemporaryDirectory() as tmp:
        def put(name, text):
            path = os.path.join(tmp, name)
            with open(path, "w") as f:
                f.write(text)
            return path
        yield tmp, put


@contextlib.contextmanager
def no_root():
    """PyROOT reported as not importable, whatever the environment."""
    saved = sys.modules.get("ROOT", ...), _root._ROOT
    sys.modules["ROOT"] = None
    _root._ROOT = None
    try:
        yield
    finally:
        if saved[0] is ...:
            sys.modules.pop("ROOT", None)
        else:
            sys.modules["ROOT"] = saved[0]
        _root._ROOT = saved[1]


def run(*args):
    res = CliRunner(mix_stderr=False).invoke(extract.cli, list(args)) \
        if _mix_stderr_supported() else CliRunner().invoke(extract.cli, list(args))
    return res


def _mix_stderr_supported():
    import inspect
    return "mix_stderr" in inspect.signature(CliRunner.__init__).parameters


def err_of(res):
    """stderr of a run (click 8.1 separates it with mix_stderr=False; 8.2+ always)."""
    try:
        return res.stderr
    except ValueError:
        return res.output


def load(path):
    with open(path) as f:
        return json.load(f)


def kinds(ir):
    out = {}
    for v in ir["volumes"]:
        out.setdefault(v["kind"], []).append(v["name"])
    return out


def test_gdml_to_ir_and_summary_table():
    with workdir() as (tmp, put):
        gdml, out = put("g.gdml", _GDML), os.path.join(tmp, "ir.json")
        res = run(gdml, "--out", out)
        assert res.exit_code == 0, res.output
        ir = load(out)
        k = kinds(ir)
        assert k["anode"] == ["volAnodePlate"]
        assert k["opdet"] == ["volArapuca_0", "volArapuca_1"]
        assert k["cathode"] == ["volCathodeGrid"]
        assert ir["detector"] == "vd" and ir.get("frame") is None
        # the summary table is on stdout, one row per placement, values in cm
        for text in ("kind", "name", "dx", "volAnodePlate", "volArapuca_1", "volCathodeGrid"):
            assert text in res.stdout, text
        row = next(r for r in res.stdout.splitlines() if "volAnodePlate" in r)
        assert row.split()[:2] == ["anode", "volAnodePlate"]
        assert row.split()[2:8] == ["-100.00", "0.00", "100.00", "2.00", "100.00", "100.00"]
        err = err_of(res)
        assert "[info] detected VD (Vertical Drift), perl GDML  (profile vd-perl)" in err
        assert "[info] 2/2 opdet placements" in err
        assert f"[info] wrote 4 volume entries to {out}" in err


def test_material_and_density_in_meta():
    """meta.material / meta.density [g/cm3] of every geometry entry that has them."""
    with workdir() as (tmp, put):
        gdml, out = put("g.gdml", _GDML), os.path.join(tmp, "ir.json")
        res = run(gdml, "--out", out, "--with-passive", "--struct-pattern", "volExtra")
        assert res.exit_code == 0, res.output
        meta = {v["name"]: v["meta"] for v in load(out)["volumes"]}
    assert meta["volAnodePlate"]["material"] == "Alu"
    assert meta["volAnodePlate"]["density"] == 2.699
    assert meta["volSteelShell"]["material"] == "Steel"
    assert abs(meta["volSteelShell"]["density"] - 7.85) < 1e-12      # 7850 kg/m3
    assert meta["volCathodeGrid"]["material"] == "G4_Cu"             # undefined: no density
    assert "density" not in meta["volCathodeGrid"]
    assert "material" not in meta["volArapuca_0"]                    # no materialref
    assert "density" not in meta["volArapuca_0"]
    assert meta["volExtra"]["material"] == "Odd" and "density" not in meta["volExtra"]
    assert "density unit 'furlong' not understood" in err_of(res)


def test_density_units():
    from duneradgeo._gdml import density_factor
    assert density_factor("") == 1.0 and density_factor("g/cm3") == 1.0
    assert density_factor("kg/m3") == 1e-3 and density_factor("mg/cm3") == 1e-3
    assert density_factor("g / mm3") == 1e3
    assert density_factor("furlong") is None and density_factor("kg") is None


def test_short_out_flag_and_help():
    with workdir() as (tmp, put):
        gdml, out = put("g.gdml", _GDML), os.path.join(tmp, "ir.json")
        assert run(gdml, "-o", out, "-v").exit_code == 0
        assert os.path.isfile(out)
    for flag in ("--help", "-h"):
        res = run(flag)
        assert res.exit_code == 0
        for text in ("--out", "--anode-pattern", "--fieldshaper-pattern", "--struct",
                     "--with-passive", "--frame", "--containment", "--npx",
                     "--gdml-backend", "--fcl-cache", "--no-fcl-cache", "--verbose",
                     "radgeo extract geometry.gdml --out ir.json",
                     "From a FCL file", "Output: a JSON IR"):
            assert text in res.output, text


def test_unknown_profile_warns():
    with workdir() as (tmp, put):
        gdml = put("g.gdml", _GDML.replace("volAnodePlate", "volAnodeX"))
        res = run(gdml, "-o", os.path.join(tmp, "ir.json"), "--anode-pattern", "^volAnodeX$")
        assert res.exit_code == 0, res.output
        assert "[warn] unknown detector type; using generic patterns" in err_of(res)


def test_frame():
    with workdir() as (tmp, put):
        gdml, out = put("g.gdml", _GDML), os.path.join(tmp, "ir.json")
        res = run(gdml, "-o", out, "--frame", "volEnclosureTPC")
        assert res.exit_code == 0, res.output
        ir = load(out)
        assert ir["frame"]["volume"] == "volEnclosureTPC"
        assert ir["frame"]["pos"] == [0.0, 0.0, 100.0]       # cm in the world
        anode = next(v for v in ir["volumes"] if v["kind"] == "anode")
        assert anode["pos"] == [-100.0, 0.0, 0.0]
        assert "[info] frame volEnclosureTPC at x=0.00 y=0.00 z=100.00" in err_of(res)
        # "world" is the default frame
        assert run(gdml, "-o", out, "--frame", "world").exit_code == 0
        assert load(out).get("frame") is None
        # a volume that is not placed exactly once
        res = run(gdml, "-o", out, "--frame", "volArapuca_.*")
        assert res.exit_code == 1
        assert "[error] --frame volArapuca_.*: expected exactly one placement" in err_of(res)


def test_with_passive():
    with workdir() as (tmp, put):
        gdml, out = put("g.gdml", _GDML), os.path.join(tmp, "ir.json")
        assert run(gdml, "-o", out).exit_code == 0
        assert "struct" not in kinds(load(out))
        res = run(gdml, "-o", out, "--with-passive")
        assert res.exit_code == 0, res.output
        assert kinds(load(out))["struct"] == ["volSteelShell"]
        assert "volSteelShell" in res.stdout


def test_variadic_pattern_options():
    with workdir() as (tmp, put):
        gdml, out = put("g.gdml", _GDML), os.path.join(tmp, "ir.json")
        res = run(gdml, "-o", out,
                  "--anode-pattern", "^volAnodePlate$", "^volNothing$",
                  "--opdet-pattern", "^volArapuca_0$",
                  "--cathode-pattern", "^volCathodeGrid$",
                  "--fieldshaper-pattern", "^volNothing$",
                  "--struct-pattern", "^volExtra$", "^volNothing2$")
        assert res.exit_code == 0, res.output
        k = kinds(load(out))
        assert k["anode"] == ["volAnodePlate"]
        assert k["opdet"] == ["volArapuca_0"]            # only one of the two
        assert k["struct"] == ["volExtra"]
        assert "fieldshaper" not in k
        # --struct: exact names, anchored, overriding --struct-pattern
        res = run(gdml, "-o", out, "--struct", "volExtra", "volSteelShell",
                  "--struct-pattern", "^volCathodeGrid$")
        assert res.exit_code == 0, res.output
        assert kinds(load(out))["struct"] == ["volExtra", "volSteelShell"]
        # the pattern options before the positional input
        res = run("-o", out, "--opdet-pattern", "^volArapuca_1$", "--", gdml)
        assert res.exit_code == 0, res.output
        assert kinds(load(out))["opdet"] == ["volArapuca_1"]


def test_nothing_extracted():
    with workdir() as (tmp, put):
        gdml, out = put("g.gdml", _GDML), os.path.join(tmp, "ir.json")
        res = run(gdml, "-o", out, "--anode-pattern", "^no$", "--opdet-pattern", "^no$",
                  "--cathode-pattern", "^no$")
        assert res.exit_code == 1
        assert "[error] nothing extracted" in err_of(res)
        assert not os.path.exists(out)


def test_missing_input_and_usage_errors():
    with workdir() as (tmp, put):
        out = os.path.join(tmp, "ir.json")
        res = run(os.path.join(tmp, "nope.gdml"), "-o", out)
        assert res.exit_code == 1
        assert "[error] " in err_of(res) and "not found on disk or in $FW_SEARCH_PATH" in err_of(res)
        assert "Traceback" not in err_of(res)
        res = run(os.path.join(tmp, "nope.fcl"), "-o", out)
        assert res.exit_code == 1 and "in $FHICL_FILE_PATH" in err_of(res)
        gdml = put("g.gdml", _GDML)
        res = run(gdml, "-o", out, "--gdml-backend", "xml")        # bad choice: usage error
        assert res.exit_code == 2
        assert run(gdml).exit_code == 2                           # --out is required
        assert run("--out", out).exit_code == 2                   # input is required


def test_fcl_radio_volumes_with_gdml():
    with workdir() as (tmp, put):
        fcl, gdml = put("job.fcl", _FCL), put("g.gdml", _GDML)
        out = os.path.join(tmp, "ir.json")
        res = run(fcl, "-o", out, "--gdml", gdml, "--no-fcl-cache")
        assert res.exit_code == 0, res.output
        ir = load(out)
        k = kinds(ir)
        assert k["radio"] == ["bulk"] and "anode" in k
        assert list(ir["producers"]) == ["bulk"]
        radio = next(v for v in ir["volumes"] if v["kind"] == "radio")
        assert radio["pos"] == [0.0, 0.0, 25.0] and radio["ext"] == [20.0, 40.0, 50.0]
        err = err_of(res)
        assert "[info] parsing FCL job.fcl" in err
        assert "[info] 1 radio volume entries" in err
        assert "and 1 producer configs" in err
        # --no-geometry: radio only, and --frame then cannot work
        res = run(fcl, "-o", out, "--no-geometry", "--no-fcl-cache")
        assert res.exit_code == 0 and set(kinds(load(out))) == {"radio"}
        res = run(fcl, "-o", out, "--no-geometry", "--no-fcl-cache", "--frame", "volX")
        assert res.exit_code == 1 and "[error] --frame needs the GDML" in err_of(res)
        # --prefix / --select filter the producers
        res = run(fcl, "-o", out, "--no-geometry", "--no-fcl-cache", "--prefix", "zzz")
        assert res.exit_code == 1 and "nothing extracted" in err_of(res)
        res = run(fcl, "-o", out, "--no-geometry", "--no-fcl-cache",
                  "--select", "^bu", "--containment", "0.9", "--npx", "200")
        assert res.exit_code == 0 and kinds(load(out))["radio"] == ["bulk"]
        res = run(fcl, "-o", out, "--no-geometry", "--no-fcl-cache", "--all-configured")
        assert res.exit_code == 0
        # --gdml that does not exist
        res = run(fcl, "-o", out, "--gdml", os.path.join(tmp, "no.gdml"), "--no-fcl-cache")
        assert res.exit_code == 1 and "[error] --gdml " in err_of(res)


def test_root_missing_is_a_clean_error():
    with workdir() as (tmp, put):
        fcl = put("job.fcl", _FCL_DISTRIB)
        out = os.path.join(tmp, "ir.json")
        with no_root():
            res = run(fcl, "-o", out, "--no-geometry", "--no-fcl-cache")
        assert res.exit_code == 1, res.output
        err = err_of(res)
        assert "[error] bulk: " in err and "PyROOT cannot be imported" in err
        assert "Traceback" not in err
        assert res.exception is None or isinstance(res.exception, SystemExit)
        assert not os.path.exists(out)


def test_root_backend_missing_is_a_clean_error():
    with workdir() as (tmp, put):
        gdml = put("g.gdml", _GDML)
        with no_root():
            res = run(gdml, "-o", os.path.join(tmp, "ir.json"), "--gdml-backend", "root")
        assert res.exit_code == 1, res.output
        assert "[error] " in err_of(res) and "ROOT" in err_of(res)
        assert isinstance(res.exception, SystemExit)


if __name__ == "__main__":
    run_tests(globals())

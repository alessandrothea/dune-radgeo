"""Checks of bki-bkg (duneradgeo.bkg): the background model, select
expressions, regions, suggestions and the commands' output.

Run with pytest, or directly: python tests/test_bkg.py.  The IR is a small
in-memory IRCollection; the logic needs neither ROOT nor rich, the rendering
checks are skipped when rich is missing.
"""
import io
import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _runner import run_tests  # noqa: E402

from duneradgeo import bkg  # noqa: E402
from duneradgeo.ir import IRCollection, VolumeEntry  # noqa: E402

try:
    import rich.console
except ImportError:
    rich = None


def raises(exc_type, fn, *args, **kw) -> str:
    try:
        fn(*args, **kw)
    except exc_type as exc:
        return str(exc)
    raise AssertionError(f"{fn.__name__}{args} did not raise {exc_type.__name__}")


def radio(name, lo, hi, **meta) -> VolumeEntry:
    lo, hi = np.array(lo, float), np.array(hi, float)
    return VolumeEntry(kind="radio", name=name, pos=(lo + hi) / 2, ext=hi - lo,
                       rot=np.eye(3), meta=meta)


def make_ir() -> IRCollection:
    """AnodeK40: 2 boxes on the plane x=325; CathodeK40: 3 boxes on x=-328;
    BulkAr39: the cryostat; MeshPo218: a layer at x=-320.9 with a distrib_x."""
    vols = []
    for i, (z0, z1) in enumerate([(0, 100), (100, 200)]):
        vols.append(radio("AnodeK40", (323, -50, z0), (327, 50, z1),
                          source="volume_gen:volAnodePlate", nuclide="K40",
                          module_type="Decay0Gen", BqPercc=0.4, surface="plane_x@325",
                          node="volAnodePlate", index=i))
    for i, (z0, z1) in enumerate([(0, 100), (100, 200), (200, 300)]):
        vols.append(radio("CathodeK40", (-330, -50, z0), (-326, 50, z1),
                          source="volume_gen:volCathodeGrid", nuclide="K40",
                          module_type="Decay0Gen", BqPercc=0.01, surface="plane_x@-328",
                          node="volCathodeGrid", index=i))
    vols.append(radio("BulkAr39", (-425, -750, -100), (425, 750, 2200),
                      source="volume_rand:volCryostat", nuclide="Ar39",
                      module_type="Decay0Gen", BqPercc=0.0014, surface="bulk", index=0))
    vols.append(radio("MeshPo218", (-323, -750, -100), (-318.8, 750, 2200),
                      source="volume_rand:volCryostat", nuclide="Po218",
                      module_type="Decay0Gen", rate=1.2, surface="plane_x@-320.9",
                      distrib_x="TMath::Gaus(x,-320,0.5)", index=0))
    vols.append(VolumeEntry(kind="anode", name="volAnodePlate", pos=np.zeros(3),
                            ext=np.ones(3), rot=np.eye(3), meta={}))
    producers = {
        "AnodeK40": {"module_type": "Decay0Gen", "isotope": "K40", "BqPercc": 0.4,
                     "material": "vm2000", "volume_gen": "volAnodePlate"},
        "CathodeK40": {"module_type": "Decay0Gen", "isotope": "K40", "BqPercc": 0.01,
                       "material": "G10", "volume_gen": "volCathodeGrid"},
        "BulkAr39": {"module_type": "Decay0Gen", "isotope": "Ar39", "BqPercc": 0.0014,
                     "material": "LAr", "volume_rand": "volCryostat"},
        "MeshPo218": {"module_type": "Decay0Gen", "isotope": "Po218", "rate": 1.2,
                      "material": "LAr", "volume_rand": "volCryostat",
                      "decay_chain": {"isotope_0": "Po218", "isotope_1": "Pb214"},
                      "distrib_x": "TMath::Gaus(x,-320,0.5)"},
    }
    return IRCollection(source="x.fcl", detector="vd", volumes=vols, producers=producers)


def names(bgs) -> list:
    return [b.name for b in bgs]


def sel(*exprs) -> list:
    bgs = bkg.backgrounds(make_ir())
    es = [bkg.parse_expr(e) for e in exprs]
    bkg.check_fields(es, bgs)
    return names(bkg.select(bgs, es))


# --- model -----------------------------------------------------------------

def test_backgrounds():
    bgs = bkg.backgrounds(make_ir())
    assert names(bgs) == ["AnodeK40", "BulkAr39", "CathodeK40", "MeshPo218"]
    assert [len(b.volumes) for b in bgs] == [2, 1, 3, 1]       # geometry is not a background


def test_fields_merge_parameters_and_volume_metadata():
    b = {b.name: b for b in bkg.backgrounds(make_ir())}["MeshPo218"]
    f = b.fields
    assert f["name"] == ["MeshPo218"] and f["nvolumes"] == [1]
    assert f["isotope"] == ["Po218"] and f["nuclide"] == ["Po218"]
    assert f["decay_chain"] == ["Po218", "Pb214"]              # a dict contributes its values
    assert f["surface"] == ["plane_x@-320.9"] and f["rate"] == [1.2]
    assert f["module_type"] == ["Decay0Gen"]                   # in both, once


def test_summary():
    b = {b.name: b for b in bkg.backgrounds(make_ir())}
    s = bkg.summary(b["CathodeK40"])
    assert (s["nuclide"], s["activity"], s["volumes"], s["material"]) == ("K40", "0.01", 3, "G10")
    assert s["surfaces"] == [{"surface": "plane_x@-328", "volumes": 3}]
    m = bkg.summary(b["MeshPo218"])
    assert m["distrib"] == ["x"] and m["activity_field"] == "rate"


# --- select ----------------------------------------------------------------

def test_parse_expr():
    e = bkg.parse_expr("nuclide=Po218")
    assert (e.field, e.op, e.value) == ("nuclide", "=", "Po218") and str(e) == "nuclide=Po218"
    assert bkg.parse_expr(" BqPercc >= 0.5 ").number == 0.5
    assert bkg.parse_expr("distrib_x").op is None
    assert "bad expression" in raises(ValueError, bkg.parse_expr, "=x")
    assert "needs a number" in raises(ValueError, bkg.parse_expr, "BqPercc<abc")
    assert "bad regex" in raises(ValueError, bkg.parse_expr, "name=(")


def test_select_regex():
    assert sel("nuclide=K40") == ["AnodeK40", "CathodeK40"]
    assert sel("name=^Cathode") == ["CathodeK40"]
    assert sel("material=vm2000|G10") == ["AnodeK40", "CathodeK40"]
    assert sel("decay_chain=Pb214") == ["MeshPo218"]             # a value inside a dict
    assert sel("surface=plane_x@-3") == ["CathodeK40", "MeshPo218"]


def test_select_negation_numeric_exists_and_and():
    assert sel("nuclide!=K40") == ["BulkAr39", "MeshPo218"]
    assert sel("node!=Anode") == ["BulkAr39", "CathodeK40", "MeshPo218"]   # absent counts
    assert sel("BqPercc>0.1") == ["AnodeK40"]
    assert sel("BqPercc<=0.01") == ["BulkAr39", "CathodeK40"]
    assert sel("rate>1") == ["MeshPo218"]                        # absent field: no match
    assert sel("nvolumes>=3") == ["CathodeK40"]
    assert sel("distrib_x") == ["MeshPo218"]
    assert sel("nuclide=K40", "BqPercc<0.1") == ["CathodeK40"]
    assert sel("nuclide=nothing") == []


def test_unknown_field_suggests():
    bgs = bkg.backgrounds(make_ir())
    msg = raises(ValueError, bkg.check_fields, [bkg.parse_expr("nucide=K40")], bgs)
    assert "unknown field 'nucide'" in msg and "nuclide" in msg
    assert "known fields" in raises(ValueError, bkg.check_fields,
                                    [bkg.parse_expr("qqqqqq=1")], bgs)


# --- region ----------------------------------------------------------------

def hits(**kw) -> dict:
    r = bkg.make_region(**kw)
    return {b.name: idx for b, idx in bkg.in_region(bkg.backgrounds(make_ir()), r)}


def test_region_slice_touching_counts():
    h = hits(slices={"x": -328})
    assert sorted(h) == ["BulkAr39", "CathodeK40"]             # mesh layer is at -323..-318.8
    assert h["CathodeK40"] == [0, 1, 2]
    assert sorted(hits(slices={"x": -323})) == ["BulkAr39", "MeshPo218"]   # touching counts
    assert sorted(hits(slices={"x": -320}, thickness=1)) == ["BulkAr39", "MeshPo218"]


def test_region_ranges_and_box():
    assert sorted(hits(ranges={"x": (300, 400)})) == ["AnodeK40", "BulkAr39"]
    h = hits(ranges={"x": (300, 400), "z": (150, 160)})
    assert h["AnodeK40"] == [1]                                  # only the second box reaches z=150
    h = hits(box=(-400, -300, -10, 10, 250, 260))
    assert sorted(h) == ["BulkAr39", "CathodeK40", "MeshPo218"] and h["CathodeK40"] == [2]
    assert hits(box=(500, 600, 0, 1, 0, 1)) == {}


def test_region_errors_and_text():
    assert "empty range" in raises(ValueError, bkg.make_region, ranges={"x": (5, 1)})
    assert "do not intersect" in raises(ValueError, bkg.make_region,
                                        ranges={"x": (0, 1)}, slices={"x": 5})
    assert bkg.make_region(slices={"x": -328}).text() == "x = -328"
    assert bkg.make_region(ranges={"z": (0, 5)}).text() == "z in [0, 5]"
    assert bkg.make_region().text() == "everywhere"
    assert "bad --slice" in raises(ValueError, bkg._parse_slices, ["x"])
    assert bkg._parse_slices(["X=1.5", "y = -2"]) == {"x": 1.5, "y": -2.0}


def test_region_rows():
    bgs = bkg.backgrounds(make_ir())
    rows = {r["name"]: r for r in bkg.region_rows(bkg.in_region(
        bgs, bkg.make_region(ranges={"x": (300, 400), "z": (150, 160)})))}
    r = rows["AnodeK40"]
    assert (r["overlapping"], r["volumes"]) == (1, 2)
    assert r["lo"] == [323.0, -50.0, 100.0] and r["hi"] == [327.0, 50.0, 200.0]


# --- suggestions -----------------------------------------------------------

def test_suggestions():
    s = bkg.suggestions(bkg.backgrounds(make_ir()))
    fields = {f["field"]: f for f in s["fields"]}
    assert "name" not in fields and "index" not in fields
    k40 = [v for v in fields["nuclide"]["values"] if v["value"] == "K40"][0]
    assert k40["backgrounds"] == 2 and k40["command"] == "select nuclide=K40"
    assert all(v["backgrounds"] >= 2 for f in s["fields"] for v in f["values"])
    num = {n["field"]: n for n in s["numeric"]}
    assert set(num) == {"BqPercc"}                 # rate has one value: nothing to split
    assert num["BqPercc"]["min"] == 0.0014 and num["BqPercc"]["max"] == 0.4
    planes = {r["slice"]: r for r in s["regions"]}
    assert planes["x=-328"]["on_plane"] == 1 and planes["x=-328"]["volumes_on_plane"] == 3
    assert planes["x=-328"]["crossing"] == 2
    assert s["extent"]["x"] == [-425.0, 425.0]
    # every suggested command is a valid selection of the stated size
    bgs = bkg.backgrounds(make_ir())
    for f in s["fields"]:
        for v in f["values"]:
            expr = v["command"].split(" ", 1)[1].strip("'")
            assert len(bkg.select(bgs, [bkg.parse_expr(expr)])) == v["backgrounds"]


def test_anchored_value_regex():
    assert bkg._plain_or_anchored("Pb21", ["Pb21", "Pb210"]) == "^Pb21$"
    assert bkg._plain_or_anchored("Bi214", ["Bi214", "K40"]) == "Bi214"
    assert bkg._plain_or_anchored("a+b", ["a+b"]) == r"a\+b"


# --- command line ----------------------------------------------------------

def run(*argv) -> tuple:
    """(exit code or None, stdout) of the command line; stderr text on exit."""
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "ir.json")
        make_ir().save(path)
        out, old = io.StringIO(), sys.stdout
        sys.stdout = out
        try:
            code = bkg.run(bkg._build_parser().parse_args([path, *argv]))
        except SystemExit as exc:
            code = exc.code
        finally:
            sys.stdout = old
    return code, out.getvalue()


def test_json_commands():
    code, out = run("--json", "list")
    rows = json.loads(out)
    assert code == 0 and [r["name"] for r in rows] == ["AnodeK40", "BulkAr39", "CathodeK40", "MeshPo218"]
    code, out = run("--json", "select", "nuclide=K40", "BqPercc<0.1")
    d = json.loads(out)
    assert d["expressions"] == ["nuclide=K40", "BqPercc<0.1"]
    assert [r["name"] for r in d["backgrounds"]] == ["CathodeK40"]
    code, out = run("--json", "show", "Mesh")
    d = json.loads(out)
    assert d[0]["name"] == "MeshPo218" and d[0]["params"]["rate"] == 1.2
    assert d[0]["volumes"][0]["lo"][0] == -323.0
    code, out = run("--json", "region", "--slice", "x=325", "--volumes")
    d = json.loads(out)
    assert d["frame"] == "world" and [r["name"] for r in d["backgrounds"]] == ["AnodeK40", "BulkAr39"]
    assert len(d["volumes"]["AnodeK40"]) == 2
    assert d["region"]["lo"] == [325.0, None, None] and d["region"]["hi"] == [325.0, None, None]
    code, out = run("--json", "suggest")
    assert "fields" in json.loads(out)


def test_command_errors_exit_cleanly():
    for argv, text in ((("select", "nucide=K40"), "unknown field"),
                       (("select", "name=("), "bad regex"),
                       (("show", "NoSuchBackground"), "no background matches"),
                       (("region", "--slice", "q=1"), "bad --slice"),
                       (("region", "--x", "5", "1"), "empty range")):
        code, _ = run("--json", *argv)
        assert isinstance(code, str) and text in code and code.startswith("[error]"), (argv, code)


# --- rendering (rich) ------------------------------------------------------

def render(fn, *a, **kw) -> str:
    buf = io.StringIO()
    console = rich.console.Console(file=buf, width=200, force_terminal=False, color_system=None)
    fn(console, *a, **kw)
    return buf.getvalue()


def test_render_tables():
    if rich is None:
        return
    bgs = bkg.backgrounds(make_ir())
    text = render(bkg.render_list, [bkg.summary(b) for b in bgs])
    assert "CathodeK40" in text and "4 background(s)" in text and "plane_x@-328 ×3" in text
    d = bkg.detail(bgs[3])
    text = render(bkg.render_show, d)
    assert "decay_chain" in text and "distrib_x" in text and "plane_x@-320.9" in text
    r = bkg.make_region(slices={"x": 325})
    text = render(bkg.render_region, r, bkg.region_rows(bkg.in_region(bgs, r)), 4, "world")
    assert "x = 325" in text and "AnodeK40" in text and "2 of 4" in text
    text = render(bkg.render_suggestions, bkg.suggestions(bgs), 4)
    assert "select nuclide=K40" in text and "region --slice x=-328" in text


if __name__ == "__main__":
    run_tests(globals())

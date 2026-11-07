"""ROOT-free checks of duneradgeo._fcl: the FHiCL parser, the #include
expansion and producer_volumes.

Run with pytest, or directly: python tests/test_fcl.py.  Needs no PyROOT:
producer_volumes is only given configurations without distrib_x/y/z, and the
geometry is a small GDML of plain boxes read by the pure-Python backend.
"""
import contextlib
import functools
import io
import math
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _runner import run_tests  # noqa: E402

from duneradgeo._fcl import FhiclParser, expand_fcl, producer_volumes  # noqa: E402
from duneradgeo._gdml import GDMLGeometry  # noqa: E402


def parse(text: str) -> dict:
    return FhiclParser(text).parse()


def raises(exc_type, fn, *args, **kw) -> str:
    """Call fn, check that it raises exc_type, return the message."""
    try:
        fn(*args, **kw)
    except exc_type as exc:
        return str(exc)
    raise AssertionError(f"{fn.__name__}{args} did not raise {exc_type.__name__}")


@contextlib.contextmanager
def environ(**values):
    """Set environment variables for the duration of the block."""
    old = {k: os.environ.get(k) for k in values}
    os.environ.update(values)
    try:
        yield
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


@contextlib.contextmanager
def quiet():
    """Capture stderr (the [info]/[warn] lines) into the yielded buffer."""
    buf = io.StringIO()
    with contextlib.redirect_stderr(buf):
        yield buf


# ---------------------------------------------------------------------------
# FhiclParser
# ---------------------------------------------------------------------------

def test_scalars():
    d = parse("""
        i: 3  neg: -12  f: 2.5  e: -1.5e-3  E: 2E+4  dot: .5  trail: 1.  plus: +4
        t: true  F: FALSE  nil: @nil
        s: "a \\"quoted\\" word\\n"  raw: 'no \\n escapes'  empty: ""
        word: bareword
    """)
    assert d["i"] == 3 and d["neg"] == -12 and d["plus"] == 4
    assert d["f"] == 2.5 and d["e"] == -1.5e-3 and d["E"] == 2e4
    assert d["dot"] == 0.5 and d["trail"] == 1.0
    assert d["t"] is True and d["F"] is False
    assert d["nil"] is None
    assert d["s"] == 'a "quoted" word\n'
    assert d["raw"] == "no \\n escapes"      # single quotes: no escape processing
    assert d["empty"] == ""
    assert d["word"] == "bareword"


def test_infinity():
    d = parse("a: infinity  b: Infinity")
    assert d["a"] == math.inf and d["b"] == math.inf


def test_signed_infinity():
    # fhicl-cpp writes infinite values as +infinity / -infinity
    d = parse("a: +infinity  b: -infinity  c: [-infinity, +infinity]")
    assert d["a"] == math.inf and d["b"] == -math.inf
    assert d["c"] == [-math.inf, math.inf]


def test_non_ascii_string():
    assert parse('s: "été → \\"x\\""')["s"] == 'été → "x"'


def test_nested_tables_and_sequences():
    d = parse("""
        a: {
          b: { c: 1  d: "x" }
          seq: [1, [2, 3], [], { e: 4 }]
          tables: [ { name: "p" rate: 1 }, { name: "q" rate: 2 } ]
          empty: {}
        }
    """)
    assert d["a"]["b"] == {"c": 1, "d": "x"}
    assert d["a"]["seq"] == [1, [2, 3], [], {"e": 4}]
    assert [t["name"] for t in d["a"]["tables"]] == ["p", "q"]
    assert d["a"]["tables"][1]["rate"] == 2
    assert d["a"]["empty"] == {}


def test_dotted_keys_and_overrides():
    p = FhiclParser("""
        a: { b: { c: 1 } d: [1, [2, 3]] }
        z: 0
        a.b.c: 5
        a.d[1][0]: 9
        a.new.deep: "n"
        z: 7
    """)
    d = p.parse()
    assert d["a"]["b"]["c"] == 5
    assert d["a"]["d"] == [1, [9, 3]]
    assert d["a"]["new"] == {"deep": "n"}
    assert d["z"] == 7
    assert p.order == ["a", "z"]     # overrides do not add top-level keys


def test_erase():
    d = parse("a: { x: 1 y: 2 }  a.x: @erase  b: 3  b: @erase")
    assert d == {"a": {"y": 2}}


def test_local_references():
    d = parse("""
        proto: { x: 1 y: { z: 2 } s: [10, 20] }
        copy: @local::proto
        copy.y.z: 3
        leaf: @local::proto.y.z
        elem: @local::proto.s[1]
    """)
    assert d["copy"] == {"x": 1, "y": {"z": 3}, "s": [10, 20]}
    assert d["proto"]["y"]["z"] == 2     # the reference is a copy
    assert d["leaf"] == 2
    assert d["elem"] == 20


def test_table_and_sequence_splices():
    d = parse("""
        p: { x: 1 y: 2 }
        q: { @table::p y: 3 w: 4 }
        s: [1, 2]
        t: [0, @sequence::s, 3]
        @table::p
    """)
    assert d["q"] == {"x": 1, "y": 3, "w": 4}
    assert d["t"] == [0, 1, 2, 3]
    assert d["x"] == 1 and d["y"] == 2    # top-level @table:: splices into the root


def test_comments():
    d = parse("""
        # a comment
        // another
        a: 1   # trailing
        b: "x # not a comment // either"   // trailing
        c: [1, # inside a sequence
            2]
    """)
    assert d == {"a": 1, "b": "x # not a comment // either", "c": [1, 2]}


def test_prolog():
    p = FhiclParser("""
        BEGIN_PROLOG
        proto: { module_type: "RadioGen" }
        END_PROLOG
        prod: @local::proto
    """)
    d = p.parse()
    assert d["prod"] == {"module_type": "RadioGen"}
    assert p.order == ["proto", "prod"]


def test_protect_annotations():
    d = parse("a @protect_ignore: 1  b @protect_error: { c: 2 }")
    assert d == {"a": 1, "b": {"c": 2}}


def test_parse_errors():
    msg = raises(SyntaxError, parse, "a 5")
    assert "expected ':' after 'a'" in msg and "'5'" in msg, msg
    msg = raises(SyntaxError, parse, ": 1")
    assert "expected key" in msg, msg
    msg = raises(SyntaxError, parse, "a: }")
    assert "unexpected token '}'" in msg, msg
    msg = raises(SyntaxError, parse, "a: 1\nb: 2\nc: $x")
    assert msg.startswith("line 3:") and "cannot tokenize" in msg and "$x" in msg, msg
    # unterminated table / sequence
    raises(SyntaxError, parse, "a: { b: 1")
    raises(SyntaxError, parse, "a: [1, 2")


def test_unknown_reference():
    try:
        parse("a: @local::nosuchkey")
    except (LookupError, SyntaxError) as exc:
        assert "nosuchkey" in str(exc), exc
    else:
        raise AssertionError("an undefined @local:: must raise")


# ---------------------------------------------------------------------------
# #include expansion (inline resolver; fhicl-dump is kept off $PATH)
# ---------------------------------------------------------------------------

def _write(path: str, text: str) -> str:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def _expand(path: str, search: str, **kw) -> str:
    with environ(PATH="", FHICL_FILE_PATH=search), quiet():
        return expand_fcl(path, **kw)


def test_include_search_path():
    with tempfile.TemporaryDirectory() as tmp:
        inc = os.path.join(tmp, "inc")
        _write(os.path.join(inc, "common.fcl"),
               'BEGIN_PROLOG\nproto: { module_type: "Decay0Gen" rate: 1 }\nEND_PROLOG\n')
        _write(os.path.join(inc, "sub", "more.fcl"),
               '#include "common.fcl"\nother: @local::proto\n')
        top = _write(os.path.join(tmp, "job", "top.fcl"),
                     '#include "sub/more.fcl"\n'
                     'physics: { producers: { p: @local::proto } }\n'
                     'physics.producers.p.rate: 2\n')
        text = _expand(top, os.pathsep.join([os.path.join(tmp, "missing"), inc]))
        assert "#include" not in text
        d = parse(text)
        assert d["physics"]["producers"]["p"] == {"module_type": "Decay0Gen", "rate": 2}
        assert d["other"]["rate"] == 1


def test_include_missing():
    with tempfile.TemporaryDirectory() as tmp:
        top = _write(os.path.join(tmp, "top.fcl"), '#include "nowhere.fcl"\na: 1\n')
        with environ(PATH="", FHICL_FILE_PATH=tmp), quiet():
            msg = raises(FileNotFoundError, expand_fcl, top)
        assert '"nowhere.fcl"' in msg and "$FHICL_FILE_PATH" in msg, msg


def test_include_search_var_and_absolute():
    with tempfile.TemporaryDirectory() as tmp:
        _write(os.path.join(tmp, "a", "x.fcl"), "x: 1\n")
        absy = _write(os.path.join(tmp, "b", "y.fcl"), "y: 2\n")
        top = _write(os.path.join(tmp, "top.fcl"),
                     f'#include "x.fcl"\n#include "{absy}"\nz: 3\n')
        with environ(PATH="", FHICL_FILE_PATH="", MY_FCL_PATH=os.path.join(tmp, "a")), \
                quiet():
            text = expand_fcl(top, search_var="MY_FCL_PATH")
        assert parse(text) == {"x": 1, "y": 2, "z": 3}


def test_include_circular():
    with tempfile.TemporaryDirectory() as tmp:
        _write(os.path.join(tmp, "a.fcl"), '#include "b.fcl"\na: 1\n')
        _write(os.path.join(tmp, "b.fcl"), '#include "a.fcl"\nb: 2\n')
        top = _write(os.path.join(tmp, "top.fcl"), '#include "a.fcl"\ntop: 0\n')
        d = parse(_expand(top, tmp))
        assert d == {"a": 1, "b": 2, "top": 0}


def test_no_include_verbatim():
    with tempfile.TemporaryDirectory() as tmp:
        body = "# only a comment\na: 1  // and #include in a comment is no directive\n"
        top = _write(os.path.join(tmp, "top.fcl"), body)
        assert _expand(top, tmp) == body


def test_include_cache():
    with tempfile.TemporaryDirectory() as tmp:
        cache = os.path.join(tmp, "cache")
        inc = _write(os.path.join(tmp, "inc.fcl"), "v: 1\n")
        top = _write(os.path.join(tmp, "top.fcl"), '#include "inc.fcl"\n')
        assert parse(_expand(top, tmp, cache_dir=cache)) == {"v": 1}
        assert len(os.listdir(cache)) == 1
        # an edited included file is not detected (documented): cached text
        _write(inc, "v: 2\n")
        assert parse(_expand(top, tmp, cache_dir=cache)) == {"v": 1}
        # editing the top-level file invalidates the entry
        _write(top, '#include "inc.fcl"\nw: 3\n')
        assert parse(_expand(top, tmp, cache_dir=cache)) == {"v": 2, "w": 3}
        # so does a different search path
        other = os.path.join(tmp, "other")
        _write(os.path.join(other, "inc.fcl"), "v: 4\n")
        assert parse(_expand(top, other, cache_dir=cache)) == {"v": 4, "w": 3}


# ---------------------------------------------------------------------------
# producer_volumes (no distrib_*: no ROOT)
# ---------------------------------------------------------------------------

# Lengths in mm in the GDML, cm in the rows.  volCryo (200 x 100 x 300 cm) is
# placed at z = +50 cm in the world; inside it, two volBox (10 x 20 x 30 cm),
# the second rotated by 90 deg about z, a thin volThin (0.2 x 100 x 100 cm)
# and a volBoxInner (same solid as volBox, a name volBox is a prefix of).
_GDML = """<?xml version="1.0"?>
<gdml>
 <define><rotation name="rz90" z="90" unit="deg"/></define>
 <solids>
  <box name="World" x="10000" y="10000" z="10000" lunit="mm"/>
  <box name="Cryo" x="2000" y="1000" z="3000" lunit="mm"/>
  <box name="Box" x="100" y="200" z="300" lunit="mm"/>
  <box name="Thin" x="2" y="1000" z="1000" lunit="mm"/>
 </solids>
 <structure>
  <volume name="volBox"><solidref ref="Box"/></volume>
  <volume name="volBoxInner"><solidref ref="Box"/></volume>
  <volume name="volThin"><solidref ref="Thin"/></volume>
  <volume name="volCryo"><solidref ref="Cryo"/>
   <physvol><volumeref ref="volBox"/>
    <position name="p1" x="100" y="0" z="0" unit="mm"/></physvol>
   <physvol><volumeref ref="volBox"/>
    <position name="p2" x="-300" y="0" z="0" unit="mm"/><rotationref ref="rz90"/></physvol>
   <physvol><volumeref ref="volThin"/>
    <position name="p3" x="500" y="0" z="0" unit="mm"/></physvol>
   <physvol><volumeref ref="volBoxInner"/>
    <position name="p4" x="0" y="0" z="1000" unit="mm"/></physvol>
  </volume>
  <volume name="volWorld"><solidref ref="World"/>
   <physvol><volumeref ref="volCryo"/>
    <position name="pc" x="0" y="0" z="500" unit="mm"/></physvol>
  </volume>
 </structure>
 <setup name="Default" version="1.0"><world ref="volWorld"/></setup>
</gdml>
"""
CRYO   = ([-100, -50, -100], [100, 50, 200])
BOX_1  = ([5, -10, 35], [15, 10, 65])
BOX_2  = ([-40, -5, 35], [-20, 5, 65])      # rotated: x and y extents swapped
THIN   = ([49.9, -50, 0], [50.1, 50, 100])
INNER  = ([-5, -10, 135], [5, 10, 165])


@functools.lru_cache(maxsize=None)
def geometry() -> GDMLGeometry:
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(os.path.join(tmp, "boxes.gdml"), _GDML)
        return GDMLGeometry(path)


def same_box(row: dict, box: tuple) -> bool:
    return (np.allclose(row["lo"], box[0], atol=1e-9)
            and np.allclose(row["hi"], box[1], atol=1e-9))


def volumes(name, cfg, geo=None):
    with quiet() as err:
        rows = producer_volumes(name, cfg, geo)
    return rows, err.getvalue()


def test_explicit_box():
    cfg = dict(module_type="RadioGen", Nuclide="K42", BqPercc=0.5, rate=2.0,
               X0=-10, X1=10, Y0=-20, Y1=20, Z0=0, Z1=50, T0=-1e6, T1=1e6)
    rows, _ = volumes("ar42", cfg)
    assert len(rows) == 1
    r = rows[0]
    assert (r["producer"], r["module_type"], r["index"], r["source"]) == \
        ("ar42", "RadioGen", 0, "explicit")
    assert r["lo"] == [-10, -20, 0] and r["hi"] == [10, 20, 50]
    assert (r["nuclide"], r["BqPercc"], r["rate"]) == ("K42", 0.5, 2.0)
    assert (r["T0"], r["T1"]) == (-1e6, 1e6)
    assert r["material"] == ".*"
    assert r["surface"] == "bulk"
    assert not any(k.startswith("distrib") for k in r)


def test_explicit_box_arrays():
    cfg = dict(module_type="RadioGen", material="LAr",
               Nuclide=["Ar39", "Kr85"], BqPercc=[1.0, 2.0],
               X0=[0, 10], X1=[1, 20], Y0=-5, Y1=5, Z0=0, Z1=[100, 200],
               T0=[0, 0], T1=[1, 2])
    rows, _ = volumes("multi", cfg)
    assert [r["index"] for r in rows] == [0, 1]
    assert [r["nuclide"] for r in rows] == ["Ar39", "Kr85"]
    assert [r["BqPercc"] for r in rows] == [1.0, 2.0]
    assert [r["T1"] for r in rows] == [1, 2]
    assert rows[0]["lo"] == [0, -5, 0] and rows[0]["hi"] == [1, 5, 100]
    assert rows[1]["lo"] == [10, -5, 0] and rows[1]["hi"] == [20, 5, 200]
    assert all(r["material"] == "LAr" for r in rows)


def test_explicit_box_errors():
    box = dict(X0=0, X1=1, Y0=0, Y1=1, Z0=0, Z1=1)
    cfg = dict(box, module_type="RadioGen")
    del cfg["Z1"]
    msg = raises(ValueError, producer_volumes, "p", cfg)
    assert "p:" in msg and "Z1" in msg, msg
    cfg = dict(box, module_type="RadioGen", X0=[0, 1, 2], X1=[1, 2])
    msg = raises(ValueError, producer_volumes, "p", cfg)
    assert "length mismatch" in msg, msg


def test_explicit_plane_surface():
    """A box thinner than 5 % of its smallest other side is a plane at its
    mid-coordinate; a thicker one is bulk."""
    cases = [
        (dict(X0=-320.1, X1=-319.9, Y0=-100, Y1=100, Z0=0, Z1=300), "plane_x@-320"),
        (dict(X0=-100, X1=100, Y0=599.5, Y1=600.5, Z0=0, Z1=300), "plane_y@600"),
        (dict(X0=-100, X1=100, Y0=-100, Y1=100, Z0=0, Z1=0.5), "plane_z@0.25"),
        (dict(X0=0, X1=10, Y0=-100, Y1=100, Z0=0, Z1=300), "bulk"),
    ]
    for box, surface in cases:
        rows, _ = volumes("p", dict(box, module_type="RadioGen", Nuclide="Rn222"))
        assert rows[0]["surface"] == surface, (box, rows[0]["surface"])


def test_decay0gen_nuclide():
    geo = geometry()
    rows, _ = volumes("p", dict(module_type="Decay0Gen", isotope="Rn222",
                                volume_rand="volCryo"), geo)
    assert rows[0]["nuclide"] == "Rn222"
    chain = {"isotope_2": "Pb214", "isotope_10": "Pb210", "isotope_1": "Bi214"}
    rows, _ = volumes("p", dict(module_type="Decay0Gen", decay_chain=chain,
                                volume_rand="volCryo"), geo)
    assert rows[0]["nuclide"] == "Bi214+Pb214+Pb210"    # numeric suffix order
    rows, _ = volumes("p", dict(module_type="Decay0Gen", isotope="Ar39",
                                X0=0, X1=1, Y0=0, Y1=1, Z0=0, Z1=1))
    assert rows[0]["nuclide"] == "Ar39"


def test_decay0gen_explicit_decay_chain():
    cfg = dict(module_type="Decay0Gen", decay_chain={"isotope_1": "Bi214", "isotope_2": "Pb214"},
               X0=0, X1=1, Y0=0, Y1=1, Z0=0, Z1=1)
    rows, _ = volumes("p", cfg)
    assert rows[0]["nuclide"] == "Bi214+Pb214"


def test_volume_rand():
    geo = geometry()
    cfg = dict(module_type="Decay0Gen", isotope="Ar39", volume_rand="volCryo",
               material="LAr", rate=3.0)
    rows, err = volumes("ar39", cfg, geo)
    assert len(rows) == 1 and err == ""
    r = rows[0]
    assert (r["index"], r["source"], r["material"], r["rate"]) == \
        (0, "volume_rand:volCryo", "LAr", 3.0)
    assert same_box(r, CRYO), (r["lo"], r["hi"])
    assert r["surface"] == "bulk"
    # the name is matched exactly: a regex character is literal
    rows, err = volumes("p", dict(cfg, volume_rand="volCry."), geo)
    assert rows[0]["lo"] is None and "not found" in err
    # no geometry: no box, no surface
    rows, err = volumes("p", cfg)
    assert rows[0]["lo"] is None and rows[0]["hi"] is None
    assert "surface" not in rows[0] and err == ""
    # unknown volume: a warning and no box
    rows, err = volumes("p", dict(cfg, volume_rand="volNone"), geo)
    assert rows[0]["lo"] is None and "surface" not in rows[0]
    assert "[warn]" in err and "volNone" in err, err


def test_volume_rand_first_placement():
    rows, _ = volumes("p", dict(module_type="Decay0Gen", volume_rand="volBox"), geometry())
    assert same_box(rows[0], BOX_1), (rows[0]["lo"], rows[0]["hi"])


def test_volume_gen():
    geo = geometry()
    cfg = dict(module_type="Decay0Gen", isotope="K40", volume_gen="^(volBox|volThin)$")
    rows, err = volumes("k40", cfg, geo)
    assert err == ""
    assert [r["index"] for r in rows] == [0, 1, 2]
    assert [r["node"] for r in rows] == ["volBox", "volBox", "volThin"]
    assert all(r["source"] == "volume_gen:^(volBox|volThin)$" for r in rows)
    assert all(r["nuclide"] == "K40" for r in rows)
    for r, box in zip(rows, (BOX_1, BOX_2, THIN)):
        assert same_box(r, box), (r["node"], r["lo"], r["hi"])
    assert [r["surface"] for r in rows] == ["bulk", "bulk", "plane_x@50"]
    # default volume_gen: every placement, the world included
    rows, _ = volumes("p", dict(module_type="Decay0Gen"), geo)
    assert all(r["source"] == "volume_gen:.*" for r in rows)
    assert sorted(r["node"] for r in rows) == sorted(
        ["volWorld", "volCryo", "volBox", "volBox", "volThin", "volBoxInner"])
    # no geometry: a single row without a box
    rows, err = volumes("p", cfg)
    assert len(rows) == 1 and rows[0]["lo"] is None and err == ""
    assert "surface" not in rows[0] and "node" not in rows[0]
    # no match: a warning and a single row without a box
    rows, err = volumes("p", dict(cfg, volume_gen="^nomatch$"), geo)
    assert len(rows) == 1 and rows[0]["lo"] is None
    assert "[warn]" in err and "nomatch" in err, err


def test_volume_gen_full_match():
    rows, _ = volumes("p", dict(module_type="Decay0Gen", volume_gen="volBox"), geometry())
    assert [r["node"] for r in rows] == ["volBox", "volBox"]


if __name__ == "__main__":
    run_tests(globals())

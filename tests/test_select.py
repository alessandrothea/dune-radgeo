"""ROOT-free checks of duneradgeo._select: the --pick / --surface /
--highlight spec parsers, validate_picks, select_volumes and the highlight
grouping.

Run with pytest, or directly: python tests/test_select.py.  The IR is a small
in-memory IRCollection built from VolumeEntry objects; no file is read.
"""
import contextlib
import io
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _runner import run_tests  # noqa: E402

from duneradgeo._select import (  # noqa: E402
    highlight_groups, match_highlight, parse_highlight_spec, parse_pick_spec,
    select_volumes, validate_picks)
from duneradgeo.ir import IRCollection, VolumeEntry  # noqa: E402


def raises(exc_type, fn, *args, **kw) -> str:
    """Call fn, check that it raises exc_type, return the message."""
    try:
        fn(*args, **kw)
    except exc_type as exc:
        return str(exc)
    raise AssertionError(f"{fn.__name__}{args} did not raise {exc_type.__name__}")


def vol(kind: str, name: str, **meta) -> VolumeEntry:
    return VolumeEntry(kind=kind, name=name, pos=np.zeros(3), ext=np.ones(3),
                       rot=np.eye(3), meta=meta)


_CRYO = "volDetEnclosure:0/volCryostat:0"


def make_ir() -> IRCollection:
    """Two anode placements, a cathode, a two-part field shaper, the
    cryostat, a geometry volume from an old IR (no path), and two radio
    producers (one box per index)."""
    vols = [
        vol("anode", "volAnodePlate", path=f"{_CRYO}/volAnodePlate:0",
            pv="volAnodePlate_pv", pvpath="volCryostat_pv/volAnodePlate_pv"),
        vol("anode", "volAnodePlate", path=f"{_CRYO}/volAnodePlate:1",
            pv="volAnodePlate_pv", pvpath="volCryostat_pv/volAnodePlate_pv"),
        vol("cathode", "volCathodeGrid", path=f"{_CRYO}/volCathodeGrid:0"),
        vol("fieldshaper", "volFieldShaper", part="FieldShaperLongtube", copy=0,
            path=f"{_CRYO}/volFieldShaper:0"),
        vol("fieldshaper", "volFieldShaper", part="FieldShaperShorttube", copy=0,
            path=f"{_CRYO}/volFieldShaper:0"),
        vol("struct", "volCryostat", path=_CRYO),
        vol("opdet", "volOldArapuca", pv="volOldArapuca_pv"),     # no path
        vol("radio", "Ar39GenInLAr", source="volume_rand:volCryostat",
            nuclide="Ar39", surface="bulk", node=None, index=0),
        vol("radio", "K42Gen", nuclide="K42", surface="cathode", index=0),
        vol("radio", "K42Gen", nuclide="K42", surface="anode", index=1),
        vol("radio", "K42Gen", nuclide="K42", surface="bulk", node="volTPC",
            index=2),
        vol("radio", "K42Gen", nuclide="K42", surface="anode", index=3),
    ]
    return IRCollection(source="test", detector="vd", volumes=vols)


def ids(vols) -> list:
    """(name, index-or-path-or-part) per volume, for readable comparisons."""
    out = []
    for v in vols:
        m = v.meta
        key = m.get("index") if v.kind == "radio" else m.get("part", m.get("path"))
        out.append((v.name, key))
    return out


def all_ids() -> list:
    return ids(make_ir().volumes)


# ---------------------------------------------------------------------------
# parse_pick_spec
# ---------------------------------------------------------------------------

def test_pick_regex_selector():
    name_re, indices, match = parse_pick_spec("^K42Gen#@cathode|anode")
    assert name_re.pattern == "^K42Gen"
    assert indices is None
    assert match.pattern == "cathode|anode"


def test_pick_empty_regex_selector():
    # '#@' alone: an empty @REGEX, which matches any present metadata field
    _, indices, match = parse_pick_spec("x#@")
    assert indices is None and match.pattern == ""


def test_pick_index_lists():
    for sel, want in (("0", {0}), ("0,2", {0, 2}), ("0,2,5-7", {0, 2, 5, 6, 7}),
                      ("4-4", {4}), (" 1 , 3 ", {1, 3}), ("1,,2,", {1, 2}),
                      ("2,1-3", {1, 2, 3})):
        name_re, indices, match = parse_pick_spec(f"^Gen#{sel}")
        assert name_re.pattern == "^Gen", sel
        assert indices == want, (sel, indices)
        assert match is None, sel


def test_pick_index_list_degenerate():
    # observed: a reversed range and a list of separators only are accepted
    # and select nothing (no error)
    assert parse_pick_spec("x#3-1")[1] == set()
    assert parse_pick_spec("x#,")[1] == set()


def test_pick_hash_split_is_last_hash():
    # rpartition: NAME_RE may contain '#', an @REGEX may not
    name_re, indices, _ = parse_pick_spec("a#b#1")
    assert name_re.pattern == "a#b" and indices == {1}
    msg = raises(ValueError, parse_pick_spec, "x#@a#b")
    assert msg == "bad index list 'b' in spec 'x#@a#b'", msg


def test_pick_bad_specs():
    for spec in ("nohash", "#3", "x#", "#", ""):
        msg = raises(ValueError, parse_pick_spec, spec)
        assert msg == f"bad spec '{spec}' (expected NAME_RE#SEL)", msg


def test_pick_bad_index_lists():
    for spec, sel in (("^volTPCPlane#1-x", "1-x"), ("x#1-", "1-"), ("x#-1", "-1"),
                      ("x#a", "a"), ("x#1-2-3", "1-2-3"), ("x#1.5", "1.5")):
        msg = raises(ValueError, parse_pick_spec, spec)
        assert msg == f"bad index list '{sel}' in spec '{spec}'", msg


def test_pick_bad_regex_is_re_error():
    # a bad NAME_RE or @REGEX raises re.error (not a ValueError);
    # validate_picks turns it into "--OPT: ..."
    for spec in ("(#@x", "x#@(", "(#1"):
        try:
            parse_pick_spec(spec)
        except re.error:
            continue
        except ValueError as exc:
            raise AssertionError(f"{spec}: ValueError {exc}")
        raise AssertionError(f"{spec}: no error")


# ---------------------------------------------------------------------------
# parse_highlight_spec
# ---------------------------------------------------------------------------

def test_highlight_plain_regex():
    label, name_re, indices, match = parse_highlight_spec("^volAnode")
    assert label == "^volAnode" and name_re.pattern == "^volAnode"
    assert indices is None and match is None


def test_highlight_label():
    for spec, want_label, want_re in (("a=^Gen", "a", "^Gen"),
                                      ("my group=^Gen", "my group", "^Gen"),
                                      ("  a  =^Gen", "a", "^Gen"),
                                      ("a=b=c", "a", "b=c")):
        label, name_re, indices, match = parse_highlight_spec(spec)
        assert (label, name_re.pattern) == (want_label, want_re), spec
        assert indices is None and match is None, spec


def test_highlight_no_label_with_regex_chars():
    # a regex character before '=' means there is no LABEL: the whole spec is
    # the NAME_RE and the label
    for spec in ("a.b=x", "^a=b", "a|b=c", "a(1)=c", "a\\d=c"):
        label, name_re, _, _ = parse_highlight_spec(spec)
        assert label == spec and name_re.pattern == spec, spec


def test_highlight_empty_after_label():
    # observed: 'x=' is not split (LABEL= needs a non-empty rest), so the
    # NAME_RE is the literal 'x=' and it matches nothing (no error)
    label, name_re, indices, match = parse_highlight_spec("x=")
    assert (label, name_re.pattern, indices, match) == ("x=", "x=", None, None)


def test_highlight_with_selector():
    label, name_re, indices, match = parse_highlight_spec("a=^Gen#0-1,3")
    assert (label, name_re.pattern, indices, match) == ("a", "^Gen", {0, 1, 3}, None)

    label, name_re, indices, match = parse_highlight_spec("^Gen#@bulk")
    assert label == "^Gen#@bulk" and name_re.pattern == "^Gen"
    assert indices is None and match.pattern == "bulk"

    # '=' after the '#' is part of the @REGEX
    label, name_re, _, match = parse_highlight_spec("lbl=x#@k=v")
    assert (label, name_re.pattern, match.pattern) == ("lbl", "x", "k=v")


def test_highlight_bad_specs():
    msg = raises(ValueError, parse_highlight_spec, "")
    assert msg == "bad highlight spec '' (expected [LABEL=]NAME_RE[#SEL])", msg
    # with a '#', the NAME_RE#SEL part is reported as a pick spec
    for spec, want in (("a=#3", "bad spec '#3' (expected NAME_RE#SEL)"),
                       ("a=^Gen#", "bad spec '^Gen#' (expected NAME_RE#SEL)"),
                       ("#3", "bad spec '#3' (expected NAME_RE#SEL)"),
                       ("a=x#1-", "bad index list '1-' in spec 'x#1-'")):
        msg = raises(ValueError, parse_highlight_spec, spec)
        assert msg == want, (spec, msg)
    for spec in ("(", "a=(", "a=x#@("):
        try:
            parse_highlight_spec(spec)
        except re.error:
            continue
        raise AssertionError(f"{spec}: no re.error")


def test_highlight_label_spaces_around_equals():
    label, name_re, _, _ = parse_highlight_spec("grp = ^K42")
    assert label == "grp"
    assert name_re.search("K42Gen"), name_re.pattern


# ---------------------------------------------------------------------------
# validate_picks
# ---------------------------------------------------------------------------

def test_validate_ok():
    validate_picks()
    validate_picks(surface=None, pick=[], highlight=None)
    validate_picks(surface=["^K42#@anode"], pick=["^Gen#0,2", "^vol#@:1$"],
                   highlight=["a=^K42#1", "^volCathode", "x="])


def test_validate_messages():
    cases = (
        (dict(pick=["nohash"]), "--pick: bad spec 'nohash' (expected NAME_RE#SEL)"),
        (dict(pick=["^volTPCPlane#1-x"]),
         "--pick: bad index list '1-x' in spec '^volTPCPlane#1-x'"),
        (dict(surface=["#3"]), "--surface: bad spec '#3' (expected NAME_RE#SEL)"),
        (dict(highlight=[""]),
         "--highlight: bad highlight spec '' (expected [LABEL=]NAME_RE[#SEL])"),
        (dict(highlight=["a=x#1-"]),
         "--highlight: bad index list '1-' in spec 'x#1-'"),
        # a bad regex: re.error becomes ValueError("--OPT: <re message>")
        (dict(pick=["(#@x"]), "--pick: missing ), unterminated subpattern at position 0"),
        (dict(surface=["x#@a["]),
         "--surface: unterminated character set at position 1"),
        (dict(highlight=["("]),
         "--highlight: missing ), unterminated subpattern at position 0"),
    )
    for kw, want in cases:
        msg = raises(ValueError, validate_picks, **kw)
        assert msg == want, (kw, msg)


def test_validate_order():
    # surface before pick before highlight, whatever the keyword order
    msg = raises(ValueError, validate_picks, highlight=[""], pick=["p"], surface=["s"])
    assert msg.startswith("--surface: bad spec 's'"), msg
    msg = raises(ValueError, validate_picks, highlight=[""], pick=["p"])
    assert msg.startswith("--pick: bad spec 'p'"), msg
    # the first bad spec of an option
    msg = raises(ValueError, validate_picks, pick=["a#1", "b", "c"])
    assert msg == "--pick: bad spec 'b' (expected NAME_RE#SEL)", msg


def test_validate_does_not_check_volumes():
    # spec syntax only: an index list for a geometry name is accepted here
    # (select_volumes rejects it once it meets a geometry volume)
    validate_picks(pick=["^volAnodePlate#0"], highlight=["^volAnodePlate#0-3"])


# ---------------------------------------------------------------------------
# select_volumes
# ---------------------------------------------------------------------------

def test_select_no_filter():
    assert ids(select_volumes(make_ir())) == all_ids()


def test_select_kinds():
    ir = make_ir()
    got = select_volumes(ir, kinds=["anode", "struct"])
    assert {v.kind for v in got} == {"anode", "struct"} and len(got) == 3
    assert select_volumes(ir, kinds=[]) == []          # [] is not None
    assert ids(select_volumes(ir, kinds=None)) == all_ids()


def test_select_name_patterns():
    ir = make_ir()
    # re.search (unanchored), patterns OR'ed, order kept
    got = ids(select_volumes(ir, name_patterns=["Cath", "^K42"]))
    assert got == [("volCathodeGrid", f"{_CRYO}/volCathodeGrid:0"),
                   ("K42Gen", 0), ("K42Gen", 1), ("K42Gen", 2), ("K42Gen", 3)], got
    assert select_volumes(ir, name_patterns=["^nothing$"]) == []
    assert ids(select_volumes(ir, name_patterns=[])) == all_ids()


def test_select_kind_patterns():
    ir = make_ir()
    # only the kinds named are filtered; the others are kept
    got = select_volumes(ir, kind_patterns={"radio": ["^K42"], "anode": ["nomatch"]})
    names = [v.name for v in got]
    assert "Ar39GenInLAr" not in names and names.count("K42Gen") == 4
    assert "volAnodePlate" not in names
    assert {"volCathodeGrid", "volFieldShaper", "volCryostat",
            "volOldArapuca"} <= set(names)


def test_select_surfaces():
    ir = make_ir()
    got = ids(select_volumes(ir, surfaces=["^K42#@anode"]))
    want = [i for i in all_ids() if i[0] != "K42Gen"] + [("K42Gen", 1), ("K42Gen", 3)]
    assert sorted(got, key=str) == sorted(want, key=str), got
    # index list; Ar39 is not named by the spec, so it is kept
    got = ids(select_volumes(ir, kinds=["radio"], surfaces=["^K42#0,2"]))
    assert got == [("Ar39GenInLAr", 0), ("K42Gen", 0), ("K42Gen", 2)], got


def test_select_surfaces_ignore_geometry():
    # --surface only filters radio volumes: an index list naming geometry
    # volumes (or volumes without a path) is not an error
    ir = make_ir()
    got = ids(select_volumes(ir, surfaces=["^vol#0", "Arapuca#@x"]))
    assert got == all_ids(), got


def test_select_picks_geometry_by_path():
    ir = make_ir()
    got = ids(select_volumes(ir, picks=["^volAnodePlate#@:1$"]))
    want = [i for i in all_ids() if i != ("volAnodePlate", f"{_CRYO}/volAnodePlate:0")]
    assert got == want, got


def test_select_picks_metadata_keys():
    ir = make_ir()
    # meta.part
    got = ids(select_volumes(ir, kinds=["fieldshaper"], picks=["Shaper#@Short"]))
    assert got == [("volFieldShaper", "FieldShaperShorttube")], got
    # meta.pv / meta.pvpath (both anodes have them)
    for pick in ("^volAnodePlate#@^volAnodePlate_pv$", "^volAnodePlate#@Cryostat_pv/"):
        got = select_volumes(ir, kinds=["anode"], picks=[pick])
        assert len(got) == 2, pick
    # radio meta.nuclide / meta.node; a None field is skipped (no 'None' match)
    got = ids(select_volumes(ir, kinds=["radio"], picks=["Gen#@^volTPC$"]))
    assert got == [("K42Gen", 2)], got
    ar39 = dict(kinds=["radio"], name_patterns=["Ar39"])
    got = ids(select_volumes(ir, picks=["Ar39#@None"], **ar39))
    assert got == [], got
    got = ids(select_volumes(ir, picks=["Ar39#@^Ar39$"], **ar39))
    assert got == [("Ar39GenInLAr", 0)], got
    # keys outside _PICK_KEYS (source, index) are not searched
    got = ids(select_volumes(ir, picks=["Ar39#@volume_rand"], **ar39))
    assert got == [], got


def test_select_picks_radio():
    ir = make_ir()
    got = ids(select_volumes(ir, kinds=["radio"], picks=["^K42Gen#1-2"]))
    assert got == [("Ar39GenInLAr", 0), ("K42Gen", 1), ("K42Gen", 2)], got
    # specs with the same NAME_RE are OR'ed
    got = ids(select_volumes(ir, kinds=["radio"], picks=["^K42Gen#0", "^K42Gen#@^bulk$"]))
    assert got == [("Ar39GenInLAr", 0), ("K42Gen", 0), ("K42Gen", 2)], got
    # one spec for each producer
    got = ids(select_volumes(ir, kinds=["radio"], picks=["^K42Gen#3", "^Ar39#5"]))
    assert got == [("K42Gen", 3)], got


def test_select_picks_unnamed_volumes_kept():
    # a pick whose NAME_RE matches no volume filters nothing (no error)
    ir = make_ir()
    assert ids(select_volumes(ir, picks=["^volTPCPlane#0"])) == all_ids()
    assert ids(select_volumes(ir, picks=["^nothing#@x"])) == all_ids()


def test_select_index_on_geometry_raises():
    msg = raises(ValueError, select_volumes, make_ir(), picks=["^volAnodePlate#0"])
    assert msg == ("index selector matches anode volume 'volAnodePlate'; "
                   "geometry volumes are picked by path (NAME_RE#@REGEX)"), msg
    # the kinds filter runs first: radgeo profiles (kinds=["radio"]) does not see it
    got = select_volumes(make_ir(), kinds=["radio"], picks=["^vol#0", "^K42#0"])
    assert ids(got) == [("Ar39GenInLAr", 0), ("K42Gen", 0)]


def test_select_no_path_raises():
    msg = raises(ValueError, select_volumes, make_ir(), picks=["Arapuca#@_pv"])
    assert msg == ("opdet volume 'volOldArapuca' has no placement path; "
                   "re-run radgeo extract to pick geometry volumes"), msg
    # radio volumes need no path
    assert select_volumes(make_ir(), kinds=["radio"], picks=["Gen#@K42"])


def test_select_bad_pick_spec():
    # select_volumes reports the bare parser message (no "--OPT:" prefix)
    msg = raises(ValueError, select_volumes, make_ir(), picks=["nohash"])
    assert msg == "bad spec 'nohash' (expected NAME_RE#SEL)", msg
    msg = raises(ValueError, select_volumes, make_ir(), surfaces=["x#1-"])
    assert msg == "bad index list '1-' in spec 'x#1-'", msg


def test_select_filters_compose():
    ir = make_ir()
    got = ids(select_volumes(ir, kinds=["radio", "anode"], name_patterns=["Gen|Anode"],
                             kind_patterns={"radio": ["K42"]},
                             surfaces=["K42#@anode"], picks=["K42#3", "Anode#@:0$"]))
    assert got == [("volAnodePlate", f"{_CRYO}/volAnodePlate:0"), ("K42Gen", 3)], got


def test_select_index_on_geometry_raises_any_order():
    for picks in (["^volCathodeGrid#0", "^volCathodeGrid#@Cathode"],   # raises
                  ["^volCathodeGrid#@Cathode", "^volCathodeGrid#0"]):  # does not
        try:
            select_volumes(make_ir(), picks=picks)
        except ValueError as exc:
            assert "geometry volumes are picked by path" in str(exc), exc
        else:
            raise AssertionError(f"picks={picks}: no ValueError")


def test_select_bad_name_regex_is_value_error():
    for kw in (dict(name_patterns=["("]), dict(kind_patterns={"radio": ["("]})):
        try:
            select_volumes(make_ir(), **kw)
        except ValueError:
            continue
        except re.error as exc:
            raise AssertionError(f"{kw}: re.error {exc}")
        raise AssertionError(f"{kw}: no error")


# ---------------------------------------------------------------------------
# match_highlight / highlight_groups
# ---------------------------------------------------------------------------

def specs(*hl) -> list:
    return [parse_highlight_spec(h) for h in hl]


def test_match_highlight():
    ir = make_ir()
    anode0, anode1, cathode = ir.volumes[:3]
    k42_0, k42_1 = ir.volumes[8:10]
    sp = specs("a=^K42#1", "^volAnode#@:1$", "Gen", "volCathode")
    assert match_highlight(k42_1, sp) == 0
    assert match_highlight(k42_0, sp) == 2          # first matching spec wins
    assert match_highlight(anode1, sp) == 1
    assert match_highlight(anode0, sp) is None      # named, but SEL misses it
    assert match_highlight(cathode, sp) == 3
    assert match_highlight(ir.volumes[5], sp) is None   # not named: no match
    assert match_highlight(cathode, []) is None


def test_match_highlight_errors_and_no_path():
    ir = make_ir()
    msg = raises(ValueError, match_highlight, ir.volumes[0], specs("a=^volAnodePlate#0-3"))
    assert msg == ("index selector matches anode volume 'volAnodePlate'; "
                   "geometry volumes are highlighted by path (NAME_RE#@REGEX)"), msg
    # no placement path: not an error, the other fields are searched
    old = ir.volumes[6]
    assert match_highlight(old, specs("Arapuca#@x")) is None
    assert match_highlight(old, specs("Arapuca#@_pv$")) == 0


def test_highlight_groups():
    ir = make_ir()
    vols = ir.volumes
    labels, groups, context = highlight_groups(
        vols, ["cat=^K42#@cathode", "^volAnodePlate", "Gen#@bulk", "rest=Gen"])
    assert labels == ["cat", "^volAnodePlate", "Gen#@bulk", "rest"]
    assert ids(groups[0]) == [("K42Gen", 0)]
    assert [v.meta["path"] for v in groups[1]] == [f"{_CRYO}/volAnodePlate:0",
                                                  f"{_CRYO}/volAnodePlate:1"]
    assert ids(groups[2]) == [("Ar39GenInLAr", 0), ("K42Gen", 2)]
    assert ids(groups[3]) == [("K42Gen", 1), ("K42Gen", 3)]
    assert [v.name for v in context] == ["volCathodeGrid", "volFieldShaper",
                                         "volFieldShaper", "volCryostat", "volOldArapuca"]
    # a partition of the input, by identity
    out = [id(v) for g in groups for v in g] + [id(v) for v in context]
    assert sorted(out) == sorted(id(v) for v in vols)


def test_highlight_groups_no_match():
    # a spec matching nothing gives an empty group and no output: the
    # "[warn] --highlight 'LABEL' matches no selected volume" line is
    # printed by plot_ir, not here
    vols = make_ir().volumes
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        labels, groups, context = highlight_groups(
            vols, ["none=^NoSuchVolume", "a=^Gen#@nomatch", "x=", "^volOldArapuca#@x"])
    assert labels == ["none", "a", "x=", "^volOldArapuca#@x"]
    assert groups == [[], [], [], []]
    assert ids(context) == all_ids()
    assert out.getvalue() == "" and err.getvalue() == ""


def test_highlight_groups_empty():
    vols = make_ir().volumes
    labels, groups, context = highlight_groups(vols, [])
    assert labels == [] and groups == [] and ids(context) == all_ids()
    assert highlight_groups([], ["a=x"]) == (["a"], [[]], [])


def test_highlight_groups_errors():
    vols = make_ir().volumes
    msg = raises(ValueError, highlight_groups, vols, ["a=#3"])
    assert msg == "bad spec '#3' (expected NAME_RE#SEL)", msg
    msg = raises(ValueError, highlight_groups, vols, ["^volCathodeGrid#0"])
    assert msg == ("index selector matches cathode volume 'volCathodeGrid'; "
                   "geometry volumes are highlighted by path (NAME_RE#@REGEX)"), msg


if __name__ == "__main__":
    run_tests(globals())

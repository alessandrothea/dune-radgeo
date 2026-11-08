"""Checks of the ``radgeo`` command group (duneradgeo.cli): the subcommands,
the ``bkg`` alias.

Run with pytest, or directly: python tests/test_radgeo_cli.py.  Needs click
and rich (skipped without them).
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
if importlib.util.find_spec("click") is None or importlib.util.find_spec("rich") is None:
    try:
        import pytest
        pytest.skip("click and rich are needed", allow_module_level=True)
    except ImportError:
        print("skipped: click and rich are needed")
        sys.exit(0)

from click.testing import CliRunner  # noqa: E402
from _runner import run_tests  # noqa: E402

from duneradgeo import cli as radgeo  # noqa: E402

SUBCOMMANDS = ("extract", "plot", "measure", "profiles", "backgrounds")


def test_group_lists_subcommands():
    r = CliRunner().invoke(radgeo.cli, ["--help"], prog_name="radgeo")
    assert r.exit_code == 0
    for name in SUBCOMMANDS:
        assert name in r.output


def test_subcommand_help_names_radgeo():
    for name in SUBCOMMANDS:
        r = CliRunner().invoke(radgeo.cli, [name, "--help"], prog_name="radgeo")
        assert r.exit_code == 0, r.output
        assert f"Usage: radgeo {name}" in r.output


def test_bkg_alias():
    r = CliRunner().invoke(radgeo.cli, ["bkg", "--help"], prog_name="radgeo")
    assert r.exit_code == 0
    assert "Usage: radgeo backgrounds" in r.output


if __name__ == "__main__":
    run_tests(globals())

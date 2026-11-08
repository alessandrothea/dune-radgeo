"""Checks of duneradgeo._cli, the click/rich helpers shared by the radgeo
commands: VariadicOption, CliError and the info/warn/error lines.

Run with pytest, or directly: python tests/test_cli.py.  Needs click and rich
(skipped without them).
"""
import importlib.util
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
if importlib.util.find_spec("click") is None or importlib.util.find_spec("rich") is None:
    try:
        import pytest
        pytest.skip("click and rich are needed", allow_module_level=True)
    except ImportError:
        print("skipped: click and rich are needed")
        sys.exit(0)

import click  # noqa: E402
from click.testing import CliRunner  # noqa: E402
from _runner import run_tests  # noqa: E402

from duneradgeo import _cli  # noqa: E402


@click.command(context_settings=_cli.CONTEXT_SETTINGS)
@click.argument("src")
@click.option("--name", "-n", cls=_cli.VariadicOption, metavar="PATTERN...")
@click.option("--kind", cls=_cli.VariadicOption, metavar="KIND...")
@click.option("--x", type=float)
@click.option("--flag", is_flag=True)
@click.option("--num", type=int)
def demo(src, name, kind, x, flag, num):
    click.echo(f"src={src} name={list(name)} kind={list(kind)} x={x} flag={flag} num={num}")


def run(*args):
    res = CliRunner().invoke(demo, list(args))
    return res.exit_code, res.output.strip()


def test_variadic_collects_until_next_option():
    assert run("a.json", "--name", "A", "B", "C", "--kind", "radio") == \
        (0, "src=a.json name=['A', 'B', 'C'] kind=['radio'] x=None flag=False num=None")


def test_variadic_stops_at_flags_and_short_options():
    code, out = run("a.json", "-n", "A", "B", "--flag")
    assert code == 0 and "name=['A', 'B']" in out and "flag=True" in out
    code, out = run("--flag", "a.json", "--name", "A")
    assert code == 0 and "src=a.json" in out and "name=['A']" in out


def test_variadic_repeats_accumulate_and_equals_form():
    code, out = run("a.json", "--name", "A", "--kind", "k", "--name", "B", "C")
    assert code == 0 and "name=['A', 'B', 'C']" in out
    code, out = run("a.json", "--name=-weird", "--flag")
    assert code == 0 and "name=['-weird']" in out and "src=a.json" in out


def test_variadic_numbers_are_values():
    code, out = run("a.json", "--name", "-5", "7", "--x", "-3")
    assert code == 0 and "name=['-5', '7']" in out and "x=-3.0" in out


def test_variadic_absent_is_empty_and_missing_value_is_an_error():
    code, out = run("a.json")
    assert code == 0 and "name=[]" in out
    code, out = run("a.json", "--name")
    assert code == 2 and "requires an argument" in out


def test_variadic_leaves_the_positional_to_the_user():
    # like argparse nargs="+": a positional after the values is swallowed,
    # so the command line puts it first (or after another option)
    code, out = run("--name", "A", "B", "a.json")
    assert code == 2 and "Missing argument" in out
    code, out = run("--name", "A", "B", "--flag", "a.json")
    assert code == 0 and "src=a.json" in out


def test_help_flags():
    code, out = run("-h")
    assert code == 0 and "--name PATTERN..." in out


@click.command()
def fail():
    raise _cli.CliError("something is wrong")


def test_cli_error_format():
    res = CliRunner().invoke(fail, [])
    assert res.exit_code == 1 and res.output.strip() == "[error] something is wrong"


@click.command()
def talk():
    _cli.info("parsing [x] …")
    _cli.warn("careful")
    _cli.error("bad")
    click.echo("stdout line")


def test_message_lines_keep_their_text():
    res = CliRunner().invoke(talk, [])
    assert res.exit_code == 0
    lines = res.output.splitlines()
    # brackets in the message are not rich markup
    assert lines == ["[info] parsing [x] …", "[warn] careful", "[error] bad", "stdout line"]


def test_console_is_wide_when_piped():
    os.environ.pop("COLUMNS", None)
    assert _cli.console().width == _cli.PIPE_WIDTH
    os.environ["COLUMNS"] = "100"
    try:
        assert _cli.console().width == 100
    finally:
        os.environ.pop("COLUMNS")


def test_messages_without_rich():
    import builtins
    real = builtins.__import__

    def no_rich(name, *a, **kw):
        if name == "rich" or name.startswith("rich."):
            raise ImportError(name)
        return real(name, *a, **kw)
    builtins.__import__ = no_rich
    try:
        res = CliRunner().invoke(talk, [])
    finally:
        builtins.__import__ = real
    assert res.exit_code == 0 and res.output.splitlines()[:3] == [
        "[info] parsing [x] …", "[warn] careful", "[error] bad"]


if __name__ == "__main__":
    run_tests(globals())

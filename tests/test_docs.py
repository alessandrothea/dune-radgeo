"""Checks that the documentation matches the commands.

* docs/bki-<command>.md contains the current ``--help`` text of the command
  (click's layout, width 80; ``bki-bkg``: the group and each subcommand), so
  the option blocks cannot drift silently;
* every relative markdown link ``[text](path#anchor)`` in README.md and
  docs/*.md points to an existing file and, with an anchor, to a heading of
  that file (GitHub anchor rules).

Run with pytest, or directly: python tests/test_docs.py.  The help check
needs click and rich (skipped without them).

To refresh a block after changing a command, print its help and paste it into
the fenced block of the docs file:

    env XDG_CACHE_HOME=/home/user/.cache python tests/test_docs.py --print bki-extract

(``--print`` with no name prints all five; XDG_CACHE_HOME is fixed so that the
``--fcl-cache`` default does not depend on the machine.)
"""
import importlib.util
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _runner import run_tests  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COMMANDS = ("bki-extract", "bki-plot", "bki-measure", "bki-distrib", "bki-bkg")
CACHE_HOME = "/home/user/.cache"       # fixes the --fcl-cache default shown

_HELP_PROGRAM = r"""
import importlib, sys
from click.testing import CliRunner
cmd = sys.argv[1]
mod = importlib.import_module("duneradgeo." + cmd.split("-", 1)[1])
runner = CliRunner()
def show(args):
    r = runner.invoke(mod.cli, args, terminal_width=80, prog_name=cmd)
    sys.stdout.write(r.output + "\n====\n")
show(["--help"])
if cmd == "bki-bkg":                      # a group: every subcommand too
    for name in mod.cli.commands:
        show(["IR_FILE", name, "--help"])
"""


def have_click_rich() -> bool:
    return all(importlib.util.find_spec(m) is not None for m in ("click", "rich"))


def help_texts(cmd: str) -> list:
    """The ``--help`` outputs of `cmd` at width 80 (several for bki-bkg)."""
    env = dict(os.environ, XDG_CACHE_HOME=CACHE_HOME)
    env["PYTHONPATH"] = os.pathsep.join(
        [os.path.join(ROOT, "src")] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    out = subprocess.run([sys.executable, "-c", _HELP_PROGRAM, cmd], env=env,
                         capture_output=True, text=True, check=True).stdout
    return [t.strip("\n") for t in out.split("\n====\n") if t.strip()]


def _norm(text: str) -> str:
    return "\n".join(line.rstrip() for line in text.splitlines())


def doc_path(cmd: str) -> str:
    return os.path.join(ROOT, "docs", cmd + ".md")


def test_help_in_docs():
    if not have_click_rich():
        try:
            import pytest
            pytest.skip("click and rich are needed")
        except ImportError:
            print("skipped: click and rich are needed")
            return
    for cmd in COMMANDS:
        with open(doc_path(cmd), encoding="utf-8") as f:
            doc = _norm(f.read())
        for text in help_texts(cmd):
            text = _norm(text)
            assert text in doc, (
                f"docs/{cmd}.md does not contain the current `{cmd} --help` text:\n"
                f"{text}\n\nRegenerate it with\n"
                f"  env XDG_CACHE_HOME={CACHE_HOME} python tests/test_docs.py --print {cmd}\n"
                f"and paste it in a fenced block of docs/{cmd}.md")


# ---- links --------------------------------------------------------------

LINK = re.compile(r"(?<!!)\[[^\]]*\]\(([^)\s]+)\)")


def markdown_files() -> list:
    docs = os.path.join(ROOT, "docs")
    files = [os.path.join(ROOT, "README.md")]
    files += sorted(os.path.join(docs, n) for n in os.listdir(docs) if n.endswith(".md"))
    return files


def github_anchor(heading: str) -> str:
    """GitHub's slug: lowercase, spaces to '-', punctuation dropped."""
    heading = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", heading)   # links: keep the text
    heading = heading.replace("`", "").replace("*", "").strip().lower()
    return re.sub(r"[^\w\- ]", "", heading).replace(" ", "-")


def anchors(path: str) -> set:
    """The heading anchors of a markdown file (fenced blocks skipped,
    repeated headings get -1, -2, ... as on GitHub)."""
    seen, out, fenced = {}, set(), False
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.lstrip().startswith("```"):
                fenced = not fenced
                continue
            m = None if fenced else re.match(r"#{1,6}\s+(.*?)\s*#*\s*$", line)
            if m:
                a = github_anchor(m.group(1))
                n = seen.get(a, 0)
                seen[a] = n + 1
                out.add(a if n == 0 else f"{a}-{n}")
    return out


def broken_links() -> list:
    bad = []
    for path in markdown_files():
        fenced = False
        with open(path, encoding="utf-8") as f:
            lines = f.read().splitlines()
        for no, line in enumerate(lines, 1):
            if line.lstrip().startswith("```"):
                fenced = not fenced
            if fenced:
                continue
            for target in LINK.findall(re.sub(r"`[^`]*`", "", line)):
                if re.match(r"[a-z][a-z0-9+.-]*:", target):      # http:, mailto:
                    continue
                file_part, _, anchor = target.partition("#")
                dest = path if not file_part else os.path.normpath(
                    os.path.join(os.path.dirname(path), file_part))
                where = f"{os.path.relpath(path, ROOT)}:{no}: ({target})"
                if not os.path.exists(dest):
                    bad.append(where + " no such file")
                elif anchor and dest.endswith(".md") and anchor not in anchors(dest):
                    bad.append(where + f" no heading #{anchor} in {os.path.relpath(dest, ROOT)}")
    return bad


def test_links_resolve():
    bad = broken_links()
    assert not bad, "broken links:\n" + "\n".join(bad)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--print":
        for name in (sys.argv[2:] or COMMANDS):
            print("\n\n".join(help_texts(name)))
            print()
    else:
        run_tests(globals())

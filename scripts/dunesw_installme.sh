#!/bin/bash
# Install duneradgeo into the active venv of a dunesw/LArSoft (UPS) environment.
#
# Usage (from the host, with dsw):
#   dsw -C <area> exec bash /path/to/dune-radgeo/scripts/dunesw_installme.sh [--no-editable] [pip args...]
# or from inside a set-up environment (dsw shell, setup_dunesw.sh) with the venv active:
#   bash scripts/dunesw_installme.sh
#
# Default: editable install (pip install -e) with --no-deps; numpy comes from the stack.
# Extra arguments are passed on to pip install.
#
# Why PYTHONPATH is cleared for pip: the UPS product python_future_six_request puts
# pip 21.2.4 on PYTHONPATH, which shadows the venv's own pip (even with python -m pip).
# pip < 21.3 can't do editable installs from pyproject.toml and fails with
# "editable mode currently requires a setuptools-based build". Only pip runs without
# PYTHONPATH; the installed tools still use the stack's numpy at run time.
set -euo pipefail

PKG=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)

if [ -z "${VIRTUAL_ENV:-}" ]; then
  echo "dunesw_installme.sh: no venv active (VIRTUAL_ENV unset)." >&2
  echo "  Create one in the area (python -m venv <area>/.venv; dsw then activates it)" >&2
  echo "  or activate it before running this script." >&2
  exit 1
fi

editable=(-e)
if [ "${1:-}" = "--no-editable" ]; then
  editable=()
  shift
fi

echo "venv:    $VIRTUAL_ENV"
echo "package: $PKG"
env -u PYTHONPATH python -m pip --version
env -u PYTHONPATH python -m pip install --no-deps "$@" "${editable[@]}" "$PKG"

# The commands need click and rich at run time; the install above is --no-deps,
# so add whichever of the two the venv's python cannot import yet.
missing=()
for mod in click rich; do
  python -c "import $mod" 2>/dev/null || missing+=("$mod")
done
if [ "${#missing[@]}" -gt 0 ]; then
  echo "installing missing run-time packages: ${missing[*]}"
  if ! env -u PYTHONPATH python -m pip install "${missing[@]}"; then
    echo "dunesw_installme.sh: could not install ${missing[*]} (no network?)." >&2
    echo "  duneradgeo is installed but its commands need them; install them by hand with" >&2
    echo "  env -u PYTHONPATH python -m pip install ${missing[*]}" >&2
    exit 1
  fi
fi
python -c "import click, rich" || { echo "dunesw_installme.sh: click/rich still not importable" >&2; exit 1; }

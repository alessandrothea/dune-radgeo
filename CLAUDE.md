# dune-radgeo (package `duneradgeo`, command `radgeo`)

Inspect DUNE background sources without a LArSoft environment: detector
geometry from GDML, radiological emission volumes from RadioGen / Decay0Gen
FHiCL, and the `distrib_x/y/z` position profiles of those producers.  Everything
flows through a JSON "IR" (`ir.py`: `IRCollection`, `VolumeEntry`).  The README
and `docs/radgeo-<command>.md` are the user documentation; this file is the
map for working on the code.

## Commands (`src/duneradgeo/cli.py`, click + rich)

| Command | Module | Does |
|---|---|---|
| `radgeo extract` | `extract.py` | GDML / FHiCL → JSON IR |
| `radgeo plot` | `plot.py` | render an IR as PDF/PNG (`notebook_plot.py` for Plotly) |
| `radgeo measure` | `measure.py`, `measurements/` | geometric measurements on an IR |
| `radgeo profiles` | `distrib.py`, `_tf1.py` | draw `distrib_*` TF1 profiles (needs PyROOT) |
| `radgeo backgrounds` (alias `bkg`) | `bkg.py` | list/select backgrounds of an IR |

Scripts: `scripts/compare_backends.py` (compare two GDML backends),
`scripts/compare_ir.py` (compare two IRs), `scripts/run_notebook.py`.

## Layout

- `_gdml.py`: `GeometryBase` (placement search, box decomposition, bounding
  boxes) and the pure-Python `GDMLGeometry` parser; `load_geometry(path,
  backend=...)`.  Tables a backend fills (`_vols`, `_prims`, `_bools`,
  `_boxes`, `_shape_bbox`, `_world`) are documented on `GeometryBase`; lengths in cm.
- GDML backends (`BACKENDS`): `python` (default, numpy only, ~0.1 s) and
  `root` (`_gdml_root.py`, TGeoManager::Import, needs ROOT with TGeo/GDML).
  A pyg4ometry backend (evaluates constants and expressions the built-in
  parser cannot) was tried and archived on the local branch
  `archive/pyg4ometry` (commit e01ae0b): too much weight for now.
- `_fcl.py`: FHiCL parser and `producer_volumes` (radio boxes, `volume_rand` /
  `volume_gen`).  `_profiles.py`: detector profiles (name patterns per kind).
  `_search.py`: FCL → GDML path search.  `_root.py`: shared lazy PyROOT import
  (`ROOTUnavailable`).  `_select.py`, `_style.py`, `_draw.py`, `_cli.py`: shared helpers.
- `docs/`: one page per command; the `--help` text is pasted into it and checked.

## Conventions worth knowing

- Physvols are placed with the inverse of the GDML rotation (`rot.T`); boolean
  operands use it directly (see `_gdml.py`).  All backends follow this.
- ROOT is imported lazily; the package must work without it.  Features that
  need it raise `ROOTUnavailable` (a `RuntimeError`), which the commands print
  as `[error] …`, exit 1.
- Dependencies stay light: `radgeo extract` needs only numpy, click, rich.
  `root` is a uv dependency *group* (Python ≥ 3.14 host venv).

## Tests

```bash
uv run pytest -q tests          # or: python tests/test_x.py (each file self-runs)
```

- `tests/test_tf1.py` needs PyROOT and fails without it.
- `tests/test_docs.py` fails when a command's `--help` changes without
  updating its docs page; regenerate with
  `env XDG_CACHE_HOME=/home/user/.cache python tests/test_docs.py --print "radgeo extract"`
  and paste the block into `docs/radgeo-extract.md`.  It also checks markdown links.
- GDML fixtures in the other tests are inline strings.

## Environment

- Host: AlmaLinux 9, Python 3.9 system; the project venv (`.venv`, `uv`) is the
  place to run things.  `lar`, `root`, `fhicl-dump` etc. must run in a dunesw
  area through the `dsw` tool (see the `dsw` skill), never on the bare host.
- `requires-python >= 3.9`: keep `from __future__ import annotations` in new modules.

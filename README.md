# duneradgeo

Inspect DUNE background sources without a LArSoft environment:

- **detector geometry** — anode (CRP/APA), optical detector, cathode,
  field-shaper, and structural volume positions from GDML files
- **radiological volumes** — emission boxes from RadioGen / Decay0Gen FHiCL
  configurations, with each producer's full FHiCL table
- **radiological position profiles** — the `distrib_x/y/z` formulas of those
  producers, drawn over the range the generator samples
- **background inspection** — list the backgrounds of an IR, select them on a
  metadata field or a region, see [bki-bkg](docs/bki-bkg.md)

Python (NumPy + Matplotlib + Plotly), no LArSoft required.  PyROOT is needed
for the `distrib_x/y/z` formulas, which are ROOT TF1 expressions and are
evaluated by ROOT (see [bki-distrib](docs/bki-distrib.md)).  Optionally, the GDML can
also be read with ROOT instead of the built-in parser (see
[GDML backends](docs/bki-extract.md#gdml-backends)).

---
## Installation

Every command line needs **click** and **rich**, besides NumPy.  The minimal
install, for a Python that provides nothing else, is

```bash
pip install --no-deps . numpy click rich
```

(`--no-deps` leaves matplotlib, plotly and the notebook tools out; add what
you need from the table below.)  The full installs:

```bash
cd dune-radgeo
uv sync                  # everything, including ROOT (pip wheel, Python 3.14)
uv sync --no-group root  # everything but ROOT
pip install .            # everything but ROOT
```

PyROOT is needed to evaluate `distrib_x/y/z` formulas (`bki-extract` on a FCL
with such producers, `bki-distrib`, `tests/test_tf1.py`) and for
`bki-extract --gdml-backend root`.  On the host it comes from the `root` pip
wheel, a uv dependency group installed by default, not a package
dependency: `pip` never installs it, so a dunesw area uses the stack's ROOT,
never shadowed by the wheel.  The group's marker limits it to Python ≥ 3.14
(the host venv).  The current wheels (0.1a12) are built without TGeo/GDML,
so for now `--gdml-backend root` works only with a full ROOT such as the
dunesw one (TF1, used for the formulas, is in every build).  The wheel
prints an "alpha stage" warning on import; duneradgeo silences it.  Without
PyROOT, the features that need it raise `ROOTUnavailable` (see
[Python API](#python-api)).

`bki-extract` of a geometry, or of producers without `distrib_*`, needs only
NumPy, click and rich.  For an extraction-only install on top of a LArSoft
Python that already provides NumPy and PyROOT, install without dependencies
and add whichever of NumPy, click and rich is missing:

```bash
pip install --no-deps . numpy click rich
```

| Package | Used by |
|---------|---------|
| `numpy` | everything |
| `click` | the command line of every `bki-*` command |
| `rich` | the messages and tables of every `bki-*` command (`bki-bkg --json` prints plain JSON) |
| `matplotlib` | `bki-plot` and `bki-distrib` PDF/PNG output, `plot_ir` |
| `plotly` | `make_figure` |
| `jupyter`, `nbconvert` | `nb_B.ipynb`, `scripts/run_notebook.py` |
| `root` (group) | `distrib_*` formulas: `bki-extract` of producers with `distrib_*`, `bki-distrib`, `tests/test_tf1.py`.  `bki-extract --gdml-backend root` and `scripts/compare_backends.py` need a ROOT with TGeo/GDML (dunesw) |

`bki-plot --list` / `-L` and `bki-distrib --list` need no matplotlib.  In a
`--no-deps` install, rendering a plot without matplotlib stops with an
`ImportError` that says what to install.

### In a dunesw area

To run `bki-extract` with the area's `$FHICL_FILE_PATH`, `$FW_SEARCH_PATH` and
`fhicl-dump`, install it into the area's venv with `scripts/dunesw_installme.sh`.
The script does an editable `--no-deps` install, and NumPy comes from the
stack. The commands need click and rich: the dunesw area has `rich` but not
`click`, so the script also installs whichever of the two the venv cannot
import. With `dsw`:

```bash
dsw -C <area> exec bash /path/to/dune-radgeo/scripts/dunesw_installme.sh
dsw -C <area> exec bki-extract <job>.fcl --out /abs/path/detector.json
```

Add `--no-editable` for a regular install. Any other arguments are passed to
pip. The area needs a `.venv`, which `dsw newarea` creates.

The script exists because a plain `pip install -e .` fails in a UPS
environment with "editable mode currently requires a setuptools-based
build". The UPS product `python_future_six_request` puts pip 21.2.4 on
`PYTHONPATH`, and that pip is used instead of the venv's own, even with
`python -m pip`. The script clears `PYTHONPATH` for pip only.

Without installing anything, you can also run from the source tree:
`dsw -C <area> exec env PYTHONPATH=/path/to/dune-radgeo/src python -m duneradgeo.extract ...`.

The area's ROOT has TGeo/GDML, so `--gdml-backend root` and
`scripts/compare_backends.py` work there (`dsw -C <area> exec bki-extract
geometry.gdml -o /abs/ir.json --gdml-backend root`).

`bki-plot` and `bki-distrib` plots need matplotlib, which the area does not
have: plot the JSON on the host with `uv run bki-plot` / `uv run bki-distrib`.
`bki-distrib --list` works in the area (it uses the area's ROOT).
If the script cannot install click (no network), do it by hand with
`env -u PYTHONPATH python -m pip install click rich` in the area.

After pulling a version that adds a command, re-run the install script: an editable install only creates the
command scripts at install time (`python -m duneradgeo.<module>` works
regardless).

---

## Quick start

### 1. Extract volumes → JSON IR

```bash
# From a GDML file (detector geometry only)
bki-extract geometry.gdml --out detector.json

# From a FHiCL job file (radio volumes + geometry via AuxDetGeometry.GDML)
bki-extract job.fcl --out detector.json
```

Bare filenames are resolved via the LArSoft search paths
([File search paths](#file-search-paths)).  `#include` expansion and its
cache, pattern overrides, `--frame`, `--with-passive` and the GDML backends are
in [bki-extract](docs/bki-extract.md).

### 2. Render static PDF

```bash
bki-plot detector.json --out detector.pdf
bki-plot detector.json --out detector.pdf --view 30 -45
```

A 3D view and three 2D views; selectors, highlighting and listings are in
[bki-plot](docs/bki-plot.md).

### 3. Radiological position profiles

```bash
bki-distrib detector.json            # detector_distrib.pdf, one page per distribution
bki-distrib detector.json --list     # the same as text
```

See [bki-distrib](docs/bki-distrib.md).

### 4. Inspect the backgrounds

```bash
bki-bkg radbkg.json list                       # every background
bki-bkg radbkg.json select nuclide=Po218       # on a metadata field
bki-bkg radbkg.json region --slice x=-328      # active on a plane
bki-bkg radbkg.json suggest                    # other selections worth trying
```

See [bki-bkg](docs/bki-bkg.md).

### 5. Interactive notebook

```bash
# Run from the command line (executes nb_B.ipynb non-interactively)
python scripts/run_notebook.py detector.json
python scripts/run_notebook.py detector.json --out executed.ipynb

# Or open in Jupyter and set IR_FILE manually
jupyter lab notebooks/nb_B.ipynb
```

The notebook reads `$BKI_IR_FILE` if set; otherwise use the `IR_FILE` variable
in the first cell.  It renders an interactive Plotly figure (3D + 3×2D) and
can also save a static matplotlib PDF.

---

## Commands

| Command | What it does | Documentation |
|---------|--------------|---------------|
| `bki-extract` | GDML / FHiCL → JSON IR | [docs/bki-extract.md](docs/bki-extract.md) |
| `bki-plot` | IR → static PDF/PNG, listings, selectors, highlighting | [docs/bki-plot.md](docs/bki-plot.md) |
| `bki-measure` | geometric measurements on an IR | [docs/bki-measure.md](docs/bki-measure.md) |
| `bki-distrib` | `distrib_x/y/z` position profiles | [docs/bki-distrib.md](docs/bki-distrib.md) |
| `bki-bkg` | list, select and locate the backgrounds of an IR | [docs/bki-bkg.md](docs/bki-bkg.md) |

Each page opens with a description and usage, then the command's `--help`
(checked against the code by `tests/test_docs.py`), then its examples and
explanations.  Sections that used to
be here:

| Section | Now in |
|---------|--------|
| Usage examples: extracting an IR | [bki-extract](docs/bki-extract.md#examples) |
| perl vs ggd VD GDML | [bki-extract](docs/bki-extract.md#perl-vs-ggd-vd-gdml) |
| GDML backends | [bki-extract](docs/bki-extract.md#gdml-backends) |
| FHiCL expansion and cache | [bki-extract](docs/bki-extract.md#fhicl-input) |
| Inspecting an IR (`-l` / `-L`) | [bki-plot](docs/bki-plot.md#inspecting-an-ir) |
| Plotting the detector, backgrounds, generation surfaces | [bki-plot](docs/bki-plot.md#plotting-the-detector) |
| `--surface` syntax | [bki-plot](docs/bki-plot.md#--surface-syntax) |
| Picking individual placements | [bki-plot](docs/bki-plot.md#picking-individual-placements) |
| Highlighting volumes | [bki-plot](docs/bki-plot.md#highlighting-volumes) |
| Adding a measurement | [bki-measure](docs/bki-measure.md#adding-a-measurement) |
| The formula and the generated distribution (Npx) | [bki-distrib](docs/bki-distrib.md#the-formula-and-the-generated-distribution-npx) |

---

## JSON IR format

`bki-extract` writes a JSON file with this structure:

```json
{
  "version": "1",
  "source": "path/to/geometry.gdml",
  "detector": "vd",
  "frame": {"volume": "volEnclosureTPC",
            "path": "volDetEnclosure:0/volCryostat:0/volEnclosureTPC:0",
            "pos": [0.08, 0.0, 3003.55], "rot": [[1,0,0],[0,1,0],[0,0,1]]},
  "producers": {
    "Ar39GenInLAr": {"module_type": "Decay0Gen", "isotope": "Ar39",
                     "material": "LAr", "volume_rand": "volCryostat",
                     "BqPercc": 0.00141}
  },
  "volumes": [
    {
      "kind": "anode",
      "name": "volAnodePlate",
      "pos":  [325.08, -505.10, 148.92],
      "ext":  [0.01, 335.40, 297.84],
      "rot":  [[1,0,0],[0,1,0],[0,0,1]],
      "meta": {"path": "volDetEnclosure:0/volCryostat:0/volAnodePlate:0"}
    },
    {
      "kind": "radio",
      "name": "K42From42ArGenInUpperMesh1x8x14",
      "pos":  [-320.917, 0.0, 1045.94],
      "ext":  [4.329, 1507.84, 2306.12],
      "rot":  [[1,0,0],[0,1,0],[0,0,1]],
      "meta": {"source": "volume_rand:volCryostat", "nuclide": "K42",
               "module_type": "Decay0Gen", "rate": 1.94372,
               "surface": "plane_x@-320.917",
               "distrib_x": "(7.2257e+4*TMath::Exp(-(x+320.)*(x+320.)/(2.*0.5*0.5)))+1.",
               "distrib_x_range": [-425.12, 425.20], "distrib_x_eff_range": [-323.082, -318.753],
               "distrib_x_formula_eff_range": [-321.704, -318.328],
               "distrib_eff_containment": 0.99, "distrib_eff_npx": 100, "index": 0}
    }
  ]
}
```

| Field  | Description |
|--------|-------------|
| `kind` | `anode`, `opdet`, `cathode`, `fieldshaper`, `struct`, or `radio` |
| `pos`  | centre `[x, y, z]` in cm, in the IR frame |
| `ext`  | full extents `[dx, dy, dz]` in cm |
| `rot`  | 3×3 rotation matrix in the IR frame |
| `meta` | kind-specific metadata: GDML placement `path` for geometry volumes (see [Picking individual placements](docs/bki-plot.md#picking-individual-placements)) and, when the GDML names its physvols (ggd), the physvol name `pv` and the physvol names along the path `pvpath` (`placeDetEnclosure/…/volTPC_PV_Bot0/volTPCPlaneU_PV`); nuclide, BqPercc, surface, index, … for radio |

For radio entries, `meta` also holds, for each axis with a `distrib_<axis>`
formula: the formula, `distrib_<axis>_range` (the full box range it is
sampled over, the source of truth), `distrib_<axis>_eff_range` and
`distrib_<axis>_formula_eff_range`, plus `distrib_eff_containment` and
`distrib_eff_npx` (see [the formula and the generated
distribution](docs/bki-distrib.md#the-formula-and-the-generated-distribution-npx)). The
effective ranges are derived by `bki-extract` with ROOT, not module
parameters. These are world coordinates. The box (`pos`/`ext`) is the full
box on the axes without a formula and the generated effective range on the
others.

`producers` is present for FCL input. It holds the complete FHiCL table of
every selected radiological producer, keyed by module label, and is stored
even when none of the producer's boxes resolve. The `name` of a radio volume
is its producer label. ±inf values are stored as the strings `"inf"` and
`"-inf"`.

`frame` is present only for IRs written with [`bki-extract --frame VOLNAME`](docs/bki-extract.md).
The IR frame is then that volume's local frame; `frame.pos` / `frame.rot`
are its world placement (`world = pos + rot @ local`).  Without `frame`,
the IR frame is the world.  `meta.path` always starts at the world volume.

Every entry is a box.  Non-box solids are approximated: a tube or torus by
its bounding box, a union by one box per operand, a subtraction or
intersection by its first (outer) operand.  A volume that splits into several
boxes gives one entry per box, all with the same `name`, plus
`meta.copy` (placement index) and `meta.part` (the primitive solid).  For
example, each VD field-shaper ring becomes 8 (or 12, for the slim rings)
tube and corner boxes.  `--list` counts placements and shows the box count
separately; `-L` prints `name#copy` and the part.

Axis convention: **x** = drift, **z** = beam; **y** is vertical in HD and
horizontal in VD (where the drift, x, is vertical).

### Radio boxes

A radio entry's box is not decomposed: it is the box larsimrad's
`BaseRadioGen` samples in, which it finds without looking at the shape
type.

- `X0..X1`, `Y0..Y1`, `Z0..Z1`: that box, as given.
- `volume_rand`, `volume_gen`: for every shape, the ROOT bounding box of the
  volume's solid (`TGeoBBox` origin and half-lengths), so a tube or torus
  gets its bounding box and a union the box of the whole composite, not of
  one operand (a subtraction: of its first operand; an intersection: the
  overlap of the operands' boxes).  That box is placed in the world through
  the mother chain, and the entry is the axis-aligned envelope of its 8
  corners.  As in `BaseRadioGen`, `volume_rand` is the exact volume name
  and gives one entry, for its first placement in depth-first order from
  the world (GDML daughter order); `volume_gen` is a regex matched against
  the full volume name (`volAnodePlate` does not select
  `volAnodePlateBottom`) and gives one entry per matching node
  (`meta.node`).
- The ROOT backend takes ROOT's own box for every shape; the built-in
  parser computes it for boxes, tubes, tori and their booleans, and warns
  about a volume whose solid it cannot box.
- Then each axis with a `distrib_<axis>` formula is narrowed to its
  effective range.

Like every entry, the box is written in the IR frame (that of the
`--frame` volume, if given); the `distrib_*` ranges stay in world
coordinates.

---

## Python API

The public entry points are the package top level (`duneradgeo.IRCollection`,
`GDMLGeometry`, `FhiclParser`, …), `duneradgeo.ir`, `plot.plot_ir`,
`notebook_plot.make_figure`, `distrib` and `measurements`.  Modules with a
leading underscore are internal: the examples use some of them, but their
names and signatures may change.

```python
from duneradgeo.ir import IRCollection
from duneradgeo.notebook_plot import make_figure
from duneradgeo.plot import plot_ir
from duneradgeo._select import select_volumes      # internal

ir = IRCollection.load("detector.json")

# filter by kind
anodes = ir.by_kind("anode")
radios = ir.by_kind("radio")

# interactive Plotly figure
fig = make_figure(ir)
fig.show()

# static matplotlib PDF (same layout as bki-plot)
plot_ir(ir, "detector.pdf")
plot_ir(ir, "detector_stretched.pdf", proportional=False)
plot_ir(ir, "detector_vd.pdf", orientation="vd")   # drift axis x vertical

# plot with selectors
plot_ir(ir, "radio_only.pdf", kinds=["radio"],
        kind_patterns={"radio": ["GenInCathode"]})

# selected generation surfaces only
plot_ir(ir, "surfaces.pdf", surfaces=["GenInCathode#0,2"])

# selected placements of a geometry volume, by GDML path
plot_ir(ir, "tpc5.pdf", picks=["^volTPC$#@volTPC:5$"])

# highlight a subset, the rest drawn as grey context
plot_ir(ir, "cathode.pdf", highlight=["cathode=^volCathodeGrid$"],
        dim=0.02, zoom_highlight=False)
make_figure(ir, highlight=["cathode=^volCathodeGrid$"]).show()

# same selection as a volume list (e.g. for make_figure)
sub = IRCollection(source=ir.source, detector=ir.detector,
                   volumes=select_volumes(ir, surfaces=["GenInCathode#@plane_x"]))
make_figure(sub).show()
```

Position profiles (`bki-distrib`; needs PyROOT):

```python
from duneradgeo.distrib import collect, fractions, plot_groups, profile, with_npx

groups = collect(select_volumes(ir, kinds=["radio"]))   # one DistribGroup per formula + range
g = groups[0]           # g.axis, g.formula, g.range_, g.eff, g.formula_eff, g.npx, g.producers
pr = profile(g)         # pr.x, pr.formula, pr.generated: both pdfs over g.range_
fractions(g)            # (generated, formula): measured fractions inside g.eff, g.formula_eff
g2 = with_npx(g, 10000) # a copy, g2.eff recomputed for a generator with SetNpx(10000)
plot_groups(groups, "profiles.pdf", log=True)           # one page per group
```

`g.dist` is the group's formula as a `Distribution` of the internal module
`duneradgeo._tf1`: one formula on one range with one Npx, which owns its
TF1 (see [the formula and the generated
distribution](docs/bki-distrib.md#the-formula-and-the-generated-distribution-npx)).
`_tf1.distribution(formula, lo, hi, npx=100)` returns the same object for
the same arguments (the 256 most recent are cached;
`distribution.cache_clear()` drops them).

```python
from duneradgeo import _tf1                      # internal

d = g.dist                       # = _tf1.distribution(g.formula, *g.range_, g.npx)
d.eval([-320.0, 0.0])            # TF1::Eval
d.formula_interval(0.99)         # shortest interval holding 99% of the formula
d.interval(0.99)                 # ... of what GetRandom generates (Npx = d.npx)
d.quantiles([0.16, 0.5, 0.84])   # TF1::GetQuantiles on the Npx table
d.cdf([-320.0]), d.pdf(pr.x)     # the generated distribution
d.sample(100000)                 # TF1::GetRandom, as the generator draws
d.with_npx(10000).interval(0.99)
_tf1.distribution("TMath::Gaus(x,100,1)", -425.12, 425.20).interval(0.99)
```

Errors:

- `ROOTUnavailable` (`duneradgeo._root`, a `RuntimeError`): PyROOT cannot
  be imported where it is needed (a `distrib_*` formula, the ROOT GDML
  backend), or it has no TGeo/GDML (`--gdml-backend root` with the pip
  wheel).
- `ValueError`: ROOT cannot parse a formula, or the formula integrates to 0
  over its range.
- `ExtractError` (`duneradgeo.extract`): extraction cannot go on
  (input or `--gdml` not found, GDML unreadable, a `--frame` volume not
  placed exactly once, a `distrib_*` producer without PyROOT, nothing
  extracted).  The extraction steps (the internal `_gdml_volumes`,
  `_fcl_radio_volumes`, `_to_frame`) raise it instead of exiting, so they
  can be called from Python; `bki-extract` prints it as `[error] …` and
  exits with status 1.

Internal modules shared by the commands:

| Module | Contents |
|--------|----------|
| `duneradgeo._select` | the selectors of `bki-plot`, `bki-distrib` and `make_figure`: `select_volumes`, `parse_pick_spec`, `parse_highlight_spec`, `validate_picks`, `highlight_groups` |
| `duneradgeo._style` | colours and palettes (`KIND_COLOUR`, `STRUCT_PALETTE`, `SURFACE_PALETTE`, `CONTEXT_COL`, …) and `import_mpl()`, the matplotlib import with an install hint |
| `duneradgeo._draw` | `aabb(pos, ext, rot)`: the (lo, hi) envelope of a rotated box's 8 corners; `box_corners_world` |
| `duneradgeo._tf1` | `Distribution`, `distribution()`, `GENERATOR_NPX` (100) |
| `duneradgeo._root` | `import_root(feature)`, the shared PyROOT import (batch mode, wheel warning silenced); `ROOTUnavailable` |

Parsing utilities:

```python
from duneradgeo import FhiclParser, GDMLGeometry, load_geometry
from duneradgeo._fcl import expand_fcl                     # internal

geo = load_geometry("geometry.gdml")             # GDMLGeometry; backend="root" for ROOT's reader
geo.find_positions([r"^volCryostat$"])           # [(name, world pos, world rot), …]
doc = FhiclParser(expand_fcl("job.fcl")).parse() # the FHiCL document as nested dicts
```

---

## Detector profiles

A profile is a detector plus a GDML flavour: `vd-perl` and `vd-ggd` are the
same VD detector, described by perl- and ggd-generated GDML, and both write
`"detector": "vd"` to the IR.  Auto-detected from logical volume names (`vd-ggd` is recognised by
`volSupportEnc` / `volEnclosureCryostat` and checked before `vd-perl`):

| Key  | Geometry | Anode | Opdet | Cathode | Field shaper | Default struct |
|------|----------|-------|-------|---------|--------------|----------------|
| `vd-perl` | Vertical Drift (FD2), perl GDML | `volAnodePlate`, `volAnodePlateBottom` | `volArapuca*` | `volCathodeGrid` (`volGroundGrid` in v4) | `volFieldShaper`, `volFieldShaperSlim` | `volTPCActive`, `volGaseousArgon`, `volTPC`, `volCryostat` |
| `vd-ggd` | Vertical Drift (FD2), ggd GDML | as `vd-perl` | as `vd-perl` | `volCathodeGrid` | as `vd-perl` | as `vd-perl` |
| `hd` | Horizontal Drift / ProtoDUNE-HD | `volAPAFrameYSide` | `volArapuca_N` | — | — | `volTPCActiveInner`, `volGaseousArgon`, `volTPC`, `volCryostat` |

Override with `--anode-pattern` / `--opdet-pattern` / `--cathode-pattern` /
`--fieldshaper-pattern` / `--struct-pattern`.

`--with-passive` adds each profile's passive material as `struct` volumes:

| Key  | Passive volumes |
|------|-----------------|
| `vd-perl` | `volSteelShell`, `volFoamPadding`, `volSteelSupport` |
| `vd-ggd` | `volSteelShell`, `fShellLog`, `WoodLog`, `FoamLog`, `ShellOutLog` (cryostat and insulation); `IBeamTop/Bot/Side`, `BeltUni`, `BeltUniTop`, `BeltHoleUni` (support); `ShieldBlockContainerLog`, `boxshapeVol`, `boxshapeLeadVol` (shielding); `volShotbox`, `volGrout`, `volConcrete`, `volRadioRockShell` (cavern) |
| `hd` | — |

The VD cathode solid is a box with 16 holes cut out (a GDML subtraction);
it is drawn as its outer box.  The cathode-mounted photon detectors
(`volArapucaCathode`) stay in the `opdet` kind.

---

## File search paths

| File type      | Environment variable  |
|----------------|-----------------------|
| FHiCL (`.fcl`) | `$FHICL_FILE_PATH`    |
| GDML geometry  | `$FW_SEARCH_PATH`     |

Both mirror LArSoft's `cet::search_path` behaviour (first-match-wins).

---

## Tests

| File | Checks | Needs |
|------|--------|-------|
| `tests/test_tf1.py` | `duneradgeo._tf1`: the generated distribution against a `lar` Decay0Gen reference, the formula, caching and TF1 ownership | PyROOT |
| `tests/test_fcl.py` | the FHiCL parser, `#include` expansion and `producer_volumes` (on a small GDML of boxes) | nothing (ROOT-free) |
| `tests/test_formula_lines.py` | the formula text layout of the `bki-distrib` pages | nothing (ROOT-free) |
| `tests/test_select.py` | the `--pick` / `--surface` / `--highlight` specs and `select_volumes` | nothing (ROOT-free) |
| `tests/test_cli.py`, `tests/test_extract_cli.py`, `tests/test_plot_cli.py`, `tests/test_measure_cli.py`, `tests/test_distrib_cli.py`, `tests/test_scripts_cli.py` | the click/rich command lines of the `bki-*` commands and scripts | click and rich (skipped without them) |
| `tests/test_docs.py` | the `--help` text of each command is in its `docs/` page; the links between README and `docs/` resolve | the help check needs click and rich |
| `tests/test_bkg.py` | `bki-bkg`: backgrounds, select expressions, regions, suggestions, `--json` output | nothing; the table rendering needs `rich` (skipped without it) |

Each file runs on its own or under pytest (not a dependency).  On the host:

```bash
.venv/bin/python tests/test_fcl.py        # prints ok per test
.venv/bin/python -m pytest tests/          # if pytest is installed
```

In a dunesw area (after `scripts/dunesw_installme.sh`), with the area's
Python and ROOT:

```bash
dsw -C <area> py tests/test_tf1.py
```

A known issue can be recorded as a test of the correct behaviour marked
`@xfail("reason")` (from `tests/_runner.py`); none is open at present.  Run
directly, such tests are skipped and listed as `xfail`; under pytest they
run as expected failures.

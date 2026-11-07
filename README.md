# duneradgeo

Inspect DUNE background sources without a LArSoft environment:

- **detector geometry** — anode (CRP/APA), optical detector, cathode,
  field-shaper, and structural volume positions from GDML files
- **radiological volumes** — emission boxes from RadioGen / Decay0Gen FHiCL
  configurations, with each producer's full FHiCL table
- **radiological position profiles** — the `distrib_x/y/z` formulas of those
  producers, drawn over the range the generator samples

Python (NumPy + Matplotlib + Plotly), no LArSoft required.  PyROOT is needed
for the `distrib_x/y/z` formulas, which are ROOT TF1 expressions and are
evaluated by ROOT (see [bki-distrib](#bki-distrib)).  Optionally, the GDML can
also be read with ROOT instead of the built-in parser (see
[GDML backends](#gdml-backends)).

---

## Installation

```bash
cd dune-radgeo
uv sync                  # everything, including ROOT (pip wheel, Python 3.14)
uv sync --no-group root  # everything but ROOT
pip install .            # everything but ROOT
```

ROOT is needed to evaluate `distrib_x/y/z` formulas (`bki-extract` on a FCL
with such producers, `bki-distrib`) and for `bki-extract --gdml-backend
root`.  It is a uv dependency group installed by default, not a package
dependency, so `pip` never installs it: in a dunesw area the stack's ROOT is used and is not
shadowed by the pip wheel.  The `root` pip wheels exist for CPython 3.14
only (the group is skipped on older Pythons), and the current ones
(0.1a12) are built without TGeo/GDML, so for now `--gdml-backend root`
works only with a full ROOT such as the dunesw one (TF1, used for the
formulas, is in every build).  The pip wheel prints an "alpha stage" warning
on import; duneradgeo silences it.

`bki-extract` itself needs only NumPy, plus PyROOT for producers with
`distrib_x/y/z`.  For an extraction-only install, e.g. on top of a LArSoft
Python that already provides NumPy and PyROOT, skip the dependencies:

```bash
pip install --no-deps .          # add numpy if the environment lacks it
```

| Package | Used by |
|---------|---------|
| `numpy` | everything (`bki-extract` of geometry and producers without `distrib_*`, `bki-plot --list` / `-L` need nothing else) |
| `matplotlib` | `bki-plot` and `bki-distrib` PDF/PNG output, `plot_ir` |
| `plotly` | `make_figure` |
| `jupyter`, `nbconvert` | `nb_B.ipynb`, `scripts/run_notebook.py` |
| `root` (group) | `distrib_*` formulas: `bki-extract` of producers with `distrib_*`, `bki-distrib`, `tests/test_tf1.py`; `bki-extract --gdml-backend root`, `scripts/compare_backends.py` (need TGeo/GDML) |

In a `--no-deps` install, rendering a plot without matplotlib stops with an
`ImportError` that says what to install.

### In a dunesw area

To run `bki-extract` with the area's `$FHICL_FILE_PATH`, `$FW_SEARCH_PATH` and
`fhicl-dump`, install it into the area's venv with `scripts/dunesw_installme.sh`.
The script does an editable `--no-deps` install, and NumPy comes from the
stack. With `dsw`:

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
`bki-distrib --list` works in the area (it uses the area's ROOT). After pulling a version that adds a
command, re-run the install script: an editable install only creates the
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

# Bare filenames are resolved via the LArSoft search paths:
#   FHiCL → $FHICL_FILE_PATH
#   GDML  → $FW_SEARCH_PATH
bki-extract prodbackground_radio.fcl --out detector.json
```

FHiCL files with `#include` directives are expanded automatically:
`fhicl-dump` is used when available; otherwise includes are resolved inline
via `$FHICL_FILE_PATH`. The expanded text is cached in
`~/.cache/bki-extract/fcl/` (or `$XDG_CACHE_HOME/bki-extract/fcl/`) and
reused on later runs. The cache is keyed on the top-level file's path, size
and mtime plus `$FHICL_FILE_PATH`; edits to *included* files are not
detected, so pass `--no-fcl-cache` after changing them. `--fcl-cache DIR`
moves the cache elsewhere.

Override patterns or supply the GDML explicitly:

```bash
bki-extract geometry.gdml --out ir.json \
    --anode-pattern '^volAnodePlate$' \
    --opdet-pattern '^volArapuca' \
    --struct volTPCActive volGaseousArgon volTPC volCryostat

bki-extract job.fcl --gdml /path/to/geometry.gdml --out ir.json
```

### 2. Render static PDF

```bash
bki-plot detector.json --out detector.pdf
bki-plot detector.json --out detector.pdf --view 30 -45
```

Layout: 3D ortho view (proportional axes) + Z–Y beam face + Z–X top view +
X–Y side view.  The axis ranges enclose every drawn volume (plus a 3%
margin), so use the selectors to zoom in.  The 2D views are drawn to scale
by default and share the same axis ranges; `--no-proportional` stretches each one to fill its panel,
which helps when the detector is much longer along one axis than another.

`--orientation` sets which axis points up.  Each panel keeps its role; only
the axes change:

| Panel        | `hd` (default)        | `vd` (drift vertical) |
|--------------|-----------------------|-----------------------|
| top-left     | 3D, y up              | 3D, x up              |
| top-right    | Z–Y  beam face        | Z–X  beam face        |
| bottom-left  | Z–X  top view         | Z–Y  top view         |
| bottom-right | X–Y  side view        | Y–X  side view        |

### 3. Radiological position profiles

```bash
bki-distrib detector.json            # detector_distrib.pdf, one page per profile
bki-distrib detector.json --list     # the same as text
```

See [bki-distrib](#bki-distrib).

### 4. Interactive notebook

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

## bki-extract options

```
bki-extract <input> --out <ir.json>

  input               GDML geometry file or expanded FHiCL job file

  --gdml PATH         override the GDML path (when input is FCL)
  --gdml-backend {python,root}
                      GDML reader: built-in parser (default) or ROOT's
                      TGeoManager::Import (see GDML backends)
  --no-geometry       skip GDML parsing (radio volume_rand/gen stay unresolved)

  --anode-pattern P   regex(es) for anode volume names   (auto-detected)
  --opdet-pattern P   regex(es) for opdet volume names   (auto-detected)
  --cathode-pattern P regex(es) for cathode volume names (auto-detected, VD only)
  --fieldshaper-pattern P
                      regex(es) for field-shaper volume names (auto-detected, VD only)
  --struct-pattern P  regex(es) for structural volumes   (auto-detected)
  --struct VOLNAME    exact structural volume names (anchored, overrides --struct-pattern)
  --with-passive      also extract the profile's passive material (steel shell,
                      insulation, support structure, shielding, cavern) as struct
  --frame VOLNAME     express positions in the frame of this singly placed volume
                      (e.g. volCryostat, volEnclosureTPC) instead of the world;
                      its world placement is stored in the IR as "frame"

  --prefix PREFIX     only producers whose label starts with PREFIX
  --select REGEX      regex filter on producer labels
  --all-configured    include producers not on any trigger path
  --containment F     distrib_x/y/z containment fraction for the effective
                      ranges, stored as meta.distrib_eff_containment (default 0.99)
  --npx N             TF1 Npx of the generator: the effective range is that of
                      the distribution TF1::GetRandom draws from an N-point
                      table, stored as meta.distrib_eff_npx (default 100, ROOT's
                      default, which larsimrad BaseRadioGen keeps)
  --fcl-cache DIR     cache directory for expanded FCL files
                      (default: ~/.cache/bki-extract/fcl)
  --no-fcl-cache      always re-run the FCL expansion
  --verbose           extra debug output
```

---

## bki-plot options

```
bki-plot <ir.json> [--out <file.pdf>]

  --list, -l          print volume kinds and names grouped by kind, with one
                      line per radio generation surface, then exit
  --list-detail, -L   print every placement with its x/y/z ranges (IR frame:
                      the world unless extracted with --frame) and its path
                      (pvpath when the GDML names its physvols), then exit
                      (both listings respect all selectors below)

  --view ELEV AZIM    3D view angles in degrees (default: 20 -60)
  --title TEXT        plot title
  --[no-]proportional 2D views to scale (default) or stretched to fill
                      their panels
  --orientation O     hd: drift axis x horizontal, y vertical (default)
                      vd: x vertical in every view
                      auto: vd for VD IRs, hd otherwise

volume selectors (applied in order):
  --no-detector       exclude anode / opdet / cathode / fieldshaper / struct volumes
  --no-radio          exclude radio volumes
  --kind K [K …]      keep only these kinds: anode opdet cathode fieldshaper
                      struct radio
  --anode P [P …]     regex filter on anode volume names
  --opdet P [P …]     regex filter on opdet volume names
  --cathode P [P …]   regex filter on cathode volume names
  --fieldshaper P [P …]
                      regex filter on field-shaper volume names
  --struct P [P …]    regex filter on struct volume names
  --radio P [P …]     regex filter on radio producer names
  --name P [P …]      cross-kind name regex
  --surface S [S …]   radio generation-surface selector
  --pick S [S …]      placement selector for any kind, by GDML placement
                      path (applied last)

highlighting (applied after the selectors):
  --highlight S [S …] [LABEL=]NAME_RE[#SEL]: draw matching volumes on top, one
                      colour and legend entry per spec, the rest as grey context
  --dim F             face opacity of the context volumes (default: 0.02)
  --zoom-highlight    fit the axes to the highlighted volumes
```

See [Usage examples](#usage-examples) below.

---

## bki-measure

Geometric measurements on an IR (lengths in cm, in the IR frame — the world
unless the IR was extracted with `--frame`; the measured distances do not
depend on it).

```
bki-measure <ir.json> [MEASUREMENT …]

  MEASUREMENT         measurements to run (default: all)
  --list              list the available measurements and exit
  --detail, -d        also print each measurement's per-row table
  --json FILE         write the results (summary, values, rows) as JSON
```

```bash
bki-measure geo.json
# == drift-gaps: cathode ↔ anode distance along x (drift)  [cm]
#   top     gap   650.080  (min 650.080, max 650.080, 80 anodes)   centre-to-centre   652.085
#           volTPCActive x-extent   650.000  (320 volumes)
#   bottom  gap   650.080  (min 650.080, max 650.080, 160 anodes)   centre-to-centre   652.085
#           volTPCActive x-extent   650.000  (320 volumes)

bki-measure geo.json drift-gaps -d --json gaps.json   # per anode–cathode pair
```

| Measurement  | What it measures |
|--------------|------------------|
| `anode-cryostat-gap` | Along x, from each anode to the cryostat face behind it: top anodes to the cryostat top face, bottom anodes to the bottom face (faces of the `volCryostat` box, i.e. the inner cryostat surface; the top one encloses the gas ullage).  `gap_inner` is measured from the anode face looking into the drift volume, `gap_outer` from the face looking at the cryostat.  Top / bottom is split at the cathode plane (cryostat mid-height without cathodes).  On top, the liquid surface (bottom of `volGaseousArgon`) is shown for reference.  Needs `anode` and `struct volCryostat`. |
| `drift-gaps` | Pairs every anode with the cathode module it faces (largest y–z overlap) and measures along x: `gap` between the facing surfaces (anode bottom − cathode top on the top side, cathode bottom − anode top on the bottom side) and `centre_dist` between centres.  Anodes above the cathode are "top", below "bottom"; per side it reports mean/min/max, flags pairs more than 0.01 cm from the median, and shows the `volTPCActive` x-extent as a cross-check.  Needs `anode` and `cathode` entries in the IR. |
| `opdet-fc-gap` | Wall optical detectors (opdet boxes whose centre is outside the field cage's y–z envelope; the cathode-mounted ones are inside and skipped) ↔ closest field-cage surface: for each detector, the nearest field-shaper part and the shortest box-to-box distance, summarised per wall (`y-`, `y+`, `z-`, `z+`) with the closest tube types.  Tubes are stored as bounding boxes, so the gap to a facing straight tube is exact and to a corner a lower bound.  Needs `opdet` and `fieldshaper` entries in the IR. |
| `fc-cryostat-gap` | The four field-cage side walls (y−, y+, z−, z+; top and bottom are open) ↔ cryostat walls (faces of `volCryostat`).  Each field-shaper piece is assigned to the wall it is nearest; per wall, `gap_min` is the distance from the cryostat wall to the closest piece (the standard tubes) and `gap_max` to the most recessed one (e.g. the slim tubes of the photon-detector window).  `--detail` lists the gaps per wall and piece type.  Needs `fieldshaper` and `struct volCryostat`. |

### Adding a measurement

Measurements live in the `duneradgeo.measurements` package, one module each:

```
src/duneradgeo/
  measure.py                  bki-measure CLI (output, --detail, --json)
  measurements/
    __init__.py               auto-imports every module below
    _core.py                  Result, @measurement registry, box helpers, X / YZ axes
    anode_cryostat_gap.py     anode-cryostat-gap
    drift_gaps.py             drift-gaps
    fc_cryostat_gap.py        fc-cryostat-gap
    opdet_fc_gap.py           opdet-fc-gap
```

To add one, create a module in `src/duneradgeo/measurements/` (name not
starting with `_`) with a function that takes an `IRCollection`, returns a
`Result` and is registered with `@measurement`.  It is imported
automatically, runs by default (in module-name order) and can be selected by
name; the CLI, `--detail` and `--json` need no changes.

```python
# src/duneradgeo/measurements/cathode_thickness.py
"""cathode-thickness: x-extent of the cathode modules."""
from __future__ import annotations

from ..ir import IRCollection
from ._core import X, Result, aabb, measurement, tags


@measurement("cathode-thickness", "x-extent of the cathode modules")
def cathode_thickness(ir: IRCollection) -> Result:
    res = Result("cathode-thickness", "cathode module thickness along x",
                 columns=["cathode", "dx"])
    cathodes = ir.by_kind("cathode")
    for tag, v in zip(tags(cathodes), cathodes):
        lo, hi = aabb(v)                        # world-frame box corners
        res.rows.append(dict(cathode=tag, dx=float(hi[X] - lo[X])))
    dx = [r["dx"] for r in res.rows]
    res.values["dx"] = max(dx) if dx else None  # headline numbers (→ JSON)
    res.summary.append(f"dx {res.values['dx']}")  # lines printed by default
    return res
```

`Result` fields: `summary` (printed lines), `values` (headline numbers,
`None` when not applicable), `columns` + `rows` (the `--detail` table).
Helpers: `aabb(v)` world-frame box corners, `tags(vols)` stable `name#k`
labels, `overlap_area(...)` projected overlap of two boxes,
`box_distances(...)` shortest distance from one box to many,
`struct_box(ir, pattern)` bounding box of matching struct volumes (e.g. the cryostat).  From Python,
`duneradgeo.measurements.run(ir, ["drift-gaps"])` returns the `Result` list.

---

## bki-distrib

Draws the position profiles of the radiological producers that set
`distrib_x`, `distrib_y` or `distrib_z`. In larsimrad's `BaseRadioGen`, each
of these is a 1D TF1 formula in the variable `x`. The generator samples it
with `TF1::GetRandom` over the box range on that axis: `X0..X1` etc., the
`volume_rand` bounding box, or each `volume_gen` node's box. Axes without a
formula are sampled uniformly.

The formulas are ROOT TFormula expressions, so they are evaluated by ROOT
(`duneradgeo._tf1`), never translated. Two distributions are described:

- **the formula**: `TF1::Eval`, normalised over the range;
- **as generated**: `GetRandom` draws from a cumulative table of `Npx`
  points with a parabolic interpolation inside each bin. `BaseRadioGen` does
  not call `SetNpx`, so ROOT's default of 100 applies: over the 850 cm of
  `volCryostat` the bins are 8.5 cm wide, and a narrower feature is not
  reproduced. The distribution is read from ROOT's own table
  (`TF1::GetQuantiles`, at most `Npx / 2` probabilities per call so that ROOT
  keeps the `Npx`-point table), without sampling.

For the upper-mesh producers (a Gaussian of σ = 0.5 cm at x = −320 cm plus
a constant), the generated distribution is a ramp from −323.08 to −318.75 cm
instead of the Gaussian. A `lar` job with Decay0Gen confirms this: its 2.4M
decay vertices match the `GetRandom` table to 0.01 cm (`tests/test_tf1.py`
holds the reference quantiles). Smooth profiles, such as the bulk LAr ones,
are not affected.

`bki-extract` records, per axis with a formula, the full range
(`distrib_<axis>_range`), the shortest interval holding `--containment` of the
generated distribution (`distrib_<axis>_eff_range`, the effective range, which
becomes the box on that axis) and the same for the formula
(`distrib_<axis>_formula_eff_range`), plus `distrib_eff_containment` and `distrib_eff_npx`.

`bki-distrib` draws, over the full range, the formula (black) and the
generated distribution (blue), shades the generated effective range and
marks the formula's (dotted), one distribution per page of the PDF. When the
effective ranges are much narrower than the full range, a zoomed view on
them is added below. Producers that share a formula and range share a page,
which lists them. Coordinates are world coordinates, which is what the
generator uses, even for IRs written with `--frame`.

```
bki-distrib <ir.json> [options]

  --out, -o FILE      output (default: <stem>_distrib.pdf). A PDF has one page per
                      distribution; other formats give one file per
                      distribution, <stem>_NN<suffix>
  --list, -l          print the distributions, both effective ranges with their
                      measured containment, and the producers, then exit
                      (needs ROOT, not matplotlib)
  --radio PATTERN …   regex filter on radio producer names
  --pick SPEC …       PRODUCER_RE#SEL, as for bki-plot --surface
  --axis {x,y,z} …    only these axes
  --npx N             Npx of the generated distribution (default: distrib_eff_npx from
                      the IR); another value recomputes its effective range,
                      e.g. --npx 10000 for a generator that sets SetNpx(10000)
  --log               logarithmic y axis
  --no-zoom           no zoomed view on the effective ranges (default: added
                      below the plot when they are narrower than half the
                      full range; window = their span ± 25%)
  --title TEXT        figure title
```

```bash
bki-distrib radbkg.json --log
bki-distrib radbkg.json --list --radio 'UpperMesh'
# [0] distrib_x  range [-425.120, 425.200]  containment 0.99  npx 100
#     eff_range (generated) [-323.082, -318.753]  (measured 0.9900)
#     formula eff range     [-321.704, -318.328]  (measured 0.9900)
#     surface: plane_x@-320.917
#     formula: (7.2257e+4*TMath::Exp(-(x+320.)*(x+320.)/(2.*0.5*0.5)))+1.
#     K42From42ArGenInUpperMesh1x8x14  (K42)
#     ...
```

Each page has a formula box with two captioned monospace blocks: the FHiCL
parameter exactly as ROOT evaluates it (`distrib_x: "…"`, wrapped between
top-level terms, nothing added or dropped), and the same formula compacted
for reading only (`f(x) = …`, one term per line, `TMath::` dropped, `Abs(a)`
as `|a|`, `Power(a,n)` as `a^n`, `*` as `·`). A side panel has the full and
both effective ranges, the containment (stored, and measured for each
distribution), `Npx`, the generator volume (`meta.source`) and the producers
with their nuclides.

Not modelled: the generator also redraws points that fall outside the
volumes allowed by `volume_gen` / `material`. With `material: LAr` in the
`volCryostat` box, nothing is generated in the gas above the liquid (x >
325.09 cm in the 1x8x14 VD geometry), and the profile is renormalised over
the liquid.

IR files must carry the keys above: re-run `bki-extract` for older ones.

---

## Usage examples

### Extracting an IR

```bash
# geometry only, detector profile auto-detected from the GDML
bki-extract geometry.gdml -o geo.json

# full background job: radio producers on the trigger paths + geometry
bki-extract prodbackground_radio.fcl -o radio.json

# only the producers whose label starts with "GenIn" ...
bki-extract prodbackground_radio.fcl -o cathode.json --prefix GenIn

# ... or matching a regex, including producers not on any trigger path
bki-extract prodbackground_radio.fcl -o rn.json --select 'Rn222|Po210' --all-configured

# use a different geometry than the one named in AuxDetGeometry.GDML
bki-extract prodbackground_radio.fcl -o radio.json --gdml /path/to/geometry.gdml

# narrower effective ranges for producers with a distrib_x/y/z profile (95% containment)
bki-extract prodbackground_radio.fcl -o radio.json --containment 0.95

# an included .fcl was edited: bypass the expansion cache
bki-extract prodbackground_radio.fcl -o radio.json --no-fcl-cache
```

### Refactored vs ggd-generated VD geometries

Both GDML flavours use the same detector volume names (the ggd one is
detected as profile `vd-ggd`, written to the IR as detector `vd`), and each
geometry volume carries its placement `path` and, for ggd, its physvol name
`pv`.  The ggd geometry nests the cryostat in `volSupportEnc` /
`volEnclosureCryostat`, so world coordinates differ by a global offset:
extract both in a common frame and compare them entry by entry.

```bash
bki-extract dunevd10kt_3view_30deg_v7_refactored_2x8x40.gdml -o ref.json --frame volEnclosureTPC
bki-extract dunevd10kt_v7_full10kt_ggd.gdml                 -o ggd.json --frame volEnclosureTPC
python scripts/compare_ir.py ref.json ggd.json -v     # exit status 1 if any group differs

# ggd only: cryostat shell, insulation, I-beams, belts, shielding, cavern
bki-extract dunevd10kt_v7_full10kt_ggd.gdml -o ggd_passive.json --with-passive
```

### GDML backends

`--gdml-backend` picks how the GDML is read:

| Backend | Reader | Needs |
|---------|--------|-------|
| `python` (default) | built-in XML parser (`GDMLGeometry`) | NumPy |
| `root` | `TGeoManager::Import` (`RootGDMLGeometry`) | ROOT with TGeo/GDML (dunesw) |

Both fill the same tables (volumes and placements, box/tube/torus bounding
boxes, boolean operands) and share the search, so placement paths, physvol
names, frames and box decomposition are identical; only the parsing
differs.  ROOT evaluates GDML expressions and units itself and computes
each solid's bounding box.  Solids other than box, tube, torus and their
booleans are skipped by both.  Reading with ROOT takes a few seconds
instead of a fraction of one.

`scripts/compare_backends.py` loads one or more GDMLs with both backends
and compares volumes, solid tables and every profile placement (paired by
path), exiting with status 1 on a difference:

```bash
dsw -C <area> py /path/to/dune-radgeo/scripts/compare_backends.py geometry.gdml [more.gdml] [-v]

# or through the IR
dsw -C <area> exec bki-extract geometry.gdml -o /abs/py.json --gdml-backend python
dsw -C <area> exec bki-extract geometry.gdml -o /abs/rt.json --gdml-backend root
dsw -C <area> exec python /path/to/dune-radgeo/scripts/compare_ir.py /abs/py.json /abs/rt.json
```

The two backends agree exactly on the VD refactored (2x8x40, 1x8x14 v5
and v7a), VD ggd (with and without wires) and HD refactored (1x2x6)
geometries.

### Inspecting an IR

```bash
# kinds, names and placement counts; radio producers list each surface
bki-plot radio.json --list

# per-placement world-frame ranges, restricted to structural volumes
bki-plot radio.json -L --kind struct

# where exactly are the cathode generation surfaces?
bki-plot radio.json -L --kind radio --radio GenInCathode
```

Every selector below also works with `--list` / `-L`, so a selection can be
checked before plotting it.

### Plotting the detector

```bash
# everything, output next to the IR (radio.pdf)
bki-plot radio.json

# detector geometry only, custom title and view angle
bki-plot radio.json --no-radio --title "FD2 VD geometry" --view 30 -45 -o geometry.pdf

# anodes and active volume only
bki-plot radio.json --kind anode struct --struct '^volTPCActive$' -o tpc.pdf

# drift geometry: anode and cathode planes
bki-plot geo.json --kind anode cathode --orientation auto -o drift.pdf

# field cage only, slim rings (photon-detector window) only
bki-plot geo.json --kind fieldshaper --fieldshaper 'Slim$' -o fc_slim.pdf

# Vertical Drift: drift axis x drawn vertically ("auto" picks it for VD IRs)
bki-plot radio.json --orientation vd -o radio_vd.pdf
bki-plot radio.json --orientation auto

# 2D views stretched to fill their panels instead of drawn to scale
bki-plot radio.json --no-proportional -o radio_stretched.pdf

# a PNG instead of a PDF (format follows the extension)
bki-plot radio.json -o radio.png
```

### Plotting backgrounds

```bash
# radio volumes only
bki-plot radio.json --no-detector -o radio_only.pdf

# selected producers, drawn inside the detector outline
bki-plot radio.json --radio 'GenInCathode|GenInAnode' -o cathode_anode.pdf

# selected producers + only the cryostat for context
bki-plot radio.json --kind radio struct --struct volCryostat \
    --radio '^Ar39' -o ar39.pdf
```

### Plotting individual generation surfaces

```bash
# 1. see which surfaces each background defines
bki-plot radio.json --list --kind radio
#   radio  (5 placements)
#     GenInCathode  (3)
#       #0    plane_x@-300           Rn222            explicit
#       #1    plane_x@0              Rn222            explicit
#       #2    plane_x@300            Rn222            explicit
#     Ar39GenInLAr  (2)
#       #0    bulk                   Ar39             volFC_A
#       #1    bulk                   Ar39             volFC_B

# 2. by index: surfaces 0 and 2 of the cathode producer
#    (all other backgrounds are still drawn in full)
bki-plot radio.json --surface 'GenInCathode#0,2' -o cathode_02.pdf

# index ranges
bki-plot radio.json --surface 'GenInCathode#1-2'

# by regex on surface label / geometry node / nuclide
bki-plot radio.json --surface 'GenInCathode#@plane_x@-'   # surfaces at negative x
bki-plot radio.json --surface 'Ar39GenInLAr#@volFC_B'     # one geometry node

# several specs: each narrows only the producers it names
bki-plot radio.json --surface 'GenInCathode#0' 'Ar39GenInLAr#@bulk'

# only the chosen surfaces, nothing else
bki-plot radio.json --no-detector --radio GenInCathode --surface 'GenInCathode#1'
```

### `--surface` syntax

A background (RadioGen / Decay0Gen producer) can define several generation
surfaces; each becomes one `radio` entry sharing the producer name and
distinguished by `meta.index`.  `--surface` picks among them with specs of
the form `PRODUCER_RE#SEL`:

| `SEL`        | Keeps surfaces whose …                                         |
|--------------|----------------------------------------------------------------|
| `0,2,5-7`    | `meta.index` is in the list                                    |
| `@REGEX`     | `meta.surface` (`bulk`, `plane_x@…`), `meta.node` or `meta.nuclide` matches |

- `meta.surface` is derived by `bki-extract` from the box shape, not read
  from the FHiCL: `plane_<axis>@<centre>` when the box's thinnest extent is
  under 5 % of the next one, `bulk` otherwise. It is computed after the box
  is shrunk to the effective range, so a peaked `distrib_x` source (e.g. a
  mesh sampled in the `volCryostat` box) is labelled a plane.
- Only producers matched by `PRODUCER_RE` (re.search) are narrowed; other
  radio producers are left as they are.  Several specs for the same producer
  are OR-ed.
- Use `--list` to see each producer's surfaces (index, surface label,
  nuclide, geometry node or source).
- With `--surface`, every selected surface gets its own colour and legend
  entry (`producer#index`); with more than 20 surfaces the colouring falls
  back to one colour per producer.
- `meta.node` (the matched geometry node of a `volume_gen` producer) is only
  present in IR files written by a recent `bki-extract`.

### Picking individual placements

Many geometry volumes share a name (112 `volTPC`, 28 `volAnodePlate`, …).
`bki-extract` records where each placement sits in the GDML tree as
`meta.path`, one `volume:k` step per level below the world volume, where `k`
counts the placements of that volume within its parent in GDML order:

```
volDetEnclosure:0/volCryostat:0/volEnclosureTPC:0/volTPC:5/volTPCActive:0
```

The path depends only on the GDML, not on the extraction patterns, so it
identifies the same placement every time the GDML is extracted, and an edit
to one branch of the tree leaves the paths elsewhere unchanged.  `-L` prints
it in the last column — or, for GDMLs that name their physvols (ggd), the
physvol-name path `meta.pvpath` instead, e.g.
`…/placeTPCEnclosure/volTPC_PV_Bot0/volTPCPlaneU_PV`.  Selectors match both.

```bash
# 1. find the placements
bki-plot geo.json -L --kind anode
#   0  anode  volAnodePlate  ...  volDetEnclosure:0/volCryostat:0/volEnclosureTPC:0/volAnodePlate:0
#   1  anode  volAnodePlate  ...  volDetEnclosure:0/volCryostat:0/volEnclosureTPC:0/volAnodePlate:1

# 2. anode plates 0 and 2 only (everything else still drawn)
bki-plot geo.json --pick 'volAnodePlate#@volAnodePlate:[02]$'

# TPC 5 and its active volume
bki-plot geo.json --kind struct --pick '^volTPC$#@volTPC:5$' '^volTPCActive$#@/volTPC:5/'

# arapucas mounted directly in the cryostat (not in the TPC enclosure)
bki-plot geo.json --kind opdet --pick 'volArapuca#@volCryostat:0/volArapuca:'
```

`--pick` takes `NAME_RE#@REGEX` specs, like `--surface`, but for every kind:

- `REGEX` is matched (re.search) against `meta.path` and `meta.part` (the
  primitive solid of a multi-part volume); on radio entries `--pick` works
  exactly like `--surface`, including index lists.
- Index lists (`#0,2`) are for radio entries only; on a geometry volume they
  are an error.
- `NAME_RE` is a `re.search` too: `volTPC` also matches `volTPCActive`, and
  that volume would then be narrowed by the same `REGEX`.  Anchor it
  (`^volTPC$`) to target one volume.
- Volumes not matched by any `NAME_RE` are untouched; several specs for the
  same volume are OR-ed.
- `REGEX` also matches `meta.pv` and `meta.pvpath`, the physvol names, for
  GDMLs that name their physvols (ggd).
- IR files written before `meta.path` existed cannot be picked from; re-run
  `bki-extract`.

### Highlighting volumes

`--pick` and the other selectors decide what is drawn (and the axis
ranges).  `--highlight` then picks out a subset of the drawn volumes: they
are drawn on top in a strong colour, one colour and legend entry per spec,
and everything else is drawn as faint grey context.

```bash
# ggd geometry with the wire planes extracted: where do the bottom-CRM
# readout planes sit?  (they should be at the bottom anode, x ≈ -652)
bki-extract dunevd10kt_v7_full10kt_ggd.gdml -o planes.json --frame volEnclosureTPC \
    --struct-pattern '^volTPCPlane[UVZ]$' '^volTPC$'
bki-plot planes.json --orientation vd \
    --highlight 'bottom planes=^volTPCPlane#@volTPC_PV_Bot' \
                'top planes=^volTPCPlane#@volTPC_PV_Top' 'cathode=^volCathodeGrid$'

# same check on a refactored GDML (unnamed physvols): bottom CRMs are the odd copies
bki-plot planes_ref.json --orientation vd \
    --highlight 'bottom planes=^volTPCPlane#@volTPC:\d*[13579]/'

# wall-mounted photon detectors only, axes fitted to them
bki-plot geo.json --kind opdet cathode --zoom-highlight \
    --highlight 'wall=volArapuca$#@volCryostat:0/volArapuca:'

# check a spec before plotting: highlighted placements are marked *k
bki-plot planes.json -l --highlight 'bottom planes=^volTPCPlane#@volTPC_PV_Bot'
bki-plot planes.json -L --kind struct --highlight '^volTPCPlaneU$#@_PV_Bot'
```

`--highlight` takes `[LABEL=]NAME_RE[#SEL]` specs:

- `NAME_RE` and `#SEL` work as for `--pick`, but `#SEL` is optional: a bare
  `NAME_RE` highlights every placement of the matching volumes.  Unlike
  `--pick`, volumes not named by any spec are not highlighted.
- `LABEL=` (letters, digits, spaces, `-`, `_`, …; no regex characters) names
  the group in the legend; without it the spec itself is used.
- A volume matched by several specs belongs to the first one.
- Highlighted boxes thinner than 0.4 % of an axis span in a 2D view (wire
  planes, X-ARAPUCAs seen edge-on) are drawn as lines so they stay visible.
- With `--highlight`, the per-kind colours and the per-surface colouring of
  `--surface` are replaced by the spec colours and grey.

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
| `meta` | kind-specific metadata: GDML placement `path` for geometry volumes (see [Picking individual placements](#picking-individual-placements)) and, when the GDML names its physvols (ggd), the physvol name `pv` and the physvol names along the path `pvpath` (`placeDetEnclosure/…/volTPC_PV_Bot0/volTPCPlaneU_PV`); nuclide, BqPercc, surface, index, … for radio |

For radio entries, `meta` also holds, for each axis with a `distrib_<axis>`
formula: the formula, `distrib_<axis>_range` (the full box range it is sampled over,
the source of truth), `distrib_<axis>_eff_range` (the effective range: the shortest
interval holding `distrib_eff_containment` of the distribution `TF1::GetRandom`
generates with `distrib_eff_npx` table points, which becomes the box extent on that
axis) and `distrib_<axis>_formula_eff_range` (the same for the formula pdf). The
effective ranges are derived by `bki-extract` with ROOT, not module
parameters. These are world coordinates. The box (`pos`/`ext`) is the full
box on the axes without a formula and the effective range on the others.

`producers` is present for FCL input. It holds the complete FHiCL table of
every selected radiological producer, keyed by module label, and is stored
even when none of the producer's boxes resolve. The `name` of a radio volume
is its producer label. ±inf values are stored as the strings `"inf"` and
`"-inf"`.

`frame` is present only for IRs written with `bki-extract --frame VOLNAME`.
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

---

## Python API

```python
from duneradgeo.ir import IRCollection
from duneradgeo.notebook_plot import make_figure
from duneradgeo.plot import plot_ir

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
from duneradgeo.ir import IRCollection
from duneradgeo.plot import select_volumes
sub = IRCollection(source=ir.source, detector=ir.detector,
                   volumes=select_volumes(ir, surfaces=["GenInCathode#@plane_x"]))
make_figure(sub).show()
```

Position profiles (`bki-distrib`):

```python
from duneradgeo.distrib import collect, profile, with_npx, plot_groups
from duneradgeo.plot import select_volumes

groups = collect(select_volumes(ir, kinds=["radio"]))   # one group per formula + range
pr = profile(groups[0])          # pr.x, pr.formula, pr.generated (pdfs), pr.frac_gen
plot_groups(groups, "profiles.pdf", log=True)           # one page per group

# the formulas directly, evaluated by ROOT
from duneradgeo import _tf1
f = groups[0].formula; lo, hi = groups[0].rng
_tf1.evaluate(f, [-320.0, 0.0])                  # TF1::Eval
_tf1.generated_interval(f, lo, hi, 0.99)         # 99% of what GetRandom generates (Npx=100)
_tf1.quantiles(f, lo, hi, [0.16, 0.5, 0.84])     # TF1::GetQuantiles on the Npx table
_tf1.sample(f, lo, hi, 100000)                   # TF1::GetRandom
```

Parsing utilities:

```python
from duneradgeo import GDMLGeometry, FhiclParser
from duneradgeo._fcl import expand_fcl
```

---

## Detector profiles

Auto-detected from logical volume names (`vd-ggd` is recognised by
`volSupportEnc` / `volEnclosureCryostat` and checked before `vd`):

| Key  | Geometry | Anode | Opdet | Cathode | Field shaper | Default struct |
|------|----------|-------|-------|---------|--------------|----------------|
| `vd` | Vertical Drift (FD2) | `volAnodePlate`, `volAnodePlateBottom` | `volArapuca*` | `volCathodeGrid` (`volGroundGrid` in v4) | `volFieldShaper`, `volFieldShaperSlim` | `volTPCActive`, `volGaseousArgon`, `volTPC`, `volCryostat` |
| `vd-ggd` | Vertical Drift (FD2), ggd-generated GDML; written to the IR as `vd` | as `vd` | as `vd` | `volCathodeGrid` | as `vd` | as `vd` |
| `hd` | Horizontal Drift / ProtoDUNE-HD | `volAPAFrameYSide` | `volArapuca_N` | — | — | `volTPCActiveInner`, `volGaseousArgon`, `volTPC`, `volCryostat` |

Override with `--anode-pattern` / `--opdet-pattern` / `--cathode-pattern` /
`--fieldshaper-pattern` / `--struct-pattern`.

`--with-passive` adds each profile's passive material as `struct` volumes:

| Key  | Passive volumes |
|------|-----------------|
| `vd` | `volSteelShell`, `volFoamPadding`, `volSteelSupport` |
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

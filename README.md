# duneradgeo

Inspect DUNE background sources without a LArSoft environment:

- **detector geometry** — anode (CRP/APA), optical detector, cathode,
  field-shaper, and structural volume positions from GDML files
- **radiological volumes** — emission boxes from RadioGen / Decay0Gen FHiCL
  configurations

Pure Python (NumPy + Matplotlib + Plotly), no PyROOT or LArSoft required.

---

## Installation

```bash
cd dune-radgeo
uv sync                  # or: pip install .   — everything
```

`bki-extract` itself needs only NumPy.  For an extraction-only install,
e.g. on top of a LArSoft Python that already provides NumPy, skip the
dependencies:

```bash
pip install --no-deps .          # add numpy if the environment lacks it
```

| Package | Used by |
|---------|---------|
| `numpy` | everything (`bki-extract`, `bki-plot --list` / `-L` need nothing else) |
| `matplotlib` | `bki-plot` PDF/PNG output, `plot_ir` |
| `plotly` | `make_figure` |
| `jupyter`, `nbconvert` | `nb_B.ipynb`, `scripts/run_notebook.py` |

In a `--no-deps` install, rendering a plot without matplotlib stops with an
`ImportError` that says what to install.

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

### 3. Interactive notebook

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
  --no-geometry       skip GDML parsing (radio volume_rand/gen stay unresolved)

  --anode-pattern P   regex(es) for anode volume names   (auto-detected)
  --opdet-pattern P   regex(es) for opdet volume names   (auto-detected)
  --cathode-pattern P regex(es) for cathode volume names (auto-detected, VD only)
  --fieldshaper-pattern P
                      regex(es) for field-shaper volume names (auto-detected, VD only)
  --struct-pattern P  regex(es) for structural volumes   (auto-detected)
  --struct VOLNAME    exact structural volume names (anchored, overrides --struct-pattern)

  --prefix PREFIX     only producers whose label starts with PREFIX
  --select REGEX      regex filter on producer labels
  --all-configured    include producers not on any trigger path
  --containment F     distrib_x containment fraction (default 0.99)
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
  --list-detail, -L   print every placement with world-frame x/y/z ranges,
                      then exit
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
  --surface S [S …]   radio generation-surface selector (applied last)
```

See [Usage examples](#usage-examples) below.

---

## bki-measure

Geometric measurements on an IR (lengths in cm, world frame).

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

# tighter x-extent for producers with a distrib_x profile (95% containment)
bki-extract prodbackground_radio.fcl -o radio.json --containment 0.95

# an included .fcl was edited: bypass the expansion cache
bki-extract prodbackground_radio.fcl -o radio.json --no-fcl-cache
```

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

---

## JSON IR format

`bki-extract` writes a JSON file with this structure:

```json
{
  "version": "1",
  "source": "path/to/geometry.gdml",
  "detector": "vd",
  "volumes": [
    {
      "kind": "anode",
      "name": "volAnodePlate",
      "pos":  [325.08, -505.10, 148.92],
      "ext":  [0.01, 335.40, 297.84],
      "rot":  [[1,0,0],[0,1,0],[0,0,1]],
      "meta": {}
    }
  ]
}
```

| Field  | Description |
|--------|-------------|
| `kind` | `anode`, `opdet`, `cathode`, `fieldshaper`, `struct`, or `radio` |
| `pos`  | world-frame centre `[x, y, z]` in cm |
| `ext`  | full extents `[dx, dy, dz]` in cm |
| `rot`  | 3×3 world-frame rotation matrix |
| `meta` | kind-specific metadata (nuclide, BqPercc, surface, …) |

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

# same selection as a volume list (e.g. for make_figure)
from duneradgeo.ir import IRCollection
from duneradgeo.plot import select_volumes
sub = IRCollection(source=ir.source, detector=ir.detector,
                   volumes=select_volumes(ir, surfaces=["GenInCathode#@plane_x"]))
make_figure(sub).show()
```

Parsing utilities:

```python
from duneradgeo import GDMLGeometry, FhiclParser
from duneradgeo._fcl import expand_fcl
```

---

## Detector profiles

Auto-detected from logical volume names:

| Key  | Geometry | Anode | Opdet | Cathode | Field shaper | Default struct |
|------|----------|-------|-------|---------|--------------|----------------|
| `vd` | Vertical Drift (FD2) | `volAnodePlate`, `volAnodePlateBottom` | `volArapuca*` | `volCathodeGrid` (`volGroundGrid` in v4) | `volFieldShaper`, `volFieldShaperSlim` | `volTPCActive`, `volGaseousArgon`, `volTPC`, `volCryostat` |
| `hd` | Horizontal Drift / ProtoDUNE-HD | `volAPAFrameYSide` | `volArapuca_N` | — | — | `volTPCActiveInner`, `volGaseousArgon`, `volTPC`, `volCryostat` |

Override with `--anode-pattern` / `--opdet-pattern` / `--cathode-pattern` /
`--fieldshaper-pattern` / `--struct-pattern`.

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

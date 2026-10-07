# duneradgeo

Inspect DUNE background sources without a LArSoft environment:

- **detector geometry** — anode (CRP/APA), optical detector, and structural volume
  positions from GDML files
- **radiological volumes** — emission boxes from RadioGen / Decay0Gen FHiCL
  configurations

Pure Python (NumPy + Matplotlib + Plotly), no PyROOT or LArSoft required.

---

## Installation

```bash
cd duneradgeo
uv sync          # or: pip install -e .
```

Core dependencies: `numpy`, `matplotlib`, `plotly`.
Notebook execution also requires `jupyter` and `nbconvert` (included by default).

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

volume selectors (applied in order):
  --no-detector       exclude anode / opdet / struct volumes
  --no-radio          exclude radio volumes
  --kind K [K …]      keep only these kinds: anode opdet struct radio
  --anode P [P …]     regex filter on anode volume names
  --opdet P [P …]     regex filter on opdet volume names
  --struct P [P …]    regex filter on struct volume names
  --radio P [P …]     regex filter on radio producer names
  --name P [P …]      cross-kind name regex
  --surface S [S …]   radio generation-surface selector (applied last)
```

See [Usage examples](#usage-examples) below.

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
| `kind` | `anode`, `opdet`, `struct`, or `radio` |
| `pos`  | world-frame centre `[x, y, z]` in cm |
| `ext`  | full extents `[dx, dy, dz]` in cm |
| `rot`  | 3×3 world-frame rotation matrix |
| `meta` | kind-specific metadata (nuclide, BqPercc, surface, …) |

Axis convention: **x** = drift, **y** = vertical, **z** = beam.

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

| Key  | Geometry | Anode | Opdet | Default struct |
|------|----------|-------|-------|----------------|
| `vd` | Vertical Drift (FD2) | `volAnodePlate`, `volAnodePlateBottom` | `volArapuca*` | `volTPCActive`, `volGaseousArgon`, `volTPC`, `volCryostat` |
| `hd` | Horizontal Drift / ProtoDUNE-HD | `volAPAFrameYSide` | `volArapuca_N` | `volTPCActiveInner`, `volGaseousArgon`, `volTPC`, `volCryostat` |

Override with `--anode-pattern` / `--opdet-pattern` / `--struct-pattern`.

---

## File search paths

| File type      | Environment variable  |
|----------------|-----------------------|
| FHiCL (`.fcl`) | `$FHICL_FILE_PATH`    |
| GDML geometry  | `$FW_SEARCH_PATH`     |

Both mirror LArSoft's `cet::search_path` behaviour (first-match-wins).

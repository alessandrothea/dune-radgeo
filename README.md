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
via `$FHICL_FILE_PATH`.

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
X–Y side view.

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
  --verbose           extra debug output
```

---

## bki-plot options

```
bki-plot <ir.json> [--out <file.pdf>]

  --list, -l          print volume kinds and names grouped by kind, then exit
  --list-detail, -L   print every placement with world-frame x/y/z ranges,
                      then exit (respects all selectors below)

  --view ELEV AZIM    3D view angles in degrees (default: 20 -60)
  --title TEXT        plot title

volume selectors (applied in order):
  --no-detector       exclude anode / opdet / struct volumes
  --no-radio          exclude radio volumes
  --kind K [K …]      keep only these kinds: anode opdet struct radio
  --anode P [P …]     regex filter on anode volume names
  --opdet P [P …]     regex filter on opdet volume names
  --struct P [P …]    regex filter on struct volume names
  --radio P [P …]     regex filter on radio producer names
  --name P [P …]      cross-kind name regex (applied last)
```

Examples:

```bash
# detector geometry only
bki-plot ir.json --no-radio --out geometry.pdf

# radio volumes only, specific producers
bki-plot ir.json --no-detector --radio 'GenInCathode|GenInAnode' --out cathode_anode.pdf

# list all available volumes
bki-plot ir.json --list

# detailed table with x/y/z ranges for structural volumes only
bki-plot ir.json -L --kind struct
```

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

# plot with selectors
plot_ir(ir, "radio_only.pdf", kinds=["radio"],
        kind_patterns={"radio": ["GenInCathode"]})
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

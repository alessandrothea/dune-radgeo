# radgeo extract

`radgeo extract` reads a GDML geometry or a FHiCL job file and writes the JSON IR
([format](../README.md#json-ir-format)) that the other commands read: the
detector volumes (anode, optical detector, cathode, field shaper, structure) and,
from a FHiCL job, the boxes in which its RadioGen / Decay0Gen producers
generate, with each producer's full FHiCL table.  Further commands:
[radgeo plot](radgeo-plot.md), [radgeo profiles](radgeo-profiles.md),
[radgeo measure](radgeo-measure.md), [radgeo backgrounds](radgeo-backgrounds.md).

## Usage

```bash
radgeo extract geometry.gdml --out detector.json     # detector geometry only
radgeo extract job.fcl --out detector.json           # radio volumes + the job's geometry
radgeo extract job.fcl --gdml /path/to/geometry.gdml --out ir.json
```

Bare filenames are resolved via the LArSoft search paths (FHiCL:
`$FHICL_FILE_PATH`, GDML: `$FW_SEARCH_PATH`, see
[File search paths](../README.md#file-search-paths)).  The command needs
click and rich, plus PyROOT for producers with `distrib_x/y/z`
([Installation](../README.md#installation)).

## Options

```
Usage: radgeo extract [OPTIONS] INPUT

  Usage:
    radgeo extract geometry.gdml --out ir.json
    radgeo extract job.fcl       --out ir.json
    radgeo extract geometry.gdml --out ir.json \
        --anode-pattern 'volAnodePlate' --opdet-pattern '^volArapuca'
    radgeo extract geometry.gdml --out ir.json --frame volEnclosureTPC

  From a FCL file (#includes expanded by fhicl-dump, else via $FHICL_FILE_PATH):
    * RadioGen / Decay0Gen producers → kind="radio" volumes, and each
      producer's full FHiCL table → the IR's top-level "producers";
    * the GDML named in AuxDetGeometry.GDML → the detector volumes
      (anode / opdet / cathode / fieldshaper / struct) and the radio boxes;
    * a distrib_x/y/z formula (needs PyROOT) narrows the box on that axis to
      the effective range of what TF1::GetRandom generates (--npx,
      --containment; see duneradgeo._tf1).

  --frame expresses every volume in the frame of one singly placed volume and
  records its world placement in the IR: the perl and the ggd VD GDMLs nest the
  cryostat differently (a global offset), so their IRs in the frame of
  volEnclosureTPC overlay.  --with-passive adds the profile's passive material
  (steel shell, insulation, support, shielding, cavern) as struct volumes.
  --gdml-backend root reads the GDML with ROOT (TGeoManager::Import), pyg4ometry
  with the pyg4ometry package; scripts/compare_ir.py compares the IRs of two
  backends.

  Non-box solids become boxes: tubes and tori their bounding box, unions one box
  per primitive operand (e.g. each straight tube and corner of a field-shaper
  ring), subtractions and intersections their first operand.  A radio
  volume_rand / volume_gen box is the ROOT bounding box of the volume's solid (a
  union: the whole composite), as BaseRadioGen samples it.

  Output: a JSON IR, read by radgeo plot, radgeo profiles and the nb_B notebook.

Options:
  -o, --out OUT                   output JSON file (IR)  [required]
  --anode-pattern ANODE_PATTERN...
                                  regexes of the anode volume names (default:
                                  the profile's)
  --opdet-pattern OPDET_PATTERN...
                                  regexes of the optical detector volume names
                                  (default: the profile's)
  --cathode-pattern CATHODE_PATTERN...
                                  regexes of the cathode volume names (default:
                                  the profile's)
  --fieldshaper-pattern FIELDSHAPER_PATTERN...
                                  regexes of the field shaper volume names
                                  (default: the profile's)
  --struct-pattern STRUCT_PATTERN...
                                  regexes of the structural volume names
                                  (default: the profile's)
  --struct VOLNAME...             exact structural volume names (anchored,
                                  overrides --struct-pattern)
  --with-passive                  also extract the profile's passive material
                                  (steel shell, insulation, support, shielding,
                                  cavern) as struct
  --frame VOLNAME                 express positions in the frame of this (singly
                                  placed) volume, e.g. volCryostat or
                                  volEnclosureTPC, instead of the world
  --prefix TEXT                   only producers whose label starts with this
  --select TEXT                   regex filter on producer labels
  --all-configured                include producers not on any trigger path
  --containment FLOAT             distrib_x/y/z containment fraction for the
                                  effective range  [default: 0.99]
  --npx INTEGER                   TF1 Npx of the generator: the effective range
                                  is that of the distribution TF1::GetRandom
                                  draws from an Npx-point table; larsimrad
                                  BaseRadioGen keeps ROOT's default  [default:
                                  100]
  --no-geometry                   skip GDML parsing even when a FCL is given
  --gdml GDML                     override GDML path (when input is FCL)
  --gdml-backend [python|root|pyg4ometry]
                                  GDML reader: the pure-Python parser, ROOT's
                                  TGeoManager::Import (needs ROOT with
                                  TGeo/GDML, e.g. in the dunesw container), or
                                  pyg4ometry (optional extra)  [default: python]
  --fcl-cache DIR                 directory where expanded FCL files are cached
                                  and reused across runs  [default:
                                  /home/user/.cache/radgeo/fcl]
  --no-fcl-cache                  always re-run the FCL expansion, ignoring the
                                  cache
  -v, --verbose
  -h, --help                      Show this message and exit.
```

`--containment` and `--npx` define the effective ranges of the producers with
a `distrib_x/y/z` profile, stored as `meta.distrib_eff_containment` and
`meta.distrib_eff_npx`; see
[the formula and the generated distribution](radgeo-profiles.md#the-formula-and-the-generated-distribution-npx).
The detector profile (and with it the default patterns) is auto-detected, see
[Detector profiles](../README.md#detector-profiles).

## FHiCL input

FHiCL files with `#include` directives are expanded automatically:
`fhicl-dump` is used when available; otherwise includes are resolved inline
via `$FHICL_FILE_PATH`. The expanded text is cached in
`~/.cache/radgeo/fcl/` (or `$XDG_CACHE_HOME/radgeo/fcl/`) and
reused on later runs. The cache is keyed on the top-level file's path, size
and mtime plus `$FHICL_FILE_PATH`; edits to *included* files are not
detected, so pass `--no-fcl-cache` after changing them. `--fcl-cache DIR`
moves the cache elsewhere.

Override patterns or supply the GDML explicitly:

```bash
radgeo extract geometry.gdml --out ir.json \
    --anode-pattern '^volAnodePlate$' \
    --opdet-pattern '^volArapuca' \
    --struct volTPCActive volGaseousArgon volTPC volCryostat

radgeo extract job.fcl --gdml /path/to/geometry.gdml --out ir.json
```

## Examples

```bash
# geometry only, detector profile auto-detected from the GDML
radgeo extract geometry.gdml -o geo.json

# full background job: radio producers on the trigger paths + geometry
radgeo extract prodbackground_radio.fcl -o radio.json

# only the producers whose label starts with "GenIn" ...
radgeo extract prodbackground_radio.fcl -o cathode.json --prefix GenIn

# ... or matching a regex, including producers not on any trigger path
radgeo extract prodbackground_radio.fcl -o rn.json --select 'Rn222|Po210' --all-configured

# use a different geometry than the one named in AuxDetGeometry.GDML
radgeo extract prodbackground_radio.fcl -o radio.json --gdml /path/to/geometry.gdml

# narrower effective ranges for producers with a distrib_x/y/z profile (95% containment)
radgeo extract prodbackground_radio.fcl -o radio.json --containment 0.95

# an included .fcl was edited: bypass the expansion cache
radgeo extract prodbackground_radio.fcl -o radio.json --no-fcl-cache
```

## perl vs ggd VD GDML

Both GDML flavours use the same detector volume names (the ggd one is
detected as profile `vd-ggd`, the perl one as `vd-perl`; both are written to
the IR as detector `vd`), and each
geometry volume carries its placement `path` and, for ggd, its physvol name
`pv`.  The ggd geometry nests the cryostat in `volSupportEnc` /
`volEnclosureCryostat`, so world coordinates differ by a global offset:
extract both in a common frame and compare them entry by entry.

```bash
radgeo extract dunevd10kt_3view_30deg_v7_refactored_2x8x40.gdml -o ref.json --frame volEnclosureTPC
radgeo extract dunevd10kt_v7_full10kt_ggd.gdml                 -o ggd.json --frame volEnclosureTPC
python scripts/compare_ir.py ref.json ggd.json -v     # exit status 1 if any group differs

# ggd only: cryostat shell, insulation, I-beams, belts, shielding, cavern
radgeo extract dunevd10kt_v7_full10kt_ggd.gdml -o ggd_passive.json --with-passive
```

## GDML backends

`--gdml-backend` picks how the GDML is read:

| Backend | Reader | Needs |
|---------|--------|-------|
| `python` (default) | built-in XML parser (`GDMLGeometry`) | NumPy |
| `root` | `TGeoManager::Import` (`RootGDMLGeometry`) | ROOT with TGeo/GDML (dunesw) |
| `pyg4ometry` | `pyg4ometry.gdml.Reader` (`Pyg4GDMLGeometry`) | the `pyg4ometry` extra |

Both fill the same tables (volumes and placements, box/tube/torus
decomposition, boolean operands) and share the search, so placement paths, physvol
names, frames and box decomposition are identical; only the parsing
differs.  ROOT evaluates GDML expressions and units itself and computes
each solid's bounding box; the built-in parser builds the same boxes,
booleans included, as ROOT does (see [radio boxes](../README.md#radio-boxes)).  Detector
volumes whose solids are not boxes, tubes, tori or their booleans are
skipped by both; for radio boxes the ROOT backend still uses ROOT's
bounding box, while the built-in parser warns and skips them.
Reading with ROOT takes a few seconds instead of a fraction of one.

`pyg4ometry` is an opt-in extra (`pip install ".[pyg4ometry]"`, or `uv sync
--extra pyg4ometry`), not installed by default: it brings scipy, pandas, sympy
and VTK, about 1 GB.  It evaluates what the built-in parser cannot: `<constant>`
and `<variable>` definitions, expressions and units in any attribute, loops,
scales and `<!ENTITY>` includes.  It takes about 8 s on the VD 1×8×14 geometry
(the built-in parser 0.1 s) and agrees with the built-in parser on the three VD GDMLs in
`examples/gdml/` (`compare_backends.py --backends python pyg4ometry`).  Notes:

- pyg4ometry imports `vtk` although the reader never uses it; when `vtk` is not
  installed a placeholder module is supplied for the import, so a VTK-free
  install works.
- A solid that is not a box, tube, torus or boolean gets the bounding box of its
  pyg4ometry mesh (as in the ROOT backend), so radio boxes work for it; it still
  yields no box decomposition.
- pyg4ometry names an unnamed physvol `<volume>_PV` itself, so the physvol
  names (IR `meta.pv`) are read from the GDML's XML.  If that cannot be parsed
  (`<!ENTITY>` includes), a physvol named exactly `<volume>_PV` is reported as
  unnamed.
- Assemblies, replicas, divisions and parameterised volumes are not followed:
  they are counted and reported on stderr.

`scripts/compare_backends.py` loads one or more GDMLs with both backends
and compares volumes, solid tables and every profile placement (paired by
path), exiting with status 1 on a difference:

```bash
dsw -C <area> py /path/to/dune-radgeo/scripts/compare_backends.py geometry.gdml [more.gdml] [-v]

# any two backends: --backends REF OTHER (default python root)
python /path/to/dune-radgeo/scripts/compare_backends.py geometry.gdml --backends python pyg4ometry

# or through the IR
dsw -C <area> exec radgeo extract geometry.gdml -o /abs/py.json --gdml-backend python
dsw -C <area> exec radgeo extract geometry.gdml -o /abs/rt.json --gdml-backend root
dsw -C <area> exec python /path/to/dune-radgeo/scripts/compare_ir.py /abs/py.json /abs/rt.json
```

The two backends agree exactly on the VD perl GDMLs (refactored 2x8x40,
1x8x14 v5 and v7a), the VD ggd one (with and without wires) and the HD
refactored 1x2x6 one.
[Back to the README](../README.md)

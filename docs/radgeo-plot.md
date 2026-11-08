# radgeo plot

`radgeo plot` renders a JSON IR ([radgeo extract](radgeo-extract.md)) as a matplotlib
PDF (or PNG): a 3D orthographic view and three 2D views of the volumes,
coloured by kind, with selectors to restrict what is drawn, a highlight mode
and listings (`-l` / `-L`) to check a selection before plotting it.

## Usage

```bash
radgeo plot detector.json --out detector.pdf
radgeo plot detector.json --out detector.pdf --view 30 -45
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

## Options

```
Usage: radgeo plot [OPTIONS] IR_FILE

  Render a JSON IR as a matplotlib PDF.

  Produces the same 3D + 3x2D layout as the original plot_gdml_detectors.py:

                   --orientation hd (default)   --orientation vd (x vertical)
    top-left       3D ortho view                3D ortho view
    top-right      Z-Y  beam face               Z-X  beam face
    bot-left       Z-X  top view                Z-Y  top view
    bot-right      X-Y  side view               Y-X  side view

  Usage:
    radgeo plot detector.json --out detector.pdf
    radgeo plot detector.json --out detector.pdf --view 20 -60
    radgeo plot detector.json --highlight 'bottom planes=^volTPCPlane#@volTPC_PV_Bot'

Options:
  -l, --list                      Print volume kinds and names grouped by kind
                                  (radio producers with one row per surface),
                                  then exit (respects selectors).
  -L, --list-detail               Print every placement with its x/y/z ranges in
                                  the IR frame and its path (physvol names when
                                  the GDML has them), then exit (respects
                                  selectors).
  -o, --out TEXT                  Output plot file (default: <stem>.pdf).
  --view ELEV AZIM                3D view elevation and azimuth in degrees.
                                  [default: 20.0, -60.0]
  --title TEXT                    Plot title (default: derived from IR source).
  --proportional / --no-proportional
                                  Draw the 2D views to scale (1 cm == 1 cm on
                                  both axes); --no-proportional stretches each
                                  view to fill its panel (default:
                                  proportional).
  --orientation [hd|vd|auto]      hd: drift axis x horizontal, y vertical; vd: x
                                  vertical in every view; auto: vd for VD IRs,
                                  hd otherwise (default: hd).
  --no-detector                   Exclude detector volumes (anode, opdet,
                                  cathode, fieldshaper, struct).
  --no-radio                      Exclude radiological background volumes
                                  (radio).
  --kind KIND...                  Show only these kinds: anode opdet cathode
                                  fieldshaper struct radio (applied after --no-
                                  detector / --no-radio).
  --name PATTERN...               Show only volumes whose name matches at least
                                  one regex across all kinds (re.search).
  --anode PATTERN...              Regex filter on anode volume names.
  --opdet PATTERN...              Regex filter on opdet volume names.
  --cathode PATTERN...            Regex filter on cathode volume names.
  --fieldshaper PATTERN...        Regex filter on field-shaper volume names.
  --struct PATTERN...             Regex filter on struct volume names.
  --radio PATTERN...              Regex filter on radio producer names.
  --surface SPEC...               Radio generation-surface selector
                                  PRODUCER_RE#SEL; SEL is an index list
                                  (0,2,5-7) or @REGEX matched against surface
                                  label, node and nuclide.  Producers not
                                  matched by any spec are left unfiltered;
                                  selected surfaces are coloured individually.
  --pick SPEC...                  Placement selector NAME_RE#@REGEX for any
                                  kind: keep placements of matching volumes
                                  whose path (shown by -L, e.g.
                                  .../volTPC:37/volTPCPlaneV:0) or part matches
                                  REGEX; for radio entries it works as
                                  --surface.  Volumes not matched by any NAME_RE
                                  are left unfiltered.
  --highlight SPEC...             [LABEL=]NAME_RE[#SEL]: draw the matching
                                  volumes on top, one colour and legend entry
                                  per spec, and the rest as grey context.  SEL
                                  as for --pick (optional; @REGEX also matches
                                  the physvol name pv and pvpath).  Applied
                                  after the selectors.
  --dim F                         Face opacity of the context volumes (with
                                  --highlight).  [default: 0.02]
  --zoom-highlight                Fit the axes to the highlighted volumes.
  -h, --help                      Show this message and exit.
```

The selectors (`--kind` to `--pick`) are applied in the order of the
listing; `--highlight` is applied after them.  The sections below show them in
use.  From Python, the same selection is available through `plot_ir` and
`select_volumes` ([Python API](../README.md#python-api)).

## Inspecting an IR

```bash
# kinds, names and placement counts; radio producers list each surface
radgeo plot radio.json --list

# per-placement IR-frame ranges, restricted to structural volumes
radgeo plot radio.json -L --kind struct

# where exactly are the cathode generation surfaces?
radgeo plot radio.json -L --kind radio --radio GenInCathode
```

Every selector of this page also works with `--list` / `-L`, so a selection can be
checked before plotting it.

## Plotting the detector

```bash
# everything, output next to the IR (radio.pdf)
radgeo plot radio.json

# detector geometry only, custom title and view angle
radgeo plot radio.json --no-radio --title "FD2 VD geometry" --view 30 -45 -o geometry.pdf

# anodes and active volume only
radgeo plot radio.json --kind anode struct --struct '^volTPCActive$' -o tpc.pdf

# drift geometry: anode and cathode planes
radgeo plot geo.json --kind anode cathode --orientation auto -o drift.pdf

# field cage only, slim rings (photon-detector window) only
radgeo plot geo.json --kind fieldshaper --fieldshaper 'Slim$' -o fc_slim.pdf

# Vertical Drift: drift axis x drawn vertically ("auto" picks it for VD IRs)
radgeo plot radio.json --orientation vd -o radio_vd.pdf
radgeo plot radio.json --orientation auto

# 2D views stretched to fill their panels instead of drawn to scale
radgeo plot radio.json --no-proportional -o radio_stretched.pdf

# a PNG instead of a PDF (format follows the extension)
radgeo plot radio.json -o radio.png
```

## Plotting backgrounds

```bash
# radio volumes only
radgeo plot radio.json --no-detector -o radio_only.pdf

# selected producers, drawn inside the detector outline
radgeo plot radio.json --radio 'GenInCathode|GenInAnode' -o cathode_anode.pdf

# selected producers + only the cryostat for context
radgeo plot radio.json --kind radio struct --struct volCryostat \
    --radio '^Ar39' -o ar39.pdf
```

## Plotting individual generation surfaces

```bash
# 1. see which surfaces each background defines
radgeo plot radio.json --list --kind radio
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
radgeo plot radio.json --surface 'GenInCathode#0,2' -o cathode_02.pdf

# index ranges
radgeo plot radio.json --surface 'GenInCathode#1-2'

# by regex on surface label / geometry node / nuclide
radgeo plot radio.json --surface 'GenInCathode#@plane_x@-'   # surfaces at negative x
radgeo plot radio.json --surface 'Ar39GenInLAr#@volFC_B'     # one geometry node

# several specs: each narrows only the producers it names
radgeo plot radio.json --surface 'GenInCathode#0' 'Ar39GenInLAr#@bulk'

# only the chosen surfaces, nothing else
radgeo plot radio.json --no-detector --radio GenInCathode --surface 'GenInCathode#1'
```

## `--surface` syntax

A background (RadioGen / Decay0Gen producer) can define several generation
surfaces; each becomes one `radio` entry sharing the producer name and
distinguished by `meta.index`.  `--surface` picks among them with specs of
the form `PRODUCER_RE#SEL`:

| `SEL`        | Keeps surfaces whose …                                         |
|--------------|----------------------------------------------------------------|
| `0,2,5-7`    | `meta.index` is in the list                                    |
| `@REGEX`     | `meta.surface` (`bulk`, `plane_x@…`), `meta.node` or `meta.nuclide` matches |

- `meta.surface` is derived by `radgeo extract` from the box shape, not read
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

## Picking individual placements

Many geometry volumes share a name (112 `volTPC`, 28 `volAnodePlate`, …).
`radgeo extract` records where each placement sits in the GDML tree as
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
radgeo plot geo.json -L --kind anode
#   0  anode  volAnodePlate  ...  volDetEnclosure:0/volCryostat:0/volEnclosureTPC:0/volAnodePlate:0
#   1  anode  volAnodePlate  ...  volDetEnclosure:0/volCryostat:0/volEnclosureTPC:0/volAnodePlate:1

# 2. anode plates 0 and 2 only (everything else still drawn)
radgeo plot geo.json --pick 'volAnodePlate#@volAnodePlate:[02]$'

# TPC 5 and its active volume
radgeo plot geo.json --kind struct --pick '^volTPC$#@volTPC:5$' '^volTPCActive$#@/volTPC:5/'

# arapucas mounted directly in the cryostat (not in the TPC enclosure)
radgeo plot geo.json --kind opdet --pick 'volArapuca#@volCryostat:0/volArapuca:'
```

`--pick` takes `NAME_RE#@REGEX` specs, like `--surface`, but for every kind:

- `REGEX` is matched (re.search) against `meta.path`, `meta.part` (the
  primitive solid of a multi-part volume) and, for GDMLs that name their
  physvols (ggd), `meta.pv` and `meta.pvpath`; on radio entries `--pick`
  works exactly like `--surface`, including index lists.
- Index lists (`#0,2`) are for radio entries only; on a geometry volume they
  are an error.
- `NAME_RE` is a `re.search` too: `volTPC` also matches `volTPCActive`, and
  that volume would then be narrowed by the same `REGEX`.  Anchor it
  (`^volTPC$`) to target one volume.
- Volumes not matched by any `NAME_RE` are untouched; several specs for the
  same volume are OR-ed.
- IR files written before `meta.path` existed cannot be picked from; re-run
  `radgeo extract`.

## Highlighting volumes

`--pick` and the other selectors decide what is drawn (and the axis
ranges).  `--highlight` then picks out a subset of the drawn volumes: they
are drawn on top in a strong colour, one colour and legend entry per spec,
and everything else is drawn as faint grey context.

```bash
# ggd geometry with the wire planes extracted: where do the bottom-CRM
# readout planes sit?  (they should be at the bottom anode, x ≈ -652)
radgeo extract dunevd10kt_v7_full10kt_ggd.gdml -o planes.json --frame volEnclosureTPC \
    --struct-pattern '^volTPCPlane[UVZ]$' '^volTPC$'
radgeo plot planes.json --orientation vd \
    --highlight 'bottom planes=^volTPCPlane#@volTPC_PV_Bot' \
                'top planes=^volTPCPlane#@volTPC_PV_Top' 'cathode=^volCathodeGrid$'

# same check on a perl GDML (unnamed physvols): bottom CRMs are the odd copies
radgeo plot planes_ref.json --orientation vd \
    --highlight 'bottom planes=^volTPCPlane#@volTPC:\d*[13579]/'

# wall-mounted photon detectors only, axes fitted to them
radgeo plot geo.json --kind opdet cathode --zoom-highlight \
    --highlight 'wall=volArapuca$#@volCryostat:0/volArapuca:'

# check a spec before plotting: highlighted placements are marked *k
radgeo plot planes.json -l --highlight 'bottom planes=^volTPCPlane#@volTPC_PV_Bot'
radgeo plot planes.json -L --kind struct --highlight '^volTPCPlaneU$#@_PV_Bot'
```

`--highlight` takes `[LABEL=]NAME_RE[#SEL]` specs:

- `NAME_RE` and `#SEL` work as for `--pick`, but `#SEL` is optional: a bare
  `NAME_RE` highlights every placement of the matching volumes.  Unlike
  `--pick`, volumes not named by any spec are not highlighted.
- `LABEL=` (letters, digits, spaces, `-`, `_`, …; no regex characters) names
  the group in the legend; without it the spec itself is used.  Spaces
  around the `=` are ignored (`grp = ^K42` is `grp=^K42`).
- A volume matched by several specs belongs to the first one.
- Highlighted boxes thinner than 0.4 % of an axis span in a 2D view (wire
  planes, X-ARAPUCAs seen edge-on) are drawn as lines so they stay visible.
- With `--highlight`, the per-kind colours and the per-surface colouring of
  `--surface` are replaced by the spec colours and grey.
[Back to the README](../README.md)

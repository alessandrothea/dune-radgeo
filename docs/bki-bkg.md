# bki-bkg

`bki-bkg` inspects the radiological backgrounds of an IR. A background is one
radiological producer: the radio volumes of the IR that share its name (its
sampling boxes) and its FHiCL parameters (`producers`). Tables are drawn with
[rich](https://rich.readthedocs.io); `--json`, given after the command (or
before the IR file), prints the same data as JSON instead, for scripts (and
needs no `rich`). The command line uses [click](https://click.palletsprojects.com),
so `bki-bkg IR_FILE COMMAND --help` describes each command (the help texts are
reproduced below).  The IR comes from [bki-extract](bki-extract.md).

## Usage

```bash
bki-bkg radbkg.json list
bki-bkg radbkg.json select nuclide=K40 'BqPercc<0.1'      # both must hold
bki-bkg radbkg.json select 'material=G10|vm2000' --json
bki-bkg radbkg.json region --x -330 -300 --z 0 500 --volumes
bki-bkg radbkg.json show Mesh --limit 0
```

## Options

```
Usage: bki-bkg [OPTIONS] IR_FILE COMMAND [ARGS]...

  Inspect the radiological backgrounds of a JSON IR (bki-extract).

  A background is one radiological producer: its radio volumes (sampling boxes)
  and its FHiCL parameters.

  Examples:
    bki-bkg ir.json list
    bki-bkg ir.json show Po218 Cathode
    bki-bkg ir.json select nuclide=Po218 'BqPercc<0.01'
    bki-bkg ir.json region --slice x=-328
    bki-bkg ir.json region --x -330 -300 --z 0 500 --volumes
    bki-bkg ir.json select material=LAr --json

Options:
  --json      Print JSON instead of tables (also accepted after the command;
              needs no rich).
  -h, --help  Show this message and exit.

Commands:
  list     Every background: module, nuclide, activity, material, sampling,...
  region   The backgrounds active in a region.
  select   The backgrounds for which every EXPR holds.
  show     The FHiCL parameters, surfaces and volumes of the backgrounds...
  suggest  Selections and regions that split the backgrounds, with the...
```

### `list`

```
Usage: bki-bkg IR_FILE list [OPTIONS]

  Every background: module, nuclide, activity, material, sampling, volumes.

Options:
  --sort [name|volumes]  Order of the rows (volumes: most first).  [default:
                         name]
  --json                 Print JSON instead of tables.
  -h, --help             Show this message and exit.
```

### `show`

```
Usage: bki-bkg IR_FILE show [OPTIONS] NAME_RE...

  The FHiCL parameters, surfaces and volumes of the backgrounds whose name
  matches a regex (re.search).

Options:
  --limit INTEGER  Volumes shown per background (0: all).  [default: 25]
  --json           Print JSON instead of tables.
  -h, --help       Show this message and exit.
```

### `select`

```
Usage: bki-bkg IR_FILE select [OPTIONS] EXPR...

  The backgrounds for which every EXPR holds.

  EXPR is
    FIELD=REGEX    re.search on any value of the field
    FIELD!=REGEX   no value matches (also true if the field is absent)
    FIELD<N        also <=, >, >=: any numeric value
    FIELD          the background has the field

  A FIELD is a FHiCL parameter of the producer or a metadata key of its volumes
  (nuclide, surface, node, source, module_type, BqPercc, T0, ...), plus name and
  nvolumes; `suggest` lists them.

Options:
  --json      Print JSON instead of tables.
  -h, --help  Show this message and exit.
```

### `region`

```
Usage: bki-bkg IR_FILE region [OPTIONS]

  The backgrounds active in a region.

  The region is the intersection of the constraints, in the frame stored in the
  IR; a background is active when the box of at least one of its volumes
  overlaps it (touching counts).

Options:
  --x LO HI
  --y LO HI
  --z LO HI
  --box X0 X1 Y0 Y1 Z0 Z1  A 3D box.
  --slice AXIS=VALUE       The plane AXIS = VALUE, e.g. x=-328; repeat for a
                           line or a point.
  --thickness T            Thickness of the --slice planes.  [default: 0.0]
  --volumes                Also list the overlapping volumes.
  --json                   Print JSON instead of tables.
  -h, --help               Show this message and exit.
```

### `suggest`

```
Usage: bki-bkg IR_FILE suggest [OPTIONS]

  Selections and regions that split the backgrounds, with the commands to run
  them.

Options:
  --json      Print JSON instead of tables.
  -h, --help  Show this message and exit.
```

## What the commands do

**`select`.** An `EXPR` is `FIELD=REGEX` (`re.search` on any value of the
field), `FIELD!=REGEX` (no value matches; also true when the background lacks
the field), `FIELD<N`, `<=`, `>`, `>=` (any numeric value), or a bare `FIELD`
(the background has it). A `FIELD` is a FHiCL parameter of the producer
(`isotope`, `material`, `volume_gen`, `decay_chain`, `BqPercc`, ...) or a
metadata key of its volumes (`nuclide`, `surface`, `node`, `source`,
`module_type`, `T0`, ...), plus `name` and `nvolumes`. A list or a table
contributes its items. An unknown field is an error that names the closest
known ones.

**`region`.** The region is the intersection of the constraints, in the frame
stored in the IR (world unless extracted with `--frame`). A background is
active when the axis-aligned box of at least one of its volumes overlaps the
region; touching counts. `--slice x=-328` is the plane `x = -328` (a
`--thickness` gives it a width). The table gives, per background, how many of
its volumes overlap and the extent of those; `--volumes` lists them.

**`suggest`.** The values shared by several backgrounds for each field, with
the `select` command and the number of backgrounds; thresholds for the numeric
fields; the words of the names (`Mesh`, `Cathode`, ...); and the planes the
volumes lie on, with the `region --slice` command, how many backgrounds lie on
the plane and how many merely cross it (bulk backgrounds cross every plane in
the cryostat).

[Back to the README](../README.md)

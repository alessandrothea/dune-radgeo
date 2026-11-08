# radgeo measure

Geometric measurements on an IR (lengths in cm, in the IR frame -- the world
unless the IR was extracted with `--frame`; the measured distances do not
depend on it).  The IR comes from [radgeo extract](radgeo-extract.md).

## Usage

```bash
radgeo measure geo.json drift-gaps        # (without a name: every measurement, one panel each)
# ╭─ drift-gaps: cathode ↔ anode distance along x (drift)  [cm] ──────────────────────────────╮
# │ top     gap   651.080  (min 651.080, max 651.080, 28 anodes)   centre-to-centre   653.085 │
# │         volTPCActive x-extent   650.000  (112 volumes)                                    │
# │ bottom  no anodes                                                                         │
# ╰───────────────────────────────────────────────────────────────────────────────────────────╯

radgeo measure geo.json drift-gaps -d --json gaps.json   # per anode–cathode pair
```

## Options

```
Usage: radgeo measure [OPTIONS] [IR_FILE] [MEASUREMENT]...

  Geometric measurements on a JSON IR.

  Usage:
    radgeo measure detector.json                    # every measurement, summaries
    radgeo measure detector.json drift-gaps -d      # one measurement, with per-row table
    radgeo measure detector.json --json out.json    # also write the results as JSON
    radgeo measure --list                           # available measurements

  Lengths are in cm, in the IR frame (x = drift, z = beam).

  The measurements live in the duneradgeo.measurements package, one module each;
  adding a module there adds a measurement (see its docstring).

Options:
  --list        List the available measurements and exit.
  -d, --detail  Also print each measurement's per-row table.
  --json FILE   Write the results (summary, values, rows) as JSON.
  -h, --help    Show this message and exit.
```

| Measurement  | What it measures | Needs |
|--------------|------------------|-------|
| `anode-cryostat-gap` | along x, each anode ↔ the `volCryostat` face behind it (top / bottom, split at the cathode plane); `gap_inner` / `gap_outer` from the anode's drift-side / cryostat-side face | `anode`, `struct volCryostat` |
| `drift-gaps` | along x, each anode ↔ the cathode module it faces: `gap` between facing surfaces, `centre_dist`; per side, with outliers flagged | `anode`, `cathode` |
| `opdet-fc-gap` | each wall optical detector ↔ the closest field-shaper part (box-to-box distance), per wall | `opdet`, `fieldshaper` |
| `fc-cryostat-gap` | the four field-cage side walls ↔ the `volCryostat` walls: `gap_min` (closest piece), `gap_max` (most recessed piece) | `fieldshaper`, `struct volCryostat` |

The full definitions are the docstrings of the measurement functions in
`src/duneradgeo/measurements/`.

## Adding a measurement

Measurements live in the `duneradgeo.measurements` package, one module each:

```
src/duneradgeo/
  measure.py                  radgeo measure CLI (output, --detail, --json)
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
        lo, hi = aabb(v)                        # IR-frame box corners
        res.rows.append(dict(cathode=tag, dx=float(hi[X] - lo[X])))
    dx = [r["dx"] for r in res.rows]
    res.values["dx"] = max(dx) if dx else None  # headline numbers (→ JSON)
    res.summary.append(f"dx {res.values['dx']}")  # lines printed by default
    return res
```

`Result` fields: `summary` (printed lines), `values` (headline numbers,
`None` when not applicable), `columns` + `rows` (the `--detail` table).
Helpers: `aabb(v)` the IR-frame (lo, hi) corners of a volume's box (via
`duneradgeo._draw.aabb`), `tags(vols)` stable `name#k` labels,
`overlap_area(...)` projected overlap of two boxes, `box_distances(...)`
shortest distance from one box to many, `struct_box(ir, pattern)` bounding
box of matching struct volumes (e.g. the cryostat).  From Python,
`duneradgeo.measurements.run(ir, ["drift-gaps"])` returns the `Result` list.

[Back to the README](../README.md)

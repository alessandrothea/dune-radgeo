# radgeo profiles

Draws the position profiles of the radiological producers that set
`distrib_x`, `distrib_y` or `distrib_z`. In larsimrad's `BaseRadioGen`, each
of these is a 1D TF1 formula in the variable `x`. The generator samples it
with `TF1::GetRandom` over the box range on that axis: `X0..X1` etc., the
`volume_rand` bounding box, or each `volume_gen` node's box (see
[radio boxes](../README.md#radio-boxes)). Axes without a formula are sampled uniformly.
The IR comes from [radgeo extract](radgeo-extract.md).

## Usage

```bash
radgeo profiles radbkg.json --log
radgeo profiles radbkg.json --list --radio 'UpperMesh'
# ╭─ [0] UpperMesh1x8x14 · distrib_x ─────────────────────────────────────────
# │ range [-425.120, 425.200]  containment 0.99  npx 100
# │ eff_range (generated) [-323.082, -318.753]  (measured 0.9900)
# │ formula eff range     [-321.704, -318.328]  (measured 0.9900)
# │ surface: plane_x@-320.917
# │ formula: (7.2257e+4*TMath::Exp(-(x+320.)*(x+320.)/(2.*0.5*0.5)))+1.
# │
# │  producer                                  nuclide             placements
# │  K42From42ArGenInUpperMesh1x8x14           K42                          1
# │  ...   (panel abbreviated)
```

## Options

```
Usage: radgeo profiles [OPTIONS] IR_FILE

  radgeo profiles  —  draw the distrib_x/y/z position profiles stored in a JSON
  IR

  RadioGen/Decay0Gen producers with a distrib_<axis> parameter draw that
  coordinate from a TF1 formula in x over the box range (TF1::GetRandom,
  larsimrad BaseRadioGen).  Each page shows the formula (TF1::Eval) and the
  distribution as generated, with the effective range of each marked (the
  generated one, shaded, is the box radgeo extract writes). GetRandom samples
  from an Npx-point table (default 100), so features narrower than range/Npx are
  lost; see the duneradgeo._tf1 docstring. Producers sharing a formula and range
  share one plot: a page of the PDF, or one file per plot for other formats.
  Needs PyROOT.  Not modelled: the redraw of points outside the volume_gen /
  material volumes.  Coordinates are world coordinates, as the generator samples
  them, also for IRs written with radgeo extract --frame.

  Usage:
    radgeo profiles detector.json                     # -> detector_distrib.pdf
    radgeo profiles detector.json --out distrib.png --log   # distrib_0.png, …
    radgeo profiles detector.json --list              # text only, no matplotlib
    radgeo profiles detector.json --radio 'Mesh|Cathode'

  Reads the distrib_* ranges, containment and Npx that radgeo extract writes.

Options:
  -o, --out FILE      Output file (default: <stem>_distrib.pdf): a PDF gets one
                      page per distribution, other formats one file per
                      distribution, <stem>_NN<suffix>.
  -l, --list          Print the distributions, effective ranges and producers,
                      then exit (needs ROOT, not matplotlib).
  --radio PATTERN...  Regex filter on radio producer names.
  --pick SPEC...      PRODUCER_RE#SEL as in radgeo plot --surface (index list or
                      @REGEX on surface, node, nuclide).
  --axis {x,y,z}...   Only these axes (default: all present).
  --npx N             TF1 Npx for the generated distribution (default: the
                      distrib_eff_npx stored in the IR); another value
                      recomputes the generated effective range.
  --log               Logarithmic y axis.
  --no-zoom           No zoomed view on the effective range (by default it is
                      added below the plot when the effective range is narrower
                      than half the full range).
  --title TEXT        FHiCL shown in the LArSoft configuration section (default:
                      the source FHiCL name).
  -h, --help          Show this message and exit.
```

`--pick` takes the same specs as [`radgeo plot --surface`](radgeo-plot.md#--surface-syntax).

## The formula and the generated distribution (Npx)

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

`radgeo extract` records, per axis with a formula, the full range
(`distrib_<axis>_range`), the shortest interval holding `--containment` of the
generated distribution (`distrib_<axis>_eff_range`, the effective range, which
becomes the box on that axis) and the same for the formula
(`distrib_<axis>_formula_eff_range`), plus `distrib_eff_containment` and `distrib_eff_npx`.

Not modelled: the generator also redraws points that fall outside the
volumes allowed by `volume_gen` / `material`. With `material: LAr` in the
`volCryostat` box, nothing is generated in the gas above the liquid (x >
325.09 cm in the 1x8x14 VD geometry), and the profile is renormalised over
the liquid.

## Pages

`radgeo profiles` draws, over the full range, the formula (black) and the
generated distribution (blue), shades the generated effective range and
marks the formula's (dotted), one distribution per page of the PDF. When the
effective ranges are much narrower than the full range, a zoomed view on
them is added below. Producers that share a formula and range share a page,
which lists them. Each page is headed by the producers' common name and
their nuclides (`UpperMesh1x8x14 · x profile, 7 producers (K42, Po218, …)`;
a single producer gives its own name). The side panel has a Location
section (shape: a thin layer at x = c or the surface label; the range
holding the containment fraction; the sampling volume), a Sampling section,
the producers, and a LArSoft configuration section naming the source FHiCL;
parameter names are in bold. `--list` prints the same name on each
`[n]` line. Coordinates are world coordinates, which is what the generator
uses, even for IRs written with `--frame`.

Each page has a formula box with two captioned monospace blocks: the FHiCL
parameter exactly as ROOT evaluates it (`distrib_x: "…"`, wrapped between
top-level terms, nothing added or dropped), and the same formula compacted
for reading only (`f(x) = …`, one term per line, `TMath::` dropped, `Abs(a)`
as `|a|`, `Power(a,n)` as `a^n`, `*` as `·`). A side panel has the full and
both effective ranges, the containment (stored, and measured for each
distribution), `Npx`, the generator volume (`meta.source`) and the producers
with their nuclides.

IR files must carry the `distrib_*` keys above: re-run `radgeo extract` for older ones.
[Back to the README](../README.md)

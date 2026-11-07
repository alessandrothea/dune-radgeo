"""Colours, palettes and the matplotlib import shared by the plotting modules.

Importing this module does not import matplotlib; :func:`import_mpl` does.
"""
from __future__ import annotations

# ---------------------------------------------------------------------------
# Per-kind colours (bki-plot; notebook_plot uses them too)
# ---------------------------------------------------------------------------

ANODE_COL   = "royalblue"
OPDET_COL   = "darkorange"
CATHODE_COL = "teal"
FSHAPER_COL = "slategray"
RADIO_COL   = "crimson"
FACE_ALPHA  = 0.25

# struct volumes are coloured per name (STRUCT_PALETTE), not per kind
KIND_COLOUR: dict[str, str] = {
    "anode":       ANODE_COL,
    "opdet":       OPDET_COL,
    "cathode":     CATHODE_COL,
    "fieldshaper": FSHAPER_COL,
    "radio":       RADIO_COL,
}

# struct volume name → (colour, face alpha)
STRUCT_PALETTE: dict[str, tuple] = {
    "volTPCActive":      ("#4CAF50", 0.08),   # green
    "volTPCActiveInner": ("#4CAF50", 0.08),
    "volGaseousArgon":   ("#90CAF9", 0.18),   # light blue
    "volTPC":            ("#FFF176", 0.06),   # pale yellow
    "volCryostat":       ("#CE93D8", 0.05),   # pale purple
}
STRUCT_DEFAULT = ("#BDBDBD", 0.08)

# Per-surface colours used when individual radio surfaces are selected;
# above this many groups the colouring collapses to one colour per producer.
SURFACE_PALETTE = ("tab:red", "tab:blue", "tab:green", "tab:purple",
                   "tab:brown", "tab:pink", "tab:olive", "tab:cyan",
                   "tab:orange", "tab:gray")
MAX_SURFACE_GROUPS = 20

# --highlight: highlighted groups take SURFACE_PALETTE colours, everything
# else is drawn as grey context with face opacity `dim`.
CONTEXT_COL   = "#9E9E9E"
# low: faces of overlapping context boxes stack up along the view direction
DIM_DEFAULT   = 0.02
HL_FACE_ALPHA = 0.45
HL_LINEWIDTH  = 1.2
# a highlighted box whose 2D projection is thinner than this fraction of the
# axis span is drawn as a line, so e.g. 0.02 cm wire planes stay visible
THIN_FRAC     = 0.004


# ---------------------------------------------------------------------------
# matplotlib
# ---------------------------------------------------------------------------

def import_mpl(note: str = ""):
    """Import matplotlib with the non-interactive Agg backend; return pyplot.

    Raises ImportError with an install hint (plus `note`, e.g. which
    options work without matplotlib) when matplotlib is missing."""
    try:
        import matplotlib
    except ImportError as exc:
        raise ImportError("plotting needs matplotlib (pip install matplotlib)"
                          + (f"; {note}" if note else "")) from exc
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt

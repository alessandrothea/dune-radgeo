"""FHiCL parser and radiological-volume extractor (no PyROOT required)."""
from __future__ import annotations

import copy
import math
import os
import re
import sys

import numpy as np

# ---------------------------------------------------------------------------
# #include preprocessor
# ---------------------------------------------------------------------------

_INCLUDE_RE = re.compile(r'^\s*#include\s+"([^"]+)"', re.MULTILINE)


def _has_includes(text: str) -> bool:
    return bool(_INCLUDE_RE.search(text))


def default_fcl_cache_dir() -> str:
    """``$XDG_CACHE_HOME/bki-extract/fcl`` (default ``~/.cache/bki-extract/fcl``)."""
    base = os.environ.get("XDG_CACHE_HOME") or os.path.expanduser("~/.cache")
    return os.path.join(base, "bki-extract", "fcl")


def _cache_file(path: str, search_var: str, cache_dir: str) -> str:
    """
    Cache file name for the expanded form of `path`.

    The key covers the resolved path, its size/mtime and the include search
    path, so editing or moving the top-level file or changing
    $FHICL_FILE_PATH invalidates the entry.  Changes to *included* files are
    not detected: use ``--no-fcl-cache`` (or delete the entry) in that case.
    """
    import hashlib
    st  = os.stat(path)
    key = "\0".join([
        os.path.realpath(path),
        str(st.st_size),
        str(st.st_mtime_ns),
        os.environ.get(search_var, ""),
    ])
    h = hashlib.sha1(key.encode()).hexdigest()[:16]
    stem = os.path.splitext(os.path.basename(path))[0]
    return os.path.join(cache_dir, f"{stem}.{h}.expanded.fcl")


def expand_fcl(path: str, search_var: str = "FHICL_FILE_PATH",
               cache_dir: str | None = None) -> str:
    """
    Return the fully expanded text of a FHiCL file.

    Strategy (in order):
      1. Run ``fhicl-dump <path>`` if it is available on $PATH — this is the
         authoritative LArSoft expander and handles all FHiCL edge cases.
      2. Fall back to our inline ``#include`` resolver (pure Python, uses
         $FHICL_FILE_PATH) when fhicl-dump is not found.

    Files without any ``#include`` are returned as-is without spawning a
    subprocess.

    When `cache_dir` is given, the expanded text is stored there (see
    :func:`_cache_file`) and reused on later runs so the expensive
    expansion step is skipped.
    """
    with open(path) as fh:
        text = fh.read()

    if not _has_includes(text):
        return text   # nothing to expand

    bname = os.path.basename(path)

    # -- cached expansion from a previous run --------------------------------
    cpath = _cache_file(path, search_var, cache_dir) if cache_dir else None
    if cpath and os.path.isfile(cpath):
        print(f"[info] using cached expansion of '{bname}': {cpath}", file=sys.stderr)
        with open(cpath) as fh:
            return fh.read()

    expanded = _expand_fcl_uncached(path, text, bname, search_var)

    if cpath:
        try:
            os.makedirs(cache_dir, exist_ok=True)
            tmp = cpath + ".tmp"
            with open(tmp, "w") as fh:
                fh.write(expanded)
            os.replace(tmp, cpath)
            print(f"[info] cached expansion → {cpath}", file=sys.stderr)
        except OSError as exc:
            print(f"[warn] could not write FCL cache ({exc})", file=sys.stderr)
    return expanded


def _expand_fcl_uncached(path: str, text: str, bname: str, search_var: str) -> str:
    # -- try fhicl-dump first ------------------------------------------------
    import shutil
    import subprocess
    if shutil.which("fhicl-dump"):
        print(f"[info] expanding '{bname}' with fhicl-dump …", file=sys.stderr)
        try:
            result = subprocess.run(
                ["fhicl-dump", path],
                capture_output=True, text=True, check=True,
            )
            print(f"[info] expansion done", file=sys.stderr)
            return result.stdout
        except subprocess.CalledProcessError as exc:
            print(f"[warn] fhicl-dump failed ({exc}); falling back to inline resolver",
                  file=sys.stderr)

    # -- inline fallback -----------------------------------------------------
    print(f"[info] expanding '{bname}' via ${search_var} …", file=sys.stderr)
    expanded = _preprocess_fcl(text, search_var=search_var)
    print(f"[info] expansion done", file=sys.stderr)
    return expanded


def _preprocess_fcl(text: str,
                    search_var: str = "FHICL_FILE_PATH",
                    _seen: frozenset = frozenset()) -> str:
    """
    Recursively inline #include "filename.fcl" directives via $FHICL_FILE_PATH.
    Used as a fallback when fhicl-dump is not available.
    """
    out = []
    for line in text.splitlines():
        m = _INCLUDE_RE.match(line)
        if not m:
            out.append(line)
            continue
        fname = m.group(1)
        found = _resolve_include(fname, search_var)
        if found is None:
            raise FileNotFoundError(
                f"#include \"{fname}\" not found in ${search_var}")
        key = os.path.realpath(found)
        if key in _seen:
            continue   # circular include guard
        print(f"[info] include '{fname}' → {found}", file=sys.stderr)
        with open(found) as fh:
            sub = _preprocess_fcl(fh.read(), search_var, _seen | {key})
        out.append(sub)
    return "\n".join(out)


def _resolve_include(fname: str, var: str) -> str | None:
    """First <dir>/<fname> that exists across the colon-separated $PATH."""
    if os.path.isabs(fname) and os.path.isfile(fname):
        return fname
    for d in (d for d in os.environ.get(var, "").split(":") if d):
        cand = os.path.join(d, fname)
        if os.path.isfile(cand):
            return cand
    return None

# ---------------------------------------------------------------------------
# Minimal FHiCL tokenizer / parser
# ---------------------------------------------------------------------------

ERASE = object()

_TOKEN_RE = re.compile(r"""
    (?P<ws>\s+)
  | (?P<comment>(\#|//)[^\n]*)
  | (?P<string>"(?:\\.|[^"\\])*"|'[^']*')
  | (?P<ref>@(?:local|table|sequence|id)::[A-Za-z_][\w.\[\]]*)
  | (?P<protect>@protect_ignore:|@protect_error:)
  | (?P<erase>@erase)
  | (?P<nil>@nil)
  | (?P<number>[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?(?![\w.]))
  | (?P<ident>[A-Za-z_][\w]*(?:\.[A-Za-z_]\w*|\[\d+\])*)
  | (?P<punct>[{}\[\]:,])
""", re.VERBOSE)


def _tokenize(text: str) -> list:
    pos, out = 0, []
    while pos < len(text):
        m = _TOKEN_RE.match(text, pos)
        if not m:
            line = text.count("\n", 0, pos) + 1
            raise SyntaxError(f"line {line}: cannot tokenize near {text[pos:pos+40]!r}")
        pos = m.end()
        kind = m.lastgroup
        if kind == "ws":
            continue
        if kind == "comment":
            continue
        out.append((kind, m.group()))
    return out


def _split_key(key: str) -> list:
    parts = []
    for name, idx in re.findall(r"([A-Za-z_]\w*)|\[(\d+)\]", key):
        parts.append(name if name else int(idx))
    return parts


def _lookup(root: dict, key: str):
    node = root
    for p in _split_key(key):
        node = node[p]
    return copy.deepcopy(node)


def _assign(target: dict, key: str, value):
    parts = _split_key(key)
    node  = target
    for p in parts[:-1]:
        if isinstance(p, str) and p not in node:
            node[p] = {}
        node = node[p]
    last = parts[-1]
    if value is ERASE:
        if isinstance(node, dict):
            node.pop(last, None)
    else:
        node[last] = value


class FhiclParser:
    """Parse a FHiCL document; text must already be expanded (no #includes)."""

    def __init__(self, text: str):
        self.toks  = _tokenize(text)
        self.i     = 0
        self.root  = {}
        self.order = []   # top-level keys in definition order

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else (None, None)

    def _next(self):
        t = self.peek()
        self.i += 1
        return t

    def _expect(self, val: str):
        k, v = self._next()
        if v != val:
            raise SyntaxError(f"expected {val!r}, got {v!r}")

    def parse(self) -> dict:
        while self.peek()[0] is not None:
            k, v = self.peek()
            if v in ("BEGIN_PROLOG", "END_PROLOG"):
                self._next()
                continue
            self._statement(self.root, top=True)
        return self.root

    def _statement(self, table: dict, top: bool = False):
        k, v = self._next()
        if k == "ref" and v.startswith("@table::"):
            table.update(_lookup(self.root, v.split("::", 1)[1]))
            return
        if k != "ident":
            raise SyntaxError(f"expected key, got {v!r}")
        kk, vv = self._next()
        if not (vv == ":" or kk == "protect"):
            raise SyntaxError(f"expected ':' after {v!r}, got {vv!r}")
        value = self._value()
        _assign(table, v, value)
        if top:
            base = _split_key(v)[0]
            if base not in self.order:
                self.order.append(base)

    def _value(self):
        k, v = self._next()
        if k == "erase":   return ERASE
        if k == "nil":     return None
        if k == "number":  return float(v)
        if k == "string":
            return (bytes(v[1:-1], "utf-8").decode("unicode_escape")
                    if v[0] == '"' else v[1:-1])
        if k == "ident":
            low = v.lower()
            if low in ("true", "false"):    return low == "true"
            if low in ("infinity", "+infinity"): return math.inf
            if low == "-infinity":           return -math.inf
            return v
        if k == "ref":
            return _lookup(self.root, v.split("::", 1)[1])
        if v == "{":
            tbl = {}
            while self.peek()[1] != "}":
                self._statement(tbl)
            self._next()
            return tbl
        if v == "[":
            seq = []
            while self.peek()[1] != "]":
                if (self.peek()[0] == "ref"
                        and self.peek()[1].startswith("@sequence::")):
                    seq.extend(_lookup(self.root, self._next()[1].split("::", 1)[1]))
                else:
                    seq.append(self._value())
                if self.peek()[1] == ",":
                    self._next()
            self._next()
            return seq
        raise SyntaxError(f"unexpected token {v!r}")


# ---------------------------------------------------------------------------
# distrib_x helper
# ---------------------------------------------------------------------------

_TF1_MAP = [
    ("TMath::Exp",   "np.exp"),
    ("TMath::Abs",   "np.abs"),
    ("TMath::Power", "np.power"),
    ("TMath::Sqrt",  "np.sqrt"),
    ("TMath::Log",   "np.log"),
    ("TMath::Gaus",  "_gaus"),
]


def _gaus(x, mean=0.0, sigma=1.0):
    return np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def x_interval(formula: str, x0: float, x1: float,
               containment: float, npts: int = 400001) -> tuple:
    """Shortest [a,b] ⊂ [x0,x1] containing `containment` fraction of the pdf."""
    expr = formula
    for a, b in _TF1_MAP:
        expr = expr.replace(a, b)
    x = np.linspace(x0, x1, npts)
    f = eval(expr, {"np": np, "_gaus": _gaus, "x": x})   # noqa: S307
    f = np.clip(np.broadcast_to(f, x.shape).astype(float), 0, None)
    c = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])
    c /= c[-1]
    j  = np.searchsorted(c, c + containment)
    ok = j < len(x)
    i  = np.nonzero(ok)[0]
    k  = np.argmin(x[j[ok]] - x[i])
    return float(x[i[k]]), float(x[j[ok][k]])


# ---------------------------------------------------------------------------
# GDML-backed geometry resolution for volume_rand / volume_gen
# ---------------------------------------------------------------------------

def _world_aabb(geo, patterns: list[str]) -> tuple | None:
    """
    Axis-aligned bounding box (lo, hi) in world coordinates for all placements
    of volumes matching `patterns`.  Returns None if nothing is found.
    """
    from ._draw import box_corners_world

    hits = geo.find_positions(patterns)
    if not hits:
        return None
    lo = np.full(3,  np.inf)
    hi = np.full(3, -np.inf)
    found = 0
    for vname, pos, rot in hits:
        ext = geo.box_extents(vname)
        if ext is None:
            continue
        corners = box_corners_world(pos, ext, rot)
        lo = np.minimum(lo, corners.min(axis=0))
        hi = np.maximum(hi, corners.max(axis=0))
        found += 1
    return (lo.tolist(), hi.tolist()) if found else None


def _world_aabb_per_hit(geo, patterns: list[str]) -> list[tuple]:
    """
    One (vol_name, lo, hi) per individual placement — used for volume_gen
    when the caller wants one row per placed volume.
    """
    from ._draw import box_corners_world

    out = []
    for vname, pos, rot in geo.find_positions(patterns):
        ext = geo.box_extents(vname)
        if ext is None:
            continue
        corners = box_corners_world(pos, ext, rot)
        out.append((vname, corners.min(axis=0).tolist(),
                    corners.max(axis=0).tolist()))
    return out


# ---------------------------------------------------------------------------
# FHiCL → radio volume rows
# ---------------------------------------------------------------------------

BOX_KEYS = ("X0", "X1", "Y0", "Y1", "Z0", "Z1")
RADIO_MODULES = ("RadioGen", "Decay0Gen")


def _as_list(v, n: int) -> list:
    if isinstance(v, list):
        if len(v) != n:
            raise ValueError(f"array length mismatch ({len(v)} vs {n})")
        return v
    return [v] * n


def _surface_label(lo: list, hi: list) -> str:
    ext = [h - lo_i for lo_i, h in zip(lo, hi)]
    ax  = int(np.argmin(ext))
    other = sorted(e for i, e in enumerate(ext) if i != ax)
    if ext[ax] < 0.05 * other[0]:
        return f"plane_{'xyz'[ax]}@{0.5 * (lo[ax] + hi[ax]):g}"
    return "bulk"


def producer_volumes(name: str, cfg: dict, geo=None,
                     xrange=None, containment: float = 0.99) -> list[dict]:
    """
    Extract emission volume rows from a single RadioGen/Decay0Gen config dict.

    Parameters
    ----------
    geo : GDMLGeometry | None
        When provided, used to resolve ``volume_rand`` and ``volume_gen``
        entries to world-frame bounding boxes.

    Returns a list of dicts with keys:
      producer, module_type, index, source, nuclide, material,
      lo [x,y,z], hi [x,y,z], surface, T0, T1, BqPercc, rate, distrib_x, x_eff
    """
    mt   = cfg.get("module_type")
    rows = []
    base = dict(producer=name, module_type=mt,
                material=cfg.get("material", ".*"),
                rate=cfg.get("rate"), BqPercc=cfg.get("BqPercc"))

    if any(k in cfg for k in BOX_KEYS):
        missing = [k for k in BOX_KEYS if k not in cfg]
        if missing:
            raise ValueError(f"{name}: incomplete box, missing {missing}")
        n    = max(len(cfg[k]) if isinstance(cfg[k], list) else 1 for k in BOX_KEYS)
        cols = {k: _as_list(cfg[k], n) for k in BOX_KEYS + ("T0", "T1") if k in cfg}
        nucl = _as_list(cfg.get("Nuclide", cfg.get("isotope", "")), n)
        bq   = _as_list(cfg.get("BqPercc"), n) if "BqPercc" in cfg else [None] * n
        for i in range(n):
            lo = [cols["X0"][i], cols["Y0"][i], cols["Z0"][i]]
            hi = [cols["X1"][i], cols["Y1"][i], cols["Z1"][i]]
            rows.append(dict(base, index=i, source="explicit",
                             nuclide=nucl[i], BqPercc=bq[i],
                             lo=lo, hi=hi,
                             T0=cols.get("T0", [None]*n)[i],
                             T1=cols.get("T1", [None]*n)[i]))
    elif "volume_rand" in cfg:
        vr = cfg["volume_rand"]
        lo = hi = None
        if geo is not None:
            bbox = _world_aabb(geo, [f"^{re.escape(vr)}$"])
            if bbox is not None:
                lo, hi = bbox
            else:
                print(f"[warn] {name}: volume_rand '{vr}' not found in geometry",
                      file=sys.stderr)
        elif xrange:
            lo = [xrange[0], None, None]
            hi = [xrange[1], None, None]
        rows.append(dict(base, index=0, source=f"volume_rand:{vr}", lo=lo, hi=hi))
    else:
        vg = cfg.get("volume_gen", ".*")
        if geo is not None:
            hits = _world_aabb_per_hit(geo, [vg])
            if hits:
                for i, (vname, lo, hi) in enumerate(hits):
                    rows.append(dict(base, index=i,
                                     source=f"volume_gen:{vg}",
                                     node=vname, lo=lo, hi=hi))
            else:
                print(f"[warn] {name}: volume_gen '{vg}' not found in geometry",
                      file=sys.stderr)
                rows.append(dict(base, index=0, source=f"volume_gen:{vg}",
                                 lo=None, hi=None))
        else:
            rows.append(dict(base, index=0, source=f"volume_gen:{vg}",
                             lo=None, hi=None))

    if mt != "RadioGen":
        iso = cfg.get("isotope")
        if iso is None and isinstance(cfg.get("decay_chain"), dict):
            iso = "+".join(cfg["decay_chain"][k] for k in
                           sorted(cfg["decay_chain"], key=lambda s: int(s.split("_")[-1])))
        for r in rows:
            r.setdefault("nuclide", iso)

    if "distrib_x" in cfg:
        for r in rows:
            r["distrib_x"] = cfg["distrib_x"]
            if r["lo"] is not None and r["lo"][0] is not None:
                a, b = x_interval(cfg["distrib_x"], r["lo"][0], r["hi"][0], containment)
                r["x_eff"] = [a, b]
                r["lo"] = [a] + r["lo"][1:]
                r["hi"] = [b] + r["hi"][1:]

    for r in rows:
        if r["lo"] is not None and None not in r["lo"]:
            r["surface"] = _surface_label(r["lo"], r["hi"])
    return rows


def _scheduled_labels(physics: dict) -> list[str]:
    if "trigger_paths" in physics:
        paths = list(physics["trigger_paths"])
    else:
        ends  = set(physics.get("end_paths", []))
        paths = [k for k, v in physics.items()
                 if isinstance(v, list)
                 and k not in ("trigger_paths", "end_paths")
                 and k not in ends]
    labels = []
    for path in paths:
        for lab in physics.get(path, []):
            if isinstance(lab, str) and lab not in labels:
                labels.append(lab)
    return labels


def select_producers(doc: dict, top_order: list,
                     prefix: str = "",
                     module_types: tuple = RADIO_MODULES,
                     all_configured: bool = False,
                     select_re: str | None = None) -> tuple:
    """
    Return (table, [labels]) for the radiological producers in `doc`.

    Handles both full job configs (physics.producers) and prolog-only files.
    """
    types   = set(module_types)
    physics = doc.get("physics")
    if isinstance(physics, dict) and isinstance(physics.get("producers"), dict):
        table = physics["producers"]
        if all_configured:
            labels = list(table)
        else:
            labels    = _scheduled_labels(physics)
            unknown   = [l for l in labels
                         if l not in table
                         and l not in doc.get("outputs", {})
                         and l not in physics.get("filters", {})
                         and l not in physics.get("analyzers", {})]
            for l in unknown:
                print(f"[warn] scheduled label '{l}' has no configuration", file=sys.stderr)
            labels       = [l for l in labels if l in table]
            unscheduled  = [l for l in table
                            if l not in labels
                            and table[l].get("module_type") in types]
            if unscheduled:
                print(f"[info] {len(unscheduled)} unscheduled radio producer(s) skipped "
                      f"(use --all-configured): {', '.join(unscheduled)}", file=sys.stderr)
        src = "physics.producers"
    else:
        table, labels, src = doc, list(top_order), "top-level tables"

    sel = re.compile(select_re) if select_re else None
    labels = [l for l in labels
              if isinstance(table.get(l), dict)
              and table[l].get("module_type") in types
              and l.startswith(prefix)
              and (sel is None or sel.search(l))]
    print(f"[info] {len(labels)} radiological producer(s) selected from {src}", file=sys.stderr)
    return table, labels

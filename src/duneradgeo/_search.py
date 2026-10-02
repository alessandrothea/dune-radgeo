"""Search-path resolution and GDML name extraction from FHiCL documents."""
from __future__ import annotations

import os
import sys
from pathlib import Path


def resolve_search_path(fname: str, var: str = "FW_SEARCH_PATH") -> Path | None:
    """
    Resolve a bare filename via $FW_SEARCH_PATH (first-match-wins, mirrors
    cet::search_path).  Returns None if not found.
    """
    p = Path(fname)
    if p.is_absolute():
        return p if p.is_file() else None
    search = os.environ.get(var, "")
    for d in (d for d in search.split(":") if d):
        cand = Path(d) / fname
        if cand.is_file():
            print(f"[info] resolved '{fname}' → {cand}", file=sys.stderr)
            return cand
    return None


def _find_tables(node, key: str, path: str = ""):
    """Yield (dotted_path, table) for every nested table named `key`."""
    if isinstance(node, dict):
        for k, v in node.items():
            p = f"{path}.{k}" if path else k
            if k == key and isinstance(v, dict):
                yield p, v
            yield from _find_tables(v, key, p)
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from _find_tables(v, key, f"{path}[{i}]")


def gdml_from_fcl(doc: dict) -> str | None:
    """Return the GDML filename declared in AuxDetGeometry.GDML (or None)."""
    names: dict[str, list] = {}
    for path, tbl in _find_tables(doc, "AuxDetGeometry"):
        if isinstance(tbl.get("GDML"), str):
            names.setdefault(tbl["GDML"], []).append(path)
    if not names:
        return None
    if len(names) > 1:
        detail = "; ".join(f"{g} ({', '.join(p)})" for g, p in names.items())
        raise ValueError(f"conflicting AuxDetGeometry.GDML values: {detail}")
    gdml, paths = next(iter(names.items()))
    print(f"[info] GDML '{gdml}' from {', '.join(paths)}", file=sys.stderr)
    return gdml

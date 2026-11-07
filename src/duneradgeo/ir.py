"""Intermediate representation: VolumeEntry + IRCollection (JSON I/O)."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class VolumeEntry:
    """One placed detector or radiological volume in the IR frame
    (the world, unless IRCollection.frame names another volume)."""
    kind: str          # "anode" | "opdet" | "cathode" | "fieldshaper" | "struct" | "radio"
    name: str          # logical-volume or producer name
    pos:  np.ndarray   # (3,) centre in the IR frame [cm]
    ext:  np.ndarray   # (3,) full extents dx,dy,dz [cm]
    rot:  np.ndarray   # (3,3) rotation matrix in the IR frame
    meta: dict         # kind-specific metadata; multi-part solids: part, copy

    # ------------------------------------------------------------------
    def to_dict(self) -> dict:
        return {
            "kind": self.kind,
            "name": self.name,
            "pos":  self.pos.tolist(),
            "ext":  self.ext.tolist(),
            "rot":  self.rot.tolist(),
            "meta": self.meta,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "VolumeEntry":
        return cls(
            kind=d["kind"],
            name=d["name"],
            pos=np.array(d["pos"], dtype=float),
            ext=np.array(d["ext"], dtype=float),
            rot=np.array(d["rot"], dtype=float),
            meta=d.get("meta", {}),
        )


@dataclass
class IRCollection:
    """Container for all extracted volumes, with JSON save/load."""
    source:   str                  = ""
    detector: str                  = ""   # "vd" | "hd" | ""
    volumes:  list[VolumeEntry]    = field(default_factory=list)
    version:  str                  = "1"
    # None: world coordinates.  Otherwise {"volume", "path", "pos", "rot"}:
    # the volume whose frame pos/rot are given in, and its world placement
    # (world = pos + rot @ local), so IRs of differently nested GDMLs overlay.
    frame:    dict | None          = None
    # Full FHiCL configuration of each radiological producer, keyed by module
    # label (FCL input only).  Radio volume rows refer to it by name.
    producers: dict[str, dict]     = field(default_factory=dict)

    # ------------------------------------------------------------------
    def save(self, path: str | Path) -> None:
        d = {
            "version":  self.version,
            "source":   self.source,
            "detector": self.detector,
        }
        if self.frame is not None:
            d["frame"] = self.frame
        if self.producers:
            d["producers"] = self.producers
        d["volumes"] = [v.to_dict() for v in self.volumes]
        Path(path).write_text(json.dumps(d, indent=2))

    @classmethod
    def load(cls, path: str | Path) -> "IRCollection":
        d = json.loads(Path(path).read_text())
        return cls(
            version=d.get("version", "1"),
            source=d.get("source", ""),
            detector=d.get("detector", ""),
            volumes=[VolumeEntry.from_dict(v) for v in d.get("volumes", [])],
            frame=d.get("frame"),
            producers=d.get("producers", {}),
        )

    # -- convenience filters -------------------------------------------
    def by_kind(self, *kinds: str) -> list[VolumeEntry]:
        return [v for v in self.volumes if v.kind in kinds]

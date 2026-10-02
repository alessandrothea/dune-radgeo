"""Intermediate representation: VolumeEntry + IRCollection (JSON I/O)."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class VolumeEntry:
    """One placed detector or radiological volume in world coordinates."""
    kind: str          # "anode" | "opdet" | "struct" | "radio"
    name: str          # logical-volume or producer name
    pos:  np.ndarray   # (3,) world-frame centre [cm]
    ext:  np.ndarray   # (3,) full extents dx,dy,dz [cm]
    rot:  np.ndarray   # (3,3) world-frame rotation matrix
    meta: dict         # kind-specific metadata

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

    # ------------------------------------------------------------------
    def save(self, path: str | Path) -> None:
        d = {
            "version":  self.version,
            "source":   self.source,
            "detector": self.detector,
            "volumes":  [v.to_dict() for v in self.volumes],
        }
        Path(path).write_text(json.dumps(d, indent=2))

    @classmethod
    def load(cls, path: str | Path) -> "IRCollection":
        d = json.loads(Path(path).read_text())
        return cls(
            version=d.get("version", "1"),
            source=d.get("source", ""),
            detector=d.get("detector", ""),
            volumes=[VolumeEntry.from_dict(v) for v in d.get("volumes", [])],
        )

    # -- convenience filters -------------------------------------------
    def by_kind(self, *kinds: str) -> list[VolumeEntry]:
        return [v for v in self.volumes if v.kind in kinds]

"""
Geometric measurements on an IR, run by bki-measure.

Each measurement lives in its own module in this package and registers
itself with ``@measurement("name", "one-line help")``.  Every module here
whose name does not start with ``_`` is imported automatically (in name
order, which is also the default run order), so adding a measurement means
adding one file; nothing else needs to change.
"""
import importlib
import pkgutil

from ._core import (MEASUREMENTS, X, YZ, Measurement, Result, aabb,
                    box_distances, measurement, overlap_area, run,
                    struct_box, tags)

for _mod in sorted(pkgutil.iter_modules(__path__), key=lambda m: m.name):
    if not _mod.name.startswith("_"):
        importlib.import_module(f"{__name__}.{_mod.name}")

__all__ = ["MEASUREMENTS", "Measurement", "Result", "X", "YZ", "aabb",
           "box_distances", "measurement", "overlap_area", "run",
           "struct_box", "tags"]

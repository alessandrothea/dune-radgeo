"""The PyROOT import shared by the modules that need ROOT (_tf1, _gdml_root).

ROOT is imported on first use, never at module import, so the package works
without it wherever no ROOT feature is used."""
from __future__ import annotations

import warnings

_ROOT = None


class ROOTUnavailable(RuntimeError):
    """PyROOT cannot be imported, or lacks a feature that is needed."""


def import_root(feature: str):
    """The ROOT module, imported once and set to batch mode.

    `feature` names what needs ROOT, for the message of the ROOTUnavailable
    raised when PyROOT cannot be imported."""
    global _ROOT
    if _ROOT is None:
        with warnings.catch_warnings():
            # the root pip wheel warns that it is an alpha release on import
            warnings.filterwarnings("ignore", message=r"(?s).*alpha stage")
            try:
                import ROOT
            except ImportError as exc:
                raise ROOTUnavailable(
                    f"{feature} needs ROOT, but PyROOT cannot be imported: use "
                    "the dune-radgeo venv (uv sync, dependency group 'root') "
                    "or a dunesw area") from exc
        ROOT.gROOT.SetBatch(True)
        _ROOT = ROOT
    return _ROOT

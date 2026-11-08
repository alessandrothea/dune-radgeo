"""ROOT TF1 evaluation of the RadioGen/Decay0Gen distrib_x/y/z formulas.

The code's one full account of how the generator samples these formulas (the
other modules point here; docs/radgeo-profiles.md is the user-facing account).

larsimrad's BaseRadioGen samples each axis independently.  distrib_<axis> is
a 1D formula in the variable x, built as ``TF1(name, formula, lo, hi)`` over
the box range on that axis (X0..X1, the volume_rand bounding box, or a
volume_gen node's box) and drawn with ``GetRandom(lo, hi)``.  The formulas are
ROOT TFormula expressions, so they are evaluated by ROOT itself, never
translated.

GetRandom draws from a table, not from the formula: it integrates the formula
over Npx equal bins of [lo, hi], approximates the cumulative inside each bin
by a parabola, and inverts that.  BaseRadioGen never calls SetNpx, so ROOT's
default Npx = 100 applies (:data:`GENERATOR_NPX`): over the 850 cm of
volCryostat the bins are 8.5 cm wide.  Each bin keeps its probability, but
the generated density is linear inside it, so a feature narrower than
(hi - lo) / Npx is not reproduced: the upper-mesh K42 Gaussian (sigma 0.5 cm
at x = -320 cm on a constant) is generated as a ramp across its 8.5 cm bin, as
the decay vertices of a lar job confirm (tests/test_tf1.py).  Smooth profiles
are not affected.

The generated distribution is read from ROOT's own table, without sampling:
:meth:`Distribution.quantiles` calls TF1::GetQuantiles, which builds the same
table but with max(Npx, 2 n) points for n probabilities, so it is asked at
most Npx / 2 probabilities per call to stay on GetRandom's Npx-point table.
The CDF tabulated from it at 20001 probabilities gives
:meth:`~Distribution.cdf`, :meth:`~Distribution.pdf` and
:meth:`~Distribution.interval`; :meth:`~Distribution.sample` calls GetRandom
itself.  The formula is evaluated with TF1::Eval
(:meth:`~Distribution.eval`), clipped at 0 and normalised over [lo, hi] on a
dense 400001-point grid (TF1::Integral misses peaks much narrower than the
range): :meth:`~Distribution.formula_pdf`,
:meth:`~Distribution.formula_interval`, :meth:`~Distribution.formula_fraction`.

An effective range is the shortest interval holding a fraction `containment`
of a distribution.  radgeo extract writes both, per axis: that of the generated
distribution (``distrib_<axis>_eff_range``, which becomes the box on that
axis) and that of the formula (``distrib_<axis>_formula_eff_range``).  For the
mesh Gaussian at 0.99 they are [-323.08, -318.75] and [-321.70, -318.33] cm;
with Npx = 10000 (:meth:`~Distribution.with_npx`) the generated one matches
the formula's.

A :class:`Distribution` is one formula on one range with one Npx, best
obtained from the cached factory :func:`distribution`.

Not modelled: the generator also redraws points that land outside the
volumes allowed by volume_gen / material (e.g. the gas above the liquid).

Needs PyROOT: the ``root`` pip wheel (uv dependency group) on the host, or the
UPS ROOT of a dunesw area.  TF1 is in ROOT's core, so any build will do.
"""
from __future__ import annotations

import functools
import itertools
from functools import cached_property

import numpy as np

from ._root import import_root

GENERATOR_NPX = 100      # ROOT's TF1 default, kept by BaseRadioGen (no SetNpx)

_FORMULA_GRID = 400001   # points of the dense formula grid: TF1::Integral misses
                         # peaks much narrower than the range
_CDF_TABLE = 20001       # probabilities at which the GetRandom CDF is tabulated
_SEED = 12345            # gRandom seed of Distribution.sample

_names = itertools.count()     # TF1 names: unique for the life of the process
_declared = False
_HELPERS = r"""
#include "TF1.h"
namespace duneradgeo_tf1 {
void eval(TF1& f, const double* x, double* y, long n) {
  for (long i = 0; i < n; ++i) y[i] = f.Eval(x[i]);
}
void sample(TF1& f, double lo, double hi, double* out, long n) {
  for (long i = 0; i < n; ++i) out[i] = f.GetRandom(lo, hi);
}
}
"""


def _root():
    """ROOT, with the C++ loops of _HELPERS declared."""
    global _declared
    ROOT = import_root("evaluating the distrib_x/y/z formulas")
    if not _declared:
        ROOT.gInterpreter.Declare(_HELPERS)
        _declared = True
    return ROOT


def _trapezoids(f: np.ndarray, x: np.ndarray) -> np.ndarray:
    """The trapezoid rule's terms for f over x: their sum is the integral,
    their cumulative sum the running integral (np.trapz is gone in numpy 2,
    np.trapezoid absent in 1.x)."""
    return 0.5 * (f[1:] + f[:-1]) * np.diff(x)


def _integrate(f: np.ndarray, x: np.ndarray) -> float:
    return float(np.sum(_trapezoids(f, x)))


def _shortest(xs: np.ndarray, cdf: np.ndarray, containment: float) -> tuple:
    """Shortest [a, b] with cdf(b) - cdf(a) >= containment, on a tabulated cdf."""
    j  = np.searchsorted(cdf, cdf + containment)
    ok = j < len(xs)
    i  = np.nonzero(ok)[0]
    if not len(i):
        return float(xs[0]), float(xs[-1])
    k  = np.argmin(xs[j[ok]] - xs[i])
    return float(xs[i[k]]), float(xs[j[ok][k]])


class Distribution:
    """distrib_<axis> `formula` on [lo, hi]: the formula itself, and what
    TF1::GetRandom generates from it with an `npx`-point table.

    Owns its TF1 (kept out of ROOT's global list of functions, so it is
    deleted with this object and its name never clashes).  Raises ValueError
    when ROOT cannot parse the formula, duneradgeo._root.ROOTUnavailable
    without PyROOT."""

    def __init__(self, formula: str, lo: float, hi: float, npx: int = GENERATOR_NPX):
        ROOT = _root()
        self.formula = str(formula)
        self.lo, self.hi, self.npx = float(lo), float(hi), int(npx)
        f = ROOT.TF1(f"duneradgeo_tf1_{next(_names)}", self.formula, self.lo, self.hi)
        f.AddToGlobalList(False)      # owned by Python, not by gROOT
        if not f.IsValid():
            raise ValueError(f"ROOT cannot parse the formula {formula!r}")
        f.SetNpx(self.npx)            # explicit, so the result does not depend on ROOT's default
        self.tf1 = f

    def __repr__(self) -> str:
        return (f"Distribution({self.formula!r}, {self.lo!r}, {self.hi!r}, "
                f"npx={self.npx})")

    def with_npx(self, npx: int) -> "Distribution":
        """The same formula and range with another Npx (cached)."""
        return distribution(self.formula, self.lo, self.hi, npx)

    # ---- the formula ------------------------------------------------------

    def eval(self, x) -> np.ndarray:
        """TF1::Eval at every point of `x` (raw values, not clipped; TF1::Eval
        ignores the range)."""
        x = np.ascontiguousarray(x, dtype=float)
        y = np.empty_like(x)
        if x.size:
            _root().duneradgeo_tf1.eval(self.tf1, x, y, x.size)
        return y

    def _formula_grid(self) -> tuple:
        """(x, f, norm): the formula on the dense grid over [lo, hi], negative
        values clipped to 0, and its integral (trapezoid sum).  Not cached
        (6 MB); the norm is."""
        x = np.linspace(self.lo, self.hi, _FORMULA_GRID)
        f = np.clip(self.eval(x), 0, None)
        norm = _integrate(f, x)
        if not norm > 0:
            raise ValueError("formula integrates to 0 over the range")
        return x, f, norm

    @cached_property
    def norm(self) -> float:
        """Integral of the formula (clipped at 0) over [lo, hi] on the dense
        grid.  Raises ValueError when it is 0."""
        return self._formula_grid()[2]

    def formula_pdf(self, x) -> np.ndarray:
        """The formula at `x`, negative values clipped to 0, normalised to 1
        over [lo, hi]."""
        return np.clip(self.eval(x), 0, None) / self.norm

    def formula_interval(self, containment: float) -> tuple:
        """Shortest [a, b] ⊂ [lo, hi] holding `containment` of the formula pdf."""
        x, f, norm = self._formula_grid()
        c = np.concatenate([[0.0], np.cumsum(_trapezoids(f / norm, x))])
        return _shortest(x, c / c[-1], containment)

    def formula_fraction(self, a: float, b: float) -> float:
        """Fraction of the formula pdf (normalised over [lo, hi]) inside
        [a, b], integrated on its own grid over [a, b]."""
        xa = np.linspace(a, b, _FORMULA_GRID)
        return _integrate(np.clip(self.eval(xa), 0, None), xa) / self.norm

    # ---- the generated distribution (TF1::GetRandom) ----------------------

    def quantiles(self, p) -> np.ndarray:
        """Quantiles of the GetRandom(lo, hi) distribution at `p` (increasing
        in [0, 1]), from TF1::GetQuantiles: deterministic, no sampling."""
        # at most Npx / 2 probabilities per call: see the module docstring
        p = np.ascontiguousarray(p, dtype=float)
        q = np.empty_like(p)
        step = max(1, self.npx // 2)
        for i in range(0, p.size, step):
            pc = np.ascontiguousarray(p[i:i + step])
            qc = np.empty_like(pc)
            self.tf1.GetQuantiles(pc.size, qc, pc)
            q[i:i + step] = qc
        return q

    @cached_property
    def _cdf_table(self) -> tuple:
        """(q, p): the GetRandom CDF tabulated at _CDF_TABLE probabilities."""
        p = np.linspace(0.0, 1.0, _CDF_TABLE)
        q = np.maximum.accumulate(self.quantiles(p))   # guard rounding
        return q, p

    def cdf(self, x) -> np.ndarray:
        """CDF of the GetRandom distribution at `x`."""
        q, p = self._cdf_table
        return np.interp(np.asarray(x, dtype=float), q, p, left=0.0, right=1.0)

    def pdf(self, x) -> np.ndarray:
        """Density of the GetRandom distribution at `x` [1/len], the derivative
        of the tabulated CDF.  Normalised over [lo, hi], like :meth:`formula_pdf`."""
        x = np.asarray(x, dtype=float)
        xe = np.concatenate([[x[0]], 0.5 * (x[1:] + x[:-1]), [x[-1]]])
        c = self.cdf(xe)
        w = np.diff(xe)
        return np.divide(np.diff(c), w, out=np.zeros_like(w), where=w > 0)

    def interval(self, containment: float) -> tuple:
        """Shortest [a, b] holding `containment` of the GetRandom distribution."""
        q, p = self._cdf_table
        return _shortest(q, p, containment)

    def sample(self, n: int, seed: int = _SEED) -> np.ndarray:
        """`n` positions drawn as the generator does, GetRandom(lo, hi), after
        gRandom.SetSeed(seed)."""
        ROOT = _root()
        ROOT.gRandom.SetSeed(seed)
        out = np.empty(int(n))
        ROOT.duneradgeo_tf1.sample(self.tf1, self.lo, self.hi, out, out.size)
        return out


@functools.lru_cache(maxsize=256)
def _cached(formula: str, lo: float, hi: float, npx: int) -> Distribution:
    return Distribution(formula, lo, hi, npx)


def distribution(formula: str, lo: float, hi: float,
                 npx: int = GENERATOR_NPX) -> Distribution:
    """The :class:`Distribution` of `formula` on [lo, hi] with `npx` points,
    shared by every caller with the same arguments (the 256 most recent are
    kept; ``distribution.cache_clear()`` drops them)."""
    return _cached(str(formula), float(lo), float(hi), int(npx))


distribution.cache_clear = _cached.cache_clear
distribution.cache_info = _cached.cache_info

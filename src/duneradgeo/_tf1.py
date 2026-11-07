"""ROOT TF1 evaluation of the RadioGen/Decay0Gen distrib_x/y/z formulas.

The formulas are ROOT TFormula expressions, so they are evaluated by ROOT
itself, never translated.  The TF1 is built as larsimrad's BaseRadioGen does,
``TF1(name, formula, lo, hi)``, and the generator draws positions with
``GetRandom(lo, hi)``.  GetRandom samples from a cumulative table of ``Npx``
points (BaseRadioGen keeps ROOT's default, 100) with a parabolic interpolation
inside each bin, so for a feature narrower than (hi - lo) / Npx the generated
distribution differs from the formula.  Two descriptions are provided:

* the formula: :func:`evaluate`, :func:`formula_interval`;
* the generated distribution: :func:`quantiles` (TF1::GetQuantiles, ROOT's
  own table, deterministic), :func:`generated_cdf`, :func:`generated_pdf`,
  :func:`generated_interval`, and :func:`sample` (TF1::GetRandom).

Not modelled: the generator also redraws points that land outside the
volumes allowed by volume_gen / material (e.g. the gas above the liquid).

Needs PyROOT: the ``root`` pip wheel (uv dependency group) on the host, or the
UPS ROOT of a dunesw area.  TF1 is in ROOT's core, so any build will do.
"""
from __future__ import annotations

import warnings

import numpy as np

GENERATOR_NPX = 100      # ROOT's TF1 default, kept by BaseRadioGen (no SetNpx)

_ROOT = None
_CACHE: dict = {}
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


class ROOTUnavailable(RuntimeError):
    """PyROOT cannot be imported: distrib_* formulas cannot be evaluated."""


def _root():
    global _ROOT
    if _ROOT is None:
        with warnings.catch_warnings():
            # the root pip wheel warns that it is an alpha release on import
            warnings.filterwarnings("ignore", message=r"(?s).*alpha stage")
            try:
                import ROOT
            except ImportError as exc:
                raise ROOTUnavailable(
                    "distrib_x/y/z formulas are evaluated with ROOT, but PyROOT "
                    "cannot be imported: use the dune-radgeo venv (uv sync, "
                    "dependency group 'root') or a dunesw area") from exc
        ROOT.gROOT.SetBatch(True)
        ROOT.gInterpreter.Declare(_HELPERS)
        _ROOT = ROOT
    return _ROOT


def tf1(formula: str, lo: float, hi: float, npx: int = GENERATOR_NPX):
    """The TF1 for `formula` on [lo, hi] with `npx` points (cached).

    Raises ValueError when ROOT cannot parse the formula."""
    key = (formula, float(lo), float(hi), int(npx))
    f = _CACHE.get(key)
    if f is None:
        ROOT = _root()
        f = ROOT.TF1(f"duneradgeo_tf1_{len(_CACHE)}", formula, float(lo), float(hi))
        if not f.IsValid():
            raise ValueError(f"ROOT cannot parse the formula {formula!r}")
        f.SetNpx(int(npx))     # explicit, so the result does not depend on ROOT's default
        _CACHE[key] = f
    return f


# ---- the formula ----------------------------------------------------------

def evaluate(formula: str, x) -> np.ndarray:
    """TF1::Eval of `formula` at every point of `x` (raw values, not clipped)."""
    x = np.ascontiguousarray(x, dtype=float)
    y = np.empty_like(x)
    if x.size:
        lo, hi = float(x.min()), float(x.max())      # TF1::Eval ignores the range
        f = tf1(formula, lo, hi if hi > lo else lo + 1.0)
        _root().duneradgeo_tf1.eval(f, x, y, x.size)
    return y


def _integrate(f: np.ndarray, x: np.ndarray) -> float:
    """Trapezoid rule (np.trapz is gone in numpy 2, np.trapezoid absent in 1.x)."""
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(x)))


def formula_pdf(formula: str, lo: float, hi: float, npts: int = 400001) -> tuple:
    """(x, pdf, norm): the formula on a dense grid over [lo, hi], negative
    values clipped to 0, normalised to 1 over the range.  The integral is a
    dense-grid trapezoid sum: TF1::Integral misses peaks much narrower than
    the range."""
    x = np.linspace(lo, hi, npts)
    f = np.clip(evaluate(formula, x), 0, None)
    norm = _integrate(f, x)
    if not norm > 0:
        raise ValueError("formula integrates to 0 over the range")
    return x, f / norm, norm


def _shortest(xs: np.ndarray, cdf: np.ndarray, containment: float) -> tuple:
    """Shortest [a, b] with cdf(b) - cdf(a) >= containment, on a tabulated cdf."""
    j  = np.searchsorted(cdf, cdf + containment)
    ok = j < len(xs)
    i  = np.nonzero(ok)[0]
    if not len(i):
        return float(xs[0]), float(xs[-1])
    k  = np.argmin(xs[j[ok]] - xs[i])
    return float(xs[i[k]]), float(xs[j[ok][k]])


def formula_interval(formula: str, lo: float, hi: float, containment: float,
                     npts: int = 400001) -> tuple:
    """Shortest [a, b] ⊂ [lo, hi] holding `containment` of the formula pdf."""
    x, pdf, _ = formula_pdf(formula, lo, hi, npts)
    c = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(x))])
    return _shortest(x, c / c[-1], containment)


def formula_fraction(formula: str, lo: float, hi: float, a: float, b: float,
                     npts: int = 400001) -> float:
    """Fraction of the formula pdf (normalised over [lo, hi]) inside [a, b],
    integrated on its own grid over [a, b]."""
    _, _, norm = formula_pdf(formula, lo, hi, npts)
    xa = np.linspace(a, b, npts)
    return _integrate(np.clip(evaluate(formula, xa), 0, None), xa) / norm


# ---- the generated distribution (TF1::GetRandom) --------------------------

def quantiles(formula: str, lo: float, hi: float, p,
              npx: int = GENERATOR_NPX) -> np.ndarray:
    """Quantiles of the GetRandom(lo, hi) distribution, from ROOT's own CDF
    table (TF1::GetQuantiles): deterministic, no sampling.  `p` increasing
    in [0, 1]."""
    # GetQuantiles builds its table with max(Npx, 2 * len(p)) points, so ask for
    # at most Npx / 2 probabilities per call to stay on GetRandom's Npx table.
    p = np.ascontiguousarray(p, dtype=float)
    q = np.empty_like(p)
    f = tf1(formula, lo, hi, npx)
    step = max(1, int(npx) // 2)
    for i in range(0, p.size, step):
        pc = np.ascontiguousarray(p[i:i + step])
        qc = np.empty_like(pc)
        f.GetQuantiles(pc.size, qc, pc)
        q[i:i + step] = qc
    return q


_CDF: dict = {}


def _cdf_table(formula: str, lo: float, hi: float, npx: int, n: int = 20001) -> tuple:
    """(q, p): the GetRandom CDF tabulated at n probabilities (cached)."""
    key = (formula, float(lo), float(hi), int(npx), n)
    if key not in _CDF:
        p = np.linspace(0.0, 1.0, n)
        q = np.maximum.accumulate(quantiles(formula, lo, hi, p, npx))   # guard rounding
        _CDF[key] = (q, p)
    return _CDF[key]


def generated_cdf(formula: str, lo: float, hi: float, x,
                  npx: int = GENERATOR_NPX) -> np.ndarray:
    """CDF of the GetRandom distribution at `x`."""
    q, p = _cdf_table(formula, lo, hi, npx)
    return np.interp(np.asarray(x, dtype=float), q, p, left=0.0, right=1.0)


def generated_pdf(formula: str, lo: float, hi: float, x,
                  npx: int = GENERATOR_NPX) -> np.ndarray:
    """Density of the GetRandom distribution at `x` [1/len], the derivative of
    the tabulated CDF.  Normalised over [lo, hi], like :func:`formula_pdf`."""
    x = np.asarray(x, dtype=float)
    xe = np.concatenate([[x[0]], 0.5 * (x[1:] + x[:-1]), [x[-1]]])
    c = generated_cdf(formula, lo, hi, xe, npx)
    w = np.diff(xe)
    return np.divide(np.diff(c), w, out=np.zeros_like(w), where=w > 0)


def generated_interval(formula: str, lo: float, hi: float, containment: float,
                       npx: int = GENERATOR_NPX) -> tuple:
    """Shortest [a, b] holding `containment` of the GetRandom distribution."""
    q, p = _cdf_table(formula, lo, hi, npx)
    return _shortest(q, p, containment)


def sample(formula: str, lo: float, hi: float, n: int,
           npx: int = GENERATOR_NPX, seed: int = 12345) -> np.ndarray:
    """`n` positions drawn as the generator does, GetRandom(lo, hi)."""
    ROOT = _root()
    ROOT.gRandom.SetSeed(seed)
    out = np.empty(int(n))
    ROOT.duneradgeo_tf1.sample(tf1(formula, lo, hi, npx), float(lo), float(hi), out, out.size)
    return out

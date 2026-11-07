"""Checks of the ROOT TF1 evaluation (duneradgeo._tf1).

Run with pytest, or directly: python tests/test_tf1.py (host venv or a dunesw
area via `dsw py`).  Needs PyROOT.

Reference: the upper-mesh K42 producer of
prodbackground_radiological_decay0_dunevd10kt_1x8x14_3view_30deg.fcl
(dunesw v10_26_00d05), sampled over the volCryostat x range.  LAR_QUANTILES are
the 0.5/16/50/84 % quantiles of the x of 2.37M decay vertices written by
Decay0Gen in a lar job (ROOT 6.28/12).  GetRandom's Npx = 100 table turns the
0.5 cm Gaussian into a ramp over its 8.5 cm bin [-323.08, -314.58] (see
duneradgeo._tf1).
"""
import numpy as np

from duneradgeo import _tf1

MESH = "(7.2257e+4*TMath::Exp(-(x+320.)*(x+320.)/(2.*0.5*0.5)))+1."
BULK = ("10255.2+(41.3788*TMath::Abs(325.-x))-(0.135462*(325.-x)*(325.-x))"
        "+(0.000203379*TMath::Power(TMath::Abs(325.-x),3))"
        "-(1.14356e-007*TMath::Power(TMath::Abs(325.-x),4))")
LO, HI = -425.12, 425.20
P = [0.005, 0.16, 0.50, 0.84]
LAR_QUANTILES = [-323.07, -322.60, -321.42, -319.82]


def setup_function(_fn=None):
    """Each test builds its own Distributions."""
    _tf1.distribution.cache_clear()


def _d(formula, npx=_tf1.GENERATOR_NPX):
    return _tf1.distribution(formula, LO, HI, npx)


def test_formula_values():
    """TF1::Eval agrees with the closed forms."""
    x = np.linspace(LO, HI, 2001)
    mesh = 7.2257e4 * np.exp(-(x + 320.0) ** 2 / (2 * 0.5 * 0.5)) + 1.0
    u = np.abs(325.0 - x)
    bulk = (10255.2 + 41.3788 * u - 0.135462 * u**2 + 0.000203379 * u**3
            - 1.14356e-7 * u**4)
    assert np.allclose(_d(MESH).eval(x), mesh, rtol=1e-12)
    assert np.allclose(_d(BULK).eval(x), bulk, rtol=1e-12)


def test_invalid_formula():
    try:
        _tf1.distribution("10*nosuchfunction(x)", LO, HI)
    except ValueError:
        return
    raise AssertionError("an unparsable formula must raise ValueError")


def test_narrow_peak_accepted():
    """bki-distrib validates a formula by its norm on the dense formula grid:
    a peak much narrower than the range is valid, although it vanishes at lo,
    mid and hi (a 3-point integral of it is 0); its effective range comes
    from the same grid."""
    d = _d("TMath::Gaus(x,100,1)")
    assert d.norm > 0
    form = d.formula_interval(0.99)
    assert np.allclose(form, (97.40, 102.55), atol=0.01), form


def test_zero_norm_rejected():
    """A formula that parses but vanishes on the range has no pdf."""
    d = _d("0*x")
    try:
        d.norm
    except ValueError:
        return
    raise AssertionError("a formula integrating to 0 must raise ValueError")


def test_generated_matches_lar():
    """GetRandom's table (Npx = 100) reproduces the lar decay vertices."""
    q = _d(MESH).quantiles(P)
    assert np.allclose(q, LAR_QUANTILES, atol=0.02), q


def test_quantile_batching():
    """Many probabilities give the same quantiles as a few (the Npx / 2
    batching of Distribution.quantiles)."""
    many = _d(MESH).quantiles(np.linspace(0, 1, 2001))
    few = _d(MESH).quantiles(P)
    assert np.allclose(np.interp(P, np.linspace(0, 1, 2001), many), few, atol=0.01)


def test_quantiles_match_sampling():
    d = _d(MESH)
    s = d.sample(1_000_000)
    assert np.allclose(np.percentile(s, [100 * p for p in P]),
                       d.quantiles(P), atol=0.02)


def test_effective_ranges():
    d = _d(MESH)
    gen = d.interval(0.99)
    form = d.formula_interval(0.99)
    assert np.allclose(gen, (-323.082, -318.753), atol=0.005), gen
    assert np.allclose(form, (-321.704, -318.328), atol=0.005), form
    c = d.cdf(list(gen))
    assert abs(c[1] - c[0] - 0.99) < 1e-3
    assert abs(d.formula_fraction(*form) - 0.99) < 1e-3
    # a fine table reproduces the formula
    fine = d.with_npx(10000)
    assert fine.npx == 10000 and fine.tf1.GetNpx() == 10000
    assert np.allclose(fine.interval(0.99), form, atol=0.05), fine.interval(0.99)


def test_smooth_profile_unaffected():
    """For a smooth profile the Npx = 100 table makes no difference."""
    d = _d(BULK)
    gen = d.interval(0.99)
    form = d.formula_interval(0.99)
    assert np.allclose(gen, form, atol=0.01), (gen, form)


def test_cache_and_ownership():
    """distribution() shares one object per (formula, range, npx), however
    the arguments are spelled; each owns a TF1 with a unique name that is
    not in ROOT's global list, so evicted objects leave nothing behind."""
    import ROOT
    a = _tf1.distribution(MESH, LO, HI)
    assert _tf1.distribution(MESH, np.float64(LO), HI, npx=100) is a
    assert _tf1.distribution(MESH, LO, HI, 200) is not a
    names = {_d(f"{k}+x").tf1.GetName() for k in range(5)}
    assert len(names) == 5
    assert not any(ROOT.gROOT.GetListOfFunctions().FindObject(n) for n in names)
    _tf1.distribution.cache_clear()
    assert _tf1.distribution(MESH, LO, HI) is not a
    assert np.array_equal(a.quantiles(P), _d(MESH).quantiles(P))


if __name__ == "__main__":
    import ROOT
    print("ROOT", ROOT.gROOT.GetVersion())
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            setup_function(fn)
            fn()
            print("ok  ", name)

"""Checks of the ROOT TF1 evaluation (duneradgeo._tf1).

Run with pytest, or directly: python tests/test_tf1.py (host venv or a dunesw
area via `dsw py`).  Needs PyROOT.

Reference: the upper-mesh K42 producer of
prodbackground_radiological_decay0_dunevd10kt_1x8x14_3view_30deg.fcl
(dunesw v10_26_00d05), sampled over the volCryostat x range.  LAR_QUANTILES are
the 0.5/16/50/84 % quantiles of the x of 2.37M decay vertices written by
Decay0Gen in a lar job (ROOT 6.28/12): BaseRadioGen draws with
TF1::GetRandom and ROOT's default Npx = 100, which turns the 0.5 cm Gaussian
into a ramp over the 8.5 cm table bin [-323.08, -314.58].
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


def test_formula_values():
    """TF1::Eval agrees with the closed forms."""
    x = np.linspace(LO, HI, 2001)
    mesh = 7.2257e4 * np.exp(-(x + 320.0) ** 2 / (2 * 0.5 * 0.5)) + 1.0
    u = np.abs(325.0 - x)
    bulk = (10255.2 + 41.3788 * u - 0.135462 * u**2 + 0.000203379 * u**3
            - 1.14356e-7 * u**4)
    assert np.allclose(_tf1.evaluate(MESH, x), mesh, rtol=1e-12)
    assert np.allclose(_tf1.evaluate(BULK, x), bulk, rtol=1e-12)


def test_invalid_formula():
    try:
        _tf1.tf1("10*nosuchfunction(x)", LO, HI)
    except ValueError:
        return
    raise AssertionError("an unparsable formula must raise ValueError")


def test_narrow_peak_accepted():
    """bki-distrib validates a formula by building its TF1: a peak much
    narrower than the range is valid, although it vanishes at lo, mid and hi
    (a 3-point integral of it is 0); its effective range comes from the dense
    formula grid."""
    narrow = "TMath::Gaus(x,100,1)"
    _tf1.tf1(narrow, LO, HI)
    form = _tf1.formula_interval(narrow, LO, HI, 0.99)
    assert np.allclose(form, (97.40, 102.55), atol=0.01), form


def test_generated_matches_lar():
    """GetRandom's table (Npx = 100) reproduces the lar decay vertices."""
    q = _tf1.quantiles(MESH, LO, HI, P)
    assert np.allclose(q, LAR_QUANTILES, atol=0.02), q


def test_quantile_batching():
    """Many probabilities give the same quantiles as a few: GetQuantiles
    would otherwise switch to a finer table (max(Npx, 2 * n) points)."""
    many = _tf1.quantiles(MESH, LO, HI, np.linspace(0, 1, 2001))
    few = _tf1.quantiles(MESH, LO, HI, P)
    assert np.allclose(np.interp(P, np.linspace(0, 1, 2001), many), few, atol=0.01)


def test_quantiles_match_sampling():
    s = _tf1.sample(MESH, LO, HI, 1_000_000)
    assert np.allclose(np.percentile(s, [100 * p for p in P]),
                       _tf1.quantiles(MESH, LO, HI, P), atol=0.02)


def test_effective_ranges():
    gen = _tf1.generated_interval(MESH, LO, HI, 0.99)
    form = _tf1.formula_interval(MESH, LO, HI, 0.99)
    assert np.allclose(gen, (-323.082, -318.753), atol=0.005), gen
    assert np.allclose(form, (-321.704, -318.328), atol=0.005), form
    c = _tf1.generated_cdf(MESH, LO, HI, list(gen))
    assert abs(c[1] - c[0] - 0.99) < 1e-3
    assert abs(_tf1.formula_fraction(MESH, LO, HI, *form) - 0.99) < 1e-3
    # a fine table reproduces the formula
    fine = _tf1.generated_interval(MESH, LO, HI, 0.99, npx=10000)
    assert np.allclose(fine, form, atol=0.05), fine


def test_smooth_profile_unaffected():
    """For a smooth profile the Npx = 100 table makes no difference."""
    gen = _tf1.generated_interval(BULK, LO, HI, 0.99)
    form = _tf1.formula_interval(BULK, LO, HI, 0.99)
    assert np.allclose(gen, form, atol=0.01), (gen, form)


if __name__ == "__main__":
    import ROOT
    print("ROOT", ROOT.gROOT.GetVersion())
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok  ", name)

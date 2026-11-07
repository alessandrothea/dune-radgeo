"""Helpers shared by the ROOT-free tests, which run with or without pytest.

``xfail(reason)`` marks a test of the correct behaviour that the current code
does not have (pytest.mark.xfail when pytest is importable); ``run_tests``
runs the ``test_*`` functions of a module for ``python tests/test_x.py`` and
skips the marked ones.
"""
import sys

try:
    import pytest
except ImportError:          # the direct runner does not need it
    pytest = None


def xfail(reason: str):
    def mark(fn):
        fn.xfail_reason = reason
        if pytest is not None:
            fn = pytest.mark.xfail(reason=reason, strict=False)(fn)
        return fn
    return mark


def run_tests(namespace: dict) -> None:
    """Run every ``test_*`` callable of `namespace` (a module's globals())."""
    failed = 0
    for name, fn in list(namespace.items()):
        if not (name.startswith("test_") and callable(fn)):
            continue
        reason = getattr(fn, "xfail_reason", None)
        if reason is not None:
            print("xfail", name, "(skipped:", reason + ")")
            continue
        try:
            fn()
        except Exception as exc:   # report every failure, then exit non-zero
            failed += 1
            print("FAIL ", name, f"{type(exc).__name__}: {exc}")
        else:
            print("ok   ", name)
    if failed:
        sys.exit(f"{failed} test(s) failed")

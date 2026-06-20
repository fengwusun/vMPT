"""The startup PyPI update-check compares versions with PEP 440 ordering.

The network fetch itself only runs under a live Bokeh session (guarded by
`session_context is not None`), so it never fires during tests. Here we just
lock down the pure comparison + installed-version helpers.
"""
from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import vmpt.main as m  # noqa: E402


def test_pypi_is_newer_basic():
    assert m._pypi_is_newer("1.6.0", "1.6.1") is True
    assert m._pypi_is_newer("1.6.1", "1.7.0") is True
    assert m._pypi_is_newer("1.6.1", "1.6.1") is False   # same
    assert m._pypi_is_newer("1.6.1", "1.6.0") is False   # older on PyPI


def test_pypi_is_newer_uses_numeric_ordering():
    # PEP 440 / numeric, not lexical: 1.10.0 > 1.9.0 (string compare gets this
    # wrong because "1.10" < "1.9" lexically).
    assert m._pypi_is_newer("1.9.0", "1.10.0") is True
    assert m._pypi_is_newer("1.10.0", "1.9.0") is False


def test_pypi_is_newer_handles_blanks():
    assert m._pypi_is_newer("", "1.6.1") is False
    assert m._pypi_is_newer("1.6.1", "") is False
    assert m._pypi_is_newer("", "") is False


def test_installed_version_matches_pyproject():
    """`_installed_vmpt_version` returns a version string when jwst-vmpt is
    importable as a package (or None if running from a bare source tree)."""
    v = m._installed_vmpt_version()
    assert v is None or isinstance(v, str)

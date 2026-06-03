"""Sphinx config for the vMPT documentation site.

Hosted at https://vmpt.readthedocs.io/. Built by Read the Docs from
``.readthedocs.yaml`` in the repo root; build settings (Python
version, OS, install steps) live there.

The site uses Markdown source via ``myst-parser`` so existing
top-level ``README.md`` and ``CHANGELOG.md`` can be ``{include}``-ed
verbatim — no need to maintain a separate RST copy.
"""

from __future__ import annotations

import importlib.metadata as _md
import os
import sys
from unittest.mock import MagicMock

# Make the in-tree ``vmpt`` package importable so autodoc can read
# the live source rather than the (possibly stale) installed wheel.
sys.path.insert(0, os.path.abspath(".."))


# `vmpt.coords` does
#     _msa_ap = pysiaf.Siaf("NIRSpec")["NRS_FULL_MSA"]
#     MSA_V2_REF = float(_msa_ap.V2Ref)
# at module import time. Under plain ``autodoc_mock_imports``,
# ``_msa_ap.V2Ref`` is itself a MagicMock — `float()` then raises
# TypeError and autodoc fails to import any module that transitively
# pulls in `vmpt.coords`. We inject a tighter pysiaf mock that
# returns the JWST-PRD fiducial numbers for the three attributes
# coords.py reads. Production code stays untouched.
class _SiafApertureMock:
    V2Ref = 378.842
    V3Ref = -428.575
    V3IdlYAngle = 138.5745697


class _SiafInstanceMock:
    def __getitem__(self, _key):
        return _SiafApertureMock()


_pysiaf_mock = MagicMock()
_pysiaf_mock.Siaf = MagicMock(return_value=_SiafInstanceMock())
sys.modules["pysiaf"] = _pysiaf_mock


# -- Project information -----------------------------------------------

project = "vMPT"
author = "Fengwu Sun"
copyright = "2026, Fengwu Sun"

# Pull the version from the installed package so RTD never gets out
# of sync with pyproject.toml.
try:
    release = _md.version("jwst-vmpt")
except _md.PackageNotFoundError:
    release = "0.0.0+unknown"
version = release


# -- General configuration ---------------------------------------------

extensions = [
    "myst_parser",          # Markdown → Sphinx (lets us reuse README.md)
    "sphinx.ext.autodoc",   # autogenerate from docstrings
    "sphinx.ext.napoleon",  # Google / NumPy-style docstrings
    "sphinx.ext.viewcode",  # "view source" links next to each symbol
    "sphinx.ext.intersphinx",
    "sphinx_design",        # `{grid}` / `{grid-item-card}` shorthand
]

# Cross-reference standard libraries used in our docstrings.
intersphinx_mapping = {
    "python":  ("https://docs.python.org/3", None),
    "numpy":   ("https://numpy.org/doc/stable/", None),
    "scipy":   ("https://docs.scipy.org/doc/scipy/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
}

# Files we accept as source. Lets us mix .md and .rst freely.
source_suffix = {
    ".md":  "markdown",
    ".rst": "restructuredtext",
}

master_doc = "index"

# Skip the build cruft.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Some of vMPT's deps are heavy / native (jwst, jwst_gtvt, pysiaf)
# and may not install cleanly inside RTD's container. Mocking them
# lets autodoc still import vmpt's source modules and pull
# docstrings even when the imports themselves can't resolve.
autodoc_mock_imports = [
    "astropy",
    "bokeh",
    "jwst",
    "jwst_gtvt",
    "pandas",
    "PIL",
    # `pysiaf` intentionally NOT mocked here — we provide a tighter
    # custom mock above so `vmpt.coords`'s float() conversions work.
    "scipy",
]

# MyST options. `colon_fence` lets us write
#   :::{include} ../README.md
# instead of needing the triple-backtick variant.
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "linkify",
    "smartquotes",
    "substitution",
]
myst_heading_anchors = 3  # auto-id H1..H3 so cross-page links work


# -- HTML output -------------------------------------------------------

html_theme = "furo"
html_static_path = ["_static"]
html_title = f"vMPT {release}"
html_logo = "../docs/logo.svg"
html_favicon = "../vmpt/static/favicon.svg"

# Furo-specific knobs. The blue accent matches the canvas
# "tip card" colour family used in the live app, so the docs site
# and the running tool feel like the same product.
html_theme_options = {
    "source_repository": "https://github.com/fengwusun/vMPT/",
    "source_branch":     "main",
    "source_directory":  "docs/",
    "light_css_variables": {
        "color-brand-primary":  "#1a3b66",
        "color-brand-content":  "#1a3b66",
    },
    "dark_css_variables": {
        "color-brand-primary":  "#5db0ff",
        "color-brand-content":  "#5db0ff",
    },
}

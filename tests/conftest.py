"""Shared pytest setup for the vMPT suite.

Disable the startup CRDS operability auto-update during tests: importing
``vmpt.main`` runs ``ensure_current_operability()``, which would otherwise
hit the network (and add latency) on the first import. Tests use whatever
operability reference is already in the local CRDS cache.
"""
import os

os.environ.setdefault("VMPT_OPERABILITY_AUTOUPDATE", "0")

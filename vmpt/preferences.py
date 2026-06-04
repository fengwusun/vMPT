"""Per-user display preferences for vMPT (v1.3.4+).

Saved to ``~/.vmpt/preferences.json`` so changes the user makes in
the Settings tab (canvas size, overlay alpha + stroke per layer,
slitlet size, layer visibility, stats-bar / catalog-hover ordering,
help-panel visibility) survive across sessions of the Bokeh app.

The file is plain JSON with a flat dict at the top level. Missing
keys are filled in from hard-coded defaults on load, so an
upgrade that adds a new persisted setting doesn't crash when an
older prefs file is read. Conversely, unknown keys are ignored
silently — making it easy to roll back without manual cleanup.

The module deliberately has NO Bokeh dependency. ``vmpt/main.py``
calls :func:`load_preferences` once at startup and
:func:`save_preferences` on every widget change. There's also a
``reset_preferences`` helper that wipes the file entirely (the
"Reset to defaults" button calls this then re-applies hard-coded
defaults to every widget).
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any


_PREFS_DIR_NAME = ".vmpt"
_PREFS_FILE_NAME = "preferences.json"


def _resolve_prefs_path() -> Path:
    """Return the absolute path of the preferences file.

    Defaults to ``~/.vmpt/preferences.json``. The home directory
    is resolved via ``Path.home()`` which honours the standard
    ``$HOME`` env var. Tests / CI can override via the
    ``VMPT_PREFS_PATH`` env var.
    """
    override = os.environ.get("VMPT_PREFS_PATH")
    if override:
        return Path(override).expanduser()
    return Path.home() / _PREFS_DIR_NAME / _PREFS_FILE_NAME


def load_preferences() -> dict[str, Any]:
    """Read the user's preferences file.

    Returns an empty dict when the file doesn't exist, when it's
    unreadable, or when it contains malformed JSON. Never raises —
    a corrupt prefs file shouldn't break app startup.
    """
    path = _resolve_prefs_path()
    if not path.exists():
        return {}
    try:
        raw = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return {}
    if not raw.strip():
        return {}
    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        return {}
    if not isinstance(data, dict):
        return {}
    return data


def save_preferences(prefs: dict[str, Any]) -> None:
    """Write the preferences dict to disk.

    Creates the parent directory if it doesn't exist. Writes via a
    temp file + rename so we never leave a half-written prefs file
    if the process is killed mid-write. Silently swallows
    OS errors — if the prefs directory isn't writable the user
    just loses their preferences across sessions; nothing else
    breaks.
    """
    path = _resolve_prefs_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
    except OSError:
        return
    tmp = path.with_suffix(path.suffix + ".tmp")
    try:
        tmp.write_text(
            json.dumps(prefs, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        os.replace(tmp, path)
    except OSError:
        # Best effort: try to clean up the half-written tmp file.
        try:
            tmp.unlink()
        except OSError:
            pass


def reset_preferences() -> None:
    """Remove the preferences file. Idempotent — no-op when the
    file isn't there. Never raises."""
    path = _resolve_prefs_path()
    try:
        path.unlink()
    except FileNotFoundError:
        return
    except OSError:
        return

"""Regression: late-wired button handlers must actually receive clicks.

`reset_prefs_btn.on_click(_on_reset_prefs)` is wired near the very end of
`vmpt.main` — long after the Settings-tab layout root was added to the
document. Bokeh only enrols a model in its document's event dispatcher
(`Document.callbacks._subscribed_models`) at *attach* time (via
`_attach_document` -> `_update_event_callbacks`). A handler wired
post-attach is recorded on the model's `_event_callbacks` but the model is
never subscribed, so the Bokeh server silently drops its ButtonClick events
and the button "does nothing" when clicked.

`vmpt.main` calls `_resubscribe_late_event_handlers()` right after such late
wiring to repair this. These tests guard both the helper (in isolation) and
the actual app buttons that depend on it.
"""
from __future__ import annotations

import sys
from pathlib import Path

from bokeh.document import Document
from bokeh.events import ButtonClick
from bokeh.models import Button, Column

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import vmpt.main as m  # noqa: E402


def _button_click_subscribers(doc: Document) -> set[str]:
    refs = doc.callbacks._subscribed_models.get("button_click", set())
    return {r().id for r in refs if r() is not None}


def test_resubscribe_helper_enrols_post_attach_handler():
    """Isolated bug -> fix: a handler wired after the model is attached is
    not subscribed; the helper enrols it so the server event path reaches
    it."""
    doc = Document()
    btn = Button(label="x")
    doc.add_root(Column(children=[btn]))      # attach FIRST...
    fired: list[int] = []
    btn.on_click(lambda: fired.append(1))     # ...then wire (post-attach)

    # Reproduce the bug: the model is not in the document's dispatcher.
    assert btn.id not in _button_click_subscribers(doc)

    m._resubscribe_late_event_handlers(btn)   # the repair the app applies

    assert btn.id in _button_click_subscribers(doc)
    # And the exact event the Bokeh server raises on a click now lands.
    doc.callbacks.trigger_event(ButtonClick(model=btn))
    assert fired == [1]


def test_resubscribe_helper_is_safe_when_unattached():
    """No document yet -> no-op, no crash (the eventual attach subscribes)."""
    orphan = Button(label="orphan")
    orphan.on_click(lambda: None)
    m._resubscribe_late_event_handlers(orphan)   # must not raise
    assert orphan.document is None


def test_reset_and_help_buttons_are_subscribed_in_app():
    """The real app's late-wired buttons are enrolled in the dispatcher.

    Before the fix `reset_prefs_btn` was absent here (its only handler was
    wired post-attach), so its clicks never reached `_on_reset_prefs`.
    """
    doc = m.reset_prefs_btn.document
    assert doc is not None, "reset button should be attached to the document"

    subs = _button_click_subscribers(doc)
    assert m.reset_prefs_btn.id in subs, "Reset-display button not subscribed"
    assert m.help_toggle_btn.id in subs, "Help-toggle button not subscribed"

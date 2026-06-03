"""Console entry point for the ``vmpt`` command.

Two subcommands:

  vmpt                    Start the Bokeh server on the bundled app.
  vmpt examples download  Fetch example_a370 / example_r0600 from
                          a GitHub release tarball into the current
                          working directory.

The default subcommand forwards the same flags as ``./run.sh``:
``--port``, ``--fits``, ``--jpg``, ``--wcs``, ``--catalog``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

# Where the in-package examples-bundle is fetched from when the user
# runs ``vmpt examples download``. The GitHub release page hosts a
# tarball asset of the same name; ``v1.2.1`` is set explicitly so
# upgrading vMPT doesn't invalidate cached example downloads.
EXAMPLES_RELEASE_TAG = "v1.3.0"
EXAMPLES_TARBALL_URL = (
    "https://github.com/fengwusun/vMPT/releases/download/"
    f"{EXAMPLES_RELEASE_TAG}/vmpt-examples.tar.gz"
)


def _app_dir() -> str:
    """Absolute path to the Bokeh app directory inside the installed
    package. Bokeh serves *this* directory; it contains main.py."""
    return os.path.dirname(os.path.abspath(__file__))


def _serve(argv: list[str]) -> int:
    """Forward to ``bokeh serve <pkg> --show`` with vMPT's tuned
    defaults (large WebSocket frame so big uploads don't truncate)."""
    # Split flags into "consumed by bokeh" and "forwarded to the app
    # via --args". `--port` is for bokeh; everything else (--fits,
    # --jpg, --wcs, --catalog) is for vmpt/main.py's autoload.
    bokeh_flags: list[str] = []
    app_args: list[str] = []
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok == "--port" and i + 1 < len(argv):
            bokeh_flags += [tok, argv[i + 1]]
            i += 2
        elif tok in ("--fits", "--jpg", "--wcs", "--catalog") and i + 1 < len(argv):
            app_args += [tok, argv[i + 1]]
            i += 2
        elif tok in ("-h", "--help"):
            _print_help()
            return 0
        else:
            # Unknown — forward to bokeh so e.g. --log-level Trace works.
            bokeh_flags.append(tok)
            i += 1

    cmd: list[str] = [
        "bokeh", "serve", _app_dir(),
        "--websocket-max-message-size", "524288000",
        "--show",
        *bokeh_flags,
    ]
    if app_args:
        cmd += ["--args", *app_args]

    # Defer to the system `bokeh` rather than calling the bootstrap
    # entry point — that way Ctrl+C cleanly terminates the server and
    # we don't need to wrap signal handling here.
    try:
        return subprocess.call(cmd)
    except FileNotFoundError:
        sys.stderr.write(
            "error: `bokeh` not found on PATH. "
            "Install with `pip install bokeh` "
            "or re-install vMPT into an active environment.\n"
        )
        return 2


def _examples_download(target_dir: str | None = None) -> int:
    """Fetch the example_a370 / example_r0600 tarball from the
    matching GitHub release into the current directory (or an
    explicit ``target_dir`` if given)."""
    dest = Path(target_dir or os.getcwd()).resolve()
    dest.mkdir(parents=True, exist_ok=True)
    print(f"vmpt: downloading examples for {EXAMPLES_RELEASE_TAG}…")
    print(f"      source: {EXAMPLES_TARBALL_URL}")
    print(f"      dest:   {dest}")

    with tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        try:
            urllib.request.urlretrieve(EXAMPLES_TARBALL_URL, tmp_path)
        except Exception as e:  # network / 404 / …
            sys.stderr.write(
                f"error: download failed: {e}\n"
                f"If you're running a pre-release build, the asset "
                f"may not exist yet — grab the examples from the "
                f"vMPT repo's example_a370/ and example_r0600/ "
                f"directories instead.\n"
            )
            return 3
        with tarfile.open(tmp_path, "r:gz") as tf:
            tf.extractall(dest)
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
    print("vmpt: examples downloaded. Try `vmpt --jpg "
          "example_r0600/JWST_F090W_F200W_F444W.jpg "
          "--wcs example_r0600/wcs.fits "
          "--catalog example_r0600/v01_fsun.cat`.")
    return 0


def _print_help() -> None:
    print(
        "Usage:\n"
        "  vmpt [--port N] [--fits PATH] [--jpg PATH --wcs PATH] \\\n"
        "       [--catalog PATH]…\n"
        "  vmpt examples download [DIR]\n"
        "  vmpt --help\n"
        "\n"
        "Start the vMPT Bokeh server (default) or fetch the example\n"
        "datasets. The serve subcommand accepts the same flags as\n"
        "./run.sh; --jpg requires --wcs, --fits is mutually exclusive\n"
        "with --jpg/--wcs, and --catalog may be repeated to stack\n"
        "multiple catalogs."
    )


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "examples":
        sub = argv[1] if len(argv) > 1 else ""
        if sub == "download":
            return _examples_download(argv[2] if len(argv) > 2 else None)
        sys.stderr.write(
            "error: unknown examples subcommand. "
            "Try `vmpt examples download [DIR]`.\n"
        )
        return 2
    return _serve(argv)


if __name__ == "__main__":
    sys.exit(main())

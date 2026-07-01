#!/bin/bash
# Start the NIRSpec MSA planner.
# Bumps WebSocket message size to 500 MB so file uploads don't truncate.
# For routine use, prefer the "path" text inputs in the sidebar.
#
# Optional arguments:
#   --port N          Bokeh server port (default 5006).
#   --fits PATH       Open this FITS image (mutually exclusive with --jpg/--wcs).
#   --jpg PATH        Open this JPG/PNG (must be paired with --wcs).
#   --wcs PATH        WCS sidecar FITS for the JPG (must be paired with --jpg).
#   --catalog PATH    Load this catalog file (CSV / ASCII / FITS).
#                     May be repeated to layer multiple catalogs.
#   --addon PATH      Load a DS9 region (.reg) / contour (.ctr/.con) overlay.
#                     May be repeated.
#   --v3pa DEG        Initial V3 PA in degrees (e.g. --v3pa 209).
#   --apa DEG         Initial NIRSpec aperture PA (APA = V3 PA + V3IdlYAngle);
#                     --v3pa wins if both are given.
#
# Examples:
#   ./run.sh
#   ./run.sh --port 5010
#   ./run.sh --fits example_a370/a370_f182m_f200w_f210m.fits
#   ./run.sh --jpg image.jpg --wcs wcs.fits --catalog targets.csv --v3pa 209
#   ./run.sh --catalog a.csv --catalog b.csv

# Conda's activate hooks reference unset env vars (e.g. GFORTRAN on
# macOS), so we activate FIRST without `-u`, then turn on strict mode
# for the rest of the script.
set -eo pipefail

cd "$(dirname "$0")"

# Activate stenv if not already in it
if [[ "${CONDA_DEFAULT_ENV:-}" != "stenv" ]]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
    conda activate stenv
fi

set -u

PORT=5006
FITS_PATH=""
JPG_PATH=""
WCS_PATH=""
CATALOG_PATHS=()
ADDON_PATHS=()
V3PA=""
APA=""

usage() {
    cat <<'EOF' >&2
Usage: ./run.sh [--port N] [--fits PATH] [--jpg PATH --wcs PATH] [--catalog PATH]...
  --port N          Bokeh server port (default 5006).
  --fits PATH       FITS image (with WCS in its header).
  --jpg PATH        JPG/PNG image; REQUIRES --wcs alongside.
  --wcs PATH        Sidecar FITS holding the WCS for --jpg.
  --catalog PATH    Catalog (CSV / ASCII / FITS). Can be repeated.
  --addon PATH      DS9 region (.reg) / contour (.ctr/.con) overlay. Repeatable.
  --v3pa DEG        Initial V3 PA in degrees (e.g. --v3pa 209).
  --apa DEG         Initial NIRSpec aperture PA (APA = V3 PA + V3IdlYAngle).
                    --v3pa wins if both are given.
  -h, --help        Show this message.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --port)
            [[ $# -ge 2 ]] || { echo "run.sh: --port requires a number" >&2; exit 2; }
            PORT="$2"; shift 2 ;;
        --fits)
            [[ $# -ge 2 ]] || { echo "run.sh: --fits requires a path" >&2; exit 2; }
            FITS_PATH="$2"; shift 2 ;;
        --jpg)
            [[ $# -ge 2 ]] || { echo "run.sh: --jpg requires a path" >&2; exit 2; }
            JPG_PATH="$2"; shift 2 ;;
        --wcs)
            [[ $# -ge 2 ]] || { echo "run.sh: --wcs requires a path" >&2; exit 2; }
            WCS_PATH="$2"; shift 2 ;;
        --catalog)
            [[ $# -ge 2 ]] || { echo "run.sh: --catalog requires a path" >&2; exit 2; }
            CATALOG_PATHS+=("$2"); shift 2 ;;
        --addon)
            [[ $# -ge 2 ]] || { echo "run.sh: --addon requires a path" >&2; exit 2; }
            ADDON_PATHS+=("$2"); shift 2 ;;
        --v3pa)
            [[ $# -ge 2 ]] || { echo "run.sh: --v3pa requires a number (deg)" >&2; exit 2; }
            V3PA="$2"; shift 2 ;;
        --apa)
            [[ $# -ge 2 ]] || { echo "run.sh: --apa requires a number (deg)" >&2; exit 2; }
            APA="$2"; shift 2 ;;
        -h|--help)
            usage; exit 0 ;;
        --)
            shift; break ;;
        *)
            echo "run.sh: unknown argument $1" >&2
            usage; exit 2 ;;
    esac
done

# Sanity-check the port.
if ! [[ "$PORT" =~ ^[0-9]+$ ]] || (( PORT < 1 || PORT > 65535 )); then
    echo "run.sh: --port must be an integer in [1, 65535], got '$PORT'" >&2
    exit 2
fi

# Validate combinations.
if [[ -n "$JPG_PATH" && -z "$WCS_PATH" ]] || \
   [[ -z "$JPG_PATH" && -n "$WCS_PATH" ]]; then
    echo "run.sh: --jpg and --wcs must be specified together." >&2
    exit 2
fi
if [[ -n "$FITS_PATH" && (-n "$JPG_PATH" || -n "$WCS_PATH") ]]; then
    echo "run.sh: --fits is mutually exclusive with --jpg / --wcs." >&2
    exit 2
fi

# Resolve paths to absolute so the app doesn't have to guess the cwd.
abspath() { python3 -c 'import os, sys; print(os.path.abspath(sys.argv[1]))' "$1"; }

[[ -n "$FITS_PATH" ]] && FITS_PATH=$(abspath "$FITS_PATH")
[[ -n "$JPG_PATH" ]]  && JPG_PATH=$(abspath "$JPG_PATH")
[[ -n "$WCS_PATH" ]]  && WCS_PATH=$(abspath "$WCS_PATH")
RESOLVED_CATALOGS=()
for c in "${CATALOG_PATHS[@]:-}"; do
    [[ -n "$c" ]] && RESOLVED_CATALOGS+=("$(abspath "$c")")
done
RESOLVED_ADDONS=()
for a in "${ADDON_PATHS[@]:-}"; do
    [[ -n "$a" ]] && RESOLVED_ADDONS+=("$(abspath "$a")")
done

# Build the `--args` list for Bokeh. Each --args value lands in
# `sys.argv` inside vmpt/main.py at startup.
APP_ARGS=()
[[ -n "$FITS_PATH" ]]   && APP_ARGS+=(--fits "$FITS_PATH")
[[ -n "$JPG_PATH" ]]    && APP_ARGS+=(--jpg "$JPG_PATH")
[[ -n "$WCS_PATH" ]]    && APP_ARGS+=(--wcs "$WCS_PATH")
# `${arr[@]:-}` yields one empty-string element for an empty array (needed so
# `set -u` doesn't trip on bash 3.2 / macOS) — so guard on non-empty to avoid
# appending a spurious `--catalog ""` / `--addon ""`.
for c in "${RESOLVED_CATALOGS[@]:-}"; do
    [[ -n "$c" ]] && APP_ARGS+=(--catalog "$c")
done
for a in "${RESOLVED_ADDONS[@]:-}"; do
    [[ -n "$a" ]] && APP_ARGS+=(--addon "$a")
done
[[ -n "$V3PA" ]] && APP_ARGS+=(--v3pa "$V3PA")
[[ -n "$APA" ]]  && APP_ARGS+=(--apa "$APA")

CMD=(bokeh serve vmpt/
     --port "$PORT"
     --websocket-max-message-size 524288000
     --show)
if [[ ${#APP_ARGS[@]} -gt 0 ]]; then
    CMD+=(--args "${APP_ARGS[@]}")
fi

exec "${CMD[@]}"

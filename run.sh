#!/bin/bash
# Start the NIRSpec MSA planner.
# Bumps WebSocket message size to 500 MB so file uploads don't truncate.
# For routine use, prefer the "path" text inputs in the sidebar.

set -euo pipefail

cd "$(dirname "$0")"

# Activate stenv if not already in it
if [[ "${CONDA_DEFAULT_ENV:-}" != "stenv" ]]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
    conda activate stenv
fi

exec bokeh serve app/ \
    --port 5006 \
    --websocket-max-message-size 524288000 \
    --show

#!/usr/bin/env bash
# Time- and x/z-averaged wall-normal profiles -> z.openfoam_output/openfoam_averaged.csv
# Settings (fields, averaging window 100-400, directions) are at the top of paraview_openfoam.py.
set -euo pipefail
cd "${0%/*}"

# Path to ParaView's pvpython: override with   PV_PYTHON=/path/to/pvpython ./b.post-processing.sh
PV_PYTHON="${PV_PYTHON:-$(command -v pvpython || echo /shared/ParaView/bin/pvpython)}"

if [[ ! -x "$PV_PYTHON" ]]; then
  echo "ERROR: pvpython not found: $PV_PYTHON (set PV_PYTHON=/path/to/pvpython)" >&2
  exit 1
fi

[[ -f foam.foam ]] || : > foam.foam
echo "Running in: $(pwd)"
exec "$PV_PYTHON" paraview_openfoam.py

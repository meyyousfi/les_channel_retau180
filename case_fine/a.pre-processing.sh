#!/bin/bash
# Pre-processing for the FINE (144^3) mesh.
# Initial condition: the final (latest) solution of case_medium, interpolated onto this
# mesh with mapFields. case_medium must therefore have been run AND reconstructed first.
cd "${0%/*}" || exit 1
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions

SOURCE_CASE=../case_medium

# Latest reconstructed time directory in the source case (ignores 0 and 0.orig)
latest=$(ls -1 "$SOURCE_CASE" 2>/dev/null | grep -E '^[0-9]+(\.[0-9]+)?(e[+-]?[0-9]+)?$' \
         | grep -vx '0' | sort -g | tail -1)
if [ -z "$latest" ]; then
    echo "ERROR: no reconstructed results found in $SOURCE_CASE."
    echo "       Run that case first, then run 'reconstructPar' inside it."
    exit 1
fi
echo "Mapping initial fields from $SOURCE_CASE, time $latest"

restore0Dir                      # 0.orig -> 0 (placeholder fields with correct BCs)
runApplication blockMesh
tail log.blockMesh
runApplication checkMesh
runApplication mapFields "$SOURCE_CASE" -sourceTime "$latest" -consistent
tail log.mapFields
runApplication decomposePar
tail log.decomposePar

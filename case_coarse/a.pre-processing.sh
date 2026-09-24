#!/bin/bash
# Pre-processing for the COARSE mesh (64^3).
# Initial condition: synthetic turbulent field stored in 0.orig/initial*Profile
# (generated with tools/matlab/turbulent_field_generator.m for a 64^3 mesh).
cd "${0%/*}" || exit 1
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions

restore0Dir                      # copy 0.orig -> 0
runApplication blockMesh
tail log.blockMesh
runApplication checkMesh
runApplication decomposePar
tail log.decomposePar

#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory
set -x

export EXTBLOCKMESH_CODE=$PWD

wmake MeshSmoother
wmake
wmake hexMeshSmoother

# ----------------------------------------------------------------- end	of-file

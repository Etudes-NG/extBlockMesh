#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory
set -x

export EXTBLOCKMESH_CODE=$PWD

echo $EXTBLOCKMESH_CODE

wmake blockMeshSmoother
wmake


# ----------------------------------------------------------------- end-of-file

#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory
set -x

wclean blockMeshSmoother
wclean


# ----------------------------------------------------------------- end-of-file

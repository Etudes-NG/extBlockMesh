#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory
set -x

wclean MeshSmoother
wclean


# ---------------------------------------------------------------- end-of-file

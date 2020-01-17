#!/bin/bash

set -e
set -u
set -x

# SLEPc also does not support out-of-tree builds. It would actually have
# surprised me if that library did.
builddir="/qbigwork2/$USER/slepc"
pushd "$builddir"

# Although PETSc has a certain `maint` branch which has the releases, SLEPc
# does not. So we need to check out some specific version within git.
git co v3.12.1

export CC=mpicc
export CFLAGS="-O3 -march=sandybridge"

export PETSC_ARCH=idontcare
export PETSC_DIR="/qbigwork2/$USER/petsc"
export SLEPC_DIR="/qbigwork2/$USER/slepc"

/qbigwork2/$USER/slepc/configure \
    --prefix="/qbigwork2/$USER/local"

make -j $(nproc)

make install

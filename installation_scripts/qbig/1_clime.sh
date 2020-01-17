#!/bin/bash

set -e
set -u
set -x

mkdir -p "/qbigwork2/$USER/local"

builddir="/qbigwork2/$USER/Build/clime"
mkdir -p "$builddir"
pushd "$builddir"

export CC=mpicc
export CFLAGS="-O3 -march=sandybridge"

/qbigwork2/$USER/c-lime/configure \
    --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
    --prefix="/qbigwork2/$USER/local"

make -j $(nproc)

make install

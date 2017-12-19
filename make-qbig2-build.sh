#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
# Licensed under The MIT/Expat License

# Sets up a debug and a release build outside of the source code with the paths
# that are needed on QBiG. The debug build will automatically set the `-g` flag
# for debugging, the release version will be performance and has `-O2 -DNDEBUG`
# set.

set -e
set -u
set -x

sourcedir=$(pwd)

for buildtype in release; do
    builddir=../LapH_EigSys-$buildtype
    rm -rf "$builddir"
    mkdir "$builddir"
    pushd "$builddir"

    cmake \
        "$sourcedir" \
        -DEIGEN3_INCLUDE_DIRS='/hadron/helmes/libraries/eigen-3.3.4' \
        -DPETSC_INCLUDE_DIRS="${PETSC_DIR}/include;${PETSC_DIR}/${PETSC_ARCH}/include"\
        -DPETSC_LIBRARIES="-L ${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc" \
        -DSLEPC_INCLUDE_DIRS="${SLEPC_DIR}/include;${SLEPC_DIR}/${SLEPC_ARCH}/include" \
        -DSLEPC_LIBRARIES="-L ${SLEPC_DIR}/${SLEPC_ARCH}/lib -lslepc"\
        -DLIME_INCLUDE_DIRS='/hadron/helmes/libraries/lime-1.3.2/include' \
        -DLIME_LIBRARIES='-L /hadron/helmes/libraries/lime-1.3.2/lib -llime' \
        -DCMAKE_CXX_COMPILER='mpic++' \
        -DCMAKE_BUILD_TYPE=$buildtype

    make -j $(nproc)
    popd
done

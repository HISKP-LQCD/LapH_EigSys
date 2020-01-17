#!/bin/bash

# The LapH eigensystem code has many different options for building:
#
# - Use the `CMakeLists.txt` and call `cmake` yourself.
# - Use the `make-qbig2-build.sh` which sounds about right. One just needs to
#   add a little extra stuff.
# - Follow the documentation and use the `src/Makefile` for compilation.
# - You would also use `src/Makefile_jureca` …
# - … or `Makefile_qbig` if you want to.

set -e
set -u
set -x

builddir="/qbigwork2/$USER/Build/laph_eigsys"
sourcedir="/qbigwork2/$USER//LapH_EigSys"
mkdir -p "$builddir"
pushd "$builddir"

export CC=mpicc
export CXX=mpicxx
export CFLAGS="-O3 -march=sandybridge"
export CXXFLAGS='-O3 -march=sandybridge'

export PETSC_ARCH=idontcare
export PETSC_DIR="/qbigwork2/$USER/petsc"
export SLEPC_DIR="/qbigwork2/$USER/slepc"

# We can't use the `make-qbig2-build.sh` because it has various hard-coded
# paths. Well shit, so we just copy the bulk from there and call `cmake`
# ourselves.
cmake \
    "$sourcedir" \
    -DPETSC_INCLUDE_DIRS="${PETSC_DIR}/include;${PETSC_DIR}/${PETSC_ARCH}/include"\
    -DPETSC_LIBRARIES="-L ${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc" \
    -DSLEPC_INCLUDE_DIRS="${SLEPC_DIR}/include;${SLEPC_DIR}/${PETSC_ARCH}/include" \
    -DSLEPC_LIBRARIES="-L ${SLEPC_DIR}/${PETSC_ARCH}/lib -lslepc"\
    -DLIME_INCLUDE_DIRS="/qbigwork2/$USER/local/include" \
    -DLIME_LIBRARIES="-L /qbigwork2/$USER/local/lib -llime" \
    -DCMAKE_CXX_COMPILER='mpicxx' \
    -DCMAKE_BUILD_TYPE=Release

make -j $(nproc)

make install

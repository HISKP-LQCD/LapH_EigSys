#!/bin/bash

# PETSc is a fun library to compile. You see the `configure` file in the
# repository and think “oh, it uses GNU Autotools” and then think that there
# should not be a `configure` shipped in the git repo but rather a
# `configure.in` and `configure.ac` such that `autoconf` can generate it. But
# if you take a look, it is Python script. Yeah, exotic build system, this is
# going to be a blast!
#
# In order to find out which options it can take, use `./configure --help |
# less`, the output is sheer endless. That is another sure indicator that the
# project is overly complex.

set -e
set -u
set -x

builddir="/qbigwork2/$USER/petsc"
pushd "$builddir"

# Then there is this `petsc-arch` thingy. Usually you just have Debug and
# Release builds, but of course that would be overly simplistic for a
# everything-and-the-kitchen-sink library. Usually a library supports out of
# tree builds and then one can configure it differently in the various build
# directories. But for PETSc you cannot call the `configure` script from
# outside the directory because it loads some other Python modules which are
# not in the `PYTHON_PATH`. They have these “build architectures” which you can
# apparently freely name and then it will create a subdirectory. I call that
# bullshit.
export PETSC_ARCH=idontcare

# And for some reason it needs to have its source tree directory as an
# environment variable because the working directory is not enough. If you do
# not specify this, then the `make install` in the end will not do anything.
export PETSC_DIR="$builddir"

# It seems to be a C library but you can also compile it using a C++ compiler.
# I don't even want to know why this is a special thing and one does not just
# switch out `gcc` for `g++`, but at this point I don't even care.
#
# The configure uses booleans but their values are 0 and 1. Fucking weak C
# types.
/qbigwork2/$USER/petsc/configure \
    --with-clanguage=C++ \
    --with-64-bit-indices=1 \
    --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
    --with-debugging=0 \
    COPTFLAGS='-O3 -march=sandybridge' \
    CXXOPTFLAGS='-O3 -march=sandybridge' \
    FOPTFLAGS='-O3 -march=sandybridge' \
    --prefix="/qbigwork2/$USER/local"

make -j $(nproc) all

make install

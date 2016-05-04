LapH_EigSys
===================

Serial version of the eigensystem calculation in the stochastic Laplacian Heaviside Smearing

Dependencies: - c-lime by USQCD (http://usqcd-software.github.io/c-lime/)
              - Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
              - PETSc by the Argonne National Laboratory (http://www.mcs.anl.gov/petsc/)
              - SLEPc an extension to PETSc under maintenance of the Univeristat
                Politecnica de Valencia (http://slepc.upv.es/)

Preparations:
-------------
A PETSc example configure script can be found in install/conf_scripts. For the
package to work complex scalars, lapack and blas are needed. The C language should be C++ and the g++ compiler and
linker should be of version later than 4.6.

Compiling:
----------
Once the dependencies are built the
Makefile in ./src can be invoked to compile the executable ev_ts. The include of
mpiuni avoids certain issues with PETSc looking for mpi headers which are not
needed. 
Cleaning up works by calling the target

make myclean

Execution:
----------
For execution make sure parameters.txt is adjusted accordingly and gauge
configurations are at hand. Normal execution is done by invoking

./ev_ts =CONF=

where =CONF= should be a valid configuration number. The according
parameters.txt is read from the execution directory. The program output contains
small sanity checks as the $\tr(V^dag * V}$ that should equal the number of eigenvectors.
Eigenvectors, eigenvalues and phases are saved in binary format where the eigenvectors
are phasefixed. 

MPI Compatibility:
-------------------
It is now possible to use a version parallelized with MPI. Therefore PETSc and SLEPc have to be compiled with MPI-Support. In the programm itself every rank has an own EPS-Instance, ensured via the usage of PETSC_COMM_SELF instead of PETSC_COMM_WORLD. At the moment the ranks play as well the role of choosing the timeslice, therefore MPI_Size() has to be the time extent of the lattice

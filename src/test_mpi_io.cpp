#include <cstdlib>
#include <string>
#include <mpi.h>
#include "read_write_mpi.h"

int main(int argc, char **argv) {
  int mpistat = 0;
  MPI::Init();
  int mpi_size, mpi_rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  MpiIO save_eigsys("test.h5",mpi_rank);

  //save_eigsys.write_rec();
  save_eigsys.finalize(mpi_rank);
  MPI::Finalize();
  return 0;
}


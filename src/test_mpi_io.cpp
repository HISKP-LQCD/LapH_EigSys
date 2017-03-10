#include <cstdlib>
#include <iostream>
#include <string>
#include <Eigen/Core>
#include <mpi.h>
#include "par_io.h"
#include "read_write_mpi.h"

int main(int argc, char **argv) {
  int mpistat = 0;
  MPI::Init();
  int mpi_size = 0, mpi_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  //get number of configuration from last argument to main
  int config = atoi( argv[ (argc-1) ] );
  --argc;

  // handle input file
  IO* pars = IO::getInstance();
  pars->set_values("parameters.txt");
  if(mpi_rank == 0) {
    pars->print_summary();
    printf("calculating config %d\n", config);
  }

  // calculate the timeslices to work on
  // is done on every communicator because I don't
  // know how else to do
  const int L0 = pars -> get_int("LT");
  int tstart = 0;
  int tend = 0;
  int todo = 0;
  int *tstarts = new int[mpi_size];
  int *tends = new int[mpi_size];
  int *todos = new int[mpi_size];
  if (mpi_rank == 0) {
    // calculate chunks minimal chunk size
    int tmp = L0/mpi_size;
    // fill array
    for(int i = 0; i < mpi_size; ++i)
      todos[i] = tmp;
    // increase the first chunks until all
    // time slices are distributed
    for(int i = 0; i < L0%mpi_size; ++i)
      ++todos[i];
    // fill tstarts and tend;
    tstarts[0] = 0;
    tends[0] = todos[0] - 1;
    for(int i = 1; i < mpi_size; ++i) {
      tstarts[i] = tends[i-1] + 1;
      tends[i] = tstarts[i] + todos[i] - 1;
    }
  }
  //print information about every process and its associated timeslices
  for(int i = 0; i < mpi_size; ++i){
    std::cout << todos[i] << std::endl;
  }
  MpiIO save_eigsys("test.h5",mpi_rank);
  MPI_Scatter(tstarts, 1, MPI_INT, &tstart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(tends, 1, MPI_INT, &tend, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(todos, 1, MPI_INT, &todo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  delete [] tstarts;
  delete [] tends;
  delete [] todos;
  //loop over timeslices of a configuration
  for(int ts = 0; ts < todo; ++ts) {
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(pars->get_int("MAT_ENTRIES"),pars->get_int("NEV"));   
    save_eigsys.write_ds(ts+tstart,V);
    std::vector<double> eval (pars -> get_int("NEV"));
    for (auto& i:eval) i = 0;
    save_eigsys.write_ds(ts+tstart, "ev", eval);
    std::vector<double> phs (pars -> get_int("NEV"));
    for (auto& i:phs) i = 0;
    save_eigsys.write_ds(ts+tstart, "phs", phs);
  }
  save_eigsys.finalize(mpi_rank);
  MPI::Finalize();
  return 0;
}


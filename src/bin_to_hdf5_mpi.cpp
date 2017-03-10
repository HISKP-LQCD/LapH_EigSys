// Convert a read in Eigensystem to an hdf5 file with the following structure
//
// The input parameters should just be taken over from the parameters structure
// /meta/parameters.txt
// Some information about generation, a timestamp, a user id etc.
// /meta/descriptor
// the eigenvectors for one configuration each timeslice is one object, important
// for mpi?
// /data/evecs/evec_[000-096]
// ... same for the eigenvalues ...
// /data/evals/eval_[000-096]
// ... and for the phases
// /data/phases/phase_[000-096]
//
// Systemincludes
#include <cstdlib>
#include <iostream>
//#include <Eigen/Core>
#include <mpi.h>
#include <string>
#include <hdf5.h>

// own includes
#include "read_write.h"
#include "par_io.h"
int main(int argc, char **argv){
  //__Initialize MPI__
  int mpistat = 0;
  MPI_Init(&argc,&argv);
  int world = 0, rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Info info  = MPI_INFO_NULL;
  // Get configuration argument
  int conf = atoi( argv[ (argc-1) ] );
  --argc;
  // Get parameters for read in
  IO* pars = IO::getInstance();
  pars -> set_values("./parameters.txt");
  printf("calculating config %d\n", conf);

  const int L0 = pars -> get_int("LT");
  int tstart = 0;
  int tend = 0;
  int todo = 0;
  int *tstarts = new int[world];
  int *tends = new int[world];
  int *todos = new int[world];
  if (rank == 0) {
    // calculate chunks minimal chunk size
    int tmp = L0/world;
    // fill array
    for(int i = 0; i < world; ++i)
      todos[i] = tmp;
    // increase the first chunks until all
    // time slices are distributed
    for(int i = 0; i < L0%world; ++i)
      ++todos[i];
    // fill tstarts and tend;
    tstarts[0] = 0;
    tends[0] = todos[0] - 1;
    for(int i = 1; i < world; ++i) {
      tstarts[i] = tends[i-1] + 1;
      tends[i] = tstarts[i] + todos[i] - 1;
      std::cout << "Tstart: " << tstarts[i] << " Tend: " << tends[i] << std::endl;
    }
  }
  MPI_Scatter(tstarts, 1, MPI_INT, &tstart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(tends, 1, MPI_INT, &tend, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(todos, 1, MPI_INT, &todo, 1, MPI_INT, 0, MPI_COMM_WORLD);

  delete [] tstarts;
  delete [] tends;
  delete [] todos;
  //calculation parameters from infile
  const int MAT_ENTRIES = pars -> get_int("MAT_ENTRIES");
  const int nev = pars -> get_int("NEV");
  pars->print_summary();
  std::vector<double> eval;
  std::vector<double> phase;

  // Define a new compund type for complex data
  hid_t h5_comp = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<double>));
  herr_t err = H5Tinsert(h5_comp,"real",0,H5T_NATIVE_DOUBLE);
  err = H5Tinsert(h5_comp,"imag", sizeof(double), H5T_NATIVE_DOUBLE);

  // Typedef for rowmajor complex eigenmatrix
  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RM_MatrixXcd;
  // Rewrite this thing in standard C to facilitate MPI functionality
  // initialize hdf5
//std::string filename = "eigensys_"+std::to_string(conf)+".h5";
  char filename [150];
  sprintf(filename,"eigensys_%04d.h5",conf);
  /* 
   * Set up file access property list with parallel I/O access
   */
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
              H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

  /*
   * Create a new file collectively and release property list identifier.
   */
  hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
            H5Pclose(plist_id);
  
	
	// Create a group named "/MygGroup" in the file
  hid_t group_id_info = H5Gcreate(file_id, "/info", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t group_id_data = H5Gcreate(file_id, "/data", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t group_id_evecs = H5Gcreate(file_id, "/data/evecs", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t group_id_evals = H5Gcreate(file_id, "/data/evals", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t group_id_phase = H5Gcreate(file_id, "/data/phase", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);


  for(int ts = 0; ts < L0; ++ts) {
    Eigen::MatrixXcd esys(MAT_ENTRIES,nev);
    // eigenvectors
    // read
    read_evectors_bin_ts("eigenvectors", conf, ts, nev, esys);
	// only debug info
    std::cout << "entry (32453,50)" <<  esys(32453,50) << std::endl; 
    // The eigenvectors should be saved as a Mat_entries x nev dimensional array
    // of complex numbers
    // determine dataspace
	char dset_name[80];
	sprintf(dset_name, "ts_%03d",ts);
    /*
     * Create property list for collective dataset write (same for all datasets).
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    const int rank = 2;
    hsize_t dims[rank];
    dims[0] = MAT_ENTRIES;
    dims[1] = nev;
	// write eigenvectors
    //convert eigensystem to 3d-array
    hid_t dataspace = H5Screate_simple(rank,dims,NULL);
    hid_t evecs_ts = H5Dcreate(group_id_evecs, dset_name, h5_comp, dataspace,
							   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // write
    RM_MatrixXcd write_esys = esys;
    err = H5Dwrite(evecs_ts, h5_comp, H5S_ALL, H5S_ALL, plist_id, write_esys.data());
  	H5Dclose(evecs_ts);
  	H5Sclose(dataspace);

    // eigenvalues
    // read
    read_eigenvalues_bin("eigenvalues", conf, ts, nev, eval);
    rank = 1;
    hsize_t dim[rank];
    dim[0] = nev;
    hid_t dataspace = H5Screate_simple(rank,dims,NULL);
    hid_t evals_ts = H5Dcreate(group_id_evals, dset_name, H5T_NATIVE_DOUBLE, dataspace,
							   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // write
    /*
     * Create property list for collective dataset write.
     */
    err = H5Dwrite(evals_ts, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, evals.data());
  	H5Dclose(evals_ts_ts);
  	H5Sclose(dataspace);

    // phases
    // read
    read_eigenvalues_bin("phases", conf, ts, nev, phase);
    // write, sane rank as for eigenvalues
    hid_t dataspace = H5Screate_simple(rank,dims,NULL);
    hid_t phase_ts = H5Dcreate(group_id_phase, dset_name, H5T_NATIVE_DOUBLE, dataspace,
							   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // write
    /*
     * Create property list for collective dataset write.
     */
    err = H5Dwrite(phase_ts, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, evals.data());
  	H5Dclose(phase_ts);
  	H5Sclose(dataspace);
  }
  /*
   * Close/release resources.
   */
  H5Gclose(group_id_info);
  H5Gclose(group_id_data);
  H5Gclose(group_id_evecs);
  H5Gclose(group_id_evals);
  H5Gclose(group_id_phase);
  H5Pclose(plist_id);
  H5Fclose(file_id);
  MPI_Finalize();
  return 0;
}

// source functions for parallel IO. Doc in read_write_mpi.h

#include "read_write_mpi.h"

  void MpiIO::setup(const int mpi_rank, const MPI_Info info){

    // if not master process, do nothing
    if (mpi_rank == 0){

      // property list for file creation parallel I/O access
      hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

      // create file with specified list
      // convert filename
      const char* _fn = filename.c_str();
      file_id = H5Fcreate(_fn, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
      H5Pclose(plist_id);
      // create the necessary groups
      // info for metadata
      gr_info = H5Gcreate(file_id, "/params", H5P_DEFAULT, H5P_DEFAULT,
                                H5P_DEFAULT);

      // group holding datasets
      gr_data = H5Gcreate(file_id, "/data", H5P_DEFAULT, H5P_DEFAULT,
                                H5P_DEFAULT);
      //inside groups create groups for datasets
      gr_evecs = H5Gcreate(file_id, "data/evecs", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      gr_evals = H5Gcreate(file_id, "/data/evals", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      gr_phase = H5Gcreate(file_id, "/data/phase", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
       
       
    }
    else {
      std::cout << "MpiIO: mpi rank is not 0, no file setup." << std::endl;
    }
  }


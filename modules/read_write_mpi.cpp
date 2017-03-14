// source functions for parallel IO. Doc in read_write_mpi.h

#include "read_write_mpi.h"
#include "par_io.h"

static IO* const pars=IO::getInstance();
void MpiIO::setup(const int mpi_rank, const MPI_Info info){

  // if not master process, do nothing
  //if (mpi_rank == 0){

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
    std::cout<< "Set up HDF5 File on process: " << mpi_rank << std::endl; 
  //}
  //else {
  //  std::cout << "MpiIO: mpi rank is not 0, no file setup." << std::endl;
  //}
}

void MpiIO::write_ds(const size_t ts, Eigen::MatrixXcd& V){
  //Create property list for collective dataset write (same for all datasets).
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  const int rank = 2;
  hsize_t dims[rank];
  dims[0] = pars -> get_int("MAT_ENTRIES");
  dims[1] = pars -> get_int("NEV");
  // Define a new compund type for complex data
  // TODO: Code Doubling, typedef that somehow
  hid_t h5_comp = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<double>));
  herr_t err = H5Tinsert(h5_comp,"real",0,H5T_NATIVE_DOUBLE);
  err = H5Tinsert(h5_comp,"imag", sizeof(double), H5T_NATIVE_DOUBLE);

  //write eigenvectors
  char _dn [6];
  sprintf(_dn,"ts_%03d",ts);
  //convert eigensystem to 3d-array
  hid_t dataspace = H5Screate_simple(rank,dims,NULL);
  hid_t evecs_ts = H5Dcreate(gr_evecs, _dn, h5_comp, dataspace,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // write
  RM_MatrixXcd write_esys = V;
  err = H5Dwrite(evecs_ts, h5_comp, H5S_ALL, H5S_ALL, plist_id, write_esys.data());
  H5Dclose(evecs_ts);
  H5Sclose(dataspace);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  std::cout << "written timeslice " << ts << " from rank:" << mpi_rank << " into file: " << this->file_id <<  std::endl; 
}

void MpiIO::write_ds(const size_t ts, const std::string id, const std::vector<double>& ev){
  // Create property list for collective dataset write (same for all datasets).
  // get group identifier depending on string id
  hid_t gr_id;
  if (id.compare("ev") == 0){gr_id = this -> gr_evals;}
  else if (id.compare("phs") == 0){gr_id = this -> gr_phase;}
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  char _dn [6];
  sprintf(_dn,"ts_%03d",ts);
  const int rank = 1;
  hsize_t dim[rank];
  dim[0] = pars -> get_int("NEV");
  hid_t dataspace = H5Screate_simple(rank,dim,NULL);
  hid_t evals_ts = H5Dcreate(gr_id, _dn, H5T_NATIVE_DOUBLE, dataspace,
						   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // write
  herr_t err = H5Dwrite(evals_ts, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, ev.data());
  H5Dclose(evals_ts);
  H5Sclose(dataspace);
}


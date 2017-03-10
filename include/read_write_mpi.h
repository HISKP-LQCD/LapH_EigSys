#ifndef READ_WRITE_MPI_H_
#define READ_WRITE_MPI_H_

#include <array>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <hdf5.h>
#include <mpi.h>
#include <vector>

//#include "hdf5_typedefs.h"

//! Parallel IO for eigensystems
/*! Class to setup, read and write eigensystems to a binary file. Uses the hdf5
 * interface.
 */
 

class MpiIO {

  private:
  //! Filename 
  /*! std::string used to name the file 
   */
  std::string filename;
  //! hdf5 File identifier
  hid_t file_id;
  //! hdf5 group identifier for metadata
  hid_t gr_info;
  //! hdf5 group identifier for Data
  hid_t gr_data;
  //! hdf5 group identifier for datasets
  hid_t gr_evecs;
  hid_t gr_evals;
  hid_t gr_phase;


  protected:
  public:
  //! Default ctor
  MpiIO();
  //! Overloaded ctor taking filename
  MpiIO(const std::string name){
    filename = name;
  }; 
  //! Overloaded ctor taking filename
  MpiIO(const std::string name, const int mpi_rank){
    filename = name;
    setup(mpi_rank, MPI_INFO_NULL);
  }; 

  //! Set up file
  /*! Set up a file for writing datasets to it in parallel. Should only be
   * called ones. To ensure that make this only accessible from MPI process 0.
   * Arguments:
   * ----------
   *
   */
  void setup(const int mpi_rank, const MPI_Info info);

  inline void finalize(const int mpi_rank){
    if (mpi_rank == 0){
      // close all open groups
      H5Gclose(gr_info);
      H5Gclose(gr_evecs);
      H5Gclose(gr_evals);
      H5Gclose(gr_phase);
      H5Gclose(gr_data);
      // close open file
      H5Fclose(file_id);
    }
    else {
      std::cout << "MpiIO: mpi rank is not 0, no file finalize." << std::endl;
    }
  };
  //! Write a record to a hdf5 file in parallel
  /*! The record is set up in the function, the datatype is specified on input
  */
  //void write_rec(const to_write, const rec_id);
  //! read a record from a hdf5 file in parallel
  //void read_rec()


};
#endif // READ_WRITE_MPI_H_

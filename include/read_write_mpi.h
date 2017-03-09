#ifndef READ_WRITE_MPI_H_
#define READ_WRITE_MPI_H_
#include <array>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <hdf5>
#include <vector>
#include <Eigen/Eigen>

//#include "sourceshape.h"
#include "par_io.h"
#include "structs.h"
#include "timeslice.h"
#include "variables.h"
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
  

  protected:
  public:
  //! Default ctor
  MpiIO();
  //! Overloaded ctor taking filename
  MpiIO(std::string name){
    filename = name;
  }  
  //! Default dtor
  ~MpiIO();

  //! write a record to a hdf5 file in parallel
  /*!
  */
  void write_rec();
  //! read a record from a hdf5 file in parallel
  void read_rec()


};
//**************************** Output to files ********************************/
///***************************Input from files**********************************/
//
////Read in Eigenvectors from one Timeslice in binary format to V
//void read_evectors_bin_ts(const char* prefix, const int config_i, const int t,
//    const int nb_ev, Eigen::MatrixXcd& V);
//
////Read in eigenvalues from ascii file to std::array
//void read_eigenvalues_ascii( const char* prefix,const int config_i, const int t,
//    const int nb_ev, std::vector<double>& ev);
//
////Read in eigenvalues from binary to std::vector
//void read_eigenvalues_bin( const char* prefix, const int config_i, const int t,
//    const int nb_ev, std::vector<double>& ev);
//
////Reads in Array of gauge-trafo matrices from binary file to Array of
////Eigen::3cd matrices
//void read_gauge_matrices (const char* prefix, Eigen::Matrix3cd* G);
//
////Reads in Sourceshape from binary
////void read_sourceshape_bin(const char* filename, std::vector<shp> sourceshape);
///***************************Output to files***********************************/
////Write Eigenvectors from one timeslice in binary format to file
//void write_eig_sys_bin(const char* prefix, const int config_i, const int t, const int nb_ev, Eigen::MatrixXcd& V);
////Write Results for source shape to ASCII-file
////void write_sourceshape_ascii(const char* prefix, const int config,
////    const int tslice, const int nb_ev, const std::vector<std::pair<double,double> >& results);
//
////Write Results for source shape to binary file
//void write_sourceshape_bin(const char* prefix, const int config,
//    const int tslice, const int nb_ev, const std::vector<std::pair<double,double> >& results);
//
////Write eigenvalues from std::vector to binary
//void write_eigenvalues_bin( const char* prefix, const int config_i, const int t,
//    const int nb_ev, std::vector<double>& ev);
////Write gauge trafo matrices to binary file
//void write_gauge_matrices(const char* prefix, Eigen::Matrix3cd* G);
///*
////Write Results for source shape to binary file
//void write_sourceshape_bin();
//*/
////write gauge link matrices of one timeslice to binary file
//void write_link_matrices_ts(const char* prefix);
#endif // READ_WRITE_MPI_H_

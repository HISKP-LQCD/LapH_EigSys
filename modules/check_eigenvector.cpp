//calculate action of unboosted laplace on k-th e_vector of Timeslice t

#include <array>
#include <cstdlib>
#include <complex>
#include <fstream>
#include <iostream>
#include <Eigen/Eigen>
#include "read_write.h"

//Function declarations:
void check_e_vector( const Eigen::MatrixXcd& V, const std::array<double, 120>& values,
    const int nev );
/*void read_eigenvectors_from_binary_ts(const int config_i, const int t,
    Eigen::MatrixXcd& V);
void read_eigenvalues_from_ascii(const int config_i, const int t,
    std::array<double, 120>& ev); */

/*****************************************************************************/
//Main
int main(int argc, char* argv[]) {
  const int dim_row = V3 * 3;
  int config = atoi(argv[1]);
  int timeslice = atoi(argv[2]);
  int eigenvalue = atoi(argv[3]);

  //memory for timeslice 
  Eigen::MatrixXcd V;
  V = Eigen::MatrixXcd::Zero(dim_row, 120);

  //memory for eigenvalues
  std::array< double, 120 > e_values;

  //get stuff to objects
  read_evectors_bin_ts("eigenvectors", config, timeslice, V);
  read_eigenvalues_from_binary("eigenvectors", config, timeslice, e_values);

  //check correctness
  check_e_vector(V, e_values, eigenvalue);

  return 0;
}
//End Main
/*****************************************************************************/

//Function definitions
//Checks whether Eigenvector nev was calculated correctly
void check_e_vector( const Eigen::MatrixXcd& V,
    const std::array<double,120>& values, const int nev ) {
  
  //instantiate vectors
  Eigen::VectorXcd lhs(3 * V3);
  lhs = Eigen::VectorXcd::Zero(3 * V3);
  //Laplace times vector in terms of Eigen
  for ( int k = 0; k < V3; ++k) {
    lhs.segment(3 * k, 3) = -(eigen_timeslice[k][0] 
        * V.block(3, 1, 3 * up_3d[k][0], nev) 
        + (eigen_timeslice[ down_3d[k][0] ][0].adjoint()) 
        * V.block(3, 1, 3 * down_3d[k][0], nev) 
        + eigen_timeslice[k][1] * V.block(3, 1, 3 * up_3d[k][1], nev)
        + (eigen_timeslice[ down_3d[k][1] ][1].adjoint()) 
        * V.block(3, 1, 3 * down_3d[k][1], nev)
        + eigen_timeslice[k][2] * V.block(3, 1, 3 * up_3d[k][2], nev)
        + (eigen_timeslice[ down_3d[k][2] ][2].adjoint()) 
        * V.block(3, 1, 3 * down_3d[k][2], nev) 
        - 6.0 * V.block(3, 1, 3 * k , nev) );
  }

  lhs /= values.at(nev);
  double dp = std::abs(lhs.dot( V.col(nev) )); 
  if( dp <= 0.0001 ) 
    std::cout << "Eigenvector recovered successfully" << std::endl;
  else std::cout << "Eigenvectors do not match" << std::endl;

}

/*
//Reads in Eigenvectors from one Timeslice in binary format to V
void read_eigenvectors_from_binary_ts(const int config_i, const int t,
    Eigen::MatrixXcd& V) {

  const int dim_row = 3 * V3;
  //buffer for read in
  std::complex<double>* eigen_vec = new std::complex<double>[dim_row];
    //setting up file
    char filename[200];
  sprintf(filename, "eigenvectors.%04d.%03d", config_i, t);
    std::ifstream infile(filename, std::ifstream::binary);
  for (int nev = 0; nev < 120; ++nev) {
    infile.read( (char*) eigen_vec, sizeof(*eigen_vec));
    V.col(nev) = Eigen::Map<Eigen::Matrix< std::complex<double>, dim_row, 1> >(eigen_vec);
  }
  //clean up
  infile.close();
  delete[] eigen_vec;

}

//Reads in eigenvalues from ascii file to std::array
void read_eigenvalues_from_ascii(const int config_i, const int t,
    std::array<double,120>& ev) {
  
  //buffer for read in
  char* value;
  //setting up file
  char filename[200];
  sprintf(filename, "eigenvalues.%04d.%03d", config_i, t);
  std::ifstream infile(filename);
  for (int nev = 0; nev < 120; ++nev) {
    infile.read( value, sizeof(double));
    ev.at(nev) = *value;
  }
  infile.close();
} 
*/
